#!/usr/bin/env python3
"""
benchmarks/posada_crandall_2001/run_benchmark.py
=================================================
Self-contained, fully reproducible replication of the Posada & Crandall (2001, PNAS)
simulation benchmark evaluating RhizAeon across the 32 standard conditions:
  - Simulations I (Power Grid): 4 theta levels x 5 rho levels = 20 conditions (alpha = None)
  - Simulations II (False Positive / Heterotachy Grid): 4 theta levels x 3 alpha levels = 12 conditions (rho = 0.0)

Evaluates:
  1. Alignment-level detection (Power on rho > 0, False Positive Rate on rho = 0)
  2. Informative SNP Jaccard Index (J_SNP)
  3. Directed Triad Fidelity with phylogenetic partial credit (S_triad)
  4. 1-to-1 Bipartite Hungarian matching
  5. Reticulation Precision (PPV), Recall (TPR), and Composite Fidelity Index (R_F1)
  6. 4D evaluation profile: S_R = (TPR, PPV, S_triad, J_SNP)
"""

import os
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
import json
import argparse
import tempfile
import subprocess
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Any, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.optimize import linear_sum_assignment

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RHIZAEON_BIN = str(PROJECT_ROOT / "target" / "release" / "rhizaeon")
MASTER_SEED = 20011120  # Posada & Crandall PNAS publication date

NT_CHARS = ["A", "C", "G", "T"]
# HKY85 equilibrium frequencies: A=0.40, C=0.20, G=0.30, T=0.10 (Posada & Crandall 2001 Table 4)
PI = np.array([0.40, 0.20, 0.30, 0.10], dtype=np.float64)
KAPPA = 2.0


def build_hky_rate_matrix(pi: np.ndarray = PI, kappa: float = KAPPA) -> np.ndarray:
    """Constructs normalized continuous-time HKY generator Q matrix."""
    Q = np.zeros((4, 4), dtype=np.float64)
    # State mapping: 0: A, 1: C, 2: G, 3: T
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            is_transition = (i in (0, 2) and j in (0, 2)) or (i in (1, 3) and j in (1, 3))
            mult = kappa if is_transition else 1.0
            Q[i, j] = mult * pi[j]
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])
    average_rate = -np.sum(pi * np.diag(Q))
    Q /= average_rate
    return Q


GLOBAL_HKY_Q = build_hky_rate_matrix()


@dataclass
class PosadaScenario:
    scenario_id: str
    category: str       # 'power' or 'false_positive'
    theta: float        # Population mutation parameter (10, 50, 100, 200)
    rho: float          # Population recombination parameter (0, 1, 4, 16, 64)
    alpha: Optional[float] # Gamma shape parameter (None for uniform, 2.0, 0.5, 0.05)
    n_taxa: int = 10
    length_nt: int = 1000
    n_reps: int = 100


def build_scenarios(n_reps: int = 100) -> List[PosadaScenario]:
    scenarios = []
    # Simulations I: Power Grid (alpha = None)
    for theta in [10.0, 50.0, 100.0, 200.0]:
        for rho in [0.0, 1.0, 4.0, 16.0, 64.0]:
            sid = f"power_theta{int(theta)}_rho{int(rho)}"
            scenarios.append(
                PosadaScenario(
                    scenario_id=sid,
                    category="power",
                    theta=theta,
                    rho=rho,
                    alpha=None,
                    n_taxa=10,
                    length_nt=1000,
                    n_reps=n_reps
                )
            )
    # Simulations II: False Positive Grid (rho = 0.0, rate heterogeneity)
    for theta in [10.0, 50.0, 100.0, 200.0]:
        for alpha in [2.0, 0.5, 0.05]:
            alpha_str = f"{alpha}".replace(".", "p")
            sid = f"fp_theta{int(theta)}_alpha{alpha_str}"
            scenarios.append(
                PosadaScenario(
                    scenario_id=sid,
                    category="false_positive",
                    theta=theta,
                    rho=0.0,
                    alpha=alpha,
                    n_taxa=10,
                    length_nt=1000,
                    n_reps=n_reps
                )
            )
    return scenarios


def simulate_pure_coalescent_tree(n: int = 10, seed: Optional[int] = None):
    if seed is not None:
        np.random.seed(seed)
    active = list(range(n))
    node_time = {i: 0.0 for i in range(n)}
    children = {}
    branch_len = {}
    next_node = n
    t = 0.0
    while len(active) > 1:
        k = len(active)
        rate = k * (k - 1) / 2.0
        dt = np.random.exponential(1.0 / rate)
        t += dt
        pair_idx = np.random.choice(len(active), size=2, replace=False)
        u, v = active[pair_idx[0]], active[pair_idx[1]]
        node_time[next_node] = t
        children[next_node] = [u, v]
        branch_len[u] = t - node_time[u]
        branch_len[v] = t - node_time[v]
        active.remove(u)
        active.remove(v)
        active.append(next_node)
        next_node += 1
    root = active[0]
    return root, children, branch_len


def simulate_hudson_arg(n: int = 10, L: int = 1000, rho: float = 4.0, seed: Optional[int] = None):
    if rho <= 0.0:
        root, ch, bl = simulate_pure_coalescent_tree(n=n, seed=seed)
        return [(0, L, root, ch, bl)]
    if seed is not None:
        np.random.seed(seed)
    lineages = {i: np.ones(L, dtype=bool) for i in range(n)}
    node_time = {i: 0.0 for i in range(n)}
    active_node = {i: i for i in range(n)}
    site_parents = {s: {} for s in range(L)}
    next_node = n
    site_k = np.full(L, n, dtype=np.int32)
    t = 0.0
    max_steps = 50000
    step = 0
    while np.any(site_k > 1) and step < max_steps:
        step += 1
        g_vals = {}
        for lid, mask in lineages.items():
            true_idx = np.where(mask)[0]
            g_vals[lid] = int(true_idx[-1] - true_idx[0]) if len(true_idx) > 1 else 0
        G = sum(g_vals.values())
        k = len(lineages)
        if k <= 1 and G == 0:
            break
        lambda_coal = k * (k - 1) / 2.0
        lambda_rec = (rho * G) / (2.0 * L) if (L > 0 and rho > 0 and G > 0) else 0.0
        tot = lambda_coal + lambda_rec
        dt = np.random.exponential(1.0 / tot)
        t += dt
        if np.random.rand() < (lambda_coal / tot):
            lids = list(lineages.keys())
            idx = np.random.choice(len(lids), size=2, replace=False)
            u, v = lids[idx[0]], lids[idx[1]]
            m_u, m_v = lineages[u], lineages[v]
            overlap = m_u & m_v
            new_node = next_node
            next_node += 1
            node_time[new_node] = t
            node_u, node_v = active_node[u], active_node[v]
            for s in np.where(overlap)[0]:
                site_parents[s][node_u] = (new_node, t - node_time[node_u])
                site_parents[s][node_v] = (new_node, t - node_time[node_v])
            site_k[overlap] -= 1
            for s in np.where(m_u & ~overlap)[0]:
                site_parents[s][node_u] = (new_node, t - node_time[node_u])
            for s in np.where(m_v & ~overlap)[0]:
                site_parents[s][node_v] = (new_node, t - node_time[node_v])
            del lineages[u]
            del lineages[v]
            del active_node[u]
            del active_node[v]
            union_active = (m_u | m_v) & (site_k > 1)
            if np.any(union_active):
                lineages[new_node] = union_active
                active_node[new_node] = new_node
        else:
            eligible = [lid for lid in lineages if g_vals[lid] > 0]
            w = [g_vals[lid] for lid in eligible]
            lid = eligible[np.random.choice(len(eligible), p=np.array(w) / sum(w))]
            mask = lineages[lid]
            true_idx = np.where(mask)[0]
            b = np.random.randint(true_idx[0], true_idx[-1] + 1)
            m_l = mask.copy()
            m_l[b + 1:] = False
            m_r = mask.copy()
            m_r[:b + 1] = False
            old_node = active_node[lid]
            del lineages[lid]
            del active_node[lid]
            if np.any(m_l):
                nl = next_node
                next_node += 1
                node_time[nl] = t
                lineages[nl] = m_l
                active_node[nl] = nl
                for s in np.where(m_l)[0]:
                    site_parents[s][old_node] = (nl, t - node_time[old_node])
            if np.any(m_r):
                nr = next_node
                next_node += 1
                node_time[nr] = t
                lineages[nr] = m_r
                active_node[nr] = nr
                for s in np.where(m_r)[0]:
                    site_parents[s][old_node] = (nr, t - node_time[old_node])

    def tree_key(s):
        return tuple(sorted((c, p[0], round(p[1], 6)) for c, p in site_parents[s].items()))

    partitions = []
    curr_start = 0
    curr_key = tree_key(0)
    for s in range(1, L):
        k = tree_key(s)
        if k != curr_key:
            partitions.append((curr_start, s))
            curr_start = s
            curr_key = k
    partitions.append((curr_start, L))

    result = []
    for start, end in partitions:
        sp = site_parents[start]
        children = {}
        branch_len = {}
        all_children = set(sp.keys())
        all_parents = set(p[0] for p in sp.values())
        root_candidates = all_parents - all_children
        root = list(root_candidates)[0] if root_candidates else (max(all_parents) if all_parents else 0)
        for c, (p, bl) in sp.items():
            children.setdefault(p, []).append(c)
            branch_len[c] = bl
        result.append((start, end, root, children, branch_len))
    return result


def evolve_sequences(partitions, n=10, L=1000, theta=50.0, alpha=None, seed=None):
    if seed is not None:
        np.random.seed(seed)
    mu_base = theta / (2.0 * L)
    if alpha is not None and alpha > 0:
        site_rates = np.random.gamma(shape=alpha, scale=1.0 / alpha, size=L)
    else:
        site_rates = np.ones(L, dtype=np.float64)
    taxa = [f"Taxon_{i}" for i in range(n)]
    mat = np.zeros((n, L), dtype=np.int8)

    for start, end, root, children, branch_len in partitions:
        seg_len = end - start
        if seg_len <= 0:
            continue
        seg_rates = site_rates[start:end]
        root_seq = np.random.choice(4, size=seg_len, p=PI)
        node_seqs = {root: root_seq}
        stack = [root]
        while stack:
            curr = stack.pop()
            for ch in children.get(curr, []):
                bl_raw = branch_len.get(ch, 0.0)
                curr_seq = node_seqs[curr]
                ch_seq = np.zeros(seg_len, dtype=np.int8)
                if alpha is None:
                    eff_bl = bl_raw * mu_base
                    P = np.tile(PI, (4, 1)) if eff_bl > 50.0 else expm(GLOBAL_HKY_Q * eff_bl)
                    P = np.nan_to_num(P, nan=0.25)
                    P = np.clip(P, 0.0, 1.0)
                    row_sums = np.sum(P, axis=1, keepdims=True)
                    row_sums[row_sums == 0] = 1.0
                    P /= row_sums
                    for s in range(4):
                        mask = (curr_seq == s)
                        if np.any(mask):
                            ch_seq[mask] = np.random.choice(4, size=np.sum(mask), p=P[s])
                else:
                    eff_bls = bl_raw * mu_base * seg_rates
                    unique_rates, inv = np.unique(np.round(eff_bls, 4), return_inverse=True)
                    for idx_u, r in enumerate(unique_rates):
                        site_mask = (inv == idx_u)
                        if not np.any(site_mask):
                            continue
                        P = np.tile(PI, (4, 1)) if r > 50.0 else expm(GLOBAL_HKY_Q * r)
                        P = np.nan_to_num(P, nan=0.25)
                        P = np.clip(P, 0.0, 1.0)
                        row_sums = np.sum(P, axis=1, keepdims=True)
                        row_sums[row_sums == 0] = 1.0
                        P /= row_sums
                        sub_curr = curr_seq[site_mask]
                        sub_ch = np.zeros(len(sub_curr), dtype=np.int8)
                        for s in range(4):
                            m = (sub_curr == s)
                            if np.any(m):
                                sub_ch[m] = np.random.choice(4, size=np.sum(m), p=P[s])
                        ch_seq[site_mask] = sub_ch
                node_seqs[ch] = ch_seq
                stack.append(ch)
        for i in range(n):
            mat[i, start:end] = node_seqs[i]

    true_bps = [p[1] for p in partitions[:-1]]
    metadata = {
        "n": n, "L": L, "theta": theta, "alpha": alpha,
        "num_partitions": len(partitions),
        "true_breakpoints": true_bps,
        "is_recombinant": len(partitions) > 1
    }
    return mat, taxa, metadata


def compute_patristic_distances(n, root, children, branch_len):
    parent_map = {}
    for p, ch_list in children.items():
        for ch in ch_list:
            parent_map[ch] = (p, branch_len.get(ch, 0.0))
    leaf_paths = {}
    for i in range(n):
        path = [(i, 0.0)]
        curr = i
        while curr in parent_map:
            p, bl = parent_map[curr]
            path.append((p, bl))
            curr = p
        leaf_paths[i] = path
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        p_i = [node for node, bl in leaf_paths[i]]
        bl_i = [bl for node, bl in leaf_paths[i]]
        cum_i = np.cumsum(bl_i)
        node_to_dist_i = {node: dist for node, dist in zip(p_i, cum_i)}
        for j in range(i + 1, n):
            p_j = [node for node, bl in leaf_paths[j]]
            bl_j = [bl for node, bl in leaf_paths[j]]
            cum_j = np.cumsum(bl_j)
            mrca = None
            dist_j = 0.0
            for node, dist in zip(p_j, cum_j):
                if node in node_to_dist_i:
                    mrca = node
                    dist_j = dist
                    break
            dist_i = node_to_dist_i[mrca]
            d_tot = dist_i + dist_j
            D[i, j] = d_tot
            D[j, i] = d_tot
    return D


def extract_ground_truth_events(partitions, n, L, mat, taxa):
    if len(partitions) <= 1:
        return []
    sister_map = np.zeros((n, L), dtype=np.int32)
    for start, end, root, children, branch_len in partitions:
        D = compute_patristic_distances(n, root, children, branch_len)
        for i in range(n):
            d_i = D[i].copy()
            d_i[i] = np.inf
            sister_map[i, start:end] = np.argmin(d_i)

    events = []
    for i in range(n):
        sisters, counts = np.unique(sister_map[i], return_counts=True)
        home = sisters[np.argmax(counts)]
        curr_donor = None
        start_u = None
        for u in range(L):
            s = sister_map[i, u]
            if s != home:
                if curr_donor is None:
                    curr_donor = s
                    start_u = u
                elif s != curr_donor:
                    events.append((i, home, curr_donor, start_u, u))
                    curr_donor = s
                    start_u = u
            else:
                if curr_donor is not None:
                    events.append((i, home, curr_donor, start_u, u))
                    curr_donor = None
                    start_u = None
        if curr_donor is not None:
            events.append((i, home, curr_donor, start_u, L))

    filtered = []
    for c, h, d, u1, u2 in events:
        snps = np.where(mat[h, u1:u2] != mat[d, u1:u2])[0]
        if len(snps) >= 1:
            filtered.append({
                'child': taxa[c], 'home': taxa[h], 'donor': taxa[d],
                'u1': u1, 'u2': u2, 'snps_in_tract': len(snps)
            })
    return filtered


def compute_reticulation_metrics(true_events, pred_events, mat, taxa_to_idx, D_tree, L=1000):
    K = len(true_events)
    K_hat = len(pred_events)

    if K == 0:
        if K_hat == 0:
            return {'detected': False, 'K': 0, 'K_hat': 0, 'rf1': 1.0, 'ppv': 1.0, 'tpr': 1.0, 's_triad': 1.0, 'j_snp': 1.0}
        else:
            return {'detected': True, 'K': 0, 'K_hat': K_hat, 'rf1': 0.0, 'ppv': 0.0, 'tpr': 0.0, 's_triad': 0.0, 'j_snp': 0.0}

    if K_hat == 0:
        return {'detected': False, 'K': K, 'K_hat': 0, 'rf1': 0.0, 'ppv': 0.0, 'tpr': 0.0, 's_triad': 0.0, 'j_snp': 0.0}

    sigma_clade = max(1e-4, np.mean(D_tree) / 2.0)

    def w(name_a, name_b):
        if name_a not in taxa_to_idx or name_b not in taxa_to_idx:
            return 1.0 if name_a == name_b else 0.0
        i_a = taxa_to_idx[name_a]
        i_b = taxa_to_idx[name_b]
        return float(np.exp(-D_tree[i_a, i_b] / sigma_clade))

    affinity_matrix = np.zeros((K, K_hat), dtype=np.float64)
    triad_matrix = np.zeros((K, K_hat), dtype=np.float64)
    jsnp_matrix = np.zeros((K, K_hat), dtype=np.float64)
    cost_matrix = np.ones((K, K_hat), dtype=np.float64)

    for i, te in enumerate(true_events):
        h_idx = taxa_to_idx[te['home']]
        d_idx = taxa_to_idx[te['donor']]
        all_snps = set(np.where(mat[h_idx] != mat[d_idx])[0])
        s_true = set(s for s in all_snps if te['u1'] <= s < te['u2'])

        for j, pe in enumerate(pred_events):
            s_inferred = set(s for s in all_snps if pe['u1'] <= s < pe['u2'])
            union = s_true | s_inferred
            inter = s_true & s_inferred
            j_snp = len(inter) / len(union) if union else (1.0 if not s_true and not s_inferred else 0.0)

            w_c = w(te['child'], pe['candidate_name'])
            w_h = w(te['home'], pe['home_name'])
            w_d = w(te['donor'], pe['donor_name'])
            s_triad = w_c * ((w_h + w_d) / 2.0)

            aff = s_triad * j_snp
            affinity_matrix[i, j] = aff
            triad_matrix[i, j] = s_triad
            jsnp_matrix[i, j] = j_snp
            cost_matrix[i, j] = 1.0 - aff

    r_ind, c_ind = linear_sum_assignment(cost_matrix)
    matched_affs = [affinity_matrix[r, c] for r, c in zip(r_ind, c_ind)]
    matched_triads = [triad_matrix[r, c] for r, c in zip(r_ind, c_ind)]
    matched_jsnps = [jsnp_matrix[r, c] for r, c in zip(r_ind, c_ind)]

    tot_aff = sum(matched_affs)
    ppv = tot_aff / K_hat
    tpr = tot_aff / K
    rf1 = (2.0 * ppv * tpr) / (ppv + tpr) if (ppv + tpr) > 0 else 0.0

    mean_triad = float(np.mean(matched_triads)) if matched_triads else 0.0
    mean_jsnp = float(np.mean(matched_jsnps)) if matched_jsnps else 0.0

    return {
        'detected': True, 'K': K, 'K_hat': K_hat, 'rf1': rf1, 'ppv': ppv, 'tpr': tpr,
        's_triad': mean_triad, 'j_snp': mean_jsnp
    }


def evaluate_single_replicate(scenario: PosadaScenario, sc_idx: int, rep_idx: int) -> Dict[str, Any]:
    seed = MASTER_SEED + sc_idx * 10000 + rep_idx
    # 1. Simulate ARG & sequences
    parts = simulate_hudson_arg(n=scenario.n_taxa, L=scenario.length_nt, rho=scenario.rho, seed=seed)
    evo_seed = seed + 1000003
    mat, taxa, meta = evolve_sequences(
        partitions=parts,
        n=scenario.n_taxa,
        L=scenario.length_nt,
        theta=scenario.theta,
        alpha=scenario.alpha,
        seed=evo_seed
    )
    taxa_to_idx = {t: i for i, t in enumerate(taxa)}

    # Background tree for patristic distance
    root, ch, bl = parts[0][2], parts[0][3], parts[0][4]
    D_tree = compute_patristic_distances(scenario.n_taxa, root, ch, bl)

    # Ground truth events
    true_events = extract_ground_truth_events(parts, scenario.n_taxa, scenario.length_nt, mat, taxa) if scenario.rho > 0 else []

    # 2. Run RhizAeon
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f_fa:
        for i, t in enumerate(taxa):
            seq = "".join(NT_CHARS[c] for c in mat[i])
            f_fa.write(f">{t}\n{seq}\n")
        fa_path = f_fa.name

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f_json:
        json_path = f_json.name

    t0 = time.perf_counter()
    res_code = subprocess.run([RHIZAEON_BIN, "-i", fa_path, "-o", json_path], capture_output=True, text=True)
    t_rhiz_ms = (time.perf_counter() - t0) * 1000.0

    pred_events = []
    if os.path.exists(json_path) and os.path.getsize(json_path) > 0:
        try:
            with open(json_path, "r") as f:
                d = json.load(f)
                pred_events = d.get("events", [])
        except Exception:
            pass

    # Clean up temp files
    try:
        os.remove(fa_path)
        os.remove(json_path)
    except Exception:
        pass

    # 3. Compute metrics
    metrics = compute_reticulation_metrics(true_events, pred_events, mat, taxa_to_idx, D_tree, L=scenario.length_nt)
    detected = len(pred_events) > 0

    return {
        "scenario_id": scenario.scenario_id,
        "category": scenario.category,
        "theta": scenario.theta,
        "rho": scenario.rho,
        "alpha": scenario.alpha if scenario.alpha is not None else -1.0,
        "rep_idx": rep_idx,
        "seed": seed,
        "num_partitions": len(parts),
        "detected": detected,
        "K_true": len(true_events),
        "K_hat": len(pred_events),
        "rf1": metrics["rf1"],
        "ppv": metrics["ppv"],
        "tpr": metrics["tpr"],
        "s_triad": metrics["s_triad"],
        "j_snp": metrics["j_snp"],
        "runtime_ms": t_rhiz_ms
    }


def main():
    parser = argparse.ArgumentParser(description="Posada & Crandall (2001) Replication for RhizAeon")
    parser.add_argument("--reps", type=int, default=100, help="Number of replicates per scenario [default: 100]")
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1), help="Number of parallel workers")
    parser.add_argument("--scenarios", type=str, default="all", help="Subset of scenarios ('all', 'power', 'fp')")
    args = parser.parse_args()

    scenarios = build_scenarios(n_reps=args.reps)
    if args.scenarios == "power":
        scenarios = [s for s in scenarios if s.category == "power"]
    elif args.scenarios == "fp":
        scenarios = [s for s in scenarios if s.category == "false_positive"]

    print("=" * 80)
    print("  POSADA & CRANDALL (2001, PNAS) REPLICATION BENCHMARK")
    print(f"  Scenarios: {len(scenarios)} conditions | Replicates: {args.reps} | Total trials: {len(scenarios) * args.reps}")
    print(f"  Workers: {args.workers} | Engine: {RHIZAEON_BIN}")
    print("=" * 80)

    tasks = []
    for sc_idx, sc in enumerate(scenarios):
        for rep in range(args.reps):
            tasks.append((sc, sc_idx, rep))

    t_start = time.perf_counter()
    raw_results = []
    completed = 0
    total = len(tasks)

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(evaluate_single_replicate, sc, sc_idx, rep): (sc.scenario_id, rep) for sc, sc_idx, rep in tasks}
        for future in as_completed(futures):
            res = future.result()
            raw_results.append(res)
            completed += 1
            if completed % 100 == 0 or completed == total:
                elapsed = time.perf_counter() - t_start
                rate = completed / elapsed
                print(f"  Progress: {completed:4d}/{total:4d} ({completed/total*100:5.1f}%) | {rate:5.1f} trials/s | Elapsed: {elapsed:5.1f}s")

    df_raw = pd.DataFrame(raw_results)
    raw_csv = SCRIPT_DIR / "posada_raw_results.csv"
    df_raw.to_csv(raw_csv, index=False)
    print(f"\n[DONE] Saved raw results ({len(df_raw)} trials) to {raw_csv}")

    # Summary table aggregation
    summary_rows = []
    for sc in scenarios:
        sub = df_raw[df_raw["scenario_id"] == sc.scenario_id]
        n_reps = len(sub)
        det_rate = sub["detected"].mean() * 100.0
        summary_rows.append({
            "scenario_id": sc.scenario_id,
            "category": sc.category,
            "theta": sc.theta,
            "rho": sc.rho,
            "alpha": sc.alpha if sc.alpha is not None else np.nan,
            "n_reps": n_reps,
            "mean_partitions": sub["num_partitions"].mean(),
            "mean_K_true": sub["K_true"].mean(),
            "mean_K_hat": sub["K_hat"].mean(),
            "detection_pct": round(det_rate, 2),
            "mean_rf1": round(sub["rf1"].mean(), 4),
            "mean_ppv": round(sub["ppv"].mean(), 4),
            "mean_tpr": round(sub["tpr"].mean(), 4),
            "mean_s_triad": round(sub["s_triad"].mean(), 4),
            "mean_j_snp": round(sub["j_snp"].mean(), 4),
            "mean_runtime_ms": round(sub["runtime_ms"].mean(), 2)
        })

    df_summary = pd.DataFrame(summary_rows)
    summary_csv = SCRIPT_DIR / "posada_summary.csv"
    df_summary.to_csv(summary_csv, index=False)
    print(f"[DONE] Saved summary table ({len(df_summary)} scenarios) to {summary_csv}")
    print("=" * 80)


if __name__ == "__main__":
    main()
