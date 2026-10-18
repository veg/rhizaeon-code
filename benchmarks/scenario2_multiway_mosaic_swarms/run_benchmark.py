#!/usr/bin/env python3
"""
Scenario 1.2: Multi-Way Mosaic Swarms Benchmark Harness.
Compares RhizAeon vs 3SEQ on complex multi-parental recombination topologies:
  1. Dual-Cassette: P1 background with 2 distinct cassettes from P2 and P3.
  2. Nested Russian Doll: P1 background with outer P2 cassette enclosing inner P3 cassette.
  3. Triple Crossover Cascade: Sequential serial transitions P1 -> P2 -> P3 across the chromosome.

Evaluated across:
  - Divergences d in {0.02, 0.05, 0.10}
  - Descendant recombinant swarm sizes K in {1, 5, 10}
  - 3 Background clades of 15 taxa each (45 background taxa)
  - Sequence length L = 6,000 nt
  - 20 Replicates per condition (540 trials per tool, 1080 evaluations total)
"""

import os
import sys
import json
import time
import tempfile
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd
from scipy.linalg import expm

NT_CHARS = ["A", "C", "G", "T"]
PI = np.array([0.30, 0.20, 0.20, 0.30], dtype=np.float64)
KAPPA = 2.5

def build_hky_rate_matrix(pi: np.ndarray = PI, kappa: float = KAPPA) -> np.ndarray:
    Q = np.zeros((4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            is_ts = (i in (0, 2) and j in (0, 2)) or (i in (1, 3) and j in (1, 3))
            mult = kappa if is_ts else 1.0
            Q[i, j] = mult * pi[j]
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])
    average_rate = -np.sum(pi * np.diag(Q))
    Q /= average_rate
    return Q

GLOBAL_HKY_Q = build_hky_rate_matrix()

def evolve_sequence(seq: np.ndarray, branch_length: float, rng: np.random.Generator) -> np.ndarray:
    if branch_length <= 0.0:
        return seq.copy()
    P_mat = expm(GLOBAL_HKY_Q * branch_length)
    new_seq = np.empty_like(seq)
    for state in range(4):
        mask = (seq == state)
        count = int(np.sum(mask))
        if count > 0:
            probs = np.maximum(P_mat[state, :], 0.0)
            probs /= np.sum(probs)
            new_seq[mask] = rng.choice(4, size=count, p=probs)
    return new_seq

def generate_multiway_alignment(
    topology: str = "dual_cassette",
    d_parent: float = 0.05,
    K: int = 5,
    n_p1: int = 15,
    n_p2: int = 15,
    n_p3: int = 15,
    L: int = 6000,
    drift: float = 0.002,
    seed: int = 42,
):
    rng = np.random.Generator(np.random.PCG64(seed))
    
    anc_root = rng.choice(4, size=L, p=PI)
    p1_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p2_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p3_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    
    taxa = []
    seqs = []
    
    # Background clades
    for i in range(n_p1):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, drift, rng))
    for i in range(n_p2):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, drift, rng))
    for i in range(n_p3):
        taxa.append(f"P3_{i}")
        seqs.append(evolve_sequence(p3_anc, drift, rng))
        
    rec_anc = p1_anc.copy()
    
    if topology == "dual_cassette":
        # P1 background with Cassette 1 [1200..2400] from P2 and Cassette 2 [3600..4800] from P3
        rec_anc[1200 - 1 : 2400] = p2_anc[1200 - 1 : 2400]
        rec_anc[3600 - 1 : 4800] = p3_anc[3600 - 1 : 4800]
        true_tracts = [
            {"home_clade": "P1", "donor_clade": "P2", "u1": 1200, "u2": 2400},
            {"home_clade": "P1", "donor_clade": "P3", "u1": 3600, "u2": 4800},
        ]
        k_expected_events = 2
    elif topology == "nested_russian_doll":
        # P1 background with Outer Cassette [1500..4500] from P2,
        # inside which [2500..3500] is replaced by P3.
        rec_anc[1500 - 1 : 4500] = p2_anc[1500 - 1 : 4500]
        rec_anc[2500 - 1 : 3500] = p3_anc[2500 - 1 : 3500]
        true_tracts = [
            {"home_clade": "P1", "donor_clade": "P2", "u1": 1500, "u2": 2500},
            {"home_clade": "P1", "donor_clade": "P3", "u1": 2500, "u2": 3500},
            {"home_clade": "P1", "donor_clade": "P2", "u1": 3500, "u2": 4500},
        ]
        k_expected_events = 2 # 2 distinct donor lineages (P2 and P3)
    elif topology == "triple_crossover":
        # P1: 1..2000, P2: 2001..4000, P3: 4001..6000
        rec_anc[2000 : 4000] = p2_anc[2000 : 4000]
        rec_anc[4000 : 6000] = p3_anc[4000 : 6000]
        true_tracts = [
            {"home_clade": "P1", "donor_clade": "P2", "u1": 2001, "u2": 4000},
            {"home_clade": "P2", "donor_clade": "P3", "u1": 4001, "u2": 6000},
        ]
        k_expected_events = 2 # 2 crossover transitions
    else:
        raise ValueError(f"Unknown topology: {topology}")
        
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat, true_tracts, k_expected_events

def compute_jaccard(true_u1: int, true_u2: int, pred_u1: int, pred_u2: int) -> float:
    start = max(true_u1, pred_u1)
    end = min(true_u2, pred_u2)
    inter = max(0, end - start + 1)
    union = (true_u2 - true_u1 + 1) + (pred_u2 - pred_u1 + 1) - inter
    return inter / union if union > 0 else 0.0

def parse_3seq_bps(bp_str: str):
    import re
    if not bp_str:
        return (0, 0)
    if "&" in bp_str:
        parts = bp_str.split("&")
        m1 = re.findall(r"\d+", parts[0])
        m2 = re.findall(r"\d+", parts[1])
        if len(m1) >= 2 and len(m2) >= 2:
            return ((int(m1[0]) + int(m1[1])) // 2, (int(m2[0]) + int(m2[1])) // 2)
        elif len(m1) >= 1 and len(m2) >= 1:
            return (int(m1[0]), int(m2[0]))
    else:
        m = re.findall(r"\d+", bp_str)
        if len(m) >= 2:
            return (int(m[0]), int(m[1]))
    return (0, 0)

def consolidate_3seq_rows(rows, taxa, mat, bp_tol=200):
    N, L = mat.shape
    name_to_idx = {t: i for i, t in enumerate(taxa)}
    dist_mat = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        dist_mat[i] = np.mean(mat[i] != mat, axis=1)
    triu = dist_mat[np.triu_indices(N, k=1)]
    clade_rad = max(0.01, 0.35 * float(np.median(triu))) if len(triu) > 0 else 0.01

    def same_clade(a, b):
        if a == b:
            return True
        ia, ib = name_to_idx.get(a), name_to_idx.get(b)
        if ia is None or ib is None:
            return False
        return dist_mat[ia, ib] <= clade_rad

    sorted_rows = sorted(rows, key=lambda r: (r["log10p"], r["dsp"]))
    accepted = []
    for r in sorted_rows:
        is_mirror = False
        for acc in accepted:
            overlap = min(r["u2"], acc["u2"]) - max(r["u1"], acc["u1"])
            shared_bp = abs(r["u1"] - acc["u1"]) <= bp_tol or abs(r["u2"] - acc["u2"]) <= bp_tol or overlap > 0
            if shared_bp:
                cites_acc_rec = any(same_clade(r["pa"], iso) or same_clade(r["pb"], iso) for iso in acc["isolates"])
                child_is_acc_parent = same_clade(r["child"], acc["pa"]) or same_clade(r["child"], acc["pb"])
                if cites_acc_rec and child_is_acc_parent:
                    is_mirror = True
                    break
        if is_mirror:
            continue

        merged = False
        for acc in accepted:
            same_child_clade = any(same_clade(r["child"], iso) for iso in acc["isolates"])
            same_parents = (same_clade(r["pa"], acc["pa"]) and same_clade(r["pb"], acc["pb"])) or (
                same_clade(r["pa"], acc["pb"]) and same_clade(r["pb"], acc["pa"])
            )
            overlap = min(r["u2"], acc["u2"]) - max(r["u1"], acc["u1"])
            min_len = max(1, min(r["u2"] - r["u1"] + 1, acc["u2"] - acc["u1"] + 1))
            shared_bp = (overlap >= 0.5 * min_len) or (abs(r["u1"] - acc["u1"]) <= bp_tol and abs(r["u2"] - acc["u2"]) <= bp_tol)
            if (same_child_clade or same_parents) and shared_bp:
                acc["isolates"].add(r["child"])
                merged = True
                break
        if not merged:
            accepted.append({
                "pa": r["pa"],
                "pb": r["pb"],
                "child": r["child"],
                "u1": r["u1"],
                "u2": r["u2"],
                "log10p": r["log10p"],
                "isolates": {r["child"]},
            })
    return accepted

def run_single_trial(args):
    topology, d_parent, K, rep, rhiz_bin, three_bin = args
    top_code = {"dual_cassette": 1, "nested_russian_doll": 2, "triple_crossover": 3}[topology]
    seed = 200000 + top_code * 10000 + int(d_parent * 1000) * 100 + K * 10 + rep
    n_p1, n_p2, n_p3 = 15, 15, 15
    L = 6000
    true_recs = set(f"Rec_{k}" for k in range(K))
    bg_taxa = set(f"P1_{i}" for i in range(n_p1)).union(
        f"P2_{j}" for j in range(n_p2)
    ).union(f"P3_{m}" for m in range(n_p3))
    
    taxa, mat, true_tracts, k_expected_events = generate_multiway_alignment(
        topology=topology, d_parent=d_parent, K=K, n_p1=n_p1, n_p2=n_p2, n_p3=n_p3, L=L, seed=seed
    )
    
    pid = os.getpid()
    tmp_fa = tempfile.NamedTemporaryFile(suffix=f"_{pid}_{topology}_{rep}.fasta", delete=False)
    for name, s in zip(taxa, mat):
        s_str = "".join(NT_CHARS[c] for c in s)
        tmp_fa.write(f">{name}\n{s_str}\n".encode())
    tmp_fa.close()
    
    results = []
    
    # ---------------- 1. RhizAeon ----------------
    out_json = tmp_fa.name + ".json"
    t0 = time.perf_counter()
    r_proc = subprocess.run(
        [rhiz_bin, "--input", tmp_fa.name, "--output", out_json],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        timeout=180,
    )
    t_rhiz_ms = (time.perf_counter() - t0) * 1000.0
    
    rhiz_unique_events = 0
    rhiz_raw_alerts = 0
    rhiz_recs_called = set()
    rhiz_bg_called = set()
    rhiz_mean_jaccard = 0.0
    rhiz_detected_tracts_count = 0
    rhiz_donors_detected = set()
    
    if r_proc.returncode == 0 and os.path.exists(out_json):
        try:
            with open(out_json) as f:
                r_res = json.load(f)
            events = r_res.get("events", [])
            rhiz_unique_events = len(events)
            
            jaccard_scores = []
            for ev in events:
                rhiz_raw_alerts += len(ev.get("isolates", []))
                for iso in ev.get("isolates", []):
                    if iso in true_recs:
                        rhiz_recs_called.add(iso)
                    elif iso in bg_taxa:
                        rhiz_bg_called.add(iso)
                
                # Check donor clade
                d_name = ev.get("donor_name", "")
                if d_name.startswith("P2"):
                    rhiz_donors_detected.add("P2")
                elif d_name.startswith("P3"):
                    rhiz_donors_detected.add("P3")
                    
                # Compute best Jaccard against true tracts
                best_j = 0.0
                for tt in true_tracts:
                    j = compute_jaccard(tt["u1"], tt["u2"], ev.get("u1", 0), ev.get("u2", 0))
                    if j > best_j:
                        best_j = j
                if best_j > 0.0:
                    jaccard_scores.append(best_j)
                    if best_j >= 0.50:
                        rhiz_detected_tracts_count += 1
                        
            rhiz_mean_jaccard = float(np.mean(jaccard_scores)) if jaccard_scores else 0.0
        except Exception:
            pass
        finally:
            if os.path.exists(out_json):
                os.remove(out_json)
                
    rhiz_recall = len(rhiz_recs_called) / K
    rhiz_precision = (
        len(rhiz_recs_called) / (len(rhiz_recs_called) + len(rhiz_bg_called))
        if (len(rhiz_recs_called) + len(rhiz_bg_called)) > 0
        else 1.0
    )
    rhiz_bg_fpr = len(rhiz_bg_called) / len(bg_taxa)
    
    results.append({
        "tool": "RhizAeon",
        "topology": topology,
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "expected_events": k_expected_events,
        "unique_events": rhiz_unique_events,
        "raw_alerts": rhiz_raw_alerts,
        "unique_rec_taxa": len(rhiz_recs_called),
        "recall": rhiz_recall,
        "precision": rhiz_precision,
        "bg_fp_count": len(rhiz_bg_called),
        "bg_fpr": rhiz_bg_fpr,
        "mean_jaccard": rhiz_mean_jaccard,
        "donors_resolved": len(rhiz_donors_detected),
        "time_ms": t_rhiz_ms,
    })
    
    # ---------------- 2. 3SEQ (Raw & Consolidated) ----------------
    run_id = f"3seq_{pid}_{top_code}_{int(d_parent*100)}_{K}_{rep}"
    for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped"]:
        p = Path(f"{run_id}{ext}")
        if p.exists():
            p.unlink()
            
    t0 = time.perf_counter()
    try:
        subprocess.run(
            [three_bin, "-full", tmp_fa.name, "-id", run_id],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=180,
        )
    except subprocess.TimeoutExpired:
        pass
    t_3seq_ms = (time.perf_counter() - t0) * 1000.0
    
    rec_csv = f"{run_id}.3s.rec.csv"
    three_triplets = 0
    three_recs_called = set()
    three_bg_called = set()
    three_jaccards = []
    three_donors_detected = set()
    raw_rows = []
    
    if os.path.exists(rec_csv):
        try:
            with open(rec_csv) as f:
                for line in f:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) >= 10 and parts[0] != "P_name":
                        try:
                            p_val = float(parts[9])
                            if p_val <= 0.05:
                                three_triplets += 1
                                child = parts[2]
                                pa = parts[0]
                                pb = parts[1]
                                log10p = float(parts[8]) if parts[8] != "-inf" else -999.0
                                if child in true_recs:
                                    three_recs_called.add(child)
                                    for p_tax in (pa, pb):
                                        if p_tax.startswith("P2"):
                                            three_donors_detected.add("P2")
                                        elif p_tax.startswith("P3"):
                                            three_donors_detected.add("P3")
                                elif child in bg_taxa:
                                    three_bg_called.add(child)
                                    
                                u1_3s, u2_3s = 0, 0
                                if len(parts) >= 13:
                                    bp_str = parts[12].strip()
                                    u1_3s, u2_3s = parse_3seq_bps(bp_str)
                                    if u1_3s > 0 and u2_3s > u1_3s:
                                        best_j = 0.0
                                        for tt in true_tracts:
                                            j = compute_jaccard(tt["u1"], tt["u2"], u1_3s, u2_3s)
                                            if j > best_j:
                                                best_j = j
                                        if best_j > 0.0:
                                            three_jaccards.append(best_j)
                                raw_rows.append({
                                    "pa": pa,
                                    "pb": pb,
                                    "child": child,
                                    "log10p": log10p,
                                    "dsp": p_val,
                                    "u1": u1_3s,
                                    "u2": u2_3s,
                                })
                        except Exception:
                            pass
        except Exception:
            pass
        finally:
            for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped"]:
                p = Path(f"{run_id}{ext}")
                if p.exists():
                    p.unlink()
                    
    three_recall = len(three_recs_called) / K
    three_precision = (
        len(three_recs_called) / (len(three_recs_called) + len(three_bg_called))
        if (len(three_recs_called) + len(three_bg_called)) > 0
        else 1.0
    )
    three_bg_fpr = len(three_bg_called) / len(bg_taxa)
    three_mean_jaccard = float(np.mean(three_jaccards)) if three_jaccards else 0.0
    
    results.append({
        "tool": "3SEQ",
        "topology": topology,
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "expected_events": k_expected_events,
        "unique_events": len(three_recs_called) + len(three_bg_called),
        "raw_alerts": three_triplets,
        "unique_rec_taxa": len(three_recs_called),
        "recall": three_recall,
        "precision": three_precision,
        "bg_fp_count": len(three_bg_called),
        "bg_fpr": three_bg_fpr,
        "mean_jaccard": three_mean_jaccard,
        "donors_resolved": len(three_donors_detected),
        "time_ms": t_3seq_ms,
    })

    # Consolidated 3SEQ ("Synthesis")
    cons_events = consolidate_3seq_rows(raw_rows, taxa, mat) if raw_rows else []
    cons_recs_called = set()
    cons_bg_called = set()
    cons_jaccards = []
    cons_donors_detected = set()
    for cev in cons_events:
        for iso in cev["isolates"]:
            if iso in true_recs:
                cons_recs_called.add(iso)
                for p_tax in (cev["pa"], cev["pb"]):
                    if p_tax.startswith("P2"):
                        cons_donors_detected.add("P2")
                    elif p_tax.startswith("P3"):
                        cons_donors_detected.add("P3")
            elif iso in bg_taxa:
                cons_bg_called.add(iso)
        if cev["u1"] > 0 and cev["u2"] > cev["u1"]:
            best_j = 0.0
            for tt in true_tracts:
                j = compute_jaccard(tt["u1"], tt["u2"], cev["u1"], cev["u2"])
                if j > best_j:
                    best_j = j
            if best_j > 0.0:
                cons_jaccards.append(best_j)

    cons_recall = len(cons_recs_called) / K
    cons_precision = (
        len(cons_recs_called) / (len(cons_recs_called) + len(cons_bg_called))
        if (len(cons_recs_called) + len(cons_bg_called)) > 0
        else 1.0
    )
    cons_bg_fpr = len(cons_bg_called) / len(bg_taxa)
    cons_mean_jaccard = float(np.mean(cons_jaccards)) if cons_jaccards else 0.0

    results.append({
        "tool": "3SEQ (Consolidated)",
        "topology": topology,
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "expected_events": k_expected_events,
        "unique_events": len(cons_events),
        "raw_alerts": len(cons_recs_called) + len(cons_bg_called),
        "unique_rec_taxa": len(cons_recs_called),
        "recall": cons_recall,
        "precision": cons_precision,
        "bg_fp_count": len(cons_bg_called),
        "bg_fpr": cons_bg_fpr,
        "mean_jaccard": cons_mean_jaccard,
        "donors_resolved": len(cons_donors_detected),
        "time_ms": t_3seq_ms,
    })
    
    if os.path.exists(tmp_fa.name):
        os.remove(tmp_fa.name)
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Run Scenario 1.2 Multi-Way Mosaic Swarms Benchmark")
    parser.add_argument("--replicates", type=int, default=20, help="Replicates per condition (default 20)")
    parser.add_argument("--threads", type=int, default=min(8, cpu_count()), help="Multiprocessing workers")
    args = parser.parse_args()
    
    rhiz_bin = "/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon"
    three_bin = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
    
    if not os.path.exists(rhiz_bin):
        print(f"Error: RhizAeon binary not found at {rhiz_bin}")
        sys.exit(1)
    if not os.path.exists(three_bin):
        print(f"Error: 3SEQ binary not found at {three_bin}")
        sys.exit(1)
        
    topologies = ["dual_cassette", "nested_russian_doll", "triple_crossover"]
    divergences = [0.02, 0.05, 0.10]
    ks = [1, 5, 10]
    
    tasks = []
    for top in topologies:
        for d in divergences:
            for k in ks:
                for rep in range(args.replicates):
                    tasks.append((top, d, k, rep, rhiz_bin, three_bin))
                    
    total_trials = len(tasks)
    print("=" * 80)
    print("SCENARIO 1.2: MULTI-WAY MOSAIC SWARMS BENCHMARK")
    print(f"  Topologies ({len(topologies)}): {', '.join(topologies)}")
    print(f"  Divergence levels ({len(divergences)}): {divergences}")
    print(f"  Swarm sizes K ({len(ks)}): {ks}")
    print(f"  Replicates per condition: {args.replicates}")
    print(f"  Total trials: {total_trials} (x2 tools = {total_trials * 2} runs)")
    print(f"  Workers: {args.threads}")
    print("=" * 80)
    
    t_start = time.perf_counter()
    all_results = []
    completed = 0
    
    with Pool(processes=args.threads) as pool:
        for res_list in pool.imap_unordered(run_single_trial, tasks, chunksize=1):
            all_results.extend(res_list)
            completed += 1
            if completed % 25 == 0 or completed == total_trials:
                elapsed = time.perf_counter() - t_start
                rate = completed / elapsed if elapsed > 0 else 0.0
                eta = (total_trials - completed) / rate if rate > 0 else 0.0
                print(f"Progress: {completed}/{total_trials} trials complete ({completed/total_trials*100:.1f}%) | Elapsed: {elapsed:.1f}s | Rate: {rate:.2f} trials/s | ETA: {eta:.1f}s")
                sys.stdout.flush()
                
    total_wall_time = time.perf_counter() - t_start
    print("=" * 80)
    print(f"Benchmark finished in {total_wall_time:.1f}s ({total_wall_time/60.0:.2f} min). Processing summary statistics...")
    
    df_raw = pd.DataFrame(all_results)
    out_dir = Path("/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms")
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = out_dir / "multiway_raw_results.csv"
    df_raw.to_csv(raw_csv, index=False)
    print(f"Saved raw trial results to {raw_csv}")
    
    # Aggregation
    summary = df_raw.groupby(["topology", "d_parent", "K", "tool"]).agg(
        total_runs=("rep", "count"),
        mean_recall=("recall", "mean"),
        mean_precision=("precision", "mean"),
        mean_bg_fpr=("bg_fpr", "mean"),
        mean_unique_events=("unique_events", "mean"),
        mean_raw_alerts=("raw_alerts", "mean"),
        mean_jaccard=("mean_jaccard", "mean"),
        mean_donors_resolved=("donors_resolved", "mean"),
        mean_time_ms=("time_ms", "mean"),
    ).reset_index()
    
    summary_csv = out_dir / "multiway_summary.csv"
    summary.to_csv(summary_csv, index=False)
    print(f"Saved summary table to {summary_csv}")
    print("\nSummary Results Preview:")
    print(summary.to_string())

if __name__ == "__main__":
    main()
