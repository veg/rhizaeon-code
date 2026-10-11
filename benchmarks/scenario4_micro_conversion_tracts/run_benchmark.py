#!/usr/bin/env python3
"""
Scenario 1.4: Bacterial Micro-Conversion & Scale-Space Resolution Suite Benchmark Harness.
Compares RhizAeon vs 3SEQ across micro-conversion tract lengths:
  L_tract in {50, 100, 200, 400, 800} nt (in L = 5,000 nt genome).

Evaluated across:
  - Divergences d_parent in {0.02, 0.05, 0.10}
  - Descendant recombinant swarm sizes K in {2, 5, 10}
  - Background clades: P1 (15 taxa), P2 (15 taxa) -> 30 background taxa
  - Total cohort size N = 30 + K in {32, 35, 40} taxa
  - 20 Replicates per condition (45 conditions * 20 reps = 900 trials per tool, 1800 total evaluations)

Ground Truth:
  K_true = 1 ancestral micro-conversion event
  True tract = [2500 - L_tract//2, 2500 + L_tract//2 - 1]
  True recombinant roster = {Rec_0, ..., Rec_{K-1}}
  True background non-recombinants = {P1_0..P1_14, P2_0..P2_14}
"""

import os
import sys
import re
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

def generate_micro_conversion_alignment(
    tract_len: int = 100,
    d_parent: float = 0.05,
    K: int = 5,
    n_p1: int = 15,
    n_p2: int = 15,
    L: int = 5000,
    drift: float = 0.002,
    seed: int = 42,
):
    rng = np.random.Generator(np.random.PCG64(seed))
    
    # 1. Ancestral root
    anc_root = rng.choice(4, size=L, p=PI)
    
    # 2. Parental ancestral nodes
    p1_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p2_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    
    taxa = []
    seqs = []
    
    # 3. Background P1 clade
    for i in range(n_p1):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, drift, rng))
        
    # 4. Background P2 clade
    for i in range(n_p2):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, drift, rng))
        
    # 5. Micro-conversion coordinates
    u1 = 2500 - tract_len // 2
    u2 = u1 + tract_len - 1
    
    mut_in_tract = int(np.sum(p1_anc[u1 - 1 : u2] != p2_anc[u1 - 1 : u2]))
    total_mut = int(np.sum(p1_anc != p2_anc))
    
    # 6. Ancestral recombinant node A_rec (P1 backbone with P2 micro-conversion)
    rec_anc = p1_anc.copy()
    rec_anc[u1 - 1 : u2] = p2_anc[u1 - 1 : u2]
    
    # 7. Sprout K descendant recombinant isolates with transmission drift
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat, u1, u2, mut_in_tract, total_mut

def compute_jaccard(true_u1: int, true_u2: int, pred_u1: int, pred_u2: int) -> float:
    start = max(true_u1, pred_u1)
    end = min(true_u2, pred_u2)
    inter = max(0, end - start + 1)
    union = (true_u2 - true_u1 + 1) + (pred_u2 - pred_u1 + 1) - inter
    return inter / union if union > 0 else 0.0

def parse_3seq_bps(bp_str: str):
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

def run_single_trial(args):
    tract_len, d_parent, K, rep, rhiz_bin, three_bin = args
    seed = 400000 + tract_len * 1000 + int(d_parent * 1000) * 50 + K * 5 + rep
    n_p1, n_p2 = 15, 15
    L = 5000
    true_recs = set(f"Rec_{k}" for k in range(K))
    bg_taxa = set(f"P1_{i}" for i in range(n_p1)).union(f"P2_{j}" for j in range(n_p2))
    
    taxa, mat, true_u1, true_u2, mut_in_tract, total_mut = generate_micro_conversion_alignment(
        tract_len=tract_len, d_parent=d_parent, K=K, n_p1=n_p1, n_p2=n_p2, L=L, seed=seed
    )
    
    # Write temporary fasta
    pid = os.getpid()
    tmp_fa = tempfile.NamedTemporaryFile(suffix=f"_mc_{pid}_{rep}.fasta", delete=False)
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
    )
    t_rhiz_ms = (time.perf_counter() - t0) * 1000.0
    
    rhiz_unique_events = 0
    rhiz_raw_alerts = 0
    rhiz_recs_called = set()
    rhiz_bg_called = set()
    rhiz_best_jaccard = 0.0
    rhiz_best_bp_error = float(L)
    
    if r_proc.returncode == 0 and os.path.exists(out_json):
        try:
            with open(out_json) as f:
                r_res = json.load(f)
            events = r_res.get("events", [])
            rhiz_unique_events = len(events)
            
            for ev in events:
                rhiz_raw_alerts += len(ev.get("isolates", []))
                for iso in ev.get("isolates", []):
                    if iso in true_recs:
                        rhiz_recs_called.add(iso)
                    elif iso in bg_taxa:
                        rhiz_bg_called.add(iso)
                pred_u1 = ev.get("u1", 0)
                pred_u2 = ev.get("u2", 0)
                jac = compute_jaccard(true_u1, true_u2, pred_u1, pred_u2)
                bp_err = abs(pred_u1 - true_u1) + abs(pred_u2 - true_u2)
                if jac > rhiz_best_jaccard:
                    rhiz_best_jaccard = jac
                    rhiz_best_bp_error = bp_err
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
        "tract_len": tract_len,
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "mut_in_tract": mut_in_tract,
        "total_mut": total_mut,
        "unique_events": rhiz_unique_events,
        "raw_alerts": rhiz_raw_alerts,
        "unique_rec_taxa": len(rhiz_recs_called),
        "recall": rhiz_recall,
        "precision": rhiz_precision,
        "bg_fp_count": len(rhiz_bg_called),
        "bg_fpr": rhiz_bg_fpr,
        "jaccard": rhiz_best_jaccard,
        "bp_error": rhiz_best_bp_error if rhiz_unique_events > 0 else float(L),
        "time_ms": t_rhiz_ms,
    })
    
    # ---------------- 2. 3SEQ ----------------
    run_id = f"3seq_{pid}_{tract_len}_{int(d_parent*100)}_{K}_{rep}"
    t0 = time.perf_counter()
    subprocess.run(
        [three_bin, "-full", tmp_fa.name, "-id", run_id],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    t_3seq_ms = (time.perf_counter() - t0) * 1000.0
    
    rec_csv = f"{run_id}.3s.rec.csv"
    three_triplets = 0
    three_recs_called = set()
    three_bg_called = set()
    three_best_jaccard = 0.0
    three_best_bp_error = float(L)
    
    if os.path.exists(rec_csv):
        try:
            with open(rec_csv) as f:
                for line in f:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) >= 10:
                        try:
                            p_val = float(parts[9])
                            if p_val <= 0.05:
                                three_triplets += 1
                                child = parts[2]
                                if child in true_recs:
                                    three_recs_called.add(child)
                                elif child in bg_taxa:
                                    three_bg_called.add(child)
                                if len(parts) >= 13:
                                    bp_str = parts[12].strip()
                                    u1_3s, u2_3s = parse_3seq_bps(bp_str)
                                    if u1_3s > 0 and u2_3s > u1_3s:
                                        jac = compute_jaccard(true_u1, true_u2, u1_3s, u2_3s)
                                        bp_err = abs(u1_3s - true_u1) + abs(u2_3s - true_u2)
                                        if jac > three_best_jaccard:
                                            three_best_jaccard = jac
                                            three_best_bp_error = bp_err
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
    
    results.append({
        "tool": "3SEQ",
        "tract_len": tract_len,
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "mut_in_tract": mut_in_tract,
        "total_mut": total_mut,
        "unique_events": len(three_recs_called) + len(three_bg_called),
        "raw_alerts": three_triplets,
        "unique_rec_taxa": len(three_recs_called),
        "recall": three_recall,
        "precision": three_precision,
        "bg_fp_count": len(three_bg_called),
        "bg_fpr": three_bg_fpr,
        "jaccard": three_best_jaccard,
        "bp_error": three_best_bp_error if three_triplets > 0 else float(L),
        "time_ms": t_3seq_ms,
    })
    
    if os.path.exists(tmp_fa.name):
        os.remove(tmp_fa.name)
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Run Scenario 1.4 Bacterial Micro-Conversion Benchmark")
    parser.add_argument("--replicates", type=int, default=20, help="Replicates per condition (default 20)")
    parser.add_argument("--threads", type=int, default=min(8, cpu_count()), help="Multiprocessing workers")
    parser.add_argument("--rhiz-bin", type=str, default="/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon", help="Path to rhizaeon binary")
    parser.add_argument("--three-bin", type=str, default="/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq", help="Path to 3seq binary")
    parser.add_argument("--out-dir", type=str, default="/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario4_micro_conversion_tracts", help="Directory for output CSVs")
    args = parser.parse_args()
    
    rhiz_bin = args.rhiz_bin
    three_bin = args.three_bin
    
    assert os.path.exists(rhiz_bin), f"RhizAeon binary not found at {rhiz_bin}"
    assert os.path.exists(three_bin), f"3SEQ binary not found at {three_bin}"
    
    tract_lengths = [50, 100, 200, 400, 800]
    divergences = [0.02, 0.05, 0.10]
    cohort_sizes = [2, 5, 10]
    
    jobs = []
    for t_len in tract_lengths:
        for d in divergences:
            for k in cohort_sizes:
                for rep in range(args.replicates):
                    jobs.append((t_len, d, k, rep, rhiz_bin, three_bin))
                    
    total_jobs = len(jobs)
    print(f"Executing Scenario 1.4 Benchmark: {total_jobs} trials across {args.threads} workers...")
    
    t0 = time.time()
    all_results = []
    
    with Pool(processes=args.threads) as pool:
        for i, res_pair in enumerate(pool.imap_unordered(run_single_trial, jobs)):
            all_results.extend(res_pair)
            if (i + 1) % 50 == 0 or (i + 1) == total_jobs:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                rem = (total_jobs - (i + 1)) / rate if rate > 0 else 0
                print(f"  Progress: {i+1}/{total_jobs} trials completed ({rate:.1f} trials/s, ~{rem:.0f}s remaining)...")
                
    df_raw = pd.DataFrame(all_results)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = out_dir / "micro_conversion_raw_results.csv"
    df_raw.to_csv(raw_csv, index=False)
    print(f"\nRaw results saved to {raw_csv}")
    
    # Aggregate summary
    summary_cols = ["tract_len", "d_parent", "K", "tool"]
    metric_cols = ["total_taxa", "mut_in_tract", "unique_events", "raw_alerts", "recall", "precision", "bg_fpr", "jaccard", "bp_error", "time_ms"]
    
    df_summary = df_raw.groupby(summary_cols)[metric_cols].mean().reset_index()
    summary_csv = out_dir / "micro_conversion_summary.csv"
    df_summary.to_csv(summary_csv, index=False)
    print(f"Summary table saved to {summary_csv}")
    
    # Print formatted table
    print("\n" + "=" * 125)
    print(f"{'Tract':>5} {'d':>5} {'K':>3} {'N':>3} {'Mut':>4} | {'3S_K':>5} {'3S_Rec':>7} {'3S_FPR':>7} {'3S_IoU':>6} {'3S_ms':>6} | {'Rh_K':>5} {'Rh_Rec':>7} {'Rh_FPR':>7} {'Rh_IoU':>6} {'Rh_ms':>6}")
    print("=" * 125)
    
    # Pivot for side-by-side comparison
    piv = df_summary.pivot_table(index=["tract_len", "d_parent", "K", "total_taxa", "mut_in_tract"], columns="tool")
    
    for (t_len, d, k, n_taxa, muts), row in piv.iterrows():
        print(f"{t_len:5d} {d:5.2f} {k:3d} {int(n_taxa):3d} {muts:4.1f} | "
              f"{row[('unique_events', '3SEQ')]:5.1f} {row[('recall', '3SEQ')]*100:6.1f}% {row[('bg_fpr', '3SEQ')]*100:6.1f}% {row[('jaccard', '3SEQ')]:6.3f} {row[('time_ms', '3SEQ')]:6.0f} | "
              f"{row[('unique_events', 'RhizAeon')]:5.2f} {row[('recall', 'RhizAeon')]*100:6.1f}% {row[('bg_fpr', 'RhizAeon')]*100:6.1f}% {row[('jaccard', 'RhizAeon')]:6.3f} {row[('time_ms', 'RhizAeon')]:6.0f}")
    print("=" * 125)
    print(f"Total benchmark elapsed time: {time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
