#!/usr/bin/env python3
"""
Scenario 1.1: Epidemic Clonal Expansion Benchmark Harness.
Compares RhizAeon vs 3SEQ on a single ancestral recombination event
that undergoes epidemic clonal expansion into K descendant isolates
with independent transmission drift.

Ground Truth:
  K_true = 1 ancestral recombination event
  True tract = [1500, 3000] (1,501 nt cassette in 5,000 nt genome)
  True recombinant roster = {Rec_0, ..., Rec_{K-1}}
  True background non-recombinants = {P1_0, ..., P1_{N1-1}, P2_0, ..., P2_{N2-1}}
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

def generate_epidemic_alignment(
    d_parent: float = 0.05,
    K: int = 10,
    n_p1: int = 20,
    n_p2: int = 20,
    L: int = 5000,
    u1: int = 1500,
    u2: int = 3000,
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
        
    # 5. Ancestral recombinant node A_rec
    rec_anc = p1_anc.copy()
    rec_anc[u1 - 1 : u2] = p2_anc[u1 - 1 : u2]
    
    # 6. Sprout K descendant recombinant isolates with transmission drift
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat

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

def run_single_trial(args):
    d_parent, K, rep, rhiz_bin, three_bin = args
    seed = 100000 + int(d_parent * 1000) * 1000 + K * 50 + rep
    n_p1, n_p2 = 20, 20
    L = 5000
    true_u1, true_u2 = 1500, 3000
    true_recs = set(f"Rec_{k}" for k in range(K))
    bg_taxa = set(f"P1_{i}" for i in range(n_p1)).union(f"P2_{j}" for j in range(n_p2))
    
    taxa, mat = generate_epidemic_alignment(
        d_parent=d_parent, K=K, n_p1=n_p1, n_p2=n_p2, L=L, u1=true_u1, u2=true_u2, seed=seed
    )
    
    # Write temporary fasta
    pid = os.getpid()
    tmp_fa = tempfile.NamedTemporaryFile(suffix=f"_{pid}_{rep}.fasta", delete=False)
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
    rhiz_jaccard = 0.0
    
    if r_proc.returncode == 0 and os.path.exists(out_json):
        with open(out_json) as f:
            r_res = json.load(f)
        events = r_res.get("events", [])
        rhiz_unique_events = len(events)
        
        # In RhizAeon, each event has an isolates roster
        best_jaccard = 0.0
        for ev in events:
            rhiz_raw_alerts += len(ev.get("isolates", []))
            for iso in ev.get("isolates", []):
                if iso in true_recs:
                    rhiz_recs_called.add(iso)
                elif iso in bg_taxa:
                    rhiz_bg_called.add(iso)
            jac = compute_jaccard(true_u1, true_u2, ev.get("u1", 0), ev.get("u2", 0))
            if jac > best_jaccard:
                best_jaccard = jac
        rhiz_jaccard = best_jaccard
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
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "unique_events": rhiz_unique_events,
        "raw_alerts": rhiz_raw_alerts,
        "unique_rec_taxa": len(rhiz_recs_called),
        "recall": rhiz_recall,
        "precision": rhiz_precision,
        "bg_fp_count": len(rhiz_bg_called),
        "bg_fpr": rhiz_bg_fpr,
        "jaccard": rhiz_jaccard,
        "time_ms": t_rhiz_ms,
    })
    
    # ---------------- 2. 3SEQ ----------------
    run_id = f"3seq_{pid}_{d_parent}_{K}_{rep}"
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
    
    if os.path.exists(rec_csv):
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
                                    if jac > three_best_jaccard:
                                        three_best_jaccard = jac
                    except Exception:
                        pass
                        
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
        "d_parent": d_parent,
        "K": K,
        "rep": rep,
        "total_taxa": len(taxa),
        "unique_events": len(three_recs_called) + len(three_bg_called), # 3SEQ cannot collapse; each child is a separate event
        "raw_alerts": three_triplets,
        "unique_rec_taxa": len(three_recs_called),
        "recall": three_recall,
        "precision": three_precision,
        "bg_fp_count": len(three_bg_called),
        "bg_fpr": three_bg_fpr,
        "jaccard": three_best_jaccard,
        "time_ms": t_3seq_ms,
    })
    
    if os.path.exists(tmp_fa.name):
        os.remove(tmp_fa.name)
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Run Scenario 1.1 Epidemic Clonal Expansion Benchmark")
    parser.add_argument("--replicates", type=int, default=20, help="Replicates per condition (default 20)")
    parser.add_argument("--threads", type=int, default=min(8, cpu_count()), help="Multiprocessing workers")
    args = parser.parse_args()
    
    rhiz_bin = "/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon"
    three_bin = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
    
    assert os.path.exists(rhiz_bin), f"RhizAeon binary not found at {rhiz_bin}"
    assert os.path.exists(three_bin), f"3SEQ binary not found at {three_bin}"
    
    divergences = [0.02, 0.05, 0.10]
    cohort_sizes = [2, 5, 8, 10, 15]
    
    jobs = []
    for d in divergences:
        for k in cohort_sizes:
            for rep in range(args.replicates):
                jobs.append((d, k, rep, rhiz_bin, three_bin))
                
    total_jobs = len(jobs)
    print(f"Executing Scenario 1.1 Benchmark: {total_jobs} trials across {args.threads} workers...")
    
    t0 = time.time()
    all_results = []
    with Pool(args.threads) as pool:
        for idx, res_list in enumerate(pool.imap_unordered(run_single_trial, jobs), 1):
            all_results.extend(res_list)
            if idx % 25 == 0 or idx == total_jobs:
                elapsed = time.time() - t0
                print(f"  Processed {idx}/{total_jobs} alignments ({idx/total_jobs*100:.1f}%) in {elapsed:.1f}s")
                
    # Save raw results
    out_dir = Path("/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion")
    raw_csv = out_dir / "epidemic_raw_results.csv"
    
    import csv
    fieldnames = list(all_results[0].keys())
    with open(raw_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nRaw results exported to: {raw_csv}")
    
    # Compute summary table
    summary_csv = out_dir / "epidemic_summary.csv"
    import pandas as pd
    df = pd.DataFrame(all_results)
    summary = df.groupby(["tool", "d_parent", "K"]).agg(
        total_taxa=("total_taxa", "mean"),
        unique_events_mean=("unique_events", "mean"),
        raw_alerts_mean=("raw_alerts", "mean"),
        recall_mean=("recall", "mean"),
        precision_mean=("precision", "mean"),
        bg_fp_mean=("bg_fp_count", "mean"),
        bg_fpr_mean=("bg_fpr", "mean"),
        jaccard_mean=("jaccard", "mean"),
        time_ms_mean=("time_ms", "mean"),
    ).reset_index()
    
    summary.to_csv(summary_csv, index=False)
    print(f"Summary table exported to: {summary_csv}")
    
    # Print formatted markdown table
    print("\n" + "=" * 110)
    print("SCENARIO 1.1: EPIDEMIC CLONAL EXPANSION SUMMARY (RhizAeon vs 3SEQ)")
    print("=" * 110)
    print(f"{'Tool':10s} | {'d_par':5s} | {'K':3s} | {'Events (True=1)':15s} | {'Raw Alerts':11s} | {'Recall':7s} | {'Precision':9s} | {'Parent FPR':10s} | {'Jaccard':7s} | {'Time (ms)':9s}")
    print("-" * 110)
    for _, row in summary.iterrows():
        print(
            f"{row['tool']:10s} | {row['d_parent']:5.2f} | {int(row['K']):3d} | "
            f"{row['unique_events_mean']:15.2f} | {row['raw_alerts_mean']:11.1f} | "
            f"{row['recall_mean']:7.2f} | {row['precision_mean']:9.2f} | "
            f"{row['bg_fpr_mean']:10.2f} | {row['jaccard_mean']:7.3f} | {row['time_ms_mean']:9.1f}"
        )
    print("=" * 110)

if __name__ == "__main__":
    main()
