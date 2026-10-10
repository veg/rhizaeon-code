#!/usr/bin/env python3
"""
benchmarks/posada_crandall_2001/run_benchmark_3seq.py
=====================================================
Head-to-head 3SEQ evaluation across the exact same 32 Posada & Crandall (2001)
simulation scenarios using identical random seeds, alignments, and evaluation metrics.
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
from typing import Dict, Any, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from run_benchmark import (
    PosadaScenario, build_scenarios, simulate_hudson_arg, evolve_sequences,
    extract_ground_truth_events, compute_patristic_distances,
    compute_reticulation_metrics, MASTER_SEED, NT_CHARS
)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
THREE_SEQ_BIN = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
PTABLE_PATH = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq_ptable_250"


def evaluate_single_replicate_3seq(scenario: PosadaScenario, sc_idx: int, rep_idx: int) -> Dict[str, Any]:
    seed = MASTER_SEED + sc_idx * 10000 + rep_idx
    # 1. Simulate ARG & sequences (identical seed and generator)
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

    root, ch, bl = parts[0][2], parts[0][3], parts[0][4]
    D_tree = compute_patristic_distances(scenario.n_taxa, root, ch, bl)
    true_events = extract_ground_truth_events(parts, scenario.n_taxa, scenario.length_nt, mat, taxa) if scenario.rho > 0 else []

    # 2. Write FASTA
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f_fa:
        for i, t in enumerate(taxa):
            f_fa.write(">" + t + "\n" + "".join(NT_CHARS[c] for c in mat[i]) + "\n")
        fa_path = f_fa.name

    # 3. Run 3SEQ
    t0 = time.perf_counter()
    pred_events = []
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [THREE_SEQ_BIN, "-full", fa_path, "-ptable", PTABLE_PATH, "-id", "run"]
        subprocess.run(cmd, cwd=tmpdir, capture_output=True, text=True)
        t_3seq_ms = (time.perf_counter() - t0) * 1000.0
        csv_file = os.path.join(tmpdir, "run.3s.rec.csv")
        if os.path.exists(csv_file):
            try:
                with open(csv_file, "r") as f:
                    lines = [l.strip() for l in f if l.strip()]
                for line in lines[1:]:
                    parts_rec = [p.strip() for p in line.split(",")]
                    if len(parts_rec) >= 13:
                        dad = parts_rec[0]
                        mum = parts_rec[1]
                        child = parts_rec[2]
                        try:
                            ds_p = float(parts_rec[9])
                        except Exception:
                            try:
                                ds_p = float(parts_rec[6])
                            except Exception:
                                ds_p = 1.0
                        if ds_p <= 0.05:
                            for bp_str in parts_rec[12:]:
                                if " & " in bp_str:
                                    left_str, right_str = bp_str.split(" & ")
                                    try:
                                        l_min, l_max = map(int, left_str.strip().split("-"))
                                        r_min, r_max = map(int, right_str.strip().split("-"))
                                        pred_events.append({
                                            "candidate_name": child,
                                            "donor_name": dad,
                                            "home_name": mum,
                                            "u1": max(0, int((l_min + l_max) / 2.0) - 1),
                                            "u2": min(scenario.length_nt, int((r_min + r_max) / 2.0)),
                                            "p_val": ds_p
                                        })
                                    except Exception:
                                        continue
            except Exception:
                pass

    try:
        os.remove(fa_path)
    except Exception:
        pass

    # 4. Compute metrics (same metric evaluator)
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
        "runtime_ms": t_3seq_ms
    }


def main():
    parser = argparse.ArgumentParser(description="Posada & Crandall (2001) Replication for 3SEQ")
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
    print("  POSADA & CRANDALL (2001, PNAS) REPLICATION BENCHMARK -- 3SEQ")
    print(f"  Scenarios: {len(scenarios)} conditions | Replicates: {args.reps} | Total trials: {len(scenarios) * args.reps}")
    print(f"  Workers: {args.workers} | Engine: {THREE_SEQ_BIN}")
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
        futures = {executor.submit(evaluate_single_replicate_3seq, sc, sc_idx, rep): (sc.scenario_id, rep) for sc, sc_idx, rep in tasks}
        for future in as_completed(futures):
            res = future.result()
            raw_results.append(res)
            completed += 1
            if completed % 100 == 0 or completed == total:
                elapsed = time.perf_counter() - t_start
                rate = completed / elapsed
                print(f"  Progress: {completed:4d}/{total:4d} ({completed/total*100:5.1f}%) | {rate:5.1f} trials/s | Elapsed: {elapsed:5.1f}s")

    df_raw = pd.DataFrame(raw_results)
    raw_csv = SCRIPT_DIR / "posada_raw_results_3seq.csv"
    df_raw.to_csv(raw_csv, index=False)
    print(f"\n[DONE] Saved raw 3SEQ results ({len(df_raw)} trials) to {raw_csv}")

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
    summary_csv = SCRIPT_DIR / "posada_summary_3seq.csv"
    df_summary.to_csv(summary_csv, index=False)
    print(f"[DONE] Saved summary table ({len(df_summary)} scenarios) to {summary_csv}")
    print("=" * 80)


if __name__ == "__main__":
    main()
