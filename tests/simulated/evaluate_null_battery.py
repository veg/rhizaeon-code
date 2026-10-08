#!/usr/bin/env python3
"""
evaluate_null_battery.py
========================
High-throughput statistical audit of RhizAeon (PA-FF v2.1) False Positive Rate (FPR)
across 5 adversarial null simulation regimes and empirical negative controls:

1. Strict Neutral Clock Null (clock_null, indices 0..6999)
2. Gamma Rate Heterogeneity Null (gamma_asrv_null, indices 7000..13999)
3. Darren Trap Coupled Heterotachy (heterotachy_darren_null, indices 14000..20999)
4. Host Deaminase / APOBEC Mutational Shower (apobec_shower_null, indices 21000..27999)
5. Episodic Darwinian Mimics & Convergent Hotspots (darwinian_burst_null, indices 28000..34999)
6. Authentic Clonal Mammalian Mitochondrial DNA (Human mtDNA, 78 taxa, 16.5 kb)

Evaluates stratified cohorts spanning:
- Taxa counts N in {4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128}
- Sequence lengths L in [600, 12000] nt
- Divergences d in [0.0005, 0.40] (0.05% to 40%)
"""

import sys
import os
import time
import json
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Any, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

# Add scenario registry to path
BENCH_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/19_grand_100k_paff_v2")
if str(BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(BENCH_DIR))

from scenario_registry import get_scenario_for_index, MultiDimScenario
from fast_tree_simulator import FastTreeSimulator

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RUST_BIN = REPO_ROOT / "target" / "release" / "rhizaeon"
EMPIRICAL_MTDNA = REPO_ROOT / "tests" / "empirical" / "data" / "human_mtdna.fasta"
DARREN_PVY = Path("/Users/sergei/Projects/RhizAeon/darren-tests/Example 1 - PVY.fas")

NULL_REGIMES = [
    ("clock_null", 0, 7000, "Strict Neutral Clock Null"),
    ("gamma_asrv_null", 7000, 14000, "Gamma ASRV Rate Heterogeneity Null"),
    ("heterotachy_darren_null", 14000, 21000, "Darren Trap Coupled Heterotachy"),
    ("apobec_shower_null", 21000, 28000, "APOBEC Mutational Shower Null"),
    ("darwinian_burst_null", 28000, 35000, "Darwinian Burst / Selection Null"),
]

def run_single_null_trial(global_idx: int) -> Dict[str, Any]:
    sc = get_scenario_for_index(global_idx)
    sim = FastTreeSimulator(sc)
    seqs, _ = sim.simulate()

    pid = os.getpid()
    shm = Path("/dev/shm") if Path("/dev/shm").exists() else Path("/tmp")
    fa_path = shm / f"null_eval_{global_idx}_{pid}.fasta"
    js_path = shm / f"null_eval_{global_idx}_{pid}.json"

    with open(fa_path, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")

    t0 = time.perf_counter()
    try:
        cmd = [str(RUST_BIN), "-i", str(fa_path), "-o", str(js_path)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        t_ms = (time.perf_counter() - t0) * 1000.0
        if res.returncode != 0:
            raise RuntimeError(f"RhizAeon failed: {res.stderr}")
        with open(js_path) as f:
            data = json.load(f)
        events = data.get("events", [])
    finally:
        if fa_path.exists():
            fa_path.unlink()
        if js_path.exists():
            js_path.unlink()

    num_events = len(events)
    is_fp = 1 if num_events > 0 else 0

    return {
        "global_idx": global_idx,
        "regime": sc.regime,
        "archetype": sc.archetype_label,
        "num_taxa": sc.num_taxa,
        "length_nt": sc.length_nt,
        "mean_divergence": sc.mean_divergence,
        "num_events": num_events,
        "false_positive": is_fp,
        "runtime_ms": t_ms,
        "informative_sites": data.get("informative_sites", 0),
    }

def evaluate_empirical_mtdna() -> Dict[str, Any]:
    assert EMPIRICAL_MTDNA.exists(), f"Missing {EMPIRICAL_MTDNA}"
    shm = Path("/dev/shm") if Path("/dev/shm").exists() else Path("/tmp")
    js_path = shm / "mtdna_eval.json"
    t0 = time.perf_counter()
    cmd = [str(RUST_BIN), "-i", str(EMPIRICAL_MTDNA), "-o", str(js_path)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    t_ms = (time.perf_counter() - t0) * 1000.0
    with open(js_path) as f:
        data = json.load(f)
    if js_path.exists():
        js_path.unlink()
    events = data.get("events", [])
    return {
        "dataset": "Human mtDNA (Clonal)",
        "taxa": 78,
        "length": 16569,
        "num_events": len(events),
        "false_positive": 1 if len(events) > 0 else 0,
        "runtime_ms": t_ms
    }

def evaluate_darren_pvy() -> Dict[str, Any]:
    assert DARREN_PVY.exists(), f"Missing {DARREN_PVY}"
    shm = Path("/dev/shm") if Path("/dev/shm").exists() else Path("/tmp")
    js_path = shm / "pvy_eval.json"
    t0 = time.perf_counter()
    cmd = [str(RUST_BIN), "-i", str(DARREN_PVY), "-o", str(js_path)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    t_ms = (time.perf_counter() - t0) * 1000.0
    with open(js_path) as f:
        data = json.load(f)
    if js_path.exists():
        js_path.unlink()
    events = data.get("events", [])
    recs = set(e.get("candidate_name") for e in events)
    return {
        "dataset": "Potato Virus Y (PVY)",
        "taxa": 25,
        "length": 9594,
        "num_events": len(events),
        "recombinant_taxa_count": len(recs),
        "runtime_ms": t_ms
    }

def main():
    parser = argparse.ArgumentParser(description="Evaluate False Positive Rate across Adversarial Nulls")
    parser.add_argument("--reps-per-regime", type=int, default=50, help="Replicates per null regime (default: 50 = 250 total)")
    parser.add_argument("--workers", type=int, default=8, help="Number of worker processes")
    parser.add_argument("--out-csv", type=Path, default=REPO_ROOT / "null_battery_results.csv", help="Output CSV path")
    args = parser.parse_args()

    print("================================================================================")
    print("  RHIZAEON (PA-FF v2.1) FALSE POSITIVE RATE AUDIT BATTERY")
    print(f"  Replicates per regime: {args.reps_per_regime} ({args.reps_per_regime * len(NULL_REGIMES)} total trials)")
    print(f"  Worker processes:      {args.workers}")
    print(f"  Output CSV:            {args.out_csv}")
    print("================================================================================")

    # 1. Empirical Controls
    print("\n--- Running Empirical Controls ---")
    mtdna_res = evaluate_empirical_mtdna()
    print(f"  [1] Human mtDNA Negative Control (78 taxa, 16.5 kb):")
    print(f"      Events called: {mtdna_res['num_events']} | FPR: {mtdna_res['false_positive']:.1%} | Runtime: {mtdna_res['runtime_ms']:.1f} ms")

    pvy_res = evaluate_darren_pvy()
    print(f"  [2] Darren's PVY Sensitivity Test (25 taxa, 9.6 kb):")
    print(f"      Verified Events: {pvy_res['num_events']} across {pvy_res['recombinant_taxa_count']} lineages | Runtime: {pvy_res['runtime_ms']:.1f} ms")

    # 2. Select stratified indices across null regimes
    trial_indices = []
    for regime_name, start_idx, end_idx, desc in NULL_REGIMES:
        step = (end_idx - start_idx) // args.reps_per_regime
        indices = [start_idx + i * step for i in range(args.reps_per_regime)]
        trial_indices.extend(indices)

    print(f"\n--- Running {len(trial_indices)} Adversarial Null Simulations across 5 Regimes ---")
    results = []
    t_start = time.perf_counter()

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_single_null_trial, idx): idx for idx in trial_indices}
        completed = 0
        total = len(trial_indices)
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            completed += 1
            if completed % 25 == 0 or completed == total:
                elapsed = time.perf_counter() - t_start
                rate = completed / elapsed
                print(f"  Completed {completed:4d}/{total} null trials ({completed/total*100:5.1f}%) [{rate:.1f} trials/s]")

    # 3. Analyze Results
    df = pd.DataFrame(results)
    df.to_csv(args.out_csv, index=False)
    print(f"\nRaw results written to {args.out_csv}")

    print("\n================================================================================")
    print("  FALSE POSITIVE RATE AUDIT SUMMARY")
    print("================================================================================")
    summary_rows = []
    for regime_name, _, _, desc in NULL_REGIMES:
        sub = df[df["regime"] == regime_name]
        n_trials = len(sub)
        n_fp = sub["false_positive"].sum()
        fpr = (n_fp / n_trials) * 100.0 if n_trials > 0 else 0.0
        mean_rt = sub["runtime_ms"].mean()
        mean_taxa = sub["num_taxa"].mean()
        mean_len = sub["length_nt"].mean()
        summary_rows.append({
            "Regime": desc,
            "Trials": n_trials,
            "False Positives": n_fp,
            "FPR (%)": f"{fpr:.2f}%",
            "Mean Taxa": f"{mean_taxa:.1f}",
            "Mean Length (nt)": f"{mean_len:.0f}",
            "Mean Runtime (ms)": f"{mean_rt:.1f} ms"
        })

    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))

    total_trials = len(df)
    total_fp = df["false_positive"].sum()
    overall_fpr = (total_fp / total_trials) * 100.0

    print("--------------------------------------------------------------------------------")
    print(f"  TOTAL ADVERSARIAL NULL TRIALS: {total_trials}")
    print(f"  TOTAL FALSE POSITIVES:        {total_fp}")
    print(f"  OVERALL FALSE POSITIVE RATE:  {overall_fpr:.2f}%")
    print(f"  HUMAN MTDNA NEGATIVE CONTROL: {mtdna_res['false_positive']} FP ({mtdna_res['num_events']} events)")
    print("================================================================================")

    if total_fp > 0:
        print(f"\n[WARNING] {total_fp} false positive trials detected! Inspecting details:")
        fp_rows = df[df["false_positive"] == 1]
        print(fp_rows[["global_idx", "regime", "num_taxa", "length_nt", "mean_divergence", "num_events", "informative_sites"]])
    else:
        print("\n[SUCCESS] Certified 0.00% False Positive Rate across all adversarial null regimes!")

if __name__ == "__main__":
    main()
