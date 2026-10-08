#!/usr/bin/env python3
"""
evaluate_large_cohort_nulls.py
==============================
Stress test RhizAeon (PA-FF v2.1) on large-cohort adversarial nulls:
Taxa counts: N in {24, 32, 48, 64, 96, 128}
Regimes:
  - Clock Null
  - Gamma ASRV Null
  - Darren Heterotachy Trap (private accelerated tract)
  - APOBEC Shower (clustered mutational burst)
  - Darwinian Selection Burst
"""

import sys
import os
import time
import json
import copy
import subprocess
from pathlib import Path
from typing import Dict, Any, List
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

BENCH_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/19_grand_100k_paff_v2")
if str(BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(BENCH_DIR))

from scenario_registry import get_scenario_for_index
from fast_tree_simulator import FastTreeSimulator

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RUST_BIN = REPO_ROOT / "target" / "release" / "rhizaeon"

LARGE_TAXA_COUNTS = [24, 32, 48, 64, 96, 128]
NULL_BASE_INDICES = {
    "clock_null": 100,
    "gamma_asrv_null": 7100,
    "heterotachy_darren_null": 14100,
    "apobec_shower_null": 21100,
    "darwinian_burst_null": 28100,
}

def run_large_null_trial(regime: str, n_taxa: int, rep: int) -> Dict[str, Any]:
    base_idx = NULL_BASE_INDICES[regime] + rep * 7
    sc = copy.deepcopy(get_scenario_for_index(base_idx))
    sc.num_taxa = n_taxa
    sc.scenario_id = f"large_{regime}_N{n_taxa}_rep{rep}"

    sim = FastTreeSimulator(sc)
    seqs, _ = sim.simulate()

    pid = os.getpid()
    shm = Path("/dev/shm") if Path("/dev/shm").exists() else Path("/tmp")
    fa_path = shm / f"large_null_{regime}_{n_taxa}_{rep}_{pid}.fasta"
    js_path = shm / f"large_null_{regime}_{n_taxa}_{rep}_{pid}.json"

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
        "regime": regime,
        "num_taxa": n_taxa,
        "length_nt": sc.length_nt,
        "rep": rep,
        "num_events": num_events,
        "false_positive": is_fp,
        "runtime_ms": t_ms,
        "events_summary": json.dumps([{"cand": e.get("candidate_name"), "home": e.get("home_name"), "donor": e.get("donor_name"), "p": e.get("p_fisher")} for e in events])
    }

def main():
    print("================================================================================")
    print("  RHIZAEON LARGE-COHORT ADVERSARIAL NULL AUDIT")
    print(f"  Taxa counts: {LARGE_TAXA_COUNTS}")
    print(f"  Regimes:     {list(NULL_BASE_INDICES.keys())}")
    print(f"  Replicates:  5 per (regime, N) combination (150 total deep-cohort trials)")
    print("================================================================================")

    tasks = []
    for regime in NULL_BASE_INDICES:
        for n in LARGE_TAXA_COUNTS:
            for rep in range(5):
                tasks.append((regime, n, rep))

    results = []
    t_start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(run_large_null_trial, r, n, rep): (r, n, rep) for r, n, rep in tasks}
        completed = 0
        total = len(tasks)
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            completed += 1
            if completed % 25 == 0 or completed == total:
                elapsed = time.perf_counter() - t_start
                rate = completed / elapsed
                print(f"  Completed {completed:3d}/{total} ({completed/total*100:5.1f}%) [{rate:.1f} trials/s]")

    df = pd.DataFrame(results)
    out_csv = REPO_ROOT / "large_cohort_null_results.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nSaved raw results to {out_csv}")

    print("\n================================================================================")
    print("  LARGE-COHORT NULL AUDIT BREAKDOWN BY TAXA COUNT N")
    print("================================================================================")
    by_n = df.groupby("num_taxa").agg(
        trials=("false_positive", "count"),
        false_positives=("false_positive", "sum"),
        mean_runtime_ms=("runtime_ms", "mean")
    ).reset_index()
    by_n["fpr"] = (by_n["false_positives"] / by_n["trials"] * 100.0).map("{:.2f}%".format)
    by_n["mean_runtime_ms"] = by_n["mean_runtime_ms"].map("{:.1f} ms".format)
    print(by_n.to_string(index=False))

    print("\n================================================================================")
    print("  LARGE-COHORT NULL AUDIT BREAKDOWN BY ADVERSARIAL REGIME")
    print("================================================================================")
    by_reg = df.groupby("regime").agg(
        trials=("false_positive", "count"),
        false_positives=("false_positive", "sum"),
        mean_runtime_ms=("runtime_ms", "mean")
    ).reset_index()
    by_reg["fpr"] = (by_reg["false_positives"] / by_reg["trials"] * 100.0).map("{:.2f}%".format)
    by_reg["mean_runtime_ms"] = by_reg["mean_runtime_ms"].map("{:.1f} ms".format)
    print(by_reg.to_string(index=False))

    total_fp = df["false_positive"].sum()
    print("--------------------------------------------------------------------------------")
    print(f"  TOTAL LARGE COHORT NULL TRIALS: {len(df)}")
    print(f"  TOTAL FALSE POSITIVES:          {total_fp}")
    print(f"  OVERALL FPR (N=24..128):        {(total_fp / len(df) * 100.0):.2f}%")
    print("================================================================================")

if __name__ == "__main__":
    main()
