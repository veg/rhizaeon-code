#!/usr/bin/env python3
"""
evaluate_rust_battery.py
========================
High-throughput evaluation harness for Native Rust RhizAeon (PA-FF v2.1)
across the 100,000-simulation Multi-Dimensional Scaling Recombination Benchmark.

Evaluates:
  1. Recombination detection power / sensitivity.
  2. False positive rate (FPR) across adversarial null regimes.
  3. Correct parental lineage attribution (Home, Donor, Joint).
  4. Breakpoint resolution and missed breakpoints (false negatives).
"""

import sys
import os
import time
import json
import sqlite3
import argparse
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Set
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.optimize import linear_sum_assignment

# Paths
GRAND_100K_DIR = Path("/data/users/sergei/Projects/19_grand_100k_paff_v2")
if str(GRAND_100K_DIR) not in sys.path:
    sys.path.insert(0, str(GRAND_100K_DIR))

from scenario_registry import get_scenario_for_index, MultiDimScenario
from fast_tree_simulator import FastTreeSimulator

RUST_BIN = "/data/users/sergei/Projects/RhizAeon/rhizaeon-rust/target/release/rhizaeon"

REGIME_SLICES = [
    ("clock_null", 0, 7000, True),
    ("gamma_asrv_null", 7000, 14000, True),
    ("heterotachy_darren_null", 14000, 21000, True),
    ("apobec_shower_null", 21000, 28000, True),
    ("darwinian_burst_null", 28000, 35000, True),
    ("single_crossover", 35000, 48000, False),
    ("cassette_swap", 48000, 63000, False),
    ("micro_conversion", 63000, 76000, False),
    ("ghost_donor_introgression", 76000, 88000, False),
    ("multiway_mosaic", 88000, 100000, False),
]

DB_SCHEMA = """
CREATE TABLE IF NOT EXISTS rust_benchmark_results (
    global_idx INTEGER PRIMARY KEY,
    scenario_id TEXT,
    partition TEXT,
    regime TEXT,
    archetype_label TEXT,
    is_null INTEGER,
    num_taxa INTEGER,
    length_nt INTEGER,
    mean_divergence REAL,
    detected INTEGER,
    true_positive INTEGER,
    false_positive INTEGER,
    num_events INTEGER,
    inferred_recombinant TEXT,
    inferred_home TEXT,
    inferred_donor TEXT,
    correct_recombinant_exact INTEGER,
    correct_recombinant_sister INTEGER,
    correct_recombinant_in_set INTEGER,
    correct_home INTEGER,
    correct_donor INTEGER,
    correct_both_parents INTEGER,
    true_breakpoints TEXT,
    inferred_breakpoints TEXT,
    num_true_bps INTEGER,
    num_inferred_bps INTEGER,
    matched_bps_50 INTEGER,
    missed_bps_50 INTEGER,
    hallucinated_bps_50 INTEGER,
    matched_bps_100 INTEGER,
    missed_bps_100 INTEGER,
    hallucinated_bps_100 INTEGER,
    bp_error_mean REAL,
    bp_error_median REAL,
    bp_all_within_35 INTEGER,
    bp_all_within_50 INTEGER,
    bp_all_within_100 INTEGER,
    informative_sites INTEGER,
    runtime_rust_ms REAL,
    runtime_sim_ms REAL,
    runtime_total_ms REAL,
    seed INTEGER,
    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_regime ON rust_benchmark_results(regime);
CREATE INDEX IF NOT EXISTS idx_is_null ON rust_benchmark_results(is_null);
CREATE INDEX IF NOT EXISTS idx_taxa ON rust_benchmark_results(num_taxa);
CREATE INDEX IF NOT EXISTS idx_div ON rust_benchmark_results(mean_divergence);
"""

def evaluate_single_trial(global_idx: int) -> Dict[str, Any]:
    """Simulates scenario and evaluates using native Rust RhizAeon."""
    t_start_total = time.perf_counter()
    sc = get_scenario_for_index(global_idx)

    # 1. Simulate alignment
    t_sim_0 = time.perf_counter()
    sim = FastTreeSimulator(sc)
    seqs, meta = sim.simulate()
    t_sim_ms = (time.perf_counter() - t_sim_0) * 1000.0

    pid = os.getpid()
    shm_base = Path("/dev/shm") if Path("/dev/shm").exists() else Path("/tmp")
    fasta_path = shm_base / f"trial_{global_idx}_{pid}.fasta"
    json_path = shm_base / f"trial_{global_idx}_{pid}.json"

    with open(fasta_path, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")

    # 2. Run Rust RhizAeon
    t_rust_0 = time.perf_counter()
    try:
        cmd = [RUST_BIN, "-i", str(fasta_path), "-o", str(json_path)]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        t_rust_ms = (time.perf_counter() - t_rust_0) * 1000.0
        with open(json_path, "r") as f:
            res_json = json.load(f)
    except Exception as e:
        t_rust_ms = (time.perf_counter() - t_rust_0) * 1000.0
        res_json = {"events": [], "informative_sites": 0, "run_time_ms": t_rust_ms}
    finally:
        if fasta_path.exists():
            fasta_path.unlink()
        if json_path.exists():
            json_path.unlink()

    # 3. Analyze detection
    events = res_json.get("events", [])
    detected = len(events) > 0
    is_null = sc.is_null
    tp = 1 if (detected and not is_null) else 0
    fp = 1 if (detected and is_null) else 0

    # 4. Extract inferred breakpoints and candidates
    L = sc.length_nt
    raw_bps: List[float] = []
    inferred_candidates: Set[str] = set()

    for ev in events:
        cand = ev.get("candidate_name", "")
        if cand:
            inferred_candidates.add(cand)
        is_cross = ev.get("is_crossover", False)
        u1_disc = ev.get("u1", 0)
        u2_disc = ev.get("u2", 0)
        u1_cont = ev.get("u1_continuous", float(u1_disc))
        u2_cont = ev.get("u2_continuous", float(u2_disc))

        if is_cross:
            if u1_disc <= 50:
                raw_bps.append(u2_cont)
            elif u2_disc >= L - 50:
                raw_bps.append(u1_cont)
            else:
                if u1_disc < (L - u2_disc):
                    raw_bps.append(u2_cont)
                else:
                    raw_bps.append(u1_cont)
        else:
            raw_bps.append(u1_cont)
            raw_bps.append(u2_cont)

    # Deduplicate within 5 nt
    all_inferred_bps: List[float] = []
    for b in sorted(raw_bps):
        if not all_inferred_bps or min(abs(b - ub) for ub in all_inferred_bps) > 5.0:
            all_inferred_bps.append(b)

    # 5. Parental lineage attribution
    true_recs = sc.recombinant_taxa
    true_parents = sc.parent_taxa
    true_bps = sc.true_breakpoints

    correct_rec_exact = 0
    correct_rec_sister = 0
    correct_rec_in_set = 0
    correct_home = 0
    correct_donor = 0
    correct_both = 0

    top_rec = events[0].get("candidate_name", "-") if events else "-"
    top_home = events[0].get("home_name", "-") if events else "-"
    top_donor = ("Ghost" if events[0].get("is_ghost_donor", False) else events[0].get("donor_name", "-")) if events else "-"

    if not is_null and detected and true_recs:
        target_rec = true_recs[0]
        expected_pars = true_parents.get(target_rec, [])
        is_true_ghost = any("Ghost" in p for p in expected_pars)
        expected_home = expected_pars[0] if len(expected_pars) > 0 else "P1"
        expected_donor = expected_pars[1] if len(expected_pars) > 1 else ("Ghost" if is_true_ghost else "P2")

        if top_rec == target_rec:
            correct_rec_exact = 1
        if target_rec in inferred_candidates:
            correct_rec_in_set = 1
        if top_rec == target_rec or (top_home == target_rec and top_rec == expected_home):
            correct_rec_sister = 1

        rec_ev = next((e for e in events if e.get("candidate_name") == target_rec), events[0])
        ev_home = rec_ev.get("home_name", "")
        ev_is_ghost = rec_ev.get("is_ghost_donor", False)
        ev_donor = "Ghost" if ev_is_ghost else rec_ev.get("donor_name", "")

        # Home check (direct or reciprocal)
        if ev_home == expected_home or (top_rec == expected_home and top_home == target_rec):
            correct_home = 1

        # Donor check
        if is_true_ghost:
            if ev_is_ghost or "Ghost" in ev_donor or (top_donor == "Ghost"):
                correct_donor = 1
        else:
            if (ev_donor in expected_pars and not ev_is_ghost and ev_donor != ev_home) or (top_donor in expected_pars and not events[0].get("is_ghost_donor", False)):
                correct_donor = 1

        if (ev_home == expected_home and ((is_true_ghost and (ev_is_ghost or "Ghost" in ev_donor)) or (not is_true_ghost and ev_donor in expected_pars and not ev_is_ghost and ev_donor != ev_home))):
            correct_both = 1


    # 6. Breakpoint Bipartite Matching & Error
    num_true_bps = len(true_bps)
    num_inf_bps = len(all_inferred_bps)

    matched_50 = 0
    matched_100 = 0
    bp_error_mean = np.nan
    bp_error_median = np.nan
    bp_all_35 = 0
    bp_all_50 = 0
    bp_all_100 = 0

    if not is_null and num_true_bps > 0:
        if num_inf_bps > 0:
            cost = np.abs(np.array(true_bps)[:, None] - np.array(all_inferred_bps)[None, :])
            r_ind, c_ind = linear_sum_assignment(cost)
            errors = [float(cost[r, c]) for r, c in zip(r_ind, c_ind)]
            bp_error_mean = float(np.mean(errors))
            bp_error_median = float(np.median(errors))

            matched_50 = sum(1 for e in errors if e <= 50.0)
            matched_100 = sum(1 for e in errors if e <= 100.0)

            if num_inf_bps == num_true_bps:
                bp_all_35 = 1 if all(e <= 35.0 for e in errors) else 0
                bp_all_50 = 1 if all(e <= 50.0 for e in errors) else 0
                bp_all_100 = 1 if all(e <= 100.0 for e in errors) else 0

        missed_50 = num_true_bps - matched_50
        missed_100 = num_true_bps - matched_100
        hallucinated_50 = max(0, num_inf_bps - matched_50)
        hallucinated_100 = max(0, num_inf_bps - matched_100)
    else:
        missed_50 = 0
        missed_100 = 0
        hallucinated_50 = num_inf_bps
        hallucinated_100 = num_inf_bps

    t_total_ms = (time.perf_counter() - t_start_total) * 1000.0

    return {
        "global_idx": global_idx,
        "scenario_id": sc.scenario_id,
        "partition": sc.partition,
        "regime": sc.regime,
        "archetype_label": sc.archetype_label,
        "is_null": 1 if is_null else 0,
        "num_taxa": sc.num_taxa,
        "length_nt": sc.length_nt,
        "mean_divergence": sc.mean_divergence,
        "detected": 1 if detected else 0,
        "true_positive": tp,
        "false_positive": fp,
        "num_events": len(events),
        "inferred_recombinant": top_rec,
        "inferred_home": top_home,
        "inferred_donor": top_donor,
        "correct_recombinant_exact": correct_rec_exact,
        "correct_recombinant_sister": correct_rec_sister,
        "correct_recombinant_in_set": correct_rec_in_set,
        "correct_home": correct_home,
        "correct_donor": correct_donor,
        "correct_both_parents": correct_both,
        "true_breakpoints": json.dumps(true_bps),
        "inferred_breakpoints": json.dumps(all_inferred_bps),
        "num_true_bps": num_true_bps,
        "num_inferred_bps": num_inf_bps,
        "matched_bps_50": matched_50,
        "missed_bps_50": missed_50,
        "hallucinated_bps_50": hallucinated_50,
        "matched_bps_100": matched_100,
        "missed_bps_100": missed_100,
        "hallucinated_bps_100": hallucinated_100,
        "bp_error_mean": bp_error_mean if not np.isnan(bp_error_mean) else None,
        "bp_error_median": bp_error_median if not np.isnan(bp_error_median) else None,
        "bp_all_within_35": bp_all_35,
        "bp_all_within_50": bp_all_50,
        "bp_all_within_100": bp_all_100,
        "informative_sites": res_json.get("informative_sites", 0),
        "runtime_rust_ms": t_rust_ms,
        "runtime_sim_ms": t_sim_ms,
        "runtime_total_ms": t_total_ms,
        "seed": sc.seed
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate Native Rust RhizAeon on Grand 100k Benchmark Battery")
    parser.add_argument("--samples-per-regime", type=int, default=500, help="Number of trials per regime (default: 500 = 5,000 total)")
    parser.add_argument("--workers", type=int, default=48, help="Number of worker processes (default: 48)")
    parser.add_argument("--db-path", type=Path, default=Path("rust_benchmark_5000.db"), help="SQLite output path")
    parser.add_argument("--csv-summary", type=Path, default=Path("rust_benchmark_summary.csv"), help="Summary CSV output path")
    parser.add_argument("--regimes", nargs="+", default=None, help="Optional subset of regimes to run")
    args = parser.parse_args()

    print("================================================================================")
    print("  EVALUATING NATIVE RUST RHIZAEON (PA-FF v2.1) BENCHMARK BATTERY")
    print(f"  Samples per regime: {args.samples_per_regime}")
    print(f"  Worker processes:   {args.workers}")
    print(f"  Database output:    {args.db_path}")
    print("================================================================================")

    # 1. Initialize SQLite database
    conn = sqlite3.connect(str(args.db_path), timeout=60.0)
    conn.executescript(DB_SCHEMA)
    conn.commit()

    # Find already completed indices
    cursor = conn.cursor()
    cursor.execute("SELECT global_idx FROM rust_benchmark_results")
    completed = {row[0] for row in cursor.fetchall()}
    print(f"  Found {len(completed)} already completed trials in {args.db_path}")

    # 2. Construct stratified indices
    selected_indices: List[int] = []
    for regime_name, start, end, is_null in REGIME_SLICES:
        if args.regimes and regime_name not in args.regimes:
            continue
        step = (end - start) / float(args.samples_per_regime)
        regime_idxs = [int(start + i * step) for i in range(args.samples_per_regime)]
        selected_indices.extend(regime_idxs)

    indices_to_run = [idx for idx in selected_indices if idx not in completed]
    print(f"  Total planned trials: {len(selected_indices)} | Remaining to run: {len(indices_to_run)}")

    if not indices_to_run:
        print("  All planned trials already completed!")
    else:
        t0 = time.time()
        batch_size = 50
        batch: List[Dict[str, Any]] = []
        completed_count = len(completed)

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(evaluate_single_trial, idx): idx for idx in indices_to_run}
            
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    res = future.result()
                    batch.append(res)
                    completed_count += 1
                except Exception as e:
                    print(f"[Error on trial {idx}]: {e}", flush=True)

                if len(batch) >= batch_size:
                    cols = list(batch[0].keys())
                    placeholders = ", ".join(["?" for _ in cols])
                    col_names = ", ".join(cols)
                    insert_sql = f"INSERT OR REPLACE INTO rust_benchmark_results ({col_names}) VALUES ({placeholders})"
                    rows = [[b[c] for c in cols] for b in batch]
                    conn.executemany(insert_sql, rows)
                    conn.commit()
                    batch = []

                    elapsed = time.time() - t0
                    speed = (completed_count - len(completed)) / elapsed if elapsed > 0 else 0
                    eta_sec = (len(indices_to_run) - (completed_count - len(completed))) / speed if speed > 0 else 0
                    print(f"  [Progress] {completed_count}/{len(selected_indices)} trials ({completed_count/len(selected_indices)*100:.1f}%) | {speed:.1f} trials/s | ETA: {eta_sec/60.0:.1f} min", flush=True)

        if batch:
            cols = list(batch[0].keys())
            placeholders = ", ".join(["?" for _ in cols])
            col_names = ", ".join(cols)
            insert_sql = f"INSERT OR REPLACE INTO rust_benchmark_results ({col_names}) VALUES ({placeholders})"
            rows = [[b[c] for c in cols] for b in batch]
            conn.executemany(insert_sql, rows)
            conn.commit()

        total_elapsed = time.time() - t0
        print(f"\n  Completed all simulations in {total_elapsed:.1f}s ({len(indices_to_run)/total_elapsed:.1f} trials/sec)")

    # 3. Generate Detailed Aggregations and Summary CSV
    print("\n================================================================================")
    print("  COMPUTING DETAILED BENCHMARK SUMMARY & METRICS")
    print("================================================================================")

    summary_query = """
    SELECT 
        regime,
        is_null,
        count(*) as total_trials,
        sum(detected) as detected_trials,
        round(avg(detected) * 100.0, 2) as power_or_fpr_pct,
        round(avg(correct_recombinant_exact) * 100.0, 2) as rec_exact_accuracy_pct,
        round(avg(correct_recombinant_sister) * 100.0, 2) as rec_sister_accuracy_pct,
        round(avg(correct_recombinant_in_set) * 100.0, 2) as rec_in_set_accuracy_pct,
        round(avg(correct_home) * 100.0, 2) as home_accuracy_pct,
        round(avg(correct_donor) * 100.0, 2) as donor_accuracy_pct,
        round(avg(correct_both_parents) * 100.0, 2) as joint_parent_accuracy_pct,
        sum(num_true_bps) as total_true_bps,
        sum(num_inferred_bps) as total_inferred_bps,
        sum(matched_bps_50) as matched_bps_50,
        sum(missed_bps_50) as missed_bps_50,
        round(cast(sum(missed_bps_50) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_50_pct,
        sum(matched_bps_100) as matched_bps_100,
        sum(missed_bps_100) as missed_bps_100,
        round(cast(sum(missed_bps_100) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_100_pct,
        round(avg(bp_error_mean), 2) as mean_bp_error_nt,
        round(avg(bp_all_within_35) * 100.0, 2) as bp_within_35_pct,
        round(avg(bp_all_within_50) * 100.0, 2) as bp_within_50_pct,
        round(avg(bp_all_within_100) * 100.0, 2) as bp_within_100_pct,
        round(avg(runtime_rust_ms), 2) as mean_rust_ms,
        round(avg(runtime_total_ms), 2) as mean_total_ms
    FROM rust_benchmark_results
    GROUP BY regime, is_null
    ORDER BY is_null DESC, min(global_idx) ASC;
    """

    cursor.execute(summary_query)
    columns = [d[0] for d in cursor.description]
    summary_rows = cursor.fetchall()

    import csv
    with open(args.csv_summary, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(summary_rows)

    print(f"  Summary written to: {args.csv_summary}")

    conn.close()

if __name__ == "__main__":
    main()
