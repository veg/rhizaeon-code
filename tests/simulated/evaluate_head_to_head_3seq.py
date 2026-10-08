#!/usr/bin/env python3
"""
evaluate_rhizaeon_vs_3seq.py
============================
Rigorous head-to-head empirical comparison between:
  1. RhizAeon PA-FF v2.1 (Native Rust Engine)
  2. 3SEQ (Boni et al., 2007)

Evaluates 11 representative benchmark datasets spanning:
  - 4 Adversarial Null Regimes (Clock, Heterotachy, APOBEC, Darwinian Burst)
  - 2 Single Crossovers (Low & High Divergence)
  - 2 Cassette Swaps (Standard & Long-span)
  - 1 Micro-conversion (Short tract 120 nt)
  - 1 Ghost Donor Introgression (Unsampled donor lineage)
  - 1 Multi-Way Mosaicism (3 distinct parental lineages, 4 breakpoints)
"""

import sys
import os
import time
import json
import subprocess
import re
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional, Set
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

# Paths
BENCH_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/19_grand_100k_paff_v2")
if str(BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(BENCH_DIR))

from scenario_registry import get_scenario_for_index, MultiDimScenario
from fast_tree_simulator import FastTreeSimulator

RUST_BIN = "/Users/sergei/Projects/RhizAeon/rhizaeon-rust/target/release/rhizaeon"
THREE_SEQ_BIN = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"

TRIALS = [
    (28, "clock_null", "Clock Null"),
    (14050, "heterotachy_darren_null", "Heterotachy Darren Null"),
    (21050, "apobec_shower_null", "APOBEC Shower Null"),
    (28050, "darwinian_burst_null", "Darwinian Burst Null"),
    (35000, "single_crossover", "Single Crossover (Low Div)"),
    (35100, "single_crossover", "Single Crossover (High Div)"),
    (48010, "cassette_swap", "Cassette Swap (Standard)"),
    (48005, "cassette_swap", "Cassette Swap (Long Span)"),
    (63031, "micro_conversion", "Micro-conversion (120 nt)"),
    (76010, "ghost_donor_introgression", "Ghost Donor Introgression"),
    (88015, "multiway_mosaic", "Multi-Way Mosaic (3 Parents)"),
]

def parse_3seq_breakpoints(bp_str: str, seq_len: int) -> List[float]:
    """Parses 3SEQ breakpoint range string into representative point coordinates."""
    if not bp_str or bp_str.strip() == "":
        return []
    
    # 3SEQ can report multiple candidate pairs separated by commas
    # Take the first / best candidate pair
    primary_candidate = bp_str.split(",")[0].strip()
    parts = [p.strip() for p in primary_candidate.split("&") if p.strip()]
    
    parsed_ranges: List[Tuple[float, float]] = []
    for part in parts:
        m = re.match(r"(\d+)-(\d+)", part)
        if m:
            start_val = float(m.group(1))
            end_val = float(m.group(2))
            parsed_ranges.append((start_val, end_val))
            
    if not parsed_ranges:
        return []
        
    inferred_bps: List[float] = []
    if len(parsed_ranges) == 1:
        # Single range
        mid = (parsed_ranges[0][0] + parsed_ranges[0][1]) / 2.0
        inferred_bps.append(mid)
    elif len(parsed_ranges) >= 2:
        r1 = parsed_ranges[0]
        r2 = parsed_ranges[1]
        mid1 = (r1[0] + r1[1]) / 2.0
        mid2 = (r2[0] + r2[1]) / 2.0
        
        # Check boundary status
        is_r1_start = (r1[0] <= 100 or mid1 <= 100)
        is_r2_end = (r2[1] >= seq_len - 100 or mid2 >= seq_len - 100)
        
        if is_r1_start and not is_r2_end:
            # Single crossover starting at sequence beginning and ending at r2
            inferred_bps.append(mid2)
        elif not is_r1_start and is_r2_end:
            # Single crossover starting at r1 and ending at sequence end
            inferred_bps.append(mid1)
        elif is_r1_start and is_r2_end:
            # Spans whole sequence or degenerate
            inferred_bps.append(mid1)
            inferred_bps.append(mid2)
        else:
            # Internal cassette swap
            inferred_bps.append(mid1)
            inferred_bps.append(mid2)
            
    return sorted(inferred_bps)

def match_breakpoints(true_bps: List[int], inf_bps: List[float]) -> Tuple[float, int, int]:
    """Computes bipartite matching MAE, matched within 50 nt, and missed count."""
    if not true_bps:
        return (0.0, 0, 0)
    if not inf_bps:
        return (float("nan"), 0, len(true_bps))
        
    cost = np.abs(np.array(true_bps)[:, None] - np.array(inf_bps)[None, :])
    r_ind, c_ind = linear_sum_assignment(cost)
    errors = [float(cost[r, c]) for r, c in zip(r_ind, c_ind)]
    
    mae = float(np.mean(errors))
    matched_50 = sum(1 for e in errors if e <= 50.0)
    missed = len(true_bps) - matched_50
    return (mae, matched_50, missed)

def evaluate_trial(global_idx: int, label: str) -> Dict[str, Any]:
    sc = get_scenario_for_index(global_idx)
    sim = FastTreeSimulator(sc)
    seqs, meta = sim.simulate()
    
    fasta_path = Path(f"/tmp/bench_trial_{global_idx}.fasta")
    json_path = Path(f"/tmp/bench_trial_{global_idx}_rhiz.json")
    prefix_3seq = f"bench_trial_{global_idx}_3seq"
    csv_3seq = Path(f"/tmp/{prefix_3seq}.3s.rec.csv")
    
    with open(fasta_path, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")
            
    # Ground truth
    is_null = sc.is_null
    true_recs = sc.recombinant_taxa
    target_rec = true_recs[0] if true_recs else "-"
    expected_pars = sc.parent_taxa.get(target_rec, []) if true_recs else []
    is_true_ghost = any("Ghost" in p for p in expected_pars)
    expected_home = expected_pars[0] if len(expected_pars) > 0 else "-"
    expected_donor = expected_pars[1] if len(expected_pars) > 1 else ("Ghost" if is_true_ghost else "-")
    true_bps = sorted(sc.true_breakpoints)
    L = sc.length_nt
    N = sc.num_taxa
    
    # ----------------------------------------------------
    # 1. RUN RHIZAEON
    # ----------------------------------------------------
    cmd_rhiz = [RUST_BIN, "-i", str(fasta_path), "-o", str(json_path)]
    t0 = time.perf_counter()
    res_rhiz = subprocess.run(cmd_rhiz, capture_output=True, text=True)
    t_rhiz_ms = (time.perf_counter() - t0) * 1000.0
    
    rhiz_events = []
    if json_path.exists():
        try:
            with open(json_path) as f:
                rhiz_json = json.load(f)
                rhiz_events = rhiz_json.get("events", [])
        except Exception:
            pass
            
    rhiz_detected = len(rhiz_events) > 0
    rhiz_rec = rhiz_events[0].get("candidate_name", "-") if rhiz_detected else "-"
    rhiz_home = rhiz_events[0].get("home_name", "-") if rhiz_detected else "-"
    rhiz_is_ghost = rhiz_events[0].get("is_ghost_donor", False) if rhiz_detected else False
    rhiz_donor = ("Ghost" if rhiz_is_ghost else rhiz_events[0].get("donor_name", "-")) if rhiz_detected else "-"
    
    # Extract RhizAeon breakpoints
    rhiz_raw_bps: List[float] = []
    for ev in rhiz_events:
        is_cross = ev.get("is_crossover", False)
        u1_disc = ev.get("u1", 0)
        u2_disc = ev.get("u2", 0)
        u1_cont = ev.get("u1_continuous", float(u1_disc))
        u2_cont = ev.get("u2_continuous", float(u2_disc))
        
        if is_cross:
            if u1_disc <= 50:
                rhiz_raw_bps.append(u2_cont)
            elif u2_disc >= L - 50:
                rhiz_raw_bps.append(u1_cont)
            else:
                if u1_disc < (L - u2_disc):
                    rhiz_raw_bps.append(u2_cont)
                else:
                    rhiz_raw_bps.append(u1_cont)
        else:
            rhiz_raw_bps.append(u1_cont)
            rhiz_raw_bps.append(u2_cont)
            
    rhiz_bps: List[float] = []
    for b in sorted(rhiz_raw_bps):
        if not rhiz_bps or min(abs(b - ub) for ub in rhiz_bps) > 5.0:
            rhiz_bps.append(b)
            
    rhiz_mae, rhiz_matched_50, rhiz_missed = match_breakpoints(true_bps, rhiz_bps)
    
    rhiz_rec_correct = int(rhiz_rec == target_rec) if not is_null and rhiz_detected else 0
    rhiz_home_correct = int(rhiz_home == expected_home or (rhiz_rec == expected_home and rhiz_home == target_rec)) if not is_null and rhiz_detected else 0
    if is_true_ghost:
        rhiz_donor_correct = int(rhiz_is_ghost or rhiz_donor == "Ghost") if not is_null and rhiz_detected else 0
    else:
        rhiz_donor_correct = int(rhiz_donor in expected_pars and not rhiz_is_ghost and rhiz_donor != rhiz_home) if not is_null and rhiz_detected else 0
    rhiz_parents_correct = int(rhiz_home_correct and rhiz_donor_correct)
    
    # ----------------------------------------------------
    # 2. RUN 3SEQ
    # ----------------------------------------------------
    cmd_3seq = [THREE_SEQ_BIN, "-full", str(fasta_path), "-id", prefix_3seq]
    t0 = time.perf_counter()
    res_3seq = subprocess.run(cmd_3seq, cwd="/tmp", capture_output=True, text=True)
    t_3seq_ms = (time.perf_counter() - t0) * 1000.0
    
    seq3_triplets = []
    if csv_3seq.exists():
        try:
            lines = csv_3seq.read_text().strip().split("\n")
            if len(lines) > 1:
                header = [h.strip() for h in lines[0].split(",")]
                for line in lines[1:]:
                    if line.strip():
                        row = [c.strip() for c in line.split(",")]
                        if len(row) >= len(header):
                            seq3_triplets.append(dict(zip(header, row)))
        except Exception:
            pass
            
    seq3_detected = len(seq3_triplets) > 0
    seq3_triplet_count = len(seq3_triplets)
    
    # Check 3seq candidates
    # Find triplet involving target_rec, or take first triplet
    top_3seq_row = None
    if seq3_triplets:
        for row in seq3_triplets:
            if row.get("C_name") == target_rec:
                top_3seq_row = row
                break
        if not top_3seq_row:
            top_3seq_row = seq3_triplets[0]
            
    seq3_rec = top_3seq_row.get("C_name", "-") if top_3seq_row else "-"
    seq3_p1 = top_3seq_row.get("P_name", "-") if top_3seq_row else "-"
    seq3_p2 = top_3seq_row.get("Q_name", "-") if top_3seq_row else "-"
    
    # 3seq breakpoints
    seq3_bp_str = top_3seq_row.get("breakpoints", "") if top_3seq_row else ""
    seq3_bps = parse_3seq_breakpoints(seq3_bp_str, L)
    seq3_mae, seq3_matched_50, seq3_missed = match_breakpoints(true_bps, seq3_bps)
    
    seq3_rec_correct = int(seq3_rec == target_rec) if not is_null and seq3_detected else 0
    seq3_parents_in_set = 0
    if not is_null and seq3_detected:
        if is_true_ghost:
            # 3seq cannot call Ghost directly; check if at least one parent matches expected_home
            seq3_parents_in_set = int(seq3_p1 == expected_home or seq3_p2 == expected_home)
        else:
            # Both parents should match expected_pars
            seq3_parents_in_set = int((seq3_p1 in expected_pars) and (seq3_p2 in expected_pars) and (seq3_p1 != seq3_p2))
            
    # Cleanup tmp files
    for p in [fasta_path, json_path, csv_3seq, Path(f"/tmp/{prefix_3seq}.3s.log"), Path(f"/tmp/{prefix_3seq}.3s.pvalHist"), Path(f"/tmp/{prefix_3seq}.3s.longRec")]:
        if p.exists():
            try:
                p.unlink()
            except Exception:
                pass
                
    return {
        "global_idx": global_idx,
        "regime": sc.regime,
        "label": label,
        "num_taxa": N,
        "length_nt": L,
        "divergence": sc.mean_divergence,
        "is_null": int(is_null),
        "true_rec": target_rec,
        "true_home": expected_home,
        "true_donor": expected_donor,
        "true_bps": ";".join(str(b) for b in true_bps),
        # RhizAeon
        "rhiz_detected": int(rhiz_detected),
        "rhiz_rec": rhiz_rec,
        "rhiz_home": rhiz_home,
        "rhiz_donor": rhiz_donor,
        "rhiz_bps": ";".join(f"{b:0.1f}" for b in rhiz_bps),
        "rhiz_rec_correct": rhiz_rec_correct,
        "rhiz_parents_correct": rhiz_parents_correct,
        "rhiz_bp_mae": rhiz_mae,
        "rhiz_matched_50": rhiz_matched_50,
        "rhiz_missed_bps": rhiz_missed,
        "rhiz_time_ms": t_rhiz_ms,
        # 3SEQ
        "seq3_detected": int(seq3_detected),
        "seq3_triplets": seq3_triplet_count,
        "seq3_rec": seq3_rec,
        "seq3_p1": seq3_p1,
        "seq3_p2": seq3_p2,
        "seq3_bps": ";".join(f"{b:0.1f}" for b in seq3_bps),
        "seq3_rec_correct": seq3_rec_correct,
        "seq3_parents_correct": seq3_parents_in_set,
        "seq3_bp_mae": seq3_mae,
        "seq3_matched_50": seq3_matched_50,
        "seq3_missed_bps": seq3_missed,
        "seq3_time_ms": t_3seq_ms,
        "speedup": t_3seq_ms / max(0.1, t_rhiz_ms),
    }

def main():
    print("=" * 80)
    print("  RHIZAEON (PA-FF v2.1) vs 3SEQ: HEAD-TO-HEAD BENCHMARK EVALUATION")
    print("=" * 80)
    
    results = []
    for idx, regime, label in TRIALS:
        print(f"\n---> Evaluating Trial {idx}: {label} ({regime}) ...", flush=True)
        res = evaluate_trial(idx, label)
        results.append(res)
        
        print(f"     Status: True Null={bool(res['is_null'])}, N={res['num_taxa']}, L={res['length_nt']}, d={res['divergence']:0.4f}")
        print(f"     RhizAeon: Called H1={bool(res['rhiz_detected'])}, Rec={res['rhiz_rec']}, Home={res['rhiz_home']}, Donor={res['rhiz_donor']}, BPs={res['rhiz_bps']}, MAE={res['rhiz_bp_mae']:0.1f}nt, Time={res['rhiz_time_ms']:0.1f}ms")
        print(f"     3SEQ:     Called H1={bool(res['seq3_detected'])}, Triplets={res['seq3_triplets']}, Rec={res['seq3_rec']}, P1={res['seq3_p1']}, P2={res['seq3_p2']}, BPs={res['seq3_bps']}, MAE={res['seq3_bp_mae']:0.1f}nt, Time={res['seq3_time_ms']:0.1f}ms")
        print(f"     Speedup:  RhizAeon is {res['speedup']:0.2f}x faster than 3SEQ")

    df = pd.DataFrame(results)
    out_csv = Path("/Users/sergei/Projects/TOGA_MEME/recombination/viz/rhizaeon_vs_3seq_head_to_head.csv")
    df.to_csv(out_csv, index=False)
    print(f"\n[DONE] Saved summary CSV to: {out_csv}")

if __name__ == "__main__":
    main()
