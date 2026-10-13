#!/usr/bin/env python3
"""
evaluate_darren_arg.py
Evaluates competing evolutionary hypotheses and Ancestral Recombination Graph (ARG)
on Darren Martin's 61-taxon Sarbecovirus alignment.
"""

import os
import sys
import json
import time
import csv
import numpy as np
from Bio import AlignIO

PROJECT_ROOT = "/Users/sergei/Projects/rhizaeon-code"
sys.path.insert(0, PROJECT_ROOT)

from tools.history_validation.history_engine import (
    ReticulationEvent,
    AncestralConsolidator,
    RecombinationHistory,
    HistoryScorer
)

DARREN_FASTA = os.path.join(PROJECT_ROOT, "empirical_datasets/darren_sarbeco/Sarbeco alignment.fas")
DARREN_JSON = os.path.join(PROJECT_ROOT, "empirical_datasets/darren_sarbeco/darren_sarbeco_result.json")
TABLE_S3 = os.path.join(PROJECT_ROOT, "empirical_datasets/01_sarbecovirus_origins_lytras2022/Table_S3.csv")

def run():
    print("Loading alignment...")
    aln = AlignIO.read(DARREN_FASTA, "fasta")
    N = len(aln)
    L = aln.get_alignment_length()
    taxa = [s.id for s in aln]
    print(f"Alignment loaded: {N} taxa, {L} nt")

    # Map Wuhan-Hu-1 reference coordinates to Darren alignment coordinates
    d_wuhan = [s for s in aln if 'wuhan' in s.id.lower()][0]
    w_seq = str(d_wuhan.seq).upper()
    wuhan_to_aln = {}
    w_idx = 0
    for aln_idx, ch in enumerate(w_seq, 1):
        if ch != '-':
            w_idx += 1
            wuhan_to_aln[w_idx] = aln_idx

    # Load GARD breakpoints from Table S3
    gard_breaks = [0]
    with open(TABLE_S3) as f:
        reader = csv.DictReader(f)
        for r in reader:
            sc2_end = int(r["SC2_end"])
            aln_pos = wuhan_to_aln.get(sc2_end, sc2_end)
            if aln_pos < L:
                gard_breaks.append(aln_pos)
    gard_breaks.append(L)
    gard_breaks = sorted(list(set(gard_breaks)))
    gard_partitions = [(gard_breaks[i], gard_breaks[i+1]) for i in range(len(gard_breaks)-1)]

    # Load RhizAeon events and consolidate
    with open(DARREN_JSON) as f:
        d_res = json.load(f)
    raw_events = d_res.get("events", [])
    consolidator = AncestralConsolidator(jaccard_thresh=0.60, max_boundary_dist=100)
    cons_events = consolidator.consolidate(raw_events)

    # RhizAeon primary chromosomal blocks
    # Breakpoints at major structural junctions in Darren's alignment:
    # 0, 3804 (nsp3), 8387 (nsp4/5), 13036 (RdRp entry), 20331 (Spike upstream),
    # 21866 (Spike start), 24753 (S1/S2 furin), 25925 (Spike stop), 28960 (ORF8/N), L
    rhiz_breakpoints = [0, 3804, 8387, 13036, 20331, 21866, 24753, 25925, 28960, L]
    rhiz_partitions = [(rhiz_breakpoints[i], rhiz_breakpoints[i+1]) for i in range(len(rhiz_breakpoints)-1)]

    # Construct hypotheses
    hypotheses = []

    # H0: Clonal Null
    h0 = RecombinationHistory(
        history_id="H0_Clonal_Null",
        name="Genome-Wide Clonal Null",
        description="Single unpartitioned phylogenetic tree across the entire 31,214 nt genome.",
        partitions=[(0, L)],
        reticulations=[]
    )
    hypotheses.append(h0)

    # H_GARD: 22 unpolarized partitions
    h_gard = RecombinationHistory(
        history_id="H_GARD_Partitions",
        name="GARD 22-Partition Model",
        description="22 unpolarized chromosomal partitions based on Lytras et al. Table S3 mapped to Darren's alignment.",
        partitions=gard_partitions,
        reticulations=[]
    )
    hypotheses.append(h_gard)

    # H_RhizAeon: Consolidated ARG
    h_rhiz = RecombinationHistory(
        history_id="H_RhizAeon_ARG",
        name="RhizAeon Consolidated ARG",
        description="Directed ARG: 9 consolidated chromosomal blocks coupled to 41 ancestral reticulations.",
        partitions=rhiz_partitions,
        reticulations=cons_events,
        k_params=len(rhiz_partitions) * (2 * N - 3) + len(cons_events) * 2 + 4
    )
    hypotheses.append(h_rhiz)

    # Evaluate
    scorer = HistoryScorer(aln)
    results = []
    print("\nEvaluating hypotheses...")
    for h in hypotheses:
        t0 = time.time()
        print(f"  Evaluating {h.history_id} ({h.name})...")
        scored = scorer.score_history(h)
        scored["evaluation_time_sec"] = round(time.time() - t0, 2)
        results.append(scored)

    min_aicc = min(r["phylogenetics"]["aicc"] for r in results)
    min_bic = min(r["phylogenetics"]["bic"] for r in results)
    min_strain = min(r["manifold_strain"]["residual_strain_ratio"] for r in results)

    for r in results:
        r["phylogenetics"]["delta_aicc"] = float(r["phylogenetics"]["aicc"] - min_aicc)
        r["phylogenetics"]["delta_bic"] = float(r["phylogenetics"]["bic"] - min_bic)
        r["manifold_strain"]["strain_excess"] = float(r["manifold_strain"]["residual_strain_ratio"] - min_strain)

    print("\n" + "="*115)
    print(f"{'Model ID':<25} | {'ln L':<12} | {'k':<6} | {'AICc':<12} | {'Delta BIC':<12} | {'Strain %':<9} | {'Explained %':<11} | {'Spectral SNR':<12}")
    print("="*115)
    for r in results:
        phy = r["phylogenetics"]
        st = r["manifold_strain"]
        print(f"{r['history_id']:<25} | {phy['total_lnL']:<12.1f} | {phy['k_parameters']:<6d} | {phy['aicc']:<12.1f} | {phy['delta_bic']:<12.1f} | {st['residual_strain_ratio']*100:<9.2f} | {st['strain_explained_pct']:<11.2f} | {st['spectral_snr']:<12.2f}")
    print("="*115)

    sscr = results[2]["synapomorphic_concordance"]
    print(f"\nRhizAeon ARG SSCR Validation:")
    print(f"  Reticulations evaluated: {sscr['events_evaluated']}")
    print(f"  Mean in-tract SSCR:      {sscr['mean_sscr_in']:.3f}")
    print(f"  Mean flanking SSCR:      {sscr['mean_sscr_flank']:.3f}")
    print(f"  Mean Delta SSCR:         {sscr['mean_delta_sscr']:.3f}")

    out_file = os.path.join(PROJECT_ROOT, "empirical_datasets/darren_sarbeco/darren_history_validation_results.json")
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved complete validation results to: {out_file}")

if __name__ == "__main__":
    run()
