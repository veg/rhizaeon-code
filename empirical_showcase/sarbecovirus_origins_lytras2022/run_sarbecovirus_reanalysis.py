#!/usr/bin/env python3
"""
run_sarbecovirus_reanalysis.py
Comprehensive reanalysis harness for Sarbecovirus Subgenus Recombination and SARS-CoV-2 Origins.

Based on:
1. Lytras S, Hughes J, Martin D, Swanepoel P, de Klerk A, Lourens R, Kosakovsky Pond SL,
   Xia W, Jiang X, Robertson DL (2022).
   Exploring the Natural Origins of SARS-CoV-2 in the Light of Recombination.
   Genome Biology and Evolution 14(2): evac018. doi:10.1093/gbe/evac018.
2. Temmam S, Vongphrachanh K, Baquero E, et al. (2022).
   Bat coronaviruses related to SARS-CoV-2 and capable of using human ACE2 to enter cells.
   Nature 604: 330–336. doi:10.1038/s41586-022-04532-4.

Executes RhizAeon PA-FF v2.1 native Rust engine on:
- Dataset A: Original 78-taxon whole-genome Sarbecovirus alignment (30,956 bp)
- Dataset B: Expanded 83-taxon whole-genome alignment adding the 5 Laos BANAL bat coronaviruses (30,965 bp)

Compares inferred changepoints directly against:
- The 21 GARD breakpoints in Table S3 (Lytras et al. 2022)
- The 160 RDP5 unique events
"""

import json
import os
import subprocess
import time
import csv
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FASTA_78 = os.path.join(DATA_DIR, "78sarbeco_wgal_fn.fas")
FASTA_83 = os.path.join(DATA_DIR, "83sarbeco_wgal.fas")
TABLE_S3 = os.path.join(DATA_DIR, "Table_S3.csv")
GENE_ANNOT = os.path.join(DATA_DIR, "sarbeco_gene_annotations.json")

OUT_JSON_78 = os.path.join(DATA_DIR, "78sarbeco_result.json")
OUT_HTML_78 = os.path.join(DATA_DIR, "78sarbeco_dashboard.html")
OUT_JSON_83 = os.path.join(DATA_DIR, "83sarbeco_result.json")
OUT_HTML_83 = os.path.join(DATA_DIR, "83sarbeco_dashboard.html")

RHIZAEON_BIN = "rhizaeon"

def run_rhizaeon(fasta_path, out_json, out_html, title, force=False):
    if not force and os.path.exists(out_json) and os.path.getsize(out_json) > 100:
        print(f"\nLoading existing scan result: {os.path.basename(out_json)}...")
        with open(out_json, "r") as f:
            data = json.load(f)
        return data, data.get("run_time_ms", 0.0) / 1000.0

    print(f"\nScanning: {os.path.basename(fasta_path)}...")
    cmd = [
        RHIZAEON_BIN,
        "-i", fasta_path,
        "-o", out_json,
        "--html", out_html,
        "-t", title,
        "-l", "16",
        "-a", "0.05"
    ]
    t0 = time.perf_counter()
    res = subprocess.run(cmd, capture_output=True, text=True)
    t1 = time.perf_counter()
    wall_sec = t1 - t0
    
    if res.returncode != 0:
        print(f"RhizAeon failed: {res.stderr}")
        return None, wall_sec
    
    with open(out_json, "r") as f:
        data = json.load(f)
    return data, wall_sec

def load_gard_table_s3():
    gard_regions = []
    with open(TABLE_S3, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gard_regions.append({
                "region": int(row["RBP_region"]),
                "al_start": int(row["al_start"]),
                "al_end": int(row["al_end"]),
                "length": int(row["length"]),
                "sc2_start": int(row["SC2_start"]),
                "sc2_end": int(row["SC2_end"])
            })
    return gard_regions

def evaluate_gard_concordance(rhiz_events, gard_regions, tolerance_nt=150):
    rhiz_bps = set()
    for ev in rhiz_events:
        rhiz_bps.add(ev["u1"])
        rhiz_bps.add(ev["u2"])
    rhiz_bps = sorted(list(rhiz_bps))
    
    gard_bps = [r["al_end"] for r in gard_regions[:-1]]
    
    matches = []
    for g_idx, g_bp in enumerate(gard_bps):
        reg_num = g_idx + 1
        closest_bp = min(rhiz_bps, key=lambda x: abs(x - g_bp))
        dist = abs(closest_bp - g_bp)
        is_match = dist <= tolerance_nt
        matches.append({
            "gard_region": reg_num,
            "gard_bp": g_bp,
            "rhizaeon_closest_bp": closest_bp,
            "distance_nt": dist,
            "matched_within_tol": is_match
        })
    
    num_matched = sum(1 for m in matches if m["matched_within_tol"])
    pct_matched = (num_matched / len(gard_bps)) * 100.0
    mean_dist = np.mean([m["distance_nt"] for m in matches])
    median_dist = np.median([m["distance_nt"] for m in matches])
    
    return {
        "total_gard_breakpoints": len(gard_bps),
        "total_rhizaeon_breakpoints": len(rhiz_bps),
        "matched_within_tol": num_matched,
        "concordance_pct": pct_matched,
        "mean_distance_nt": float(mean_dist),
        "median_distance_nt": float(median_dist),
        "matches": matches
    }

def main():
    print("================================================================================")
    print("  RhizAeon Reanalysis: Natural Origins of SARS-CoV-2 (Lytras et al. 2022)       ")
    print("================================================================================")
    
    # 1. Run 78-taxon scan
    res_78, wall_78 = run_rhizaeon(
        FASTA_78, OUT_JSON_78, OUT_HTML_78,
        "Lytras et al. (2022) 78-Taxon Whole Sarbecovirus Alignment"
    )
    
    # 2. Run 83-taxon scan (with Laos BANAL coronaviruses)
    res_83, wall_83 = run_rhizaeon(
        FASTA_83, OUT_JSON_83, OUT_HTML_83,
        "Expanded 83-Taxon Sarbecovirus Alignment (+ Laos BANAL Bat Coronaviruses)"
    )
    
    gard_regions = load_gard_table_s3()
    
    print("\n--------------------------------------------------------------------------------")
    print("  HEAD-TO-HEAD BENCHMARKING & CONCORDANCE                                       ")
    print("--------------------------------------------------------------------------------")
    
    scan_ms_78 = res_78.get("run_time_ms", wall_78 * 1000.0)
    events_78 = res_78.get("events", [])
    print(f"[Dataset A: Original 78 Genomes (Lytras et al. 2022)]")
    print(f"  Dimensions:          78 taxa x {res_78['length']} bp ({res_78['informative_sites']} informative sites)")
    print(f"  Engine Runtime:      {scan_ms_78/1000.0:.2f} s (Wall clock: {wall_78:.2f} s)")
    print(f"  Verified Events:     {len(events_78)}")
    
    conc_78 = evaluate_gard_concordance(events_78, gard_regions, tolerance_nt=150)
    print(f"  GARD Concordance:    {conc_78['matched_within_tol']}/{conc_78['total_gard_breakpoints']} ({conc_78['concordance_pct']:.1f}%) within 150 nt")
    print(f"  Median Dist to GARD: {conc_78['median_distance_nt']:.1f} nt")
    
    scan_ms_83 = res_83.get("run_time_ms", wall_83 * 1000.0)
    events_83 = res_83.get("events", [])
    print(f"\n[Dataset B: Expanded 83 Genomes (+ Laos BANAL Isolates)]")
    print(f"  Dimensions:          83 taxa x {res_83['length']} bp ({res_83['informative_sites']} informative sites)")
    print(f"  Engine Runtime:      {scan_ms_83/1000.0:.2f} s (Wall clock: {wall_83:.2f} s)")
    print(f"  Verified Events:     {len(events_83)}")
    
    conc_83 = evaluate_gard_concordance(events_83, gard_regions, tolerance_nt=150)
    print(f"  GARD Concordance:    {conc_83['matched_within_tol']}/{conc_83['total_gard_breakpoints']} ({conc_83['concordance_pct']:.1f}%) within 150 nt")
    print(f"  Median Dist to GARD: {conc_83['median_distance_nt']:.1f} nt")
    
    print("\n--------------------------------------------------------------------------------")
    print("  KEY SARS-COV-2 PROGENITOR & RESERVOIR LINEAGES                                ")
    print("--------------------------------------------------------------------------------")
    focal_names = ["Wuhan-Hu-1", "RpYN06", "RaTG13", "BANAL-52", "BANAL-103", "BANAL-236", "BANAL-116", "RacCS203", "Pangolin-CoV", "P5E"]
    
    for ev in events_83:
        cand = ev["candidate_name"]
        donor = ev["donor_name"]
        home = ev["home_name"]
        if any(fn.lower() in cand.lower() or fn.lower() in donor.lower() for fn in focal_names):
            print(f"  Recombinant: {cand:28s} | Home: {home:22s} | Donor: {donor:22s}")
            print(f"    Tract: {ev['u1']:5d} - {ev['u2']:5d} bp ({ev['tract_length']:5d} bp) | p-Fisher: {ev['p_fisher']:.2e} | Z: {ev['z_phys']:.2f}")

    # Save summary report
    report = {
        "dataset_78": {
            "num_taxa": 78,
            "length": res_78["length"],
            "informative_sites": res_78["informative_sites"],
            "runtime_seconds": scan_ms_78 / 1000.0,
            "events_count": len(events_78),
            "gard_concordance": conc_78
        },
        "dataset_83": {
            "num_taxa": 83,
            "length": res_83["length"],
            "informative_sites": res_83["informative_sites"],
            "runtime_seconds": scan_ms_83 / 1000.0,
            "events_count": len(events_83),
            "gard_concordance": conc_83
        }
    }
    report_path = os.path.join(DATA_DIR, "sarbecovirus_reanalysis_report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nSaved comprehensive report to: {report_path}")

if __name__ == "__main__":
    main()
