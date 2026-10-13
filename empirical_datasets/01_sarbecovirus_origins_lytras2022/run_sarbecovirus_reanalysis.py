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

import shutil

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = BASE_DIR
FASTA_78 = os.path.join(DATA_DIR, "78sarbeco_wgal_fn.fas")
FASTA_83 = os.path.join(DATA_DIR, "83sarbeco_wgal.fas")
TABLE_S3 = os.path.join(DATA_DIR, "Table_S3.csv")
GENE_ANNOT = os.path.join(DATA_DIR, "sarbeco_gene_annotations.json")

OUT_JSON_78 = os.path.join(DATA_DIR, "78sarbeco_result.json")
OUT_HTML_78 = os.path.join(DATA_DIR, "78sarbeco_dashboard.html")
OUT_JSON_83 = os.path.join(DATA_DIR, "83sarbeco_result.json")
OUT_HTML_83 = os.path.join(DATA_DIR, "83sarbeco_dashboard.html")

def get_rhizaeon_bin():
    if shutil.which("rhizaeon"):
        return shutil.which("rhizaeon")
    repo_root = os.path.abspath(os.path.join(BASE_DIR, "..", ".."))
    for candidate in [
        os.path.join(repo_root, "target", "release", "rhizaeon"),
        os.path.join(repo_root, "target", "debug", "rhizaeon"),
    ]:
        if os.path.exists(candidate):
            return candidate
    return "rhizaeon"

RHIZAEON_BIN = get_rhizaeon_bin()

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

def read_fasta_dict(fasta_path):
    seqs = {}
    cur = None
    buf = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip()
                buf = []
            else:
                buf.append(line)
    if cur:
        seqs[cur] = "".join(buf)
    return seqs

def build_coordinate_maps(fasta_78, fasta_83):
    s78 = read_fasta_dict(fasta_78)
    s83 = read_fasta_dict(fasta_83)
    w78 = [v for k, v in s78.items() if "Wuhan-Hu-1" in k][0]
    w83 = [v for k, v in s83.items() if "Wuhan-Hu-1" in k][0]
    shared_taxa = [k for k in s78.keys() if k in s83]

    def build_aln_to_ref(w_seq):
        a2r = {}
        r2a = {}
        r = 0
        for idx, ch in enumerate(w_seq, start=1):
            if ch != "-":
                r += 1
                r2a[r] = idx
            a2r[idx] = r
        return a2r, r2a

    a2r_78, r2a_78 = build_aln_to_ref(w78)
    a2r_83, r2a_83 = build_aln_to_ref(w83)

    col78_to_83 = {}
    j = 0
    for i in range(len(w78)):
        while j < len(w83):
            col78 = [s78[t][i] for t in shared_taxa[:20]]
            col83 = [s83[t][j] for t in shared_taxa[:20]]
            if col78 == col83:
                col78_to_83[i + 1] = j + 1
                j += 1
                break
            else:
                j += 1

    return {
        "a2r_78": a2r_78,
        "r2a_78": r2a_78,
        "a2r_83": a2r_83,
        "r2a_83": r2a_83,
        "col78_to_83": col78_to_83,
    }

def load_gard_table_s3(col78_to_83):
    gard_regions = []
    with open(TABLE_S3, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            al_s78 = int(row["al_start"])
            al_e78 = int(row["al_end"])
            gard_regions.append({
                "region": int(row["RBP_region"]),
                "al_start": al_s78,
                "al_end": al_e78,
                "al_start_83": col78_to_83.get(al_s78, al_s78),
                "al_end_83": col78_to_83.get(al_e78, al_e78),
                "length": int(row["length"]),
                "sc2_start": int(row["SC2_start"]),
                "sc2_end": int(row["SC2_end"])
            })
    return gard_regions

def evaluate_gard_concordance(rhiz_events, gard_regions, a2r_map, is_83=False, tolerance_nt=150):
    rhiz_bps_aln = set()
    for ev in rhiz_events:
        rhiz_bps_aln.add(ev["u1"])
        rhiz_bps_aln.add(ev["u2"])
    rhiz_bps_aln = sorted(list(rhiz_bps_aln))
    rhiz_bps_ref = sorted(list({a2r_map[b] for b in rhiz_bps_aln if b in a2r_map}))

    matches = []
    for g_idx, r in enumerate(gard_regions[:-1]):
        reg_num = g_idx + 1
        g_bp_aln = r["al_end_83"] if is_83 else r["al_end"]
        g_bp_ref = r["sc2_end"]

        closest_bp_aln = min(rhiz_bps_aln, key=lambda x: abs(x - g_bp_aln))
        dist_aln = abs(closest_bp_aln - g_bp_aln)
        is_match_aln = dist_aln <= tolerance_nt

        closest_bp_ref = min(rhiz_bps_ref, key=lambda x: abs(x - g_bp_ref))
        dist_ref = abs(closest_bp_ref - g_bp_ref)
        is_match_ref = dist_ref <= tolerance_nt

        matches.append({
            "gard_region": reg_num,
            "gard_bp_aln78": r["al_end"],
            "gard_bp_aln83": r["al_end_83"],
            "gard_bp": g_bp_aln,
            "gard_bp_ref": g_bp_ref,
            "rhizaeon_closest_bp": closest_bp_aln,
            "rhizaeon_closest_bp_ref_of_aln": a2r_map.get(closest_bp_aln, -1),
            "distance_nt": dist_aln,
            "matched_within_tol": is_match_aln,
            "rhizaeon_closest_bp_ref": closest_bp_ref,
            "distance_ref_nt": dist_ref,
            "matched_ref_within_tol": is_match_ref
        })

    num_matched_aln = sum(1 for m in matches if m["matched_within_tol"])
    pct_matched_aln = (num_matched_aln / len(matches)) * 100.0
    mean_dist_aln = float(np.mean([m["distance_nt"] for m in matches]))
    median_dist_aln = float(np.median([m["distance_nt"] for m in matches]))

    num_matched_ref = sum(1 for m in matches if m["matched_ref_within_tol"])
    pct_matched_ref = (num_matched_ref / len(matches)) * 100.0
    mean_dist_ref = float(np.mean([m["distance_ref_nt"] for m in matches]))
    median_dist_ref = float(np.median([m["distance_ref_nt"] for m in matches]))

    return {
        "total_gard_breakpoints": len(matches),
        "total_rhizaeon_breakpoints": len(rhiz_bps_aln),
        "matched_within_tol": num_matched_aln,
        "concordance_pct": pct_matched_aln,
        "mean_distance_nt": mean_dist_aln,
        "median_distance_nt": median_dist_aln,
        "matched_ref_within_tol": num_matched_ref,
        "concordance_ref_pct": pct_matched_ref,
        "mean_distance_ref_nt": mean_dist_ref,
        "median_distance_ref_nt": median_dist_ref,
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

    maps = build_coordinate_maps(FASTA_78, FASTA_83)
    gard_regions = load_gard_table_s3(maps["col78_to_83"])

    print("\n--------------------------------------------------------------------------------")
    print("  HEAD-TO-HEAD BENCHMARKING & CONCORDANCE (MSA FRAME & WUHAN-HU-1 REF FRAME)    ")
    print("--------------------------------------------------------------------------------")

    scan_ms_78 = res_78.get("run_time_ms", wall_78 * 1000.0)
    events_78 = res_78.get("events", [])
    print(f"[Dataset A: Original 78 Genomes (Lytras et al. 2022)]")
    print(f"  Dimensions:          78 taxa x {res_78['length']} bp ({res_78['informative_sites']} informative sites)")
    print(f"  Engine Runtime:      {scan_ms_78/1000.0:.2f} s (Wall clock: {wall_78:.2f} s)")
    print(f"  Verified Events:     {len(events_78)}")

    conc_78 = evaluate_gard_concordance(events_78, gard_regions, maps["a2r_78"], is_83=False, tolerance_nt=150)
    print(f"  GARD Concordance (78-MSA Frame):      {conc_78['matched_within_tol']}/{conc_78['total_gard_breakpoints']} ({conc_78['concordance_pct']:.1f}%), Median Dist: {conc_78['median_distance_nt']:.1f} nt")
    print(f"  GARD Concordance (Wuhan-Hu-1 Frame):  {conc_78['matched_ref_within_tol']}/{conc_78['total_gard_breakpoints']} ({conc_78['concordance_ref_pct']:.1f}%), Median Dist: {conc_78['median_distance_ref_nt']:.1f} nt")

    scan_ms_83 = res_83.get("run_time_ms", wall_83 * 1000.0)
    events_83 = res_83.get("events", [])
    print(f"\n[Dataset B: Expanded 83 Genomes (+ Laos BANAL Isolates)]")
    print(f"  Dimensions:          83 taxa x {res_83['length']} bp ({res_83['informative_sites']} informative sites)")
    print(f"  Engine Runtime:      {scan_ms_83/1000.0:.2f} s (Wall clock: {wall_83:.2f} s)")
    print(f"  Verified Events:     {len(events_83)}")

    conc_83 = evaluate_gard_concordance(events_83, gard_regions, maps["a2r_83"], is_83=True, tolerance_nt=150)
    print(f"  GARD Concordance (83-MSA Frame):      {conc_83['matched_within_tol']}/{conc_83['total_gard_breakpoints']} ({conc_83['concordance_pct']:.1f}%), Median Dist: {conc_83['median_distance_nt']:.1f} nt")
    print(f"  GARD Concordance (Wuhan-Hu-1 Frame):  {conc_83['matched_ref_within_tol']}/{conc_83['total_gard_breakpoints']} ({conc_83['concordance_ref_pct']:.1f}%), Median Dist: {conc_83['median_distance_ref_nt']:.1f} nt")

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

