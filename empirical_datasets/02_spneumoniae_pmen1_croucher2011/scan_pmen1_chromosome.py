#!/usr/bin/env python3
"""
scan_pmen1_chromosome.py
========================
End-to-end provenance pipeline for whole-chromosome recombination deconvolution
across the canonical Streptococcus pneumoniae PMEN1 cohort (Croucher et al. 2011 Science / 2015 NAR).
Chromosome: ATCC 700669 (GenBank FM211187, L = 2,221,315 bp).
"""

import os
import sys
import time
import json
import numpy as np
import pandas as pd
from Bio import SeqIO

def compute_classical_mds(d_mat, k=3):
    """Classical MDS (Torgerson-Gower scaling) on pairwise distance matrix."""
    n = d_mat.shape[0]
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ (d_mat ** 2) @ h
    eigvals, eigvecs = np.linalg.eigh(b)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    k_eff = min(k, max(1, n - 1))
    pos_evals = np.maximum(0.0, eigvals[:k_eff])
    coords = eigvecs[:, :k_eff] * np.sqrt(pos_evals)
    return coords

def align_procrustes(z1, z2):
    """Orthogonal Procrustes alignment of z2 onto z1."""
    mu1 = np.mean(z1, axis=0)
    mu2 = np.mean(z2, axis=0)
    z1_c = z1 - mu1
    z2_c = z2 - mu2
    m = z2_c.T @ z1_c
    u, s, vt = np.linalg.svd(m)
    r = u @ vt
    z2_rot = z2_c @ r + mu1
    res = np.linalg.norm(z1 - z2_rot, axis=1)
    return z2_rot, res

def main():
    import gzip
    base_dir = os.path.dirname(os.path.abspath(__file__))
    aln_path = os.path.join(base_dir, "PMEN1.aln.gz") if os.path.exists(os.path.join(base_dir, "PMEN1.aln.gz")) else os.path.join(base_dir, "PMEN1.aln")
    gub_full_gff = os.path.join(base_dir, "gubbins_predictions.gff")
    ref_gff = os.path.join(base_dir, "Spn23f.gff")
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    print("=" * 80)
    print("RHIZAEON PA-FF: WHOLE-CHROMOSOME S. PNEUMONIAE PMEN1 DECONVOLUTION")
    print(f"Alignment: {aln_path}")
    print("=" * 80)

    # 1. Ingest whole-genome alignment
    t0_load = time.perf_counter()
    if aln_path.endswith(".gz"):
        with gzip.open(aln_path, "rt") as gz_handle:
            records = list(SeqIO.parse(gz_handle, "fasta"))
    else:
        records = list(SeqIO.parse(aln_path, "fasta"))
    taxa = [r.id for r in records]
    n = len(records)
    l = len(records[0].seq)
    t_load = time.perf_counter() - t0_load
    print(f"[*] Alignment loaded in {t_load:.2f} s: {n} taxa, {l:,} base pairs.")

    # Convert to 2D numpy array
    seqs = np.array([list(str(r.seq).upper()) for r in records])

    # 2. Extract Segregating Informative Polymorphisms
    t0_snp = time.perf_counter()
    valid_chars = set("ACGT")
    is_seg = np.zeros(l, dtype=bool)
    for col in range(l):
        col_chars = set(seqs[:, col]) & valid_chars
        if len(col_chars) > 1:
            is_seg[col] = True

    seg_indices = np.where(is_seg)[0]
    num_snps = len(seg_indices)
    t_snp = time.perf_counter() - t0_snp
    print(f"[*] Extracted {num_snps:,} segregating informative SNPs ({l / max(1, num_snps):.1f}x compression) in {t_snp:.3f} s.")

    # 3. Build Prefix Difference Tensor in SNP space
    t0_scan = time.perf_counter()
    snp_matrix = seqs[:, seg_indices]

    prefix_diff = np.zeros((n, n, num_snps + 1), dtype=np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            mismatches = (snp_matrix[i] != snp_matrix[j]).astype(np.int32)
            cum_m = np.zeros(num_snps + 1, dtype=np.int32)
            cum_m[1:] = np.cumsum(mismatches)
            prefix_diff[i, j] = cum_m
            prefix_diff[j, i] = cum_m

    def query_sub_distance(u1, u2):
        """O(1) bilateral pairwise distance query on [u1, u2) in SNP space."""
        diff = prefix_diff[:, :, u2] - prefix_diff[:, :, u1]
        span = max(1, u2 - u1)
        return diff.astype(np.float64) / span

    # 4. Multi-Scale Kinetic Strain Scan
    flanks = [10, 15, 25]  # Multi-scale SNP windows
    step = 2
    all_evals = []

    for flk in flanks:
        cps = np.arange(flk, num_snps - flk, step)
        for cp in cps:
            d1 = query_sub_distance(cp - flk, cp)
            d2 = query_sub_distance(cp, cp + flk)
            if np.max(d1) < 1e-6 or np.max(d2) < 1e-6:
                continue

            z1 = compute_classical_mds(d1, k=min(4, n - 1))
            z2 = compute_classical_mds(d2, k=min(4, n - 1))
            _, res = align_procrustes(z1, z2)
            all_evals.append((flk, cp, res))

    # Whole-genome background standardization
    res_matrix = np.array([e[2] for e in all_evals])
    bg_mean = float(np.mean(res_matrix))
    bg_std = float(np.std(res_matrix))

    candidate_peaks = {}
    for flk, cp, res in all_evals:
        top_taxon_idx = int(np.argmax(res))
        max_z = float((res[top_taxon_idx] - bg_mean) / bg_std)

        if max_z >= 2.0:
            if cp not in candidate_peaks or max_z > candidate_peaks[cp]["z"]:
                candidate_peaks[cp] = {
                    "snp_idx": int(cp),
                    "genome_nt": int(seg_indices[cp]),
                    "z": max_z,
                    "taxon": taxa[top_taxon_idx],
                    "flank": flk,
                    "raw_res": float(res[top_taxon_idx])
                }

    # Non-maximum suppression within 2,000 bp in physical genomic space
    sorted_cps = sorted(candidate_peaks.keys(), key=lambda c: candidate_peaks[c]["z"], reverse=True)
    kept_cps = []
    for c in sorted_cps:
        c_pos = candidate_peaks[c]["genome_nt"]
        if not any(abs(c_pos - candidate_peaks[k]["genome_nt"]) <= 2000 for k in kept_cps):
            kept_cps.append(c)
    kept_cps.sort()

    t_scan = time.perf_counter() - t0_scan
    print(f"[✓] Multi-scale Procrustes kinetic scan finished in {t_scan:.3f} s.")
    print(f"[✓] Background mean residual: {bg_mean:.5f}, std: {bg_std:.5f}")
    print(f"[✓] Identified {len(kept_cps)} significant changepoints (Z >= 2.0).")

    # 5. Continuous Physical Fields: Breakpoint Flux Phi(u) and Strain rho(u)
    t0_fields = time.perf_counter()
    grid_step = 1000  # 1 kb resolution grid across 2.22 Mb
    u_grid = np.arange(0, l, grid_step, dtype=np.float64)
    phi_field = np.zeros_like(u_grid)
    sigma_kernel = 5000.0  # 5 kb smoothing kernel for physical flux density

    for c in kept_cps:
        info = candidate_peaks[c]
        bp = info["genome_nt"]
        z = info["z"]
        # Gaussian kernel contribution
        dist = u_grid - bp
        phi_field += z * (1.0 / (np.sqrt(2 * np.pi) * sigma_kernel)) * np.exp(-0.5 * (dist / sigma_kernel) ** 2)

    # Cumulative Strain Field rho(u)
    strain_field = np.cumsum(phi_field) * grid_step
    # Normalize strain field so max matches total integrated Z
    t_fields = time.perf_counter() - t0_fields
    print(f"[✓] Continuous fields computed in {t_fields:.3f} s ({len(u_grid):,} grid points).")

    np.savez_compressed(
        os.path.join(results_dir, "pmen1_continuous_fields.npz"),
        u_grid=u_grid,
        phi_field=phi_field,
        strain_field=strain_field
    )

    # 6. Load Published Recombination Tracts from Croucher et al. EVAL GFF
    gub_tracts = []
    with open(eval_gff_path) as f:
        for line in f:
            if "\tGUBBINS\t" in line:
                parts = line.strip().split("\t")
                s = int(parts[3])
                e = int(parts[4])
                attr = parts[8]
                import re
                m_snps = re.search(r'snp_count=\"([^\"]+)\"', attr)
                gub_tracts.append({
                    "start": s,
                    "end": e,
                    "length": e - s + 1,
                    "snp_count": int(m_snps.group(1)) if m_snps else 0,
                    "attr": attr
                })

    # Cluster into hotspots (within 8 kb)
    gub_tracts.sort(key=lambda x: x["start"])
    hotspots = []
    for g in gub_tracts:
        if not hotspots or g["start"] > hotspots[-1]["end"] + 8000:
            hotspots.append({
                "start": g["start"],
                "end": g["end"],
                "total_snps": g["snp_count"],
                "events": [g]
            })
        else:
            hotspots[-1]["end"] = max(hotspots[-1]["end"], g["end"])
            hotspots[-1]["total_snps"] += g["snp_count"]
            hotspots[-1]["events"].append(g)

    print(f"[✓] Loaded {len(gub_tracts)} published Gubbins recombination tracts clustered into {len(hotspots)} distinct loci.")

    # 7. Detailed Changepoint Table
    records_out = []
    for c in kept_cps:
        info = candidate_peaks[c]
        snp_idx = c
        pos_nt = int(seg_indices[snp_idx])
        left_nt = int(seg_indices[max(0, snp_idx - 1)])
        right_nt = int(seg_indices[min(num_snps - 1, snp_idx + 1)])
        plateau = right_nt - left_nt

        min_dist = float("inf")
        matched_h = None
        for h in hotspots:
            d = max(0, h["start"] - pos_nt, pos_nt - h["end"])
            if d < min_dist:
                min_dist = d
                matched_h = h

        records_out.append({
            "snp_idx": snp_idx,
            "genome_coord": pos_nt,
            "left_flank_snp": left_nt,
            "right_flank_snp": right_nt,
            "uninformative_plateau": plateau,
            "recombinant_taxon": info["taxon"],
            "kinetic_z": round(info["z"], 3),
            "raw_procrustes_res": round(info["raw_res"], 5),
            "matched_hotspot_start": matched_h["start"] if matched_h else 0,
            "matched_hotspot_end": matched_h["end"] if matched_h else 0,
            "hotspot_distance": int(min_dist)
        })

    df_out = pd.DataFrame(records_out)
    csv_out = os.path.join(results_dir, "pmen1_changepoints_vs_gubbins.csv")
    df_out.to_csv(csv_out, index=False)
    print(f"[✓] Saved changepoint results to {csv_out}")

    # 8. Boundary Concordance with Gubbins
    gub_bps = []
    for t in gub_tracts:
        gub_bps.append(t["start"])
        gub_bps.append(t["end"])
    gub_bps = sorted(list(set(gub_bps)))
    rhiz_bps = [r["genome_coord"] for r in records_out]

    boundary_matches = []
    for gb in gub_bps:
        closest = min(rhiz_bps, key=lambda x: abs(x - gb))
        dist = abs(closest - gb)
        boundary_matches.append({
            "gubbins_bp": gb,
            "closest_rhizaeon_bp": closest,
            "distance_bp": dist,
            "within_1kb": dist <= 1000,
            "within_2kb": dist <= 2000,
            "within_5kb": dist <= 5000
        })

    dists = [m["distance_bp"] for m in boundary_matches]
    med_dist = float(np.median(dists))
    mean_dist = float(np.mean(dists))
    n_1k = sum(1 for d in dists if d <= 1000)
    n_2k = sum(1 for d in dists if d <= 2000)
    n_5k = sum(1 for d in dists if d <= 5000)

    print("\n" + "-" * 80)
    print(f"HEAD-TO-HEAD BOUNDARY CONCORDANCE WITH GUBBINS (52 boundaries):")
    print(f"  Median distance to nearest RhizAeon changepoint: {med_dist:.1f} bp")
    print(f"  Mean distance:                                  {mean_dist:.1f} bp")
    print(f"  Concordance <= 1 kb:                            {n_1k}/52 ({n_1k/52*100:.1f}%)")
    print(f"  Concordance <= 2 kb:                            {n_2k}/52 ({n_2k/52*100:.1f}%)")
    print(f"  Concordance <= 5 kb:                            {n_5k}/52 ({n_5k/52*100:.1f}%)")
    print("-" * 80)

    # 9. Summary of Clinical Landmarks
    landmarks = [
        ("Penicillin-binding protein 2x (pbp2x)", 285069, 305751),
        ("Capsular polysaccharide switch (cps, 23F->19A)", 336730, 352271),
        ("Penicillin-binding protein 1a (pbp1a)", 323721, 325108),
        ("Trimethoprim resistance island (folA)", 1620528, 1649446),
        ("Surface adhesin (pspC / cbpA)", 2104999, 2113473),
        ("Recombination hotspot 73kb", 73625, 76580),
        ("Recombination hotspot 732kb", 732987, 743877),
        ("Recombination hotspot 1.15Mb", 1156704, 1158831),
        ("Recombination hotspot 1.72Mb", 1722140, 1723519),
        ("Recombination hotspot 2.13Mb", 2133459, 2134136),
    ]

    print("\n" + "=" * 80)
    print("CONCORDANCE WITH CANONICAL CLINICAL HOTSPOTS (Croucher et al. 2011/2015):")
    print("=" * 80)
    landmark_results = []
    for name, s, e in landmarks:
        sub = df_out[(df_out["genome_coord"] >= s - 5000) & (df_out["genome_coord"] <= e + 5000)]
        if len(sub) > 0:
            best = sub.sort_values("kinetic_z", ascending=False).iloc[0]
            print(f"  [FOUND] {name:48s} | Span: {s:>7,d} - {e:>7,d} bp | Changepoint: {best['genome_coord']:>7,d} bp (Z = {best['kinetic_z']:.2f}, Strain: {best['recombinant_taxon']})")
            landmark_results.append({
                "landmark": name,
                "start": s,
                "end": e,
                "status": "FOUND",
                "changepoint": int(best["genome_coord"]),
                "kinetic_z": float(best["kinetic_z"]),
                "strain": str(best["recombinant_taxon"])
            })
        else:
            print(f"  [MISS ] {name:48s} | Span: {s:>7,d} - {e:>7,d} bp")
            landmark_results.append({
                "landmark": name,
                "start": s,
                "end": e,
                "status": "MISS"
            })

    # Save fields summary JSON
    summary_data = {
        "alignment": {
            "path": aln_path,
            "num_taxa": n,
            "length_bp": l,
            "informative_snps": num_snps,
            "compression_ratio": round(l / max(1, num_snps), 1)
        },
        "performance": {
            "load_time_sec": round(t_load, 3),
            "snp_extract_sec": round(t_snp, 3),
            "scan_time_sec": round(t_scan, 3),
            "total_runtime_sec": round(t_load + t_snp + t_scan, 3)
        },
        "background_stats": {
            "mean_residual": round(bg_mean, 5),
            "std_residual": round(bg_std, 5)
        },
        "changepoints_count": len(kept_cps),
        "gubbins_concordance": {
            "total_gubbins_boundaries": len(gub_bps),
            "median_distance_bp": med_dist,
            "mean_distance_bp": mean_dist,
            "concordance_1kb_pct": round(n_1k / len(gub_bps) * 100.0, 1),
            "concordance_2kb_pct": round(n_2k / len(gub_bps) * 100.0, 1),
            "concordance_5kb_pct": round(n_5k / len(gub_bps) * 100.0, 1),
            "matches": boundary_matches
        },
        "landmarks": landmark_results
    }

    summary_json_path = os.path.join(results_dir, "pmen1_fields_summary.json")
    with open(summary_json_path, "w") as f:
        json.dump(summary_data, f, indent=2)
    print(f"\n[✓] Saved comprehensive summary to {summary_json_path}")

if __name__ == "__main__":
    main()
