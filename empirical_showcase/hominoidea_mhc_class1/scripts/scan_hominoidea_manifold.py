#!/usr/bin/env python3
"""
scan_hominoidea_manifold.py
Whole-gene continuous physical field deconvolution across the
Exhaustive Superfamily Hominoidea MHC Class I dataset (N = 379 taxa, L = 1,128 bp).
"""

import os
import json
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.ndimage import gaussian_filter1d

BASE_DIR = "/Users/sergei/Projects/rhizaeon-code/empirical_showcase/hominoidea_mhc_class1"
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

ALN_PATH = os.path.join(DATA_DIR, "hominoidea_class1_aligned.fasta")

def cmdscale(D, k=4):
    n = D.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    pos = evals > 1e-9
    k_eff = min(k, np.sum(pos))
    if k_eff == 0:
        return np.zeros((n, k))
    coords = evecs[:, :k_eff] * np.sqrt(evals[:k_eff])
    if k_eff < k:
        coords = np.pad(coords, ((0, 0), (0, k - k_eff)))
    return coords

def align_procrustes(z1, z2):
    mu1 = np.mean(z1, axis=0)
    mu2 = np.mean(z2, axis=0)
    z1_c = z1 - mu1
    z2_c = z2 - mu2
    u, _, vt = np.linalg.svd(z2_c.T @ z1_c)
    r = u @ vt
    z2_rot = z2_c @ r + mu1
    res = np.linalg.norm(z1 - z2_rot, axis=1)
    return z2_rot, res

def main():
    t0 = time.time()
    records = list(SeqIO.parse(ALN_PATH, "fasta"))
    taxa = [r.id for r in records]
    N = len(taxa)
    L = len(records[0].seq)
    print(f"[1] Loaded alignment: N = {N} taxa, L = {L} bp")

    # Extract Segregating Informative Sites
    seq_matrix = np.array([list(str(r.seq).upper()) for r in records])
    valid_chars = set("ACGT")
    is_seg = np.zeros(L, dtype=bool)
    for col in range(L):
        chars = set(seq_matrix[:, col]) & valid_chars
        if len(chars) > 1:
            is_seg[col] = True

    seg_indices = np.where(is_seg)[0]
    M = len(seg_indices)
    print(f"[2] Extracted M = {M} segregating informative sites ({100.0 * M / L:.1f}% of alignment)")

    # Build Prefix Hamming Difference Tensor
    compressed_seqs = seq_matrix[:, seg_indices]
    prefix_diff = np.zeros((N, N, M + 1), dtype=np.int32)
    for i in range(N):
        for j in range(i + 1, N):
            mismatches = (compressed_seqs[i] != compressed_seqs[j]) & (compressed_seqs[i] != '-') & (compressed_seqs[j] != '-')
            cum = np.zeros(M + 1, dtype=np.int32)
            cum[1:] = np.cumsum(mismatches)
            prefix_diff[i, j] = cum
            prefix_diff[j, i] = cum

    def query_sub_distance(u1, u2):
        diff = prefix_diff[:, :, u2] - prefix_diff[:, :, u1]
        span = max(1, u2 - u1)
        return diff.astype(np.float64) / span

    # Multi-scale bilateral Procrustes kinetic scan across dyadic scales
    flanks = [15, 30]
    all_evals = []
    print("[3] Performing multi-scale bilateral Procrustes kinetic scan...")
    for flk in flanks:
        for cp in range(flk, M - flk, 2):
            d1 = query_sub_distance(cp - flk, cp)
            d2 = query_sub_distance(cp, cp + flk)
            if np.max(d1) < 1e-6 or np.max(d2) < 1e-6:
                continue
            z1 = cmdscale(d1, k=4)
            z2 = cmdscale(d2, k=4)
            _, res = align_procrustes(z1, z2)
            all_evals.append((flk, cp, res))

    res_matrix = np.array([e[2] for e in all_evals])  # Shape: (num_evals, N)
    bg_mean = float(np.mean(res_matrix))
    bg_std = float(np.std(res_matrix))
    print(f"[4] Background calibration: Mean residual = {bg_mean:.4f}, Std = {bg_std:.4f}")

    # Maximum kinetic disparity across taxa at each SNP position
    snp_max_res = np.zeros(M)
    snp_eval_counts = np.zeros(M)
    for (flk, cp, res) in all_evals:
        max_r = np.max(res)
        if max_r > snp_max_res[cp]:
            snp_max_res[cp] = max_r
        snp_eval_counts[cp] += 1

    snp_z_scores = np.maximum(0.0, (snp_max_res - bg_mean) / bg_std)

    # Project onto physical continuous genomic coordinates u in [0, L-1]
    phi_field = np.zeros(L, dtype=np.float64)
    phi_field[seg_indices] = snp_z_scores
    # Gentle smoothing for visualization
    phi_smooth = gaussian_filter1d(phi_field, sigma=3.0)
    strain_field = np.cumsum(phi_smooth)

    # Taxon-resolved reticulations
    taxon_peaks = []
    for i, t in enumerate(taxa):
        r_i = res_matrix[:, i]
        max_idx = np.argmax(r_i)
        max_r = r_i[max_idx]
        max_z = (max_r - bg_mean) / bg_std
        flk, cp, _ = all_evals[max_idx]
        bp = seg_indices[cp]
        
        # Determine structural domain
        if bp < 72:
            dom = "Exon 1 (Signal)"
        elif bp < 342:
            dom = "Exon 2 (alpha-1 ARS)"
        elif bp < 618:
            dom = "Exon 3 (alpha-2 ARS)"
        elif bp < 894:
            dom = "Exon 4 (alpha-3 Stalk)"
        elif bp < 1011:
            dom = "Exon 5 (Transmembrane)"
        else:
            dom = "Exons 6-7 (Cytoplasmic)"

        taxon_peaks.append({
            "taxon": t,
            "max_z": float(max_z),
            "snp_idx": int(cp),
            "bp": int(bp),
            "flank": int(flk),
            "domain": dom
        })

    df_taxa = pd.DataFrame(taxon_peaks)
    df_sig = df_taxa[df_taxa["max_z"] >= 3.0].sort_values("max_z", ascending=False)
    print(f"[5] Identified {len(df_sig)} / {N} taxa with significant reticulation (Z >= 3.0)")

    # Save data
    np.savez_compressed(
        os.path.join(RESULTS_DIR, "hominoidea_continuous_fields.npz"),
        u_grid=np.arange(L),
        phi_field=phi_smooth,
        strain_field=strain_field,
        seg_indices=seg_indices,
        snp_z_scores=snp_z_scores
    )
    df_taxa.to_csv(os.path.join(RESULTS_DIR, "hominoidea_taxon_peaks.csv"), index=False)
    df_sig.to_csv(os.path.join(RESULTS_DIR, "hominoidea_significant_reticulations.csv"), index=False)

    # Domain summary
    dom_counts = df_sig["domain"].value_counts().to_dict()

    dur = time.time() - t0
    summary = {
        "num_taxa": N,
        "length": L,
        "informative_snps": M,
        "bg_mean": bg_mean,
        "bg_std": bg_std,
        "num_significant_taxa_z3": len(df_sig),
        "runtime_seconds": dur,
        "domain_distribution": dom_counts,
        "b46_changepoint": df_taxa[df_taxa["taxon"].str.contains("HLA_B_46_01")].to_dict(orient="records")[0],
        "pal_changepoint": df_taxa[df_taxa["taxon"].str.contains("Patr_AL_01_01_01_01")].to_dict(orient="records")[0]
    }

    with open(os.path.join(RESULTS_DIR, "hominoidea_fields_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[✓] Analysis successfully completed in {dur:.2f} s. Summary written to JSON.")

if __name__ == "__main__":
    main()
