#!/usr/bin/env python3
"""
analyze_grand_sweep.py
Comprehensive statistical analysis and operational phase diagram generator
for the RhizAeon 50k-100k continuous parameter sweep.
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def analyze_sweep(csv_path: str, out_prefix: str = "grand_sweep"):
    print(f"Loading sweep data from: {csv_path}...")
    df = pd.read_csv(csv_path)
    total_trials = len(df)
    print(f"Loaded {total_trials:,} total simulation trials.")

    # Split into Recombinant and Null
    recs = df[df["is_recombinant_true"] == 1].copy()
    nulls = df[df["is_recombinant_true"] == 0].copy()

    n_recs = len(recs)
    n_nulls = len(nulls)
    print(f"  Recombinant trials: {n_recs:,}")
    print(f"  Null trials:        {n_nulls:,}")

    # =========================================================================
    # 1. Null Partition Statistics (Empirical False Positive Calibration)
    # =========================================================================
    null_fps = (nulls["detected"] > 0).sum()
    null_fpr = null_fps / n_nulls if n_nulls > 0 else 0.0
    bg_taxa_fps = nulls["bg_fp_count"].sum()
    total_null_taxa = nulls["num_taxa"].sum()
    taxon_level_fpr = bg_taxa_fps / total_null_taxa if total_null_taxa > 0 else 0.0

    print("\n--- NULL PARTITION (FALSE POSITIVE RATE CALIBRATION) ---")
    print(f"  Alignment-level FPR: {null_fpr * 100:.4f}% ({null_fps} / {n_nulls:,})")
    print(f"  Taxon-level FPR:     {taxon_level_fpr * 100:.5f}% ({bg_taxa_fps} / {total_null_taxa:,})")

    # Stratified by Heterotachy Alpha
    strict_clock_nulls = nulls[nulls["heterotachy_alpha"] >= 50.0]
    gamma_nulls = nulls[nulls["heterotachy_alpha"] < 50.0]
    sc_fpr = (strict_clock_nulls["detected"] > 0).sum() / len(strict_clock_nulls) if len(strict_clock_nulls) > 0 else 0.0
    gamma_fpr = (gamma_nulls["detected"] > 0).sum() / len(gamma_nulls) if len(gamma_nulls) > 0 else 0.0
    print(f"  Strict clock null FPR: {sc_fpr * 100:.4f}%")
    print(f"  Gamma heterotachy FPR: {gamma_fpr * 100:.4f}%")

    # =========================================================================
    # 2. Recombinant Partition Statistics
    # =========================================================================
    overall_power = (recs["detected"] > 0).sum() / n_recs if n_recs > 0 else 0.0
    mean_recall = recs["recall"].mean()
    mean_precision = recs["precision"].mean()
    
    # Detected subset statistics
    detected_recs = recs[recs["detected"] > 0]
    det_jaccard = detected_recs["jaccard"].mean()
    det_bp1_err = detected_recs["bp1_error"].median()
    det_bp2_err = detected_recs["bp2_error"].median()

    print("\n--- RECOMBINANT PARTITION (OVERALL POWER & ACCURACY) ---")
    print(f"  Overall Empirical Power: {overall_power * 100:.2f}%")
    print(f"  Mean Taxon Recall:       {mean_recall * 100:.2f}%")
    print(f"  Mean Taxon Precision:    {mean_precision * 100:.2f}%")
    print(f"  Mean Tract Jaccard (when detected): {det_jaccard:.4f}")
    print(f"  Median Breakpoint 1 Error: {det_bp1_err:.1f} nt")
    print(f"  Median Breakpoint 2 Error: {det_bp2_err:.1f} nt")

    # =========================================================================
    # 3. Information-Theoretic Resolution Envelope
    # =========================================================================
    print("\n--- INFORMATION-THEORETIC RESOLUTION (MUTATIONS IN TRACT) ---")
    recs["mut_bin"] = pd.cut(
        recs["mut_in_tract"],
        bins=[-1, 2, 5, 10, 15, 25, 50, 100, 10000],
        labels=["0-2", "3-5", "6-10", "11-15", "16-25", "26-50", "51-100", ">100"]
    )
    mut_summary = recs.groupby("mut_bin", observed=True).agg(
        n_trials=("detected", "count"),
        power=("detected", "mean"),
        mean_jaccard=("jaccard", "mean")
    )
    for idx, row in mut_summary.iterrows():
        print(f"  M_tract {idx:8s}: Trials={int(row['n_trials']):5d} | Power={row['power']*100:5.1f}% | Jaccard={row['mean_jaccard']:.3f}")

    # =========================================================================
    # 4. Generate Publication Figures
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
    plt.subplots_adjust(wspace=0.28, hspace=0.32)

    # --- Plot A: 2D Operational Phase Diagram (d_parent vs tract_len / length) ---
    ax = axes[0, 0]
    recs["tract_ratio"] = recs["tract_len"] / recs["length"]
    
    # 2D Grid Binning
    n_bins_x, n_bins_y = 35, 35
    x_edges = np.geomspace(0.005, 0.25, n_bins_x + 1)
    y_edges = np.geomspace(0.008, 0.50, n_bins_y + 1)
    
    grid_power = np.zeros((n_bins_y, n_bins_x))
    grid_counts = np.zeros((n_bins_y, n_bins_x))

    x_idx = np.digitize(recs["d_parent"], x_edges) - 1
    y_idx = np.digitize(recs["tract_ratio"], y_edges) - 1

    for i in range(len(recs)):
        xi = x_idx[i]
        yi = y_idx[i]
        if 0 <= xi < n_bins_x and 0 <= yi < n_bins_y:
            grid_power[yi, xi] += recs["detected"].iloc[i]
            grid_counts[yi, xi] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        grid_rate = np.where(grid_counts >= 3, grid_power / grid_counts, np.nan)

    # Smooth for clean isoclines
    smoothed_rate = gaussian_filter(np.nan_to_num(grid_rate, nan=0.0), sigma=1.0)
    
    x_centers = np.sqrt(x_edges[:-1] * x_edges[1:])
    y_centers = np.sqrt(y_edges[:-1] * y_edges[1:])
    X, Y = np.meshgrid(x_centers, y_centers)

    c = ax.pcolormesh(x_edges, y_edges, grid_rate, cmap="Spectral_r", vmin=0.0, vmax=1.0, shading="flat")
    cb = fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Detection Power $P(\\mathrm{Detect})$", fontsize=10)

    # Contours
    cs = ax.contour(X, Y, smoothed_rate, levels=[0.50, 0.90, 0.99], colors=["white", "yellow", "cyan"], linewidths=1.8)
    ax.clabel(cs, fmt=lambda v: f"{int(v*100)}%", fontsize=9, inline=True)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Parental Divergence $d_{\\mathrm{parent}}$ (subs/site)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Relative Tract Length $L_{\\mathrm{tract}} / L$", fontsize=11, fontweight="bold")
    ax.set_title("A. Continuous Operational Phase Diagram", fontsize=12, fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.4, color="gray")

    # --- Plot B: Information-Theoretic Resolution Curve ---
    ax = axes[0, 1]
    m_bins = np.geomspace(1, 500, 30)
    recs["mut_bin_num"] = pd.cut(recs["mut_in_tract"], bins=m_bins)
    mut_power = recs.groupby("mut_bin_num", observed=True).agg(
        m_mid=("mut_in_tract", "median"),
        power=("detected", "mean"),
        count=("detected", "count")
    ).dropna()

    ax.plot(mut_power["m_mid"], mut_power["power"] * 100, marker="o", color="#1f77b4", lw=2.2, label="Empirical Power")
    ax.axhline(50, color="gray", ls="--", alpha=0.7, label="50% Threshold")
    ax.axhline(90, color="orange", ls="--", alpha=0.7, label="90% Threshold")
    ax.axhline(99, color="green", ls="--", alpha=0.7, label="99% Threshold")

    ax.set_xscale("log")
    ax.set_xlabel("Informative Mutations in Recombinant Tract ($M_{\\mathrm{tract}}$)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Empirical Detection Power (%)", fontsize=11, fontweight="bold")
    ax.set_title("B. Information-Theoretic Resolution Limit", fontsize=12, fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.set_ylim(-2, 104)
    ax.legend(loc="lower right", frameon=True)

    # --- Plot C: Breakpoint Localization Error ---
    ax = axes[1, 0]
    det_sub = recs[(recs["detected"] > 0) & (recs["mut_in_tract"] >= 5)].copy()
    m_err_bins = np.geomspace(5, 500, 20)
    det_sub["mut_bin_err"] = pd.cut(det_sub["mut_in_tract"], bins=m_err_bins)
    bp_err_df = det_sub.groupby("mut_bin_err", observed=True).agg(
        m_mid=("mut_in_tract", "median"),
        bp1_med=("bp1_error", "median"),
        bp1_q25=("bp1_error", lambda s: s.quantile(0.25)),
        bp1_q75=("bp1_error", lambda s: s.quantile(0.75)),
        bp2_med=("bp2_error", "median"),
        count=("bp1_error", "count")
    ).dropna()

    ax.plot(bp_err_df["m_mid"], bp_err_df["bp1_med"], marker="s", color="#d62728", lw=2, label="Breakpoint 1 Error")
    ax.fill_between(bp_err_df["m_mid"], bp_err_df["bp1_q25"], bp_err_df["bp1_q75"], color="#d62728", alpha=0.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Informative Mutations in Recombinant Tract ($M_{\\mathrm{tract}}$)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Breakpoint Spatial Error $|\\hat{u} - u^*|$ (nt)", fontsize=11, fontweight="bold")
    ax.set_title("C. Sub-Gene Breakpoint Spatial Resolution", fontsize=12, fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(loc="upper right", frameon=True)

    # --- Plot D: Computational Throughput & Scaling vs Problem Size ---
    ax = axes[1, 1]
    df["problem_size"] = df["num_taxa"] * df["length"]
    size_bins = np.geomspace(df["problem_size"].min(), df["problem_size"].max(), 25)
    df["size_bin"] = pd.cut(df["problem_size"], bins=size_bins)
    perf_df = df.groupby("size_bin", observed=True).agg(
        size_mid=("problem_size", "median"),
        time_med=("run_time_ms", "median"),
        time_q25=("run_time_ms", lambda s: s.quantile(0.25)),
        time_q75=("run_time_ms", lambda s: s.quantile(0.75)),
        count=("run_time_ms", "count")
    ).dropna()

    ax.plot(perf_df["size_mid"] / 1000.0, perf_df["time_med"], marker="^", color="#2ca02c", lw=2.2, label="RhizAeon Empirical Runtime")
    ax.fill_between(perf_df["size_mid"] / 1000.0, perf_df["time_q25"], perf_df["time_q75"], color="#2ca02c", alpha=0.2)

    # Theoretical O(N^3 L^2) traditional scaling curve for comparison
    ref_x = perf_df["size_mid"].values / 1000.0
    ref_y = 0.05 * (ref_x / ref_x[0]) ** 1.3
    ax.plot(ref_x, ref_y, ls="--", color="black", alpha=0.7, label="Empirical Scaling $\\mathcal{O}((N M U)^{\\approx 1.3})$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Alignment Problem Size ($N \\times L / 10^3$ matrix cells)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Execution Time per Trial (ms)", fontsize=11, fontweight="bold")
    ax.set_title("D. Computational Scaling ($N \\leq 250, L \\leq 30\\,\\mathrm{kb}$)", fontsize=12, fontweight="bold", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(loc="upper left", frameon=True)

    out_fig = f"{out_prefix}_phase_diagram.png"
    plt.tight_layout()
    plt.savefig(out_fig, dpi=300)
    print(f"\nPublication figure saved to: {out_fig}")

    # Output detailed summary metrics to text/LaTeX table
    summary_txt = f"{out_prefix}_summary.txt"
    with open(summary_txt, "w") as f:
        f.write("================================================================================\n")
        f.write("               RHIZAEON GRAND-SCALE CONTINUOUS SWEEP SUMMARY                    \n")
        f.write("================================================================================\n")
        f.write(f"Total Trials:              {total_trials:,}\n")
        f.write(f"Null Trials:              {n_nulls:,}\n")
        f.write(f"Null False Positive Rate: {null_fpr * 100:.4f}% ({null_fps}/{n_nulls})\n")
        f.write(f"Taxon-level False Alarm:  {taxon_level_fpr * 100:.5f}%\n")
        f.write(f"Recombinant Trials:       {n_recs:,}\n")
        f.write(f"Overall Empirical Power:  {overall_power * 100:.2f}%\n")
        f.write(f"Mean Taxon Precision:     {mean_precision * 100:.2f}%\n")
        f.write(f"Mean Taxon Recall:        {mean_recall * 100:.2f}%\n")
        f.write(f"Median BP1 Error:         {det_bp1_err:.1f} nt\n")
        f.write(f"Median BP2 Error:         {det_bp2_err:.1f} nt\n")
        f.write(f"Median Runtime per Trial: {df['run_time_ms'].median():.2f} ms\n")
        f.write("================================================================================\n")
    print(f"Summary report written to: {summary_txt}")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "/Users/sergei/Projects/rhizaeon-code/grand_sweep_50k.csv"
    prefix = sys.argv[2] if len(sys.argv) > 2 else "grand_sweep_50k"
    analyze_sweep(path, prefix)
