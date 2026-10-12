#!/usr/bin/env python3
"""
plot_operational_envelope.py
Production figure generator for the RhizAeon 1,000,000-genome continuous
simulation sweep operational envelope.

Generates:
1. fig_operational_envelope_3panel.png (and .pdf):
   - Panel A: 2D Operational Phase Envelope (Parental Divergence vs Relative Tract Length)
   - Panel B: Information-Theoretic Observability Frontier (Power vs Informative Mutations)
   - Panel C: Sub-Gene Breakpoint Spatial Resolution (Median Error & IQR vs Informative Mutations)
2. fig_operational_envelope_2panel.png (and .pdf):
   - Panel A: 2D Operational Phase Envelope
   - Panel B: Dual-Axis Resolution & Localization Limits (Power % and BP Error nt vs M_tract)
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

# Set publication typography standards (BioVis-Expert)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9.5,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

def generate_envelope_figures(csv_path: str, out_dir: str = "."):
    os.makedirs(out_dir, exist_ok=True)
    print(f"Loading 1M sweep data from: {csv_path}...")
    df = pd.read_csv(csv_path)
    total_trials = len(df)
    print(f"Loaded {total_trials:,} trials.")

    recs = df[df["is_recombinant_true"] == 1].copy()
    nulls = df[df["is_recombinant_true"] == 0].copy()
    n_recs = len(recs)
    n_nulls = len(nulls)
    print(f"  Recombinant trials: {n_recs:,}")
    print(f"  Null trials:        {n_nulls:,}")

    null_fps = (nulls["detected"] > 0).sum()
    null_fpr = null_fps / n_nulls if n_nulls > 0 else 0.0
    print(f"  Null FPR:           {null_fpr * 100:.4f}% ({null_fps} / {n_nulls:,})")

    recs["tract_ratio"] = recs["tract_len"] / recs["length"]

    # =========================================================================
    # PRECOMPUTE DATA FOR PANELS
    # =========================================================================

    # --- 1. 2D Phase Diagram Grid ---
    n_bins_x, n_bins_y = 35, 35
    x_edges = np.geomspace(0.005, 0.35, n_bins_x + 1)
    y_edges = np.geomspace(0.010, 0.45, n_bins_y + 1)

    grid_power = np.zeros((n_bins_y, n_bins_x))
    grid_counts = np.zeros((n_bins_y, n_bins_x))

    x_idx = np.digitize(recs["d_parent"], x_edges) - 1
    y_idx = np.digitize(recs["tract_ratio"], y_edges) - 1

    valid_mask = (x_idx >= 0) & (x_idx < n_bins_x) & (y_idx >= 0) & (y_idx < n_bins_y)
    valid_x = x_idx[valid_mask]
    valid_y = y_idx[valid_mask]
    valid_det = recs["detected"].values[valid_mask]

    np.add.at(grid_power, (valid_y, valid_x), valid_det)
    np.add.at(grid_counts, (valid_y, valid_x), 1)

    with np.errstate(divide="ignore", invalid="ignore"):
        grid_rate = np.where(grid_counts >= 5, grid_power / grid_counts, np.nan)

    # Normalized Gaussian convolution to prevent boundary edge curl
    mask = ~np.isnan(grid_rate)
    grid_filled = np.nan_to_num(grid_rate, nan=0.0)
    smoothed_num = gaussian_filter(grid_filled, sigma=0.85, mode="nearest")
    smoothed_den = gaussian_filter(mask.astype(float), sigma=0.85, mode="nearest")
    with np.errstate(divide="ignore", invalid="ignore"):
        smoothed_rate = np.where(smoothed_den > 0.05, smoothed_num / smoothed_den, np.nan)

    x_centers = np.sqrt(x_edges[:-1] * x_edges[1:])
    y_centers = np.sqrt(y_edges[:-1] * y_edges[1:])
    X, Y = np.meshgrid(x_centers, y_centers)

    # --- 2. 1D Resolution Curve (Power vs M_tract) ---
    m_bins = np.geomspace(1, 500, 32)
    recs["mut_bin_num"] = pd.cut(recs["mut_in_tract"], bins=m_bins)
    mut_power = recs.groupby("mut_bin_num", observed=True).agg(
        m_mid=("mut_in_tract", "median"),
        power=("detected", "mean"),
        count=("detected", "count")
    ).dropna()

    # --- 3. Breakpoint Spatial Error vs M_tract ---
    det_sub = recs[(recs["detected"] > 0) & (recs["mut_in_tract"] >= 5)].copy()
    m_err_bins = np.geomspace(5, 500, 22)
    det_sub["mut_bin_err"] = pd.cut(det_sub["mut_in_tract"], bins=m_err_bins)
    bp_err_df = det_sub.groupby("mut_bin_err", observed=True).agg(
        m_mid=("mut_in_tract", "median"),
        bp1_med=("bp1_error", "median"),
        bp1_q25=("bp1_error", lambda s: s.quantile(0.25)),
        bp1_q75=("bp1_error", lambda s: s.quantile(0.75)),
        bp2_med=("bp2_error", "median"),
        bp2_q25=("bp2_error", lambda s: s.quantile(0.25)),
        bp2_q75=("bp2_error", lambda s: s.quantile(0.75)),
        count=("bp1_error", "count")
    ).dropna()

    # =========================================================================
    # FIGURE 1: 3-PANEL OPERATIONAL FLIGHT ENVELOPE (1x3 HORIZONTAL BANNER)
    # =========================================================================
    print("Rendering 3-Panel Operational Flight Envelope...")
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), dpi=300)
    plt.subplots_adjust(left=0.06, right=0.96, bottom=0.15, top=0.88, wspace=0.32)

    # --- Panel A: 2D Continuous Operational Phase Envelope ---
    ax = axes[0]
    c = ax.pcolormesh(x_edges, y_edges, grid_rate, cmap="viridis", vmin=0.0, vmax=1.0, shading="flat")
    cb = fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Detection Probability $P(\\mathrm{Detect})$", fontsize=9, fontweight="bold")
    cb.ax.tick_params(labelsize=8)

    cs = ax.contour(X, Y, smoothed_rate, levels=[0.50, 0.90], colors=["#ffffff", "#ffd700"], linewidths=[2.0, 2.0])
    ax.clabel(cs, fmt={0.50: "50% Power", 0.90: "90% Power"}, fontsize=8.5, inline=True)

    # Stall limit hyperbola: M = d * L_tract = 3 SNPs at median L = 5000 nt -> L_tract/L = 3 / (5000 * d)
    d_vals = np.geomspace(0.005, 0.35, 100)
    median_L = 5000.0
    stall_ratio = 3.0 / (median_L * d_vals)
    valid_stall = (stall_ratio >= 0.010) & (stall_ratio <= 0.45)
    ax.plot(d_vals[valid_stall], stall_ratio[valid_stall], ls=":", color="#e74c3c", lw=1.8, label="Stall Horizon ($M < 3$ SNPs)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, 0.35)
    ax.set_ylim(0.010, 0.45)
    ax.set_xlabel("Parental Divergence $d_{\\mathrm{parent}}$ (subs/site)")
    ax.set_ylabel("Relative Tract Length $L_{\\mathrm{tract}} / L$")
    ax.set_title("A. Continuous Operational Flight Envelope", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.3, color="gray")

    # --- Panel B: Information-Theoretic Resolution Limit ---
    ax = axes[1]
    ax.plot(mut_power["m_mid"], mut_power["power"] * 100, marker="o", markersize=4.5,
            color="#1a5276", lw=2.2, label="Empirical Detection Power", zorder=4)
    ax.axhline(50, color="#7f8c8d", ls="--", lw=1.2, alpha=0.85, label="50% Operational Ceiling")
    ax.axhline(90, color="#d35400", ls="--", lw=1.2, alpha=0.85, label="90% Saturation Horizon")

    # Null FPR floor
    ax.axhspan(0, null_fpr * 100, color="#27ae60", alpha=0.25, zorder=2)
    ax.axhline(null_fpr * 100, color="#27ae60", ls="-", lw=1.0, alpha=0.9,
               label=f"Non-Rec Null FPR: {null_fpr*100:.2f}%")

    # Highlight stall vs saturation
    ax.axvspan(1, 2.5, color="#e74c3c", alpha=0.08, zorder=1)
    ax.text(1.2, 8, "Stall Horizon\n(Power < 1%)", fontsize=7.5, color="#c0392b", fontweight="bold")

    ax.set_xscale("log")
    ax.set_xlim(1.0, 500)
    ax.set_ylim(-2, 104)
    ax.set_xlabel("Informative Mutations in Tract ($M_{\\mathrm{tract}}$)")
    ax.set_ylabel("Empirical Detection Power (%)")
    ax.set_title("B. Informational Observability Horizon", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.35, color="gray")
    ax.legend(loc="lower right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#dcdde1")

    # --- Panel C: Sub-Gene Breakpoint Spatial Resolution ---
    ax = axes[2]
    ax.plot(bp_err_df["m_mid"], bp_err_df["bp1_med"], marker="s", markersize=4.5,
            color="#c0392b", lw=2.0, label="Median Breakpoint Error", zorder=4)
    ax.fill_between(bp_err_df["m_mid"], bp_err_df["bp1_q25"], bp_err_df["bp1_q75"],
                    color="#c0392b", alpha=0.18, label="Interquartile Range (25th–75th %ile)", zorder=3)

    # Reference line for sub-gene margin (15 nt)
    ax.axhline(15, color="#2c3e50", ls=":", lw=1.3, alpha=0.7, label="Sub-15 nt Gene Margin")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(4.5, 500)
    ax.set_ylim(2, 130)
    ax.set_xlabel("Informative Mutations in Tract ($M_{\\mathrm{tract}}$)")
    ax.set_ylabel("Spatial Breakpoint Error $|\\hat{u} - u^*|$ (nt)")
    ax.set_title("C. Sub-Gene Breakpoint Precision", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.35, color="gray")
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#dcdde1")

    # Save 3-panel figures
    fig3_png = os.path.join(out_dir, "fig_operational_envelope_3panel.png")
    fig3_pdf = os.path.join(out_dir, "fig_operational_envelope_3panel.pdf")
    plt.savefig(fig3_png, dpi=300)
    plt.savefig(fig3_pdf)
    plt.close()
    print(f"✓ Saved: {fig3_png}")
    print(f"✓ Saved: {fig3_pdf}")

    # =========================================================================
    # FIGURE 2: 2-PANEL COMPACT ENVELOPE (1x2 SIDE-BY-SIDE)
    # =========================================================================
    print("Rendering 2-Panel Compact Operational Envelope...")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.92, bottom=0.15, top=0.88, wspace=0.34)

    # --- Panel A: 2D Operational Phase Diagram ---
    ax = axes[0]
    c = ax.pcolormesh(x_edges, y_edges, grid_rate, cmap="viridis", vmin=0.0, vmax=1.0, shading="flat")
    cb = fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Detection Probability $P(\\mathrm{Detect})$", fontsize=9, fontweight="bold")
    cb.ax.tick_params(labelsize=8)

    cs = ax.contour(X, Y, smoothed_rate, levels=[0.50, 0.90], colors=["#ffffff", "#ffd700"], linewidths=[2.0, 2.0])
    ax.clabel(cs, fmt={0.50: "50% Power", 0.90: "90% Power"}, fontsize=8.5, inline=True)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, 0.35)
    ax.set_ylim(0.010, 0.45)
    ax.set_xlabel("Parental Divergence $d_{\\mathrm{parent}}$ (subs/site)")
    ax.set_ylabel("Relative Tract Length $L_{\\mathrm{tract}} / L$")
    ax.set_title("A. Continuous Operational Phase Envelope", loc="left")
    ax.grid(True, which="both", ls=":", alpha=0.3, color="gray")

    # --- Panel B: Dual-Axis Power & Spatial Precision ---
    ax1 = axes[1]
    ax2 = ax1.twinx()

    # Left axis: Power
    l1 = ax1.plot(mut_power["m_mid"], mut_power["power"] * 100, marker="o", markersize=4.5,
                  color="#1a5276", lw=2.2, label="Detection Power (%)", zorder=4)
    ax1.axhline(50, color="#7f8c8d", ls="--", lw=1.0, alpha=0.7)
    ax1.axhline(90, color="#d35400", ls="--", lw=1.0, alpha=0.7)
    ax1.set_xscale("log")
    ax1.set_xlim(1.0, 500)
    ax1.set_ylim(-2, 104)
    ax1.set_xlabel("Informative Mutations in Tract ($M_{\\mathrm{tract}}$)")
    ax1.set_ylabel("Empirical Detection Power (%)", color="#1a5276")
    ax1.tick_params(axis="y", labelcolor="#1a5276")

    # Right axis: Breakpoint Error
    l2 = ax2.plot(bp_err_df["m_mid"], bp_err_df["bp1_med"], marker="s", markersize=4.0,
                  color="#c0392b", lw=1.8, label="Median BP Error (nt)", zorder=3)
    ax2.fill_between(bp_err_df["m_mid"], bp_err_df["bp1_q25"], bp_err_df["bp1_q75"],
                     color="#c0392b", alpha=0.15, zorder=2)
    ax2.set_yscale("log")
    ax2.set_ylim(2, 100)
    ax2.set_ylabel("Breakpoint Spatial Error $|\\hat{u} - u^*|$ (nt)", color="#c0392b")
    ax2.tick_params(axis="y", labelcolor="#c0392b")
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color("#c0392b")
    ax2.spines["left"].set_color("#1a5276")

    ax1.set_title("B. Observability & Sub-Gene Precision Frontier", loc="left")
    ax1.grid(True, which="both", ls=":", alpha=0.3, color="gray")

    # Unified Legend
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="center right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#dcdde1")

    fig2_png = os.path.join(out_dir, "fig_operational_envelope_2panel.png")
    fig2_pdf = os.path.join(out_dir, "fig_operational_envelope_2panel.pdf")
    plt.savefig(fig2_png, dpi=300)
    plt.savefig(fig2_pdf)
    plt.close()
    print(f"✓ Saved: {fig2_png}")
    print(f"✓ Saved: {fig2_pdf}")

    print("\nAll figures generated successfully.")

if __name__ == "__main__":
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "/home/sergei/grand_sweep_1m/grand_sweep_1m.csv"
    output_directory = sys.argv[2] if len(sys.argv) > 2 else "."
    generate_envelope_figures(csv_file, output_directory)
