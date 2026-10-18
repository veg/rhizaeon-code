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
    m_bins = np.geomspace(1, 1000, 36)
    recs["mut_bin_num"] = pd.cut(recs["mut_in_tract"], bins=m_bins)
    mut_power = recs.groupby("mut_bin_num", observed=True).agg(
        m_mid=("mut_in_tract", "median"),
        power=("detected", "mean"),
        count=("detected", "count")
    ).dropna()

    # --- 3. Breakpoint Spatial Error vs M_tract (Matched Events: Jaccard > 0) ---
    det_sub = recs[(recs["detected"] > 0) & (recs["jaccard"] > 0.0) & (recs["mut_in_tract"] >= 5)].copy()
    m_err_bins = np.geomspace(5, 1000, 26)
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
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8), dpi=300)
    plt.subplots_adjust(left=0.06, right=0.96, bottom=0.15, top=0.88, wspace=0.32)

    # --- Panel A: 2D Continuous Operational Flight Envelope ---
    ax = axes[0]
    # Smooth continuous filled contours
    cf_levels = np.linspace(0.0, 1.0, 21)
    cf = ax.contourf(X, Y, smoothed_rate, levels=cf_levels, cmap="viridis", alpha=0.92)
    cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Detection Probability $P(\\mathrm{Detect})$", fontsize=9.5, fontweight="bold")
    cb.ax.tick_params(labelsize=8)

    # Iso-power contour curves (50% Operational Boundary & 90% Certified Threshold)
    cs_50 = ax.contour(X, Y, smoothed_rate, levels=[0.50], colors=["#ffffff"], linewidths=[2.0], linestyles=["-"])
    cs_90 = ax.contour(X, Y, smoothed_rate, levels=[0.90], colors=["#ffd700"], linewidths=[2.2], linestyles=["-"])

    cl_50 = ax.clabel(cs_50, fmt={0.50: "50% Boundary"}, fontsize=8.0, inline=True, manual=[(0.022, 0.14)], colors="#ffffff")
    cl_90 = ax.clabel(cs_90, fmt={0.90: "90% Certified"}, fontsize=8.0, inline=True, manual=[(0.22, 0.065)], colors="#7d6608")
    for txt in cl_50 + cl_90:
        txt.set_fontweight("bold")

    # Stall limit hyperbola: M = d * L_tract = 3 SNPs at median L = 5000 nt -> L_tract/L = 3 / (5000 * d)
    d_vals = np.geomspace(0.005, 0.35, 150)
    median_L = 5000.0
    stall_ratio = 3.0 / (median_L * d_vals)
    valid_stall = (stall_ratio >= 0.010) & (stall_ratio <= 0.45)
    
    # Shade Stall Zone with diagonal hatching
    ax.fill_between(d_vals[valid_stall], 0.010, stall_ratio[valid_stall],
                    facecolor="#e74c3c", alpha=0.22, hatch="//", edgecolor="#c0392b", zorder=3)
    ax.plot(d_vals[valid_stall], stall_ratio[valid_stall], ls="--", color="#e74c3c", lw=1.8,
            label="Stall Horizon ($M < 3$ SNPs)", zorder=4)
    ax.text(0.0065, 0.013, "Informational Stall\n($M < 3$ SNPs)",
            fontsize=7.5, color="#ffffff", fontweight="bold", zorder=5)

    # Biological Landmarks
    landmarks = [
        {"name": "S. pneumoniae PMEN1", "d": 0.08, "r": 0.20, "marker": "s", "color": "#1b9e77", "dx": -6, "dy": 8, "ha": "center"},
        {"name": "HIV-1 CRF01_AE", "d": 0.12, "r": 0.16, "marker": "D", "color": "#e7298a", "dx": -4, "dy": -16, "ha": "center"},
        {"name": "Sarbecovirus Spike RBD", "d": 0.04, "r": 0.035, "marker": "^", "color": "#d95f02", "dx": 8, "dy": -4, "ha": "left"},
        {"name": "Primate MHC (Micro-conv.)", "d": 0.05, "r": 0.07, "marker": "o", "color": "#7570b3", "dx": 4, "dy": 8, "ha": "center"},
    ]
    for lm in landmarks:
        ax.scatter(lm["d"], lm["r"], marker=lm["marker"], s=60, color=lm["color"],
                   edgecolors="white", linewidths=1.2, zorder=6)
        ax.annotate(lm["name"], xy=(lm["d"], lm["r"]), xytext=(lm["dx"], lm["dy"]), textcoords="offset points",
                    fontsize=7.3, fontweight="bold", color="#111111", ha=lm["ha"],
                    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#cccccc", lw=0.6, alpha=0.92), zorder=7)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, 0.35)
    ax.set_ylim(0.010, 0.45)
    ax.set_xlabel("Parental Divergence $d_{\\mathrm{parent}}$ (subs/site)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Relative Tract Length $L_{\\mathrm{tract}} / L$", fontsize=10, fontweight="bold")
    ax.set_title("A. Aeronautical Operating Flight Envelope", loc="left", fontsize=11, fontweight="bold")
    ax.grid(True, which="both", ls=":", alpha=0.25, color="gray")

    # --- Panel B: Informational Observability Frontier ---
    ax = axes[1]
    sat_power = recs[recs["mut_in_tract"] > 100]["detected"].mean() * 100.0
    stall_pwr = recs[recs["mut_in_tract"] <= 2]["detected"].mean() * 100.0
    ax.plot(mut_power["m_mid"], mut_power["power"] * 100, marker="o", markersize=4.0,
            color="#1a5276", lw=2.2, label="Empirical Power", zorder=4)
    ax.axhline(50, color="#7f8c8d", ls="--", lw=1.2, alpha=0.85, label="50% Operational Boundary")
    ax.axhline(90, color="#d35400", ls="--", lw=1.2, alpha=0.85, label="90% Certified Threshold")
    ax.axhline(sat_power, color="#27ae60", ls="--", lw=1.4, alpha=0.90, label=f"{sat_power:.1f}% Saturation Plateau ($M > 100$)")

    # Null FPR floor
    ax.axhspan(0, null_fpr * 100, color="#16a085", alpha=0.20, zorder=2)
    ax.axhline(null_fpr * 100, color="#16a085", ls="-", lw=1.1, alpha=0.95,
               label=f"Non-Rec Null FPR: {null_fpr*100:.2f}%")

    # Highlight stall horizon
    ax.axvspan(1, 2.5, color="#e74c3c", alpha=0.08, zorder=1)
    ax.text(1.15, 25, f"Stall Horizon\n(Power $\\approx$ {stall_pwr:.1f}%)", fontsize=7.5, color="#c0392b", fontweight="bold")

    ax.set_xscale("log")
    ax.set_xlim(1.0, 1000)
    ax.set_ylim(-2, 106)
    ax.set_xlabel("Informative Mutations in Tract ($M_{\\mathrm{tract}}$)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Empirical Detection Power (%)", fontsize=10, fontweight="bold")
    ax.set_title("B. Informational Observability Frontier", loc="left", fontsize=11, fontweight="bold")
    ax.grid(True, which="both", ls=":", alpha=0.30, color="gray")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, facecolor="white", edgecolor="#cccccc", fontsize=7.8)

    # --- Panel C: Sub-Gene Breakpoint Spatial Resolution ---
    ax = axes[2]
    ax.plot(bp_err_df["m_mid"], bp_err_df["bp1_med"], marker="s", markersize=4.0,
            color="#c0392b", lw=2.0, label="Median Breakpoint Error", zorder=4)
    ax.fill_between(bp_err_df["m_mid"], bp_err_df["bp1_q25"], bp_err_df["bp1_q75"],
                    color="#c0392b", alpha=0.18, label="Interquartile Range (25th–75th %ile)", zorder=3)

    # Reference line for sub-gene margin (15 nt)
    ax.axhline(15, color="#2c3e50", ls=":", lw=1.3, alpha=0.75, label="Sub-15 nt Gene Margin")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(4.5, 1000)
    ax.set_ylim(1.5, 130)
    ax.set_xlabel("Informative Mutations in Tract ($M_{\\mathrm{tract}}$)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Spatial Breakpoint Error $|\\hat{u} - u^*|$ (nt)", fontsize=10, fontweight="bold")
    ax.set_title("C. Sub-Gene Breakpoint Precision", loc="left", fontsize=11, fontweight="bold")
    ax.grid(True, which="both", ls=":", alpha=0.30, color="gray")
    ax.legend(loc="upper right", frameon=True, framealpha=0.92, facecolor="white", edgecolor="#cccccc", fontsize=8.0)

    # Save 3-panel figures (both as operational envelope and as phase diagram replacement)
    fig3_png = os.path.join(out_dir, "fig_operational_envelope_3panel.png")
    fig3_pdf = os.path.join(out_dir, "fig_operational_envelope_3panel.pdf")
    plt.savefig(fig3_png, dpi=300)
    plt.savefig(fig3_pdf)
    print(f"✓ Saved: {fig3_png}")
    print(f"✓ Saved: {fig3_pdf}")

    # Mirror as fig_grand_sweep_1m_phase_diagram for drop-in manuscript integration
    phase_png = os.path.join(out_dir, "fig_grand_sweep_1m_phase_diagram.png")
    phase_pdf = os.path.join(out_dir, "fig_grand_sweep_1m_phase_diagram.pdf")
    plt.savefig(phase_png, dpi=300)
    plt.savefig(phase_pdf)
    plt.close()
    print(f"✓ Mirrored: {phase_png}")
    print(f"✓ Mirrored: {phase_pdf}")

    # =========================================================================
    # FIGURE 2: 2-PANEL COMPACT ENVELOPE (1x2 SIDE-BY-SIDE)
    # =========================================================================
    print("Rendering 2-Panel Compact Operational Envelope...")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.92, bottom=0.15, top=0.88, wspace=0.34)

    # --- Panel A: 2D Operational Phase Diagram ---
    ax = axes[0]
    cf = ax.contourf(X, Y, smoothed_rate, levels=cf_levels, cmap="viridis", alpha=0.92)
    cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Detection Probability $P(\\mathrm{Detect})$", fontsize=9, fontweight="bold")
    cb.ax.tick_params(labelsize=8)

    cs_50 = ax.contour(X, Y, smoothed_rate, levels=[0.50], colors=["#ffffff"], linewidths=[1.8], linestyles=["-"])
    cs_90 = ax.contour(X, Y, smoothed_rate, levels=[0.90], colors=["#ffd700"], linewidths=[2.0], linestyles=["-"])
    cl_50 = ax.clabel(cs_50, fmt={0.50: "50% Boundary"}, fontsize=7.5, inline=True, manual=[(0.022, 0.14)], colors="#ffffff")
    cl_90 = ax.clabel(cs_90, fmt={0.90: "90% Certified"}, fontsize=7.5, inline=True, manual=[(0.22, 0.065)], colors="#7d6608")
    for txt in cl_50 + cl_90:
        txt.set_fontweight("bold")

    # Stall limit hyperbola
    ax.fill_between(d_vals[valid_stall], 0.010, stall_ratio[valid_stall],
                    facecolor="#e74c3c", alpha=0.22, hatch="//", edgecolor="#c0392b", zorder=3)
    ax.plot(d_vals[valid_stall], stall_ratio[valid_stall], ls="--", color="#e74c3c", lw=1.5, zorder=4)
    ax.text(0.0065, 0.013, "Stall ($M < 3$)", fontsize=7.0, color="#ffffff", fontweight="bold", zorder=5)

    # Biological Landmarks
    for lm in landmarks:
        ax.scatter(lm["d"], lm["r"], marker=lm["marker"], s=45, color=lm["color"],
                   edgecolors="white", linewidths=1.0, zorder=6)
        ax.annotate(lm["name"], xy=(lm["d"], lm["r"]), xytext=(lm["dx"], lm["dy"]), textcoords="offset points",
                    fontsize=6.8, fontweight="bold", color="#111111", ha=lm["ha"],
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc", lw=0.5, alpha=0.90), zorder=7)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, 0.35)
    ax.set_ylim(0.010, 0.45)
    ax.set_xlabel("Parental Divergence $d_{\\mathrm{parent}}$ (subs/site)", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("Relative Tract Length $L_{\\mathrm{tract}} / L$", fontsize=9.5, fontweight="bold")
    ax.set_title("A. Continuous Operational Phase Envelope", loc="left", fontsize=10.5, fontweight="bold")
    ax.grid(True, which="both", ls=":", alpha=0.25, color="gray")

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
