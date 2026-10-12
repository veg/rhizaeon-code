#!/usr/bin/env python3
"""
generate_fig_posada_crandall.py

Generates a publication-grade, 3-panel figure visualizing the full 3,200-trial
Posada & Crandall (2001) benchmark comparing RhizAeon vs. 3SEQ.

Adheres strictly to MBE/Nature BioVis guidelines:
- Zero narrative text annotations inside plot areas.
- Legends placed in dedicated non-overlapping zones.
- Expected event cardinality E[R] under Hudson's coalescent plotted in Panel B.

Panels:
  (A) Recombination Power vs. Recombination Intensity (rho) across Mutation Rates (theta)
  (B) Event Cardinality: Expected E[R] vs. Inferred Events (K_hat)
  (C) False Positive Rate Control under Gamma Site-Rate Heterogeneity (alpha)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Nature / MBE styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 8.5,
    'axes.labelsize': 9.0,
    'axes.titlesize': 9.5,
    'xtick.labelsize': 8.0,
    'ytick.labelsize': 8.0,
    'legend.fontsize': 7.5,
    'figure.titlesize': 10.5,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.35,
})

CSV_PATH = "/Users/sergei/Projects/rhizaeon-code/benchmarks/posada_crandall_2001/posada_comparison_rhizaeon_vs_3seq.csv"
OUT_PDF = "/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_posada_crandall_benchmark.pdf"
OUT_PNG = "/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_posada_crandall_benchmark.png"

df = pd.read_csv(CSV_PATH)

power_df = df[df['category_rhiz'] == 'power'].copy()
fp_df = df[df['category_rhiz'] == 'false_positive'].copy()

# Double column figure dimensions: 7.2 in wide, 2.7 in tall
fig, axes = plt.subplots(1, 3, figsize=(10.5, 2.95), dpi=300)
fig.subplots_adjust(wspace=0.32, bottom=0.18, top=0.76, left=0.065, right=0.985)

# Color palette for theta
colors = {
    10.0: '#64748b',   # Slate
    50.0: '#0ea5e9',   # Sky blue
    100.0: '#2563eb',  # Royal blue
    200.0: '#1e3a8a',  # Navy
}

rhos = [0.0, 1.0, 4.0, 16.0, 64.0]
rho_labels = ['0', '1', '4', '16', '64']
x_indices = np.arange(len(rhos))

# Hudson's theoretical expectation for n=10 lineages: E[R] = (rho / 2) * sum_{i=1}^9 (1/i)
H_9 = sum(1.0 / i for i in range(1, 10))
expected_R = [0.5 * r * H_9 for r in rhos]

# ----------------------------------------------------
# Panel A: Alignment-Level Detection Power
# ----------------------------------------------------
ax0 = axes[0]
ax0.set_title(r"A. Alignment Power vs. $\rho$", fontweight='bold', loc='left', pad=6)
ax0.set_xlabel(r"Recombination Parameter ($\rho = 4N_e r$)")
ax0.set_ylabel("Detection Power (%)")
ax0.set_xticks(x_indices)
ax0.set_xticklabels(rho_labels)
ax0.set_ylim(-2, 104)
ax0.grid(True, linestyle='--', alpha=0.35)

for theta in [10.0, 50.0, 100.0, 200.0]:
    sub = power_df[power_df['theta_rhiz'] == theta].sort_values('rho_rhiz')
    ax0.plot(x_indices, sub['detection_pct_rhiz'], marker='o', markersize=4.2,
             linewidth=1.6, color=colors[theta], zorder=4)
    ax0.plot(x_indices, sub['detection_pct_3seq'], marker='s', markersize=3.8,
             linewidth=1.2, linestyle=':', color=colors[theta], alpha=0.85, zorder=3)

# ----------------------------------------------------
# Panel B: Expected vs. Inferred Event Cardinality
# ----------------------------------------------------
ax1 = axes[1]
ax1.set_title(r"B. Event Cardinality ($\hat{K}$ vs. $\rho$)", fontweight='bold', loc='left', pad=6)
ax1.set_xlabel(r"Recombination Parameter ($\rho = 4N_e r$)")
ax1.set_ylabel(r"Mean Inferred Events ($\hat{K}$)")
ax1.set_xticks(x_indices)
ax1.set_xticklabels(rho_labels)
ax1.grid(True, linestyle='--', alpha=0.35)

# Expected E[R] under Hudson's ARG
ax1.plot(x_indices, expected_R, color='#0f172a', linestyle='--', marker='D', markersize=4.2,
         linewidth=1.5, zorder=5, label=r"Expected $E[R]$")

for theta in [50.0, 100.0, 200.0]:
    sub = power_df[power_df['theta_rhiz'] == theta].sort_values('rho_rhiz')
    ax1.plot(x_indices, sub['mean_K_hat_3seq'], marker='s', markersize=3.8,
             linewidth=1.2, linestyle=':', color=colors[theta], alpha=0.85, zorder=3)
    ax1.plot(x_indices, sub['mean_K_hat_rhiz'], marker='o', markersize=4.2,
             linewidth=1.6, color=colors[theta], zorder=4)

ax1.set_yscale('symlog', linthresh=1.0)
ax1.set_yticks([0, 1, 5, 20, 100])
ax1.set_yticklabels(['0', '1', '5', '20', '100'])
ax1.set_ylim(-0.15, 125)

# ----------------------------------------------------
# Panel C: False Positive Rate under Rate Variation
# ----------------------------------------------------
ax2 = axes[2]
ax2.set_title("C. Null FPR under Heterotachy", fontweight='bold', loc='left', pad=6)
ax2.set_xlabel(r"Gamma Rate Heterogeneity ($\alpha$)")
ax2.set_ylabel("False Positive Rate (%)")

alphas = ['2.0\n(Mild)', '0.5\n(Moderate)', '0.05\n(Severe)']
alpha_vals = [2.0, 0.5, 0.05]
bar_width = 0.32
x_pos = np.arange(len(alphas))

rhiz_fp_means = [fp_df[fp_df['alpha_rhiz'] == a]['detection_pct_rhiz'].mean() for a in alpha_vals]
s3_fp_means = [fp_df[fp_df['alpha_rhiz'] == a]['detection_pct_3seq'].mean() for a in alpha_vals]

ax2.bar(x_pos - bar_width/2, rhiz_fp_means, bar_width, label='RhizAeon (Mean 0.58%)',
        color='#2563eb', edgecolor='#1d4ed8', alpha=0.9, zorder=3)
ax2.bar(x_pos + bar_width/2, s3_fp_means, bar_width, label='3SEQ (Mean 0.17%)',
        color='#cbd5e1', edgecolor='#94a3b8', alpha=0.9, zorder=3)

ax2.axhline(5.0, color='#dc2626', linestyle='--', linewidth=1.0, alpha=0.85,
            label=r"Nominal $\alpha = 5\%$", zorder=2)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(alphas)
ax2.set_ylim(0, 6.8)
ax2.grid(True, axis='y', linestyle='--', alpha=0.35)

# Place Panel C legend inside clear vertical band between y=1.2 and y=4.5
ax2.legend(loc='center right', frameon=True, framealpha=0.92, fontsize=7.2, edgecolor='#e2e8f0')

# Shared Legend on top of Panels A & B
handles_common = [
    plt.Line2D([0], [0], color='#0f172a', lw=1.5, linestyle='--', marker='D', markersize=4, label=r"Expected $E[R]$"),
    plt.Line2D([0], [0], color='#2563eb', lw=1.6, marker='o', markersize=4, label='RhizAeon'),
    plt.Line2D([0], [0], color='#2563eb', lw=1.2, linestyle=':', marker='s', markersize=3.8, label='3SEQ'),
    plt.Line2D([0], [0], color=colors[10.0], lw=2.2, label=r"$\theta=10$"),
    plt.Line2D([0], [0], color=colors[50.0], lw=2.2, label=r"$\theta=50$"),
    plt.Line2D([0], [0], color=colors[100.0], lw=2.2, label=r"$\theta=100$"),
    plt.Line2D([0], [0], color=colors[200.0], lw=2.2, label=r"$\theta=200$"),
]

fig.legend(handles=handles_common, loc='upper left', bbox_to_anchor=(0.065, 0.98),
           ncol=7, frameon=False, fontsize=7.5, columnspacing=1.1, handletextpad=0.4)

plt.savefig(OUT_PDF, bbox_inches='tight')
plt.savefig(OUT_PNG, bbox_inches='tight', dpi=300)
print(f"Successfully generated:\n  {OUT_PDF}\n  {OUT_PNG}")
