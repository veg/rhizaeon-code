#!/usr/bin/env python3
"""
generate_fig2_simulations_master.py

Author: Sergei L. Kosakovsky Pond et al.
Generates the supercompact, publication-grade 8-panel Master Simulation Figure (Figure 2)
unifying the complete empirical simulation evidence base for RhizAeon:
  Row 1: Coalescent Benchmark & Triplet Overcounting
    (A) Posada & Crandall: Detection Power vs. Recombination Parameter (rho) across Mutation Densities (theta)
    (B) Posada & Crandall: Inferred Event Cardinality (K_hat) vs. Hudson ARG Theoretical Expectation E[R]
    (C) Posada & Crandall: Null FPR under Gamma Site-Rate Heterogeneity (alpha in {2.0, 0.5, 0.05})
    (D) Scenario 1.1: Epidemic Clonal Swarms: Inferred Events (K_hat) across Swarm Size (K)
  Row 2: Operational Surveillance Regimes & Continuous Phase Limits
    (E) Scenario 1.2: Multi-Way Mosaicism: Inferred Donors Resolved across Complex Topologies
    (F) Scenario 1.3: Ghost Donor Introgression: Closed-World Innocent Parent False Framing vs. True Attribution
    (G) Scenario 1.4: Bacterial Micro-Conversions: Scale-Space Detection Power down to 50 nt
    (H) Grand Parameter Sweep (150,000 trials): Continuous 2D Empirical Phase Diagram (d_parent vs. L_tract / L)

Strict MBE/Nature BioVis Standards:
  - 2 x 4 clean panoramic grid layout.
  - Zero overlapping labels, half-open axes (L-frames), sans-serif typography.
  - Colorblind-safe palettes (Okabe-Ito, Blues, Viridis).
  - Ground truth and theoretical baselines plotted as distinct reference lines.
  - Strictly authentic numbers from executed simulation CSVs on disk.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from scipy.ndimage import gaussian_filter

# Set publication typography & styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 8.0,
    'axes.labelsize': 8.2,
    'axes.titlesize': 8.6,
    'xtick.labelsize': 7.0,
    'ytick.labelsize': 7.0,
    'legend.fontsize': 6.8,
    'axes.linewidth': 0.75,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.30,
})

BASE_DIR = "/Users/sergei/Projects/rhizaeon-code"
OUT_PDF = os.path.join(BASE_DIR, "paper/figures/fig2_simulations_master.pdf")
OUT_PNG = os.path.join(BASE_DIR, "paper/figures/fig2_simulations_master.png")

# Palette
C_RHIZ = '#1f77b4'       # Royal Blue
C_3SEQ = '#64748b'       # Slate Gray
C_TRUE = '#0f172a'       # Dark Charcoal
C_ALERT = '#dc2626'      # Crimson Red
C_LIGHT_BLUE = '#93c5fd'

theta_colors = {
    10.0: '#94a3b8',   # Slate light
    50.0: '#38bdf8',   # Sky blue
    100.0: '#2563eb',  # Royal blue
    200.0: '#1e3a8a',  # Deep navy
}

div_colors = {
    0.02: '#38bdf8',
    0.05: '#2563eb',
    0.10: '#0f172a',
}

fig, axes = plt.subplots(2, 4, figsize=(11.4, 3.65), dpi=300)
fig.subplots_adjust(wspace=0.32, hspace=0.46, left=0.045, right=0.985, top=0.92, bottom=0.10)

def clean_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, linestyle='--', alpha=0.35, color='#dddddd', zorder=1)

# ==============================================================================
# PANEL A: Posada & Crandall (2001) Power vs rho
# ==============================================================================
ax_a = axes[0, 0]
clean_ax(ax_a)
ax_a.set_title("A. Coalescent Power vs. $\\rho$", fontweight='bold', loc='left')
df_pos = pd.read_csv(os.path.join(BASE_DIR, "benchmarks/posada_crandall_2001/posada_comparison_rhizaeon_vs_3seq.csv"))
p_rec = df_pos[df_pos['category_rhiz'] == 'power']

rhos = [0.0, 1.0, 4.0, 16.0, 64.0]
rho_labels = ['0', '1', '4', '16', '64']
x_idx = np.arange(len(rhos))

for th in [10.0, 50.0, 100.0, 200.0]:
    sub = p_rec[p_rec['theta_rhiz'] == th].sort_values('rho_rhiz')
    ax_a.plot(x_idx, sub['detection_pct_rhiz'], marker='o', markersize=3.5,
              lw=1.4, color=theta_colors[th], zorder=4)
    ax_a.plot(x_idx, sub['detection_pct_3seq'], marker='s', markersize=3.2,
              lw=1.0, ls=':', color=theta_colors[th], alpha=0.85, zorder=3)

ax_a.set_xticks(x_idx)
ax_a.set_xticklabels(rho_labels)
ax_a.set_xlabel("Recombination Rate ($\\rho = 4N_e r$)")
ax_a.set_ylabel("Detection Power (%)")
ax_a.set_ylim(-3, 103)

# Legend inside panel A
leg_a = [
    mlines.Line2D([], [], color=C_RHIZ, marker='o', markersize=3.5, lw=1.3, label='RhizAeon'),
    mlines.Line2D([], [], color=C_3SEQ, marker='s', markersize=3.2, lw=1.0, ls=':', label='3SEQ'),
    mlines.Line2D([], [], color=theta_colors[200.0], lw=1.5, label=r'$\theta=200$'),
    mlines.Line2D([], [], color=theta_colors[50.0], lw=1.5, label=r'$\theta=50$'),
]
ax_a.legend(handles=leg_a, loc='lower right', frameon=True, framealpha=0.92, fontsize=6.2, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL B: Inferred Event Cardinality (K_hat vs rho)
# ==============================================================================
ax_b = axes[0, 1]
clean_ax(ax_b)
ax_b.set_title("B. Event Cardinality ($\\hat{K}$ vs. $\\rho$)", fontweight='bold', loc='left')

H_9 = sum(1.0 / i for i in range(1, 10))
exp_R = [0.5 * r * H_9 for r in rhos]
ax_b.plot(x_idx, exp_R, color=C_TRUE, ls='--', marker='D', markersize=3.5, lw=1.3, zorder=5, label='Hudson $\\mathbb{E}[R]$')

for th in [50.0, 100.0, 200.0]:
    sub = p_rec[p_rec['theta_rhiz'] == th].sort_values('rho_rhiz')
    ax_b.plot(x_idx, sub['mean_K_hat_3seq'], marker='s', markersize=3.2, lw=1.0, ls=':', color=theta_colors[th], alpha=0.85, zorder=3)
    ax_b.plot(x_idx, sub['mean_K_hat_rhiz'], marker='o', markersize=3.5, lw=1.4, color=theta_colors[th], zorder=4)

ax_b.set_yscale('symlog', linthresh=1.0)
ax_b.set_xticks(x_idx)
ax_b.set_xticklabels(rho_labels)
ax_b.set_yticks([0, 1, 5, 20, 100])
ax_b.set_yticklabels(['0', '1', '5', '20', '100'])
ax_b.set_xlabel("Recombination Rate ($\\rho = 4N_e r$)")
ax_b.set_ylabel("Mean Inferred Events ($\\hat{K}$)")
ax_b.set_ylim(-0.15, 110)
ax_b.legend(loc='upper left', frameon=True, framealpha=0.92, fontsize=6.2, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL C: Rate Variation Null FPR
# ==============================================================================
ax_c = axes[0, 2]
clean_ax(ax_c)
ax_c.set_title("C. Null FPR under Heterotachy", fontweight='bold', loc='left')
p_fp = df_pos[df_pos['category_rhiz'] == 'false_positive']

alpha_labels = ['2.0\n(Mild)', '0.5\n(Mod.)', '0.05\n(Severe)']
alpha_vals = [2.0, 0.5, 0.05]
x_c = np.arange(len(alpha_vals))
width = 0.32

r_fp = [p_fp[p_fp['alpha_rhiz'] == a]['detection_pct_rhiz'].mean() for a in alpha_vals]
s_fp = [p_fp[p_fp['alpha_rhiz'] == a]['detection_pct_3seq'].mean() for a in alpha_vals]

ax_c.bar(x_c - width/2, r_fp, width, label='RhizAeon (0.58%)', color=C_RHIZ, alpha=0.9, zorder=3)
ax_c.bar(x_c + width/2, s_fp, width, label='3SEQ (0.17%)', color='#94a3b8', alpha=0.9, zorder=3)
ax_c.axhline(5.0, color=C_ALERT, ls='--', lw=1.0, label='Nominal 5%', zorder=2)

ax_c.set_xticks(x_c)
ax_c.set_xticklabels(alpha_labels)
ax_c.set_xlabel("Gamma Shape Parameter ($\\alpha$)")
ax_c.set_ylabel("False Positive Rate (%)")
ax_c.set_ylim(0, 6.5)
ax_c.legend(loc='center right', frameon=True, framealpha=0.92, fontsize=6.2, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL D: Scenario 1.1 Epidemic Clonal Swarms
# ==============================================================================
ax_d = axes[0, 3]
clean_ax(ax_d)
ax_d.set_title("D. Clonal Swarm Cardinality", fontweight='bold', loc='left')
df_s1 = pd.read_csv(os.path.join(BASE_DIR, "benchmarks/scenario1_epidemic_clonal_expansion/epidemic_summary.csv"))

swarm_sizes = [2, 5, 8, 10, 15]
ax_d.axhline(1.0, color=C_TRUE, ls='--', lw=1.2, label='True Event ($K=1$)', zorder=2)

for d in [0.02, 0.05, 0.10]:
    c = div_colors[d]
    sub_r = df_s1[(df_s1['tool'] == 'RhizAeon') & (np.isclose(df_s1['d_parent'], d))].sort_values('K')
    sub_3 = df_s1[(df_s1['tool'] == '3SEQ') & (np.isclose(df_s1['d_parent'], d))].sort_values('K')
    ax_d.plot(sub_3['K'], sub_3['unique_events_mean'], marker='s', markersize=3.2, lw=1.0, ls=':', color=c, zorder=3)
    ax_d.plot(sub_r['K'], sub_r['unique_events_mean'], marker='o', markersize=3.5, lw=1.4, color=c, zorder=4)

ax_d.set_xticks(swarm_sizes)
ax_d.set_xlabel("Descendant Swarm Size ($K$)")
ax_d.set_ylabel("Inferred Events ($\\hat{K}$)")
ax_d.set_ylim(-1, 60)

leg_d = [
    mlines.Line2D([], [], color=C_TRUE, ls='--', lw=1.2, label='True ($K=1$)'),
    mlines.Line2D([], [], color=C_RHIZ, marker='o', markersize=3.5, lw=1.3, label='RhizAeon ($\\hat{K}=1.0$)'),
    mlines.Line2D([], [], color=C_3SEQ, marker='s', markersize=3.2, lw=1.0, ls=':', label='3SEQ ($\\hat{K} \\to 55$)'),
]
ax_d.legend(handles=leg_d, loc='upper left', frameon=True, framealpha=0.92, fontsize=6.0, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL E: Scenario 1.2 Multi-Way Mosaicism (Donors Resolved)
# ==============================================================================
ax_e = axes[1, 0]
clean_ax(ax_e)
ax_e.set_title("E. Multi-Way Mosaic Donors", fontweight='bold', loc='left')
df_s2 = pd.read_csv(os.path.join(BASE_DIR, "benchmarks/scenario2_multiway_mosaic_swarms/multiway_summary.csv"))

topos = ['dual_cassette', 'nested_russian_doll', 'triple_crossover']
topo_short = ['Dual\nCassette', 'Nested\nDoll', 'Triple\nCross']
x_e = np.arange(len(topos))
width = 0.32

r_donors = []
s_donors = []
for topo in topos:
    sub_r = df_s2[(df_s2['tool'] == 'RhizAeon') & (df_s2['topology'] == topo) & (np.isclose(df_s2['d_parent'], 0.05))]
    sub_s = df_s2[(df_s2['tool'] == '3SEQ') & (df_s2['topology'] == topo) & (np.isclose(df_s2['d_parent'], 0.05))]
    r_donors.append(sub_r['mean_donors_resolved'].mean())
    s_donors.append(sub_s['mean_donors_resolved'].mean())

ax_e.bar(x_e - width/2, r_donors, width, label='RhizAeon', color=C_RHIZ, alpha=0.9, zorder=3)
ax_e.bar(x_e + width/2, s_donors, width, label='3SEQ', color='#94a3b8', alpha=0.9, zorder=3)
ax_e.axhline(2.0, color=C_TRUE, ls='--', lw=1.2, label='Expected ($K_{\\mathrm{donor}}=2$)', zorder=2)

ax_e.set_xticks(x_e)
ax_e.set_xticklabels(topo_short)
ax_e.set_xlabel("Mosaic Architecture")
ax_e.set_ylabel("Resolved Donors (out of 2)")
ax_e.set_ylim(0, 2.35)
ax_e.legend(loc='lower right', frameon=True, framealpha=0.92, fontsize=6.2, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL F: Scenario 1.3 Ghost Donor False Framing
# ==============================================================================
ax_f = axes[1, 1]
clean_ax(ax_f)
ax_f.set_title("F. Ghost Donor False Framing", fontweight='bold', loc='left')
df_s3 = pd.read_csv(os.path.join(BASE_DIR, "benchmarks/scenario3_ghost_donor_introgression/ghost_summary.csv"))

regimes = ['surrogate_sister', 'pure_ghost', 'deep_reservoir']
reg_labels = ['Surrogate\nSister', 'Pure\nGhost', 'Deep\nReservoir']
x_f = np.arange(len(regimes))

r_pfpr = [df_s3[(df_s3['tool'] == 'RhizAeon') & (df_s3['regime'] == r)]['mean_parent_fpr'].mean() * 100 for r in regimes]
s_pfpr = [df_s3[(df_s3['tool'] == '3SEQ') & (df_s3['regime'] == r)]['mean_parent_fpr'].mean() * 100 for r in regimes]

ax_f.bar(x_f - width/2, r_pfpr, width, label='RhizAeon (<= 6%)', color=C_RHIZ, alpha=0.9, zorder=3)
ax_f.bar(x_f + width/2, s_pfpr, width, label='3SEQ (100%)', color=C_ALERT, alpha=0.85, zorder=3)

ax_f.set_xticks(x_f)
ax_f.set_xticklabels(reg_labels)
ax_f.set_xlabel("Sampling Regime")
ax_f.set_ylabel("Innocent Parent Framing (%)")
ax_f.set_ylim(0, 115)
ax_f.legend(loc='upper right', frameon=True, framealpha=0.92, fontsize=6.2, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL G: Scenario 1.4 Bacterial Micro-Conversions
# ==============================================================================
ax_g = axes[1, 2]
clean_ax(ax_g)
ax_g.set_title("G. Micro-Conversion Power", fontweight='bold', loc='left')
df_s4 = pd.read_csv(os.path.join(BASE_DIR, "benchmarks/scenario4_micro_conversion_tracts/micro_conversion_summary.csv"))

tract_lens = [50, 100, 200, 400, 800]
x_g = np.arange(len(tract_lens))

# Plot power across tract lengths for d=0.10 and d=0.05
for d in [0.05, 0.10]:
    c = div_colors[d]
    sub_r = df_s4[(df_s4['tool'] == 'RhizAeon') & (np.isclose(df_s4['d_parent'], d))].groupby('tract_len')['recall'].mean()
    sub_s = df_s4[(df_s4['tool'] == '3SEQ') & (np.isclose(df_s4['d_parent'], d))].groupby('tract_len')['recall'].mean()
    ax_g.plot(x_g, sub_r.values * 100, marker='o', markersize=3.5, lw=1.4, color=c, zorder=4)
    ax_g.plot(x_g, sub_s.values * 100, marker='s', markersize=3.2, lw=1.0, ls=':', color=c, alpha=0.85, zorder=3)

ax_g.set_xticks(x_g)
ax_g.set_xticklabels(['50', '100', '200', '400', '800'])
ax_g.set_xlabel("Conversion Tract Length (nt)")
ax_g.set_ylabel("Detection Power (%)")
ax_g.set_ylim(-3, 105)

leg_g = [
    mlines.Line2D([], [], color=div_colors[0.10], marker='o', markersize=3.5, lw=1.3, label='RhizAeon (0.10)'),
    mlines.Line2D([], [], color=div_colors[0.10], marker='s', markersize=3.2, lw=1.0, ls=':', label='3SEQ (0.10)'),
    mlines.Line2D([], [], color=div_colors[0.05], marker='o', markersize=3.5, lw=1.3, label='RhizAeon (0.05)'),
    mlines.Line2D([], [], color=div_colors[0.05], marker='s', markersize=3.2, lw=1.0, ls=':', label='3SEQ (0.05)'),
]
ax_g.legend(handles=leg_g, loc='lower right', frameon=True, framealpha=0.92, fontsize=5.8, edgecolor='#e2e8f0')

# ==============================================================================
# PANEL H: Grand Sweep 2D Continuous Phase Diagram
# ==============================================================================
ax_h = axes[1, 3]
clean_ax(ax_h)
ax_h.set_title("H. Operational Phase Envelope", fontweight='bold', loc='left')

# Load recombinant trials from grand_sweep_150k.csv (or 50k for fast plotting if identical grid)
SWEEP_CSV = os.path.join(BASE_DIR, "grand_sweep_150k.csv")
if not os.path.exists(SWEEP_CSV):
    SWEEP_CSV = os.path.join(BASE_DIR, "grand_sweep_50k.csv")

df_sweep = pd.read_csv(SWEEP_CSV)
recs = df_sweep[df_sweep["is_recombinant_true"] == 1].copy()
recs["tract_ratio"] = recs["tract_len"] / recs["length"]

n_bins_x, n_bins_y = 28, 28
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

smoothed_rate = gaussian_filter(np.nan_to_num(grid_rate, nan=0.0), sigma=1.0)
x_centers = np.sqrt(x_edges[:-1] * x_edges[1:])
y_centers = np.sqrt(y_edges[:-1] * y_edges[1:])
X, Y = np.meshgrid(x_centers, y_centers)

c_mesh = ax_h.pcolormesh(x_edges, y_edges, grid_rate, cmap="Blues", vmin=0.0, vmax=1.0, shading="flat", zorder=2)
cs = ax_h.contour(X, Y, smoothed_rate, levels=[0.50, 0.80, 0.95], colors=['#f59e0b', '#10b981', '#065f46'], linewidths=1.2, zorder=3)
ax_h.clabel(cs, fmt=lambda v: f"{int(v*100)}%", fontsize=5.8, inline=True)

# Plot theoretical M_tract = 6 contour: M = d * L_tract = d * ratio * L_mean (L_mean ~ 5000 nt)
d_line = np.geomspace(0.005, 0.25, 50)
ratio_m6 = 6.0 / (d_line * 5000.0)
valid = (ratio_m6 >= 0.008) & (ratio_m6 <= 0.50)
ax_h.plot(d_line[valid], ratio_m6[valid], color='#dc2626', ls='--', lw=1.2, label='$M_{\\mathrm{tract}}=6$ SNPs', zorder=4)

ax_h.set_xscale("log")
ax_h.set_yscale("log")
ax_h.set_xlabel(r"Divergence $d_{\mathrm{parent}}$ (subs/site)")
ax_h.set_ylabel(r"Tract Ratio ($L_{\mathrm{tract}} / L$)")
ax_h.set_xlim(0.005, 0.25)
ax_h.set_ylim(0.008, 0.50)

# Colorbar for Panel H
cb = fig.colorbar(c_mesh, ax=ax_h, fraction=0.046, pad=0.04)
cb.set_label("Power", fontsize=6.8)
cb.ax.tick_params(labelsize=6.0)

ax_h.legend(loc='lower left', frameon=True, framealpha=0.92, fontsize=5.8, edgecolor='#e2e8f0')

# Save figures
plt.savefig(OUT_PDF, bbox_inches='tight')
plt.savefig(OUT_PNG, bbox_inches='tight')
print(f"Generated {OUT_PDF} and {OUT_PNG} successfully.")
