#!/usr/bin/env python3
"""
Publication-Quality Row 2: Scenario 1.2 Multi-Way Mosaic Swarms
Minimal iconic schematic matching Row 1 rhythm (minimal text) + High-Density Empirical Metrics.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

# Load authentic benchmark data
raw_df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms/multiway_raw_results.csv')
sum_df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms/multiway_summary.csv')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

# 11.8 x 3.1 in, 300 DPI: Widescreen row matching Row 1 exactly
fig = plt.figure(figsize=(11.8, 3.1), dpi=300)
gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 0.95, 0.95, 0.95], wspace=0.34,
                      left=0.035, right=0.985, top=0.78, bottom=0.18)

C_P1 = '#1f77b4'       # Parent 1 Home (Blue)
C_P2 = '#d95f02'       # Parent 2 Donor (Orange)
C_P3 = '#2ca02c'       # Parent 3 Donor (Green)
C_REC = '#0088cc'      # Recombinant (Cyan/Blue)
colors_div = {
    0.02: '#4dbbd5',   # light cyan-blue
    0.05: '#1f77b4',   # mid blue
    0.10: '#002b66',   # dark navy
}

# ==============================================================================
# Helper for Clean Axes (Matching Row 1 Exactly)
# ==============================================================================
def clean_axis(ax, title, xlabel, ylabel, ylim=None):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, linestyle='--', alpha=0.35, color='#dddddd', zorder=1)
    ax.set_title(title, fontsize=8.5, fontweight='bold', pad=6)
    ax.set_xlabel(xlabel, fontsize=7.5, labelpad=4)
    ax.set_ylabel(ylabel, fontsize=7.5, labelpad=4)
    ax.tick_params(labelsize=6.8, length=3.0, width=0.7, pad=2)
    if ylim is not None:
        ax.set_ylim(ylim)

# ==============================================================================
# 1. Column 1: Minimal Iconic Schematic (Matching Row 1 Rhythm, Minimal Text)
# ==============================================================================
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.axis('off')

# Title (Zero narrative subtitle, exactly like Row 1)
ax0.text(0.0, 1.22, 'B  Multi-Way Mosaic Swarms', fontsize=9.5, fontweight='bold', va='top')

# Upper Portion: Reticulation Tree & Swarm Expansion Fan
# Lineage branches
ax0.plot([0.04, 0.16], [0.92, 0.92], color=C_P1, lw=1.6)
ax0.text(0.05, 0.94, r'$P_1$', fontsize=7.2, fontweight='bold', color=C_P1, va='bottom')

ax0.plot([0.04, 0.16], [0.77, 0.77], color=C_P2, lw=1.6)
ax0.text(0.05, 0.79, r'$P_2$', fontsize=7.2, fontweight='bold', color=C_P2, va='bottom')

ax0.plot([0.04, 0.16], [0.62, 0.62], color=C_P3, lw=1.6)
ax0.text(0.05, 0.64, r'$P_3$', fontsize=7.2, fontweight='bold', color=C_P3, va='bottom')

# Background taxon count label (matching N_bg=40 in 1.1)
ax0.text(0.20, 0.92, r'$N_{\mathrm{bg}}=45$', va='center', fontsize=6.8, color='#555555')

# Dual transfers into R_anc at (0.44, 0.77)
ax0.plot([0.16, 0.44], [0.92, 0.77], color=C_P1, lw=1.4)
ax0.annotate('', xy=(0.44, 0.77), xytext=(0.16, 0.77),
             arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.3, ls='--'))
ax0.annotate('', xy=(0.44, 0.77), xytext=(0.16, 0.62),
             arrowprops=dict(arrowstyle="->", color=C_P3, lw=1.3, ls='--'))

# Reticulation Node Dot
ax0.plot(0.44, 0.77, 'o', color='#222222', markersize=4.0, zorder=5)
ax0.text(0.42, 0.83, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

# Clonal expansion fan to contemporary isolates
fan_y = [0.91, 0.82, 0.72, 0.63]
for yt in fan_y:
    ax0.plot([0.44, 0.65], [0.77, yt], color=C_REC, lw=1.1)
    ax0.plot(0.65, yt, 'o', color=C_REC, markersize=3.0)

# Minimal swarm label (identical to R_1 ... R_K in 1.1)
ax0.text(0.70, 0.77, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# --- Lower Portion: 3 Minimal Chromosome Tracks ---
h_bar = 0.075
w_total = 0.92
x0 = 0.04

# 1. Dual Cassette
y1 = 0.36
ax0.text(x0, y1 + h_bar + 0.015, 'Dual Cassette', fontsize=6.5, color='#333333', fontweight='bold', va='bottom')
w_seg = w_total * 0.20
ax0.add_patch(patches.Rectangle((x0, y1), w_seg, h_bar, facecolor=C_P1))
ax0.add_patch(patches.Rectangle((x0 + w_seg, y1), w_seg, h_bar, facecolor=C_P2))
ax0.add_patch(patches.Rectangle((x0 + 2*w_seg, y1), w_seg, h_bar, facecolor=C_P1))
ax0.add_patch(patches.Rectangle((x0 + 3*w_seg, y1), w_seg, h_bar, facecolor=C_P3))
ax0.add_patch(patches.Rectangle((x0 + 4*w_seg, y1), w_seg, h_bar, facecolor=C_P1))
ax0.add_patch(patches.Rectangle((x0, y1), w_total, h_bar, fill=False, edgecolor='#222222', lw=0.6))
ax0.text(x0 + 0.5*w_seg, y1 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + 1.5*w_seg, y1 + h_bar/2, r'$P_2$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + 2.5*w_seg, y1 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + 3.5*w_seg, y1 + h_bar/2, r'$P_3$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + 4.5*w_seg, y1 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')

# 2. Nested Russian Doll
y2 = 0.19
ax0.text(x0, y2 + h_bar + 0.015, 'Nested Doll', fontsize=6.5, color='#333333', fontweight='bold', va='bottom')
w_p1 = w_total * 0.24
w_p2 = w_total * 0.16
w_p3 = w_total * 0.20
ax0.add_patch(patches.Rectangle((x0, y2), w_p1, h_bar, facecolor=C_P1))
ax0.add_patch(patches.Rectangle((x0 + w_p1, y2), w_p2, h_bar, facecolor=C_P2))
ax0.add_patch(patches.Rectangle((x0 + w_p1 + w_p2, y2), w_p3, h_bar, facecolor=C_P3))
ax0.add_patch(patches.Rectangle((x0 + w_p1 + w_p2 + w_p3, y2), w_p2, h_bar, facecolor=C_P2))
ax0.add_patch(patches.Rectangle((x0 + w_p1 + 2*w_p2 + w_p3, y2), w_p1, h_bar, facecolor=C_P1))
ax0.add_patch(patches.Rectangle((x0, y2), w_total, h_bar, fill=False, edgecolor='#222222', lw=0.6))
ax0.text(x0 + w_p1/2, y2 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + w_p1 + w_p2/2, y2 + h_bar/2, r'$P_2$', color='white', fontsize=4.8, fontweight='bold', ha='center', va='center')
ax0.text(x0 + w_p1 + w_p2 + w_p3/2, y2 + h_bar/2, r'$P_3$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + w_p1 + 1.5*w_p2 + w_p3, y2 + h_bar/2, r'$P_2$', color='white', fontsize=4.8, fontweight='bold', ha='center', va='center')
ax0.text(x0 + w_p1 + 2*w_p2 + 1.5*w_p3, y2 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')

# 3. Triple Crossover
y3 = 0.02
ax0.text(x0, y3 + h_bar + 0.015, 'Triple Crossover', fontsize=6.5, color='#333333', fontweight='bold', va='bottom')
w_c1 = w_total * 0.33
w_c2 = w_total * 0.33
w_c3 = w_total * 0.34
ax0.add_patch(patches.Rectangle((x0, y3), w_c1, h_bar, facecolor=C_P1))
ax0.add_patch(patches.Rectangle((x0 + w_c1, y3), w_c2, h_bar, facecolor=C_P2))
ax0.add_patch(patches.Rectangle((x0 + w_c1 + w_c2, y3), w_c3, h_bar, facecolor=C_P3))
ax0.add_patch(patches.Rectangle((x0, y3), w_total, h_bar, fill=False, edgecolor='#222222', lw=0.6))
ax0.text(x0 + w_c1/2, y3 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + 1.5*w_c1, y3 + h_bar/2, r'$P_2$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax0.text(x0 + 2.5*w_c1, y3 + h_bar/2, r'$P_3$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')


# ==============================================================================
# 2. Column 2: Donors Correctly Resolved (Solid Fills Only, Zero Hatching)
# ==============================================================================
ax1 = fig.add_subplot(gs[0, 1])
topos = ['dual_cassette', 'nested_russian_doll', 'triple_crossover']
topo_names = ['Dual\nCassette', 'Nested\nDoll', 'Triple\nCrossover']
x_t = np.arange(len(topos))
bar_w = 0.35

ax1.axhline(2.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

r_donors_mean = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['topology'] == t)]['donors_resolved'].mean() for t in topos]
r_donors_std  = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['topology'] == t)]['donors_resolved'].std() for t in topos]
s_donors_mean = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['topology'] == t)]['donors_resolved'].mean() for t in topos]
s_donors_std  = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['topology'] == t)]['donors_resolved'].std() for t in topos]

# Pure Solid Fills (ZERO HATCHING)
# RhizAeon: Solid Deep Navy (#002b66)
# 3SEQ: Solid Sky Blue (#85b3d1) with crisp dark border
ax1.bar(x_t - bar_w/2, r_donors_mean, bar_w, yerr=r_donors_std, capsize=3,
        label='RhizAeon', color='#002b66', edgecolor='none', zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))
ax1.bar(x_t + bar_w/2, s_donors_mean, bar_w, yerr=s_donors_std, capsize=3,
        label='3SEQ', color='#85b3d1', edgecolor='#002b66', lw=0.9, zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))

clean_axis(ax1, 'Donors Correctly Resolved', 'Mosaic Topology', 'Donors Resolved', ylim=(0, 2.45))
ax1.set_xticks(x_t)
ax1.set_xticklabels(topo_names, fontsize=6.8)


# ==============================================================================
# 3. Column 3: Breakpoint Overlap Jaccard (Solid Fills Only, Zero Hatching)
# ==============================================================================
ax2 = fig.add_subplot(gs[0, 2])
ax2.axhline(1.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

r_jacc_mean = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['topology'] == t)]['mean_jaccard'].mean() for t in topos]
r_jacc_std  = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['topology'] == t)]['mean_jaccard'].std() for t in topos]
s_jacc_mean = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['topology'] == t)]['mean_jaccard'].mean() for t in topos]
s_jacc_std  = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['topology'] == t)]['mean_jaccard'].std() for t in topos]

# Pure Solid Fills (ZERO HATCHING)
ax2.bar(x_t - bar_w/2, r_jacc_mean, bar_w, yerr=r_jacc_std, capsize=3,
        color='#002b66', edgecolor='none', zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))
ax2.bar(x_t + bar_w/2, s_jacc_mean, bar_w, yerr=s_jacc_std, capsize=3,
        color='#85b3d1', edgecolor='#002b66', lw=0.9, zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))

clean_axis(ax2, r'Breakpoint Overlap ($J_{\mathrm{SNP}}$)', 'Mosaic Topology', 'Jaccard Index', ylim=(0, 1.15))
ax2.set_xticks(x_t)
ax2.set_xticklabels(topo_names, fontsize=6.8)


# ==============================================================================
# 4. Column 4: Inferred Events (\hat{K}) vs Swarm Size (K) with Authentic IQR Bands
# ==============================================================================
ax3 = fig.add_subplot(gs[0, 3])
clean_axis(ax3, r'Inferred Events ($\hat{K}$)', r'Descendant Swarm Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-2, 60))
ax3.axhline(2.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

swarm_sizes = [1, 5, 10]
ax3.set_xticks(swarm_sizes)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_r = sum_df[(sum_df['tool'] == 'RhizAeon') & (np.isclose(sum_df['d_parent'], d))].groupby('K')['mean_unique_events'].mean().reset_index()
    sub_3 = sum_df[(sum_df['tool'] == '3SEQ') & (np.isclose(sum_df['d_parent'], d))].groupby('K')['mean_unique_events'].mean().reset_index()
    
    # 3SEQ IQR dispersion
    raw_3 = raw_df[(raw_df['tool'] == '3SEQ') & (np.isclose(raw_df['d_parent'], d))]
    q25 = raw_3.groupby('K')['unique_events'].quantile(0.25)
    q75 = raw_3.groupby('K')['unique_events'].quantile(0.75)
    ax3.fill_between(sub_3['K'], q25, q75, color=c, alpha=0.15, zorder=2)
    
    # Curves
    ax3.plot(sub_3['K'], sub_3['mean_unique_events'], marker='s', mfc='none', mec=c, mew=1.3,
             ms=4.8, ls=':', lw=1.2, color=c, zorder=3)
    ax3.plot(sub_r['K'], sub_r['mean_unique_events'], marker='o', mfc=c, mec='white', mew=0.5,
             ms=5.0, ls='-', lw=1.6, color=c, zorder=4)


# ==============================================================================
# Unified Top-Hoisted Legend (Identical to Row 1)
# ==============================================================================
legend_handles = [
    mlines.Line2D([], [], color='#002b66', marker='o', lw=1.5, markersize=4.5, label='RhizAeon'),
    mlines.Line2D([], [], color='#002b66', marker='s', lw=1.3, linestyle=':', markersize=4.8,
                  markerfacecolor='none', markeredgewidth=1.3, label='3SEQ'),
    mlines.Line2D([], [], color='#222222', linestyle='--', lw=1.2, label='Ground Truth'),
    patches.Patch(facecolor=colors_div[0.02], label=r'$d = 0.02$'),
    patches.Patch(facecolor=colors_div[0.05], label=r'$d = 0.05$'),
    patches.Patch(facecolor=colors_div[0.10], label=r'$d = 0.10$'),
]

fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.58, 0.98),
           ncol=6, frameon=False, fontsize=7.2, handletextpad=0.4, columnspacing=1.3)

# Save to local viz directory and manuscript figures directory
out_png = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_panel_row2_scenario12.png'
out_pdf = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_panel_row2_scenario12.pdf'
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)

paper_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_panel_row2_scenario12.png'
paper_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_panel_row2_scenario12.pdf'
plt.savefig(paper_png, dpi=300)
plt.savefig(paper_pdf)
plt.close()
print("Row 2 successfully updated with minimal text schematic.")
