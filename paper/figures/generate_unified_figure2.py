#!/usr/bin/env python3
"""
Master Composer for Figure 2: Specialized Genomic & Epidemiological Simulation Scenarios.
Unifies Scenarios 1.1 - 1.4 into a cohesive 4-row publication-quality figure.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

# ==============================================================================
# Global Matplotlib Typography & Style
# ==============================================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

C_P1 = '#1f77b4'       # Home / Parent 1 (Blue)
C_P2 = '#d95f02'       # Donor / Parent 2 (Orange)
C_P3 = '#2ca02c'       # Donor / Parent 3 (Green)
C_SURR = '#2ca02c'     # Surrogate Sister (Green)
C_GHOST = '#7570b3'    # Unsampled Ghost Donor (Purple / Dashed)
C_REC = '#0088cc'      # Recombinant swarm (Cyan/Blue)

colors_div = {
    0.02: '#4dbbd5',   # light cyan-blue
    0.05: '#1f77b4',   # mid blue
    0.10: '#002b66',   # dark navy
}

# 11.8 x 12.0 in, 300 DPI: Unified 4-row Master Figure
fig = plt.figure(figsize=(11.8, 12.6), dpi=300)
gs = fig.add_gridspec(4, 4, width_ratios=[1.15, 0.95, 0.95, 0.95],
                      height_ratios=[1.0, 1.0, 1.0, 1.0],
                      wspace=0.34, hspace=0.42,
                      left=0.035, right=0.985, top=0.938, bottom=0.045)

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
# ROW 1 (A): SCENARIO 1.1 EPIDEMIC CLONAL SWARMS
# ==============================================================================
df1_raw = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/epidemic_raw_results.csv')
df1_sum = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/epidemic_summary.csv')

# --- Row 1, Col 1: Minimal Tree Schematic ---
ax1_0 = fig.add_subplot(gs[0, 0])
ax1_0.set_xlim(0, 1)
ax1_0.set_ylim(0, 1)
ax1_0.axis('off')
ax1_0.text(0.0, 1.22, 'A  Epidemic Clonal Clusters', fontsize=9.5, fontweight='bold', va='top')

# P1 Home Clade (Top)
ax1_0.plot([0.05, 0.16], [0.84, 0.84], color=C_P1, lw=1.6)
ax1_0.plot([0.16, 0.16], [0.72, 0.94], color=C_P1, lw=1.6)
ax1_0.plot([0.16, 0.36], [0.94, 0.94], color=C_P1, lw=1.3)
ax1_0.plot([0.16, 0.36], [0.72, 0.72], color=C_P1, lw=1.5)
ax1_0.text(0.39, 0.94, r'$N_{\mathrm{bg}}=40$', va='center', fontsize=6.8, color='#555555')
ax1_0.text(0.20, 0.77, r'$P_1$', va='bottom', fontsize=7.8, fontweight='bold', color=C_P1)

# P2 Donor Clade (Bottom)
ax1_0.plot([0.05, 0.16], [0.44, 0.44], color=C_P2, lw=1.6)
ax1_0.plot([0.16, 0.16], [0.32, 0.54], color=C_P2, lw=1.6)
ax1_0.plot([0.16, 0.36], [0.54, 0.54], color=C_P2, lw=1.5)
ax1_0.plot([0.16, 0.36], [0.32, 0.32], color=C_P2, lw=1.3)
ax1_0.text(0.20, 0.49, r'$P_2$', va='top', fontsize=7.8, fontweight='bold', color=C_P2)

# Reticulation transfer to R_anc
ax1_0.plot([0.36, 0.58], [0.72, 0.60], color=C_P1, lw=1.4)
ax1_0.annotate('', xy=(0.58, 0.60), xytext=(0.36, 0.54),
               arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.4, ls='--'))
ax1_0.plot(0.58, 0.60, 'o', color='#222222', markersize=4.2, zorder=5)
ax1_0.text(0.57, 0.67, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

# Clonal fan to contemporary isolates
for yt in [0.80, 0.67, 0.53, 0.40]:
    ax1_0.plot([0.58, 0.76], [0.60, yt], color=C_REC, lw=1.2)
    ax1_0.plot(0.76, yt, 'o', color=C_REC, markersize=3.2)
ax1_0.text(0.80, 0.60, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# Minimal Chromosome Track
bar_y = 0.08
bar_h = 0.13
ax1_0.add_patch(patches.Rectangle((0.05, bar_y), 0.27, bar_h, facecolor=C_P1, edgecolor='none'))
ax1_0.add_patch(patches.Rectangle((0.32, bar_y), 0.27, bar_h, facecolor=C_P2, edgecolor='none'))
ax1_0.add_patch(patches.Rectangle((0.59, bar_y), 0.36, bar_h, facecolor=C_P1, edgecolor='none'))
ax1_0.add_patch(patches.Rectangle((0.05, bar_y), 0.90, bar_h, fill=False, edgecolor='#222222', lw=0.7))
ax1_0.text(0.05 + 0.135, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')
ax1_0.text(0.32 + 0.135, bar_y + bar_h/2, r'$P_2$ (1.5 kb)', color='white', fontsize=6.2, fontweight='bold', ha='center', va='center')
ax1_0.text(0.59 + 0.18, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

ax1_0.plot([0.05, 0.05], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax1_0.plot([0.32, 0.32], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax1_0.plot([0.59, 0.59], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax1_0.plot([0.95, 0.95], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax1_0.text(0.05, bar_y - 0.038, '0', fontsize=5.5, ha='left', va='top', color='#555555')
ax1_0.text(0.32, bar_y - 0.038, '1.5k', fontsize=5.5, ha='center', va='top', color='#555555')
ax1_0.text(0.59, bar_y - 0.038, '3.0k', fontsize=5.5, ha='center', va='top', color='#555555')
ax1_0.text(0.95, bar_y - 0.038, '5.0 kb', fontsize=5.5, ha='right', va='top', color='#555555')

# --- Row 1, Col 2: Isolate Precision (%) with IQR ---
ax1_1 = fig.add_subplot(gs[0, 1])
clean_axis(ax1_1, 'Isolate Precision (%)', r'Descendant Cluster Size ($K$)', 'Precision (%)', ylim=(-3, 105))
ax1_1.axhline(100.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)
ks_1 = [2, 5, 8, 10, 15]
ax1_1.set_xticks(ks_1)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_r = df1_sum[(df1_sum['tool'] == 'RhizAeon') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    sub_3c = df1_sum[(df1_sum['tool'] == '3SEQ (Consolidated)') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    sub_3 = df1_sum[(df1_sum['tool'] == '3SEQ') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    
    raw_3 = df1_raw[(df1_raw['tool'] == '3SEQ') & (np.isclose(df1_raw['d_parent'], d))]
    q25 = raw_3.groupby('K')['precision'].quantile(0.25) * 100
    q75 = raw_3.groupby('K')['precision'].quantile(0.75) * 100
    ax1_1.fill_between(sub_3['K'], q25, q75, color=c, alpha=0.12, zorder=2)
    
    ax1_1.plot(sub_3['K'], sub_3['precision_mean'] * 100.0, marker='s', mfc='none', mec=c, mew=1.2,
               ms=4.4, ls=':', lw=1.1, color=c, zorder=3)
    ax1_1.plot(sub_3c['K'], sub_3c['precision_mean'] * 100.0, marker='^', mfc='white', mec=c, mew=1.2,
               ms=4.6, ls='--', lw=1.3, color=c, zorder=4)
    ax1_1.plot(sub_r['K'], sub_r['precision_mean'] * 100.0, marker='o', mfc=c, mec='white', mew=0.5,
               ms=4.8, ls='-', lw=1.6, color=c, zorder=5)

# --- Row 1, Col 3: Runtime Scalability (ms) ---
ax1_2 = fig.add_subplot(gs[0, 2])
clean_axis(ax1_2, 'Runtime Scalability', r'Descendant Cluster Size ($K$)', 'Execution Time (ms)', ylim=(0, 11500))
ax1_2.set_xticks(ks_1)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_r = df1_sum[(df1_sum['tool'] == 'RhizAeon') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    sub_3 = df1_sum[(df1_sum['tool'] == '3SEQ') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    
    ax1_2.plot(sub_3['K'], sub_3['time_ms_mean'], marker='s', mfc='none', mec=c, mew=1.2,
               ms=4.4, ls=':', lw=1.2, color=c, zorder=3)
    ax1_2.plot(sub_r['K'], sub_r['time_ms_mean'], marker='o', mfc=c, mec='white', mew=0.5,
               ms=4.8, ls='-', lw=1.6, color=c, zorder=4)

# --- Row 1, Col 4: Inferred Events (\hat{K}) vs K with IQR ---
ax1_3 = fig.add_subplot(gs[0, 3])
clean_axis(ax1_3, r'Inferred Events ($\hat{K}$)', r'Descendant Cluster Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-2, 60))
ax1_3.axhline(1.0, color='#222222', linestyle='-', linewidth=1.4, zorder=2)
ax1_3.set_xticks(ks_1)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_r = df1_sum[(df1_sum['tool'] == 'RhizAeon') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    sub_3c = df1_sum[(df1_sum['tool'] == '3SEQ (Consolidated)') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    sub_3 = df1_sum[(df1_sum['tool'] == '3SEQ') & (np.isclose(df1_sum['d_parent'], d))].sort_values('K')
    
    raw_3 = df1_raw[(df1_raw['tool'] == '3SEQ') & (np.isclose(df1_raw['d_parent'], d))]
    q25 = raw_3.groupby('K')['unique_events'].quantile(0.25)
    q75 = raw_3.groupby('K')['unique_events'].quantile(0.75)
    ax1_3.fill_between(sub_3['K'], q25, q75, color=c, alpha=0.12, zorder=2)
    
    ax1_3.plot(sub_3['K'], sub_3['unique_events_mean'], marker='s', mfc='none', mec=c, mew=1.2,
               ms=4.4, ls=':', lw=1.1, color=c, zorder=3)
    ax1_3.plot(sub_3c['K'], sub_3c['unique_events_mean'], marker='^', mfc='white', mec=c, mew=1.2,
               ms=4.6, ls='--', lw=1.3, color=c, zorder=4)
    ax1_3.plot(sub_r['K'], sub_r['unique_events_mean'], marker='o', mfc=c, mec='white', mew=0.5,
               ms=4.8, ls='-', lw=1.6, color=c, zorder=5)


# ==============================================================================
# ROW 2 (B): SCENARIO 1.2 MULTI-WAY MOSAIC CLUSTERS
# ==============================================================================
df2_raw = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms/multiway_raw_results.csv')
df2_sum = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms/multiway_summary.csv')

ax2_0 = fig.add_subplot(gs[1, 0])
ax2_0.set_xlim(0, 1)
ax2_0.set_ylim(0, 1)
ax2_0.axis('off')
ax2_0.text(0.0, 1.22, 'B  Multi-Way Mosaic Clusters', fontsize=9.5, fontweight='bold', va='top')

# Lineage branches
ax2_0.plot([0.04, 0.16], [0.92, 0.92], color=C_P1, lw=1.6)
ax2_0.text(0.05, 0.94, r'$P_1$', fontsize=7.2, fontweight='bold', color=C_P1, va='bottom')
ax2_0.plot([0.04, 0.16], [0.77, 0.77], color=C_P2, lw=1.6)
ax2_0.text(0.05, 0.79, r'$P_2$', fontsize=7.2, fontweight='bold', color=C_P2, va='bottom')
ax2_0.plot([0.04, 0.16], [0.62, 0.62], color=C_P3, lw=1.6)
ax2_0.text(0.05, 0.64, r'$P_3$', fontsize=7.2, fontweight='bold', color=C_P3, va='bottom')
ax2_0.text(0.20, 0.92, r'$N_{\mathrm{bg}}=45$', va='center', fontsize=6.8, color='#555555')

# Dual transfers into R_anc at (0.44, 0.77)
ax2_0.plot([0.16, 0.44], [0.92, 0.77], color=C_P1, lw=1.4)
ax2_0.annotate('', xy=(0.44, 0.77), xytext=(0.16, 0.77),
               arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.3, ls='--'))
ax2_0.annotate('', xy=(0.44, 0.77), xytext=(0.16, 0.62),
               arrowprops=dict(arrowstyle="->", color=C_P3, lw=1.3, ls='--'))
ax2_0.plot(0.44, 0.77, 'o', color='#222222', markersize=4.0, zorder=5)
ax2_0.text(0.42, 0.83, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

fan_y = [0.91, 0.82, 0.72, 0.63]
for yt in fan_y:
    ax2_0.plot([0.44, 0.65], [0.77, yt], color=C_REC, lw=1.1)
    ax2_0.plot(0.65, yt, 'o', color=C_REC, markersize=3.0)
ax2_0.text(0.70, 0.77, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# 3 Minimal Chromosome Tracks
h_bar = 0.075
w_total = 0.92
x0 = 0.04

# 1. Dual Cassette
y1 = 0.36
ax2_0.text(x0, y1 + h_bar + 0.015, 'Dual Cassette', fontsize=6.5, color='#333333', fontweight='bold', va='bottom')
w_seg = w_total * 0.20
ax2_0.add_patch(patches.Rectangle((x0, y1), w_seg, h_bar, facecolor=C_P1))
ax2_0.add_patch(patches.Rectangle((x0 + w_seg, y1), w_seg, h_bar, facecolor=C_P2))
ax2_0.add_patch(patches.Rectangle((x0 + 2*w_seg, y1), w_seg, h_bar, facecolor=C_P1))
ax2_0.add_patch(patches.Rectangle((x0 + 3*w_seg, y1), w_seg, h_bar, facecolor=C_P3))
ax2_0.add_patch(patches.Rectangle((x0 + 4*w_seg, y1), w_seg, h_bar, facecolor=C_P1))
ax2_0.add_patch(patches.Rectangle((x0, y1), w_total, h_bar, fill=False, edgecolor='#222222', lw=0.6))
ax2_0.text(x0 + 0.5*w_seg, y1 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + 1.5*w_seg, y1 + h_bar/2, r'$P_2$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + 2.5*w_seg, y1 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + 3.5*w_seg, y1 + h_bar/2, r'$P_3$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + 4.5*w_seg, y1 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')

# 2. Nested Russian Doll
y2 = 0.19
ax2_0.text(x0, y2 + h_bar + 0.015, 'Nested Doll', fontsize=6.5, color='#333333', fontweight='bold', va='bottom')
w_p1 = w_total * 0.24
w_p2 = w_total * 0.16
w_p3 = w_total * 0.20
ax2_0.add_patch(patches.Rectangle((x0, y2), w_p1, h_bar, facecolor=C_P1))
ax2_0.add_patch(patches.Rectangle((x0 + w_p1, y2), w_p2, h_bar, facecolor=C_P2))
ax2_0.add_patch(patches.Rectangle((x0 + w_p1 + w_p2, y2), w_p3, h_bar, facecolor=C_P3))
ax2_0.add_patch(patches.Rectangle((x0 + w_p1 + w_p2 + w_p3, y2), w_p2, h_bar, facecolor=C_P2))
ax2_0.add_patch(patches.Rectangle((x0 + w_p1 + 2*w_p2 + w_p3, y2), w_p1, h_bar, facecolor=C_P1))
ax2_0.add_patch(patches.Rectangle((x0, y2), w_total, h_bar, fill=False, edgecolor='#222222', lw=0.6))
ax2_0.text(x0 + w_p1/2, y2 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + w_p1 + w_p2/2, y2 + h_bar/2, r'$P_2$', color='white', fontsize=4.8, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + w_p1 + w_p2 + w_p3/2, y2 + h_bar/2, r'$P_3$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + w_p1 + 1.5*w_p2 + w_p3, y2 + h_bar/2, r'$P_2$', color='white', fontsize=4.8, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + w_p1 + 2*w_p2 + 1.5*w_p3, y2 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')

# 3. Triple Crossover
y3 = 0.02
ax2_0.text(x0, y3 + h_bar + 0.015, 'Triple Crossover', fontsize=6.5, color='#333333', fontweight='bold', va='bottom')
w_c1 = w_total * 0.33
w_c2 = w_total * 0.33
w_c3 = w_total * 0.34
ax2_0.add_patch(patches.Rectangle((x0, y3), w_c1, h_bar, facecolor=C_P1))
ax2_0.add_patch(patches.Rectangle((x0 + w_c1, y3), w_c2, h_bar, facecolor=C_P2))
ax2_0.add_patch(patches.Rectangle((x0 + w_c1 + w_c2, y3), w_c3, h_bar, facecolor=C_P3))
ax2_0.add_patch(patches.Rectangle((x0, y3), w_total, h_bar, fill=False, edgecolor='#222222', lw=0.6))
ax2_0.text(x0 + w_c1/2, y3 + h_bar/2, r'$P_1$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + 1.5*w_c1, y3 + h_bar/2, r'$P_2$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')
ax2_0.text(x0 + 2.5*w_c1, y3 + h_bar/2, r'$P_3$', color='white', fontsize=5.2, fontweight='bold', ha='center', va='center')

# --- Row 2, Col 2: Donors Correctly Resolved ---
ax2_1 = fig.add_subplot(gs[1, 1])
topos = ['dual_cassette', 'nested_russian_doll', 'triple_crossover']
topo_names = ['Dual\nCassette', 'Nested\nDoll', 'Triple\nCrossover']
x_t = np.arange(len(topos))
bar_w = 0.25

ax2_1.axhline(2.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)
r_donors_mean  = [df2_raw[(df2_raw['tool'] == 'RhizAeon') & (df2_raw['topology'] == t)]['donors_resolved'].mean() for t in topos]
r_donors_std   = [df2_raw[(df2_raw['tool'] == 'RhizAeon') & (df2_raw['topology'] == t)]['donors_resolved'].std() for t in topos]
sc_donors_mean = [df2_raw[(df2_raw['tool'] == '3SEQ (Consolidated)') & (df2_raw['topology'] == t)]['donors_resolved'].mean() for t in topos]
sc_donors_std  = [df2_raw[(df2_raw['tool'] == '3SEQ (Consolidated)') & (df2_raw['topology'] == t)]['donors_resolved'].std() for t in topos]
s_donors_mean  = [df2_raw[(df2_raw['tool'] == '3SEQ') & (df2_raw['topology'] == t)]['donors_resolved'].mean() for t in topos]
s_donors_std   = [df2_raw[(df2_raw['tool'] == '3SEQ') & (df2_raw['topology'] == t)]['donors_resolved'].std() for t in topos]

ax2_1.bar(x_t - bar_w, r_donors_mean, bar_w, yerr=r_donors_std, capsize=2.5,
          color='#002b66', edgecolor='none', zorder=4,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax2_1.bar(x_t, sc_donors_mean, bar_w, yerr=sc_donors_std, capsize=2.5,
          color='#1f77b4', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax2_1.bar(x_t + bar_w, s_donors_mean, bar_w, yerr=s_donors_std, capsize=2.5,
          color='#85b3d1', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
clean_axis(ax2_1, 'Donors Correctly Resolved', 'Mosaic Topology', 'Donors Resolved', ylim=(0, 2.45))
ax2_1.set_xticks(x_t)
ax2_1.set_xticklabels(topo_names, fontsize=6.8)

# --- Row 2, Col 3: Breakpoint Overlap Jaccard (J_SNP) ---
ax2_2 = fig.add_subplot(gs[1, 2])
ax2_2.axhline(1.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)
r_jacc_mean  = [df2_raw[(df2_raw['tool'] == 'RhizAeon') & (df2_raw['topology'] == t)]['mean_jaccard'].mean() for t in topos]
r_jacc_std   = [df2_raw[(df2_raw['tool'] == 'RhizAeon') & (df2_raw['topology'] == t)]['mean_jaccard'].std() for t in topos]
sc_jacc_mean = [df2_raw[(df2_raw['tool'] == '3SEQ (Consolidated)') & (df2_raw['topology'] == t)]['mean_jaccard'].mean() for t in topos]
sc_jacc_std  = [df2_raw[(df2_raw['tool'] == '3SEQ (Consolidated)') & (df2_raw['topology'] == t)]['mean_jaccard'].std() for t in topos]
s_jacc_mean  = [df2_raw[(df2_raw['tool'] == '3SEQ') & (df2_raw['topology'] == t)]['mean_jaccard'].mean() for t in topos]
s_jacc_std   = [df2_raw[(df2_raw['tool'] == '3SEQ') & (df2_raw['topology'] == t)]['mean_jaccard'].std() for t in topos]

ax2_2.bar(x_t - bar_w, r_jacc_mean, bar_w, yerr=r_jacc_std, capsize=2.5,
          color='#002b66', edgecolor='none', zorder=4,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax2_2.bar(x_t, sc_jacc_mean, bar_w, yerr=sc_jacc_std, capsize=2.5,
          color='#1f77b4', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax2_2.bar(x_t + bar_w, s_jacc_mean, bar_w, yerr=s_jacc_std, capsize=2.5,
          color='#85b3d1', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
clean_axis(ax2_2, r'Breakpoint Overlap ($J_{\mathrm{SNP}}$)', 'Mosaic Topology', 'Jaccard Index', ylim=(0, 1.15))
ax2_2.set_xticks(x_t)
ax2_2.set_xticklabels(topo_names, fontsize=6.8)

# --- Row 2, Col 4: Inferred Events (\hat{K}) vs K with IQR ---
ax2_3 = fig.add_subplot(gs[1, 3])
clean_axis(ax2_3, r'Inferred Events ($\hat{K}$)', r'Descendant Cluster Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-2, 60))
ax2_3.axhline(2.0, color='#222222', linestyle='-', linewidth=1.4, zorder=2)
swarm_sizes_2 = [1, 5, 10]
ax2_3.set_xticks(swarm_sizes_2)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_r  = df2_sum[(df2_sum['tool'] == 'RhizAeon') & (np.isclose(df2_sum['d_parent'], d))].groupby('K')['mean_unique_events'].mean().reset_index()
    sub_3c = df2_sum[(df2_sum['tool'] == '3SEQ (Consolidated)') & (np.isclose(df2_sum['d_parent'], d))].groupby('K')['mean_unique_events'].mean().reset_index()
    sub_3  = df2_sum[(df2_sum['tool'] == '3SEQ') & (np.isclose(df2_sum['d_parent'], d))].groupby('K')['mean_unique_events'].mean().reset_index()
    
    raw_3 = df2_raw[(df2_raw['tool'] == '3SEQ') & (np.isclose(df2_raw['d_parent'], d))]
    q25 = raw_3.groupby('K')['unique_events'].quantile(0.25)
    q75 = raw_3.groupby('K')['unique_events'].quantile(0.75)
    ax2_3.fill_between(sub_3['K'], q25, q75, color=c, alpha=0.12, zorder=2)
    ax2_3.plot(sub_3['K'], sub_3['mean_unique_events'], marker='s', mfc='none', mec=c, mew=1.2,
               ms=4.4, ls=':', lw=1.1, color=c, zorder=3)
    ax2_3.plot(sub_3c['K'], sub_3c['mean_unique_events'], marker='^', mfc='white', mec=c, mew=1.2,
               ms=4.6, ls='--', lw=1.3, color=c, zorder=4)
    ax2_3.plot(sub_r['K'], sub_r['mean_unique_events'], marker='o', mfc=c, mec='white', mew=0.5,
               ms=4.8, ls='-', lw=1.6, color=c, zorder=5)


# ==============================================================================
# ROW 3 (C): SCENARIO 1.3 GHOST DONOR INTROGRESSION
# ==============================================================================
df3_raw = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario3_ghost_donor_introgression/ghost_raw_results.csv')
df3_sum = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario3_ghost_donor_introgression/ghost_summary.csv')

ax3_0 = fig.add_subplot(gs[2, 0])
ax3_0.set_xlim(0, 1)
ax3_0.set_ylim(0, 1)
ax3_0.axis('off')
ax3_0.text(0.0, 1.22, 'C  Ghost Donor Introgression', fontsize=9.5, fontweight='bold', va='top')

# Upper Portion: Reticulation Tree with Unsampled Ghost Donor & Surrogate
ax3_0.plot([0.04, 0.16], [0.92, 0.92], color=C_P1, lw=1.6)
ax3_0.text(0.05, 0.94, r'$P_1$', fontsize=7.2, fontweight='bold', color=C_P1, va='bottom')
ax3_0.plot([0.04, 0.16], [0.77, 0.77], color=C_P2, lw=1.6)
ax3_0.text(0.05, 0.79, r'$P_2$', fontsize=7.2, fontweight='bold', color=C_P2, va='bottom')
ax3_0.plot([0.04, 0.16], [0.60, 0.60], color=C_GHOST, lw=1.6, ls='--')
ax3_0.text(0.05, 0.62, r'$D_{\mathrm{ghost}}$', fontsize=7.0, fontweight='bold', color=C_GHOST, va='bottom')

# Surrogate Sister branch splitting from Ghost lineage
ax3_0.plot([0.16, 0.16], [0.60, 0.48], color=C_GHOST, lw=1.3, ls='--')
ax3_0.plot([0.16, 0.28], [0.48, 0.48], color=C_SURR, lw=1.5)
ax3_0.text(0.30, 0.48, r'$\mathrm{Surr}$', fontsize=7.0, fontweight='bold', color=C_SURR, va='center')
ax3_0.text(0.20, 0.92, r'$N_{\mathrm{bg}}=30\text{--}40$', va='center', fontsize=6.8, color='#555555')

# Ghost Introgression transfer to R_anc
ax3_0.plot([0.16, 0.44], [0.92, 0.77], color=C_P1, lw=1.4)
ax3_0.annotate('', xy=(0.44, 0.77), xytext=(0.16, 0.60),
               arrowprops=dict(arrowstyle="->", color=C_GHOST, lw=1.3, ls='--'))
ax3_0.plot(0.44, 0.77, 'o', color='#222222', markersize=4.0, zorder=5)
ax3_0.text(0.42, 0.83, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

fan_y = [0.91, 0.82, 0.72, 0.63]
for yt in fan_y:
    ax3_0.plot([0.44, 0.65], [0.77, yt], color=C_REC, lw=1.1)
    ax3_0.plot(0.65, yt, 'o', color=C_REC, markersize=3.0)
ax3_0.text(0.70, 0.77, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# Minimal Chromosome Track
bar_y = 0.08
bar_h = 0.13
w_bar = 0.90
x_bar = 0.05
ax3_0.add_patch(patches.Rectangle((x_bar, bar_y), 0.30 * w_bar, bar_h, facecolor=C_P1, edgecolor='none'))
ax3_0.text(x_bar + 0.15 * w_bar, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')
ax3_0.add_patch(patches.Rectangle((x_bar + 0.30 * w_bar, bar_y), 0.30 * w_bar, bar_h, facecolor=C_GHOST, edgecolor='none'))
ax3_0.text(x_bar + 0.45 * w_bar, bar_y + bar_h/2, r'$D_{\mathrm{ghost}}$ (1.5 kb)', color='white', fontsize=6.2, fontweight='bold', ha='center', va='center')
ax3_0.add_patch(patches.Rectangle((x_bar + 0.60 * w_bar, bar_y), 0.40 * w_bar, bar_h, facecolor=C_P1, edgecolor='none'))
ax3_0.text(x_bar + 0.80 * w_bar, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')
ax3_0.add_patch(patches.Rectangle((x_bar, bar_y), w_bar, bar_h, fill=False, edgecolor='#222222', lw=0.7))

ax3_0.plot([x_bar, x_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax3_0.plot([x_bar + 0.30 * w_bar, x_bar + 0.30 * w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax3_0.plot([x_bar + 0.60 * w_bar, x_bar + 0.60 * w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax3_0.plot([x_bar + w_bar, x_bar + w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax3_0.text(x_bar, bar_y - 0.038, '0', fontsize=5.5, ha='left', va='top', color='#555555')
ax3_0.text(x_bar + 0.30 * w_bar, bar_y - 0.038, '1.5k', fontsize=5.5, ha='center', va='top', color='#555555')
ax3_0.text(x_bar + 0.60 * w_bar, bar_y - 0.038, '3.0k', fontsize=5.5, ha='center', va='top', color='#555555')
ax3_0.text(x_bar + w_bar, bar_y - 0.038, '5.0 kb', fontsize=5.5, ha='right', va='top', color='#555555')

# --- Row 3, Col 2: Background Parent Misclassification (%) ---
ax3_1 = fig.add_subplot(gs[2, 1])
clean_axis(ax3_1, 'Background Parent FPR (%)', 'Ghost Donor Regime', 'False Alarms (%)', ylim=(-3, 115))
ax3_1.axhline(0.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

regimes_3 = ['surrogate_sister', 'pure_ghost', 'deep_reservoir']
reg_labels_3 = ['Surrogate\nSister', 'Pure\nGhost', 'Deep\nReservoir']
x_r3 = np.arange(len(regimes_3))
bar_w = 0.25

r_pfpr_mean  = [df3_raw[(df3_raw['tool'] == 'RhizAeon') & (df3_raw['regime'] == r)]['parent_fpr'].mean() * 100 for r in regimes_3]
r_pfpr_std   = [df3_raw[(df3_raw['tool'] == 'RhizAeon') & (df3_raw['regime'] == r)]['parent_fpr'].std() * 100 for r in regimes_3]
sc_pfpr_mean = [df3_raw[(df3_raw['tool'] == '3SEQ (Consolidated)') & (df3_raw['regime'] == r)]['parent_fpr'].mean() * 100 for r in regimes_3]
sc_pfpr_std  = [df3_raw[(df3_raw['tool'] == '3SEQ (Consolidated)') & (df3_raw['regime'] == r)]['parent_fpr'].std() * 100 for r in regimes_3]
s_pfpr_mean  = [df3_raw[(df3_raw['tool'] == '3SEQ') & (df3_raw['regime'] == r)]['parent_fpr'].mean() * 100 for r in regimes_3]
s_pfpr_std   = [df3_raw[(df3_raw['tool'] == '3SEQ') & (df3_raw['regime'] == r)]['parent_fpr'].std() * 100 for r in regimes_3]

ax3_1.bar(x_r3 - bar_w, r_pfpr_mean, bar_w, yerr=r_pfpr_std, capsize=2.5,
          color='#002b66', edgecolor='none', zorder=4,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax3_1.bar(x_r3, sc_pfpr_mean, bar_w, yerr=sc_pfpr_std, capsize=2.5,
          color='#1f77b4', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax3_1.bar(x_r3 + bar_w, s_pfpr_mean, bar_w, yerr=s_pfpr_std, capsize=2.5,
          color='#85b3d1', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax3_1.set_xticks(x_r3)
ax3_1.set_xticklabels(reg_labels_3, fontsize=6.8)

# --- Row 3, Col 3: Donor Attribution Accuracy (%) ---
ax3_2 = fig.add_subplot(gs[2, 2])
clean_axis(ax3_2, 'Surrogate Donor Acc. (%)', 'Parent Divergence (d)', 'Attribution Acc. (%)', ylim=(-3, 115))
ax3_2.axhline(100.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)
divs = [0.02, 0.05, 0.10]
div_labels = ['0.02', '0.05', '0.10']
x_d3 = np.arange(len(divs))

r_dacc_mean  = [df3_raw[(df3_raw['tool'] == 'RhizAeon') & (df3_raw['regime'] == 'surrogate_sister') & (np.isclose(df3_raw['d_parent'], d))]['donor_acc'].mean() * 100 for d in divs]
r_dacc_std   = [df3_raw[(df3_raw['tool'] == 'RhizAeon') & (df3_raw['regime'] == 'surrogate_sister') & (np.isclose(df3_raw['d_parent'], d))]['donor_acc'].std() * 100 for d in divs]
sc_dacc_mean = [df3_raw[(df3_raw['tool'] == '3SEQ (Consolidated)') & (df3_raw['regime'] == 'surrogate_sister') & (np.isclose(df3_raw['d_parent'], d))]['donor_acc'].mean() * 100 for d in divs]
sc_dacc_std  = [df3_raw[(df3_raw['tool'] == '3SEQ (Consolidated)') & (df3_raw['regime'] == 'surrogate_sister') & (np.isclose(df3_raw['d_parent'], d))]['donor_acc'].std() * 100 for d in divs]
s_dacc_mean  = [df3_raw[(df3_raw['tool'] == '3SEQ') & (df3_raw['regime'] == 'surrogate_sister') & (np.isclose(df3_raw['d_parent'], d))]['donor_acc'].mean() * 100 for d in divs]
s_dacc_std   = [df3_raw[(df3_raw['tool'] == '3SEQ') & (df3_raw['regime'] == 'surrogate_sister') & (np.isclose(df3_raw['d_parent'], d))]['donor_acc'].std() * 100 for d in divs]

ax3_2.bar(x_d3 - bar_w, r_dacc_mean, bar_w, yerr=r_dacc_std, capsize=2.5,
          color='#002b66', edgecolor='none', zorder=4,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax3_2.bar(x_d3, sc_dacc_mean, bar_w, yerr=sc_dacc_std, capsize=2.5,
          color='#1f77b4', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax3_2.bar(x_d3 + bar_w, s_dacc_mean, bar_w, yerr=s_dacc_std, capsize=2.5,
          color='#85b3d1', edgecolor='#002b66', lw=0.8, zorder=3,
          error_kw=dict(lw=0.9, capthick=0.9, ecolor='#002b66'))
ax3_2.set_xticks(x_d3)
ax3_2.set_xticklabels(div_labels, fontsize=6.8)

# --- Row 3, Col 4: Inferred Events (\hat{K}) vs K with IQR ---
ax3_3 = fig.add_subplot(gs[2, 3])
clean_axis(ax3_3, r'Inferred Events ($\hat{K}$)', r'Descendant Cluster Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-1.5, 45))
ax3_3.axhline(1.0, color='#222222', linestyle='-', linewidth=1.4, zorder=2)
ks_3 = [2, 5, 10]
ax3_3.set_xticks(ks_3)
sub_surr_raw = df3_raw[df3_raw['regime'] == 'surrogate_sister']

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    raw_3 = sub_surr_raw[(sub_surr_raw['tool'] == '3SEQ') & (np.isclose(sub_surr_raw['d_parent'], d))]
    t_medians = [np.median(raw_3[raw_3['K'] == k]['unique_events'].values) for k in ks_3]
    t_q25 = [np.percentile(raw_3[raw_3['K'] == k]['unique_events'].values, 25) for k in ks_3]
    t_q75 = [np.percentile(raw_3[raw_3['K'] == k]['unique_events'].values, 75) for k in ks_3]
    ax3_3.plot(ks_3, t_medians, 's:', color=c, lw=1.1, markersize=4.4, markerfacecolor='none', markeredgewidth=1.2, zorder=3)
    ax3_3.fill_between(ks_3, t_q25, t_q75, color=c, alpha=0.12, zorder=2)

    raw_3c = sub_surr_raw[(sub_surr_raw['tool'] == '3SEQ (Consolidated)') & (np.isclose(sub_surr_raw['d_parent'], d))]
    tc_medians = [np.median(raw_3c[raw_3c['K'] == k]['unique_events'].values) for k in ks_3]
    ax3_3.plot(ks_3, tc_medians, '^--', color=c, lw=1.3, markersize=4.6, markerfacecolor='white', markeredgewidth=1.2, zorder=4)

    raw_r = sub_surr_raw[(sub_surr_raw['tool'] == 'RhizAeon') & (np.isclose(sub_surr_raw['d_parent'], d))]
    r_medians = [np.median(raw_r[raw_r['K'] == k]['unique_events'].values) for k in ks_3]
    r_q25 = [np.percentile(raw_r[raw_r['K'] == k]['unique_events'].values, 25) for k in ks_3]
    r_q75 = [np.percentile(raw_r[raw_r['K'] == k]['unique_events'].values, 75) for k in ks_3]
    ax3_3.plot(ks_3, r_medians, 'o-', color=c, lw=1.5, markersize=4.4, zorder=5)
    ax3_3.fill_between(ks_3, r_q25, r_q75, color=c, alpha=0.18, zorder=3)


# ==============================================================================
# ROW 4 (D): SCENARIO 1.4 BACTERIAL MICRO-CONVERSIONS
# ==============================================================================
df4_raw = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario4_micro_conversion_tracts/micro_conversion_raw_results.csv')
df4_sum = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario4_micro_conversion_tracts/micro_conversion_summary.csv')

ax4_0 = fig.add_subplot(gs[3, 0])
ax4_0.set_xlim(0, 1)
ax4_0.set_ylim(0, 1)
ax4_0.axis('off')
ax4_0.text(0.0, 1.22, 'D  Bacterial Micro-Conversions', fontsize=9.5, fontweight='bold', va='top')

# Tree
ax4_0.plot([0.04, 0.16], [0.88, 0.88], color=C_P1, lw=1.6)
ax4_0.plot([0.16, 0.16], [0.76, 0.98], color=C_P1, lw=1.6)
ax4_0.plot([0.16, 0.36], [0.98, 0.98], color=C_P1, lw=1.3)
ax4_0.plot([0.16, 0.36], [0.76, 0.76], color=C_P1, lw=1.5)
ax4_0.text(0.39, 0.98, r'$N_{\mathrm{bg}}=30$', va='center', fontsize=6.8, color='#555555')
ax4_0.text(0.20, 0.81, r'$P_1$', va='bottom', fontsize=7.8, fontweight='bold', color=C_P1)

ax4_0.plot([0.04, 0.16], [0.48, 0.48], color=C_P2, lw=1.6)
ax4_0.plot([0.16, 0.16], [0.36, 0.58], color=C_P2, lw=1.6)
ax4_0.plot([0.16, 0.36], [0.58, 0.58], color=C_P2, lw=1.5)
ax4_0.plot([0.16, 0.36], [0.36, 0.36], color=C_P2, lw=1.3)
ax4_0.text(0.20, 0.53, r'$P_2$', va='top', fontsize=7.8, fontweight='bold', color=C_P2)

ax4_0.plot([0.36, 0.58], [0.76, 0.64], color=C_P1, lw=1.4)
ax4_0.annotate('', xy=(0.58, 0.64), xytext=(0.36, 0.58),
               arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.4, ls='--'))
ax4_0.plot(0.58, 0.64, 'o', color='#222222', markersize=4.0, zorder=5)
ax4_0.text(0.57, 0.71, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

for yt in [0.82, 0.70, 0.58, 0.46]:
    ax4_0.plot([0.58, 0.76], [0.64, yt], color=C_REC, lw=1.2)
    ax4_0.plot(0.76, yt, 'o', color=C_REC, markersize=3.0)
ax4_0.text(0.80, 0.64, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# Minimal Chromosome Track with Micro-Conversion Tract
bar_y = 0.08
bar_h = 0.13
w_bar = 0.90
x_bar = 0.05
w1 = 0.42 * w_bar
w2 = 0.10 * w_bar
w3 = 0.48 * w_bar

ax4_0.add_patch(patches.Rectangle((x_bar, bar_y), w1, bar_h, facecolor=C_P1, edgecolor='none'))
ax4_0.text(x_bar + w1/2, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')
ax4_0.add_patch(patches.Rectangle((x_bar + w1, bar_y), w2, bar_h, facecolor=C_P2, edgecolor='none'))
ax4_0.text(x_bar + w1 + w2/2, bar_y + bar_h + 0.02, r'$L_{\mathrm{tract}}$', color=C_P2, fontsize=6.0, fontweight='bold', ha='center', va='bottom')
ax4_0.add_patch(patches.Rectangle((x_bar + w1 + w2, bar_y), w3, bar_h, facecolor=C_P1, edgecolor='none'))
ax4_0.text(x_bar + w1 + w2 + w3/2, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')
ax4_0.add_patch(patches.Rectangle((x_bar, bar_y), w_bar, bar_h, fill=False, edgecolor='#222222', lw=0.7))

ax4_0.plot([x_bar, x_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax4_0.plot([x_bar + w1, x_bar + w1], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax4_0.plot([x_bar + w_bar, x_bar + w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax4_0.text(x_bar, bar_y - 0.038, '0', fontsize=5.5, ha='left', va='top', color='#555555')
ax4_0.text(x_bar + w1, bar_y - 0.038, '2.1k', fontsize=5.5, ha='center', va='top', color='#555555')
ax4_0.text(x_bar + w_bar, bar_y - 0.038, '5.0 kb', fontsize=5.5, ha='right', va='top', color='#555555')

# --- Row 4, Col 2: Detection Power vs Tract Length ---
ax4_1 = fig.add_subplot(gs[3, 1])
clean_axis(ax4_1, 'Detection Power vs Tract Length', 'Tract Length (nt)', 'Sensitivity (%)', ylim=(-3, 105))
ax4_1.axhline(100.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)
tract_lens = [50, 100, 200, 400, 800]
x_t4 = np.arange(len(tract_lens))
ax4_1.set_xticks(x_t4)
ax4_1.set_xticklabels([str(t) for t in tract_lens], fontsize=6.8)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub = df4_sum[np.isclose(df4_sum['d_parent'], d)]
    s_power  = [sub[(sub['tool'] == '3SEQ') & (sub['tract_len'] == t)]['recall'].mean() * 100 for t in tract_lens]
    sc_power = [sub[(sub['tool'] == '3SEQ (Consolidated)') & (sub['tract_len'] == t)]['recall'].mean() * 100 for t in tract_lens]
    r_power  = [sub[(sub['tool'] == 'RhizAeon') & (sub['tract_len'] == t)]['recall'].mean() * 100 for t in tract_lens]
    ax4_1.plot(x_t4, s_power, 's:', color=c, lw=1.1, markersize=4.4, markerfacecolor='none', markeredgewidth=1.2, zorder=3)
    ax4_1.plot(x_t4, sc_power, '^--', color=c, lw=1.3, markersize=4.6, markerfacecolor='white', markeredgewidth=1.2, zorder=4)
    ax4_1.plot(x_t4, r_power, 'o-', color=c, lw=1.5, markersize=4.5, zorder=5)

# --- Row 4, Col 3: Background Taxon FPR (%) ---
ax4_2 = fig.add_subplot(gs[3, 2])
clean_axis(ax4_2, 'Background Taxon FPR (%)', 'Tract Length (nt)', 'False Alarms (%)', ylim=(-3, 105))
ax4_2.axhline(0.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)
ax4_2.set_xticks(x_t4)
ax4_2.set_xticklabels([str(t) for t in tract_lens], fontsize=6.8)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub = df4_sum[np.isclose(df4_sum['d_parent'], d)]
    s_fpr  = [sub[(sub['tool'] == '3SEQ') & (sub['tract_len'] == t)]['bg_fpr'].mean() * 100 for t in tract_lens]
    sc_fpr = [sub[(sub['tool'] == '3SEQ (Consolidated)') & (sub['tract_len'] == t)]['bg_fpr'].mean() * 100 for t in tract_lens]
    r_fpr  = [sub[(sub['tool'] == 'RhizAeon') & (sub['tract_len'] == t)]['bg_fpr'].mean() * 100 for t in tract_lens]
    ax4_2.plot(x_t4, s_fpr, 's:', color=c, lw=1.1, markersize=4.4, markerfacecolor='none', markeredgewidth=1.2, zorder=3)
    ax4_2.plot(x_t4, sc_fpr, '^--', color=c, lw=1.3, markersize=4.6, markerfacecolor='white', markeredgewidth=1.2, zorder=4)
    ax4_2.plot(x_t4, r_fpr, 'o-', color=c, lw=1.5, markersize=4.5, zorder=5)

# --- Row 4, Col 4: Inferred Events (\hat{K}) vs K with IQR ---
ax4_3 = fig.add_subplot(gs[3, 3])
clean_axis(ax4_3, r'Inferred Events ($\hat{K}$)', r'Descendant Cluster Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-1.5, 45))
ax4_3.axhline(1.0, color='#222222', linestyle='-', linewidth=1.4, zorder=2)
ks_4 = [2, 5, 10]
ax4_3.set_xticks(ks_4)
sub_detectable = df4_raw[df4_raw['tract_len'] >= 200]

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    raw_3 = sub_detectable[(sub_detectable['tool'] == '3SEQ') & (np.isclose(sub_detectable['d_parent'], d))]
    t_medians = [np.median(raw_3[raw_3['K'] == k]['unique_events'].values) for k in ks_4]
    t_q25 = [np.percentile(raw_3[raw_3['K'] == k]['unique_events'].values, 25) for k in ks_4]
    t_q75 = [np.percentile(raw_3[raw_3['K'] == k]['unique_events'].values, 75) for k in ks_4]
    ax4_3.plot(ks_4, t_medians, 's:', color=c, lw=1.1, markersize=4.4, markerfacecolor='none', markeredgewidth=1.2, zorder=3)
    ax4_3.fill_between(ks_4, t_q25, t_q75, color=c, alpha=0.12, zorder=2)

    raw_3c = sub_detectable[(sub_detectable['tool'] == '3SEQ (Consolidated)') & (np.isclose(sub_detectable['d_parent'], d))]
    tc_medians = [np.median(raw_3c[raw_3c['K'] == k]['unique_events'].values) for k in ks_4]
    ax4_3.plot(ks_4, tc_medians, '^--', color=c, lw=1.3, markersize=4.6, markerfacecolor='white', markeredgewidth=1.2, zorder=4)

    raw_r = sub_detectable[(sub_detectable['tool'] == 'RhizAeon') & (np.isclose(sub_detectable['d_parent'], d))]
    r_medians = [np.median(raw_r[raw_r['K'] == k]['unique_events'].values) for k in ks_4]
    r_q25 = [np.percentile(raw_r[raw_r['K'] == k]['unique_events'].values, 25) for k in ks_4]
    r_q75 = [np.percentile(raw_r[raw_r['K'] == k]['unique_events'].values, 75) for k in ks_4]
    ax4_3.plot(ks_4, r_medians, 'o-', color=c, lw=1.5, markersize=4.4, zorder=5)
    ax4_3.fill_between(ks_4, r_q25, r_q75, color=c, alpha=0.18, zorder=3)


# ==============================================================================
# Unified Top-Hoisted Master Legend
# ==============================================================================
legend_handles = [
    mlines.Line2D([], [], color='#002b66', marker='o', lw=1.5, markersize=4.5, label='RhizAeon'),
    mlines.Line2D([], [], color='#002b66', marker='^', lw=1.3, linestyle='--', markersize=4.6,
                  markerfacecolor='white', markeredgewidth=1.2, label='3SEQ (Consolidated)'),
    mlines.Line2D([], [], color='#002b66', marker='s', lw=1.1, linestyle=':', markersize=4.4,
                  markerfacecolor='none', markeredgewidth=1.2, label='3SEQ (Raw)'),
    mlines.Line2D([], [], color='#222222', linestyle='--', lw=1.2, label='Ground Truth'),
    patches.Patch(facecolor=colors_div[0.02], label=r'$d = 0.02$'),
    patches.Patch(facecolor=colors_div[0.05], label=r'$d = 0.05$'),
    patches.Patch(facecolor=colors_div[0.10], label=r'$d = 0.10$'),
]

fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.56, 0.985),
           ncol=7, frameon=False, fontsize=7.2, handletextpad=0.35, columnspacing=1.2)

# Save master outputs to both directories
out_master_png = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_specialized_scenarios.png'
out_master_pdf = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_specialized_scenarios.pdf'
plt.savefig(out_master_png, dpi=300)
plt.savefig(out_master_pdf)

paper_master_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_specialized_scenarios.png'
paper_master_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_specialized_scenarios.pdf'
plt.savefig(paper_master_png, dpi=300)
plt.savefig(paper_master_pdf)
plt.close()
print("Master Figure 2 successfully generated and saved.")
