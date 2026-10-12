#!/usr/bin/env python3
"""
Publication-Quality Row 3: Scenario 1.3 Ghost Donor Introgression
Minimal iconic schematic matching Rows 1 & 2 rhythm + High-Density Empirical Metrics.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

# Load authentic benchmark data
raw_df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario3_ghost_donor_introgression/ghost_raw_results.csv')
sum_df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario3_ghost_donor_introgression/ghost_summary.csv')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

# 11.8 x 3.1 in, 300 DPI: Widescreen row matching Rows 1 & 2 exactly
fig = plt.figure(figsize=(11.8, 3.1), dpi=300)
gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 0.95, 0.95, 0.95], wspace=0.34,
                      left=0.035, right=0.985, top=0.78, bottom=0.18)

C_P1 = '#1f77b4'       # Parent 1 Home (Blue)
C_P2 = '#d95f02'       # Parent 2 Background (Orange)
C_SURR = '#2ca02c'     # Surrogate Sister (Green)
C_GHOST = '#7570b3'    # Unsampled Ghost Donor (Purple / Dashed)
C_REC = '#0088cc'      # Recombinant swarm (Cyan/Blue)
colors_div = {
    0.02: '#4dbbd5',   # light cyan-blue
    0.05: '#1f77b4',   # mid blue
    0.10: '#002b66',   # dark navy
}

# ==============================================================================
# Helper for Clean Axes (Matching Rows 1 & 2 Exactly)
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
# 1. Column 1: Minimal Iconic Schematic (Matching Rows 1 & 2 Rhythm, Minimal Text)
# ==============================================================================
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.axis('off')

# Title (Zero narrative subtitle, exactly like Rows 1 & 2)
ax0.text(0.0, 1.22, 'C  Ghost Donor Introgression', fontsize=9.5, fontweight='bold', va='top')

# --- Upper Portion: Reticulation Tree with Unsampled Ghost Donor & Surrogate ---
# P1 Home Clade (Top)
ax0.plot([0.04, 0.16], [0.92, 0.92], color=C_P1, lw=1.6)
ax0.text(0.05, 0.94, r'$P_1$', fontsize=7.2, fontweight='bold', color=C_P1, va='bottom')

# P2 Innocent Background Clade (Middle)
ax0.plot([0.04, 0.16], [0.77, 0.77], color=C_P2, lw=1.6)
ax0.text(0.05, 0.79, r'$P_2$', fontsize=7.2, fontweight='bold', color=C_P2, va='bottom')

# Ghost Donor Lineage (Bottom - Dashed)
ax0.plot([0.04, 0.16], [0.60, 0.60], color=C_GHOST, lw=1.6, ls='--')
ax0.text(0.05, 0.62, r'$D_{\mathrm{ghost}}$', fontsize=7.0, fontweight='bold', color=C_GHOST, va='bottom')

# Surrogate Sister branch splitting from Ghost lineage (branching downward to avoid transfer arrow)
ax0.plot([0.16, 0.16], [0.60, 0.48], color=C_GHOST, lw=1.3, ls='--')
ax0.plot([0.16, 0.28], [0.48, 0.48], color=C_SURR, lw=1.5)
ax0.text(0.30, 0.48, r'$\mathrm{Surr}$', fontsize=7.0, fontweight='bold', color=C_SURR, va='center')

# Background taxon count
ax0.text(0.20, 0.92, r'$N_{\mathrm{bg}}=30\text{--}40$', va='center', fontsize=6.8, color='#555555')

# Ghost Introgression transfer to R_anc
ax0.plot([0.16, 0.44], [0.92, 0.77], color=C_P1, lw=1.4)
ax0.annotate('', xy=(0.44, 0.77), xytext=(0.16, 0.60),
             arrowprops=dict(arrowstyle="->", color=C_GHOST, lw=1.3, ls='--'))

# Reticulation Node Dot
ax0.plot(0.44, 0.77, 'o', color='#222222', markersize=4.0, zorder=5)
ax0.text(0.42, 0.83, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

# Clonal expansion fan to contemporary isolates
fan_y = [0.91, 0.82, 0.72, 0.63]
for yt in fan_y:
    ax0.plot([0.44, 0.65], [0.77, yt], color=C_REC, lw=1.1)
    ax0.plot(0.65, yt, 'o', color=C_REC, markersize=3.0)

# Minimal swarm label
ax0.text(0.70, 0.77, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# --- Lower Portion: Minimal Chromosome Track (Matching Row 1 Exactly) ---
bar_y = 0.08
bar_h = 0.13
w_bar = 0.90
x_bar = 0.05

# 0 to 1.5 kb (P1 Home) - 30%
rect1 = patches.Rectangle((x_bar, bar_y), 0.30 * w_bar, bar_h, facecolor=C_P1, edgecolor='none')
ax0.add_patch(rect1)
ax0.text(x_bar + 0.15 * w_bar, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

# 1.5 to 3.0 kb (Ghost donor tract) - 30%
rect2 = patches.Rectangle((x_bar + 0.30 * w_bar, bar_y), 0.30 * w_bar, bar_h, facecolor=C_GHOST, edgecolor='none')
ax0.add_patch(rect2)
ax0.text(x_bar + 0.45 * w_bar, bar_y + bar_h/2, r'$D_{\mathrm{ghost}}$ (1.5 kb)', color='white', fontsize=6.2, fontweight='bold', ha='center', va='center')

# 3.0 to 5.0 kb (P1 Home) - 40%
rect3 = patches.Rectangle((x_bar + 0.60 * w_bar, bar_y), 0.40 * w_bar, bar_h, facecolor=C_P1, edgecolor='none')
ax0.add_patch(rect3)
ax0.text(x_bar + 0.80 * w_bar, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

# Outer border
track_box = patches.Rectangle((x_bar, bar_y), w_bar, bar_h, fill=False, edgecolor='#222222', lw=0.7)
ax0.add_patch(track_box)

# Coordinate ticks below bar (identical to Row 1)
ax0.plot([x_bar, x_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([x_bar + 0.30 * w_bar, x_bar + 0.30 * w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([x_bar + 0.60 * w_bar, x_bar + 0.60 * w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([x_bar + w_bar, x_bar + w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.text(x_bar, bar_y - 0.038, '0', fontsize=5.5, ha='left', va='top', color='#555555')
ax0.text(x_bar + 0.30 * w_bar, bar_y - 0.038, '1.5k', fontsize=5.5, ha='center', va='top', color='#555555')
ax0.text(x_bar + 0.60 * w_bar, bar_y - 0.038, '3.0k', fontsize=5.5, ha='center', va='top', color='#555555')
ax0.text(x_bar + w_bar, bar_y - 0.038, '5.0 kb', fontsize=5.5, ha='right', va='top', color='#555555')


# ==============================================================================
# 2. Column 2: False Framing of Innocent Taxa (%) across Ghost Regimes (NO HATCHING)
# ==============================================================================
ax1 = fig.add_subplot(gs[0, 1])
clean_axis(ax1, 'False Parent Framing (%)', 'Ghost Donor Regime', 'False Alarms (%)', ylim=(-3, 115))
ax1.axhline(0.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

regimes = ['surrogate_sister', 'pure_ghost', 'deep_reservoir']
reg_labels = ['Surrogate\nSister', 'Pure\nGhost', 'Deep\nReservoir']
x_r = np.arange(len(regimes))
bar_w = 0.35

r_pfpr_mean = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['regime'] == r)]['parent_fpr'].mean() * 100 for r in regimes]
r_pfpr_std  = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['regime'] == r)]['parent_fpr'].std() * 100 for r in regimes]
s_pfpr_mean = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['regime'] == r)]['parent_fpr'].mean() * 100 for r in regimes]
s_pfpr_std  = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['regime'] == r)]['parent_fpr'].std() * 100 for r in regimes]

# Pure Solid Fills (ZERO HATCHING)
ax1.bar(x_r - bar_w/2, r_pfpr_mean, bar_w, yerr=r_pfpr_std, capsize=3,
        label='RhizAeon', color='#002b66', edgecolor='none', zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))
ax1.bar(x_r + bar_w/2, s_pfpr_mean, bar_w, yerr=s_pfpr_std, capsize=3,
        label='3SEQ', color='#85b3d1', edgecolor='#002b66', lw=0.9, zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))

ax1.set_xticks(x_r)
ax1.set_xticklabels(reg_labels, fontsize=6.8)


# ==============================================================================
# 3. Column 3: Donor Attribution Accuracy (%) in Surrogate Sister Regime (NO HATCHING)
# ==============================================================================
ax2 = fig.add_subplot(gs[0, 2])
clean_axis(ax2, 'Donor Attribution Accuracy (%)', 'Parent Divergence (d)', 'Attribution Acc. (%)', ylim=(-3, 115))
ax2.axhline(100.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

divs = [0.02, 0.05, 0.10]
div_labels = ['0.02', '0.05', '0.10']
x_d = np.arange(len(divs))

r_dacc_mean = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['regime'] == 'surrogate_sister') & (np.isclose(raw_df['d_parent'], d))]['donor_acc'].mean() * 100 for d in divs]
r_dacc_std  = [raw_df[(raw_df['tool'] == 'RhizAeon') & (raw_df['regime'] == 'surrogate_sister') & (np.isclose(raw_df['d_parent'], d))]['donor_acc'].std() * 100 for d in divs]
s_dacc_mean = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['regime'] == 'surrogate_sister') & (np.isclose(raw_df['d_parent'], d))]['donor_acc'].mean() * 100 for d in divs]
s_dacc_std  = [raw_df[(raw_df['tool'] == '3SEQ') & (raw_df['regime'] == 'surrogate_sister') & (np.isclose(raw_df['d_parent'], d))]['donor_acc'].std() * 100 for d in divs]

# Pure Solid Fills (ZERO HATCHING)
ax2.bar(x_d - bar_w/2, r_dacc_mean, bar_w, yerr=r_dacc_std, capsize=3,
        color='#002b66', edgecolor='none', zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))
ax2.bar(x_d + bar_w/2, s_dacc_mean, bar_w, yerr=s_dacc_std, capsize=3,
        color='#85b3d1', edgecolor='#002b66', lw=0.9, zorder=3,
        error_kw=dict(lw=1.0, capthick=1.0, ecolor='#002b66'))

ax2.set_xticks(x_d)
ax2.set_xticklabels(div_labels, fontsize=6.8)


# ==============================================================================
# 4. Column 4: Inferred Events (\hat{K}) vs Swarm Size (K) in Surrogate Sister
# ==============================================================================
ax3 = fig.add_subplot(gs[0, 3])
clean_axis(ax3, r'Inferred Events ($\hat{K}$)', r'Descendant Swarm Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-1.5, 45))
ax3.axhline(1.0, color='#222222', linestyle='-', linewidth=1.4, zorder=2) # True Ground Truth = 1.0

ks = [2, 5, 10]
ax3.set_xticks(ks)

sub_surr_raw = raw_df[raw_df['regime'] == 'surrogate_sister']
sub_surr_sum = sum_df[sum_df['regime'] == 'surrogate_sister']

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    
    # 3SEQ Event Inflation IQR
    raw_3 = sub_surr_raw[(sub_surr_raw['tool'] == '3SEQ') & (np.isclose(sub_surr_raw['d_parent'], d))]
    t_medians = []
    t_q25 = []
    t_q75 = []
    for k in ks:
        vals = raw_3[raw_3['K'] == k]['unique_events'].values
        t_medians.append(np.median(vals))
        t_q25.append(np.percentile(vals, 25))
        t_q75.append(np.percentile(vals, 75))
        
    ax3.plot(ks, t_medians, 's:', color=c, lw=1.3, markersize=4.8, markerfacecolor='none', markeredgewidth=1.3, zorder=4)
    ax3.fill_between(ks, t_q25, t_q75, color=c, alpha=0.15, zorder=2)

    # RhizAeon Clonal Resolution
    raw_r = sub_surr_raw[(sub_surr_raw['tool'] == 'RhizAeon') & (np.isclose(sub_surr_raw['d_parent'], d))]
    r_medians = []
    r_q25 = []
    r_q75 = []
    for k in ks:
        vals = raw_r[raw_r['K'] == k]['unique_events'].values
        r_medians.append(np.median(vals))
        r_q25.append(np.percentile(vals, 25))
        r_q75.append(np.percentile(vals, 75))
        
    ax3.plot(ks, r_medians, 'o-', color=c, lw=1.4, markersize=4.2, zorder=5)
    ax3.fill_between(ks, r_q25, r_q75, color=c, alpha=0.20, zorder=3)


# ==============================================================================
# Unified Top-Hoisted Legend (Identical to Rows 1 & 2)
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
out_png = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_panel_row3_scenario13.png'
out_pdf = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_panel_row3_scenario13.pdf'
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)

paper_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_panel_row3_scenario13.png'
paper_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_panel_row3_scenario13.pdf'
plt.savefig(paper_png, dpi=300)
plt.savefig(paper_pdf)
plt.close()
print("Row 3 successfully generated and saved.")
