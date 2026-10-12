#!/usr/bin/env python3
"""
Publication-Quality Row 4: Scenario 1.4 Bacterial Micro-Conversions
Minimal iconic schematic matching Rows 1-3 rhythm + High-Density Empirical Metrics.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

# Load authentic benchmark data
raw_df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario4_micro_conversion_tracts/micro_conversion_raw_results.csv')
sum_df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario4_micro_conversion_tracts/micro_conversion_summary.csv')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

# 11.8 x 3.1 in, 300 DPI: Widescreen row matching Rows 1-3 exactly
fig = plt.figure(figsize=(11.8, 3.1), dpi=300)
gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 0.95, 0.95, 0.95], wspace=0.34,
                      left=0.035, right=0.985, top=0.78, bottom=0.18)

C_P1 = '#1f77b4'       # Parent 1 Home (Blue)
C_P2 = '#d95f02'       # Parent 2 Donor (Orange)
C_REC = '#0088cc'      # Recombinant swarm (Cyan/Blue)
colors_div = {
    0.02: '#4dbbd5',   # light cyan-blue
    0.05: '#1f77b4',   # mid blue
    0.10: '#002b66',   # dark navy
}

# ==============================================================================
# Helper for Clean Axes (Matching Rows 1-3 Exactly)
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
# 1. Column 1: Minimal Iconic Schematic (Matching Rows 1-3 Rhythm, Minimal Text)
# ==============================================================================
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.axis('off')

# Title (Zero narrative subtitle, exactly like Rows 1-3)
ax0.text(0.0, 1.22, 'D  Bacterial Micro-Conversions', fontsize=9.5, fontweight='bold', va='top')

# Upper Portion: Reticulation Tree with Micro-Tract Transfer
# P1 Home Clade (Top)
ax0.plot([0.04, 0.16], [0.88, 0.88], color=C_P1, lw=1.6)
ax0.plot([0.16, 0.16], [0.76, 0.98], color=C_P1, lw=1.6)
ax0.plot([0.16, 0.36], [0.98, 0.98], color=C_P1, lw=1.3)
ax0.plot([0.16, 0.36], [0.76, 0.76], color=C_P1, lw=1.5)
ax0.text(0.39, 0.98, r'$N_{\mathrm{bg}}=30$', va='center', fontsize=6.8, color='#555555')
ax0.text(0.20, 0.81, r'$P_1$', va='bottom', fontsize=7.8, fontweight='bold', color=C_P1)

# P2 Donor Clade (Bottom)
ax0.plot([0.04, 0.16], [0.48, 0.48], color=C_P2, lw=1.6)
ax0.plot([0.16, 0.16], [0.36, 0.58], color=C_P2, lw=1.6)
ax0.plot([0.16, 0.36], [0.58, 0.58], color=C_P2, lw=1.5)
ax0.plot([0.16, 0.36], [0.36, 0.36], color=C_P2, lw=1.3)
ax0.text(0.20, 0.53, r'$P_2$', va='top', fontsize=7.8, fontweight='bold', color=C_P2)

# Reticulation transfer to R_anc
ax0.plot([0.36, 0.58], [0.76, 0.64], color=C_P1, lw=1.4)
ax0.annotate('', xy=(0.58, 0.64), xytext=(0.36, 0.58),
             arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.4, ls='--'))
ax0.plot(0.58, 0.64, 'o', color='#222222', markersize=4.0, zorder=5)
ax0.text(0.57, 0.71, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

# Clonal fan to contemporary isolates
for yt in [0.82, 0.70, 0.58, 0.46]:
    ax0.plot([0.58, 0.76], [0.64, yt], color=C_REC, lw=1.2)
    ax0.plot(0.76, yt, 'o', color=C_REC, markersize=3.0)

ax0.text(0.80, 0.64, r'$R_1 \dots R_K$', fontsize=7.2, fontweight='bold', color=C_REC, va='center')

# Lower Portion: Minimal Chromosome Track with Micro-Conversion Tract
bar_y = 0.08
bar_h = 0.13
w_bar = 0.90
x_bar = 0.05

# 0 to 2100 nt (42% of 5kb) - P1 Home
w1 = 0.42 * w_bar
ax0.add_patch(patches.Rectangle((x_bar, bar_y), w1, bar_h, facecolor=C_P1, edgecolor='none'))
ax0.text(x_bar + w1/2, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

# Micro-tract: 2100 to 2500 nt (e.g. 400 nt, 8% of 5kb) - P2 Donor
w2 = 0.10 * w_bar
ax0.add_patch(patches.Rectangle((x_bar + w1, bar_y), w2, bar_h, facecolor=C_P2, edgecolor='none'))
ax0.text(x_bar + w1 + w2/2, bar_y + bar_h + 0.02, r'$L_{\mathrm{tract}}$', color=C_P2, fontsize=6.0, fontweight='bold', ha='center', va='bottom')

# 2500 to 5000 nt (48% of 5kb) - P1 Home
w3 = 0.48 * w_bar
ax0.add_patch(patches.Rectangle((x_bar + w1 + w2, bar_y), w3, bar_h, facecolor=C_P1, edgecolor='none'))
ax0.text(x_bar + w1 + w2 + w3/2, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

# Track border
track_box = patches.Rectangle((x_bar, bar_y), w_bar, bar_h, fill=False, edgecolor='#222222', lw=0.7)
ax0.add_patch(track_box)

# Coordinate ticks below bar (matching Rows 1-3 style)
ax0.plot([x_bar, x_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([x_bar + w1, x_bar + w1], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([x_bar + w_bar, x_bar + w_bar], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.text(x_bar, bar_y - 0.038, '0', fontsize=5.5, ha='left', va='top', color='#555555')
ax0.text(x_bar + w1, bar_y - 0.038, '2.1k', fontsize=5.5, ha='center', va='top', color='#555555')
ax0.text(x_bar + w_bar, bar_y - 0.038, '5.0 kb', fontsize=5.5, ha='right', va='top', color='#555555')


# ==============================================================================
# 2. Column 2: Detection Power (%) vs Tract Length (L_tract)
# ==============================================================================
ax1 = fig.add_subplot(gs[0, 1])
clean_axis(ax1, 'Detection Power vs Tract Length', 'Tract Length (nt)', 'Sensitivity (%)', ylim=(-3, 105))
ax1.axhline(100.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

tract_lens = [50, 100, 200, 400, 800]
x_t = np.arange(len(tract_lens))
ax1.set_xticks(x_t)
ax1.set_xticklabels([str(t) for t in tract_lens], fontsize=6.8)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub = sum_df[np.isclose(sum_df['d_parent'], d)]
    
    # 3SEQ Power
    s_power = [sub[(sub['tool'] == '3SEQ') & (sub['tract_len'] == t)]['recall'].mean() * 100 for t in tract_lens]
    ax1.plot(x_t, s_power, 's:', color=c, lw=1.3, markersize=4.8, markerfacecolor='none', markeredgewidth=1.3, zorder=3)
    
    # RhizAeon Power
    r_power = [sub[(sub['tool'] == 'RhizAeon') & (sub['tract_len'] == t)]['recall'].mean() * 100 for t in tract_lens]
    ax1.plot(x_t, r_power, 'o-', color=c, lw=1.5, markersize=4.5, zorder=4)


# ==============================================================================
# 3. Column 3: Background Parent False Alarms (%) vs Tract Length (L_tract)
# ==============================================================================
ax2 = fig.add_subplot(gs[0, 2])
clean_axis(ax2, 'Background Taxon FPR (%)', 'Tract Length (nt)', 'False Alarms (%)', ylim=(-3, 105))
ax2.axhline(0.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

ax2.set_xticks(x_t)
ax2.set_xticklabels([str(t) for t in tract_lens], fontsize=6.8)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub = sum_df[np.isclose(sum_df['d_parent'], d)]
    
    # 3SEQ BG FPR
    s_fpr = [sub[(sub['tool'] == '3SEQ') & (sub['tract_len'] == t)]['bg_fpr'].mean() * 100 for t in tract_lens]
    ax2.plot(x_t, s_fpr, 's:', color=c, lw=1.3, markersize=4.8, markerfacecolor='none', markeredgewidth=1.3, zorder=3)
    
    # RhizAeon BG FPR
    r_fpr = [sub[(sub['tool'] == 'RhizAeon') & (sub['tract_len'] == t)]['bg_fpr'].mean() * 100 for t in tract_lens]
    ax2.plot(x_t, r_fpr, 'o-', color=c, lw=1.5, markersize=4.5, zorder=4)


# ==============================================================================
# 4. Column 4: Inferred Events (\hat{K}) vs Swarm Size (K) at L_tract = 400 nt
# ==============================================================================
ax3 = fig.add_subplot(gs[0, 3])
clean_axis(ax3, r'Inferred Events ($\hat{K}$)', r'Descendant Swarm Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-1.5, 45))
ax3.axhline(1.0, color='#222222', linestyle='-', linewidth=1.4, zorder=2) # True Ground Truth = 1.0

ks = [2, 5, 10]
ax3.set_xticks(ks)

sub_detectable = raw_df[raw_df['tract_len'] >= 200]

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    
    # 3SEQ Event Inflation IQR
    raw_3 = sub_detectable[(sub_detectable['tool'] == '3SEQ') & (np.isclose(sub_detectable['d_parent'], d))]
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
    raw_r = sub_detectable[(sub_detectable['tool'] == 'RhizAeon') & (np.isclose(sub_detectable['d_parent'], d))]
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
# Unified Top-Hoisted Legend (Identical to Rows 1-3)
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
out_png = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_panel_row4_scenario14.png'
out_pdf = '/Users/sergei/Projects/TOGA_MEME/recombination/viz/fig_panel_row4_scenario14.pdf'
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)

paper_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_panel_row4_scenario14.png'
paper_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_panel_row4_scenario14.pdf'
plt.savefig(paper_png, dpi=300)
plt.savefig(paper_pdf)
plt.close()
print("Row 4 successfully generated and saved.")
