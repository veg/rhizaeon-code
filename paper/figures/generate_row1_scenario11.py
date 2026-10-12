import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

# Load authentic benchmark data
df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/epidemic_summary.csv')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

# Calibrated dimensions: 11.8 x 3.0 in, 300 DPI
fig = plt.figure(figsize=(11.8, 3.0), dpi=300)
gs = fig.add_gridspec(1, 4, width_ratios=[1.22, 0.92, 0.92, 0.92], wspace=0.32,
                      left=0.035, right=0.985, top=0.78, bottom=0.18)

C_P1 = '#1f77b4'       # Home / Parent 1 (Blue)
C_P2 = '#d95f02'       # Donor / Parent 2 (Orange)
C_REC = '#0088cc'      # Recombinant swarm (Cyan/Blue)
colors_div = {
    0.02: '#4dbbd5',   # light cyan-blue
    0.05: '#1f77b4',   # mid blue
    0.10: '#002b66',   # dark navy
}

# ==============================================================================
# 1. Minimal & Clear Schematic: Epidemic Reticulation & Clonal Expansion
# ==============================================================================
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.axis('off')

# Panel Title & minimal metadata
ax0.text(0.0, 1.22, 'A  Scenario 1.1: Epidemic Clonal Swarms', fontsize=9.5, fontweight='bold', va='top')
ax0.text(0.0, 1.07, r'Ancestral reticulation with clonal expansion ($N = 40 + K$, $L = 5\ \mathrm{kb}$)',
         fontsize=6.8, color='#555555', va='top')

# Minimal Tree Schematic
# P1 Home Clade (Top)
ax0.plot([0.03, 0.12], [0.82, 0.82], color=C_P1, lw=1.6)
ax0.plot([0.12, 0.12], [0.70, 0.94], color=C_P1, lw=1.6)
ax0.plot([0.12, 0.28], [0.94, 0.94], color=C_P1, lw=1.3)
ax0.plot([0.12, 0.32], [0.70, 0.70], color=C_P1, lw=1.5)
ax0.text(0.30, 0.94, r'Background ($N_{\mathrm{bg}}=40$)', va='center', fontsize=6.6, color='#444444')
ax0.text(0.14, 0.74, r'$P_1$ (Home)', va='bottom', fontsize=7.2, fontweight='bold', color=C_P1)

# P2 Donor Clade (Bottom)
ax0.plot([0.03, 0.12], [0.42, 0.42], color=C_P2, lw=1.6)
ax0.plot([0.12, 0.12], [0.30, 0.54], color=C_P2, lw=1.6)
ax0.plot([0.12, 0.32], [0.54, 0.54], color=C_P2, lw=1.5)
ax0.plot([0.12, 0.28], [0.30, 0.30], color=C_P2, lw=1.3)
ax0.text(0.14, 0.58, r'$P_2$ (Donor)', va='bottom', fontsize=7.2, fontweight='bold', color=C_P2)
ax0.text(0.30, 0.30, 'Donor Relatives', va='center', fontsize=6.6, color='#444444')

# Transfer to Ancestral Node
ax0.plot([0.32, 0.56], [0.70, 0.58], color=C_P1, lw=1.4)
ax0.annotate('', xy=(0.56, 0.58), xytext=(0.32, 0.54),
             arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.4, ls='--'))
ax0.text(0.42, 0.47, 'Transfer', color=C_P2, fontsize=6.3, fontweight='bold', ha='center', va='top')

# Reticulation Node Dot
ax0.plot(0.56, 0.58, 'o', color='#222222', markersize=4.2, zorder=5)
ax0.text(0.54, 0.65, r'$R_{\mathrm{anc}}$', fontsize=7.2, fontweight='bold', ha='right')

# Clonal Expansion Fan
y_tips = [0.80, 0.65, 0.51, 0.36]
for yt in y_tips:
    ax0.plot([0.56, 0.76], [0.58, yt], color=C_REC, lw=1.2)
    ax0.plot(0.76, yt, 'o', color=C_REC, markersize=3.2)

ax0.text(0.80, 0.66, 'Clonal Swarm', fontsize=6.8, fontweight='bold', color=C_REC, va='center')
ax0.text(0.80, 0.49, r'$R_1 \dots R_K$ ($K \leq 15$)', fontsize=6.4, color=C_REC, va='center')

# Minimal Chromosome Track
ax0.text(0.0, 0.21, r'Chromosome: Single Cassette ($K_{\mathrm{true}} = 1$ across all isolates)',
         fontsize=6.4, color='#333333', va='bottom')

bar_y = 0.05
bar_h = 0.12
# Segment 1: P1 (Home) 0 to 1500 nt (30%)
rect1 = patches.Rectangle((0.0, bar_y), 0.30, bar_h, facecolor=C_P1, edgecolor='none')
ax0.add_patch(rect1)
ax0.text(0.15, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

# Segment 2: P2 (Donor) 1500 to 3000 nt (30%)
rect2 = patches.Rectangle((0.30, bar_y), 0.30, bar_h, facecolor=C_P2, edgecolor='none')
ax0.add_patch(rect2)
ax0.text(0.45, bar_y + bar_h/2, r'$P_2$ (1.5 kb)', color='white', fontsize=6.2, fontweight='bold', ha='center', va='center')

# Segment 3: P1 (Home) 3000 to 5000 nt (40%)
rect3 = patches.Rectangle((0.60, bar_y), 0.40, bar_h, facecolor=C_P1, edgecolor='none')
ax0.add_patch(rect3)
ax0.text(0.80, bar_y + bar_h/2, r'$P_1$', color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

# Track border
track_box = patches.Rectangle((0.0, bar_y), 1.0, bar_h, fill=False, edgecolor='#222222', lw=0.7)
ax0.add_patch(track_box)

# Coordinate ticks below bar
ax0.plot([0.0, 0.0], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([0.30, 0.30], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([0.60, 0.60], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.plot([1.0, 1.0], [bar_y, bar_y - 0.025], color='#333333', lw=0.6)
ax0.text(0.0, bar_y - 0.038, '1 nt', fontsize=5.5, ha='left', va='top', color='#555555')
ax0.text(0.30, bar_y - 0.038, '1.5k', fontsize=5.5, ha='center', va='top', color='#555555')
ax0.text(0.60, bar_y - 0.038, '3.0k', fontsize=5.5, ha='center', va='top', color='#555555')
ax0.text(1.0, bar_y - 0.038, '5.0 kb', fontsize=5.5, ha='right', va='top', color='#555555')


# ==============================================================================
# Helper for Clean Axes
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
    ax.set_xticks([2, 5, 8, 10, 15])


# ==============================================================================
# 2. Metric 1: Alignment-Level Power (%) — Pure data curves (NO callouts)
# ==============================================================================
ax1 = fig.add_subplot(gs[0, 1])
clean_axis(ax1, 'Alignment-Level Power', r'Descendant Swarm Size ($K$)', 'Sensitivity (%)', ylim=(-3, 105))
ax1.axhline(100.0, color='#666666', linestyle='--', linewidth=0.8, alpha=0.5, zorder=2)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_rhiz = df[(df['tool'] == 'RhizAeon') & (np.isclose(df['d_parent'], d))].sort_values('K')
    sub_3seq = df[(df['tool'] == '3SEQ') & (np.isclose(df['d_parent'], d))].sort_values('K')
    
    # 3SEQ (open square)
    ax1.plot(sub_3seq['K'], sub_3seq['recall_mean'] * 100.0,
             marker='s', markerfacecolor='none', markeredgecolor=c, markeredgewidth=1.3,
             markersize=4.8, linestyle=':', linewidth=1.2, color=c, zorder=3)
    # RhizAeon (filled circle)
    ax1.plot(sub_rhiz['K'], sub_rhiz['recall_mean'] * 100.0,
             marker='o', markerfacecolor=c, markeredgecolor='white', markeredgewidth=0.5,
             markersize=5.0, linestyle='-', linewidth=1.6, color=c, zorder=4)


# ==============================================================================
# 3. Metric 2: Inferred Event Cardinality (\hat{K}) — Pure data curves (NO callouts)
# ==============================================================================
ax2 = fig.add_subplot(gs[0, 2])
clean_axis(ax2, r'Inferred Event Cardinality ($\hat{K}$)', r'Descendant Swarm Size ($K$)', r'Inferred Events ($\hat{K}$)', ylim=(-2, 60))
ax2.axhline(1.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_rhiz = df[(df['tool'] == 'RhizAeon') & (np.isclose(df['d_parent'], d))].sort_values('K')
    sub_3seq = df[(df['tool'] == '3SEQ') & (np.isclose(df['d_parent'], d))].sort_values('K')
    
    # 3SEQ (open square)
    ax2.plot(sub_3seq['K'], sub_3seq['unique_events_mean'],
             marker='s', markerfacecolor='none', markeredgecolor=c, markeredgewidth=1.3,
             markersize=4.8, linestyle=':', linewidth=1.2, color=c, zorder=3)
    # RhizAeon (filled circle)
    ax2.plot(sub_rhiz['K'], sub_rhiz['unique_events_mean'],
             marker='o', markerfacecolor=c, markeredgecolor='white', markeredgewidth=0.5,
             markersize=5.0, linestyle='-', linewidth=1.6, color=c, zorder=4)


# ==============================================================================
# 4. Metric 3: Background False Alarms (%) — Pure data curves (NO callouts)
# ==============================================================================
ax3 = fig.add_subplot(gs[0, 3])
clean_axis(ax3, 'Innocent Background False Alarms', r'Descendant Swarm Size ($K$)', 'Background FPR (%)', ylim=(-3, 105))
ax3.axhline(0.0, color='#222222', linestyle='--', linewidth=1.1, zorder=2)

for d in [0.02, 0.05, 0.10]:
    c = colors_div[d]
    sub_rhiz = df[(df['tool'] == 'RhizAeon') & (np.isclose(df['d_parent'], d))].sort_values('K')
    sub_3seq = df[(df['tool'] == '3SEQ') & (np.isclose(df['d_parent'], d))].sort_values('K')
    
    # 3SEQ (open square)
    ax3.plot(sub_3seq['K'], sub_3seq['bg_fpr_mean'] * 100.0,
             marker='s', markerfacecolor='none', markeredgecolor=c, markeredgewidth=1.3,
             markersize=4.8, linestyle=':', linewidth=1.2, color=c, zorder=3)
    # RhizAeon (filled circle)
    ax3.plot(sub_rhiz['K'], sub_rhiz['bg_fpr_mean'] * 100.0,
             marker='o', markerfacecolor=c, markeredgecolor='white', markeredgewidth=0.5,
             markersize=5.0, linestyle='-', linewidth=1.6, color=c, zorder=4)


# ==============================================================================
# Unified Shared Legend in Header
# ==============================================================================
h_rhiz = mlines.Line2D([], [], color='#002b66', marker='o', markerfacecolor='#002b66',
                       markeredgecolor='white', markeredgewidth=0.5, markersize=5.2, lw=1.6,
                       label='RhizAeon (Continuous Manifold)')
h_3seq = mlines.Line2D([], [], color='#002b66', marker='s', markerfacecolor='none',
                       markeredgecolor='#002b66', markeredgewidth=1.3, markersize=4.8,
                       linestyle=':', lw=1.2, label='3SEQ (Pairwise Triplets)')
h_gt = mlines.Line2D([], [], color='#222222', linestyle='--', lw=1.1,
                     label=r'Ground Truth ($K_{\mathrm{true}}=1$, $\mathrm{FPR}=0\%$)')
p_d02 = patches.Patch(facecolor=colors_div[0.02], edgecolor='none', label=r'$d = 0.02$')
p_d05 = patches.Patch(facecolor=colors_div[0.05], edgecolor='none', label=r'$d = 0.05$')
p_d10 = patches.Patch(facecolor=colors_div[0.10], edgecolor='none', label=r'$d = 0.10$')

fig.legend(handles=[h_rhiz, h_3seq, h_gt, p_d02, p_d05, p_d10],
           loc='upper center', bbox_to_anchor=(0.58, 0.99),
           ncol=6, frameon=False, fontsize=7.0, columnspacing=1.3, handletextpad=0.4)

out_png = '/Users/sergei/.gemini/antigravity-cli/brain/1dfa5602-784e-4936-ae1c-a9734d91bc95/fig_panel_row1_scenario11.png'
plt.savefig(out_png, dpi=300)
plt.close()
print("Saved final publication style Row 1 to:", out_png)
