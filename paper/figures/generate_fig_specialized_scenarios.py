"""
Generate Figure 2: Specialized Genomic and Epidemiological Regimes
Multi-panel figure covering Scenarios 1.1 to 1.4:
  Panel A: Epidemic Clonal Expansion (Scenario 1.1) - Inferred Event Cardinality vs Swarm Size
  Panel B: Multi-Way Mosaicism (Scenario 1.2) - Donors Correctly Resolved across Complex Topologies
  Panel C: Ghost Donor Introgression (Scenario 1.3) - Donor Attribution & Innocent Parent False Framing
  Panel D: Bacterial Micro-Conversions (Scenario 1.4) - Scale-Space Detection Power vs Tract Length

Strict MBE BioVis Standards:
- Clean 2x2 grid layout (publication ready, 300 DPI, vector PDF + PNG)
- No narrative callout boxes or text annotations inside plot areas
- Hoisted clean shared legend in top margin (zero data overlap)
- Ground truth plotted as reference line
- Zero synthetic mockery: strictly authentic executed simulation data
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import os

# Load authentic benchmark data
df_s1 = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/epidemic_summary.csv')
df_s2 = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms/multiway_summary.csv')
df_s3 = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario3_ghost_donor_introgression/ghost_summary.csv')
df_s4 = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario4_micro_conversion_tracts/micro_conversion_summary.csv')

# Set style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.2), dpi=300)

colors = {
    0.02: '#0099cc',  # cyan / light blue
    0.05: '#1f77b4',  # medium blue
    0.10: '#002b66',  # dark navy
}

# ==============================================================================
# Panel A: Scenario 1.1 - Epidemic Clonal Expansion
# ==============================================================================
ax_a = axes[0, 0]
swarm_sizes = [2, 5, 8, 10, 15]

ax_a.axhline(1.0, color='#222222', linestyle='--', linewidth=1.5, zorder=2)

for d in [0.02, 0.05, 0.10]:
    c = colors[d]
    sub_rhiz = df_s1[(df_s1['tool'] == 'RhizAeon') & (np.isclose(df_s1['d_parent'], d))].sort_values('K')
    sub_3seq = df_s1[(df_s1['tool'] == '3SEQ') & (np.isclose(df_s1['d_parent'], d))].sort_values('K')
    
    # 3SEQ
    ax_a.plot(sub_3seq['K'], sub_3seq['unique_events_mean'],
              marker='s', markersize=5.5, linestyle=':', linewidth=1.6, color=c, zorder=3)
    # RhizAeon
    ax_a.plot(sub_rhiz['K'], sub_rhiz['unique_events_mean'],
              marker='o', markersize=6.0, linestyle='-', linewidth=2.0, color=c, zorder=4)

ax_a.set_title('A. Epidemic Clonal Swarms (Scenario 1.1)', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_a.set_xlabel(r'Descendant Swarm Size ($K$)', fontsize=10, fontweight='bold')
ax_a.set_ylabel(r'Inferred Events ($\hat{K}$)', fontsize=10, fontweight='bold')
ax_a.set_xticks(swarm_sizes)
ax_a.set_ylim(-1, 60)
ax_a.grid(True, linestyle='--', alpha=0.4, color='#cccccc')


# ==============================================================================
# Panel B: Scenario 1.2 - Multi-Way Mosaic Swarms
# ==============================================================================
ax_b = axes[0, 1]
topologies = ['dual_cassette', 'nested_russian_doll', 'triple_crossover']
topo_labels = [
    'Dual Cassette\n($K_{\\mathrm{donor}}=2$)',
    'Nested Doll\n($K_{\\mathrm{donor}}=2$)',
    'Triple Crossover\n($K_{\\mathrm{donor}}=2$)'
]
x = np.arange(len(topologies))
width = 0.35

# Aggregate at d=0.05 across swarm sizes K in {1, 5, 10}
rhiz_donors = []
seq3_donors = []
for topo in topologies:
    sub_r = df_s2[(df_s2['tool'] == 'RhizAeon') & (df_s2['topology'] == topo) & (np.isclose(df_s2['d_parent'], 0.05))]
    sub_3 = df_s2[(df_s2['tool'] == '3SEQ') & (df_s2['topology'] == topo) & (np.isclose(df_s2['d_parent'], 0.05))]
    rhiz_donors.append(sub_r['mean_donors_resolved'].mean())
    seq3_donors.append(sub_3['mean_donors_resolved'].mean())

rects1 = ax_b.bar(x - width/2, rhiz_donors, width, label='RhizAeon', color='#1f77b4', edgecolor='#002b66', linewidth=1.0, zorder=3)
rects2 = ax_b.bar(x + width/2, seq3_donors, width, label='3SEQ', color='#aec7e8', edgecolor='#1f77b4', linewidth=1.0, hatch='//', zorder=3)
ax_b.axhline(2.0, color='#222222', linestyle='--', linewidth=1.5, zorder=2)

ax_b.set_title('B. Multi-Way Mosaic Swarms (Scenario 1.2)', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_b.set_ylabel('Donors Correctly Resolved', fontsize=10, fontweight='bold')
ax_b.set_xticks(x)
ax_b.set_xticklabels(topo_labels, fontsize=9.0)
ax_b.set_ylim(0, 2.5)
ax_b.grid(True, linestyle='--', alpha=0.4, color='#cccccc', axis='y')

# Add bar value labels
for rect in rects1:
    h = rect.get_height()
    ax_b.annotate(f'{h:.2f}', xy=(rect.get_x() + rect.get_width()/2, h),
                  xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8.5, fontweight='bold')
for rect in rects2:
    h = rect.get_height()
    ax_b.annotate(f'{h:.2f}', xy=(rect.get_x() + rect.get_width()/2, h),
                  xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8.5)


# ==============================================================================
# Panel C: Scenario 1.3 - Ghost Donor Introgression
# ==============================================================================
ax_c = axes[1, 0]
categories_c = [
    'Surrogate Donor\nAttribution Accuracy',
    'False Framing\n(Surrogate Sister)',
    'False Framing\n(Pure Ghost)',
    'False Framing\n(Deep Reservoir)'
]
xc = np.arange(len(categories_c))

# Metrics from authentic benchmark summary
# 1. Donor attribution accuracy in surrogate sister
rhiz_acc = df_s3[(df_s3['regime'] == 'surrogate_sister') & (df_s3['tool'] == 'RhizAeon')]['mean_donor_acc'].mean() * 100.0
seq3_acc = df_s3[(df_s3['regime'] == 'surrogate_sister') & (df_s3['tool'] == '3SEQ')]['mean_donor_acc'].mean() * 100.0

# 2. Parent FPR in surrogate sister
rhiz_fpr_surr = df_s3[(df_s3['regime'] == 'surrogate_sister') & (df_s3['tool'] == 'RhizAeon')]['mean_parent_fpr'].mean() * 100.0
seq3_fpr_surr = df_s3[(df_s3['regime'] == 'surrogate_sister') & (df_s3['tool'] == '3SEQ')]['mean_parent_fpr'].mean() * 100.0

# 3. Parent FPR in pure ghost
rhiz_fpr_pure = df_s3[(df_s3['regime'] == 'pure_ghost') & (df_s3['tool'] == 'RhizAeon')]['mean_parent_fpr'].mean() * 100.0
seq3_fpr_pure = df_s3[(df_s3['regime'] == 'pure_ghost') & (df_s3['tool'] == '3SEQ')]['mean_parent_fpr'].mean() * 100.0

# 4. Parent FPR in deep reservoir
rhiz_fpr_deep = df_s3[(df_s3['regime'] == 'deep_reservoir') & (df_s3['tool'] == 'RhizAeon')]['mean_parent_fpr'].mean() * 100.0
seq3_fpr_deep = df_s3[(df_s3['regime'] == 'deep_reservoir') & (df_s3['tool'] == '3SEQ')]['mean_parent_fpr'].mean() * 100.0

rhiz_vals_c = [rhiz_acc, rhiz_fpr_surr, rhiz_fpr_pure, rhiz_fpr_deep]
seq3_vals_c = [seq3_acc, seq3_fpr_surr, seq3_fpr_pure, seq3_fpr_deep]

rects_c1 = ax_c.bar(xc - width/2, rhiz_vals_c, width, label='RhizAeon', color='#1f77b4', edgecolor='#002b66', linewidth=1.0, zorder=3)
rects_c2 = ax_c.bar(xc + width/2, seq3_vals_c, width, label='3SEQ', color='#aec7e8', edgecolor='#1f77b4', linewidth=1.0, hatch='//', zorder=3)

# Subtle vertical divider separating attribution accuracy from false framing
ax_c.axvline(0.5, color='#999999', linestyle=':', linewidth=1.0, zorder=2)

ax_c.set_title('C. Ghost Donor Introgression (Scenario 1.3)', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_c.set_ylabel('Percentage (%)', fontsize=10, fontweight='bold')
ax_c.set_xticks(xc)
ax_c.set_xticklabels(categories_c, fontsize=8.2)
ax_c.set_ylim(0, 118)
ax_c.grid(True, linestyle='--', alpha=0.4, color='#cccccc', axis='y')

for rect in rects_c1:
    h = rect.get_height()
    ax_c.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h),
                  xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8.0, fontweight='bold')
for rect in rects_c2:
    h = rect.get_height()
    ax_c.annotate(f'{h:.1f}%', xy=(rect.get_x() + rect.get_width()/2, h),
                  xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8.0)


# ==============================================================================
# Panel D: Scenario 1.4 - Bacterial Micro-Conversions
# ==============================================================================
ax_d = axes[1, 1]
tract_lengths = [50, 100, 200, 400, 800]

ax_d.set_xscale('log')

for d in [0.02, 0.05, 0.10]:
    c = colors[d]
    sub = df_s4[np.isclose(df_s4['d_parent'], d)]
    r_recall = [sub[(sub['tool'] == 'RhizAeon') & (sub['tract_len'] == tl)]['recall'].mean() * 100.0 for tl in tract_lengths]
    s_recall = [sub[(sub['tool'] == '3SEQ') & (sub['tract_len'] == tl)]['recall'].mean() * 100.0 for tl in tract_lengths]
    
    # 3SEQ
    ax_d.plot(tract_lengths, s_recall, marker='s', markersize=5.5, linestyle=':', linewidth=1.6, color=c, zorder=3)
    # RhizAeon
    ax_d.plot(tract_lengths, r_recall, marker='o', markersize=6.0, linestyle='-', linewidth=2.0, color=c, zorder=4)

ax_d.set_title('D. Bacterial Micro-Conversions (Scenario 1.4)', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_d.set_xlabel('Conversion Tract Length (nt)', fontsize=10, fontweight='bold')
ax_d.set_ylabel('Detection Power / Recall (%)', fontsize=10, fontweight='bold')
ax_d.set_xticks(tract_lengths)
ax_d.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax_d.set_ylim(-5, 105)
ax_d.grid(True, linestyle='--', alpha=0.4, color='#cccccc')


# ==============================================================================
# Unified Shared Legend in Top Margin
# ==============================================================================
gt_handle = mlines.Line2D([], [], color='#222222', linestyle='--', linewidth=1.5,
                          label=r'Ground Truth')
rhiz_line = mlines.Line2D([], [], color='#1f77b4', marker='o', markersize=6.0, linestyle='-', linewidth=2.0,
                          label='RhizAeon')
seq3_line = mlines.Line2D([], [], color='#1f77b4', marker='s', markersize=5.5, linestyle=':', linewidth=1.6,
                          label='3SEQ')

d02_patch = mpatches.Patch(color=colors[0.02], label=r'$d_{\mathrm{parent}} = 0.02$')
d05_patch = mpatches.Patch(color=colors[0.05], label=r'$d_{\mathrm{parent}} = 0.05$')
d10_patch = mpatches.Patch(color=colors[0.10], label=r'$d_{\mathrm{parent}} = 0.10$')

fig.legend(handles=[gt_handle, rhiz_line, seq3_line, d02_patch, d05_patch, d10_patch],
           loc='upper center', bbox_to_anchor=(0.5, 0.99),
           ncol=6, frameon=False, fontsize=9.5, columnspacing=1.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])

output_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_specialized_scenarios.pdf'
output_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_specialized_scenarios.png'
plt.savefig(output_pdf, bbox_inches='tight')
plt.savefig(output_png, bbox_inches='tight')
print(f"Generated {output_pdf} and {output_png}")
