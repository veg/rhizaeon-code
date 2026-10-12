"""
Generate Figure 2: Scenario 1.1 Epidemic Clonal Expansion Benchmark
Compares RhizAeon vs 3SEQ on:
  Panel A: Inferred Event Cardinality (\hat{K}) vs Descendant Swarm Size (K)
  Panel B: Background Parent False Positive Rate (%) vs Descendant Swarm Size (K)
  Panel C: Recombinant Isolate Precision (%) vs Descendant Swarm Size (K)

Strict BioVis guidelines:
- Publication ready (300 DPI, vector PDF + PNG)
- No narrative callout boxes / text boxes inside plot area
- Hoisted clean shared legend in top margin (zero data overlap)
- Ground truth plotted as reference line
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

# Load authentic benchmark data
df = pd.read_csv('/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/epidemic_summary.csv')

# Set style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), dpi=300)

# Color scheme for divergence levels d_parent
colors = {
    0.02: '#0099cc',  # cyan / light blue
    0.05: '#1f77b4',  # medium blue
    0.10: '#002b66',  # dark navy
}

swarm_sizes = [2, 5, 8, 10, 15]

# --- Panel A: Inferred Event Cardinality ---
ax_a = axes[0]
ax_a.axhline(1.0, color='#222222', linestyle='--', linewidth=1.5, zorder=2)

for d in [0.02, 0.05, 0.10]:
    c = colors[d]
    sub_rhiz = df[(df['tool'] == 'RhizAeon') & (np.isclose(df['d_parent'], d))].sort_values('K')
    sub_3seq = df[(df['tool'] == '3SEQ') & (np.isclose(df['d_parent'], d))].sort_values('K')
    
    # 3SEQ
    ax_a.plot(sub_3seq['K'], sub_3seq['unique_events_mean'],
              marker='s', markersize=5.5, linestyle=':', linewidth=1.6, color=c, zorder=3)
    # RhizAeon
    ax_a.plot(sub_rhiz['K'], sub_rhiz['unique_events_mean'],
              marker='o', markersize=6.0, linestyle='-', linewidth=2.0, color=c, zorder=4)

ax_a.set_title('A. Inferred Event Cardinality ($\hat{K}$)', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_a.set_xlabel('Descendant Swarm Size ($K$)', fontsize=10, fontweight='bold')
ax_a.set_ylabel('Inferred Events ($\hat{K}$)', fontsize=10, fontweight='bold')
ax_a.set_xticks(swarm_sizes)
ax_a.set_ylim(-1, 60)
ax_a.grid(True, linestyle='--', alpha=0.4, color='#cccccc')

# --- Panel B: Background Parent False Positive Rate (%) ---
ax_b = axes[1]
ax_b.axhline(0.0, color='#888888', linestyle='-', linewidth=0.8, alpha=0.5, zorder=1)

for d in [0.02, 0.05, 0.10]:
    c = colors[d]
    sub_rhiz = df[(df['tool'] == 'RhizAeon') & (np.isclose(df['d_parent'], d))].sort_values('K')
    sub_3seq = df[(df['tool'] == '3SEQ') & (np.isclose(df['d_parent'], d))].sort_values('K')
    
    # 3SEQ
    ax_b.plot(sub_3seq['K'], sub_3seq['bg_fpr_mean'] * 100.0,
              marker='s', markersize=5.5, linestyle=':', linewidth=1.6, color=c, zorder=3)
    # RhizAeon
    ax_b.plot(sub_rhiz['K'], sub_rhiz['bg_fpr_mean'] * 100.0,
              marker='o', markersize=6.0, linestyle='-', linewidth=2.0, color=c, zorder=4)

ax_b.set_title('B. Background Parent False Positive Rate', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_b.set_xlabel('Descendant Swarm Size ($K$)', fontsize=10, fontweight='bold')
ax_b.set_ylabel('Parent False Positive Rate (%)', fontsize=10, fontweight='bold')
ax_b.set_xticks(swarm_sizes)
ax_b.set_ylim(-3, 105)
ax_b.grid(True, linestyle='--', alpha=0.4, color='#cccccc')

# --- Panel C: Recombinant Isolate Precision (%) ---
ax_c = axes[2]

for d in [0.02, 0.05, 0.10]:
    c = colors[d]
    sub_rhiz = df[(df['tool'] == 'RhizAeon') & (np.isclose(df['d_parent'], d))].sort_values('K')
    sub_3seq = df[(df['tool'] == '3SEQ') & (np.isclose(df['d_parent'], d))].sort_values('K')
    
    # 3SEQ
    ax_c.plot(sub_3seq['K'], sub_3seq['precision_mean'] * 100.0,
              marker='s', markersize=5.5, linestyle=':', linewidth=1.6, color=c, zorder=3)
    # RhizAeon
    ax_c.plot(sub_rhiz['K'], sub_rhiz['precision_mean'] * 100.0,
              marker='o', markersize=6.0, linestyle='-', linewidth=2.0, color=c, zorder=4)

ax_c.set_title('C. Recombinant Isolate Precision', fontsize=11, fontweight='bold', loc='left', pad=10)
ax_c.set_xlabel('Descendant Swarm Size ($K$)', fontsize=10, fontweight='bold')
ax_c.set_ylabel('Isolate Precision (%)', fontsize=10, fontweight='bold')
ax_c.set_xticks(swarm_sizes)
ax_c.set_ylim(-3, 105)
ax_c.grid(True, linestyle='--', alpha=0.4, color='#cccccc')

# Shared legend in top margin
gt_handle = mlines.Line2D([], [], color='#222222', linestyle='--', linewidth=1.5,
                          label=r'Ground Truth ($K_{\mathrm{true}}=1$)')
rhiz_handle = mlines.Line2D([], [], color='#1f77b4', marker='o', markersize=6.0, linestyle='-', linewidth=2.0,
                            label=r'RhizAeon (Continuous Roster)')
seq3_handle = mlines.Line2D([], [], color='#1f77b4', marker='s', markersize=5.5, linestyle=':', linewidth=1.6,
                            label=r'3SEQ (Triplet Scanning)')

d02_handle = mlines.Line2D([], [], color=colors[0.02], linewidth=3.0, label=r'$d_{\mathrm{parent}} = 0.02$')
d05_handle = mlines.Line2D([], [], color=colors[0.05], linewidth=3.0, label=r'$d_{\mathrm{parent}} = 0.05$')
d10_handle = mlines.Line2D([], [], color=colors[0.10], linewidth=3.0, label=r'$d_{\mathrm{parent}} = 0.10$')

fig.legend(handles=[gt_handle, rhiz_handle, seq3_handle, d02_handle, d05_handle, d10_handle],
           loc='upper center', bbox_to_anchor=(0.5, 1.05),
           ncol=6, frameon=False, fontsize=9.0, columnspacing=1.5)

plt.tight_layout(rect=[0, 0, 1, 0.93])

output_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_scenario1_epidemic_swarms.pdf'
output_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_scenario1_epidemic_swarms.png'
plt.savefig(output_pdf, bbox_inches='tight')
plt.savefig(output_png, bbox_inches='tight')
print(f"Generated {output_pdf} and {output_png}")
