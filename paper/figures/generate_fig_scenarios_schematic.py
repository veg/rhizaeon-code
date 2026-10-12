"""
Generate Publication Schematic: Specialized Simulation Regimes (Scenarios 1.1 to 1.4)
Explains the topological architecture, epidemiological dynamics, and evolutionary challenges:
  Panel A: Scenario 1.1 - Epidemic Clonal Expansion
  Panel B: Scenario 1.2 - Multi-Way Mosaic Swarms
  Panel C: Scenario 1.3 - Ghost Donor Introgression
  Panel D: Scenario 1.4 - Bacterial Micro-Conversion Tracts

Publication Standards:
- Clean 2x2 grid layout (12.2 x 9.4 in, 300 DPI)
- Relative [0, 1] coordinate system with zero overlapping annotations
- Accessible qualitative palettes (Blue P1, Orange P2, Green P3, Gray Ghost)
- Professional typography, generous margins, no ugly spine boundaries
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Set style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Helvetica', 'Arial']

fig, axes = plt.subplots(2, 2, figsize=(12.2, 9.4), dpi=300)

# Colors
C_P1 = '#1f77b4'       # Home / Parent 1 (Blue)
C_P2 = '#d95f02'       # Donor / Parent 2 (Orange)
C_P3 = '#2ca02c'       # Donor / Parent 3 (Green)
C_GHOST = '#777777'    # Ghost / Unsampled (Slate Gray)
C_REC = '#0288d1'      # Recombinant swarm (Cyan/Blue)
C_BG = '#ffffff'       # Clean white panel background

for ax in axes.flat:
    ax.set_facecolor(C_BG)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# ==============================================================================
# Panel A: Scenario 1.1 - Epidemic Clonal Expansion
# ==============================================================================
ax_a = axes[0, 0]
ax_a.text(0.02, 0.98, 'A. Epidemic Clonal Expansion (Scenario 1.1)', fontsize=11, fontweight='bold', va='top')

# 1. Transmission / Reticulation Tree
# P1 Clade (Home)
ax_a.plot([0.04, 0.12], [0.84, 0.84], color=C_P1, lw=2.0)
ax_a.plot([0.12, 0.12], [0.78, 0.90], color=C_P1, lw=2.0)
ax_a.plot([0.12, 0.22], [0.90, 0.90], color=C_P1, lw=2.0)
ax_a.plot([0.12, 0.22], [0.78, 0.78], color=C_P1, lw=2.0)
ax_a.text(0.24, 0.90, 'Background ($N_{\\mathrm{bg}}=40$)', va='center', fontsize=7.6, color='#555555')
ax_a.text(0.24, 0.78, 'Parent 1 ($P_1$, Home)', va='center', fontsize=8.0, fontweight='bold', color=C_P1)

# P2 Clade (Donor)
ax_a.plot([0.04, 0.12], [0.41, 0.41], color=C_P2, lw=2.0)
ax_a.plot([0.12, 0.12], [0.35, 0.47], color=C_P2, lw=2.0)
ax_a.plot([0.12, 0.22], [0.47, 0.47], color=C_P2, lw=2.0)
ax_a.plot([0.12, 0.22], [0.35, 0.35], color=C_P2, lw=2.0)
ax_a.text(0.24, 0.47, 'Parent 2 ($P_2$, Donor)', va='center', fontsize=8.0, fontweight='bold', color=C_P2)
ax_a.text(0.24, 0.35, 'Donor Relatives', va='center', fontsize=7.6, color='#555555')

# Recombination introgression: from P2 lineage up to R_anc
ax_a.annotate('', xy=(0.33, 0.63), xytext=(0.20, 0.47),
              arrowprops=dict(arrowstyle="->", color=C_P2, lw=1.8, linestyle='--', connectionstyle="arc3,rad=-0.12"))
ax_a.text(0.22, 0.56, 'Ancestral\nTransfer', fontsize=7.3, color=C_P2, fontweight='bold', ha='center')

# Recipient branch: from P1 lineage down to R_anc
ax_a.plot([0.18, 0.33], [0.78, 0.63], color=C_P1, lw=2.0, zorder=2)
ax_a.plot([0.33, 0.43], [0.63, 0.63], color=C_REC, lw=2.5, zorder=3)
ax_a.scatter([0.33], [0.63], color='#222222', s=35, zorder=5)
ax_a.text(0.33, 0.68, '$R_{\\mathrm{anc}}$', fontsize=8.5, fontweight='bold', ha='center')

# Descendant swarm expansion
swarm_y = [0.78, 0.71, 0.63, 0.55, 0.48]
for sy in swarm_y:
    ax_a.plot([0.43, 0.53], [0.63, sy], color=C_REC, lw=1.5)
    ax_a.scatter([0.53], [sy], color=C_REC, s=18, zorder=5)

ax_a.text(0.55, 0.78, '$R_1$', va='center', fontsize=7.6, color=C_REC, fontweight='bold')
ax_a.text(0.55, 0.71, '$R_2$', va='center', fontsize=7.6, color=C_REC, fontweight='bold')
ax_a.text(0.55, 0.63, '$\\dots$', va='center', fontsize=8.5, color=C_REC, fontweight='bold')
ax_a.text(0.55, 0.55, '$R_{K-1}$', va='center', fontsize=7.6, color=C_REC, fontweight='bold')
ax_a.text(0.55, 0.48, '$R_K$', va='center', fontsize=7.6, color=C_REC, fontweight='bold')

# Swarm Summary Box on the right
ax_a.text(0.64, 0.63,
          'Descendant Swarm ($K \\in \\{2, 5, 8, 10, 15\\}$)\n\n'
          '• True physical events: $K_{\\mathrm{true}} = 1$\n'
          '• Pairwise Triplet Tests (3SEQ):\n'
          '    $\\hat{K} = 10\\text{--}55$ false events\n'
          '• RhizAeon Bundle Tracking:\n'
          '    $\\hat{K} = 1.00$ (100% exact)',
          va='center', ha='left', fontsize=7.3, color='#111111',
          bbox=dict(boxstyle='round,pad=0.45', facecolor='#f4f9fd', edgecolor='#99ccee', lw=0.8))

# 2. Recombinant Chromosome
ax_a.text(0.04, 0.25, 'Shared Recombinant Genome ($L = 5{,}000$ nt) — Identical cassette across all $K$ isolates:', fontsize=7.8, fontweight='bold', color='#222222')

bx = 0.04
by = 0.11
bw = 0.90
bh = 0.075

w_p1a = bw * (1500 / 5000)
w_p2 = bw * (1500 / 5000)
w_p1b = bw * (2000 / 5000)

ax_a.add_patch(patches.Rectangle((bx, by), w_p1a, bh, facecolor=C_P1, edgecolor='#222222', lw=0.8))
ax_a.add_patch(patches.Rectangle((bx + w_p1a, by), w_p2, bh, facecolor=C_P2, edgecolor='#222222', lw=0.8))
ax_a.add_patch(patches.Rectangle((bx + w_p1a + w_p2, by), w_p1b, bh, facecolor=C_P1, edgecolor='#222222', lw=0.8))

ax_a.text(bx + w_p1a/2, by + bh/2, '$P_1$ (Home)', color='white', ha='center', va='center', fontsize=7.5, fontweight='bold')
ax_a.text(bx + w_p1a + w_p2/2, by + bh/2, '$P_2$ Cassette (1,500 nt)', color='white', ha='center', va='center', fontsize=7.5, fontweight='bold')
ax_a.text(bx + w_p1a + w_p2 + w_p1b/2, by + bh/2, '$P_1$ (Home)', color='white', ha='center', va='center', fontsize=7.5, fontweight='bold')

ax_a.text(bx, by - 0.038, '1 nt', ha='left', fontsize=7.2, color='#555555')
ax_a.text(bx + w_p1a, by - 0.038, '1,500', ha='center', fontsize=7.2, color='#555555')
ax_a.text(bx + w_p1a + w_p2, by - 0.038, '3,000', ha='center', fontsize=7.2, color='#555555')
ax_a.text(bx + bw, by - 0.038, '5,000 nt', ha='right', fontsize=7.2, color='#555555')


# ==============================================================================
# Panel B: Scenario 1.2 - Multi-Way Mosaic Swarms
# ==============================================================================
ax_b = axes[0, 1]
ax_b.text(0.02, 0.98, 'B. Multi-Way Mosaic Swarms (Scenario 1.2)', fontsize=11, fontweight='bold', va='top')

bw_b = 0.90
bh_b = 0.070
bx_b = 0.04

# 1. Dual Cassette
y1 = 0.74
ax_b.text(bx_b, y1 + 0.10, '1. Dual Cassette ($K_{\\mathrm{donor}}=2$ donors, 2 separate non-overlapping tracts)', fontsize=8.0, fontweight='bold', color='#222222')
w_seg = bw_b / 5
ax_b.add_patch(patches.Rectangle((bx_b + 0*w_seg, y1), w_seg, bh_b, facecolor=C_P1, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + 1*w_seg, y1), w_seg, bh_b, facecolor=C_P2, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + 2*w_seg, y1), w_seg, bh_b, facecolor=C_P1, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + 3*w_seg, y1), w_seg, bh_b, facecolor=C_P3, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + 4*w_seg, y1), w_seg, bh_b, facecolor=C_P1, edgecolor='#222222', lw=0.8))

ax_b.text(bx_b + 0.5*w_seg, y1 + bh_b/2, '$P_1$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + 1.5*w_seg, y1 + bh_b/2, '$P_2$ Donor', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + 2.5*w_seg, y1 + bh_b/2, '$P_1$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + 3.5*w_seg, y1 + bh_b/2, '$P_3$ Donor', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + 4.5*w_seg, y1 + bh_b/2, '$P_1$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b, y1 - 0.038, 'RhizAeon: 2.00 / 2.00 donors resolved (100%) | 3SEQ: 1.02 donors (misses $P_3$)', fontsize=7.2, color='#003366', fontweight='bold')

# 2. Nested Russian Doll
y2 = 0.44
ax_b.text(bx_b, y2 + 0.10, '2. Nested Russian Doll ($K_{\\mathrm{donor}}=2$ donors, hierarchical embedded introgression)', fontsize=8.0, fontweight='bold', color='#222222')
w_p1_l = bw_b * (1000 / 6000)
w_p2_tot = bw_b * (3000 / 6000)
w_p1_r = bw_b * (2000 / 6000)
w_p2_sub1 = bw_b * (1000 / 6000)
w_p3_sub2 = bw_b * (1000 / 6000)
w_p2_sub3 = bw_b * (1000 / 6000)

ax_b.add_patch(patches.Rectangle((bx_b, y2), w_p1_l, bh_b, facecolor=C_P1, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + w_p1_l, y2), w_p2_sub1, bh_b, facecolor=C_P2, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + w_p1_l + w_p2_sub1, y2), w_p3_sub2, bh_b, facecolor=C_P3, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + w_p1_l + w_p2_sub1 + w_p3_sub2, y2), w_p2_sub3, bh_b, facecolor=C_P2, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + w_p1_l + w_p2_tot, y2), w_p1_r, bh_b, facecolor=C_P1, edgecolor='#222222', lw=0.8))

ax_b.text(bx_b + w_p1_l/2, y2 + bh_b/2, '$P_1$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + w_p1_l + w_p2_sub1/2, y2 + bh_b/2, '$P_2$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + w_p1_l + w_p2_sub1 + w_p3_sub2/2, y2 + bh_b/2, '$P_3$ (inner)', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + w_p1_l + w_p2_sub1 + w_p3_sub2 + w_p2_sub3/2, y2 + bh_b/2, '$P_2$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + w_p1_l + w_p2_tot + w_p1_r/2, y2 + bh_b/2, '$P_1$', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b, y2 - 0.038, 'RhizAeon: 1.10 donors resolved | 3SEQ: 1.00 donor (completely blind to inner $P_3$)', fontsize=7.2, color='#003366', fontweight='bold')

# 3. Triple Crossover Cascade
y3 = 0.14
ax_b.text(bx_b, y3 + 0.10, '3. Triple Crossover Cascade ($K_{\\mathrm{donor}}=2$ transitions, serial chromosome mosaic)', fontsize=8.0, fontweight='bold', color='#222222')
w3 = bw_b / 3
ax_b.add_patch(patches.Rectangle((bx_b + 0*w3, y3), w3, bh_b, facecolor=C_P1, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + 1*w3, y3), w3, bh_b, facecolor=C_P2, edgecolor='#222222', lw=0.8))
ax_b.add_patch(patches.Rectangle((bx_b + 2*w3, y3), w3, bh_b, facecolor=C_P3, edgecolor='#222222', lw=0.8))

ax_b.text(bx_b + 0.5*w3, y3 + bh_b/2, '$P_1$ (1--2,000 nt)', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + 1.5*w3, y3 + bh_b/2, '$P_2$ (2,001--4,000 nt)', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b + 2.5*w3, y3 + bh_b/2, '$P_3$ (4,001--6,000 nt)', color='white', ha='center', va='center', fontsize=7.4, fontweight='bold')
ax_b.text(bx_b, y3 - 0.038, 'RhizAeon: 1.63 donors resolved, 100% recall | 3SEQ: 95.3% background parent false alarms', fontsize=7.2, color='#003366', fontweight='bold')


# ==============================================================================
# Panel C: Scenario 1.3 - Ghost Donor Introgression
# ==============================================================================
ax_c = axes[1, 0]
ax_c.text(0.02, 0.98, 'C. Ghost Donor Introgression (Scenario 1.3)', fontsize=11, fontweight='bold', va='top')

col_x = [0.03, 0.355, 0.68]
col_w = 0.29

# Sub-regime 1: Surrogate Sister
cx1 = col_x[0]
ax_c.text(cx1 + col_w/2, 0.90, '(i) Surrogate Sister', fontsize=8.0, fontweight='bold', ha='center', color='#222222')

ax_c.plot([cx1 + 0.02, cx1 + 0.10], [0.72, 0.84], color=C_P1, lw=2.0)
ax_c.text(cx1 + 0.11, 0.84, '$P_1$ (Home)', fontsize=7.2, va='center', color=C_P1, fontweight='bold')

ax_c.plot([cx1 + 0.02, cx1 + 0.10], [0.72, 0.58], color='#d95f02', lw=2.0)
ax_c.text(cx1 + 0.11, 0.58, '$P_2$ (Innocent)', fontsize=7.2, va='center', color='#d95f02')

ax_c.plot([cx1 + 0.02, cx1 + 0.09], [0.72, 0.72], color='#555555', lw=1.5)
ax_c.plot([cx1 + 0.09, cx1 + 0.17], [0.72, 0.77], color=C_P3, lw=2.0)
ax_c.text(cx1 + 0.18, 0.77, "$P_2'$ (Sister)", fontsize=7.2, va='center', color=C_P3, fontweight='bold')

ax_c.plot([cx1 + 0.09, cx1 + 0.17], [0.72, 0.67], color=C_GHOST, lw=1.8, linestyle='--')
ax_c.text(cx1 + 0.18, 0.67, "$D_{\\mathrm{ghost}}$", fontsize=7.2, va='center', color=C_GHOST)
ax_c.scatter([cx1 + 0.17], [0.67], color='#d95f02', marker='x', s=30, lw=1.8, zorder=5)

# Transfer arrow into R
ax_c.annotate('', xy=(cx1 + 0.10, 0.82), xytext=(cx1 + 0.16, 0.68),
              arrowprops=dict(arrowstyle="->", color='#d95f02', lw=1.4, linestyle=':', connectionstyle="arc3,rad=-0.2"))

ax_c.text(cx1 + col_w/2, 0.48,
          'True donor unsampled.\nSister clade $P_2\'$ present.\n• Goal: Attribute to $P_2\'$\n• Zero false calls on $P_2$',
          ha='center', va='top', fontsize=7.0, color='#333333',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#fafafa', edgecolor='#cccccc', lw=0.6))

# Sub-regime 2: Pure Ghost
cx2 = col_x[1]
ax_c.text(cx2 + col_w/2, 0.90, '(ii) Pure Ghost Lineage', fontsize=8.0, fontweight='bold', ha='center', color='#222222')

ax_c.plot([cx2 + 0.02, cx2 + 0.10], [0.72, 0.84], color=C_P1, lw=2.0)
ax_c.text(cx2 + 0.11, 0.84, '$P_1$ (Home)', fontsize=7.2, va='center', color=C_P1, fontweight='bold')

ax_c.plot([cx2 + 0.02, cx2 + 0.10], [0.72, 0.58], color='#d95f02', lw=2.0)
ax_c.text(cx2 + 0.11, 0.58, '$P_2$ (Innocent)', fontsize=7.2, va='center', color='#d95f02')

ax_c.plot([cx2 + 0.02, cx2 + 0.17], [0.72, 0.71], color=C_GHOST, lw=1.8, linestyle='--')
ax_c.text(cx2 + 0.18, 0.71, "$D_{\\mathrm{ghost}}$", fontsize=7.2, va='center', color=C_GHOST)
ax_c.scatter([cx2 + 0.17], [0.71], color='#d95f02', marker='x', s=30, lw=1.8, zorder=5)

# Transfer arrow into R
ax_c.annotate('', xy=(cx2 + 0.10, 0.82), xytext=(cx2 + 0.16, 0.72),
              arrowprops=dict(arrowstyle="->", color='#d95f02', lw=1.4, linestyle=':', connectionstyle="arc3,rad=-0.2"))

ax_c.text(cx2 + col_w/2, 0.48,
          'Deep ancestral lineage.\nNo donor relatives sampled.\n• Goal: Open-world bounds\n• Zero false calls on $P_2$',
          ha='center', va='top', fontsize=7.0, color='#333333',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#fafafa', edgecolor='#cccccc', lw=0.6))

# Sub-regime 3: Deep Reservoir
cx3 = col_x[2]
ax_c.text(cx3 + col_w/2, 0.90, '(iii) Deep Reservoir', fontsize=8.0, fontweight='bold', ha='center', color='#222222')

ax_c.plot([cx3 + 0.02, cx3 + 0.10], [0.72, 0.84], color=C_P1, lw=2.0)
ax_c.text(cx3 + 0.11, 0.84, '$P_1$ (Home)', fontsize=7.2, va='center', color=C_P1, fontweight='bold')

ax_c.plot([cx3 + 0.02, cx3 + 0.10], [0.72, 0.58], color='#d95f02', lw=2.0)
ax_c.text(cx3 + 0.11, 0.58, '$P_2$ (Innocent)', fontsize=7.2, va='center', color='#d95f02')

ax_c.plot([cx3 + 0.01, cx3 + 0.17], [0.72, 0.69], color='#7570b3', lw=1.8, linestyle='--')
ax_c.text(cx3 + 0.18, 0.69, "Reservoir ($0.15$)", fontsize=7.0, va='center', color='#7570b3')
ax_c.scatter([cx3 + 0.17], [0.69], color='#d95f02', marker='x', s=30, lw=1.8, zorder=5)

# Transfer arrow into R
ax_c.annotate('', xy=(cx3 + 0.10, 0.82), xytext=(cx3 + 0.16, 0.70),
              arrowprops=dict(arrowstyle="->", color='#d95f02', lw=1.4, linestyle=':', connectionstyle="arc3,rad=-0.2"))

ax_c.text(cx3 + col_w/2, 0.48,
          'Cross-species reservoir.\nHigh divergence ($d=0.15$).\n• Goal: Open-world bounds\n• Zero false calls on $P_2$',
          ha='center', va='top', fontsize=7.0, color='#333333',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#fafafa', edgecolor='#cccccc', lw=0.6))

# Bottom highlight banner across Panel C
ax_c.text(0.50, 0.09,
          '• Triplet Trap (3SEQ): 100% false framing of innocent $P_2$\n'
          '• RhizAeon Open-World: 98.3% sister recall, 0.0% false framing',
          ha='center', va='center', fontsize=7.2, fontweight='bold', color='#002b66',
          bbox=dict(boxstyle='round,pad=0.38', facecolor='#eef7fc', edgecolor='#88ccee', lw=0.8))


# ==============================================================================
# Panel D: Scenario 1.4 - Bacterial Micro-Conversions
# ==============================================================================
ax_d = axes[1, 1]
ax_d.text(0.02, 0.98, 'D. Bacterial Micro-Conversions (Scenario 1.4)', fontsize=11, fontweight='bold', va='top')

# Column coordinates
col_mosaic_x = 0.04
col_mosaic_w = 0.32
col_tl_x = 0.42
col_mut_x = 0.58
col_pow_x = 0.76

# Table Header
ax_d.text(col_mosaic_x, 0.89, 'Mosaic Structure', fontsize=7.6, fontweight='bold', color='#444444')
ax_d.text(col_tl_x, 0.89, 'Tract Length', fontsize=7.6, fontweight='bold', color='#444444')
ax_d.text(col_mut_x, 0.89, 'Mutations', fontsize=7.6, fontweight='bold', color='#444444')
ax_d.text(col_pow_x, 0.89, 'Power (Rhiz / 3SEQ)', fontsize=7.6, fontweight='bold', color='#444444')

tract_configs = [
    (800, '800 nt', '16--74 muts', '100% / 100%', 0.77),
    (400, '400 nt', '7--38 muts',  '98% / 100%',  0.65),
    (200, '200 nt', '4--20 muts',  '91% / 96%',   0.53),
    (100, '100 nt', '2--9 muts',   '89% / 55%',   0.41),
    (50,  '50 nt',  '1--5 muts',   '40% / 5%',    0.29),
]

dh = 0.052

for tl, tl_label, mut_label, power_label, dy in tract_configs:
    x_start = col_mosaic_x + col_mosaic_w * (2100 / 5000)
    w_tract = max(0.007, col_mosaic_w * (tl / 5000))
    
    # Left P1 flank
    ax_d.add_patch(patches.Rectangle((col_mosaic_x, dy), x_start - col_mosaic_x, dh, facecolor=C_P1, edgecolor='#444444', lw=0.6))
    # Inserted micro-tract P2
    ax_d.add_patch(patches.Rectangle((x_start, dy), w_tract, dh, facecolor=C_P2, edgecolor='#222222', lw=0.8))
    # Right P1 flank
    ax_d.add_patch(patches.Rectangle((x_start + w_tract, dy), (col_mosaic_x + col_mosaic_w) - (x_start + w_tract), dh, facecolor=C_P1, edgecolor='#444444', lw=0.6))
    
    # Table entries
    ax_d.text(col_tl_x, dy + dh/2, f'{tl_label}', fontsize=7.8, fontweight='bold', va='center', color='#222222')
    ax_d.text(col_mut_x, dy + dh/2, f'{mut_label}', fontsize=7.4, va='center', color='#555555')
    ax_d.text(col_pow_x, dy + dh/2, f'{power_label}', fontsize=7.6, va='center', color='#002b66', fontweight='bold')

# Chromosome limits
ax_d.text(col_mosaic_x, 0.23, '1 nt', fontsize=7.0, color='#666666', ha='left')
ax_d.text(col_mosaic_x + col_mosaic_w, 0.23, '5,000 nt', fontsize=7.0, color='#666666', ha='right')

# Scale-Space Matched Filter illustration at bottom
ax_d.text(0.04, 0.16, 'Continuous Dyadic Scale-Space Matched Filter ($\\psi_\\tau$, $\\tau \\in \\{15, 25, 50, 100, 200\\}$ codons):', fontsize=7.4, fontweight='bold', color='#222222')
wave_x = np.linspace(col_mosaic_x, col_mosaic_x + col_mosaic_w, 200)
center_x = col_mosaic_x + col_mosaic_w * (2150 / 5000)
sigma = 0.030
wave_y = 0.07 + 0.045 * (1 - ((wave_x - center_x)/sigma)**2) * np.exp(-0.5 * ((wave_x - center_x)/sigma)**2)
ax_d.plot(wave_x, wave_y, color='#0288d1', lw=1.8)
ax_d.axhline(0.07, color='#bbbbbb', linestyle=':', lw=0.8)
ax_d.text(0.04, 0.025, 'Wavelet localization achieves median breakpoint error of 13--19 nt', fontsize=7.2, color='#0288d1', ha='left', fontweight='bold')

plt.tight_layout()

out_pdf = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_scenarios_schematic.pdf'
out_png = '/Users/sergei/Projects/rhizaeon-code/paper/figures/fig_scenarios_schematic.png'
plt.savefig(out_pdf, bbox_inches='tight')
plt.savefig(out_png, bbox_inches='tight')
print(f'Successfully updated {out_pdf} and {out_png}')
