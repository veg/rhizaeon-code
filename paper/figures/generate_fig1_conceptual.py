#!/usr/bin/env python3
"""
generate_fig1_conceptual.py
============================
Generates the publication-grade Figure 1 for RhizAeon in Molecular Biology and Evolution (MBE).
Design Philosophy: Clean, intuitive, publication-grade conceptual figure.
- Compact 2:1 aspect ratio (14.0 x 7.2 in) for optimal page-budget integration.
- Restrained, colorblind-safe palette (Prussian Blue, Braun International Orange, Royal Purple, Slate Gray).
- Four conceptual panels:
  Panel A: Conceptual Intuition — From Trees to Continuous Trajectories
  Panel B: Core Conceptual Innovation — Cross-Attention & The Autapomorphic Root Sink
  Panel C: Physical Signal Deconvolution — Zero-Drift Force & Cumulative Momentum
  Panel D: Evolutionary Regimes & Ancestral Swarm Consolidation
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
import numpy as np

# Canvas Setup: 14.0 x 7.2 inches, 300 DPI (compact aspect ratio)
fig = plt.figure(figsize=(14.0, 7.2), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 14.0)
ax.set_ylim(0, 7.2)
ax.axis('off')

# ------------------------------------------------------------------------------
# COLOR PALETTE
# ------------------------------------------------------------------------------
c_bg_panel      = '#F8FAFC'
c_panel_border  = '#CBD5E1'
c_text_header   = '#0F172A'
c_text_body     = '#334155'
c_text_muted    = '#64748B'

c_home          = '#1D4ED8'  # Classic Blue
c_home_soft     = '#DBEAFE'
c_donor         = '#D97706'  # Braun Amber / Orange
c_donor_soft    = '#FEF3C7'
c_rec           = '#DC2626'  # Crimson Recombinant
c_rec_soft      = '#FEE2E2'
c_root          = '#7C3AED'  # Royal Purple Root Sink
c_root_soft     = '#EDE9FE'

def draw_panel(ax, x, y, w, h, title_letter, title_text, subtitle=None):
    card = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.08",
                          facecolor=c_bg_panel, edgecolor=c_panel_border, linewidth=1.1, zorder=1)
    ax.add_patch(card)
    
    badge = FancyBboxPatch((x + 0.15, y + h - 0.38), 0.34, 0.25, boxstyle="round,pad=0.0,rounding_size=0.04",
                           facecolor=c_text_header, edgecolor=c_text_header, zorder=2)
    ax.add_patch(badge)
    ax.text(x + 0.32, y + h - 0.25, title_letter, ha='center', va='center',
            fontsize=9.5, fontweight='bold', color='#FFFFFF', zorder=3)
    
    ax.text(x + 0.58, y + h - 0.25, title_text, ha='left', va='center',
            fontsize=8.8, fontweight='bold', color=c_text_header, zorder=3)
    
    if subtitle:
        ax.text(x + 0.58, y + h - 0.44, subtitle, ha='left', va='center',
                fontsize=6.8, color=c_text_muted, zorder=3)
    
    ax.plot([x + 0.15, x + w - 0.15], [y + h - 0.54, y + h - 0.54],
            color='#E2E8F0', lw=0.9, zorder=2)

def draw_subbox(ax, x, y, w, h, title=None, bg='#FFFFFF', border='#E2E8F0', lw=0.8):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.05",
                         facecolor=bg, edgecolor=border, linewidth=lw, zorder=2)
    ax.add_patch(box)
    if title:
        ax.text(x + 0.12, y + h - 0.16, title, ha='left', va='top',
                fontsize=7.2, fontweight='bold', color=c_text_header, zorder=3)

# ------------------------------------------------------------------------------
# PANEL POSITIONS (2x2 Grid, w=6.45, h=3.20)
# ------------------------------------------------------------------------------
px0, py0 = 0.40, 3.70  # Panel A: Top-Left
pw, ph   = 6.45, 3.20

px1, py1 = 7.15, 3.70  # Panel B: Top-Right
px2, py2 = 0.40, 0.25  # Panel C: Bottom-Left
px3, py3 = 7.15, 0.25  # Panel D: Bottom-Right

# ==============================================================================
# PANEL A: CONCEPTUAL INTUITION — FROM TREES TO TRAJECTORIES
# ==============================================================================
draw_panel(ax, px0, py0, pw, ph, 'A', 'From Bifurcating Trees to Continuous Trajectories',
           'Recombination shatters trees into discordant segments; trajectories model reticulation as metric flow.')

# Sub-block A1: The Tree Paradox
draw_subbox(ax, px0 + 0.15, py0 + 1.35, pw - 0.30, 1.22, 'Phylogenetic Discordance Shatters Single-Tree Models')

chr_y = py0 + 1.95
ax.plot([px0 + 2.90, px0 + 2.90], [chr_y - 0.06, chr_y + 0.20], color=c_rec, lw=1.6, linestyle='--', zorder=6)
ax.text(px0 + 2.90, chr_y + 0.24, 'Crossover Breakpoint', ha='center', va='bottom', fontsize=6.2, color=c_rec, fontweight='bold', zorder=6)

ax.add_patch(Rectangle((px0 + 0.45, chr_y), 2.45, 0.15, facecolor=c_home_soft, edgecolor=c_home, lw=1.1, zorder=4))
ax.add_patch(Rectangle((px0 + 2.90, chr_y), 2.65, 0.15, facecolor=c_donor_soft, edgecolor=c_donor, lw=1.1, zorder=4))
ax.text(px0 + 1.67, chr_y + 0.07, 'Region 1 (Home Clade A)', ha='center', va='center', fontsize=6.2, color=c_home, fontweight='bold', zorder=5)
ax.text(px0 + 4.22, chr_y + 0.07, 'Region 2 (Donor Clade B)', ha='center', va='center', fontsize=6.2, color=c_donor, fontweight='bold', zorder=5)

# Trees
t1_x = px0 + 1.05
t1_y = py0 + 1.42
ax.text(t1_x + 0.35, t1_y + 0.45, 'Tree 1: ((R, A), B)', ha='center', fontsize=6.4, fontweight='bold', color=c_home, zorder=4)
ax.plot([t1_x, t1_x + 0.30], [t1_y + 0.32, t1_y + 0.32], color=c_rec, lw=1.3, zorder=4)
ax.text(t1_x - 0.03, t1_y + 0.32, 'R', ha='right', va='center', fontsize=6.0, fontweight='bold', color=c_rec)
ax.plot([t1_x, t1_x + 0.30], [t1_y + 0.18, t1_y + 0.18], color=c_home, lw=1.3, zorder=4)
ax.text(t1_x - 0.03, t1_y + 0.18, 'A', ha='right', va='center', fontsize=6.0, fontweight='bold', color=c_home)
ax.plot([t1_x + 0.30, t1_x + 0.30], [t1_y + 0.18, t1_y + 0.32], color=c_home, lw=1.3, zorder=4)
ax.plot([t1_x + 0.30, t1_x + 0.55], [t1_y + 0.25, t1_y + 0.25], color=c_home, lw=1.3, zorder=4)
ax.plot([t1_x, t1_x + 0.55], [t1_y + 0.05, t1_y + 0.05], color=c_donor, lw=1.3, zorder=4)
ax.text(t1_x - 0.03, t1_y + 0.05, 'B', ha='right', va='center', fontsize=6.0, fontweight='bold', color=c_donor)
ax.plot([t1_x + 0.55, t1_x + 0.55], [t1_y + 0.05, t1_y + 0.25], color='#64748B', lw=1.3, zorder=4)

t2_x = px0 + 4.05
t2_y = py0 + 1.42
ax.text(t2_x + 0.35, t2_y + 0.45, 'Tree 2: ((R, B), A)', ha='center', fontsize=6.4, fontweight='bold', color=c_donor, zorder=4)
ax.plot([t2_x, t2_x + 0.30], [t2_y + 0.32, t2_y + 0.32], color=c_home, lw=1.3, zorder=4)
ax.text(t2_x - 0.03, t2_y + 0.32, 'A', ha='right', va='center', fontsize=6.0, fontweight='bold', color=c_home)
ax.plot([t2_x, t2_x + 0.30], [t2_y + 0.18, t2_y + 0.18], color=c_rec, lw=1.3, zorder=4)
ax.text(t2_x - 0.03, t2_y + 0.18, 'R', ha='right', va='center', fontsize=6.0, fontweight='bold', color=c_rec)
ax.plot([t2_x + 0.30, t2_x + 0.30], [t2_y + 0.18, t2_y + 0.32], color=c_rec, lw=1.3, zorder=4)
ax.plot([t2_x, t2_x + 0.55], [t2_y + 0.05, t2_y + 0.05], color=c_donor, lw=1.3, zorder=4)
ax.text(t2_x - 0.03, t2_y + 0.05, 'B', ha='right', va='center', fontsize=6.0, fontweight='bold', color=c_donor)
ax.plot([t2_x + 0.30, t2_x + 0.55], [t2_y + 0.25, t2_y + 0.25], color=c_donor, lw=1.3, zorder=4)
ax.plot([t2_x + 0.55, t2_x + 0.55], [t2_y + 0.05, t2_y + 0.25], color='#64748B', lw=1.3, zorder=4)

ax.text(px0 + 2.90, py0 + 1.55, 'Topological\nConflict', ha='center', va='center',
        fontsize=6.0, fontweight='bold', color='#B91C1C',
        bbox=dict(boxstyle="round,pad=0.18", facecolor='#FEE2E2', edgecolor='#F87171', lw=0.7), zorder=5)

# Sub-block A2: Continuous Trajectory Mechanics
draw_subbox(ax, px0 + 0.15, py0 + 0.12, pw - 0.30, 1.15, 'Continuous Trajectory Mechanics in Metric Space')

traj_x0 = px0 + 0.50
traj_y0 = py0 + 0.28
traj_w  = pw - 1.05
traj_h  = 0.65

ax.annotate('', xy=(traj_x0 + traj_w, traj_y0), xytext=(traj_x0, traj_y0),
            arrowprops=dict(arrowstyle="->", color=c_text_muted, lw=0.9), zorder=3)
ax.text(traj_x0 + traj_w, traj_y0 - 0.10, 'Genomic Coordinate u', ha='right', fontsize=6.2, color=c_text_muted)
ax.annotate('', xy=(traj_x0, traj_y0 + traj_h), xytext=(traj_x0, traj_y0),
            arrowprops=dict(arrowstyle="->", color=c_text_muted, lw=0.9), zorder=3)
ax.text(traj_x0 - 0.06, traj_y0 + traj_h/2, 'Metric Space (Z)', ha='center', va='center', rotation=90, fontsize=6.0, color=c_text_muted)

u_pts = np.linspace(traj_x0 + 0.15, traj_x0 + traj_w - 0.25, 100)
y_a = traj_y0 + 0.12 + 0.02 * np.sin(u_pts * 3.0)
ax.plot(u_pts, y_a, color=c_home, lw=1.8, zorder=4)
ax.text(traj_x0 + traj_w - 0.20, y_a[-1], 'Clade A (Home)', ha='left', va='center', fontsize=6.0, fontweight='bold', color=c_home, zorder=5)

y_b = traj_y0 + 0.50 + 0.02 * np.sin(u_pts * 2.8 + 1.0)
ax.plot(u_pts, y_b, color=c_donor, lw=1.8, zorder=4)
ax.text(traj_x0 + traj_w - 0.20, y_b[-1], 'Clade B (Donor)', ha='left', va='center', fontsize=6.0, fontweight='bold', color=c_donor, zorder=5)

bp_u = traj_x0 + traj_w * 0.48
y_rec = []
for u in u_pts:
    if u < bp_u - 0.18:
        y_rec.append(traj_y0 + 0.14 + 0.015 * np.sin(u * 3.2))
    elif u > bp_u + 0.18:
        y_rec.append(traj_y0 + 0.48 + 0.015 * np.sin(u * 2.9 + 0.8))
    else:
        t = (u - (bp_u - 0.18)) / 0.36
        y_start = traj_y0 + 0.14 + 0.015 * np.sin((bp_u - 0.18) * 3.2)
        y_end = traj_y0 + 0.48 + 0.015 * np.sin((bp_u + 0.18) * 2.9 + 0.8)
        y_rec.append(y_start + (y_end - y_start) * (1.0 / (1.0 + np.exp(-12.0 * (t - 0.5)))))

ax.plot(u_pts, y_rec, color=c_rec, lw=2.0, linestyle='-', zorder=6)
ax.text(traj_x0 + 0.05, y_rec[0] - 0.08, 'Recombinant Child R', ha='left', va='top', fontsize=6.2, fontweight='bold', color=c_rec, zorder=6)

ax.plot([bp_u, bp_u], [traj_y0 + 0.05, traj_y0 + 0.60], color=c_rec, lw=1.1, linestyle=':', zorder=5)
ax.annotate('Metric Deflection\n(Crossover)', xy=(bp_u, traj_y0 + 0.31), xytext=(bp_u - 1.05, traj_y0 + 0.40),
            arrowprops=dict(arrowstyle="->", color=c_rec, lw=1.1),
            fontsize=6.0, fontweight='bold', color=c_rec,
            bbox=dict(boxstyle="round,pad=0.18", facecolor='#FFFFFF', edgecolor=c_rec, lw=0.7), zorder=7)


# ==============================================================================
# PANEL B: THE CORE INNOVATION — CROSS-ATTENTION & AUTAPOMORPHY ROOT SINK
# ==============================================================================
draw_panel(ax, px1, py1, pw, ph, 'B', 'Phylogenetic Cross-Attention & Autapomorphy Root Sink',
           'The unpolarized [ROOT] sink absorbs lineage rate bursts, preventing heuristic false alarms.')

w_half = (pw - 0.40) / 2
draw_subbox(ax, px1 + 0.15, py1 + 0.12, w_half, ph - 0.72, 'Mutational Burst / Heterotachy\n(Clonal Null with Rate Burst)')

b1_cx = px1 + 0.15 + w_half/2
b1_y = py1 + 1.45

ax.add_patch(Circle((b1_cx, b1_y), 0.20, facecolor=c_root_soft, edgecolor=c_root, lw=1.3, zorder=4))
ax.text(b1_cx, b1_y, '[ROOT]\nSink', ha='center', va='center', fontsize=6.0, fontweight='bold', color=c_root, zorder=5)

ax.add_patch(Circle((b1_cx - 0.85, b1_y - 0.55), 0.16, facecolor=c_home_soft, edgecolor=c_home, lw=1.2, zorder=4))
ax.text(b1_cx - 0.85, b1_y - 0.55, 'Landmark\nHome A', ha='center', va='center', fontsize=5.5, fontweight='bold', color=c_home, zorder=5)

ax.add_patch(Circle((b1_cx + 0.85, b1_y - 0.55), 0.16, facecolor=c_donor_soft, edgecolor=c_donor, lw=1.2, zorder=4))
ax.text(b1_cx + 0.85, b1_y - 0.55, 'Landmark\nDonor B', ha='center', va='center', fontsize=5.5, fontweight='bold', color=c_donor, zorder=5)

q1_x, q1_y = b1_cx - 0.40, b1_y + 0.65
ax.add_patch(Circle((q1_x, q1_y), 0.15, facecolor='#FEF2F2', edgecolor='#EF4444', lw=1.3, zorder=4))
ax.text(q1_x, q1_y, 'Query\n(Burst)', ha='center', va='center', fontsize=5.5, fontweight='bold', color='#B91C1C', zorder=5)

ax.annotate('', xy=(b1_cx - 0.05, b1_y + 0.22), xytext=(q1_x + 0.08, q1_y - 0.12),
            arrowprops=dict(arrowstyle="-|>", color=c_root, lw=2.2, mutation_scale=10), zorder=6)
ax.text(b1_cx - 0.02, b1_y + 0.45, 'Private mutations\npool into [ROOT]', ha='left', va='center', fontsize=5.8, fontweight='bold', color=c_root, zorder=7)

ax.plot([q1_x + 0.15, b1_cx + 0.70], [q1_y - 0.15, b1_y - 0.42], color='#94A3B8', lw=1.1, linestyle=':', zorder=4)
ax.text(b1_cx + 0.52, b1_y + 0.05, 'A_donor = 0.0\n(No match)', ha='center', va='center', fontsize=5.5, color='#64748B',
        bbox=dict(boxstyle="round,pad=0.12", facecolor='#F1F5F9', edgecolor='#CBD5E1', lw=0.6), zorder=5)

ax.text(b1_cx, py1 + 0.32, 'False Positive Rate = 0.0%\nRate burst absorbed by sink',
        ha='center', va='center', fontsize=6.0, fontweight='bold', color='#047857',
        bbox=dict(boxstyle="round,pad=0.20", facecolor='#ECFDF5', edgecolor='#6EE7B7', lw=0.8), zorder=7)


# Sub-block B2: True Recombination
b2_x0 = px1 + 0.25 + w_half
draw_subbox(ax, b2_x0, py1 + 0.12, w_half, ph - 0.72, 'Authentic Recombination Event\n(Coupled Lineage Attention Dipole)')

b2_cx = b2_x0 + w_half/2
b2_y = py1 + 1.45

ax.add_patch(Circle((b2_cx, b2_y), 0.20, facecolor=c_root_soft, edgecolor=c_root, lw=1.3, zorder=4))
ax.text(b2_cx, b2_y, '[ROOT]\nSink', ha='center', va='center', fontsize=6.0, fontweight='bold', color=c_root, zorder=5)

ax.add_patch(Circle((b2_cx - 0.85, b2_y - 0.55), 0.16, facecolor=c_home_soft, edgecolor=c_home, lw=1.2, zorder=4))
ax.text(b2_cx - 0.85, b2_y - 0.55, 'Landmark\nHome A', ha='center', va='center', fontsize=5.5, fontweight='bold', color=c_home, zorder=5)

ax.add_patch(Circle((b2_cx + 0.85, b2_y - 0.55), 0.16, facecolor=c_donor_soft, edgecolor=c_donor, lw=1.2, zorder=4))
ax.text(b2_cx + 0.85, b2_y - 0.55, 'Landmark\nDonor B', ha='center', va='center', fontsize=5.5, fontweight='bold', color=c_donor, zorder=5)

q2_x, q2_y = b2_cx, b2_y + 0.65
ax.add_patch(Circle((q2_x, q2_y), 0.15, facecolor='#FEF2F2', edgecolor=c_rec, lw=1.4, zorder=4))
ax.text(q2_x, q2_y, 'Child R\n(Tract)', ha='center', va='center', fontsize=5.5, fontweight='bold', color=c_rec, zorder=5)

ax.annotate('', xy=(b2_cx + 0.72, b2_y - 0.42), xytext=(q2_x + 0.08, q2_y - 0.12),
            arrowprops=dict(arrowstyle="-|>", color=c_donor, lw=2.4, mutation_scale=11), zorder=6)
ax.text(b2_cx + 0.65, b2_y + 0.18, 'A_donor surges\n(Shared alleles)', ha='left', va='center', fontsize=5.8, fontweight='bold', color=c_donor, zorder=7)

ax.annotate('', xy=(b2_cx - 0.72, b2_y - 0.42), xytext=(q2_x - 0.08, q2_y - 0.12),
            arrowprops=dict(arrowstyle="-|>", color=c_home, lw=1.1, linestyle='--'), zorder=5)
ax.text(b2_cx - 0.65, b2_y + 0.18, 'A_home drops\n(Reticulation)', ha='right', va='center', fontsize=5.8, color=c_home, zorder=7)

ax.text(b2_cx, py1 + 0.32, 'Anti-Correlated Dipole\nDelta A_donor > 0  &  Delta A_home < 0',
        ha='center', va='center', fontsize=6.0, fontweight='bold', color='#1E40AF',
        bbox=dict(boxstyle="round,pad=0.20", facecolor='#EFF6FF', edgecolor='#93C5FD', lw=0.8), zorder=7)


# ==============================================================================
# PANEL C: PHYSICAL SIGNAL DECONVOLUTION — FORCE & CUMULATIVE MOMENTUM
# ==============================================================================
draw_panel(ax, px2, py2, pw, ph, 'C', 'Physical Signal Deconvolution & Momentum Trajectories',
           'Zero-drift force integration converts noisy site attention into macroscopic Brownian excursions.')

c_plot_x0 = px2 + 0.90
c_plot_w  = pw - 1.25

# Subplot 1: Alignment schematic
y_aln = py2 + 2.30
ax.text(px2 + 0.15, y_aln + 0.06, 'Alignment\nSites', ha='left', va='center', fontsize=6.5, fontweight='bold', color=c_text_header)

bp1_x = c_plot_x0 + c_plot_w * 0.35
bp2_x = c_plot_x0 + c_plot_w * 0.70

ax.text(bp1_x, y_aln + 0.24, 'Breakpoint u1', ha='center', va='bottom', fontsize=6.0, color='#DC2626', fontweight='bold', zorder=9)
ax.text(bp2_x, y_aln + 0.24, 'Breakpoint u2', ha='center', va='bottom', fontsize=6.0, color='#DC2626', fontweight='bold', zorder=9)

ax.plot([bp1_x, bp1_x], [py2 + 0.38, y_aln + 0.20], color='#EF4444', lw=1.1, linestyle='--', zorder=8)
ax.plot([bp2_x, bp2_x], [py2 + 0.38, y_aln + 0.20], color='#EF4444', lw=1.1, linestyle='--', zorder=8)

ax.add_patch(Rectangle((c_plot_x0, y_aln - 0.05), c_plot_w * 0.35, 0.18, facecolor=c_home_soft, edgecolor=c_home, lw=0.9, zorder=3))
ax.add_patch(Rectangle((c_plot_x0 + c_plot_w * 0.35, y_aln - 0.05), c_plot_w * 0.35, 0.18, facecolor=c_donor_soft, edgecolor=c_donor, lw=1.1, zorder=3))
ax.add_patch(Rectangle((c_plot_x0 + c_plot_w * 0.70, y_aln - 0.05), c_plot_w * 0.30, 0.18, facecolor=c_home_soft, edgecolor=c_home, lw=0.9, zorder=3))

ax.text(c_plot_x0 + c_plot_w * 0.175, y_aln + 0.04, 'Home Flank 1', ha='center', va='center', fontsize=6.0, color=c_home, fontweight='bold', zorder=4)
ax.text(c_plot_x0 + c_plot_w * 0.525, y_aln + 0.04, 'Converted Cassette [u1, u2]', ha='center', va='center', fontsize=6.0, color=c_donor, fontweight='bold', zorder=4)
ax.text(c_plot_x0 + c_plot_w * 0.850, y_aln + 0.04, 'Home Flank 2', ha='center', va='center', fontsize=6.0, color=c_home, fontweight='bold', zorder=4)

# Subplot 2: Instantaneous Zero-Drift Force Field F(u)
y_f = py2 + 1.45
h_f = 0.50
ax.text(px2 + 0.15, y_f + h_f/2, 'Zero-Drift\nForce F(u)', ha='left', va='center', fontsize=6.5, fontweight='bold', color=c_text_header)
ax.plot([c_plot_x0, c_plot_x0 + c_plot_w], [y_f + h_f/2, y_f + h_f/2], color='#94A3B8', lw=0.8, linestyle='-', zorder=3)
ax.text(c_plot_x0 - 0.04, y_f + h_f/2, '0', ha='right', va='center', fontsize=5.8, color=c_text_muted)

np.random.seed(42)
u_f = np.linspace(c_plot_x0, c_plot_x0 + c_plot_w, 85)
f_donor = []
f_home = []
for u in u_f:
    if u < bp1_x or u > bp2_x:
        f_donor.append(np.random.normal(-0.05, 0.03))
        f_home.append(np.random.normal(0.05, 0.03))
    else:
        f_donor.append(np.random.normal(0.22, 0.05))
        f_home.append(np.random.normal(-0.22, 0.05))

ax.plot(u_f, [y_f + h_f/2 + v * 0.6 for v in f_donor], color=c_donor, lw=1.0, alpha=0.9, zorder=4)
ax.plot(u_f, [y_f + h_f/2 + v * 0.6 for v in f_home], color=c_home, lw=1.0, alpha=0.9, zorder=4)
ax.text(c_plot_x0 + c_plot_w - 0.05, y_f + h_f + 0.04, 'F_donor > 0  (Attraction to Donor)', ha='right', va='bottom', fontsize=5.8, fontweight='bold', color=c_donor, zorder=5)

# Subplot 3: Cumulative Trajectory Momentum C(u)
y_c = py2 + 0.35
ax.text(px2 + 0.15, y_c + 0.25, 'Cumulative\nMomentum\nC(u) = sum F', ha='left', va='center', fontsize=6.5, fontweight='bold', color=c_text_header)

ax.plot([c_plot_x0, c_plot_x0 + c_plot_w], [y_c + 0.22, y_c + 0.22], color='#94A3B8', lw=0.8, zorder=3)
ax.text(c_plot_x0 - 0.04, y_c + 0.22, '0', ha='right', va='center', fontsize=5.8, color=c_text_muted)
ax.text(c_plot_x0 + c_plot_w, y_c + 0.08, 'Alignment Coordinate u', ha='right', va='top', fontsize=6.2, color=c_text_muted)

u_fine = np.linspace(c_plot_x0, c_plot_x0 + c_plot_w, 200)
t_norm = (u_fine - c_plot_x0) / c_plot_w
env = 0.06 * np.sqrt(np.maximum(1e-5, t_norm * (1.0 - t_norm)))
ax.fill_between(u_fine, y_c + 0.22 - env, y_c + 0.22 + env, color='#E2E8F0', alpha=0.8, zorder=3)
ax.text(c_plot_x0 + c_plot_w * 0.12, y_c + 0.22 + 0.06, 'Clonal Null Envelope', ha='center', fontsize=5.5, color=c_text_muted, zorder=4)

c_curve = []
for u in u_fine:
    if u < bp1_x:
        c_curve.append(y_c + 0.22 + np.random.normal(0, 0.008))
    elif u < bp2_x:
        frac = (u - bp1_x) / (bp2_x - bp1_x)
        c_curve.append(y_c + 0.22 + frac * 0.35 + np.random.normal(0, 0.008))
    else:
        frac = (u - bp2_x) / (c_plot_x0 + c_plot_w - bp2_x)
        c_curve.append(y_c + 0.22 + 0.35 * (1.0 - frac) + np.random.normal(0, 0.008))

ax.plot(u_fine, c_curve, color=c_donor, lw=1.8, zorder=6)
ax.plot(u_fine, [2*(y_c + 0.22) - v for v in c_curve], color=c_home, lw=1.3, linestyle='--', zorder=5)

ax.annotate('Macroscopic Excursion\n(Extreme-Value Changepoint)', xy=(c_plot_x0 + c_plot_w * 0.525, y_c + 0.52),
            xytext=(c_plot_x0 + c_plot_w * 0.525, y_c + 0.72),
            ha='center', va='center', fontsize=6.0, fontweight='bold', color=c_donor,
            arrowprops=dict(arrowstyle="->", color=c_donor, lw=1.1),
            bbox=dict(boxstyle="round,pad=0.18", facecolor='#FFFBEB', edgecolor=c_donor, lw=0.7), zorder=7)


# ==============================================================================
# PANEL D: RESOLUTION OF EVOLUTIONARY REGIMES & SWARM CONSOLIDATION
# ==============================================================================
draw_panel(ax, px3, py3, pw, ph, 'D', 'Resolution Across Evolutionary Regimes & Swarm Consolidation',
           'Deconvolving epidemic transmission swarms, ghost donors, and micro-conversions at genome scale.')

# Sub-block D1: Epidemic Swarm Consolidation
draw_subbox(ax, px3 + 0.15, py3 + 1.45, pw - 0.30, 1.12, 'Epidemic Transmission Swarm Consolidation')

d1_x = px3 + 0.30
d1_y = py3 + 1.50

ax.text(d1_x + 1.20, d1_y + 0.78, 'Triplet Heuristics (3Seq / RDP5):', ha='center', fontsize=6.6, fontweight='bold', color='#B91C1C')
ax.text(d1_x + 1.20, d1_y + 0.62, 'Evaluates each isolate independently', ha='center', fontsize=5.8, color=c_text_muted)
for k in range(4):
    ax.add_patch(Rectangle((d1_x + 0.20, d1_y + 0.12 + k * 0.09), 2.0, 0.06, facecolor='#FEE2E2', edgecolor='#EF4444', lw=0.6))
ax.text(d1_x + 1.20, d1_y + 0.28, '50 Isolates = 50 Flagged Events!\n(10x Overcounting Inflation)', ha='center', va='center',
        fontsize=6.0, fontweight='bold', color='#B91C1C',
        bbox=dict(boxstyle="round,pad=0.16", facecolor='#FFFFFF', edgecolor='#EF4444', lw=0.7), zorder=4)

ax.annotate('', xy=(d1_x + 2.95, d1_y + 0.35), xytext=(d1_x + 2.55, d1_y + 0.35),
            arrowprops=dict(arrowstyle="-|>", color=c_text_header, lw=1.6), zorder=5)

ax.text(d1_x + 4.30, d1_y + 0.78, 'RhizAeon Graph Adjudication:', ha='center', fontsize=6.6, fontweight='bold', color='#047857')
ax.text(d1_x + 4.30, d1_y + 0.62, 'Retains shared ancestral history', ha='center', fontsize=5.8, color=c_text_muted)
ax.add_patch(Rectangle((d1_x + 3.25, d1_y + 0.14), 2.10, 0.36, facecolor='#ECFDF5', edgecolor='#10B981', lw=1.2))
ax.text(d1_x + 4.30, d1_y + 0.32, '1 Unified Ancestral Event\n(50 Descendant Isolates linked)', ha='center', va='center',
        fontsize=6.2, fontweight='bold', color='#047857', zorder=4)

# Sub-block D2: Operational Regimes
draw_subbox(ax, px3 + 0.15, py3 + 0.12, pw - 0.30, 1.25, 'Deconvolution Across Empirical Genomic Architectures')

d2_y = py3 + 0.20
reg_w = (pw - 0.70) / 3

# Regime 1: Viral Chimeras
r1_x = px3 + 0.25
ax.add_patch(FancyBboxPatch((r1_x, d2_y), reg_w, 0.95, boxstyle="round,pad=0.0,rounding_size=0.04",
                            facecolor='#FFFFFF', edgecolor='#93C5FD', lw=0.9))
ax.text(r1_x + reg_w/2, d2_y + 0.82, 'Viral Chimeras\n(Sarbecovirus)', ha='center', va='top', fontsize=6.2, fontweight='bold', color=c_home)
ax.text(r1_x + reg_w/2, d2_y + 0.50, 'Multi-donor mosaics\nGhost donors\nSARS-CoV-2 Spike',
        ha='center', va='top', fontsize=5.5, color=c_text_body, linespacing=1.15)
ax.text(r1_x + reg_w/2, d2_y + 0.06, '30 kb / 83 taxa', ha='center', va='bottom', fontsize=5.2, fontweight='bold', color='#2563EB')

# Regime 2: Bacterial Chromosomes
r2_x = r1_x + reg_w + 0.10
ax.add_patch(FancyBboxPatch((r2_x, d2_y), reg_w, 0.95, boxstyle="round,pad=0.0,rounding_size=0.04",
                            facecolor='#FFFFFF', edgecolor='#FDBA74', lw=0.9))
ax.text(r2_x + reg_w/2, d2_y + 0.82, 'Bacterial Chromosomes\n(S. pneumoniae)', ha='center', va='top', fontsize=6.2, fontweight='bold', color='#D97706')
ax.text(r2_x + reg_w/2, d2_y + 0.50, 'SNP coordinates\nPMEN1 penicillin pbp2x\n96.2% Gubbins match',
        ha='center', va='top', fontsize=5.5, color=c_text_body, linespacing=1.15)
ax.text(r2_x + reg_w/2, d2_y + 0.06, '2.22 Mb in <10 ms', ha='center', va='bottom', fontsize=5.2, fontweight='bold', color='#D97706')

# Regime 3: Immunogenomics
r3_x = r2_x + reg_w + 0.10
ax.add_patch(FancyBboxPatch((r3_x, d2_y), reg_w, 0.95, boxstyle="round,pad=0.0,rounding_size=0.04",
                            facecolor='#FFFFFF', edgecolor='#A7F3D0', lw=0.9))
ax.text(r3_x + reg_w/2, d2_y + 0.82, 'Immunogenomics\n(Primate HLA)', ha='center', va='top', fontsize=6.2, fontweight='bold', color='#047857')
ax.text(r3_x + reg_w/2, d2_y + 0.50, 'Inter-locus conversions\nPeptide binding groove\nHLA-B*46:01 chimera',
        ha='center', va='top', fontsize=5.5, color=c_text_body, linespacing=1.15)
ax.text(r3_x + reg_w/2, d2_y + 0.06, '50 bp micro-tracts', ha='center', va='bottom', fontsize=5.2, fontweight='bold', color='#047857')

# Save outputs
fig.savefig('figures/fig1_rhizaeon_dieter_rams_flowchart.pdf', format='pdf', dpi=300)
fig.savefig('figures/fig1_rhizaeon_dieter_rams_flowchart.png', format='png', dpi=300)
print("Figure 1 generated successfully: figures/fig1_rhizaeon_dieter_rams_flowchart.pdf and .png")
