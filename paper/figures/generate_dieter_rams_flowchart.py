#!/usr/bin/env python3
"""
generate_dieter_rams_flowchart.py
=================================
Generates the publication-grade Figure 1 for RhizAeon in Molecular Biology and Evolution (MBE).
Design Philosophy: Dieter Rams ("Weniger, aber besser" / "Less, but better")
- High-contrast, restrained, functional palette (Slate, Prussian Navy, Braun Orange, Emerald, Charcoal).
- Strict visual hierarchy, zero chart-junk, precision geometric alignment.
- Illustrates the complete 5-stage architectural pipeline from raw sequence alignment to polarized ARG.
"""

from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
import numpy as np

# Canvas Setup: 16.0 x 8.8 inches, 300 DPI
fig = plt.figure(figsize=(16.0, 8.8), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 16.0)
ax.set_ylim(0, 8.8)
ax.axis('off')

# ------------------------------------------------------------------------------
# 1. DIETER RAMS COLOR PALETTE
# ------------------------------------------------------------------------------
c_banner_bg     = '#111827'  # Deepest slate-gray / charcoal
c_banner_text   = '#FFFFFF'
c_banner_sub    = '#94A3B8'

# Main Panel Cards
c_card_bg       = '#F8FAFC'  # Very soft neutral gray-white
c_card_border   = '#CBD5E1'  # Subtle precision border
c_card_header   = '#0F172A'  # Primary title color
c_card_sub      = '#475569'  # Secondary label

# Inner Functional Blocks
c_block_bg      = '#FFFFFF'
c_block_border  = '#E2E8F0'
c_block_text    = '#1E293B'
c_block_detail  = '#475569'

# Accent Colors (Restrained, Functional)
c_orange_accent = '#EA580C'  # Braun International Orange: Recombinant / Breakpoint / Crossover
c_orange_soft   = '#FFF7ED'
c_orange_border = '#FDBA74'

c_blue_accent   = '#2563EB'  # Classic Tech Blue: Metric Space / Trajectory / Manifold
c_blue_soft     = '#EFF6FF'
c_blue_border   = '#93C5FD'

c_teal_accent   = '#0D9488'  # Precision Emerald / Teal: Polishing / ARG / Output
c_teal_soft     = '#F0FDFA'
c_teal_border   = '#5EEAD4'

c_purple_accent = '#7C3AED'  # Deep Violet: Unpolarized Root Token / Force Field
c_purple_soft   = '#F5F3FF'
c_purple_border = '#C4B5FD'

# Helper Functions
def draw_card(ax, x, y, w, h, bg=c_card_bg, border=c_card_border, radius=0.08, lw=1.2, zorder=1):
    patch = FancyBboxPatch((x, y), w, h,
                           boxstyle=f"round,pad=0.0,rounding_size={radius}",
                           facecolor=bg, edgecolor=border, linewidth=lw, zorder=zorder)
    ax.add_patch(patch)
    return patch

def draw_block(ax, x, y, w, h, title, subtitle=None, bg=c_block_bg, border=c_block_border,
               title_color=c_block_text, title_size=8.0, sub_color=c_block_detail, sub_size=6.6,
               radius=0.05, lw=1.0, zorder=2, title_weight='bold'):
    draw_card(ax, x, y, w, h, bg, border, radius=radius, lw=lw, zorder=zorder)
    if subtitle is None:
        ax.text(x + w/2, y + h/2, title, ha='center', va='center', fontsize=title_size,
                fontweight=title_weight, color=title_color, zorder=zorder+1)
    else:
        # Title anchored cleanly at the top
        ax.text(x + w/2, y + h - 0.20, title, ha='center', va='top', fontsize=title_size,
                fontweight=title_weight, color=title_color, zorder=zorder+1)
        # Subtitle / multi-line details anchored below title with generous breathing room
        ax.text(x + w/2, y + h - 0.44, subtitle, ha='center', va='top', fontsize=sub_size,
                color=sub_color, zorder=zorder+1, multialignment='center', linespacing=1.22)

def draw_badge(ax, x, y, w, h, text, bg=c_blue_soft, border=c_blue_border, fg=c_blue_accent,
               fontsize=7.2, zorder=3):
    draw_card(ax, x, y, w, h, bg, border, radius=0.04, lw=0.9, zorder=zorder)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fontsize,
            fontweight='bold', color=fg, zorder=zorder+1)

def draw_arrow(ax, x1, y1, x2, y2, color='#64748B', lw=1.6, zorder=5):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                mutation_scale=11, shrinkA=0, shrinkB=0),
                zorder=zorder)

# ------------------------------------------------------------------------------
# 2. TOP BANNER (DIETER RAMS TYPOGRAPHIC MASTHEAD)
# ------------------------------------------------------------------------------
draw_card(ax, 0.40, 7.95, 15.20, 0.62, bg=c_banner_bg, border=c_banner_bg, radius=0.06, zorder=2)

ax.text(0.65, 8.26, "RhizAeon: Physical Force Field Deconvolution of Sequence Manifolds",
        fontsize=11.2, fontweight='bold', color=c_banner_text, va='center', zorder=3)
ax.text(15.35, 8.26, "Continuous Metric Geometry  |  O(1) Prefix Tensors  |  Linear O(N·L) Complexity  |  Native Rust & WASM",
        fontsize=8.0, fontweight='bold', color=c_banner_sub, va='center', ha='right', zorder=3)

# ------------------------------------------------------------------------------
# 3. 5-COLUMN HORIZONTAL SWIMLANE LAYOUT (STRICT UNIFORM GRID)
# ------------------------------------------------------------------------------
col_w = 2.82
col_gap = 0.26
col_x = [0.40 + i * (col_w + col_gap) for i in range(5)]
base_y = 0.45
base_h = 7.25

# Standardized vertical grid coordinates for all columns
# Column total height = 7.25 (y from 0.45 to 7.70)
# Header: 7.15 to 7.70
# Block 1: 5.68 to 7.03 (h = 1.35)
# Block 2: 4.21 to 5.53 (h = 1.32)
# Block 3: 2.74 to 4.06 (h = 1.32)
# Block 4: 1.25 to 2.57 (h = 1.32)
# Badge:   0.62 to 1.07 (h = 0.45)

b1_y, b1_h = 5.68, 1.35
b2_y, b2_h = 4.21, 1.32
b3_y, b3_h = 2.74, 1.32
b4_y, b4_h = 1.25, 1.32
bdg_y, bdg_h = 0.62, 0.45

col_titles = [
    ("A", "Constant-Time Metric Engine", "Prefix Distance Tensors & SNP Indexing"),
    ("B", "Continuous Metric Manifolds", "Classical MDS & Trajectory Mechanics"),
    ("C", "Force Field Deconvolution", "Restoring Potentials & Root Token Sink"),
    ("D", "Polishing & ARG Consolidation", "Profile ML & Procrustes Swarm Grouping"),
    ("E", "Resolved Outputs & Surveillance", "Polarized ARG & Real-Time Dashboard")
]

for i in range(5):
    cx = col_x[i]
    draw_card(ax, cx, base_y, col_w, base_h, bg=c_card_bg, border=c_card_border, radius=0.08, zorder=1)
    
    # Column Letter Circle Badge
    letter, title, subtitle = col_titles[i]
    draw_card(ax, cx + 0.15, base_y + base_h - 0.45, 0.32, 0.32, bg=c_banner_bg, border=c_banner_bg, radius=0.16, zorder=3)
    ax.text(cx + 0.31, base_y + base_h - 0.29, letter, color='#FFFFFF', fontsize=9.5, fontweight='bold', ha='center', va='center', zorder=4)
    
    # Title & Subtitle
    ax.text(cx + 0.55, base_y + base_h - 0.22, title, fontsize=9.0, fontweight='bold', color=c_card_header, va='center', zorder=3)
    ax.text(cx + 0.55, base_y + base_h - 0.42, subtitle, fontsize=7.0, color=c_card_sub, va='center', zorder=3)
    
    # Header divider rule
    ax.plot([cx + 0.15, cx + col_w - 0.15], [base_y + base_h - 0.56, base_y + base_h - 0.56], color=c_card_border, lw=0.8, zorder=2)

# Inter-Column Global Pipeline Arrows (Connecting Columns at Mid-Heights)
arrow_ys = [4.87, 4.87, 4.87, 4.87]
for i in range(4):
    x_start = col_x[i] + col_w
    x_end = col_x[i+1]
    draw_arrow(ax, x_start + 0.04, arrow_ys[i], x_end - 0.04, arrow_ys[i], color='#94A3B8', lw=2.0)

# ==============================================================================
# COLUMN A: CONSTANT-TIME METRIC ENGINE
# ==============================================================================
xa = col_x[0]

# Block A1: Multiple Sequence Alignment Input
draw_block(ax, xa + 0.15, b1_y, col_w - 0.30, b1_h,
           "1. Multiple Sequence Alignment",
           "FASTA Matrix X in {A, C, G, T, -}^(N x L)\nPathogen Cohort: N Taxa x L Coordinates\nDirect Ingestion: Zero Pre-Screen Filtering",
           bg='#FFFFFF', border=c_card_border)

# Block A2: Segregating SNP Compression Engine
draw_block(ax, xa + 0.15, b2_y, col_w - 0.30, b2_h,
           "2. SNP-Compressed Coordinates",
           "Extract Segregating Polymorphic Sites |S|\nO(log |S|) Interval Binary Search Mapping\n220x Memory Reduction (88 GB -> 400 MB)",
           bg=c_blue_soft, border=c_blue_border, title_color=c_blue_accent)

# Block A3: 3D Cumulative Prefix Distance Tensor
draw_block(ax, xa + 0.15, b3_y, col_w - 0.30, b3_h,
           "3. Prefix Distance Tensor",
           "Precompute T in R^(N x N x L) in Single Pass\nT_diff(i, j, u) = sum_{k=0}^u D_ij(k)\nExact Linear O(L·N^2) Pre-Accumulation",
           bg='#FFFFFF', border=c_card_border)

# Block A4: O(1) Interval Distance Engine
draw_block(ax, xa + 0.15, b4_y, col_w - 0.30, b4_h,
           "4. Exact O(1) Sub-Interval Query",
           "Arbitrary Window D(s1, s2) = T(s2) - T(s1)\nPairwise TN93 or Hamming Distances in O(1)\nEliminates Repeated Window Re-Scans\nEvaluates Cutpoints in Microseconds",
           bg=c_teal_soft, border=c_teal_border, title_color=c_teal_accent)

# Connectors inside Col A
draw_arrow(ax, xa + col_w/2, b1_y, xa + col_w/2, b2_y + b2_h, color='#94A3B8')
draw_arrow(ax, xa + col_w/2, b2_y, xa + col_w/2, b3_y + b3_h, color='#94A3B8')
draw_arrow(ax, xa + col_w/2, b3_y, xa + col_w/2, b4_y + b4_h, color='#94A3B8')

# Badge in Col A
draw_badge(ax, xa + 0.20, bdg_y, col_w - 0.40, bdg_h, "O(1) Interval Retrieval", bg=c_teal_soft, border=c_teal_border, fg=c_teal_accent)

# ==============================================================================
# COLUMN B: CONTINUOUS METRIC MANIFOLDS
# ==============================================================================
xb = col_x[1]

# Block B1: Classical MDS & Buneman Metric Space
draw_block(ax, xb + 0.15, b1_y, col_w - 0.30, b1_h,
           "1. Canonical Metric Basis (R^k)",
           "Doubly Centered Gram Matrix B = -0.5 H D^2 H\nEigendecomposition: Leading k Axes (k = 3-4)\nBuneman Quartet Four-Point Additive Space\nStationary Global Frame (Zero Procrustes)",
           bg='#FFFFFF', border=c_card_border)

# Block B2: Clonal Corridors vs Recombinant Deflection (Schematic mini-diagram)
draw_card(ax, xb + 0.15, b2_y, col_w - 0.30, b2_h, bg='#FFFFFF', border=c_card_border, radius=0.05, zorder=2)
ax.text(xb + col_w/2, b2_y + b2_h - 0.18, "2. Sequence Trajectory Projection", ha='center', va='top', fontsize=8.0, fontweight='bold', color=c_card_header, zorder=3)
ax.text(xb + col_w/2, b2_y + b2_h - 0.38, "Closed-form Gower Out-of-Sample Mapping", ha='center', va='top', fontsize=6.6, color=c_card_sub, zorder=3)

# Mini trajectory visual inside card B2
tx = np.linspace(xb + 0.35, xb + col_w - 0.35, 60)
p1_y_line = b2_y + 0.58
p2_y_line = b2_y + 0.20

# Corridor 1 (Parent A)
ax.plot(tx, np.full_like(tx, p1_y_line), color=c_blue_accent, lw=1.8, alpha=0.35, zorder=3)
ax.text(xb + 0.38, p1_y_line + 0.05, "Parent A Corridor", fontsize=6.2, fontweight='bold', color=c_blue_accent, zorder=4)
# Corridor 2 (Parent B)
ax.plot(tx, np.full_like(tx, p2_y_line), color=c_teal_accent, lw=1.8, alpha=0.35, zorder=3)
ax.text(xb + 0.38, p2_y_line + 0.05, "Parent B Corridor", fontsize=6.2, fontweight='bold', color=c_teal_accent, zorder=4)

# Recombinant Trajectory jumping from A to B
split_idx = 28
rx1 = tx[:split_idx]
ry1 = np.full_like(rx1, p1_y_line)
rx2 = tx[split_idx:split_idx+8]
ry2 = np.linspace(p1_y_line, p2_y_line, len(rx2))
rx3 = tx[split_idx+8:]
ry3 = np.full_like(rx3, p2_y_line)
rx_full = np.concatenate([rx1, rx2, rx3])
ry_full = np.concatenate([ry1, ry2, ry3])

ax.plot(rx_full, ry_full, color=c_orange_accent, lw=2.0, zorder=5)
# Breakpoint marker
bp_x = tx[split_idx+4]
bp_y = (p1_y_line + p2_y_line) / 2.0
ax.plot(bp_x, bp_y, 'o', color=c_orange_accent, markersize=4.5, zorder=6)
ax.text(bp_x + 0.08, bp_y + 0.02, "Crossover s*", fontsize=6.2, fontweight='bold', color=c_orange_accent, zorder=6)

# Block B3: Kinetic Energy & Velocity Vectors
draw_block(ax, xb + 0.15, b3_y, col_w - 0.30, b3_h,
           "3. Kinetic Velocity & Energy",
           "Velocity v_i(s) = [y(s+Delta) - y(s-Delta)] / 2Delta\nInstantaneous Kinetic Energy K_i(s) = ||v_i(s)||^2\nClonal Lineages: v_i(s) ~ 0 (Stationary)\nRecombinants: Sharp Velocity Spikes",
           bg=c_blue_soft, border=c_blue_border, title_color=c_blue_accent)

# Block B4: Breakpoint Flux Density Field
draw_block(ax, xb + 0.15, b4_y, col_w - 0.30, b4_h,
           "4. Breakpoint Flux Density Field",
           "Phi_bp(s) = sum_e Z_e · K_sigma(s - s_e)\nIntegrates Local Turbulence Across Cohort\nDetects Recombination Hotspots (Spike NTD)\nIdentifies Catalytic Coldspots (RdRp Core)",
           bg='#FFFFFF', border=c_card_border)

# Connectors inside Col B
draw_arrow(ax, xb + col_w/2, b1_y, xb + col_w/2, b2_y + b2_h, color='#94A3B8')
draw_arrow(ax, xb + col_w/2, b2_y, xb + col_w/2, b3_y + b3_h, color='#94A3B8')
draw_arrow(ax, xb + col_w/2, b3_y, xb + col_w/2, b4_y + b4_h, color='#94A3B8')

# Badge in Col B
draw_badge(ax, xb + 0.20, bdg_y, col_w - 0.40, bdg_h, "Continuous Trajectories", bg=c_blue_soft, border=c_blue_border, fg=c_blue_accent)

# ==============================================================================
# COLUMN C: FORCE FIELD DECONVOLUTION
# ==============================================================================
xc = col_x[2]

# Block C1: Hookean Restoring Potential & Velocity Reversal
draw_block(ax, xc + 0.15, b1_y, col_w - 0.30, b1_h,
           "1. Dynamic Velocity Reversal",
           "Restoring Potential Zero-Crossing:\nDelta f = f_right(s) - f_left(s)\nStandardized Z-Score: Z_flip = |Delta f · d_hat| / SE\nDistinguishes Active Reticulation from Drift",
           bg='#FFFFFF', border=c_card_border)

# Block C2: The Unpolarized Root Token Sink (Crucial differentiator!)
draw_block(ax, xc + 0.15, b2_y, col_w - 0.30, b2_h,
           "2. Ancestral Root Token Sink",
           "Token [ROOT] Anchored at Origin m_0 = [0, 0, 0, 0]\nAbsorbs Private Mutations & Heterotachy Spikes\nAPOBEC Hypermutation Naturally Pools into Root\nProvably Suppresses Clonal False Positives",
           bg=c_purple_soft, border=c_purple_border, title_color=c_purple_accent)

# Block C3: Dyadic Scale-Space Matched Filter
draw_block(ax, xc + 0.15, b3_y, col_w - 0.30, b3_h,
           "3. Dyadic Scale-Space Filter",
           "Bipolar Wavelet Convolution: dy * h_tau\nDyadic Scales tau in {15, 25, 50, 100, 200, 500}\nLocalizes Micro-Conversions Down to 15 Codons\nTaylor Running Drawdown Bounded Peaks",
           bg='#FFFFFF', border=c_card_border)

# Block C4: Directional Chord Collinearity
draw_block(ax, xc + 0.15, b4_y, col_w - 0.30, b4_h,
           "4. Directional Collinearity Filter",
           "Four-Point Metric Polarization: cos(theta) >= 0.40\nRejects Orthogonal Noise & Asymmetric Drift\nGuarantees True Anti-Parallel Bipolarity\nEliminates False Calls on Long Lineages",
           bg=c_orange_soft, border=c_orange_border, title_color=c_orange_accent)

# Connectors inside Col C
draw_arrow(ax, xc + col_w/2, b1_y, xc + col_w/2, b2_y + b2_h, color='#94A3B8')
draw_arrow(ax, xc + col_w/2, b2_y, xc + col_w/2, b3_y + b3_h, color='#94A3B8')
draw_arrow(ax, xc + col_w/2, b3_y, xc + col_w/2, b4_y + b4_h, color='#94A3B8')

# Badge in Col C
draw_badge(ax, xc + 0.20, bdg_y, col_w - 0.40, bdg_h, "False Positive Immunity (0.415%)", bg=c_purple_soft, border=c_purple_border, fg=c_purple_accent)

# ==============================================================================
# COLUMN D: POLISHING & ARG CONSOLIDATION
# ==============================================================================
xd = col_x[3]

# Block D1: Single-Base Profile Likelihood Polisher
draw_block(ax, xd + 0.15, b1_y, col_w - 0.30, b1_h,
           "1. Profile Likelihood Polishing",
           "Single-Nucleotide Likelihood Profile L(x | D, H)\nEvaluates Exact Transition Probabilities < 0.1 ms\nTightly Bounds Physical Plateaus [x_left, x_right]\nExposes Inter-Mutation Uncertainty Intervals",
           bg=c_teal_soft, border=c_teal_border, title_color=c_teal_accent)

# Block D2: Orthogonal Procrustes Swarm Grouping
draw_block(ax, xd + 0.15, b2_y, col_w - 0.30, b2_h,
           "2. Procrustes Swarm Grouping",
           "Shared Ancestral Recombination Grouping\nAligns In-Tract Subspaces via Orthogonal Rotations\nCollapses K Descendant Isolates into Single Event\nPrevents Triplet Heuristic 5-10x Event Inflation",
           bg='#FFFFFF', border=c_card_border)

# Block D3: Ghost Donor & Unsampled Lineage Isolation
draw_block(ax, xd + 0.15, b3_y, col_w - 0.30, b3_h,
           "3. Ghost Donor Isolation",
           "Studentized Coordinate Dislocation Z > 3.0\nSolves Triplet Closed-World Entrapment\nDetects Recombination When Both Parents Unsampled\nPrevents Misidentifying Recombinants as Parents",
           bg=c_orange_soft, border=c_orange_border, title_color=c_orange_accent)

# Block D4: RMT Bulk Spectral Edge Halting
draw_block(ax, xd + 0.15, b4_y, col_w - 0.30, b4_h,
           "4. Marchenko-Pastur Halting",
           "Random Matrix Theory Bulk Edge Stopping Rule\nSingular Value Ratio s1 / s_bulk >= SNR_cutoff\nAutomated Termination: Deflates Peeling Modes\nPrevents Overfitting to Stochastic Noise",
           bg='#FFFFFF', border=c_card_border)

# Connectors inside Col D
draw_arrow(ax, xd + col_w/2, b1_y, xd + col_w/2, b2_y + b2_h, color='#94A3B8')
draw_arrow(ax, xd + col_w/2, b2_y, xd + col_w/2, b3_y + b3_h, color='#94A3B8')
draw_arrow(ax, xd + col_w/2, b3_y, xd + col_w/2, b4_y + b4_h, color='#94A3B8')

# Badge in Col D
draw_badge(ax, xd + 0.20, bdg_y, col_w - 0.40, bdg_h, "Exact Spatial Resolution (16 nt)", bg=c_teal_soft, border=c_teal_border, fg=c_teal_accent)

# ==============================================================================
# COLUMN E: RESOLVED OUTPUTS & SURVEILLANCE
# ==============================================================================
xe = col_x[4]

# Block E1: Polarized Ancestral Recombination Graph (ARG)
draw_card(ax, xe + 0.15, b1_y, col_w - 0.30, b1_h, bg='#FFFFFF', border=c_card_border, radius=0.05, zorder=2)
ax.text(xe + col_w/2, b1_y + b1_h - 0.18, "1. Polarized ARG Network", ha='center', va='top', fontsize=8.0, fontweight='bold', color=c_card_header, zorder=3)
ax.text(xe + col_w/2, b1_y + b1_h - 0.38, "Reconciles Conflicting Lineages", ha='center', va='top', fontsize=6.6, color=c_card_sub, zorder=3)

# Mini ARG Graphic inside Block E1
n1 = (xe + 0.52, b1_y + 0.58)
n2 = (xe + 1.41, b1_y + 0.58)
n3 = (xe + 2.30, b1_y + 0.58)
r1 = (xe + 0.96, b1_y + 0.20)

ax.plot([n1[0], r1[0]], [n1[1], r1[1]], color=c_blue_accent, lw=1.6, zorder=4)
ax.plot([n2[0], r1[0]], [n2[1], r1[1]], color=c_orange_accent, lw=1.6, linestyle='--', zorder=4)
ax.plot([n2[0], n2[0]], [n2[1], b1_y + 0.18], color='#CBD5E1', lw=1.2, zorder=4)
ax.plot([n3[0], n3[0]], [n3[1], b1_y + 0.18], color='#CBD5E1', lw=1.2, zorder=4)

ax.plot(n1[0], n1[1], 'o', color=c_blue_accent, markersize=4.5, zorder=5)
ax.plot(n2[0], n2[1], 'o', color=c_orange_accent, markersize=4.5, zorder=5)
ax.plot(n3[0], n3[1], 'o', color='#94A3B8', markersize=4.5, zorder=5)
ax.plot(r1[0], r1[1], 's', color=c_orange_accent, markersize=5.5, zorder=5)

ax.text(n1[0], n1[1]+0.06, "Parent A", fontsize=5.8, ha='center', color=c_blue_accent, fontweight='bold', zorder=6)
ax.text(n2[0], n2[1]+0.06, "Donor B", fontsize=5.8, ha='center', color=c_orange_accent, fontweight='bold', zorder=6)
ax.text(r1[0], r1[1]-0.12, "Recombinant R", fontsize=5.8, ha='center', color=c_orange_accent, fontweight='bold', zorder=6)

# Block E2: Machine-Readable Mosaic Formula
draw_block(ax, xe + 0.15, b2_y, col_w - 0.30, b2_h,
           "2. Mosaic Nomenclature",
           "Standardized Machine-Readable Syntax:\nHG:Rec[P1(1..u1) + P2(u1..u2) + P1(u2..L)]\nDeterministic Lineage & Cassette Labeling\nReplaces Ad-Hoc Public Health Aliasing",
           bg=c_blue_soft, border=c_blue_border, title_color=c_blue_accent)

# Block E3: Native Rust CLI & WebAssembly Dashboard
draw_block(ax, xe + 0.15, b3_y, col_w - 0.30, b3_h,
           "3. Interactive WebAssembly App",
           "Native Rust Binary (CLI & Library)\nClient-Side WebAssembly (WASM) & WebGPU\nZero Server Uploads: 100% In-Browser Privacy\nSynchronized Inspection, MDS Plots & NJ Trees",
           bg=c_teal_soft, border=c_teal_border, title_color=c_teal_accent)

# Block E4: Surveillance Latency & Throughput
draw_block(ax, xe + 0.15, b4_y, col_w - 0.30, b4_h,
           "4. Real-Time Surveillance Latency",
           "> 30 Genomes/s Throughput on Laptops\nSub-Millisecond Single-Genome Screening\nScales to Hundreds of Thousands of Isolates\nOperates without Supercomputers or Clusters",
           bg='#FFFFFF', border=c_card_border)

# Connectors inside Col E
draw_arrow(ax, xe + col_w/2, b1_y, xe + col_w/2, b2_y + b2_h, color='#94A3B8')
draw_arrow(ax, xe + col_w/2, b2_y, xe + col_w/2, b3_y + b3_h, color='#94A3B8')
draw_arrow(ax, xe + col_w/2, b3_y, xe + col_w/2, b4_y + b4_h, color='#94A3B8')

# Badge in Col E
draw_badge(ax, xe + 0.20, bdg_y, col_w - 0.40, bdg_h, "Real-Time Desktop Surveillance", bg=c_teal_soft, border=c_teal_border, fg=c_teal_accent)

# ------------------------------------------------------------------------------
# 4. SAVE OUTPUTS (VECTOR PDF & 300 DPI RASTER PNG)
# ------------------------------------------------------------------------------
out_dir = Path("/Users/sergei/Projects/rhizaeon-code/paper/figures")
out_dir.mkdir(parents=True, exist_ok=True)
pdf_path = out_dir / "fig1_rhizaeon_dieter_rams_flowchart.pdf"
png_path = out_dir / "fig1_rhizaeon_dieter_rams_flowchart.png"

plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.04)
plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight', pad_inches=0.04)
plt.close(fig)

print(f"[OK] Successfully generated Dieter Rams Flowchart:")
print(f"     PDF: {pdf_path}")
print(f"     PNG: {png_path}")
