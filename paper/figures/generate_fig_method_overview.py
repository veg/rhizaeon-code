"""
paper/figures/generate_fig_method_overview.py
=============================================
Generates the publication-ready conceptual method flowchart and architectural
pipeline figure (Figure 1) for Section 2 of the RhizAeon manuscript.

Design Mandates:
- 50% LESS TEXT: 1 clear, punchy conceptual bullet per section (zero math equations).
- LARGE, PROMINENT, HIGH-CONTRAST FONTS:
    * Main Title: 22.0 pt bold
    * Subtitle: 13.5 pt italic
    * Card Banners: 13.5 pt bold (centered, single string, generous margin)
    * Section Headers: 13.5 pt bold
    * Body Bullets: 12.5 pt (strict max 28-30 chars/line)
    * Axis Labels: 11.0 pt bold
    * Ticks: 9.5 pt
    * Badges & Annotations: 10.2-11.0 pt bold
- ZERO OVERFLOW & ZERO COLLISIONS:
    * All headings and bullets wrapped strictly <= 30 chars per line.
    * Generous margins (>=0.40 in) inside cards on left, right, top, and bottom.
    * Subplot Y-axes positioned with >=0.20 in clearance from card borders.
    * Landmark genome bars have clean labels ("Landmark 1", "Seq 2", etc.).
    * Stage 4 Genome Track has visually proportional cassette (1400 nt / 3000 nt)
      with clean short labels ("Home", "Donor", "Home") and well-separated breakpoint ticks.
    * Stage 4 Event Roster Badge is a centered 2-line bounded card with ample padding.
- Output: 300 DPI vector PDF and high-res PNG.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Output paths
script_dir = Path(__file__).resolve().parent
pdf_path = script_dir / "fig_method_overview.pdf"
png_path = script_dir / "fig_method_overview.png"

# Global typography and styling
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['text.color'] = '#1A1A1A'
plt.rcParams['axes.edgecolor'] = '#CCCCCC'

fig_w = 16.5
fig_h = 10.5
fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis('off')

# Color palette
C_STAGE1 = "#0072B2"  # Blue
C_STAGE2 = "#D55E00"  # Vermilion/Orange
C_STAGE3 = "#009E73"  # Green
C_STAGE4 = "#782869"  # Deep Purple
C_ROOT   = "#555555"  # Dark Gray
C_BG1    = "#F4F8FB"  # Light Blue Tint
C_BG2    = "#FDF8F2"  # Light Orange Tint
C_BG3    = "#F2F9F5"  # Light Green Tint
C_BG4    = "#FAF4FA"  # Light Purple Tint

def draw_card(ax, x, y, w, h, bg_color, border_color, banner_text, tag_color):
    """Draws a structured column card with a centered, bold header banner."""
    # Outer card
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.20",
        facecolor=bg_color, edgecolor=border_color, linewidth=1.8, zorder=1
    )
    ax.add_patch(rect)
    
    # Header banner (centered text to completely prevent any text overlap)
    header_w = w * 0.94
    header_h = 0.60
    header_x = x + (w - header_w) / 2
    header_y = y + h - header_h - 0.16
    banner = patches.FancyBboxPatch(
        (header_x, header_y), header_w, header_h,
        boxstyle="round,pad=0.01,rounding_size=0.14",
        facecolor=tag_color, edgecolor='none', zorder=2
    )
    ax.add_patch(banner)
    
    ax.text(header_x + header_w * 0.5, header_y + header_h * 0.5, banner_text,
            fontsize=13.5, fontweight='bold', color='white', va='center', ha='center', zorder=3)

# Column coordinates
col_w = 3.70
gap = 0.38
start_x = 0.40
card_y = 0.40
card_h = 8.90

x1 = start_x
x2 = x1 + col_w + gap
x3 = x2 + col_w + gap
x4 = x3 + col_w + gap

# Main title and subtitle
ax.text(fig_w * 0.5, 10.10, "RhizAeon Computational Architecture and Recombination Inference Pipeline",
        fontsize=21.5, fontweight='bold', ha='center', va='center', color='#111111')
ax.text(fig_w * 0.5, 9.68, "From Multiple Sequence Alignment to Exact Characterization of Mosaic Reticulations in Linear Time",
        fontsize=13.2, style='italic', ha='center', va='center', color='#555555')


# ==============================================================================
# CARD 1: STAGE 1 - METRIC MANIFOLD EMBEDDING
# ==============================================================================
draw_card(ax, x1, card_y, col_w, card_h, C_BG1, C_STAGE1,
          "STAGE 1: Metric Embedding", C_STAGE1)

# Section A: Landmarks
y_cursor = card_y + card_h - 1.02
ax.text(x1 + 0.28, y_cursor, "A. Landmark Sampling",
        fontsize=13.5, fontweight='bold', color=C_STAGE1, va='top')
y_cursor -= 0.34
ax.text(x1 + 0.30, y_cursor, "• Minimax sampling selects\n  diverse reference genomes.",
        fontsize=12.5, color='#222222', va='top')

# Subplot 1A: Landmark Genome Bars
plot1a_top = y_cursor - 0.56
plot1a_h = 1.65
ax_sub1 = fig.add_axes([(x1 + 0.35)/fig_w, (plot1a_top - plot1a_h)/fig_h, 3.00/fig_w, plot1a_h/fig_h])
ax_sub1.set_xlim(0, 10)
ax_sub1.set_ylim(0, 5)
ax_sub1.axis('off')

np.random.seed(42)
labels_1a = ["Landmark 1", "Seq 2", "Landmark 2", "Seq 4", "Landmark 3"]
for row in range(5):
    y_row = 4 - row * 0.88
    is_landmark = row in [0, 2, 4]
    row_col = C_STAGE1 if is_landmark else '#777777'
    ax_sub1.text(0.1, y_row + 0.12, labels_1a[row],
                 fontsize=10.0, fontweight='bold' if is_landmark else 'normal', color=row_col)
    rect_bar = patches.Rectangle((4.5, y_row), 5.2, 0.48, facecolor='#EBF3FA' if is_landmark else '#F4F4F4',
                                 edgecolor=row_col, linewidth=1.3)
    ax_sub1.add_patch(rect_bar)
    for tick in np.linspace(4.7, 9.4, 7):
        if np.random.rand() > 0.4:
            ax_sub1.plot([tick, tick], [y_row, y_row + 0.48], color='#D55E00' if is_landmark else '#AAAAAA', lw=1.5)

# Section B: Continuous Metric Trajectories
y_cursor = plot1a_top - plot1a_h - 0.38
ax.text(x1 + 0.28, y_cursor, "B. Metric Trajectories",
        fontsize=13.5, fontweight='bold', color=C_STAGE1, va='top')
y_cursor -= 0.34
ax.text(x1 + 0.30, y_cursor, "• MDS & Gower projection track\n  continuous spatial trajectories.",
        fontsize=12.5, color='#222222', va='top')

# Subplot 1B: Continuous Trajectory in Metric Space
plot1b_top = y_cursor - 0.56
plot1b_h = 2.40
ax_sub2 = fig.add_axes([(x1 + 0.62)/fig_w, (plot1b_top - plot1b_h)/fig_h, 2.62/fig_w, plot1b_h/fig_h])
ax_sub2.set_xlim(-1.2, 1.2)
ax_sub2.set_ylim(-0.95, 0.95)
ax_sub2.set_facecolor('#FCFDFE')
for spine in ax_sub2.spines.values():
    spine.set_edgecolor('#CCCCCC')
    spine.set_linewidth(1.0)
ax_sub2.axhline(0, color='#E0E0E0', lw=0.9, ls='--')
ax_sub2.axvline(0, color='#E0E0E0', lw=0.9, ls='--')

ax_sub2.scatter([0.65], [0.40], color=C_STAGE1, s=120, zorder=4)
ax_sub2.scatter([-0.65], [-0.40], color='#D55E00', s=120, zorder=4)
ax_sub2.text(0.65, 0.62, 'Donor Lineage', fontsize=10.5, ha='center', fontweight='bold', color=C_STAGE1)
ax_sub2.text(-0.65, -0.68, 'Home Lineage', fontsize=10.5, ha='center', fontweight='bold', color='#D55E00')

t_vals = np.linspace(-0.65, 0.65, 45)
arc_y = -0.40 + 0.80 / (1 + np.exp(-10 * t_vals)) + 0.04 * np.sin(10 * t_vals)
ax_sub2.plot(t_vals, arc_y, color='#009E73', lw=2.8, zorder=3)
ax_sub2.scatter(t_vals[::8], arc_y[::8], color='#009E73', s=36, zorder=5)
ax_sub2.text(0.0, -0.16, 'Recombinant Trajectory', fontsize=10.0, ha='center', color='#009E73', fontweight='bold')
ax_sub2.set_xlabel('Phylogenetic Split Axis 1', fontsize=11.0, fontweight='bold', labelpad=3)
ax_sub2.set_ylabel('Split Axis 2', fontsize=11.0, fontweight='bold', labelpad=3)
ax_sub2.tick_params(labelsize=9.5)


# ==============================================================================
# CARD 2: STAGE 2 - DIRECTED ATTENTION & ROOT SINK DECOUPLING
# ==============================================================================
draw_card(ax, x2, card_y, col_w, card_h, C_BG2, C_STAGE2,
          "STAGE 2: Directed Attention", C_STAGE2)

# Section A: Tree-RoPE & Prior
y_cursor = card_y + card_h - 1.02
ax.text(x2 + 0.28, y_cursor, "A. Tree-RoPE & Prior",
        fontsize=13.5, fontweight='bold', color=C_STAGE2, va='top')
y_cursor -= 0.34
ax.text(x2 + 0.30, y_cursor, "• Tree-RoPE encodes distances;\n  Markov prior elevates key SNPs.",
        fontsize=12.5, color='#222222', va='top')

# Subplot 2A: Tree-RoPE rotary dials
plot2a_top = y_cursor - 0.56
plot2a_h = 1.65
ax_sub3 = fig.add_axes([(x2 + 0.35)/fig_w, (plot2a_top - plot2a_h)/fig_h, 3.00/fig_w, plot2a_h/fig_h])
ax_sub3.set_xlim(-1.5, 3.5)
ax_sub3.set_ylim(-1.35, 1.35)
ax_sub3.axis('off')

for idx, (cx, label, theta, col, subtext) in enumerate([
    (-0.3, "Sister Clade", 0.2, '#0072B2', "In-Phase (Aligned)"),
    (2.3, "Divergent Clade", 1.8, '#D55E00', "Out of Phase")
]):
    circle = plt.Circle((cx, 0), 0.85, facecolor='#FAF7F0', edgecolor=col, lw=1.8)
    ax_sub3.add_patch(circle)
    ax_sub3.plot([cx, cx + 0.85*np.cos(theta)], [0, 0.85*np.sin(theta)], color=col, lw=2.6)
    ax_sub3.plot([cx, cx + 0.85], [0, 0], color='#888888', lw=1.2, ls=':')
    ax_sub3.scatter([cx + 0.85*np.cos(theta)], [0.85*np.sin(theta)], color=col, s=48, zorder=5)
    ax_sub3.text(cx, -1.22, label, fontsize=10.5, ha='center', fontweight='bold', color=col)
    ax_sub3.text(cx, 1.05, subtext, fontsize=9.5, ha='center', color='#444444', fontweight='bold')

# Section B: Ancestral Root Sink
y_cursor = plot2a_top - plot2a_h - 0.38
ax.text(x2 + 0.28, y_cursor, "B. Ancestral Root Sink",
        fontsize=13.5, fontweight='bold', color=C_STAGE2, va='top')
y_cursor -= 0.34
ax.text(x2 + 0.30, y_cursor, "• Origin token absorbs bursts,\n  strictly eliminating false donors.",
        fontsize=12.5, color='#222222', va='top')

# Subplot 2B: Attention Allocation Stacked Bar
plot2b_top = y_cursor - 0.56
plot2b_h = 2.15
ax_sub4 = fig.add_axes([(x2 + 0.62)/fig_w, (plot2b_top - plot2b_h)/fig_h, 2.62/fig_w, plot2b_h/fig_h])
sites = np.arange(1, 15)
att_home = np.array([0.7, 0.65, 0.68, 0.15, 0.12, 0.10, 0.14, 0.12, 0.65, 0.70, 0.20, 0.68, 0.72, 0.70])
att_donor = np.array([0.15, 0.20, 0.18, 0.75, 0.80, 0.82, 0.78, 0.80, 0.20, 0.15, 0.15, 0.18, 0.14, 0.15])
att_root = 1.0 - att_home - att_donor
att_root[10] = 0.65
att_home[10] = 0.20
att_donor[10] = 0.15

ax_sub4.bar(sites, att_donor, color='#0072B2', label='Donor', width=0.75)
ax_sub4.bar(sites, att_home, bottom=att_donor, color='#D55E00', label='Home', width=0.75)
ax_sub4.bar(sites, att_root, bottom=att_donor + att_home, color='#555555', label='Root Sink', width=0.75)

ax_sub4.set_ylim(0, 1.05)
ax_sub4.set_xlabel('Nucleotide Position', fontsize=11.0, fontweight='bold', labelpad=3)
ax_sub4.set_ylabel('Attention Share', fontsize=11.0, fontweight='bold', labelpad=3)
ax_sub4.tick_params(labelsize=9.5)
ax_sub4.legend(loc='lower center', bbox_to_anchor=(0.5, -0.48), ncol=3, frameon=False, fontsize=9.2)


# ==============================================================================
# CARD 3: STAGE 3 - THREE-TIER HIERARCHICAL STATISTICAL INFERENCE
# ==============================================================================
draw_card(ax, x3, card_y, col_w, card_h, C_BG3, C_STAGE3,
          "STAGE 3: Three-Tier Inference", C_STAGE3)

# Tier 1
y_cursor = card_y + card_h - 1.02
ax.text(x3 + 0.28, y_cursor, "Tier 1: Single Crossovers",
        fontsize=13.5, fontweight='bold', color=C_STAGE3, va='top')
y_cursor -= 0.32
ax.text(x3 + 0.30, y_cursor, "• Centered Brownian bridge with\n  analytical boundary test.",
        fontsize=12.2, color='#222222', va='top')

# Subplot 3A
plot3a_top = y_cursor - 0.46
plot3a_h = 1.25
ax_sub5 = fig.add_axes([(x3 + 0.62)/fig_w, (plot3a_top - plot3a_h)/fig_h, 2.62/fig_w, plot3a_h/fig_h])
u_axis = np.linspace(0, 1, 100)
bridge_arch = 2.4 * np.sin(np.pi * u_axis) + 0.15 * np.random.randn(100)
bridge_arch[0] = 0; bridge_arch[-1] = 0
ax_sub5.plot(u_axis, bridge_arch, color='#009E73', lw=2.2)
ax_sub5.axhline(1.36, color='#CC79A7', ls='--', lw=1.5)
ax_sub5.text(0.50, 1.48, 'Significance Boundary', fontsize=9.0, color='#CC79A7', fontweight='bold', ha='center')
ax_sub5.axhline(0, color='#AAAAAA', lw=0.8)
ax_sub5.tick_params(labelsize=9.5)
ax_sub5.set_ylabel('Bridge Excursion', fontsize=10.0, fontweight='bold', labelpad=2)

# Tier 2
y_cursor = plot3a_top - plot3a_h - 0.34
ax.text(x3 + 0.28, y_cursor, "Tier 2: Short Conversions",
        fontsize=13.5, fontweight='bold', color=C_STAGE3, va='top')
y_cursor -= 0.32
ax.text(x3 + 0.30, y_cursor, "• Multi-scale wavelets and\n  drawdown isolate short tracts.",
        fontsize=12.2, color='#222222', va='top')

# Subplot 3B
plot3b_top = y_cursor - 0.46
plot3b_h = 1.25
ax_sub6 = fig.add_axes([(x3 + 0.62)/fig_w, (plot3b_top - plot3b_h)/fig_h, 2.62/fig_w, plot3b_h/fig_h])
u_axis2 = np.linspace(0, 1000, 150)
impulse = np.exp(-((u_axis2 - 400)/30)**2) - np.exp(-((u_axis2 - 600)/30)**2)
dd_curve = np.zeros_like(u_axis2)
dd_curve[u_axis2 >= 600] = 1.8 * (1.0 - np.exp(-(u_axis2[u_axis2 >= 600] - 600)/60))
ax_sub6.plot(u_axis2, impulse, color='#E69F00', lw=2.2, label='Wavelet Impulse')
ax_sub6.plot(u_axis2, dd_curve, color='#D55E00', lw=2.2, ls=':', label='Running Drawdown')
ax_sub6.tick_params(labelsize=9.5)
ax_sub6.set_ylabel('Filter Response', fontsize=10.0, fontweight='bold', labelpad=2)
ax_sub6.legend(loc='upper right', frameon=False, fontsize=8.8)

# Tier 3
y_cursor = plot3b_top - plot3b_h - 0.34
ax.text(x3 + 0.28, y_cursor, "Tier 3: Multi-Parent Mosaics",
        fontsize=13.5, fontweight='bold', color=C_STAGE3, va='top')
y_cursor -= 0.32
ax.text(x3 + 0.30, y_cursor, "• Force field + SVD mode peeling\n  stops at random noise edge.",
        fontsize=12.2, color='#222222', va='top')

# Subplot 3C
plot3c_top = y_cursor - 0.46
plot3c_h = 1.25
ax_sub7 = fig.add_axes([(x3 + 0.62)/fig_w, (plot3c_top - plot3c_h)/fig_h, 2.62/fig_w, plot3c_h/fig_h])
modes = np.arange(1, 16)
s_vals = 45.0 * np.exp(-0.35 * modes) + 3.5
ax_sub7.plot(modes, s_vals, 'o-', color='#0072B2', ms=6.0, lw=2.0)
ax_sub7.axhline(4.62, color='#D55E00', ls='--', lw=1.5)
ax_sub7.text(4, 21, 'Reticulation Modes', fontsize=9.5, color='#0072B2', fontweight='bold')
ax_sub7.text(9.5, 7.8, 'Noise Floor', fontsize=9.5, color='#D55E00', fontweight='bold')
ax_sub7.tick_params(labelsize=9.5)
ax_sub7.set_xlabel('Spatial Mode Index', fontsize=10.0, fontweight='bold', labelpad=2)
ax_sub7.set_ylabel('Singular Value', fontsize=10.0, fontweight='bold', labelpad=2)


# ==============================================================================
# CARD 4: STAGE 4 - RESOLUTION LIMITS & RECOMBINATION INVENTORY
# ==============================================================================
draw_card(ax, x4, card_y, col_w, card_h, C_BG4, C_STAGE4,
          "STAGE 4: Event Resolution", C_STAGE4)

# Section A: Resolution Limits
y_cursor = card_y + card_h - 1.02
ax.text(x4 + 0.28, y_cursor, "A. Informational Resolution",
        fontsize=13.5, fontweight='bold', color=C_STAGE4, va='top')
y_cursor -= 0.34
ax.text(x4 + 0.30, y_cursor, "• Profile likelihood measures\n  uncertainty at SNP bounds.",
        fontsize=12.5, color='#222222', va='top')

# Subplot 4A: Profile Likelihood Plateau
plot4a_top = y_cursor - 0.56
plot4a_h = 1.95
ax_sub8 = fig.add_axes([(x4 + 0.62)/fig_w, (plot4a_top - plot4a_h)/fig_h, 2.62/fig_w, plot4a_h/fig_h])
pos_grid = np.linspace(800, 950, 100)
ll_profile = -10.0 - 0.05 * np.maximum(0, 860 - pos_grid)**1.4 - 0.05 * np.maximum(0, pos_grid - 890)**1.4
ax_sub8.plot(pos_grid, ll_profile, color='#782869', lw=2.6)
ax_sub8.axvspan(860, 890, color='#E69F00', alpha=0.30)
ax_sub8.plot([875, 875], [-25, -10], color='#D55E00', ls='--', lw=1.6)
ax_sub8.scatter([875], [-10], color='#D55E00', s=60, zorder=5)
ax_sub8.text(875, -13.5, 'ML Midpoint', fontsize=10.2, ha='center', fontweight='bold', color='#D55E00',
             bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#D55E00', lw=1.2, alpha=0.95))
ax_sub8.text(875, -23.0, 'Uncertainty Plateau', fontsize=9.2, ha='center', fontweight='bold', color='#8B5A00')
ax_sub8.set_xlabel('Genomic Position (nt)', fontsize=11.0, fontweight='bold', labelpad=2)
ax_sub8.set_ylabel('Profile Likelihood', fontsize=11.0, fontweight='bold', labelpad=2)
ax_sub8.tick_params(labelsize=9.5)

# Section B: Event Inventory
y_cursor = plot4a_top - plot4a_h - 0.38
ax.text(x4 + 0.28, y_cursor, "B. Recombination Inventory",
        fontsize=13.5, fontweight='bold', color=C_STAGE4, va='top')
y_cursor -= 0.34
ax.text(x4 + 0.30, y_cursor, "• Continuous polarity separates\n  crossovers from cassettes.",
        fontsize=12.5, color='#222222', va='top')

# Subplot 4B: Resolved Mosaic Genome Track
plot4b_top = y_cursor - 0.46
plot4b_h = 1.30
ax_sub9 = fig.add_axes([(x4 + 0.40)/fig_w, (plot4b_top - plot4b_h)/fig_h, 2.90/fig_w, plot4b_h/fig_h])
ax_sub9.set_xlim(0, 3000)
ax_sub9.set_ylim(-0.4, 2.3)
ax_sub9.axis('off')

ax_sub9.text(0, 2.05, 'Resolved Mosaic Genome Track:', fontsize=10.5, fontweight='bold', color='#333333')

# Visually clear proportions: Left Home [0, 800], Donor Cassette [800, 2200] (1400 nt wide), Right Home [2200, 3000]
ax_sub9.add_patch(patches.Rectangle((0, 0.70), 800, 0.90, facecolor='#D55E00', edgecolor='white', lw=2.0))
ax_sub9.text(400, 1.15, 'Home', color='white', fontsize=11.0, ha='center', va='center', fontweight='bold')

ax_sub9.add_patch(patches.Rectangle((800, 0.70), 1400, 0.90, facecolor='#0072B2', edgecolor='white', lw=2.0))
ax_sub9.text(1500, 1.15, 'Donor', color='white', fontsize=11.0, ha='center', va='center', fontweight='bold')

ax_sub9.add_patch(patches.Rectangle((2200, 0.70), 800, 0.90, facecolor='#D55E00', edgecolor='white', lw=2.0))
ax_sub9.text(2600, 1.15, 'Home', color='white', fontsize=11.0, ha='center', va='center', fontweight='bold')

# Breakpoint lines & coordinate callouts with wide horizontal separation (1400 nt = 1.4 inches)
ax_sub9.plot([800, 800], [0.45, 1.85], color='#222222', lw=1.6)
ax_sub9.plot([2200, 2200], [0.45, 1.85], color='#222222', lw=1.6)
ax_sub9.text(800, 0.08, 'Entry: 953 nt', fontsize=9.8, ha='center', fontweight='bold', color='#222222')
ax_sub9.text(2200, 0.08, 'Exit: 1232 nt', fontsize=9.8, ha='center', fontweight='bold', color='#222222')

# Stage 4 Event Roster Callout Box (Centered, bounded, zero overflow)
badge_w = 3.20
badge_h = 0.76
badge_x = x4 + (col_w - badge_w) / 2
badge_y = card_y + 0.30

badge_patch = patches.FancyBboxPatch(
    (badge_x, badge_y), badge_w, badge_h,
    boxstyle="round,pad=0.01,rounding_size=0.12",
    facecolor='#EDE6F0', edgecolor='#782869', linewidth=1.5, zorder=2
)
ax.add_patch(badge_patch)

ax.text(badge_x + badge_w * 0.5, badge_y + badge_h * 0.68,
        "Cassette Event: [953 – 1232 nt]",
        fontsize=10.5, fontweight='bold', color='#782869', ha='center', va='center', zorder=3)
ax.text(badge_x + badge_w * 0.5, badge_y + badge_h * 0.26,
        "Donor: Clade B   |   Home: Clade C",
        fontsize=9.8, fontweight='bold', color='#4A154B', ha='center', va='center', zorder=3)


# ==============================================================================
# CONNECTING FLOW ARROWS BETWEEN CARDS
# ==============================================================================
arrow_y = card_y + card_h * 0.50
for arr_x in [x1 + col_w, x2 + col_w, x3 + col_w]:
    ax.annotate(
        "", xy=(arr_x + gap - 0.05, arrow_y), xytext=(arr_x + 0.05, arrow_y),
        arrowprops=dict(
            arrowstyle="->,head_width=0.48,head_length=0.65",
            lw=3.4, color="#555555"
        ),
        zorder=10
    )

# Save figures
fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
plt.close(fig)

print(f"[✓] Large-font, zero-overflow flowchart successfully regenerated at:\n    PDF: {pdf_path}\n    PNG: {png_path}")
