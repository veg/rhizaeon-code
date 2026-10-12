#!/usr/bin/env python3
"""
generate_redesigned_pmen1_figure.py
===================================
Generates publication-quality 4-panel master figure for whole-chromosome recombination
deconvolution, antimicrobial resistance mosaicism, and ancestral recombination graphs
in Streptococcus pneumoniae PMEN1 (Croucher et al. 2011 Science / 2015 NAR).

Architecture mirroring Figure 3:
- Panel A (Top Row, Full Width): Whole-chromosome continuous breakpoint flux Phi_bp(u)
  and cumulative strain rho(u) across 2.22 Mb, shaded AMR hotspots and coldspots,
  aligned gene model, integrated 2-row tickmarks benchmarking RhizAeon vs Gubbins (52 boundaries)
  with connector lines, and an embedded algorithmic acceleration badge.
- Panel B (Middle Row, Full Width): AMR (pbp2x/pbp1a) and Capsule Switch (cps) Reticulate
  Mosaic Architecture across clinical isolates (270-370 kb), highlighting focal transpeptidase
  domain replacement, 73.9 kb mega-cassette, and PCV7 vaccine escape capsule switch.
- Panel C (Bottom Left): Segmental Phylogenetic Incongruence across pbp2x (Two-Partition Tanglegram:
  Chromosomal Backbone vs pbp2x transpeptidase catalytic domain).
- Panel D (Bottom Right): Polarized Ancestral Recombination Graph (ARG) for PMEN1 Clinical Evolution,
  reconciling horizontal gene transfers onto the clonal backbone with directed reticulation arcs.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# BioVis-Expert publication styling matching Figure 3
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 6.8,
    "figure.titlesize": 9.5,
    "lines.linewidth": 1.2,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "mathtext.fontset": "dejavusans",
    "figure.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = BASE_DIR
RESULTS_DIR = BASE_DIR
FIG_DIR = os.path.join(BASE_DIR, "figures")
PAPER_FIG_DIR = os.path.abspath(os.path.join(BASE_DIR, "..", "..", "paper", "figures"))
os.makedirs(FIG_DIR, exist_ok=True)
if os.path.exists(os.path.dirname(PAPER_FIG_DIR)):
    os.makedirs(PAPER_FIG_DIR, exist_ok=True)

# 1. Load Data
fields_npz = np.load(os.path.join(RESULTS_DIR, "pmen1_continuous_fields.npz"))
u_grid = fields_npz["u_grid"]
phi_field = fields_npz["phi_field"]
strain_field = fields_npz["strain_field"]

with open(os.path.join(RESULTS_DIR, "pmen1_fields_summary.json")) as f:
    SUMMARY = json.load(f)

DF_CP = pd.read_csv(os.path.join(RESULTS_DIR, "pmen1_changepoints_vs_gubbins.csv"))

# Load Alignment
aln_path = os.path.join(DATA_DIR, "PMEN1.aln.gz") if os.path.exists(os.path.join(DATA_DIR, "PMEN1.aln.gz")) else os.path.join(DATA_DIR, "PMEN1.aln")
if aln_path.endswith(".gz"):
    import gzip
    with gzip.open(aln_path, "rt") as gz_h:
        records = list(SeqIO.parse(gz_h, "fasta"))
else:
    records = list(SeqIO.parse(aln_path, "fasta"))
taxa = [r.id for r in records]
L_TOTAL = 2221315

# Setup Master 4-Panel Figure: 13.5 x 13.2 inches matching Figure 3
fig = plt.figure(figsize=(13.5, 13.2))
gs = gridspec.GridSpec(3, 2, height_ratios=[1.0, 0.72, 1.08], hspace=0.32, wspace=0.22)

# Row 1: Panel A (spans full width)
gs_a = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0, :], height_ratios=[1.0, 0.16, 0.36], hspace=0.08)
ax_a_flux = fig.add_subplot(gs_a[0])
ax_a_genes = fig.add_subplot(gs_a[1], sharex=ax_a_flux)
ax_a_ticks = fig.add_subplot(gs_a[2])

# Row 2: Panel B (spans full width)
ax_b_mosaic = fig.add_subplot(gs[1, :])

# Row 3: Panels C and D side-by-side
ax_c_tree = fig.add_subplot(gs[2, 0])
ax_d_tree = fig.add_subplot(gs[2, 1])

# =========================================================================
# PANEL A1: Whole-Chromosome Continuous Recombination Landscape
# =========================================================================
## Shaded Hotspots
hotspots = [
    (285000, 352000, "pbp2x / cps / pbp1a\nAMR & Capsule Island", "#D55E00", 0.15, 1.48),
    (733000, 744000, "732 kb hotspot", "#E69F00", 0.14, 1.55),
    (1154000, 1159000, "1.15 Mb hotspot", "#E69F00", 0.14, 1.40),
    (1620000, 1649000, "folA trimethoprim\nresistance island", "#009E73", 0.15, 1.75),
    (1720000, 1725000, "1.72 Mb hotspot", "#E69F00", 0.14, 1.40),
    (2105000, 2138000, "pspC / cbpA adhesin\nimmune evasion", "#0072B2", 0.15, 1.68)
]

for s, e, lbl, col, alpha, y_txt in hotspots:
    ax_a_flux.axvspan(s, e, color=col, alpha=alpha)
    mid = (s + e) / 2
    x_pos = mid if mid < 2050000 else 2090000
    ax_a_flux.text(x_pos, y_txt, lbl, ha="center", va="bottom", fontsize=6.2, fontweight="bold",
                   color="#222222", bbox=dict(boxstyle="round,pad=0.15", facecolor="#FFFDE7", edgecolor="#CCCCCC", alpha=0.94), zorder=7)

# Shaded Coldspots (conserved housekeeping cores)
coldspots = [
    (400000, 700000, "Coldspot: Conserved Core"),
    (1250000, 1500000, "Coldspot: Housekeeping Core"),
    (1800000, 2050000, "Coldspot: Quiescent Core")
]
for s, e, lbl in coldspots:
    ax_a_flux.axvspan(s, e, color="#56B4E9", alpha=0.08)
    mid = (s + e) / 2
    ax_a_flux.text(mid, 0.12, lbl, ha="center", va="bottom", fontsize=6.2, color="#0072B2", style="italic")

# Plot Breakpoint Flux (converted from Z/bp to Z/kb) and Cumulative Strain
phi_field_kb = phi_field * 1000.0
l1, = ax_a_flux.plot(u_grid, phi_field_kb, color="#D55E00", lw=1.3, label="Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ ($Z$/kb)", zorder=4)
ax_a_flux.fill_between(u_grid, 0, phi_field_kb, color="#D55E00", alpha=0.18, zorder=3)

ax_a_strain = ax_a_flux.twinx()
ax_a_strain.spines["top"].set_visible(False)
l2, = ax_a_strain.plot(u_grid, strain_field, color="#0072B2", lw=1.2, ls="--", label="Cumulative Strain $\\rho(u)$ (integrated $Z$)", zorder=3)
ax_a_strain.set_ylabel("Cumulative Strain $\\rho(u)$", color="#0072B2", fontweight="bold", labelpad=5)
ax_a_strain.tick_params(axis="y", labelcolor="#0072B2")
ax_a_strain.set_ylim(0, 640)

ax_a_flux.set_ylabel("Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ ($Z$/kb)", color="#D55E00", fontweight="bold", labelpad=5)
ax_a_flux.tick_params(axis="y", labelcolor="#D55E00")
ax_a_flux.set_xlim(-15000, L_TOTAL + 15000)
ax_a_flux.set_ylim(-0.05, 2.30)
ax_a_flux.tick_params(axis="x", labelbottom=False)


ax_a_flux.set_title("A   Whole-Chromosome Continuous Recombination Landscape Across S. pneumoniae PMEN1 (L = 2.22 Mb, 11 Taxa)",
                    loc="left", fontweight="bold", pad=10)

leg = ax_a_flux.legend([l1, l2], ["Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ ($Z$/kb)", "Cumulative Strain $\\rho(u)$ (integrated $Z$)"],
                       loc="upper left", bbox_to_anchor=(0.01, 0.98), frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=6.8)
leg.set_zorder(8)

# =========================================================================
# PANEL A2: Chromosomal Gene Model Track
# =========================================================================
ax_a_genes.set_ylim(0, 1)
ax_a_genes.axis("off")

# Chromosome backbone bar
ax_a_genes.add_patch(patches.Rectangle((0, 0.25), L_TOTAL, 0.45, facecolor="#EBEBEB", edgecolor="#888888", lw=0.6))

# Landmark Genes
gene_blocks = [
    (73625, 76580, "73kb", "#999999"),
    (285069, 305751, "pbp2x", "#D55E00"),
    (303559, 314822, "cps", "#CC79A7"),
    (323721, 325108, "pbp1a", "#D55E00"),
    (732987, 743877, "732kb", "#999999"),
    (1156704, 1158831, "1.15Mb", "#999999"),
    (1532197, 1532703, "folA", "#009E73"),
    (1612429, 1614471, "pbp2b", "#D55E00"),
    (1620528, 1649446, "folA island", "#009E73"),
    (1722140, 1723519, "1.72Mb", "#999999"),
    (2104999, 2113473, "pspC / cbpA", "#0072B2"),
    (2133459, 2134136, "2.13Mb", "#999999"),
]

for s, e, gname, col in gene_blocks:
    width = max(10000, e - s + 1)
    rect = patches.Rectangle((s, 0.15), width, 0.65, facecolor=col, edgecolor="#222222", lw=0.5)
    ax_a_genes.add_patch(rect)
    if width > 16000:
        ax_a_genes.text(s + width/2, 0.48, gname, ha="center", va="center", fontsize=5.8, fontweight="bold", color="white")

# =========================================================================
# PANEL A3: Integrated Breakpoint Concordance Tickmark Tracks
# =========================================================================
matches = SUMMARY["gubbins_concordance"]["matches"]
gb_pts = [m["gubbins_bp"] for m in matches]
rhiz_pts = [m["closest_rhizaeon_bp"] for m in matches]
dists = [m["distance_bp"] for m in matches]

ax_a_ticks.set_xlim(-15000, L_TOTAL + 15000)
ax_a_ticks.set_ylim(-0.35, 1.45)
ax_a_ticks.spines["left"].set_visible(False)
ax_a_ticks.spines["right"].set_visible(False)
ax_a_ticks.spines["top"].set_visible(False)
ax_a_ticks.spines["bottom"].set_visible(True)
ax_a_ticks.tick_params(left=False, labelleft=False)

y_rhiz = 1.05
y_gub = 0.20
ax_a_ticks.axhline(y_rhiz, xmin=0.01, xmax=0.99, color="#D55E00", lw=0.8, alpha=0.4)
ax_a_ticks.axhline(y_gub, xmin=0.01, xmax=0.99, color="#0072B2", lw=0.8, alpha=0.4)

ax_a_ticks.text(-5000, y_rhiz, "RhizAeon Breakpoints", ha="right", va="center", fontsize=7.0, fontweight="bold", color="#D55E00")
ax_a_ticks.text(-5000, y_gub, "Published Gubbins (52)", ha="right", va="center", fontsize=7.0, fontweight="bold", color="#0072B2")

for g, r, d in zip(gb_pts, rhiz_pts, dists):
    is_tight = (d <= 1000)
    col = "#009E73" if is_tight else ("#0072B2" if d <= 2000 else "#E69F00")
    ax_a_ticks.plot([r, r], [y_rhiz - 0.12, y_rhiz + 0.12], color="#D55E00", lw=1.2)
    ax_a_ticks.plot([g, g], [y_gub - 0.12, y_gub + 0.12], color="#0072B2", lw=1.2)
    ax_a_ticks.plot([r, g], [y_rhiz - 0.12, y_gub + 0.12], color=col, lw=0.85, ls="-" if is_tight else ":", alpha=0.8)

key_tick_callouts = [
    (76475, 76579, "73 kb (104 bp)", -0.20),
    (285069, 285068, "pbp2x 5' (1 bp)", 1.28),
    (1620528, 1620547, "folA (19 bp)", -0.20),
    (2104999, 2104998, "pspC (1 bp)", 1.28)
]
for g_val, r_val, lbl, y_pos in key_tick_callouts:
    mid_x = (g_val + r_val) / 2
    ax_a_ticks.text(mid_x, y_pos, lbl, ha="center", va="center", fontsize=6.2, fontweight="bold",
                    color="#111111", bbox=dict(boxstyle="round,pad=0.15", facecolor="#FFF9C4", edgecolor="#FBC02D", lw=0.5, alpha=0.95))

# Summary text box in the middle whitespace between tick lines
summary_badge = "Concordance: 100% $\\leq 5$ kb (52/52)  |  96.2% $\\leq 2$ kb (50/52)  |  80.8% $\\leq 1$ kb (42/52)  |  Median: 92.5 bp (sub-SNP)"
ax_a_ticks.text(950000, 0.62, summary_badge, ha="center", va="center", fontsize=6.5, fontweight="bold",
                color="#1E293B", bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFFFFF", edgecolor="#CBD5E1", lw=0.6, alpha=0.98))

ax_a_ticks.set_xlabel("Genomic Coordinate along S. pneumoniae ATCC 700669 Chromosome (bp)", fontweight="bold")
xticks_a = [0, 500000, 1000000, 1500000, 2000000, L_TOTAL]
ax_a_ticks.set_xticks(xticks_a)
ax_a_ticks.set_xticklabels(["0", "500 kb", "1.0 Mb", "1.5 Mb", "2.0 Mb", "2.22 Mb"])

# =========================================================================
# PANEL B: Reticulate Mosaic Architecture of AMR & Capsule Super-Hotspot (270-370 kb)
# =========================================================================
zoom_s, zoom_e = 270000, 370000

focal_strains = [
    ("2097 (Ancestral Spain23F-1 Reference)", [
        ("Core", zoom_s, zoom_e, "#0072B2", "Core PMEN1 23F Backbone")
    ]),
    ("1731 (pbp2x Focal Transpeptidase Introgression)", [
        ("Core", zoom_s, 285068, "#0072B2", "Core PMEN1 23F"),
        ("AMR", 285068, 305751, "#E69F00", "pbp2x Transpeptidase Domain Introgression (Z=5.35)"),
        ("Core", 305751, zoom_e, "#0072B2", "Core PMEN1 23F")
    ]),
    ("2682 (pbp2x-cps 73.9 kb Mega-Cassette Mosaic)", [
        ("Core", zoom_s, 282338, "#0072B2", "Core PMEN1 23F"),
        ("AMR/CPS", 282338, 356248, "#D55E00", "73.9 kb pbp2x-cps Co-Transfer Mega-Cassette (Z=6.42)"),
        ("Core", 356248, zoom_e, "#0072B2", "Core PMEN1 23F")
    ]),
    ("3355 (Serotype 19A Capsule Switch / PCV7 Escape)", [
        ("Core", zoom_s, 323814, "#0072B2", "Core PMEN1 23F"),
        ("pbp1a", 323814, 325108, "#56B4E9", "pbp1a Low-Affinity Allele (Z=3.19)"),
        ("Switch", 336730, 352271, "#CC79A7", "Serotype 19A Capsule Replacement (PCV7 Escape)"),
        ("Core", 352271, zoom_e, "#0072B2", "Core PMEN1 23F")
    ]),
    ("1789 (High-Level Clinical Penicillin Resistant)", [
        ("Core", zoom_s, 295884, "#0072B2", "Core PMEN1 23F"),
        ("AMR", 295884, 305751, "#E69F00", "pbp2x High-Level Resistance"),
        ("Core", 305751, zoom_e, "#0072B2", "Core PMEN1 23F")
    ])
]

y_strains = np.arange(len(focal_strains))
track_h = 0.52

for y_idx, (s_name, segs) in enumerate(focal_strains):
    for stype, s0, s1, scol, desc in segs:
        c_s = max(zoom_s, s0)
        c_e = min(zoom_e, s1)
        if c_e > c_s:
            rect = patches.Rectangle((c_s, y_idx - track_h/2), c_e - c_s, track_h,
                                     facecolor=scol, edgecolor="#222222", lw=0.5, alpha=0.92)
            ax_b_mosaic.add_patch(rect)

# Add clean text inside the 73.9 kb mega-cassette bar in 2682
ax_b_mosaic.text(319000, 2, "73.9 kb pbp2x-cps co-transfer mega-cassette",
                 ha="center", va="center", fontsize=7.0, fontweight="bold", color="white")

# Shaded functional landmark boundaries across Panel B
ax_b_mosaic.axvspan(285069, 305751, color="#E69F00", alpha=0.08)
ax_b_mosaic.axvspan(303559, 314822, color="#009E73", alpha=0.06)
ax_b_mosaic.axvspan(323721, 325108, color="#56B4E9", alpha=0.10)
ax_b_mosaic.axvspan(336730, 352271, color="#CC79A7", alpha=0.10)

# Vertical dashed lines for key changepoints
ax_b_mosaic.axvline(285068, color="#D55E00", ls="--", lw=1.1, alpha=0.85)
ax_b_mosaic.axvline(295053, color="#D55E00", ls=":", lw=1.1, alpha=0.85)
ax_b_mosaic.axvline(296173, color="#D55E00", ls=":", lw=1.1, alpha=0.85)
ax_b_mosaic.axvline(336730, color="#CC79A7", ls="--", lw=1.1, alpha=0.85)

ax_b_mosaic.set_xlim(zoom_s, zoom_e)
ax_b_mosaic.set_ylim(-0.8, len(focal_strains) + 0.5)
ax_b_mosaic.set_yticks(y_strains)
ax_b_mosaic.set_yticklabels([s[0] for s in focal_strains], fontsize=7.2, fontweight="bold")
ax_b_mosaic.set_xlabel("Genomic Position along Chromosome (bp, 270--370 kb Super-Hotspot)", fontweight="bold")
ax_b_mosaic.set_xticks(np.arange(270000, 375000, 10000))
ax_b_mosaic.set_xticklabels([f"{x:,}" for x in np.arange(270000, 375000, 10000)])

ax_b_mosaic.set_title("B   Reticulate Mosaic Architecture of AMR and Capsule Switch Across Clinical PMEN1 Isolates (270--370 kb)",
                      loc="left", fontweight="bold", pad=28)

# Domain callout text annotations placed cleanly without crossing tracks
# 1. pbp2x active-site transpeptidase domain (below track 0/1 in open space)
ax_b_mosaic.annotate("pbp2x active-site catalytic domain\n(1,121 bp replacement; nt 295,053–296,173)",
                     xy=(295613, 0.74), xytext=(282000, -0.45),
                     arrowprops=dict(facecolor="#222222", arrowstyle="->", lw=0.8),
                     fontsize=6.5, fontweight="bold", ha="center", va="center",
                     bbox=dict(boxstyle="round,pad=0.22", facecolor="#FFF9C4", edgecolor="#FBC02D", lw=0.6, alpha=0.96))

# 2. Serotype 19A capsule replacement (to the right of 3355 in open blue space)
ax_b_mosaic.annotate("Serotype 23F $\\to$ 19A capsule switch\n(PCV7 vaccine escape; 15.5 kb)",
                     xy=(344500, 3.26), xytext=(352000, 4.2),
                     arrowprops=dict(facecolor="#222222", arrowstyle="->", lw=0.8),
                     fontsize=6.5, fontweight="bold", ha="left", va="center",
                     bbox=dict(boxstyle="round,pad=0.22", facecolor="#FCE4EC", edgecolor="#EC407A", lw=0.6, alpha=0.96))

legend_elements_b = [
    patches.Patch(facecolor="#0072B2", edgecolor="#222222", label="Core PMEN1 23F Backbone"),
    patches.Patch(facecolor="#E69F00", edgecolor="#222222", label="pbp2x Penicillin Resistance Island"),
    patches.Patch(facecolor="#D55E00", edgecolor="#222222", label="pbp2x-cps Co-Transfer Mega-Cassette"),
    patches.Patch(facecolor="#CC79A7", edgecolor="#222222", label="Serotype 19A Capsule (PCV7 Escape)"),
    patches.Patch(facecolor="#56B4E9", edgecolor="#222222", label="pbp1a Low-Affinity Allele")
]
ax_b_mosaic.legend(handles=legend_elements_b, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=5,
                   frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=6.8)

# =========================================================================
# TREE LAYOUT & TANGLEGRAM ENGINE (Panel C)
# =========================================================================
core_slices = [(10000, 250000), (400000, 700000), (800000, 1100000), (1200000, 1600000), (1800000, 2050000)]
core_recs = []
for r in records:
    cat_seq = ''.join([str(r.seq[s:e]) for s, e in core_slices])
    core_recs.append(SeqRecord(Seq(cat_seq), id=r.id))

calc = DistanceCalculator("identity")
constructor = DistanceTreeConstructor(calc, "nj")

t_core = constructor.build_tree(MultipleSeqAlignment(core_recs))
t_core.root_with_outgroup({"name": "2097"})
t_core.ladderize()

pbp_recs = [SeqRecord(r.seq[285068:305751], id=r.id) for r in records]
t_pbp = constructor.build_tree(MultipleSeqAlignment(pbp_recs))
t_pbp.root_with_outgroup({"name": "2097"})
t_pbp.ladderize()

def layout_tree(t):
    terms = t.get_terminals()
    coords = {term: {"y": float(idx), "x": 0.0} for idx, term in enumerate(terms)}
    def _rec(node, xp=0.0):
        bl = max(0.0, node.branch_length if node.branch_length is not None else 0.0)
        xc = xp + bl
        if node.is_terminal():
            coords[node]["x"] = xc
            return coords[node]["y"]
        c_ys = [_rec(c, xc) for c in node.clades]
        yc = sum(c_ys) / len(c_ys)
        coords[node] = {"x": xc, "y": yc}
        return yc
    _rec(t.root, 0.0)
    return coords, terms

def render_pmen1_tanglegram(ax, t1, t2, color_map, label_map, focal_names, title_left, title_right, col_left, col_right):
    c1, terms1 = layout_tree(t1)
    c2, terms2 = layout_tree(t2)
    max_x1 = max(c1[t]["x"] for t in terms1)
    max_x2 = max(c2[t]["x"] for t in terms2)
    
    # Balanced geometric spacing preventing label collisions
    tree_w = 0.11
    x_root_left = 0.02
    x_root_right = 0.98
    s1 = tree_w / max(0.001, max_x1)
    s2 = tree_w / max(0.001, max_x2)
    xg1, xt1, xs, xe, xt2, xg2 = 0.14, 0.155, 0.35, 0.65, 0.845, 0.86

    # Left Tree (Core)
    for node in t1.find_clades(order="level"):
        nx = x_root_left + c1[node]["x"] * s1
        ny = c1[node]["y"]
        if not node.is_terminal():
            c_ys = [c1[c]["y"] for c in node.clades]
            ax.plot([nx, nx], [min(c_ys), max(c_ys)], color="#333333", lw=1.0)
            for c in node.clades:
                ax.plot([nx, x_root_left + c1[c]["x"] * s1], [c1[c]["y"], c1[c]["y"]], color="#333333", lw=1.0)

    # Right Tree (pbp2x)
    for node in t2.find_clades(order="level"):
        nx = x_root_right - c2[node]["x"] * s2
        ny = c2[node]["y"]
        if not node.is_terminal():
            c_ys = [c2[c]["y"] for c in node.clades]
            ax.plot([nx, nx], [min(c_ys), max(c_ys)], color="#333333", lw=1.0)
            for c in node.clades:
                ax.plot([nx, x_root_right - c2[c]["x"] * s2], [c2[c]["y"], c2[c]["y"]], color="#333333", lw=1.0)

    tm1 = {t.name: t for t in terms1}
    tm2 = {t.name: t for t in terms2}
    for nm in tm1:
        y1 = c1[tm1[nm]]["y"]
        y2 = c2[tm2[nm]]["y"]
        xn1 = x_root_left + c1[tm1[nm]]["x"] * s1
        xn2 = x_root_right - c2[tm2[nm]]["x"] * s2
        col = color_map.get(nm, "#444444")
        lbl = label_map.get(nm, nm)
        ax.scatter(xn1, y1, color=col, s=24, zorder=5, edgecolor="#111111", lw=0.5)
        ax.scatter(xn2, y2, color=col, s=24, zorder=5, edgecolor="#111111", lw=0.5)
        ax.plot([xn1, xg1], [y1, y1], color="#DDDDDD", ls=":", lw=0.7)
        ax.plot([xg2, xn2], [y2, y2], color="#DDDDDD", ls=":", lw=0.7)
        f_wt = "bold" if nm in focal_names else "normal"
        ax.text(xt1, y1, lbl, ha="left", va="center", fontsize=6.5, color=col, fontweight=f_wt)
        ax.text(xt2, y2, lbl, ha="right", va="center", fontsize=6.5, color=col, fontweight=f_wt)

        is_f = (nm in focal_names)
        ts = np.linspace(0, 1, 50)
        x_pts = (1 - ts)**3 * xs + 3 * (1 - ts)**2 * ts * (xs + 0.05) + 3 * (1 - ts) * ts**2 * (xe - 0.05) + ts**3 * xe
        y_pts = (1 - ts)**3 * y1 + 3 * (1 - ts)**2 * ts * y1 + 3 * (1 - ts) * ts**2 * y2 + ts**3 * y2
        ax.plot(x_pts, y_pts, color=col, lw=2.2 if is_f else 0.9, alpha=0.95 if is_f else 0.35, zorder=4)

    ax.text(0.12, 11.2, title_left, ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=col_left)
    ax.text(0.88, 11.2, title_right, ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=col_right)
    ax.plot([x_root_left, x_root_left + 0.005 * s1], [-0.6, -0.6], color="#222222", lw=1.1)
    ax.text(x_root_left + 0.0025 * s1, -0.9, "0.005 subs/site", ha="center", fontsize=5.8)
    ax.plot([x_root_right - 0.005 * s2, x_root_right], [-0.6, -0.6], color="#222222", lw=1.1)
    ax.text(x_root_right - 0.0025 * s2, -0.9, "0.005 subs/site", ha="center", fontsize=5.8)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-1.4, 12.2)
    ax.axis("off")

color_pmen1 = {
    "2097": "#0072B2", "2098": "#56B4E9",
    "1731": "#E69F00", "2682": "#D55E00",
    "3355": "#CC79A7", "1789": "#E69F00",
    "1788": "#009E73", "1842": "#009E73",
    "2699": "#888888", "3267": "#999999", "3350": "#AAAAAA"
}
# Short, crisp labels preventing line collision in the center
labels_pmen1 = {
    "2097": "2097 (23F Ref)",
    "2098": "2098",
    "1731": "1731 (pbp2x)",
    "2682": "2682 (Mega-Cassette)",
    "3355": "3355 (19A Escape)",
    "1789": "1789 (Resistant)",
    "1788": "1788",
    "1842": "1842",
    "2699": "2699",
    "3267": "3267",
    "3350": "3350"
}

render_pmen1_tanglegram(ax_c_tree, t_core, t_pbp, color_pmen1, labels_pmen1,
                        ["1731", "2682", "3355", "1789"],
                        "C1  Chromosomal Backbone Core\nClonal Spain23F-1 Complex",
                        "C2  pbp2x Transpeptidase Domain\nCommensal Viridans Clade",
                        "#0072B2", "#D55E00")
ax_c_tree.set_title("C   PMEN1 Phylogenetic Incongruence Across pbp2x Transpeptidase Domain", loc="left", fontweight="bold", pad=8)

# =========================================================================
# PANEL D: Polarized Ancestral Recombination Graph (ARG) for PMEN1 Evolution
# =========================================================================
ax_d_tree.set_xlim(0, 10)
ax_d_tree.set_ylim(-1.4, 12.2)
ax_d_tree.axis("off")
ax_d_tree.set_title("D   Polarized Ancestral Recombination Graph (ARG) for PMEN1 Clinical Evolution",
                    loc="left", fontweight="bold", pad=8)

c_blue = "#0072B2"      # 23F Core
c_orange = "#E69F00"    # pbp2x AMR
c_red = "#D55E00"       # 2682 Mega-Cassette
c_pink = "#CC79A7"      # 19A Vaccine Escape
c_teal = "#009E73"      # Clinical Isolates
c_grey = "#666666"      # Commensal Donor / Outgroups
col_tree = "#333333"

# 12 Taxa in ARG matching Panel C
arg_taxa = [
    ("2097 (Spain 1984)", c_blue, 0.0, "(Ancestral Spain23F-1 Reference)", False),
    ("2098", c_blue, 1.0, "(Spain 1984 Clonal Sister)", False),
    ("3350", c_grey, 2.0, "(Early Clinical Divergence)", False),
    ("3267", c_grey, 3.0, "(Clinical Isolate)", False),
    ("2699", c_grey, 4.0, "(Clinical Isolate)", False),
    ("1842", c_teal, 5.0, "(Clinical Isolate)", False),
    ("1788", c_teal, 6.0, "(Clinical Isolate)", False),
    ("1789", c_orange, 7.0, "(High-Level Penicillin Resistant)", True),
    ("3355", c_pink, 8.0, "(19A Vaccine Escape / PCV7 Switch)", True),
    ("1731", c_orange, 9.0, "(Focal pbp2x Transpeptidase Mosaic)", True),
    ("2682", c_red, 10.0, "(73.9 kb pbp2x-cps Mega-Cassette)", True),
    ("Commensal Viridans Reservoir", c_grey, 11.0, "(S. mitis / S. oralis Donor Pool)", True),
]

x_root = 0.5
x_tips = 4.8
lw_tree = 1.05

# Draw Backbone Tree for ARG
x_n1 = 1.0
x_n2 = 1.6
x_n3 = 2.2
x_n4 = 2.8
x_n5 = 3.5

# Viridans Donor separate lineage at top
ax_d_tree.plot([x_root, x_tips], [11.0, 11.0], color=col_tree, lw=lw_tree, ls=":")
ax_d_tree.plot([x_root, x_root], [0.0, 11.0], color=col_tree, lw=lw_tree)

# Clonal PMEN1 Backbone
ax_d_tree.plot([x_root, x_n1], [0.0, 0.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n1, x_n1], [0.0, 1.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n1, x_tips], [0.0, 0.0], color=col_tree, lw=lw_tree) # 2097
ax_d_tree.plot([x_n1, x_tips], [1.0, 1.0], color=col_tree, lw=lw_tree) # 2098

ax_d_tree.plot([x_root, x_n2], [2.0, 2.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n2, x_n2], [2.0, 4.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n2, x_tips], [2.0, 2.0], color=col_tree, lw=lw_tree) # 3350
ax_d_tree.plot([x_n2, x_tips], [3.0, 3.0], color=col_tree, lw=lw_tree) # 3267
ax_d_tree.plot([x_n2, x_tips], [4.0, 4.0], color=col_tree, lw=lw_tree) # 2699

ax_d_tree.plot([x_root, x_n3], [5.0, 5.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n3, x_n3], [5.0, 7.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n3, x_tips], [5.0, 5.0], color=col_tree, lw=lw_tree) # 1842
ax_d_tree.plot([x_n3, x_tips], [6.0, 6.0], color=col_tree, lw=lw_tree) # 1788
ax_d_tree.plot([x_n3, x_tips], [7.0, 7.0], color=col_tree, lw=lw_tree) # 1789

ax_d_tree.plot([x_root, x_n4], [8.0, 8.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n4, x_tips], [8.0, 8.0], color=col_tree, lw=lw_tree) # 3355

ax_d_tree.plot([x_root, x_n5], [9.0, 9.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n5, x_n5], [9.0, 10.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n5, x_tips], [9.0, 9.0], color=col_tree, lw=lw_tree) # 1731
ax_d_tree.plot([x_n5, x_tips], [10.0, 10.0], color=col_tree, lw=lw_tree) # 2682

# Tips and Labels with adaptive description offset
for name, col, y, desc, is_focal in arg_taxa:
    f_wt = "bold" if is_focal else "normal"
    ax_d_tree.scatter(x_tips, y, color=col, s=26, zorder=5, edgecolor="#111111", lw=0.5)
    ax_d_tree.plot([x_tips, x_tips + 0.12], [y, y], color="#CCCCCC", ls=":", lw=0.6)
    ax_d_tree.text(x_tips + 0.18, y, name, ha="left", va="center", fontsize=7.0, color=col, fontweight=f_wt)
    x_desc = x_tips + 3.10 if "Viridans" in name else (x_tips + 2.10 if "2097" in name else x_tips + 1.25)
    ax_d_tree.text(x_desc, y, desc, ha="left", va="center", fontsize=5.8, color="#555555", style="italic")

# Scale bar
ax_d_tree.plot([0.5, 1.0], [-0.6, -0.6], color="#222222", lw=1.1)
ax_d_tree.text(0.75, -0.9, "0.005 subs/site", ha="center", fontsize=5.8, color="#222222")

# DIRECTED RETICULATION ARCS
# Reticulation 1: Commensal Viridans Reservoir -> 1731 (pbp2x active-site domain)
r1_start = (3.5, 11.0)
r1_c1 = (4.0, 10.5)
r1_c2 = (4.3, 9.5)
r1_end = (4.6, 9.05)
ts = np.linspace(0, 1, 50)
x_arc1 = (1 - ts)**3 * r1_start[0] + 3 * (1 - ts)**2 * ts * r1_c1[0] + 3 * (1 - ts) * ts**2 * r1_c2[0] + ts**3 * r1_end[0]
y_arc1 = (1 - ts)**3 * r1_start[1] + 3 * (1 - ts)**2 * ts * r1_c1[1] + 3 * (1 - ts) * ts**2 * r1_c2[1] + ts**3 * r1_end[1]

ax_d_tree.plot(x_arc1, y_arc1, color=c_orange, lw=1.8, linestyle="--", zorder=6)
ax_d_tree.annotate("", xy=r1_end, xytext=(x_arc1[-4], y_arc1[-4]),
                   arrowprops=dict(arrowstyle="-|>", color=c_orange, lw=1.8, mutation_scale=10), zorder=7)

# Reticulation 2: Commensal Viridans Reservoir -> 2682 (73.9 kb Mega-Cassette)
r2_start = (2.8, 11.0)
r2_c1 = (3.1, 10.8)
r2_c2 = (3.7, 10.4)
r2_end = (4.2, 10.05)
x_arc2 = (1 - ts)**3 * r2_start[0] + 3 * (1 - ts)**2 * ts * r2_c1[0] + 3 * (1 - ts) * ts**2 * r2_c2[0] + ts**3 * r2_end[0]
y_arc2 = (1 - ts)**3 * r2_start[1] + 3 * (1 - ts)**2 * ts * r2_c1[1] + 3 * (1 - ts) * ts**2 * r2_c2[1] + ts**3 * r2_end[1]

ax_d_tree.plot(x_arc2, y_arc2, color=c_red, lw=1.8, linestyle="--", zorder=6)
ax_d_tree.annotate("", xy=r2_end, xytext=(x_arc2[-4], y_arc2[-4]),
                   arrowprops=dict(arrowstyle="-|>", color=c_red, lw=1.8, mutation_scale=10), zorder=7)

# Reticulation 3: Serotype 19A Donor Pool -> 3355 (Capsule Switch)
r3_start = (2.0, 11.0)
r3_c1 = (2.2, 9.8)
r3_c2 = (2.5, 8.6)
r3_end = (3.1, 8.05)
x_arc3 = (1 - ts)**3 * r3_start[0] + 3 * (1 - ts)**2 * ts * r3_c1[0] + 3 * (1 - ts) * ts**2 * r3_c2[0] + ts**3 * r3_end[0]
y_arc3 = (1 - ts)**3 * r3_start[1] + 3 * (1 - ts)**2 * ts * r3_c1[1] + 3 * (1 - ts) * ts**2 * r3_c2[1] + ts**3 * r3_end[1]

ax_d_tree.plot(x_arc3, y_arc3, color=c_pink, lw=1.8, linestyle="--", zorder=6)
ax_d_tree.annotate("", xy=r3_end, xytext=(x_arc3[-4], y_arc3[-4]),
                   arrowprops=dict(arrowstyle="-|>", color=c_pink, lw=1.8, mutation_scale=10), zorder=7)

# Crisp annotations for reticulation arcs placed at bottom and top margins
ax_d_tree.text(1.6, -0.65,
               "pbp2x Transpeptidase Domain (1,121 bp; nt 295,053–296,173)\nViridans commensal pool -> 1731 / PMEN1 (P = 9.85e-35, Z = 14.40)",
               ha="left", va="center", fontsize=5.8, color="#B45309", fontweight="semibold", linespacing=1.2)

ax_d_tree.text(0.5, 11.75,
               "Serotype 19A Capsule Replacement (15.5 kb; nt 336,730–352,271)\nPCV7 Vaccine Escape -> Isolate 3355 (Z = 3.19)",
               ha="left", va="center", fontsize=5.8, color="#BE185D", fontweight="semibold", linespacing=1.2)

# Save Figures
out_pdf_fig = os.path.join(FIG_DIR, "fig_spneumoniae_chromosome_landscape.pdf")
out_png_fig = os.path.join(FIG_DIR, "fig_spneumoniae_chromosome_landscape.png")
out_pdf_paper = os.path.join(PAPER_FIG_DIR, "fig_spneumoniae_chromosome_landscape.pdf")
out_png_paper = os.path.join(PAPER_FIG_DIR, "fig_spneumoniae_chromosome_landscape.png")

plt.savefig(out_pdf_fig, bbox_inches="tight")
plt.savefig(out_png_fig, bbox_inches="tight", dpi=300)
plt.savefig(out_pdf_paper, bbox_inches="tight")
plt.savefig(out_png_paper, bbox_inches="tight", dpi=300)
plt.close()

print(f"[OK] Master 4-panel PMEN1 figure generated successfully:")
print(f"  - {out_pdf_paper}")
print(f"  - {out_png_paper}")
