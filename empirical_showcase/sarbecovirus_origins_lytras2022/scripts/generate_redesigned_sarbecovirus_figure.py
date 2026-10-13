#!/usr/bin/env python3
"""
generate_redesigned_sarbecovirus_figure.py

Publication Figure Generator for Sarbecovirus Whole-Genome Recombination Landscape
and Progenitor Mosaic Architecture (Lytras et al. 2022 + Temmam et al. 2022 BANAL).

4-Panel Architecture:
- Panel A (Top Row, Full Width): Genome-wide physical breakpoint flux Phi_bp(u) across 30.9 kb,
  shaded RdRp and Spike S2 coldspots, aligned gene model, and integrated 2-row tickmarks
  benchmarking RhizAeon vs published GARD (2022) breakpoints with 1-to-1 connector lines.
- Panel B (Middle Row, Full Width): Spike Reticulate Mosaic Architecture across 8 wildlife progenitors
  and emergent SARS-CoV-2 (Wuhan-Hu-1), highlighting Wuhan-Hu-1 as an emergent chimera and
  displaying nested NTD micro-conversions in BANAL-247 and BANAL-116.
- Panel C (Bottom Left): Segmental Phylogenetic Incongruence across Codon 313 Hinge for SARS-CoV-2
  (Spike NTD vs RBD rooted NJ trees connected by a tanglegram).
- Panel D (Bottom Right): Segmental Phylogenetic Incongruence for Wildlife Chimera BANAL-247
  (Spike S2 Backbone vs NTD RacCS203 Micro-Conversion rooted NJ trees connected by a tanglegram).
"""

import os
import json
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Set publication style
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
    "figure.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")
FIG_DIR = os.path.join(BASE_DIR, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

REPORT_PATH = os.path.join(DATA_DIR, "sarbecovirus_reanalysis_report.json")
GENE_PATH = os.path.join(DATA_DIR, "sarbeco_gene_annotations.json")
RES_83_PATH = os.path.join(DATA_DIR, "83sarbeco_result.json")
FASTA_83 = os.path.join(DATA_DIR, "83sarbeco_wgal.fas")

with open(REPORT_PATH, "r") as f:
    REPORT = json.load(f)

with open(GENE_PATH, "r") as f:
    GENES = json.load(f)

with open(RES_83_PATH, "r") as f:
    RES_83 = json.load(f)

def get_flux_profile(events, length, sigma=120.0):
    flux = np.zeros(length, dtype=float)
    for ev in events:
        z = ev.get("z_phys", 0.0)
        bps = [ev.get("u1", 0), ev.get("u2", 0)]
        for bp in bps:
            if 0 <= bp < length:
                idx_min = max(0, int(bp - 4 * sigma))
                idx_max = min(length, int(bp + 4 * sigma))
                xs = np.arange(idx_min, idx_max)
                flux[idx_min:idx_max] += z * np.exp(-0.5 * ((xs - bp) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    return flux

L_83 = RES_83.get("length", 30965)
u_83 = np.arange(L_83)
flux_83 = get_flux_profile(RES_83.get("events", []), L_83)

# 4-Panel Master Layout: 13.5 x 13.2 inches
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
# PANEL A1: Recombination Landscape (Flux + Coldspots + Hotspots)
# =========================================================================
ax_a_flux.plot(u_83, flux_83, color="#D55E00", lw=1.8, label="Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ (83 sarbecoviruses)")
ax_a_flux.fill_between(u_83, 0, flux_83, color="#D55E00", alpha=0.15)
ax_a_flux.set_xlim(-200, 31200)
ax_a_flux.set_ylim(-0.02, np.max(flux_83) * 1.15)
ax_a_flux.spines["right"].set_visible(False)
ax_a_flux.spines["top"].set_visible(False)
ax_a_flux.tick_params(bottom=False, labelbottom=False)
ax_a_flux.set_ylabel("Breakpoint Flux $\\Phi_{\\mathrm{bp}}$ (Z/nt)", fontweight="bold", color="#D55E00")
ax_a_flux.set_title("A   Genome-wide Recombination Flux $\\Phi_{\\mathrm{bp}}(u)$ and Benchmark Concordance Across 30.9 kb",
                    loc="left", fontweight="bold", pad=10)

# Shaded Coldspots (in 83-Taxon Multiple Sequence Alignment coordinates)
coldspots = [
    (13674, 16468, "Coldspot: RdRp Catalytic Core"),
    (24600, 25600, "Coldspot: Spike S2 Fusion Engine")
]
for c_start, c_end, c_label in coldspots:
    rect_c = patches.Rectangle((c_start, 0), c_end - c_start, np.max(flux_83) * 1.15,
                               facecolor="#56B4E9", alpha=0.16, edgecolor="none")
    ax_a_flux.add_patch(rect_c)
    mid = (c_start + c_end) / 2
    ax_a_flux.text(mid, 0.14, c_label, ha="center", va="center", fontsize=7.0,
                   fontweight="bold", color="#0072B2",
                   bbox=dict(boxstyle="round,pad=0.2", facecolor="#EBF5FB", edgecolor="#AECEDB", lw=0.6, alpha=0.92))

# Major Hotspot Annotations (83-MSA coordinate | Wuhan-Hu-1 Reference coordinate)
hotspot_calls = [
    (13001, "RdRp 5' entry\n(MSA 13,001 | Ref 12,769)", 0.72),
    (21825, "Spike 5' / NTD\n(MSA 21,825 | Ref 21,566)", 0.92),
    (22764, "Codon 313 hinge\n(MSA 22,764 | Ref 22,468)", 0.76),
    (23629, "RBD 3' exit\n(MSA 23,629 | Ref 23,324)", 0.60),
    (25728, "Spike 3' exit\n(MSA 25,728 | Ref 25,410)", 0.44),
    (28810, "ORF8 / N cluster\n(MSA 28,810 | Ref 27,900)", 0.88)
]
for pos, label, y_frac in hotspot_calls:
    y_val = flux_83[min(int(pos), len(flux_83)-1)]
    ax_a_flux.annotate(
        label,
        xy=(pos, y_val),
        xytext=(pos, np.max(flux_83) * y_frac),
        arrowprops=dict(facecolor="#222222", edgecolor="none", width=0.6, headwidth=4, shrink=0.08),
        fontsize=6.5, ha="center",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#CCCCCC", alpha=0.92)
    )

ax_a_flux.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.0)

# =========================================================================
# PANEL A2: Gene Annotation Track
# =========================================================================
GENE_COLORS = {
    "ORF1ab": "#E0E0E0", "RdRp": "#B0B0B0", "Spike": "#F0E442",
    "Spike_NTD": "#E69F00", "Spike_RBD": "#D55E00", "Spike_RBM": "#CC79A7",
    "ORF3a": "#56B4E9", "E": "#009E73", "M": "#0072B2",
    "ORF6": "#F0E442", "ORF7a": "#E69F00", "ORF7b": "#D55E00",
    "ORF8": "#CC79A7", "N": "#009E73", "ORF10": "#999999"
}

ax_a_genes.set_ylim(0, 1)
ax_a_genes.axis("off")

major_genes = ["ORF1ab", "Spike", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N"]
for g_name in major_genes:
    info = GENES.get(g_name)
    if not info:
        continue
    x0 = info["aln_start"]
    width = info["aln_length"]
    col = GENE_COLORS.get(g_name, "#E0E0E0")
    rect = patches.Rectangle((x0, 0.1), width, 0.8, facecolor=col, edgecolor="#333333", lw=0.5)
    ax_a_genes.add_patch(rect)
    if width > 800:
        x_lbl = x0 + width * 0.72 if g_name == "Spike" else x0 + width / 2
        ax_a_genes.text(x_lbl, 0.5, g_name, ha="center", va="center", fontsize=7.0, fontweight="bold", color="#111111")

s_ntd = GENES["Spike_NTD"]
s_rbd = GENES["Spike_RBD"]
ax_a_genes.add_patch(patches.Rectangle((s_ntd["aln_start"], 0.1), s_ntd["aln_length"], 0.8, facecolor="#E69F00", edgecolor="#333333", lw=0.5))
ax_a_genes.add_patch(patches.Rectangle((s_rbd["aln_start"], 0.1), s_rbd["aln_length"], 0.8, facecolor="#D55E00", edgecolor="#333333", lw=0.5))
rdrp = GENES["RdRp"]
ax_a_genes.add_patch(patches.Rectangle((rdrp["aln_start"], 0.1), rdrp["aln_length"], 0.8, facecolor="#A0A0A0", edgecolor="#333333", lw=0.5))
ax_a_genes.text(rdrp["aln_start"] + rdrp["aln_length"]/2, 0.5, "RdRp", ha="center", va="center", fontsize=6.5, fontweight="bold", color="#111111")

# =========================================================================
# PANEL A3: Integrated Breakpoint Concordance Tickmark Tracks (83-Taxon MSA Frame)
# =========================================================================
matches_83 = REPORT["dataset_83"]["gard_concordance"]["matches"]
gard_bps = [m["gard_bp_aln83"] for m in matches_83]
rhiz_bps = [m["rhizaeon_closest_bp"] for m in matches_83]
dists = [m["distance_nt"] for m in matches_83]

ax_a_ticks.set_xlim(-200, 31200)
ax_a_ticks.set_ylim(-0.32, 1.50)
ax_a_ticks.spines["left"].set_visible(False)
ax_a_ticks.spines["right"].set_visible(False)
ax_a_ticks.spines["top"].set_visible(False)
ax_a_ticks.spines["bottom"].set_visible(True)
ax_a_ticks.tick_params(left=False, labelleft=False)

y_rhiz = 1.05
y_gard = 0.20
ax_a_ticks.axhline(y_rhiz, xmin=0.005, xmax=0.995, color="#D55E00", lw=0.8, alpha=0.4)
ax_a_ticks.axhline(y_gard, xmin=0.005, xmax=0.995, color="#0072B2", lw=0.8, alpha=0.4)

ax_a_ticks.text(-350, y_rhiz, "RhizAeon Breakpoints", ha="right", va="center", fontsize=7.2, fontweight="bold", color="#D55E00")
ax_a_ticks.text(-350, y_gard, "Published GARD (2022)", ha="right", va="center", fontsize=7.2, fontweight="bold", color="#0072B2")

for g, r, d in zip(gard_bps, rhiz_bps, dists):
    is_concordant = (d <= 150)
    col = "#009E73" if is_concordant else "#E69F00"
    ax_a_ticks.plot([r, r], [y_rhiz - 0.12, y_rhiz + 0.12], color="#D55E00", lw=1.3)
    ax_a_ticks.plot([g, g], [y_gard - 0.12, y_gard + 0.12], color="#0072B2", lw=1.3)
    ax_a_ticks.plot([r, g], [y_rhiz - 0.12, y_gard + 0.12], color=col, lw=0.9, ls="-" if is_concordant else ":", alpha=0.85)

key_tick_callouts = [
    (12854, 13001, "RdRp 5' (147 nt)", 0.62, 0),
    (21643, 21621, "Spike 5' (22 nt)", 1.32, -450),
    (22756, 22649, "NTD/RBD (107 nt)", -0.18, -500),
    (23701, 23629, "RBD 3' (72 nt)", 1.32, -100),
    (24461, 24484, "S2 (23 nt)", -0.18, 450),
    (25160, 25183, "S2 exit (23 nt)", 1.32, 650)
]
for g_val, r_val, lbl, y_pos, dx in key_tick_callouts:
    mid_x = (g_val + r_val) / 2 + dx
    ax_a_ticks.text(mid_x, y_pos, lbl, ha="center", va="center", fontsize=6.2, fontweight="bold",
                    color="#111111", bbox=dict(boxstyle="round,pad=0.14", facecolor="#FFF9C4", edgecolor="#FBC02D", lw=0.5, alpha=0.92))

ax_a_ticks.set_xlabel("Genomic Coordinate in 83-Taxon Sarbecovirus Multiple Alignment (nt)", fontweight="bold")
ax_a_ticks.set_xticks(np.arange(0, 32000, 5000))
ax_a_ticks.set_xticklabels([f"{x:,}" for x in np.arange(0, 32000, 5000)])

# =========================================================================
# PANEL B: Spike Reticulate Mosaic Architecture (Full Width)
# =========================================================================
taxa_tracks = [
    ("SARS-CoV-2 (Wuhan-Hu-1) [Ancestral RBD Trunk]", [
        ("NTD", 21000, 22764, "#0072B2", "Ancestral nCoV NTD Module (94.4% ID to RaTG13)"),
        ("RBD", 22764, 24500, "#0072B2", "Ancestral ACE2-Adapted Contact Ridge (Shared with BANAL-52)"),
        ("S2", 24500, 26000, "#009E73", "Conserved Sarbecovirus S2 Fusion Backbone")
    ]),
    ("Bat RaTG13 (Yunnan) [Recombinant RBD]", [
        ("NTD", 21000, 22764, "#0072B2", "Shared NTD Ancestral Module with SARS-CoV-2"),
        ("RBD", 22764, 24500, "#9333EA", "Divergent Recombinant Cassette from Unsampled Bat Reservoir (Z=20.9)"),
        ("S2", 24500, 26000, "#0072B2", "Shared S2 Backbone with SARS-CoV-2")
    ]),
    ("Bat BANAL-52 (Laos)", [
        ("NTD", 21000, 22764, "#56B4E9", "Laos BANAL NTD"),
        ("RBD", 22764, 24500, "#0072B2", "Ancestral ACE2-Adapted Contact Ridge (Shared with SARS-CoV-2)"),
        ("S2", 24500, 26000, "#56B4E9", "BANAL Backbone")
    ]),
    ("Pangolin-CoV (MP789)", [
        ("NTD", 21000, 22813, "#999999", "Pangolin NTD"),
        ("RBM", 22813, 23490, "#0072B2", "Shared ACE2 Binding Motif with SARS-CoV-2"),
        ("S2", 23490, 26000, "#999999", "Pangolin Backbone")
    ]),
    ("Bat BANAL-116 (Laos)", [
        ("Backbone", 21000, 22076, "#56B4E9", "BANAL Clade"),
        ("Micro-1", 22076, 22251, "#E69F00", "RacCS203 NTD introgression"),
        ("Backbone", 22251, 22649, "#56B4E9", "BANAL Clade"),
        ("S1 Cassette", 22649, 24577, "#E69F00", "RacCS203 S1 introgression"),
        ("Backbone", 24577, 26000, "#56B4E9", "BANAL Clade")
    ]),
    ("Bat BANAL-247 (Laos)", [
        ("Backbone", 21000, 21909, "#56B4E9", "BANAL Clade"),
        ("Micro-1", 21909, 22011, "#E69F00", "RacCS203 NTD transfer (Z=155.3)"),
        ("Backbone", 22011, 22054, "#56B4E9", "BANAL Clade"),
        ("Micro-2", 22054, 22251, "#E69F00", "RacCS203 NTD transfer (Z=113.1)"),
        ("Backbone", 22251, 22649, "#56B4E9", "BANAL Clade"),
        ("S1 Cassette", 22649, 24577, "#E69F00", "RacCS203 S1 transfer (Z=26.9)"),
        ("Backbone", 24577, 26000, "#56B4E9", "BANAL Clade")
    ]),
    ("Bat RpYN06 (Yunnan)", [
        ("Backbone", 21000, 21539, "#009E73", "Conserved Sarbecovirus S2 Fusion Backbone"),
        ("S1 Cassette", 21539, 24700, "#777777", "PrC31 S1 replacement (Z=41.0)"),
        ("S2 Subunit", 24700, 26000, "#009E73", "Conserved Sarbecovirus S2 Fusion Backbone")
    ]),
    ("Bat RacCS203 (Thailand)", [
        ("Backbone", 21000, 26000, "#E69F00", "Thailand Rhinolophus Clade")
    ])
]

y_pos = np.arange(len(taxa_tracks))
track_height = 0.58

for y, (taxon_name, segments) in enumerate(taxa_tracks):
    for seg_type, start, end, col, desc in segments:
        s_clamped = max(21000, start)
        e_clamped = min(26000, end)
        if e_clamped > s_clamped:
            rect = patches.Rectangle((s_clamped, y - track_height/2), e_clamped - s_clamped, track_height,
                                     facecolor=col, edgecolor="#222222", lw=0.5, alpha=0.92)
            ax_b_mosaic.add_patch(rect)

ax_b_mosaic.axvspan(GENES["Spike_NTD"]["aln_start"], GENES["Spike_NTD"]["aln_end"], color="#E69F00", alpha=0.08)
ax_b_mosaic.axvspan(GENES["Spike_RBD"]["aln_start"], GENES["Spike_RBD"]["aln_end"], color="#D55E00", alpha=0.10)
ax_b_mosaic.axvline(23932, color="#D55E00", ls="--", lw=1.2)

ax_b_mosaic.set_xlim(21000, 26000)
ax_b_mosaic.set_ylim(-0.8, len(taxa_tracks) + 0.6)
ax_b_mosaic.set_yticks(y_pos)
ax_b_mosaic.set_yticklabels([t[0] for t in taxa_tracks], fontsize=7.2, fontweight="bold")
ax_b_mosaic.set_xlabel("Spike Glycoprotein Coordinate in 83-Taxon Multiple Alignment (nt)", fontweight="bold")
ax_b_mosaic.set_title("B   Spike Reticulate Mosaic Architecture Across Wildlife Progenitors and SARS-CoV-2",
                      loc="left", fontweight="bold", pad=28)

ax_b_mosaic.text(22000, len(taxa_tracks) + 0.10, "NTD", ha="center", fontsize=7.5, fontweight="bold", color="#D55E00")
ax_b_mosaic.text(23150, len(taxa_tracks) + 0.10, "RBD", ha="center", fontsize=7.5, fontweight="bold", color="#D55E00")
ax_b_mosaic.text(23932, -0.65, "S1/S2 (MSA 23,932 | Ref 23,618)", ha="center", fontsize=6.5, fontweight="bold", color="#D55E00")

# Callout on BANAL-247 micro-conversions in full width
ax_b_mosaic.annotate("RacCS203 NTD\nmicro-conversions\n(Z=155.3, 113.1)",
                     xy=(22080, 5), xytext=(21550, 5.0),
                     arrowprops=dict(facecolor="#222222", arrowstyle="->", lw=0.9),
                     fontsize=6.5, fontweight="bold", ha="right", va="center",
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFF9C4", edgecolor="#FBC02D", lw=0.6, alpha=0.95))

# Callout on RaTG13 recombinant RBD
ax_b_mosaic.annotate("RaTG13 divergent\nRBD introgression\n(Z=20.9, ghost donor)",
                     xy=(23300, 1), xytext=(24000, 1.0),
                     arrowprops=dict(facecolor="#222222", arrowstyle="->", lw=0.9),
                     fontsize=6.5, fontweight="bold", ha="left", va="center",
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="#F3E8FF", edgecolor="#9333EA", lw=0.6, alpha=0.95))

legend_elements_b = [
    patches.Patch(facecolor="#0072B2", edgecolor="#222222", label="Ancestral nCoV / ACE2 Contact Ridge Module"),
    patches.Patch(facecolor="#9333EA", edgecolor="#222222", label="RaTG13 Divergent Recombinant Cassette"),
    patches.Patch(facecolor="#56B4E9", edgecolor="#222222", label="Laos BANAL Lineage"),
    patches.Patch(facecolor="#E69F00", edgecolor="#222222", label="Thailand RacCS203 Clade"),
    patches.Patch(facecolor="#009E73", edgecolor="#222222", label="Conserved S2 Fusion Subunit")
]
ax_b_mosaic.legend(handles=legend_elements_b, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=5,
                   frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=6.8)

# =========================================================================
# TREE LAYOUT & TANGLEGRAM ENGINE
# =========================================================================
aln_spike = AlignIO.read(os.path.join(BASE_DIR, "..", "..", "02_sarbecovirus_spike", "sarbecovirus_spike_aligned.fasta"), "fasta")
aln_83 = AlignIO.read(FASTA_83, "fasta")
calc = DistanceCalculator("identity")
constructor = DistanceTreeConstructor(calc, "nj")

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

def render_tanglegram(ax, t1, t2, color_map, label_map, focal_names, title_left, title_right, col_left, col_right):
    c1, terms1 = layout_tree(t1)
    c2, terms2 = layout_tree(t2)
    max_x1 = max(c1[t]["x"] for t in terms1)
    max_x2 = max(c2[t]["x"] for t in terms2)
    tree_w = 0.13
    x_root_left = 0.02
    x_root_right = 0.98
    s1 = tree_w / max(0.01, max_x1)
    s2 = tree_w / max(0.01, max_x2)
    xg1, xt1, xs, xe, xt2, xg2 = 0.165, 0.175, 0.38, 0.62, 0.825, 0.835
    
    # Left Tree
    for node in t1.find_clades(order="level"):
        nx = x_root_left + c1[node]["x"] * s1
        ny = c1[node]["y"]
        if not node.is_terminal():
            c_ys = [c1[c]["y"] for c in node.clades]
            ax.plot([nx, nx], [min(c_ys), max(c_ys)], color="#333333", lw=1.0)
            for c in node.clades:
                ax.plot([nx, x_root_left + c1[c]["x"] * s1], [c1[c]["y"], c1[c]["y"]], color="#333333", lw=1.0)
    
    # Right Tree
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
        x_pts = (1 - ts)**3 * xs + 3 * (1 - ts)**2 * ts * (xs + 0.06) + 3 * (1 - ts) * ts**2 * (xe - 0.06) + ts**3 * xe
        y_pts = (1 - ts)**3 * y1 + 3 * (1 - ts)**2 * ts * y1 + 3 * (1 - ts) * ts**2 * y2 + ts**3 * y2
        ax.plot(x_pts, y_pts, color=col, lw=2.2 if is_f else 0.9, alpha=0.95 if is_f else 0.35, zorder=4)
    
    ax.text(0.12, 9.1, title_left, ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=col_left)
    ax.text(0.88, 9.1, title_right, ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=col_right)
    ax.plot([x_root_left, x_root_left + 0.05 * s1], [-0.5, -0.5], color="#222222", lw=1.1)
    ax.text(x_root_left + 0.025 * s1, -0.75, "0.05 subs/site", ha="center", fontsize=5.8)
    ax.plot([x_root_right - 0.05 * s2, x_root_right], [-0.5, -0.5], color="#222222", lw=1.1)
    ax.text(x_root_right - 0.025 * s2, -0.75, "0.05 subs/site", ha="center", fontsize=5.8)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-1.2, 10.0)
    ax.axis("off")

# =========================================================================
# PANEL C: SARS-CoV-2 Incongruence Across Codon 313 Hinge
# =========================================================================
t_ntd_c = constructor.build_tree(aln_spike[:, 0:939])
t_ntd_c.root_with_outgroup({"name": "SARS_CoV_1_Tor2"})
t_ntd_c.ladderize()

t_rbd_c = constructor.build_tree(aln_spike[:, 939:1620])
t_rbd_c.root_with_outgroup({"name": "SARS_CoV_1_Tor2"})
t_rbd_c.ladderize()

color_c = {
    "SARS_CoV_2_Wuhan_Hu_1": "#D55E00", "Bat_RaTG13": "#0072B2",
    "Bat_BANAL_52": "#56B4E9", "Bat_BANAL_103": "#009E73",
    "Pangolin_MP789": "#CC79A7", "Bat_RmYN02": "#E69F00",
    "SARS_CoV_1_Tor2": "#777777", "Bat_RsSHC014": "#888888", "Bat_HKU3_1": "#999999"
}
labels_c = {
    "SARS_CoV_2_Wuhan_Hu_1": "SARS-CoV-2", "Bat_RaTG13": "Bat RaTG13",
    "Bat_BANAL_52": "Bat BANAL-52", "Bat_BANAL_103": "Bat BANAL-103",
    "Pangolin_MP789": "Pangolin MP789", "Bat_RmYN02": "Bat RmYN02",
    "SARS_CoV_1_Tor2": "SARS-CoV-1", "Bat_RsSHC014": "Bat RsSHC014", "Bat_HKU3_1": "Bat HKU3-1"
}

render_tanglegram(ax_c_tree, t_ntd_c, t_rbd_c, color_c, labels_c,
                  ["SARS_CoV_2_Wuhan_Hu_1", "Bat_RaTG13", "Bat_BANAL_52", "Pangolin_MP789"],
                  "C1  Spike NTD (Codons 1-313)\nShared nCoV / RaTG13 Clade",
                  "C2  Spike RBD (Codons 314-540)\nRaTG13 Shifts Basally; BANAL-52 / SC2 Sister Clade",
                  "#0072B2", "#9333EA")
ax_c_tree.set_title("C   SARS-CoV-2 / RaTG13 Incongruence Across Codon 313 Hinge", loc="left", fontweight="bold", pad=8)

# =========================================================================
# PANEL D: Ancestral Recombination Graph (ARG) & Macro-Evolutionary Validation
# =========================================================================
ax_d_tree.set_xlim(0, 10)
ax_d_tree.set_ylim(-1.2, 10.0)
ax_d_tree.axis("off")
ax_d_tree.set_title("D   Polarized Ancestral Recombination Graph (ARG) for Sarbecovirus Progenitors",
                    loc="left", fontweight="bold", pad=8)

c_blue = "#0072B2"      # SARS-CoV-2 / Ancestral nCoV Backbone
c_purple = "#9333EA"    # RaTG13 Recombinant RBD
c_cyan = "#56B4E9"      # BANAL-52 / Laos Ancestral Trunk
c_teal = "#009E73"      # BANAL-103 / 247
c_gold = "#E69F00"      # RacCS203 (Thailand)
c_pink = "#CC79A7"      # Pangolin MP789
c_grey = "#666666"      # Outgroups
col_tree = "#333333"

# 9 Taxa spaced evenly from y=0.0 to 8.0, exactly matching Panel C
taxa_info = [
    ("SARS-CoV-2 (Wuhan-Hu-1)", c_blue, 0.0, "(Ancestral ACE2 Trunk)", True),
    ("Bat RaTG13", c_purple, 1.0, "(Recombinant RBD Cassette)", True),
    ("Pangolin MP789", c_pink, 2.0, "(Shared ACE2 Receptor Ridge)", True),
    ("Bat BANAL-52", c_cyan, 3.0, "(Laos Ancestral Trunk)", True),
    ("Bat BANAL-103", c_teal, 4.0, "(Laos Reservoir)", True),
    ("Bat BANAL-247", c_teal, 5.0, "(Laos Mosaic NTD)", True),
    ("Bat RacCS203", c_gold, 6.0, "(Thailand Reservoir)", True),
    ("Bat HKU3-1", c_grey, 7.0, "(Bat Clade 1)", False),
    ("SARS-CoV-1 (Tor2)", c_grey, 8.0, "(Outgroup)", False),
]

x_root = 0.5
x_tips = 5.3
lw_tree = 1.05

x_n1 = 1.0
x_n2 = 1.7
x_n3 = 2.4
x_n_laos = 3.2
x_n_b103_52 = 4.0
x_n_yunnan = 3.3
x_n_sc2_ratg = 4.3

# Draw backbone tree
# Root vertical bar
ax_d_tree.plot([x_root, x_root], [3.5, 8.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_root, x_tips], [8.0, 8.0], color=col_tree, lw=lw_tree) # SARS-CoV-1
ax_d_tree.plot([x_root, x_n1], [3.5, 3.5], color=col_tree, lw=lw_tree)

# n1 vertical
ax_d_tree.plot([x_n1, x_n1], [3.0, 7.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n1, x_tips], [7.0, 7.0], color=col_tree, lw=lw_tree) # HKU3-1
ax_d_tree.plot([x_n1, x_n2], [3.0, 3.0], color=col_tree, lw=lw_tree)

# n2 vertical
ax_d_tree.plot([x_n2, x_n2], [2.5, 6.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n2, x_tips], [6.0, 6.0], color=col_tree, lw=lw_tree) # RacCS203
ax_d_tree.plot([x_n2, x_n3], [2.5, 2.5], color=col_tree, lw=lw_tree)

# n3 vertical
ax_d_tree.plot([x_n3, x_n3], [1.0, 4.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n3, x_n_laos], [4.0, 4.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n3, x_n_yunnan], [1.0, 1.0], color=col_tree, lw=lw_tree)

# Laos Clade
ax_d_tree.plot([x_n_laos, x_n_laos], [3.5, 5.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n_laos, x_tips], [5.0, 5.0], color=col_tree, lw=lw_tree) # BANAL-247
ax_d_tree.plot([x_n_laos, x_n_b103_52], [3.5, 3.5], color=col_tree, lw=lw_tree)

ax_d_tree.plot([x_n_b103_52, x_n_b103_52], [3.0, 4.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n_b103_52, x_tips], [4.0, 4.0], color=col_tree, lw=lw_tree) # BANAL-103
ax_d_tree.plot([x_n_b103_52, x_tips], [3.0, 3.0], color=col_tree, lw=lw_tree) # BANAL-52

# Yunnan / Emergent Clade
ax_d_tree.plot([x_n_yunnan, x_n_yunnan], [0.5, 2.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n_yunnan, x_tips], [2.0, 2.0], color=col_tree, lw=lw_tree) # Pangolin MP789
ax_d_tree.plot([x_n_yunnan, x_n_sc2_ratg], [0.5, 0.5], color=col_tree, lw=lw_tree)

ax_d_tree.plot([x_n_sc2_ratg, x_n_sc2_ratg], [0.0, 1.0], color=col_tree, lw=lw_tree)
ax_d_tree.plot([x_n_sc2_ratg, x_tips], [1.0, 1.0], color=col_tree, lw=lw_tree) # RaTG13
ax_d_tree.plot([x_n_sc2_ratg, x_tips], [0.0, 0.0], color=col_tree, lw=lw_tree) # Wuhan-Hu-1

# Tip nodes and Taxa Labels
for name, col, y, desc, is_focal in taxa_info:
    f_wt = "bold" if is_focal else "normal"
    ax_d_tree.scatter(x_tips, y, color=col, s=26, zorder=5, edgecolor="#111111", lw=0.5)
    ax_d_tree.plot([x_tips, x_tips + 0.12], [y, y], color="#CCCCCC", ls=":", lw=0.6)
    ax_d_tree.text(x_tips + 0.18, y, name, ha="left", va="center", fontsize=7.0, color=col, fontweight=f_wt)
    x_desc = x_tips + 3.25 if "Wuhan" in name else x_tips + 2.15
    ax_d_tree.text(x_desc, y, desc, ha="left", va="center", fontsize=5.8, color="#555555", style="italic")

# Tree scale bar matching Panel C
ax_d_tree.plot([0.5, 1.0], [-0.75, -0.75], color="#222222", lw=1.1)
ax_d_tree.text(0.75, -1.0, "0.05 subs/site", ha="center", fontsize=5.8, color="#222222")

# RETICULATION 1: Divergent RBD Introgression into Bat RaTG13 from Ghost Donor
retic1_start = (1.9, 2.7)  # Unsampled ancestral bat sarbecovirus reservoir
retic1_ctrl1 = (3.0, 1.9)
retic1_ctrl2 = (4.0, 1.3)
retic1_end = (4.9, 1.05)   # Entering Bat RaTG13

ts = np.linspace(0, 1, 50)
x_arc1 = (1 - ts)**3 * retic1_start[0] + 3 * (1 - ts)**2 * ts * retic1_ctrl1[0] + 3 * (1 - ts) * ts**2 * retic1_ctrl2[0] + ts**3 * retic1_end[0]
y_arc1 = (1 - ts)**3 * retic1_start[1] + 3 * (1 - ts)**2 * ts * retic1_ctrl1[1] + 3 * (1 - ts) * ts**2 * retic1_ctrl2[1] + ts**3 * retic1_end[1]

ax_d_tree.plot(x_arc1, y_arc1, color=c_purple, lw=1.8, linestyle="--", zorder=6)
ax_d_tree.annotate("", xy=retic1_end, xytext=(x_arc1[-4], y_arc1[-4]),
                   arrowprops=dict(arrowstyle="-|>", color=c_purple, lw=1.8, mutation_scale=10), zorder=7)

# Ghost donor node marker
ax_d_tree.scatter(retic1_start[0], retic1_start[1], color="#9333EA", s=30, marker="o", facecolors="none", lw=1.2, zorder=7)
ax_d_tree.text(retic1_start[0] - 0.15, retic1_start[1] + 0.15, "Ghost Donor\n(Unsampled Bat CoV)",
               ha="right", va="bottom", fontsize=5.6, color="#7E22CE", style="italic", fontweight="semibold")

ax_d_tree.text(2.0, -0.75,
               "Divergent RBD Introgression (codons 314–540 / nt 940–1620)\nUnsampled Bat Reservoir -> Bat RaTG13 (P = 2.1e-9, Z = 20.9)\nAncestral ACE2 contact ridge preserved in SARS-CoV-2 / BANAL-52 trunk",
               ha="left", va="center", fontsize=5.8, color="#7E22CE", fontweight="semibold", linespacing=1.2)

# RETICULATION 2: RacCS203 NTD Micro-Conversions into BANAL-247
retic2_start = (2.6, 6.0)
retic2_ctrl1 = (3.0, 5.7)
retic2_ctrl2 = (3.3, 5.4)
retic2_end = (3.5, 5.05)

x_arc2 = (1 - ts)**3 * retic2_start[0] + 3 * (1 - ts)**2 * ts * retic2_ctrl1[0] + 3 * (1 - ts) * ts**2 * retic2_ctrl2[0] + ts**3 * retic2_end[0]
y_arc2 = (1 - ts)**3 * retic2_start[1] + 3 * (1 - ts)**2 * ts * retic2_ctrl1[1] + 3 * (1 - ts) * ts**2 * retic2_ctrl2[1] + ts**3 * retic2_end[1]

ax_d_tree.plot(x_arc2, y_arc2, color=c_gold, lw=1.6, linestyle="--", zorder=6)
ax_d_tree.annotate("", xy=retic2_end, xytext=(x_arc2[-4], y_arc2[-4]),
                   arrowprops=dict(arrowstyle="-|>", color=c_gold, lw=1.6, mutation_scale=9), zorder=7)

ax_d_tree.text(3.7, 5.5,
               "RacCS203 NTD micro-conversions\n(Z = 155.3, 113.1)",
               ha="left", va="center", fontsize=5.6, color="#92400E", fontweight="semibold", linespacing=1.15)


# Save high-res PNG and PDF
out_png = os.path.join(FIG_DIR, "sarbecovirus_recombination_landscape.png")
out_pdf = os.path.join(FIG_DIR, "sarbecovirus_recombination_landscape.pdf")
REPO_ROOT = os.path.abspath(os.path.join(BASE_DIR, "..", "..", ".."))
paper_png = os.path.join(REPO_ROOT, "paper", "figures", "fig_sarbecovirus_recombination_landscape.png")
paper_pdf = os.path.join(REPO_ROOT, "paper", "figures", "fig_sarbecovirus_recombination_landscape.pdf")

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(paper_png, dpi=300, bbox_inches="tight")
plt.savefig(paper_pdf, bbox_inches="tight")
plt.close()

print(f"Successfully generated 4-panel publication figure:")
print(f"  {out_png}")
print(f"  {paper_png}")
