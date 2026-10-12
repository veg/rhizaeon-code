#!/usr/bin/env python3
"""
plot_pmen1_chromosome_landscape.py
==================================
Generates publication-quality 4-panel figure for whole-chromosome recombination
deconvolution in Streptococcus pneumoniae PMEN1 (Croucher et al. 2011 Science / 2015 NAR).
Chromosome: ATCC 700669 (GenBank FM211187, L = 2,221,315 bp).
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
from matplotlib.gridspec import GridSpec
from Bio import SeqIO

# BioVis-Expert styling standards
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.labelsize": 9.2,
    "axes.titlesize": 10.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 7.5,
    "figure.titlesize": 11.5,
    "lines.linewidth": 1.2,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "mathtext.fontset": "dejavusans",
    "figure.dpi": 300
})

BASE_DIR = "/Users/sergei/Projects/rhizaeon-code/empirical_showcase/spneumoniae_pmen1_croucher2011"
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIG_DIR = os.path.join(BASE_DIR, "figures")
PAPER_FIG_DIR = "/Users/sergei/Projects/rhizaeon-code/paper/figures"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(PAPER_FIG_DIR, exist_ok=True)

# 1. Load Data
fields_npz = np.load(os.path.join(RESULTS_DIR, "pmen1_continuous_fields.npz"))
u_grid = fields_npz["u_grid"]
phi_field = fields_npz["phi_field"]
strain_field = fields_npz["strain_field"]

with open(os.path.join(RESULTS_DIR, "pmen1_fields_summary.json")) as f:
    SUMMARY = json.load(f)

DF_CP = pd.read_csv(os.path.join(RESULTS_DIR, "pmen1_changepoints_vs_gubbins.csv"))

# Load Alignment for Panel C
aln_path = os.path.join(DATA_DIR, "PMEN1.aln")
records = list(SeqIO.parse(aln_path, "fasta"))
taxa = [r.id for r in records]
seqs = np.array([list(str(r.seq).upper()) for r in records])
valid_chars = set("ACGT")
is_seg = np.zeros(seqs.shape[1], dtype=bool)
for col in range(seqs.shape[1]):
    chars = set(seqs[:, col]) & valid_chars
    if len(chars) > 1:
        is_seg[col] = True
seg_indices = np.where(is_seg)[0]

# 2. Setup Figure Layout
fig = plt.figure(figsize=(11.5, 9.8))

# Outer GridSpec: Top row (Panel A + gene track) vs Bottom section (Panels B, C, D)
gs_outer = GridSpec(2, 1, figure=fig, height_ratios=[1.10, 2.0], hspace=0.45)

# Top Sub-GridSpec for Panel A
gs_top = gs_outer[0].subgridspec(2, 1, height_ratios=[2.0, 0.42], hspace=0.10)
ax_a = fig.add_subplot(gs_top[0])
ax_a_genes = fig.add_subplot(gs_top[1])

# Bottom Sub-GridSpec for Panels B, C, D
# Left column: Panel B (Concordance) + Panel D (Latency)
# Right column: Panel C top (SNP density) + Panel C bottom (Mosaics)
gs_bot = gs_outer[1].subgridspec(2, 2, width_ratios=[0.92, 1.08], height_ratios=[1.0, 1.0], hspace=0.52, wspace=0.30)

ax_b = fig.add_subplot(gs_bot[0, 0])      # Panel B: Concordance Scatter
ax_d = fig.add_subplot(gs_bot[1, 0])      # Panel D: Computational Latency
ax_c = fig.add_subplot(gs_bot[0, 1])      # Panel C: SNP Density Zoom (270-370 kb)
ax_c_sub = fig.add_subplot(gs_bot[1, 1])  # Panel C Sub: Mosaic Structure (270-370 kb)

# =========================================================================
# PANEL A: Whole-Chromosome Landscape (2.22 Mb)
# =========================================================================
ax_a_strain = ax_a.twinx()
ax_a_strain.spines["top"].set_visible(False)
ax_a.spines["top"].set_visible(False)

# Breakpoint Flux (Primary, left, Vermilion)
line_flux, = ax_a.plot(u_grid, phi_field * 1000.0, color="#D55E00", lw=1.2, label="Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ (Z / kb)")
ax_a.fill_between(u_grid, 0, phi_field * 1000.0, color="#D55E00", alpha=0.20)

# Cumulative Strain (Secondary, right, Blue)
line_strain, = ax_a_strain.plot(u_grid, strain_field, color="#0072B2", lw=1.2, ls="--", label="Cumulative Strain $\\rho(u)$ (integrated Z)")
ax_a_strain.fill_between(u_grid, 0, strain_field, color="#0072B2", alpha=0.07)

ax_a.set_xlim(0, 2221315)
ax_a.set_ylim(0, np.max(phi_field * 1000.0) * 1.35)
ax_a_strain.set_ylim(0, np.max(strain_field) * 1.15)

ax_a.set_ylabel("Breakpoint Flux $\\Phi_{\\mathrm{bp}}$ (Z / kb)", color="#D55E00", fontweight="bold")
ax_a_strain.set_ylabel("Cumulative Strain $\\rho(u)$ (Z)", color="#0072B2", fontweight="bold")
ax_a.tick_params(axis="y", labelcolor="#D55E00")
ax_a_strain.tick_params(axis="y", labelcolor="#0072B2")
ax_a.tick_params(axis="x", labelbottom=False)

# Major Hotspot Callouts on Panel A
hotspots_annot = [
    (294000, "pbp2x / cps island\n(285-352 kb)", 0.78, 294000),
    (740000, "732 kb hotspot\n(733-744 kb)", 0.72, 740000),
    (1155000, "1.15 Mb hotspot\n(1,154-1,159 kb)", 0.62, 1155000),
    (1635000, "folA trimethoprim\n(1,620-1,649 kb)", 0.76, 1635000),
    (1722000, "1.72 Mb hotspot", 0.52, 1722000),
    (2135000, "pspC & 2.13 Mb\n(2,105-2,138 kb)", 0.78, 2040000)
]

for pos, label, y_frac, text_x in hotspots_annot:
    idx_g = min(int(pos / 1000.0), len(phi_field) - 1)
    y_val = phi_field[idx_g] * 1000.0
    ax_a.annotate(
        label,
        xy=(pos, y_val),
        xytext=(text_x, ax_a.get_ylim()[1] * y_frac),
        arrowprops=dict(facecolor="#333333", edgecolor="none", width=0.6, headwidth=4, shrink=0.08),
        fontsize=6.5,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#CCCCCC", alpha=0.92)
    )

ax_a.legend([line_flux, line_strain], [line_flux.get_label(), line_strain.get_label()],
            loc="upper left", frameon=True, facecolor="white", framealpha=0.95, edgecolor="#DDDDDD")
ax_a.set_title("A   Whole-Chromosome Continuous Recombination Landscape Across S. pneumoniae PMEN1 (L = 2.22 Mb, 11 Taxa)",
               loc="left", fontweight="bold", pad=8)

# Dedicated Gene Model Track for Panel A
ax_a_genes.set_ylim(0, 1)
ax_a_genes.set_xlim(0, 2221315)
ax_a_genes.spines["left"].set_visible(False)
ax_a_genes.spines["bottom"].set_visible(True)
ax_a_genes.tick_params(axis="y", left=False, labelleft=False)

# Chromosome backbone bar
ax_a_genes.add_patch(patches.Rectangle((0, 0.35), 2221315, 0.3, facecolor="#EBEBEB", edgecolor="#888888", lw=0.6))

# Landmark Genes
gene_blocks = [
    (73625, 76580, "73kb", "#999999"),
    (285069, 305751, "pbp2x", "#D55E00"),
    (303559, 314822, "cps (23F->19A)", "#CC79A7"),
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
    rect = patches.Rectangle((s, 0.15), width, 0.7, facecolor=col, edgecolor="#222222", lw=0.5)
    ax_a_genes.add_patch(rect)
    if width > 15000:
        ax_a_genes.text(s + width/2, 0.5, gname, ha="center", va="center", fontsize=6, fontweight="bold", color="white")

ax_a_genes.set_xlabel("Genomic Coordinate along S. pneumoniae ATCC 700669 Chromosome (bp)", fontweight="bold", labelpad=5)
xticks_a = [0, 500000, 1000000, 1500000, 2000000, 2221315]
ax_a_genes.set_xticks(xticks_a)
ax_a_genes.set_xticklabels(["0", "500 kb", "1.0 Mb", "1.5 Mb", "2.0 Mb", "2.22 Mb"])


# =========================================================================
# PANEL B: Head-to-Head Concordance with Gubbins (52 Boundaries)
# =========================================================================
matches = SUMMARY["gubbins_concordance"]["matches"]
gb_pts = [m["gubbins_bp"] for m in matches]
rhiz_pts = [m["closest_rhizaeon_bp"] for m in matches]
dists = [m["distance_bp"] for m in matches]

ax_b.fill_between([0, 2221315], [-2000, 2221315 - 2000], [2000, 2221315 + 2000],
                  color="#E8E8E8", alpha=0.7, label="$\\pm 2.0$ kb corridor")
ax_b.plot([0, 2221315], [0, 2221315], color="#888888", ls="--", lw=1.0, label="1:1 identity line")

# Scatter points color-coded by distance
for gb, rb, d in zip(gb_pts, rhiz_pts, dists):
    if d <= 1000:
        col = "#009E73"  # Green <= 1 kb
    elif d <= 2000:
        col = "#0072B2"  # Blue <= 2 kb
    else:
        col = "#D55E00"  # Orange <= 5 kb
    ax_b.scatter(gb, rb, color=col, s=26, zorder=4, edgecolor="#222222", lw=0.5)

# Concordance Callouts for key genes
b_callouts = [
    (285069, 285068, "pbp2x (1 bp)", (-450000, 190000)),
    (336730, 340596, "cps 5' (3.8 kb)", (-460000, 320000)),
    (1620528, 1620547, "folA (19 bp)", (-420000, 180000)),
    (2104999, 2104998, "pspC (1 bp)", (-380000, -160000)),
    (73625, 76579, "73 kb (2.9 kb)", (320000, -220000)),
    (2134136, 2138664, "2.13 Mb (4.5 kb)", (-240000, -320000))
]
for gb, rb, txt, offset in b_callouts:
    ax_b.annotate(
        txt,
        xy=(gb, rb),
        xytext=(gb + offset[0], rb + offset[1]),
        arrowprops=dict(arrowstyle="->", color="#333333", lw=0.6),
        fontsize=6.2,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="#FFFDE7", edgecolor="#CCCCCC", alpha=0.92)
    )

ax_b.set_xlim(0, 2221315)
ax_b.set_ylim(0, 2221315)
ax_b.set_aspect("equal")
ax_b.set_xlabel("Published Gubbins Boundaries (bp)", fontweight="bold", labelpad=3)
ax_b.set_ylabel("RhizAeon Changepoints (bp)", fontweight="bold")
ax_b.set_xticks([0, 500000, 1000000, 1500000, 2000000])
ax_b.set_yticks([0, 500000, 1000000, 1500000, 2000000])
ax_b.set_xticklabels(["0", "0.5", "1.0", "1.5", "2.0 Mb"])
ax_b.set_yticklabels(["0", "0.5", "1.0", "1.5", "2.0 Mb"])

summary_txt = (
    "100.0% within 5 kb (52/52)\n"
    "80.8% within 1 kb (42/52)\n"
    "Median: 92.5 bp (sub-SNP)\n"
    "Mean: 460.1 bp"
)
ax_b.text(0.04, 0.65, summary_txt, transform=ax_b.transAxes, fontsize=6.5,
          bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#BBBBBB", alpha=0.95))

ax_b.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#DDDDDD", fontsize=6.5)
ax_b.set_title("B   Concordance with 52 Gubbins Boundaries\n     (Median Deviation: 92.5 bp; 100% $\\leq 5$ kb)",
               loc="left", fontweight="bold", pad=6)


# =========================================================================
# PANEL D: Computational Latency & Algorithmic Scaling Comparison
# =========================================================================
benchmarks = [
    ("ClonalFrame\n(2011)", 2538000.0, "#999999", "705 hr\n(~1 month)"),
    ("Gubbins\n(2015)", 5400.0, "#56B4E9", "1.5 hr\n(90 min)"),
    ("RhizAeon\n(Python)", 0.114, "#0072B2", "114 ms\n(0.11 s)"),
    ("RhizAeon\n(Rust Core)", 0.127, "#D55E00", "127 ms\n(0.13 s)")
]

b_names = [b[0] for b in benchmarks]
b_times = [b[1] for b in benchmarks]
b_colors = [b[2] for b in benchmarks]

bars = ax_d.bar(range(len(b_names)), b_times, color=b_colors, edgecolor="#222222", lw=0.6, width=0.52)
ax_d.set_yscale("log")
ax_d.set_ylim(0.01, 1e8)
ax_d.set_xlim(-0.6, 3.6)
ax_d.set_xticks(range(len(b_names)))
ax_d.set_xticklabels(b_names, fontsize=7.2, fontweight="bold")
ax_d.set_ylabel("Execution Time (seconds, log scale)", fontweight="bold")
ax_d.set_title("D   Computational Latency & Scaling\n     (>30,000x Speedup Over Gubbins; >18,000,000x Over ClonalFrame)",
               loc="left", fontweight="bold", pad=10)

# Speedup Callouts
ax_d.text(0, 2538000 * 1.8, "705 hr\n(~1 mo)", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#555555")
ax_d.text(1, 5400 * 1.8, "1.5 hr\n(5,400 s)", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#0072B2")
ax_d.text(2, 0.114 * 2.2, "114 ms", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#0072B2")
ax_d.text(3, 0.127 * 2.2, "127 ms", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#D55E00")

# Horizontal Reference Lines inside ax_d
ax_d.axhline(1.0, color="#AAAAAA", ls=":", lw=0.7)
ax_d.text(1.6, 1.0 * 1.25, "1 second", va="bottom", ha="left", fontsize=6.0, color="#777777")
ax_d.axhline(60.0, color="#AAAAAA", ls=":", lw=0.7)
ax_d.text(1.6, 60.0 * 1.25, "1 minute", va="bottom", ha="left", fontsize=6.0, color="#777777")
ax_d.axhline(3600.0, color="#AAAAAA", ls=":", lw=0.7)
ax_d.text(1.6, 3600.0 * 1.25, "1 hour", va="bottom", ha="left", fontsize=6.0, color="#777777")
ax_d.axhline(86400.0, color="#AAAAAA", ls=":", lw=0.7)
ax_d.text(1.6, 86400.0 * 1.25, "1 day", va="bottom", ha="left", fontsize=6.0, color="#777777")
ax_d.axhline(2592000.0, color="#AAAAAA", ls=":", lw=0.7)
ax_d.text(0.35, 2592000.0 * 1.25, "1 month", va="bottom", ha="left", fontsize=6.0, color="#777777")


# =========================================================================
# PANEL C: Micro-Resolution Zoom on AMR & Capsule Switch (270 - 370 kb)
# =========================================================================
zoom_s, zoom_e = 270000, 370000
mask_snps = (seg_indices >= zoom_s) & (seg_indices <= zoom_e)
zoom_snps = seg_indices[mask_snps]

bins = np.arange(zoom_s, zoom_e + 1000, 1000)
snp_counts, _ = np.histogram(zoom_snps, bins=bins)
bin_centers = (bins[:-1] + bins[1:]) / 2.0

ax_c.plot(bin_centers, snp_counts, color="#333333", lw=1.0, label="Segregating SNPs / kb")
ax_c.fill_between(bin_centers, 0, snp_counts, color="#BBBBBB", alpha=0.35)

# Shaded gene regions
ax_c.axvspan(285069, 305751, color="#D55E00", alpha=0.18, label="pbp2x island (285-306 kb)")
ax_c.axvspan(303559, 314822, color="#CC79A7", alpha=0.22, label="cps operon (304-315 kb)")
ax_c.axvspan(336730, 352271, color="#56B4E9", alpha=0.25, label="Capsule switch 23F->19A")
ax_c.axvspan(323721, 325108, color="#009E73", alpha=0.18, label="pbp1a (324-325 kb)")

# Vertical dashed lines for changepoints
cps_in_zoom = DF_CP[(DF_CP["genome_coord"] >= zoom_s) & (DF_CP["genome_coord"] <= zoom_e)]
for _, row in cps_in_zoom.iterrows():
    cp_nt = row["genome_coord"]
    z_val = row["kinetic_z"]
    if z_val >= 4.0:
        ax_c.axvline(cp_nt, color="#D55E00", ls=":", lw=1.1, alpha=0.85)

ax_c.set_xlim(zoom_s, zoom_e)
ax_c.set_ylim(0, max(snp_counts) * 1.25)
ax_c.set_ylabel("SNP Density (SNPs / kb)", fontweight="bold")
ax_c.tick_params(axis="x", labelbottom=False)
ax_c.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#DDDDDD", fontsize=6.5)
ax_c.set_title("C   AMR (pbp2x/pbp1a) & Capsule Switch (cps) Super-Hotspot (270--370 kb)\n     (Segregating SNP Density and Physical Changepoints)",
               loc="left", fontweight="bold", pad=6)

# Panel C Sub: Mosaic Architecture across Key Strains
focal_strains = [
    ("2097 (23F Backbone)", [
        ("Core", zoom_s, zoom_e, "#0072B2", "Core 23F Lineage")
    ]),
    ("1731 (pbp2x Mosaic)", [
        ("Core", zoom_s, 285068, "#0072B2", "Core 23F"),
        ("AMR", 285068, 305751, "#D55E00", "pbp2x Transpeptidase Introgression (Z=5.35)"),
        ("Core", 305751, zoom_e, "#0072B2", "Core 23F")
    ]),
    ("2682 (pbp2x/cps Mosaic)", [
        ("Core", zoom_s, 282338, "#0072B2", "Core 23F"),
        ("AMR/CPS", 282338, 356248, "#D55E00", "pbp2x-cps Co-Transfer (73.9 kb, Z=6.42)"),
        ("Core", 356248, zoom_e, "#0072B2", "Core 23F")
    ]),
    ("3355 (19A Vaccine Escape)", [
        ("Core", zoom_s, 323814, "#0072B2", "Core 23F"),
        ("pbp1a", 323814, 325108, "#009E73", "pbp1a Low-Affinity Allele (Z=3.19)"),
        ("Switch", 336730, 352271, "#CC79A7", "Serotype 19A Capsule Replacement (PCV7 Escape)"),
        ("Core", 352271, zoom_e, "#0072B2", "Core 23F")
    ]),
    ("1789 (Clinical Resistant)", [
        ("Core", zoom_s, 295884, "#0072B2", "Core 23F"),
        ("AMR", 295884, 305751, "#D55E00", "pbp2x High-Level Resistance"),
        ("Core", 305751, zoom_e, "#0072B2", "Core 23F")
    ])
]

y_strains = np.arange(len(focal_strains))
track_h = 0.55

for y_idx, (s_name, segs) in enumerate(focal_strains):
    for stype, s0, s1, scol, desc in segs:
        c_s = max(zoom_s, s0)
        c_e = min(zoom_e, s1)
        if c_e > c_s:
            rect = patches.Rectangle((c_s, y_idx - track_h/2), c_e - c_s, track_h,
                                     facecolor=scol, edgecolor="#222222", lw=0.5, alpha=0.92)
            ax_c_sub.add_patch(rect)

ax_c_sub.set_xlim(zoom_s, zoom_e)
ax_c_sub.set_ylim(-0.6, len(focal_strains) - 0.4)
ax_c_sub.set_yticks(y_strains)
ax_c_sub.set_yticklabels([s[0] for s in focal_strains], fontsize=7.2, fontweight="bold")
ax_c_sub.set_xlabel("Genomic Position (bp, 270 - 370 kb Window)", fontweight="bold")
ax_c_sub.set_xticks(np.arange(270000, 375000, 20000))
ax_c_sub.set_xticklabels([f"{x:,}" for x in np.arange(270000, 375000, 20000)])

patches_c = [
    patches.Patch(facecolor="#0072B2", edgecolor="#222222", label="Core PMEN1 23F Backbone"),
    patches.Patch(facecolor="#D55E00", edgecolor="#222222", label="pbp2x Penicillin Resistance Island"),
    patches.Patch(facecolor="#CC79A7", edgecolor="#222222", label="19A Capsule Switch (PCV7 Escape)"),
    patches.Patch(facecolor="#009E73", edgecolor="#222222", label="pbp1a Low-Affinity Allele")
]
ax_c_sub.legend(handles=patches_c, loc="upper center", bbox_to_anchor=(0.5, -0.28),
                frameon=True, facecolor="white", framealpha=0.95, edgecolor="#DDDDDD", ncol=2, fontsize=6.5)

# Save high-res PNG and PDF
png_local = os.path.join(FIG_DIR, "fig_spneumoniae_chromosome_landscape.png")
pdf_local = os.path.join(FIG_DIR, "fig_spneumoniae_chromosome_landscape.pdf")
png_paper = os.path.join(PAPER_FIG_DIR, "fig_spneumoniae_chromosome_landscape.png")
pdf_paper = os.path.join(PAPER_FIG_DIR, "fig_spneumoniae_chromosome_landscape.pdf")

plt.savefig(png_local, dpi=300, bbox_inches="tight")
plt.savefig(pdf_local, bbox_inches="tight")
plt.savefig(png_paper, dpi=300, bbox_inches="tight")
plt.savefig(pdf_paper, bbox_inches="tight")
plt.close()

print(f"[✓] Publication figure successfully generated:")
print(f"  Local PNG: {png_local}")
print(f"  Local PDF: {pdf_local}")
print(f"  Paper PNG: {png_paper}")
print(f"  Paper PDF: {pdf_paper}")

if __name__ == "__main__":
    pass
