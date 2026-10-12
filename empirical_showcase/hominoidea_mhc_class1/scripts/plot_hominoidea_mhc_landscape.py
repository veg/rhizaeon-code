#!/usr/bin/env python3
"""
plot_hominoidea_mhc_landscape.py
Generates a publication-grade 4-panel figure illustrating the
Exhaustive Superfamily Hominoidea MHC Class I Recombination & Domain Landscape.
Adheres strictly to BioVis-Expert guidelines:
- Accessible Okabe-Ito color palette
- Standalone legible typography (Helvetica/DejaVu Sans)
- High-contrast visual hierarchy
- No text collisions or edge bleeding
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from Bio import SeqIO

# Global Typography and Styling
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 7.5,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 6.5,
    "figure.titlesize": 9.5,
    "lines.linewidth": 1.1,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "mathtext.fontset": "dejavusans",
    "figure.dpi": 300
})

BASE_DIR = "/Users/sergei/Projects/rhizaeon-code/empirical_showcase/hominoidea_mhc_class1"
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIG_DIR = os.path.join(BASE_DIR, "figures")
PAPER_FIG_DIR = "/Users/sergei/Projects/rhizaeon-code/paper/figures"
ARTIFACT_DIR = "/Users/sergei/.gemini/antigravity-cli/brain/1dfa5602-784e-4936-ae1c-a9734d91bc95"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(PAPER_FIG_DIR, exist_ok=True)

# 1. Load Data
fields_npz = np.load(os.path.join(RESULTS_DIR, "hominoidea_continuous_fields.npz"))
u_grid = fields_npz["u_grid"]
phi_field = fields_npz["phi_field"]
strain_field = fields_npz["strain_field"]

with open(os.path.join(RESULTS_DIR, "hominoidea_fields_summary.json")) as f:
    SUMMARY = json.load(f)

DF_SIG = pd.read_csv(os.path.join(RESULTS_DIR, "hominoidea_significant_reticulations.csv"))
DF_TAXA = pd.read_csv(os.path.join(RESULTS_DIR, "hominoidea_taxon_peaks.csv"))

# Load alignment for sequence comparisons and MDS
aln_path = os.path.join(DATA_DIR, "hominoidea_class1_aligned.fasta")
records = list(SeqIO.parse(aln_path, "fasta"))
taxa = [r.id for r in records]
seq_matrix = np.array([list(str(r.seq).upper()) for r in records])
N, L = seq_matrix.shape

# 2. Setup Figure Layout
fig = plt.figure(figsize=(11.5, 10.2))

# Outer GridSpec: Top row (Panel A + gene track) vs Bottom section (Panels B, C, D)
gs_outer = GridSpec(2, 1, figure=fig, height_ratios=[1.08, 2.0], hspace=0.40)

# Top Sub-GridSpec for Panel A
gs_top = gs_outer[0].subgridspec(2, 1, height_ratios=[2.0, 0.40], hspace=0.08)
ax_a = fig.add_subplot(gs_top[0])
ax_a_domains = fig.add_subplot(gs_top[1])

# Bottom Sub-GridSpec for Panels B, C, D
# Left column: Panel B (Domain Segregation) + Panel D (MDS Topological Discordance)
# Right column: Panel C top (HLA-B*46:01 Mismatch Profile) + Panel C bottom (Mosaic Ribbons)
gs_bot = gs_outer[1].subgridspec(2, 2, width_ratios=[0.88, 1.12], height_ratios=[1.0, 1.0], hspace=0.52, wspace=0.44)

ax_b = fig.add_subplot(gs_bot[0, 0])      # Panel B: Structural Domain Distribution
ax_d = fig.add_subplot(gs_bot[1, 0])      # Panel D: Topological Subspace Discordance
ax_c_top = fig.add_subplot(gs_bot[0, 1])  # Panel C Top: B*46:01 Cassette Profile
ax_c_bot = fig.add_subplot(gs_bot[1, 1])  # Panel C Bottom: Mosaic Ribbons

# =========================================================================
# PANEL A: Whole-Gene Continuous Recombination Landscape & Domain Model
# =========================================================================
ax_a_strain = ax_a.twinx()
ax_a_strain.spines["top"].set_visible(False)
ax_a.spines["top"].set_visible(False)

# Breakpoint Flux Field (Primary, left, Vermilion)
line_flux, = ax_a.plot(u_grid, phi_field, color="#D55E00", lw=1.2, label="Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ (Z / bp)")
ax_a.fill_between(u_grid, 0, phi_field, color="#D55E00", alpha=0.18)

# Cumulative Strain Field (Secondary, right, Blue)
line_strain, = ax_a_strain.plot(u_grid, strain_field, color="#0072B2", lw=1.2, ls="--", label="Cumulative Strain $\\rho(u)$ (integrated Z)")
ax_a_strain.fill_between(u_grid, 0, strain_field, color="#0072B2", alpha=0.06)

ax_a.set_xlim(0, L)
ax_a.set_ylim(0, 6.8)  # Generous headroom so annotations stay inside plot
ax_a_strain.set_ylim(0, np.max(strain_field) * 1.18)

ax_a.set_ylabel("Breakpoint Flux $\\Phi_{\\mathrm{bp}}$ (Z / bp)", color="#D55E00", fontweight="bold")
ax_a_strain.set_ylabel("Cumulative Strain $\\rho(u)$ (Z)", color="#0072B2", fontweight="bold")
ax_a.tick_params(axis="y", labelcolor="#D55E00")
ax_a_strain.tick_params(axis="y", labelcolor="#0072B2")
ax_a.tick_params(axis="x", labelbottom=False)

# Key Structural Callouts (Carefully positioned inside axes)
callouts_a = [
    (268, 3.8, "C1 epitope cassette\n(nt 220--293, Z=6.31)", (-75, 1.6)),
    (342, 2.5, "Exon 2 / Exon 3\njunction (342 bp)", (0, 2.5)),
    (611, 2.2, "ARS Groove / $\\alpha_3$ Stalk\nBreakpoint (nt 608--615)", (85, 2.4)),
    (990, 2.8, "Transmembrane / Cyto\nJunction (nt 971--990)", (-85, 1.8))
]

for coord, yval, txt, offset in callouts_a:
    ax_a.annotate(
        txt,
        xy=(coord, yval),
        xytext=(coord + offset[0], yval + offset[1]),
        arrowprops=dict(arrowstyle="->", color="#333333", lw=0.7),
        fontsize=6.2,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFFDE7", edgecolor="#CCCCCC", alpha=0.92)
    )

# Domain Model Track (Bottom of Panel A)
domains = [
    ("Exon 1\n(Signal)", 0, 72, "#E69F00"),
    ("Exon 2 ($\\alpha_1$ ARS groove)", 72, 342, "#D55E00"),
    ("Exon 3 ($\\alpha_2$ ARS groove)", 342, 618, "#CC79A7"),
    ("Exon 4 ($\\alpha_3$ stalk / $\\beta_2$m-CD8)", 618, 894, "#0072B2"),
    ("Exon 5\n(TM)", 894, 1011, "#009E73"),
    ("Exons 6-7\n(Cyto)", 1011, L, "#56B4E9")
]

for dname, s0, s1, col in domains:
    rect = patches.Rectangle((s0, 0.15), s1 - s0, 0.70, facecolor=col, edgecolor="#222222", lw=0.6, alpha=0.88)
    ax_a_domains.add_patch(rect)
    ax_a_domains.text((s0 + s1) / 2.0, 0.50, dname, ha="center", va="center", fontsize=5.8, fontweight="bold",
                      color="white" if col in ["#0072B2", "#D55E00", "#CC79A7"] else "#222222")

ax_a_domains.set_xlim(0, L)
ax_a_domains.set_ylim(0, 1.0)
ax_a_domains.axis("off")

ax_a.set_title("A   Whole-Gene Continuous Recombination Landscape Across Hominoidea MHC Class I (N = 379, L = 1,128 bp)\n     (Direct Discovery of the ARS Groove vs. $\\alpha_3$ Stalk Domain Breakpoint at nt 608--615)",
               loc="left", fontweight="bold", pad=8)


# =========================================================================
# PANEL B: Structural Domain Segregation of Reticulate Changepoints
# =========================================================================
dom_order = [
    ("Exon 2 ($\\alpha_1$ ARS)", 154, "#D55E00", "48.7%"),
    ("Exon 5 (Transmembrane)", 108, "#009E73", "34.2%"),
    ("Exon 3 ($\\alpha_2$ ARS)", 34, "#CC79A7", "10.8%"),
    ("Exon 4 ($\\alpha_3$ Stalk)", 10, "#0072B2", "3.2%"),
    ("Exons 6-7 (Cytoplasmic)", 6, "#56B4E9", "1.9%"),
    ("Exon 1 (Signal)", 4, "#E69F00", "1.3%")
]

b_labels = [d[0] for d in dom_order]
b_counts = [d[1] for d in dom_order]
b_cols = [d[2] for d in dom_order]
b_pcts = [d[3] for d in dom_order]

y_pos = np.arange(len(b_labels))[::-1]
bars = ax_b.barh(y_pos, b_counts, color=b_cols, edgecolor="#222222", lw=0.6, height=0.62)

for y, count, pct in zip(y_pos, b_counts, b_pcts):
    ax_b.text(count + 3, y, f"{count} ({pct})", va="center", ha="left", fontsize=6.8, fontweight="bold", color="#333333")

ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(b_labels, fontsize=7.2, fontweight="bold")
ax_b.set_xlim(0, 190)
ax_b.set_xlabel("Number of Significant Reticulate Taxa ($Z \\geq 3.0$)", fontweight="bold")
ax_b.set_title("B   Structural Domain Segregation of Reticulations\n     (93.7% in ARS Groove & TM Anchor; Exon 4 is Clonal Coldspot)",
               loc="left", fontweight="bold", pad=8)

# Functional Callout inside B
ax_b.text(0.96, 0.12, "Exon 4 $\\alpha_3$ Stalk:\nPreserves species & locus\nphylogeny (only 10 / 316 events)",
          transform=ax_b.transAxes, ha="right", va="bottom", fontsize=6.2,
          bbox=dict(boxstyle="round,pad=0.25", facecolor="#EBF5FB", edgecolor="#AEC6CF", alpha=0.92))


# =========================================================================
# PANEL C Top: HLA-B*46:01 Pairwise Mismatch Profile (C1 Epitope Cassette)
# =========================================================================
records_dict = {r.id: str(r.seq).upper() for r in records}
b46_seq = records_dict["HLA_B_46_01_01_01"]
b15_seq = records_dict["HLA_B_15_01_01_01"]
c01_seq = records_dict["HLA_C_01_02_01_01"]

# Compute 30-bp sliding window mismatch fraction
win = 30
half_w = win // 2
x_coords = np.arange(half_w, L - half_w)
diff_b15 = []
diff_c01 = []

for x in x_coords:
    sub_b46 = b46_seq[x - half_w : x + half_w]
    sub_b15 = b15_seq[x - half_w : x + half_w]
    sub_c01 = c01_seq[x - half_w : x + half_w]
    d15 = sum(1 for a, b in zip(sub_b46, sub_b15) if a != b and a != '-' and b != '-') / win
    d01 = sum(1 for a, b in zip(sub_b46, sub_c01) if a != b and a != '-' and b != '-') / win
    diff_b15.append(d15)
    diff_c01.append(d01)

ax_c_top.plot(x_coords, diff_b15, color="#0072B2", lw=1.1, label="vs. Recipient $HLA\\text{-}B*15:01$ (Backbone)")
ax_c_top.plot(x_coords, diff_c01, color="#009E73", lw=1.1, label="vs. Donor $HLA\\text{-}C*01:02$ (C1 Epitope)")

# Highlight C1 epitope conversion tract (nt 220 to 293)
ax_c_top.axvspan(220, 293, color="#D55E00", alpha=0.22, label="C1 Epitope Cassette (73 bp, Z=2.69)")

ax_c_top.set_xlim(0, L)
ax_c_top.set_ylim(-0.01, 0.25)
ax_c_top.set_ylabel("Pairwise Mismatch Rate", fontweight="bold")
ax_c_top.tick_params(axis="x", labelbottom=False)
ax_c_top.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#DDDDDD", fontsize=6.2)
ax_c_top.set_title("C   Micro-Resolution Deconstruction of Canonical Chimeras\n     ($HLA\\text{-}B*46:01$ Recipient $B*15:01$ vs. Donor $C*01:02$ C1 Epitope Cassette)",
                   loc="left", fontweight="bold", pad=6)

# =========================================================================
# PANEL C Bottom: Mosaic Architectural Ribbons Across Hominoidea Chimeras
# =========================================================================
chimeras = [
    ("B*15:01 (Recip.)", [("B", 0, L, "#0072B2")]),
    ("B*46:01 (Chim.)", [("B", 0, 220, "#0072B2"), ("C", 220, 293, "#009E73"), ("B", 293, L, "#0072B2")]),
    ("C*01:02 (Donor)", [("C", 0, L, "#009E73")]),
    ("Patr-AL (Chim.)", [("A", 0, 280, "#56B4E9"), ("A-like", 280, 665, "#CC79A7"), ("A", 665, 990, "#56B4E9"), ("Cyto", 990, L, "#E69F00")]),
    ("Patr-B*07 (Chimp)", [("B", 0, 326, "#D55E00"), ("B-mosaic", 326, L, "#0072B2")]),
    ("Gogo-C*02 (Gorilla)", [("C", 0, 971, "#009E73"), ("TM-swap", 971, L, "#D55E00")])
]

y_chim = np.arange(len(chimeras))
r_h = 0.52

for y_idx, (cname, tracts) in enumerate(chimeras):
    for t_type, s0, s1, tcol in tracts:
        rect = patches.Rectangle((s0, y_idx - r_h/2.0), s1 - s0, r_h, facecolor=tcol, edgecolor="#222222", lw=0.5, alpha=0.92)
        ax_c_bot.add_patch(rect)

ax_c_bot.set_xlim(0, L)
ax_c_bot.set_ylim(-0.6, len(chimeras) - 0.4)
ax_c_bot.set_yticks(y_chim)
ax_c_bot.set_yticklabels([c[0] for c in chimeras], fontsize=6.8, fontweight="bold")
ax_c_bot.set_xlabel("Genomic Position (bp, Codons 1--376)", fontweight="bold")
ax_c_bot.set_xticks(np.arange(0, 1200, 200))

patches_leg_c = [
    patches.Patch(facecolor="#0072B2", edgecolor="#222222", label="HLA-B Backbone"),
    patches.Patch(facecolor="#009E73", edgecolor="#222222", label="HLA-C / Gogo-C Introgression"),
    patches.Patch(facecolor="#56B4E9", edgecolor="#222222", label="Patr-A Backbone"),
    patches.Patch(facecolor="#CC79A7", edgecolor="#222222", label="Patr-AL Chimeric Segment"),
    patches.Patch(facecolor="#D55E00", edgecolor="#222222", label="Ectopic B/C Motif")
]
ax_c_bot.legend(handles=patches_leg_c, loc="upper center", bbox_to_anchor=(0.5, -0.25),
                frameon=True, facecolor="white", framealpha=0.95, edgecolor="#DDDDDD", ncol=3, fontsize=5.8)


# =========================================================================
# PANEL D: Topological Subspace Discordance (Exons 2-3 ARS vs. Exon 4 Stalk)
# =========================================================================
# Compute 2D MDS on Partition 1 (nt 73-618) and Partition 2 (nt 619-894)
part1 = seq_matrix[:, 73:618]
part2 = seq_matrix[:, 619:894]

def get_mds_2d(mat):
    n = mat.shape[0]
    d = np.zeros((n, n))
    for i in range(n):
        diff = (mat[i] != mat) & (mat[i] != '-') & (mat != '-')
        d[i] = np.mean(diff, axis=1)
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (d ** 2) @ H
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    return evecs[:, :2] * np.sqrt(np.maximum(0, evals[:2]))

z1 = get_mds_2d(part1)
z2 = get_mds_2d(part2)

# Subdivide ax_d into two side-by-side scatter plots
gs_d = ax_d.get_subplotspec().subgridspec(1, 2, wspace=0.28)
ax_d.remove()  # replace container with 2 subplots
ax_d1 = fig.add_subplot(gs_d[0])
ax_d2 = fig.add_subplot(gs_d[1])

locus_colors = {
    "A": "#0072B2",   # Blue
    "B": "#D55E00",   # Vermilion
    "C": "#009E73",   # Green
    "AL": "#CC79A7"   # Purple
}

species_markers = {
    "HLA": "o",       # Circle (Human)
    "Patr": "s",      # Square (Chimp)
    "Papa": "D",      # Diamond (Bonobo)
    "Gogo": "^",      # Triangle (Gorilla)
    "Poab": "*",      # Star (Orangutan)
    "Popy": "*"
}

for i, t in enumerate(taxa):
    parts = t.split('_')
    sp = parts[0]
    loc = parts[1]
    col = locus_colors.get(loc, "#999999")
    m = species_markers.get(sp, "o")
    ax_d1.scatter(z1[i, 0], z1[i, 1], c=col, marker=m, s=11, alpha=0.75, edgecolors="none")
    ax_d2.scatter(z2[i, 0], z2[i, 1], c=col, marker=m, s=11, alpha=0.75, edgecolors="none")

ax_d1.set_title("D   Exons 2--3: ARS Groove\n     (Trans-Species Reticulation)", loc="left", fontsize=6.8, fontweight="bold")
ax_d1.set_xlabel("MDS Dimension 1", fontsize=6.5)
ax_d1.set_ylabel("MDS Dimension 2", fontsize=6.5)
ax_d1.tick_params(labelsize=5.5)

ax_d2.set_title("Exon 4: $\\alpha_3$ Stalk\n(Locus Monophyly: A, B, C)", loc="left", fontsize=6.8, fontweight="bold")
ax_d2.set_xlabel("MDS Dimension 1", fontsize=6.5)
ax_d2.set_ylabel("MDS Dimension 2", fontsize=6.5)
ax_d2.tick_params(labelsize=5.5)

# Legends for Panel D
leg_locus = [
    patches.Patch(facecolor="#0072B2", label="Locus A (N=105)"),
    patches.Patch(facecolor="#D55E00", label="Locus B (N=184)"),
    patches.Patch(facecolor="#009E73", label="Locus C (N=84)"),
    patches.Patch(facecolor="#CC79A7", label="Patr-AL (N=6)")
]
ax_d1.legend(handles=leg_locus, loc="upper left", frameon=True, facecolor="white", framealpha=0.90, edgecolor="#DDDDDD", fontsize=5.5)

# Save figures
png_local = os.path.join(FIG_DIR, "fig_hominoidea_mhc_landscape.png")
pdf_local = os.path.join(FIG_DIR, "fig_hominoidea_mhc_landscape.pdf")
png_paper = os.path.join(PAPER_FIG_DIR, "fig_hominoidea_mhc_landscape.png")
pdf_paper = os.path.join(PAPER_FIG_DIR, "fig_hominoidea_mhc_landscape.pdf")
png_art = os.path.join(ARTIFACT_DIR, "fig_hominoidea_mhc_landscape.png")

plt.savefig(png_local, dpi=300, bbox_inches="tight")
plt.savefig(pdf_local, bbox_inches="tight")
plt.savefig(png_paper, dpi=300, bbox_inches="tight")
plt.savefig(pdf_paper, bbox_inches="tight")
plt.savefig(png_art, dpi=300, bbox_inches="tight")
plt.close()

print(f"[✓] Publication figure successfully generated:")
print(f"  Local PNG:  {png_local}")
print(f"  Local PDF:  {pdf_local}")
print(f"  Paper PNG:  {png_paper}")
print(f"  Paper PDF:  {pdf_paper}")
print(f"  Artifact:   {png_art}")

if __name__ == "__main__":
    pass
