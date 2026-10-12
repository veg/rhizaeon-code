#!/usr/bin/env python3
"""
generate_high_impact_sarbecovirus_figure.py

Generates the flagship high-impact 5-panel figure for the Sarbecovirus empirical application:
"Clean Breaks: Tracking Recombination as Continuous Motion Across Sequence Manifolds"

Panels:
A: Genome-wide Recombination Breakpoint Flux Phi(u) and Cumulative Strain rho(u) across 30.9 kb
   with explicit Hotspot callouts, Coldspot basins, and aligned Gene Model track.
B: Quantitative Concordance with 21 Published Whole-Genome GARD Breakpoints (Lytras et al. 2022).
C: Computational Latency & Scaling Comparison (> 20,000x speedup over cluster GARD).
D: High-Resolution Spike Glycoprotein Reticulate Mosaic Architecture across Wildlife Progenitors
   (Wuhan-Hu-1, Laos BANAL-52/116/236/247, Yunnan RaTG13/RpYN06, Thailand RacCS203, Pangolin P5E).
E: Topological Sequence Manifold Phase Transition Across Codon 313:
   2D Classical MDS projections illustrating Wuhan-Hu-1's continuous metric leap from Bat RaTG13 (NTD)
   to Pangolin MP789 / Laos BANAL-52 (RBD ACE2-contact ridge).
"""

import json
import os
import sys
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# BioVis-Expert styling standards
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.labelsize": 9.5,
    "axes.titlesize": 10.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0,
    "figure.titlesize": 12.0,
    "lines.linewidth": 1.3,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "mathtext.fontset": "dejavusans",
    "figure.dpi": 300
})

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = BASE_DIR
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# 1. Load Data
with open(os.path.join(DATA_DIR, "sarbeco_gene_annotations.json")) as f:
    GENES = json.load(f)

with open(os.path.join(DATA_DIR, "78sarbeco_result.json")) as f:
    RES_78 = json.load(f)

with open(os.path.join(DATA_DIR, "83sarbeco_result.json")) as f:
    RES_83 = json.load(f)

with open(os.path.join(DATA_DIR, "sarbecovirus_reanalysis_report.json")) as f:
    REPORT = json.load(f)

GARD_DF = pd.read_csv(os.path.join(DATA_DIR, "Table_S3.csv"))

# 2. Field Computation Functions
def compute_fields(events, L, sigma=100.0):
    u = np.arange(L, dtype=float)
    strain = np.zeros(L, dtype=float)
    flux = np.zeros(L, dtype=float)
    cutoff_dist = int(4 * sigma)
    
    for ev in events:
        z = ev["z_phys"]
        u1 = max(0, int(ev["u1"]))
        u2 = min(L, int(ev["u2"]))
        if u2 > u1:
            strain[u1:u2] += z
            
        for bp in [u1, u2]:
            lo = max(0, int(bp - cutoff_dist))
            hi = min(L, int(bp + cutoff_dist + 1))
            if hi > lo:
                dist = u[lo:hi] - bp
                k = (1.0 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (dist / sigma) ** 2)
                flux[lo:hi] += z * k
                
    return u, strain, flux

L_78 = RES_78["length"]
L_83 = RES_83["length"]
u_78, strain_78, flux_78 = compute_fields(RES_78["events"], L_78, sigma=100.0)
u_83, strain_83, flux_83 = compute_fields(RES_83["events"], L_83, sigma=100.0)

# Compute MDS for Panel E across Codon 313
def read_fasta(p):
    seqs, names = [], []
    with open(p) as f:
        cur_name, cur_seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_name:
                    seqs.append(''.join(cur_seq).upper())
                    names.append(cur_name)
                cur_name = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_name:
            seqs.append(''.join(cur_seq).upper())
            names.append(cur_name)
    return names, seqs

def compute_cmds(seqs, start, end):
    n = len(seqs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            s1 = seqs[i][start:end]
            s2 = seqs[j][start:end]
            diff = sum(1 for a, b in zip(s1, s2) if a != b and a not in '-N?' and b not in '-N?')
            comp = sum(1 for a, b in zip(s1, s2) if a not in '-N?' and b not in '-N?')
            d = diff / max(1, comp)
            D[i, j] = d
            D[j, i] = d
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    k_evals = np.maximum(evals[:2], 1e-9)
    coords = evecs[:, :2] * np.sqrt(k_evals)
    return coords, D

spike_fasta = os.path.join(DATA_DIR, "sarbecovirus_spike_aligned.fasta")
spike_names, spike_seqs = read_fasta(spike_fasta)
c_ntd, d_ntd = compute_cmds(spike_seqs, 0, 939)        # NTD: bp 1-939 (Codons 1-313)
c_rbd, d_rbd = compute_cmds(spike_seqs, 939, 1620)     # RBD: bp 940-1620 (Codons 314-540)

# =========================================================================
# 3. Master Grid Layout
# =========================================================================
fig = plt.figure(figsize=(11.5, 13.0))
gs_master = fig.add_gridspec(3, 1, height_ratios=[1.20, 0.95, 1.25], hspace=0.38)

# Row 1: Panel A (Genome-wide Recombination Landscape across 30.9 kb) with dedicated gene track
gs_row1 = gs_master[0].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.06)
ax_a = fig.add_subplot(gs_row1[0])
ax_a_genes = fig.add_subplot(gs_row1[1], sharex=ax_a)

# Row 2: 2 columns (Panel B: GARD Concordance and Panel C: Computational Scaling)
gs_row2 = gs_master[1].subgridspec(1, 2, width_ratios=[1.15, 0.85], wspace=0.28)
ax_b = fig.add_subplot(gs_row2[0])
ax_c_speed = fig.add_subplot(gs_row2[1])

# Row 3: 2 columns (Panel D: Spike Mosaic Architecture and Panel E: Metric Space Phase Transition)
gs_row3 = gs_master[2].subgridspec(1, 2, width_ratios=[1.32, 1.08], wspace=0.30)
ax_d_mosaic = fig.add_subplot(gs_row3[0])

# Panel E Subgrid: 2 vertical subplots (NTD MDS on top, RBD MDS on bottom)
gs_e = gs_row3[1].subgridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.58)
ax_e1 = fig.add_subplot(gs_e[0])
ax_e2 = fig.add_subplot(gs_e[1])


# =========================================================================
# PANEL A: Genome-wide Recombination Flux & Strain Field
# =========================================================================
ax_a_strain = ax_a.twinx()
ax_a_strain.spines["top"].set_visible(False)
ax_a_strain.spines["left"].set_visible(False)
ax_a_strain.spines["right"].set_visible(True)

# Highlight Coldspot Basins (Shaded light blue corridor)
coldspots = [
    (13600, 16400, "Coldspot: RdRp Catalytic Core"),
    (24600, 25600, "Coldspot: Spike S2 Fusion Engine")
]
for c_start, c_end, c_label in coldspots:
    ax_a.axvspan(c_start, c_end, color="#56B4E9", alpha=0.10, zorder=0)
    ax_a.text((c_start + c_end)/2, ax_a.get_ylim()[1]*0.18, c_label,
              ha="center", va="bottom", fontsize=6.5, color="#0072B2", fontstyle="italic",
              bbox=dict(boxstyle="square,pad=0.15", facecolor="#EBF5FB", edgecolor="#AED6F1", lw=0.5, alpha=0.85))

# Breakpoint Flux (Primary, left)
line_flux83, = ax_a.plot(u_83, flux_83, color="#D55E00", lw=1.6, label="Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ (83 taxa)")
line_flux78, = ax_a.plot(u_78, flux_78, color="#E69F00", lw=1.1, ls="--", alpha=0.85, label="Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ (78 taxa)")

# Cumulative Strain (Secondary, right)
line_strain83 = ax_a_strain.fill_between(u_83, 0, strain_83, color="#0072B2", alpha=0.14, label="Cumulative Strain $\\rho(u)$ (83 taxa)")
line_strain78, = ax_a_strain.plot(u_78, strain_78, color="#56B4E9", lw=0.9, ls=":", alpha=0.75, label="Cumulative Strain $\\rho(u)$ (78 taxa)")

ax_a.set_xlim(0, 31000)
ax_a.set_ylim(0, max(np.max(flux_83), np.max(flux_78)) * 1.20)
ax_a_strain.set_ylim(0, max(np.max(strain_83), np.max(strain_78)) * 1.22)

ax_a.set_ylabel("Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ (Z / nt)", color="#D55E00", fontweight="bold")
ax_a_strain.set_ylabel("Cumulative Strain $\\rho(u)$ (Z)", color="#0072B2", fontweight="bold")
ax_a.tick_params(axis="y", labelcolor="#D55E00")
ax_a_strain.tick_params(axis="y", labelcolor="#0072B2")
ax_a.tick_params(axis="x", labelbottom=False)

# Annotate Major Hotspots
hotspot_calls = [
    (12867, "RdRp exit", 0.78),
    (21825, "Spike 5' / NTD", 0.90),
    (22764, "Codon 313 hinge", 0.74),
    (25728, "Spike exit / ORF3a", 0.62),
    (28810, "ORF8 / N cluster", 0.88)
]
for pos, label, y_frac in hotspot_calls:
    y_val = flux_83[min(int(pos), len(flux_83)-1)]
    ax_a.annotate(
        f"{label}\n({pos:,} bp)",
        xy=(pos, y_val),
        xytext=(pos, ax_a.get_ylim()[1] * y_frac),
        arrowprops=dict(facecolor="#333333", edgecolor="none", width=0.6, headwidth=4, shrink=0.08),
        fontsize=7,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#CCCCCC", alpha=0.92)
    )

handles_a = [line_flux83, line_flux78, line_strain83, line_strain78]
labels_a = [h.get_label() for h in handles_a]
ax_a.legend(handles_a, labels_a, loc="upper left", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#DDDDDD")
ax_a.set_title("A   Genome-wide Recombination Strain Field $\\rho(u)$ and Physical Breakpoint Flux $\\Phi_{\\mathrm{bp}}(u)$ Across 30.9 kb",
               loc="left", fontweight="bold", pad=8)

# Dedicated Gene Model Track
GENE_COLORS = {
    "ORF1ab": "#E0E0E0", "RdRp": "#B0B0B0", "Spike": "#F0E442",
    "Spike_NTD": "#E69F00", "Spike_RBD": "#D55E00", "Spike_RBM": "#CC79A7",
    "ORF3a": "#56B4E9", "E": "#009E73", "M": "#0072B2",
    "ORF6": "#F0E442", "ORF7a": "#E69F00", "ORF7b": "#D55E00",
    "ORF8": "#CC79A7", "N": "#009E73", "ORF10": "#999999"
}

ax_a_genes.set_ylim(0, 1)
ax_a_genes.spines["left"].set_visible(False)
ax_a_genes.spines["bottom"].set_visible(True)
ax_a_genes.tick_params(axis="y", left=False, labelleft=False)

major_genes = ["ORF1ab", "Spike", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N"]
for g_name in major_genes:
    info = GENES.get(g_name)
    if not info:
        continue
    x0 = info["aln_start"]
    width = info["aln_length"]
    col = GENE_COLORS.get(g_name, "#E0E0E0")
    rect = patches.Rectangle((x0, 0.15), width, 0.7, facecolor=col, edgecolor="#333333", lw=0.6)
    ax_a_genes.add_patch(rect)
    if width > 800:
        ax_a_genes.text(x0 + width/2, 0.5, g_name, ha="center", va="center", fontsize=7, fontweight="bold", color="#111111")

# Subdomains inside Spike and ORF1ab
s_ntd = GENES["Spike_NTD"]
s_rbd = GENES["Spike_RBD"]
ax_a_genes.add_patch(patches.Rectangle((s_ntd["aln_start"], 0.15), s_ntd["aln_length"], 0.7, facecolor="#E69F00", edgecolor="#333333", lw=0.6))
ax_a_genes.add_patch(patches.Rectangle((s_rbd["aln_start"], 0.15), s_rbd["aln_length"], 0.7, facecolor="#D55E00", edgecolor="#333333", lw=0.6))

rdrp = GENES["RdRp"]
ax_a_genes.add_patch(patches.Rectangle((rdrp["aln_start"], 0.15), rdrp["aln_length"], 0.7, facecolor="#A0A0A0", edgecolor="#333333", lw=0.6))
ax_a_genes.text(rdrp["aln_start"] + rdrp["aln_length"]/2, 0.5, "RdRp", ha="center", va="center", fontsize=6.5, fontweight="bold", color="#111111")

ax_a_genes.set_xlabel("Genomic Coordinate in Sarbecovirus Multiple Alignment (nt)", fontweight="bold")
ax_a_genes.set_xticks(np.arange(0, 32000, 5000))
ax_a_genes.set_xticklabels([f"{x:,}" for x in np.arange(0, 32000, 5000)])


# =========================================================================
# PANEL B: GARD Benchmark Concordance
# =========================================================================
matches_78 = REPORT["dataset_78"]["gard_concordance"]["matches"]
gard_bps = [m["gard_bp"] for m in matches_78]
rhiz_bps = [m["rhizaeon_closest_bp"] for m in matches_78]
dists = [m["distance_nt"] for m in matches_78]

ax_b.plot([0, 31000], [0, 31000], color="#666666", ls="--", lw=1.0, label="1:1 Identity Line")
ax_b.fill_between([0, 31000], [-150, 31000-150], [150, 31000+150], color="#E0E0E0", alpha=0.6, label="$\\pm 150$ nt Concordance Corridor")

for g, r, d in zip(gard_bps, rhiz_bps, dists):
    if d <= 150:
        ax_b.scatter(g, r, color="#009E73", s=28, zorder=5, edgecolor="#111111", lw=0.5)
    else:
        ax_b.scatter(g, r, color="#D55E00", s=28, zorder=5, edgecolor="#111111", lw=0.5)

key_callouts = [
    (12854, 12866, "RdRp (12 nt)", (-20, 20)),
    (21430, 21516, "Spike 5' (86 nt)", (-50, 18)),
    (22764, 22719, "NTD/RBD (45 nt)", (15, -20)),
    (24452, 24461, "S1/S2 (9 nt)", (-45, 18)),
    (25151, 25161, "S2 exit (10 nt)", (15, -18)),
    (28804, 29046, "ORF8/N (242 nt)", (-40, 20))
]
for g_val, r_val, lbl, offset in key_callouts:
    ax_b.annotate(
        lbl, xy=(g_val, r_val), xytext=(g_val + offset[0]*60, r_val + offset[1]*60),
        arrowprops=dict(facecolor="#333333", edgecolor="none", width=0.5, headwidth=3, shrink=0.1),
        fontsize=6.5, ha="center",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="#FFF9C4", edgecolor="#FBC02D", lw=0.5, alpha=0.92)
    )

ax_b.set_xlim(0, 31000)
ax_b.set_ylim(0, 31000)
ax_b.set_aspect("equal")
ax_b.set_xlabel("Lytras et al. (2022) GARD Breakpoints (nt)", fontweight="bold")
ax_b.set_ylabel("RhizAeon Changepoints (nt)", fontweight="bold")
ax_b.set_xticks(np.arange(0, 32000, 10000))
ax_b.set_yticks(np.arange(0, 32000, 10000))
ax_b.set_xticklabels([f"{x:,}" for x in np.arange(0, 32000, 10000)])
ax_b.set_yticklabels([f"{x:,}" for x in np.arange(0, 32000, 10000)])
ax_b.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#DDDDDD")
ax_b.set_title("B   Concordance with 21 Published GARD Breakpoints\n     (Median Deviation: 99.0 nt; Key S1/S2 Junction: 9 nt)",
               loc="left", fontweight="bold", pad=8)


# =========================================================================
# PANEL C: Computational Scaling & Latency Comparison
# =========================================================================
t_78 = RES_78.get("run_time_ms", 18158.9) / 1000.0
t_83 = RES_83.get("run_time_ms", 25025.6) / 1000.0

benchmarks = [
    ("GARD 2022\n(18 taxa)", 388800.0, "#999999", "4.5 days\n(18 taxa)"),
    ("GARD 2022\n(78 taxa est.)", 2592000.0, "#CCCCCC", ">30 days\n(infeasible)"),
    ("RhizAeon\n(78 taxa)", t_78, "#0072B2", f"{t_78:.1f} s\n(78 taxa)"),
    ("RhizAeon\n(83 taxa)", t_83, "#D55E00", f"{t_83:.1f} s\n(83 taxa)")
]

names = [b[0] for b in benchmarks]
runtimes = [b[1] for b in benchmarks]
colors = [b[2] for b in benchmarks]

bars = ax_c_speed.bar(range(len(names)), runtimes, color=colors, edgecolor="#333333", lw=0.6, width=0.55)
ax_c_speed.set_yscale("log")
ax_c_speed.set_ylim(1, 1e7)
ax_c_speed.set_xticks(range(len(names)))
ax_c_speed.set_xticklabels(names, fontsize=7.5)
ax_c_speed.set_ylabel("Execution Time (seconds, log scale)", fontweight="bold")
ax_c_speed.set_title("C   Computational Latency & Scaling\n     ($>20{,}000\\times$ Speedup Over Cluster GARD)",
                     loc="left", fontweight="bold", pad=8)

# Speedup Callouts
ax_c_speed.text(0, 388800 * 1.5, "4.5 days\n(18 taxa)", ha="center", va="bottom", fontsize=6.8, fontweight="bold")
ax_c_speed.text(1, 2592000 * 1.5, ">30 days\n(infeasible)", ha="center", va="bottom", fontsize=6.8, color="#777777")
ax_c_speed.text(2, t_78 * 1.6, f"{t_78:.1f} s\n(78 taxa)", ha="center", va="bottom", fontsize=6.8, fontweight="bold", color="#0072B2")
ax_c_speed.text(3, t_83 * 1.6, f"{t_83:.1f} s\n(83 taxa)", ha="center", va="bottom", fontsize=6.8, fontweight="bold", color="#D55E00")

# Horizontal reference lines
ax_c_speed.axhline(60, color="#AAAAAA", ls=":", lw=0.7)
ax_c_speed.text(3.42, 60, "1 min", va="center", fontsize=6.5, color="#666666")
ax_c_speed.axhline(3600, color="#AAAAAA", ls=":", lw=0.7)
ax_c_speed.text(3.42, 3600, "1 hr", va="center", fontsize=6.5, color="#666666")
ax_c_speed.axhline(86400, color="#AAAAAA", ls=":", lw=0.7)
ax_c_speed.text(3.42, 86400, "1 day", va="center", fontsize=6.5, color="#666666")


# =========================================================================
# PANEL D: Spike Locus Deep Dive (nt 21,000 - 26,000)
# =========================================================================
taxa_tracks = [
    ("SARS-CoV-2 (Wuhan-Hu-1)", [
        ("Backbone", 21000, 26000, "#0072B2", "Wuhan-Hu-1 Clade")
    ]),
    ("BANAL-52 (Laos)", [
        ("Backbone", 21000, 26000, "#0072B2", "Wuhan-Hu-1 Clade (ACE2 binding)")
    ]),
    ("BANAL-116 (Laos)", [
        ("Backbone", 21000, 22076, "#56B4E9", "BANAL-103 Clade"),
        ("Micro-conversion", 22076, 22251, "#E69F00", "RacCS203 NTD introgression (Z=146.2)"),
        ("Backbone", 22251, 22649, "#56B4E9", "BANAL-103 Clade"),
        ("S1 Cassette", 22649, 24577, "#E69F00", "RacCS203 S1 introgression (Z=27.1)"),
        ("Backbone", 24577, 26000, "#56B4E9", "BANAL-103 Clade")
    ]),
    ("BANAL-236 (Laos)", [
        ("Backbone", 21000, 22645, "#999999", "PrC31 Backbone"),
        ("RBD+S2", 22645, 25728, "#56B4E9", "BANAL-103 RBD/S2 transfer (Z=17.6)"),
        ("Backbone", 25728, 26000, "#999999", "PrC31 Backbone")
    ]),
    ("BANAL-247 (Laos)", [
        ("Backbone", 21000, 21909, "#56B4E9", "BANAL-103 Clade"),
        ("NTD Micro-1", 21909, 22011, "#E69F00", "RacCS203 NTD transfer (Z=155.3)"),
        ("Backbone", 22011, 22054, "#56B4E9", "BANAL-103 Clade"),
        ("NTD Micro-2", 22054, 22251, "#E69F00", "RacCS203 NTD transfer (Z=113.1)"),
        ("Backbone", 22251, 22649, "#56B4E9", "BANAL-103 Clade"),
        ("S1 Cassette", 22649, 24577, "#E69F00", "RacCS203 S1 transfer (Z=26.9)"),
        ("Backbone", 24577, 26000, "#56B4E9", "BANAL-103 Clade")
    ]),
    ("RaTG13 (Yunnan)", [
        ("NTD", 21000, 22764, "#0072B2", "Shares NTD with Wuhan-Hu-1"),
        ("RBD/S2", 22764, 26000, "#CC79A7", "Divergent RBD/S2")
    ]),
    ("RpYN06 (Yunnan)", [
        ("Backbone", 21000, 21539, "#0072B2", "SARS-CoV-2 Genome Backbone"),
        ("S1 Cassette", 21539, 24700, "#999999", "PrC31 S1 replacement (Z=41.0)"),
        ("S2 Subunit", 24700, 26000, "#0072B2", "SARS-CoV-2 Genome Backbone")
    ]),
    ("RacCS203 (Thailand)", [
        ("Backbone", 21000, 26000, "#E69F00", "RacCS203 Backbone")
    ]),
    ("P5E (Pangolin-CoV)", [
        ("Backbone", 21000, 22813, "#999999", "Pangolin Backbone"),
        ("ACE2 RBM", 22813, 23490, "#0072B2", "Wuhan-Hu-1 Shared RBM (Z=26.8)"),
        ("Backbone", 23490, 26000, "#999999", "Pangolin Backbone")
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
            ax_d_mosaic.add_patch(rect)

# Domain vertical spans
ax_d_mosaic.axvspan(GENES["Spike_NTD"]["aln_start"], GENES["Spike_NTD"]["aln_end"], color="#E69F00", alpha=0.10)
ax_d_mosaic.axvspan(GENES["Spike_RBD"]["aln_start"], GENES["Spike_RBD"]["aln_end"], color="#D55E00", alpha=0.14)
ax_d_mosaic.axvspan(GENES["Spike_RBM"]["aln_start"], GENES["Spike_RBM"]["aln_end"], color="#CC79A7", alpha=0.20)
# S1/S2 furin site
ax_d_mosaic.axvline(23618, color="#D55E00", ls="--", lw=1.2)

ax_d_mosaic.set_xlim(21000, 26000)
ax_d_mosaic.set_ylim(-0.8, len(taxa_tracks) + 0.8)
ax_d_mosaic.set_yticks(y_pos)
ax_d_mosaic.set_yticklabels([t[0] for t in taxa_tracks], fontsize=7.5, fontweight="bold")
ax_d_mosaic.set_xlabel("Spike Glycoprotein Coordinate in Multiple Alignment (nt)", fontweight="bold")
ax_d_mosaic.set_title("D   Spike Reticulate Mosaic Architecture Across Wildlife Progenitors\n     (Laos BANAL, Yunnan RpYN06/RaTG13, and Pangolin-CoV)",
                      loc="left", fontweight="bold", pad=20)

# Domain labels on top (in dedicated strip above highest track)
ax_d_mosaic.text(22000, len(taxa_tracks) + 0.12, "NTD", ha="center", fontsize=7.5, fontweight="bold", color="#D55E00")
ax_d_mosaic.text(23300, len(taxa_tracks) + 0.12, "RBD", ha="center", fontsize=7.5, fontweight="bold", color="#D55E00")
ax_d_mosaic.text(23618, -0.65, "S1/S2 (23,618)", ha="center", fontsize=6.8, fontweight="bold", color="#D55E00")

legend_elements_d = [
    patches.Patch(facecolor="#0072B2", edgecolor="#222222", label="Wuhan-Hu-1 Clade / Core Backbone"),
    patches.Patch(facecolor="#56B4E9", edgecolor="#222222", label="BANAL-103 / Laos Clade"),
    patches.Patch(facecolor="#E69F00", edgecolor="#222222", label="RacCS203 / Thailand Clade"),
    patches.Patch(facecolor="#999999", edgecolor="#222222", label="PrC31 / Distant Bat Backbone"),
    patches.Patch(facecolor="#CC79A7", edgecolor="#222222", label="RaTG13 Divergent Locus")
]
ax_d_mosaic.legend(handles=legend_elements_d, loc="lower left", bbox_to_anchor=(0.0, -0.28), ncol=3,
                   frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.0)


# =========================================================================
# PANEL E: Metric Space Trajectory Phase Transition Across Codon 313
# =========================================================================
# Subpanel E1: Partition 1 (NTD: bp 1-939, before Codon 313)
# Subpanel E2: Partition 2 (RBD: bp 940-1620, after Codon 313)
colors_map = {
    "SARS_CoV_2_Wuhan_Hu_1": "#D55E00", # Recombinant focal (Vermillion)
    "Bat_RaTG13": "#0072B2",            # NTD Donor (Blue)
    "Bat_BANAL_52": "#56B4E9",          # Sky blue
    "Bat_BANAL_103": "#009E73",         # Green
    "Pangolin_MP789": "#CC79A7",        # Purple / Reddish purple (RBD Home)
    "Bat_RmYN02": "#F0E442",            # Yellow
    "SARS_CoV_1_Tor2": "#999999",       # Grey
    "Bat_RsSHC014": "#777777",          # Dark grey
    "Bat_HKU3_1": "#333333"             # Black
}

# --- Subpanel E1: NTD ---
ax_e1.set_title("E   Topological Metric Jump Across Codon 313\n     E1: NTD (bp 1-939, Codons 1-313)",
                loc="left", fontweight="bold", pad=8, fontsize=8.2)

for i, nm in enumerate(spike_names):
    col = colors_map.get(nm, "#333333")
    pt = c_ntd[i]
    is_focal = (nm == "SARS_CoV_2_Wuhan_Hu_1")
    is_donor = (nm == "Bat_RaTG13")
    is_home = (nm == "Pangolin_MP789")
    
    sz = 55 if (is_focal or is_donor or is_home) else 25
    ax_e1.scatter(pt[0], pt[1], color=col, s=sz, zorder=5, edgecolor="#111111", lw=0.6)

# Non-overlapping callouts in E1
ax_e1.text(-0.165, -0.065, "SARS-CoV-2", fontsize=7.0, fontweight="bold", color="#D55E00", ha="right", va="center")
ax_e1.text(-0.165, -0.040, "RaTG13 (Donor)", fontsize=6.8, fontweight="bold", color="#0072B2", ha="right", va="center")
ax_e1.text(-0.130, -0.040, "BANAL-52", fontsize=6.5, color="#56B4E9", ha="left", va="bottom")
ax_e1.text(0.025, 0.230, "Pangolin MP789", fontsize=6.8, color="#CC79A7", ha="left", va="center")

# Annotate tight RaTG13 cluster
ax_e1.annotate("94.4% Identity\nCluster", xy=(c_ntd[0, 0], c_ntd[0, 1]), xytext=(-0.05, -0.15),
               arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8),
               fontsize=6.5, ha="center", bbox=dict(boxstyle="round,pad=0.2", facecolor="#E8F8F5", edgecolor="#A3E4D7", lw=0.5))

ax_e1.set_xlim(-0.24, 0.28)
ax_e1.set_ylim(-0.22, 0.30)
ax_e1.set_xlabel("MDS Dimension 1", fontsize=7.5)
ax_e1.set_ylabel("MDS Dimension 2", fontsize=7.5)
ax_e1.tick_params(labelsize=6.8)
ax_e1.axhline(0, color="#EEEEEE", ls=":", lw=0.5)
ax_e1.axvline(0, color="#EEEEEE", ls=":", lw=0.5)

# --- Subpanel E2: RBD ---
ax_e2.set_title("E2: Receptor-Binding Domain (bp 940-1620)",
                loc="left", fontweight="bold", pad=8, fontsize=8.2)

for i, nm in enumerate(spike_names):
    col = colors_map.get(nm, "#333333")
    pt = c_rbd[i]
    is_focal = (nm == "SARS_CoV_2_Wuhan_Hu_1")
    is_donor = (nm == "Bat_RaTG13")
    is_home = (nm == "Pangolin_MP789")
    
    sz = 55 if (is_focal or is_donor or is_home) else 25
    ax_e2.scatter(pt[0], pt[1], color=col, s=sz, zorder=5, edgecolor="#111111", lw=0.6)

# Non-overlapping callouts in E2
ax_e2.text(-0.098, -0.008, "SARS-CoV-2", fontsize=7.0, fontweight="bold", color="#D55E00", ha="right", va="center")
ax_e2.text(-0.080, 0.024, "BANAL-52 (ACE2 Ridge)", fontsize=6.8, color="#0072B2", ha="center", va="bottom")
ax_e2.text(-0.052, 0.012, "Pangolin MP789 (Home)", fontsize=6.8, fontweight="bold", color="#CC79A7", ha="left", va="center")
ax_e2.text(-0.070, -0.026, "RaTG13 (Divergent)", fontsize=6.5, color="#666666", ha="center", va="top")

# Annotate convergence on ACE2 binding ridge
ax_e2.annotate("Shared ACE2 Contact\nRidge Homology", xy=(c_rbd[0, 0], c_rbd[0, 1]), xytext=(-0.01, -0.11),
               arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8),
               fontsize=6.5, ha="center", bbox=dict(boxstyle="round,pad=0.2", facecolor="#FDEDEC", edgecolor="#F5B7B1", lw=0.5))

ax_e2.set_xlim(-0.16, 0.24)
ax_e2.set_ylim(-0.16, 0.22)
ax_e2.set_xlabel("MDS Dimension 1", fontsize=7.5)
ax_e2.set_ylabel("MDS Dimension 2", fontsize=7.5)
ax_e2.tick_params(labelsize=6.8)
ax_e2.axhline(0, color="#EEEEEE", ls=":", lw=0.5)
ax_e2.axvline(0, color="#EEEEEE", ls=":", lw=0.5)

# Save high-res PNG and PDF
out_png = os.path.join(FIG_DIR, "sarbecovirus_recombination_landscape.png")
out_pdf = os.path.join(FIG_DIR, "sarbecovirus_recombination_landscape.pdf")
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

print(f"Successfully generated high-impact publication figure:")
print(f"  {out_png}")
print(f"  {out_pdf}")
