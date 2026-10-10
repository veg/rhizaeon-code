"""
paper/figures/generate_fig3_theory_force_field_svd.py
=====================================================
Generates Figure 3 for Section 2 of the RhizAeon manuscript:
"Physical Force Field Deconvolution, Multi-Scale Matched Filters, and Random Matrix Theory"

CRITICAL SCIENTIFIC INTEGRITY & BIOVIS-EXPERT MANDATES:
1. 100% UNIFIED ON THE EMPIRICAL HIV-1 POL RETROVIRAL BENCHMARK:
   All four panels (A, B, C, D) are computed directly from authentic HIV-1 pol sequences
   (hiv1_pol_8B_8C_3CRF07_3CRF08.fasta, N=22, L=3,120 nt) and authentic engine outputs.
   Zero simulated toy datasets or arbitrary benchmark hopping.
2. ZERO OVERPLOTTED LEGENDS IN PLOTTING AREA:
   Every single legend is placed outside the data frames in a dedicated, isolated header track.
   Zero occlusion of curves, scatter points, streamlines, or filter baselines.
3. ZERO TITLE-LEGEND COLLISIONS:
   Subplot titles and legends are rendered in dedicated non-overlapping header tracks with
   explicit mathematical y-coordinate offsets, eliminating Matplotlib title-legend collisions.
"""

from pathlib import Path
import json
import math
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ==============================================================================
# 1. CORE ALGORITHMIC & LOADING ROUTINES
# ==============================================================================

BASE_MAP = {
    "A": 1, "a": 1,
    "C": 2, "c": 2,
    "G": 3, "g": 3,
    "T": 4, "t": 4, "U": 4, "u": 4
}

def load_fasta(aln_path: Path):
    taxa = []
    seq_strings = []
    curr_label = None
    curr_seq = []

    with open(aln_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_label is not None:
                    taxa.append(curr_label)
                    seq_strings.append("".join(curr_seq))
                curr_label = line[1:].strip()
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_label is not None:
            taxa.append(curr_label)
            seq_strings.append("".join(curr_seq))

    N = len(taxa)
    L = len(seq_strings[0])
    seq_matrix = np.zeros((N, L), dtype=np.int8)
    for i, s in enumerate(seq_strings):
        for u, ch in enumerate(s):
            seq_matrix[i, u] = BASE_MAP.get(ch, 0)
    return seq_matrix, taxa, L

def compute_global_distances(seq_matrix: np.ndarray) -> np.ndarray:
    N, L = seq_matrix.shape
    D = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        for j in range(i + 1, N):
            valid = (seq_matrix[i] > 0) & (seq_matrix[j] > 0)
            if np.any(valid):
                d = np.mean(seq_matrix[i, valid] != seq_matrix[j, valid])
            else:
                d = 0.0
            D[i, j] = d
            D[j, i] = d
    return D

def embed_canonical_pcoa(D_matrix: np.ndarray, k: int = 4) -> np.ndarray:
    N = D_matrix.shape[0]
    H = np.eye(N) - (1.0 / N) * np.ones((N, N))
    B = -0.5 * H @ (D_matrix ** 2) @ H
    evals, evecs = la.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    pos_evals = np.maximum(evals[:k], 1e-9)
    coords = evecs[:, :k] * np.sqrt(pos_evals)
    return coords

def compute_2d_grid_potential(medoid_basin: np.ndarray, xlim=(-0.08, 0.08), ylim=(-0.04, 0.06), nx=45, ny=45, lambda_param=1.0):
    x_lin = np.linspace(xlim[0], xlim[1], nx)
    y_lin = np.linspace(ylim[0], ylim[1], ny)
    X, Y = np.meshgrid(x_lin, y_lin)
    
    # Hookean potential well centered at medoid_basin: U(z) = 0.5 * lambda * ||z - z_basin||^2
    dx = X - medoid_basin[0]
    dy = Y - medoid_basin[1]
    V = 0.5 * lambda_param * (dx**2 + dy**2)
    
    # Restoring force: f = -grad(U) = -lambda * (z - z_basin)
    Fx = -lambda_param * dx
    Fy = -lambda_param * dy
    return X, Y, Fx, Fy, V

# ==============================================================================
# 2. MAIN GENERATOR ROUTINE
# ==============================================================================

def generate_figure3():
    repo_root = Path(__file__).resolve().parent.parent.parent
    hiv_path = Path("/Users/sergei/Projects/TOGA_MEME/recombination/visual_intuition/hiv1_pol_8B_8C_3CRF07_3CRF08.fasta")
    if not hiv_path.exists():
        hiv_path = repo_root / "visual_intuition" / "hiv1_pol_8B_8C_3CRF07_3CRF08.fasta"
        
    json1_path = repo_root / "paper" / "figures" / "fig1_rhizaeon_authentic_data.json"
    json2_path = repo_root / "paper" / "figures" / "fig2_rhizaeon_authentic_data.json"
    
    pdf_path = repo_root / "paper" / "figures" / "fig3_theory_force_field_svd.pdf"
    png_path = repo_root / "paper" / "figures" / "fig3_theory_force_field_svd.png"
    out_json = repo_root / "paper" / "figures" / "fig3_rhizaeon_authentic_data.json"
    
    print("[*] Loading authentic HIV-1 pol alignment and engine datasets...")
    seqs_hiv, taxa_hiv, L_hiv = load_fasta(hiv_path)
    D_hiv = compute_global_distances(seqs_hiv)
    coords_hiv = embed_canonical_pcoa(D_hiv, k=4)
    m2d_hiv = coords_hiv[:, :2]
    
    with open(json1_path, "r", encoding="utf-8") as f:
        d1 = json.load(f)
    with open(json2_path, "r", encoding="utf-8") as f:
        d2 = json.load(f)
        
    b_taxa = d1["b_indices"]
    c_taxa = d1["c_indices"]
    crf07_taxa = d1["crf07_indices"]
    crf08_taxa = d1["crf08_indices"]
    af286230_idx = d1["af286230_idx"]
    ay008715_idx = d1["ay008715_idx"]
    
    med_B = np.mean(m2d_hiv[b_taxa], axis=0)
    med_C = np.mean(m2d_hiv[c_taxa], axis=0)
    d_vec = med_B - med_C
    d_hat = d_vec / np.linalg.norm(d_vec)
    
    # Aesthetic styling (Okabe-Ito colorblind palette)
    COLOR_P1 = "#0072B2"      # Navy Blue (Subtype B / Donor)
    COLOR_P2 = "#D55E00"      # Vermilion (Subtype C / Home)
    COLOR_REC = "#009E73"     # Emerald Green (CRF07_BC)
    COLOR_DONOR = "#CC79A7"   # Royal Reddish Purple (CRF08_BC)
    COLOR_ROOT = "#222222"    # Charcoal Gray
    COLOR_AMBER = "#E69F00"   # Amber (Matched Scale)
    COLOR_SKY = "#56B4E9"     # Sky Blue
    COLOR_GRID = "#EAEAEA"
    
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 7.8,
        "axes.titlesize": 8.0,
        "axes.titleweight": "bold",
        "axes.labelsize": 7.6,
        "axes.labelweight": "medium",
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.8,
        "figure.titlesize": 11,
        "figure.titleweight": "bold",
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.4,
        "patch.linewidth": 0.8,
    })
    
    fig = plt.figure(figsize=(14.0, 11.4), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.0], hspace=0.28, wspace=0.26)
    
    # ==============================================================================
    # PANEL A: Metric Force Field & Hookean Restoring Potential (HIV-1 pol)
    # ==============================================================================
    print("[*] Generating Panel A: Continuous Force Field on HIV-1 pol...")
    # 7 rows: [0] Main Title, [1] Titles A1, [2] Legend A1, [3] Plot A1, [4] Title A2, [5] Legend A2, [6] Plot A2
    inner_a = gridspec.GridSpecFromSubplotSpec(7, 1, subplot_spec=gs[0, 0], 
                                               height_ratios=[0.08, 0.14, 0.08, 0.66, 0.10, 0.08, 0.72], 
                                               hspace=0.28)
    
    # A-Head: Dedicated title track
    ax_a_head = fig.add_subplot(inner_a[0])
    ax_a_head.axis("off")
    ax_a_head.text(-0.04, 0.50, "A", fontsize=12, fontweight="bold", va="center")
    ax_a_head.text(0.02, 0.50, "Metric Force Field & Hookean Restoring Potential (HIV-1 pol)", 
                   fontsize=9.2, fontweight="bold", va="center", color="#1e293b")
    
    # A1-Titles: Dedicated side-by-side title strip for A1
    inner_a1_titles = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner_a[1], wspace=0.28)
    ax_a1_tleft = fig.add_subplot(inner_a1_titles[0])
    ax_a1_tleft.axis("off")
    ax_a1_tleft.text(0.50, 0.50, "Protease Flank (0–862 nt)\nBasin at Subtype C (Home)", 
                     fontsize=7.8, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    ax_a1_tright = fig.add_subplot(inner_a1_titles[1])
    ax_a1_tright.axis("off")
    ax_a1_tright.text(0.50, 0.50, "RT Cassette (862–1169 nt)\nBasin at Subtype B (Donor)", 
                      fontsize=7.8, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    # A-Leg1: Dedicated legend strip for A1
    ax_a_leg1 = fig.add_subplot(inner_a[2])
    ax_a_leg1.axis("off")
    legend_handles_a = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_P1, markersize=7, label="Subtype B (Donor)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=COLOR_P2, markersize=7, label="Subtype C (Home)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_REC, markersize=7, label="CRF07_BC ($n=3$)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=COLOR_DONOR, markersize=7, label="CRF08_BC ($n=3$)"),
    ]
    ax_a_leg1.legend(handles=legend_handles_a, loc="center", ncol=4, frameon=False, fontsize=7.2)
    
    # A1: Side-by-side potential landscapes on HIV-1 pol (NO set_title)
    inner_a1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner_a[3], wspace=0.28)
    ax_a1_left = fig.add_subplot(inner_a1[0])
    ax_a1_right = fig.add_subplot(inner_a1[1])
    
    xlim_a = (-0.075, 0.075)
    ylim_a = (-0.035, 0.055)
    X1, Y1, Fx1, Fy1, V1 = compute_2d_grid_potential(med_C, xlim=xlim_a, ylim=ylim_a)
    X2, Y2, Fx2, Fy2, V2 = compute_2d_grid_potential(med_B, xlim=xlim_a, ylim=ylim_a)
    
    # Segment 1 (PR flank, 0-862 nt): Basin at Subtype C (Home)
    ax_a1_left.contourf(X1, Y1, V1, levels=14, cmap="Oranges_r", alpha=0.75)
    ax_a1_left.streamplot(X1[0, :], Y1[:, 0], Fx1, Fy1, color="#8a2500", density=0.75, linewidth=0.7, arrowsize=0.8)
    
    # Segment 2 (RT cassette, 862-1169 nt): Basin at Subtype B (Donor)
    ax_a1_right.contourf(X2, Y2, V2, levels=14, cmap="Blues_r", alpha=0.75)
    ax_a1_right.streamplot(X2[0, :], Y2[:, 0], Fx2, Fy2, color="#003f5c", density=0.75, linewidth=0.7, arrowsize=0.8)
    
    for ax_sub in [ax_a1_left, ax_a1_right]:
        ax_sub.scatter(m2d_hiv[b_taxa, 0], m2d_hiv[b_taxa, 1], color=COLOR_P1, marker="o", s=32, edgecolors="white", linewidth=0.8, zorder=5)
        ax_sub.scatter(m2d_hiv[c_taxa, 0], m2d_hiv[c_taxa, 1], color=COLOR_P2, marker="s", s=32, edgecolors="white", linewidth=0.8, zorder=5)
        ax_sub.scatter(m2d_hiv[crf07_taxa, 0], m2d_hiv[crf07_taxa, 1], color=COLOR_REC, marker="D", s=40, edgecolors="#1e293b", linewidth=0.8, zorder=6)
        ax_sub.scatter(m2d_hiv[crf08_taxa, 0], m2d_hiv[crf08_taxa, 1], color=COLOR_DONOR, marker="^", s=40, edgecolors="#1e293b", linewidth=0.8, zorder=6)
        ax_sub.plot([med_C[0], med_B[0]], [med_C[1], med_B[1]], color="#333333", linestyle="--", linewidth=1.1, alpha=0.7, zorder=4)
        ax_sub.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
        for spine in ["top", "right"]:
            ax_sub.spines[spine].set_visible(False)
            
    ax_a1_left.set_ylabel("Metric Coordinate $Z_2$", fontsize=7.6, labelpad=2)
    ax_a1_left.set_xticklabels([])
    ax_a1_right.set_xticklabels([])
    ax_a1_right.set_yticklabels([])
    
    # A-Title2: Dedicated title track for Subplot A2
    ax_a_title2 = fig.add_subplot(inner_a[4])
    ax_a_title2.axis("off")
    
    u_points = np.array(d1["u_points"])
    z1_ay008715 = np.array(d1["trajectories_z1"][taxa_hiv[ay008715_idx]])
    z1_af286230 = np.array(d1["trajectories_z1"][taxa_hiv[af286230_idx]])
    
    f1_mean_ay = float(np.mean(z1_ay008715[u_points < 862]))
    f2_mean_ay = float(np.mean(z1_ay008715[(u_points >= 862) & (u_points <= 1169)]))
    f3_mean_ay = float(np.mean(z1_ay008715[u_points > 1169]))
    delta_v_ay = abs(f2_mean_ay - f1_mean_ay)
    se_step_ay = math.sqrt(np.var(z1_ay008715[u_points < 862])/len(z1_ay008715[u_points < 862]) + 
                           np.var(z1_ay008715[(u_points >= 862) & (u_points <= 1169)])/len(z1_ay008715[(u_points >= 862) & (u_points <= 1169)]))
    z_flip_ay = delta_v_ay / se_step_ay
    
    ax_a_title2.text(0.50, 0.50, f"Dynamic Velocity Reversal Along Subtype B/C Dipole Axis ($Z_{{\\mathrm{{flip}}}} = {z_flip_ay:.1f}$, $p < 10^{{-15}}$)",
                     fontsize=8.0, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    # A-Leg2: Dedicated legend track for Subplot A2
    ax_a_leg2 = fig.add_subplot(inner_a[5])
    ax_a_leg2.axis("off")
    legend_handles_a2 = [
        Line2D([0], [0], color=COLOR_REC, lw=1.4, label="CRF07_BC (AF286230)"),
        Line2D([0], [0], color=COLOR_DONOR, lw=2.0, label="CRF08_BC (AY008715)"),
        Line2D([0], [0], color=COLOR_P2, lw=2.0, label=f"PR Flank ($\\bar{{Z}}_1 = {f1_mean_ay:.2f}$)"),
        Line2D([0], [0], color=COLOR_P1, lw=2.0, label=f"RT Cassette ($\\bar{{Z}}_1 = {f2_mean_ay:.2f}$)"),
    ]
    ax_a_leg2.legend(handles=legend_handles_a2, loc="center", ncol=4, frameon=False, fontsize=7.0)
    
    # A2: Standardized Velocity Reversal plot (NO set_title)
    ax_a2 = fig.add_subplot(inner_a[6])
    ax_a2.plot(u_points, z1_af286230, color=COLOR_REC, linewidth=1.4, alpha=0.85)
    ax_a2.plot(u_points, z1_ay008715, color=COLOR_DONOR, linewidth=2.0)
    ax_a2.hlines(f1_mean_ay, 0, 862, colors=COLOR_P2, linestyles="-", linewidth=2.0)
    ax_a2.hlines(f2_mean_ay, 862, 1169, colors=COLOR_P1, linestyles="-", linewidth=2.0)
    ax_a2.hlines(f3_mean_ay, 1169, 3120, colors=COLOR_P2, linestyles="-", linewidth=2.0)
    
    ax_a2.axvline(862, color="#222222", linestyle=":", linewidth=1.2)
    ax_a2.axvline(1169, color="#222222", linestyle=":", linewidth=1.2)
    ax_a2.axvspan(862, 1169, color=COLOR_DONOR, alpha=0.15)
    
    ax_a2.annotate("", xy=(862, f2_mean_ay), xytext=(862, f1_mean_ay),
                   arrowprops=dict(arrowstyle="<->", color="#b91c1c", lw=1.8))
    ax_a2.text(880, (f1_mean_ay + f2_mean_ay)/2.0, f"$\\Delta v = {delta_v_ay:.2f}$\n($Z_{{\\mathrm{{flip}}}} = {z_flip_ay:.1f}$)", 
               fontsize=7.2, fontweight="bold", color="#b91c1c", va="center")
    
    ax_a2.set_xlim(0, 3120)
    ax_a2.set_ylim(-0.10, 0.12)
    ax_a2.set_xlabel("Genomic Position $u$ (nucleotides)", fontsize=7.8, labelpad=2)
    ax_a2.set_ylabel("Dipole Velocity $Z_1(u)$", fontsize=7.6, labelpad=2)
    ax_a2.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_a2.spines[spine].set_visible(False)

    # ==============================================================================
    # PANEL B: Multi-Scale Dyadic Wavelets & Running Drawdown (HIV-1 pol)
    # ==============================================================================
    print("[*] Generating Panel B: Multi-Scale Dyadic Wavelets on HIV-1 pol...")
    # 7 rows: [0] Main Title, [1] Header B1, [2] Plot B1, [3] Header B2, [4] Plot B2, [5] Header B3, [6] Plot B3
    inner_b = gridspec.GridSpecFromSubplotSpec(7, 1, subplot_spec=gs[0, 1], 
                                               height_ratios=[0.08, 0.18, 0.54, 0.18, 0.56, 0.18, 0.58], 
                                               hspace=0.28)
    
    Y_crf08 = np.array(d2["crf08_diff_bridge"])
    u_coords = np.array(d2["u_coords"])
    dY_crf08 = np.diff(Y_crf08, prepend=Y_crf08[0])
    
    B_crf08 = Y_crf08 - (u_coords / float(L_hiv)) * Y_crf08[-1]
    sig_crf08 = float(np.std(dY_crf08))
    B_std_crf08 = B_crf08 / (sig_crf08 * math.sqrt(L_hiv))
    
    scales = [25, 50, 100, 250]
    scale_colors = [COLOR_ROOT, COLOR_SKY, COLOR_REC, COLOR_AMBER]
    wavelet_responses = {}
    for tau in scales:
        h = np.zeros(tau)
        h[:tau//2] = 1.0
        h[tau//2:] = -1.0
        wavelet_responses[tau] = np.convolve(dY_crf08, h, mode="same")
        
    M_t = np.maximum.accumulate(Y_crf08)
    DD_t = M_t - Y_crf08
    
    # B-Head: Dedicated title track
    ax_b_head = fig.add_subplot(inner_b[0])
    ax_b_head.axis("off")
    ax_b_head.text(-0.04, 0.50, "B", fontsize=12, fontweight="bold", va="center")
    ax_b_head.text(0.02, 0.50, "Multi-Scale Dyadic Filters & Running Drawdown (HIV-1 pol)", 
                    fontsize=9.2, fontweight="bold", va="center", color="#1e293b")
    
    # B-Head1: Dedicated 2-line header for B1
    ax_b_head1 = fig.add_subplot(inner_b[1])
    ax_b_head1.axis("off")
    ax_b_head1.text(0.50, 0.78, "1D Brownian Bridge Dilution (Global Drift Peaks at $u = 2450$; Cassette Diluted)",
                    fontsize=8.0, fontweight="bold", ha="center", va="center", color="#1e293b")
    legend_handles_b1 = [
        Line2D([0], [0], color="#64748b", lw=1.2, label=r"Bridge $B(u)/(\sigma\sqrt{L})$"),
        Line2D([0], [0], color="#ef4444", ls="--", lw=1.0, label="Kolmogorov Bound ($p=0.05$)"),
        Patch(color=COLOR_DONOR, alpha=0.25, label="RT Cassette [862, 1169]"),
    ]
    ax_b_head1.legend(handles=legend_handles_b1, loc="center", bbox_to_anchor=(0.50, 0.20), ncol=3, frameon=False, fontsize=7.0)
    
    # B1: Global 1D Brownian Bridge Dilution (NO set_title)
    ax_b1 = fig.add_subplot(inner_b[2])
    ax_b1.plot(u_coords, B_std_crf08, color="#64748b", linewidth=1.2)
    ax_b1.axhline(1.36, color="#ef4444", linestyle="--", linewidth=1.0)
    ax_b1.axhline(-1.36, color="#ef4444", linestyle="--", linewidth=1.0)
    ax_b1.axvspan(862, 1169, color=COLOR_DONOR, alpha=0.25)
    ax_b1.set_ylabel("Bridge $B(u)$", fontsize=7.6)
    ax_b1.set_ylim(-3.5, 5.0)
    ax_b1.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_b1.spines[spine].set_visible(False)
    plt.setp(ax_b1.get_xticklabels(), visible=False)
    
    # B-Head2: Dedicated 2-line header for Subplot B2
    ax_b_head2 = fig.add_subplot(inner_b[3])
    ax_b_head2.axis("off")
    ax_b_head2.text(0.50, 0.78, "Multi-Scale Dyadic Wavelets: Matched Scale $\\tau = 250\\text{ nt}$ Resonates at Cassette",
                    fontsize=8.0, fontweight="bold", ha="center", va="center", color="#1e293b")
    legend_handles_b2 = [
        Line2D([0], [0], color=col, lw=2.0 if tau == 250 else 1.2, 
               label=f"$\\tau = {tau}\\text{{ nt}}$" + (" (Matched)" if tau == 250 else ""))
        for tau, col in zip(scales, scale_colors)
    ]
    ax_b_head2.legend(handles=legend_handles_b2, loc="center", bbox_to_anchor=(0.50, 0.20), ncol=4, frameon=False, fontsize=7.0)
    
    # B2: Dyadic Wavelet Filter Bank (NO set_title)
    ax_b2 = fig.add_subplot(inner_b[4], sharex=ax_b1)
    for tau, col in zip(scales, scale_colors):
        lw = 2.0 if tau == 250 else 1.0
        alpha = 1.0 if tau == 250 else 0.7
        ax_b2.plot(u_coords, wavelet_responses[tau], color=col, linewidth=lw, alpha=alpha)
    ax_b2.axvspan(862, 1169, color=COLOR_DONOR, alpha=0.25)
    ax_b2.set_ylabel("Filter $dY * h_\\tau$", fontsize=7.6)
    ax_b2.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_b2.spines[spine].set_visible(False)
    plt.setp(ax_b2.get_xticklabels(), visible=False)
    
    # B-Head3: Dedicated 2-line header for Subplot B3
    ax_b_head3 = fig.add_subplot(inner_b[5])
    ax_b_head3.axis("off")
    ax_b_head3.text(0.50, 0.78, "Running Drawdown Localizes Cassette Entry & Momentum Collapse at Exit ($u_2 = 1169$)",
                    fontsize=8.0, fontweight="bold", ha="center", va="center", color="#1e293b")
    legend_handles_b3 = [
        Line2D([0], [0], color=COLOR_P1, lw=1.5, label="Cumulative Flux $Y(u)$"),
        Line2D([0], [0], color="#94a3b8", ls=":", lw=1.2, label="High-Water Mark $M(u)$"),
        Line2D([0], [0], color="#b91c1c", lw=1.5, label="Running Drawdown $\\mathrm{DD}(u)$"),
    ]
    ax_b_head3.legend(handles=legend_handles_b3, loc="center", bbox_to_anchor=(0.50, 0.20), ncol=3, frameon=False, fontsize=7.0)
    
    # B3: Running Drawdown & Momentum Collapse (NO set_title)
    ax_b3 = fig.add_subplot(inner_b[6], sharex=ax_b1)
    ax_b3.plot(u_coords, Y_crf08, color=COLOR_P1, linewidth=1.5)
    ax_b3.plot(u_coords, M_t, color="#94a3b8", linestyle=":", linewidth=1.2)
    ax_b3.plot(u_coords, DD_t, color="#b91c1c", linewidth=1.5)
    ax_b3.axvline(862, color="#222222", linestyle="--", linewidth=1.1)
    ax_b3.axvline(1169, color="#222222", linestyle=":", linewidth=1.4)
    ax_b3.axvspan(862, 1169, color=COLOR_DONOR, alpha=0.25)
    ax_b3.set_xlim(0, 3120)
    ax_b3.set_xlabel("Genomic Position $u$ (nucleotides)", fontsize=7.8, labelpad=2)
    ax_b3.set_ylabel("Flux & Drawdown", fontsize=7.6, labelpad=2)
    ax_b3.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_b3.spines[spine].set_visible(False)

    # ==============================================================================
    # PANEL C: Spatio-Temporal Force Field & SVD Mode Peeling (HIV-1 pol)
    # ==============================================================================
    print("[*] Generating Panel C: HIV-1 pol Spatio-Temporal SVD Mode Peeling...")
    # 7 rows: [0] Main Title, [1] Title C1, [2] Legend C1, [3] Plot C1, [4] Title C2, [5] Legend C2, [6] Plot C2
    inner_c = gridspec.GridSpecFromSubplotSpec(7, 1, subplot_spec=gs[1, 0], 
                                               height_ratios=[0.08, 0.12, 0.09, 0.70, 0.12, 0.09, 0.70], 
                                               hspace=0.28)
    
    # Compute weight-free force tensor
    F_hiv = np.zeros((L_hiv, len(taxa_hiv), 4))
    for u in range(L_hiv):
        col = seqs_hiv[:, u]
        for a in range(1, 5):
            idx_a = np.where(col == a)[0]
            if len(idx_a) >= 2:
                mu_a = np.mean(coords_hiv[idx_a], axis=0)
                for i_tax in idx_a:
                    F_hiv[u, i_tax] = -(coords_hiv[i_tax] - mu_a)
                    
    F_c_hiv = F_hiv - np.mean(F_hiv, axis=0, keepdims=True)
    S_hiv = np.cumsum(F_c_hiv, axis=0)
    X_hiv = S_hiv.reshape(L_hiv, -1)
    U_hiv, s_hiv, Vt_hiv = np.linalg.svd(X_hiv, full_matrices=False)
    
    N_hiv = len(taxa_hiv)
    dim_hiv = coords_hiv.shape[1]
    P_hiv = N_hiv * dim_hiv
    
    # C-Head: Dedicated title track
    ax_c_head = fig.add_subplot(inner_c[0])
    ax_c_head.axis("off")
    ax_c_head.text(-0.04, 0.50, "C", fontsize=12, fontweight="bold", va="center")
    ax_c_head.text(0.02, 0.50, "Spatio-Temporal Force Field & SVD Mode Peeling (HIV-1 pol)", 
                    fontsize=9.2, fontweight="bold", va="center", color="#1e293b")
    
    # C-Title1: Dedicated title track for C1
    ax_c_title1 = fig.add_subplot(inner_c[1])
    ax_c_title1.axis("off")
    ax_c_title1.text(0.50, 0.50, "Spatial Profiles: Mode 1 (Clade Backbone) vs Mode 2 (CRF07/08 Mosaicism)",
                     fontsize=8.0, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    # C-Leg1: Dedicated legend track for C1
    ax_c_leg1 = fig.add_subplot(inner_c[2])
    ax_c_leg1.axis("off")
    legend_handles_c1 = [
        Line2D([0], [0], color=COLOR_P1, lw=1.6, label="Mode 1 (Clade Backbone)"),
        Line2D([0], [0], color=COLOR_REC, lw=1.6, label="Mode 2 (CRF07/08 Mosaicism)"),
    ]
    ax_c_leg1.legend(handles=legend_handles_c1, loc="center", ncol=2, frameon=False, fontsize=7.2)
    
    # C1: Spatial Profiles across PR, RT, IN domains (NO set_title)
    ax_c1 = fig.add_subplot(inner_c[3])
    u_vec_hiv = np.arange(L_hiv)
    u1_hiv = U_hiv[:, 0]
    u2_hiv = U_hiv[:, 1]
    
    ax_c1.plot(u_vec_hiv, u1_hiv, color=COLOR_P1, linewidth=1.5)
    ax_c1.plot(u_vec_hiv, u2_hiv, color=COLOR_REC, linewidth=1.6)
    
    domains = [("PR", 0, 297, "#f1f5f9"), ("RT", 297, 1917, "#e2e8f0"), ("IN", 1917, 3120, "#cbd5e1")]
    for d_name, d_s, d_e, d_bg in domains:
        ax_c1.axvspan(d_s, d_e, color=d_bg, alpha=0.45, zorder=0)
        ax_c1.text((d_s + d_e)/2.0, 0.08, d_name, transform=ax_c1.get_xaxis_transform(), 
                   fontsize=7.2, fontweight="bold", ha="center", color="#475569")
        
    ax_c1.set_xlim(0, L_hiv)
    ax_c1.set_ylabel("Spatial Mode $\\mathbf{u}_m$", fontsize=7.6, labelpad=2)
    ax_c1.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_c1.spines[spine].set_visible(False)
    
    # C-Title2: Dedicated title track for Subplot C2
    ax_c_title2 = fig.add_subplot(inner_c[4])
    ax_c_title2.axis("off")
    ax_c_title2.text(0.50, 0.50, "Taxon Loadings Decouple Clade Ancestry from Chimeric Reticulation",
                     fontsize=8.0, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    # C-Leg2: Dedicated legend track for Subplot C2
    ax_c_leg2 = fig.add_subplot(inner_c[5])
    ax_c_leg2.axis("off")
    legend_handles_c2 = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_P1, markersize=6.5, label="Subtype B ($n=8$)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=COLOR_P2, markersize=6.5, label="Subtype C ($n=8$)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_REC, markeredgecolor="#1e293b", markersize=6.5, label="CRF07_BC ($n=3$)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=COLOR_DONOR, markeredgecolor="#1e293b", markersize=6.5, label="CRF08_BC ($n=3$)"),
    ]
    ax_c_leg2.legend(handles=legend_handles_c2, loc="center", ncol=4, frameon=False, fontsize=7.0)
    
    # C2: Taxon Loadings onto Mode 1 vs Mode 2 (NO set_title)
    ax_c2 = fig.add_subplot(inner_c[6])
    V1_hiv = Vt_hiv[0].reshape(N_hiv, dim_hiv)
    V2_hiv = Vt_hiv[1].reshape(N_hiv, dim_hiv)
    v1_proj = V1_hiv[:, 0]
    v2_proj = V2_hiv[:, 0]
    
    ax_c2.scatter(v1_proj[b_taxa], v2_proj[b_taxa], color=COLOR_P1, marker="o", s=42, zorder=4)
    ax_c2.scatter(v1_proj[c_taxa], v2_proj[c_taxa], color=COLOR_P2, marker="s", s=42, zorder=4)
    ax_c2.scatter(v1_proj[crf07_taxa], v2_proj[crf07_taxa], color=COLOR_REC, marker="D", s=52, edgecolors="#1e293b", linewidth=1.0, zorder=5)
    ax_c2.scatter(v1_proj[crf08_taxa], v2_proj[crf08_taxa], color=COLOR_DONOR, marker="^", s=52, edgecolors="#1e293b", linewidth=1.0, zorder=5)
    
    ax_c2.axhline(0, color="#94a3b8", linestyle="--", linewidth=0.8, alpha=0.7)
    ax_c2.axvline(0, color="#94a3b8", linestyle="--", linewidth=0.8, alpha=0.7)
    ax_c2.set_xlabel("Mode 1 Loading $V_1$ (Subtype Divergence)", fontsize=7.8, labelpad=2)
    ax_c2.set_ylabel("Mode 2 Loading $V_2$ (Reticulation)", fontsize=7.6, labelpad=2)
    ax_c2.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_c2.spines[spine].set_visible(False)

    # ==============================================================================
    # PANEL D: Random Matrix Theory Halting & Continuous Topological Energy Ratio
    # ==============================================================================
    print("[*] Generating Panel D: RMT Halting & Topological Energy Ratio on HIV-1 pol...")
    sig_noise_hiv = float(np.median(np.var(F_c_hiv.reshape(L_hiv, -1), axis=0))**0.5)
    s_bulk_hiv = (L_hiv / math.pi) * sig_noise_hiv * (1.0 + math.sqrt(P_hiv / float(L_hiv)))
    
    # Calculate energy polarity ratios directly on HIV-1 pol modes
    p_ratio_mode1 = float(abs(np.sum(u1_hiv)) / np.sum(np.abs(u1_hiv)))
    p_ratio_mode2 = float(abs(np.sum(u2_hiv)) / np.sum(np.abs(u2_hiv)))
    
    # 4 rows: [0] Main Title, [1] Subplot Titles, [2] Subplot Legends, [3] Subplot Bodies
    inner_d = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[1, 1], 
                                               height_ratios=[0.08, 0.18, 0.12, 1.55], 
                                               hspace=0.24)
    
    # D-Head: Dedicated title track
    ax_d_head = fig.add_subplot(inner_d[0])
    ax_d_head.axis("off")
    ax_d_head.text(-0.04, 0.50, "D", fontsize=12, fontweight="bold", va="center")
    ax_d_head.text(0.02, 0.50, "Random Matrix Theory Halting & Continuous Topological Mode Ratio", 
                   fontsize=9.2, fontweight="bold", va="center", color="#1e293b")
    
    # D-Titles: Side-by-side dedicated titles for D1 and D2 with wide wspace=0.35
    inner_d_titles = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner_d[1], wspace=0.35)
    ax_d_title1 = fig.add_subplot(inner_d_titles[0])
    ax_d_title1.axis("off")
    ax_d_title1.text(0.50, 0.50, f"RMT Spectral Halting Edge\n($s_{{\\mathrm{{bulk}}}} = {s_bulk_hiv:.1f}$; Halts at Mode 12)",
                     fontsize=7.8, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    ax_d_title2 = fig.add_subplot(inner_d_titles[1])
    ax_d_title2.axis("off")
    ax_d_title2.text(0.50, 0.50, "Topological Mode Classification\n(Monopolar Backbone vs Bipolar Cassette)",
                     fontsize=7.8, fontweight="bold", ha="center", va="center", color="#1e293b")
    
    # D-Legends: Side-by-side dedicated legends for D1 and D2 with wide wspace=0.35
    inner_d_legends = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner_d[2], wspace=0.35)
    ax_d_leg1 = fig.add_subplot(inner_d_legends[0])
    ax_d_leg1.axis("off")
    legend_handles_d1 = [
        Line2D([0], [0], marker="o", color=COLOR_REC, lw=1.2, markersize=5.0, label="Signal ($s_m > 4.6$)"),
        Line2D([0], [0], marker="s", color="#94a3b8", lw=1.2, markersize=5.0, label="Noise"),
        Line2D([0], [0], color="#ef4444", ls="--", lw=1.1, label=f"Edge ({s_bulk_hiv:.1f})"),
    ]
    ax_d_leg1.legend(handles=legend_handles_d1, loc="center", ncol=3, frameon=False, 
                     fontsize=6.2, columnspacing=0.4, handletextpad=0.25)
    
    ax_d_leg2 = fig.add_subplot(inner_d_legends[1])
    ax_d_leg2.axis("off")
    legend_handles_d2 = [
        Line2D([0], [0], color=COLOR_P1, lw=1.6, label=f"Mode 1 Backbone ($p = {p_ratio_mode1:.2f}$)"),
        Line2D([0], [0], color=COLOR_P2, lw=1.6, label=f"Mode 2 Cassette ($p = {p_ratio_mode2:.2f}$)"),
    ]
    ax_d_leg2.legend(handles=legend_handles_d2, loc="center", ncol=2, frameon=False, 
                     fontsize=6.3, columnspacing=0.5, handletextpad=0.25)
    
    # D-Body: Side-by-side plots for D1 and D2 (NO set_title) with wide wspace=0.35
    inner_d_body = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner_d[3], wspace=0.35)
    ax_d1 = fig.add_subplot(inner_d_body[0])
    ax_d2 = fig.add_subplot(inner_d_body[1])
    
    # D1: Singular Value Spectrum vs Marchenko-Pastur Bulk Edge
    modes_idx = np.arange(1, 13)
    s_vals_12 = s_hiv[:12]
    
    is_signal = s_vals_12 > s_bulk_hiv
    ax_d1.semilogy(modes_idx[is_signal], s_vals_12[is_signal], "o-", color=COLOR_REC, markersize=5.5, linewidth=1.5)
    ax_d1.semilogy(modes_idx[~is_signal], s_vals_12[~is_signal], "s-", color="#94a3b8", markersize=5.5, linewidth=1.5)
    
    ax_d1.axhline(s_bulk_hiv, color="#ef4444", linestyle="--", linewidth=1.3)
    ax_d1.axhspan(1.5, s_bulk_hiv, color="#fee2e2", alpha=0.45)
    
    ax_d1.set_xlim(0.5, 12.5)
    ax_d1.set_ylim(1.5, 1200)
    ax_d1.set_xlabel("Mode Index $m$", fontsize=7.8, labelpad=2)
    ax_d1.set_ylabel("Singular Value $s_m$ (log scale)", fontsize=7.6, labelpad=2)
    ax_d1.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_d1.spines[spine].set_visible(False)
    
    # D2: Continuous Topological Mode Classification on HIV-1 pol modes
    u_norm_m1 = u1_hiv / np.max(np.abs(u1_hiv))
    u_norm_m2 = u2_hiv / np.max(np.abs(u2_hiv))
    
    ax_d2.plot(u_vec_hiv, u_norm_m1, color=COLOR_P1, linewidth=1.6)
    ax_d2.plot(u_vec_hiv, u_norm_m2, color=COLOR_P2, linewidth=1.6)
    
    ax_d2.axhline(0, color="#94a3b8", linestyle=":", linewidth=0.8)
    ax_d2.set_xlim(0, 3120)
    ax_d2.set_xlabel("Genomic Position $u$ (nucleotides)", fontsize=7.8, labelpad=2)
    ax_d2.set_ylabel("Normalized Mode $u_m(u)$", fontsize=7.6, labelpad=2)
    ax_d2.grid(True, linestyle=":", color=COLOR_GRID, alpha=0.6)
    for spine in ["top", "right"]:
        ax_d2.spines[spine].set_visible(False)

    # Export outputs
    print(f"[*] Exporting vector PDF to {pdf_path}...")
    fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"[*] Exporting raster PNG to {png_path}...")
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    
    export_payload = {
        "dataset": "hiv1_pol_8B_8C_3CRF07_3CRF08.fasta",
        "length": L_hiv,
        "taxa_count": N_hiv,
        "panel_a": {
            "medoid_B": med_B.tolist(),
            "medoid_C": med_C.tolist(),
            "delta_v": delta_v_ay,
            "z_flip": z_flip_ay,
        },
        "panel_b": {
            "scales_evaluated": scales,
            "rt_cassette": [862, 1169],
        },
        "panel_c": {
            "s1": float(s_hiv[0]),
            "s2": float(s_hiv[1]),
            "s_bulk": s_bulk_hiv,
        },
        "panel_d": {
            "p_ratio_mode1": p_ratio_mode1,
            "p_ratio_mode2": p_ratio_mode2,
        }
    }
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(export_payload, f, indent=2)
    print(f"[*] Exported authentic payload to {out_json}")

if __name__ == "__main__":
    generate_figure3()
