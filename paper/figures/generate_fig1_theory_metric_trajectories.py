"""
paper/figures/generate_fig1_theory_metric_trajectories.py
=========================================================
Generates Figure 1 for Section 2 of the RhizAeon manuscript:
"Continuous Sequence Metric Manifolds Resolve Competing Retroviral Chimeras in Linear Time"

CRITICAL SCIENTIFIC INTEGRITY MANDATE:
- All data, coordinates, trajectories, corridors, breakpoints, and tangent flow vectors
  are loaded directly from the authentic outputs of the native Rust RhizAeon engine
  (`crates/rhizaeon-cli/src/bin/export_figure1_data.rs` / `fig1_rhizaeon_authentic_data.json`).
- ZERO ad-hoc sliding-window re-implementations.
- ZERO synthetic placeholders or parametric mockups.

BioVis-Expert Guidelines Implemented:
- Strict Okabe-Ito colorblind-accessible palette:
    Subtype B: Navy Blue (#0072B2)
    Subtype C: Vermilion (#D55E00)
    CRF07_BC: Emerald Green (#009E73)
    CRF08_BC: Royal Reddish-Purple (#782869 / #CC79A7)
- 4-Panel Progressive Macro-Layout (2x2 balanced grid):
    (A) Canonical Sequence Metric Space (Buneman Coordinates / Classical MDS)
    (B) Clonal Stationarity and Coalescent Wobble (Non-Recombinant Reference Corridors)
    (C) Dynamic Recombinant Trajectories and Resolved Breakpoints Across Gene Domains
    (D) Metric Phase Space Dynamics (Continuous 2D Orbital Loops)
- Despined L-frames, explicit % variance explained, aligned gene domain track, high data-to-ink ratio.
- High-res vector PDF and 300 DPI raster PNG export.
"""

from pathlib import Path
import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter1d

# Paths
script_dir = Path(__file__).resolve().parent
repo_root = script_dir.parent.parent
json_path = script_dir / "fig1_rhizaeon_authentic_data.json"
pdf_path = script_dir / "fig1_theory_metric_trajectories.pdf"
png_path = script_dir / "fig1_theory_metric_trajectories.png"

# Ensure data exists; if not, execute native Rust exporter
if not json_path.exists():
    print(f"[*] Authentic data JSON missing at {json_path}. Compiling and running Rust exporter...")
    cmd = ["cargo", "run", "--bin", "export_figure1_data"]
    subprocess.run(cmd, cwd=str(repo_root), check=True)

print(f"[*] Loading authentic RhizAeon engine data from {json_path}...")
with open(json_path, "r", encoding="utf-8") as f:
    data = json.load(f)

# Extract core metadata and arrays
L = data["alignment_length"]
N = data["num_taxa"]
taxa = data["taxa"]
b_indices = data["b_indices"]
c_indices = data["c_indices"]
crf07_indices = data["crf07_indices"]
crf08_indices = data["crf08_indices"]
af286230_idx = data["af286230_idx"]
ay008715_idx = data["ay008715_idx"]
af_name = taxa[af286230_idx]
ay_name = taxa[ay008715_idx]

var_exp = data["variance_explained"]
u_points = np.array(data["u_points"])
num_pts = len(u_points)

# Global MDS coordinates
coords = np.zeros((N, 2), dtype=np.float64)
for i, t in enumerate(taxa):
    coords[i] = data["global_coords"][t]

# Trajectories
trajs_z1 = np.zeros((N, num_pts), dtype=np.float64)
trajs_z2 = np.zeros((N, num_pts), dtype=np.float64)
for i, t in enumerate(taxa):
    trajs_z1[i] = data["trajectories_z1"][t]
    trajs_z2[i] = data["trajectories_z2"][t]

# Corridors
corr = data["reference_corridors"]
mean_b = np.array(corr["subtype_b_mean"])
std_b = np.array(corr["subtype_b_std"])
min_b = np.array(corr["subtype_b_min"])
max_b = np.array(corr["subtype_b_max"])
b_overall_mean = corr["subtype_b_overall_mean"]
b_flutter = corr["subtype_b_overall_flutter"]

mean_c = np.array(corr["subtype_c_mean"])
std_c = np.array(corr["subtype_c_std"])
min_c = np.array(corr["subtype_c_min"])
max_c = np.array(corr["subtype_c_max"])
c_overall_mean = corr["subtype_c_overall_mean"]
c_flutter = corr["subtype_c_overall_flutter"]

# Breakpoints
bps = data["resolved_breakpoints"]
crf08_entry = bps["crf08_entry_nt"]
crf07_entry = bps["crf07_entry_nt"]
stagger_offset = bps["stagger_offset_nt"]
shared_exit = bps["shared_exit_nt"]

# Orbits
orbits = data["phase_space_orbits"]
af_orbit = orbits[af_name]
ay_orbit = orbits[ay_name]

print(f"[*] Authentic data loaded: N={N}, L={L}, VarExp=[{var_exp[0]:.1f}%, {var_exp[1]:.1f}%]")
print(f"[*] Corridors: B mean={b_overall_mean:+.3f} (sigma={b_flutter:.3f}), C mean={c_overall_mean:+.3f} (sigma={c_flutter:.3f})")
print(f"[*] RT Breakpoints: CRF08={crf08_entry} nt, CRF07={crf07_entry} nt, Exit={shared_exit} nt")

# -------------------------------------------------------------
# Figure Architecture & Styling (BioVis-Expert Standards)
# -------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.edgecolor": "#334155",
    "axes.linewidth": 0.8,
    "grid.color": "#f1f5f9",
    "grid.linestyle": "-",
    "grid.linewidth": 0.5,
})

# Okabe-Ito Colorblind-Safe Palette
blue_main = "#0072B2"       # Subtype B
blue_light = "#56B4E9"
vermilion_main = "#D55E00"  # Subtype C
vermilion_light = "#E69F00"
crf07_color = "#009E73"     # CRF07_BC (Emerald Bluish-Green)
crf07_dark = "#006644"
crf08_color = "#CC79A7"     # CRF08_BC (Reddish Purple)
crf08_dark = "#782869"

fig = plt.figure(figsize=(15.2, 11.8), dpi=300)

outer_grid = gridspec.GridSpec(
    2, 2,
    height_ratios=[1.0, 1.05],
    width_ratios=[1.0, 1.0],
    hspace=0.28,
    wspace=0.22,
    left=0.065,
    right=0.975,
    bottom=0.06,
    top=0.95
)

# =============================================================
# PANEL A: Global Sequence Metric Space (Buneman Coordinates)
# =============================================================
ax_a = fig.add_subplot(outer_grid[0, 0])
ax_a.grid(True, zorder=0)
ax_a.axvline(0, color="#94a3b8", linestyle="--", linewidth=0.8, zorder=1)
ax_a.axhline(0, color="#94a3b8", linestyle=":", linewidth=0.6, zorder=1)

# Plot Subtype B references
ax_a.scatter(
    coords[b_indices, 0], coords[b_indices, 1],
    color=blue_main, s=55, marker="o", edgecolors="#0f172a", lw=0.8,
    label=f"Subtype B References ($n=8$)", zorder=3
)
# Plot Subtype C references
ax_a.scatter(
    coords[c_indices, 0], coords[c_indices, 1],
    color=vermilion_main, s=55, marker="o", edgecolors="#0f172a", lw=0.8,
    label=f"Subtype C References ($n=8$)", zorder=3
)
# Plot CRF07_BC isolates
ax_a.scatter(
    coords[crf07_indices, 0], coords[crf07_indices, 1],
    color=crf07_color, s=75, marker="D", edgecolors="#0f172a", lw=0.9,
    label=f"CRF07_BC Chimeras ($n=3$)", zorder=4
)
# Plot CRF08_BC isolates
ax_a.scatter(
    coords[crf08_indices, 0], coords[crf08_indices, 1],
    color=crf08_dark, s=70, marker="s", edgecolors="#0f172a", lw=0.9,
    label=f"CRF08_BC Chimeras ($n=3$)", zorder=4
)

# Text Callouts & Cluster Labels
ax_a.text(
    0.050, -0.016, "Subtype B Cluster\n($Z_1 > 0$, Pure Clade)",
    ha="center", fontsize=7.5, fontweight="bold", color=blue_main,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#eff6ff", edgecolor=blue_main, lw=0.6, alpha=0.9)
)
ax_a.text(
    -0.050, -0.016, "Subtype C Cluster\n($Z_1 < 0$, Pure Clade)",
    ha="center", fontsize=7.5, fontweight="bold", color=vermilion_main,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#fff7ed", edgecolor=vermilion_main, lw=0.6, alpha=0.9)
)
ax_a.text(
    -0.034, 0.018, "Recombinants Compressed at\nIntermediate Coordinates ($Z_1 \\approx -0.03$ to $-0.04$)",
    ha="center", fontsize=7.2, fontweight="semibold", color="#1e293b",
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#f8fafc", edgecolor="#94a3b8", lw=0.6, alpha=0.9)
)

ax_a.set_title("A   Global Buneman Metric Space ($N=22$ sequences)", loc="left", fontsize=10.5, fontweight="bold", color="#0f172a", pad=8)
ax_a.set_xlabel(f"Primary Canonical Axis $Z_1$ ({var_exp[0]:.1f}% variance)", fontsize=9.2, fontweight="semibold")
ax_a.set_ylabel(f"Secondary Canonical Axis $Z_2$ ({var_exp[1]:.1f}% variance)", fontsize=9.2, fontweight="semibold")
ax_a.tick_params(labelsize=8)
ax_a.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.2)
ax_a.set_ylim(-0.026, 0.028)

# =============================================================
# PANEL B: Clonal Stationarity and Coalescent Wobble
# =============================================================
ax_b = fig.add_subplot(outer_grid[0, 1])
ax_b.grid(True, zorder=0)
ax_b.axhline(0, color="#64748b", linestyle="--", linewidth=0.9, zorder=1)

# Shaded reference corridors (empirical min-max bounds)
ax_b.fill_between(u_points, min_b, max_b, color=blue_main, alpha=0.15, zorder=1, label="Subtype B Empirical Envelope")
ax_b.fill_between(u_points, min_c, max_c, color=vermilion_main, alpha=0.15, zorder=1, label="Subtype C Empirical Envelope")

# Individual non-recombinant reference trajectories from authentic Rust engine
for idx in b_indices:
    ax_b.plot(u_points, trajs_z1[idx], color=blue_main, alpha=0.25, lw=0.8, zorder=2)
for idx in c_indices:
    ax_b.plot(u_points, trajs_z1[idx], color=vermilion_main, alpha=0.25, lw=0.8, zorder=2)

# Clade mean trajectories
ax_b.plot(u_points, mean_b, color=blue_main, lw=2.2, label="Subtype B Mean ($n=8$)", zorder=3)
ax_b.plot(u_points, mean_c, color=vermilion_main, lw=2.2, label="Subtype C Mean ($n=8$)", zorder=3)

# Annotations emphasizing stationarity
ax_b.text(
    1560, 0.090, f"Subtype B Clonal Corridor ($Z_1 > 0$)\nStationary Mean = {b_overall_mean:+.3f}, Flutter $\\sigma = {b_flutter:.3f}$",
    ha="center", fontsize=7.5, fontweight="bold", color=blue_main,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#eff6ff", edgecolor=blue_main, lw=0.6, alpha=0.9)
)
ax_b.text(
    1560, -0.090, f"Subtype C Clonal Corridor ($Z_1 < 0$)\nStationary Mean = {c_overall_mean:+.3f}, Flutter $\\sigma = {c_flutter:.3f}$",
    ha="center", fontsize=7.5, fontweight="bold", color=vermilion_main,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#fff7ed", edgecolor=vermilion_main, lw=0.6, alpha=0.9)
)
ax_b.text(
    2750, 0.008, "Wobble Never Crosses Clades ($Z_1 = 0$)",
    ha="center", fontsize=7.2, fontweight="bold", color="#475569"
)

ax_b.set_title("B   Clonal Stationarity and Coalescent Wobble ($16$ Non-Recombinants)", loc="left", fontsize=10.5, fontweight="bold", color="#0f172a", pad=8)
ax_b.set_xlabel("Genomic Coordinate $u$ along HIV-1 pol (nt)", fontsize=9.2, fontweight="semibold")
ax_b.set_ylabel("Instantaneous Trajectory $Z_1(u)$", fontsize=9.2, fontweight="semibold")
ax_b.set_xlim(0, L)
ax_b.set_ylim(-0.135, 0.135)
ax_b.tick_params(labelsize=8)
ax_b.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.0, ncol=2)

# =============================================================
# PANEL C: Dynamic Recombinants + Seamless Gene Domain Track
# =============================================================
inner_c = gridspec.GridSpecFromSubplotSpec(
    2, 1,
    subplot_spec=outer_grid[1, 0],
    height_ratios=[1.0, 0.08],
    hspace=0.03
)

ax_c = fig.add_subplot(inner_c[0])
ax_c.grid(True, zorder=0)
ax_c.axhline(0, color="#64748b", linestyle="--", linewidth=0.9, zorder=1)

# Background reference corridors
ax_c.fill_between(u_points, min_b, max_b, color=blue_main, alpha=0.10, zorder=1)
ax_c.fill_between(u_points, min_c, max_c, color=vermilion_main, alpha=0.10, zorder=1)
ax_c.plot(u_points, mean_b, color=blue_main, lw=1.0, alpha=0.45, linestyle=":", zorder=2)
ax_c.plot(u_points, mean_c, color=vermilion_main, lw=1.0, alpha=0.45, linestyle=":", zorder=2)

# Shaded cassette spans based on authentic breakpoints
ax_c.axvspan(0, 465, color="#f1f5f9", alpha=0.6, zorder=0)  # PR Discordance
ax_c.axvspan(crf08_entry, shared_exit, color=crf08_color, alpha=0.12, zorder=0)
ax_c.axvspan(crf07_entry, shared_exit, color=crf07_color, alpha=0.12, zorder=0)

# Plot CRF07_BC isolates
for idx in crf07_indices:
    if idx == af286230_idx:
        continue
    ax_c.plot(u_points, trajs_z1[idx], color=crf07_color, alpha=0.45, lw=1.1, ls="--", zorder=3)
ax_c.plot(
    u_points, trajs_z1[af286230_idx],
    color=crf07_color, lw=2.4, ls="-",
    label="CRF07_BC (AF286230 Prototype, $n=3$)", zorder=5
)

# Plot CRF08_BC isolates
for idx in crf08_indices:
    if idx == ay008715_idx:
        continue
    ax_c.plot(u_points, trajs_z1[idx], color=crf08_dark, alpha=0.45, lw=1.1, ls=":", zorder=3)
ax_c.plot(
    u_points, trajs_z1[ay008715_idx],
    color=crf08_dark, lw=2.4, ls="-",
    label="CRF08_BC (AY008715 Prototype, $n=3$)", zorder=5
)

# Breakpoint Annotations directly from authentic breakpoints
ax_c.annotate(
    f"CRF08 Entry into B ({crf08_entry} nt)\n(Enters ~{stagger_offset} nt earlier)",
    xy=(crf08_entry, 0.0), xytext=(520, 0.082),
    arrowprops=dict(facecolor=crf08_dark, edgecolor=crf08_dark, shrink=0.08, width=0.8, headwidth=4),
    fontsize=7.2, fontweight="bold", color=crf08_dark
)
ax_c.annotate(
    f"CRF07 Entry into B ({crf07_entry} nt)\n(Shifted 3' downstream)",
    xy=(crf07_entry, 0.0), xytext=(1220, 0.082),
    arrowprops=dict(facecolor=crf07_dark, edgecolor=crf07_dark, shrink=0.08, width=0.8, headwidth=4),
    fontsize=7.2, fontweight="bold", color=crf07_dark
)
ax_c.annotate(
    "PR Discordance:\nCRF07 in B vs. CRF08 in C",
    xy=(280, 0.0), xytext=(90, -0.015),
    arrowprops=dict(facecolor="#0f172a", edgecolor="#0f172a", shrink=0.08, width=0.8, headwidth=4),
    fontsize=7.2, fontweight="bold", color="#0f172a"
)

ax_c.set_title("C   Dynamic Recombinant Trajectories ($Z_1(u)$ vs. Parental Corridors)", loc="left", fontsize=10.5, fontweight="bold", color="#0f172a", pad=8)
ax_c.set_ylabel("Instantaneous Trajectory $Z_1(u)$", fontsize=9.2, fontweight="semibold")
ax_c.set_xlim(0, L)
ax_c.set_ylim(-0.135, 0.135)
ax_c.tick_params(labelsize=8)
ax_c.set_xticklabels([])
ax_c.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=7.0)

# Seamless Gene Domain Track
ax_track = fig.add_subplot(inner_c[1])
ax_track.set_xlim(0, L)
ax_track.set_ylim(0, 1)
ax_track.grid(False)

domains = [
    ("p6*", 0, 168, "#94a3b8"),
    ("PR", 168, 465, "#64748b"),
    ("Reverse Transcriptase (Polymerase)", 465, 1785, "#3b82f6"),
    ("RNase H", 1785, 2145, "#0ea5e9"),
    ("Integrase (IN)", 2145, 3120, "#6366f1"),
]

for name, start, end, col in domains:
    rect = Rectangle((start, 0.05), end - start, 0.90, facecolor=col, edgecolor="#1e293b", lw=0.7, alpha=0.90)
    ax_track.add_patch(rect)
    mid = (start + end) / 2
    ax_track.text(mid, 0.50, name, ha="center", va="center", fontsize=7.0, fontweight="bold", color="white")

ax_track.spines["top"].set_visible(False)
ax_track.spines["left"].set_visible(False)
ax_track.spines["right"].set_visible(False)
ax_track.set_yticks([])
ax_track.set_xlabel("Genomic Coordinate $u$ along HIV-1 pol (nt)", fontsize=9.0, fontweight="semibold", labelpad=2)
ax_track.tick_params(axis="x", labelsize=8)

# =============================================================
# PANEL D: Metric Phase Space Dynamics (2D Orbital Loops)
# =============================================================
ax_d = fig.add_subplot(outer_grid[1, 1])
ax_d.grid(True, zorder=0)
ax_d.axvline(0, color="#64748b", linestyle="--", linewidth=0.8, alpha=0.6, zorder=1)
ax_d.axhline(0, color="#64748b", linestyle=":", linewidth=0.6, alpha=0.6, zorder=1)

# Reference attractor basins (stationary scatter)
ax_d.scatter(coords[b_indices, 0], coords[b_indices, 1], color=blue_main, s=40, alpha=0.75, edgecolors="#0f172a", lw=0.6, zorder=2)
ax_d.scatter(coords[c_indices, 0], coords[c_indices, 1], color=vermilion_main, s=40, alpha=0.75, edgecolors="#0f172a", lw=0.6, zorder=2)

# Raw authentic trajectory points from Rust engine
raw_07_x = np.array(af_orbit["z1"])
raw_07_y = np.array(af_orbit["z2"])
raw_08_x = np.array(ay_orbit["z1"])
raw_08_y = np.array(ay_orbit["z2"])

# Light gaussian filter to reveal macroscopic dynamical orbit
smooth_07_x = gaussian_filter1d(raw_07_x, sigma=1.8)
smooth_07_y = gaussian_filter1d(raw_07_y, sigma=1.8)
smooth_08_x = gaussian_filter1d(raw_08_x, sigma=1.8)
smooth_08_y = gaussian_filter1d(raw_08_y, sigma=1.8)

ax_d.plot(raw_07_x, raw_07_y, color=crf07_color, lw=0.7, alpha=0.35, zorder=3)
ax_d.plot(raw_08_x, raw_08_y, color=crf08_dark, lw=0.7, alpha=0.35, zorder=3)

ax_d.plot(smooth_07_x, smooth_07_y, color=crf07_color, lw=2.4, alpha=0.95, zorder=4, label="CRF07_BC Phase Orbit")
ax_d.plot(smooth_08_x, smooth_08_y, color=crf08_dark, lw=2.4, alpha=0.95, zorder=4, label="CRF08_BC Phase Orbit")

# Tangent flow vectors computed directly by Rust engine
for arrow in af_orbit["tangent_arrows"]:
    ax_d.annotate(
        "", xy=(arrow["x"] + 0.007 * arrow["dx"], arrow["y"] + 0.007 * arrow["dy"]),
        xytext=(arrow["x"], arrow["y"]),
        arrowprops=dict(arrowstyle="->", color=crf07_dark, lw=1.6, mutation_scale=10),
        zorder=5
    )

for arrow in ay_orbit["tangent_arrows"]:
    ax_d.annotate(
        "", xy=(arrow["x"] + 0.007 * arrow["dx"], arrow["y"] + 0.007 * arrow["dy"]),
        xytext=(arrow["x"], arrow["y"]),
        arrowprops=dict(arrowstyle="->", color=crf08_dark, lw=1.6, mutation_scale=10),
        zorder=5
    )

# Start (5') and End (3') markers from authentic coordinates
ax_d.scatter(af_orbit["start_5p"][0], af_orbit["start_5p"][1], marker="o", color=crf07_color, s=75, edgecolors="#0f172a", lw=1.2, zorder=6, label="CRF07 5' Start")
ax_d.scatter(ay_orbit["start_5p"][0], ay_orbit["start_5p"][1], marker="o", color=crf08_dark, s=75, edgecolors="#0f172a", lw=1.2, zorder=6, label="CRF08 5' Start")
ax_d.scatter(af_orbit["end_3p"][0], af_orbit["end_3p"][1], marker="s", color=crf07_color, s=65, edgecolors="#0f172a", lw=1.2, zorder=6, label="CRF07 3' End")
ax_d.scatter(ay_orbit["end_3p"][0], ay_orbit["end_3p"][1], marker="s", color=crf08_dark, s=65, edgecolors="#0f172a", lw=1.2, zorder=6, label="CRF08 3' End")

# Basin labels
ax_d.text(
    0.055, 0.046, "Subtype B Basin\n($Z_1 > 0$)", ha="center", fontsize=7.5, fontweight="bold", color=blue_main,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#eff6ff", edgecolor=blue_main, lw=0.6, alpha=0.9)
)
ax_d.text(
    -0.065, 0.046, "Subtype C Basin\n($Z_1 < 0$)", ha="center", fontsize=7.5, fontweight="bold", color=vermilion_main,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="#fff7ed", edgecolor=vermilion_main, lw=0.6, alpha=0.9)
)

ax_d.set_title("D   2D Metric Phase Space Dynamics $(Z_1 \\times Z_2)$", loc="left", fontsize=10.5, fontweight="bold", color="#0f172a", pad=8)
ax_d.set_xlabel(f"Primary Canonical Axis $Z_1$ ({var_exp[0]:.1f}% variance)", fontsize=9.2, fontweight="semibold")
ax_d.set_ylabel(f"Secondary Canonical Axis $Z_2$ ({var_exp[1]:.1f}% variance)", fontsize=9.2, fontweight="semibold")
ax_d.set_ylim(-0.060, 0.065)
ax_d.tick_params(labelsize=8)
ax_d.legend(loc="lower left", frameon=True, framealpha=0.92, fontsize=6.8, ncol=2)

# Export high-res PDF and PNG
print(f"[*] Saving vector PDF to {pdf_path}...")
fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
print(f"[*] Saving 300 DPI raster PNG to {png_path}...")
fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
plt.close(fig)

print("[+] Figure 1 successfully generated from authentic RhizAeon engine data.")
