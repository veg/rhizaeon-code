import pandas as pd
import numpy as np

df = pd.read_csv("/home/sergei/grand_sweep_1m/grand_sweep_1m.csv")
recs = df[df["is_recombinant_true"] == 1]

bins = [(0, 2, "0--2"), (3, 5, "3--5"), (6, 10, "6--10"), (11, 15, "11--15"), (16, 25, "16--25"), (26, 50, "26--50"), (51, 100, "51--100"), (101, 1000000, ">100")]
print("=== M_TRACT BINS (MATCHED BP1_ERROR) ===")
for lo, hi, label in bins:
    sub = recs[(recs["mut_in_tract"] >= lo) & (recs["mut_in_tract"] <= hi)]
    det = sub[sub["detected"] == 1]
    m_det = det[det["jaccard"] > 0]
    bp_err = m_det["bp1_error"].median() if len(m_det) > 0 else np.nan
    pwr = sub["detected"].mean() * 100.0
    mj = sub["jaccard"].mean()
    print(f"Tract SNPs $M \\in [{label}]$ & $L_{{\\mathrm{{tract}}}} \\in [45, 0.45 L]\\text{{ nt}}$ & {len(sub):,} & {len(det):,} & {pwr:.1f}\\% & {mj:.3f} & {bp_err:.1f} nt \\\\")

sub500 = recs[recs["mut_in_tract"] >= 500]
m500 = sub500[(sub500["detected"] == 1) & (sub500["jaccard"] > 0)]
print("M >= 500:", len(sub500), int(sub500["detected"].sum()), f"{sub500['detected'].mean()*100:.2f}%", f"med_bp1={m500['bp1_error'].median():.1f} nt")

print("=== COMPLEXITY & TOPOLOGY ===")
for k, name, pspace in [
    (1, "Single Event ($K_{\\mathrm{true}} = 1$)", "Crossover / Cassette"),
    (2, "Dual Events ($K_{\\mathrm{true}} = 2$)", "Two Reticulations"),
    (3, "Triple Events ($K_{\\mathrm{true}} = 3$)", "Three Reticulations"),
    (4, "Quad Events ($K_{\\mathrm{true}} = 4$)", "Complex Mosaicism"),
]:
    sub = recs[recs["num_events_true"] == k]
    det = sub[sub["detected"] == 1]
    m_det = det[det["mean_jaccard"] > 0]
    pwr = sub["detected"].mean() * 100.0
    mj = sub["mean_jaccard"].mean()
    bp = m_det["bp1_error"].median()
    ep = m_det["event_precision"].mean() * 100.0
    print(f"{name} & {pspace} & {len(sub):,} & {len(det):,} & {pwr:.1f}\\% & {mj:.3f} & {bp:.1f} nt \\\\  % event_prec={ep:.1f}%")

for name, pspace, sub in [
    ("Terminal Crossovers", "5' or 3' Exchange", recs[recs["has_crossover"] == 1]),
    ("Internal Cassettes Only", "Intercalary Swap", recs[recs["has_crossover"] == 0]),
]:
    det = sub[sub["detected"] == 1]
    m_det = det[det["mean_jaccard"] > 0]
    pwr = sub["detected"].mean() * 100.0
    mj = sub["mean_jaccard"].mean()
    bp = m_det["bp1_error"].median()
    med_j = m_det["mean_jaccard"].median()
    print(f"{name} & {pspace} & {len(sub):,} & {len(det):,} & {pwr:.1f}\\% & {mj:.3f} & {bp:.1f} nt \\\\  % matched_med_j={med_j:.3f}")
