import pandas as pd

df = pd.read_csv("benchmarks/scenario4_micro_conversion_tracts/micro_conversion_summary.csv")
piv = df.pivot_table(index=["tract_len", "d_parent", "K", "total_taxa", "mut_in_tract"], columns="tool")

out = []
out.append(r"\begin{table*}[t!]")
out.append(r"\centering")
out.append(r"\scriptsize")
out.append(r"\setlength{\tabcolsep}{2.5pt}")
out.append(r"\caption{\textbf{Scenario 1.4: Bacterial Micro-Conversion and Scale-Space Resolution Benchmark.} Performance of \RhizAeon versus 3SEQ across 900 authentic simulation trials ($L = 5{,}000\text{ nt}$, 20 replicates across 45 parameter regimes) spanning micro-conversion tract lengths $L_{\text{tract}} \in \{50, 100, 200, 400, 800\}$ nt ($1.0\%\text{--}16.0\%$ of genome length). Each tract length is evaluated across three parental divergences ($d_{\text{parent}} \in \{0.02, 0.05, 0.10\}$) and three descendant recombinant swarm sizes ($K \in \{2, 5, 10\}$ isolates) with 30 non-recombinant background parental taxa (15 in $P_1$, 15 in $P_2$; total $N \in \{32, 35, 40\}$). Ground truth consists of $K_{\text{true}} = 1.0$ ancestral micro-conversion event. Columns report parental divergence ($d_{\text{parent}}$), descendant isolates ($K$), cohort size ($N$), expected mutations in tract (Mut), and for each engine: inferred unique events ($\hat{K}$), member isolate recall (\%), background parental false positive rate ($\text{FPR}_{\text{parent}}$, \%), spatial Jaccard overlap ($J$), and execution runtime in milliseconds (ms).}")
out.append(r"\label{tab:micro_conversion_benchmark}")
out.append(r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}} ccc c rrrrr rrrrr @{}}")
out.append(r"\toprule")
out.append(r"& & & & \multicolumn{5}{c}{\textbf{3SEQ (Triplet Mosaic Scanner)}} & \multicolumn{5}{c}{\textbf{\RhizAeon (Continuous Manifold Tracking)}} \\")
out.append(r"\cmidrule(lr){5-9} \cmidrule(lr){10-14}")
out.append(r"$d_{\text{parent}}$ & $K$ & $N$ & Mut & $\hat{K}$ & Recall & $\text{FPR}_{\text{parent}}$ & Jaccard & ms & $\hat{K}$ & Recall & $\text{FPR}_{\text{parent}}$ & Jaccard & ms \\")
out.append(r"\midrule")

current_tract = None
for (t_len, d, k, n_taxa, muts), row in piv.iterrows():
    if t_len != current_tract:
        if current_tract is not None:
            out.append(r"\midrule")
        current_tract = t_len
        pct = (t_len / 5000.0) * 100.0
        out.append(f"\multicolumn{{14}}{{l}}{{\textbf{{Tract Length $L_{{\text{{tract}}}} = {t_len}\text{{ nt}}$ ({pct:.1f}\% Genome Length)}}}} \\")
    
    k_3s = row[("unique_events", "3SEQ")]
    rec_3s = row[("recall", "3SEQ")] * 100.0
    fpr_3s = row[("bg_fpr", "3SEQ")] * 100.0
    jac_3s = row[("jaccard", "3SEQ")]
    ms_3s = row[("time_ms", "3SEQ")]
    
    k_rh = row[("unique_events", "RhizAeon")]
    rec_rh = row[("recall", "RhizAeon")] * 100.0
    fpr_rh = row[("bg_fpr", "RhizAeon")] * 100.0
    jac_rh = row[("jaccard", "RhizAeon")]
    ms_rh = row[("time_ms", "RhizAeon")]
    
    line = (f"{d:.2f} & {k:2d} & {int(n_taxa):2d} & {muts:4.1f} & "
            f"{k_3s:5.2f} & {rec_3s:5.1f}\% & {fpr_3s:5.1f}\% & {jac_3s:5.3f} & {ms_3s:5.0f} & "
            f"{k_rh:5.2f} & {rec_rh:5.1f}\% & {fpr_rh:5.1f}\% & {jac_rh:5.3f} & {ms_rh:5.0f} \\")
    out.append(line)

out.append(r"\bottomrule")
out.append(r"\end{tabular*}")
out.append(r"\end{table*}")

with open("paper/sections/table6_micro_conversion.tex", "w") as f_out:
    f_out.write("\n".join(out) + "\n")
print("Table written to paper/sections/table6_micro_conversion.tex")
