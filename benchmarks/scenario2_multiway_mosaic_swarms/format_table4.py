#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

def main():
    p = Path("/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario2_multiway_mosaic_swarms/multiway_summary.csv")
    if not p.exists():
        print("Summary CSV not found.")
        return
    df = pd.read_csv(p)
    
    topologies = [
        ("dual_cassette", "Dual-Cassette ($P_1 \\to P_2 \\to P_1 \\to P_3 \\to P_1$)"),
        ("nested_russian_doll", "Nested Russian Doll ($P_1 \\to P_2 \\to P_3 \\to P_2 \\to P_1$)"),
        ("triple_crossover", "Triple Crossover Cascade ($P_1 \\to P_2 \\to P_3$)"),
    ]
    divergences = [0.02, 0.05, 0.10]
    ks = [1, 5, 10]
    
    lines = []
    for top_key, top_title in topologies:
        lines.append(f"\\multicolumn{{12}}{{l}}{{\\textbf{{{top_title}}}}} \\\\")
        for d_idx, d in enumerate(divergences):
            for k in ks:
                n = 45 + k
                r_3seq = df[(df["topology"] == top_key) & (df["d_parent"] == d) & (df["K"] == k) & (df["tool"] == "3SEQ")]
                r_rhiz = df[(df["topology"] == top_key) & (df["d_parent"] == d) & (df["K"] == k) & (df["tool"] == "RhizAeon")]
                
                if r_3seq.empty or r_rhiz.empty:
                    continue
                s_row = r_3seq.iloc[0]
                z_row = r_rhiz.iloc[0]
                
                s_k = f"{s_row['mean_unique_events']:4.1f}"
                s_fpr = f"{s_row['mean_bg_fpr'] * 100:4.1f}\\%"
                s_time = f"{int(round(s_row['mean_time_ms'])):,}"
                
                z_k = f"{z_row['mean_unique_events']:4.2f}"
                z_rec = f"{z_row['mean_recall'] * 100:4.1f}\\%"
                z_fpr = f"{z_row['mean_bg_fpr'] * 100:4.1f}\\%"
                z_prec = f"{z_row['mean_precision'] * 100:4.1f}\\%"
                z_jacc = f"{z_row['mean_jaccard']:0.3f}"
                z_time = f"{int(round(z_row['mean_time_ms'])):,}"
                
                d_str = f"{d:0.2f}"
                k_str = f"{k:2d}"
                n_str = f"{n:2d}"
                
                line = f"{d_str} & {k_str} & {n_str} & {s_k} & {s_fpr} & {s_time:>5} & {z_k} & {z_rec} & {z_fpr} & {z_prec} & {z_jacc} & {z_time:>4} \\\\"
                lines.append(line)
            if d_idx < len(divergences) - 1:
                lines.append("\\addlinespace")
        if top_key != topologies[-1][0]:
            lines.append("\\midrule")
            
    out_tex = "\n".join(lines)
    print(out_tex)
    
    # Also print overall stats
    print("\n% OVERALL AGGREGATE STATS:")
    rhiz_df = df[df["tool"] == "RhizAeon"]
    three_df = df[df["tool"] == "3SEQ"]
    print(f"% 3SEQ Mean BG FPR: {three_df['mean_bg_fpr'].mean() * 100:.2f}%")
    print(f"% 3SEQ Mean Precision: {three_df['mean_precision'].mean() * 100:.2f}%")
    print(f"% 3SEQ Mean Unique Events: {three_df['mean_unique_events'].mean():.2f}")
    print(f"% 3SEQ Mean Time: {three_df['mean_time_ms'].mean():.1f} ms")
    print(f"% RhizAeon Mean Recall: {rhiz_df['mean_recall'].mean() * 100:.2f}%")
    print(f"% RhizAeon Mean BG FPR: {rhiz_df['mean_bg_fpr'].mean() * 100:.2f}%")
    print(f"% RhizAeon Mean Precision: {rhiz_df['mean_precision'].mean() * 100:.2f}%")
    print(f"% RhizAeon Mean Unique Events: {rhiz_df['mean_unique_events'].mean():.2f}")
    print(f"% RhizAeon Mean Jaccard: {rhiz_df['mean_jaccard'].mean():.3f}")
    print(f"% RhizAeon Mean Time: {rhiz_df['mean_time_ms'].mean():.1f} ms")

if __name__ == "__main__":
    main()
