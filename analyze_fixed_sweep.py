import pandas as pd
import numpy as np
import sys

def analyze(csv_path):
    print(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)
    print(f"Total trials loaded: {len(df)}")

    null_df = df[df['is_recombinant_true'] == 0]
    rec_df = df[df['is_recombinant_true'] == 1]
    print(f"Null trials: {len(null_df)}, Recombinant trials: {len(rec_df)}")

    # Null FPR Analysis
    print("\n" + "="*70)
    print("NULL PARTITION: EMPIRICAL FALSE POSITIVE RATE")
    print("="*70)
    total_null = len(null_df)
    null_fps = (null_df['detected'] > 0).sum()
    print(f"Overall Null FPR: {null_fps}/{total_null} ({null_fps/total_null*100:.3f}%)")

    clock_null = null_df[null_df['heterotachy_alpha'] >= 50]
    mild_null = null_df[(null_df['heterotachy_alpha'] >= 1.0) & (null_df['heterotachy_alpha'] < 50)]
    severe_null = null_df[null_df['heterotachy_alpha'] < 1.0]

    for name, sub in [("Strict Clock (alpha >= 50)", clock_null),
                      ("Mild Heterotachy (1.0 <= alpha < 50)", mild_null),
                      ("Severe Heterotachy (alpha < 1.0)", severe_null)]:
        if len(sub) > 0:
            fps = (sub['detected'] > 0).sum()
            print(f"  {name:38}: {fps:4}/{len(sub):4} ({fps/len(sub)*100:.3f}%)")

    # Recombinant Power by M_tract
    print("\n" + "="*70)
    print("RECOMBINANT PARTITION: ALIGNMENT RECALL & TAXON RECALL BY M_TRACT")
    print("="*70)
    bins = [0, 2, 5, 10, 15, 25, 50, 100, 100000]
    labels = ["[0, 2] (Floor)", "[3, 5]", "[6, 10]", "[11, 15]", "[16, 25]", "[26, 50]", "[51, 100]", "> 100 (Saturation)"]
    rec_df = rec_df.copy()
    rec_df['m_bin'] = pd.cut(rec_df['mut_in_tract'], bins=bins, labels=labels, right=True, include_lowest=True)

    header = f"{'Stratum':20} | {'Trials':6} | {'Align Recall':12} | {'Taxon Recall':12} | {'Precision':10} | {'Med BP1':8} | {'Med BP2':8} | {'Mean Jaccard':12}"
    print(header)
    print("-" * len(header))

    for label in labels:
        sub = rec_df[rec_df['m_bin'] == label]
        if len(sub) == 0:
            continue
        trials = len(sub)
        aln_recall = (sub['detected'] > 0).mean() * 100.0
        taxon_recall = sub['recall'].mean() * 100.0
        precision = sub[sub['detected'] > 0]['precision'].mean() * 100.0 if (sub['detected'] > 0).sum() > 0 else 0.0
        det_sub = sub[sub['detected'] > 0]
        med_bp1 = det_sub['bp1_error'].median() if len(det_sub) > 0 else 0.0
        med_bp2 = det_sub['bp2_error'].median() if len(det_sub) > 0 else 0.0
        mean_jaccard = det_sub['jaccard'].mean() if len(det_sub) > 0 else 0.0

        print(f"{label:20} | {trials:6} | {aln_recall:11.2f}% | {taxon_recall:11.2f}% | {precision:9.2f}% | {med_bp1:7.1f}nt | {med_bp2:7.1f}nt | {mean_jaccard:11.3f}")

    # Power for M >= 50 by cohort size N
    print("\n" + "="*70)
    print("POWER FOR HIGH-INFORMATION STRATUM (M_tract >= 50) BY COHORT SIZE N")
    print("="*70)
    high_m = rec_df[rec_df['mut_in_tract'] >= 50].copy()
    n_bins = [0, 8, 15, 30, 60, 300]
    n_labels = ["N in [4, 8]", "N in [9, 15]", "N in [16, 30]", "N in [31, 60]", "N > 60"]
    high_m['n_bin'] = pd.cut(high_m['num_taxa'], bins=n_bins, labels=n_labels, right=True, include_lowest=True)

    header_n = f"{'Cohort Size Stratum':20} | {'Trials':6} | {'Align Recall':12} | {'Taxon Recall':12} | {'Precision':10}"
    print(header_n)
    print("-" * len(header_n))
    for label in n_labels:
        sub = high_m[high_m['n_bin'] == label]
        if len(sub) == 0:
            continue
        trials = len(sub)
        aln_recall = (sub['detected'] > 0).mean() * 100.0
        taxon_recall = sub['recall'].mean() * 100.0
        precision = sub[sub['detected'] > 0]['precision'].mean() * 100.0 if (sub['detected'] > 0).sum() > 0 else 0.0
        print(f"{label:20} | {trials:6} | {aln_recall:11.2f}% | {taxon_recall:11.2f}% | {precision:9.2f}%")

    # High-information & high-divergence check (M >= 50, d >= 0.02)
    clean_stratum = rec_df[(rec_df['mut_in_tract'] >= 50) & (rec_df['d_parent'] >= 0.02)]
    print("\n" + "="*70)
    print("HIGH-INFORMATION & MODERATE DIVERGENCE (M_tract >= 50, d_parent >= 0.02)")
    print("="*70)
    if len(clean_stratum) > 0:
        det = (clean_stratum['detected'] > 0).sum()
        tot = len(clean_stratum)
        print(f"Alignment Recall: {det}/{tot} ({det/tot*100:.2f}%)")
        print(f"Mean Taxon Recall: {clean_stratum['recall'].mean()*100:.2f}%")
        print(f"Mean Precision:    {clean_stratum['precision'].mean()*100:.2f}%")

if __name__ == '__main__':
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'test_sweep_2k_fixed.csv'
    analyze(csv_file)
