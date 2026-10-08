#!/usr/bin/env python3
"""
generate_benchmark_dossier.py
=============================
Generates comprehensive quantitative reports, stratified analysis tables,
and publication-quality figures from the Native Rust RhizAeon simulation battery.
"""

import sys
import os
import sqlite3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def generate_dossier(db_path: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    conn = sqlite3.connect(db_path)

    # 1. Overall Regime Summary
    summary_query = """
    SELECT 
        regime,
        is_null,
        count(*) as total_trials,
        sum(detected) as detected_trials,
        round(avg(detected) * 100.0, 2) as power_or_fpr_pct,
        round(avg(correct_recombinant_exact) * 100.0, 2) as rec_exact_accuracy_pct,
        round(avg(correct_recombinant_sister) * 100.0, 2) as rec_sister_accuracy_pct,
        round(avg(correct_recombinant_in_set) * 100.0, 2) as rec_in_set_accuracy_pct,
        round(avg(correct_home) * 100.0, 2) as home_accuracy_pct,
        round(avg(correct_donor) * 100.0, 2) as donor_accuracy_pct,
        round(avg(correct_both_parents) * 100.0, 2) as joint_parent_accuracy_pct,
        sum(num_true_bps) as total_true_bps,
        sum(num_inferred_bps) as total_inferred_bps,
        sum(matched_bps_50) as matched_bps_50,
        sum(missed_bps_50) as missed_bps_50,
        round(cast(sum(missed_bps_50) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_50_pct,
        sum(matched_bps_100) as matched_bps_100,
        sum(missed_bps_100) as missed_bps_100,
        round(cast(sum(missed_bps_100) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_100_pct,
        round(avg(bp_error_mean), 2) as mean_bp_error_nt,
        round(avg(bp_error_median), 2) as median_bp_error_nt,
        round(avg(bp_all_within_35) * 100.0, 2) as bp_within_35_pct,
        round(avg(bp_all_within_50) * 100.0, 2) as bp_within_50_pct,
        round(avg(bp_all_within_100) * 100.0, 2) as bp_within_100_pct,
        round(avg(runtime_rust_ms), 2) as mean_rust_ms,
        round(avg(runtime_total_ms), 2) as mean_total_ms
    FROM rust_benchmark_results
    GROUP BY regime, is_null
    ORDER BY is_null DESC, min(global_idx) ASC;
    """
    df_regime = pd.read_sql_query(summary_query, conn)
    df_regime.to_csv(os.path.join(output_dir, "summary_by_regime.csv"), index=False)

    # 2. Stratified by Taxa Bucket (N)
    taxa_query = """
    SELECT 
        CASE 
            WHEN num_taxa <= 6 THEN 'N=4-6 (Quartets/Sextets)'
            WHEN num_taxa <= 16 THEN 'N=8-16 (Medium Clades)'
            WHEN num_taxa <= 48 THEN 'N=24-48 (Large Cohorts)'
            ELSE 'N=64-128 (Deep Surveillance)'
        END as taxa_bucket,
        is_null,
        count(*) as total_trials,
        round(avg(detected) * 100.0, 2) as power_or_fpr_pct,
        round(avg(correct_recombinant_exact) * 100.0, 2) as rec_exact_pct,
        round(avg(correct_home) * 100.0, 2) as home_pct,
        round(avg(correct_donor) * 100.0, 2) as donor_pct,
        round(avg(correct_both_parents) * 100.0, 2) as joint_parent_pct,
        round(cast(sum(missed_bps_50) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_50_pct,
        round(avg(bp_error_mean), 2) as mean_bp_error_nt,
        round(avg(runtime_rust_ms), 2) as mean_rust_ms
    FROM rust_benchmark_results
    GROUP BY taxa_bucket, is_null
    ORDER BY is_null DESC, min(num_taxa) ASC;
    """
    df_taxa = pd.read_sql_query(taxa_query, conn)
    df_taxa.to_csv(os.path.join(output_dir, "summary_by_taxa.csv"), index=False)

    # 3. Stratified by Sequence Length L
    len_query = """
    SELECT 
        length_nt,
        is_null,
        count(*) as total_trials,
        round(avg(detected) * 100.0, 2) as power_or_fpr_pct,
        round(avg(correct_recombinant_exact) * 100.0, 2) as rec_exact_pct,
        round(avg(correct_both_parents) * 100.0, 2) as joint_parent_pct,
        round(cast(sum(missed_bps_50) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_50_pct,
        round(avg(bp_error_mean), 2) as mean_bp_error_nt,
        round(avg(runtime_rust_ms), 2) as mean_rust_ms
    FROM rust_benchmark_results
    GROUP BY length_nt, is_null
    ORDER BY is_null DESC, length_nt ASC;
    """
    df_len = pd.read_sql_query(len_query, conn)
    df_len.to_csv(os.path.join(output_dir, "summary_by_length.csv"), index=False)

    # 4. Stratified by Divergence Tier
    div_query = """
    SELECT 
        CASE 
            WHEN mean_divergence < 0.005 THEN '1. Outbreak/Sparse (d < 0.005)'
            WHEN mean_divergence < 0.025 THEN '2. Low Divergence (0.005 <= d < 0.025)'
            WHEN mean_divergence < 0.080 THEN '3. Moderate Divergence (0.025 <= d < 0.08)'
            WHEN mean_divergence < 0.200 THEN '4. Deep Divergence (0.08 <= d < 0.20)'
            ELSE '5. Saturation/Noise (d >= 0.20)'
        END as divergence_tier,
        is_null,
        count(*) as total_trials,
        round(avg(detected) * 100.0, 2) as power_or_fpr_pct,
        round(avg(correct_recombinant_exact) * 100.0, 2) as rec_exact_pct,
        round(avg(correct_both_parents) * 100.0, 2) as joint_parent_pct,
        round(cast(sum(missed_bps_50) as float) / nullif(sum(num_true_bps), 0) * 100.0, 2) as missed_bps_rate_50_pct,
        round(avg(bp_error_mean), 2) as mean_bp_error_nt,
        round(avg(runtime_rust_ms), 2) as mean_rust_ms
    FROM rust_benchmark_results
    GROUP BY divergence_tier, is_null
    ORDER BY is_null DESC, divergence_tier ASC;
    """
    df_div = pd.read_sql_query(div_query, conn)
    df_div.to_csv(os.path.join(output_dir, "summary_by_divergence.csv"), index=False)

    # 5. Extract all spatial errors for detected recombinant trials
    err_query = "SELECT bp_error_mean FROM rust_benchmark_results WHERE is_null=0 AND detected=1 AND bp_error_mean IS NOT NULL"
    errors = [r[0] for r in conn.execute(err_query).fetchall()]

    conn.close()

    # -------------------------------------------------------------
    # Plot 1: Performance Dashboard (Power, FPR, Attribution, Missed BPs)
    # -------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    plt.subplots_adjust(hspace=0.35, wspace=0.25)

    # Panel A: Power and FPR across regimes
    null_df = df_regime[df_regime['is_null'] == 1]
    rec_df = df_regime[df_regime['is_null'] == 0]

    all_regimes = list(null_df['regime']) + list(rec_df['regime'])
    rates = list(null_df['power_or_fpr_pct']) + list(rec_df['power_or_fpr_pct'])
    colors = ['#d9534f' if 'null' in r else '#2e6da4' for r in all_regimes]
    labels = [r.replace('_null', ' (H0)').replace('_', ' ').title() for r in all_regimes]

    y_pos = np.arange(len(all_regimes))
    axes[0, 0].barh(y_pos, rates, color=colors, alpha=0.85, edgecolor='black', linewidth=0.8)
    axes[0, 0].axvline(5.0, color='red', linestyle='--', linewidth=1.2, label='Nominal FPR alpha=5%')
    axes[0, 0].set_yticks(y_pos)
    axes[0, 0].set_yticklabels(labels, fontsize=10)
    axes[0, 0].set_xlabel('Empirical Rate (%)', fontsize=11, fontweight='bold')
    axes[0, 0].set_title('A. Recombination Detection: Power (H1) vs FPR (H0)', fontsize=12, fontweight='bold')
    axes[0, 0].set_xlim(0, 105)
    axes[0, 0].grid(axis='x', linestyle=':', alpha=0.6)
    for i, v in enumerate(rates):
        axes[0, 0].text(v + 1.5, i, f"{v:.1f}%", va='center', fontsize=9, fontweight='bold')

    # Panel B: Lineage Attribution in Recombinant Regimes
    bar_width = 0.22
    y_rec = np.arange(len(rec_df))
    r_labels = [r.replace('_', ' ').title() for r in rec_df['regime']]

    axes[0, 1].barh(y_rec - bar_width, rec_df['rec_exact_accuracy_pct'], height=bar_width, color='#337ab7', edgecolor='black', label='Recombinant (Exact)')
    axes[0, 1].barh(y_rec, rec_df['home_accuracy_pct'], height=bar_width, color='#5cb85c', edgecolor='black', label='Home Parent')
    axes[0, 1].barh(y_rec + bar_width, rec_df['donor_accuracy_pct'], height=bar_width, color='#f0ad4e', edgecolor='black', label='Donor Parent')

    axes[0, 1].set_yticks(y_rec)
    axes[0, 1].set_yticklabels(r_labels, fontsize=10)
    axes[0, 1].set_xlabel('Attribution Accuracy (%)', fontsize=11, fontweight='bold')
    axes[0, 1].set_title('B. Parental Lineage Attribution Accuracy (H1)', fontsize=12, fontweight='bold')
    axes[0, 1].set_xlim(0, 105)
    axes[0, 1].legend(loc='lower right', frameon=True, fontsize=9)
    axes[0, 1].grid(axis='x', linestyle=':', alpha=0.6)

    # Panel C: Breakpoint Localization Error CDF
    if errors:
        sorted_errs = np.sort(errors)
        cdf = np.arange(1, len(sorted_errs) + 1) / len(sorted_errs) * 100.0
        axes[1, 0].plot(sorted_errs, cdf, color='#2e6da4', linewidth=2.5, label='Empirical CDF')
        axes[1, 0].axvline(35.0, color='#5cb85c', linestyle='--', linewidth=1.5, label='<= 35 nt resolution')
        axes[1, 0].axvline(50.0, color='#f0ad4e', linestyle='--', linewidth=1.5, label='<= 50 nt resolution')
        axes[1, 0].axvline(100.0, color='#d9534f', linestyle='--', linewidth=1.5, label='<= 100 nt resolution')
        axes[1, 0].set_xscale('log')
        axes[1, 0].set_xlabel('Breakpoint Absolute Error |u_inf - u_true| (nt, log scale)', fontsize=11, fontweight='bold')
        axes[1, 0].set_ylabel('Cumulative Percentage (%)', fontsize=11, fontweight='bold')
        axes[1, 0].set_title('C. Breakpoint Spatial Error Cumulative Distribution', fontsize=12, fontweight='bold')
        axes[1, 0].grid(True, which='both', linestyle=':', alpha=0.6)
        axes[1, 0].legend(loc='lower right', frameon=True, fontsize=9)

    # Panel D: Missed Breakpoint Rates (False Negatives)
    axes[1, 1].barh(y_rec - bar_width/2, rec_df['missed_bps_rate_50_pct'], height=bar_width, color='#d9534f', edgecolor='black', label='Missed at <= 50 nt')
    axes[1, 1].barh(y_rec + bar_width/2, rec_df['missed_bps_rate_100_pct'], height=bar_width, color='#f0ad4e', edgecolor='black', label='Missed at <= 100 nt')
    axes[1, 1].set_yticks(y_rec)
    axes[1, 1].set_yticklabels(r_labels, fontsize=10)
    axes[1, 1].set_xlabel('Missed Breakpoint Rate (%) [False Negatives]', fontsize=11, fontweight='bold')
    axes[1, 1].set_title('D. Breakpoint Omission / False Negative Rate', fontsize=12, fontweight='bold')
    axes[1, 1].set_xlim(0, 105)
    axes[1, 1].legend(loc='lower left', frameon=True, fontsize=9)
    axes[1, 1].grid(axis='x', linestyle=':', alpha=0.6)

    fig_path = os.path.join(output_dir, "fig_rust_benchmark_dashboard.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved to: {fig_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: generate_benchmark_dossier.py <db_path> <output_dir>")
        sys.exit(1)
    generate_dossier(sys.argv[1], sys.argv[2])
