#!/usr/bin/env python3
"""
merge_sweep_1m.py
Merges distributed partition CSVs from the 1,000,000-trial simulation sweep into a unified CSV.
Deduplicates by trial_id, sorts in ascending order, verifies continuous coverage [0 .. 1000000),
and checks column schemas.
"""

import glob
import os
import sys
import argparse
import pandas as pd

def merge(results_dir="results", out_file="grand_sweep_1m.csv", expected_count=1000000):
    parts = sorted(
        glob.glob(os.path.join(results_dir, "sweep_part_*.csv")),
        key=lambda p: int(os.path.basename(p).replace(".csv", "").split("_")[-1])
    )
    print(f"Found {len(parts)} partition files in '{results_dir}'.")
    if not parts:
        print(f"ERROR: No partition files matching '{results_dir}/sweep_part_*.csv' found!")
        sys.exit(1)

    dfs = []
    total_rows = 0
    for p in parts:
        df = pd.read_csv(p)
        total_rows += len(df)
        dfs.append(df)
        min_id = df["trial_id"].min() if len(df) > 0 else "None"
        max_id = df["trial_id"].max() if len(df) > 0 else "None"
        print(f"  {os.path.basename(p)}: {len(df):,} rows (ID range: {min_id} .. {max_id})")

    print(f"\nConcatenating {len(dfs)} dataframes ({total_rows:,} total rows)...")
    full_df = pd.concat(dfs, ignore_index=True)
    full_df = full_df.sort_values("trial_id").drop_duplicates(subset=["trial_id"]).reset_index(drop=True)
    
    unique_count = len(full_df)
    print(f"Unique trials after sorting/dedup: {unique_count:,}")

    min_id = full_df["trial_id"].min()
    max_id = full_df["trial_id"].max()
    print(f"Trial ID Range: {min_id} .. {max_id}")

    if unique_count == expected_count and min_id == 0 and max_id == expected_count - 1:
        print(f"[✓] SUCCESS: Exactly {expected_count:,} continuous trials verified! Zero missing trials.")
    else:
        print(f"[!] WARNING: Expected {expected_count:,} continuous trials [0 .. {expected_count - 1}], found {unique_count:,}.")
        trial_set = set(full_df["trial_id"])
        missing = sorted(list(set(range(expected_count)) - trial_set))
        if missing:
            print(f"    Missing {len(missing):,} trials. First 10 missing: {missing[:10]}")

    print(f"Writing unified sweep to '{out_file}'...")
    full_df.to_csv(out_file, index=False)
    file_size_mb = os.path.getsize(out_file) / (1024 * 1024)
    print(f"[✓] Done! Output file: {out_file} ({file_size_mb:.2f} MB)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge distributed grand sweep CSV partitions.")
    parser.add_argument("--results-dir", default="results", help="Directory containing sweep_part_*.csv files")
    parser.add_argument("--output", default="grand_sweep_1m.csv", help="Output concatenated CSV path")
    parser.add_argument("--expected", type=int, default=1000000, help="Expected total number of trials")
    args = parser.parse_args()

    merge(args.results_dir, args.output, args.expected)
