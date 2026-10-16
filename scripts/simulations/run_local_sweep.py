#!/usr/bin/env python3
"""
run_local_sweep.py
Runs a local multi-threaded synthetic simulation parameter sweep using the native RhizAeon grand_sweep binary.
"""

import os
import sys
import argparse
import subprocess
import time
import shutil

def find_binary():
    # 1. In PATH
    b = shutil.which("grand_sweep")
    if b:
        return b

    # 2. Built binary inside a HyphAeon monorepo checkout (HYPHAEON_ROOT)
    hypha = os.environ.get("HYPHAEON_ROOT")
    if hypha:
        for candidate in [
            os.path.join(hypha, "rust", "target", "release", "grand_sweep"),
            os.path.join(hypha, "rust", "target", "debug", "grand_sweep"),
        ]:
            if os.path.exists(candidate):
                return candidate
    return None

def main():
    parser = argparse.ArgumentParser(description="Execute a local synthetic parameter sweep with grand_sweep.")
    parser.add_argument("--trials", type=int, default=1000, help="Number of trials to simulate (default: 1000)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4, help="Number of worker threads")
    parser.add_argument("--start-trial", type=int, default=0, help="Starting trial ID")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed")
    parser.add_argument("--output", default="local_sweep.csv", help="Output CSV path")
    args = parser.parse_args()

    bin_path = find_binary()
    if not bin_path:
        print("ERROR: grand_sweep binary not found! Build it in the HyphAeon monorepo with:")
        print("  cargo build -p rhizaeon-tools --release")
        print("  (then put grand_sweep on PATH or set HYPHAEON_ROOT to the checkout)")
        sys.exit(1)

    print("=" * 70)
    print("RHIZAEON SYNTHETIC SIMULATION LOCAL SWEEP")
    print(f"Binary:       {bin_path}")
    print(f"Trials:       {args.trials:,}")
    print(f"Threads:      {args.threads}")
    print(f"Start Trial:  {args.start_trial}")
    print(f"Seed:         {args.seed}")
    print(f"Output:       {args.output}")
    print("=" * 70)

    cmd = [
        bin_path,
        "--start-trial", str(args.start_trial),
        "--trials", str(args.trials),
        "--threads", str(args.threads),
        "--output", args.output,
        "--seed", str(args.seed),
    ]

    t0 = time.perf_counter()
    res = subprocess.run(cmd)
    t1 = time.perf_counter()

    if res.returncode == 0:
        elapsed = t1 - t0
        rate = args.trials / max(0.001, elapsed)
        print("=" * 70)
        print(f"[✓] Sweep completed successfully in {elapsed:.2f} s ({rate:.1f} trials/s)")
        if os.path.exists(args.output):
            print(f"    Saved: {args.output} ({os.path.getsize(args.output)/(1024*1024):.2f} MB)")
        print("=" * 70)
    else:
        print(f"[!] Simulation failed with exit code {res.returncode}")
        sys.exit(res.returncode)

if __name__ == "__main__":
    main()
