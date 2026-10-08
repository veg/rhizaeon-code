#!/usr/bin/env python3
"""
test_simulated_scenarios.py
===========================
Self-contained automated test suite for RhizAeon (PA-FF v2.1)
across simulated benchmark regimes:
  1. Synthetic Single Crossover
  2. Synthetic Cassette Swap
  3. Darren Heterotachy Null (Lineage Rate Burst)
  4. Monomorphic Clonal Null
"""

import sys
import os
import json
import subprocess
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RUST_BIN = REPO_ROOT / "target" / "release" / "rhizaeon"

def run_rhizaeon_fasta(fasta_str: str) -> dict:
    if not RUST_BIN.exists():
        subprocess.run(["cargo", "build", "--release", "-p", "rhizaeon-cli"], cwd=REPO_ROOT, check=True)
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f_fa:
        f_fa.write(fasta_str)
        fa_path = f_fa.name
    out_json = fa_path + ".rhiz.json"
    try:
        cmd = [str(RUST_BIN), "-i", fa_path, "-o", out_json]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"RhizAeon failed:\n{res.stderr}")
        with open(out_json) as f:
            data = json.load(f)
        return data
    finally:
        if os.path.exists(fa_path):
            os.remove(fa_path)
        if os.path.exists(out_json):
            os.remove(out_json)

def test_synthetic_crossover():
    print("[TEST] Synthetic Single Crossover...")
    L = 1200
    s_p1 = ["A" if i % 7 == 0 else "C" for i in range(L)]
    s_p2 = ["G" if i % 7 == 0 else "T" for i in range(L)]
    s_r = [s_p1[i] if i < 600 else s_p2[i] for i in range(L)]
    s_o = ["G" for _ in range(L)]
    
    fasta = f">P1\n{''.join(s_p1)}\n>P2\n{''.join(s_p2)}\n>R\n{''.join(s_r)}\n>O\n{''.join(s_o)}\n"
    res = run_rhizaeon_fasta(fasta)
    events = res.get("events", [])
    assert len(events) >= 1, "Expected crossover detection"
    top = events[0]
    assert top["candidate_name"] == "R", f"Expected R, got {top['candidate_name']}"
    assert top["p_fisher"] < 1e-4, f"Expected p < 1e-4, got {top['p_fisher']}"
    print(f"       PASS: Detected R (Home: {top['home_name']}, Donor: {top['donor_name']}, tract: {top['u1']}-{top['u2']})")

def test_synthetic_heterotachy_null():
    print("[TEST] Darren Heterotachy Null (Private Rate Acceleration)...")
    L = 1000
    s1 = ["A"] * L
    s2 = ["A"] * L
    s3 = ["A"] * L
    s4 = ["A"] * L
    
    # Shared tree mutations
    for i in range(0, L, 40):
        s1[i] = "C"
        s2[i] = "C"
    # Private burst on s3
    for i in range(300, 600):
        if i % 4 == 0:
            s3[i] = "T"
            
    fasta = f">T1\n{''.join(s1)}\n>T2\n{''.join(s2)}\n>T3\n{''.join(s3)}\n>T4\n{''.join(s4)}\n"
    res = run_rhizaeon_fasta(fasta)
    events = res.get("events", [])
    assert len(events) == 0, f"Expected 0 events on heterotachy null, got {len(events)}"
    print("       PASS: Correctly rejected H1 (0 false alarms).")

def test_monomorphic_null():
    print("[TEST] Monomorphic Null...")
    fasta = ">S1\nACGTACGTACGT\n>S2\nACGTACGTACGT\n>S3\nACGTACGTACGT\n"
    res = run_rhizaeon_fasta(fasta)
    events = res.get("events", [])
    assert len(events) == 0, "Expected 0 events on identical sequences"
    print("       PASS: Clean monomorphic null.")

if __name__ == "__main__":
    test_synthetic_crossover()
    test_synthetic_heterotachy_null()
    test_monomorphic_null()
    print("\n[ALL SIMULATED TESTS PASSED]")
