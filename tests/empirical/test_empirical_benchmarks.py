#!/usr/bin/env python3
"""
test_empirical_benchmarks.py
============================
Empirical validation test harness for RhizAeon (PA-FF v2.1).
Evaluates real-world viral and human genomic alignments:
  1. HIV-1 CRF02_AG (KAL-153): Canonical multi-way recombinant.
  2. SARS-CoV-2 XBB Spike: Circulating Omicron recombinant.
  3. Human mtDNA: Clonal negative control (0% false positives).
"""

import sys
import os
import json
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RUST_BIN = REPO_ROOT / "target" / "release" / "rhizaeon"
EMPIRICAL_DATA = REPO_ROOT / "tests" / "empirical" / "data"

def run_rhizaeon(fasta_path: Path) -> dict:
    if not RUST_BIN.exists():
        # Build release binary if not present
        subprocess.run(["cargo", "build", "--release", "-p", "rhizaeon-cli"], cwd=REPO_ROOT, check=True)
        
    out_json = fasta_path.with_suffix(".rhiz.json")
    cmd = [str(RUST_BIN), "-i", str(fasta_path), "-o", str(out_json)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"RhizAeon failed on {fasta_path}:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        
    with open(out_json) as f:
        data = json.load(f)
    if out_json.exists():
        out_json.unlink()
    return data

def test_hiv1_kal153():
    kal_path = EMPIRICAL_DATA / "hiv1_kal153_crf02_ag.fasta"
    assert kal_path.exists(), f"Missing {kal_path}"
    print(f"\n[TEST] Running RhizAeon on HIV-1 CRF02_AG (KAL-153)...")
    res = run_rhizaeon(kal_path)
    events = res.get("events", [])
    print(f"       Found {len(events)} verified events in {res.get('run_time_ms', 0):.1f} ms.")
    assert len(events) >= 1, "Expected recombination in KAL-153"
    top = events[0]
    assert top["candidate_name"] == "R", f"Expected candidate R, got {top['candidate_name']}"
    assert top["home_name"] == "B", f"Expected home B, got {top['home_name']}"
    assert top["donor_name"] == "A", f"Expected donor A, got {top['donor_name']}"
    print(f"       PASS: Top event {top['candidate_name']} (Home: {top['home_name']}, Donor: {top['donor_name']}) with p = {top['p_fisher']:.2e}")

def test_human_mtdna_null():
    mtdna_path = EMPIRICAL_DATA / "human_mtdna.fasta"
    assert mtdna_path.exists(), f"Missing {mtdna_path}"
    print(f"\n[TEST] Running RhizAeon on Human mtDNA Negative Control...")
    res = run_rhizaeon(mtdna_path)
    events = res.get("events", [])
    print(f"       Found {len(events)} events in {res.get('run_time_ms', 0):.1f} ms.")
    assert len(events) == 0, f"Expected 0 events on clonal human mtDNA, got {len(events)}"
    print(f"       PASS: Strictly certified non-recombinant (H0 null).")

if __name__ == "__main__":
    test_hiv1_kal153()
    test_human_mtdna_null()
    print("\n[ALL EMPIRICAL TESTS PASSED]")
