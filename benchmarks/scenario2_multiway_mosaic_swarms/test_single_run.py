#!/usr/bin/env python3
"""
Test script for Scenario 1.2: Multi-Way Mosaic Swarms.
Generates an alignment under Dual-Cassette, Nested Russian Doll, or Triple Crossover Cascade,
and inspects the outputs of RhizAeon and 3SEQ.
"""

import os
import sys
import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy.linalg import expm

NT_CHARS = ["A", "C", "G", "T"]
PI = np.array([0.30, 0.20, 0.20, 0.30], dtype=np.float64)
KAPPA = 2.5

def build_hky_rate_matrix(pi: np.ndarray = PI, kappa: float = KAPPA) -> np.ndarray:
    Q = np.zeros((4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            is_ts = (i in (0, 2) and j in (0, 2)) or (i in (1, 3) and j in (1, 3))
            mult = kappa if is_ts else 1.0
            Q[i, j] = mult * pi[j]
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])
    average_rate = -np.sum(pi * np.diag(Q))
    Q /= average_rate
    return Q

GLOBAL_HKY_Q = build_hky_rate_matrix()

def evolve_sequence(seq: np.ndarray, branch_length: float, rng: np.random.Generator) -> np.ndarray:
    if branch_length <= 0.0:
        return seq.copy()
    P_mat = expm(GLOBAL_HKY_Q * branch_length)
    new_seq = np.empty_like(seq)
    for state in range(4):
        mask = (seq == state)
        count = int(np.sum(mask))
        if count > 0:
            probs = np.maximum(P_mat[state, :], 0.0)
            probs /= np.sum(probs)
            new_seq[mask] = rng.choice(4, size=count, p=probs)
    return new_seq

def generate_multiway_alignment(
    topology: str = "dual_cassette",
    d_parent: float = 0.05,
    K: int = 5,
    n_p1: int = 15,
    n_p2: int = 15,
    n_p3: int = 15,
    L: int = 6000,
    drift: float = 0.002,
    seed: int = 42,
):
    rng = np.random.Generator(np.random.PCG64(seed))
    
    # 1. Root and 3 parental ancestors equidistant from root
    anc_root = rng.choice(4, size=L, p=PI)
    p1_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p2_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p3_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    
    taxa = []
    seqs = []
    
    # Background clades
    for i in range(n_p1):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, drift, rng))
    for i in range(n_p2):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, drift, rng))
    for i in range(n_p3):
        taxa.append(f"P3_{i}")
        seqs.append(evolve_sequence(p3_anc, drift, rng))
        
    # Recombinant ancestor
    rec_anc = p1_anc.copy()
    
    if topology == "dual_cassette":
        # P1 background with Cassette 1 [1200..2400] from P2 and Cassette 2 [3600..4800] from P3
        # 1-indexed coordinates: [1200, 2400] and [3600, 4800]
        rec_anc[1200 - 1 : 2400] = p2_anc[1200 - 1 : 2400]
        rec_anc[3600 - 1 : 4800] = p3_anc[3600 - 1 : 4800]
        true_tracts = [
            {"home": "P1", "donor": "P2", "u1": 1200, "u2": 2400, "type": "cassette"},
            {"home": "P1", "donor": "P3", "u1": 3600, "u2": 4800, "type": "cassette"},
        ]
    elif topology == "nested_russian_doll":
        # P1 background with Outer Cassette [1500..4500] from P2,
        # inside which [2500..3500] is replaced by P3.
        # So P1: 1..1499, P2: 1500..2499, P3: 2500..3500, P2: 3501..4500, P1: 4501..6000
        rec_anc[1500 - 1 : 4500] = p2_anc[1500 - 1 : 4500]
        rec_anc[2500 - 1 : 3500] = p3_anc[2500 - 1 : 3500]
        true_tracts = [
            {"home": "P1", "donor": "P2", "u1": 1500, "u2": 4500, "type": "outer_cassette"},
            {"home": "P2", "donor": "P3", "u1": 2500, "u2": 3500, "type": "inner_cassette"},
        ]
    elif topology == "triple_crossover":
        # P1: 1..2000, P2: 2001..4000, P3: 4001..6000
        rec_anc[2000 : 4000] = p2_anc[2000 : 4000]
        rec_anc[4000 : 6000] = p3_anc[4000 : 6000]
        true_tracts = [
            {"home": "P1", "donor": "P2", "u1": 2001, "u2": 4000, "type": "segment2"},
            {"home": "P2", "donor": "P3", "u1": 4001, "u2": 6000, "type": "segment3"},
        ]
    else:
        raise ValueError(f"Unknown topology: {topology}")
        
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat, true_tracts

def test_topology(topology: str, d_parent: float = 0.05, K: int = 5):
    print("=" * 80)
    print(f"Testing Topology: {topology} (d={d_parent}, K={K})")
    print("=" * 80)
    taxa, mat, true_tracts = generate_multiway_alignment(topology=topology, d_parent=d_parent, K=K)
    
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp_fa:
        for name, s in zip(taxa, mat):
            s_str = "".join(NT_CHARS[c] for c in s)
            tmp_fa.write(f">{name}\n{s_str}\n".encode())
        tmp_name = tmp_fa.name
        
    rhiz_bin = "/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon"
    out_json = tmp_name + ".json"
    
    # 1. Run RhizAeon
    print("\n--- Running RhizAeon ---")
    res = subprocess.run([rhiz_bin, "--input", tmp_name, "--output", out_json], capture_output=True, text=True)
    if res.returncode == 0 and os.path.exists(out_json):
        with open(out_json) as f:
            r_data = json.load(f)
        events = r_data.get("events", [])
        print(f"RhizAeon detected {len(events)} unique events:")
        for ev in events:
            print(f"  Event #{ev.get('event_id')}: Cand={ev.get('candidate_name')}, Home={ev.get('home_name')}, Donor={ev.get('donor_name')}, bp=[{ev.get('u1')}..{ev.get('u2')}], tract_len={ev.get('tract_length')}, isolates={ev.get('isolates')}")
        os.remove(out_json)
    else:
        print("RhizAeon failed:", res.stderr)
        
    # 2. Run 3SEQ
    print("\n--- Running 3SEQ ---")
    three_bin = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
    run_id = f"test_3s_{topology}_{os.getpid()}"
    for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped"]:
        p = Path(f"{run_id}{ext}")
        if p.exists():
            p.unlink()
    subprocess.run([three_bin, "-full", tmp_name, "-id", run_id], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rec_csv = f"{run_id}.3s.rec.csv"
    if os.path.exists(rec_csv):
        print(f"3SEQ generated {rec_csv}:")
        count = 0
        with open(rec_csv) as f:
            for line in f:
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 10:
                    try:
                        p_val = float(parts[9])
                        if p_val <= 0.05:
                            count += 1
                            if count <= 10:
                                print(f"  3SEQ Triplet #{count}: P_A={parts[0]}, P_B={parts[1]}, Child={parts[2]}, p={p_val}, bp={parts[12] if len(parts) >= 13 else 'N/A'}")
                    except Exception:
                        pass
        print(f"Total significant 3SEQ triplet alerts: {count}")
        for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped"]:
            p = f"{run_id}{ext}"
            if os.path.exists(p):
                os.remove(p)
    else:
        print("3SEQ produced no rec.csv")
        
    if os.path.exists(tmp_name):
        os.remove(tmp_name)

if __name__ == "__main__":
    for top in ["dual_cassette", "nested_russian_doll", "triple_crossover"]:
        test_topology(top, d_parent=0.05, K=5)
