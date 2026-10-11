#!/usr/bin/env python3
"""
Diagnostic test for Scenario 1.3: Ghost Donor Introgression.
Tests RhizAeon vs 3SEQ on an alignment where the donor is completely unsampled.
"""

import os
import sys
import json
import time
import tempfile
import subprocess
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

def generate_ghost_alignment(
    d_parent: float = 0.05,
    d_ghost: float = 0.08,
    K: int = 5,
    n_p1: int = 20,
    n_p2: int = 20,
    L: int = 5000,
    u1: int = 1500,
    u2: int = 3000,
    drift: float = 0.002,
    seed: int = 42,
):
    rng = np.random.Generator(np.random.PCG64(seed))
    
    anc_root = rng.choice(4, size=L, p=PI)
    p1_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p2_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    
    # Ghost donor lineage diverges from root or P1
    ghost_anc = evolve_sequence(anc_root, d_ghost, rng)
    
    taxa = []
    seqs = []
    
    # Background P1 clade
    for i in range(n_p1):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, drift, rng))
        
    # Background P2 clade (unrelated background)
    for i in range(n_p2):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, drift, rng))
        
    # Recombinant inherits u1..u2 from ghost_anc, rest from p1_anc
    rec_anc = p1_anc.copy()
    rec_anc[u1 - 1 : u2] = ghost_anc[u1 - 1 : u2]
    
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat

def main():
    taxa, mat = generate_ghost_alignment(d_parent=0.05, d_ghost=0.08, K=5, n_p1=20, n_p2=20, seed=42)
    print(f"Alignment generated: {len(taxa)} taxa, {mat.shape[1]} nt.")
    
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tf:
        fa_path = tf.name
        for t, s in zip(taxa, mat):
            s_str = "".join(NT_CHARS[c] for c in s)
            tf.write(f">{t}\n{s_str}\n".encode())
            
    rhiz_bin = "/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon"
    three_bin = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
    out_json = fa_path + ".json"
    
    print("\n--- Running RhizAeon ---")
    t0 = time.perf_counter()
    r = subprocess.run([rhiz_bin, "-i", fa_path, "-o", out_json], capture_output=True, text=True)
    t_rhiz = (time.perf_counter() - t0) * 1000.0
    print(f"RhizAeon exit code: {r.returncode} in {t_rhiz:.1f} ms")
    if os.path.exists(out_json):
        with open(out_json) as f:
            j = json.load(f)
        events = j.get("events", [])
        print(f"Detected {len(events)} events:")
        for ev in events:
            print(f"  Event #{ev.get('event_id')}: cand={ev.get('candidate_name')}, home={ev.get('home_name')}, donor={ev.get('donor_name')}, is_ghost={ev.get('is_ghost_donor')}, tract=[{ev.get('u1')}, {ev.get('u2')}], isolates={len(ev.get('isolates', []))}")
            
    print("\n--- Running 3SEQ ---")
    t0 = time.perf_counter()
    r_3seq = subprocess.run([three_bin, "-f", fa_path, "-id", fa_path + "_3s"], capture_output=True, text=True)
    t_3seq = (time.perf_counter() - t0) * 1000.0
    print(f"3SEQ exit code: {r_3seq.returncode} in {t_3seq:.1f} ms")
    csv_path = fa_path + "_3s.3s.rec.csv"
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            lines = [l.strip() for l in f if l.strip()]
        print(f"3SEQ generated {len(lines)-1} alert lines.")
        for l in lines[1:10]:
            print(f"  {l}")

if __name__ == "__main__":
    main()
