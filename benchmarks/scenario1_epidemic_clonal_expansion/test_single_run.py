#!/usr/bin/env python3
"""
Test single run of Epidemic Clonal Expansion simulation.
"""
import os
import sys
import json
import time
import tempfile
import subprocess
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
            is_transition = (i in (0, 2) and j in (0, 2)) or (i in (1, 3) and j in (1, 3))
            mult = kappa if is_transition else 1.0
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
            probs = P_mat[state, :]
            probs = np.maximum(probs, 0.0)
            probs /= np.sum(probs)
            new_seq[mask] = rng.choice(4, size=count, p=probs)
    return new_seq

def generate_epidemic_alignment(
    d_parent: float = 0.05,
    K: int = 25,
    L: int = 5000,
    u1: int = 1500,
    u2: int = 3000,
    seed: int = 42
):
    rng = np.random.Generator(np.random.PCG64(seed))
    
    # 1. Root sequence
    anc_root = rng.choice(4, size=L, p=PI)
    
    # 2. Parental ancestors separated by d_parent (split branch length d_parent / 2 each)
    p1_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p2_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    
    taxa = []
    seqs = []
    
    # 3. Background P1 taxa (10 isolates) with small intra-clade drift (0.002)
    for i in range(10):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, 0.002, rng))
        
    # 4. Background P2 taxa (10 isolates) with small intra-clade drift (0.002)
    for i in range(10):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, 0.002, rng))
        
    # 5. Ancestral recombinant node A_rec
    # Cassette from P2 (u1 to u2, 1-indexed, so u1-1 to u2 0-indexed)
    rec_anc = p1_anc.copy()
    rec_anc[u1-1:u2] = p2_anc[u1-1:u2]
    
    # 6. Sprout K descendant recombinant isolates with terminal drift (0.002)
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, 0.002, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat

def main():
    taxa, mat = generate_epidemic_alignment(d_parent=0.05, K=25, seed=123)
    print(f"Generated alignment with {len(taxa)} taxa, {mat.shape[1]} nt.")
    
    # Write fasta
    tmp_fa = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    for name, s in zip(taxa, mat):
        s_str = "".join(NT_CHARS[c] for c in s)
        tmp_fa.write(f">{name}\n{s_str}\n".encode())
    tmp_fa.close()
    
    rhizaeon_bin = "/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon"
    out_json = tmp_fa.name + ".json"
    t0 = time.perf_counter()
    subprocess.run([rhizaeon_bin, "--input", tmp_fa.name, "--output", out_json], check=True, stdout=subprocess.DEVNULL)
    t_rhiz = (time.perf_counter() - t0) * 1000.0
    
    with open(out_json) as f:
        res = json.load(f)
    print(f"RhizAeon completed in {t_rhiz:.1f} ms.")
    print(f"RhizAeon detected {len(res['events'])} consolidated events.")
    for ev in res['events']:
        print(f"  Event #{ev['event_id']}: candidate={ev['candidate_name']}, home={ev['home_name']}, donor={ev['donor_name']}, tract=[{ev['u1']}, {ev['u2']}], isolates={len(ev['isolates'])}: {ev['isolates'][:5]}...")

    # Run 3SEQ
    three_seq_bin = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
    run_id = f"test_{os.getpid()}"
    t0 = time.perf_counter()
    cmd = [three_seq_bin, "-full", tmp_fa.name, "-id", run_id]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    t_3seq = (time.perf_counter() - t0) * 1000.0
    
    rec_csv = f"{run_id}.3s.rec.csv"
    sig_triplets = 0
    rec_children = set()
    if os.path.exists(rec_csv):
        with open(rec_csv) as f:
            for line in f:
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 10:
                    try:
                        p_val = float(parts[9])
                        if p_val <= 0.05:
                            sig_triplets += 1
                            rec_children.add(parts[2])
                    except Exception:
                        pass
        # Clean up
        for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped"]:
            p = Path(f"{run_id}{ext}")
            if p.exists():
                p.unlink()

    print(f"3SEQ completed in {t_3seq:.1f} ms.")
    print(f"3SEQ raw significant triplets: {sig_triplets}")
    print(f"3SEQ unique recombinant children: {len(rec_children)}: {sorted(list(rec_children))[:5]}...")

    os.remove(tmp_fa.name)
    os.remove(out_json)

if __name__ == "__main__":
    main()
