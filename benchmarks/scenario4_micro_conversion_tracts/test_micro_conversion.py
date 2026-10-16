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

def generate_micro_conversion_alignment(
    tract_len: int = 100,
    d_parent: float = 0.05,
    K: int = 5,
    n_p1: int = 15,
    n_p2: int = 15,
    L: int = 5000,
    drift: float = 0.002,
    seed: int = 42,
):
    rng = np.random.Generator(np.random.PCG64(seed))
    
    anc_root = rng.choice(4, size=L, p=PI)
    p1_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    p2_anc = evolve_sequence(anc_root, d_parent / 2.0, rng)
    
    taxa = []
    seqs = []
    
    for i in range(n_p1):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, drift, rng))
        
    for i in range(n_p2):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, drift, rng))
        
    u1 = 2500 - tract_len // 2
    u2 = u1 + tract_len - 1
    
    mut_in_tract = int(np.sum(p1_anc[u1 - 1 : u2] != p2_anc[u1 - 1 : u2]))
    total_mut = int(np.sum(p1_anc != p2_anc))
    
    rec_anc = p1_anc.copy()
    rec_anc[u1 - 1 : u2] = p2_anc[u1 - 1 : u2]
    
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat, u1, u2, mut_in_tract, total_mut

def compute_jaccard(true_u1: int, true_u2: int, pred_u1: int, pred_u2: int) -> float:
    start = max(true_u1, pred_u1)
    end = min(true_u2, pred_u2)
    inter = max(0, end - start + 1)
    union = (true_u2 - true_u1 + 1) + (pred_u2 - pred_u1 + 1) - inter
    return inter / union if union > 0 else 0.0

def test_evaluation():
    rhiz_bin = "rhizaeon"
    three_bin = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"
    
    for tract_len in [100, 200, 400, 800]:
        for d in [0.05, 0.10]:
            taxa, mat, u1, u2, mut_in_tract, total_mut = generate_micro_conversion_alignment(
                tract_len=tract_len, d_parent=d, K=5, seed=123
            )
            
            with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp_fa:
                for name, s in zip(taxa, mat):
                    s_str = "".join(NT_CHARS[c] for c in s)
                    tmp_fa.write(f">{name}\n{s_str}\n".encode())
                tmp_fa_path = tmp_fa.name
                
            out_json = tmp_fa_path + ".json"
            t0 = time.perf_counter()
            r_proc = subprocess.run([rhiz_bin, "--input", tmp_fa_path, "--output", out_json], capture_output=True, text=True)
            t_rhiz = (time.perf_counter() - t0) * 1000
            
            events = []
            if os.path.exists(out_json):
                with open(out_json) as f:
                    events = json.load(f).get("events", [])
                os.remove(out_json)
                
            # 3SEQ
            run_id = f"test_3s_{tract_len}_{int(d*100)}"
            t0 = time.perf_counter()
            subprocess.run([three_bin, "-full", tmp_fa_path, "-id", run_id], capture_output=True, text=True)
            t_3seq = (time.perf_counter() - t0) * 1000
            
            three_triplets = 0
            three_children = set()
            rec_csv = f"{run_id}.3s.rec.csv"
            if os.path.exists(rec_csv):
                with open(rec_csv) as f:
                    for line in f:
                        parts = [p.strip() for p in line.split(",")]
                        if len(parts) >= 10:
                            try:
                                if float(parts[9]) <= 0.05:
                                    three_triplets += 1
                                    three_children.add(parts[2])
                            except Exception:
                                pass
                for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped"]:
                    if os.path.exists(f"{run_id}{ext}"):
                        os.remove(f"{run_id}{ext}")
                        
            os.remove(tmp_fa_path)
            
            # Summarize RhizAeon
            rhiz_called_recs = set()
            rhiz_called_bg = set()
            rhiz_jaccards = []
            rhiz_bounds = []
            for ev in events:
                rhiz_bounds.append((ev.get("u1"), ev.get("u2")))
                for iso in ev.get("isolates", []):
                    if iso.startswith("Rec_"):
                        rhiz_called_recs.add(iso)
                    else:
                        rhiz_called_bg.add(iso)
                rhiz_jaccards.append(compute_jaccard(u1, u2, ev.get("u1"), ev.get("u2")))
                
            print(f"=== Tract {tract_len:3d} nt, d={d:.2f} (Mut in tract: {mut_in_tract:2d}, True: [{u1}, {u2}]) ===")
            print(f"  RhizAeon: {len(events)} events, Recs: {len(rhiz_called_recs)}/5, BG FP: {len(rhiz_called_bg)}, Bounds: {rhiz_bounds}, Jaccard: {max(rhiz_jaccards) if rhiz_jaccards else 0.0:.3f}, Time: {t_rhiz:.1f}ms")
            print(f"  3SEQ:     {three_triplets} triplets, Children: {len(three_children)}, Time: {t_3seq:.1f}ms")

if __name__ == "__main__":
    test_evaluation()
