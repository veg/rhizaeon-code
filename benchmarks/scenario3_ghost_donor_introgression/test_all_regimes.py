import os, sys, json, tempfile, subprocess, time
import numpy as np
from scipy.linalg import expm
from test_single_run import build_hky_rate_matrix, evolve_sequence, PI, NT_CHARS

def generate_scenario3_alignment(
    regime: str = "pure_ghost",
    d_parent: float = 0.05,
    K: int = 5,
    n_p1: int = 15,
    n_p2: int = 15,
    n_surr: int = 10,
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
    
    taxa = []
    seqs = []
    
    # P1 background
    for i in range(n_p1):
        taxa.append(f"P1_{i}")
        seqs.append(evolve_sequence(p1_anc, drift, rng))
        
    # P2 background
    for i in range(n_p2):
        taxa.append(f"P2_{i}")
        seqs.append(evolve_sequence(p2_anc, drift, rng))
        
    if regime == "pure_ghost":
        # Unsampled ghost donor diverging from root
        d_ghost = d_parent + 0.03
        ghost_anc = evolve_sequence(anc_root, d_ghost, rng)
        rec_donor = ghost_anc
        expected_donor_type = "ghost"
    elif regime == "surrogate_sister":
        # Donor clade ancestor diverged from root
        d_donor_clade = d_parent
        donor_clade_anc = evolve_sequence(anc_root, d_donor_clade, rng)
        delta_surr = 0.03
        ghost_anc = evolve_sequence(donor_clade_anc, delta_surr / 2.0, rng)
        surr_anc = evolve_sequence(donor_clade_anc, delta_surr / 2.0, rng)
        # Sample surrogate taxa
        for i in range(n_surr):
            taxa.append(f"Surr_{i}")
            seqs.append(evolve_sequence(surr_anc, drift, rng))
        rec_donor = ghost_anc
        expected_donor_type = "surrogate_or_ghost"
    elif regime == "deep_reservoir":
        # Deep divergent reservoir introgression
        d_reservoir = 0.15
        ghost_anc = evolve_sequence(anc_root, d_reservoir, rng)
        rec_donor = ghost_anc
        expected_donor_type = "ghost"
    else:
        raise ValueError(f"Unknown regime: {regime}")
        
    # Recombinant
    rec_anc = p1_anc.copy()
    rec_anc[u1 - 1 : u2] = rec_donor[u1 - 1 : u2]
    
    for k in range(K):
        taxa.append(f"Rec_{k}")
        seqs.append(evolve_sequence(rec_anc, drift, rng))
        
    mat = np.array(seqs, dtype=np.uint8)
    return taxa, mat, expected_donor_type

def test_single(regime, d_parent, K, seed=123):
    taxa, mat, exp_type = generate_scenario3_alignment(
        regime=regime, d_parent=d_parent, K=K, seed=seed
    )
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tf:
        fa_path = tf.name
        for t, s in zip(taxa, mat):
            s_str = "".join(NT_CHARS[c] for c in s)
            tf.write(f">{t}\n{s_str}\n".encode())
            
    rhiz_bin = "rhizaeon"
    out_json = fa_path + ".json"
    r = subprocess.run([rhiz_bin, "-i", fa_path, "-o", out_json], capture_output=True, text=True)
    events = []
    if os.path.exists(out_json):
        with open(out_json) as f:
            j = json.load(f)
        events = j.get("events", [])
        os.remove(out_json)
    os.remove(fa_path)
    
    print(f"[{regime:<16} | d={d_parent:.2f} | K={K:>2}] -> Events: {len(events)}")
    for ev in events:
        print(f"    cand={ev.get('candidate_name')}, home={ev.get('home_name')}, donor={ev.get('donor_name')}, is_ghost={ev.get('is_ghost_donor')}, tract=[{ev.get('u1')}, {ev.get('u2')}], isolates={len(ev.get('isolates', []))}")

if __name__ == "__main__":
    print("Testing regimes across divergences and swarm sizes:")
    for reg in ["pure_ghost", "surrogate_sister", "deep_reservoir"]:
        for d in [0.02, 0.05, 0.10]:
            for k in [2, 5, 10]:
                test_single(reg, d, k)
    