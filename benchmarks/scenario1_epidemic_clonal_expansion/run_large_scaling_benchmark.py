#!/usr/bin/env python3
"""
Large-Alignment Runtime & Accuracy Scaling Benchmark (N = 100 to 1,500 sequences).
Compares RhizAeon vs 3SEQ on epidemic recombinant cohorts (L = 5,000 nt, d_parent = 0.05).
"""
import os
import sys
import time
import json
import tempfile
import subprocess
from pathlib import Path

sys.path.insert(0, "/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion")
from run_benchmark import generate_epidemic_alignment, NT_CHARS, compute_jaccard

RHIZ_BIN = "/Users/sergei/Development/testing/HyphAeon/rust/target/release/rhizaeon"
THREE_BIN = "/Users/sergei/Projects/TOGA_MEME/recombination/benchmarks/tools/3seq/cmake-build/3seq"

def run_scaling():
    configs = [
        (50, 10, 20),
        (100, 10, 45),
        (150, 15, 67),
        (250, 20, 115),
        (500, 50, 225),
        (1000, 100, 450),
        (1500, 150, 675),
    ]
    out_path = "/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/large_alignment_scaling.json"
    existing_3seq = {}
    if os.path.exists(out_path):
        try:
            with open(out_path) as f:
                for item in json.load(f):
                    existing_3seq[item["N"]] = item
        except Exception:
            pass
    rows = []
    for N_total, K, n_p in configs:
        taxa, mat = generate_epidemic_alignment(
            d_parent=0.05, K=K, n_p1=n_p, n_p2=n_p, L=5000, u1=1500, u2=3000, drift=0.005, seed=42 + N_total
        )
        tmp_fa = tempfile.NamedTemporaryFile(suffix=f"_{N_total}.fasta", delete=False)
        for name, s in zip(taxa, mat):
            seq_str = "".join(NT_CHARS[c] for c in s)
            tmp_fa.write(f">{name}\n{seq_str}\n".encode())
        tmp_fa.close()

        # Run RhizAeon
        out_json = tmp_fa.name + ".json"
        t0 = time.perf_counter()
        subprocess.run([RHIZ_BIN, "--input", tmp_fa.name, "--output", out_json], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        t_rhiz = time.perf_counter() - t0
        with open(out_json) as f:
            r_res = json.load(f)
        os.remove(out_json)

        evs = r_res.get("events", [])
        true_recs = {f"Rec_{k}" for k in range(K)}
        called_recs = set()
        called_bg = set()
        for ev in evs:
            for iso in ev.get("isolates", []):
                if iso in true_recs:
                    called_recs.add(iso)
                else:
                    called_bg.add(iso)
        recall = len(called_recs) / K
        prec = len(called_recs) / max(1, len(called_recs) + len(called_bg))

        # Run 3SEQ (reuse already-measured 3SEQ timing on identical seed if present; once 3SEQ exceeds 240s at smaller N, its O(N^3) loop cannot finish in <240s at larger N)
        run_id = f"scale3s_{N_total}"
        t0 = time.perf_counter()
        timed_out = False
        triplets_3seq = 0
        if N_total in existing_3seq:
            t_3seq = existing_3seq[N_total]["3seq_sec"]
            timed_out = existing_3seq[N_total]["3seq_timed_out"]
            triplets_3seq = existing_3seq[N_total]["3seq_triplets"]
        elif any(r.get("3seq_timed_out") for r in rows):
            t_3seq = 240.0
            timed_out = True
        else:
            try:
                subprocess.run([THREE_BIN, "-full", tmp_fa.name, "-id", run_id], input=b"Y\n", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=240)
                t_3seq = time.perf_counter() - t0
            except subprocess.TimeoutExpired:
                t_3seq = 240.0
                timed_out = True

            rec_csv = f"{run_id}.3s.rec.csv"
            if os.path.exists(rec_csv):
                with open(rec_csv) as f:
                    for line in f:
                        parts = [p.strip() for p in line.split(",")]
                        if len(parts) >= 10 and parts[0] != "P_name":
                            try:
                                if float(parts[9]) <= 0.05:
                                    triplets_3seq += 1
                            except Exception:
                                pass
            for ext in [".3s.rec.csv", ".3s.pvalHist", ".3s.log", ".3s.skipped", ".3s.longRec"]:
                p = Path(f"{run_id}{ext}")
                if p.exists():
                    p.unlink()

        os.remove(tmp_fa.name)
        speedup = t_3seq / t_rhiz
        print(f"N={N_total:4d} (K={K:3d}): RhizAeon={t_rhiz:6.2f}s (events={len(evs)}, rec={recall*100:.1f}%, prec={prec*100:.1f}%) | 3SEQ={t_3seq:7.2f}s ({'TIMEOUT' if timed_out else f'triplets={triplets_3seq}'}) | Speedup={speedup:6.1f}x", flush=True)
        rows.append({
            "N": N_total, "K": K, "rhizaeon_sec": t_rhiz, "rhizaeon_events": len(evs),
            "rhizaeon_recall": recall, "rhizaeon_precision": prec,
            "3seq_sec": t_3seq, "3seq_timed_out": timed_out, "3seq_triplets": triplets_3seq,
            "speedup": speedup
        })

    out_path = "/Users/sergei/Projects/rhizaeon-code/benchmarks/scenario1_epidemic_clonal_expansion/large_alignment_scaling.json"
    with open(out_path, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"Saved {out_path}")

if __name__ == "__main__":
    run_scaling()
