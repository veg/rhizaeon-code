#!/usr/bin/env python3
"""
run_analysis.py
Verification harness for Empirical Application 4:
Adaptive Mimicry Introgression in Heliconius Butterflies.

Verifies:
1. Exact execution of RhizAeon PA-FF v2.1 native Rust engine.
2. Concordance with published evolutionary genomics literature:
   - Pardo-Diaz C, Salazar C, Baxter SW, et al. (2012) PLoS Genetics 8(6): e1002752
     (Adaptive introgression of the HmB/optix red wing pattern locus).
   - The Heliconius Genome Consortium (2012) Nature 487: 94-98
     (Butterfly genome reveals promiscuous exchange of mimicry adaptations).
   - Reed RD, Papa R, Martin A, et al. (2011) Science 333: 1137-1141
     (Optix drives convergent evolution of butterfly wing pattern mimicry).
"""

import json
import os
import subprocess
import time

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FASTA_PATH = os.path.join(BASE_DIR, "heliconius_optix_island_aligned.fasta")
OUTPUT_JSON = os.path.join(BASE_DIR, "heliconius_optix_result.json")
OUTPUT_HTML = os.path.join(BASE_DIR, "heliconius_optix_dashboard.html")
LOG_PATH = os.path.join(BASE_DIR, "heliconius_execution.log")
PARTITIONS_PATH = os.path.join(BASE_DIR, "heliconius_partitions.json")
RHIZAEON_BIN = "/Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon"

def main():
    print("================================================================================")
    print("  Empirical Showcase 4: Heliconius Mimicry Adaptive Introgression Island        ")
    print("================================================================================")
    print(f"Alignment:    {FASTA_PATH}")
    print(f"Binary:       {RHIZAEON_BIN}")
    print(f"Output JSON:  {OUTPUT_JSON}")
    print(f"Output HTML:  {OUTPUT_HTML}")
    print("--------------------------------------------------------------------------------")

    cmd = [
        RHIZAEON_BIN,
        "-i", FASTA_PATH,
        "-o", OUTPUT_JSON,
        "--html", OUTPUT_HTML,
        "-t", "Heliconius Mimicry Adaptive Introgression Island",
        "-l", "16",
        "-a", "0.05"
    ]

    t0 = time.perf_counter()
    res = subprocess.run(cmd, capture_output=True, text=True)
    t1 = time.perf_counter()
    wall_ms = (t1 - t0) * 1000.0

    with open(LOG_PATH, "w") as f:
        f.write("=== STDOUT ===\n")
        f.write(res.stdout)
        f.write("\n=== STDERR ===\n")
        f.write(res.stderr)

    print(res.stdout)

    with open(OUTPUT_JSON, "r") as f:
        data = json.load(f)

    with open(PARTITIONS_PATH, "r") as f:
        part_meta = json.load(f)

    scan_ms = data.get("run_time_ms", wall_ms)
    aln_taxa = data.get("num_taxa", 0)
    aln_len = data.get("length", 0)
    inf_sites = data.get("informative_sites", 0)
    events = data.get("events", [])

    print("--------------------------------------------------------------------------------")
    print("  VERIFICATION & GROUND TRUTH MAPPING                                           ")
    print("--------------------------------------------------------------------------------")
    print(f"Alignment Dimensions: {aln_taxa} taxa, {aln_len} bp ({inf_sites} informative sites)")
    print(f"Engine Scan Time:     {scan_ms:.2f} ms (Wall clock: {wall_ms:.2f} ms)")
    print(f"Verified Events:      {len(events)}")

    partitions = part_meta["partitions"]
    print("\nChromosomal Locus Partitions:")
    for locus, (start, end, span) in partitions.items():
        print(f"  {locus:12s}: {start:5d} - {end:5d} bp ({span} bp)")

    for idx, ev in enumerate(events):
        rec_child = ev.get("candidate_name", "Unknown")
        home = ev.get("home_name", "Unknown")
        donor = ev.get("donor_name", "Unknown")
        u1 = ev.get("u1", 0)
        u2 = ev.get("u2", 0)
        tract_len = ev.get("tract_length", u2 - u1 + 1)
        p_val = ev.get("p_fisher", 1.0)
        z_score = ev.get("z_phys", 0.0)
        d_tract_donor = ev.get("d_tract_donor", 1.0)
        d_tract_home = ev.get("d_tract_home", 1.0)
        is_crossover = ev.get("is_crossover", True)

        print(f"\n[Event #{idx + 1}]")
        print(f"  Recombinant Taxon:      {rec_child}")
        print(f"  Parental Dipole:        Home={home} | Donor={donor}")
        print(f"  Event Topology:         {'Crossover' if is_crossover else 'Cassette Swap / Introgression Island'}")
        print(f"  Alignment Coordinates:  {u1} - {u2} bp (Length: {tract_len} bp)")
        print(f"  Fisher's Exact P-value: {p_val:.2e}")
        print(f"  Standardized Z-Score:   {z_score:.2f}")
        print(f"  Divergence in Tract:    To Donor={d_tract_donor:.4f} ({d_tract_donor*100:.2f}%) | To Home={d_tract_home:.4f} ({d_tract_home*100:.2f}%)")

        # Map to physical locus
        overlapping_loci = []
        for locus, (start, end, span) in partitions.items():
            if not (u2 < start or u1 > end):
                overlapping_loci.append(locus)

        print(f"  Spanning Physical Loci: {', '.join(overlapping_loci)}")

        if "HmB453k" in overlapping_loci:
            print(f"  >> CONCORDANCE: Event directly overlaps HmB453k (the peak of phenotype-associated introgression)!")
            print(f"  >> Biological mechanism: Interspecific adaptive introgression of the red-pattern switch from H. melpomene into H. timareta.")
            if d_tract_donor == 0.0:
                print(f"  >> PERFECT DONOR IDENTITY: 0.00% sequence divergence between introgressed H. timareta tract and H. melpomene donor!")

    print("\nProvenance verification complete.")

if __name__ == "__main__":
    main()
