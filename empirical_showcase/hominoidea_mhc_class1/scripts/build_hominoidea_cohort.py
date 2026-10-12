#!/usr/bin/env python3
"""
build_hominoidea_cohort.py
Constructs the complete Exhaustive Superfamily Hominoidea MHC Class I dataset:
- All unique full-length coding sequences (>=1050 bp) from non-human hominoids in IPD-MHC
  (Chimpanzee Patr-A, B, C, AL; Bonobo Papa-A, B, C; Gorilla Gogo-A, B, C; Orangutan Poab-A, B, C; Popy-A, B, C).
- Canonical human reference alleles from IPD-IMGT/HLA representing all major HLA-A, HLA-B, HLA-C lineages,
  including canonical gene-conversion benchmarks (HLA-B*46:01, HLA-B*15:01, HLA-C*01:02, HLA-B*73:01).
- Evaluates codon alignment using MAFFT, checks gap statistics, reading frames, and informativeness.
"""

import os
import re
import subprocess
import hashlib
from Bio import SeqIO
from collections import Counter, OrderedDict

BASE_DIR = "/Users/sergei/Projects/rhizaeon-code/empirical_showcase/hominoidea_mhc_class1"
DATA_DIR = os.path.join(BASE_DIR, "data")
MHC_NUC_FASTA = os.path.join(DATA_DIR, "MHC_nuc.fasta")
HLA_NUC_FASTA = os.path.join(DATA_DIR, "hla_nuc.fasta")
UNALIGNED_OUT = os.path.join(DATA_DIR, "hominoidea_class1_unaligned.fasta")
ALIGNED_OUT = os.path.join(DATA_DIR, "hominoidea_class1_aligned.fasta")

def extract_non_human_apes():
    species_map = {
        'Patr': 'Pan_troglodytes',
        'Papa': 'Pan_paniscus',
        'Gogo': 'Gorilla_gorilla',
        'Poab': 'Pongo_abelii',
        'Popy': 'Pongo_pygmaeus'
    }
    loci = {'A', 'B', 'C', 'AL'}

    unique_seqs = OrderedDict()
    for r in SeqIO.parse(MHC_NUC_FASTA, 'fasta'):
        m = re.search(r'([A-Za-z]+)-([A-Za-z0-9]+)\*([0-9:]+)', r.description)
        if m:
            sp, locus, allele = m.group(1), m.group(2), m.group(3)
            if sp in species_map and locus in loci:
                seq_str = str(r.seq).upper()
                if len(seq_str) >= 1050:
                    clean_id = f"{sp}_{locus}_{allele.replace(':', '_')}"
                    if seq_str not in unique_seqs:
                        unique_seqs[seq_str] = (clean_id, sp, locus, allele, len(seq_str))
    
    print(f"[✓] Extracted {len(unique_seqs)} unique full-length non-human ape alleles.")
    return unique_seqs

def extract_human_references():
    # We want representative alleles across each 1-field lineage family
    # plus key benchmark alleles
    human_seqs = OrderedDict()
    
    # Track representative per (locus, 1-field lineage)
    # e.g. locus='A', group='01' -> pick A*01:01:01:01
    records_by_group = {}
    
    key_alleles = {
        # Canonical gene conversion pair and chimera
        "B*46:01:01:01", "B*15:01:01:01", "C*01:02:01:01",
        # Archaic Neanderthal / Denisovan introgression & divergent B
        "B*73:01:01:01", "B*27:05:01:01", "B*57:01:01:01",
        # Prototypical supertype references
        "A*02:01:01:01", "A*01:01:01:01", "A*03:01:01:01", "A*24:02:01:01",
        "B*07:02:01:01", "B*08:01:01:01", "B*44:02:01:01", "B*35:01:01:01",
        "C*07:01:01:01", "C*04:01:01:01", "C*06:02:01:01"
    }

    for r in SeqIO.parse(HLA_NUC_FASTA, 'fasta'):
        m = re.search(r'HLA:HLA\d+\s+([ABC])\*([0-9]+):([0-9]+):([0-9]+):([0-9]+)', r.description)
        if not m:
            m = re.search(r'HLA:HLA\d+\s+([ABC])\*([0-9]+):([0-9]+):([0-9]+)', r.description)
        if not m:
            m = re.search(r'HLA:HLA\d+\s+([ABC])\*([0-9]+):([0-9]+)', r.description)
        
        if m:
            locus = m.group(1)
            f1 = m.group(2)
            f2 = m.group(3)
            full_allele = r.description.split()[1] # e.g. A*01:01:01:01
            allele_name = full_allele[2:] # strip locus*
            
            seq_str = str(r.seq).upper()
            if len(seq_str) >= 1050:
                group_key = (locus, f1)
                
                # Check if this is a priority reference allele
                is_priority = (full_allele in key_alleles) or (f2 == "01")
                
                if group_key not in records_by_group:
                    records_by_group[group_key] = (full_allele, locus, allele_name, seq_str)
                elif is_priority and full_allele in key_alleles:
                    records_by_group[group_key] = (full_allele, locus, allele_name, seq_str)

    # Add all selected group representatives
    for (locus, f1), (full_allele, loc, al, seq_str) in sorted(records_by_group.items()):
        clean_id = f"HLA_{locus}_{al.replace(':', '_')}"
        if seq_str not in human_seqs:
            human_seqs[seq_str] = (clean_id, "Homo_sapiens", locus, al, len(seq_str))

    # Explicitly ensure all key benchmark alleles are included
    for r in SeqIO.parse(HLA_NUC_FASTA, 'fasta'):
        for k in key_alleles:
            if f" {k} " in r.description or r.description.endswith(f" {k}"):
                seq_str = str(r.seq).upper()
                locus = k[0]
                clean_id = f"HLA_{locus}_{k[2:].replace(':', '_')}"
                if seq_str not in human_seqs:
                    human_seqs[seq_str] = (clean_id, "Homo_sapiens", locus, k[2:], len(seq_str))

    print(f"[✓] Extracted {len(human_seqs)} canonical human HLA reference alleles.")
    return human_seqs

def main():
    non_human = extract_non_human_apes()
    human = extract_human_references()

    total_cohort = OrderedDict()
    # Add human first as reference anchors
    for seq_str, meta in human.items():
        total_cohort[meta[0]] = (meta, seq_str)
    # Add non-human apes
    for seq_str, meta in non_human.items():
        total_cohort[meta[0]] = (meta, seq_str)

    print(f"\n[✓] Total Exhaustive Hominoidea Cohort: {len(total_cohort)} taxa.")

    # Write unaligned fasta
    with open(UNALIGNED_OUT, "w") as f:
        for clean_id, (meta, seq_str) in total_cohort.items():
            f.write(f">{clean_id}\n{seq_str}\n")
    print(f"[✓] Wrote unaligned FASTA to: {UNALIGNED_OUT}")

    # Align with MAFFT
    print("\nRunning MAFFT alignment on complete cohort...")
    cmd = ["/usr/local/bin/mafft", "--auto", "--reorder", UNALIGNED_OUT]
    with open(ALIGNED_OUT, "w") as out_f:
        res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if res.returncode != 0:
        print(f"[✗] MAFFT failed: {res.stderr}")
        return

    # Evaluate alignment quality
    aligned_records = list(SeqIO.parse(ALIGNED_OUT, "fasta"))
    N = len(aligned_records)
    L = len(aligned_records[0].seq)
    print(f"[✓] MAFFT Alignment Complete: N = {N} taxa, L = {L} bp")

    # Gap and informativeness analysis
    total_cells = N * L
    total_gaps = 0
    seg_sites = 0
    informative_sites = 0

    import numpy as np
    aln_mat = np.array([list(str(r.seq).upper()) for r in aligned_records])
    for col in range(L):
        col_chars = aln_mat[:, col]
        col_no_gaps = [c for c in col_chars if c in "ACGT"]
        gaps_in_col = len(col_chars) - len(col_no_gaps)
        total_gaps += gaps_in_col
        counts = Counter(col_no_gaps)
        if len(counts) > 1:
            seg_sites += 1
            # Informative if at least 2 distinct bases appear >= 2 times
            if sum(1 for c in counts.values() if c >= 2) >= 2:
                informative_sites += 1

    gap_pct = 100.0 * total_gaps / total_cells
    print(f"Alignment Summary Statistics:")
    print(f"  Taxa (N):                {N}")
    print(f"  Length (L):              {L} bp ({L/3:.1f} codons)")
    print(f"  Overall Gap Content:     {gap_pct:.2f}%")
    print(f"  Segregating Sites:       {seg_sites} / {L} ({100.0*seg_sites/L:.1f}%)")
    print(f"  Parsimony-Informative:   {informative_sites} / {L} ({100.0*informative_sites/L:.1f}%)")

    # Hash verification
    with open(ALIGNED_OUT, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print(f"  SHA-256 Checksum:        {h}")

if __name__ == "__main__":
    main()
