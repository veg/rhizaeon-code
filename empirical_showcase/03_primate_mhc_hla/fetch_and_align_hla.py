#!/usr/bin/env python3
"""
fetch_and_align_hla.py
Fetches authentic Primate MHC / Human Leukocyte Antigen (HLA) Class I coding sequences
directly from the official EBI IPD-IMGT/HLA REST API, and aligns them using MAFFT.

Focuses on the canonical textbook case of ectopic (inter-locus) gene conversion:
The emergence of HLA-B*46:01:01 via an interlocus gene conversion between
an HLA-B recipient backbone (HLA-B*15:01 / B*62) and an HLA-C donor (HLA-C*01:02),
transferring a 73-bp segment in exon 2 (alpha-1 domain).

Primary Literature:
1. Barber LD, Percival L, Valiante NM, et al. (1996) J. Immunol. 156(9): 3275-3284.
2. Parham P, Adams EJ, Arnett KL (1995) Immunol. Rev. 143: 141-180.
3. Guethlein LA, Flodin LR, Adams EJ, Parham P (2002) J. Immunol. 169(1): 220-229.
"""

import urllib.request
import json
import time
import subprocess
import os
import hashlib
from Bio import SeqIO

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

ALLELE_TARGETS = [
    ("HLA_B_46_01_01", "HLA00331", "B*46:01:01:01", "Chimeric inter-locus recombinant allele"),
    ("HLA_B_15_01_01", "HLA00162", "B*15:01:01:01", "Recipient backbone parent (HLA-B*62 clade)"),
    ("HLA_C_01_02_01", "HLA00401", "C*01:02:01:01", "Donor parent (HLA-C*01 alpha-1 helix donor)"),
    ("HLA_A_02_01_01", "HLA00005", "A*02:01:01:01", "HLA-A paralog outgroup"),
    ("HLA_B_07_02_01", "HLA00138", "B*07:02:01:01", "HLA-B paralog control"),
    ("HLA_B_08_01_01", "HLA00147", "B*08:01:01:01", "HLA-B paralog control"),
    ("HLA_C_07_01_01", "HLA00411", "C*07:01:01:01", "HLA-C paralog control"),
    ("HLA_C_04_01_01", "HLA00405", "C*04:01:01:01", "HLA-C paralog control")
]

def fetch_imgt_hla_coding(label, acc):
    print(f"Fetching {label} ({acc}) from EBI IPD-IMGT/HLA API...")
    time.sleep(0.4)
    url = f"https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/{acc}"
    req = urllib.request.Request(url, headers={"User-Agent": "RhizAeon-Benchmarking/1.0 (academic research)", "Accept": "application/json"})
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read())
    
    coding = data["sequence"]["coding"]
    print(f"  Retrieved {data['name']}: {len(coding)} bp ({len(coding)//3} codons)")
    return coding

def main():
    unaligned_path = os.path.join(BASE_DIR, "hla_class1_unaligned.fasta")
    aligned_path = os.path.join(BASE_DIR, "hla_class1_aligned.fasta")

    sequences = {}
    for label, acc, name, role in ALLELE_TARGETS:
        seq = fetch_imgt_hla_coding(label, acc)
        sequences[label] = seq

    with open(unaligned_path, "w") as f:
        for label, seq in sequences.items():
            f.write(f">{label}\n{seq}\n")
    print(f"\nSaved unaligned sequences to: {unaligned_path}")

    # Align with MAFFT
    print("Running MAFFT alignment...")
    cmd = ["/usr/local/bin/mafft", "--auto", unaligned_path]
    with open(aligned_path, "w") as out_f:
        res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if res.returncode != 0:
        print(f"MAFFT failed: {res.stderr}")
        return

    aligned_records = list(SeqIO.parse(aligned_path, "fasta"))
    aln_len = len(aligned_records[0].seq)
    print(f"MAFFT alignment complete: {len(aligned_records)} taxa, {aln_len} bp")

    with open(aligned_path, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print(f"Alignment SHA-256: {h}")

if __name__ == "__main__":
    main()
