#!/usr/bin/env python3
"""
fetch_and_align_sarbecovirus.py
Fetches authentic Sarbecovirus full-length Spike glycoprotein coding sequences
directly from NCBI GenBank, extracts exact CDS regions, and performs MAFFT alignment.

All accessions are authentic primary literature genomes:
- NC_045512.2: SARS-CoV-2 reference (Wuhan-Hu-1)
- MN996532.2: Bat coronavirus RaTG13 (Zhou et al. 2020 Nature)
- MZ937000.1: Bat coronavirus BANAL-52 (Temmam et al. 2022 Nature)
- MZ937001.1: Bat coronavirus BANAL-103 (Temmam et al. 2022 Nature)
- MT121216.1: Pangolin coronavirus MP789 (Xiao et al. 2020 Nature / Lam et al. 2020)
- MT121215.1: Bat coronavirus RmYN02 (Zhou et al. 2020 Curr Biol)
- NC_004718.3: SARS-CoV-1 reference (Tor2 / Urbani)
- KC881005.1: Bat SARS-like coronavirus RsSHC014 (Ge et al. 2013 Nature)
- DQ022305.2: Bat SARS coronavirus HKU3-1 (Lau et al. 2005 PNAS)
"""

import urllib.request
import time
import io
import subprocess
import os
import hashlib
from Bio import SeqIO

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

TARGETS = [
    ("SARS_CoV_2_Wuhan_Hu_1", "NC_045512.2", "Human pandemic index reference"),
    ("Bat_RaTG13", "MN996532.2", "Yunnan bat sarbecovirus backbone sister"),
    ("Bat_BANAL_52", "MZ937000.1", "Laos bat sarbecovirus with ACE2-binding RBM"),
    ("Bat_BANAL_103", "MZ937001.1", "Laos bat sarbecovirus with ACE2-binding RBM"),
    ("Pangolin_MP789", "MT121216.1", "Guangdong pangolin-CoV sharing ACE2 contact ridge"),
    ("Bat_RmYN02", "MT121215.1", "Yunnan bat sarbecovirus with natural S1/S2 insertions"),
    ("SARS_CoV_1_Tor2", "NC_004718.3", "Human SARS-CoV-1 pandemic reference"),
    ("Bat_RsSHC014", "KC881005.1", "Bat SARS-like coronavirus (clade 1)"),
    ("Bat_HKU3_1", "DQ022305.2", "Bat SARS coronavirus outgroup (clade 2)")
]

def fetch_spike_cds(label, acc):
    print(f"Fetching {label} ({acc}) from NCBI Entrez...")
    time.sleep(0.5)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "RhizAeon-Benchmarking/1.0 (academic research)"})
    with urllib.request.urlopen(req) as resp:
        text = resp.read().decode("utf-8")
    
    rec = SeqIO.read(io.StringIO(text), "genbank")
    for feat in rec.features:
        if feat.type == "CDS":
            prod = feat.qualifiers.get("product", [""])[0].lower()
            gene = feat.qualifiers.get("gene", [""])[0].lower()
            note = feat.qualifiers.get("note", [""])[0].lower()
            if "spike" in prod or "spike" in note or gene == "s" or prod == "s protein":
                s_seq = str(feat.extract(rec.seq))
                print(f"  Extracted Spike CDS: {len(s_seq)} bp ({len(s_seq)//3} aa)")
                return s_seq
    raise RuntimeError(f"Spike CDS not found for {acc}")

def main():
    unaligned_path = os.path.join(BASE_DIR, "sarbecovirus_spike_unaligned.fasta")
    aligned_path = os.path.join(BASE_DIR, "sarbecovirus_spike_aligned.fasta")

    sequences = {}
    for label, acc, note in TARGETS:
        seq = fetch_spike_cds(label, acc)
        sequences[label] = seq

    with open(unaligned_path, "w") as f:
        for label, seq in sequences.items():
            f.write(f">{label}\n{seq}\n")
    print(f"\nSaved unaligned sequences to: {unaligned_path}")

    # Run MAFFT alignment
    print("Running MAFFT alignment...")
    cmd = ["/usr/local/bin/mafft", "--auto", unaligned_path]
    with open(aligned_path, "w") as out_f:
        res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    
    if res.returncode != 0:
        print(f"MAFFT failed: {res.stderr}")
        return

    # Check aligned dimensions
    aligned_records = list(SeqIO.parse(aligned_path, "fasta"))
    aln_len = len(aligned_records[0].seq)
    print(f"MAFFT alignment complete: {len(aligned_records)} taxa, {aln_len} bp")

    # Compute sha256
    with open(aligned_path, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print(f"Alignment SHA-256: {h}")

if __name__ == "__main__":
    main()
