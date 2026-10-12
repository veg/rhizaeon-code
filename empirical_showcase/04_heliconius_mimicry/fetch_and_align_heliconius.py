#!/usr/bin/env python3
"""
fetch_and_align_heliconius.py
Fetches authentic Heliconius butterfly sequences across the Chromosome 18
optix / HmB adaptive mimicry introgression island directly from NCBI GenBank.

Source Study:
Pardo-Diaz C, Salazar C, Baxter SW, Merot C, Linares M, Jiggins CD, McMillan WO (2012).
Adaptive Introgression across Species Boundaries in Heliconius Butterflies.
PLoS Genetics 8(6): e1002752. doi:10.1371/journal.pgen.1002752. PMID: 22737081.

Chromosome 18 Loci Ordered Along the Physical Introgression Island:
1. rabkinesin  (~350 kb) - Upstream flanking control locus
2. HmB00012    (~400 kb) - Upstream transition marker
3. HmB449k     ( 449 kb) - Introgression boundary marker
4. HmB453k     ( 453 kb) - Peak phenotype-associated cis-regulatory element
5. optix       ( 515 kb) - Red wing-patterning master transcription factor
6. HmB520k     ( 520 kb) - Downstream flanking control marker
"""

import os
import time
import urllib.request
import subprocess
import json
import hashlib
from Bio import SeqIO
from io import StringIO

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

VOUCHERS = [
    # Introgressed Rayed Mimics (Recombinant)
    ("H_tim_tim_9227", "9227", "Heliconius timareta timareta", "Peruvian rayed mimic"),
    ("H_tim_tim_9185", "9185", "Heliconius timareta timareta", "Peruvian rayed mimic"),
    ("H_tim_tim_9182", "9182", "Heliconius timareta timareta", "Peruvian rayed mimic"),
    # Co-mimetic Donors (Red-pattern source)
    ("H_mel_mal_2312", "2312", "Heliconius melpomene malleti", "Amazonian rayed co-mimic donor"),
    ("H_mel_mal_1261", "1261", "Heliconius melpomene malleti", "Amazonian rayed co-mimic donor"),
    ("H_mel_agl_1385", "1385", "Heliconius melpomene aglaope", "Amazonian rayed co-mimic donor"),
    ("H_mel_ros_1023", "1023", "Heliconius melpomene rosina", "Central American postman race"),
    # Non-mimetic Recipient Clade (Home species)
    ("H_cyd_cor_2182", "2182", "Heliconius cydno cordula", "White/yellow sister"),
    ("H_cyd_cor_2253", "2253", "Heliconius cydno cordula", "White/yellow sister"),
    ("H_cyd_chi_2424", "2424", "Heliconius cydno chioneus", "White/yellow sister"),
    ("H_cyd_wey_1911", "1911", "Heliconius cydno weymeri", "Polymorphic sister"),
    # Silvaniform Outgroup
    ("H_num_2387",     "2387", "Heliconius numata", "Silvaniform outgroup")
]

LOCUS_ORDER = ["rabkinesin", "HmB00012", "HmB449k", "HmB453k", "optix", "HmB520k"]

# Exact GenBank accession mapping verified from Pardo-Diaz et al. (2012)
ACCESSIONS = {
    "2312": {
        "rabkinesin": "JX004173",
        "HmB00012":   "JX005704",
        "HmB449k":    "JX005487",
        "HmB453k":    "JX005268",
        "optix":      "JX004389",
        "HmB520k":    "JX005049"
    },
    "1261": {
        "rabkinesin": "JX004171",
        "HmB00012":   "JX005702",
        "HmB449k":    "JX005485",
        "HmB453k":    "JX005266",
        "optix":      "JX004416",
        "HmB520k":    "JX005047"
    },
    "1385": {
        "rabkinesin": "JX004161",
        "HmB00012":   "JX005692",
        "HmB449k":    "JX005475",
        "HmB453k":    "JX005256",
        "optix":      "JX004379",
        "HmB520k":    "JX005037"
    },
    "1023": {
        "rabkinesin": "JX004137",
        "HmB00012":   "JX005668",
        "HmB449k":    "JX005451",
        "HmB453k":    "JX005232",
        "optix":      "JX004355",
        "HmB520k":    "JX005013"
    },
    "9227": {
        "rabkinesin": "JX004097",
        "HmB00012":   "JX005628",
        "HmB449k":    "JX005411",
        "HmB453k":    "JX005192",
        "optix":      "JX004315",
        "HmB520k":    "JX004973"
    },
    "9182": {
        "rabkinesin": "JX004089",
        "HmB00012":   "JX005620",
        "HmB449k":    "JX005403",
        "HmB453k":    "JX005184",
        "optix":      "JX004307",
        "HmB520k":    "JX004965"
    },
    "9185": {
        "rabkinesin": "JX004077",
        "HmB00012":   "JX005608",
        "HmB449k":    "JX005391",
        "HmB453k":    "JX005172",
        "optix":      "JX004295",
        "HmB520k":    "JX004953"
    },
    "1911": {
        "rabkinesin": "JX004035",
        "HmB00012":   "JX005566",
        "HmB449k":    "JX005349",
        "HmB453k":    "JX005130",
        "optix":      "JX004254",
        "HmB520k":    "JX004911"
    },
    "2424": {
        "rabkinesin": "JX004011",
        "HmB00012":   "JX005542",
        "HmB449k":    "JX005325",
        "HmB453k":    "JX005106",
        "optix":      "JX004230",
        "HmB520k":    "JX004887"
    },
    "2253": {
        "rabkinesin": "JX004003",
        "HmB00012":   "JX005534",
        "HmB449k":    "JX005317",
        "HmB453k":    "JX005098",
        "optix":      "JX004222",
        "HmB520k":    "JX004879"
    },
    "2182": {
        "rabkinesin": "JX004001",
        "HmB00012":   "JX005532",
        "HmB449k":    "JX005315",
        "HmB453k":    "JX005096",
        "optix":      "JX004220",
        "HmB520k":    "JX004877"
    },
    "2387": {
        "rabkinesin": "JX003980",
        "HmB00012":   "JX005511",
        "HmB449k":    "JX005294",
        "HmB453k":    "JX005075",
        "optix":      "JX004199",
        "HmB520k":    "JX004856"
    }
}

def fetch_all_fasta():
    all_accs = []
    for v_id, l_map in ACCESSIONS.items():
        for locus, acc in l_map.items():
            all_accs.append(acc)
    
    print(f"Total unique accessions to fetch: {len(set(all_accs))}")
    acc_list_str = ",".join(all_accs)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc_list_str}&rettype=fasta&retmode=text"
    
    req = urllib.request.Request(url, headers={"User-Agent": "RhizAeon-Benchmarking/1.0 (academic research)"})
    with urllib.request.urlopen(req) as resp:
        fasta_text = resp.read().decode()
    
    seqs_by_acc = {}
    for rec in SeqIO.parse(StringIO(fasta_text), "fasta"):
        # The ID in FASTA is like JX004173.1
        acc_base = rec.id.split(".")[0]
        seqs_by_acc[acc_base] = str(rec.seq).upper()
    
    print(f"Successfully retrieved {len(seqs_by_acc)} sequences from NCBI GenBank.")
    return seqs_by_acc

def align_locus(locus_name, seqs_dict):
    """Align a single locus with MAFFT."""
    unaligned_file = os.path.join(BASE_DIR, f"temp_{locus_name}_unaligned.fasta")
    aligned_file = os.path.join(BASE_DIR, f"temp_{locus_name}_aligned.fasta")
    
    with open(unaligned_file, "w") as f:
        for label, v_id, sp, role in VOUCHERS:
            seq = seqs_dict[label]
            f.write(f">{label}\n{seq}\n")
    
    cmd = ["/usr/local/bin/mafft", "--auto", "--quiet", unaligned_file]
    with open(aligned_file, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)
    
    aligned_seqs = {}
    for rec in SeqIO.parse(aligned_file, "fasta"):
        aligned_seqs[rec.id] = str(rec.seq)
    
    # Clean up temp files
    os.remove(unaligned_file)
    os.remove(aligned_file)
    
    return aligned_seqs

def main():
    print("================================================================================")
    print("  Fetching & Aligning Heliconius Mimicry Introgression Island (Chromosome 18)   ")
    print("================================================================================")
    
    seqs_by_acc = fetch_all_fasta()
    
    aligned_by_locus = {}
    partition_coords = {}
    current_pos = 1
    
    for locus in LOCUS_ORDER:
        print(f"Aligning locus: {locus}...")
        locus_seqs = {}
        for label, v_id, sp, role in VOUCHERS:
            acc = ACCESSIONS[v_id][locus]
            locus_seqs[label] = seqs_by_acc[acc]
        
        aligned = align_locus(locus, locus_seqs)
        aligned_by_locus[locus] = aligned
        
        locus_len = len(list(aligned.values())[0])
        start_pos = current_pos
        end_pos = current_pos + locus_len - 1
        partition_coords[locus] = (start_pos, end_pos, locus_len)
        current_pos = end_pos + 1
        print(f"  {locus:12s}: {locus_len} bp (Coordinates {start_pos} - {end_pos})")
    
    # Concatenate in chromosomal order
    concatenated = {label: "" for label, _, _, _ in VOUCHERS}
    for locus in LOCUS_ORDER:
        aligned = aligned_by_locus[locus]
        for label in concatenated:
            concatenated[label] += aligned[label]
    
    total_len = len(list(concatenated.values())[0])
    print(f"\nConcatenation complete: {len(concatenated)} taxa, {total_len} bp total alignment length.")
    
    # Write concatenated alignment
    out_fasta = os.path.join(BASE_DIR, "heliconius_optix_island_aligned.fasta")
    with open(out_fasta, "w") as f:
        for label, v_id, sp, role in VOUCHERS:
            f.write(f">{label}\n{concatenated[label]}\n")
    
    with open(out_fasta, "rb") as f:
        sha256 = hashlib.sha256(f.read()).hexdigest()
    
    print(f"Saved alignment to: {out_fasta}")
    print(f"SHA-256: {sha256}")
    
    # Save partition metadata
    meta = {
        "study": "Pardo-Diaz et al. (2012) PLoS Genetics 8(6): e1002752",
        "locus_order": LOCUS_ORDER,
        "partitions": partition_coords,
        "total_length": total_len,
        "taxa": VOUCHERS,
        "sha256": sha256
    }
    meta_path = os.path.join(BASE_DIR, "heliconius_partitions.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved partition metadata to: {meta_path}")

if __name__ == "__main__":
    main()
