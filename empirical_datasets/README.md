# RhizAeon Empirical Datasets

This directory contains the three authentic empirical datasets and reproduction pipelines used in the RhizAeon manuscript (*"Continuous Sequence Manifolds Disentangle Reticulate Evolution across Genomes"*):

1. **`01_sarbecovirus_origins_lytras2022/`** — **Sarbecovirus Subgenus Recombination and SARS-CoV-2 Origins**
   - **Scale:** $N = 83$ taxa, $L = 30{,}965$ bp whole-genome alignment (plus historical 78-taxon $L = 30{,}956$ bp alignment).
   - **Biological Focus:** Fine-scale reticulation across the sarbecovirus lineage, resolving the mosaic origins of SARS-CoV-2 (Wuhan-Hu-1) and bat progenitors (BANAL-52, RaTG13, RpYN06) across the Spike glycoprotein receptor-binding domain (RBD) and ACE2 contact ridge.
   - **Benchmark Comparison:** Direct head-to-head comparison against 21 whole-genome Genetic Algorithm for Recombination Detection (GARD) breakpoints from Lytras et al. (2022) *Genome Biology and Evolution*.

2. **`02_spneumoniae_pmen1_croucher2011/`** — **Whole-Chromosome Antimicrobial Resistance & Capsule Switching in *Streptococcus pneumoniae* PMEN1**
   - **Scale:** $N = 250$ clinical isolates, $L = 2{,}221{,}315$ bp whole-chromosome alignment.
   - **Biological Focus:** Macro-scale bacterial horizontal gene transfer, penicillin-binding protein transpeptidase domain replacements (*pbp2x*, *pbp1a*), trimethoprim resistance (*folA*), and the serotype 23F $\to$ 19A capsular polysaccharide (*cps*) vaccine-escape switch.
   - **Benchmark Comparison:** Direct concordance against 52 published recombination boundaries inferred by Gubbins (Croucher et al. 2011 *Science*, 2015 *Nucleic Acids Research*).

3. **`03_hominoidea_mhc_class1/`** — **Primate MHC Class I / HLA Inter-Locus Ectopic Gene Conversion**
   - **Scale:** $N = 379$ hominoid alleles (*Homo sapiens*, *Pan troglodytes*, *Pan paniscus*, *Gorilla gorilla*, *Pongo abelii*), $L = 1{,}128$ bp full-length coding sequence across Exons 1–8.
   - **Biological Focus:** Micro-scale ectopic (non-allelic) gene conversion tracts, resolving the canonical 73-bp C1 functional epitope transfer from *HLA-C\*01:02* into the *HLA-B\*15:01* backbone generating *HLA-B\*46:01*, alongside chimpanzee inter-locus chimeras (*Patr-AL\*01*).

---

## Directory Organization

```
empirical_datasets/
├── README.md                                # Master empirical index and guide
├── 01_sarbecovirus_origins_lytras2022/      # Sarbecovirus 30.9 kb whole genomes
│   ├── 83sarbeco_wgal.fas                   # 83-taxon alignment with Laos BANAL viruses
│   ├── 78sarbeco_wgal_fn.fas                # 78-taxon historical alignment
│   ├── sarbecovirus_spike_aligned.fasta     # Spike gene alignment (N=83, L=3,822 bp)
│   ├── Table_S3.csv                         # 21 GARD breakpoints from Lytras et al. (2022)
│   ├── sarbeco_gene_annotations.json        # Curated ORF coordinates (ORF1ab, S, E, M, N)
│   ├── run_sarbecovirus_reanalysis.py       # RhizAeon scan and GARD comparison harness
│   ├── generate_figure.py                   # Generates 5-panel Fig. 3 landscape
│   └── README.md                            # Detailed documentation and walkthrough
├── 02_spneumoniae_pmen1_croucher2011/       # S. pneumoniae 2.22 Mb whole chromosomes
│   ├── PMEN1.aln.gz                         # Gzipped 250-isolate whole-chromosome alignment
│   ├── Spn23f.gff                           # ATCC 700669 reference genome annotation
│   ├── gubbins_predictions.gff              # Published Gubbins recombination boundaries
│   ├── pmen1_metadata.csv                   # Clinical isolate metadata and collection dates
│   ├── scan_pmen1_chromosome.py             # Whole-chromosome continuous field scanner
│   ├── generate_figure.py                   # Generates 4-panel Fig. 4 AMR/ARG landscape
│   └── README.md                            # Detailed documentation and walkthrough
└── 03_hominoidea_mhc_class1/                # Primate MHC Class I 1.1 kb alleles
    ├── hominoidea_class1_aligned.fasta      # 379-allele aligned coding sequences
    ├── scan_hominoidea_manifold.py          # Continuous physical field deconvolution
    ├── generate_figure.py                   # Generates 4-panel Fig. 5 MHC/ARG landscape
    └── README.md                            # Detailed documentation and walkthrough
```

---

## Quick-Start Reproduction Guide

Each dataset directory is self-contained. You can reproduce the analysis and figures in two modes:

### 1. Generating Figures from Precomputed Continuous Fields (< 5 seconds)
To immediately reproduce the publication figures:
```bash
# Sarbecovirus (Figure 3)
cd empirical_datasets/01_sarbecovirus_origins_lytras2022
python3 generate_figure.py

# Streptococcus pneumoniae PMEN1 (Figure 4)
cd ../02_spneumoniae_pmen1_croucher2011
python3 generate_figure.py

# Hominoidea MHC Class I (Figure 5)
cd ../03_hominoidea_mhc_class1
python3 generate_figure.py
```

### 2. Rerunning Scans from Raw Alignments
To execute the complete physical field scan from raw FASTA sequences:
```bash
# Install the RhizAeon engine (PyPI wheel with native extension)
pip install rhizaeon

# Sarbecovirus scan (~1.5 seconds)
cd empirical_datasets/01_sarbecovirus_origins_lytras2022
python3 run_sarbecovirus_reanalysis.py

# S. pneumoniae PMEN1 chromosome scan (~45 seconds)
cd ../02_spneumoniae_pmen1_croucher2011
python3 scan_pmen1_chromosome.py

# Hominoidea MHC Class I scan (~2 seconds)
cd ../03_hominoidea_mhc_class1
python3 scan_hominoidea_manifold.py
```

---

## References

1. **Lytras S, Hughes J, Martin D, Swanepoel P, de Klerk A, Lourens R, Kosakovsky Pond SL, Xia W, Jiang X, Robertson DL** (2022). Exploring the Natural Origins of SARS-CoV-2 in the Light of Recombination. *Genome Biology and Evolution*, 14(2): evac018. [doi:10.1093/gbe/evac018](https://doi.org/10.1093/gbe/evac018).
2. **Temmam S, Vongphrachanh K, Baquero E, et al.** (2022). Bat coronaviruses related to SARS-CoV-2 and capable of using human ACE2 to enter cells. *Nature*, 604: 330–336. [doi:10.1038/s41586-022-04532-4](https://doi.org/10.1038/s41586-022-04532-4).
3. **Croucher NJ, Harris SR, Fraser C, et al.** (2011). Rapid pneumococcal evolution in response to clinical interventions. *Science*, 331(6016): 430–434. [doi:10.1126/science.1198545](https://doi.org/10.1126/science.1198545).
4. **Croucher NJ, Page AJ, Connor TR, et al.** (2015). Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Research*, 43(3): e15. [doi:10.1093/nar/gku1196](https://doi.org/10.1093/nar/gku1196).
5. **Robinson J, Barker DJ, Georgiou X, et al.** (2020). IPD-IMGT/HLA Database. *Nucleic Acids Research*, 48(D1): D948–D955. [doi:10.1093/nar/gkz950](https://doi.org/10.1093/nar/gkz950).
