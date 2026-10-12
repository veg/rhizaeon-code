# Empirical Dataset 2: Whole-Chromosome Antimicrobial Resistance & Capsule Switching in *Streptococcus pneumoniae* PMEN1

This directory contains the data, annotations, and analysis harnesses for Empirical Application 2 of RhizAeon:
**Whole-Chromosome Continuous Field Deconvolution across the Pandemic Multidrug-Resistant PMEN1 Lineage**.

---

## 1. Biological Background & Context

*Streptococcus pneumoniae* is a major human pathogen responsible for pneumonia, meningitis, and septicemia. The PMEN1 clone (Spain23F-1, sequence type ST81) is a globally disseminated multidrug-resistant pandemic lineage that emerged in Europe in the 1980s. Its clinical evolutionary trajectory is driven by extensive horizontal gene transfer (HGT) via natural genetic transformation:

1. **Antimicrobial Resistance (AMR):** High-level resistance to beta-lactam antibiotics (penicillin) was acquired through inter-species homologous recombination of mosaic *pbp2x* and *pbp1a* alleles (penicillin-binding protein transpeptidase domains) from commensal viridans streptococci (*S. mitis*, *S. oralis*). Resistance to trimethoprim (*folA*) and macrolides/tetracyclines (*tetM*, *ermB*) was similarly integrated via mobile cassettes.
2. **Capsule Switching (PCV7 Vaccine Escape):** Following the introduction of the 7-valent pneumococcal conjugate vaccine (PCV7) targeting the ancestral 23F capsular serotype, PMEN1 isolates underwent a massive 15.5-kb recombination replacement of the capsular polysaccharide synthesis (*cps*) locus, switching from serotype 23F to serotype 19A and completely evading vaccine-induced immunity (Croucher et al. 2011 *Science*).

In this application, RhizAeon scales to whole bacterial chromosomes ($L = 2{,}221{,}315$ bp across $N = 250$ clinical isolates), extracting continuous breakpoint flux $\Phi_{\mathrm{bp}}(u)$ and localized kinetic strain $\rho(u)$ without partitioning windows or heuristics.

---

## 2. Directory Contents

| File | Description | Size |
|---|---|---|
| `PMEN1.aln.gz` | Gzipped whole-chromosome multiple sequence alignment of 250 PMEN1 clinical isolates ($L = 2{,}221{,}315$ bp) mapped to ATCC 700669. | 6.4 MB |
| `Spn23f.gff` | Standard GFF3 gene annotations for the *S. pneumoniae* ATCC 700669 reference chromosome (GenBank FM211187). | 560 KB |
| `gubbins_predictions.gff` | The 52 published recombination event boundaries inferred by Gubbins (Croucher et al. 2011, 2015). | 25 KB |
| `pmen1_metadata.csv` | Clinical isolate metadata: isolation year, geographic origin, MLST sequence type, serotype, and AMR phenotypes. | 18 KB |
| `scan_pmen1_chromosome.py` | Provenance script running the whole-chromosome continuous physical field deconvolution pipeline. | 15 KB |
| `generate_figure.py` | BioVis-compliant plotting script generating the master 4-panel publication figure (Figure 4). | 30 KB |
| `pmen1_continuous_fields.npz` | Precomputed continuous breakpoint flux $\Phi(u)$ and strain $\rho(u)$ across the 2.22-Mb chromosome. | 39 KB |
| `pmen1_fields_summary.json` | JSON summary containing detected changepoints, kinetic Z-scores, and runtime statistics. | 13 KB |
| `pmen1_changepoints_vs_gubbins.csv` | Head-to-head boundary concordance table comparing RhizAeon vs Gubbins. | 9.2 KB |
| `rhizaeon_pmen1_summary.json` | High-level summary of significant reticulation events and detected clinical landmarks. | 2.5 KB |

---

## 3. Algorithmic Execution & Reproduction

### Step 1: Execute Whole-Chromosome Scan
Run the end-to-end continuous field deconvolution script:
```bash
python3 scan_pmen1_chromosome.py
```

**Computational Highlights:**
- **Segregating SNP Extraction:** Out of 2,221,315 genomic positions, RhizAeon identifies 38,162 segregating informative SNPs ($58.2\times$ matrix compression) in **0.42 seconds**.
- **Whole-Chromosome Deconvolution:** Constructs the 3-axis Gower metric space, tracks the $250 \times 2{,}221{,}315$ continuous trajectory tensor, and computes spatial force fields across the entire 2.22-Mb chromosome in **42.6 seconds**.
- **Gubbins Benchmarking:** Gubbins required $\approx 14.5$ hours on the same multicore host, demonstrating a **$>1{,}200\times$ speedup** for RhizAeon.

**Head-to-Head Boundary Concordance with Gubbins (52 Boundaries):**
- Median distance to nearest RhizAeon changepoint: **342 bp**
- Concordance $\le 1\text{ kb}$: **88.5%** (46/52)
- Concordance $\le 2\text{ kb}$: **96.2%** (50/52)
- Concordance $\le 5\text{ kb}$: **100.0%** (52/52)

### Step 2: Generate Publication Figure (Figure 4)
```bash
python3 generate_figure.py
```
This produces:
- `figures/fig_spneumoniae_chromosome_landscape.png` (300 DPI high-resolution raster)
- `figures/fig_spneumoniae_chromosome_landscape.pdf` (Vector format)

---

## 4. Key Clinical Landmarks Identified

RhizAeon precisely localizes all major historical recombination events documented in PMEN1 evolution:

| Landmark | Genomic Span (bp) | Kinetic Z-Score | Biological Impact |
|---|---|---|---|
| *pbp2x* (Penicillin resistance) | 285,069 – 305,751 | **$Z = 14.40$** ($P = 9.85 \times 10^{-35}$) | Focal 1,121-bp replacement of the transpeptidase catalytic domain from commensal viridans streptococci. |
| *cps* (Capsule 23F $\to$ 19A switch) | 336,730 – 352,271 | **$Z = 3.19$** | 15.5-kb horizontal transfer conferring PCV7 vaccine escape in clinical isolate 3355. |
| *pbp1a* (High-level beta-lactam res.) | 323,721 – 325,108 | **$Z = 11.23$** | Transpeptidase domain remodeling synergistic with *pbp2x*. |
| *folA* (Trimethoprim resistance) | 1,620,528 – 1,649,446 | **$Z = 9.87$** | Dihydrofolate reductase replacement cassette. |
| *pspC* / *cbpA* (Adhesin diversity) | 2,104,999 – 2,113,473 | **$Z = 7.64$** | Surface choline-binding protein domain shuffling. |
| Hotspot 73 kb Mega-Cassette | 73,625 – 76,580 | **$Z = 8.12$** | Large chromosomal integration island. |

---

## 5. References

- **Croucher NJ, Harris SR, Fraser C, et al.** (2011). Rapid pneumococcal evolution in response to clinical interventions. *Science*, 331(6016): 430–434. [doi:10.1126/science.1198545](https://doi.org/10.1126/science.1198545).
- **Croucher NJ, Page AJ, Connor TR, et al.** (2015). Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Research*, 43(3): e15. [doi:10.1093/nar/gku1196](https://doi.org/10.1093/nar/gku1196).
