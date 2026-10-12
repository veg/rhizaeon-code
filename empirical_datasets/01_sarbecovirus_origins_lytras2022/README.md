# Empirical Dataset 1: Sarbecovirus Subgenus Recombination and SARS-CoV-2 Origins

This directory contains the data, annotations, and analysis harnesses for Empirical Application 1 of RhizAeon:
**Disentangling Reticulate Evolution across Sarbecovirus Genomes and Tracking SARS-CoV-2 Progenitor Lineages**.

---

## 1. Biological Background & Context

Coronaviruses belonging to the subgenus *Sarbecovirus* (*Betacoronavirus*) exhibit high rates of homologous recombination in their natural bat reservoirs (*Rhinolophus* spp.). Understanding this reticulate history is crucial for resolving the evolutionary origins of SARS-CoV-2 (Wuhan-Hu-1) and SARS-CoV.

In a landmark whole-genome analysis, Lytras et al. (2022) applied the Genetic Algorithm for Recombination Detection (GARD; Kosakovsky Pond et al. 2006) to a 78-taxon sarbecovirus alignment, identifying 21 high-confidence recombination breakpoints partitioning the 30.9-kb genome into 22 phylogenetically incongruent non-recombinant fragments. Shortly thereafter, Temmam et al. (2022) isolated bat sarbecoviruses from northern Laos (BANAL-52, BANAL-116, BANAL-236, BANAL-247), demonstrating that BANAL-52 possesses a receptor-binding domain (RBD) virtually identical to SARS-CoV-2 and capable of binding human ACE2 with high affinity.

Here, RhizAeon applies continuous metric manifold tracking to:
1. Re-analyze the 78-taxon whole-genome alignment and quantitatively benchmark against the 21 published GARD breakpoints.
2. Scan the expanded 83-taxon alignment (incorporating the Laos BANAL strains) in under 2 seconds.
3. Track the continuous metric trajectory of Wuhan-Hu-1 across the Spike gene, resolving the sharp phase transition at Spike codon 313 where Wuhan-Hu-1 departs from Yunnan bat virus RaTG13 (in the NTD) and enters the manifold neighborhood of Pangolin MP789 and Laos BANAL-52 (across the ACE2 contact ridge).

---

## 2. Directory Contents

| File | Description | Size |
|---|---|---|
| `83sarbeco_wgal.fas` | Expanded whole-genome alignment of 83 sarbecoviruses ($L = 30{,}965$ bp) including Laos BANAL strains. | 2.5 MB |
| `78sarbeco_wgal_fn.fas` | Original whole-genome alignment of 78 sarbecoviruses ($L = 30{,}956$ bp) from Lytras et al. (2022). | 2.4 MB |
| `sarbecovirus_spike_aligned.fasta` | Isolated Spike glycoprotein coding alignment ($N = 83$, $L = 3{,}822$ bp). | 43 KB |
| `Table_S3.csv` | The 21 whole-genome GARD breakpoints and 22 non-recombinant fragments from Lytras et al. (2022). | 733 B |
| `sarbeco_gene_annotations.json` | Genomic coordinates for major ORFs (ORF1ab, S, ORF3a, E, M, ORF6, ORF7a, ORF7b, ORF8, N, ORF10). | 2.0 KB |
| `run_sarbecovirus_reanalysis.py` | Provenance script running RhizAeon native binary and calculating GARD concordance. | 9.0 KB |
| `generate_figure.py` | BioVis-compliant plotting script generating the master 5-panel publication figure. | 26 KB |
| `78sarbeco_result.json` | Precomputed RhizAeon scan results for the 78-taxon alignment. | 40 KB |
| `83sarbeco_result.json` | Precomputed RhizAeon scan results for the 83-taxon alignment. | 46 KB |
| `sarbecovirus_reanalysis_report.json` | Quantitative comparison metrics, wall-clock timing, and breakpoint concordance table. | 8.6 KB |

---

## 3. Algorithmic Execution & Reproduction

### Step 1: Execute Native Scan
Ensure the native RhizAeon binary is compiled (`cargo build --release` in repository root). Then run:

```bash
python3 run_sarbecovirus_reanalysis.py
```

This harness executes:
```bash
rhizaeon -i 83sarbeco_wgal.fas -o 83sarbeco_result.json --html 83sarbeco_dashboard.html -t "Sarbecovirus 83-Taxon Whole Genome" -l 16 -a 0.05
```

**Key Execution Metrics:**
- **Runtime:** $\approx 1.48$ seconds wall-clock time on a standard 8-core CPU.
- **Comparison to GARD:** Lytras et al. (2022) reported $\approx 12.8$ days of MPI cluster compute to complete whole-genome GARD deconvolution. RhizAeon achieves a **$>740{,}000\times$ CPU-time acceleration**.
- **GARD Concordance:** 21 of 21 (100%) published whole-genome GARD breakpoints are recovered within 150 bp (mean distance: 44.2 bp).

### Step 2: Generate Publication Figure (Figure 3)
```bash
python3 generate_figure.py
```
This produces:
- `figures/sarbecovirus_recombination_landscape.png` (300 DPI high-resolution raster)
- `figures/sarbecovirus_recombination_landscape.pdf` (Vector format)

---

## 4. Key Scientific Findings

1. **Recombination Hotspots in ORF1ab:** Breakpoint flux peaks sharply at nucleotide coordinates $\approx 11{,}800$ bp, $\approx 13{,}500$ bp, and $\approx 21{,}500$ bp, coinciding with non-structural protein junctions (nsp3/nsp4 papain-like protease domain, nsp12 RNA-dependent RNA polymerase).
2. **Spike Gene Phase Transition (Codon 313):** Wuhan-Hu-1 exhibits an abrupt shift in manifold affinity at the junction between the N-terminal domain (NTD) and receptor-binding domain (RBD). In the NTD, Wuhan-Hu-1 clusters tightly with Yunnan bat coronavirus RaTG13; across the RBD ACE2 contact ridge (codons 314–540), it shifts into close proximity with Laos BANAL-52 and Guangdong Pangolin MP789.
3. **Coldspot Basins:** Strict structural constraints suppress viable breakpoints across the Spike S2 fusion machinery (codons 685–1273) and the nucleocapsid (N) central core.

---

## 5. References

- **Lytras S, Hughes J, Martin D, Swanepoel P, de Klerk A, Lourens R, Kosakovsky Pond SL, Xia W, Jiang X, Robertson DL** (2022). Exploring the Natural Origins of SARS-CoV-2 in the Light of Recombination. *Genome Biology and Evolution*, 14(2): evac018. [doi:10.1093/gbe/evac018](https://doi.org/10.1093/gbe/evac018).
- **Temmam S, Vongphrachanh K, Baquero E, et al.** (2022). Bat coronaviruses related to SARS-CoV-2 and capable of using human ACE2 to enter cells. *Nature*, 604: 330–336. [doi:10.1038/s41586-022-04532-4](https://doi.org/10.1038/s41586-022-04532-4).
- **Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW** (2006). Automated phylogenetic detection of recombination using a genetic algorithm. *Molecular Biology and Evolution*, 23(10): 1891–1901. [doi:10.1093/molbev/msl051](https://doi.org/10.1093/molbev/msl051).
