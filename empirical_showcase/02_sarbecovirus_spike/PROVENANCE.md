# Empirical Showcase 2: Natural Sarbecovirus Spike Receptor-Binding Domain Chimerism

## 1. Scientific Overview & Primary Literature Citations

* **Biological System**: Natural Sarbecovirus Genus (Coronaviridae: Betacoronavirus: Sarbecovirus).
* **Genomic Locus**: Full-length Spike Glycoprotein Coding Sequence ($S$ gene, $L = 3{,}843\text{ bp}$).
* **Evolutionary Question**: What is the natural recombinant mosaic origin of the SARS-CoV-2 Spike protein? Resolving the structural boundary between the N-terminal Domain (NTD) ancestral bat backbone and the Receptor-Binding Domain (RBD) introgression sharing human ACE2 receptor contact homology.

### Primary Publications
1. **Zhou P, Yang XL, Wang XG, et al.** (2020). A pneumonia outbreak associated with a new coronavirus of probable bat origin. *Nature*, 579(7798): 270–273. [doi:10.1038/s41586-020-2012-7](https://doi.org/10.1038/s41586-020-2012-7). PMID: 32015507.
2. **Xiao K, Zhai J, Feng Y, et al.** (2020). Isolation of SARS-CoV-2-related coronavirus from Malayan pangolins. *Nature*, 583(7815): 286–289. [doi:10.1038/s41586-020-2313-x](https://doi.org/10.1038/s41586-020-2313-x). PMID: 32380510.
3. **Lam TT, Jia N, Zhang YW, et al.** (2020). Identifying SARS-CoV-2-related coronaviruses in Malayan pangolins. *Nature*, 583(7815): 282–285. [doi:10.1038/s41586-020-2169-0](https://doi.org/10.1038/s41586-020-2169-0). PMID: 32218527.
4. **Temmam S, Vongphrachanh K, Baquero E, et al.** (2022). Bat coronaviruses related to SARS-CoV-2 and capable of using human ACE2 to enter cells. *Nature*, 604(7905): 330–336. [doi:10.1038/s41586-022-04532-4](https://doi.org/10.1038/s41586-022-04532-4). PMID: 35172323.
5. **Boni MF, Lemey P, Jiang X, et al.** (2020). Evolutionary origins of the SARS-CoV-2 sarbecovirus lineage responsible for the COVID-19 pandemic. *Nature Microbiology*, 5(11): 1408–1417. [doi:10.1038/s41564-020-0771-4](https://doi.org/10.1038/s41564-020-0771-4). PMID: 32724171.

---

## 2. Input Alignment Metadata & Accessions

All sequences were retrieved directly from NCBI GenBank via Entrez E-utilities:

| Taxon Label | GenBank Accession | Host Species | Geographic Origin | Clinical / Evolutionary Role | Spike CDS Length |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **SARS_CoV_2_Wuhan_Hu_1** | [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) | *Homo sapiens* | Wuhan, China (2019) | Pandemic index reference | 3,822 bp (1,274 aa) |
| **Bat_RaTG13** | [`MN996532.2`](https://www.ncbi.nlm.nih.gov/nuccore/MN996532.2) | *Rhinolophus affinis* | Yunnan, China (2013) | 96% genome-wide sister | 3,810 bp (1,270 aa) |
| **Bat_BANAL_52** | [`MZ937000.1`](https://www.ncbi.nlm.nih.gov/nuccore/MZ937000.1) | *Rhinolophus malayanus*| Feuang, Laos (2020) | ACE2-binding RBM homolog | 3,810 bp (1,270 aa) |
| **Bat_BANAL_103** | [`MZ937001.1`](https://www.ncbi.nlm.nih.gov/nuccore/MZ937001.1) | *Rhinolophus pusillus* | Fueng, Laos (2020) | ACE2-binding RBM homolog | 3,798 bp (1,266 aa) |
| **Pangolin_MP789** | [`MT121216.1`](https://www.ncbi.nlm.nih.gov/nuccore/MT121216.1) | *Manis javanica* | Guangdong, China (2019)| ACE2-contact ridge donor | 3,798 bp (1,266 aa) |
| **Bat_RmYN02** | [`MT121215.1`](https://www.ncbi.nlm.nih.gov/nuccore/MT121215.1) | *Rhinolophus malayanus*| Mengla, China (2019) | Natural S1/S2 insertions | 3,822 bp (1,274 aa) |
| **SARS_CoV_1_Tor2** | [`NC_004718.3`](https://www.ncbi.nlm.nih.gov/nuccore/NC_004718.3) | *Homo sapiens* | Toronto, Canada (2003) | 2002-2003 pandemic reference | 3,768 bp (1,256 aa) |
| **Bat_RsSHC014** | [`KC881005.1`](https://www.ncbi.nlm.nih.gov/nuccore/KC881005.1) | *Rhinolophus sinicus* | Kunming, China (2011) | Bat SARS-like coronavirus | 3,771 bp (1,257 aa) |
| **Bat_HKU3_1** | [`DQ022305.2`](https://www.ncbi.nlm.nih.gov/nuccore/DQ022305.2) | *Rhinolophus sinicus* | Hong Kong, China (2005)| Clade 2 sarbecovirus outgroup| 3,729 bp (1,243 aa) |

* **Aligned File**: `sarbecovirus_spike_aligned.fasta`
* **Alignment Method**: MAFFT v7.485 (`--auto`)
* **SHA-256 Hash**: `5e2107346e6db6c95fe5f6a62bca14498cc9720ab5603731a8ac1cf8beacb2c4`
* **Alignment Length**: $3{,}843\text{ bp}$ ($1{,}153\text{ informative sites}$)

---

## 3. Fresh Execution & Verified Results

* **Execution Command**:
  ```bash
  /Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon \
    -i sarbecovirus_spike_aligned.fasta \
    -o sarbecovirus_spike_result.json \
    --html sarbecovirus_spike_dashboard.html \
    -t "Natural Sarbecovirus Spike RBD Chimerism" \
    -l 16 -a 0.05
  ```
* **Engine Execution Time**: **`15.84 ms`** (Wall-clock: $28.3\text{ ms}$).
* **Verified Recombination Events**: 1

### Event Details
* **Recombinant Child**: `SARS_CoV_2_Wuhan_Hu_1`
* **Donor Parent**: `Bat_RaTG13` (Yunnan bat virus)
* **Home Parent**: `Pangolin_MP789` (Guangdong pangolin virus)
* **Alignment Coordinates**: `1 – 939 bp` (Length: $939\text{ bp}$)
* **Protein Structural Boundary**: Codon **`313`** (spanning the N-terminal Domain / Receptor-Binding Domain boundary)
* **Statistical Support**:
  * Fisher's Exact Test: **$p = 4.07 \times 10^{-22}$**
  * Standardized Physical $Z$-Score: **$Z = 8.23$**
  * Informative Substitutions: 275 SNPs
  * Divergence in NTD tract ($1\text{--}939\text{ bp}$):
    * Divergence to `Bat_RaTG13`: **`5.56%`**
    * Divergence to `Pangolin_MP789`: **`34.11%`**
  * Divergence downstream ($940\text{--}3{,}843\text{ bp}$, spanning RBD/RBM):
    * Divergence to `Pangolin_MP789`: drops to **`9.75%`**!

---

## 4. Ground-Truth Concordance

* The detected breakpoint at nucleotide 939 (codon 313) falls precisely at the structural junction between the N-terminal Domain (NTD, residues 14–305) and the Receptor-Binding Domain (RBD, residues 319–541).
* In the NTD, SARS-CoV-2 shares high sequence identity with Bat RaTG13 ($94.4\%$). Immediately at codon 313, the trajectory crosses over into the Pangolin-CoV / BANAL-52 coordinate subspace, sharing identical ACE2-binding residues.
* **Execution Time**: **`15.84 milliseconds`** for the complete 3.84 kb alignment across the Sarbecovirus genus.
