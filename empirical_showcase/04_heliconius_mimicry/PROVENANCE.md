# Empirical Showcase 4: Heliconius Mimicry Adaptive Introgression Island

## 1. Scientific Overview & Primary Literature Citations

* **Biological System**: *Heliconius* Butterflies (*Nymphalidae: Heliconiinae: Heliconiini*).
* **Genomic Locus**: Chromosome 18 *HmB* / *optix* Adaptive Introgression Island ($L = 3{,}677\text{ bp}$, physical span $\approx 170\text{ kb}$ across 6 linked loci).
* **Evolutionary Question**: Do complex adaptive phenotypes arise by convergent de novo mutation or via interspecific introgressive hybridization? Resolving the adaptive introgression of the red wing-patterning switch between *Heliconius melpomene* and *Heliconius timareta*, where reticulate gene flow transfers a functional cis-regulatory cassette across strict species boundaries.

### Primary Publications
1. **Pardo-Diaz C, Salazar C, Baxter SW, Merot C, Linares M, Jiggins CD, McMillan WO** (2012). Adaptive introgression across species boundaries in *Heliconius* butterflies. *PLoS Genetics*, 8(6): e1002752. [doi:10.1371/journal.pgen.1002752](https://doi.org/10.1371/journal.pgen.1002752). PMID: 22737081.
2. **The Heliconius Genome Consortium** (2012). Butterfly genome reveals promiscuous exchange of mimicry adaptations among species. *Nature*, 487(7405): 94–98. [doi:10.1038/nature11041](https://doi.org/10.1038/nature11041). PMID: 22722851.
3. **Reed RD, Papa R, Martin A, et al.** (2011). *optix* drives the repeated convergent evolution of butterfly wing pattern mimicry. *Science*, 333(6046): 1137–1141. [doi:10.1126/science.1208227](https://doi.org/10.1126/science.1208227). PMID: 21852495.

---

## 2. Input Alignment Metadata & Accessions

The alignment concatenates 6 Sanger-sequenced loci along Chromosome 18 in physical order across 12 authentic vouchers from Pardo-Diaz et al. (2012) retrieved directly from NCBI GenBank (`JX003980–JX005837`):

### 2.1 Taxon Vouchers
| Taxon Label | Voucher ID | Species / Subspecies | Geographic Origin | Phenotype / Evolutionary Role |
| :--- | :--- | :--- | :--- | :--- |
| **`H_tim_tim_9185`** | `9227` | *H. timareta timareta* | San Martin, Peru | Peruvian rayed mimic (**introgressed recombinant**) |
| **`H_tim_tim_9182`** | `9182` | *H. timareta timareta* | San Martin, Peru | Peruvian rayed mimic (**introgressed recombinant**) |
| **`H_tim_tim_9227`** | `9227` | *H. timareta timareta* | San Martin, Peru | Peruvian rayed mimic (**introgressed recombinant**) |
| **`H_mel_mal_2312`** | `2312` | *H. melpomene malleti* | Pastaza, Ecuador | Amazonian rayed co-mimic (**donor parent**) |
| **`H_mel_mal_1261`** | `1261` | *H. melpomene malleti* | Pastaza, Ecuador | Amazonian rayed co-mimic (**donor parent**) |
| **`H_mel_agl_1385`** | `1385` | *H. melpomene aglaope* | San Martin, Peru | Amazonian rayed co-mimic (**donor parent**) |
| **`H_mel_ros_1023`** | `1023` | *H. melpomene rosina* | Gamboa, Panama | Central American postman morph |
| **`H_cyd_cor_2182`** | `2182` | *H. cydno cordula* | Tachira, Venezuela | White/yellow sister (**home recipient clade**) |
| **`H_cyd_cor_2253`** | `2253` | *H. cydno cordula* | Tachira, Venezuela | White/yellow sister (**home recipient clade**) |
| **`H_cyd_chi_2424`** | `2424` | *H. cydno chioneus* | Colon, Panama | White/yellow sister (**home recipient clade**) |
| **`H_cyd_wey_1911`** | `1911` | *H. cydno weymeri* | Valle del Cauca, Colombia | Polymorphic sister (**home recipient clade**) |
| **`H_num_2387`** | `2387` | *H. numata* | San Martin, Peru | Silvaniform outgroup |

### 2.2 Physical Locus Architecture on Chromosome 18
| Locus Name | Chromosome 18 Coordinate | Alignment Coordinates | Partition Span | Biological Role |
| :--- | :--- | :--- | :--- | :--- |
| **`rabkinesin`** | $\approx 350\text{ kb}$ | `1 – 750 bp` | $750\text{ bp}$ | Upstream flanking control locus |
| **`HmB00012`** | $\approx 400\text{ kb}$ | `751 – 987 bp` | $237\text{ bp}$ | Upstream transition marker |
| **`HmB449k`** | $449\text{ kb}$ | `988 – 1,763 bp` | $776\text{ bp}$ | Introgression island 5' boundary marker |
| **`HmB453k`** | $453\text{ kb}$ | `1,764 – 2,650 bp` | $887\text{ bp}$ | **Peak phenotype-associated cis-regulatory switch** |
| **`optix`** | $515\text{ kb}$ | `2,651 – 3,381 bp` | $731\text{ bp}$ | Master red-pattern transcription factor |
| **`HmB520k`** | $520\text{ kb}$ | `3,382 – 3,677 bp` | $296\text{ bp}$ | Downstream flanking control marker |

* **Aligned File**: `heliconius_optix_island_aligned.fasta`
* **Alignment Method**: MAFFT v7.485 (`--auto`)
* **SHA-256 Hash**: `1b6695c9341602ccec83f41a2c3fbfc4096edf880eadd71513e396ecdf5e4022`
* **Alignment Dimensions**: 12 taxa, $3{,}677\text{ bp}$ ($126\text{ informative sites}$)

---

## 3. Fresh Execution & Verified Results

* **Execution Command**:
  ```bash
  /Users/sergei/Projects/rhizaeon-code/target/release/rhizaeon \
    -i heliconius_optix_island_aligned.fasta \
    -o heliconius_optix_result.json \
    --html heliconius_optix_dashboard.html \
    -t "Heliconius Mimicry Adaptive Introgression Island" \
    -l 16 -a 0.05
  ```
* **Engine Execution Time**: **`52.12 ms`** (Wall-clock: $77.5\text{ ms}$).
* **Verified Recombination Events**: 2

### Primary Event Details: Interspecific Adaptive Introgression
* **Recombinant Child**: `H_tim_tim_9185` (*Heliconius timareta timareta*)
* **Donor Parent**: `H_mel_mal_2312` (*Heliconius melpomene malleti*)
* **Home Parent**: `H_mel_agl_1385`
* **Event Topology**: **Cassette Swap / Introgression Island** (`is_crossover: false`)
* **Alignment Coordinates**: **`2,170 – 2,650 bp`** (Tract length: $481\text{ bp}$)
* **Statistical Support**:
  * Fisher's Exact Test: **$p = 6.00 \times 10^{-7}$**
  * Standardized Physical $Z$-Score: **$Z = 26.68$**
  * Informative Substitutions: 18 SNPs
* **Sequence Divergence in Tract**:
  * Divergence to donor `H_mel_mal_2312`: **`0.0000` (0.00% — IDENTICAL!)**
  * Divergence to home background: **`12.86%`**

---

## 4. Ground-Truth Biological Concordance

1. **Exact Localization to Peak Regulatory Switch (`HmB453k`)**:
   * The detected introgression cassette spans coordinates `2,170 – 2,650 bp`, which terminates at the exact physical boundary of the `HmB453k` partition (`1,764 – 2,650 bp`).
   * In Pardo-Diaz et al. (2012, Table 1), `HmB453k` was identified as the single marker across the entire *Heliconius* genome with the highest percentage of variation explained by colour pattern ($62.42\%$), rather than species phylogeny.
2. **Bitwise Zero Divergence to Amazonian Co-Mimic Donor**:
   * Within the $481\text{ bp}$ tract, `H_tim_tim_9185` is **100% identical** ($d = 0.000$) to the Amazonian co-mimic `H_mel_mal_2312`, confirming direct, recent adaptive introgression between species.
3. **Execution Latency**:
   * RhizAeon parsed and scanned the 12-taxon, 3.68-kb multi-locus island across 6 physical markers in **`52.12 milliseconds`**.
