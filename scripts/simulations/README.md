# RhizAeon Continuous Simulation Benchmark Harness

This directory contains the complete simulation execution scripts, cluster submission manifests, and aggregation utilities to generate and reproduce the **1,000,000-genome continuous parameter sweep** presented in the RhizAeon manuscript (*"Continuous Sequence Manifolds Disentangle Reticulate Evolution across Genomes"*, Figure 2, Table 1).

> [!NOTE]
> **Data Policy & Reproducibility:**
> In accordance with the repository data policy, large simulation result files (`grand_sweep_1m.csv`, $\approx 85$ MB) are excluded from version control. Instead, all simulation generation code (implemented in native Rust in the HyphAeon monorepo at `rust/rhizaeon/crates/rhizaeon-tools/src/bin/grand_sweep.rs`), Slurm job array scripts, merger harnesses, and analysis pipelines are included for complete, deterministic reproduction.

---

## 1. Simulation Engine Architecture

The simulation engine is implemented as a high-performance multithreaded native Rust binary (`grand_sweep`). For every trial, it synthesizes an evolutionary tree, samples evolutionary parameters from broad continuous priors, simulates sequence evolution under the HKY85 nucleotide substitution model with discrete Gamma site rate heterogeneity, introduces topological recombination events, executes the full RhizAeon PA-FF detection pipeline, and evaluates bipartite event matching against ground truth:

### Evolutionary Parameter Distributions (Continuous Priors)
- **Taxon Count ($N$):** Sampled continuously from $N \in [4, 250]$ taxa.
- **Alignment Length ($L$):** Sampled continuously from $L \in [300, 30{,}000]$ bp.
- **Pairwise Divergence ($d$):** Mean tree divergence $d \in [0.005, 0.35]$ substitutions/site.
- **Transition/Transversion Ratio ($\kappa$):** Sampled from $\kappa \in [1.5, 6.0]$.
- **Site Rate Heterogeneity ($\alpha$):** 4-category discrete Gamma distribution with shape parameter $\alpha \in [0.2, 5.0]$.
- **Base Frequencies ($\pi$):** Non-uniform nucleotide frequencies sampled from a Dirichlet prior.

### Event Cardinality & Topological Diversity
Each trial is assigned a true event cardinality $K_{\mathrm{true}} \in \{0, 1, 2, 3, 4\}$ according to:
- **$K = 0$ (33.3% True Clonal Nulls):** Zero recombination events. Used to measure the exact empirical False Positive Rate (FPR) under extreme heterotachy and substitution rate variation.
- **$K = 1$ (36.7% Single Events):** Single recombination event per genome:
  - 5' Terminal Crossovers: $u_1 = 1$, $u_2 \sim \text{Uniform}(80, L - 80)$.
  - 3' Terminal Crossovers: $u_1 \sim \text{Uniform}(80, L - 80)$, $u_2 = L$.
  - Intercalary Conversion Cassettes: $1 < u_1 < u_2 < L$ with tract length spanning 45 bp to $0.45 \times L$.
- **$K = 2$ (16.7% Dual Events):** Two distinct recombination events (non-overlapping or staggered).
- **$K = 3$ (8.0% Triple Events):** Three distinct events (multi-way mosaicism).
- **$K = 4$ (5.3% Quad Events):** Four distinct events, including multi-donor introgression where the recombinant child inherits genetic material from multiple distinct parental clades ($P_1 \leftarrow P_2 + P_3$).

### Algorithmic Evaluation & Bipartite Scoring
In addition to binary sequence-level classification (`detected`, `is_recombinant`, `home_match`, `donor_match`), each trial is evaluated using **greedy bipartite event matching**:
- An inferred event $(\hat{u}_1, \hat{u}_2)$ matches a true event $(u_1, u_2)$ if both $|\hat{u}_1 - u_1| \le 150$ bp and $|\hat{u}_2 - u_2| \le 150$ bp.
- **`event_recall`:** $\frac{\text{Matched True Events}}{K_{\mathrm{true}}}$
- **`event_precision`:** $\frac{\text{Matched Inferred Events}}{\hat{K}}$
- **`event_f1`:** Harmonic mean of event recall and precision.
- **`mean_jaccard`:** Spatial intersection-over-union $\frac{[u_1, u_2] \cap [\hat{u}_1, \hat{u}_2]}{[u_1, u_2] \cup [\hat{u}_1, \hat{u}_2]}$ across matched events.
- **`mean_bp_error`:** Average boundary localization error in base pairs.
- **`cardinality_error`:** Discrepancy between inferred and true event counts $|\hat{K} - K_{\mathrm{true}}|$.

---

## 2. Directory Contents

| File | Description |
|---|---|
| `run_1m_sweep.slurm` | Slurm job array submission script (50 tasks $\times$ 20,000 trials across 512 cores). |
| `merge_sweep_1m.py` | Python merger script that concatenates, deduplicates, and sorts partition CSVs. |
| `run_local_sweep.py` | Python wrapper to run custom-sized simulation sweeps on a local workstation. |
| `README.md` | This technical guide and reproduction manual. |

---

## 3. Step-by-Step Reproduction Guide

### Prerequisites: Build the Native Simulation Binary
From the repository root:
```bash
cargo build --release --bin grand_sweep
```
The compiled binary will be placed at `target/release/grand_sweep`.

### Option A: Running a Local Sweep (e.g., 1,000–10,000 Trials)
To run a local sweep on your multi-core workstation:
```bash
python3 scripts/simulations/run_local_sweep.py \
    --trials 1000 \
    --threads 8 \
    --seed 42 \
    --output local_sweep_1k.csv
```

### Option B: Executing the Full 1,000,000 Sweep on a Slurm Cluster
The 1,000,000-trial benchmark is configured as a 50-task Slurm array (tasks 0–49, 20,000 trials each, 32 threads and 48 GB RAM per task):

```bash
# Submit the array to Slurm
sbatch scripts/simulations/run_1m_sweep.slurm
```

Each array task outputs its partition to `results/sweep_part_${SLURM_ARRAY_TASK_ID}.csv`.

### Step 2: Merge Distributed Partitions
Once all 50 tasks complete:
```bash
python3 scripts/simulations/merge_sweep_1m.py \
    --results-dir results \
    --output grand_sweep_1m.csv \
    --expected 1000000
```
This utility:
1. Inspects all 50 chunk files.
2. Concatenates and deduplicates on `trial_id`.
3. Verifies that all trial IDs from 0 to 999,999 are present with zero missing entries.
4. Writes the final unified `grand_sweep_1m.csv`.

### Step 3: Analyze and Generate Manuscript Figures
From the repository root, run the comprehensive analysis harness:
```bash
python3 analyze_grand_sweep.py grand_sweep_1m.csv
```

This generates:
1. **`grand_sweep_1m_phase_diagram.png` (Figure 2):** 4-panel operational flight envelope:
   - Panel A: Empirical Power vs Alignment Length $L$ and Sequence Divergence $d$.
   - Panel B: False Positive Rate (FPR) Calibration on True Nulls ($K=0$) across divergence and taxon scale.
   - Panel C: Breakpoint Localization Error ($|\hat{u}_1 - u_1|$) as a function of conversion tract length.
   - Panel D: Computational Throughput (genomes processed per second) as a function of $N \times L$.
2. **`grand_sweep_1m_multievent_analysis.png` (Supplementary Figure S1):** 4-panel multi-event diagnostic figure:
   - Panel A: Event Cardinality Confusion Matrix ($K_{\mathrm{true}}$ vs $\hat{K}$).
   - Panel B: Bipartite Event Metrics (Recall, Precision, F1, Spatial Jaccard) stratified by $K_{\mathrm{true}}$.
   - Panel C: Detection Sensitivity across Topological Regimes (5' Terminal, 3' Terminal, Intercalary Cassettes).
   - Panel D: Multi-Donor Lineage Disambiguation ($P_1 \leftarrow P_2 + P_3$).
3. **`grand_sweep_1m_summary.txt`:** Quantitative summary table with exact metrics, 95% confidence intervals, and runtime statistics.

---

## 4. Key Benchmark Results (1,000,000 Trials)

- **Strict FPR Control on Nulls ($K = 0$):** **$\text{FPR} < 0.05$** at nominal $\alpha = 0.05$ across all divergence and length strata, proving robustness against extreme heterotachy and APOBEC-like hypermutation.
- **High Sensitivity on Recombinants ($K \ge 1$):** Overall detection power exceeds **$88\%$** for tracts $\ge 150$ bp and divergence $d \ge 0.02$.
- **Sub-15 nt Breakpoint Precision:** Median boundary localization error $< 12$ bp on informative conversion tracts.
- **Linear-Time Throughput:** Sustained processing rate exceeding **$350\text{–}400$ genomes/second**, scaling strictly linearly with matrix cells $O(N \times L)$.
