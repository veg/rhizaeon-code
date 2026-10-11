# Posada & Crandall (2001) Sample Replicates Dossier

This directory contains 10 representative replicate alignments (`.fas`) saved during the full Posada & Crandall (2001) benchmark simulation suite for review.

## Summary of Saved Replicates

| # | Alignment File | Category | $\theta$ | $\rho$ | $\alpha$ | Partitions | Ground Truth Events | RhizAeon Detection |
|---|----------------|----------|-----------|---------|-----------|------------|---------------------|--------------------|
| 1 | [`fp_theta100_alpha0p5_rep00.fas`](fp_theta100_alpha0p5_rep00.fas) | `false_positive` | 100.0 | 0.0 | 0.5 | 1 | 0 events | Certified Null (0 events) |
| 2 | [`fp_theta50_alpha0p05_rep00.fas`](fp_theta50_alpha0p05_rep00.fas) | `false_positive` | 50.0 | 0.0 | 0.05 | 1 | 0 events | Certified Null (0 events) |
| 3 | [`power_theta10_rho0_rep00.fas`](power_theta10_rho0_rep00.fas) | `power` | 10.0 | 0.0 | None | 1 | 0 events | Certified Null (0 events) |
| 4 | [`power_theta10_rho16_rep00.fas`](power_theta10_rho16_rep00.fas) | `power` | 10.0 | 16.0 | None | 51 | 12 events | Certified Null (0 events) |
| 5 | [`power_theta200_rho16_rep00.fas`](power_theta200_rho16_rep00.fas) | `power` | 200.0 | 16.0 | None | 37 | 19 events | **2 Event(s)** |
| 6 | [`power_theta50_rho0_rep00.fas`](power_theta50_rho0_rep00.fas) | `power` | 50.0 | 0.0 | None | 1 | 0 events | Certified Null (0 events) |
| 7 | [`power_theta50_rho16_rep00.fas`](power_theta50_rho16_rep00.fas) | `power` | 50.0 | 16.0 | None | 36 | 9 events | **1 Event(s)** |
| 8 | [`power_theta50_rho1_rep00.fas`](power_theta50_rho1_rep00.fas) | `power` | 50.0 | 1.0 | None | 2 | 0 events | Certified Null (0 events) |
| 9 | [`power_theta50_rho4_rep00.fas`](power_theta50_rho4_rep00.fas) | `power` | 50.0 | 4.0 | None | 14 | 4 events | Certified Null (0 events) |
| 10 | [`power_theta50_rho64_rep00.fas`](power_theta50_rho64_rep00.fas) | `power` | 50.0 | 64.0 | None | 154 | 52 events | **1 Event(s)** |

---
## Detailed Per-Replicate Breakdown

### 1. `fp_theta100_alpha0p5_rep00.fas`
- **Scenario**: `fp_theta100_alpha0p5` (Category: `false_positive`)
- **Parameters**: $\theta = 100.0$, $\rho = 0.0$, $\alpha = 0.5$, $N = 10$, $L = 1000$ bp
- **Partitions**: 1 coalescent segments
- **True Breakpoints**: `[]`
- **Ground Truth Events (0)**:
  - None (Strict clonal non-recombinant null)
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 2. `fp_theta50_alpha0p05_rep00.fas`
- **Scenario**: `fp_theta50_alpha0p05` (Category: `false_positive`)
- **Parameters**: $\theta = 50.0$, $\rho = 0.0$, $\alpha = 0.05$, $N = 10$, $L = 1000$ bp
- **Partitions**: 1 coalescent segments
- **True Breakpoints**: `[]`
- **Ground Truth Events (0)**:
  - None (Strict clonal non-recombinant null)
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 3. `power_theta10_rho0_rep00.fas`
- **Scenario**: `power_theta10_rho0` (Category: `power`)
- **Parameters**: $\theta = 10.0$, $\rho = 0.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 1 coalescent segments
- **True Breakpoints**: `[]`
- **Ground Truth Events (0)**:
  - None (Strict clonal non-recombinant null)
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 4. `power_theta10_rho16_rep00.fas`
- **Scenario**: `power_theta10_rho16` (Category: `power`)
- **Parameters**: $\theta = 10.0$, $\rho = 16.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 51 coalescent segments
- **True Breakpoints**: `[10, 12, 22, 32, 40, 67, 78, 119, 136, 168, 198, 217, 252, 255, 274]...`
- **Ground Truth Events (12)**:
  - Child: `Taxon_0`, Home: `Taxon_3`, Donor: `Taxon_4`, Tract: `[78, 119)` (1 segregating SNPs)
  - Child: `Taxon_2`, Home: `Taxon_0`, Donor: `Taxon_3`, Tract: `[0, 136)` (1 segregating SNPs)
  - Child: `Taxon_2`, Home: `Taxon_0`, Donor: `Taxon_9`, Tract: `[800, 1000)` (4 segregating SNPs)
  - Child: `Taxon_3`, Home: `Taxon_0`, Donor: `Taxon_2`, Tract: `[0, 136)` (1 segregating SNPs)
  - Child: `Taxon_4`, Home: `Taxon_5`, Donor: `Taxon_7`, Tract: `[0, 373)` (5 segregating SNPs)
  - *(and 7 additional small/silent micro-tracts)*
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 5. `power_theta200_rho16_rep00.fas`
- **Scenario**: `power_theta200_rho16` (Category: `power`)
- **Parameters**: $\theta = 200.0$, $\rho = 16.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 37 coalescent segments
- **True Breakpoints**: `[12, 26, 27, 36, 58, 94, 175, 198, 206, 239, 361, 367, 373, 403, 414]...`
- **Ground Truth Events (19)**:
  - Child: `Taxon_0`, Home: `Taxon_7`, Donor: `Taxon_4`, Tract: `[403, 502)` (23 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_7`, Donor: `Taxon_5`, Tract: `[502, 677)` (41 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_7`, Donor: `Taxon_6`, Tract: `[677, 854)` (23 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_7`, Donor: `Taxon_9`, Tract: `[854, 1000)` (4 segregating SNPs)
  - Child: `Taxon_2`, Home: `Taxon_4`, Donor: `Taxon_5`, Tract: `[206, 361)` (10 segregating SNPs)
  - *(and 14 additional small/silent micro-tracts)*
- **RhizAeon Prediction (2)**:
  - Recombinant: `Taxon_2`, Home: `Taxon_8`, Donor: `Taxon_3`, Tract: `687-1000` bp (p-value: `N/A`, Z: `0.00`)
  - Recombinant: `Taxon_0`, Home: `Taxon_7`, Donor: `Taxon_2`, Tract: `375-639` bp (p-value: `N/A`, Z: `0.00`)

### 6. `power_theta50_rho0_rep00.fas`
- **Scenario**: `power_theta50_rho0` (Category: `power`)
- **Parameters**: $\theta = 50.0$, $\rho = 0.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 1 coalescent segments
- **True Breakpoints**: `[]`
- **Ground Truth Events (0)**:
  - None (Strict clonal non-recombinant null)
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 7. `power_theta50_rho16_rep00.fas`
- **Scenario**: `power_theta50_rho16` (Category: `power`)
- **Parameters**: $\theta = 50.0$, $\rho = 16.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 36 coalescent segments
- **True Breakpoints**: `[2, 92, 113, 127, 130, 152, 167, 195, 196, 244, 251, 256, 260, 264, 277]...`
- **Ground Truth Events (9)**:
  - Child: `Taxon_0`, Home: `Taxon_1`, Donor: `Taxon_5`, Tract: `[283, 371)` (1 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_1`, Donor: `Taxon_7`, Tract: `[371, 644)` (3 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_1`, Donor: `Taxon_3`, Tract: `[717, 784)` (1 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_1`, Donor: `Taxon_6`, Tract: `[784, 1000)` (12 segregating SNPs)
  - Child: `Taxon_2`, Home: `Taxon_1`, Donor: `Taxon_7`, Tract: `[0, 371)` (2 segregating SNPs)
  - *(and 4 additional small/silent micro-tracts)*
- **RhizAeon Prediction (1)**:
  - Recombinant: `Taxon_7`, Home: `Taxon_0`, Donor: `Taxon_5`, Tract: `651-1000` bp (p-value: `N/A`, Z: `0.00`)

### 8. `power_theta50_rho1_rep00.fas`
- **Scenario**: `power_theta50_rho1` (Category: `power`)
- **Parameters**: $\theta = 50.0$, $\rho = 1.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 2 coalescent segments
- **True Breakpoints**: `[56]`
- **Ground Truth Events (0)**:
  - None (Strict clonal non-recombinant null)
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 9. `power_theta50_rho4_rep00.fas`
- **Scenario**: `power_theta50_rho4` (Category: `power`)
- **Parameters**: $\theta = 50.0$, $\rho = 4.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 14 coalescent segments
- **True Breakpoints**: `[10, 151, 179, 213, 422, 465, 477, 650, 687, 793, 893, 923, 928]`
- **Ground Truth Events (4)**:
  - Child: `Taxon_3`, Home: `Taxon_4`, Donor: `Taxon_6`, Tract: `[893, 1000)` (7 segregating SNPs)
  - Child: `Taxon_4`, Home: `Taxon_3`, Donor: `Taxon_1`, Tract: `[893, 1000)` (1 segregating SNPs)
  - Child: `Taxon_6`, Home: `Taxon_2`, Donor: `Taxon_3`, Tract: `[893, 1000)` (1 segregating SNPs)
  - Child: `Taxon_7`, Home: `Taxon_1`, Donor: `Taxon_0`, Tract: `[650, 1000)` (19 segregating SNPs)
- **RhizAeon Prediction (0)**:
  - Zero false alarms called (H0 Null preserved)

### 10. `power_theta50_rho64_rep00.fas`
- **Scenario**: `power_theta50_rho64` (Category: `power`)
- **Parameters**: $\theta = 50.0$, $\rho = 64.0$, $\alpha = None$, $N = 10$, $L = 1000$ bp
- **Partitions**: 154 coalescent segments
- **True Breakpoints**: `[12, 19, 22, 30, 43, 56, 108, 139, 145, 151, 154, 159, 162, 163, 186]...`
- **Ground Truth Events (52)**:
  - Child: `Taxon_0`, Home: `Taxon_2`, Donor: `Taxon_8`, Tract: `[217, 245)` (4 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_2`, Donor: `Taxon_1`, Tract: `[315, 370)` (8 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_2`, Donor: `Taxon_8`, Tract: `[370, 404)` (2 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_2`, Donor: `Taxon_3`, Tract: `[744, 785)` (1 segregating SNPs)
  - Child: `Taxon_0`, Home: `Taxon_2`, Donor: `Taxon_3`, Tract: `[838, 963)` (4 segregating SNPs)
  - *(and 47 additional small/silent micro-tracts)*
- **RhizAeon Prediction (1)**:
  - Recombinant: `Taxon_9`, Home: `Taxon_7`, Donor: `Taxon_4`, Tract: `541-1000` bp (p-value: `N/A`, Z: `0.00`)
