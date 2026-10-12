use clap::Parser;
use rayon::prelude::*;
use rhizaeon_core::{Alignment, RhizAeonEngine};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

/// Exact nucleotide frequencies (HKY85)
const PI: [f64; 4] = [0.30, 0.20, 0.20, 0.30];
const KAPPA: f64 = 2.5;

/// High-speed 64-bit pseudo-random number generator (XorShift64Star)
struct FastRng {
    state: u64,
}

impl FastRng {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x853c49e6748fea9b } else { seed },
        }
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545f4914f6cdd1d)
    }

    #[inline(always)]
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    #[inline(always)]
    fn gen_range_f64(&mut self, low: f64, high: f64) -> f64 {
        low + self.next_f64() * (high - low)
    }

    #[inline(always)]
    fn gen_range_usize(&mut self, low: usize, high: usize) -> usize {
        if low >= high {
            return low;
        }
        low + (self.next_u64() as usize % (high - low + 1))
    }
}

/// Marsaglia and Tsang (2000) standard normal variate
#[inline(always)]
fn sample_standard_normal(rng: &mut FastRng) -> f64 {
    let u1 = rng.next_f64().max(1e-15);
    let u2 = rng.next_f64();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// Marsaglia and Tsang (2000) exact Gamma distribution sampler
fn sample_gamma(alpha: f64, rng: &mut FastRng) -> f64 {
    if alpha < 1.0 {
        let g = sample_gamma(alpha + 1.0, rng);
        let u = rng.next_f64().max(1e-15);
        return g * u.powf(1.0 / alpha);
    }
    let d = alpha - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    loop {
        let z = sample_standard_normal(rng);
        let v = 1.0 + c * z;
        if v <= 0.0 {
            continue;
        }
        let v = v * v * v;
        let u = rng.next_f64();
        if u < 1.0 - 0.0331 * z * z * z * z {
            return d * v;
        }
        if u.ln() < 0.5 * z * z + d * (1.0 - v + v.ln()) {
            return d * v;
        }
    }
}

/// Exact analytical HKY85 transition matrix for branch length `t`
#[inline(always)]
fn compute_hky85_p_mat(t: f64, kappa: f64, pi: &[f64; 4]) -> [[f64; 4]; 4] {
    let pi_r = pi[0] + pi[2]; // Purines A, G
    let pi_y = pi[1] + pi[3]; // Pyrimidines C, T

    // Scaling coefficient beta
    let beta_coeff = pi_r * pi_y + kappa * (pi[0] * pi[2] + pi[1] * pi[3]);
    let u = t / (2.0 * beta_coeff);

    let exp_u = (-u).exp();
    let exp_u_r = (-(1.0 + (kappa - 1.0) * pi_r) * u).exp();
    let exp_u_y = (-(1.0 + (kappa - 1.0) * pi_y) * u).exp();

    let mut p = [[0.0; 4]; 4];

    for i in 0..4 {
        let is_purine = i == 0 || i == 2;
        let pi_g = if is_purine { pi_r } else { pi_y };
        let exp_u_g = if is_purine { exp_u_r } else { exp_u_y };

        let mut sum_non_diag = 0.0;
        for j in 0..4 {
            if i == j {
                continue;
            }
            let is_ts = (is_purine && (j == 0 || j == 2)) || (!is_purine && (j == 1 || j == 3));
            let p_ij = if !is_ts {
                // Transversion
                pi[j] * (1.0 - exp_u)
            } else {
                // Transition
                pi[j] + pi[j] * (1.0 / pi_g - 1.0) * exp_u - (pi[j] / pi_g) * exp_u_g
            };
            let clamped = p_ij.max(0.0);
            p[i][j] = clamped;
            sum_non_diag += clamped;
        }
        p[i][i] = (1.0 - sum_non_diag).max(0.0);
    }
    p
}

/// Evolve a nucleotide sequence along a branch with length `branch_len` and per-site rates
fn evolve_sequence(
    src: &[u8],
    branch_len: f64,
    site_rates: Option<&[f64]>,
    rng: &mut FastRng,
) -> Vec<u8> {
    let mut dest = Vec::with_capacity(src.len());
    let default_p_mat = if site_rates.is_none() {
        Some(compute_hky85_p_mat(branch_len, KAPPA, &PI))
    } else {
        None
    };

    for (u, &parent_state) in src.iter().enumerate() {
        let p_mat = match site_rates {
            Some(rates) => compute_hky85_p_mat(branch_len * rates[u], KAPPA, &PI),
            None => default_p_mat.unwrap(),
        };

        let row = &p_mat[parent_state as usize];
        let rand_val = rng.next_f64();

        let new_state = if rand_val < row[0] {
            0u8
        } else if rand_val < row[0] + row[1] {
            1u8
        } else if rand_val < row[0] + row[1] + row[2] {
            2u8
        } else {
            3u8
        };
        dest.push(new_state);
    }
    dest
}

#[derive(Debug)]
struct TrialResult {
    trial_id: usize,
    is_recombinant_true: u8,
    length: usize,
    num_taxa: usize,
    d_parent: f64,
    tract_len: usize,
    k_true: usize,
    u1_true: usize,
    u2_true: usize,
    mut_in_tract: usize,
    mut_total: usize,
    heterotachy_alpha: f64,
    informative_sites: usize,
    detected: u8,
    unique_events_inferred: usize,
    recall: f64,
    precision: f64,
    bg_fp_count: usize,
    bg_fpr: f64,
    jaccard: f64,
    bp1_error: usize,
    bp2_error: usize,
    run_time_ms: f64,
}

impl TrialResult {
    fn to_csv_row(&self) -> String {
        format!(
            "{},{},{},{},{:.5},{},{},{},{},{},{},{:.4},{},{},{},{:.4},{:.4},{},{:.5},{:.4},{},{},{:.2}",
            self.trial_id,
            self.is_recombinant_true,
            self.length,
            self.num_taxa,
            self.d_parent,
            self.tract_len,
            self.k_true,
            self.u1_true,
            self.u2_true,
            self.mut_in_tract,
            self.mut_total,
            self.heterotachy_alpha,
            self.informative_sites,
            self.detected,
            self.unique_events_inferred,
            self.recall,
            self.precision,
            self.bg_fp_count,
            self.bg_fpr,
            self.jaccard,
            self.bp1_error,
            self.bp2_error,
            self.run_time_ms
        )
    }
}

fn compute_jaccard(t_u1: usize, t_u2: usize, p_u1: usize, p_u2: usize) -> f64 {
    if p_u1 == 0 || p_u2 == 0 || p_u1 > p_u2 {
        return 0.0;
    }
    let start = t_u1.max(p_u1);
    let end = t_u2.min(p_u2);
    if end < start {
        return 0.0;
    }
    let inter = (end - start + 1) as f64;
    let union = ((t_u2 - t_u1 + 1) + (p_u2 - p_u1 + 1)) as f64 - inter;
    if union > 0.0 {
        inter / union
    } else {
        0.0
    }
}

/// Executes a single simulation trial completely in memory
fn run_single_simulation(trial_id: usize, seed_base: u64, engine: &RhizAeonEngine) -> TrialResult {
    let mut rng = FastRng::new(seed_base + (trial_id as u64).wrapping_mul(7919) + 31);

    // 1. Draw parameters from continuous phase space
    let is_recombinant = trial_id % 2 == 1;

    // Length: log-uniform in [1,000, 30,000]
    let log_l = rng.gen_range_f64((1000.0f64).ln(), (30000.0f64).ln());
    let length = ((log_l.exp() / 50.0).round() * 50.0) as usize;

    // Taxa count: 40% small [4, 30], 40% medium [31, 100], 20% large [101, 250]
    let r_n = rng.next_f64();
    let num_taxa = if r_n < 0.40 {
        rng.gen_range_usize(4, 30)
    } else if r_n < 0.80 {
        rng.gen_range_usize(31, 100)
    } else {
        rng.gen_range_usize(101, 250)
    };

    // Divergence: log-uniform in [0.005, 0.25]
    let log_d = rng.gen_range_f64((0.005f64).ln(), (0.25f64).ln());
    let d_parent = log_d.exp();

    // Heterotachy Gamma alpha: 25% strict clock (alpha = 1000.0), 75% Gamma [0.05, 5.0]
    let heterotachy_alpha = if rng.next_f64() < 0.25 {
        1000.0
    } else {
        let log_a = rng.gen_range_f64((0.05f64).ln(), (5.0f64).ln());
        log_a.exp()
    };

    let transmission_drift = rng.gen_range_f64(0.0005, 0.0035);

    // Site rates
    let site_rates: Option<Vec<f64>> = if heterotachy_alpha < 50.0 {
        let mut rates = Vec::with_capacity(length);
        for _ in 0..length {
            let r = sample_gamma(heterotachy_alpha, &mut rng) / heterotachy_alpha;
            rates.push(r.clamp(0.01, 50.0));
        }
        Some(rates)
    } else {
        None
    };

    // 2. Evolve Ancestral Root Sequence
    let mut anc_root = Vec::with_capacity(length);
    for _ in 0..length {
        let p = rng.next_f64();
        let s = if p < PI[0] {
            0u8
        } else if p < PI[0] + PI[1] {
            1u8
        } else if p < PI[0] + PI[1] + PI[2] {
            2u8
        } else {
            3u8
        };
        anc_root.push(s);
    }

    // Parental node ancestors
    let p1_anc = evolve_sequence(
        &anc_root,
        d_parent / 2.0,
        site_rates.as_deref(),
        &mut rng,
    );
    let p2_anc = evolve_sequence(
        &anc_root,
        d_parent / 2.0,
        site_rates.as_deref(),
        &mut rng,
    );

    // Total mutations between P1 and P2
    let mut mut_total = 0;
    for u in 0..length {
        if p1_anc[u] != p2_anc[u] {
            mut_total += 1;
        }
    }

    // Demography partition
    let (k_true, tract_len, u1_true, u2_true, mut_in_tract) = if is_recombinant {
        // Tract length ratio: log-uniform in [0.008, 0.50]
        let log_tr = rng.gen_range_f64((0.008f64).ln(), (0.50f64).ln());
        let t_ratio = log_tr.exp();
        let mut t_len = (length as f64 * t_ratio).round() as usize;
        t_len = t_len.clamp(30, length.saturating_sub(100));

        let u1 = rng.gen_range_usize(50, length.saturating_sub(t_len + 50));
        let u2 = u1 + t_len - 1;

        let mut mut_tr = 0;
        for u in (u1 - 1)..u2 {
            if p1_anc[u] != p2_anc[u] {
                mut_tr += 1;
            }
        }

        // Swarm size K
        let max_k = ((num_taxa as f64 * 0.35).floor() as usize).max(1);
        let k = rng.gen_range_usize(1, max_k);

        (k, t_len, u1, u2, mut_tr)
    } else {
        (0, 0, 0, 0, 0)
    };

    let n_bg_total = num_taxa - k_true;
    let n_p1 = (n_bg_total / 2).max(1);
    let n_p2 = n_bg_total.saturating_sub(n_p1).max(1);

    let mut taxa_names = Vec::with_capacity(num_taxa);
    let mut full_matrix = Vec::with_capacity(num_taxa * length);

    // Generate P1 background taxa
    for i in 0..n_p1 {
        taxa_names.push(format!("P1_{}", i));
        let leaf = evolve_sequence(
            &p1_anc,
            transmission_drift,
            site_rates.as_deref(),
            &mut rng,
        );
        for s in leaf {
            full_matrix.push(s + 1); // 1-indexed encoding
        }
    }

    // Generate P2 background taxa
    for i in 0..n_p2 {
        taxa_names.push(format!("P2_{}", i));
        let leaf = evolve_sequence(
            &p2_anc,
            transmission_drift,
            site_rates.as_deref(),
            &mut rng,
        );
        for s in leaf {
            full_matrix.push(s + 1);
        }
    }

    // Generate Recombinant taxa if applicable
    if is_recombinant && k_true > 0 {
        let mut rec_anc = p1_anc.clone();
        for u in (u1_true - 1)..u2_true {
            rec_anc[u] = p2_anc[u];
        }

        for k in 0..k_true {
            taxa_names.push(format!("Rec_{}", k));
            let leaf = evolve_sequence(
                &rec_anc,
                transmission_drift,
                site_rates.as_deref(),
                &mut rng,
            );
            for s in leaf {
                full_matrix.push(s + 1);
            }
        }
    }

    let aln = Alignment::new(taxa_names, length, full_matrix);

    // 3. Run RhizAeon Engine
    let t0 = Instant::now();
    let scan_res = engine.scan(&aln);
    let run_time_ms = t0.elapsed().as_secs_f64() * 1000.0;

    // 4. Score Metrics
    let mut rec_taxa_called = std::collections::HashSet::new();
    let mut bg_taxa_called = std::collections::HashSet::new();
    let mut best_jaccard = 0.0f64;
    let mut best_bp1_err = length;
    let mut best_bp2_err = length;

    let unique_events = scan_res.events.len();
    let detected = if unique_events > 0 { 1u8 } else { 0u8 };

    for ev in &scan_res.events {
        for iso in &ev.isolates {
            if iso.starts_with("Rec_") {
                rec_taxa_called.insert(iso.clone());
            } else {
                bg_taxa_called.insert(iso.clone());
            }
        }

        if is_recombinant {
            let j = compute_jaccard(u1_true, u2_true, ev.u1, ev.u2);
            if j > best_jaccard {
                best_jaccard = j;
                best_bp1_err = ev.u1.abs_diff(u1_true);
                best_bp2_err = ev.u2.abs_diff(u2_true);
            }
        }
    }

    let recall = if is_recombinant && k_true > 0 {
        rec_taxa_called.len() as f64 / k_true as f64
    } else {
        if detected == 0 {
            1.0
        } else {
            0.0
        }
    };

    let total_called = rec_taxa_called.len() + bg_taxa_called.len();
    let precision = if is_recombinant {
        if total_called > 0 {
            rec_taxa_called.len() as f64 / total_called as f64
        } else {
            1.0
        }
    } else {
        if detected == 0 {
            1.0
        } else {
            0.0
        }
    };

    let bg_fp_count = bg_taxa_called.len();
    let bg_fpr = if n_bg_total > 0 {
        bg_fp_count as f64 / n_bg_total as f64
    } else {
        0.0
    };

    let jaccard = if is_recombinant {
        best_jaccard
    } else {
        if detected == 0 {
            1.0
        } else {
            0.0
        }
    };

    TrialResult {
        trial_id,
        is_recombinant_true: if is_recombinant { 1 } else { 0 },
        length,
        num_taxa,
        d_parent,
        tract_len,
        k_true,
        u1_true,
        u2_true,
        mut_in_tract,
        mut_total,
        heterotachy_alpha,
        informative_sites: scan_res.informative_sites,
        detected,
        unique_events_inferred: unique_events,
        recall,
        precision,
        bg_fp_count,
        bg_fpr,
        jaccard,
        bp1_error: if is_recombinant && detected > 0 {
            best_bp1_err
        } else {
            0
        },
        bp2_error: if is_recombinant && detected > 0 {
            best_bp2_err
        } else {
            0
        },
        run_time_ms,
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "grand_sweep",
    author = "Sergei L. Kosakovsky Pond & Antigravity AI",
    version = "1.0",
    about = "Grand-Scale Parameter Sweep Harness for RhizAeon (50k-100k In-Memory Simulations)"
)]
struct Args {
    #[arg(short, long, default_value_t = 100000)]
    trials: usize,

    #[arg(short, long, default_value = "grand_sweep_results.csv")]
    output: String,

    #[arg(short = 'j', long, default_value_t = 0)]
    threads: usize,

    #[arg(short, long, default_value_t = 42)]
    seed: u64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap_or(());
    }

    let num_threads = rayon::current_num_threads();

    println!("================================================================================");
    println!("     RhizAeon Grand-Scale Parameter Sweep Harness (Continuous Phase Space)      ");
    println!("================================================================================");
    println!("Total Trials:         {}", args.trials);
    println!("Worker Threads:       {}", num_threads);
    println!("Output File:          {}", args.output);
    println!("Random Seed:          {}", args.seed);
    println!("--------------------------------------------------------------------------------");

    let engine = RhizAeonEngine::new();
    let file = File::create(&args.output)?;
    let writer = Mutex::new(BufWriter::with_capacity(4 * 1024 * 1024, file));

    // Write CSV header
    {
        let mut w = writer.lock().unwrap();
        writeln!(
            w,
            "trial_id,is_recombinant_true,length,num_taxa,d_parent,tract_len,k_true,u1_true,u2_true,mut_in_tract,mut_total,heterotachy_alpha,informative_sites,detected,unique_events_inferred,recall,precision,bg_fp_count,bg_fpr,jaccard,bp1_error,bp2_error,run_time_ms"
        )?;
    }

    let progress_counter = AtomicUsize::new(0);
    let total_recombinant_detected = AtomicUsize::new(0);
    let total_null_clean = AtomicUsize::new(0);
    let t_start = Instant::now();

    let trials: Vec<usize> = (0..args.trials).collect();

    // Process in parallel chunks to stream output efficiently across all cores
    trials.par_chunks(16).for_each(|chunk| {
        let mut local_csv = String::with_capacity(chunk.len() * 140);
        let mut rec_det_chunk = 0;
        let mut null_clean_chunk = 0;

        for &trial_id in chunk {
            let res = run_single_simulation(trial_id, args.seed, &engine);
            local_csv.push_str(&res.to_csv_row());
            local_csv.push('\n');

            if res.is_recombinant_true == 1 && res.detected == 1 {
                rec_det_chunk += 1;
            } else if res.is_recombinant_true == 0 && res.bg_fp_count == 0 {
                null_clean_chunk += 1;
            }
        }

        // Flush chunk to buffered writer
        {
            let mut w = writer.lock().unwrap();
            w.write_all(local_csv.as_bytes()).unwrap();
        }

        total_recombinant_detected.fetch_add(rec_det_chunk, Ordering::Relaxed);
        total_null_clean.fetch_add(null_clean_chunk, Ordering::Relaxed);

        let completed = progress_counter.fetch_add(chunk.len(), Ordering::Relaxed) + chunk.len();
        if completed % 1000 == 0 || completed == args.trials {
            let elapsed = t_start.elapsed().as_secs_f64();
            let rate = completed as f64 / elapsed;
            let remaining = (args.trials - completed) as f64 / rate;
            println!(
                "  Progress: {:6}/{} ({:5.1}%) | Rate: {:6.1} trials/s | Elapsed: {:5.1}s | ETA: {:5.1}s",
                completed,
                args.trials,
                (completed as f64 / args.trials as f64) * 100.0,
                rate,
                elapsed,
                remaining
            );
        }
    });

    let total_elapsed = t_start.elapsed().as_secs_f64();
    let total_rec_trials = args.trials / 2;
    let total_null_trials = args.trials - total_rec_trials;
    let rec_det = total_recombinant_detected.load(Ordering::Relaxed);
    let null_clean = total_null_clean.load(Ordering::Relaxed);

    println!("================================================================================");
    println!("                         GRAND SWEEP COMPLETE                                   ");
    println!("================================================================================");
    println!("Total Execution Time:      {:.2} seconds ({:.1} minutes)", total_elapsed, total_elapsed / 60.0);
    println!("Average Throughput:        {:.1} trials/second", args.trials as f64 / total_elapsed);
    println!("Recombinant Trials:        {} | Detected: {} ({:.2}%)", total_rec_trials, rec_det, (rec_det as f64 / total_rec_trials as f64) * 100.0);
    println!("Non-Recombinant Nulls:     {} | Zero FPR: {} ({:.4}%)", total_null_trials, null_clean, (null_clean as f64 / total_null_trials as f64) * 100.0);
    println!("Results saved to:          {}", args.output);
    println!("================================================================================");

    Ok(())
}
