use rhizaeon_core::{Alignment, RhizAeonEngine};
use std::time::Instant;

// Include the simulation functions from grand_sweep
const PI: [f64; 4] = [0.30, 0.20, 0.20, 0.30];
const KAPPA: f64 = 2.5;

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

#[inline(always)]
fn sample_standard_normal(rng: &mut FastRng) -> f64 {
    let u1 = rng.next_f64().max(1e-15);
    let u2 = rng.next_f64();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

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

#[inline(always)]
fn compute_hky85_p_mat(t: f64, kappa: f64, pi: &[f64; 4]) -> [[f64; 4]; 4] {
    let pi_r = pi[0] + pi[2];
    let pi_y = pi[1] + pi[3];
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
                pi[j] * (1.0 - exp_u)
            } else {
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

fn generate_trial_alignment(trial_id: usize, seed_base: u64) -> (Alignment, usize, usize, usize, usize, usize) {
    let mut rng = FastRng::new(seed_base + (trial_id as u64).wrapping_mul(7919) + 31);
    let is_recombinant = trial_id % 2 == 1;

    let log_l = rng.gen_range_f64((1000.0f64).ln(), (30000.0f64).ln());
    let length = ((log_l.exp() / 50.0).round() * 50.0) as usize;

    let r_n = rng.next_f64();
    let num_taxa = if r_n < 0.40 {
        rng.gen_range_usize(4, 30)
    } else if r_n < 0.80 {
        rng.gen_range_usize(31, 100)
    } else {
        rng.gen_range_usize(101, 250)
    };

    let log_d = rng.gen_range_f64((0.005f64).ln(), (0.25f64).ln());
    let d_parent = log_d.exp();

    let heterotachy_alpha = if rng.next_f64() < 0.25 {
        1000.0
    } else {
        let log_a = rng.gen_range_f64((0.05f64).ln(), (5.0f64).ln());
        log_a.exp()
    };

    let transmission_drift = rng.gen_range_f64(0.0005, 0.0035);

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

    let p1_anc = evolve_sequence(&anc_root, d_parent / 2.0, site_rates.as_deref(), &mut rng);
    let p2_anc = evolve_sequence(&anc_root, d_parent / 2.0, site_rates.as_deref(), &mut rng);

    let mut mut_total = 0;
    for u in 0..length {
        if p1_anc[u] != p2_anc[u] {
            mut_total += 1;
        }
    }

    let (k_true, tract_len, u1_true, u2_true, mut_tr) = if is_recombinant {
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

    for i in 0..n_p1 {
        taxa_names.push(format!("P1_{}", i));
        let leaf = evolve_sequence(&p1_anc, transmission_drift, site_rates.as_deref(), &mut rng);
        for s in leaf {
            full_matrix.push(s + 1);
        }
    }

    for i in 0..n_p2 {
        taxa_names.push(format!("P2_{}", i));
        let leaf = evolve_sequence(&p2_anc, transmission_drift, site_rates.as_deref(), &mut rng);
        for s in leaf {
            full_matrix.push(s + 1);
        }
    }

    if is_recombinant && k_true > 0 {
        let mut rec_anc = p1_anc.clone();
        for u in (u1_true - 1)..u2_true {
            rec_anc[u] = p2_anc[u];
        }

        for k in 0..k_true {
            taxa_names.push(format!("Rec_{}", k));
            let leaf = evolve_sequence(&rec_anc, transmission_drift, site_rates.as_deref(), &mut rng);
            for s in leaf {
                full_matrix.push(s + 1);
            }
        }
    }

    (
        Alignment::new(taxa_names, length, full_matrix),
        k_true,
        u1_true,
        u2_true,
        mut_tr,
        mut_total,
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seed_base = 42;
    let missed_trials = [
        109, 15739, 14997, 17349, 17367, 17385, 17021, 4885, 17123, 5055,
        17995, 4627, 14613, 2775, 2987, 18455, 18605, 9475, 5655, 43995,
        7185, 7249, 44107, 36059, 28569, 40969, 23217, 41059, 44369, 28895,
        41187, 38257, 44543, 28991, 40051, 32247, 30637, 36675, 30751, 36689,
        47805, 38585, 39615, 30777, 41559, 41571, 10405, 47955, 30327, 37109,
        35777, 35809, 37367, 9091, 32845, 37249, 48729, 43031, 11187, 11963,
        33305, 43567, 11701
    ];

    println!("========================================================");
    println!("Testing all {} previously missed trials with M >= 500:", missed_trials.len());
    println!("========================================================");

    let engine_default = RhizAeonEngine::new();
    let mut engine_raw = RhizAeonEngine::new();
    engine_raw.compress_snps = Some(false);

    let mut def_detected = 0;
    let mut raw_detected = 0;

    for &tid in &missed_trials {
        let (aln_t, k_t, u1_t, u2_t, m_tr, _) = generate_trial_alignment(tid, seed_base);
        let res_def = engine_default.scan(&aln_t);
        let res_uncomp = engine_raw.scan(&aln_t);

        let def_det = res_def.events.len() > 0;
        let raw_det = res_uncomp.events.len() > 0;

        if def_det { def_detected += 1; }
        if raw_det { raw_detected += 1; }

        println!("Trial {:5} (N={:2}, L={:5}, M_tr={:4}, True=[{:5}, {:5}]): Default={} | Raw={}",
            tid, aln_t.num_taxa, aln_t.length, m_tr, u1_t, u2_t,
            if def_det { "DET" } else { "MIS" },
            if raw_det { "DET" } else { "MIS" }
        );
    }

    println!("========================================================");
    println!("Summary over {} previously missed trials:", missed_trials.len());
    println!("  Default Engine Detected:     {}/{}", def_detected, missed_trials.len());
    println!("  Uncompressed Raw Alignment:  {}/{}", raw_detected, missed_trials.len());
    println!("========================================================");

    println!("\n========================================================");
    println!("Detailed Diagnostic Trace for Trial 17021:");
    println!("========================================================");
    let tid = 17021;
    let (aln, k_t, u1_t, u2_t, m_tr, m_tot) = generate_trial_alignment(tid, seed_base);
    println!("Trial {}: N={}, L={}, k_true={}, True=[{}, {}], M_tr={}, M_tot={}",
        tid, aln.num_taxa, aln.length, k_t, u1_t, u2_t, m_tr, m_tot);
    for (i, name) in aln.taxa.iter().enumerate() {
        println!("  Taxon {}: {}", i, name);
    }

    let mut engine = RhizAeonEngine::new();
    engine.compress_snps = Some(false);

    let (aln_trimmed, u_offset) = rhizaeon_core::fasta::trim_coverage_envelope(&aln);
    let n = aln_trimmed.num_taxa;
    let l = aln_trimmed.length;
    let m_inf = rhizaeon_core::fasta::count_parsimony_informative_sites(&aln_trimmed);
    let m_min = rhizaeon_core::fasta::compute_min_informative_sites(n);
    println!("m_inf = {}, m_min = {}", m_inf, m_min);

    let effective_k = engine.target_landmarks.max(n / 2).min(32).min(n);
    let landmarks = rhizaeon_core::landmarks::select_landmarks(&aln_trimmed, effective_k);
    println!("Selected {} landmarks: {:?}", landmarks.num_landmarks, landmarks.indices);

    let rope = rhizaeon_core::tree_rope::TreeRoPE::new(engine.dim);
    let m_matrix = rope.compute_query_landmark_rotary_factors(&aln_trimmed, &landmarks);
    let prior = rhizaeon_core::prior::PhyloBias::new(engine.epsilon_0, engine.decay_lambda);
    let b_prior = prior.compute_prior_matrix(&aln_trimmed, &landmarks);
    let root_seq = rhizaeon_core::soft_root::compute_consensus_soft_root(&aln_trimmed, None);

    let attention = rhizaeon_core::attention::PhyloAttentionEngine::new(engine.dim, engine.root_sink_factor);
    let (a_tensor, f_tensor, c_tensor) = attention.compute_attention_and_force_tensors(
        &aln_trimmed, &landmarks, &m_matrix, &b_prior, &root_seq,
    );

    let screener = rhizaeon_core::screener::LocalFluxScreener::new();
    let screening = screener.screen(&f_tensor, &c_tensor, &aln_trimmed, &landmarks);
    println!("Screening is_h0 = {}", screening.is_non_recombinant_h0);
    println!("Screening candidates ({}): {:?}", screening.candidate_recombinant_indices.len(), screening.candidate_recombinant_indices);
    for i in 0..n {
        println!("  Taxon {:2} ({}): v_cont={:.4}, v_root={:.4}, r_root={:.4}, gamma_root={:.4}",
            i, aln.taxa[i], screening.v_cont[i], screening.v_root[i], screening.r_root[i], screening.gamma_root[i]);
    }
    println!("Screening rate bursts: {:?}", screening.rate_burst_indices);
    println!("Screening ghost candidates: {:?}", screening.ghost_candidate_indices);

    let resolver = rhizaeon_core::attribution::TrajectoryParentalResolver::new();
    let polisher = rhizaeon_core::polishing::CumulativeTrajectoryPolisher::new();
    let due_diligence = rhizaeon_core::due_diligence::SequenceDueDiligence::new(engine.poisson_floor, engine.alpha, 100_000);

    let num_channels = landmarks.num_landmarks + 1;
    for &cand_idx in &screening.candidate_recombinant_indices {
        let parental_pairs = resolver.resolve_all_reticulation_channels(
            cand_idx, &c_tensor, &a_tensor, l, n, &landmarks, &aln_trimmed, &screening.ghost_candidate_indices,
        );
        println!("Candidate {} ({}): {} reticulation pairs found", cand_idx, aln.taxa[cand_idx], parental_pairs.len());
        for (p_i, pair) in parental_pairs.iter().enumerate() {
            println!("  Pair {}: home_ch={}, donor_ch={}",
                p_i, pair.home_channel, pair.donor_channel);
            let tracts = polisher.deconvolve_candidate(
                cand_idx, pair, &c_tensor, l, n, num_channels,
            );
            println!("    Deconvolved {} tracts", tracts.len());
            for tr in &tracts {
                let mut w_tr = vec![0.0; num_channels];
                let (s_i, e_i) = (tr.u1_discrete - 1, tr.u2_discrete - 1);
                for ch in 0..num_channels {
                    let c_u1 = if s_i > 0 { c_tensor[(s_i - 1) * n * num_channels + cand_idx * num_channels + ch] } else { 0.0 };
                    let c_u2 = c_tensor[e_i * n * num_channels + cand_idx * num_channels + ch];
                    w_tr[ch] = c_u2 - c_u1;
                }
                println!("        w_tr: root(ch0)={:.2}, max_cont={:.2}",
                    w_tr[0], w_tr[1..].iter().cloned().fold(f64::NEG_INFINITY, f64::max));
                for ch in 0..num_channels {
                    let name = if ch == 0 { "ROOT".to_string() } else { landmarks.taxa[ch - 1].clone() };
                    if w_tr[ch] > 10.0 {
                        println!("          ch {:2} ({}): w_tr={:.2}", ch, name, w_tr[ch]);
                    }
                }
                let local_parental = resolver.resolve_reticulation_channels_localized(
                    cand_idx, tr.u1_discrete, tr.u2_discrete, &c_tensor, l, n, &landmarks, &aln_trimmed,
                    Some(tr.home_channel), Some(tr.donor_channel), &screening.candidate_recombinant_indices,
                    &screening.ghost_candidate_indices,
                );
                println!("        Local parental: is_ghost={}, home_glob={:?}, donor_glob={:?}, home_name={}, donor_name={}",
                    local_parental.is_ghost_donor, local_parental.home_global_idx, local_parental.donor_global_idx,
                    local_parental.home_name, local_parental.donor_name);

                let tract_opt = due_diligence.verify_tract(
                    tr, &local_parental, &aln_trimmed, &landmarks, &screening.rate_burst_indices, &screening.candidate_recombinant_indices,
                );
                match tract_opt {
                    Some(ev) => println!("        Verified: {} (p_fisher={:.4e})", ev.is_verified, ev.p_fisher),
                    None => println!("        Due diligence REJECTED tract!"),
                }
            }
        }
    }

    Ok(())
}

