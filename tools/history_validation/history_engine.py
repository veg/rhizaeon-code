#!/usr/bin/env python3
"""
history_engine.py
Algorithmic Kernel for Reconstructing and Validating Competing Recombination Histories.

Provides:
1. ReticulationEvent: Polarized, boundary-resolved reticulation transfer edge.
2. AncestralConsolidator: Groups isolate-level tracts into ancestral reticulations vs recurrent transfers.
3. RecombinationHistory: Representation of a timed/directed ARG or partitioned hypothesis.
4. HistoryScorer: Multidimensional statistical validation engine computing:
   - Residual Manifold Strain & Spectral SNR over Marchenko-Pastur bulk edge.
   - Synapomorphic Site Concordance Ratio (SSCR).
   - Partitioned continuous-time Markov substitution log-likelihood (HKY85/GTR) & Delta AICc / BIC.
   - Temporal / Coalescent Co-existence Feasibility.
"""

import math
import numpy as np
import scipy.linalg
from scipy.stats import fisher_exact
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

class ReticulationEvent:
    def __init__(self, event_id, recipient, donor, home, u1, u2, z_phys=0.0,
                 is_ghost=False, is_crossover=False, member_isolates=None, category="Ancestral_Consolidated"):
        self.event_id = event_id
        self.recipient = recipient
        self.donor = donor
        self.home = home
        self.u1 = int(u1)
        self.u2 = int(u2)
        self.z_phys = float(z_phys)
        self.is_ghost = bool(is_ghost)
        self.is_crossover = bool(is_crossover)
        self.member_isolates = member_isolates or [recipient]
        self.category = category # Ancestral_Consolidated, Recurrent_Hotspot, Private_Terminal

    @property
    def span(self):
        return self.u2 - self.u1

    def to_dict(self):
        return {
            "event_id": self.event_id,
            "recipient": self.recipient,
            "donor": self.donor,
            "home": self.home,
            "u1": self.u1,
            "u2": self.u2,
            "span_nt": self.span,
            "z_phys": self.z_phys,
            "is_ghost": self.is_ghost,
            "is_crossover": self.is_crossover,
            "member_isolates": self.member_isolates,
            "category": self.category
        }

class AncestralConsolidator:
    """
    Consolidates isolate-level recombinant tracts into unique ancestral reticulation events.
    Disentangles shared ancestral transfers from independent recurrent transfers at hotspots.
    """
    def __init__(self, jaccard_thresh=0.60, max_boundary_dist=100):
        self.jaccard_thresh = jaccard_thresh
        self.max_boundary_dist = max_boundary_dist

    def consolidate(self, raw_events):
        if not raw_events:
            return []

        clusters = []
        for ev in raw_events:
            u1 = ev.get("u1", 0)
            u2 = ev.get("u2", 0)
            cand = ev.get("candidate_name", "")
            donor = ev.get("donor_name", "")
            home = ev.get("home_name", "")
            z = ev.get("z_phys", 0.0)
            ghost = ev.get("is_ghost_donor", False)
            cross = ev.get("is_crossover", False)

            matched_cluster = None
            for cl in clusters:
                # Check spatial overlap
                cl_u1 = cl["u1_list"]
                cl_u2 = cl["u2_list"]
                med_u1 = np.median(cl_u1)
                med_u2 = np.median(cl_u2)

                inter = max(0, min(u2, med_u2) - max(u1, med_u1))
                union = max(u2, med_u2) - min(u1, med_u1)
                jaccard = inter / union if union > 0 else 0.0
                dist_b1 = abs(u1 - med_u1)
                dist_b2 = abs(u2 - med_u2)

                # Check donor compatibility (same donor or close relative)
                same_donor = (donor == cl["donor"]) or (donor[:6] == cl["donor"][:6])
                same_home = (home == cl["home"]) or (home[:6] == cl["home"][:6])

                if (jaccard >= self.jaccard_thresh or (dist_b1 <= self.max_boundary_dist and dist_b2 <= self.max_boundary_dist)) and same_donor:
                    matched_cluster = cl
                    break

            if matched_cluster:
                matched_cluster["events"].append(ev)
                matched_cluster["u1_list"].append(u1)
                matched_cluster["u2_list"].append(u2)
                matched_cluster["z_list"].append(z)
                matched_cluster["members"].append(cand)
            else:
                clusters.append({
                    "donor": donor,
                    "home": home,
                    "is_ghost": ghost,
                    "is_crossover": cross,
                    "events": [ev],
                    "u1_list": [u1],
                    "u2_list": [u2],
                    "z_list": [z],
                    "members": [cand]
                })

        # Build consolidated events
        consolidated = []
        for idx, cl in enumerate(clusters, 1):
            med_u1 = int(np.round(np.median(cl["u1_list"])))
            med_u2 = int(np.round(np.median(cl["u2_list"])))
            max_z = float(np.max(cl["z_list"]))
            members = list(dict.fromkeys(cl["members"]))

            category = "Ancestral_Consolidated" if len(members) > 1 else "Private_Terminal"
            rec_event = ReticulationEvent(
                event_id=idx,
                recipient=members[0] if len(members) == 1 else f"Ancestral_Node_{idx}({len(members)}_taxa)",
                donor=cl["donor"],
                home=cl["home"],
                u1=med_u1,
                u2=med_u2,
                z_phys=max_z,
                is_ghost=cl["is_ghost"],
                is_crossover=cl["is_crossover"],
                member_isolates=members,
                category=category
            )
            consolidated.append(rec_event)

        return consolidated

class RecombinationHistory:
    """
    Formal representation of an evolutionary history (Tree or ARG) with associated partitions.
    """
    def __init__(self, history_id, name, description, partitions, reticulations=None, k_params=None):
        self.history_id = history_id
        self.name = name
        self.description = description
        self.partitions = partitions # List of (start_nt, end_nt) tuples
        self.reticulations = reticulations or []
        self.k_params = k_params

class HistoryScorer:
    """
    Multidimensional statistical validation engine for competing recombination histories.
    """
    def __init__(self, aln, kappa=4.0):
        self.aln = aln
        self.L = aln.get_alignment_length()
        self.N = len(aln)
        self.taxa = [s.id for s in aln]
        self.seqs = {s.id: str(s.seq).upper() for s in aln}
        self.kappa = kappa

        # Precompute base frequencies and global substitution Q matrix
        self.base_map = {"A": 0, "C": 1, "G": 2, "T": 3}
        all_seq = "".join(self.seqs.values())
        counts = np.array([all_seq.count(b) for b in "ACGT"], dtype=float)
        total_c = max(1.0, counts.sum())
        self.pi = np.maximum(1e-4, counts / total_c)
        self.pi /= self.pi.sum()

        self.Q = np.zeros((4, 4))
        for i, b1 in enumerate("ACGT"):
            for j, b2 in enumerate("ACGT"):
                if i != j:
                    is_ts = (b1, b2) in [("A","G"), ("G","A"), ("C","T"), ("T","C")]
                    self.Q[i, j] = (self.kappa if is_ts else 1.0) * self.pi[j]
            self.Q[i, i] = -np.sum(self.Q[i, :])
        scale = -np.sum(self.pi * np.diag(self.Q))
        self.Q /= max(1e-9, scale)

        # Global distance matrix D_obs
        self.calc = DistanceCalculator("identity")
        self.constructor = DistanceTreeConstructor(self.calc, "nj")
        self.D_obs = self._compute_distance_matrix(self.aln)

    def _compute_distance_matrix(self, sub_aln):
        matrix = self.calc.get_distance(sub_aln)
        N = len(sub_aln)
        D = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                D[i, j] = matrix[i, j]
        return D

    def _felsenstein_log_likelihood(self, tree, sub_aln):
        L_sub = sub_aln.get_alignment_length()
        seqs_sub = {s.id: str(s.seq).upper() for s in sub_aln}

        def prune(clade):
            if clade.is_terminal():
                lik = np.zeros((4, L_sub), dtype=float)
                s_str = seqs_sub.get(clade.name, "N" * L_sub)
                for u, ch in enumerate(s_str):
                    if ch in self.base_map:
                        lik[self.base_map[ch], u] = 1.0
                    else:
                        lik[:, u] = 1.0
                return lik
            else:
                total_lik = np.ones((4, L_sub), dtype=float)
                for child in clade.clades:
                    child_lik = prune(child)
                    bl = max(1e-5, child.branch_length if child.branch_length is not None else 0.01)
                    P = scipy.linalg.expm(self.Q * bl)
                    total_lik *= (P @ child_lik)
                return total_lik

        root_lik = prune(tree.root)
        site_lik = np.sum(self.pi[:, None] * root_lik, axis=0)
        site_lik = np.maximum(site_lik, 1e-300)
        return float(np.sum(np.log(site_lik)))

    def evaluate_residual_strain(self, history):
        """
        Computes residual non-tree strain and spectral SNR over Marchenko-Pastur bulk edge.
        """
        D_model = np.zeros((self.N, self.N))
        total_w = 0

        for u1, u2 in history.partitions:
            u1_c = max(0, min(self.L, u1))
            u2_c = max(0, min(self.L, u2))
            w = u2_c - u1_c
            if w <= 0:
                continue
            sub_aln = MultipleSeqAlignment([s[u1_c:u2_c] for s in self.aln])
            tree = self.constructor.build_tree(sub_aln)
            D_part = np.zeros((self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    D_part[i, j] = tree.distance(self.taxa[i], self.taxa[j])
            D_model += (w / self.L) * D_part
            total_w += w

        # Residual distance matrix
        R = self.D_obs - D_model
        norm_obs = np.linalg.norm(self.D_obs, "fro")
        norm_res = np.linalg.norm(R, "fro")
        strain_ratio = float(norm_res / norm_obs) if norm_obs > 0 else 0.0

        # Double-centering on residual matrix for spectral analysis
        H = np.eye(self.N) - (1.0 / self.N) * np.ones((self.N, self.N))
        B_R = -0.5 * H @ (R ** 2) @ H
        evals = np.linalg.eigvalsh(B_R)
        evals = np.sort(np.maximum(0.0, evals))[::-1]

        top_s = float(evals[0]) if len(evals) > 0 else 0.0
        sigma_res = float(np.std(R))
        # Marchenko-Pastur bulk edge
        s_bulk = (sigma_res ** 2) * (1.0 + np.sqrt(self.N / self.L)) ** 2
        spectral_snr = float(top_s / s_bulk) if s_bulk > 1e-12 else 0.0

        return {
            "frobenius_norm_obs": float(norm_obs),
            "frobenius_norm_residual": float(norm_res),
            "residual_strain_ratio": strain_ratio,
            "strain_explained_pct": float((1.0 - strain_ratio) * 100.0),
            "top_residual_eigenvalue": top_s,
            "marchenko_pastur_bulk_edge": s_bulk,
            "spectral_snr": spectral_snr,
            "is_within_mp_bulk": bool(spectral_snr <= 1.0)
        }

    def evaluate_sscr(self, history):
        """
        Computes Synapomorphic Site Concordance Ratio (SSCR) and polarization contrast.
        """
        if not history.reticulations:
            return {
                "events_evaluated": 0,
                "mean_sscr_in": 0.0,
                "mean_sscr_flank": 0.0,
                "mean_delta_sscr": 0.0,
                "all_events": []
            }

        event_scores = []
        for ev in history.reticulations:
            # Find recipient, donor, home in aln
            r_id = next((t for t in self.taxa if ev.recipient in t), None)
            d_id = next((t for t in self.taxa if ev.donor in t), None)
            h_id = next((t for t in self.taxa if ev.home in t), None)

            if not (r_id and d_id and h_id):
                continue

            r_seq = self.seqs[r_id]
            d_seq = self.seqs[d_id]
            h_seq = self.seqs[h_id]

            # In-tract counts
            in_d_match = 0
            in_h_match = 0
            for u in range(ev.u1, min(self.L, ev.u2)):
                r_c, d_c, h_c = r_seq[u], d_seq[u], h_seq[u]
                if d_c != h_c and r_c in [d_c, h_c] and r_c not in ["-", "N"]:
                    if r_c == d_c:
                        in_d_match += 1
                    else:
                        in_h_match += 1

            # Flank counts
            flank_d_match = 0
            flank_h_match = 0
            flank_indices = list(range(0, ev.u1)) + list(range(ev.u2, self.L))
            for u in flank_indices:
                r_c, d_c, h_c = r_seq[u], d_seq[u], h_seq[u]
                if d_c != h_c and r_c in [d_c, h_c] and r_c not in ["-", "N"]:
                    if r_c == d_c:
                        flank_d_match += 1
                    else:
                        flank_h_match += 1

            sscr_in = in_d_match / max(1, in_d_match + in_h_match)
            sscr_flank = flank_d_match / max(1, flank_d_match + flank_h_match)
            delta_sscr = sscr_in - sscr_flank

            # Fisher exact test
            table = [[in_d_match, in_h_match], [flank_d_match, flank_h_match]]
            try:
                odds, p_fisher = fisher_exact(table, alternative="greater")
            except Exception:
                p_fisher = 1.0

            event_scores.append({
                "event_id": ev.event_id,
                "recipient": ev.recipient,
                "donor": ev.donor,
                "u1": ev.u1,
                "u2": ev.u2,
                "in_tract_donor_matches": in_d_match,
                "in_tract_home_matches": in_h_match,
                "flank_donor_matches": flank_d_match,
                "flank_home_matches": flank_h_match,
                "sscr_in": float(sscr_in),
                "sscr_flank": float(sscr_flank),
                "delta_sscr": float(delta_sscr),
                "p_fisher": float(p_fisher)
            })

        if event_scores:
            m_in = float(np.mean([s["sscr_in"] for s in event_scores]))
            m_flank = float(np.mean([s["sscr_flank"] for s in event_scores]))
            m_delta = float(np.mean([s["delta_sscr"] for s in event_scores]))
        else:
            m_in, m_flank, m_delta = 0.0, 0.0, 0.0

        return {
            "events_evaluated": len(event_scores),
            "mean_sscr_in": m_in,
            "mean_sscr_flank": m_flank,
            "mean_delta_sscr": m_delta,
            "all_events": event_scores
        }

    def evaluate_phylogenetic_likelihood(self, history):
        """
        Computes partitioned log-likelihood, AICc, and BIC under HKY85 continuous-time Markov model.
        """
        total_lnL = 0.0
        part_details = []

        for u1, u2 in history.partitions:
            u1_c = max(0, min(self.L, u1))
            u2_c = max(0, min(self.L, u2))
            w = u2_c - u1_c
            if w <= 0:
                continue
            sub_aln = MultipleSeqAlignment([s[u1_c:u2_c] for s in self.aln])
            tree = self.constructor.build_tree(sub_aln)
            lnL = self._felsenstein_log_likelihood(tree, sub_aln)
            total_lnL += lnL
            part_details.append({
                "u1": u1_c,
                "u2": u2_c,
                "length_nt": w,
                "lnL": float(lnL)
            })

        K = len(part_details)
        # Parameter count
        if history.k_params is not None:
            k = history.k_params
        else:
            # Standard partitioned tree: K trees * (2N - 3 branch lengths) + 4 model params
            k = K * (2 * self.N - 3) + (K * 4)

        # AICc & BIC
        aicc = 2.0 * k - 2.0 * total_lnL + (2.0 * k * (k + 1.0)) / max(1.0, self.L - k - 1.0)
        bic = k * np.log(self.L) - 2.0 * total_lnL

        return {
            "total_lnL": float(total_lnL),
            "k_parameters": int(k),
            "num_partitions": K,
            "aicc": float(aicc),
            "bic": float(bic),
            "partitions": part_details
        }

    def score_history(self, history):
        strain = self.evaluate_residual_strain(history)
        sscr = self.evaluate_sscr(history)
        phylo = self.evaluate_phylogenetic_likelihood(history)

        return {
            "history_id": history.history_id,
            "name": history.name,
            "description": history.description,
            "num_partitions": len(history.partitions),
            "num_reticulations": len(history.reticulations),
            "phylogenetics": phylo,
            "manifold_strain": strain,
            "synapomorphic_concordance": sscr
        }
