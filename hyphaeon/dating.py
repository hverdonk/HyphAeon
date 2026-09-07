"""
hyphaeon/dating.py
------------------
Heterochronous Molecular Clock Calibration and Ancestor Dating (t_MRCA)
for Pathogen Genomics.

Methods:
1. Strict in-frame coding alignment validation (L_nt % 3 == 0, triplet-gap check, stop codon audit).
2. Flexible timestamp ingestion (FASTA headers, CSV/TSV metadata, Nextstrain Auspice JSON v2).
3. Root-to-tip divergence computation:
   - Tree-based: Patristic distance traversal with heuristic root search (TempEst R^2 maximization).
   - Tree-free: Direct pairwise distance estimation (TN93) and ancestral consensus anchoring.
4. Estimators:
   - Centered Root-to-Tip Ordinary Least Squares (OLS / TempEst emulation with delta-method & bootstrap CIs).
   - HyphAeon Attention-Derived Phylogenetic Generalized Least Squares (PGLS) via A_fused covariance.
   - Non-Linear Restricted Cubic Spline Clock (2 DF) and Power-Law Clock with hypothesis testing.
5. Historical outlier scoring & blind tip dating (e.g. dating the 1959 ZR59 archival isolate).
6. Publication-grade diagnostic visualization (PDF and PNG).
"""

import os
import sys
import json
import time
import math
import re
import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any

import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.stats as stats
import scipy.optimize as optimize
import torch

try:
    from Bio import Phylo
    HAS_BIOPHYLO = True
except ImportError:
    HAS_BIOPHYLO = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from .dataset import (
    parse_alignment_sequences,
    compute_tn93_distance_matrix,
    load_alignment_and_tree,
    GENETIC_CODE,
    CODON_TO_AA,
)
from .inference import (
    load_model,
    get_device,
    prepare_alignment,
)
from scipy.spatial.distance import pdist, squareform
from .splits import (
    extract_cross_taxa_attentions_and_embeddings,
    compute_fused_affinity_matrix,
)
from .temporal import parse_date_to_decimal, parse_dates_from_auspice_json, extract_date_from_string
from .io import ensure_parent_directory, write_json, write_csv


def compute_neural_covariance_kernel(
    cross_attn: np.ndarray,
    taxa_repr: Optional[np.ndarray] = None,
    mds_coords: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Computes a strictly positive semi-definite (PSD) phylogenetic correlation matrix
    fusing transformer cross-taxa attention profile similarity and 128D continuous sequence representations.
    Features are centered across taxa to eliminate representation anisotropy (the cone effect),
    ensuring proper phylogenetic decoupling between distinct clades and realistic confidence intervals.
    Guarantees unit diagonal C(i, i) = 1.0 and zero negative eigenvalues.
    """
    n = cross_attn.shape[0]

    # 1. Cross-Taxa Attention Profile Correlation Matrix (Centered)
    # Each row a_i is taxon i's attention distribution across all taxa.
    # Centering removes the background entropy baseline across taxa:
    a_centered = cross_attn - np.mean(cross_attn, axis=0, keepdims=True)
    a_cov = a_centered @ a_centered.T
    a_var = np.diag(a_cov)
    a_std = np.sqrt(np.maximum(a_var, 1e-12))[:, None]
    a_denom = a_std @ a_std.T
    K_attn = np.divide(a_cov, a_denom, where=(a_denom > 1e-12), out=np.eye(n))
    np.fill_diagonal(K_attn, 1.0)

    # 2. Continuous Latent Sequence Embedding Correlation Matrix (Centered)
    if taxa_repr is not None and len(taxa_repr) == n:
        # Centering removes the dominant common activation vector across taxa:
        z_centered = taxa_repr - np.mean(taxa_repr, axis=0, keepdims=True)
        z_cov = z_centered @ z_centered.T
        z_var = np.diag(z_cov)
        z_std = np.sqrt(np.maximum(z_var, 1e-12))[:, None]
        z_denom = z_std @ z_std.T
        K_emb = np.divide(z_cov, z_denom, where=(z_denom > 1e-12), out=np.eye(n))
        np.fill_diagonal(K_emb, 1.0)
        K_neural = 0.50 * K_attn + 0.50 * K_emb
    else:
        K_neural = K_attn

    np.fill_diagonal(K_neural, 1.0)
    return K_neural


def compute_attention_covariance_kernel(
    cross_attn: np.ndarray,
    mds_coords: Optional[np.ndarray] = None,
    taxa_repr: Optional[np.ndarray] = None
) -> np.ndarray:
    """Backward-compatible alias for compute_neural_covariance_kernel."""
    return compute_neural_covariance_kernel(cross_attn, taxa_repr=taxa_repr, mds_coords=mds_coords)



# =========================================================================
# 1. In-Frame Coding Alignment Validation
# =========================================================================

def verify_coding_alignment(
    seq_dict: Dict[str, str],
    allow_stop_codons: bool = True
) -> Tuple[int, int]:
    """
    Strictly verifies that sequences form an aligned, in-frame coding dataset.

    Requirements:
    - Non-empty alignment.
    - All sequences must possess identical aligned lengths.
    - Aligned sequence length must be a multiple of 3 (L_nt % 3 == 0).
    - Checks for internal stop codons (TAA, TAG, TGA).

    Returns:
        (num_taxa, num_codons)
    Raises:
        ValueError: If sequence lengths differ, are non-coding, or violate codon constraints.
    """
    if not seq_dict:
        raise ValueError("Provided alignment is empty.")

    taxa = list(seq_dict.keys())
    first_taxon = taxa[0]
    first_seq = seq_dict[first_taxon].upper().replace('U', 'T')
    l_nt = len(first_seq)

    if l_nt == 0:
        raise ValueError(f"Sequence '{first_taxon}' has length 0.")

    if l_nt % 3 != 0:
        raise ValueError(
            f"HyphAeon is a codon-level foundation model and strictly requires in-frame coding sequences. "
            f"Sequence '{first_taxon}' has length {l_nt} nt ({l_nt % 3} remainder modulo 3). "
            f"Please verify open reading frames and remove non-coding flanking regions or frameshift indels."
        )

    num_codons = l_nt // 3

    # Check length uniformity across all taxa
    mismatched = []
    for t in taxa:
        seq_len = len(seq_dict[t])
        if seq_len != l_nt:
            mismatched.append((t, seq_len))
            if len(mismatched) >= 5:
                break

    if mismatched:
        details = ", ".join([f"'{t}': {l} nt" for t, l in mismatched])
        raise ValueError(
            f"Alignment sequences are not uniformly aligned to length {l_nt} nt. "
            f"Mismatched examples: {details}."
        )

    # Stop codon inspection
    stop_codons = {'TAA', 'TAG', 'TGA'}
    taxa_with_stops = []
    total_stops = 0

    for t, s in seq_dict.items():
        s_clean = s.upper().replace('U', 'T')
        for c_idx in range(num_codons - 1):
            codon = s_clean[c_idx * 3:(c_idx + 1) * 3]
            if codon in stop_codons:
                total_stops += 1
                if t not in taxa_with_stops:
                    taxa_with_stops.append(t)

    if taxa_with_stops:
        pct_affected = 100.0 * len(taxa_with_stops) / len(taxa)
        msg = (
            f"[!] Notice: Detected {total_stops} internal stop codon(s) across "
            f"{len(taxa_with_stops)}/{len(taxa)} taxa ({pct_affected:.1f}%). "
            f"HyphAeon automatically tokenizes stop codons to token 64 ('*')."
        )
        if not allow_stop_codons:
            raise ValueError(f"{msg} Set allow_stop_codons=True to proceed anyway.")
        else:
            print(msg)

    return len(taxa), num_codons


def generate_consensus_sequence(seq_dict: Dict[str, str], taxa: Optional[List[str]] = None) -> str:
    """Computes the majority-rule nucleotide consensus sequence across specified taxa."""
    if taxa is None:
        taxa = list(seq_dict.keys())
    if not taxa:
        raise ValueError("Cannot compute consensus of 0 sequences.")

    seq_len = len(seq_dict[taxa[0]])
    consensus_chars = []
    for pos in range(seq_len):
        counts: Dict[str, int] = {}
        for t in taxa:
            char = seq_dict[t][pos].upper()
            if char not in ['-', '?', 'N']:
                counts[char] = counts.get(char, 0) + 1
        if counts:
            best_char = max(counts.items(), key=lambda x: x[1])[0]
        else:
            best_char = '-'
        consensus_chars.append(best_char)
    return "".join(consensus_chars)


def generate_time_decay_consensus_sequence(
    seq_dict: Dict[str, str],
    dates_map: Dict[str, float],
    taxa: Optional[List[str]] = None,
    gamma: Optional[float] = None,
    half_life: Optional[float] = None,
) -> Tuple[str, float]:
    """
    Computes a time-decay weighted nucleotide consensus sequence:
        w_i \propto \exp(-\gamma * (t_i - t_min))

    Downweights modern, densely-sampled contemporary isolates and upweights ancestral/early isolates,
    producing an optimal tree-free reference root anchor resistant to temporal sampling bias.

    Returns:
        (consensus_sequence, effective_gamma)
    """
    if taxa is None:
        taxa = list(seq_dict.keys())

    valid_taxa = [t for t in taxa if t in dates_map and not np.isnan(dates_map[t])]
    if not valid_taxa:
        # Fallback to unweighted consensus if dates are unavailable
        return generate_consensus_sequence(seq_dict, taxa), 0.0

    times = np.array([dates_map[t] for t in valid_taxa], dtype=np.float64)
    t_min = float(np.min(times))
    t_max = float(np.max(times))
    delta_t = t_max - t_min

    if half_life is not None and half_life > 0:
        eff_gamma = float(np.log(2.0) / half_life)
    elif gamma is not None and gamma > 0:
        eff_gamma = float(gamma)
    else:
        # Adaptive default: if span > 0, set gamma = 0.05 or adapt
        if delta_t > 0:
            eff_gamma = 0.05 if (0.05 * delta_t >= 1.0) else float(2.0 / delta_t)
        else:
            eff_gamma = 0.0

    weights = np.exp(-eff_gamma * (times - t_min))
    w_sum = np.sum(weights)
    if w_sum > 0:
        weights /= w_sum
    else:
        weights = np.ones(len(valid_taxa)) / len(valid_taxa)

    seq_len = len(seq_dict[valid_taxa[0]])
    consensus_chars = []
    for pos in range(seq_len):
        char_weights: Dict[str, float] = {}
        for i, t in enumerate(valid_taxa):
            c = seq_dict[t][pos].upper()
            if c not in ['-', '?', 'N']:
                char_weights[c] = char_weights.get(c, 0.0) + float(weights[i])
        if char_weights:
            best_c = max(char_weights.items(), key=lambda x: x[1])[0]
        else:
            best_c = '-'
        consensus_chars.append(best_c)

    return "".join(consensus_chars), eff_gamma


# =========================================================================
# 2. Timestamp Extraction & Parsing
# =========================================================================

def parse_header_timestamp(name: str) -> float:
    """
    Parses timestamps from FASTA sequence headers across common phylogenetic formats:
    - ISO full date: `2021-04-15` or `2021/04/15`
    - Decimal year: `2021.25`, `1959.5`
    - Year-month: `2021-04`
    - Delimiter prefixed: `|2021-04-15`, `_1986.5`, `/1993.5`
    - Bette Korber / HIV LANL format: e.g. `B86US.SFMHS18` -> 1986.5, `Z59ZR.ZHU` -> 1959.5
    - Longitudinal intra-host WPI: e.g. `16WPI` -> 16.0
    """
    if not name:
        return np.nan

    # 1. Check for archival 1959 ZR59 anchor
    if 'Z59' in name or 'ZR59' in name or '1959' in name:
        return 1959.5

    # 2. Try standard ISO/decimal extraction from temporal module
    d = extract_date_from_string(name)
    if not np.isnan(d):
        return d

    # 3. Korber HIV-1 isolate pattern: [Subtype][2-digit year][Country].[Strain]
    # Examples: B86US.SFMHS18, F93BE_VI850, C86ET.ETH2220, A85UG.U455, D84ZR.84ZR085
    m_korber = re.search(r'^[A-Za-z](\d{2})[A-Za-z]{2}[._]', name)
    if m_korber:
        yr_short = int(m_korber.group(1))
        full_yr = 1900 + yr_short if yr_short >= 30 else 2000 + yr_short
        return float(full_yr) + 0.5

    # 4. Nextstrain / LANL pipe format with year: Ref_B_FR_83_HXB2 or Ref_C_ET_86_ETH2220
    m_pipe = re.search(r'_(?:[A-Z]{2})_(\d{2})_', name)
    if m_pipe:
        yr_short = int(m_pipe.group(1))
        full_yr = 1900 + yr_short if yr_short >= 30 else 2000 + yr_short
        return float(full_yr) + 0.5

    # 5. Longitudinal intra-host WPI (Weeks Post Infection)
    m_wpi = re.search(r'(\d+(?:\.\d+)?)\s*(?:WPI|wpi)', name)
    if m_wpi:
        return float(m_wpi.group(1))

    # 6. Longitudinal DPI (Days Post Infection)
    m_dpi = re.search(r'(\d+(?:\.\d+)?)\s*(?:DPI|dpi)', name)
    if m_dpi:
        return float(m_dpi.group(1))

    return np.nan


def parse_sample_dates(
    taxa: List[str],
    dates_source: Optional[Union[str, Path, Dict[str, float]]] = None,
    date_col: Optional[str] = None,
    strain_col: Optional[str] = None,
    date_regex: Optional[str] = None
) -> Tuple[Dict[str, float], List[str]]:
    """
    Ingests and parses sampling dates for a list of taxa.

    Supports:
    - Pre-parsed dictionary {taxon: date}
    - Nextstrain Auspice JSON v2
    - CSV or TSV metadata table
    - FASTA header auto-extraction
    """
    dates_map: Dict[str, float] = {}

    # Case 1: Direct dictionary
    if isinstance(dates_source, dict):
        for t in taxa:
            if t in dates_source:
                val = parse_date_to_decimal(dates_source[t])
                if not np.isnan(val):
                    dates_map[t] = val

    # Case 2: External file (JSON, CSV, TSV)
    elif dates_source is not None:
        source_path = Path(dates_source)
        if not source_path.exists():
            raise FileNotFoundError(f"Dates source file not found: {source_path}")

        if source_path.suffix.lower() == '.json':
            with open(source_path, 'r', encoding='utf-8') as f:
                raw_json = json.load(f)
            if 'tree' in raw_json:
                auspice_dates = parse_dates_from_auspice_json(source_path)
                dates_map.update(auspice_dates)
            else:
                for k, v in raw_json.items():
                    val = parse_date_to_decimal(v)
                    if not np.isnan(val):
                        dates_map[k] = val
                    elif isinstance(v, dict) and 'year' in v:
                        dates_map[k] = float(v['year'])

        elif source_path.suffix.lower() in ['.csv', '.tsv', '.txt']:
            sep = '\t' if source_path.suffix.lower() in ['.tsv', '.txt'] else ','
            df = pd.read_csv(source_path, sep=sep)

            if not strain_col:
                cand_strains = ['strain', 'taxon', 'taxa', 'name', 'id', 'accession', 'sequence']
                for c in df.columns:
                    if c.lower() in cand_strains:
                        strain_col = c
                        break
                if not strain_col:
                    strain_col = df.columns[0]

            if not date_col:
                cand_dates = ['date', 'year', 'time', 'num_date', 'collection_date', 'sampling_date']
                for c in df.columns:
                    if c.lower() in cand_dates:
                        date_col = c
                        break
                if not date_col:
                    date_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

            for _, row in df.iterrows():
                strain_val = str(row[strain_col]).strip()
                date_val = parse_date_to_decimal(row[date_col])
                if not np.isnan(date_val):
                    dates_map[strain_val] = date_val

    # Case 3: Fallback to sequence header extraction
    for t in taxa:
        if t not in dates_map or np.isnan(dates_map[t]):
            if date_regex:
                m = re.search(date_regex, t)
                if m:
                    extracted = parse_date_to_decimal(m.group(1))
                    if not np.isnan(extracted):
                        dates_map[t] = extracted
            if t not in dates_map or np.isnan(dates_map[t]):
                val = parse_header_timestamp(t)
                if not np.isnan(val):
                    dates_map[t] = val

    missing_taxa = [t for t in taxa if t not in dates_map or np.isnan(dates_map[t])]
    return dates_map, missing_taxa


# =========================================================================
# 3. Root-to-Tip Divergence Calculation & Heuristic Rooting
# =========================================================================

def extract_tree_root_to_tip(
    tree_path: str,
    taxa: List[str],
    dates_map: Dict[str, float],
    root_taxon: Optional[str] = None,
    optimize_root: bool = True
) -> Tuple[Dict[str, float], str]:
    """
    Computes patristic root-to-tip distances from a Newick phylogenetic tree.
    If optimize_root=True, searches candidate rooting nodes to maximize
    the correlation R^2 with tip dates (TempEst / Path-O-Gen emulation).
    """
    if not HAS_BIOPHYLO:
        raise ImportError("Bio.Phylo is required to extract distances from phylogenetic trees.")

    tree = Phylo.read(tree_path, 'newick')
    taxa_set = set(taxa)

    # Explicit user rooting
    if root_taxon:
        matching_terminals = [t for t in tree.get_terminals() if t.name == root_taxon]
        if matching_terminals:
            tree.root_with_outgroup(matching_terminals[0])
            dists = {tip.name: tree.distance(tip) for tip in tree.get_terminals() if tip.name in taxa_set}
            return dists, f"user_root_{root_taxon}"

    # Heuristic Root Optimization (TempEst / Path-O-Gen emulation)
    if optimize_root:
        non_terminals = tree.get_nonterminals()
        best_r2 = -1.0
        best_node = None
        best_dists: Dict[str, float] = {}

        # Subsample candidate nodes for large trees to avoid N^2 tree traversals
        if len(non_terminals) > 60:
            step = max(1, len(non_terminals) // 60)
            candidate_nodes = non_terminals[::step]
        else:
            candidate_nodes = non_terminals

        for node in candidate_nodes:
            try:
                tree.root_with_outgroup(node)
            except Exception:
                continue

            dists = {tip.name: tree.distance(tip) for tip in tree.get_terminals() if tip.name in taxa_set}
            xs, ys = [], []
            for name, d in dists.items():
                if name in dates_map and not np.isnan(dates_map[name]):
                    xs.append(dates_map[name])
                    ys.append(d)

            if len(xs) >= 5 and np.std(ys) > 1e-7 and np.std(xs) > 1e-7:
                r_val = float(np.corrcoef(xs, ys)[0, 1])
                slope = float(np.polyfit(xs, ys, 1)[0])
                if slope > 0 and (r_val ** 2) > best_r2:
                    best_r2 = r_val ** 2
                    best_node = node
                    best_dists = dists

        if best_node is not None:
            tree.root_with_outgroup(best_node)
            return best_dists, "optimized_internal_root"

    # Default: Use original root
    dists = {tip.name: tree.distance(tip) for tip in tree.get_terminals() if tip.name in taxa_set}
    return dists, "original_tree_root"


def compute_tree_free_divergences(
    seq_dict: Dict[str, str],
    dated_taxa: List[str],
    dates_map: Dict[str, float],
    root_taxon: Optional[str] = None,
    decay_gamma: Optional[float] = None,
    decay_half_life: Optional[float] = None,
) -> Tuple[np.ndarray, str]:
    """
    Computes tree-free root-to-tip divergence directly from pairwise TN93 distances.
    Anchors root at:
    1. Specified root_taxon (if found in alignment).
    2. Explicit unweighted consensus if requested ('unweighted_consensus' / 'modal_consensus').
    3. Earliest sampled taxon / cohort if requested ('earliest' / 'earliest_cohort').
    4. Time-Decay Weighted Consensus (default & recommended for tree-free dating).
    """
    all_taxa = list(seq_dict.keys())

    # Case 1: User explicitly specified an existing taxon as root (e.g. outgroup or specific strain)
    if root_taxon and root_taxon in seq_dict:
        eval_taxa = [root_taxon] + [t for t in dated_taxa if t != root_taxon]
        dist_mat = compute_tn93_distance_matrix(seq_dict, eval_taxa)
        # Distance from root (index 0) to each dated taxon
        div_dict = {eval_taxa[i]: dist_mat[0, i] for i in range(1, len(eval_taxa))}
        divergences = np.array([div_dict[t] for t in dated_taxa if t != root_taxon], dtype=np.float64)
        return divergences, f"explicit_root_{root_taxon}"

    # Case 2: User requested unweighted modal consensus
    if root_taxon and root_taxon.lower() in ['unweighted_consensus', 'flat_consensus', 'modal_consensus']:
        con_seq = generate_consensus_sequence(seq_dict, dated_taxa)
        aug_dict = dict(seq_dict)
        aug_dict['__SYNTHETIC_CONSENSUS__'] = con_seq
        eval_taxa = ['__SYNTHETIC_CONSENSUS__'] + dated_taxa
        dist_mat = compute_tn93_distance_matrix(aug_dict, eval_taxa)
        divergences = np.array([dist_mat[0, i + 1] for i in range(len(dated_taxa))], dtype=np.float64)
        return divergences, "unweighted_modal_consensus_root"

    # Case 3: Anchor on earliest sampled cohort
    if root_taxon and root_taxon.lower() in ['earliest', 'earliest_taxon', 'earliest_cohort']:
        valid_dates = [(t, dates_map[t]) for t in dated_taxa if t in dates_map and not np.isnan(dates_map[t])]
        valid_dates.sort(key=lambda x: x[1])
        min_date = valid_dates[0][1]
        earliest_taxa = [t for t, d in valid_dates if abs(d - min_date) < 1e-4]

        dist_mat = compute_tn93_distance_matrix(seq_dict, dated_taxa)
        taxa_idx = {t: i for i, t in enumerate(dated_taxa)}
        earliest_indices = [taxa_idx[t] for t in earliest_taxa]

        if len(earliest_indices) == 1:
            root_idx = earliest_indices[0]
            divergences = dist_mat[root_idx, :].copy()
            root_desc = f"earliest_taxon_{dated_taxa[root_idx]}"
        else:
            divergences = np.mean(dist_mat[earliest_indices, :], axis=0)
            root_desc = f"earliest_cohort_n{len(earliest_indices)}"

        return divergences, root_desc

    # Case 4 (Default & Recommended for Tree-Free): Time-Decay Weighted Consensus
    decay_seq, eff_gamma = generate_time_decay_consensus_sequence(
        seq_dict, dates_map, dated_taxa, gamma=decay_gamma, half_life=decay_half_life
    )
    aug_dict = dict(seq_dict)
    aug_dict['__TIME_DECAY_ROOT__'] = decay_seq
    eval_taxa = ['__TIME_DECAY_ROOT__'] + dated_taxa
    dist_mat = compute_tn93_distance_matrix(aug_dict, eval_taxa)
    divergences = np.array([dist_mat[0, i + 1] for i in range(len(dated_taxa))], dtype=np.float64)
    root_desc = f"time_decay_consensus_root (γ={eff_gamma:.4f})"
    return divergences, root_desc


# =========================================================================
# 4. Dating Estimators: OLS, Attention PGLS, Latent Manifold Collapse
# =========================================================================

def run_ols_dating(
    times: np.ndarray,
    dists: np.ndarray,
    t_ref: Optional[float] = None,
    n_boot: int = 1000,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Fits Centered Root-to-Tip Ordinary Least Squares (OLS) regression:
        d_i = mu * (t_i - t_ref) + d_0 + epsilon_i

    Estimated ancestor date:
        t_MRCA = t_ref - (d_0 / mu)

    Reference centering (t_ref = mean(t)) guarantees orthogonal predictors,
    preventing delta-method covariance blowup when dates are in calendar years.
    """
    n = len(times)
    if n < 3:
        raise ValueError(f"At least 3 observations are required for OLS dating (got N={n}).")

    if t_ref is None:
        t_ref = float(np.mean(times))

    x = times - t_ref
    X = np.column_stack([x, np.ones(n)])

    beta_ols, residuals, rank, s = la.lstsq(X, dists)
    mu_ols = float(beta_ols[0])
    d0_ols = float(beta_ols[1])

    if abs(mu_ols) < 1e-15:
        mu_ols = 1e-15

    t_mrca = float(t_ref - (d0_ols / mu_ols))

    # Residual variance and covariance matrix
    res = dists - X @ beta_ols
    sigma2 = float(np.sum(res ** 2) / max(1, n - 2))
    cov_beta = sigma2 * la.inv(X.T @ X)

    se_mu = float(np.sqrt(max(1e-15, cov_beta[0, 0])))
    se_d0 = float(np.sqrt(max(1e-15, cov_beta[1, 1])))

    # Delta method for SE(t_MRCA)
    grad = np.array([d0_ols / (mu_ols ** 2), -1.0 / mu_ols])
    var_mrca = float(grad @ cov_beta @ grad)
    se_mrca = float(np.sqrt(max(0.0, var_mrca)))

    t_crit = float(stats.t.ppf(0.975, df=max(1, n - 2)))
    ci_analytical = [t_mrca - t_crit * se_mrca, t_mrca + t_crit * se_mrca]

    # Non-parametric Bootstrap for empirical 95% CI
    ci_bootstrap = None
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        boot_mrcas = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, size=n)
            X_b = X[idx]
            d_b = dists[idx]
            try:
                b_b, _, _, _ = la.lstsq(X_b, d_b)
                m_b, c_b = b_b[0], b_b[1]
                if abs(m_b) > 1e-8:
                    boot_mrcas.append(t_ref - (c_b / m_b))
            except Exception:
                continue
        if len(boot_mrcas) >= 50:
            ci_bootstrap = [
                float(np.percentile(boot_mrcas, 2.5)),
                float(np.percentile(boot_mrcas, 97.5))
            ]

    # Correlation and R^2
    r_val = float(np.corrcoef(times, dists)[0, 1]) if np.std(times) > 1e-8 and np.std(dists) > 1e-8 else 0.0
    r2 = r_val ** 2
    f_stat = (r2 / (1.0 - r2 + 1e-12)) * (n - 2) if r2 < 1.0 else 999.0
    p_val = float(1.0 - stats.f.cdf(f_stat, 1, max(1, n - 2)))

    return {
        'method': 'OLS',
        'mu': mu_ols,
        'd0': d0_ols,
        't_ref': t_ref,
        't_mrca': t_mrca,
        'se_mu': se_mu,
        'se_d0': se_d0,
        'se_mrca': se_mrca,
        'ci_analytical': ci_analytical,
        'ci_bootstrap': ci_bootstrap,
        'ci_mrca': ci_bootstrap if ci_bootstrap is not None else ci_analytical,
        'r': r_val,
        'r2': r2,
        'p_value': p_val,
        'sigma2': sigma2,
        'rmse': float(np.sqrt(np.mean(res ** 2))),
        'residuals': res,
        'fitted': X @ beta_ols,
        'n': n
    }


def run_pgls_dating(
    times: np.ndarray,
    dists: np.ndarray,
    cov_matrix: np.ndarray,
    ridge: float = 0.05,
    pagel_lambda: Optional[float] = None,
    t_ref: Optional[float] = None,
    n_boot: int = 1000,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Fits Centered Phylogenetic Generalized Least Squares (PGLS) regression:
        d = X * beta + epsilon,   Cov(epsilon) = sigma^2 * Sigma

    where Sigma is HyphAeon's neural phylogenetic covariance matrix.
    Directly incorporates evolutionary covariance without tree reconstruction.
    """
    n = len(times)
    if n < 3:
        raise ValueError(f"At least 3 observations are required for PGLS dating (got N={n}).")

    if t_ref is None:
        t_ref = float(np.mean(times))

    x = times - t_ref
    X = np.column_stack([x, np.ones(n)])

    w_raw, v = la.eigh(cov_matrix)
    w_pos = np.maximum(w_raw, 0.0)

    # Either Pagel's lambda covariance: C = lambda * K + (1 - lambda) * I
    # or additive ridge covariance: C = K + ridge * I
    if pagel_lambda is not None:
        eff_lam = float(np.clip(pagel_lambda, 0.001, 0.999))
        w_c = eff_lam * w_pos + (1.0 - eff_lam)
    else:
        w_c = w_pos + ridge

    inv_w = 1.0 / np.maximum(w_c, 1e-12)
    C_inv = v @ np.diag(inv_w) @ v.T

    # GLS solution: beta = (X^T C^-1 X)^-1 X^T C^-1 d
    Xt_Cinv = X.T @ C_inv
    Xt_Cinv_X = Xt_Cinv @ X
    beta_gls = la.solve(Xt_Cinv_X, Xt_Cinv @ dists)

    mu_gls = float(beta_gls[0])
    d0_gls = float(beta_gls[1])

    if abs(mu_gls) < 1e-15:
        mu_gls = 1e-15

    t_mrca = float(t_ref - (d0_gls / mu_gls))

    residuals = dists - X @ beta_gls
    sigma2_gls = float((residuals.T @ C_inv @ residuals) / max(1, n - 2))
    cov_beta = sigma2_gls * la.inv(Xt_Cinv_X)

    se_mu = float(np.sqrt(max(1e-15, cov_beta[0, 0])))
    se_d0 = float(np.sqrt(max(1e-15, cov_beta[1, 1])))

    # Delta method for SE(t_MRCA)
    grad = np.array([d0_gls / (mu_gls ** 2), -1.0 / mu_gls])
    var_mrca = float(grad @ cov_beta @ grad)
    se_mrca = float(np.sqrt(max(0.0, var_mrca)))

    t_crit = float(stats.t.ppf(0.975, df=max(1, n - 2)))
    ci_lower = float(t_mrca - t_crit * se_mrca)
    ci_upper = min(float(np.min(times)), float(t_mrca + t_crit * se_mrca))
    ci_analytical = [ci_lower, ci_upper]

    # Generalized R^2 (Buse 1973)
    one_Cinv_one = float(np.ones(n).T @ C_inv @ np.ones(n))
    weighted_mean = float(np.ones(n).T @ C_inv @ dists) / max(1e-12, one_Cinv_one)
    tot_residuals = dists - weighted_mean
    ss_tot = float(tot_residuals.T @ C_inv @ tot_residuals)
    ss_res = float(residuals.T @ C_inv @ residuals)
    r2_gls = float(max(0.0, 1.0 - (ss_res / max(1e-12, ss_tot))))

    return {
        'method': 'PGLS',
        'mu': mu_gls,
        'd0': d0_gls,
        't_ref': t_ref,
        't_mrca': t_mrca,
        'se_mu': se_mu,
        'se_d0': se_d0,
        'se_mrca': se_mrca,
        'ci_analytical': ci_analytical,
        'ci_mrca': ci_analytical,
        'r2': r2_gls,
        'ridge': ridge if pagel_lambda is None else (1.0 - eff_lam),
        'pagel_lambda': pagel_lambda,
        'sigma2': sigma2_gls,
        'rmse': float(np.sqrt(np.mean(residuals ** 2))),
        'residuals': residuals,
        'fitted': X @ beta_gls,
        'times': times,
        'n': n
    }


def estimate_reml_pagel_lambda(
    times: np.ndarray,
    dists: np.ndarray,
    cov_matrix: np.ndarray
) -> Dict[str, Any]:
    """
    Estimates the phylogenetic signal / shrinkage parameter lambda in [0, 1] (Pagel's lambda)
    by maximizing the exact profile Restricted Maximum Likelihood (REML).

    Covariance model:
        C(lambda) = lambda * cov_matrix + (1 - lambda) * I

    where lambda = 1 represents full neural phylogenetic covariance and lambda = 0 represents
    independent tip variance. Computed in O(N) using spectral projection.
    """
    n = len(times)
    t_ref = float(np.mean(times))
    x = times - t_ref
    X = np.column_stack([x, np.ones(n)])

    w_K, V = la.eigh(cov_matrix)
    w_K = np.maximum(w_K, 0.0)

    # Pre-project design matrix and responses onto eigenvectors
    Z = V.T @ X      # (N, 2)
    u = V.T @ dists  # (N,)

    def neg_reml(lam):
        w_c = lam * w_K + (1.0 - lam)
        inv_w = 1.0 / np.maximum(w_c, 1e-12)

        Z_scaled = Z * inv_w[:, None]
        Xt_Cinv_X = Z.T @ Z_scaled
        Xt_Cinv_d = Z_scaled.T @ u

        try:
            beta = la.solve(Xt_Cinv_X, Xt_Cinv_d)
        except Exception:
            return 1e9

        res_ss = float(np.sum((u ** 2) * inv_w) - beta.T @ Xt_Cinv_d)
        if res_ss <= 0:
            return 1e9

        sigma2 = res_ss / max(1, n - 2)
        log_det_C = float(np.sum(np.log(np.maximum(w_c, 1e-12))))
        _, log_det_XtCinvX = np.linalg.slogdet(Xt_Cinv_X)

        minus_2_logL = (n - 2) * np.log(sigma2) + log_det_C + log_det_XtCinvX
        return minus_2_logL

    res_opt = optimize.minimize_scalar(neg_reml, bounds=(0.001, 0.999), method='bounded')
    opt_lambda = float(res_opt.x) if res_opt.success else 0.95

    return {
        'best_lambda': opt_lambda,
        'status': 'OPTIMAL_REML'
    }


def tune_ridge_for_pgls(
    times: np.ndarray,
    dists: np.ndarray,
    cov_matrix: np.ndarray,
    method: str = 'reml'
) -> Dict[str, Any]:
    """Backward-compatible wrapper mapping to estimate_reml_pagel_lambda."""
    reml_res = estimate_reml_pagel_lambda(times, dists, cov_matrix)
    return {
        'best_lambda': reml_res['best_lambda'],
        'status': reml_res['status']
    }


def run_powerlaw_clock_dating(
    times: np.ndarray,
    dists: np.ndarray,
    cov_matrix: Optional[np.ndarray] = None,
    ridge: float = 0.05,
    n_boot: int = 500,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Fits a Time-Dependent Rate (TDR) Power-Law Molecular Clock:
        d(t) = k * (t - t_MRCA)^theta,   t > t_MRCA

    Captures sub-linear rate deceleration (theta < 1.0) caused by long-term purifying
    selection and mutational saturation (Aiewsakun & Katzourakis 2015, Membrebe et al. 2019).
    When theta == 1.0, reduces to standard linear clock regression d(t) = mu * (t - t_MRCA).

    Performs nested F-test and AIC comparison against the linear null model.
    """
    n = len(times)
    if n < 4:
        raise ValueError(f"At least 4 observations are required for non-linear power-law dating (got N={n}).")

    t_min = float(np.min(times))
    t_max = float(np.max(times))
    delta_t = max(1e-4, t_max - t_min)

    # Covariance weighting with single-call spectral projection & inversion:
    if cov_matrix is not None:
        w_raw, v = la.eigh(cov_matrix)
        w_c = np.maximum(w_raw, 0.0) + ridge
        C_inv = v @ np.diag(1.0 / w_c) @ v.T
    else:
        C_inv = np.eye(n)

    # 1. Fit Linear Null Model
    x_mean = float(np.mean(times))
    d_mean = float(np.mean(dists))
    X_lin = np.column_stack([times - x_mean, np.ones(n)])
    Xt_Cinv = X_lin.T @ C_inv
    beta_lin = la.solve(Xt_Cinv @ X_lin, Xt_Cinv @ dists)
    mu_lin = float(beta_lin[0])
    d0_lin = float(beta_lin[1])
    fitted_lin = X_lin @ beta_lin
    res_lin = dists - fitted_lin
    rss_lin = float(res_lin.T @ C_inv @ res_lin)
    t0_lin = float(x_mean - d0_lin / max(1e-12, mu_lin))
    aic_lin = float(n * np.log(max(1e-12, rss_lin / n)) + 2 * 2)

    # 2. Fit Power-Law Model d(t) = k * (t - t0)^theta
    def objective(params):
        t0, k, theta = params
        if t0 >= t_min - 0.001:
            return 1e9 + float((t0 - t_min) ** 2)
        dt = times - t0
        if np.any(dt <= 0.0) or k <= 0.0 or theta <= 0.0:
            return 1e9
        pred = k * (dt ** theta)
        res = dists - pred
        return float(res.T @ C_inv @ res)

    best_res = None
    best_val = 1e12

    # Multi-start initializations over t0 and theta
    t0_candidates = [
        t_min - 0.1 * delta_t,
        t_min - 0.3 * delta_t,
        t_min - 0.6 * delta_t,
        t_min - 1.2 * delta_t,
        t_min - 2.5 * delta_t,
        t_min - 5.0 * delta_t,
    ]
    if t0_lin < t_min:
        t0_candidates.append(t0_lin)

    for t0_cand in t0_candidates:
        if t0_cand >= t_min:
            continue
        for th_cand in [0.4, 0.7, 1.0, 1.3]:
            dt_mean = max(1e-4, x_mean - t0_cand)
            k_cand = float(max(1e-6, d_mean / (dt_mean ** th_cand)))
            try:
                opt = optimize.minimize(
                    objective,
                    [t0_cand, k_cand, th_cand],
                    bounds=[
                        (t_min - 50.0 * delta_t, t_min - 0.001),
                        (1e-7, 100.0),
                        (0.05, 3.0)
                    ],
                    method='L-BFGS-B'
                )
                if opt.fun < best_val:
                    best_val = opt.fun
                    best_res = opt
            except Exception:
                continue

    if best_res is not None and best_res.success:
        t0_nl, k_nl, th_nl = [float(x) for x in best_res.x]
        rss_nl = float(best_val)
    else:
        # Fallback to linear
        t0_nl = t0_lin
        th_nl = 1.0
        k_nl = mu_lin
        rss_nl = rss_lin

    pred_nl = k_nl * (np.maximum(1e-6, times - t0_nl) ** th_nl)
    res_nl = dists - pred_nl

    # 3. Model Comparison Metrics (F-test and AIC)
    df_lin = n - 2
    df_nl = n - 3
    diff_rss = max(0.0, rss_lin - rss_nl)
    f_stat = float((diff_rss / 1.0) / max(1e-12, rss_nl / max(1, df_nl)))
    p_f_test = float(stats.f.sf(f_stat, 1, max(1, df_nl)))

    aic_nl = float(n * np.log(max(1e-12, rss_nl / n)) + 2 * 3)
    delta_aic = float(aic_lin - aic_nl)

    # Automatic selection decision:
    # Requires p < 0.05, delta_AIC >= 2.0, and meaningful curvature deviation (|theta - 1.0| >= 0.03)
    is_nonlinear_preferred = bool(
        (p_f_test < 0.05) and (delta_aic >= 2.0) and (abs(th_nl - 1.0) >= 0.03)
    )

    # Instantaneous Rates
    r_ancestral = float(k_nl * th_nl * (max(1e-6, t_min - t0_nl) ** (th_nl - 1.0)))
    r_recent = float(k_nl * th_nl * (max(1e-6, t_max - t0_nl) ** (th_nl - 1.0)))
    r_mean = float((k_nl * (max(1e-6, t_max - t0_nl) ** th_nl) - k_nl * (max(1e-6, t_min - t0_nl) ** th_nl)) / delta_t)

    # Bootstrap CIs for non-linear parameters
    rng = np.random.RandomState(seed)
    boot_t0 = []
    boot_th = []
    boot_k = []
    if n_boot > 0:
        for _ in range(n_boot):
            b_idx = rng.choice(n, size=n, replace=True)
            b_t = times[b_idx]
            b_d = dists[b_idx]
            if len(np.unique(b_t)) < 3:
                continue
            b_tmin = float(np.min(b_t))
            b_dt = max(1e-4, float(np.max(b_t)) - b_tmin)
            try:
                def b_obj(p):
                    if p[0] >= b_tmin - 0.001:
                        return 1e9
                    dt = b_t - p[0]
                    if np.any(dt <= 0) or p[1] <= 0 or p[2] <= 0:
                        return 1e9
                    return float(np.sum((b_d - p[1] * (dt ** p[2])) ** 2))

                b_opt = optimize.minimize(
                    b_obj, [t0_nl, k_nl, th_nl],
                    bounds=[(b_tmin - 50.0 * b_dt, b_tmin - 0.001), (1e-7, 100.0), (0.05, 3.0)],
                    method='L-BFGS-B'
                )
                if b_opt.success:
                    boot_t0.append(float(b_opt.x[0]))
                    boot_k.append(float(b_opt.x[1]))
                    boot_th.append(float(b_opt.x[2]))
            except Exception:
                pass

    if len(boot_t0) >= 20:
        ci_t0 = [float(np.percentile(boot_t0, 2.5)), float(np.percentile(boot_t0, 97.5))]
        ci_th = [float(np.percentile(boot_th, 2.5)), float(np.percentile(boot_th, 97.5))]
        ci_k = [float(np.percentile(boot_k, 2.5)), float(np.percentile(boot_k, 97.5))]
    else:
        ci_t0 = [t0_nl, t0_nl]
        ci_th = [th_nl, th_nl]
        ci_k = [k_nl, k_nl]

    # Generalized R^2
    one_Cinv_one = float(np.ones(n).T @ C_inv @ np.ones(n))
    weighted_mean = float(np.ones(n).T @ C_inv @ dists) / max(1e-12, one_Cinv_one)
    tot_residuals = dists - weighted_mean
    ss_tot = float(tot_residuals.T @ C_inv @ tot_residuals)
    r2_nl = float(max(0.0, 1.0 - (rss_nl / max(1e-12, ss_tot))))

    return {
        'method': 'POWER_LAW',
        't_mrca': t0_nl,
        'ci_mrca': ci_t0,
        'theta': th_nl,
        'ci_theta': ci_th,
        'k': k_nl,
        'ci_k': ci_k,
        'rate_recent': r_recent,
        'rate_ancestral': r_ancestral,
        'rate_mean': r_mean,
        'rss': rss_nl,
        'aic': aic_nl,
        'rss_linear': rss_lin,
        'aic_linear': aic_lin,
        'delta_aic': delta_aic,
        'f_stat': f_stat,
        'p_f_test': p_f_test,
        'is_nonlinear_preferred': is_nonlinear_preferred,
        'r2': r2_nl,
        'rmse': float(np.sqrt(np.mean(res_nl ** 2))),
        'fitted': pred_nl,
        'residuals': res_nl,
        'n': n
    }


def compute_rcs_basis(
    x: np.ndarray,
    knots: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes Harrell's Restricted Cubic Spline (RCS) basis matrix and its first derivative.
    For k knots, produces (k - 2) non-linear columns.
    When knots[0] = t_min, all non-linear basis columns and derivatives
    vanish identically for t <= t_min, guaranteeing strictly linear ancestral extrapolation.
    """
    knots = np.asarray(knots, dtype=float)
    k = len(knots)
    if k < 3:
        raise ValueError("At least 3 knots are required for restricted cubic splines.")
    t1, tk_1, tk = knots[0], knots[-2], knots[-1]
    denom = (tk - t1) ** 2

    cols, d_cols = [], []
    for j in range(k - 2):
        tj = knots[j]
        term1 = np.maximum(0.0, x - tj) ** 3
        term2 = ((tk - tj) / (tk - tk_1)) * (np.maximum(0.0, x - tk_1) ** 3)
        term3 = ((tk_1 - tj) / (tk - tk_1)) * (np.maximum(0.0, x - tk) ** 3)
        cols.append((term1 - term2 + term3) / denom)

        d1 = 3.0 * (np.maximum(0.0, x - tj) ** 2)
        d2 = ((tk - tj) / (tk - tk_1)) * (3.0 * (np.maximum(0.0, x - tk_1) ** 2))
        d3 = ((tk_1 - tj) / (tk - tk_1)) * (3.0 * (np.maximum(0.0, x - tk) ** 2))
        d_cols.append((d1 - d2 + d3) / denom)

    B = np.column_stack(cols) if cols else np.empty((len(x), 0))
    dB = np.column_stack(d_cols) if d_cols else np.empty((len(x), 0))
    return B, dB


def run_restricted_spline_clock_dating(
    times: np.ndarray,
    dists: np.ndarray,
    cov_matrix: Optional[np.ndarray] = None,
    ridge: float = 0.05,
    n_boot: int = 500,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Fits a Degrees-of-Freedom (DF) Restricted Natural Cubic Spline Molecular Clock.

    Solves the boundary collapse pathology of unconstrained power-law models by
    restricting the curve to be STRICTLY LINEAR beyond the boundary knots:
        For t <= t_min: d(t) = beta_0 + beta_1 * t

    Guarantees well-behaved, singularity-free ancestral extrapolation to the MRCA:
        t_MRCA = -beta_0 / beta_1

    While allowing non-linear curvature (flexibility) across the observation window [t_min, t_max]:
        d(t) = beta_0 + beta_1 * t + beta_2 * X_2(t)

    Performs an exact nested F-test and Delta-AIC comparison against the linear null model.
    """
    n = len(times)
    if n < 5:
        raise ValueError(f"At least 5 observations are required for restricted spline dating (got N={n}).")

    t_min = float(np.min(times))
    t_max = float(np.max(times))

    # 3 knots: t_min, median, 90th percentile
    knots = np.array([t_min, float(np.median(times)), float(np.percentile(times, 90))])
    B, dB = compute_rcs_basis(times, knots)

    # Covariance weighting with single-call spectral projection & inversion:
    if cov_matrix is not None:
        w_raw, v = la.eigh(cov_matrix)
        w_c = np.maximum(w_raw, 0.0) + ridge
        C_inv = v @ np.diag(1.0 / w_c) @ v.T
    else:
        C_inv = np.eye(n)

    # 1. Fit Linear Null Model: d(t) = beta_0 + beta_1 * t
    X_lin = np.column_stack([np.ones(n), times])
    Xt_Cinv_lin = X_lin.T @ C_inv
    beta_lin = la.solve(Xt_Cinv_lin @ X_lin, Xt_Cinv_lin @ dists)
    pred_lin = X_lin @ beta_lin
    res_lin = dists - pred_lin
    rss_lin = float(res_lin.T @ C_inv @ res_lin)
    aic_lin = float(n * np.log(max(1e-12, rss_lin / n)) + 2 * 2)

    # 2. Fit Restricted Spline Model: d(t) = beta_0 + beta_1 * t + beta_2 * X_2(t)
    X_sp = np.column_stack([np.ones(n), times, B])
    Xt_Cinv_sp = X_sp.T @ C_inv
    beta_sp = la.solve(Xt_Cinv_sp @ X_sp, Xt_Cinv_sp @ dists)
    pred_sp = X_sp @ beta_sp
    res_sp = dists - pred_sp
    rss_sp = float(res_sp.T @ C_inv @ res_sp)
    aic_sp = float(n * np.log(max(1e-12, rss_sp / n)) + 2 * 3)
    delta_aic = float(aic_lin - aic_sp)

    # 3. Model Comparison Metrics (Nested F-test & AIC)
    df_sp = n - 3
    diff_rss = max(0.0, rss_lin - rss_sp)
    f_stat = float((diff_rss / 1.0) / max(1e-12, rss_sp / max(1, df_sp)))
    p_f_test = float(stats.f.sf(f_stat, 1, max(1, df_sp)))

    # Ancestral MRCA (strictly linear for t <= t_min)
    mu_ancestral = float(beta_sp[1])
    if mu_ancestral > 1e-9:
        t0_sp = float(-beta_sp[0] / mu_ancestral)
    else:
        # If ancestral slope is non-positive, backward extrapolation is undefined
        t0_sp = float(-beta_lin[0] / max(1e-12, beta_lin[1])) if beta_lin[1] > 1e-9 else float('nan')

    # Instantaneous Rates
    dB_max = float(dB[-1, 0]) if len(dB) > 0 else 0.0
    mu_recent = float(beta_sp[1] + beta_sp[2] * dB_max)
    rate_ratio = float(mu_recent / mu_ancestral) if mu_ancestral > 1e-9 else 1.0

    # Automatic selection rule:
    # Requires statistical significance (p < 0.05), positive model evidence (delta_AIC >= 2.0),
    # positive ancestral rate, and biologically meaningful rate variation (|rate_ratio - 1.0| >= 0.15).
    is_nonlinear_preferred = bool(
        p_f_test < 0.05 and delta_aic >= 2.0 and mu_ancestral > 0 and abs(rate_ratio - 1.0) >= 0.15
    )

    # Bootstrap 95% Confidence Intervals
    boot_t0 = []
    boot_mu_anc = []
    boot_mu_rec = []
    boot_beta2 = []
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        for _ in range(n_boot):
            b_idx = rng.choice(n, size=n, replace=True)
            b_t = times[b_idx]
            b_d = dists[b_idx]
            if len(np.unique(b_t)) < 4:
                continue
            b_B, b_dB = compute_rcs_basis(b_t, knots)
            b_X = np.column_stack([np.ones(n), b_t, b_B])
            try:
                b_beta = la.lstsq(b_X, b_d, rcond=None)[0]
                if b_beta[1] > 1e-9:
                    boot_t0.append(float(-b_beta[0] / b_beta[1]))
                    boot_mu_anc.append(float(b_beta[1]))
                    b_rec = float(b_beta[1] + b_beta[2] * (b_dB[-1, 0] if len(b_dB) > 0 else 0.0))
                    boot_mu_rec.append(b_rec)
                    boot_beta2.append(float(b_beta[2]))
            except Exception:
                pass

    if len(boot_t0) >= 20:
        ci_t0 = [float(np.percentile(boot_t0, 2.5)), float(np.percentile(boot_t0, 97.5))]
        ci_mu_anc = [float(np.percentile(boot_mu_anc, 2.5)), float(np.percentile(boot_mu_anc, 97.5))]
        ci_mu_rec = [float(np.percentile(boot_mu_rec, 2.5)), float(np.percentile(boot_mu_rec, 97.5))]
        ci_beta2 = [float(np.percentile(boot_beta2, 2.5)), float(np.percentile(boot_beta2, 97.5))]
    else:
        ci_t0 = [t0_sp, t0_sp]
        ci_mu_anc = [mu_ancestral, mu_ancestral]
        ci_mu_rec = [mu_recent, mu_recent]
        ci_beta2 = [float(beta_sp[2]), float(beta_sp[2])]

    # Generalized R^2
    one_Cinv_one = float(np.ones(n).T @ C_inv @ np.ones(n))
    weighted_mean = float(np.ones(n).T @ C_inv @ dists) / max(1e-12, one_Cinv_one)
    tot_residuals = dists - weighted_mean
    ss_tot = float(tot_residuals.T @ C_inv @ tot_residuals)
    r2_sp = float(max(0.0, 1.0 - (rss_sp / max(1e-12, ss_tot))))

    return {
        'method': 'RESTRICTED_SPLINE',
        't_mrca': t0_sp,
        'ci_mrca': ci_t0,
        'rate_ancestral': mu_ancestral,
        'ci_rate_ancestral': ci_mu_anc,
        'rate_recent': mu_recent,
        'ci_rate_recent': ci_mu_rec,
        'rate_ratio': rate_ratio,
        'beta_0': float(beta_sp[0]),
        'beta_1': float(beta_sp[1]),
        'beta_2': float(beta_sp[2]),
        'ci_beta_2': ci_beta2,
        'knots': knots.tolist(),
        'rss': rss_sp,
        'aic': aic_sp,
        'rss_linear': rss_lin,
        'aic_linear': aic_lin,
        'delta_aic': delta_aic,
        'f_stat': f_stat,
        'p_f_test': p_f_test,
        'is_nonlinear_preferred': is_nonlinear_preferred,
        'r2': r2_sp,
        'rmse': float(np.sqrt(np.mean(res_sp ** 2))),
        'fitted': pred_sp,
        'residuals': res_sp,
        'n': n
    }


# =========================================================================
# 5. Diagnostic Visualization & Figure Generation
# =========================================================================

def plot_mrca_dating(
    dating_results: Dict[str, Any],
    output_path: Union[str, Path],
    title: Optional[str] = None
):
    """Generates a publication-grade diagnostic PDF and PNG figure."""
    if not HAS_MATPLOTLIB:
        print("[!] Matplotlib not available; skipping diagnostic plot.")
        return

    fig = plt.figure(figsize=(13, 5.5), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1.0])

    times = dating_results['times']
    dists = dating_results['dists']
    ols = dating_results['ols']
    pgls = dating_results.get('pgls')
    spline = dating_results.get('spline')
    power = dating_results.get('power')

    # Panel A: Root-to-Tip Molecular Clock Regression
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(times, dists, color='#1f77b4', s=42, alpha=0.75, edgecolors='black', linewidth=0.5, label='Dated Strains', zorder=3)

    t_min = float(np.min(times))
    t_max = float(np.max(times))
    mrca_candidates = [ols['t_mrca']]
    if pgls:
        mrca_candidates.append(pgls['t_mrca'])
    if spline:
        mrca_candidates.append(spline['t_mrca'])
    if power:
        mrca_candidates.append(power['t_mrca'])

    plot_left = min(t_min - (t_max - t_min) * 0.35, min(mrca_candidates) - (t_max - t_min) * 0.1)
    plot_right = t_max + (t_max - t_min) * 0.05
    x_grid = np.linspace(plot_left, plot_right, 200)

    # OLS fitted line
    y_ols = ols['mu'] * (x_grid - ols['t_ref']) + ols['d0']
    ax1.plot(x_grid, y_ols, color='#e63946', linestyle='--', linewidth=2.0,
             label=f"OLS (t_MRCA={ols['t_mrca']:.1f}, μ={ols['mu']:.5f})", zorder=4)

    # PGLS fitted line
    if pgls:
        y_pgls = pgls['mu'] * (x_grid - pgls['t_ref']) + pgls['d0']
        ax1.plot(x_grid, y_pgls, color='#7b2cbf', linestyle='-', linewidth=2.5,
                 label=f"HyphAeon PGLS (t_MRCA={pgls['t_mrca']:.1f}, μ={pgls['mu']:.5f})", zorder=5)

    # Restricted Spline fitted curve
    if spline:
        knots_arr = np.array(spline['knots'])
        b0 = spline['beta_0']
        b1 = spline['beta_1']
        b2 = spline['beta_2']
        x_spline = np.linspace(max(plot_left, spline['t_mrca']), plot_right, 250)
        B_grid, _ = compute_rcs_basis(x_spline, knots_arr)
        y_spline = b0 + b1 * x_spline + (b2 * B_grid[:, 0] if B_grid.shape[1] > 0 else 0.0)
        lbl_spline = f"Restricted Spline (t_MRCA={spline['t_mrca']:.1f}, μ_anc={spline['rate_ancestral']:.5f})"
        color_s = '#2a9d8f' if spline.get('is_nonlinear_preferred') else '#f4a261'
        style_s = '-' if spline.get('is_nonlinear_preferred') else ':'
        ax1.plot(x_spline, y_spline, color=color_s, linestyle=style_s, linewidth=2.4, label=lbl_spline, zorder=6)

    # Power-law fitted curve
    if power:
        x_power = np.linspace(max(plot_left, power['t_mrca'] + 1e-4), plot_right, 200)
        y_power = power['k'] * (np.maximum(0.0, x_power - power['t_mrca']) ** power['theta'])
        lbl_power = f"Power-Law (t_MRCA={power['t_mrca']:.1f}, θ={power['theta']:.3f})"
        color_p = '#2a9d8f' if power.get('is_nonlinear_preferred') else '#f4a261'
        style_p = '-' if power.get('is_nonlinear_preferred') else ':'
        ax1.plot(x_power, y_power, color=color_p, linestyle=style_p, linewidth=2.2, label=lbl_power, zorder=6)

    # MRCA markers and CI error bars at distance = 0
    ax1.axhline(0, color='gray', linestyle=':', linewidth=0.8, zorder=1)
    e_ols_l = max(0.0, ols['t_mrca'] - min(ols['ci_mrca'][0], ols['ci_mrca'][1]))
    e_ols_r = max(0.0, max(ols['ci_mrca'][0], ols['ci_mrca'][1]) - ols['t_mrca'])
    ax1.errorbar([ols['t_mrca']], [0], xerr=[[e_ols_l], [e_ols_r]],
                 fmt='s', color='#e63946', markersize=6, capsize=4, capthick=1.5, zorder=6)
    if pgls:
        e_pgls_l = max(0.0, pgls['t_mrca'] - min(pgls['ci_mrca'][0], pgls['ci_mrca'][1]))
        e_pgls_r = max(0.0, max(pgls['ci_mrca'][0], pgls['ci_mrca'][1]) - pgls['t_mrca'])
        ax1.errorbar([pgls['t_mrca']], [-0.002], xerr=[[e_pgls_l], [e_pgls_r]],
                     fmt='D', color='#7b2cbf', markersize=6, capsize=4, capthick=1.5, zorder=6)
    if spline and spline.get('is_nonlinear_preferred'):
        e_spl_l = max(0.0, spline['t_mrca'] - min(spline['ci_mrca'][0], spline['ci_mrca'][1]))
        e_spl_r = max(0.0, max(spline['ci_mrca'][0], spline['ci_mrca'][1]) - spline['t_mrca'])
        ax1.errorbar([spline['t_mrca']], [-0.003], xerr=[[e_spl_l], [e_spl_r]],
                     fmt='^', color='#2a9d8f', markersize=6, capsize=4, capthick=1.5, zorder=6)

    ax1.set_xlim(plot_left, plot_right)
    ax1.set_xlabel("Sampling Date / Time", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Root-to-Tip Divergence (subs/site)", fontsize=11, fontweight='bold')
    ax1.set_title("(A) Heterochronous Molecular Clock Regression", fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left', frameon=True, fontsize=8.5)
    ax1.grid(True, linestyle=':', alpha=0.4)

    # Panel B: Residual Error Diagnostics
    ax2 = fig.add_subplot(gs[1])
    res_ols = ols['residuals']
    ax2.scatter(times, res_ols, color='#e63946', alpha=0.7, s=40, edgecolors='black', linewidth=0.5, label='OLS Residuals', zorder=3)
    if pgls:
        res_pgls = pgls['residuals']
        t_pgls = pgls.get('times', times)
        if len(t_pgls) != len(res_pgls):
            t_pgls = times[:len(res_pgls)]
        ax2.scatter(t_pgls, res_pgls, color='#7b2cbf', alpha=0.7, s=40, marker='^', edgecolors='black', linewidth=0.5, label='HyphAeon PGLS Residuals', zorder=4)
    ax2.axhline(0, color='black', linestyle='--', linewidth=1.2)
    ax2.set_xlabel("Sampling Date / Time", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Residual Divergence (d - d_pred)", fontsize=11, fontweight='bold')
    ax2.set_title("(B) Residual Error Diagnostics", fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', frameon=True, fontsize=8.5)
    ax2.grid(True, linestyle=':', alpha=0.4)

    plt.suptitle(title or "HyphAeon Molecular Clock Calibration & Ancestor Dating", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()

    out_p = Path(output_path)
    ensure_parent_directory(out_p)
    plt.savefig(out_p, dpi=300)
    if out_p.suffix.lower() != '.png':
        plt.savefig(out_p.with_suffix('.png'), dpi=300)
    plt.close(fig)
    print(f"[✓] Diagnostic plot generated: {out_p}")


# =========================================================================
# 6. Master MRCA Dating Pipeline
# =========================================================================

def run_mrca_dating(
    alignment_path: Union[str, Path],
    tree_path: Optional[Union[str, Path]] = None,
    dates_source: Optional[Union[str, Path, Dict[str, float]]] = None,
    date_col: Optional[str] = None,
    strain_col: Optional[str] = None,
    date_regex: Optional[str] = None,
    root_taxon: Optional[str] = None,
    decay_gamma: Optional[float] = None,
    decay_half_life: Optional[float] = None,
    optimize_root: bool = True,
    use_tn93: bool = False,
    method: str = "all",
    clock_model: str = "auto",
    ridge: Union[float, str] = "auto",
    tune_ridge: bool = False,
    n_bootstrap: int = 1000,
    model: Optional[torch.nn.Module] = None,
    weights: Optional[str] = None,
    variant: Optional[str] = None,
    device: Optional[torch.device] = None,
    batch_size: Optional[int] = None,
    max_species: Optional[int] = None,
    allow_stop_codons: bool = True,
    output_prefix: Optional[str] = None,
    plot: bool = False
) -> Dict[str, Any]:
    """
    Executes end-to-end molecular clock calibration and MRCA dating on time-stamped sequences.

    Returns:
        Structured dictionary with model parameters, confidence intervals,
        per-taxon predictions, and model diagnostics.
    """
    t0 = time.time()
    align_p = Path(alignment_path)
    if not align_p.exists():
        raise FileNotFoundError(f"Alignment file not found: {align_p}")

    # 1. Parse and Validate In-Frame Coding Alignment
    seq_dict = parse_alignment_sequences(str(align_p))
    n_taxa_raw, n_codons = verify_coding_alignment(seq_dict, allow_stop_codons=allow_stop_codons)
    print(f"[*] Alignment verified: {n_taxa_raw} taxa, {n_codons} codons ({n_codons * 3} nt in-frame).")

    # 2. Extract and Parse Sampling Timestamps
    taxa_all = list(seq_dict.keys())
    dates_map, missing = parse_sample_dates(
        taxa_all,
        dates_source=dates_source,
        date_col=date_col,
        strain_col=strain_col,
        date_regex=date_regex
    )
    
    # Filter to dated taxa (excluding root_taxon if root_taxon is an ancestral reference without timestamp)
    dated_taxa = [t for t in taxa_all if t in dates_map and not np.isnan(dates_map[t]) and (root_taxon is None or t != root_taxon)]
    n_dated = len(dated_taxa)
    print(f"[*] Timestamps mapped: {n_dated}/{len(taxa_all)} taxa successfully dated.")
    if missing and len(missing) <= 10:
        print(f"    Notice: {len(missing)} taxa omitted due to missing timestamps: {missing}")
    elif missing:
        print(f"    Notice: {len(missing)} taxa omitted due to missing timestamps.")

    if n_dated < 3:
        raise ValueError(
            f"Fewer than 3 taxa could be mapped to valid timestamps ({n_dated} dated). "
            f"Please check your metadata file (--dates) or FASTA header format."
        )

    # 3. Compute Patristic or Tree-Free Divergences
    has_tree = (tree_path is not None and Path(tree_path).exists() and not use_tn93)
    if has_tree:
        print(f"[*] Computing patristic tree distances from: {tree_path}...")
        tree_dists, root_desc = extract_tree_root_to_tip(
            str(tree_path), dated_taxa, dates_map,
            root_taxon=root_taxon, optimize_root=optimize_root
        )
        taxa = [t for t in dated_taxa if t in tree_dists]
        times = np.array([dates_map[t] for t in taxa], dtype=np.float64)
        dists = np.array([tree_dists[t] for t in taxa], dtype=np.float64)
    else:
        print(f"[*] Tree skipped: Estimating tree-free pairwise distances via TN93...")
        dists, root_desc = compute_tree_free_divergences(
            seq_dict, dated_taxa, dates_map,
            root_taxon=root_taxon, decay_gamma=decay_gamma, decay_half_life=decay_half_life
        )
        if root_taxon and root_taxon in seq_dict:
            taxa = [t for t in dated_taxa if t != root_taxon]
        else:
            taxa = list(dated_taxa)
        times = np.array([dates_map[t] for t in taxa], dtype=np.float64)

    print(f"[*] Root configuration: {root_desc} (Timespan: {np.min(times):.1f} - {np.max(times):.1f})")

    # 4. Fit Standard OLS
    ols_res = run_ols_dating(times, dists, n_boot=n_bootstrap)
    print(f"[✓] OLS Molecular Clock: t_MRCA = {ols_res['t_mrca']:.2f} [{ols_res['ci_mrca'][0]:.1f}, {ols_res['ci_mrca'][1]:.1f}], μ = {ols_res['mu']:.6f} subs/site/yr (R^2 = {ols_res['r2']:.3f})")

    # 5. HyphAeon Neural Attention PGLS (if requested)
    pgls_res = None
    cov_matrix = None
    opt_lambda = 0.95

    run_neural = method in ["all", "pgls"]
    if run_neural:
        if device is None:
            device = get_device()
        if model is None:
            print(f"[*] Loading HyphAeon transformer backbone on {device}...")
            model = load_model(weights=weights, variant=variant, device=device)

        # Load alignment into model tensors
        c, a, d, z, inv, aln_taxa, L, tree_cache = prepare_alignment(
            str(align_p),
            str(tree_path) if has_tree else None,
            model=model,
            device=device,
            max_species=max_species,
            prune_duplicates=False,
            use_tn93=(not has_tree)
        )

        t_fwd = time.time()
        print(f"[*] Extracting cross-taxa attention and 128D continuous representations...")
        msa_codons = c.to(device)
        msa_aas = a.to(device)

        cross_attn, taxa_repr = extract_cross_taxa_attentions_and_embeddings(
            model, msa_codons, msa_aas, tree_cache, device=device
        )
        K_neural = compute_neural_covariance_kernel(cross_attn, taxa_repr)
        print(f"[✓] Forward pass complete in {time.time() - t_fwd:.2f}s! Extracted {taxa_repr.shape[0]} taxa neural representations.")

        # Align taxa order with dated taxa
        aln_taxa_map = {t: i for i, t in enumerate(aln_taxa)}
        sub_indices = [aln_taxa_map[t] for t in taxa if t in aln_taxa_map]
        sub_taxa = [taxa[i] for i, t in enumerate(taxa) if t in aln_taxa_map]
        sub_times = times[[i for i, t in enumerate(taxa) if t in aln_taxa_map]]
        sub_dists = dists[[i for i, t in enumerate(taxa) if t in aln_taxa_map]]

        cov_matrix = K_neural[sub_indices, :][:, sub_indices]
        z_sub = taxa_repr[sub_indices]

        effective_ridge = 0.05
        if isinstance(ridge, (int, float)):
            effective_ridge = float(ridge)
        elif str(ridge).lower() in ["auto", "reml"]:
            opt_reml = estimate_reml_pagel_lambda(sub_times, sub_dists, cov_matrix)
            # Map Pagel's lambda to complementary nugget ridge: ridge = (1 - lambda)
            effective_ridge = float(np.clip(1.0 - opt_reml['best_lambda'], 0.01, 0.20))
            print(f"[*] REML Estimated Phylogenetic Signal: Pagel's λ* = {opt_reml['best_lambda']:.4f} (nugget ridge = {effective_ridge:.4f})")

        # Fit PGLS
        if method in ["all", "pgls"]:
            pgls_res = run_pgls_dating(sub_times, sub_dists, cov_matrix, ridge=effective_ridge, n_boot=n_bootstrap)
            print(f"[✓] HyphAeon PGLS Clock: t_MRCA = {pgls_res['t_mrca']:.2f} [{pgls_res['ci_mrca'][0]:.1f}, {pgls_res['ci_mrca'][1]:.1f}], μ = {pgls_res['mu']:.6f} subs/site/yr (R^2_gls = {pgls_res['r2']:.3f})")

    # 5b. Non-Linear Clock Models: Restricted Natural Spline & Power-Law
    spline_res = None
    power_res = None

    if clock_model in ["auto", "spline"]:
        fit_times = sub_times if (run_neural and cov_matrix is not None) else times
        fit_dists = sub_dists if (run_neural and cov_matrix is not None) else dists
        spline_cov = cov_matrix if (run_neural and cov_matrix is not None) else None
        try:
            spline_res = run_restricted_spline_clock_dating(
                fit_times, fit_dists, cov_matrix=spline_cov, ridge=effective_ridge, n_boot=min(500, n_bootstrap)
            )
            ratio_sym = "acceleration" if spline_res['rate_ratio'] > 1.0 else "deceleration"
            if np.isnan(spline_res['t_mrca']) or spline_res['rate_ancestral'] <= 0:
                t0_str = "n/a (ancestral rate <= 0)"
            else:
                t0_str = f"{spline_res['t_mrca']:.2f} [{spline_res['ci_mrca'][0]:.1f}, {spline_res['ci_mrca'][1]:.1f}]"
            print(f"[✓] Restricted Spline Clock: t_MRCA = {t0_str}, μ_anc = {spline_res['rate_ancestral']:.6f}, μ_rec = {spline_res['rate_recent']:.6f} ({ratio_sym} {spline_res['rate_ratio']:.2f}x) (R^2 = {spline_res['r2']:.3f}, ΔAIC = {spline_res['delta_aic']:+.1f}, p = {spline_res['p_f_test']:.4f})")
        except Exception as e:
            print(f"[!] Notice: Restricted spline fitting fell back to linear ({e})")

    elif clock_model == "power":
        fit_times = sub_times if (run_neural and cov_matrix is not None) else times
        fit_dists = sub_dists if (run_neural and cov_matrix is not None) else dists
        power_cov = opt_lambda * cov_matrix + (1.0 - opt_lambda) * np.eye(len(fit_times)) if (run_neural and cov_matrix is not None) else None
        try:
            power_res = run_powerlaw_clock_dating(
                fit_times, fit_dists, cov_matrix=power_cov, ridge=0.01, n_boot=min(500, n_bootstrap)
            )
            print(f"[✓] Power-Law Clock: t_MRCA = {power_res['t_mrca']:.2f} [{power_res['ci_mrca'][0]:.1f}, {power_res['ci_mrca'][1]:.1f}], θ = {power_res['theta']:.3f} [{power_res['ci_theta'][0]:.3f}, {power_res['ci_theta'][1]:.3f}], μ_mean = {power_res['rate_mean']:.6f} subs/site/yr (R^2 = {power_res['r2']:.3f}, ΔAIC = {power_res['delta_aic']:+.1f}, p = {power_res['p_f_test']:.4f})")
        except Exception as e:
            print(f"[!] Notice: Non-linear power-law fitting fell back to linear ({e})")

    # Model Selection
    selected_clock = "Linear"
    if clock_model == "spline" and spline_res is not None:
        active_model = spline_res
        selected_clock = "Restricted Spline (forced)"
    elif clock_model == "power" and power_res is not None:
        active_model = power_res
        selected_clock = "Power-Law (forced)"
    elif clock_model == "linear":
        active_model = pgls_res if pgls_res is not None else ols_res
        selected_clock = "Linear (forced)"
    else:  # auto
        if spline_res is not None and spline_res['is_nonlinear_preferred']:
            active_model = spline_res
            ratio_str = f"acceleration ({spline_res['rate_ratio']:.2f}x)" if spline_res['rate_ratio'] > 1.0 else f"deceleration ({spline_res['rate_ratio']:.2f}x)"
            selected_clock = f"Restricted Spline (rate {ratio_str} detected: F={spline_res['f_stat']:.2f}, p={spline_res['p_f_test']:.4f}, ΔAIC={spline_res['delta_aic']:+.1f})"
        else:
            active_model = pgls_res if pgls_res is not None else ols_res
            sp_p = f"p={spline_res['p_f_test']:.4f}" if spline_res else "p=n/a"
            selected_clock = f"Linear (parsimonious linear clock preferred; {sp_p})"

    print(f"[✓] Clock Model Selection: {selected_clock}")

    # 6. Per-Taxon Residuals and Predictions
    if active_model.get('method') == 'RESTRICTED_SPLINE':
        b0 = active_model['beta_0']
        b1 = active_model['beta_1']
        b2 = active_model['beta_2']
        knots_arr = np.array(active_model['knots'])
        t_kn0 = knots_arr[0]
        d_kn0 = b0 + b1 * t_kn0

        B_all, _ = compute_rcs_basis(times, knots_arr)
        fitted_all = b0 + b1 * times + (b2 * B_all[:, 0] if B_all.shape[1] > 0 else 0.0)

        pred_dates = np.zeros(len(dists))
        for idx_d, d_val in enumerate(dists):
            if d_val <= d_kn0 or abs(b2) < 1e-12:
                pred_dates[idx_d] = (d_val - b0) / max(1e-15, b1)
            else:
                def f_diff(t_cand):
                    B_c, _ = compute_rcs_basis(np.array([t_cand]), knots_arr)
                    return float(b0 + b1 * t_cand + b2 * B_c[0, 0] - d_val)
                try:
                    t_root_sol = optimize.brentq(f_diff, t_kn0, max(times) + 100.0)
                    pred_dates[idx_d] = t_root_sol
                except Exception:
                    pred_dates[idx_d] = (d_val - b0) / max(1e-15, b1)
    elif active_model.get('method') == 'POWER_LAW':
        k_val = max(1e-12, active_model['k'])
        th_val = max(1e-4, active_model['theta'])
        t0_val = active_model['t_mrca']
        fitted_all = k_val * (np.maximum(1e-6, times - t0_val) ** th_val)
        pred_dates = t0_val + (np.maximum(0.0, dists) / k_val) ** (1.0 / th_val)
    else:
        fitted_all = active_model['d0'] + active_model['mu'] * (times - active_model['t_ref'])
        pred_dates = active_model['t_ref'] + (dists - active_model['d0']) / max(1e-15, active_model['mu'])

    residuals = dists - fitted_all
    std_res = np.std(residuals) if np.std(residuals) > 1e-12 else 1.0

    taxon_records = []
    for i, t in enumerate(taxa):
        z_score = float(residuals[i] / std_res)
        is_outlier = bool(abs(z_score) >= 2.5)
        taxon_records.append({
            'taxon': t,
            'sampling_date': float(times[i]),
            'root_divergence': float(dists[i]),
            'fitted_divergence': float(fitted_all[i]),
            'predicted_date': float(pred_dates[i]),
            'divergence_residual': float(residuals[i]),
            'temporal_residual': float(pred_dates[i] - times[i]),
            'z_score': z_score,
            'is_outlier': is_outlier
        })

    elapsed_time = time.time() - t0

    results = {
        'alignment': str(align_p),
        'tree': str(tree_path) if has_tree else None,
        'root_description': root_desc,
        'taxa_count': len(taxa),
        'timespan': [float(np.min(times)), float(np.max(times))],
        'elapsed_seconds': elapsed_time,
        'times': times,
        'dists': dists,
        'taxa': taxa,
        'ols': ols_res,
        'pgls': pgls_res,
        'spline': spline_res,
        'power': power_res,
        'clock_model': clock_model,
        'selected_clock': selected_clock,
        'taxa_records': taxon_records
    }

    # 7. Export Results
    if output_prefix:
        out_p = Path(output_prefix)
        ensure_parent_directory(out_p)

        json_data = {
            'alignment': results['alignment'],
            'tree': results['tree'],
            'root_description': results['root_description'],
            'taxa_count': results['taxa_count'],
            'timespan': results['timespan'],
            'elapsed_seconds': results['elapsed_seconds'],
            'ols': {k: v for k, v in ols_res.items() if k not in ['residuals', 'fitted']},
            'pgls': {k: v for k, v in pgls_res.items() if k not in ['residuals', 'fitted', 'times']} if pgls_res else None,
            'spline': {k: v for k, v in spline_res.items() if k not in ['residuals', 'fitted']} if spline_res else None,
            'power': {k: v for k, v in power_res.items() if k not in ['residuals', 'fitted']} if power_res else None,
            'clock_model': clock_model,
            'selected_clock': selected_clock,
            'taxa_summary': taxon_records
        }
        json_file = out_p.with_suffix('.json') if not str(out_p).endswith('.json') else out_p
        write_json(json_file, json_data)
        print(f"[✓] Saved JSON summary: {json_file}")

        csv_file = out_p.with_suffix('.csv') if not str(out_p).endswith('.csv') else out_p
        pd.DataFrame(taxon_records).to_csv(csv_file, index=False)
        print(f"[✓] Saved per-taxon CSV: {csv_file}")

    # 8. Diagnostic Plotting
    if plot:
        fig_path = f"{output_prefix}_diagnostic.pdf" if output_prefix else "mrca_dating_diagnostic.pdf"
        plot_mrca_dating(results, fig_path)

    return results
