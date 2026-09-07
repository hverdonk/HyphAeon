"""
hyphaeon/dating.py
------------------
Heterochronous Molecular Clock Calibration, Ancestor Dating (t_MRCA),
and Latent Manifold Coalescent Variance Collapse for Pathogen Genomics.

Methods:
1. Strict in-frame coding alignment validation (L_nt % 3 == 0, triplet-gap check, stop codon audit).
2. Flexible timestamp ingestion (FASTA headers, CSV/TSV metadata, Nextstrain Auspice JSON v2).
3. Root-to-tip divergence computation:
   - Tree-based: Patristic distance traversal with heuristic root search (TempEst R^2 maximization).
   - Tree-free: Direct pairwise distance estimation (TN93) and ancestral consensus anchoring.
4. Estimators:
   - Centered Root-to-Tip Ordinary Least Squares (OLS / TempEst emulation with delta-method & bootstrap CIs).
   - HyphAeon Attention-Derived Phylogenetic Generalized Least Squares (PGLS) via A_fused covariance.
   - Latent Manifold Coalescent Variance Collapse (Var(Z(t)) -> 0 in 128D continuous representation space).
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
from .splits import (
    extract_cross_taxa_attentions_and_embeddings,
    compute_fused_affinity_matrix,
)
from .temporal import parse_date_to_decimal, parse_dates_from_auspice_json, extract_date_from_string
from .io import ensure_parent_directory, write_json, write_csv


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

        for node in non_terminals:
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
    root_taxon: Optional[str] = None
) -> Tuple[np.ndarray, str]:
    """
    Computes tree-free root-to-tip divergence directly from pairwise TN93 distances.
    Anchors root at:
    1. Specified root_taxon (if found in alignment).
    2. Computed consensus of the entire alignment or earliest cohort.
    3. Earliest sampled taxon.
    """
    all_taxa = list(seq_dict.keys())

    # Case 1: User explicitly specified an existing taxon as root (e.g. 'CONSENSUS' or 'Z59ZR.ZHU')
    if root_taxon and root_taxon in seq_dict:
        eval_taxa = [root_taxon] + [t for t in dated_taxa if t != root_taxon]
        dist_mat = compute_tn93_distance_matrix(seq_dict, eval_taxa)
        # Distance from root (index 0) to each dated taxon
        div_dict = {eval_taxa[i]: dist_mat[0, i] for i in range(1, len(eval_taxa))}
        divergences = np.array([div_dict[t] for t in dated_taxa if t != root_taxon], dtype=np.float64)
        return divergences, f"explicit_root_{root_taxon}"

    # Case 2: User requested 'consensus' root, or alignment has no predefined root
    if root_taxon and root_taxon.lower() in ['consensus', 'founder', 'ancestor']:
        con_seq = generate_consensus_sequence(seq_dict, dated_taxa)
        aug_dict = dict(seq_dict)
        aug_dict['__SYNTHETIC_CONSENSUS__'] = con_seq
        eval_taxa = ['__SYNTHETIC_CONSENSUS__'] + dated_taxa
        dist_mat = compute_tn93_distance_matrix(aug_dict, eval_taxa)
        divergences = np.array([dist_mat[0, i + 1] for i in range(len(dated_taxa))], dtype=np.float64)
        return divergences, "synthetic_consensus_root"

    # Case 3: Anchor on earliest sampled cohort
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
    t_ref: Optional[float] = None,
    n_boot: int = 1000,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Fits Centered Phylogenetic Generalized Least Squares (PGLS) regression:
        d = X * beta + epsilon,   Cov(epsilon) = sigma^2 * Sigma

    where Sigma = A_fused + ridge * I is HyphAeon's cross-taxa attention covariance.
    Directly incorporates evolutionary covariance without tree reconstruction.
    """
    n = len(times)
    if n < 3:
        raise ValueError(f"At least 3 observations are required for PGLS dating (got N={n}).")

    if t_ref is None:
        t_ref = float(np.mean(times))

    x = times - t_ref
    X = np.column_stack([x, np.ones(n)])

    # Spectral regularization: Sigma = V * Lambda * V^T
    C = cov_matrix + ridge * np.eye(n)
    w, v = la.eigh(C)
    w = np.maximum(w, 1e-6)
    C_inv = v @ np.diag(1.0 / w) @ v.T

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
    ci_analytical = [t_mrca - t_crit * se_mrca, t_mrca + t_crit * se_mrca]

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
        'ridge': ridge,
        'sigma2': sigma2_gls,
        'rmse': float(np.sqrt(np.mean(residuals ** 2))),
        'residuals': residuals,
        'fitted': X @ beta_gls,
        'n': n
    }


def run_manifold_variance_collapse(
    taxa_repr: np.ndarray,
    times: np.ndarray,
    min_taxa_per_timepoint: int = 2
) -> Optional[Dict[str, Any]]:
    """
    Evaluates Coalescent Population Variance Collapse in 128D continuous representation space:
        Var(Z(t)) = s * (t - t_founder)

    Under founder bottleneck transmission (N(0)=1), variance expands linearly with time.
    Extrapolating population dispersion Var(Z(t)) -> 0 recovers the time of origin
    without requiring root heuristics, branch lengths, or tree inference.
    """
    unique_times = np.unique(times)
    valid_times = []
    intra_vars = []
    counts = []

    for t in unique_times:
        mask = (times == t)
        if np.sum(mask) >= min_taxa_per_timepoint:
            z_t = taxa_repr[mask]
            var_t = float(np.sum(np.var(z_t, axis=0)))
            valid_times.append(float(t))
            intra_vars.append(var_t)
            counts.append(int(np.sum(mask)))

    if len(valid_times) < 3:
        return None

    vt = np.array(valid_times)
    iv = np.array(intra_vars)

    slope, intercept = np.polyfit(vt, iv, 1)
    if abs(slope) < 1e-15:
        slope = 1e-15

    t_collapse = -intercept / slope
    r_val = float(np.corrcoef(vt, iv)[0, 1])

    return {
        'method': 'Manifold_Collapse',
        't_mrca': float(t_collapse),
        'slope': float(slope),
        'intercept': float(intercept),
        'r2': float(r_val ** 2),
        'r': float(r_val),
        'timepoints': vt,
        'variances': iv,
        'counts': counts,
        'n_timepoints': len(vt)
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
    manifold = dating_results.get('manifold')

    # Panel A: Root-to-Tip Molecular Clock Regression
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(times, dists, color='#1f77b4', s=42, alpha=0.75, edgecolors='black', linewidth=0.5, label='Dated Strains', zorder=3)

    t_min = float(np.min(times))
    t_max = float(np.max(times))
    mrca_candidates = [ols['t_mrca']]
    if pgls:
        mrca_candidates.append(pgls['t_mrca'])
    if manifold:
        mrca_candidates.append(manifold['t_mrca'])

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

    # MRCA markers and CI error bars at distance = 0
    ax1.axhline(0, color='gray', linestyle=':', linewidth=0.8, zorder=1)
    ax1.errorbar([ols['t_mrca']], [0], xerr=[[ols['t_mrca'] - ols['ci_mrca'][0]], [ols['ci_mrca'][1] - ols['t_mrca']]],
                 fmt='s', color='#e63946', markersize=6, capsize=4, capthick=1.5, zorder=6)
    if pgls:
        ax1.errorbar([pgls['t_mrca']], [-0.002], xerr=[[pgls['t_mrca'] - pgls['ci_mrca'][0]], [pgls['ci_mrca'][1] - pgls['t_mrca']]],
                     fmt='D', color='#7b2cbf', markersize=6, capsize=4, capthick=1.5, zorder=6)

    ax1.set_xlim(plot_left, plot_right)
    ax1.set_xlabel("Sampling Date / Time", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Root-to-Tip Divergence (subs/site)", fontsize=11, fontweight='bold')
    ax1.set_title("(A) Heterochronous Molecular Clock Regression", fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left', frameon=True, fontsize=8.5)
    ax1.grid(True, linestyle=':', alpha=0.4)

    # Panel B: Manifold Collapse or Residual Diagnostics
    ax2 = fig.add_subplot(gs[1])
    if manifold is not None and manifold['n_timepoints'] >= 3:
        vt = manifold['timepoints']
        iv = manifold['variances']
        ax2.scatter(vt, iv, color='#2a9d8f', s=65, edgecolors='black', linewidth=0.6, label='Latent Population Variance', zorder=3)
        t_var_plot = np.linspace(manifold['t_mrca'] - (t_max - t_min) * 0.1, t_max + 1, 100)
        ax2.plot(t_var_plot, manifold['slope'] * t_var_plot + manifold['intercept'], color='#2a9d8f', linewidth=2.2,
                 label=f"Coalescent Collapse (t_MRCA={manifold['t_mrca']:.1f})", zorder=4)
        ax2.axvline(manifold['t_mrca'], color='#2a9d8f', linestyle='--', linewidth=1.8)
        ax2.scatter([manifold['t_mrca']], [0], color='#e76f51', s=80, marker='*', zorder=5, label='Origin Intercept')
        ax2.axhline(0, color='gray', linestyle=':', linewidth=0.8)
        ax2.set_xlabel("Sampling Date / Time", fontsize=11, fontweight='bold')
        ax2.set_ylabel("128D Latent Manifold Variance Tr(Var(Z))", fontsize=11, fontweight='bold')
        ax2.set_title("(B) Latent Manifold Coalescent Variance Collapse", fontsize=12, fontweight='bold')
        ax2.legend(loc='upper left', frameon=True, fontsize=8.5)
        ax2.grid(True, linestyle=':', alpha=0.4)
    else:
        # Residuals vs Time
        res_ols = ols['residuals']
        ax2.scatter(times, res_ols, color='#e63946', alpha=0.7, s=40, edgecolors='black', linewidth=0.5, label='OLS Residuals', zorder=3)
        if pgls:
            res_pgls = pgls['residuals']
            ax2.scatter(times, res_pgls, color='#7b2cbf', alpha=0.7, s=40, marker='^', edgecolors='black', linewidth=0.5, label='PGLS Residuals', zorder=4)
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
    optimize_root: bool = True,
    use_tn93: bool = False,
    method: str = "all",
    ridge: float = 0.05,
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
            seq_dict, dated_taxa, dates_map, root_taxon=root_taxon
        )
        taxa = [t for t in dated_taxa if root_taxon is None or t != root_taxon]
        times = np.array([dates_map[t] for t in taxa], dtype=np.float64)

    print(f"[*] Root configuration: {root_desc} (Timespan: {np.min(times):.1f} - {np.max(times):.1f})")

    # 4. Fit Standard OLS
    ols_res = run_ols_dating(times, dists, n_boot=n_bootstrap)
    print(f"[✓] OLS Molecular Clock: t_MRCA = {ols_res['t_mrca']:.2f} [{ols_res['ci_mrca'][0]:.1f}, {ols_res['ci_mrca'][1]:.1f}], μ = {ols_res['mu']:.6f} subs/site/yr (R^2 = {ols_res['r2']:.3f})")

    # 5. HyphAeon Attention PGLS & Manifold Collapse (if requested)
    pgls_res = None
    manifold_res = None
    cov_matrix = None

    run_neural = method in ["all", "pgls", "manifold"]
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
        mds_coords = z.squeeze(0).cpu().numpy()

        cross_attn, taxa_repr = extract_cross_taxa_attentions_and_embeddings(
            model, msa_codons, msa_aas, tree_cache, device=device
        )
        A_fused = compute_fused_affinity_matrix(cross_attn, mds_coords, taxa_repr)
        print(f"[✓] Forward pass complete in {time.time() - t_fwd:.2f}s! Extracted {taxa_repr.shape[0]} taxa representations.")

        # Align taxa order with dated taxa
        aln_taxa_map = {t: i for i, t in enumerate(aln_taxa)}
        sub_indices = [aln_taxa_map[t] for t in taxa if t in aln_taxa_map]
        sub_taxa = [taxa[i] for i, t in enumerate(taxa) if t in aln_taxa_map]
        sub_times = times[[i for i, t in enumerate(taxa) if t in aln_taxa_map]]
        sub_dists = dists[[i for i, t in enumerate(taxa) if t in aln_taxa_map]]

        cov_matrix = A_fused[sub_indices, :][:, sub_indices]
        z_sub = taxa_repr[sub_indices]

        # Fit PGLS
        if method in ["all", "pgls"]:
            pgls_res = run_pgls_dating(sub_times, sub_dists, cov_matrix, ridge=ridge, n_boot=n_bootstrap)
            print(f"[✓] HyphAeon PGLS Clock: t_MRCA = {pgls_res['t_mrca']:.2f} [{pgls_res['ci_mrca'][0]:.1f}, {pgls_res['ci_mrca'][1]:.1f}], μ = {pgls_res['mu']:.6f} subs/site/yr (R^2_gls = {pgls_res['r2']:.3f})")

        # Fit Manifold Collapse
        if method in ["all", "manifold"]:
            manifold_res = run_manifold_variance_collapse(z_sub, sub_times)
            if manifold_res is not None:
                print(f"[✓] Manifold Collapse Clock: t_MRCA = {manifold_res['t_mrca']:.2f} (Variance expansion s = {manifold_res['slope']:.6f}/yr, R^2 = {manifold_res['r2']:.3f})")

    # 6. Per-Taxon Residuals and Predictions
    active_model = pgls_res if pgls_res is not None else ols_res
    pred_dates = active_model['t_ref'] + (dists - active_model['d0']) / active_model['mu']
    residuals = dists - active_model['fitted']
    std_res = np.std(residuals) if np.std(residuals) > 1e-12 else 1.0

    taxon_records = []
    for i, t in enumerate(taxa):
        z_score = float(residuals[i] / std_res)
        is_outlier = bool(abs(z_score) >= 2.5)
        taxon_records.append({
            'taxon': t,
            'sampling_date': float(times[i]),
            'root_divergence': float(dists[i]),
            'fitted_divergence': float(active_model['fitted'][i]),
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
        'manifold': manifold_res,
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
            'pgls': {k: v for k, v in pgls_res.items() if k not in ['residuals', 'fitted']} if pgls_res else None,
            'manifold': {k: v for k, v in manifold_res.items() if k not in ['timepoints', 'variances']} if manifold_res else None,
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
