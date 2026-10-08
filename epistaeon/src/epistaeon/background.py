"""
epistaeon/background.py
-----------------------
Latent embeddings for hypothetical genotypes.

A genotype is scored by writing its residues into one row of the alignment --
the focal taxon -- and reading that row's embedding out of the model. The tree
is left untouched, so the genotype inherits the focal taxon's phylogenetic
position, exactly as hyphaeon's in-silico DMS does.

The expensive part is avoided by exploiting a property of the published
checkpoint: it has no cross-site mixing (window_size = 1 and no column layers),
so a taxon's embedding is the mean over sites of per-site hidden states, and
each per-site state depends only on its own alignment column. Therefore

    z(S) = (1/L) * ( sum_s h_s(anc) - sum_k h_k(anc) + sum_k h_k(state_k) )

and every genotype's embedding follows from 1 + 2K forward passes rather than
one per genotype. `verify_additivity` must pass before that shortcut is used;
if the checkpoint ever gains cross-site layers it will fail, and the caller
should fall back to `embed_genotype_direct`.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

from aeon_core.dataset import (
    AA_MAP,
    CODON_TO_AA,
    GENETIC_CODE,
    load_alignment_and_tree,
)
from aeon_core.inference import get_device, load_model
from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings

from .lattice import GAP, Lattice

# One canonical codon per amino acid, matching hyphaeon's DMS convention.
CANONICAL_AA_TO_CODON = {
    "A": "GCC", "C": "TGC", "D": "GAC", "E": "GAG", "F": "TTC",
    "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAG", "L": "CTG",
    "M": "ATG", "N": "AAC", "P": "CCC", "Q": "CAG", "R": "CGC",
    "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
}
GAP_CODON_TOKEN = 64
GAP_AA_TOKEN = 20


@dataclass
class Context:
    """A loaded alignment, tree cache and model, plus the focal row index."""

    model: torch.nn.Module
    c_tensor: torch.Tensor
    a_tensor: torch.Tensor
    tree_cache: dict
    taxa: List[str]
    n_codons: int
    device: torch.device
    focal_index: int
    focal_name: str
    variant: Optional[str] = None

    @property
    def n_taxa(self) -> int:
        return len(self.taxa)


def load_context(
    alignment: str,
    tree: Optional[str],
    focal_taxon: Optional[str] = None,
    weights: Optional[str] = None,
    variant: Optional[str] = None,
    cpu: bool = False,
    max_species: Optional[int] = None,
    use_tn93: bool = False,
) -> Context:
    """Load alignment + tree + model and pick the row genotypes are written to."""
    c, a, d, z, _inv, taxa, L = load_alignment_and_tree(
        alignment, tree, max_species=max_species, prune_duplicates=True, use_tn93=use_tn93
    )
    device = get_device(cpu=cpu)
    model = load_model(weights=weights, variant=variant, device=device)
    tree_cache = model.precompute_tree_cache(d.to(device), z.to(device))

    idx, name = 0, taxa[0]
    if focal_taxon:
        matches = [i for i, t in enumerate(taxa) if focal_taxon.lower() in t.lower()]
        if not matches:
            raise ValueError(
                f"focal taxon {focal_taxon!r} is not in the alignment after taxon "
                f"matching and duplicate pruning; first few taxa are {taxa[:5]}"
            )
        idx, name = matches[0], taxa[matches[0]]
    return Context(
        model=model, c_tensor=c, a_tensor=a, tree_cache=tree_cache, taxa=taxa,
        n_codons=L, device=device, focal_index=idx, focal_name=name, variant=variant,
    )


def _write_residue(c: torch.Tensor, a: torch.Tensor, site: int, row: int, residue: str) -> None:
    """Write one residue (or a gap) into the focal row at one site, in place."""
    if residue == GAP:
        c[site, row, 0] = GAP_CODON_TOKEN
        a[site, row, 0] = GAP_AA_TOKEN
        return
    codon = CANONICAL_AA_TO_CODON.get(residue.upper())
    if codon is None:
        raise ValueError(f"no canonical codon for residue {residue!r}")
    c[site, row, 0] = GENETIC_CODE.get(codon, GAP_CODON_TOKEN)
    a[site, row, 0] = AA_MAP.get(CODON_TO_AA.get(codon, GAP), GAP_AA_TOKEN)


def _tensors_for(ctx: Context, lattice: Lattice, mask: int, sites: Optional[Sequence[int]] = None):
    """Copies of the alignment tensors carrying this genotype in the focal row."""
    rows = list(range(ctx.n_codons)) if sites is None else list(sites)
    c = ctx.c_tensor[rows].clone()
    a = ctx.a_tensor[rows].clone()
    where = {s: i for i, s in enumerate(rows)}
    seq = lattice.sequence(mask)
    for site, i in where.items():
        if site < len(seq):
            _write_residue(c, a, i, ctx.focal_index, seq[site])
    return c, a


def _focal_embedding(ctx: Context, c: torch.Tensor, a: torch.Tensor) -> np.ndarray:
    """Mean-over-sites embedding of the focal taxon for these tensors."""
    _attn, reprs = extract_cross_taxa_attentions_and_embeddings(
        ctx.model, c, a, ctx.tree_cache, device=ctx.device
    )
    return np.asarray(reprs[ctx.focal_index], dtype=np.float64)


def embed_genotype_direct(ctx: Context, lattice: Lattice, mask: int) -> np.ndarray:
    """One full forward pass for this genotype. The slow, always-correct path."""
    c, a = _tensors_for(ctx, lattice, mask)
    return _focal_embedding(ctx, c, a)


@dataclass
class AdditiveCache:
    """Per-site deltas that let any genotype's embedding be computed by hand."""

    ancestral_sum: np.ndarray                    # sum over all sites, ancestral background
    per_unit: Dict[int, Tuple[np.ndarray, np.ndarray]]   # unit -> (ancestral, derived) site vectors
    n_codons: int
    passes_used: int = 0
    additive: bool = True
    additivity_error: float = float("nan")

    def embedding(self, lattice: Lattice, mask: int) -> np.ndarray:
        total = self.ancestral_sum.copy()
        for u in range(lattice.n_units):
            anc_vec, der_vec = self.per_unit[u]
            if mask >> u & 1:
                total += der_vec - anc_vec
        return total / float(self.n_codons)


def build_additive_cache(ctx: Context, lattice: Lattice) -> AdditiveCache:
    """1 + 2K forward passes: the ancestral background, then each unit's states."""
    anc_mask = 0
    c, a = _tensors_for(ctx, lattice, anc_mask)
    anc_mean = _focal_embedding(ctx, c, a)
    ancestral_sum = anc_mean * float(ctx.n_codons)
    passes = 1

    per_unit: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for u in range(lattice.n_units):
        sites = sorted({sub.index for sub in lattice.unit_members(u)})
        c_a, a_a = _tensors_for(ctx, lattice, anc_mask, sites=sites)
        vec_anc = _focal_embedding(ctx, c_a, a_a) * float(len(sites))
        c_d, a_d = _tensors_for(ctx, lattice, 1 << u, sites=sites)
        vec_der = _focal_embedding(ctx, c_d, a_d) * float(len(sites))
        per_unit[u] = (vec_anc, vec_der)
        passes += 2

    return AdditiveCache(
        ancestral_sum=ancestral_sum, per_unit=per_unit,
        n_codons=ctx.n_codons, passes_used=passes,
    )


def verify_additivity(
    ctx: Context,
    lattice: Lattice,
    cache: AdditiveCache,
    mask: Optional[int] = None,
    tolerance: float = 1e-4,
) -> Tuple[bool, float]:
    """Compare a cached genotype embedding against a direct forward pass.

    This gates the fast path. It fails if the checkpoint ever mixes information
    across sites, in which case embeddings are not additive and the cache is
    invalid.
    """
    if mask is None:
        mask = 0b11 if lattice.n_units >= 2 else 1
    predicted = cache.embedding(lattice, mask)
    actual = embed_genotype_direct(ctx, lattice, mask)
    scale = max(1e-12, float(np.abs(actual).max()))
    err = float(np.abs(predicted - actual).max() / scale)
    ok = err <= tolerance
    cache.additive = ok
    cache.additivity_error = err
    return ok, err


def extant_embeddings(ctx: Context) -> np.ndarray:
    """Embeddings of every taxon under the unmodified alignment.

    Used as the reference cloud for typicality scoring, and for the isometric
    calibration below.
    """
    _attn, reprs = extract_cross_taxa_attentions_and_embeddings(
        ctx.model, ctx.c_tensor, ctx.a_tensor, ctx.tree_cache, device=ctx.device
    )
    return np.asarray(reprs, dtype=np.float64)


def isometric_scale(embeddings: np.ndarray, physical: np.ndarray) -> float:
    """Least-squares alpha mapping latent distance onto substitutions/site.

    The same calibration chronaeon performs inside its hull optimiser, but
    computed directly so it needs no tip dates.
    """
    iu = np.triu_indices(embeddings.shape[0], k=1)
    lat = np.linalg.norm(embeddings[iu[0]] - embeddings[iu[1]], axis=1)
    phys = np.asarray(physical)[iu]
    denom = float(np.sum(lat ** 2))
    if denom <= 1e-12:
        return 1.0
    return float(np.sum(phys * lat) / denom)


def root_token_embedding(ctx: Context) -> np.ndarray:
    """The model's learned [ROOT] token representation, averaged over sites.

    The fallback phi anchor when no ancestral reconstruction is available. It
    lives in the same latent space as the taxon embeddings and needs no dates,
    unlike chronaeon's convex-hull root.
    """
    ctx.model.eval()
    acc = None
    n = 0
    step = max(1, min(64, ctx.n_codons))
    with torch.no_grad():
        for start in range(0, ctx.n_codons, step):
            stop = min(ctx.n_codons, start + step)
            c = ctx.c_tensor[start:stop].to(ctx.device)
            a = ctx.a_tensor[start:stop].to(ctx.device)
            out = ctx.model.forward_cached(c, a, ctx.tree_cache, return_hidden=True)
            root_repr = out[2].detach().cpu().numpy()
            acc = root_repr.sum(axis=0) if acc is None else acc + root_repr.sum(axis=0)
            n += root_repr.shape[0]
    return np.asarray(acc, dtype=np.float64) / max(1, n)


def focal_protein(ctx: Context) -> Tuple[str, List[int]]:
    """The focal row's protein sequence, with the codon column of each residue.

    Gap columns are skipped, so the returned string is ungapped and `columns[i]`
    is the alignment column holding residue `i`.
    """
    rev = {v: k for k, v in AA_MAP.items()}
    a = ctx.a_tensor.squeeze(-1).numpy()
    letters: List[str] = []
    columns: List[int] = []
    for site in range(ctx.n_codons):
        ch = rev.get(int(a[site, ctx.focal_index]), GAP)
        if ch != GAP:
            letters.append(ch)
            columns.append(site)
    return "".join(letters), columns


def map_to_alignment(ctx: Context, query: str) -> Dict[int, int]:
    """Map positions in `query` onto alignment columns, via the focal row.

    The ancestral and derived sequences are usually a domain (the receptor LBD
    is 250 residues) while the alignment is full length, so substitution
    indices have to be translated before residues can be written into the
    alignment. Returns {query_index -> codon column} for the positions that
    align; callers should treat missing keys as unmappable rather than guess.
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    target, columns = focal_protein(ctx)
    ungapped = [(i, ch) for i, ch in enumerate(query) if ch != GAP]
    probe = "".join(ch for _i, ch in ungapped)

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"
    best = aligner.align(probe, target)[0]

    out: Dict[int, int] = {}
    for qi, ti in zip(*best.indices):
        if qi < 0 or ti < 0:
            continue
        out[ungapped[qi][0]] = columns[ti]
    return out
