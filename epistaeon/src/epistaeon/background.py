"""
epistaeon/background.py
-----------------------
Endpoints and latent embeddings for hypothetical genotypes.

Both endpoints are rows of the alignment. The derived endpoint is the
`--focal-taxon` row. The ancestral endpoint is the `--ancestor` row when one is
named; otherwise each column starts from the commonest codon of its plurality
residue among the other rows, and phi is anchored at the model's [ROOT] token,
because there is no ancestral sequence to embed.

A substitution is a change of codon, so synonymous changes are ordered
alongside the rest and applying every unit reproduces the derived row exactly,
nucleotide for nucleotide. Ordering only the amino-acid changes would leave the
end of the trajectory short of the derived sequence at every synonymous column,
which a nucleotide distance then sees.

Every intermediate genotype is scored in its own row. The genotype's codons are
taken from the two endpoint rows, appended to the alignment as an extra taxon,
and placed among the others by TN93 distances to every real sequence. Its MDS
coordinates are projected into the existing frame (Gower's out-of-sample
formula) rather than recomputed, so the real taxa, their distances and the
[ROOT] origin are identical for every genotype -- the only thing that changes
between two forward passes is the genotype itself.

The endpoints are not appended, because they are not hypothetical: the
ancestral and derived sequences are rows of the alignment already, and they are
read from there rather than duplicated beside themselves. The consequence to
report is that an endpoint is embedded among N taxa and an intermediate among
N + 1; `GenotypeRows.endpoint_frame_shift` measures what that is worth.

Because each genotype brings its own distance row, the per-site hidden states
of one genotype are not reusable for another, and every genotype costs one
full forward pass. Embeddings are memoised by mask.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch

from aeon_core.dataset import (
    CODON_TO_AA,
    compute_tn93_cross_distance_matrix,
    get_aa_token,
    get_codon_token,
    load_alignment_and_tree,
    parse_alignment_sequences,
)
from aeon_core.inference import get_device, load_model
from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings

from .lattice import GAP, Lattice, Substitution

GAP_CODON = "---"


def residue_of(codon: str) -> str:
    """Amino acid for a codon; gaps, ambiguities and stops all read as GAP."""
    aa = CODON_TO_AA.get(codon.upper(), GAP)
    return GAP if aa == "*" else aa


def codons_of(seq: str, n_codons: int) -> List[str]:
    seq = seq[: 3 * n_codons].ljust(3 * n_codons, "-")
    return [seq[3 * i: 3 * i + 3].upper() for i in range(n_codons)]


# -- choosing rows -------------------------------------------------------------

def resolve_row(query: str, seqs: Dict[str, str], role: str) -> str:
    """The name of the alignment row `query` refers to.

    Exact name first, then case-insensitive exact, then a unique
    case-insensitive substring. An endpoint supplies codons only, so it is
    looked up in the whole alignment: it need not have survived duplicate
    pruning, tree matching or --max-species to be used as an endpoint.
    """
    names = list(seqs)
    hits = [n for n in names if n == query]
    if not hits:
        hits = [n for n in names if n.lower() == query.lower()]
    if not hits:
        hits = [n for n in names if query.lower() in n.lower()]
    if not hits:
        raise SystemExit(f"{role} {query!r} matches no sequence in the alignment")
    if len(hits) > 1:
        raise SystemExit(f"{role} {query!r} is ambiguous; it matches {hits[:8]}")
    return hits[0]


# -- endpoints -------------------------------------------------------------------

@dataclass
class Endpoints:
    """The two ends of the branch, in alignment-column coordinates."""

    ancestral_codons: List[str]          # mask-0 genotype, one codon per column
    derived_codons: Dict[int, str]       # column -> codon once its unit is applied
    ancestral_protein: str               # one residue (or GAP) per column
    substitutions: List[Substitution]
    reference_name: str                  # the row whose residues number positions
    residue_number: Dict[int, int]       # column -> 1-based residue in that row
    plurality_support: Dict[int, float] = field(default_factory=dict)
    classes: Dict[int, str] = field(default_factory=dict)   # column -> kind of change

    def column_of(self) -> Dict[int, int]:
        return {n: c for c, n in self.residue_number.items()}

    def class_counts(self) -> Dict[str, int]:
        return dict(Counter(self.classes.values()))


def _plurality(codons: Sequence[str]) -> Tuple[str, str, float]:
    """(residue, codon, support) of the most common residue in a column.

    Ties go to a residue over a gap, then alphabetically, so the choice is
    deterministic. The codon is the commonest one encoding that residue.
    """
    residues = Counter(residue_of(c) for c in codons)
    best = sorted(residues.items(), key=lambda kv: (-kv[1], kv[0] == GAP, kv[0]))[0][0]
    if best == GAP:
        codon = GAP_CODON
    else:
        codon = Counter(c for c in codons if residue_of(c) == best).most_common(1)[0][0]
    return best, codon, residues[best] / max(1, len(codons))


def derive_endpoints(
    rows: Dict[str, List[str]],
    derived: str,
    ancestor: Optional[str] = None,
    columns: Optional[Iterable[int]] = None,
) -> Endpoints:
    """Substitutions between the ancestral and derived rows.

    A substitution is any column whose codon differs, synonymous changes
    included, so applying every unit reproduces the derived row exactly.

    `rows` maps taxon -> codons per column. Without an ancestor, the ancestral
    state of each column is the commonest codon of the plurality residue across
    every row except the derived one. `columns` (0-based) restricts which
    columns may differ.
    """
    der = rows[derived]
    n = len(der)
    wanted = None if columns is None else set(columns)
    others = [t for t in rows if t != derived]

    if ancestor is not None:
        anc = list(rows[ancestor])
        support: Dict[int, float] = {}
    else:
        if not others:
            raise SystemExit("inferring ancestral states needs at least one row besides the derived one")
        anc, support = [], {}
        for col in range(n):
            _res, codon, frac = _plurality([rows[t][col] for t in others])
            anc.append(codon)
            support[col] = frac

    reference = rows[ancestor] if ancestor is not None else der
    residue_number: Dict[int, int] = {}
    k = 0
    for col, codon in enumerate(reference):
        if residue_of(codon) != GAP:
            k += 1
            residue_number[col] = k

    subs: List[Substitution] = []
    derived_codons: Dict[int, str] = {}
    classes: Dict[int, str] = {}
    for col in range(n):
        # a substitution is a change of codon, not of residue: a synonymous
        # change is a real event on the branch, and ordering codon changes is
        # what makes the end of the trajectory the derived row exactly
        if anc[col] == der[col] or (wanted is not None and col not in wanted):
            continue
        a, d = residue_of(anc[col]), residue_of(der[col])
        # positions are numbered in the reference row, which has no residue at
        # a column it is itself gapped at -- name those by the column instead
        num = residue_number.get(col)
        pos = str(num) if num is not None else f"@c{col + 1}"
        if a == GAP and d == GAP:
            # neither side resolves to a residue (ambiguity codes, or a partial
            # gap beside an indel), so there is no residue to name it by
            kind = "unresolved"
            label = f"c{col + 1}:{anc[col]}>{der[col]}"
        elif d == GAP:
            kind = "indel"
            label = f"{a}{pos}del"
        elif a == GAP:
            kind = "indel"
            label = f"ins{pos}{d}"
        elif a == d:
            kind = "synonymous"
            label = f"{a}{pos}={anc[col]}>{der[col]}"
        else:
            kind = "nonsynonymous"
            label = f"{a}{pos}{d}"
        subs.append(Substitution(index=col, ancestral=a, derived=d, label=label))
        derived_codons[col] = der[col]
        classes[col] = kind

    return Endpoints(
        ancestral_codons=anc,
        derived_codons=derived_codons,
        ancestral_protein="".join(residue_of(c) for c in anc),
        substitutions=subs,
        reference_name=ancestor if ancestor is not None else derived,
        residue_number=residue_number,
        plurality_support={s.index: support[s.index] for s in subs} if support else {},
        classes=classes,
    )


def resolve_positions(tokens: Iterable, endpoints: Endpoints) -> List[int]:
    """0-based columns for user positions.

    An integer is a 1-based residue number in the reference row (the ancestor
    when given, otherwise the derived row). `c<N>` is 1-based alignment column
    N, for columns the reference row has no residue at.
    """
    col_of = endpoints.column_of()
    out: List[int] = []
    for tok in tokens:
        s = str(tok).strip()
        if s.lower().startswith("c"):
            out.append(int(s[1:]) - 1)
            continue
        num = int(s)
        if num not in col_of:
            raise SystemExit(
                f"position {num} is beyond the {len(col_of)} residues of "
                f"{endpoints.reference_name!r}"
            )
        out.append(col_of[num])
    return out


# -- MDS placement ---------------------------------------------------------------

def project_mds(base_dist: np.ndarray, base_coords: np.ndarray, new_dist: np.ndarray) -> np.ndarray:
    """Place one new point in an existing classical-MDS frame (Gower 1968).

    The new point's inner products with the centred configuration are recovered
    from its squared distances, then regressed onto the existing coordinates.
    Re-projecting a point already in the configuration returns its coordinates.
    """
    D2 = np.asarray(base_dist, dtype=np.float64) ** 2
    d2 = np.asarray(new_dist, dtype=np.float64) ** 2
    b = -0.5 * (d2 - d2.mean() - D2.mean(axis=1) + D2.mean())
    X = np.asarray(base_coords, dtype=np.float64)
    coords, *_ = np.linalg.lstsq(X, b, rcond=None)
    return coords


# -- the loaded alignment -----------------------------------------------------------

@dataclass
class Context:
    """A loaded alignment, distances and model, plus the two endpoint rows."""

    model: torch.nn.Module
    c_tensor: torch.Tensor
    a_tensor: torch.Tensor
    dist: np.ndarray                 # [N, N] distances between real taxa
    mds: np.ndarray                  # [N, 4] their MDS coordinates
    tree_cache: dict
    taxa: List[str]
    seqs: Dict[str, str]             # nucleotide sequence per taxon, 3L long
    n_codons: int
    device: torch.device
    derived_name: str
    ancestor_name: Optional[str] = None
    distance_source: str = "tree"
    variant: Optional[str] = None
    _extant: Optional[np.ndarray] = field(default=None, repr=False)

    @property
    def n_taxa(self) -> int:
        return len(self.taxa)

    def rows(self) -> Dict[str, List[str]]:
        """Codons per column for the context taxa and both endpoint rows.

        An endpoint that is not one of the context taxa is included anyway: it
        contributes its codons without joining the background, and it is left
        out of the plurality vote.
        """
        wanted = list(self.taxa)
        for name in (self.derived_name, self.ancestor_name):
            if name is not None and name not in wanted:
                wanted.append(name)
        return {t: codons_of(self.seqs[t], self.n_codons) for t in wanted}


def load_context(
    alignment: str,
    tree: Optional[str],
    derived_taxon: str,
    ancestor_taxon: Optional[str] = None,
    weights: Optional[str] = None,
    variant: Optional[str] = None,
    cpu: bool = False,
    max_species: Optional[int] = None,
    use_tn93: bool = False,
) -> Context:
    """Load alignment + distances + model and locate the endpoint rows."""
    # Downsampling is applied here rather than inside the loader, because the
    # endpoints have to survive it: the derived row is where every genotype is
    # written, and Faith's PD picks divergent taxa, which drops exactly the
    # well-sampled ones (hg38 among them) that make natural focal taxa.
    c, a, d, z, _inv, taxa, L = load_alignment_and_tree(
        alignment, tree, max_species=None, prune_duplicates=True, use_tn93=use_tn93
    )
    raw = parse_alignment_sequences(alignment)
    seqs = {t: raw[t][: 3 * L].ljust(3 * L, "-") for t in raw}
    taxa = list(taxa)
    dist = d.squeeze(0).numpy().astype(np.float64)

    d_name = resolve_row(derived_taxon, seqs, "--focal-taxon")
    a_name = None
    if ancestor_taxon is not None:
        a_name = resolve_row(ancestor_taxon, seqs, "--ancestor")
        if a_name == d_name:
            raise SystemExit(f"--ancestor and --focal-taxon both resolve to {d_name!r}")

    # An endpoint that is one of the analysed taxa is read from its own row
    # rather than appended as a copy of itself, so it is worth resolving a
    # pruned duplicate onto the identical sequence that was kept.
    if d_name not in taxa:
        same = [t for t in taxa if seqs[t] == seqs[d_name]]
        if same:
            print(f"[*] --focal-taxon {d_name!r} is identical to {same[0]!r}, which "
                  "duplicate pruning kept; using that row as the derived endpoint")
            d_name = same[0]

    keep = [t for t in (d_name, a_name) if t is not None and t in taxa]
    if max_species is not None and len(taxa) > max_species:
        from aeon_core.dataset import compute_mds_coordinates, downsample_taxa_faith_pd
        pool = [i for i, t in enumerate(taxa) if t not in keep]
        budget = max(1, int(max_species) - len(keep))
        _sub, chosen = downsample_taxa_faith_pd(
            dist[np.ix_(pool, pool)], [taxa[i] for i in pool], budget
        )
        wanted = set(chosen) | set(keep)
        idx = [i for i, t in enumerate(taxa) if t in wanted]
        taxa = [taxa[i] for i in idx]
        dist = dist[np.ix_(idx, idx)]
        c, a = c[:, idx], a[:, idx]
        z = torch.tensor(compute_mds_coordinates(dist.astype(np.float32), n_components=4),
                         dtype=torch.float32).unsqueeze(0)
        d = torch.tensor(dist, dtype=torch.float32).unsqueeze(0)
        print(f"[*] Downsampling: kept {len(taxa)} taxa, including the endpoints "
              f"({', '.join(keep)})")

    for role, name in (("--focal-taxon", d_name), ("--ancestor", a_name)):
        if name is not None and name not in taxa:
            print(f"[*] {role} {name!r} supplies codons but is not one of the "
                  f"{len(taxa)} analysed taxa, so that endpoint is appended as "
                  "its own row like an intermediate rather than read from the "
                  "alignment")

    device = get_device(cpu=cpu)
    model = load_model(weights=weights, variant=variant, device=device)
    tree_cache = model.precompute_tree_cache(d.to(device), z.to(device))
    tn93 = use_tn93 or (tree is not None and str(tree).strip().lower() in ("tn93", "none", "skip"))
    return Context(
        model=model, c_tensor=c, a_tensor=a,
        dist=dist,
        mds=z.squeeze(0).numpy().astype(np.float64),
        tree_cache=tree_cache, taxa=taxa, seqs=seqs, n_codons=L, device=device,
        derived_name=d_name, ancestor_name=a_name,
        distance_source="tn93" if tn93 else "tree", variant=variant,
    )


# -- genotypes in their own row ---------------------------------------------------------

class GenotypeRows:
    """Intermediates get their own row; the endpoints use the rows they already have.

    An intermediate genotype is a sequence the alignment does not contain, so
    it is appended as an extra row and placed by TN93 distance to every taxon,
    with its MDS coordinates projected into the alignment's own frame so the
    real taxa and the [ROOT] origin never move between genotypes.

    The two endpoints are different: they *are* taxa. The ancestral and derived
    sequences already sit in the alignment, so they are read from there rather
    than appended a second time -- no genotype duplicates a real row.

    The cost of that is a frame difference the caller must report: an endpoint
    is embedded among N taxa and an intermediate among N + 1, and an extra row
    shifts every embedding slightly. `endpoint_frame_shift` measures it.
    """

    def __init__(
        self,
        ctx: Context,
        endpoints: Endpoints,
        lattice: Lattice,
        budget: Optional[int] = None,
        use_real_row: bool = True,
    ):
        self.ctx = ctx
        self.ends = endpoints
        self.lattice = lattice
        self.budget = budget
        self.use_real_row = use_real_row
        self._anc_c = torch.tensor([get_codon_token(x) for x in endpoints.ancestral_codons])
        self._anc_a = torch.tensor([get_aa_token(x) for x in endpoints.ancestral_codons])
        self._der = {
            col: (get_codon_token(x), get_aa_token(x))
            for col, x in endpoints.derived_codons.items()
        }
        self._dist_rows: Dict[int, np.ndarray] = {}
        self._embeddings: Dict[int, np.ndarray] = {}
        self.forward_passes = 0
        self.tn93_calls = 0
        self.real_row_hits = 0

    def real_row(self, mask: int) -> Optional[int]:
        """Index of the taxon this genotype *is*, if it is one of the endpoints.

        Only the endpoint taxa count. An intermediate that happens to match
        some other taxon is still a hypothetical genotype on this branch, and
        collapsing it onto that taxon would change what is being scored.
        """
        if not self.use_real_row:
            return None
        seq = self.nucleotides(mask)
        for name in (self.ctx.ancestor_name, self.ctx.derived_name):
            if name is not None and name in self.ctx.taxa and self.ctx.seqs[name] == seq:
                return self.ctx.taxa.index(name)
        return None

    # sequence ------------------------------------------------------------------
    def codons(self, mask: int) -> List[str]:
        out = list(self.ends.ancestral_codons)
        for sub in self.lattice.applied(mask):
            out[sub.index] = self.ends.derived_codons[sub.index]
        return out

    def nucleotides(self, mask: int) -> str:
        return "".join(self.codons(mask))

    # placement -----------------------------------------------------------------
    def prefetch(self, masks: Iterable[int]) -> None:
        """TN93 distances for every listed intermediate that lacks them, in one call."""
        todo = [
            m for m in dict.fromkeys(masks)
            if m not in self._dist_rows and self.real_row(m) is None
        ]
        if not todo:
            return
        ids = [f"__epistaeon_genotype_{m}" for m in todo]
        pool = {t: self.ctx.seqs[t] for t in self.ctx.taxa}
        pool.update({gid: self.nucleotides(m) for gid, m in zip(ids, todo)})
        dist = compute_tn93_cross_distance_matrix(pool, ids, self.ctx.taxa)
        self.tn93_calls += 1
        for m, row in zip(todo, dist):
            self._dist_rows[m] = np.asarray(row, dtype=np.float64)

    def distances(self, mask: int) -> np.ndarray:
        """TN93 distance from this genotype to every taxon."""
        self.prefetch([mask])
        return self._dist_rows[mask]

    def _assemble(self, mask: int, sites: Optional[Sequence[int]] = None):
        """Alignment tensors with this genotype appended as row N, plus its cache."""
        ctx = self.ctx
        row = self.distances(mask)
        n = ctx.n_taxa
        dist = np.zeros((n + 1, n + 1), dtype=np.float32)
        dist[:n, :n] = ctx.dist
        dist[n, :n] = dist[:n, n] = row
        coords = np.vstack([ctx.mds, project_mds(ctx.dist, ctx.mds, row)]).astype(np.float32)
        tree_cache = ctx.model.precompute_tree_cache(
            torch.from_numpy(dist).unsqueeze(0).to(ctx.device),
            torch.from_numpy(coords).unsqueeze(0).to(ctx.device),
        )

        g_c, g_a = self._anc_c.clone(), self._anc_a.clone()
        for sub in self.lattice.applied(mask):
            g_c[sub.index], g_a[sub.index] = self._der[sub.index]
        cols = list(range(ctx.n_codons)) if sites is None else list(sites)
        c = torch.cat([ctx.c_tensor[cols], g_c[cols].view(-1, 1, 1)], dim=1)
        a = torch.cat([ctx.a_tensor[cols], g_a[cols].view(-1, 1, 1)], dim=1)
        return c, a, tree_cache

    def endpoint_frame_shift(self) -> float:
        """How far an extra row moves an embedding, relative to its own norm.

        Measured on a taxon that is not an endpoint, with and without one
        intermediate appended, so it isolates the N vs N + 1 difference between
        how an endpoint and an intermediate are embedded.
        """
        ctx = self.ctx
        probe = next(
            (i for i in range(ctx.n_taxa)
             if ctx.taxa[i] not in (ctx.ancestor_name, ctx.derived_name)),
            0,
        )
        mask = next(
            (m for m in (1, self.lattice.full_mask ^ 1, self.lattice.full_mask)
             if self.real_row(m) is None),
            None,
        )
        if mask is None:
            return 0.0
        base = extant_embeddings(ctx)[probe]
        c, a, tree_cache = self._assemble(mask)
        _attn, reprs = extract_cross_taxa_attentions_and_embeddings(
            ctx.model, c, a, tree_cache, device=ctx.device
        )
        self.forward_passes += 1
        return float(np.linalg.norm(np.asarray(reprs[probe]) - base) / np.linalg.norm(base))

    # model outputs --------------------------------------------------------------
    def _check_budget(self, mask: int) -> None:
        """Refuse to start an unbounded run.

        Sampling can ask for up to one genotype per step per path, so without
        this a large lattice silently commits to days of forward passes.
        """
        if self.budget is None or mask in self._embeddings:
            return
        if len(self._embeddings) >= self.budget:
            raise SystemExit(
                f"this run has embedded {len(self._embeddings)} genotypes, its "
                f"--max-genotypes budget, and the scorer is asking for another. "
                "Each genotype costs one forward pass, so lower --n-samples, "
                "group units with --blocks, or raise --max-genotypes."
            )

    def embedding(self, mask: int) -> np.ndarray:
        """Mean-over-sites embedding of this genotype."""
        if mask not in self._embeddings:
            self._check_budget(mask)
            row = self.real_row(mask)
            if row is not None:
                # this genotype is an endpoint taxon: read it from the
                # alignment rather than appending a copy of it
                self.real_row_hits += 1
                self._embeddings[mask] = extant_embeddings(self.ctx)[row]
            else:
                c, a, tree_cache = self._assemble(mask)
                _attn, reprs = extract_cross_taxa_attentions_and_embeddings(
                    self.ctx.model, c, a, tree_cache, device=self.ctx.device
                )
                self.forward_passes += 1
                self._embeddings[mask] = np.asarray(reprs[self.ctx.n_taxa], dtype=np.float64)
        return self._embeddings[mask]

    def lrt(self, mask: int, sites: Sequence[int]) -> float:
        """Summed, floored predicted LRT at these columns with this genotype present.

        An endpoint is scored on the unmodified alignment, for the same reason
        its embedding is: it is already a row there.
        """
        if self.real_row(mask) is not None:
            cols = list(sites)
            c, a = self.ctx.c_tensor[cols], self.ctx.a_tensor[cols]
            tree_cache = self.ctx.tree_cache
        else:
            c, a, tree_cache = self._assemble(mask, sites=sites)
        self.ctx.model.eval()
        with torch.no_grad():
            y, *_ = self.ctx.model.forward_cached(
                c.to(self.ctx.device), a.to(self.ctx.device), tree_cache
            )
        self.forward_passes += 1
        return float(torch.clamp(y, min=0.0).sum().item())


def extant_embeddings(ctx: Context) -> np.ndarray:
    """Embeddings of every taxon under the unmodified alignment.

    The reference cloud for typicality scoring, the isometric calibration
    below, and the host row's own embedding at the derived end of the branch.
    Computed once per context and cached, since it does not depend on any
    genotype.
    """
    if ctx._extant is None:
        _attn, reprs = extract_cross_taxa_attentions_and_embeddings(
            ctx.model, ctx.c_tensor, ctx.a_tensor, ctx.tree_cache, device=ctx.device
        )
        ctx._extant = np.asarray(reprs, dtype=np.float64)
    return ctx._extant


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

    The phi anchor when no ancestor row is given. It lives in the same latent
    space as the taxon embeddings and needs no dates, unlike chronaeon's
    convex-hull root. Computed on the unmodified alignment, so it is the same
    reference for every genotype.
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
