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

A genotype is scored by **replacing the derived taxon's sequence** with it. It
is never added as an extra row: an extra row is an extra taxon, and the model
reads taxa as evidence of a realized evolutionary process, so an appended
genotype would inform how every site is embedded with evidence for a sequence
that may never have been viable. The alignment keeps exactly its own taxa and
the lineage under study occupies the row it already has.

That makes the alignment a different alignment, so its **distance matrix is
recomputed from it** -- TN93 for the host against every other taxon -- and its
MDS frame with it. Only the host's row and column can change, because every
other pair of sequences is untouched and TN93 is a function of the pair alone.

Because each genotype brings its own matrix, the per-site hidden states of one
genotype are not reusable for another, and every genotype costs one full
forward pass. Embeddings are memoised by mask.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch

from aeon_core.dataset import (
    CODON_TO_AA,
    compute_mds_coordinates,
    compute_tn93_cross_distance_matrix,
    compute_tn93_distance_matrix,
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
    derived_index: int               # the host row every genotype replaces
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
    # Patristic and TN93 distances are both substitutions per site, but they are
    # not the same matrix. A hypothetical genotype has no branch in any tree, so
    # its row can only be computed from its sequence; and the tree here is one
    # genome-wide tree shared by every gene alignment, while TN93 measures the
    # gene in front of you. On NR3C1 that is a ~7x gap (patristic median 0.76
    # against 0.11), the two correlate at r = 0.84, and no single factor
    # reconciles them (R^2 = 0.46 for the best scalar fit, the ratio falling
    # with distance as TN93 saturates). So the host's row cannot be put on the
    # tree's scale, and a patched row would read as far closer to every taxon
    # than any real sequence is.
    if not (use_tn93 or (tree is not None
                         and str(tree).strip().lower() in ("tn93", "none", "skip"))):
        raise SystemExit(
            "`order` needs --use-tn93. Each genotype's distances are recomputed "
            "from its sequence, and a tree has no branch for a sequence that is "
            "not in it. Patristic and TN93 distances are both substitutions per "
            "site, but a shared genome-wide tree and a single gene's TN93 are "
            "not interchangeable: on NR3C1 the patristic distances run ~7x the "
            "TN93 ones and no single factor reconciles them, so a TN93 row "
            "patched into a patristic matrix would place the genotype far "
            "closer to every taxon than any real sequence is. Pass --use-tn93 "
            "(the tree is then unused) and re-run."
        )

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

    # Every genotype replaces this row's sequence, so it has to be one of the
    # analysed taxa. A pruned duplicate resolves onto the identical sequence
    # that was kept, which is the same alignment row for this purpose.
    if d_name not in taxa:
        same = [t for t in taxa if seqs[t] == seqs[d_name]]
        if same:
            print(f"[*] --focal-taxon {d_name!r} is identical to {same[0]!r}, which "
                  "duplicate pruning kept; genotypes replace that row")
            d_name = same[0]
        else:
            raise SystemExit(
                f"--focal-taxon {d_name!r} is in the alignment but not among the "
                f"{len(taxa)} analysed taxa, and it is the row every genotype "
                "replaces. It is missing from the tree, so either add it there "
                "or run with --use-tn93."
            )

    keep = [t for t in (d_name, a_name) if t is not None and t in taxa]
    if max_species is not None and len(taxa) > max_species:
        from aeon_core.dataset import downsample_taxa_faith_pd
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

    if a_name is not None and a_name not in taxa:
        print(f"[*] --ancestor {a_name!r} supplies codons but is not one of the "
              f"{len(taxa)} analysed taxa, so it is not part of the background")

    device = get_device(cpu=cpu)
    model = load_model(weights=weights, variant=variant, device=device)
    tree_cache = model.precompute_tree_cache(d.to(device), z.to(device))
    tn93 = use_tn93 or (tree is not None and str(tree).strip().lower() in ("tn93", "none", "skip"))
    return Context(
        model=model, c_tensor=c, a_tensor=a,
        dist=dist,
        mds=z.squeeze(0).numpy().astype(np.float64),
        tree_cache=tree_cache, taxa=taxa, seqs=seqs, n_codons=L, device=device,
        derived_name=d_name, derived_index=taxa.index(d_name), ancestor_name=a_name,
        distance_source="tn93" if tn93 else "tree", variant=variant,
    )


# -- genotypes in their own row ---------------------------------------------------------

class GenotypeRows:
    """Each genotype replaces the derived taxon's sequence in the alignment.

    A genotype is not a taxon. Appending it as an extra row would hand the
    model one more piece of evidence about a realized evolutionary process --
    evidence for a sequence that may never have been viable -- and that
    evidence would then inform how every site is embedded. So the alignment
    keeps exactly its own taxa, and the lineage under study occupies the row it
    already has: the derived taxon's.

    Because that row's sequence changes, the alignment the model is given is a
    different alignment, and its **distance matrix is recomputed from it** --
    TN93 for the host against every other taxon, then MDS from the result. Only
    the host's row and column can change, since every other pair of sequences
    is untouched and TN93 is pairwise, so patching those entries is identical
    to recomputing the whole matrix (`check_matrix_patch_equals_full_recompute`
    verifies it).

    At the derived end of the branch the host row holds its original sequence,
    so the alignment, its matrix and its MDS frame are the originals, and the
    genotype's embedding is the derived taxon's own.
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
        self.host = ctx.derived_index
        self._others = [i for i in range(ctx.n_taxa) if i != self.host]
        self._other_names = [ctx.taxa[i] for i in self._others]
        self._dist_rows: Dict[int, np.ndarray] = {}
        self._embeddings: Dict[int, np.ndarray] = {}
        self.forward_passes = 0
        self.tn93_calls = 0
        self.real_row_hits = 0

    def real_row(self, mask: int) -> Optional[int]:
        """The host's index when this genotype is the host's own sequence.

        Only the host row counts. A genotype matching some other taxon is
        still a hypothetical sequence in the host's place, and the alignment it
        produces is not that taxon's alignment.
        """
        if not self.use_real_row:
            return None
        if self.nucleotides(mask) == self.ctx.seqs[self.ctx.derived_name]:
            return self.host
        return None

    # sequence ------------------------------------------------------------------
    def codons(self, mask: int) -> List[str]:
        out = list(self.ends.ancestral_codons)
        for sub in self.lattice.applied(mask):
            out[sub.index] = self.ends.derived_codons[sub.index]
        return out

    def nucleotides(self, mask: int) -> str:
        return "".join(self.codons(mask))

    # the modified alignment's distances ------------------------------------------
    def prefetch(self, masks: Iterable[int]) -> None:
        """TN93 distances for every listed genotype that lacks them, in one call."""
        todo = [
            m for m in dict.fromkeys(masks)
            if m not in self._dist_rows and self.real_row(m) is None
        ]
        if not todo:
            return
        ids = [f"__epistaeon_genotype_{m}" for m in todo]
        pool = {t: self.ctx.seqs[t] for t in self._other_names}
        pool.update({gid: self.nucleotides(m) for gid, m in zip(ids, todo)})
        dist = compute_tn93_cross_distance_matrix(pool, ids, self._other_names)
        self.tn93_calls += 1
        for m, row in zip(todo, dist):
            self._dist_rows[m] = np.asarray(row, dtype=np.float64)

    def distances(self, mask: int) -> np.ndarray:
        """TN93 distance from this genotype to each of the other taxa."""
        self.prefetch([mask])
        if mask not in self._dist_rows:            # the host's own sequence
            self._dist_rows[mask] = self.ctx.dist[self.host][self._others]
        return self._dist_rows[mask]

    def matrix(self, mask: int) -> np.ndarray:
        """The modified alignment's distance matrix.

        Only the host's row and column move: every other pair of sequences is
        unchanged and TN93 is a function of the pair alone.
        """
        row = self.distances(mask)
        dist = self.ctx.dist.copy()
        dist[self.host, self._others] = row
        dist[self._others, self.host] = row
        dist[self.host, self.host] = 0.0
        return dist

    def _assemble(self, mask: int, sites: Optional[Sequence[int]] = None):
        """Tensors for the modified alignment, plus the cache for its own matrix."""
        ctx = self.ctx
        if self.real_row(mask) is not None:
            dist, coords, tree_cache = ctx.dist, ctx.mds, ctx.tree_cache
        else:
            dist = self.matrix(mask)
            # MDS is recomputed from this alignment's own matrix, not projected
            coords = compute_mds_coordinates(dist.astype(np.float32), n_components=4)
            tree_cache = ctx.model.precompute_tree_cache(
                torch.from_numpy(dist.astype(np.float32)).unsqueeze(0).to(ctx.device),
                torch.from_numpy(np.asarray(coords, dtype=np.float32)).unsqueeze(0).to(ctx.device),
            )

        g_c, g_a = self._anc_c.clone(), self._anc_a.clone()
        for sub in self.lattice.applied(mask):
            g_c[sub.index], g_a[sub.index] = self._der[sub.index]
        cols = list(range(ctx.n_codons)) if sites is None else list(sites)
        c = ctx.c_tensor[cols].clone()
        a = ctx.a_tensor[cols].clone()
        c[:, self.host, 0] = g_c[cols]
        a[:, self.host, 0] = g_a[cols]
        return c, a, tree_cache

    def check_matrix_patch_equals_full_recompute(self, mask: int) -> float:
        """Largest disagreement between the patched matrix and a full recompute.

        `matrix` rewrites only the host's row and column on the argument that
        every other pair is untouched. This recomputes the whole thing from the
        modified alignment and returns the worst absolute difference, so the
        shortcut is checked rather than asserted.
        """
        seqs = {t: self.ctx.seqs[t] for t in self._other_names}
        seqs[self.ctx.derived_name] = self.nucleotides(mask)
        order = list(self.ctx.taxa)
        full = compute_tn93_distance_matrix(seqs, order)
        self.tn93_calls += 1
        return float(np.abs(np.asarray(full, dtype=np.float64) - self.matrix(mask)).max())

    def background_frame_shift(self) -> float:
        """How far recomputing the matrix moves a bystander taxon's embedding.

        The host's sequence is what the ordering is about, but changing it
        changes the matrix and so the MDS frame, which moves every taxon a
        little. Measured on a taxon that is neither endpoint, between the
        ancestral genotype and one intermediate: if it approaches the host's
        own per-unit delta, the signal is competing with the frame.
        """
        ctx = self.ctx
        probe = next(
            (i for i in range(ctx.n_taxa)
             if ctx.taxa[i] not in (ctx.ancestor_name, ctx.derived_name)),
            None,
        )
        other = next((m for m in range(1, self.lattice.full_mask + 1)
                      if self.real_row(m) is None), None)
        if probe is None or other is None:
            return 0.0
        out = []
        for mask in (0, other):
            c, a, tree_cache = self._assemble(mask)
            _attn, reprs = extract_cross_taxa_attentions_and_embeddings(
                ctx.model, c, a, tree_cache, device=ctx.device
            )
            self.forward_passes += 1
            out.append(np.asarray(reprs[probe], dtype=np.float64))
        return float(np.linalg.norm(out[0] - out[1]) / max(1e-12, np.linalg.norm(out[0])))

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
                self._embeddings[mask] = np.asarray(reprs[self.host], dtype=np.float64)
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
