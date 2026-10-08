"""
epistaeon/lattice.py
--------------------
The mutational lattice between an ancestral and a derived sequence.

A genotype is an integer bitmask over the K substitutions that separate the two
endpoints: bit k set means substitution k has been applied. Bit 0 (mask 0) is
the ancestor, and (1 << K) - 1 is the derived sequence.

Nothing here materialises the 2^K graph. Nodes are integers and edges are
generated on demand, so the same representation serves the exact subset DP (K
small) and sampling (K large).
"""

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

GAP = "-"


@dataclass(frozen=True)
class Substitution:
    """One difference between the ancestral and derived sequence.

    `index` is the position in the ungapped ancestral sequence, 0-based.
    `ancestral` and `derived` are single-character states; either may be GAP,
    which is how an indel is represented (e.g. the receptor's S212 deletion).
    """

    index: int
    ancestral: str
    derived: str
    label: Optional[str] = None

    @property
    def is_indel(self) -> bool:
        return self.ancestral == GAP or self.derived == GAP

    def name(self, offset: int = 1) -> str:
        """Human-readable name, e.g. 'S106P' or 'S212del'.

        `offset` converts the 0-based index to the numbering used for display;
        the default of 1 gives 1-based positions. Callers working in a
        structure's numbering should pass the offset from
        data/numbering_offsets.json rather than assuming a constant.
        """
        if self.label:
            return self.label
        pos = self.index + offset
        if self.derived == GAP:
            return f"{self.ancestral}{pos}del"
        if self.ancestral == GAP:
            return f"ins{pos}{self.derived}"
        return f"{self.ancestral}{pos}{self.derived}"


def derive_substitutions(
    ancestral: str,
    derived: str,
    sites: Optional[Iterable[int]] = None,
) -> List[Substitution]:
    """Find the substitutions separating two aligned, equal-length sequences.

    If `sites` is given (0-based indices), only those positions are considered,
    which is how the --sites flag restricts K to an assayed subset.
    """
    if len(ancestral) != len(derived):
        raise ValueError(
            f"ancestral and derived must be the same length "
            f"({len(ancestral)} != {len(derived)}); align them first"
        )
    wanted = None if sites is None else set(int(s) for s in sites)
    subs: List[Substitution] = []
    for i, (a, d) in enumerate(zip(ancestral, derived)):
        if a == d:
            continue
        if wanted is not None and i not in wanted:
            continue
        subs.append(Substitution(index=i, ancestral=a, derived=d))
    return subs


@dataclass
class Lattice:
    """The set of substitutions, plus genotype <-> sequence conversion.

    `units` are what the DP orders. Normally one unit per substitution; in
    block mode a unit is a group of substitutions applied together, which is
    what keeps the exact DP feasible when K is large.
    """

    ancestral: str
    substitutions: List[Substitution]
    blocks: Optional[List[List[int]]] = None
    block_names: Optional[List[str]] = None
    _unit_members: List[List[int]] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        n_subs = len(self.substitutions)
        if self.blocks is None:
            self._unit_members = [[i] for i in range(n_subs)]
        else:
            seen: set = set()
            for group in self.blocks:
                for i in group:
                    if not 0 <= i < n_subs:
                        raise ValueError(f"block member {i} out of range (K={n_subs})")
                    if i in seen:
                        raise ValueError(f"substitution {i} appears in more than one block")
                    seen.add(i)
            missing = sorted(set(range(n_subs)) - seen)
            if missing:
                raise ValueError(
                    f"blocks must cover every substitution; missing {missing}. "
                    "Add them as singleton blocks if they should be ordered individually."
                )
            self._unit_members = [list(g) for g in self.blocks]
        if self.block_names is not None and len(self.block_names) != len(self._unit_members):
            raise ValueError("block_names length must match the number of blocks")

    # -- sizes ---------------------------------------------------------------
    @property
    def n_units(self) -> int:
        return len(self._unit_members)

    @property
    def n_substitutions(self) -> int:
        return len(self.substitutions)

    @property
    def full_mask(self) -> int:
        return (1 << self.n_units) - 1

    def unit_name(self, u: int, offset: int = 1) -> str:
        if self.block_names is not None:
            return self.block_names[u]
        members = self._unit_members[u]
        if len(members) == 1:
            return self.substitutions[members[0]].name(offset)
        return "+".join(self.substitutions[i].name(offset) for i in members)

    def unit_members(self, u: int) -> List[Substitution]:
        return [self.substitutions[i] for i in self._unit_members[u]]

    # -- genotypes -----------------------------------------------------------
    def applied(self, mask: int) -> List[Substitution]:
        """Every substitution applied in this genotype."""
        out: List[Substitution] = []
        for u in range(self.n_units):
            if mask >> u & 1:
                out.extend(self.unit_members(u))
        return out

    def sequence(self, mask: int) -> str:
        """The genotype's sequence, with GAP written as '-' (not deleted)."""
        chars = list(self.ancestral)
        for sub in self.applied(mask):
            if chars[sub.index] != sub.ancestral:
                raise ValueError(
                    f"position {sub.index} holds {chars[sub.index]!r}, expected "
                    f"{sub.ancestral!r}; ancestral sequence and substitutions disagree"
                )
            chars[sub.index] = sub.derived
        return "".join(chars)

    def available(self, mask: int) -> List[int]:
        """Units not yet applied — the out-edges of this node."""
        return [u for u in range(self.n_units) if not mask >> u & 1]

    def edges(self, mask: int) -> List[Tuple[int, int]]:
        """(unit, next_mask) for each out-edge."""
        return [(u, mask | (1 << u)) for u in self.available(mask)]

    def step_index(self, mask: int) -> int:
        """How many units are applied. In this DAG that *is* the step index."""
        return bin(mask).count("1")

    def unit_index_by_site(self) -> Dict[int, int]:
        """Map sequence position -> unit, for designating permissive sets."""
        out: Dict[int, int] = {}
        for u in range(self.n_units):
            for sub in self.unit_members(u):
                out[sub.index] = u
        return out


def blocks_from_groups(
    substitutions: Sequence[Substitution],
    groups: Dict[str, Sequence[int]],
) -> Tuple[List[List[int]], List[str]]:
    """Build blocks from named groups of sequence positions.

    `groups` maps a name to ungapped sequence indices, e.g.
    {"X": [105, 110], "Y": [28, 97, 211], "Z": [25, 104]}. Substitutions not
    named in any group become singleton blocks, so nothing is silently dropped.
    """
    by_index = {sub.index: i for i, sub in enumerate(substitutions)}
    blocks: List[List[int]] = []
    names: List[str] = []
    claimed: set = set()
    for name, indices in groups.items():
        members = []
        for idx in indices:
            if idx not in by_index:
                raise ValueError(
                    f"group {name!r} names position {idx}, which is not a substitution "
                    "between these two sequences"
                )
            members.append(by_index[idx])
            claimed.add(by_index[idx])
        if members:
            blocks.append(sorted(members))
            names.append(name)
    for i, sub in enumerate(substitutions):
        if i not in claimed:
            blocks.append([i])
            names.append(sub.name())
    return blocks, names
