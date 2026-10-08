"""
Epistaeon: mutational trajectory scoring between ancestral and derived sequences.

Two layers that must not be confused:

  scoring    supplies a background-dependent edge weight w(unit | genotype)
  trajectory infers orderings exactly over the 2^n subset lattice

The inference layer is exact but inert if the scorer is a plain exponential of
a state-function difference, because such weights give every ordering the same
path product. Breaking that is the scorer's job.
"""

from .lattice import GAP, Lattice, Substitution, blocks_from_groups, derive_substitutions
from .scoring import (
    SCORERS,
    DeltaPhiProduct,
    PottsMetropolis,
    SoftmaxRepaired,
    TypicalityGate,
    delta_phi_fn,
    pairwise_energy,
    phi_divergence,
    typicality_mahalanobis,
)
from .trajectory import EXACT_UNIT_LIMIT, TrajectoryResult, exact, infer, sample

__version__ = "0.1.0"

__all__ = [
    "Lattice", "Substitution", "derive_substitutions", "blocks_from_groups", "GAP",
    "SCORERS", "PottsMetropolis", "TypicalityGate", "SoftmaxRepaired", "DeltaPhiProduct",
    "pairwise_energy", "phi_divergence", "typicality_mahalanobis", "delta_phi_fn",
    "exact", "sample", "infer", "TrajectoryResult", "EXACT_UNIT_LIMIT",
    "__version__",
]
