"""
Epistaeon: mutational trajectory scoring between ancestral and derived sequences.

Two layers that must not be confused:

  scoring    supplies a background-dependent edge weight w(unit | genotype)
  trajectory infers orderings exactly over the 2^n subset lattice

The inference layer is exact but inert if the scorer is a plain exponential of
a state-function difference, because such weights give every ordering the same
path product. Breaking that is the scorer's job.
"""

__version__ = "0.1.0"

__all__ = [
    "Lattice", "Substitution", "derive_substitutions", "blocks_from_groups", "GAP",
    "SCORERS", "PottsMetropolis", "TypicalityGate", "SoftmaxRepaired", "DeltaPhiProduct",
    "pairwise_energy", "phi_divergence", "typicality_mahalanobis", "delta_phi_fn",
    "exact", "sample", "infer", "TrajectoryResult", "EXACT_UNIT_LIMIT",
    "__version__",
]

# The artifact verifier and evaluator intentionally run without model or
# numerical dependencies. Preserve the public API without importing inference
# modules just to inspect an already completed run.
_EXPORT_MODULES = {
    **{name: "lattice" for name in (
        "GAP", "Lattice", "Substitution", "blocks_from_groups", "derive_substitutions")},
    **{name: "scoring" for name in (
        "SCORERS", "DeltaPhiProduct", "PottsMetropolis", "SoftmaxRepaired",
        "TypicalityGate", "delta_phi_fn", "pairwise_energy", "phi_divergence",
        "typicality_mahalanobis")},
    **{name: "trajectory" for name in (
        "EXACT_UNIT_LIMIT", "TrajectoryResult", "exact", "infer", "sample")},
}


def __getattr__(name):
    if name not in _EXPORT_MODULES:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module
    value = getattr(import_module(f".{_EXPORT_MODULES[name]}", __name__), name)
    globals()[name] = value
    return value
