"""
epistaeon/scoring.py
--------------------
Edge weights for the mutational lattice.

A scorer answers one question: in log space, how favoured is applying unit `u`
to genotype `mask`? Weights need not be normalised per step -- the inference
layer normalises over whole paths via Z(M) -- so a scorer is free to return any
non-positive-infinite log weight, including -inf for a blocked step.

The one property that matters is that the weights must NOT be a plain
exponential of a state-function difference. For any state function f,
sum of (f(next) - f(cur)) along a path telescopes to f(derived) - f(ancestor),
which is identical for every ordering; the DP then correctly reports a uniform
distribution. Every scorer here therefore breaks that in an explicit way:

  potts-metropolis  min(0, dE/T)       -- the min() is the non-linearity
  typicality-gate   dV/T + log gate    -- absolute level of the new state
  softmax-repaired  per-step softmax   -- varying normalisation
  dphi-product      log(dphi/T)        -- baseline only; see the class docstring
"""

from dataclasses import dataclass
from math import inf, log
from typing import Callable, Dict, Optional, Protocol

import numpy as np

from .lattice import Lattice

StateEnergy = Callable[[int], float]


class Scorer(Protocol):
    """Anything that can weight a lattice edge in log space."""

    name: str

    def log_weight(self, mask: int, unit: int) -> float: ...


@dataclass
class PottsMetropolis:
    """Default. Pairwise-coupled energy, Metropolis acceptance.

    `energy(mask)` should be higher for more viable genotypes. With site terms
    plus pairwise couplings the increment dE depends on the background, which
    is where the epistasis comes from; the Metropolis min() then makes a
    deleterious step cost something no matter where it falls in the order.

    Honest limit: with the published checkpoint the couplings come from the
    co-selection network, which is a correlation statistic, not a fitted
    coupling. This is not learned site-site coupling.
    """

    energy: StateEnergy
    temperature: float = 0.5
    name: str = "potts-metropolis"

    def log_weight(self, mask: int, unit: int) -> float:
        d = self.energy(mask | (1 << unit)) - self.energy(mask)
        return min(0.0, d / self.temperature)


@dataclass
class TypicalityGate:
    """Score the absolute plausibility of the state being entered.

    `viability(mask)` is high when the genotype's embedding looks like a real
    member of the family (e.g. negative Mahalanobis distance to the cloud of
    extant embeddings). Edges are weighted by the change in viability plus a
    gate on the new state's absolute level, so an ordering that must pass
    through an atypical intermediate is penalised directly rather than only
    through a difference.
    """

    viability: StateEnergy
    temperature: float = 0.5
    floor: Optional[float] = None
    name: str = "typicality-gate"

    def log_weight(self, mask: int, unit: int) -> float:
        nxt = mask | (1 << unit)
        v_new = self.viability(nxt)
        if self.floor is not None and v_new < self.floor:
            return -inf                      # intermediate judged non-viable
        d = v_new - self.viability(mask)
        return d / self.temperature + v_new / self.temperature


@dataclass
class SoftmaxRepaired:
    """The white paper's rule, with both defects fixed.

    Normalises over *all* single-residue moves available at each step, not just
    the remaining historical substitutions -- otherwise the final step has only
    one option and contributes probability 1 regardless of its score. An
    optional viability gate supplies the absolute term the paper lacks.

    `alternatives(mask, unit)` returns the log scores of the other moves that
    compete with this one at this step.
    """

    delta_phi: Callable[[int, int], float]
    alternatives: Callable[[int, int], np.ndarray]
    temperature: float = 0.5
    gate: Optional[Callable[[int], bool]] = None
    name: str = "softmax-repaired"

    def log_weight(self, mask: int, unit: int) -> float:
        nxt = mask | (1 << unit)
        if self.gate is not None and not self.gate(nxt):
            return -inf
        mine = self.delta_phi(mask, unit) / self.temperature
        others = np.asarray(self.alternatives(mask, unit), dtype=float) / self.temperature
        pool = np.concatenate([[mine], others])
        m = pool.max()
        return float(mine - (m + np.log(np.exp(pool - m).sum())))


@dataclass
class DeltaPhiProduct:
    """Baseline only. Raw dphi/T as a multiplicative edge weight.

    Documented so its two pathologies are visible rather than discovered later:

    1. `temperature` is inert. Every ordering has exactly n edges, so dividing
       each by T scales every path by T^-n and cannot change the ranking.
    2. What it rewards is evenly sized steps, not plausibility: for a fixed sum
       (which is path-independent here) a product is maximised when the terms
       are equal.

    Non-positive increments have no logarithm, so they are blocked outright.
    """

    delta_phi: Callable[[int, int], float]
    temperature: float = 0.5
    name: str = "dphi-product"

    def log_weight(self, mask: int, unit: int) -> float:
        d = self.delta_phi(mask, unit) / self.temperature
        if d <= 0.0:
            return -inf
        return log(d)


def pairwise_energy(
    lattice: Lattice,
    site_terms: Dict[int, float],
    couplings: Dict[tuple, float],
    scale: float = 1.0,
) -> StateEnergy:
    """Build E(S) = sum_u f_u + sum_{u<v} J_uv over applied units.

    `site_terms` maps unit -> f, `couplings` maps (unit_i, unit_j) -> J with
    i < j. Returned energy is higher for more viable genotypes; `scale` flips
    or rescales if the supplied statistics run the other way.
    """
    n = lattice.n_units
    f = np.zeros(n)
    for u, v in site_terms.items():
        f[int(u)] = float(v)
    J = np.zeros((n, n))
    for (i, j), v in couplings.items():
        i, j = int(i), int(j)
        J[i, j] = J[j, i] = float(v)

    def energy(mask: int) -> float:
        x = np.fromiter(((mask >> u) & 1 for u in range(n)), dtype=float, count=n)
        return float(scale * (f @ x + x @ J @ x / 2.0))

    return energy


SCORERS = {
    "potts-metropolis": PottsMetropolis,
    "typicality-gate": TypicalityGate,
    "softmax-repaired": SoftmaxRepaired,
    "dphi-product": DeltaPhiProduct,
}


# ---------------------------------------------------------------------------
# State functions built from latent embeddings
# ---------------------------------------------------------------------------

def phi_divergence(
    embedding_of: Callable[[int], np.ndarray],
    anchor: np.ndarray,
    alpha: float = 1.0,
) -> StateEnergy:
    """phi(S) = alpha * || z(S) - anchor ||, calibrated divergence from the anchor.

    Note this is a state function, so differences of it telescope: on its own
    it cannot rank orderings. It is useful as a coordinate (for relative branch
    position) and as an input to scorers whose combination rule is non-linear.
    """
    anchor = np.asarray(anchor, dtype=np.float64)

    def phi(mask: int) -> float:
        return float(alpha) * float(np.linalg.norm(embedding_of(mask) - anchor))

    return phi


def typicality_mahalanobis(
    embedding_of: Callable[[int], np.ndarray],
    cloud: np.ndarray,
    shrinkage: float = 1e-3,
) -> StateEnergy:
    """Viability as negative Mahalanobis distance to the extant embedding cloud.

    High when a genotype's embedding sits where real members of the family sit.
    The covariance is shrunk toward its diagonal because the cloud (hundreds of
    taxa) is usually smaller than the latent dimension (384).
    """
    cloud = np.asarray(cloud, dtype=np.float64)
    mu = cloud.mean(axis=0)
    cov = np.cov(cloud, rowvar=False)
    cov = cov + np.eye(cov.shape[0]) * (shrinkage * float(np.trace(cov)) / cov.shape[0] + 1e-12)
    inv = np.linalg.pinv(cov)

    def viability(mask: int) -> float:
        d = embedding_of(mask) - mu
        return -float(np.sqrt(max(0.0, d @ inv @ d)))

    return viability


def delta_phi_fn(phi: StateEnergy) -> Callable[[int, int], float]:
    """Edge-wise increment of a state function."""

    def delta(mask: int, unit: int) -> float:
        return phi(mask | (1 << unit)) - phi(mask)

    return delta
