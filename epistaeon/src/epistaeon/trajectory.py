"""
epistaeon/trajectory.py
-----------------------
Inference over mutational orderings on the subset lattice.

Exact mode runs three dynamic programmes over the 2^n subset DAG, all in log
space:

  Viterbi  V(S) = max_u [ V(S-u) + log w(u | S-u) ]        -> MAP ordering
  Forward  Z(S) = sum_u   Z(S-u) * w(u | S-u)              -> total path weight
  Backward B(S) = sum_u   w(u | S) * B(S + u)              -> completions

Edge marginals Z(S)*w(u|S)*B(S + u) / Z(M) then give every ordering statistic
without enumerating n! paths.

A caution worth stating in code: this machinery is exact but *inert* if the
weights are a plain exponential of a state-function difference. In that case
every ordering has identical weight, Z(M) = n! * w, and the DP correctly
reports a uniform distribution with C = 0.5 everywhere. Discrimination is the
scorer's job, not the DP's. `uniformity` in the result quantifies this.
"""

from dataclasses import dataclass, field
from math import inf, isinf, lgamma, log
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .lattice import Lattice

# log w(unit | mask); -inf means the step is blocked
LogWeightFn = Callable[[int, int], float]

EXACT_UNIT_LIMIT = 22


def _logsumexp(values: Sequence[float]) -> float:
    finite = [v for v in values if not isinf(v) or v > 0]
    if not finite:
        return -inf
    m = max(finite)
    if isinf(m):
        return m
    return m + log(sum(np.exp(v - m) for v in finite))


@dataclass
class TrajectoryResult:
    """Ordering statistics for one ancestor -> derived comparison."""

    mode: str
    n_units: int
    unit_names: List[str]
    map_order: List[int]
    map_log_prob: float
    log_total_weight: float
    map_prob_normalised: float
    before: np.ndarray                    # C[i, j] = P(unit i before unit j)
    position: np.ndarray                  # R[i, k] = P(unit i at step k)
    accessible: bool = True
    note: str = ""
    designated: Dict[str, float] = field(default_factory=dict)
    diagnostics: Dict[str, float] = field(default_factory=dict)

    @property
    def uniform_baseline(self) -> float:
        """1/n! — what P(MAP) looks like when the scorer does not discriminate."""
        return float(np.exp(-lgamma(self.n_units + 1)))

    @property
    def map_order_names(self) -> List[str]:
        return [self.unit_names[u] for u in self.map_order]

    def summary(self) -> str:
        if not self.accessible:
            return f"no accessible path ({self.note})"
        return (
            f"MAP {' -> '.join(self.map_order_names)} | "
            f"P(MAP|endpoint) = {self.map_prob_normalised:.4g} "
            f"(uniform baseline {self.uniform_baseline:.4g})"
        )


def _designated_marginals(
    edge_mass: Dict[Tuple[int, int], float],
    lattice: Lattice,
    focal: Optional[int],
    permissive: Optional[Sequence[int]],
) -> Dict[str, float]:
    """P(focal unit occurs after at least k of a permissive set), for all k.

    This is the statistic the receptor ground truth actually specifies: group Y
    is tolerated only once groups X and Z are present.
    """
    out: Dict[str, float] = {}
    if focal is None or not permissive:
        return out
    pset = set(int(p) for p in permissive)
    for k in range(len(pset) + 1):
        mass = 0.0
        for (mask, unit), m in edge_mass.items():
            if unit != focal:
                continue
            present = sum(1 for p in pset if mask >> p & 1)
            if present >= k:
                mass += m
        out[f"P(focal after >= {k} of permissive)"] = float(mass)
    return out


def exact(
    lattice: Lattice,
    log_weight: LogWeightFn,
    focal_unit: Optional[int] = None,
    permissive_units: Optional[Sequence[int]] = None,
    display_offset: int = 1,
) -> TrajectoryResult:
    """Exact subset DP. O(n * 2^n) time, O(2^n) memory."""
    n = lattice.n_units
    if n > EXACT_UNIT_LIMIT:
        raise ValueError(
            f"exact mode needs n <= {EXACT_UNIT_LIMIT} units (got {n}); "
            f"2^{n} nodes is {2 ** n:,}. Use blocks to group substitutions, "
            "or mode='sample'."
        )
    size = 1 << n
    full = size - 1
    neg = -inf

    # Ascending integer order is a valid topological order: clearing a set bit
    # always decreases the integer, so every predecessor is already computed.
    V = np.full(size, neg)
    Z = np.full(size, neg)
    back = np.full(size, -1, dtype=np.int64)
    V[0] = 0.0
    Z[0] = 0.0

    for mask in range(1, size):
        best = neg
        best_u = -1
        terms: List[float] = []
        for u in range(n):
            if not mask >> u & 1:
                continue
            prev = mask ^ (1 << u)
            if isinf(V[prev]) and V[prev] < 0 and isinf(Z[prev]) and Z[prev] < 0:
                continue
            lw = log_weight(prev, u)
            if isinf(lw) and lw < 0:
                continue
            cand = V[prev] + lw
            if cand > best:
                best, best_u = cand, u
            if not (isinf(Z[prev]) and Z[prev] < 0):
                terms.append(Z[prev] + lw)
        V[mask] = best
        back[mask] = best_u
        Z[mask] = _logsumexp(terms) if terms else neg

    names = [lattice.unit_name(u, display_offset) for u in range(n)]
    if isinf(Z[full]) and Z[full] < 0:
        return TrajectoryResult(
            mode="exact", n_units=n, unit_names=names, map_order=[],
            map_log_prob=neg, log_total_weight=neg, map_prob_normalised=0.0,
            before=np.full((n, n), np.nan), position=np.full((n, n), np.nan),
            accessible=False,
            note="every ordering is blocked by a zero-weight step",
        )

    # backward pass
    B = np.full(size, neg)
    B[full] = 0.0
    for mask in range(size - 2, -1, -1):
        terms = []
        for u in range(n):
            if mask >> u & 1:
                continue
            nxt = mask | (1 << u)
            if isinf(B[nxt]) and B[nxt] < 0:
                continue
            lw = log_weight(mask, u)
            if isinf(lw) and lw < 0:
                continue
            terms.append(lw + B[nxt])
        B[mask] = _logsumexp(terms) if terms else neg

    # edge marginals
    logZM = Z[full]
    C = np.zeros((n, n))
    R = np.zeros((n, n))
    edge_mass: Dict[Tuple[int, int], float] = {}
    for mask in range(size):
        if isinf(Z[mask]) and Z[mask] < 0:
            continue
        step = bin(mask).count("1")
        for u in range(n):
            if mask >> u & 1:
                continue
            nxt = mask | (1 << u)
            if isinf(B[nxt]) and B[nxt] < 0:
                continue
            lw = log_weight(mask, u)
            if isinf(lw) and lw < 0:
                continue
            mass = float(np.exp(Z[mask] + lw + B[nxt] - logZM))
            if mass <= 0.0:
                continue
            edge_mass[(mask, u)] = mass
            R[u, step] += mass
            for i in range(n):
                if mask >> i & 1:
                    C[i, u] += mass

    order: List[int] = []
    mask = full
    while mask:
        u = int(back[mask])
        if u < 0:
            break
        order.append(u)
        mask ^= 1 << u
    order.reverse()

    map_logp = float(V[full])
    return TrajectoryResult(
        mode="exact", n_units=n, unit_names=names, map_order=order,
        map_log_prob=map_logp, log_total_weight=float(logZM),
        map_prob_normalised=float(np.exp(map_logp - logZM)),
        before=C, position=R,
        designated=_designated_marginals(edge_mass, lattice, focal_unit, permissive_units),
        diagnostics={"n_edges_scored": float(len(edge_mass))},
    )


def sample(
    lattice: Lattice,
    log_weight: LogWeightFn,
    n_samples: int = 20000,
    seed: int = 0,
    focal_unit: Optional[int] = None,
    permissive_units: Optional[Sequence[int]] = None,
    display_offset: int = 1,
    beam_width: int = 64,
) -> TrajectoryResult:
    """Importance-sampled ordering statistics, for n beyond the exact limit.

    Orderings are drawn sequentially with probability proportional to the step
    weights, which is *not* the path-normalised distribution. Each draw is
    therefore reweighted by the product of its step denominators, making the
    estimates unbiased for P(pi) ∝ prod w. Standard errors are reported in
    diagnostics.
    """
    rng = np.random.default_rng(seed)
    n = lattice.n_units
    names = [lattice.unit_name(u, display_offset) for u in range(n)]

    C = np.zeros((n, n))
    R = np.zeros((n, n))
    des_counts: Dict[int, float] = {}
    wsum = 0.0
    wsq = 0.0
    best_logp = -inf
    best_order: List[int] = []
    pset = set(int(p) for p in (permissive_units or []))

    for _ in range(int(n_samples)):
        mask = 0
        order: List[int] = []
        log_p = 0.0        # log prod w   (path weight)
        log_q = 0.0        # log prod q   (proposal)
        blocked = False
        for _step in range(n):
            avail = [u for u in range(n) if not mask >> u & 1]
            lws = np.array([log_weight(mask, u) for u in avail])
            finite = np.isfinite(lws)
            if not finite.any():
                blocked = True
                break
            m = lws[finite].max()
            probs = np.where(finite, np.exp(lws - m), 0.0)
            denom = probs.sum()
            probs = probs / denom
            j = int(rng.choice(len(avail), p=probs))
            u = avail[j]
            log_p += float(lws[j])
            log_q += float(np.log(probs[j]))
            if focal_unit is not None and u == focal_unit and pset:
                present = sum(1 for p in pset if mask >> p & 1)
                des_counts[present] = des_counts.get(present, 0.0)
            order.append(u)
            mask |= 1 << u
        if blocked:
            continue
        iw = float(np.exp(log_p - log_q))
        wsum += iw
        wsq += iw * iw
        if log_p > best_logp:
            best_logp, best_order = log_p, list(order)
        pos_of = {u: k for k, u in enumerate(order)}
        for i in range(n):
            R[i, pos_of[i]] += iw
            for j in range(n):
                if i != j and pos_of[i] < pos_of[j]:
                    C[i, j] += iw
        if focal_unit is not None and pset:
            before_focal = {u for u in order[: pos_of[focal_unit]]}
            present = len(pset & before_focal)
            des_counts[present] = des_counts.get(present, 0.0) + iw

    if wsum <= 0.0:
        return TrajectoryResult(
            mode="sample", n_units=n, unit_names=names, map_order=[],
            map_log_prob=-inf, log_total_weight=-inf, map_prob_normalised=0.0,
            before=np.full((n, n), np.nan), position=np.full((n, n), np.nan),
            accessible=False, note="no unblocked ordering was sampled",
        )

    C /= wsum
    R /= wsum
    designated: Dict[str, float] = {}
    if focal_unit is not None and pset:
        for k in range(len(pset) + 1):
            designated[f"P(focal after >= {k} of permissive)"] = float(
                sum(m for present, m in des_counts.items() if present >= k) / wsum
            )
    # effective sample size of the importance weights
    ess = (wsum * wsum) / wsq if wsq > 0 else 0.0
    # log Z estimate = log(mean importance weight) + log(n!) is not formed here;
    # report the sampled mean weight so P(MAP) can be scaled consistently.
    log_mean_w = float(np.log(wsum / max(1, int(n_samples))))
    return TrajectoryResult(
        mode="sample", n_units=n, unit_names=names, map_order=best_order,
        map_log_prob=float(best_logp), log_total_weight=log_mean_w,
        map_prob_normalised=float(np.exp(best_logp - log_mean_w)) if np.isfinite(best_logp) else 0.0,
        before=C, position=R, designated=designated,
        diagnostics={
            "n_samples": float(n_samples),
            "effective_sample_size": float(ess),
            "note_map_prob_is_estimate": 1.0,
        },
    )


def infer(
    lattice: Lattice,
    log_weight: LogWeightFn,
    mode: str = "auto",
    **kwargs,
) -> TrajectoryResult:
    """Run exact or sampled inference, choosing by size when mode='auto'."""
    if mode == "auto":
        mode = "exact" if lattice.n_units <= EXACT_UNIT_LIMIT else "sample"
    if mode == "exact":
        return exact(lattice, log_weight, **{
            k: v for k, v in kwargs.items()
            if k in {"focal_unit", "permissive_units", "display_offset"}
        })
    if mode == "sample":
        return sample(lattice, log_weight, **kwargs)
    raise ValueError(f"unknown mode {mode!r}; expected 'exact', 'sample' or 'auto'")
