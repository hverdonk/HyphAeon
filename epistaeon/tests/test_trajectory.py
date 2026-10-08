"""Exact DP must agree with brute-force enumeration, and must stay honest
about scorers that cannot discriminate between orderings."""

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon.lattice import Lattice, derive_substitutions  # noqa: E402
from epistaeon import trajectory as T  # noqa: E402


def _toy(n, seed=3):
    """A potts-like state function with pairwise couplings, so that
    Metropolis weights genuinely depend on the background."""
    rng = np.random.default_rng(seed)
    h = rng.normal(size=n)
    J = rng.normal(size=(n, n)) * 1.2
    J = (J + J.T) / 2
    np.fill_diagonal(J, 0.0)

    def energy(mask):
        x = np.array([(mask >> i) & 1 for i in range(n)], dtype=float)
        return float(h @ x + x @ J @ x / 2)

    return energy


def _lattice(n):
    anc = "A" * n
    der = "C" * n
    return Lattice(anc, derive_substitutions(anc, der))


def _brute_force(lat, logw):
    """Every ordering, explicitly."""
    n = lat.n_units
    out = {}
    for order in itertools.permutations(range(n)):
        mask, lp = 0, 0.0
        for u in order:
            lp += logw(mask, u)
            mask |= 1 << u
        out[order] = lp
    return out


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_exact_dp_matches_brute_force(n):
    lat = _lattice(n)
    energy = _toy(n)
    T_sel = 0.5

    def logw(mask, u):  # metropolis, in log space
        d = energy(mask | (1 << u)) - energy(mask)
        return min(0.0, d / T_sel)

    res = T.exact(lat, logw)
    paths = _brute_force(lat, logw)

    # total weight
    expected_logZ = math.log(sum(math.exp(lp) for lp in paths.values()))
    assert res.log_total_weight == pytest.approx(expected_logZ, rel=1e-9)

    # MAP ordering and its probability. Metropolis clamps every uphill step to
    # log w = 0, so several orderings can tie exactly; any argmax is correct.
    best_lp = max(paths.values())
    optima = {o for o, lp in paths.items() if abs(lp - best_lp) < 1e-12}
    assert res.map_log_prob == pytest.approx(best_lp, rel=1e-9)
    assert tuple(res.map_order) in optima
    assert res.map_prob_normalised == pytest.approx(math.exp(best_lp - expected_logZ), rel=1e-9)

    # pairwise and position marginals
    C = np.zeros((n, n))
    R = np.zeros((n, n))
    for order, lp in paths.items():
        p = math.exp(lp - expected_logZ)
        pos = {u: k for k, u in enumerate(order)}
        for i in range(n):
            R[i, pos[i]] += p
            for j in range(n):
                if i != j and pos[i] < pos[j]:
                    C[i, j] += p
    assert np.allclose(res.before, C, atol=1e-10)
    assert np.allclose(res.position, R, atol=1e-10)


@pytest.mark.parametrize("n", [3, 4, 5])
def test_marginals_are_proper_distributions(n):
    lat = _lattice(n)
    energy = _toy(n)

    def logw(mask, u):
        return min(0.0, (energy(mask | (1 << u)) - energy(mask)) / 0.5)

    res = T.exact(lat, logw)
    assert np.allclose(res.position.sum(axis=1), 1.0)   # each unit occupies one step
    assert np.allclose(res.position.sum(axis=0), 1.0)   # each step holds one unit
    off = ~np.eye(n, dtype=bool)
    assert np.allclose((res.before + res.before.T)[off], 1.0)


def test_inert_scorer_yields_uniform_distribution():
    """A plain exponential of a state-function difference cannot discriminate.

    This is the regression for the original design error: the DP is correct but
    inert, so the scorer must be the thing that breaks path-independence.
    """
    n = 4
    lat = _lattice(n)
    energy = _toy(n)

    def telescoping(mask, u):           # log exp(dE/T) = dE/T
        return (energy(mask | (1 << u)) - energy(mask)) / 0.5

    res = T.exact(lat, telescoping)
    assert res.map_prob_normalised == pytest.approx(res.uniform_baseline, rel=1e-9)
    off = ~np.eye(n, dtype=bool)
    assert np.allclose(res.before[off], 0.5, atol=1e-9)

    def metropolis(mask, u):            # non-telescoping
        return min(0.0, telescoping(mask, u))

    disc = T.exact(lat, metropolis)
    assert disc.map_prob_normalised > 2 * disc.uniform_baseline
    assert not np.allclose(disc.before[off], 0.5, atol=0.05)


def test_blocked_edges_do_not_produce_nan():
    n = 3
    lat = _lattice(n)

    def logw(mask, u):
        if u == 1 and mask == 0:        # unit 1 may never go first
            return -math.inf
        return 0.0

    res = T.exact(lat, logw)
    assert res.accessible
    assert np.isfinite(res.before).all()
    assert res.position[1, 0] == pytest.approx(0.0)


def test_all_paths_blocked_is_reported():
    lat = _lattice(3)
    res = T.exact(lat, lambda mask, u: -math.inf)
    assert not res.accessible
    assert "blocked" in res.note


def test_exact_mode_refuses_oversized_lattices():
    lat = _lattice(T.EXACT_UNIT_LIMIT + 1)
    with pytest.raises(ValueError, match="exact mode needs"):
        T.exact(lat, lambda mask, u: 0.0)


def test_sampling_approximates_exact():
    n = 5
    lat = _lattice(n)
    energy = _toy(n)

    def logw(mask, u):
        return min(0.0, (energy(mask | (1 << u)) - energy(mask)) / 0.5)

    ex = T.exact(lat, logw)
    sm = T.sample(lat, logw, n_samples=40000, seed=1)
    off = ~np.eye(n, dtype=bool)
    assert np.abs(sm.before[off] - ex.before[off]).max() < 0.05
    assert sm.diagnostics["effective_sample_size"] > 100


def test_designated_marginal_matches_brute_force():
    """P(focal after >= k of a permissive set) -- the receptor's statistic."""
    n = 4
    lat = _lattice(n)
    energy = _toy(n)
    focal, perm = 3, [0, 1]

    def logw(mask, u):
        return min(0.0, (energy(mask | (1 << u)) - energy(mask)) / 0.5)

    res = T.exact(lat, logw, focal_unit=focal, permissive_units=perm)
    paths = _brute_force(lat, logw)
    logZ = math.log(sum(math.exp(lp) for lp in paths.values()))
    for k in range(len(perm) + 1):
        expect = 0.0
        for order, lp in paths.items():
            pos = {u: i for i, u in enumerate(order)}
            present = sum(1 for p in perm if pos[p] < pos[focal])
            if present >= k:
                expect += math.exp(lp - logZ)
        assert res.designated[f"P(focal after >= {k} of permissive)"] == pytest.approx(expect, abs=1e-10)
