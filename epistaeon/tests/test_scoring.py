"""Each scorer must break path-independence, except the documented baseline."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon.lattice import Lattice, derive_substitutions  # noqa: E402
from epistaeon import scoring as S, trajectory as T  # noqa: E402


def _lat(n):
    anc, der = "A" * n, "C" * n
    return Lattice(anc, derive_substitutions(anc, der))


def _coupled_energy(lat, seed=5):
    rng = np.random.default_rng(seed)
    n = lat.n_units
    site = {u: float(v) for u, v in enumerate(rng.normal(size=n))}
    pairs = {}
    for i in range(n):
        for j in range(i + 1, n):
            pairs[(i, j)] = float(rng.normal() * 1.2)
    return S.pairwise_energy(lat, site, pairs)


def test_pairwise_energy_is_background_dependent():
    """dE for the same unit must differ by background, or there is no epistasis."""
    lat = _lat(3)
    energy = _coupled_energy(lat)
    d_alone = energy(0b001) - energy(0b000)
    d_with_other = energy(0b011) - energy(0b010)
    assert not math.isclose(d_alone, d_with_other, rel_tol=1e-6)


def test_potts_metropolis_discriminates_orderings():
    lat = _lat(4)
    sc = S.PottsMetropolis(energy=_coupled_energy(lat), temperature=0.5)
    res = T.exact(lat, sc.log_weight)
    off = ~np.eye(lat.n_units, dtype=bool)
    assert res.map_prob_normalised > res.uniform_baseline
    assert not np.allclose(res.before[off], 0.5, atol=0.05)


def test_typicality_gate_blocks_below_floor():
    lat = _lat(3)
    # unit 2 produces a state far below the floor whenever it is applied
    def viability(mask):
        return -10.0 if mask >> 2 & 1 else 0.0

    sc = S.TypicalityGate(viability=viability, floor=-1.0)
    assert sc.log_weight(0b000, 2) == -math.inf
    assert math.isfinite(sc.log_weight(0b000, 0))
    res = T.exact(lat, sc.log_weight)
    assert not res.accessible          # every ordering must apply unit 2


def test_dphi_product_temperature_is_inert():
    """Documented pathology: T cannot change the ranking."""
    lat = _lat(4)
    energy = _coupled_energy(lat)

    def dphi(mask, u):
        return abs(energy(mask | (1 << u)) - energy(mask)) + 0.1   # keep positive

    orders = []
    for temp in (0.25, 1.0, 4.0):
        sc = S.DeltaPhiProduct(delta_phi=dphi, temperature=temp)
        orders.append(tuple(T.exact(lat, sc.log_weight).map_order))
    assert len(set(orders)) == 1


def test_softmax_repaired_last_step_is_informative():
    """With alternatives in the pool, the final step is not probability 1."""
    lat = _lat(2)
    energy = _coupled_energy(lat)

    def dphi(mask, u):
        return energy(mask | (1 << u)) - energy(mask)

    def alts(mask, u):            # three competing non-historical moves
        return np.array([dphi(mask, u) - 0.5, dphi(mask, u) - 1.0, dphi(mask, u) + 0.25])

    sc = S.SoftmaxRepaired(delta_phi=dphi, alternatives=alts)
    last = sc.log_weight(0b01, 1)          # only one historical move remains
    assert last < -1e-6, "final step should still cost something"
    assert math.exp(last) < 0.99


def test_build_scorer_selects_appended_genotypes_for_softmax_repaired_only():
    import torch
    from argparse import Namespace
    from epistaeon.cli import _build_scorer
    from epistaeon.background import AppendedGenotypeRows, GenotypeRows, derive_endpoints, codons_of

    lat = _lat(2)
    rows = {"anc": codons_of("AAAAAA", 2), "der": codons_of("CCCCCC", 2)}
    ends = derive_endpoints(rows, "der", "anc")

    ctx = Namespace(
        model=None, taxa=["anc", "der"], n_codons=2, derived_name="der", ancestor_name="anc",
        c_tensor=torch.zeros((2, 2, 1), dtype=torch.long), a_tensor=torch.zeros((2, 2, 1), dtype=torch.long),
        dist=np.zeros((2, 2)), mds=np.zeros((2, 4)), tree_cache={}, seqs={"anc": "AAAAAA", "der": "CCCCCC"},
        derived_index=1, _extant=np.zeros((2, 384)),
    )
    args = Namespace(t_sel=0.5, alpha=1.0, viability_floor=None)

    g_normal = GenotypeRows.__new__(GenotypeRows)
    g_normal.ends = ends
    g_normal.budget = 10
    g_normal.embedding = lambda mask: np.array([0.1, 0.2])
    prov = {}

    sc_sr, _ = _build_scorer("softmax-repaired", ctx, lat, g_normal, "ancestral-sequence", args, prov)
    assert sc_sr.name == "softmax-repaired"

    sc_dp, _ = _build_scorer("dphi-product", ctx, lat, g_normal, "ancestral-sequence", args, prov)
    assert sc_dp.name == "dphi-product"
