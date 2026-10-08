"""
epistaeon/couplings.py
----------------------
Potts terms taken from the model, for the default scorer.

Site terms are the model's own output: the change in predicted selection
evidence (LRT) at a site when the genotype row's residue is switched from the
ancestral to the derived state. Couplings come from hyphaeon's attribution
co-selection network (CESI), which scores how much two sites' attribution
vectors move together across taxa.

The honest caveat, which belongs in any write-up of results: CESI is a
correlation statistic over attribution patterns, not a fitted coupling
constant. It supplies the background dependence that the site-independent
encoder cannot, but it is not learned site-site coupling.
"""

from typing import Dict, Tuple

from .background import Context, GenotypeRows
from .lattice import Lattice


def site_terms_from_lrt(genotypes: GenotypeRows, lattice: Lattice) -> Dict[int, float]:
    """f_u = LRT(derived state) - LRT(ancestral state) at the unit's sites.

    Both states are scored as their own genotype row, the same way the
    embeddings are, so the site term sees the genotype placed by its own
    distances rather than borrowing a real taxon's position.
    """
    out: Dict[int, float] = {}
    for u in range(lattice.n_units):
        sites = sorted({s.index for s in lattice.unit_members(u)})
        out[u] = genotypes.lrt(1 << u, sites) - genotypes.lrt(0, sites)
    return out


def couplings_from_coselection(
    ctx: Context,
    lattice: Lattice,
    max_fdr: float = 0.05,
    min_sim: float = 0.0,
) -> Tuple[Dict[Tuple[int, int], float], Dict[str, float]]:
    """J_uv from the attribution co-selection network, restricted to our units.

    Returns the couplings plus diagnostics, including how many of the unit
    pairs the network actually scored -- if that is zero there is no epistatic
    signal available and the caller should say so rather than proceed silently.
    """
    from hyphaeon.epistasis import (
        compute_branch_coselection_network,
        compute_transformer_attributions,
    )

    attributions, lrts, _pvals, consensus = compute_transformer_attributions(
        ctx.model, ctx.c_tensor, ctx.a_tensor, ctx.tree_cache, ctx.taxa,
        ctx.device, progress=False,
    )
    pairs, _graph = compute_branch_coselection_network(
        attributions, lrts, ctx.taxa, consensus,
        min_sim=min_sim, max_fdr=max_fdr,
    )
    site_to_unit = lattice.unit_index_by_site()
    cesi: Dict[Tuple[int, int], float] = {}
    for p in pairs:
        # network reports 1-based site numbers
        u = site_to_unit.get(int(p["site_u"]) - 1)
        v = site_to_unit.get(int(p["site_v"]) - 1)
        if u is None or v is None or u == v:
            continue
        key = (min(u, v), max(u, v))
        cesi[key] = max(cesi.get(key, 0.0), float(p["cesi"]))
    n_possible = lattice.n_units * (lattice.n_units - 1) // 2
    diag = {
        "coselection_pairs_total": float(len(pairs)),
        "unit_pairs_with_coupling": float(len(cesi)),
        "unit_pairs_possible": float(n_possible),
    }
    return cesi, diag
