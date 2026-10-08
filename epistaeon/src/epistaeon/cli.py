"""
epistaeon/cli.py
----------------
`epistaeon order` -- rank mutational orderings between two sequences.
"""

import argparse
import json
import sys
from typing import Dict, List, Optional

import numpy as np

from . import __version__
from .io import write_json, write_matrix_csv
from .lattice import Lattice, blocks_from_groups, derive_substitutions
from .timing import TimingNotSupported, absolute_timing, relative_positions
from .trajectory import EXACT_UNIT_LIMIT, infer


def _read_fasta_first(path: str) -> str:
    seq, started = [], False
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if started:
                break
            started = True
            continue
        seq.append(line)
    return "".join(seq)


def _build_scorer(name, ctx, lattice, cache, args, provenance):
    """Assemble the requested scorer from model-derived quantities."""
    from . import scoring as S
    from .background import extant_embeddings, isometric_scale, root_token_embedding

    def embedding_of(mask: int) -> np.ndarray:
        return cache.embedding(lattice, mask)

    cloud = extant_embeddings(ctx)
    anchor_name = "ancestral-sequence"
    anchor = embedding_of(0)
    if args.anchor == "root-token":
        anchor = root_token_embedding(ctx)
        anchor_name = "root-token"
    provenance["phi_anchor"] = anchor_name

    alpha = 1.0
    if args.alpha is not None:
        alpha = float(args.alpha)
    provenance["alpha"] = alpha
    phi = S.phi_divergence(embedding_of, anchor, alpha)
    dphi = S.delta_phi_fn(phi)

    if name == "potts-metropolis":
        from .couplings import couplings_from_coselection, site_terms_from_lrt
        site_terms = site_terms_from_lrt(ctx, lattice)
        cesi, diag = couplings_from_coselection(ctx, lattice)
        provenance.update(diag)
        if not cesi:
            print(
                "[!] the co-selection network scored no coupling between these "
                "units, so this scorer has no background dependence to use. "
                "Results will rank orderings by site terms alone.",
                file=sys.stderr,
            )
        energy = S.pairwise_energy(lattice, site_terms, cesi)
        provenance["site_terms"] = {lattice.unit_name(u): v for u, v in site_terms.items()}
        return S.PottsMetropolis(energy=energy, temperature=args.t_sel), phi

    if name == "typicality-gate":
        viability = S.typicality_mahalanobis(embedding_of, cloud)
        floor = None if args.viability_floor is None else float(args.viability_floor)
        return S.TypicalityGate(viability=viability, temperature=args.t_sel, floor=floor), phi

    if name == "softmax-repaired":
        def alternatives(mask: int, unit: int) -> np.ndarray:
            # compete the historical move against the other units still available
            return np.array([dphi(mask, u) for u in lattice.available(mask) if u != unit])
        return S.SoftmaxRepaired(delta_phi=dphi, alternatives=alternatives,
                                 temperature=args.t_sel), phi

    if name == "dphi-product":
        return S.DeltaPhiProduct(delta_phi=dphi, temperature=args.t_sel), phi

    raise SystemExit(f"unknown scorer {name!r}; choose from {sorted(S.SCORERS)}")


def ctx_focal_protein(ctx) -> str:
    """Focal row as a full-length string over alignment columns, gaps included."""
    from aeon_core.dataset import AA_MAP
    rev = {v: k for k, v in AA_MAP.items()}
    a = ctx.a_tensor.squeeze(-1).numpy()
    return "".join(rev.get(int(a[s, ctx.focal_index]), "-") for s in range(ctx.n_codons))


def _relattice(aligned_ancestral, mapped, args, original_subs, colmap):
    """Rebuild the lattice in alignment coordinates, re-applying any blocking."""
    anc = list(aligned_ancestral)
    for sub in mapped:
        anc[sub.index] = sub.ancestral          # the branch starts from the ancestral state
    blocks = names = None
    if args.blocks:
        groups = json.load(open(args.blocks))
        by_old = {s.index: i for i, s in enumerate(original_subs)}
        remap = {}
        for name, positions in groups.items():
            cols = []
            for pos in positions:
                old_index = int(pos) - 1
                if old_index in colmap:
                    cols.append(colmap[old_index])
            if cols:
                remap[name] = cols
        index_of = {s.index: i for i, s in enumerate(mapped)}
        blocks, names = [], []
        claimed = set()
        for name, cols in remap.items():
            members = [index_of[c] for c in cols if c in index_of]
            if members:
                blocks.append(sorted(members)); names.append(name); claimed.update(members)
        for i, s in enumerate(mapped):
            if i not in claimed:
                blocks.append([i]); names.append(s.name())
    return Lattice("".join(anc), list(mapped), blocks=blocks, block_names=names)


def _discrimination(res) -> Dict[str, object]:
    """Did the scorer actually separate the orderings?

    A normalised P(MAP) at the uniform baseline, with every before-probability
    at 0.5, means the edge weights telescoped: the DP ran correctly and found
    that all orderings are equally weighted. That is a property of the scorer,
    not evidence that the orderings are equally viable, so it is reported
    rather than left for the reader to infer.
    """
    import numpy as _np
    ratio = (res.map_prob_normalised / res.uniform_baseline) if res.uniform_baseline else float("nan")
    off = ~_np.eye(res.n_units, dtype=bool)
    flat = bool(res.n_units > 1 and _np.allclose(res.before[off], 0.5, atol=1e-6))
    degenerate = bool(flat and abs(ratio - 1.0) < 1e-6)
    out = {
        "map_prob_over_uniform": float(ratio),
        "before_matrix_is_flat": flat,
        "degenerate_uniform": degenerate,
    }
    if degenerate:
        out["explanation"] = (
            "every ordering carries identical weight, so this scorer provides no "
            "background dependence on these units. Sum of a state-function "
            "difference telescopes; check whether the coupling terms are empty."
        )
    return out


def cmd_order(args: argparse.Namespace) -> int:
    from .background import build_additive_cache, load_context, verify_additivity

    ancestral = _read_fasta_first(args.ancestor)
    derived = _read_fasta_first(args.descendant)
    sites = None
    if args.sites:
        sites = [int(x) - 1 for x in open(args.sites).read().split()]
    subs = derive_substitutions(ancestral, derived, sites=sites)
    if not subs:
        raise SystemExit("ancestor and descendant are identical at the sites considered")

    mode = args.mode
    ctx = load_context(
        args.alignment, args.tree, focal_taxon=args.focal_taxon,
        weights=args.weights, variant=args.variant, cpu=args.cpu,
        max_species=args.max_species, use_tn93=args.use_tn93,
    )

    # The ancestor/derived pair is usually a domain while the alignment is full
    # length, so substitution indices must be translated into alignment columns
    # before residues can be written into the focal row.
    from .background import map_to_alignment
    colmap = map_to_alignment(ctx, ancestral)
    unmapped = [s for s in subs if s.index not in colmap]
    if unmapped:
        print(
            f"[!] {len(unmapped)} of {len(subs)} substitutions could not be "
            f"placed in the alignment and are dropped: "
            f"{[s.name() for s in unmapped][:8]}",
            file=sys.stderr,
        )
    mapped = [
        type(s)(index=colmap[s.index], ancestral=s.ancestral,
                derived=s.derived, label=s.name())
        for s in subs if s.index in colmap
    ]
    if not mapped:
        raise SystemExit(
            "no substitution could be placed in the alignment; check that "
            "--ancestor matches the focal taxon's protein"
        )
    aligned_ancestral = list(ctx_focal_protein(ctx))
    lattice = _relattice(aligned_ancestral, mapped, args, subs, colmap)
    print(f"[*] placed {lattice.n_substitutions} substitutions, "
          f"{lattice.n_units} ordering units")
    provenance: Dict[str, object] = {
        "epistaeon_version": __version__,
        "model_variant": args.variant or "(explicit weights)",
        "weights": args.weights,
        "alignment": args.alignment,
        "tree": args.tree,
        "focal_taxon": ctx.focal_name,
        "n_taxa": ctx.n_taxa,
        "n_codons": ctx.n_codons,
        "scorer": args.scorer,
        "mode": mode,
        "t_sel": args.t_sel,
    }

    if mode == "auto":
        mode = "exact" if lattice.n_units <= EXACT_UNIT_LIMIT else "sample"
    if mode == "exact" and lattice.n_units > EXACT_UNIT_LIMIT:
        raise SystemExit(
            f"exact mode needs <= {EXACT_UNIT_LIMIT} units, got {lattice.n_units}. "
            "Group them with --blocks, restrict with --sites, or use --mode sample."
        )
    print(f"[*] inference mode: {mode}")
    provenance["mode"] = mode

    cache = build_additive_cache(ctx, lattice)
    ok, err = verify_additivity(ctx, lattice, cache)
    provenance["forward_passes"] = cache.passes_used
    provenance["additivity_ok"] = bool(ok)
    provenance["additivity_relative_error"] = err
    if not ok:
        raise SystemExit(
            f"additivity check failed (relative error {err:.3g}). The cached "
            "embeddings are invalid for this checkpoint -- it appears to mix "
            "information across sites. Re-run without the fast path."
        )
    # how large is a single-substitution signal, against float precision?
    z_anc = cache.embedding(lattice, 0)
    z_scale = max(1e-12, float(np.linalg.norm(z_anc)))
    deltas = [
        float(np.linalg.norm((d - a) / float(ctx.n_codons)) / z_scale)
        for a, d in cache.per_unit.values()
    ]
    provenance["per_unit_relative_delta"] = {
        "min": min(deltas), "median": float(np.median(deltas)), "max": max(deltas),
        "float32_eps": float(np.finfo(np.float32).eps),
        "_meaning": "||dz|| / ||z|| for one unit: how much of the embedding a "
                    "single substitution moves, against float32 precision",
    }

    scorer, phi = _build_scorer(args.scorer, ctx, lattice, cache, args, provenance)

    focal_unit = None
    permissive = None
    if args.focal_unit is not None:
        focal_unit = int(args.focal_unit)
    if args.permissive_units:
        permissive = [int(x) for x in args.permissive_units.split(",")]

    res = infer(
        lattice, scorer.log_weight, mode=mode,
        focal_unit=focal_unit, permissive_units=permissive,
        **({"n_samples": args.n_samples, "seed": args.seed} if mode == "sample" else {}),
    )
    print(f"[*] {res.summary()}")

    payload: Dict[str, object] = {
        "provenance": provenance,
        "substitutions": [
            {
                "unit": u,
                "name": lattice.unit_name(u),
                "members": [
                    {"index0": s.index, "ancestral": s.ancestral,
                     "derived": s.derived, "is_indel": s.is_indel}
                    for s in lattice.unit_members(u)
                ],
            }
            for u in range(lattice.n_units)
        ],
        "accessible": res.accessible,
        "note": res.note,
        "map_order": res.map_order_names,
        "map_log_prob": res.map_log_prob,
        "log_total_weight": res.log_total_weight,
        "map_prob_normalised": res.map_prob_normalised,
        "uniform_baseline": res.uniform_baseline,
        "before_matrix": {"names": res.unit_names, "matrix": res.before},
        "position_matrix": {"names": res.unit_names, "matrix": res.position},
        "designated_marginals": res.designated,
        "diagnostics": res.diagnostics,
        "discrimination": _discrimination(res),
        "caveats": [
            "P(pi | endpoint) is conditional on reaching the derived sequence by "
            "exactly these substitutions; it ranks orderings and is not an "
            "absolute accessibility.",
            "A normalised P(MAP) near the uniform baseline means the scorer is "
            "not discriminating, not that all orderings are equally viable.",
            "The checkpoint scores each codon site independently, so no part of "
            "this pipeline conditions a site on the rest of its own sequence.",
        ],
    }

    if res.accessible and res.map_order:
        phis = [phi(0)]
        mask = 0
        for u in res.map_order:
            mask |= 1 << u
            phis.append(phi(mask))
        pos = relative_positions(res.unit_names, res.map_order, phis)
        payload["relative_branch_position"] = pos.as_dict()
        payload["relative_branch_position_monotone"] = pos.monotone
        if pos.note:
            payload["relative_branch_position_note"] = pos.note
        if args.ne is not None or args.s is not None:
            try:
                payload["absolute_timing"] = absolute_timing(
                    pos, args.ne, args.s, args.generation_time
                )
            except TimingNotSupported as exc:
                payload["absolute_timing_error"] = str(exc)
                print(f"[!] absolute timing not available: {exc}", file=sys.stderr)

    write_json(f"{args.output}/order.json", payload)
    write_matrix_csv(f"{args.output}/before_matrix.csv", res.before, res.unit_names)
    write_matrix_csv(f"{args.output}/position_matrix.csv", res.position, res.unit_names)
    print(f"[✓] wrote {args.output}/order.json and the two matrices")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="epistaeon",
        description="Score mutational trajectories between an ancestral and a derived sequence",
    )
    p.add_argument("--version", action="version", version=f"epistaeon {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    o = sub.add_parser("order", help="rank orderings of the substitutions between two sequences")
    o.add_argument("--alignment", required=True)
    o.add_argument("--tree", default=None, help="Newick tree; omit with --use-tn93")
    o.add_argument("--ancestor", required=True, help="FASTA with the ancestral sequence")
    o.add_argument("--descendant", required=True, help="FASTA with the derived sequence")
    o.add_argument("--focal-taxon", default=None, help="alignment row genotypes are written into")
    o.add_argument("--scorer", default="potts-metropolis",
                   choices=["potts-metropolis", "typicality-gate", "softmax-repaired", "dphi-product"])
    o.add_argument("--mode", default="auto", choices=["auto", "exact", "sample"])
    o.add_argument("--blocks", default=None, help="JSON mapping group name -> 1-based positions")
    o.add_argument("--sites", default=None, help="file of 1-based positions to restrict to")
    o.add_argument("--anchor", default="ancestral-sequence",
                   choices=["ancestral-sequence", "root-token"])
    o.add_argument("--alpha", type=float, default=None, help="isometric scale for phi")
    o.add_argument("--t-sel", type=float, default=0.5, dest="t_sel")
    o.add_argument("--viability-floor", type=float, default=None)
    o.add_argument("--focal-unit", default=None, help="unit index for designated marginals")
    o.add_argument("--permissive-units", default=None, help="comma-separated unit indices")
    o.add_argument("--n-samples", type=int, default=20000)
    o.add_argument("--seed", type=int, default=0)
    o.add_argument("--ne", type=float, default=None, help="effective population size")
    o.add_argument("--s", type=float, default=None, help="selection coefficient")
    o.add_argument("--generation-time", type=float, default=None, help="years per generation")
    o.add_argument("--weights", default=None)
    o.add_argument("--variant", default=None)
    o.add_argument("--cpu", action="store_true")
    o.add_argument("--max-species", type=int, default=None)
    o.add_argument("--use-tn93", action="store_true")
    o.add_argument("-o", "--output", default="epistaeon_out")
    o.set_defaults(func=cmd_order)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
