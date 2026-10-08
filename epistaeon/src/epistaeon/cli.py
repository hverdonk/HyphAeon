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
from .lattice import Lattice
from .timing import TimingNotSupported, absolute_timing, relative_positions
from .trajectory import EXACT_UNIT_LIMIT, infer


def _build_scorer(name, ctx, lattice, genotypes, anchor_name, args, provenance):
    """Assemble the requested scorer from model-derived quantities."""
    from . import scoring as S
    from .background import extant_embeddings, root_token_embedding

    embedding_of = genotypes.embedding
    if anchor_name == "root-token":
        anchor = root_token_embedding(ctx)
    else:
        anchor = embedding_of(0)
    provenance["phi_anchor"] = anchor_name

    alpha = 1.0
    if args.alpha is not None:
        alpha = float(args.alpha)
    provenance["alpha"] = alpha
    phi = S.phi_divergence(embedding_of, anchor, alpha)
    dphi = S.delta_phi_fn(phi)

    if name == "potts-metropolis":
        from .couplings import couplings_from_coselection, site_terms_from_lrt
        site_terms = site_terms_from_lrt(genotypes, lattice)
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
        cloud = extant_embeddings(ctx)
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
    from .background import GenotypeRows, derive_endpoints, load_context, resolve_positions

    if args.anchor == "ancestral-sequence" and args.ancestor is None:
        raise SystemExit(
            "--anchor ancestral-sequence needs --ancestor; with no ancestor row "
            "the only available anchor is --anchor root-token"
        )
    anchor_name = args.anchor or ("ancestral-sequence" if args.ancestor else "root-token")

    mode = args.mode
    ctx = load_context(
        args.alignment, args.tree, derived_taxon=args.focal_taxon,
        ancestor_taxon=args.ancestor,
        weights=args.weights, variant=args.variant, cpu=args.cpu,
        max_species=args.max_species, use_tn93=args.use_tn93,
    )

    # Positions are numbered within the reference row, so --sites and --blocks
    # can only be resolved once the endpoints exist; derive them unrestricted
    # first, then again with the requested columns.
    rows = ctx.rows()
    ends = derive_endpoints(rows, ctx.derived_name, ctx.ancestor_name)
    if args.sites:
        columns = resolve_positions(open(args.sites).read().split(), ends)
        ends = derive_endpoints(rows, ctx.derived_name, ctx.ancestor_name, columns=columns)
    subs = ends.substitutions
    if not subs:
        raise SystemExit(
            "the ancestral and derived rows are identical at the sites considered"
        )

    blocks = names = None
    if args.blocks:
        # resolved column indices are not what the user typed, so the groups are
        # checked here rather than in blocks_from_groups, which cannot name the
        # offending position in their coordinates
        groups = json.load(open(args.blocks))
        by_column = {s.index: i for i, s in enumerate(subs)}
        blocks, names, claimed = [], [], set()
        for group, tokens in groups.items():
            members = []
            for token, col in zip(tokens, resolve_positions(tokens, ends)):
                if col not in by_column:
                    raise SystemExit(
                        f"--blocks: group {group!r} names position {token}, which is "
                        f"not one of the {len(subs)} substitutions between these endpoints"
                    )
                if by_column[col] in claimed:
                    raise SystemExit(
                        f"--blocks: position {token} appears in more than one group"
                    )
                members.append(by_column[col])
                claimed.add(by_column[col])
            if members:
                blocks.append(sorted(members))
                names.append(group)
        for i, sub in enumerate(subs):
            if i not in claimed:
                blocks.append([i])
                names.append(sub.name())
    lattice = Lattice(ends.ancestral_protein, subs, blocks=blocks, block_names=names)
    source = ctx.ancestor_name or "column plurality"
    counts = ends.class_counts()
    breakdown = ", ".join(f"{n} {k}" for k, n in sorted(counts.items()))
    print(f"[*] {lattice.n_substitutions} codon substitutions from {source} to "
          f"{ctx.derived_name} ({breakdown}), {lattice.n_units} ordering units")
    provenance: Dict[str, object] = {
        "epistaeon_version": __version__,
        "model_variant": args.variant or "(explicit weights)",
        "weights": args.weights,
        "alignment": args.alignment,
        "tree": args.tree,
        "distance_source": ctx.distance_source,
        "derived_taxon": ctx.derived_name,
        "ancestor_taxon": ctx.ancestor_name,
        "ancestral_states": (
            f"alignment row {ctx.ancestor_name!r}" if ctx.ancestor_name
            else "plurality residue per column over the other rows (inferred, "
                 "not a reconstruction)"
        ),
        "position_numbering": f"1-based residues of {ends.reference_name!r}",
        "substitution_unit": "codon change (synonymous changes included)",
        "substitution_classes": ends.class_counts(),
        "endpoints_read_from_the_alignment": {
            name: name in ctx.taxa
            for name in (ctx.ancestor_name, ctx.derived_name) if name
        },
        "n_taxa": ctx.n_taxa,
        "n_codons": ctx.n_codons,
        "scorer": args.scorer,
        "mode": mode,
        "t_sel": args.t_sel,
    }
    if ends.plurality_support:
        provenance["plurality_support"] = {
            lattice.unit_name(u): min(
                ends.plurality_support[s.index] for s in lattice.unit_members(u)
            )
            for u in range(lattice.n_units)
        }

    # Each genotype carries its own distance row, so its embedding costs a full
    # forward pass and nothing is shared between genotypes. Exact inference
    # touches all 2^n of them, which is the binding cost, not the DP itself.
    n_genotypes = 1 << lattice.n_units
    if mode == "auto":
        mode = ("exact" if lattice.n_units <= EXACT_UNIT_LIMIT
                and n_genotypes <= args.max_genotypes else "sample")
    if mode == "exact":
        if lattice.n_units > EXACT_UNIT_LIMIT:
            raise SystemExit(
                f"exact mode needs <= {EXACT_UNIT_LIMIT} units, got {lattice.n_units}. "
                "Group them with --blocks, restrict with --sites, or use --mode sample."
            )
        if n_genotypes > args.max_genotypes:
            raise SystemExit(
                f"exact mode over {lattice.n_units} units needs {n_genotypes} forward "
                f"passes, above --max-genotypes {args.max_genotypes}. Group units with "
                "--blocks, use --mode sample, or raise the budget."
            )
    print(f"[*] inference mode: {mode}")
    provenance["mode"] = mode

    genotypes = GenotypeRows(ctx, ends, lattice, budget=args.max_genotypes)
    provenance["genotype_placement"] = (
        "each intermediate is appended to the alignment as its own row, placed "
        "by TN93 distance to every taxon, with MDS coordinates projected into "
        "the existing frame; the endpoints are read from their own rows rather "
        "than appended as copies of themselves"
    )
    if mode == "exact":
        # one TN93 call for the whole lattice beats one per genotype
        print(f"[*] placing {n_genotypes} genotypes by TN93")
        genotypes.prefetch(range(n_genotypes))
    else:
        reachable = min(n_genotypes, args.n_samples * lattice.n_units)
        print(f"[*] sampling will embed up to {reachable} genotypes, one forward "
              f"pass each, against a budget of {args.max_genotypes}")

    # How large is a single-substitution signal? Both genotypes compared here
    # must be appended rows: measuring a unit against the ancestral endpoint
    # instead measures the extra row, which is the larger effect by far.
    z_scale = max(1e-12, float(np.linalg.norm(genotypes.embedding(0))))
    probed = min(lattice.n_units, 6) if lattice.n_units > 1 else 0
    deltas = []
    for u in range(probed):
        base = 1 << ((u + 1) % lattice.n_units)
        deltas.append(float(
            np.linalg.norm(genotypes.embedding(base | 1 << u) - genotypes.embedding(base))
            / z_scale
        ))
    if deltas:
        provenance["per_unit_relative_delta"] = {
            "min": min(deltas), "median": float(np.median(deltas)), "max": max(deltas),
            "units_probed": probed,
            "float32_eps": float(np.finfo(np.float32).eps),
            "_meaning": "||dz|| / ||z|| for one unit applied to a background "
                        "that already carries another unit, so both genotypes "
                        "are appended rows: how much of the embedding a single "
                        "substitution moves, against float32 precision",
        }
    endpoint_steps = [
        float(np.linalg.norm(genotypes.embedding(1 << u) - genotypes.embedding(0)) / z_scale)
        for u in range(probed or lattice.n_units)
    ]
    provenance["first_step_relative_delta"] = {
        "median": float(np.median(endpoint_steps)),
        "_meaning": "the same measure for the first step out of the ancestral "
                    "endpoint, which is read from its own row while the "
                    "intermediate is appended. The gap against "
                    "per_unit_relative_delta is the frame shift, not biology.",
    }

    scorer, phi = _build_scorer(
        args.scorer, ctx, lattice, genotypes, anchor_name, args, provenance
    )

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
    provenance["forward_passes"] = genotypes.forward_passes
    provenance["tn93_calls"] = genotypes.tn93_calls
    provenance["genotypes_read_from_a_real_row"] = genotypes.real_row_hits
    provenance["endpoint_frame_shift"] = {
        "value": genotypes.endpoint_frame_shift(),
        "_meaning": "an endpoint is embedded among N taxa and an intermediate "
                    "among N + 1; this is how far that extra row moves a "
                    "taxon's embedding, relative to its norm. Compare it with "
                    "per_unit_relative_delta before trusting a step that "
                    "crosses between the two.",
    }

    payload: Dict[str, object] = {
        "provenance": provenance,
        "substitutions": [
            {
                "unit": u,
                "name": lattice.unit_name(u),
                "members": [
                    {"column0": s.index,
                     "residue_number": ends.residue_number.get(s.index),
                     "ancestral": s.ancestral,
                     "derived": s.derived,
                     "ancestral_codon": ends.ancestral_codons[s.index],
                     "derived_codon": ends.derived_codons[s.index],
                     "class": ends.classes.get(s.index),
                     "is_indel": s.is_indel}
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
            "this pipeline conditions a site on the rest of its own sequence. "
            "What differs between two genotypes' forward passes is the "
            "genotype's own distances to the other taxa, not cross-site "
            "context within the genotype.",
        ],
    }
    if not ctx.ancestor_name:
        payload["caveats"].append(
            "No --ancestor row was given, so each ancestral state is the "
            "plurality residue in that alignment column. That is a consensus, "
            "not an ancestral reconstruction, and it can differ from the real "
            "ancestor at any column; see provenance.plurality_support."
        )

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
    o.add_argument("--focal-taxon", "--descendant", "--derived", required=True,
                   dest="focal_taxon", metavar="NAME",
                   help="alignment row holding the derived (end-point) sequence")
    o.add_argument("--ancestor", default=None, metavar="NAME",
                   help="alignment row holding the ancestral sequence; omit to take "
                        "each ancestral state from its column's plurality residue "
                        "and anchor phi at the [ROOT] embedding")
    o.add_argument("--scorer", default="potts-metropolis",
                   choices=["potts-metropolis", "typicality-gate", "softmax-repaired", "dphi-product"])
    o.add_argument("--mode", default="auto", choices=["auto", "exact", "sample"])
    o.add_argument("--blocks", default=None,
                   help="JSON mapping group name -> positions, numbered as --sites")
    o.add_argument("--sites", default=None,
                   help="file of positions to restrict to: 1-based residues of the "
                        "ancestor row (or of the derived row when there is none), "
                        "or c<N> for 1-based alignment column N")
    o.add_argument("--anchor", default=None,
                   choices=["ancestral-sequence", "root-token"],
                   help="phi anchor; defaults to the ancestral sequence when "
                        "--ancestor is given and to root-token otherwise")
    o.add_argument("--alpha", type=float, default=None, help="isometric scale for phi")
    o.add_argument("--t-sel", type=float, default=0.5, dest="t_sel")
    o.add_argument("--viability-floor", type=float, default=None)
    o.add_argument("--focal-unit", default=None, help="unit index for designated marginals")
    o.add_argument("--permissive-units", default=None, help="comma-separated unit indices")
    o.add_argument("--max-genotypes", type=int, default=4096,
                   help="forward-pass budget for exact mode, which embeds all 2^n "
                        "genotypes (default 4096)")
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
