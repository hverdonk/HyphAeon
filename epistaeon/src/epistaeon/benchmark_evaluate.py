"""Evaluate sealed predictions against historical truth without loading inference.

This module deliberately depends only on the standard library and the artifact
integrity helpers. No alignment selection, checkpoint loading or rerun is possible
through this interface. Historical branch precedence is distinct from biochemical
accessibility, which is not scored here.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path
from typing import Any

from .benchmark_artifacts import SCORERS, sha256_file, verify_seal, write_json


_PROBABILITY_TOLERANCE = 1e-6
_MAX_CANDIDATE_EVENTS = 8
_RESOLVED_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON object key: {key}")
        value[key] = item
    return value


def _reject_json_constant(constant: str) -> None:
    raise ValueError(f"Nonstandard JSON numeric constant: {constant}")


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle, object_pairs_hook=_json_object,
                          parse_constant=_reject_json_constant)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path.name}")
    return value


def _inside(path: Path, directory: Path) -> bool:
    return path == directory or directory in path.parents


def _objects(value: Any, label: str) -> list[dict[str, Any]]:
    if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
        raise ValueError(f"{label} must be a list of objects")
    return value


def _identifiers(value: Any, label: str) -> list[str]:
    if (not isinstance(value, list) or not value
            or any(not isinstance(item, str) or not item for item in value)
            or len(set(value)) != len(value)):
        raise ValueError(f"{label} must be a nonempty list of unique identifiers")
    return value


def _direction(ancestral: Any, derived: Any, label: str) -> None:
    if (not isinstance(ancestral, str) or ancestral not in _RESOLVED_AMINO_ACIDS
            or not isinstance(derived, str) or derived not in _RESOLVED_AMINO_ACIDS
            or ancestral == derived):
        raise ValueError(f"{label} must be a resolved nonsynonymous amino-acid substitution")


def _reference_order(reference: dict[str, Any]) -> dict[str, Any]:
    """Expand historical groups only; do not convert assay rules into history."""
    if reference.get("schema_version", "1.0") != "1.0":
        raise ValueError("Unsupported substitution-reference schema version")
    _objects(reference.get("sources", []), "Reference sources")
    if not isinstance(reference.get("assay_data", {}), dict):
        raise ValueError("Reference assay_data must be an object")
    reference_events = _objects(reference.get("events"), "Reference events")
    if any(not isinstance(event.get("order_scoring", True), bool) for event in reference_events):
        raise ValueError("Reference order_scoring flags must be boolean")
    targets = [event for event in reference_events if event.get("order_scoring", True)]
    if not targets:
        raise ValueError("Reference contains no order-scoring events")
    events: dict[str, dict[str, Any]] = {}
    keys: dict[tuple[int, str, str], str] = {}
    for event in targets:
        identifier = event.get("id")
        if not isinstance(identifier, str) or not identifier:
            raise ValueError("Reference events require nonempty identifiers")
        if identifier in events:
            raise ValueError(f"Duplicate reference event: {identifier}")
        columns = set()
        for mapping in _objects(event.get("structure_mapping"), "Reference structure mappings"):
            crosswalk = mapping.get("alignment_reference_columns_1based", {})
            if not isinstance(crosswalk, dict):
                raise ValueError(f"Invalid alignment crosswalk for {identifier}")
            for column in crosswalk.values():
                if column is not None:
                    if type(column) is not int or column < 1:
                        raise ValueError(f"Invalid alignment column for {identifier}")
                    columns.add(column)
        if len(columns) != 1:
            raise ValueError(f"Reference event requires one unambiguous alignment column: {identifier}")
        column = next(iter(columns))
        if isinstance(column, bool) or not isinstance(column, int) or column < 1:
            raise ValueError(f"Invalid alignment column for {identifier}")
        ancestral, derived = event.get("from"), event.get("to")
        _direction(ancestral, derived, f"Reference event {identifier}")
        key = (column - 1, ancestral, derived)
        if key in keys:
            raise ValueError(f"Reference events share a substitution key: {identifier}, {keys[key]}")
        keys[key] = identifier
        events[identifier] = {"id": identifier, "label": event.get("label", identifier),
                              "column0": column - 1, "ancestral": ancestral,
                              "derived": derived}

    groups: dict[str, list[str]] = {}
    for group in _objects(reference.get("groups", []), "Reference groups"):
        identifier = group.get("id")
        if not isinstance(identifier, str) or not identifier:
            raise ValueError("Reference groups require nonempty identifiers")
        if identifier in events or identifier in groups:
            raise ValueError(f"Duplicate reference group identifier: {identifier}")
        groups[identifier] = _identifiers(group.get("members"), "Reference group members")

    def expand(identifier: str, visited: frozenset[str] = frozenset()) -> set[str]:
        if identifier in events:
            return {identifier}
        if identifier not in groups or identifier in visited:
            raise ValueError(f"Unknown or cyclic historical group: {identifier}")
        return set().union(*(expand(member, visited | {identifier})
                             for member in groups[identifier]))

    pairs: set[tuple[str, str]] = set()
    evidence = []
    for constraint in _objects(reference.get("order_constraints", []), "Reference constraints"):
        if not isinstance(constraint.get("type"), str) or not constraint["type"]:
            raise ValueError("Reference constraints require a nonempty evidence type")
        if constraint.get("type") != "phylogenetic_branch_precedence":
            continue
        if constraint.get("logic") != "all_members_before_all_members":
            raise ValueError("Historical precedence must use all_members_before_all_members")
        before = set().union(*(expand(item) for item in _identifiers(
            constraint.get("prerequisites"), "Historical prerequisites")))
        after = set().union(*(expand(item) for item in _identifiers(
            constraint.get("dependent"), "Historical dependents")))
        pairs.update(itertools.product(before, after))
        evidence.append({key: constraint[key] for key in ("id", "source_ids", "locator", "scope")
                         if key in constraint})
    # Transitive closure makes all and only resolved historical pairs testable.
    for middle in events:
        predecessors = {first for first, last in pairs if last == middle}
        successors = {last for first, last in pairs if first == middle}
        pairs.update(itertools.product(predecessors, successors))
    if any(first == last or (last, first) in pairs for first, last in pairs):
        raise ValueError("Historical precedence contains a cycle")
    ordered_pairs = sorted(pairs)
    candidates = None
    if ordered_pairs and len(events) <= _MAX_CANDIDATE_EVENTS:
        candidates = [list(order) for order in itertools.permutations(events)
                      if all(order.index(first) < order.index(last) for first, last in ordered_pairs)]
    return {"events": events, "keys": keys, "pairs": ordered_pairs,
            "candidates": candidates, "evidence": evidence}


def _predictions(payload: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    predictions = []
    names: set[str] = set()
    units: set[int] = set()
    keys: set[tuple[int, str, str]] = set()
    columns = set()
    substitutions = _objects(payload.get("substitutions"), "Predicted substitutions")
    if len(substitutions) < 2:
        raise ValueError("Benchmark predictions require at least two substitution units")
    for unit in substitutions:
        members = _objects(unit.get("members"), "Predicted unit members")
        if len(members) != 1:
            raise ValueError("Blinded benchmark predictions must contain singleton units")
        member = members[0]
        name, index = unit.get("name"), unit.get("unit")
        column = member.get("column0")
        if (not isinstance(name, str) or not name or name in names or isinstance(index, bool)
                or not isinstance(index, int) or index < 0 or index in units
                or isinstance(column, bool) or not isinstance(column, int) or column < 0
                or column in columns):
            raise ValueError("Invalid or duplicate prediction unit")
        _direction(member.get("ancestral"), member.get("derived"), "Predicted member")
        if member.get("is_indel", False) or member.get("class", "nonsynonymous") != "nonsynonymous":
            raise ValueError("Benchmark predictions cannot contain indel or synonymous units")
        key = (column, member["ancestral"], member["derived"])
        if key in keys:
            raise ValueError("Duplicate predicted substitution key")
        names.add(name)
        units.add(index)
        columns.add(column)
        keys.add(key)
        predictions.append({"unit": index, "name": name, "column0": column,
                            "ancestral": member["ancestral"], "derived": member["derived"]})
    if units != set(range(len(predictions))):
        raise ValueError("Prediction unit indices must be contiguous from zero")
    before_matrix = payload.get("before_matrix")
    if not isinstance(before_matrix, dict):
        raise ValueError("Prediction before_matrix must be an object")
    matrix_names = before_matrix.get("names")
    if (not isinstance(matrix_names, list)
            or any(not isinstance(name, str) for name in matrix_names)
            or len(matrix_names) != len(predictions) or len(set(matrix_names)) != len(matrix_names)
            or set(matrix_names) != names):
        raise ValueError("Before-matrix names must identify every prediction unit once")
    return predictions, {name: index for index, name in enumerate(matrix_names)}


def _probability(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("Before matrix must contain finite probabilities")
    if value < -_PROBABILITY_TOLERANCE or value > 1 + _PROBABILITY_TOLERANCE:
        raise ValueError("Before-matrix probability lies outside [0, 1]")
    return min(1.0, max(0.0, float(value)))


def _evaluate_method(scorer: str, status: dict[str, Any], payload: dict[str, Any] | None,
                     truth: dict[str, Any]) -> dict[str, Any]:
    targets, pairs = truth["events"], truth["pairs"]
    base = {"scorer": scorer,
            "role": "diagnostic_baseline" if scorer == "dphi-product" else "scoring_method",
            "analysis_status": status["status"], "exit_code": status.get("exit_code"),
            "coverage": None,
            "historical_order": {"status": "not_testable", "reason": "analysis_failed"},
            "historical_candidates": {"status": "not_testable", "reason": "analysis_failed"}}
    if status["status"] != "completed":
        return base
    assert payload is not None
    predictions, matrix_index = _predictions(payload)
    matches = {}
    unmatched = []
    for prediction in predictions:
        key = (prediction["column0"], prediction["ancestral"], prediction["derived"])
        identifier = truth["keys"].get(key)
        if identifier is None:
            unmatched.append(prediction)
        else:
            matches[identifier] = prediction
    covered_pairs = [(first, last) for first, last in pairs if first in matches and last in matches]
    base["coverage"] = {
        "target_event_count": len(targets), "predicted_unit_count": len(predictions),
        "covered_target_event_count": len(matches),
        "target_event_fraction": len(matches) / len(targets),
        "ordered_pair_count": len(pairs), "covered_ordered_pair_count": len(covered_pairs),
        "ordered_pair_fraction": len(covered_pairs) / len(pairs) if pairs else None,
        "matched_events": [{"event_id": identifier, **matches[identifier]}
                           for identifier in targets if identifier in matches],
        "missing_event_ids": [identifier for identifier in targets if identifier not in matches],
        "unmatched_predictions": unmatched,
        "unmatched_predictions_are_false_positives": False,
    }
    accessible = payload.get("accessible")
    if not isinstance(accessible, bool):
        raise ValueError("Prediction accessible flag must be boolean")
    discrimination = payload.get("discrimination", {})
    provenance = payload.get("provenance", {})
    if not isinstance(discrimination, dict) or not isinstance(provenance, dict):
        raise ValueError("Prediction diagnostics and provenance must be objects")
    uniform = discrimination.get("degenerate_uniform", False)
    if not isinstance(uniform, bool):
        raise ValueError("Prediction degenerate_uniform flag must be boolean")
    matrix = payload.get("before_matrix", {}).get("matrix", [])
    map_order = payload.get("map_order", [])
    map_unidentified = uniform
    map_unidentified_reason = "uniform_predictions_have_no_identified_MAP_order" if uniform else None
    if not isinstance(matrix, list) or not isinstance(map_order, list):
        raise ValueError("Prediction matrix and MAP order must be lists")
    if not accessible and map_order:
        raise ValueError("An inaccessible prediction cannot claim a complete MAP order")
    if accessible:
        size = len(predictions)
        if len(matrix) != size or any(not isinstance(row, list) or len(row) != size for row in matrix):
            raise ValueError("Before matrix has incorrect dimensions")
        matrix = [[_probability(value) for value in row] for row in matrix]
        for i in range(size):
            for j in range(i):
                if abs(matrix[i][j] + matrix[j][i] - 1) > _PROBABILITY_TOLERANCE:
                    raise ValueError("Before-matrix opposing probabilities must sum to one")
        if (not isinstance(map_order, list) or len(map_order) != size
                or any(not isinstance(name, str) for name in map_order)
                or len(set(map_order)) != size or set(map_order) != set(matrix_index)):
            raise ValueError("Accessible prediction MAP must contain every unit exactly once")
        # Pairwise marginals alone cannot distinguish all path distributions.
        flat = size > 1 and all(
            abs(matrix[i][j] - 0.5) <= _PROBABILITY_TOLERANCE
            for i in range(size) for j in range(size) if i != j)
        probability_ratio = None
        if "map_prob_normalised" in payload or "uniform_baseline" in payload:
            map_probability = _probability(payload.get("map_prob_normalised"))
            baseline = _probability(payload.get("uniform_baseline"))
            if baseline <= 0 or not math.isclose(baseline, 1 / math.factorial(size),
                                                 rel_tol=_PROBABILITY_TOLERANCE, abs_tol=0):
                raise ValueError("Uniform baseline must equal one over the factorial of the unit count")
            probability_ratio = map_probability / baseline
        if uniform and (not flat or (probability_ratio is not None
                                     and abs(probability_ratio - 1) > _PROBABILITY_TOLERANCE)):
            raise ValueError("Uniform-distribution diagnostic contradicts the prediction probabilities")
        # Flat pairwise marginals alone do not prove a uniform path distribution:
        # higher-order path preferences can cancel in every pairwise marginal.
        # When path probabilities are absent, withhold MAP credit conservatively.
        uniform = uniform or (flat and probability_ratio is not None
                              and abs(probability_ratio - 1) <= _PROBABILITY_TOLERANCE)
        map_unidentified = uniform or (flat and probability_ratio is None)
        map_unidentified_reason = (
            "uniform_predictions_have_no_identified_MAP_order" if uniform else
            "flat_pairwise_predictions_without_path_probability_evidence" if map_unidentified else None
        )
    base["diagnostics"] = {"accessible": accessible, "degenerate_uniform": uniform,
                           "map_order_unidentified": map_unidentified,
                           "discrimination": payload.get("discrimination", {}),
                           "note": payload.get("note"),
                           "placement": provenance.get("genotype_placement")}
    historical = {"status": "not_testable", "reason": "insufficient_order_coverage",
                  "scope": "conditional_on_covered_ordered_pairs",
                  "full_order_coverage": bool(pairs) and len(covered_pairs) == len(pairs),
                  "mean_probability_correct_direction": None, "brier_score": None,
                  "uniform_brier_score": None, "brier_improvement_over_uniform": None,
                  "map_inversions": None, "map_inversion_fraction": None,
                  "map_metrics_reason": "insufficient_order_coverage", "pairs": []}
    if covered_pairs and not accessible:
        historical["reason"] = "no_accessible_path"
        historical["map_metrics_reason"] = "no_accessible_path"
    elif covered_pairs:
        rank = {name: index for index, name in enumerate(map_order)}
        entries = []
        for first, last in covered_pairs:
            first_name, last_name = matches[first]["name"], matches[last]["name"]
            probability = matrix[matrix_index[first_name]][matrix_index[last_name]]
            entries.append({"before": first, "after": last,
                            "probability_correct_direction": probability,
                            "map_inversion": None if map_unidentified else rank[first_name] > rank[last_name]})
        mean = sum(item["probability_correct_direction"] for item in entries) / len(entries)
        brier = sum((1 - item["probability_correct_direction"]) ** 2 for item in entries) / len(entries)
        historical.update(status="evaluated", reason=None,
                          mean_probability_correct_direction=mean, brier_score=brier,
                          uniform_brier_score=0.25, brier_improvement_over_uniform=0.25 - brier,
                          map_metrics_reason=map_unidentified_reason,
                          pairs=entries)
        if not map_unidentified:
            inversions = sum(item["map_inversion"] for item in entries)
            historical.update(map_inversions=inversions, map_inversion_fraction=inversions / len(entries))
    base["historical_order"] = historical

    candidates = truth["candidates"]
    candidate_comparison = {
        "status": "not_testable", "reason": "incomplete_target_coverage",
        "candidate_count": len(candidates) if candidates is not None else None,
        "projected_map_event_ids": None, "matches_historically_favored_path": None,
    }
    if len(matches) == len(targets):
        if not pairs:
            candidate_comparison["reason"] = "no_historical_order_constraints"
        elif not accessible:
            candidate_comparison["reason"] = "no_accessible_path"
        elif map_unidentified:
            candidate_comparison["reason"] = map_unidentified_reason
        elif candidates is None:
            candidate_comparison["reason"] = "reference_exceeds_candidate_enumeration_limit"
        else:
            by_name = {prediction["name"]: identifier for identifier, prediction in matches.items()}
            projected = [by_name[name] for name in map_order if name in by_name]
            candidate_comparison.update(status="evaluated", reason=None,
                                        projected_map_event_ids=projected,
                                        matches_historically_favored_path=projected in candidates)
    base["historical_candidates"] = candidate_comparison
    return base


def evaluate_run(run_dir: Path, reference_path: Path, output_dir: Path) -> dict[str, Any]:
    """Verify a complete run and write deterministic reference-aware assessments.

    Prediction and reference inputs are read-only. Existing nonempty output
    directories are rejected, and outputs cannot be placed within the sealed run.
    """
    supplied_run = Path(run_dir).absolute()
    if supplied_run.is_symlink():
        raise ValueError("Sealed run must not be a symlink")
    run_dir = supplied_run.resolve()
    reference_path, output_dir = (Path(value).resolve() for value in (reference_path, output_dir))
    if (_inside(output_dir, run_dir) or _inside(reference_path, output_dir)
            or _inside(reference_path, run_dir)):
        raise ValueError("Evaluation output must be separate from the sealed run and reference")
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise ValueError("Evaluation output directory must be absent or empty")
    verify_seal(run_dir)
    seal_hash = sha256_file(run_dir / "seal.json")
    reference_hash = sha256_file(reference_path)
    reference = _read_json(reference_path)
    protocol = _read_json(run_dir / "protocol.json")
    scorers = protocol.get("scorers")
    if protocol.get("schema_version") not in {"1.0", "2.0", "3.0"} or scorers != list(SCORERS):
        raise ValueError("Protocol must use schema 1.0, 2.0 or 3.0 and the four frozen scorers in order")
    if protocol["schema_version"] in {"2.0", "3.0"} and protocol.get("ancestral_method") != "hyphy-slac":
        raise ValueError("SLAC protocols require hyphy-slac ancestry")
    if protocol["schema_version"] == "3.0" and protocol.get("event_policy") != "panel_resolved_nonsynonymous":
        raise ValueError("Protocol 3.0 requires a manuscript substitution panel")
    alignment = reference.get("alignment", {})
    if not isinstance(alignment, dict):
        raise ValueError("Reference alignment metadata must be an object")
    reference_row = alignment.get("reference_row")
    if reference_row and protocol.get("focal_taxon") != reference_row:
        raise ValueError("Reference crosswalk row does not match the protocol focal taxon")
    truth = _reference_order(reference)
    methods = []
    for scorer in scorers:
        method_dir = run_dir / "methods" / scorer
        status = _read_json(method_dir / "status.json")
        if (status.get("scorer") != scorer or not isinstance(status.get("status"), str)
                or status["status"] not in {"completed", "failed"}):
            raise ValueError(f"Invalid analysis status for {scorer}")
        if status["status"] == "completed" and (type(status.get("exit_code")) is not int
                                                    or status["exit_code"] != 0):
            raise ValueError(f"Completed analysis must have exit code zero: {scorer}")
        payload = _read_json(method_dir / "order.json") if status["status"] == "completed" else None
        methods.append(_evaluate_method(scorer, status, payload, truth))
    # Detect an accidental reference edit during assessment; no inference is ever run.
    if sha256_file(reference_path) != reference_hash:
        raise ValueError("Reference changed during evaluation")
    verify_seal(run_dir)
    if sha256_file(run_dir / "seal.json") != seal_hash:
        raise ValueError("Run seal changed during evaluation")
    report = {
        "schema_version": "1.0", "study_id": reference.get("study_id"),
        "evaluation_type": "phylogenetic_branch_precedence",
        "analysis_protocol": {
            "schema_version": protocol["schema_version"],
            "ancestral_method": "nonfocal-plurality" if protocol["schema_version"] == "1.0" else "hyphy-slac",
            "asr_node": protocol.get("asr_node"),
        },
        "reference_sha256": reference_hash, "run_seal_sha256": seal_hash,
        "source_evidence": reference.get("sources", []),
        "empirical_order": {
            "target_events": list(truth["events"].values()),
            "ordered_pairs": [list(pair) for pair in truth["pairs"]],
            "historical_candidates": truth["candidates"],
            "constraint_evidence": truth["evidence"],
            "unresolved_pairs_are_scored": False,
        },
        "methods": methods,
        "assay_accessibility": {"status": "deferred", "reason": "not_scored_by_this_historical_order_evaluator",
                                "reference_status": reference.get("assay_data", {}).get("scoring_status")},
        "limitations": [
            "Historical branch precedence is phylogenetic evidence, not proof of biochemical necessity.",
            "Accuracy is conditional on covered directed substitutions; coverage is reported separately.",
            "Unmatched predictions are not automatically false positives.",
            "Uniform predictions receive probability-baseline scores and no MAP-order credit.",
            "These assessments do not establish checkpoint training-data provenance.",
        ],
    }
    if protocol["schema_version"] == "3.0":
        report["analysis_protocol"].update(event_policy=protocol["event_policy"],
                                            substitution_panel_sha256=protocol["input_hashes"]["substitution_panel"])
        report["limitations"].append(
            "Manuscript sites were supplied in advance; this tests ordering within that panel, not site discovery. "
            "Excluded columns are fixed at focal codons, so the starting background is not the full SLAC ancestor or the assayed ancestral construct."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "evaluation.json", report)
    for method in methods:
        write_json(output_dir / "methods" / f"{method['scorer']}.json", method)
    fields = ["scorer", "analysis_status", "target_events_covered", "target_events_total",
              "ordered_pairs_covered", "ordered_pairs_total", "historical_order_status", "reason",
              "mean_probability_correct_direction", "brier_score", "uniform_brier_score",
              "brier_improvement_over_uniform", "map_inversions", "historical_candidate_match"]
    with (output_dir / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for method in methods:
            coverage = method["coverage"] or {}
            metrics = method["historical_order"]
            writer.writerow({
                "scorer": method["scorer"], "analysis_status": method["analysis_status"],
                "target_events_covered": coverage.get("covered_target_event_count"),
                "target_events_total": coverage.get("target_event_count", len(truth["events"])),
                "ordered_pairs_covered": coverage.get("covered_ordered_pair_count"),
                "ordered_pairs_total": coverage.get("ordered_pair_count", len(truth["pairs"])),
                "historical_order_status": metrics["status"], "reason": metrics.get("reason"),
                **{key: metrics.get(key) for key in fields[8:13]},
                "historical_candidate_match": method["historical_candidates"]["matches_historically_favored_path"]
                if "matches_historically_favored_path" in method["historical_candidates"] else None,
            })
    output_hashes = {path.relative_to(output_dir).as_posix(): sha256_file(path)
                     for path in sorted(output_dir.rglob("*")) if path.is_file()}
    write_json(output_dir / "evaluation_manifest.json", {
        "schema_version": "1.0", "run_seal_sha256": seal_hash,
        "protocol_sha256": sha256_file(run_dir / "protocol.json"),
        "reference_sha256": reference_hash,
        "evaluator_code_sha256": sha256_file(Path(__file__)),
        "evaluator_dependencies_sha256": {
            "benchmark_artifacts.py": sha256_file(Path(__file__).with_name("benchmark_artifacts.py"))
        },
        "outputs": output_hashes,
    })
    return report


def cmd_benchmark_evaluate(args: Any) -> int:
    """Thin CLI adapter; it accepts no analysis settings or model inputs."""
    report = evaluate_run(Path(args.run), Path(args.reference), Path(args.output))
    print(f"Evaluated {len(report['methods'])} sealed scoring methods; wrote {args.output}")
    return 0
