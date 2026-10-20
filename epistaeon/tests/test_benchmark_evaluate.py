"""Reference-aware scoring never imports or invokes inference code."""

import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon.benchmark_artifacts import seal_run, sha256_file, write_json
from epistaeon.benchmark_evaluate import evaluate_run


SCORERS = ["potts-metropolis", "typicality-gate", "softmax-repaired", "dphi-product"]
ROOT = Path(__file__).resolve().parents[2]


def _reference():
    def event(identifier, col, ancestral, derived):
        return {"id": identifier, "from": ancestral, "to": derived, "order_scoring": True,
                "structure_mapping": [{"alignment_reference_columns_1based": {"TEST": col + 1}}]}
    return {
        "study_id": "synthetic", "alignment": {"reference_row": "focus"},
        "events": [event("A", 0, "A", "G"), event("B", 1, "S", "T"), event("C", 2, "F", "L")],
        "groups": [{"id": "early", "members": ["A", "B"]}],
        "order_constraints": [
            {"id": "early_before_C", "type": "phylogenetic_branch_precedence",
             "prerequisites": ["early"], "dependent": ["C"], "logic": "all_members_before_all_members"},
            # This deliberately contradicts A-before-C, but it is an assay rule
            # rather than phylogenetic historical precedence and must be ignored.
            {"id": "assay_only", "type": "experimental_single_mutant_exclusion",
             "prerequisites": ["C"], "dependent": ["A"], "logic": "at_least_one"},
        ],
    }


def _payload(predictions=None, probabilities=None, order=None, uniform=False, accessible=True):
    predictions = predictions or [(0, "A", "G"), (1, "S", "T"), (2, "F", "L")]
    names = [f"u{i}" for i in range(len(predictions))]
    if probabilities is None:
        probabilities = [[0 if i == j else (0.8 if i < j else 0.2)
                          for j in range(len(names))] for i in range(len(names))]
    return {
        "substitutions": [{"unit": i, "name": names[i],
                           "members": [{"column0": col, "ancestral": ancestral, "derived": derived}]}
                          for i, (col, ancestral, derived) in enumerate(predictions)],
        "accessible": accessible, "map_order": names if order is None else order,
        "before_matrix": {"names": names, "matrix": probabilities},
        "discrimination": {"degenerate_uniform": uniform},
    }


def _run(tmp_path, payload=None, reference=None, statuses=None, protocol=None):
    run_dir = tmp_path / "run"
    reference_path = tmp_path / "truth.json"
    write_json(run_dir / "protocol.json", {"schema_version": "1.0", "scorers": SCORERS,
                                           "focal_taxon": "focus", **(protocol or {})})
    write_json(reference_path, _reference() if reference is None else reference)
    for scorer in SCORERS:
        status = {"scorer": scorer, "status": "completed", "exit_code": 0}
        if statuses and scorer in statuses:
            status.update(statuses[scorer])
        method = run_dir / "methods" / scorer
        write_json(method / "status.json", status)
        if status["status"] == "completed":
            write_json(method / "order.json", _payload() if payload is None else payload)
        (method / "stdout.log").write_text("analysis log\n")
        (method / "stderr.log").write_text("")
    seal_run(run_dir)
    return run_dir, reference_path


def _evaluate(tmp_path, **kwargs):
    run_dir, reference_path = _run(tmp_path, **kwargs)
    return evaluate_run(run_dir, reference_path, tmp_path / "evaluation")


def test_known_partial_order_metrics_and_candidates(tmp_path):
    report = _evaluate(tmp_path)
    assert report["analysis_protocol"]["ancestral_method"] == "nonfocal-plurality"
    assert report["empirical_order"]["ordered_pairs"] == [["A", "C"], ["B", "C"]]
    assert report["empirical_order"]["historical_candidates"] == [["A", "B", "C"], ["B", "A", "C"]]
    assert len(report["methods"]) == 4
    for method in report["methods"]:
        assert method["coverage"]["covered_target_event_count"] == 3
        assert method["coverage"]["covered_ordered_pair_count"] == 2
        metrics = method["historical_order"]
        assert metrics["status"] == "evaluated"
        assert metrics["mean_probability_correct_direction"] == pytest.approx(0.8)
        assert metrics["brier_score"] == pytest.approx(0.04)
        assert metrics["brier_improvement_over_uniform"] == pytest.approx(0.21)
        assert metrics["map_inversions"] == 0
        assert method["historical_candidates"]["matches_historically_favored_path"] is True
    assert report["methods"][-1]["role"] == "diagnostic_baseline"
    assert report["assay_accessibility"]["status"] == "deferred"
    assert (tmp_path / "evaluation" / "comparison.csv").read_text().count("\n") == 5
    manifest = json.loads((tmp_path / "evaluation" / "evaluation_manifest.json").read_text())
    for filename, digest in manifest["outputs"].items():
        assert sha256_file(tmp_path / "evaluation" / filename) == digest


def test_slac_protocol_is_reported_separately_and_unrun_methods_remain_unscored(tmp_path):
    report = _evaluate(tmp_path, protocol={"schema_version": "2.0", "ancestral_method": "hyphy-slac", "asr_node": "Anc"},
                       statuses={scorer: {"status": "failed", "exit_code": None,
                                          "analysis_executed": False, "reason": "exact_budget_exceeded"}
                                 for scorer in SCORERS})
    assert report["analysis_protocol"] == {"schema_version": "2.0", "ancestral_method": "hyphy-slac", "asr_node": "Anc"}
    for method in report["methods"]:
        assert method["historical_order"] == {"status": "not_testable", "reason": "analysis_failed"}


def test_reversed_direction_and_unmatched_events_are_not_false_positives(tmp_path):
    report = _evaluate(tmp_path, payload=_payload([(0, "G", "A"), (1, "S", "T"), (2, "F", "L"),
                                                 (4, "Q", "N")]))
    method = report["methods"][0]
    coverage = method["coverage"]
    assert coverage["covered_target_event_count"] == 2
    assert coverage["missing_event_ids"] == ["A"]
    assert coverage["covered_ordered_pair_count"] == 1
    assert len(coverage["unmatched_predictions"]) == 2
    assert coverage["unmatched_predictions_are_false_positives"] is False
    assert method["historical_order"]["full_order_coverage"] is False
    assert method["historical_order"]["mean_probability_correct_direction"] == pytest.approx(0.8)
    assert method["historical_candidates"]["reason"] == "incomplete_target_coverage"


def test_only_unresolved_pairs_have_coverage_but_no_order_score(tmp_path):
    method = _evaluate(tmp_path, payload=_payload([(0, "A", "G"), (1, "S", "T")]))["methods"][0]
    assert method["coverage"]["covered_target_event_count"] == 2
    assert method["coverage"]["covered_ordered_pair_count"] == 0
    assert method["historical_order"]["status"] == "not_testable"
    assert method["historical_order"]["reason"] == "insufficient_order_coverage"
    assert method["historical_order"]["brier_score"] is None
    assert method["historical_order"]["map_inversions"] is None


def test_one_unit_has_event_coverage_but_no_testable_order(tmp_path):
    report = _evaluate(tmp_path, payload=_payload([(0, "A", "G")]))
    for method in report["methods"]:
        assert method["coverage"]["covered_target_event_count"] == 1
        assert method["coverage"]["covered_ordered_pair_count"] == 0
        assert method["historical_order"]["status"] == "not_testable"
        assert method["historical_order"]["reason"] == "insufficient_order_coverage"


@pytest.mark.parametrize("diagnostic_flag", [False, True])
def test_uniform_predictions_get_baseline_but_no_map_credit(tmp_path, diagnostic_flag):
    matrix = [[0 if i == j else 0.5 for j in range(3)] for i in range(3)]
    method = _evaluate(tmp_path, payload=_payload(probabilities=matrix,
                                                  uniform=diagnostic_flag))["methods"][0]
    metrics = method["historical_order"]
    assert metrics["mean_probability_correct_direction"] == 0.5
    assert metrics["brier_score"] == 0.25
    assert metrics["brier_improvement_over_uniform"] == 0
    assert metrics["map_inversions"] is None
    assert metrics["map_inversion_fraction"] is None
    assert all(pair["map_inversion"] is None for pair in metrics["pairs"])
    assert method["historical_candidates"]["matches_historically_favored_path"] is None
    expected_reason = ("uniform_predictions_have_no_identified_MAP_order" if diagnostic_flag else
                       "flat_pairwise_predictions_without_path_probability_evidence")
    assert method["historical_candidates"]["reason"] == expected_reason


def test_map_inversions_and_historical_candidate_mismatch(tmp_path):
    method = _evaluate(tmp_path, payload=_payload(order=["u2", "u1", "u0"]))["methods"][0]
    assert method["historical_order"]["map_inversions"] == 2
    assert method["historical_order"]["map_inversion_fraction"] == 1
    assert method["historical_candidates"]["matches_historically_favored_path"] is False


@pytest.mark.parametrize("covered,reason", [(True, "no_accessible_path"),
                                           (False, "insufficient_order_coverage")])
def test_blocked_paths_preserve_coverage_without_order_metrics(tmp_path, covered, reason):
    payload = _payload(accessible=False, order=[], probabilities=[])
    if not covered:
        payload = _payload([(0, "A", "G"), (1, "S", "T")], accessible=False,
                           order=[], probabilities=[])
    method = _evaluate(tmp_path, payload=payload)["methods"][0]
    assert method["coverage"]["covered_target_event_count"] == (3 if covered else 2)
    assert method["historical_order"]["status"] == "not_testable"
    assert method["historical_order"]["reason"] == reason
    assert method["historical_order"]["mean_probability_correct_direction"] is None


def test_failed_method_has_no_fabricated_coverage_or_order(tmp_path):
    report = _evaluate(tmp_path, statuses={"dphi-product": {"status": "failed", "exit_code": 2}})
    failed = report["methods"][-1]
    assert failed["analysis_status"] == "failed"
    assert failed["coverage"] is None
    assert failed["historical_order"] == {"status": "not_testable", "reason": "analysis_failed"}
    assert report["methods"][0]["historical_order"]["status"] == "evaluated"


@pytest.mark.parametrize("alteration", ["changed", "added", "removed", "symlink"])
def test_checksum_rejection_precedes_output(tmp_path, alteration):
    run_dir, reference_path = _run(tmp_path)
    path = run_dir / "methods" / SCORERS[0] / "stdout.log"
    if alteration == "changed":
        path.write_text("changed\n")
    elif alteration == "added":
        (run_dir / "unsealed.txt").write_text("extra")
    elif alteration == "removed":
        path.unlink()
    else:
        path.unlink()
        path.symlink_to(reference_path)
    with pytest.raises(ValueError):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")
    assert not (tmp_path / "evaluation").exists()


def test_changed_evaluator_truth_leaves_sealed_analysis_unchanged(tmp_path):
    run_dir, reference_path = _run(tmp_path)
    original_seal = sha256_file(run_dir / "seal.json")
    original_artifacts = {p.relative_to(run_dir): sha256_file(p)
                          for p in run_dir.rglob("*") if p.is_file()}
    first = evaluate_run(run_dir, reference_path, tmp_path / "first")
    truth = _reference()
    truth["order_constraints"][0].update(prerequisites=["C"], dependent=["early"])
    write_json(reference_path, truth)
    second = evaluate_run(run_dir, reference_path, tmp_path / "second")
    assert first["methods"][0]["historical_order"]["map_inversions"] == 0
    assert second["methods"][0]["historical_order"]["map_inversions"] == 2
    assert original_seal == sha256_file(run_dir / "seal.json")
    assert original_artifacts == {p.relative_to(run_dir): sha256_file(p)
                                  for p in run_dir.rglob("*") if p.is_file()}


def test_evaluation_is_deterministic(tmp_path):
    run_dir, reference_path = _run(tmp_path)
    first = evaluate_run(run_dir, reference_path, tmp_path / "first")
    second = evaluate_run(run_dir, reference_path, tmp_path / "second")
    assert first == second
    assert (tmp_path / "first" / "comparison.csv").read_bytes() == (tmp_path / "second" / "comparison.csv").read_bytes()
    assert (tmp_path / "first" / "evaluation_manifest.json").read_bytes() == (tmp_path / "second" / "evaluation_manifest.json").read_bytes()


def test_rejects_outputs_in_run_or_containing_reference(tmp_path):
    run_dir, reference_path = _run(tmp_path)
    for output in [run_dir / "evaluation", reference_path, tmp_path]:
        with pytest.raises(ValueError):
            evaluate_run(run_dir, reference_path, output)


def test_rejects_symlink_run(tmp_path):
    run_dir, reference_path = _run(tmp_path)
    link = tmp_path / "run-link"
    link.symlink_to(run_dir, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        evaluate_run(link, reference_path, tmp_path / "evaluation")


def test_rejects_focal_crosswalk_mismatch(tmp_path):
    truth = _reference()
    truth["alignment"]["reference_row"] = "different"
    run_dir, reference_path = _run(tmp_path, reference=truth)
    with pytest.raises(ValueError, match="focal"):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")


def test_rejects_cyclic_history(tmp_path):
    truth = _reference()
    reverse = copy.deepcopy(truth["order_constraints"][0])
    reverse.update(id="reverse", prerequisites=["C"], dependent=["early"])
    truth["order_constraints"].append(reverse)
    run_dir, reference_path = _run(tmp_path, reference=truth)
    with pytest.raises(ValueError, match="cycle"):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")


def test_sws1_plurality_acceptance_without_model_or_reference_based_selection(tmp_path):
    from epistaeon.benchmark_inputs import select_benchmark_inputs

    selection = select_benchmark_inputs(ROOT / "epistaeon/data/alignments/sws1/SWS1.codonified.fa",
                                        "Homo_sapiens")
    assert selection["n_units"] == 10
    assert selection["n_genotypes"] == 1024
    predictions = [(event["column0"], event["ancestral"], event["derived"])
                   for event in selection["selected"]]
    truth = json.loads((ROOT / "epistaeon/data/substitution_reference/sws1.json").read_text())
    truth["alignment"]["reference_row"] = "focus"  # fixture protocol row only
    report = _evaluate(tmp_path, payload=_payload(predictions), reference=truth)
    assert len(report["empirical_order"]["ordered_pairs"]) == 18
    assert len(report["empirical_order"]["historical_candidates"]) == 8
    for method in report["methods"]:
        assert method["coverage"]["target_event_count"] == 7
        assert method["coverage"]["covered_target_event_count"] == 2
        assert method["coverage"]["ordered_pair_count"] == 18
        assert method["coverage"]["covered_ordered_pair_count"] == 0
        assert method["historical_order"]["status"] == "not_testable"
        assert method["historical_order"]["reason"] == "insufficient_order_coverage"
        assert method["historical_order"]["mean_probability_correct_direction"] is None


def test_evaluator_imports_no_model_or_numerical_runtime():
    src = Path(__file__).resolve().parents[1] / "src"
    script = """
import sys
sys.path.insert(0, sys.argv[1])
import epistaeon.benchmark_evaluate
forbidden = ('torch', 'numpy', 'hyphaeon', 'epistaeon.scoring', 'epistaeon.trajectory',
             'epistaeon.benchmark_runner', 'epistaeon.benchmark_worker')
assert not any(name == prefix or name.startswith(prefix + '.')
               for name in sys.modules for prefix in forbidden), sorted(sys.modules)
"""
    completed = subprocess.run([sys.executable, "-I", "-c", script, str(src)],
                               capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr


def test_historical_pairs_include_transitive_closure(tmp_path):
    truth = _reference()
    truth["order_constraints"] = [
        {"id": "A_before_B", "type": "phylogenetic_branch_precedence",
         "prerequisites": ["A"], "dependent": ["B"], "logic": "all_members_before_all_members"},
        {"id": "B_before_C", "type": "phylogenetic_branch_precedence",
         "prerequisites": ["B"], "dependent": ["C"], "logic": "all_members_before_all_members"},
    ]
    report = _evaluate(tmp_path, reference=truth)
    assert report["empirical_order"]["ordered_pairs"] == [["A", "B"], ["A", "C"], ["B", "C"]]
    assert report["empirical_order"]["historical_candidates"] == [["A", "B", "C"]]


@pytest.mark.parametrize("problem", ["shape", "probability", "asymmetry", "duplicate_map", "nonsingleton"])
def test_rejects_malformed_completed_predictions_without_partial_outputs(tmp_path, problem):
    payload = _payload()
    if problem == "shape":
        payload["before_matrix"]["matrix"].pop()
    elif problem == "probability":
        payload["before_matrix"]["matrix"][0][1] = 1.5
    elif problem == "asymmetry":
        payload["before_matrix"]["matrix"][0][1] = 0.3
    elif problem == "duplicate_map":
        payload["map_order"] = ["u0", "u0", "u2"]
    else:
        payload["substitutions"][0]["members"].append({"column0": 3, "ancestral": "N", "derived": "H"})
    run_dir, reference_path = _run(tmp_path, payload=payload)
    with pytest.raises(ValueError):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")
    assert not (tmp_path / "evaluation").exists()


def test_no_historical_constraints_cannot_earn_a_vacuous_candidate_match(tmp_path):
    truth = _reference()
    truth["order_constraints"] = []
    method = _evaluate(tmp_path, reference=truth)["methods"][0]
    assert method["historical_order"]["status"] == "not_testable"
    assert method["historical_candidates"]["status"] == "not_testable"
    assert method["historical_candidates"]["reason"] == "no_historical_order_constraints"
    assert method["historical_candidates"]["matches_historically_favored_path"] is None


def test_output_run_containment_is_checked_after_parent_alias_resolution(tmp_path):
    run_dir, reference_path = _run(tmp_path)
    parent_alias = tmp_path / "parent-alias"
    parent_alias.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="separate"):
        evaluate_run(parent_alias / "run", reference_path, run_dir / "evaluation")
    assert not (run_dir / "evaluation").exists()


@pytest.mark.parametrize("problem", ["events_container", "order_flag", "direction", "crosswalk",
                                      "group_members", "empty_prerequisites", "constraint_type",
                                      "assay_container", "source_container"])
def test_malformed_references_fail_clearly_without_output(tmp_path, problem):
    truth = _reference()
    if problem == "events_container":
        truth["events"] = {"A": truth["events"][0]}
    elif problem == "order_flag":
        truth["events"][0]["order_scoring"] = "false"
    elif problem == "direction":
        truth["events"][0]["from"] = ["A"]
    elif problem == "crosswalk":
        truth["events"][0]["structure_mapping"][0]["alignment_reference_columns_1based"] = {"TEST": [1]}
    elif problem == "group_members":
        truth["groups"][0]["members"] = "A"
    elif problem == "empty_prerequisites":
        truth["order_constraints"][0]["prerequisites"] = []
    elif problem == "constraint_type":
        truth["order_constraints"][0]["type"] = {}
    elif problem == "assay_container":
        truth["assay_data"] = []
    else:
        truth["sources"] = "a paper"
    run_dir, reference_path = _run(tmp_path, reference=truth)
    with pytest.raises(ValueError):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")
    assert not (tmp_path / "evaluation").exists()


@pytest.mark.parametrize("problem", ["substitution_container", "missing_direction", "duplicate_column",
                                      "diagnostic_container", "uniform_diagnostic", "map_container",
                                      "matrix_names", "synonymous", "indel"])
def test_malformed_prediction_metadata_fail_clearly(tmp_path, problem):
    payload = _payload()
    if problem == "substitution_container":
        payload["substitutions"] = {}
    elif problem == "missing_direction":
        del payload["substitutions"][0]["members"][0]["ancestral"]
    elif problem == "duplicate_column":
        payload["substitutions"][1]["members"][0]["column0"] = 0
    elif problem == "diagnostic_container":
        payload["discrimination"] = []
    elif problem == "uniform_diagnostic":
        payload["discrimination"]["degenerate_uniform"] = True
    elif problem == "map_container":
        payload["map_order"] = {}
    elif problem == "matrix_names":
        payload["before_matrix"]["names"] = None
    elif problem == "synonymous":
        payload["substitutions"][0]["members"][0]["derived"] = "A"
    else:
        payload["substitutions"][0]["members"][0]["is_indel"] = True
    run_dir, reference_path = _run(tmp_path, payload=payload)
    with pytest.raises(ValueError):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")
    assert not (tmp_path / "evaluation").exists()


def test_flat_pairwise_marginals_do_not_prove_uniform_paths(tmp_path):
    matrix = [[0 if i == j else 0.5 for j in range(3)] for i in range(3)]
    payload = _payload(probabilities=matrix)
    # Equal mass on an order and its reverse gives flat pairwise marginals
    # while assigning three times the uniform probability to either MAP path.
    payload.update(map_prob_normalised=1 / 2, uniform_baseline=1 / 6)
    method = _evaluate(tmp_path, payload=payload)["methods"][0]
    assert method["historical_order"]["brier_score"] == 0.25
    assert method["historical_order"]["map_inversions"] == 0
    assert method["historical_candidates"]["matches_historically_favored_path"] is True


def test_uniform_map_probability_withholds_map_credit(tmp_path):
    matrix = [[0 if i == j else 0.5 for j in range(3)] for i in range(3)]
    payload = _payload(probabilities=matrix)
    payload.update(map_prob_normalised=1 / 6, uniform_baseline=1 / 6)
    method = _evaluate(tmp_path, payload=payload)["methods"][0]
    assert method["historical_order"]["map_inversions"] is None
    assert method["historical_candidates"]["matches_historically_favored_path"] is None


@pytest.mark.parametrize("seal", [[], "not a seal", None])
def test_malformed_seal_top_level_fails_clearly(tmp_path, seal):
    run_dir, reference_path = _run(tmp_path)
    write_json(run_dir / "seal.json", seal)
    with pytest.raises(ValueError):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")


@pytest.mark.parametrize("text", ['{"events": [], "events": []}', '{"events": NaN}'])
def test_ambiguous_or_nonstandard_json_reference_is_rejected(tmp_path, text):
    run_dir, reference_path = _run(tmp_path)
    reference_path.write_text(text)
    with pytest.raises(ValueError):
        evaluate_run(run_dir, reference_path, tmp_path / "evaluation")
