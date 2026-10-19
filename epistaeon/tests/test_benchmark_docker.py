"""Opt-in, real-checkpoint integration under the pinned offline Docker runtime.

EPISTAEON_DOCKER_TEST=1 PYTHONPATH=epistaeon/src python -m pytest \
    epistaeon/tests/test_benchmark_docker.py -q

This small synthetic alignment tests execution and repeatability, not SWS1
accuracy. A running Docker daemon and the local checkpoint are required.
"""

import json
import math
import os
from pathlib import Path
import random
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon.benchmark_artifacts import SCORERS, sha256_file, verify_seal, write_json
from epistaeon.benchmark_evaluate import evaluate_run
from epistaeon.benchmark_run import run_benchmark, BenchmarkRunError


pytestmark = pytest.mark.skipif(os.environ.get("EPISTAEON_DOCKER_TEST") != "1",
                                reason="opt-in Docker/checkpoint integration")
ROOT = Path(__file__).resolve().parents[2]


def synthetic_alignment(path):
    rng = random.Random(12)
    codons = ["GCT", "GGT", "TTT", "CTG", "ATC", "ACC", "CAA", "GAG", "AAC", "TGC"]
    base = [rng.choice(codons) for _ in range(90)]
    base[8], base[49] = "GCT", "GGT"
    rows = {}
    for taxon in range(5):
        row = list(base)
        # Changes confined to separate sites ensure the plurality remains base.
        for col in range(2 + taxon, 85, 13):
            if col not in (8, 49):
                row[col] = codons[(codons.index(row[col]) + taxon + 1) % len(codons)]
        rows[f"background_{taxon}"] = row
    focus = list(base)
    focus[8], focus[49] = "GAT", "AGT"
    rows["focus"] = focus
    path.write_text("".join(f">{name}\n{''.join(row)}\n" for name, row in rows.items()))


def assert_numerically_equal(first, second):
    if isinstance(first, dict):
        assert first.keys() == second.keys()
        for key in first:
            assert_numerically_equal(first[key], second[key])
    elif isinstance(first, list):
        assert len(first) == len(second)
        for left, right in zip(first, second):
            assert_numerically_equal(left, right)
    elif isinstance(first, float) and math.isfinite(first):
        assert first == pytest.approx(second, abs=1e-7, rel=1e-6)
    else:
        assert first == second


def test_real_four_scorer_run_repeatability_and_evaluation(tmp_path):
    alignment = tmp_path / "synthetic.fa"
    synthetic_alignment(alignment)
    weights = ROOT / "epistaeon/data/model/model.safetensors"
    assert weights.is_file(), "Integration requires the local checkpoint; no downloads or native fallback"
    tree = tmp_path / "tree.nwk"
    tree.write_text("(background_0,background_1,background_2,(background_3,background_4,focus)Anc)Root;")
    panel = tmp_path / "panel.json"
    write_json(panel, {"schema_version": "1.0", "alignment_sha256": sha256_file(alignment),
                        "focal_taxon": "focus", "columns0": [8, 49]})
    runs = [run_benchmark(alignment, "focus", weights, tmp_path / f"run_{index}", asr_tree=tree, asr_node="Anc", substitution_panel=panel)
            for index in range(2)]
    for run in runs:
        verify_seal(run)
        assert json.loads((run / "run_status.json").read_text())["status"] == "completed"
        isolation = json.loads((run / "isolation_check.json").read_text())
        assert isolation["status"] == "pass"
        assert isolation["active_network_interfaces"] == ["lo"]
        assert "substitution_panel.json" in isolation["input_files"]
        selection = json.loads((run / "selection.json").read_text())
        assert selection["n_units"] == 2
        assert selection["n_genotypes"] == 4
        assert selection["full_mask_reproduces_focal"] is True
    for scorer in SCORERS:
        payloads = [json.loads((run / "methods" / scorer / "order.json").read_text()) for run in runs]
        assert_numerically_equal(payloads[0], payloads[1])

    # Truth is first created after both runs are sealed.
    reference = tmp_path / "truth.json"
    truth = {
        "study_id": "synthetic_integration", "alignment": {"reference_row": "focus"},
        "events": [{"id": name, "from": ancestral, "to": derived,
                    "structure_mapping": [{"alignment_reference_columns_1based": {"synthetic": col}}]}
                   for name, ancestral, derived, col in [("A", "A", "D", 9), ("B", "G", "S", 50)]],
        "order_constraints": [{"id": "A_before_B", "type": "phylogenetic_branch_precedence",
                               "prerequisites": ["A"], "dependent": ["B"],
                               "logic": "all_members_before_all_members"}],
    }
    write_json(reference, truth)
    before = {p.relative_to(runs[0]).as_posix(): sha256_file(p)
              for p in runs[0].rglob("*") if p.is_file()}
    report = evaluate_run(runs[0], reference, tmp_path / "evaluation")
    assert all(method["coverage"]["covered_target_event_count"] == 2 for method in report["methods"])
    truth["order_constraints"][0].update(prerequisites=["B"], dependent=["A"])
    write_json(reference, truth)
    evaluate_run(runs[0], reference, tmp_path / "changed_truth_evaluation")
    after = {p.relative_to(runs[0]).as_posix(): sha256_file(p)
             for p in runs[0].rglob("*") if p.is_file()}
    assert before == after


def test_sws1_manuscript_panel_runs_all_four_scorers(tmp_path):
    alignment = ROOT / "epistaeon/data/alignments/sws1/SWS1.codonified.fa"
    data = ROOT / "epistaeon/validation/sws1_slac/inputs"
    output = run_benchmark(alignment, "Homo_sapiens", ROOT / "epistaeon/data/model/model.safetensors", tmp_path / "sws1_panel",
                           asr_alignment=data / "SWS1.outgroups.codon.fa", asr_tree=data / "species_scaffold.nwk",
                           asr_node="Boreoeutheria", substitution_panel=alignment.parent / "substitution_panel.json")
    verify_seal(output)
    selection = json.loads((output / "selection.json").read_text())
    assert selection["unrestricted_nonsynonymous_units"] == 28
    assert selection["n_units"] == 7 and selection["n_genotypes"] == 128
    assert selection["full_mask_reproduces_focal"] and selection["unselected_panel_columns0"] == []
    assert json.loads((output / "run_status.json").read_text())["status"] == "completed"
    for scorer in SCORERS:
        assert json.loads((output / "methods" / scorer / "selection.json").read_text()) == selection
    # The reference and its historical order are first read after sealing.
    report = evaluate_run(output, ROOT / "epistaeon/data/substitution_reference/sws1.json", tmp_path / "evaluation")
    for method in report["methods"]:
        assert method["coverage"]["covered_target_event_count"] == 7
        assert method["coverage"]["covered_ordered_pair_count"] == 18
    assert report["analysis_protocol"]["schema_version"] == "3.0"


def test_sws1_reconstruction_is_sealed_before_budget_rejection(tmp_path):
    alignment = ROOT / "epistaeon/data/alignments/sws1/SWS1.codonified.fa"
    data = ROOT / "epistaeon/validation/sws1_slac/inputs"
    output = tmp_path / "sws1"
    with pytest.raises(BenchmarkRunError, match="268435456"):
        run_benchmark(alignment, "Homo_sapiens", ROOT / "epistaeon/data/model/model.safetensors", output,
                      asr_alignment=data / "SWS1.outgroups.codon.fa", asr_tree=data / "species_scaffold.nwk",
                      asr_node="Boreoeutheria")
    verify_seal(output)
    assert json.loads((output / "run_status.json").read_text())["status"] == "preflight_failed"
    selection = json.loads((output / "selection.json").read_text())
    assert selection["n_units"] == 28 and selection["n_genotypes"] == 268435456
    assert (output / "asr/ancestor.codon.fa").is_file()
    for scorer in SCORERS:
        status = json.loads((output / "methods" / scorer / "status.json").read_text())
        assert status["analysis_executed"] is False
        assert not (output / "methods" / scorer / "order.json").exists()
    # Truth enters only after the isolated preparation output has been sealed.
    from epistaeon.benchmark_evaluate import _reference_order
    truth = _reference_order(json.loads((ROOT / "epistaeon/data/substitution_reference/sws1.json").read_text()))
    keys = {(e["column0"], e["ancestral"], e["derived"]) for e in selection["selected"]}
    assert len(truth["keys"]) == 7 and set(truth["keys"]) <= keys
    report = evaluate_run(output, ROOT / "epistaeon/data/substitution_reference/sws1.json", tmp_path / "evaluation")
    assert all(method["historical_order"]["status"] == "not_testable" for method in report["methods"])
