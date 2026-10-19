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
from epistaeon.benchmark_run import run_benchmark


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
    runs = [run_benchmark(alignment, "focus", weights, tmp_path / f"run_{index}")
            for index in range(2)]
    for run in runs:
        verify_seal(run)
        assert json.loads((run / "run_status.json").read_text())["status"] == "completed"
        isolation = json.loads((run / "isolation_check.json").read_text())
        assert isolation["status"] == "pass"
        assert isolation["active_network_interfaces"] == ["lo"]
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
