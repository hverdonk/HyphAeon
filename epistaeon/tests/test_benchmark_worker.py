"""Frozen worker contract; scorer integration uses synthetic model outputs."""

from dataclasses import dataclass
import json
from pathlib import Path
import sys
import types

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon import benchmark_worker as worker  # noqa: E402
from epistaeon.benchmark_inputs import (  # noqa: E402
    BenchmarkInputError, SCORERS, read_codon_alignment, residue_of, sha256_file,
)


@pytest.fixture
def frozen(tmp_path, monkeypatch):
    paths = {
        "alignment": str(tmp_path / "alignment.fa"),
        "weights": str(tmp_path / "model.safetensors"),
        "config": str(tmp_path / "config.json"),
    }
    Path(paths["alignment"]).write_text(
        ">focal\nCCCTGGAAG\n>other\nAGCTGCAAA\n>duplicate\nAGCTGCAAA\n"
    )
    Path(paths["weights"]).write_bytes(b"synthetic checkpoint; never loaded in this test")
    Path(paths["config"]).write_text("{}")
    monkeypatch.setattr(worker, "INPUT_PATHS", paths)
    monkeypatch.setattr(worker, "OUTPUT_ROOT", tmp_path / "outputs")
    protocol = {
        "schema_version": "1.0", "focal_taxon": "focal", "scorers": list(SCORERS),
        "seed": 0, "mode": "exact", "max_genotypes": 4096,
        "temperature": 0.5, "alpha": 1.0, "anchor": "root-token",
        "viability_floor": None, "device": "cpu", "threads": 1,
        "event_policy": "resolved_nonsynonymous", "fixed_background": "focal",
        **paths, "input_hashes": {key: sha256_file(path) for key, path in paths.items()},
    }
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(protocol))
    return protocol, path


@pytest.mark.parametrize("key, value", [
    ("reference", "/secret/truth.json"), ("sites", [1, 2]), ("mode", "sample"),
    ("seed", 2), ("device", "cuda"), ("threads", 2), ("temperature", 0.2),
    ("scorers", ["potts-metropolis"]), ("fixed_background", "plurality"),
    ("threads", 1.0), ("alpha", True),
])
def test_protocol_rejects_hints_and_nonfrozen_settings(frozen, key, value):
    protocol, _ = frozen
    protocol[key] = value
    with pytest.raises(BenchmarkInputError):
        worker.validate_protocol(protocol)


def test_hash_rejection_and_budget_failure_happen_before_model_import(frozen, monkeypatch):
    protocol, path = frozen
    monkeypatch.setattr(worker, "_configure_determinism", lambda: pytest.fail("model import reached"))
    Path(protocol["weights"]).write_bytes(b"changed")
    with pytest.raises(BenchmarkInputError, match="SHA-256"):
        worker.run_worker(path, SCORERS[0], worker.OUTPUT_ROOT / "bad")
    protocol["input_hashes"]["weights"] = sha256_file(protocol["weights"])
    protocol["max_genotypes"] = 2
    path.write_text(json.dumps(protocol))
    with pytest.raises(BenchmarkInputError, match="above max_genotypes"):
        worker.run_worker(path, SCORERS[0], worker.OUTPUT_ROOT / "bad")


def test_input_symlinks_and_escaping_outputs_fail_closed(frozen, tmp_path, monkeypatch):
    protocol, path = frozen
    original = Path(protocol["config"])
    moved = tmp_path / "actual-config.json"
    original.rename(moved)
    original.symlink_to(moved)
    with pytest.raises(BenchmarkInputError, match="regular input file"):
        worker.preflight_protocol(protocol)
    original.unlink()
    original.write_text("{}")
    monkeypatch.setattr(worker, "_configure_determinism", lambda: pytest.fail("model import reached"))
    with pytest.raises(BenchmarkInputError, match="subdirectory"):
        worker.run_worker(path, SCORERS[0], tmp_path / "outside")


def _fake_engine(monkeypatch, decreasing=False, no_couplings=False):
    """Use actual scorer classes/DP while replacing only the model boundary."""
    background = types.ModuleType("epistaeon.background")

    @dataclass
    class Endpoints:
        ancestral_codons: list
        derived_codons: dict
        ancestral_protein: str
        substitutions: list
        reference_name: str
        residue_number: dict
        plurality_support: dict
        classes: dict

        def class_counts(self):
            return {"nonsynonymous": len(self.substitutions)}

    calls = []

    def load_context(alignment, tree, **kwargs):
        calls.append(kwargs)
        rows = read_codon_alignment(alignment)
        return types.SimpleNamespace(
            rows=lambda: rows, derived_name=kwargs["derived_taxon"], ancestor_name=None,
            taxa=list(rows), n_taxa=len(rows), n_codons=len(next(iter(rows.values()))),
            distance_source="tn93",
        )

    class GenotypeRows:
        def __init__(self, ctx, ends, lattice, budget=None):
            self.ctx, self.ends, self.lattice, self.budget = ctx, ends, lattice, budget
            self.forward_passes, self.tn93_calls, self.real_row_hits = 4, 1, 1

        def prefetch(self, masks):
            assert list(masks) == list(range(1 << self.lattice.n_units))

        def embedding(self, mask):
            first, second = bool(mask & 1), bool(mask & 2)
            if decreasing:
                return np.array([4.0 - first - second, 0.0])
            return np.array([1 + first + .2 * second, .4 * second + .3 * first * second])

        def endpoint_frame_shift(self):
            return 0.0

        def background_frame_shift(self):
            return 0.0

        def check_matrix_patch_equals_full_recompute(self, mask):
            return 0.0

    class AppendedGenotypeRows(GenotypeRows):
        pass

    def forbidden(*args, **kwargs):
        pytest.fail("benchmark attempted unrestricted endpoint/site selection")

    background.Endpoints = Endpoints
    background.residue_of = residue_of
    background.GenotypeRows = GenotypeRows
    background.AppendedGenotypeRows = AppendedGenotypeRows
    background.load_context = load_context
    background.derive_endpoints = forbidden
    background.resolve_positions = forbidden
    background.root_token_embedding = lambda ctx: np.zeros(2)
    background.extant_embeddings = lambda ctx: np.array([[1, 0], [2, 1], [1.2, .5]])
    couplings = types.ModuleType("epistaeon.couplings")
    couplings.site_terms_from_lrt = lambda genotypes, lattice: {
        u: value for u, value in ((0, -1.), (1, 2.)) if u < lattice.n_units
    }
    couplings.couplings_from_coselection = lambda ctx, lattice: (
        {} if no_couplings or lattice.n_units < 2 else {(0, 1): 3.},
        {"unit_pairs_with_coupling": 0 if no_couplings or lattice.n_units < 2 else 1},
    )
    monkeypatch.setitem(sys.modules, "epistaeon.background", background)
    monkeypatch.setitem(sys.modules, "epistaeon.couplings", couplings)
    monkeypatch.setattr(worker, "_configure_determinism", lambda: None)
    return calls


@pytest.mark.parametrize("scorer", SCORERS)
def test_worker_runs_actual_scorer_and_exact_dp_with_frozen_endpoints(frozen, monkeypatch, scorer):
    protocol, path = frozen
    calls = _fake_engine(monkeypatch)
    outputs = []
    for suffix in ("first", "repeat"):
        output = worker.OUTPUT_ROOT / f"{scorer}-{suffix}"
        assert worker.run_worker(path, scorer, output) == 0
        selection = json.loads((output / "selection.json").read_text())
        assert selection["starting_codons"] == ["AGC", "TGC", "AAG"]
        assert selection["full_mask_reproduces_focal"] is True
        payload = json.loads((output / "order.json").read_text())
        assert payload["provenance"]["n_taxa"] == 3
        assert payload["provenance"]["mode"] == "exact"
        assert payload["provenance"]["diagnostic_baseline"] == (scorer == "dphi-product")
        assert len(payload["substitutions"]) == 2
        assert payload["accessible"]
        assert (output / "before_matrix.csv").is_file()
        assert (output / "position_matrix.csv").is_file()
        outputs.append(payload)
    # Output paths are absent from order payload; all numerical artifacts match.
    assert outputs[0] == outputs[1]
    assert all(call["prune_duplicates"] is False and call["max_species"] is None
               for call in calls)
    placement = outputs[0]["provenance"]["genotype_placement"]
    assert ("appended" in placement) == (scorer == "softmax-repaired")


def test_blocked_dphi_paths_are_reported_without_switching_method(frozen, monkeypatch):
    _, path = frozen
    _fake_engine(monkeypatch, decreasing=True)
    output = worker.OUTPUT_ROOT / "blocked"
    assert worker.run_worker(path, "dphi-product", output) == 0
    payload = json.loads((output / "order.json").read_text())
    assert payload["accessible"] is False
    assert payload["map_order"] == []
    assert payload["provenance"]["scorer"] == "dphi-product"


def test_uniform_potts_is_explicitly_diagnosed(frozen, monkeypatch):
    _, path = frozen
    _fake_engine(monkeypatch, no_couplings=True)
    output = worker.OUTPUT_ROOT / "uniform"
    assert worker.run_worker(path, "potts-metropolis", output) == 0
    payload = json.loads((output / "order.json").read_text())
    assert payload["discrimination"]["degenerate_uniform"] is True
    assert payload["before_matrix"]["matrix"][0][1] == pytest.approx(.5)


@pytest.mark.parametrize("scorer", SCORERS)
def test_one_unit_is_valid_for_every_existing_scorer(frozen, monkeypatch, scorer):
    protocol, path = frozen
    Path(protocol["alignment"]).write_text(">focal\nCCC\n>other\nAGC\n")
    protocol["input_hashes"]["alignment"] = sha256_file(protocol["alignment"])
    path.write_text(json.dumps(protocol))
    _fake_engine(monkeypatch)
    output = worker.OUTPUT_ROOT / f"one-unit-{scorer}"
    assert worker.run_worker(path, scorer, output) == 0
    selection = json.loads((output / "selection.json").read_text())
    assert selection["n_units"] == 1
    assert selection["n_genotypes"] == 2
    payload = json.loads((output / "order.json").read_text())
    assert payload["accessible"]
    assert payload["map_order"] == ["S1P"]
