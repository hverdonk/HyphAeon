"""Manuscript site scope is declared independently of empirical order truth."""

import json
from pathlib import Path

import pytest

from epistaeon.benchmark_inputs import (
    BenchmarkInputError, read_substitution_panel, select_reconstructed,
    sha256_file, validate_substitution_panel,
)


def panel(columns, focal="focus", digest="a" * 64):
    return {"schema_version": "1.0", "alignment_sha256": digest,
            "focal_taxon": focal, "columns0": columns}


def test_panel_filters_before_budget_and_freezes_excluded_nonsynonymous_codons():
    rows = {"focus": ["GAT"] * 28, "other": ["GCT"] * 28}
    reconstruction = {"method": "hyphy-slac", "node": "Anc", "codons": ["GCT"] * 28,
                      "amino_acid_states": ["A"] * 28}
    columns = [1, 4, 10, 12, 14, 24, 27]
    selected = select_reconstructed(rows, "focus", reconstruction, 128, panel(columns))
    assert selected["schema_version"] == "3.0"
    assert selected["unrestricted_nonsynonymous_units"] == 28
    assert selected["selected_columns0"] == columns
    assert selected["n_units"] == 7 and selected["n_genotypes"] == 128 and selected["budget_ok"]
    assert selected["full_mask_reproduces_focal"]
    for i, codon in enumerate(selected["starting_codons"]):
        assert codon == ("GCT" if i in columns else "GAT")
    assert sum(e["exclusion_reason"] == "outside_substitution_panel" for e in selected["columns"]) == 21
    failed = select_reconstructed(rows, "focus", reconstruction, 64, panel(columns))
    assert not failed["budget_ok"] and failed["n_units"] == 7  # never truncate the panel to fit


def test_panel_does_not_supply_or_repair_ancestral_residues():
    rows = {"focus": ["GAT", "GGT", "GAT", "AAG"], "other": ["GCT"] * 4}
    reconstruction = {"method": "hyphy-slac", "node": "Anc", "codons": ["GCT", "GGT", "NNN", "AAA"],
                      "amino_acid_states": ["A", "G", "X", "K"]}
    selected = select_reconstructed(rows, "focus", reconstruction, panel=panel([0, 1, 2, 3]))
    assert selected["selected_columns0"] == [0]
    assert selected["unselected_panel_columns0"] == [1, 2, 3]
    assert [e["exclusion_reason"] for e in selected["columns"]] == [None, "unchanged", "unresolved", "synonymous"]
    assert selected["starting_codons"] == ["GCT", "GGT", "GAT", "AAG"]


@pytest.mark.parametrize("change", [
    {"expected_order": [1, 0]}, {"events": []}, {"columns0": [2, 1]},
    {"columns0": [0, 0]}, {"columns0": [True]}, {"columns0": [-1]},
    {"columns0": [3]}, {"columns0": []}, {"focal_taxon": "other"},
    {"alignment_sha256": "wrong"},
])
def test_panel_rejects_truth_fields_and_invalid_scope(change):
    with pytest.raises(BenchmarkInputError):
        validate_substitution_panel(dict(panel([0]), **change), 3, "focus")


def test_panel_alignment_binding_and_duplicate_keys(tmp_path):
    alignment = tmp_path / "alignment.fa"
    alignment.write_text(">focus\nGAT\n>other\nGCT\n")
    path = tmp_path / "panel.json"
    path.write_text(json.dumps(panel([0], digest=sha256_file(alignment))))
    assert read_substitution_panel(path, alignment, "focus")["columns0"] == [0]
    alignment.write_text(alignment.read_text().replace("GCT", "GGT"))
    with pytest.raises(BenchmarkInputError, match="SHA-256 mismatch"):
        read_substitution_panel(path, alignment, "focus")
    path.write_text('{"columns0": [0], "columns0": [1]}')
    with pytest.raises(BenchmarkInputError, match="Duplicate"):
        read_substitution_panel(path, alignment, "focus")


def test_sws1_panel_matches_reviewed_source_crosswalk_without_order_information():
    root = Path(__file__).resolve().parents[1]
    alignment = root / "data/alignments/sws1/SWS1.codonified.fa"
    supplied = read_substitution_panel(alignment.parent / "substitution_panel.json", alignment, "Homo_sapiens")
    reference = json.loads((root / "data/substitution_reference/sws1.json").read_text())
    columns = sorted({m["alignment_reference_columns_1based"]["SWS1"] - 1
                      for e in reference["events"] if e.get("order_scoring", True)
                      for m in e["structure_mapping"]})
    assert supplied["columns0"] == columns == [44, 47, 50, 84, 91, 112, 116]
    assert set(supplied) == {"schema_version", "alignment_sha256", "focal_taxon", "columns0"}
    assert [f"c{col + 1}" for col in columns] == (alignment.parent / "sites_columns.txt").read_text().split()
