"""Truth-free input selection and endpoint invariants for the benchmark."""

import os
from pathlib import Path
import subprocess
import sys

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from epistaeon.benchmark_inputs import (  # noqa: E402
    BenchmarkInputError, read_codon_alignment, select_benchmark_inputs, select_from_rows,
)


def test_preflight_does_not_import_numerical_or_model_libraries():
    env = dict(os.environ, PYTHONPATH=str(SRC))
    completed = subprocess.run(
        [sys.executable, "-c", "import sys; import epistaeon.benchmark_inputs; "
         "assert 'numpy' not in sys.modules; assert 'torch' not in sys.modules"],
        env=env, capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr


def test_selects_all_nonsynonymous_and_freezes_every_other_codon():
    # Two AA changes, one synonymous change, one gap, and one ambiguous focal.
    rows = {
        "focal": ["CCC", "TGG", "AAG", "---", "NNN"],
        "other": ["AGC", "TGC", "AAA", "GGC", "AAA"],
    }
    selection = select_from_rows(rows, "focal")
    assert selection["selected_columns0"] == [0, 1]
    assert selection["n_units"] == 2
    assert selection["n_genotypes"] == 4
    assert selection["starting_codons"] == ["AGC", "TGC", "AAG", "---", "NNN"]
    built = selection["starting_codons"].copy()
    for event in selection["selected"]:
        built[event["column0"]] = event["derived_codon"]
    assert built == rows["focal"]
    assert selection["classification_counts"] == {
        "indel": 1, "nonsynonymous": 2, "synonymous": 1, "unresolved": 1,
    }


def test_plurality_residue_ties_and_first_seen_codon_ties():
    rows = {
        "focal": ["CCC", "CCC", "TGG"],
        "first": ["---", "GAC", "AGT"],
        "second": ["AGC", "GCT", "AGC"],
    }
    selection = select_from_rows(rows, "focal")
    assert selection["plurality_codons"] == ["AGC", "GCT", "AGT"]
    assert selection["selected"][2]["plurality_support"] == 1.0
    # A tied S codon follows input order, not lexical codon sorting.
    reordered = {name: rows[name] for name in ("focal", "second", "first")}
    assert select_from_rows(reordered, "focal")["plurality_codons"][2] == "AGC"


def test_all_taxa_and_duplicate_sequences_contribute_to_plurality():
    rows = {"focal": ["TGG", "CCC"], "a": ["TGC", "AGC"],
            "duplicate": ["TGC", "AGC"], "b": ["TGG", "CCC"]}
    selection = select_from_rows(rows, "focal")
    assert selection["taxa"] == list(rows)
    assert selection["n_taxa"] == 4
    assert selection["selected"][0]["plurality_support"] == pytest.approx(2 / 3)


def test_unknown_plurality_is_excluded_and_not_called_a_deletion():
    rows = {"focal": ["CCC", "TGG", "AAA"], "other": ["AGC", "TGC", "NNN"]}
    selection = select_from_rows(rows, "focal")
    assert selection["columns"][2]["class"] == "unresolved"
    assert selection["starting_codons"][2] == "AAA"


def test_one_unit_is_a_valid_exact_benchmark():
    selection = select_from_rows({"focal": ["CCC"], "other": ["AGC"]}, "focal")
    assert selection["n_units"] == 1
    assert selection["n_genotypes"] == 2
    assert selection["full_mask_reproduces_focal"]


@pytest.mark.parametrize("rows, message", [
    ({"x": ["AAA"]}, "two alignment taxa"),
    ({"x": ["AAA"], "y": ["AAA"]}, "at least one resolved"),
    ({"x": ["AAA", "CCC"], "y": ["AAA"]}, "equal nonzero"),
])
def test_input_shape_and_minimum_units_fail_closed(rows, message):
    with pytest.raises(BenchmarkInputError, match=message):
        select_from_rows(rows, "x")


def test_budget_failure_drops_no_sites_and_never_samples():
    rows = {"x": ["AAA"] * 13, "y": ["CCC"] * 13}
    with pytest.raises(BenchmarkInputError, match="8192.*no sites are dropped"):
        select_from_rows(rows, "x")
    with pytest.raises(BenchmarkInputError, match="at most 22"):
        select_from_rows({"x": ["AAA"] * 23, "y": ["CCC"] * 23}, "x", 1 << 25)


@pytest.mark.parametrize("contents, message", [
    (">a\nAAA\n>a\nCCC\n", "duplicate FASTA"),
    (">a\nAA\n>b\nCC\n", "divisible by three"),
    (">a\nAAA\n>b\nCCCCCC\n", "same nonzero length"),
    (">a\nXXX\n>b\nCCC\n", "invalid nucleotide"),
])
def test_malformed_fasta_is_not_trimmed_or_padded(tmp_path, contents, message):
    path = tmp_path / "input.fa"
    path.write_text(contents)
    with pytest.raises(BenchmarkInputError, match=message):
        read_codon_alignment(path)


def test_focal_name_is_exact(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text(">Homo_sapiens\nAAACCC\n>other\nCCCAAA\n")
    with pytest.raises(BenchmarkInputError, match="exact FASTA row"):
        select_benchmark_inputs(path, "homo")


def test_sws1_selection_is_derived_only_from_extant_alignment():
    path = Path(__file__).resolve().parents[1] / "data/alignments/sws1/SWS1.codonified.fa"
    selection = select_benchmark_inputs(path, "Homo_sapiens")
    assert selection["n_taxa"] == 15
    assert selection["n_units"] == 10
    assert selection["n_genotypes"] == 1024
    assert selection["full_mask_reproduces_focal"]
