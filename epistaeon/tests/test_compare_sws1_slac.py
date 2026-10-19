"""Coordinates and unavailable states in the evaluation-only SLAC comparison."""

import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "epistaeon/scripts/compare_sws1_slac.py"
spec = importlib.util.spec_from_file_location("compare_sws1_slac", SCRIPT)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def _fixture(tmp_path, prediction="M-FFTFTASV"):
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    (run / "results").mkdir()
    codon_for = {aa: next(codon for codon, residue in comparison.CODE.items() if residue == aa)
                 for aa in comparison.AA}

    def dna(protein):
        return "".join("---" if residue == "-" else "NNN" if residue == "X" else codon_for[residue]
                       for residue in protein)

    alignment = run / "inputs/alignment.fa"
    alignment.write_text(f">Homo_sapiens\n{dna('M-TLFLPGTV')}\n>other\n{dna('M-TLFLPGTV')}\n")
    codon = run / "results/ancestors.codon.fa"
    codon.write_text(f">Boreoeutheria\n{dna(prediction)}\n")
    protein = run / "results/ancestors.protein.fa"
    protein.write_text(f">Boreoeutheria\n{prediction}\n")
    node = run / "results/selected_node.json"
    node.write_text(json.dumps({"node_id": "Boreoeutheria"}))
    seal = run / "results/inference_seal.json"
    frozen = [alignment, codon, protein, node]
    seal.write_text(json.dumps({"files": {path.relative_to(run).as_posix(): comparison.sha256(path)
                                          for path in frozen}}))
    paper = tmp_path / "paper_reference.json"
    paper.write_text(json.dumps({
        "ancestor_protein": "AFFTFTASG", "paper_human_protein": "ATLFLPGTG",
        "terminal_mask": {"comparable_core_ranges_1based": [[2, 8]],
                          "engineered_ranges_1based": [[1, 1], [9, 9]]},
    }))
    empirical = tmp_path / "events.json"
    empirical.write_text(json.dumps({"events": [
        {"id": f"e{i}", "label": f"{ancestral}{i}{derived}", "position": i,
         "from": ancestral, "to": derived,
         "structure_mapping": [{"human_cds_position_1based": i + 1}]}
        for i, (ancestral, derived) in enumerate(zip("FFTFTAS", "TLFLPGT"), 1)
    ]}))
    return {"alignment": alignment, "ancestors_codon": codon, "ancestors_protein": protein,
            "selected_node": node, "paper_reference": paper, "substitution_reference": empirical,
            "inference_seal": seal, "run_dir": run, "output": run / "comparison"}


def test_real_s4_human_bridge_uses_core_and_retains_input_gaps():
    paper = comparison.read_json(ROOT / "epistaeon/validation/sws1_slac_reference/source_reference.json")
    rows = comparison.read_fasta(ROOT / "epistaeon/data/alignments/sws1/SWS1.codonified.fa")
    human = comparison.translate(rows["Homo_sapiens"])
    core_range, core, engineered = comparison.mask_positions(paper)
    bridge, inverse, metadata = comparison.build_human_bridge(paper["paper_human_protein"], core_range, human)
    assert len(core) == 282
    assert len(engineered) == 65
    assert metadata["human_core_range_1based"] == [28, 309]
    assert metadata["alignment_core_mapped_column_range_1based"] == [30, 311]
    for event in paper["focal_sites"]:
        position = event["paper_protein_position_1based"]
        assert bridge[position]["alignment_column_1based"] == event["benchmark_alignment_codon_column_1based"]
        assert inverse[event["historical_human_cds_position_1based"]] == position


def test_bridge_handles_alignment_insertions_without_a_constant_offset():
    bridge, _, metadata = comparison.build_human_bridge("XACDFY", [2, 5], "MAC-D-FV")
    assert [bridge[position]["alignment_column_1based"] for position in range(2, 6)] == [2, 3, 5, 7]
    assert metadata["human_core_range_1based"] == [2, 5]


def test_reference_insertion_has_no_invented_human_coordinate():
    bridge, _, metadata = comparison.build_human_bridge("XAC-DFY", [2, 6], "MACDFV")
    assert 4 not in bridge
    assert metadata["unmapped_human_gap_positions"] == [4]


@pytest.mark.parametrize("human", ["MACDFACDFV", "MACGEV"])
def test_nonunique_or_inexact_human_bridge_is_rejected(human):
    with pytest.raises(ValueError, match="exactly once"):
        comparison.build_human_bridge("XACDFY", [2, 5], human)


def test_mask_must_partition_figure_coordinates():
    paper = {"ancestor_protein": "ACDF", "paper_human_protein": "ACDF",
             "terminal_mask": {"comparable_core_ranges_1based": [[2, 3]],
                               "engineered_ranges_1based": [[1, 2], [4, 4]]}}
    with pytest.raises(ValueError, match="partition"):
        comparison.mask_positions(paper)


def test_ambiguous_translation_resolves_only_unique_amino_acids():
    assert comparison.translate_codon("ATH") == "I"
    assert comparison.translate_codon("GCN") == "A"
    assert comparison.translate_codon("NNN") == "X"
    assert comparison.translate_codon("---") == "-"
    with pytest.raises(ValueError, match="partial-gap"):
        comparison.translate_codon("A--")


def test_complete_comparison_is_reproducible_and_leaves_frozen_inputs_untouched(tmp_path):
    args = _fixture(tmp_path)
    before = {path: comparison.sha256(path) for path in args.values() if path.is_file()}
    first = comparison.compare(**args)
    assert first["focal_sites"]["matched_residue_count"] == 7
    assert first["paper_natural_core"]["identity_among_resolved"] == 1
    assert first["mask"]["excluded_engineered_positions"] == 2
    assert first["coordinate_bridge"]["alignment_core_mapped_column_range_1based"] == [3, 9]
    assert before == {path: comparison.sha256(path) for path in before}
    original = (args["output"] / "comparison.json").read_bytes()
    args["output"] = tmp_path / "repeat"
    second = comparison.compare(**args)
    assert first == second
    assert original == (args["output"] / "comparison.json").read_bytes()


@pytest.mark.parametrize("unavailable", ["X", "-"])
def test_gaps_and_ambiguity_are_unavailable_not_resolved_substitutions(tmp_path, unavailable):
    args = _fixture(tmp_path, prediction=f"M-L{unavailable}TFTASV")
    report = comparison.compare(**args)
    core = report["paper_natural_core"]
    assert core["expected_residue_count"] == 7
    assert core["resolved_reconstruction_count"] == 6
    assert core["matched_residue_count"] == 5
    assert core["substitution_disagreement_count"] == 1
    assert core["unresolved_or_unmapped_count"] == 1
    assert core["identity_among_resolved"] == pytest.approx(5 / 6)
    assert core["strict_match_fraction_all_expected_residues"] == pytest.approx(5 / 7)
    assert core["unresolved_or_unmapped"][0]["reconstructed_ancestor"] == unavailable


def test_changed_predictions_fail_before_comparison_output(tmp_path):
    args = _fixture(tmp_path)
    args["ancestors_protein"].write_text(">Boreoeutheria\nM-TLFLPGTV\n")
    with pytest.raises(ValueError, match="checksum"):
        comparison.compare(**args)
    assert not args["output"].exists()


def test_unsealed_comparison_inputs_are_rejected(tmp_path):
    args = _fixture(tmp_path)
    seal = comparison.read_json(args["inference_seal"])
    del seal["files"]["results/ancestors.protein.fa"]
    args["inference_seal"].write_text(json.dumps(seal))
    with pytest.raises(ValueError, match="not frozen"):
        comparison.compare(**args)
    assert not args["output"].exists()


def test_no_numerical_or_model_imports():
    assert "torch" not in comparison.__dict__
    assert "numpy" not in comparison.__dict__
    assert "epistaeon" not in comparison.__dict__
