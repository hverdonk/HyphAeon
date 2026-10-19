"""SLAC input contracts, node identity, raw state export and native integration."""

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from epistaeon import ancestral
from epistaeon.benchmark_inputs import read_codon_alignment, select_reconstructed

ROOT = Path(__file__).resolve().parents[1]


def dataset(tmp_path):
    scoring = tmp_path / "scoring.fa"
    scoring.write_text(">focus\nGATAGT\n>other\nGCTGGT\n")
    asr = tmp_path / "asr.fa"
    asr.write_text(scoring.read_text() + ">out1\nGCTGGT\n>out2\nGCTGGT\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(out1,out2,(focus,other)Anc)Root;")
    return scoring, asr, tree


def payload():
    return {"input": {"partition count": 1, "trees": {"0": "(out1,out2,(focus,other)Anc)Root;"}},
            "branch attributes": {"0": {"Anc": {"codon": [["GCT", "GGT"]],
                                                       "amino-acid": [["A", "G"]]}}}}


def test_explicit_tree_node_and_scoring_sequences_are_validated(tmp_path):
    alignment, asr, tree = dataset(tmp_path)
    _, rows, descendants = ancestral.validate_inputs(alignment, "focus", tree, "Anc", asr)
    assert descendants == ["focus", "other"]
    result = ancestral.extract_result(payload(), "Anc", "focus", rows, descendants)
    assert result["codons"] == ["GCT", "GGT"]
    selection = select_reconstructed(read_codon_alignment(alignment), "focus", result)
    assert selection["n_units"] == 2 and selection["n_genotypes"] == 4
    assert selection["starting_codons"] == ["GCT", "GGT"]


@pytest.mark.parametrize("newick,node,error", [
    ("(out1,out2,(focus,other)Anc)Root;", "focus", "internal"),
    ("(out1,out2,(focus,other)Anc)Root;", "missing", "internal"),
    ("((out1,out2)Anc,focus,other)Root;", "Anc", "descend"),
    ("((focus,other)Node,(out1,out2)Node)Root;", "Node", "exactly one"),
    ("((focus,other)Anc,(out1,out2)Outside)Root;", "Root", "degree-two"),
    ("(out1,out2,(focus,other,unexpected)Anc)Root;", "Anc", "match"),
])
def test_invalid_targets_fail_before_hyphy(tmp_path, newick, node, error):
    alignment, asr, tree = dataset(tmp_path)
    tree.write_text(newick)
    with pytest.raises(ancestral.AncestralError, match=error):
        ancestral.validate_inputs(alignment, "focus", tree, node, asr)


def test_asr_alignment_cannot_change_or_remove_scoring_rows(tmp_path):
    alignment, asr, tree = dataset(tmp_path)
    asr.write_text(asr.read_text().replace("GATAGT", "GCTAGT"))
    with pytest.raises(ancestral.AncestralError, match="identical"):
        ancestral.validate_inputs(alignment, "focus", tree, "Anc", asr)


def test_hyphy_output_must_preserve_target_clade():
    p = payload()
    p["input"]["trees"]["0"] = "(other,out2,(focus,out1)Anc)Root;"
    with pytest.raises(ancestral.AncestralError, match="descendant"):
        ancestral.extract_result(p, "Anc", "focus", {n: ["GCT"] * 2 for n in ["focus", "other", "out1", "out2"]}, ["focus", "other"])


def test_multistate_cells_remain_raw_and_export_as_one_unknown_site():
    p = payload()
    p["branch attributes"]["0"]["Anc"].update(codon=[["GWC", "GGT"]], **{"amino-acid": [["DV", "G"]]})
    result = ancestral.extract_result(p, "Anc", "focus", {n: ["GCT"] * 2 for n in ["focus", "other", "out1", "out2"]}, ["focus", "other"])
    assert result["protein"] == "XG"
    assert result["amino_acid_states"] == ["DV", "G"]
    result["node"] = "Anc"
    selected = select_reconstructed({"focus": ["GAT", "AGT"], "other": ["GCT", "GGT"]}, "focus", result)
    assert selected["columns"][0]["class"] == "unresolved"
    assert selected["starting_codons"][0] == "GAT"


def test_normalized_tip_names_keep_fitted_branch_lengths():
    p = payload()
    p["branch attributes"]["0"]["out1"] = {"original name": "out-1", "Global MG94xREV": 0.25}
    rows = {n: ["GCT"] * 2 for n in ["focus", "other", "out-1", "out2"]}
    result = ancestral.extract_result(p, "Anc", "focus", rows, ["focus", "other"])
    fitted = ancestral._tree(result["fitted_tree"])
    assert fitted.find_any(name="out-1").branch_length == 0.25


@pytest.mark.parametrize("version,code", [("HYPHY 2.5.101(MP)", 0), ("bad executable", 1)])
def test_incompatible_existing_hyphy_is_not_silently_replaced(monkeypatch, version, code):
    monkeypatch.setattr(ancestral.shutil, "which", lambda _: "/bin/hyphy")
    monkeypatch.setattr(ancestral.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, code, version, ""))
    with pytest.raises(ancestral.AncestralError, match="2.5.96"):
        ancestral.discover_hyphy()


def test_missing_hyphy_has_no_plurality_fallback(monkeypatch):
    monkeypatch.setattr(ancestral.shutil, "which", lambda _: None)
    with pytest.raises(ancestral.AncestralError, match="No plurality fallback"):
        ancestral.discover_hyphy()


def test_executable_override_uses_the_selected_installation(monkeypatch, tmp_path):
    executable = tmp_path / "hyphy-custom"
    seen = []
    monkeypatch.setattr(ancestral.shutil, "which", lambda name: seen.append(name) or str(executable))
    monkeypatch.setattr(ancestral.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, "HYPHY 2.5.96(MP)", ""))
    path, version = ancestral.discover_hyphy(executable)
    assert seen == [str(executable)] and path == str(executable)
    assert "2.5.96" in version


@pytest.mark.parametrize("exit_code,error", [(1, "exit status 1"), (0, "valid SLAC result")])
def test_slac_failures_preserve_logs_without_fabricating_an_ancestor(monkeypatch, tmp_path, exit_code, error):
    alignment, asr, tree = dataset(tmp_path)
    executable = tmp_path / "install/bin/hyphy"
    executable.parent.mkdir(parents=True)
    executable.write_text("fake executable")
    library = executable.parents[1] / "share/hyphy"
    library.mkdir(parents=True)
    (library / "genetic-code.def").write_text("supporting library input")
    monkeypatch.setattr(ancestral, "discover_hyphy", lambda _: (str(executable), "HYPHY 2.5.96"))
    def fail(command, **kwargs):
        kwargs["stderr"].write("SLAC diagnostic\n")
        assert f"LIBPATH={library}" in command
        assert "ENV=SetParameter(RANDOM_SEED,0,0);" in command
        return subprocess.CompletedProcess(command, exit_code)
    monkeypatch.setattr(ancestral.subprocess, "run", fail)
    output = tmp_path / "result"
    with pytest.raises(ancestral.AncestralError, match=error):
        ancestral.reconstruct(alignment, "focus", tree, "Anc", output, asr)
    assert not (output / "reconstruction.json").exists()
    assert (output / "stderr.log").read_text() == "SLAC diagnostic\n"
    assert json.loads((output / "status.json").read_text())["exit_code"] == exit_code
    assert "genetic-code.def" in json.loads((output / "protocol.json").read_text())["library_sha256"]


def test_native_cli_requires_asr_and_rejects_conflicting_ancestor():
    from epistaeon.cli import build_parser, cmd_order
    parser = build_parser()
    args = parser.parse_args(["order", "--alignment", "missing", "--focal-taxon", "focus"])
    with pytest.raises(SystemExit, match="--asr-tree"):
        cmd_order(args)
    args.ancestor, args.asr_node = "existing", "Anc"
    with pytest.raises(SystemExit, match="conflicts"):
        cmd_order(args)


def test_explicit_ancestor_bypasses_slac(monkeypatch):
    from epistaeon.cli import build_parser, cmd_order
    import epistaeon.background as background
    monkeypatch.setattr(ancestral, "reconstruct", lambda *a, **k: pytest.fail("SLAC called for explicit ancestor"))
    def stop(*a, **k):
        raise RuntimeError("explicit-ancestor path reached context loading")
    monkeypatch.setattr(background, "load_context", stop)
    args = build_parser().parse_args(["order", "--alignment", "a.fa", "--focal-taxon", "focus", "--ancestor", "ancestor"])
    with pytest.raises(RuntimeError, match="explicit-ancestor"):
        cmd_order(args)


@pytest.mark.skipif(not shutil.which("hyphy"), reason="HyPhy integration requires existing installation")
def test_real_sws1_slac_and_complete_budget_diagnostics(tmp_path):
    original = ROOT / "data/alignments/sws1/SWS1.codonified.fa"
    inputs = ROOT / "validation/sws1_slac/inputs"
    r = ancestral.reconstruct(original, "Homo_sapiens", inputs / "species_scaffold.nwk", "Boreoeutheria",
                              tmp_path / "asr", inputs / "SWS1.outgroups.codon.fa")
    s = select_reconstructed(read_codon_alignment(original), "Homo_sapiens", r)
    assert s["n_units"] == 28 and s["n_genotypes"] == 268435456
    assert not s["budget_ok"] and "at most 22" in s["budget_error"]
    assert s["full_mask_reproduces_focal"]
    # Empirical reference enters only after reconstruction and candidate selection.
    reference = json.loads((ROOT / "data/substitution_reference/sws1.json").read_text())
    from epistaeon.benchmark_evaluate import _reference_order
    truth = _reference_order(reference)
    keys = {(e["column0"], e["ancestral"], e["derived"]) for e in s["selected"]}
    assert set(truth["keys"]) <= keys
    assert len(truth["keys"]) == 7 and len(truth["pairs"]) == 18
