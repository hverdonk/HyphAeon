"""Exercise host freezing/isolation contracts without a model or Docker daemon."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon import benchmark_run as runner
from epistaeon.benchmark_artifacts import SCORERS, sha256_file, verify_seal, write_json
from epistaeon.benchmark_inputs import BenchmarkInputError, read_codon_alignment, read_substitution_panel, select_reconstructed, residue_of
from epistaeon.benchmark_worker import validate_protocol


IMAGE = "sha256:" + "a" * 64


def run_benchmark(alignment, focal_taxon, weights, output, max_genotypes=4096, **kwargs):
    tree = Path(alignment).parent / "asr_tree.nwk"
    tree.write_text("(focus,other1,other2)Anc;")
    return runner.run_benchmark(alignment, focal_taxon, weights, output, max_genotypes,
                                asr_tree=tree, asr_node="Anc", **kwargs)


def inputs(tmp_path):
    alignment = tmp_path / "alignment.fa"
    alignment.write_text(">focus\nGATAGTGCC\n>other1\nGCTGGTGCT\n>other2\nGCTGGTGCT\n")
    weights = tmp_path / "model.safetensors"
    weights.write_bytes(b"fake checkpoint: never loaded by host")
    write_json(tmp_path / "config.json", {"architecture": "test"})
    return alignment, weights


def mock_runtime(monkeypatch, failed=None):
    commands = []
    monkeypatch.setattr(runner, "require_docker", lambda: {"Server": {"Version": "test"}})

    def build(context, digest, log):
        Path(log).write_text("mock build\n")
        assert not any(Path(context).rglob("*.md"))
        return {"Id": IMAGE, "Os": "linux", "Architecture": "amd64"}

    monkeypatch.setattr(runner, "build_image", build)
    monkeypatch.setattr(runner, "check_isolation", lambda *args: {"status": "pass"})

    def run(command, **kwargs):
        if command[:3] == ["docker", "rm", "--force"]:
            return subprocess.CompletedProcess(command, 0, command[-1] + "\n", "")
        import csv
        specs = [command[i + 1] for i, item in enumerate(command) if item == "--mount"]
        output_spec = next(spec for spec in specs if "dst=/outputs" in spec)
        source = next(field[4:] for field in next(csv.reader([output_spec])) if field.startswith("src="))
        if "epistaeon.benchmark_prepare" in command:
            rows = read_codon_alignment(Path(source).parent / "inputs/alignment.fa")
            codons = rows["other1"]
            r = {"method": "hyphy-slac", "node": "Anc", "codons": codons,
                 "amino_acid_states": [residue_of(c) for c in codons], "descendant_taxa": list(rows),
                 "provenance": {"hyphy_version": "HYPHY 2.5.96(MP)"}}
            method = Path(source) / "asr"
            write_json(method / "reconstruction.json", r)
            protocol = json.loads((Path(source).parent / "inputs/protocol.json").read_text())
            panel = (read_substitution_panel(Path(source).parent / "inputs/substitution_panel.json",
                                            Path(source).parent / "inputs/alignment.fa", "focus")
                     if protocol["schema_version"] in {"3.0", "4.0"} else None)
            write_json(method / "selection.json", select_reconstructed(rows, "focus", r, protocol["max_genotypes"], panel=panel, fixed_background=protocol["fixed_background"]))
            return subprocess.CompletedProcess(command, 0)
        commands.append(command)
        scorer = command[command.index("--scorer") + 1]
        method = Path(source) / scorer
        method.mkdir()
        selection = json.loads((Path(source).parent / "selection.json").read_text())
        assert sorted(p.name for p in method.iterdir()) == []
        kwargs["stdout"].write("worker output\n")
        kwargs["stderr"].write("worker diagnostics\n")
        if scorer != failed:
            selected = selection["selected"]
            names = [entry["name"] for entry in selected]
            write_json(method / "order.json", {
                "substitutions": [dict(unit=e["unit"], name=e["name"], members=[e]) for e in selected],
                "accessible": True, "map_order": names,
                "before_matrix": {"names": names, "matrix": [[0 if i == j else .5 for j in range(len(names))] for i in range(len(names))]},
                "discrimination": {"degenerate_uniform": True},
            })
        return subprocess.CompletedProcess(command, 2 if scorer == failed else 0)

    monkeypatch.setattr(runner, "_run", run)
    return commands


def test_all_four_methods_frozen_and_sealed(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    commands = mock_runtime(monkeypatch)
    output = run_benchmark(alignment, "focus", weights, tmp_path / "run")
    seal = verify_seal(output)
    protocol = json.loads((output / "protocol.json").read_text())
    validate_protocol(protocol)
    assert protocol["input_hashes"]["alignment"] == sha256_file(alignment)
    assert protocol["runtime_image_id"] == IMAGE
    assert protocol == json.loads((output / "inputs/protocol.json").read_text())
    assert len(commands) == 4
    for scorer, command in zip(SCORERS, commands):
        method = output / "methods" / scorer
        status = json.loads((method / "status.json").read_text())
        assert status["status"] == "completed"
        assert status["command"] == command
        assert status["elapsed_seconds"] >= 0
        assert status["diagnostic_baseline"] == (scorer == "dphi-product")
        assert (method / "stdout.log").read_text() == "worker output\n"
        assert f"methods/{scorer}/order.json" in seal["files"]
    assert json.loads((output / "run_status.json").read_text())["status"] == "completed"


def test_failed_method_does_not_skip_other_methods(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    commands = mock_runtime(monkeypatch, failed="typicality-gate")
    output = run_benchmark(alignment, "focus", weights, tmp_path / "run")
    assert len(commands) == 4
    status = json.loads((output / "run_status.json").read_text())
    assert status["status"] == "completed_with_failures"
    assert status["failed_methods"] == ["typicality-gate"]
    verify_seal(output)


@pytest.mark.parametrize("background", ["focal", "slac"])
def test_panel_is_staged_hashed_and_frozen_before_all_scorers(tmp_path, monkeypatch, background):
    alignment, weights = inputs(tmp_path)
    panel_path = tmp_path / "panel.json"
    write_json(panel_path, {"schema_version": "1.0", "alignment_sha256": sha256_file(alignment),
                            "focal_taxon": "focus", "columns0": [0]})
    commands = mock_runtime(monkeypatch)
    output = run_benchmark(alignment, "focus", weights, tmp_path / "run", substitution_panel=panel_path, fixed_background=background)
    verify_seal(output)
    protocol = json.loads((output / "protocol.json").read_text())
    validate_protocol(protocol)
    assert protocol["schema_version"] == ("4.0" if background == "slac" else "3.0")
    assert protocol["input_hashes"]["substitution_panel"] == sha256_file(panel_path)
    assert len(commands) == 4
    selection = json.loads((output / "selection.json").read_text())
    assert selection["selected_columns0"] == [0]
    assert selection["unrestricted_nonsynonymous_units"] == 2
    assert selection["starting_codons"] == (["GCT", "GGT", "GCT"] if background == "slac" else ["GCT", "AGT", "GCC"])


def test_budget_failure_retains_complete_selection_and_runs_no_scorers(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    commands = mock_runtime(monkeypatch)
    with pytest.raises(runner.BenchmarkRunError, match="no sites are dropped"):
        run_benchmark(alignment, "focus", weights, tmp_path / "run", max_genotypes=2)
    assert commands == []
    output = tmp_path / "run"
    verify_seal(output)
    selection = json.loads((output / "selection.json").read_text())
    assert selection["n_units"] == 2 and not selection["budget_ok"]
    assert (output / "asr/reconstruction.json").is_file()
    assert json.loads((output / "run_status.json").read_text())["status"] == "preflight_failed"
    for scorer in SCORERS:
        assert not (output / "methods" / scorer / "order.json").exists()


def test_docker_unavailable_has_no_fallback_or_output(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    monkeypatch.setattr(runner.shutil, "which", lambda _: None)
    with pytest.raises(runner.BenchmarkRunError, match="no unisolated fallback"):
        run_benchmark(alignment, "focus", weights, tmp_path / "run")
    assert not (tmp_path / "run").exists()


def test_existing_output_is_never_reused(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    output = tmp_path / "run"
    output.mkdir()
    monkeypatch.setattr(runner, "require_docker", lambda: pytest.fail("Docker reached"))
    with pytest.raises(runner.BenchmarkRunError, match="already exists"):
        run_benchmark(alignment, "focus", weights, output)


def test_setup_failure_seals_partial_artifacts(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    mock_runtime(monkeypatch)

    def fail(*args):
        raise runner.BenchmarkRunError("deliberate image failure")

    monkeypatch.setattr(runner, "build_image", fail)
    with pytest.raises(runner.BenchmarkRunError, match="Partial artifacts are sealed"):
        run_benchmark(alignment, "focus", weights, tmp_path / "run")
    output = tmp_path / "run"
    verify_seal(output)
    assert json.loads((output / "run_status.json").read_text())["status"] == "setup_or_execution_failed"
    assert all(json.loads((output / "methods" / name / "status.json").read_text())["status"] == "failed"
               for name in SCORERS)


def _interrupt_first_worker(monkeypatch, cleanup_failure=None, already_absent=False):
    """Preserve real host orchestration while interrupting its Docker boundary."""
    mock_runtime(monkeypatch)
    ordinary_run = runner._run
    commands = []

    def interrupted(command, **kwargs):
        if "epistaeon.benchmark_prepare" in command:
            return ordinary_run(command, **kwargs)
        commands.append(command)
        if command[:2] == ["docker", "run"]:
            ordinary_run(command, **kwargs)
            raise KeyboardInterrupt()
        assert command[:3] == ["docker", "rm", "--force"]
        if isinstance(cleanup_failure, BaseException):
            raise cleanup_failure
        if cleanup_failure is not None:
            return subprocess.CompletedProcess(command, 1, "", cleanup_failure)
        if already_absent:
            return subprocess.CompletedProcess(
                command, 1, "", f"Error response from daemon: No such container: {command[-1]}"
            )
        return subprocess.CompletedProcess(command, 0, command[-1] + "\n", "")

    monkeypatch.setattr(runner, "_run", interrupted)
    return commands


@pytest.mark.parametrize("already_absent", [False, True])
def test_interruption_cleanup_preserves_logs_command_timing_and_seal(
    tmp_path, monkeypatch, already_absent,
):
    alignment, weights = inputs(tmp_path)
    commands = _interrupt_first_worker(monkeypatch, already_absent=already_absent)
    output = tmp_path / "run"
    with pytest.raises(runner.BenchmarkRunError, match="Partial artifacts are sealed"):
        run_benchmark(alignment, "focus", weights, output)
    assert len(commands) == 2  # interrupted first method, then forced cleanup
    command, cleanup = commands
    assert cleanup == ["docker", "rm", "--force", command[command.index("--name") + 1]]
    status = json.loads((output / "methods" / SCORERS[0] / "status.json").read_text())
    assert status["status"] == "failed"
    assert status["exit_code"] is None
    assert status["command"] == command
    assert status["started_at"] and status["finished_at"]
    assert status["elapsed_seconds"] >= 0
    assert "KeyboardInterrupt" in status["reason"]
    assert (output / "methods" / SCORERS[0] / "stdout.log").read_text() == "worker output\n"
    assert (output / "methods" / SCORERS[0] / "stderr.log").read_text() == "worker diagnostics\n"
    assert json.loads((output / "run_status.json").read_text())["status"] == "interrupted"
    verify_seal(output)


@pytest.mark.parametrize("cleanup_failure", [
    "Cannot connect to the Docker daemon",
    runner.BenchmarkRunError("docker executable became unavailable"),
    subprocess.TimeoutExpired(["docker", "rm"], 30),
    KeyboardInterrupt(),
    "No such container: some-other-container",
])
def test_unconfirmed_cleanup_preserves_failure_but_never_seals(
    tmp_path, monkeypatch, cleanup_failure,
):
    alignment, weights = inputs(tmp_path)
    commands = _interrupt_first_worker(monkeypatch, cleanup_failure=cleanup_failure)
    output = tmp_path / "run"
    with pytest.raises(runner.ContainerCleanupUnconfirmed, match="Artifacts are unsealed") as error:
        run_benchmark(alignment, "focus", weights, output)
    assert str(output) in str(error.value)
    assert not (output / "seal.json").exists()
    with pytest.raises(ValueError, match="Cannot read benchmark seal"):
        verify_seal(output)
    status = json.loads((output / "run_status.json").read_text())
    assert status["status"] == "cleanup_unconfirmed"
    assert status["seal_status"] == "withheld_cleanup_unconfirmed"
    assert status["container_cleanup"]["confirmed"] is False
    assert status["container_cleanup"]["container_name"] == commands[1][-1]
    assert status["container_cleanup"]["interrupted_by"] == "KeyboardInterrupt"
    method = output / "methods" / SCORERS[0]
    method_status = json.loads((method / "status.json").read_text())
    assert method_status["command"] == commands[0]
    assert method_status["exit_code"] is None
    assert "ContainerCleanupUnconfirmed" in method_status["reason"]
    assert (method / "stdout.log").read_text() == "worker output\n"
    assert (method / "stderr.log").read_text() == "worker diagnostics\n"


@pytest.mark.parametrize("already_absent", [False, True])
def test_nonzero_worker_exit_confirms_cleanup_and_preserves_normal_failure(
    tmp_path, monkeypatch, already_absent,
):
    alignment, weights = inputs(tmp_path)
    commands = mock_runtime(monkeypatch, failed=SCORERS[0])
    ordinary_run = runner._run
    cleanups = []

    def run(command, **kwargs):
        if command[:3] == ["docker", "rm", "--force"]:
            cleanups.append(command)
            if already_absent:
                return subprocess.CompletedProcess(
                    command, 1, "", f"Error response from daemon: No such container: {command[-1]}"
                )
            return subprocess.CompletedProcess(command, 0, command[-1] + "\n", "")
        return ordinary_run(command, **kwargs)

    monkeypatch.setattr(runner, "_run", run)
    output = run_benchmark(alignment, "focus", weights, tmp_path / "run")
    assert len(commands) == 4
    assert cleanups == [["docker", "rm", "--force", commands[0][commands[0].index("--name") + 1]]]
    method = output / "methods" / SCORERS[0]
    status = json.loads((method / "status.json").read_text())
    assert status["status"] == "failed"
    assert status["exit_code"] == 2
    assert status["command"] == commands[0]
    assert (method / "stderr.log").read_text() == "worker diagnostics\n"
    assert json.loads((output / "run_status.json").read_text())["status"] == "completed_with_failures"
    verify_seal(output)


def test_nonzero_docker_exit_with_unconfirmed_cleanup_remains_unsealed(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    commands = mock_runtime(monkeypatch, failed=SCORERS[0])
    ordinary_run = runner._run

    def run(command, **kwargs):
        if command[:3] == ["docker", "rm", "--force"]:
            return subprocess.CompletedProcess(command, 1, "", "Daemon connection failed")
        return ordinary_run(command, **kwargs)

    monkeypatch.setattr(runner, "_run", run)
    output = tmp_path / "run"
    with pytest.raises(runner.ContainerCleanupUnconfirmed, match="Artifacts are unsealed") as error:
        run_benchmark(alignment, "focus", weights, output)
    assert error.value.docker_exit_code == 2
    assert len(commands) == 1
    assert not (output / "seal.json").exists()
    status = json.loads((output / "run_status.json").read_text())
    assert status["status"] == "cleanup_unconfirmed"
    assert status["container_cleanup"]["interrupted_by"] == "Docker exit 2"
    assert status["container_cleanup"]["docker_exit_code"] == 2
    method = output / "methods" / SCORERS[0]
    method_status = json.loads((method / "status.json").read_text())
    assert method_status["exit_code"] == 2
    assert method_status["command"] == commands[0]
    assert (method / "stdout.log").read_text() == "worker output\n"
    with pytest.raises(ValueError, match="Cannot read benchmark seal"):
        verify_seal(output)


def test_isolation_probe_interruption_also_requires_confirmed_cleanup(tmp_path, monkeypatch):
    alignment, weights = inputs(tmp_path)
    actual_check_isolation = runner.check_isolation
    mock_runtime(monkeypatch)
    monkeypatch.setattr(runner, "check_isolation", actual_check_isolation)
    commands = []

    def interrupted(command, **kwargs):
        commands.append(command)
        if command[:2] == ["docker", "run"]:
            raise KeyboardInterrupt()
        assert command[:3] == ["docker", "rm", "--force"]
        return subprocess.CompletedProcess(command, 1, "", "Daemon unavailable")

    monkeypatch.setattr(runner, "_run", interrupted)
    output = tmp_path / "run"
    with pytest.raises(runner.ContainerCleanupUnconfirmed):
        run_benchmark(alignment, "focus", weights, output)
    assert len(commands) == 2
    assert "--scorer" not in commands[0]  # probe failed before a worker was started
    assert not (output / "seal.json").exists()
    assert json.loads((output / "run_status.json").read_text())["status"] == "cleanup_unconfirmed"


def test_allowlisted_snapshot_excludes_truth_docs_history_and_host_tools(tmp_path):
    source = tmp_path / "source"
    manifest = runner.snapshot_sources(source)
    files = set(manifest["files"])
    assert files == {p.relative_to(source).as_posix() for p in source.rglob("*") if p.is_file()}
    assert "epistaeon/benchmark_worker.py" in files
    for forbidden in ("benchmark_evaluate", "benchmark_artifacts", "benchmark_run", "reference",
                      "experimental", "README", "AGENTS", ".claude", "history", "test_"):
        assert not any(forbidden in filename for filename in files)
    context = tmp_path / "context"
    digest = runner.prepare_build_context(source, context)
    assert len(digest) == 64
    assert {p.name for p in context.iterdir()} == {"src", "Dockerfile", "requirements.lock", "runtime.json"}


def test_runtime_has_only_input_output_mounts_and_no_network(tmp_path):
    command = runner.docker_command(IMAGE, tmp_path / "inputs", tmp_path / "outputs", "python", "worker")
    assert command[command.index("--network") + 1] == "none"
    assert "--read-only" in command
    assert command[command.index("--cap-drop") + 1] == "ALL"
    assert command[command.index("--security-opt") + 1] == "no-new-privileges"
    assert command[command.index("--cpus") + 1] == "1"
    mounts = [command[i + 1] for i, item in enumerate(command) if item == "--mount"]
    assert len(mounts) == 2
    assert "dst=/inputs,readonly" in mounts[0]
    assert "dst=/outputs" in mounts[1]
    assert IMAGE in command
    assert not any("docker.sock" in value or "reference" in value for value in command)


@pytest.mark.parametrize("option", ["--reference", "--sites", "--blocks", "--ancestor", "--scorer"])
def test_runner_cli_rejects_truth_and_selection_options(option):
    from epistaeon.cli import build_parser
    with pytest.raises(SystemExit):
        build_parser().parse_args(["benchmark-run", "--alignment", "a.fa", "--focal-taxon", "focus",
                                  "--weights", "model.safetensors", "--output", "run", option, "anything"])
