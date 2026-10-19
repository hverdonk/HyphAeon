"""Host-side orchestration of an offline, reference-blind scorer comparison.

The host builds a content-addressed inference image from an explicit source
allowlist. Only four staged input files and an output directory enter workers.
There is deliberately no empirical-reference argument or native-run fallback.
"""

import csv
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from .benchmark_artifacts import (
    SCORERS, canonical_digest, seal_run, sha256_file, utc_now, write_json,
)
from .benchmark_inputs import BenchmarkInputError, select_benchmark_inputs


ROOT = Path(__file__).resolve().parents[3]
RUNTIME = ROOT / "epistaeon/benchmark_runtime"
SOURCE_ALLOWLIST = {
    "epistaeon": ("epistaeon/src/epistaeon", (
        "__init__.py", "lattice.py", "scoring.py", "trajectory.py", "timing.py",
        "background.py", "couplings.py", "io.py", "cli.py",
        "benchmark_inputs.py", "benchmark_worker.py",
    )),
    "aeon_core": ("aeon-core/src/aeon_core", (
        "__init__.py", "model.py", "weights.py", "inference.py", "splits.py",
        "dataset.py", "stats.py", "temporal.py", "io.py", "_progress.py",
    )),
    "hyphaeon": ("hyphaeon/src/hyphaeon", ("epistasis.py",)),
}
NAMESPACE_INIT = '"""Inference-only namespace; only the co-selection module is exposed."""\n'


class BenchmarkRunError(ValueError):
    pass


class ContainerCleanupUnconfirmed(BenchmarkRunError):
    """A container may still write artifacts, so those artifacts cannot be sealed."""

    def __init__(self, container_name, cleanup_detail, interrupted_by, output=None,
                 docker_exit_code=None):
        self.container_name = container_name
        self.cleanup_detail = cleanup_detail
        self.interrupted_by = interrupted_by
        self.docker_exit_code = docker_exit_code
        location = f" at {output}" if output is not None else ""
        super().__init__(
            f"Container cleanup unconfirmed for {container_name} after {interrupted_by}: "
            f"{cleanup_detail}. Artifacts are unsealed{location}; the container may still be running."
        )


def _run(command, **kwargs):
    try:
        return subprocess.run(command, **kwargs)
    except OSError as exc:
        raise BenchmarkRunError(f"Cannot execute {command[0]}: {exc}") from exc


def require_docker():
    if not shutil.which("docker"):
        raise BenchmarkRunError("Docker is required; no unisolated fallback is available")
    result = _run(["docker", "version", "--format", "{{json .}}"],
                  capture_output=True, text=True, timeout=30)
    if result.returncode:
        raise BenchmarkRunError(
            "Docker is not running or is inaccessible. Start Docker, then retry; "
            "no analysis was run and there is no unisolated fallback.\n" + result.stderr.strip()
        )
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise BenchmarkRunError("Docker returned invalid version metadata") from exc


def _regular_file(path, label):
    path = Path(path).expanduser()
    if path.is_symlink() or not path.is_file():
        raise BenchmarkRunError(f"{label} must be a regular file: {path}")
    return path.resolve()


def snapshot_sources(destination):
    """Copy only production modules; never copy an entire repository tree."""
    destination = Path(destination)
    files = {}
    for package, (directory, names) in SOURCE_ALLOWLIST.items():
        for name in names:
            origin = _regular_file(ROOT / directory / name, "Allowlisted source")
            target = destination / package / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(origin, target)
            files[target.relative_to(destination).as_posix()] = sha256_file(target)
    init = destination / "hyphaeon/__init__.py"
    init.write_text(NAMESPACE_INIT)
    files["hyphaeon/__init__.py"] = sha256_file(init)
    return {
        "schema_version": "1.0", "files": dict(sorted(files.items())),
        "source_snapshot_sha256": canonical_digest(files),
        "generated_files": {"hyphaeon/__init__.py": "Minimal namespace; no eager phenotype/disease imports"},
        "policy": "Explicit production Python allowlist. No tests, data, references, documentation or agent history.",
    }


def prepare_build_context(source_dir, destination, runtime_dir=RUNTIME):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("Dockerfile", "requirements.lock", "runtime.json"):
        shutil.copyfile(_regular_file(Path(runtime_dir) / name, "Runtime input"), destination / name)
    shutil.copytree(source_dir, destination / "src")
    entries = {p.relative_to(destination).as_posix(): sha256_file(p)
               for p in sorted(destination.rglob("*")) if p.is_file()}
    return canonical_digest(entries)


def build_image(context, context_digest, log_path):
    """Build the pinned image and return its immutable local image ID."""
    tag = f"epistaeon-benchmark:{context_digest[:24]}"
    with Path(log_path).open("w") as log:
        result = _run([
            "docker", "build", "--platform", "linux/amd64", "--progress=plain",
            "--label", f"org.epistaeon.build-context={context_digest}",
            "--tag", tag, str(context),
        ], stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise BenchmarkRunError(f"Pinned runtime image build failed; see {log_path}")
    result = _run(["docker", "image", "inspect", tag], capture_output=True, text=True, timeout=30)
    if result.returncode:
        raise BenchmarkRunError("Cannot inspect the built runtime image")
    metadata = json.loads(result.stdout)[0]
    if (metadata.get("Os") != "linux" or metadata.get("Architecture") != "amd64"
            or metadata.get("Config", {}).get("Labels", {}).get("org.epistaeon.build-context") != context_digest):
        raise BenchmarkRunError("Runtime image platform or source label does not match the frozen build")
    return metadata


def docker_command(image_id, inputs, outputs, *worker_command):
    """Construct the only runtime path; no repository, credentials or socket mount."""
    def mount(source, target, readonly=False):
        fields = ["type=bind", f"src={Path(source).resolve()}", f"dst={target}"]
        if readonly:
            fields.append("readonly")
        stream = io.StringIO()
        csv.writer(stream, lineterminator="\n").writerow(fields)
        return stream.getvalue().removesuffix("\n")

    # An identifiable container can be removed even if the Docker client is
    # interrupted before it has finished forwarding the stop signal.
    name = "epistaeon-" + canonical_digest({
        "image": image_id, "output": str(Path(outputs).resolve()),
        "worker_command": worker_command,
    })[:24]
    return [
        "docker", "run", "--rm", "--platform", "linux/amd64",
        "--name", name,
        "--network", "none", "--read-only", "--cap-drop", "ALL",
        "--security-opt", "no-new-privileges", "--pids-limit", "256",
        "--cpus", "1", "--memory", "4g",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--tmpfs", "/tmp:rw,nosuid,nodev,size=536870912",
        "--mount", mount(inputs, "/inputs", readonly=True),
        "--mount", mount(outputs, "/outputs"),
        image_id, *worker_command,
    ]


def _confirm_container_cleanup(name, interrupted_by, docker_exit_code=None):
    """A failed client is insufficient evidence that its container has stopped."""
    try:
        cleanup = _run(["docker", "rm", "--force", name], capture_output=True,
                       text=True, timeout=30)
    except BaseException as cleanup_error:
        raise ContainerCleanupUnconfirmed(
            name, f"{type(cleanup_error).__name__}: {cleanup_error}", interrupted_by,
            docker_exit_code=docker_exit_code,
        ) from cleanup_error
    detail = ((cleanup.stdout or "") + "\n" + (cleanup.stderr or "")).strip()
    already_absent = "no such container" in detail.lower() and name in detail
    if cleanup.returncode != 0 and not already_absent:
        raise ContainerCleanupUnconfirmed(
            name, detail or f"docker rm exited {cleanup.returncode}", interrupted_by,
            docker_exit_code=docker_exit_code,
        )


def _run_container(command, **kwargs):
    """Confirm termination after exceptions and every unsuccessful Docker exit."""
    name = command[command.index("--name") + 1]
    try:
        result = _run(command, **kwargs)
    except BaseException as run_error:
        _confirm_container_cleanup(name, type(run_error).__name__)
        raise
    if result.returncode != 0:
        _confirm_container_cleanup(name, f"Docker exit {result.returncode}", result.returncode)
    return result


ISOLATION_PROBE = r'''
import importlib.metadata, json, os, pathlib, shutil, socket
p = pathlib.Path('/inputs')
assert sorted(x.name for x in p.iterdir()) == ['alignment.fa', 'config.json', 'model.safetensors', 'protocol.json']
assert all(not x.is_symlink() and x.is_file() for x in p.iterdir())
for name in ('/app/AGENTS.md', '/app/src/epistaeon/AGENTS.md', '/app/data', '/app/epistaeon/data', '/workspace', '/var/run/docker.sock', '/app/src/epistaeon/benchmark_evaluate.py', '/app/src/epistaeon/benchmark_run.py'):
    assert not pathlib.Path(name).exists(), name
interfaces = {name: int(pathlib.Path('/sys/class/net', name, 'flags').read_text().strip(), 16) for _, name in socket.if_nameindex()}
active_interfaces = sorted(name for name, flags in interfaces.items() if flags & 1)
assert active_interfaces == ['lo'], interfaces
# Some Docker Desktop kernels create inactive default tunnel devices even in
# network=none namespaces. None may be up or have a route outside loopback.
assert not pathlib.Path('/proc/net/route').read_text().splitlines()[1:]
assert all(line.split()[-1] == 'lo' for line in pathlib.Path('/proc/net/ipv6_route').read_text().splitlines() if line.strip())
assert shutil.which('tn93') is None
assert importlib.metadata.version('tn93') == '1.2.2'
for name in ('/app/forbidden-write', '/inputs/forbidden-write'):
    try:
        pathlib.Path(name).write_text('must fail')
    except OSError:
        pass
    else:
        raise AssertionError('Unexpected write access: '+name)
print(json.dumps({'status':'pass','input_files':sorted(x.name for x in p.iterdir()),'network_interfaces':sorted(interfaces),'active_network_interfaces':active_interfaces,'interface_flags':interfaces,'rootfs_read_only':True,'inputs_read_only':True,'tn93':'1.2.2','packages':{x:importlib.metadata.version(x) for x in ['torch','numpy','scipy','biopython','tn93','safetensors']}}))
'''


def check_isolation(image_id, inputs, outputs):
    command = docker_command(image_id, inputs, outputs, "python", "-c", ISOLATION_PROBE)
    result = _run_container(command, capture_output=True, text=True)
    if result.returncode:
        raise BenchmarkRunError("Runtime isolation probe failed:\n" + result.stderr)
    return json.loads(result.stdout)


def _protocol(focal_taxon, hashes, max_genotypes, source_digest, image_id):
    return {
        "schema_version": "1.0", "focal_taxon": focal_taxon,
        "scorers": list(SCORERS), "seed": 0, "mode": "exact",
        "max_genotypes": max_genotypes, "temperature": 0.5, "alpha": 1.0,
        "anchor": "root-token", "viability_floor": None, "device": "cpu", "threads": 1,
        "event_policy": "resolved_nonsynonymous", "fixed_background": "focal",
        "alignment": "/inputs/alignment.fa", "weights": "/inputs/model.safetensors",
        "config": "/inputs/config.json", "input_hashes": hashes,
        "source_snapshot_sha256": source_digest, "runtime_image_id": image_id,
    }


def run_benchmark(alignment, focal_taxon, weights, output, max_genotypes=4096):
    alignment = _regular_file(alignment, "Alignment")
    weights = _regular_file(weights, "Checkpoint")
    if weights.suffix != ".safetensors":
        raise BenchmarkRunError("Benchmark checkpoints must use .safetensors")
    config = _regular_file(weights.parent / "config.json", "Checkpoint config.json")
    try:
        architecture = json.loads(config.read_text())
    except json.JSONDecodeError as exc:
        raise BenchmarkRunError("Checkpoint config.json is invalid JSON") from exc
    if not isinstance(architecture, dict) or not architecture:
        raise BenchmarkRunError("Checkpoint config.json must be a nonempty JSON object")
    selection = select_benchmark_inputs(alignment, focal_taxon, max_genotypes)
    output = Path(output).expanduser()
    if output.exists() or output.is_symlink():
        raise BenchmarkRunError("Output already exists; each frozen benchmark needs a new directory")
    docker_version = require_docker()  # fail before creating a partial run
    output = output.resolve()
    output.mkdir(parents=True)
    inputs, methods = output / "inputs", output / "methods"
    inputs.mkdir()
    methods.mkdir()
    started_at = utc_now()
    run_status = {"schema_version": "1.0", "started_at": started_at,
                  "status": "preparing", "scorers": list(SCORERS)}
    write_json(output / "run_status.json", run_status)
    failure = None
    try:
        for key, origin, target in (("alignment", alignment, "alignment.fa"),
                                    ("weights", weights, "model.safetensors"),
                                    ("config", config, "config.json")):
            shutil.copyfile(origin, inputs / target)
        hashes = {"alignment": sha256_file(inputs / "alignment.fa"),
                  "weights": sha256_file(inputs / "model.safetensors"),
                  "config": sha256_file(inputs / "config.json")}
        # Recheck the staged bytes, not just the original path before copying.
        staged_selection = select_benchmark_inputs(inputs / "alignment.fa", focal_taxon, max_genotypes)
        if staged_selection != selection:
            raise BenchmarkRunError("Alignment changed while staging the benchmark")
        write_json(output / "selection.json", selection)
        source_manifest = snapshot_sources(output / "source")
        write_json(output / "source_manifest.json", source_manifest)
        runtime_config = json.loads((RUNTIME / "runtime.json").read_text())
        write_json(output / "runtime_config.json", runtime_config)
        for name in ("Dockerfile", "requirements.lock", "requirements.in"):
            shutil.copyfile(RUNTIME / name, output / name)
        with tempfile.TemporaryDirectory(prefix="epistaeon-build-") as temporary:
            context = Path(temporary) / "context"
            digest = prepare_build_context(output / "source", context)
            image = build_image(context, digest, output / "build.log")
        write_json(output / "runtime_image.json", image)
        source_digest = source_manifest["source_snapshot_sha256"]
        protocol = _protocol(focal_taxon, hashes, max_genotypes, source_digest, image["Id"])
        write_json(output / "protocol.json", protocol)
        write_json(inputs / "protocol.json", protocol)
        write_json(output / "execution_manifest.json", {
            "schema_version": "1.0", "protocol_sha256": sha256_file(output / "protocol.json"),
            "source_snapshot_sha256": source_digest, "build_context_sha256": digest,
            "runtime_image_id": image["Id"], "runtime_config_sha256": sha256_file(output / "runtime_config.json"),
            "requirements_sha256": sha256_file(output / "requirements.lock"),
            "docker_version": docker_version, "input_hashes": hashes,
            "blinding": "Offline allowlisted inference image; alignment/checkpoint/config/protocol inputs only; no reference or evaluator access",
            "scorer_semantics": {
                "softmax-repaired": "Existing appended-row placement and remaining-unit alternatives retained",
                "dphi-product": "Diagnostic baseline; existing nonpositive-increment blocking retained",
            },
        })
        write_json(output / "isolation_check.json", check_isolation(image["Id"], inputs, methods))
        # No truth is available here. No scorer's outcome can change settings
        # or decide whether another scorer is run.
        for scorer in SCORERS:
            method = methods / scorer
            command = docker_command(
                image["Id"], inputs, methods,
                "python", "-m", "epistaeon.benchmark_worker", "--protocol", "/inputs/protocol.json",
                "--scorer", scorer, "--output", f"/outputs/{scorer}",
            )
            method_start, start = utc_now(), time.monotonic()
            print(f"[*] blinded benchmark: {scorer} ({selection['n_genotypes']} genotypes)", flush=True)
            result, method_error = None, None
            with tempfile.TemporaryDirectory(prefix="epistaeon-method-logs-") as temporary:
                stdout, stderr = Path(temporary) / "stdout.log", Path(temporary) / "stderr.log"
                try:
                    with stdout.open("w") as out, stderr.open("w") as err:
                        result = _run_container(command, stdout=out, stderr=err)
                except BaseException as exc:
                    method_error = exc
                finally:
                    method.mkdir(exist_ok=True)
                    for origin in (stdout, stderr):
                        if origin.exists():
                            shutil.copyfile(origin, method / origin.name)
                    complete = (result is not None and result.returncode == 0
                                and (method / "order.json").is_file())
                    reason = (None if complete else
                              (f"Worker interrupted: {type(method_error).__name__}: {method_error}"
                               if method_error is not None else
                               "Worker failed or did not produce order.json; inspect logs"))
                    write_json(method / "status.json", {
                        "schema_version": "1.0", "scorer": scorer,
                        "status": "completed" if complete else "failed",
                        "exit_code": (result.returncode if result is not None else
                                      getattr(method_error, "docker_exit_code", None)),
                        "started_at": method_start, "finished_at": utc_now(),
                        "elapsed_seconds": time.monotonic() - start, "command": command,
                        "diagnostic_baseline": scorer == "dphi-product", "reason": reason,
                    })
            if method_error is not None:
                raise method_error
            print(f"[{'ok' if complete else 'failed'}] {scorer}", flush=True)
        failed = [name for name in SCORERS
                  if json.loads((methods / name / "status.json").read_text())["status"] != "completed"]
        run_status.update(status="completed" if not failed else "completed_with_failures", failed_methods=failed)
    except (Exception, KeyboardInterrupt) as exc:
        failure = exc
        run_status.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "setup_or_execution_failed",
                          error=str(exc))
        if isinstance(exc, ContainerCleanupUnconfirmed):
            run_status.update(
                status="cleanup_unconfirmed", seal_status="withheld_cleanup_unconfirmed",
                container_cleanup={"confirmed": False, "container_name": exc.container_name,
                                   "detail": exc.cleanup_detail, "interrupted_by": exc.interrupted_by,
                                   "docker_exit_code": exc.docker_exit_code},
            )
        for scorer in SCORERS:
            status = methods / scorer / "status.json"
            if not status.exists():
                write_json(status, {"schema_version": "1.0", "scorer": scorer, "status": "failed",
                                    "exit_code": None, "reason": "Run did not complete: " + str(exc)})
    finally:
        run_status["finished_at"] = utc_now()
        write_json(output / "run_status.json", run_status)
        if isinstance(failure, ContainerCleanupUnconfirmed):
            print(f"[!] left {output} unsealed: container cleanup is unconfirmed", flush=True)
        else:
            seal = seal_run(output)
            print(f"[*] sealed {output} ({seal['content_digest']})", flush=True)
    if failure is not None:
        if isinstance(failure, ContainerCleanupUnconfirmed):
            raise ContainerCleanupUnconfirmed(
                failure.container_name, failure.cleanup_detail, failure.interrupted_by, output,
                docker_exit_code=failure.docker_exit_code,
            ) from failure
        raise BenchmarkRunError(f"Benchmark stopped: {failure}. Partial artifacts are sealed at {output}") from failure
    return output


def cmd_benchmark_run(args):
    try:
        output = run_benchmark(args.alignment, args.focal_taxon, args.weights,
                               args.output, args.max_genotypes)
    except (BenchmarkRunError, BenchmarkInputError, OSError, subprocess.TimeoutExpired) as exc:
        print(f"benchmark-run: {exc}", file=sys.stderr)
        return 2
    status = json.loads((output / "run_status.json").read_text())
    return 0 if status["status"] == "completed" else 1
