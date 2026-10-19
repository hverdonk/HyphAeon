"""Offline, single-method worker for a frozen blinded benchmark protocol.

This entry point has no reference/manuscript argument and never imports the
evaluator. Its only data inputs are the three fixed read-only input files.
"""

import argparse
import json
import os
from pathlib import Path
import random
import re

from .benchmark_inputs import (
    BenchmarkInputError, SCORERS, select_benchmark_inputs, sha256_file,
)


INPUT_PATHS = {
    "alignment": "/inputs/alignment.fa",
    "weights": "/inputs/model.safetensors",
    "config": "/inputs/config.json",
}
OUTPUT_ROOT = Path("/outputs")
_REQUIRED = {
    "schema_version", "focal_taxon", "scorers", "seed", "mode",
    "max_genotypes", "temperature", "alpha", "anchor", "viability_floor",
    "device", "threads", "event_policy", "fixed_background", "alignment",
    "weights", "config", "input_hashes",
}
_OPTIONAL = {"source_snapshot_sha256", "runtime_image_id"}


def validate_protocol(protocol):
    """Reject unknown fields and any departure from the benchmark profile."""
    if not isinstance(protocol, dict):
        raise BenchmarkInputError("protocol must be a JSON object")
    missing = _REQUIRED - protocol.keys()
    extra = protocol.keys() - (_REQUIRED | _OPTIONAL)
    if missing or extra:
        raise BenchmarkInputError(
            f"invalid protocol fields: missing={sorted(missing)}, unknown={sorted(extra)}"
        )
    fixed = {
        "schema_version": "1.0", "scorers": list(SCORERS), "seed": 0,
        "mode": "exact", "temperature": 0.5, "alpha": 1.0,
        "anchor": "root-token", "viability_floor": None, "device": "cpu",
        "threads": 1, "event_policy": "resolved_nonsynonymous",
        "fixed_background": "focal", **INPUT_PATHS,
    }
    for key, expected in fixed.items():
        value = protocol[key]
        wrong_type = ((type(expected) is int and type(value) is not int)
                      or (type(expected) is float and type(value) not in (int, float)))
        if value != expected or wrong_type:
            raise BenchmarkInputError(f"protocol {key!r} must be {expected!r}")
    if not isinstance(protocol["focal_taxon"], str) or not protocol["focal_taxon"].strip():
        raise BenchmarkInputError("focal_taxon must be a nonempty exact row name")
    if type(protocol["max_genotypes"]) is not int or protocol["max_genotypes"] < 1:
        raise BenchmarkInputError("max_genotypes must be a positive integer")
    hashes = protocol["input_hashes"]
    if not isinstance(hashes, dict) or set(hashes) != set(INPUT_PATHS):
        raise BenchmarkInputError("input_hashes must contain exactly alignment, weights, config")
    for key, digest in hashes.items():
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise BenchmarkInputError(f"invalid SHA-256 for {key}")
    if "source_snapshot_sha256" in protocol:
        if not re.fullmatch(r"[0-9a-f]{64}", str(protocol["source_snapshot_sha256"])):
            raise BenchmarkInputError("invalid source_snapshot_sha256")
    if "runtime_image_id" in protocol:
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", str(protocol["runtime_image_id"])):
            raise BenchmarkInputError("invalid runtime_image_id")
    return protocol


def preflight_protocol(protocol):
    """Check actual bytes and the exact-inference budget before model imports."""
    validate_protocol(protocol)
    for key, fixed_path in INPUT_PATHS.items():
        path = Path(fixed_path)
        if path.is_symlink() or not path.is_file():
            raise BenchmarkInputError(f"{key} must be a regular input file: {path}")
        if sha256_file(path) != protocol["input_hashes"][key]:
            raise BenchmarkInputError(f"{key} input SHA-256 does not match frozen protocol")
    return select_benchmark_inputs(
        protocol["alignment"], protocol["focal_taxon"], protocol["max_genotypes"]
    )


def _configure_determinism():
    # Environment is installed before numerical/model imports. Docker also sets
    # PYTHONHASHSEED before interpreter startup; setting it here covers children.
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = "1"
    os.environ["PYTHONHASHSEED"] = "0"
    random.seed(0)
    import numpy as np
    import torch
    np.random.seed(0)
    torch.manual_seed(0)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)


def order_namespace(protocol, scorer, output):
    """Construct fixed engine arguments without accepting site/order hints."""
    validate_protocol(protocol)
    if scorer not in SCORERS:
        raise BenchmarkInputError(f"unknown benchmark scorer {scorer!r}")
    return argparse.Namespace(
        alignment=protocol["alignment"], focal_taxon=protocol["focal_taxon"],
        weights=protocol["weights"], output=str(output), scorer=scorer,
        tree=None, ancestor=None, variant=None, cpu=True, max_species=None,
        use_tn93=True, anchor="root-token", mode="exact", sites=None,
        blocks=None, alpha=1.0, t_sel=0.5, viability_floor=None,
        max_genotypes=protocol["max_genotypes"], seed=0, n_samples=20000,
        focal_unit=None, permissive_units=None, ne=None, s=None,
        generation_time=None, _benchmark_profile=True,
    )


def run_worker(protocol_path, scorer, output):
    protocol = json.loads(Path(protocol_path).read_text(encoding="utf-8"))
    selection = preflight_protocol(protocol)
    output = Path(output)
    root = OUTPUT_ROOT.resolve()
    resolved = output.resolve()
    if not output.is_absolute() or resolved == root or root not in resolved.parents:
        raise BenchmarkInputError(f"worker output must be a subdirectory of {root}")
    if output.exists() and any(output.iterdir()):
        raise BenchmarkInputError("worker output directory must be empty")
    output.mkdir(parents=True, exist_ok=True)
    (output / "selection.json").write_text(
        json.dumps(selection, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _configure_determinism()
    from .cli import cmd_order
    return cmd_order(order_namespace(protocol, scorer, output))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--scorer", required=True, choices=SCORERS)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        return run_worker(args.protocol, args.scorer, args.output)
    except (BenchmarkInputError, OSError, json.JSONDecodeError) as exc:
        parser.exit(2, f"benchmark worker: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
