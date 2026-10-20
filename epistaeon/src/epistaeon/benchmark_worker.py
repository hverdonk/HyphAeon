"""Offline, single-method worker for a frozen blinded benchmark protocol.

This entry point has no reference/manuscript argument and never imports the
evaluator. Its data inputs are fixed read-only files declared by the protocol.
"""

import argparse
import json
import os
from pathlib import Path
import random
import re

from .benchmark_inputs import (
    BenchmarkInputError, SCORERS, select_benchmark_inputs, sha256_file,
    read_codon_alignment, select_reconstructed, read_substitution_panel,
)


INPUT_PATHS = {
    "alignment": "/inputs/alignment.fa",
    "weights": "/inputs/model.safetensors",
    "config": "/inputs/config.json",
}
ASR_PATHS = {"asr_alignment": "/inputs/asr_alignment.fa", "asr_tree": "/inputs/asr_tree.nwk"}
PANEL_PATH = "/inputs/substitution_panel.json"
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
    v4 = protocol.get("schema_version") == "4.0"
    v3 = protocol.get("schema_version") in {"3.0", "4.0"}
    v2 = protocol.get("schema_version") in {"2.0", "3.0", "4.0"}
    asr_keys = {"ancestral_method", "asr_alignment", "asr_tree", "asr_node"} if v2 else set()
    if v3:
        asr_keys.add("substitution_panel")
    if v4:
        asr_keys.add("endpoint_policy")
    missing = (_REQUIRED | asr_keys) - protocol.keys()
    extra = protocol.keys() - (_REQUIRED | asr_keys | _OPTIONAL | ({"reconstruction_sha256"} if v2 else set()))
    if missing or extra:
        raise BenchmarkInputError(
            f"invalid protocol fields: missing={sorted(missing)}, unknown={sorted(extra)}"
        )
    fixed = {
        "schema_version": protocol["schema_version"] if v2 else "1.0", "scorers": list(SCORERS), "seed": 0,
        "mode": "exact", "temperature": 0.5, "alpha": 1.0,
        "anchor": "root-token", "viability_floor": None, "device": "cpu",
        "threads": 1, "event_policy": "panel_resolved_nonsynonymous" if v3 else "resolved_nonsynonymous",
        "fixed_background": "slac" if v4 else "focal", **INPUT_PATHS,
    }
    paths = dict(INPUT_PATHS)
    if v2:
        paths.update(ASR_PATHS)
        fixed.update(ancestral_method="hyphy-slac", asr_alignment=paths["asr_alignment"], asr_tree=paths["asr_tree"])
        if not isinstance(protocol["asr_node"], str) or not protocol["asr_node"].strip():
            raise BenchmarkInputError("asr_node must be a nonempty internal-node name")
        if "reconstruction_sha256" in protocol and not re.fullmatch(r"[0-9a-f]{64}", str(protocol["reconstruction_sha256"])):
            raise BenchmarkInputError("invalid reconstruction SHA-256")
    if v3:
        paths["substitution_panel"] = PANEL_PATH
        fixed["substitution_panel"] = PANEL_PATH
    if v4:
        fixed["endpoint_policy"] = "panel_sites_match_focal"
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
    if not isinstance(hashes, dict) or set(hashes) != set(paths):
        raise BenchmarkInputError(f"input_hashes must contain exactly {', '.join(sorted(paths))}")
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
    paths = dict(INPUT_PATHS)
    if protocol["schema_version"] in {"2.0", "3.0", "4.0"}:
        paths.update(asr_alignment=protocol["asr_alignment"], asr_tree=protocol["asr_tree"])
    if protocol["schema_version"] in {"3.0", "4.0"}:
        paths["substitution_panel"] = PANEL_PATH
    for key, fixed_path in paths.items():
        path = Path(fixed_path)
        if path.is_symlink() or not path.is_file():
            raise BenchmarkInputError(f"{key} must be a regular input file: {path}")
        if sha256_file(path) != protocol["input_hashes"][key]:
            raise BenchmarkInputError(f"{key} input SHA-256 does not match frozen protocol")
    if protocol["schema_version"] in {"2.0", "3.0", "4.0"}:
        path = Path(INPUT_PATHS["alignment"]).parent / "reconstruction.json"
        if path.is_symlink() or not path.is_file() or sha256_file(path) != protocol.get("reconstruction_sha256"):
            raise BenchmarkInputError("SLAC reconstruction SHA-256 does not match the frozen protocol")
        reconstruction = json.loads(path.read_text())
        if reconstruction.get("node") != protocol["asr_node"]:
            raise BenchmarkInputError("SLAC reconstruction node differs from the protocol")
        panel = (read_substitution_panel(protocol["substitution_panel"], protocol["alignment"], protocol["focal_taxon"])
                 if protocol["schema_version"] in {"3.0", "4.0"} else None)
        selection = select_reconstructed(read_codon_alignment(protocol["alignment"]), protocol["focal_taxon"], reconstruction, protocol["max_genotypes"], panel=panel, fixed_background=protocol["fixed_background"])
        if not selection["budget_ok"]:
            raise BenchmarkInputError(selection["budget_error"])
        return selection
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
    """Construct fixed engine arguments; the panel carries no ordering hints."""
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
    namespace = order_namespace(protocol, scorer, output)
    if protocol["schema_version"] in {"2.0", "3.0", "4.0"}:
        namespace._benchmark_selection = selection
        namespace._benchmark_reconstruction = json.loads((Path(protocol["alignment"]).parent / "reconstruction.json").read_text())
    _configure_determinism()
    from .cli import cmd_order
    return cmd_order(namespace)


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
