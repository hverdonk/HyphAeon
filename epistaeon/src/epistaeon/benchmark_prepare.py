"""Isolated SLAC preparation; no checkpoint or numerical inference imports."""

import argparse
import json
from pathlib import Path

from .ancestral import AncestralError, reconstruct
from .benchmark_inputs import BenchmarkInputError, read_codon_alignment, select_reconstructed, sha256_file
from .benchmark_worker import validate_protocol


def prepare(protocol_path, output):
    protocol = json.loads(Path(protocol_path).read_text())
    validate_protocol(protocol)
    if protocol["schema_version"] != "2.0":
        raise BenchmarkInputError("New benchmark preparation requires SLAC protocol 2.0")
    output = Path(output)
    if not output.is_absolute() or output.is_symlink() or Path("/outputs") not in output.resolve().parents:
        raise BenchmarkInputError("Preparation output must be a real subdirectory of /outputs")
    for name, expected in protocol["input_hashes"].items():
        path = Path(protocol[name])
        if path.is_symlink() or not path.is_file() or sha256_file(path) != expected:
            raise BenchmarkInputError(f"Preparation input hash mismatch: {name}")
    reconstruction = reconstruct(protocol["alignment"], protocol["focal_taxon"], protocol["asr_tree"],
                                 protocol["asr_node"], output, protocol["asr_alignment"])
    selection = select_reconstructed(read_codon_alignment(protocol["alignment"]), protocol["focal_taxon"],
                                     reconstruction, protocol["max_genotypes"])
    (output / "selection.json").write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        return prepare(args.protocol, args.output)
    except (AncestralError, BenchmarkInputError, OSError, json.JSONDecodeError) as exc:
        parser.exit(2, f"SLAC preparation: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
