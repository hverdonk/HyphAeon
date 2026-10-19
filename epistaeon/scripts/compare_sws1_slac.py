#!/usr/bin/env python3
"""Evaluation-only comparison of frozen SLAC ASR with Yokoyama 2014 S4.

No inference or model code is imported. The human core, rather than reconstructed
ancestral residues, determines the coordinate bridge. All inference artifacts
must already be covered by the independent inference checksum manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path, PurePosixPath

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = dict(zip(("".join(codon) for codon in itertools.product(BASES, repeat=3)), AMINO_ACIDS))
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
         "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
         "H": "ACT", "V": "ACG", "N": "ACGT", "?": "ACGT"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def read_fasta(path):
    rows, name = {}, None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            names = line[1:].split()
            if not names or names[0] in rows:
                raise ValueError(f"Empty or duplicate FASTA identifier: {path}")
            name = names[0]
            rows[name] = ""
        elif name is None:
            raise ValueError(f"Sequence preceding FASTA header: {path}")
        else:
            rows[name] += "".join(line.split()).upper()
    if not rows or any(not sequence for sequence in rows.values()):
        raise ValueError(f"Empty FASTA sequence: {path}")
    return rows


def translate_codon(codon):
    if codon == "---":
        return "-"
    if len(codon) != 3 or any(base not in IUPAC for base in codon):
        raise ValueError(f"Invalid or partial-gap codon: {codon}")
    residues = {CODE["".join(concrete)] for concrete in itertools.product(*(IUPAC[base] for base in codon))}
    return next(iter(residues)) if len(residues) == 1 else "X"


def translate(sequence):
    if len(sequence) % 3:
        raise ValueError("Codon sequence length is not divisible by three")
    return "".join(translate_codon(sequence[index:index + 3]) for index in range(0, len(sequence), 3))


def verify_inference_seal(seal_path, base_dir, required_paths):
    """Verify every frozen input/result and require the actual comparison inputs."""
    base_dir = Path(base_dir).resolve()
    manifest = read_json(seal_path)
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("Inference seal must contain a nonempty files map")
    for name, expected in files.items():
        relative = PurePosixPath(name)
        if (relative.is_absolute() or ".." in relative.parts or str(relative) != name
                or not isinstance(expected, str) or len(expected) != 64
                or any(character not in "0123456789abcdef" for character in expected)):
            raise ValueError(f"Invalid inference seal entry: {name}")
        path = base_dir / name
        if any(parent.is_symlink() for parent in [path, *path.parents] if parent != base_dir
               and base_dir in parent.parents):
            raise ValueError(f"Symlink in frozen inference artifacts: {name}")
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Inference checksum verification failed: {name}")
    for path in required_paths:
        try:
            name = Path(path).resolve().relative_to(base_dir).as_posix()
        except ValueError as exc:
            raise ValueError("Comparison inference inputs must be inside the frozen run") from exc
        if name not in files:
            raise ValueError(f"Comparison input was not frozen before evaluation: {name}")
    return manifest


def mask_positions(reference):
    ancestor, human = reference["ancestor_protein"], reference["paper_human_protein"]
    if len(ancestor) != len(human) or not ancestor:
        raise ValueError("S4 ancestor and human rows must have equal nonzero lengths")
    if set(ancestor + human) - (AA | {"-"}):
        raise ValueError("S4 rows contain unresolved or invalid residues")
    mask = reference["terminal_mask"]
    core_ranges = mask["comparable_core_ranges_1based"]
    engineered_ranges = mask["engineered_ranges_1based"]
    if len(core_ranges) != 1:
        raise ValueError("SWS1 comparison expects one contiguous natural core")

    def positions(ranges):
        values = []
        for start, end in ranges:
            if type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(ancestor):
                raise ValueError("Invalid S4 mask coordinates")
            values.extend(range(start, end + 1))
        if len(values) != len(set(values)):
            raise ValueError("Overlapping S4 mask ranges")
        return set(values)

    core, engineered = positions(core_ranges), positions(engineered_ranges)
    if core & engineered or core | engineered != set(range(1, len(ancestor) + 1)):
        raise ValueError("Core and engineered masks must partition the displayed S4 row")
    return core_ranges[0], core, engineered


def build_human_bridge(paper_human, core_range, aligned_human):
    """Map the unique identical human core; never align prediction to truth."""
    start, end = core_range
    core = paper_human[start - 1:end].replace("-", "")
    natural = aligned_human.replace("-", "")
    if not core:
        raise ValueError("S4 human core is empty")
    first = natural.find(core)
    if first < 0 or natural.find(core, first + 1) >= 0:
        raise ValueError("S4 human core must occur exactly once in the input human sequence")
    human_to_column = {}
    human_position = 0
    for column, residue in enumerate(aligned_human, 1):
        if residue != "-":
            human_position += 1
            human_to_column[human_position] = column
    bridge, inverse = {}, {}
    human_position = first
    for paper_position in range(start, end + 1):
        if paper_human[paper_position - 1] != "-":
            human_position += 1
            item = {"human_cds_position_1based": human_position,
                    "alignment_column_1based": human_to_column[human_position]}
            bridge[paper_position] = item
            inverse[human_position] = paper_position
    metadata = {
        "method": "Unique exact match of the complete S4-human non-engineered core to the ungapped input Homo_sapiens; human residues then map through actual input gaps.",
        "paper_core_range_1based": [start, end],
        "human_core_range_1based": [first + 1, first + len(core)],
        "alignment_core_mapped_column_range_1based": [min(item["alignment_column_1based"] for item in bridge.values()),
                                                     max(item["alignment_column_1based"] for item in bridge.values())],
        "mapped_paper_positions": len(bridge),
        "unmapped_human_gap_positions": [position for position in range(start, end + 1) if position not in bridge],
    }
    return bridge, inverse, metadata


def summarize(rows):
    residues = [row for row in rows if row["paper_ancestor"] in AA]
    resolved = [row for row in residues if row["status"] in {"match", "substitution"}]
    matched = [row for row in residues if row["status"] == "match"]
    mismatches = [row for row in residues if row["status"] == "substitution"]
    missing = [row for row in residues if row["status"] not in {"match", "substitution"}]
    return {"expected_residue_count": len(residues), "resolved_reconstruction_count": len(resolved),
            "matched_residue_count": len(matched), "substitution_disagreement_count": len(mismatches),
            "unresolved_or_unmapped_count": len(missing),
            "identity_among_resolved": len(matched) / len(resolved) if resolved else None,
            "strict_match_fraction_all_expected_residues": len(matched) / len(residues) if residues else None,
            "reference_gap_columns_excluded_from_amino_acid_identity": sum(row["paper_ancestor"] == "-" for row in rows),
            "disagreements": mismatches, "unresolved_or_unmapped": missing}


def _status(expected, observed):
    if expected == "-":
        return "paper_gap"
    if observed not in AA:
        return "unresolved_reconstruction"
    return "match" if observed == expected else "substitution"


def _write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare(alignment, ancestors_codon, ancestors_protein, selected_node, paper_reference,
            substitution_reference, inference_seal, run_dir, output):
    inputs = {name: Path(path).resolve() for name, path in {
        "alignment": alignment, "ancestors_codon": ancestors_codon, "ancestors_protein": ancestors_protein,
        "selected_node": selected_node, "paper_reference": paper_reference,
        "substitution_reference": substitution_reference, "inference_seal": inference_seal}.items()}
    run_dir, output = Path(run_dir).resolve(), Path(output).resolve()
    if any(path == output or output in path.parents for path in inputs.values()):
        raise ValueError("Comparison output must not contain its input files")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Comparison output must be new or empty")
    frozen = [inputs[name] for name in ("alignment", "ancestors_codon", "ancestors_protein", "selected_node")]
    verify_inference_seal(inputs["inference_seal"], run_dir, frozen)
    hashes = {name: sha256(path) for name, path in inputs.items()}
    paper, empirical = read_json(inputs["paper_reference"]), read_json(inputs["substitution_reference"])
    selected = read_json(inputs["selected_node"])
    node_id = selected.get("node_id")
    if not isinstance(node_id, str) or not node_id:
        raise ValueError("selected_node.json must identify one frozen node_id")
    extant = read_fasta(inputs["alignment"])
    codon_rows, protein_rows = read_fasta(inputs["ancestors_codon"]), read_fasta(inputs["ancestors_protein"])
    if node_id not in codon_rows or node_id not in protein_rows:
        raise ValueError("Selected node is absent from ancestral FASTA files")
    if "Homo_sapiens" not in extant:
        raise ValueError("Frozen input alignment lacks Homo_sapiens")
    ancestral_codons = codon_rows[node_id]
    reconstructed = translate(ancestral_codons)
    if reconstructed != protein_rows[node_id].replace("?", "X"):
        raise ValueError("Ancestral nucleotide/protein FASTA translations disagree")
    human = translate(extant["Homo_sapiens"])
    if len(human) != len(reconstructed) or any(len(row) != 3 * len(human) for row in extant.values()):
        raise ValueError("Input alignment and ancestral sequence widths disagree")
    core_range, core, engineered = mask_positions(paper)
    bridge, inverse, bridge_metadata = build_human_bridge(paper["paper_human_protein"], core_range, human)
    core_rows = []
    for position in sorted(core):
        mapping = bridge.get(position, {})
        column = mapping.get("alignment_column_1based")
        prediction = reconstructed[column - 1] if column else None
        expected = paper["ancestor_protein"][position - 1]
        core_rows.append({"paper_protein_position_1based": position,
                          "human_cds_position_1based": mapping.get("human_cds_position_1based"),
                          "alignment_column_1based": column, "paper_ancestor": expected,
                          "paper_human": paper["paper_human_protein"][position - 1],
                          "reconstructed_ancestor": prediction,
                          "reconstructed_codon": ancestral_codons[3 * (column - 1):3 * column] if column else None,
                          "status": _status(expected, prediction) if column else "unmapped_reference_insertion"})
    by_position = {row["paper_protein_position_1based"]: row for row in core_rows}
    focal_rows = []
    for event in empirical["events"]:
        if not event.get("order_scoring", True):
            continue
        human_positions = {mapping["human_cds_position_1based"] for mapping in event["structure_mapping"]
                           if mapping.get("human_cds_position_1based") is not None}
        if len(human_positions) != 1:
            raise ValueError(f"Ambiguous human coordinate for {event['id']}")
        human_position = next(iter(human_positions))
        paper_position = inverse.get(human_position)
        if paper_position is None:
            raise ValueError(f"Focal event lies outside the comparable natural core: {event['id']}")
        row = dict(by_position[paper_position])
        if row["paper_ancestor"] != event["from"] or row["paper_human"] != event["to"]:
            raise ValueError(f"S4 transcription and empirical focal crosswalk disagree: {event['id']}")
        focal_rows.append({"substitution": event["label"], "paper_source_position_1based": event["position"], **row})
    if len(focal_rows) != 7:
        raise ValueError("SWS1 comparison requires exactly seven independently mapped focal events")
    report = {
        "schema_version": 1, "analysis": "Evaluation-only comparison of frozen HyPhy SLAC point ASR",
        "node_id": node_id, "selected_node": selected,
        "focal_sites": summarize(focal_rows), "focal_site_rows": focal_rows,
        "paper_natural_core": summarize(core_rows), "coordinate_bridge": bridge_metadata,
        "mask": {"basis": paper["terminal_mask"].get("basis"),
                 "engineered_ranges_1based": paper["terminal_mask"]["engineered_ranges_1based"],
                 "excluded_engineered_positions": len(engineered),
                 "methods_prose_verbatim_bounds": paper["terminal_mask"].get("methods_prose_verbatim_bounds")},
        "source_evidence": paper.get("source_files", []),
        "input_sha256": hashes, "comparison_script_sha256": sha256(Path(__file__)),
        "limitations": [
            "S4 displays terminal-standardized assay chimeras; compare only its non-engineered core with a natural ASR.",
            "Full natural ancestral termini and true ancestral codons are unavailable; no nucleotide-accuracy claim is possible.",
            "The rebuilt taxon panel, alignment, phylogeny and SLAC reconstruction do not reproduce the original author analysis exactly.",
            "Point-state agreement does not measure posterior uncertainty, substitution order or epistasis.",
            "Prediction gaps, ambiguous states and unmapped insertions are reported as unavailable, separately from resolved substitutions.",
        ],
    }
    verify_inference_seal(inputs["inference_seal"], run_dir, frozen)
    if hashes != {name: sha256(path) for name, path in inputs.items()}:
        raise ValueError("An input changed during comparison")
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "comparison.json", report)
    _write_tsv(output / "focal_sites.tsv", focal_rows)
    _write_tsv(output / "core_residues.tsv", core_rows)
    focal, full = report["focal_sites"], report["paper_natural_core"]
    lines = [
        "# SWS1 SLAC ancestral-state comparison", "",
        f"The frozen `{node_id}` reconstruction matches **{focal['matched_residue_count']}/{focal['expected_residue_count']} focal ancestral residues** and **{full['matched_residue_count']}/{full['expected_residue_count']} available paper-core residues**.", "",
        f"Among {full['resolved_reconstruction_count']} resolved core comparisons, there are {full['substitution_disagreement_count']} substitutions and {full['unresolved_or_unmapped_count']} unresolved or unmapped positions. Resolved amino-acid identity is {full['identity_among_resolved']:.2%}." if full['identity_among_resolved'] is not None else "No core positions have resolved comparisons.", "",
        "| Substitution | Paper ancestor | SLAC ancestor | Reconstructed codon | Alignment column | Result |",
        "|---|---|---|---|---:|---|",
    ]
    for row in focal_rows:
        lines.append(f"| {row['substitution']} | {row['paper_ancestor']} | {row['reconstructed_ancestor']} | {row['reconstructed_codon']} | {row['alignment_column_1based']} | {row['status']} |")
    lines.extend([
        "", "Coordinates were bridged by the unique identical S4-human core, independently of the reconstructed ancestor. Actual input gaps are retained. The figure's red terminal blocks are excluded; its sequential protein indices differ from the Methods prose bounds. [S4 Fig](https://ndownloader.figshare.com/files/1847637), [Yokoyama et al. 2014](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004884).", "",
        "The paper provides a chimeric assay protein, rather than a complete natural ancestor or ancestral DNA. Consequently this comparison covers the 282-residue non-engineered core and the seven focal amino acids; it cannot assess true ancestral codons or natural termini. This is a separate SLAC reconstruction using rebuilt inputs, not an exact rerun of the original reconstruction. Agreement does not validate substitution order or epistasis.", "",
        "Inference artifact checksums were verified before reading empirical truth and again before writing these outputs. No inference settings or predictions were modified. Full core discrepancies, unavailable states and coordinate mappings are in `core_residues.tsv`; hashes and detailed metrics are in `comparison.json`.", "",
    ])
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")
    write_json(output / "manifest.json", {"schema_version": 1, "input_sha256": hashes,
               "comparison_script_sha256": sha256(Path(__file__)),
               "outputs": {path.name: sha256(path) for path in sorted(output.iterdir()) if path.is_file()}})
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("alignment", "ancestors-codon", "ancestors-protein", "selected-node", "paper-reference",
                 "substitution-reference", "inference-seal", "run-dir", "output"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        report = compare(**vars(args))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"SWS1 comparison: {exc}\n")
    print(f"Wrote comparison for {report['node_id']} to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
