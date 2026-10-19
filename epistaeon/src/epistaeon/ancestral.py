"""Shared, reference-free HyPhy SLAC joint-ML ancestral reconstruction."""

import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from .benchmark_inputs import read_codon_alignment, sha256_file


class AncestralError(ValueError):
    pass


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def discover_hyphy(executable=None):
    path = shutil.which(str(executable or "hyphy"))
    if not path:
        raise AncestralError("HyPhy is required for SLAC reconstruction; install it or supply --hyphy. No plurality fallback is available")
    try:
        result = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise AncestralError(f"Cannot query the selected HyPhy executable: {exc}") from exc
    version = (result.stdout + result.stderr).strip()
    if result.returncode or not re.match(r"HYPHY 2\.5\.96(?:\(|\s|$)", version):
        raise AncestralError(f"SLAC requires HyPhy 2.5.96; found {version!r}")
    return str(Path(path).resolve()), version


def _tree(text):
    from Bio import Phylo
    try:
        return Phylo.read(io.StringIO(text), "newick")
    except Exception as exc:
        raise AncestralError(f"Invalid ASR Newick tree: {exc}") from exc


def _target(tree, node, focal, taxa):
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != set(taxa):
        raise AncestralError("ASR tree tips must match the ASR alignment exactly, with unique names")
    matches = list(tree.find_clades(name=node))
    if len(matches) != 1 or matches[0].is_terminal():
        raise AncestralError("--asr-node must identify exactly one named internal node")
    descendants = sorted(tip.name for tip in matches[0].get_terminals())
    if focal not in descendants:
        raise AncestralError("The focal taxon must descend from --asr-node")
    if matches[0] is tree.root and len(tree.root.clades) == 2:
        raise AncestralError("A degree-two ASR root is suppressed by SLAC; supply outgroups and an identifiable internal target")
    return matches[0], descendants


def validate_inputs(alignment, focal_taxon, tree, node, asr_alignment=None):
    scoring = read_codon_alignment(alignment)
    augmented = read_codon_alignment(asr_alignment or alignment)
    if focal_taxon not in scoring:
        raise AncestralError(f"Focal taxon {focal_taxon!r} is absent from the scoring alignment")
    if any(name not in augmented or augmented[name] != sequence for name, sequence in scoring.items()):
        raise AncestralError("ASR alignment must contain identical codon columns and sequences for every scoring taxon")
    if len(next(iter(augmented.values()))) != len(next(iter(scoring.values()))):
        raise AncestralError("ASR and scoring alignments must have the same codon columns")
    parsed = _tree(Path(tree).read_text())
    _, descendants = _target(parsed, node, focal_taxon, augmented)
    return scoring, augmented, descendants


def extract_result(payload, node, focal_taxon, rows, descendants):
    """Read site arrays without confusing a multi-AA cell with extra sites."""
    from Bio.Seq import Seq
    try:
        if payload["input"]["partition count"] != 1:
            raise AncestralError("SLAC reconstruction requires one alignment partition")
        tree = _tree(payload["input"]["trees"]["0"])
        attributes = payload["branch attributes"]["0"]
        # Apply fitted lengths while names still match HyPhy's attribute keys.
        for branch in tree.find_clades():
            values = attributes.get(branch.name, {})
            if "Global MG94xREV" in values:
                branch.branch_length = values["Global MG94xREV"]
        # HyPhy can normalize identifiers. Restore biological tip names before
        # validating clade identity, but retain the original node identifier.
        for tip in tree.get_terminals():
            tip.name = attributes.get(tip.name, {}).get("original name", tip.name)
        _, observed = _target(tree, node, focal_taxon, rows)
        if observed != descendants:
            raise AncestralError("SLAC changed the selected node's descendant taxa")
        codons = attributes[node]["codon"][0]
        states = attributes[node]["amino-acid"][0]
    except (KeyError, TypeError, IndexError) as exc:
        raise AncestralError("SLAC output lacks the selected internal node or its sequence arrays") from exc
    width = len(next(iter(rows.values())))
    if not isinstance(codons, list) or not isinstance(states, list) or len(codons) != width or len(states) != width:
        raise AncestralError("SLAC ancestral arrays do not match alignment columns")
    for codon, amino_acids in zip(codons, states):
        if not isinstance(codon, str) or len(codon) != 3 or set(codon) - set("ACGTRYSWKMBDHVN?-"):
            raise AncestralError("SLAC returned an invalid ancestral codon")
        if "-" in codon and codon != "---":
            raise AncestralError("SLAC returned a partial-gap codon")
        if not isinstance(amino_acids, str) or not amino_acids or set(amino_acids) - set("ACDEFGHIKLMNPQRSTVWY-X"):
            raise AncestralError("SLAC returned invalid amino-acid states")
        if not set(codon) - set("ACGT"):
            translated = str(Seq(codon).translate())
            if translated != amino_acids:
                raise AncestralError("SLAC ancestral codons and amino-acid states disagree")
    from Bio import Phylo
    stream = io.StringIO()
    Phylo.write(tree, stream, "newick")
    return {"schema_version": "1.0", "method": "hyphy-slac", "node": node,
            "descendant_taxa": observed, "codons": codons, "amino_acid_states": states,
            "protein": "".join(state if len(state) == 1 else "X" for state in states),
            "fitted_tree": stream.getvalue().strip(), "n_codons": width,
            "uncertainty": "Joint maximum-likelihood point estimate; no residue posterior probabilities"}


def reconstruct(alignment, focal_taxon, tree, node, output, asr_alignment=None, hyphy=None):
    scoring, rows, descendants = validate_inputs(alignment, focal_taxon, tree, node, asr_alignment)
    executable, version = discover_hyphy(hyphy)
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise AncestralError("ASR output must be absent or empty; frozen reconstructions are not overwritten")
    output.mkdir(parents=True, exist_ok=True)
    inputs = {"alignment": str(Path(alignment).resolve()), "asr_alignment": str(Path(asr_alignment or alignment).resolve()),
              "tree": str(Path(tree).resolve())}
    hashes = {key: sha256_file(path) for key, path in inputs.items()}
    library = Path(executable).parents[1] / "share/hyphy"
    command = [executable, "CPU=1", f"LIBPATH={library}", "ENV=SetParameter(RANDOM_SEED,0,0);",
               "slac", "--alignment", inputs["asr_alignment"],
               "--tree", inputs["tree"], "--code", "Universal", "--branches", "All",
               "--samples", "0", "--pvalue", "0.1", "--kill-zero-lengths", "No",
               "--output", str(output / "slac.json")]
    libraries = {p.relative_to(library).as_posix(): sha256_file(p) for p in sorted(library.rglob("*"))
                 if p.is_file()}
    if not libraries:
        raise AncestralError(f"Cannot fingerprint the HyPhy batch library at {library}")
    metadata = {"method": "hyphy-slac", "hyphy_version": version,
                "executable_sha256": sha256_file(executable), "library_path": str(library),
                "library_sha256": libraries, "seed": 0,
                "command": command, "input_hashes": hashes, "node": node,
                "descendant_taxa": descendants, "scoring_taxa": list(scoring)}
    write_json(output / "protocol.json", metadata)
    start = time.monotonic()
    with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
        result = subprocess.run(command, cwd=output, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
                                env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"))
    write_json(output / "status.json", {"exit_code": result.returncode, "elapsed_seconds": time.monotonic() - start})
    if result.returncode:
        raise AncestralError(f"HyPhy SLAC failed with exit status {result.returncode}; see {output}/stderr.log")
    if hashes != {key: sha256_file(path) for key, path in inputs.items()}:
        raise AncestralError("ASR input changed during reconstruction")
    try:
        payload = json.loads((output / "slac.json").read_text())
        reconstruction = extract_result(payload, node, focal_taxon, rows, descendants)
    except (OSError, json.JSONDecodeError) as exc:
        raise AncestralError("HyPhy did not produce a valid SLAC result") from exc
    reconstruction["provenance"] = metadata
    write_json(output / "reconstruction.json", reconstruction)
    (output / "ancestor.codon.fa").write_text(f">{node}\n{''.join(reconstruction['codons'])}\n")
    (output / "ancestor.protein.fa").write_text(f">{node}\n{reconstruction['protein']}\n")
    (output / "fitted_tree.nwk").write_text(reconstruction["fitted_tree"] + "\n")
    return reconstruction
