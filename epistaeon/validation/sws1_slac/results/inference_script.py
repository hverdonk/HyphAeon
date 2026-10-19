#!/usr/bin/env python3
"""Run a separate, truth-free HyPhy SLAC ancestral-reconstruction experiment.

The original 15-taxon codon alignment is preserved. Three paper-listed SWS1
outgroups identify the Boreoeutherian node; their proteins are added to the fixed
alignment and backtranslated with MAFFT's explicit residue map. Empirical
ancestral sequences are never inputs to this script.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import Bio
from Bio import Phylo, SeqIO
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "epistaeon/validation/sws1_slac"
ORIGINAL = ROOT / "epistaeon/data/alignments/sws1"
OUTGROUPS = {
    "Loxodonta_africana": "AY686753.1",
    "Notamacropus_eugenii": "AY286017.1",
    "Sminthopsis_crassicaudata": "AY442173.1",
}
EXPECTED_HYPHY = "HYPHY 2.5.96(MP)"
EXPECTED_MAFFT = "v7.526 (2024/Apr/26)"


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def fasta(path, rows):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in rows.items()))


def read_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    rows = {record.id: str(record.seq).upper() for record in records}
    if len(rows) != len(records) or not rows:
        raise ValueError(f"Empty FASTA or duplicate names: {path}")
    return rows


def clade(name, *children):
    return "(" + ",".join(children) + ")" + name


def topology():
    hominini = clade("Hominini", "Homo_sapiens", "Pan_troglodytes")
    homininae = clade("Homininae", hominini, "Gorilla_gorilla")
    hominidae = clade("Hominidae", homininae, "Pongo_abelii")
    hominoidea = clade("Hominoidea", hominidae, "Nomascus_leucogenys")
    catarrhini = clade("Catarrhini", hominoidea, "Macaca_fascicularis")
    platyrrhini = clade("Platyrrhini", "Saimiri_boliviensis", "Callithrix_jacchus")
    simiiformes = clade("Simiiformes", catarrhini, platyrrhini)
    haplorhini = clade("Haplorhini", simiiformes, "Cephalopachus_bancanus")
    primates = clade("Primates", haplorhini, "Daubentonia_madagascariensis")
    murinae = clade("Murinae", "Mus_musculus", "Rattus_norvegicus")
    rodentia = clade("Rodentia", murinae, "Sciurus_carolinensis", "Cavia_porcellus")
    euarchontoglires = clade("Euarchontoglires", primates, rodentia)
    boreoeutheria = clade("Boreoeutheria", "Bos_taurus", euarchontoglires)
    placentalia = clade("Placentalia", "Loxodonta_africana", boreoeutheria)
    # Degree-three computational root at the marsupial MRCA; the target is the
    # explicitly named INTERNAL Boreoeutheria node, not the computational root.
    return clade("ASRRoot", "Notamacropus_eugenii", "Sminthopsis_crassicaudata", placentalia) + ";\n"


def prepare(base, mafft):
    inputs = base / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    downloads = json.loads((inputs / "outgroup_downloads.json").read_text())
    pinned = {item["accession"]: item for item in downloads}
    coding, proteins, metadata = {}, {}, []
    for identifier, accession in OUTGROUPS.items():
        path = inputs / "raw" / f"{accession}.embl"
        if digest(path) != pinned[accession]["sha256"]:
            raise ValueError(f"Outgroup record changed: {path}")
        record = SeqIO.read(path, "embl")
        features = [feature for feature in record.features if feature.type == "CDS"]
        if record.id != accession or len(features) != 1:
            raise ValueError(f"Unexpected accession or CDS features: {path}")
        feature = features[0]
        if feature.qualifiers.get("codon_start", ["1"]) != ["1"]:
            raise ValueError(f"Unexpected CDS frame: {path}")
        sequence = str(feature.extract(record.seq)).upper()
        translated = str(Seq(sequence).translate())
        if not translated.endswith("*") or "*" in translated[:-1]:
            raise ValueError(f"Unexpected stop codons: {path}")
        protein = translated[:-1]
        if protein != feature.qualifiers["translation"][0]:
            raise ValueError(f"CDS disagrees with annotated translation: {path}")
        coding[identifier], proteins[identifier] = sequence[:-3], protein
        metadata.append({"taxon": identifier, "accession": accession,
                         "record_organism": record.annotations["organism"],
                         "record_sha256": digest(path), "source_protein_length": len(protein),
                         "terminal_stop_removed": True, "cds_location": str(feature.location)})
    fasta(inputs / "outgroups.cds.fa", coding)
    additions = inputs / "outgroups.protein.fa"
    fasta(additions, proteins)
    version = subprocess.run([mafft, "--version"], capture_output=True, text=True, check=True)
    mafft_version = (version.stdout + version.stderr).strip()
    if mafft_version != EXPECTED_MAFFT:
        raise ValueError(f"Expected MAFFT {EXPECTED_MAFFT}, got {mafft_version}")
    command = [mafft, "--quiet", "--amino", "--globalpair", "--thread", "1", "--inputorder",
               "--add", str(additions), "--keeplength", "--mapout",
               str(ORIGINAL / "SWS1.protein.aligned.fa")]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    (inputs / "mafft.log").write_text(result.stderr)
    protein_path = inputs / "SWS1.outgroups.protein.fa"
    protein_path.write_text(result.stdout)
    aligned = read_fasta(protein_path)
    original_codons = read_fasta(ORIGINAL / "SWS1.codonified.fa")
    original_proteins = read_fasta(ORIGINAL / "SWS1.protein.aligned.fa")
    if any(aligned[name] != seq for name, seq in original_proteins.items()):
        raise ValueError("MAFFT changed an original alignment column")
    if len(aligned) != 18 or {len(seq) for seq in aligned.values()} != {350}:
        raise ValueError("Expected 18 taxa and 350 fixed protein columns")
    codon_rows = dict(original_codons)
    for identifier in OUTGROUPS:
        codon_rows[identifier] = ["---"] * 350
    discarded = {identifier: [] for identifier in OUTGROUPS}
    mapped = {identifier: set() for identifier in OUTGROUPS}
    current = None
    for line in Path(str(additions) + ".map").read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
        elif line.strip() and not line.startswith("#"):
            aa, source, column = [field.strip() for field in line.split(",")]
            source = int(source)
            if current not in OUTGROUPS or proteins[current][source - 1] != aa:
                raise ValueError("Invalid MAFFT residue mapping")
            codon = coding[current][3 * (source - 1):3 * source]
            if column == "-":
                discarded[current].append({"source_residue": source, "amino_acid": aa, "codon": codon})
            else:
                column = int(column)
                if column in mapped[current] or aligned[current][column - 1] != aa:
                    raise ValueError("Duplicate or inconsistent mapped alignment column")
                mapped[current].add(column)
                codon_rows[current][column - 1] = codon
    for identifier in OUTGROUPS:
        codon_rows[identifier] = "".join(codon_rows[identifier])
        if str(Seq(codon_rows[identifier]).translate(gap="-")) != aligned[identifier]:
            raise ValueError(f"Outgroup backtranslation failed: {identifier}")
        if len(mapped[identifier]) + len(discarded[identifier]) != len(proteins[identifier]):
            raise ValueError(f"Incomplete source mapping: {identifier}")
    if any(codon_rows[name] != seq for name, seq in original_codons.items()):
        raise ValueError("Original codon sequence changed")
    fasta(inputs / "SWS1.outgroups.codon.fa", codon_rows)
    tree_path = inputs / "species_scaffold.nwk"
    tree_path.write_text(topology())
    tree = Phylo.read(tree_path, "newick")
    if {tip.name for tip in tree.get_terminals()} != set(codon_rows):
        raise ValueError("Tree/alignment taxon mismatch")
    target = next(tree.find_clades(name="Boreoeutheria"))
    if {tip.name for tip in target.get_terminals()} != set(original_codons):
        raise ValueError("Target node is not exactly the 15-taxon ingroup MRCA")
    dump(inputs / "preparation.json", {
        "original_taxa": list(original_codons), "outgroups": metadata,
        "original_15_rows_preserved_exactly": True, "original_350_columns_preserved_exactly": True,
        "mafft_version": mafft_version, "biopython_version": Bio.__version__, "mafft_command": command,
        "outgroup_insertions_excluded_to_preserve_reference_columns": discarded,
        "target_node": "Boreoeutheria", "tree_policy": "Fixed species scaffold; no topology search; rodent deep relationships left as a polytomy; branch lengths estimated by HyPhy",
        "computational_root": "Marsupial MRCA (degree three), distinct from the target Boreoeutheria node",
        "sources": ["https://doi.org/10.1371/journal.pgen.1004884", "https://doi.org/10.1371/journal.pgen.1001342", "https://mafft.cbrc.jp/alignment/software/addsequences.html"],
    })


def extract(base):
    results = base / "results"
    payload = json.loads((results / "slac.json").read_text())
    output_tree = Phylo.read(io.StringIO(payload["input"]["trees"]["0"]), "newick")
    ingroup = set(read_fasta(ORIGINAL / "SWS1.codonified.fa"))
    targets = list(output_tree.find_clades(name="Boreoeutheria"))
    if len(targets) != 1 or {tip.name for tip in targets[0].get_terminals()} != ingroup:
        raise ValueError("SLAC tree does not identify the target as exactly the original ingroup MRCA")
    attributes = payload["branch attributes"]["0"]
    codons, proteins = {}, {}
    for node, values in attributes.items():
        if "codon" not in values or "amino-acid" not in values:
            continue
        codon = values["codon"][0]
        protein = values["amino-acid"][0]
        codons[node] = "".join(codon) if isinstance(codon, list) else codon
        proteins[node] = "".join(protein) if isinstance(protein, list) else protein
        if len(codons[node]) != 1050 or len(proteins[node]) != 350:
            raise ValueError(f"Unexpected reconstructed sequence length: {node}")
    if "Boreoeutheria" not in codons:
        raise ValueError("HyPhy output did not retain the explicit Boreoeutheria target node")
    for node in codons:
        translated = str(Seq(codons[node]).translate(gap="-"))
        if translated != proteins[node]:
            raise ValueError(f"HyPhy ancestral codons and amino acids disagree: {node}")
    fasta(results / "ancestors.codon.fa", codons)
    fasta(results / "ancestors.protein.fa", proteins)
    fasta(results / "Boreoeutheria.codon.fa", {"Boreoeutheria": codons["Boreoeutheria"]})
    fasta(results / "Boreoeutheria.protein.fa", {"Boreoeutheria": proteins["Boreoeutheria"]})
    dump(results / "selected_node.json", {
        "node_id": "Boreoeutheria", "label": "Boreoeutheria", "method": "SLAC joint maximum-likelihood codon ancestral reconstruction",
        "identified_by": "Internal node subtending exactly the 15 original Boreoeutherian taxa in the frozen species-scaffold topology; elephant and two marsupials are outside",
        "input_alignment": "inputs/SWS1.outgroups.codon.fa", "input_tree": "inputs/species_scaffold.nwk",
        "codon_fasta": "results/Boreoeutheria.codon.fa", "protein_fasta": "results/Boreoeutheria.protein.fa",
        "codon_columns": 350, "uncertainty": "Joint ML point estimate; no per-residue posterior probabilities computed; samples=0",
    })


def run(base, hyphy):
    results = base / "results"
    results.mkdir(parents=True, exist_ok=True)
    if (results / "slac.json").exists() or (results / "inference_seal.json").exists():
        raise ValueError("Results already exist; use a new --directory to run another experiment")
    version = subprocess.run([hyphy, "--version"], capture_output=True, text=True, check=True)
    hyphy_version = (version.stdout + version.stderr).strip()
    if not hyphy_version.startswith(EXPECTED_HYPHY):
        raise ValueError(f"Expected {EXPECTED_HYPHY}, got {hyphy_version}")
    command = [hyphy, "CPU=1", "slac", "--alignment", str(base / "inputs/SWS1.outgroups.codon.fa"),
               "--tree", str(base / "inputs/species_scaffold.nwk"), "--code", "Universal",
               "--branches", "All", "--samples", "0", "--pvalue", "0.1", "--kill-zero-lengths", "No",
               "--output", str(results / "slac.json")]
    input_hashes = {str(path.relative_to(base)): digest(path) for path in sorted((base / "inputs").rglob("*")) if path.is_file()}
    hyphy_library = Path(hyphy).resolve().parents[1] / "share/hyphy"
    library_hashes = {path.relative_to(hyphy_library).as_posix(): digest(path)
                      for path in sorted(hyphy_library.rglob("*.bf"))}
    dump(base / "inference_protocol.json", {
        "schema_version": "1.0", "hyphy_version": hyphy_version, "hyphy_executable_sha256": digest(hyphy),
        "script_sha256": digest(__file__), "original_alignment_sha256": digest(ORIGINAL / "SWS1.codonified.fa"),
        "hyphy_batch_library_sha256": library_hashes,
        "input_hashes": input_hashes, "command": command, "cpu": 1, "target_node": "Boreoeutheria",
        "settings": {"genetic_code": "Universal", "reconstruction": "joint ML under SLAC's fitted MG94xREV model", "samples": 0, "pvalue": .1, "kill_zero_lengths": False},
        "empirical_reference_used_in_inference": False,
    })
    started = datetime.now(timezone.utc).isoformat()
    beginning = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (results / "stdout.log").open("w") as stdout, (results / "stderr.log").open("w") as stderr:
        completed = subprocess.run(command, cwd=results, stdout=stdout, stderr=stderr, env=env)
    dump(results / "status.json", {"started_at": started, "finished_at": datetime.now(timezone.utc).isoformat(),
                                   "elapsed_seconds": time.monotonic() - beginning, "exit_code": completed.returncode,
                                   "command": command})
    if completed.returncode:
        raise RuntimeError("SLAC failed; inspect results/stdout.log and stderr.log")
    if any(digest(base / name) != value for name, value in input_hashes.items()):
        raise ValueError("Inference inputs changed during SLAC execution")
    extract(base)
    paths = [base / "inference_protocol.json", *sorted((base / "inputs").rglob("*")), *sorted(results.rglob("*"))]
    files = {str(path.relative_to(base)): digest(path) for path in paths if path.is_file()}
    dump(results / "inference_seal.json", {"algorithm": "sha256", "files": files,
                                           "scope": "Inference protocol, all input and result files; empirical evaluation occurs separately"})
    print(f"SLAC completed and inference artifacts sealed at {results}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=BASE)
    parser.add_argument("--hyphy", default=shutil.which("hyphy"))
    parser.add_argument("--mafft", default=shutil.which("mafft"))
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--run-only", action="store_true")
    args = parser.parse_args()
    if not args.hyphy or not args.mafft:
        parser.error("HyPhy and MAFFT must be installed")
    base = args.directory.resolve()
    if not args.run_only:
        prepare(base, args.mafft)
    if not args.prepare_only:
        run(base, args.hyphy)


if __name__ == "__main__":
    main()
