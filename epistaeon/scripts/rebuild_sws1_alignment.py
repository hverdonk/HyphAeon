#!/usr/bin/env python3
"""Rebuild the extant, codon-aware SWS1 panel from saved accession records.

The normal operation is entirely offline. MAFFT aligns translated proteins;
backtranslation consumes each source codon exactly once, without nucleotide
realignment or imputation. --check rebuilds in a temporary directory and checks
the committed products; --download-missing retrieves only missing pinned inputs.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess
import tempfile
import urllib.request

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATA = Path(__file__).resolve().parents[1] / "data/alignments/sws1"
MAFFT_ARGS = ["--quiet", "--amino", "--globalpair", "--maxiterate", "1000", "--thread", "1", "--inputorder"]
TAXA = [
    ("Homo_sapiens", "Homo sapiens", "M13295", "M13299.embl"),
    ("Pan_troglodytes", "Pan troglodytes", "AF039433", "AF039435.embl"),
    ("Gorilla_gorilla", "Gorilla gorilla gorilla", "XM_004046176", "XM_004046176.1.gb"),
    ("Pongo_abelii", "Pongo abelii", "XM_002818421", "XM_002818421.1.gb"),
    ("Nomascus_leucogenys", "Nomascus leucogenys", "XM_003261297", "XM_003261297.1.gb"),
    ("Macaca_fascicularis", "Macaca fascicularis", "AF158977", "AF158977.embl"),
    ("Saimiri_boliviensis", "Saimiri sciureus", "U53875", "U53875.embl"),
    ("Callithrix_jacchus", "Callithrix jacchus", "L76201", "L76201.embl"),
    ("Cephalopachus_bancanus", "Tarsius bancanus", "AB111463", "AB111463.embl"),
    ("Daubentonia_madagascariensis", "Daubentonia madagascariensis", "EF667285", "EF667285.embl"),
    ("Bos_taurus", "Bos taurus", "U92557", "U92557.embl"),
    ("Mus_musculus", "Mus musculus", "U49720", "U49720.embl"),
    ("Rattus_norvegicus", "Rattus norvegicus", "U63972", "U63972.embl"),
    ("Sciurus_carolinensis", "Sciurus carolinensis", "DQ302163", "DQ302163.embl"),
    ("Cavia_porcellus", "Cavia porcellus", "AY552608", "AY552608.embl"),
]


def digest(data: bytes) -> str:
    return sha256(data).hexdigest()


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fasta(path: Path, rows: dict[str, str]) -> None:
    records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in rows.items()]
    SeqIO.write(records, path, "fasta")


def prepare_inputs(data: Path, download_missing: bool) -> tuple[dict, dict]:
    manifest = json.loads((data / "manifest.json").read_text())
    records = {}
    for item in manifest["raw_records"]:
        path = data / "raw" / item["file"]
        if not path.exists() and download_missing:
            path.parent.mkdir(parents=True, exist_ok=True)
            content = urllib.request.urlopen(item["url"], timeout=60).read()
            if digest(content) != item["sha256"]:
                raise ValueError(f"Downloaded record differs from pinned input: {path.name}")
            path.write_bytes(content)
        content = path.read_bytes()
        if digest(content) != item["sha256"]:
            raise ValueError(f"Raw record checksum mismatch: {path}")
        record = SeqIO.read(io.StringIO(content.decode()), item["format"])
        if record.id != item["record_id"]:
            raise ValueError(f"Record accession mismatch: {path}")
        records[item["file"]] = record
    return manifest, records


def extract(records: dict) -> tuple[dict, dict, dict, list]:
    references = {record.id: record.seq for record in records.values()}
    source_cds, coding, proteins, taxa = {}, {}, {}, []
    for identifier, paper_species, cited_accession, filename in TAXA:
        record = records[filename]
        features = [f for f in record.features if f.type == "CDS"]
        notes = []
        excluded_partial = ""
        if filename == "L76201.embl":
            if features or len(record.seq) != 1030 or not str(record.seq).startswith("ATG"):
                raise ValueError("Unexpected unannotated L76201 source")
            # This record contains a CDS fragment in frame 0. The single final
            # nucleotide cannot define a codon and is retained only in raw/.
            cds = record.seq[:1029]
            excluded_partial = str(record.seq[1029:])
            feature = None
            notes.append("No CDS feature: frame 0 starts ATG and translates 343 codons without stops; supported by full-length protein homology. One terminal T is excluded as a partial codon; unsequenced C terminus is not filled.")
        else:
            if len(features) != 1:
                raise ValueError(f"Expected one annotated CDS: {filename}")
            feature = features[0]
            if feature.qualifiers.get("codon_start", ["1"]) != ["1"]:
                raise ValueError(f"Unexpected CDS frame: {filename}")
            cds = feature.extract(record.seq, references=references)
        if len(cds) % 3:
            raise ValueError(f"Non-triplet CDS: {identifier}")
        translated = str(cds.translate(table=1))
        terminal_stop = translated.endswith("*")
        protein = translated[:-1] if terminal_stop else translated
        if "*" in protein:
            raise ValueError(f"Internal stop codon: {identifier}")
        if feature and protein != feature.qualifiers["translation"][0]:
            raise ValueError(f"Translation disagrees with source annotation: {identifier}")
        if len(protein) < 340 or len(protein) > 360:
            raise ValueError(f"Unexpected SWS1 protein length: {identifier}")
        if identifier == "Homo_sapiens":
            notes.append("Full historical CDS is explicitly annotated in M13299.1, joining M13295.1 through M13299.1. Preserve its 348-residue sequence rather than replacing it with current 345-residue UniProt P03999.")
        if identifier == "Pan_troglodytes":
            notes.append("Full CDS is explicitly annotated in AF039435.1, joining AF039433.1, AF039434.1 and AF039435.1.")
        if identifier == "Saimiri_boliviensis":
            notes.append("Paper labels U53875 as Saimiri sciureus; the current source record identifies Saimiri boliviensis. Accession sequence is preserved; taxon identifier follows the source record.")
            notes.append("Record title says partial CDS, but the annotated CDS has an initiating ATG and terminal stop, with no internal stop or ambiguous nucleotide. This discrepancy is retained in metadata.")
        if identifier == "Cephalopachus_bancanus":
            notes.append("Paper uses the historical name Tarsius bancanus; current source uses Cephalopachus bancanus.")
        coding_sequence = str(cds[:-3] if terminal_stop else cds).upper()
        ambiguous = [{"source_codon": i // 3 + 1, "codon": coding_sequence[i:i+3], "amino_acid": protein[i // 3]}
                     for i in range(0, len(coding_sequence), 3) if set(coding_sequence[i:i+3]) - set("ACGT")]
        if ambiguous:
            notes.append("Source ambiguity is retained as IUPAC nucleotides; translated unresolved codons remain X, without imputation.")
        source_cds[identifier] = str(cds).upper()
        coding[identifier] = coding_sequence
        proteins[identifier] = protein
        taxon = {"id": identifier, "paper_species": paper_species, "record_organism": record.annotations["organism"],
                 "cited_accession": cited_accession, "cds_record": record.id, "cds_record_file": filename,
                 "record_description": record.description, "source_cds_length_nt": len(cds), "protein_length_aa": len(protein),
                 "terminal_stop_removed_for_alignment": terminal_stop, "excluded_partial_codon_bases": excluded_partial,
                 "cds_feature_location": str(feature.location) if feature else None,
                 "cds_feature_qualifiers": dict(feature.qualifiers) if feature else None,
                 "ambiguous_codons": ambiguous, "notes": notes}
        taxa.append(taxon)
    return source_cds, coding, proteins, taxa


def build(data: Path, out: Path, manifest: dict, records: dict, mafft: str) -> None:
    version = subprocess.run([mafft, "--version"], capture_output=True, text=True, check=True)
    version_string = (version.stdout + version.stderr).strip()
    if version_string != manifest["software"]["mafft"]:
        raise ValueError(f"MAFFT version must be {manifest['software']['mafft']}; got {version_string}")
    if Bio.__version__ != manifest["software"]["biopython"]:
        raise ValueError(f"Biopython version must be {manifest['software']['biopython']}; got {Bio.__version__}")
    source_cds, coding, proteins, taxa = extract(records)
    out.mkdir(parents=True, exist_ok=True)
    fasta(out / "SWS1.source_cds.fa", source_cds)
    fasta(out / "SWS1.cds.fa", coding)
    fasta(out / "SWS1.protein.fa", proteins)
    aligned = subprocess.run([mafft, *MAFFT_ARGS, str(out / "SWS1.protein.fa")], capture_output=True, text=True, check=True)
    aligned_records = list(SeqIO.parse(io.StringIO(aligned.stdout), "fasta"))
    aligned_proteins = {r.id: str(r.seq).upper() for r in aligned_records}
    if list(aligned_proteins) != list(proteins):
        raise ValueError("Aligner changed row identity or order")
    aligned_codons = {}
    for identifier, protein in aligned_proteins.items():
        consumed = 0
        codons = []
        for aa in protein:
            if aa == "-":
                codons.append("---")
            else:
                codon = coding[identifier][consumed:consumed+3]
                if str(Seq(codon).translate()) != aa:
                    raise ValueError(f"Backtranslation mismatch: {identifier}")
                codons.append(codon)
                consumed += 3
        if consumed != len(coding[identifier]) or protein.replace("-", "") != proteins[identifier]:
            raise ValueError(f"Source sequence not preserved: {identifier}")
        aligned_codons[identifier] = "".join(codons)
    widths = {len(s) for s in aligned_codons.values()}
    if len(widths) != 1 or next(iter(widths)) % 3 or len(aligned_codons) != 15:
        raise ValueError("Alignment dimension mismatch")
    fasta(out / "SWS1.protein.aligned.fa", aligned_proteins)
    fasta(out / "SWS1.codonified.fa", aligned_codons)
    dump(out / "taxa.json", taxa)
    human_position = 0
    map_rows = ["alignment_codon\talignment_nt_start\thuman_residue\thuman_amino_acid\thuman_codon"]
    for i, aa in enumerate(aligned_proteins["Homo_sapiens"]):
        if aa != "-":
            human_position += 1
        map_rows.append(f"{i+1}\t{3*i+1}\t{human_position if aa != '-' else ''}\t{aa}\t{aligned_codons['Homo_sapiens'][3*i:3*i+3]}")
    (out / "human_coordinate_map.tsv").write_text("\n".join(map_rows) + "\n")
    quality = {"status": "pass_with_documented_source_limitations", "taxon_count": len(taxa),
               "extant_only": True, "ancestral_rows": 0, "alignment_width_nt": next(iter(widths)),
               "alignment_width_codons": next(iter(widths)) // 3, "human_row": "Homo_sapiens",
               "human_length_aa": len(proteins["Homo_sapiens"]),
               "source_codon_sequences_preserved_exactly": True, "triplet_gaps_only": True,
               "protein_backtranslation_exact": True, "internal_stop_codons": 0,
               "retained_ambiguous_codons": sum(len(t["ambiguous_codons"]) for t in taxa),
               "partial_rows": ["Callithrix_jacchus"], "terminal_stops_excluded": sum(t["terminal_stop_removed_for_alignment"] for t in taxa),
               "excluded_partial_codon_bases": {t["id"]: t["excluded_partial_codon_bases"] for t in taxa if t["excluded_partial_codon_bases"]},
               "alignment_gaps_codons_per_taxon": {k: v.count("-") for k, v in aligned_proteins.items()},
               "checks": ["raw accession record SHA256", "annotated CDS translation agrees with source", "no internal stops", "15 unique extant taxa", "equal alignment lengths divisible by three", "gap codons exactly ---", "ungapped codons reproduce source after documented terminal trimming", "translated codon alignment reproduces protein alignment"],
               "limitations": ["Rebuilt alignment, not the unreleased original author alignment.", "Taxa are the boreoeutherian union of both Methods accession lists; this is not asserted to reproduce the exact Figure 4 panel.", "Marmoset source lacks CDS annotation and full C terminus; excluded partial codon and missing codons are not imputed.", "Rat source has one ambiguous W nucleotide translating to X.", "Paper substitution positions require an explicit coordinate crosswalk; they are not necessarily human CDS indices."]}
    dump(out / "quality_report.json", quality)
    for filename, expected in manifest.get("products", {}).items():
        if digest((out / filename).read_bytes()) != expected:
            raise ValueError(f"Generated product differs from pinned SHA256: {filename}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--check", action="store_true", help="Rebuild offline in a temporary directory and compare products")
    parser.add_argument("--download-missing", action="store_true", help="Download missing raw records only; hashes must match pinned inputs")
    args = parser.parse_args()
    manifest, records = prepare_inputs(args.data_dir, args.download_missing)
    if args.check:
        with tempfile.TemporaryDirectory(prefix="sws1-rebuild-") as temp:
            out = Path(temp)
            build(args.data_dir, out, manifest, records, args.mafft)
            for path in sorted(out.iterdir()):
                if path.read_bytes() != (args.data_dir / path.name).read_bytes():
                    raise ValueError(f"Rebuild differs: {path.name}")
        print("PASS: 15 taxa; exact offline rebuild of all SWS1 alignment products")
    else:
        build(args.data_dir, args.data_dir, manifest, records, args.mafft)
        print("Built 15-taxon SWS1 codon alignment; see quality_report.json for source limitations")


if __name__ == "__main__":
    main()
