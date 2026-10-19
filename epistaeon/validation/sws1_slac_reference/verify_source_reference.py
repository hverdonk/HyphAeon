#!/usr/bin/env python3
"""Verify the paper-only S4 transcription; never infer ancestral states."""
from hashlib import sha256
import csv
import json
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ref = json.loads((HERE / "source_reference.json").read_text())
ancestor = ref["ancestor_protein"]
human = ref["paper_human_protein"]
pigment_g = ref["pigment_g_protein"]
assert len(ancestor) == len(human) == len(pigment_g) == 347
assert sum(a != b for a, b in zip(ancestor, human)) == 25
assert sum(a != b for a, b in zip(ancestor, pigment_g)) == 9
assert "-" not in ancestor + human + pigment_g
for entry in ref["source_files"]:
    assert sha256((ROOT / entry["local_path"]).read_bytes()).hexdigest() == entry["sha256"]
rows = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / ref["human_coordinate_bridge"]["reference_human_protein"], "fasta")}
extant_human = rows["Homo_sapiens"]
core_human = human[25:307]
assert len(core_human) == 282 and extant_human.count(core_human) == 1
assert extant_human.index(core_human) == 27
assert ref["terminal_mask"]["engineered_ranges_1based"] == [[1, 25], [308, 347]]
assert ref["terminal_mask"]["comparable_core_ranges_1based"] == [[26, 307]]
alignment = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / "epistaeon/data/alignments/sws1/SWS1.codonified.fa", "fasta")}
assert len(alignment) == 15
with (HERE / "paper_S4.core_positions.tsv").open() as stream:
    core_map = list(csv.DictReader(stream, delimiter="\t"))
assert len(core_map) == 282
for row in core_map:
    p, hp, col = (int(row[key]) for key in ["paper_protein_position_1based", "historical_human_cds_position_1based", "benchmark_alignment_codon_column_1based"])
    assert row["ancestor_aa"] == ancestor[p - 1]
    assert row["paper_human_aa"] == human[p - 1] == extant_human[hp - 1]
    assert row["pigment_g_aa"] == pigment_g[p - 1]
    codon = alignment["Homo_sapiens"][3 * (col - 1):3 * col]
    from Bio.Seq import Seq
    assert str(Seq(codon).translate()) == human[p - 1]
for site in ref["focal_sites"]:
    p = site["paper_protein_position_1based"]
    col = site["benchmark_alignment_codon_column_1based"]
    assert ancestor[p - 1] == site["ancestor_state"]
    assert human[p - 1] == site["human_state"]
    assert all(set(seq[3 * (col - 1):3 * col]) <= set("ACGT") for seq in alignment.values())
manifest = json.loads((HERE / "source_manifest.json").read_text())
for filename, expected in manifest["reference_files_sha256"].items():
    assert sha256((HERE / filename).read_bytes()).hexdigest() == expected
print("PASS: S4 347-aa source rows; 25/9 differences; 282-aa unique human bridge; 7 focal states and 15-taxon codon coverage")
