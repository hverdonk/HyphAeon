"""Check the literature reference against its sources and current model scaffolds.

Run from any directory. --write-tsv refreshes the derived flat residue list.
This checks data integrity, not biological truth or model performance.
"""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from itertools import permutations, product
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "epistaeon/data"
REFERENCE = DATA / "substitution_reference"
STUDIES = ("hemoglobin", "myoglobin", "steroid_receptor", "sws1")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reference_columns(gene):
    record = next(
        r for r in SeqIO.parse(DATA / f"alignments/toga2/{gene}.codonified.fa", "fasta")
        if r.id == "hg38"
    )
    sequence = str(record.seq)
    result = {}
    position = 0
    for start in range(0, len(sequence), 3):
        if sequence[start:start + 3] != "---":
            position += 1
            result[position] = start // 3 + 1
    return result


def check_acyclic(constraints):
    edges = defaultdict(list)
    for constraint in constraints:
        if constraint["type"] == "phylogenetic_branch_precedence":
            for parent in constraint["prerequisites"]:
                edges[parent].extend(constraint["dependent"])
    visiting, complete = set(), set()

    def visit(node):
        require(node not in visiting, f"Cyclic branch precedence at {node}")
        if node in complete:
            return
        visiting.add(node)
        for child in edges.get(node, []):
            visit(child)
        visiting.remove(node)
        complete.add(node)

    for node in list(edges):
        visit(node)


def check_sws1_panel():
    """Audit the separate measured-genotype and inferred-history references."""
    panel = json.loads((DATA / "experimental/sws1_2014_genotypes.json").read_text())
    paths = json.loads((DATA / "experimental/sws1_2014_paths.json").read_text())
    mutations = panel["mutation_bit_order"]
    records = panel["records_including_ancestral_control"]
    lookup = {r["genotype"]: r for r in records}
    require(len(records) == len(lookup) == 128, "Incomplete SWS1 genotype cube")
    require(sum(not r["functional_pigment"] for r in records) == 60,
            "SWS1 nonfunctional genotype count changed")
    for row in records:
        bits = row["genotype"]
        require(len(bits) == 7 and set(bits) <= {"0", "1"}, "Invalid SWS1 genotype")
        require(row["mutations"] == [m for m, b in zip(mutations, bits) if b == "1"],
                "SWS1 genotype/mutation mismatch")
        require(row["mask"] == sum(int(b) << i for i, b in enumerate(bits)),
                "SWS1 genotype/mask mismatch")
        require((row["lambda_max_nm"] is not None) == row["functional_pigment"],
                "Inferred nonfunctional spectrum treated as a measured SWS1 spectrum")
        if row["functional_pigment"] and bits[2] == "1":
            require(bits[0] == "1" or bits[4] == "1", "T52F rescue rule changed")
    require(lookup["0000000"]["lambda_max_nm"] == 357 and
            lookup["1111111"]["lambda_max_nm"] == 411, "SWS1 endpoints changed")
    require(not lookup["0010000"]["functional_pigment"], "T52F-alone result changed")
    accessible = set()
    gradual = 0
    for order in permutations(range(7)):
        bits = ["0"] * 7
        visited = [lookup["".join(bits)]]
        for i in order:
            bits[i] = "1"
            visited.append(lookup["".join(bits)])
        if all(row["functional_pigment"] for row in visited):
            accessible.add(tuple(mutations[i] for i in order))
            gradual += all(abs(b["lambda_max_nm"] - a["lambda_max_nm"]) < 25
                           for a, b in zip(visited, visited[1:]))
    stored = paths["all_accessible_paths_from_normalized_table"]
    require(len(stored) == len(accessible) == 1076 and gradual == 336,
            "Normalized SWS1 table no longer reproduces documented counts")
    require({tuple(p["mutation_order"]) for p in stored} == accessible,
            "Stored SWS1 path set is stale")
    expected_historical = {
        tuple(a) + ("F86L",) + tuple(b) + tuple(c)
        for a, b, c in product(permutations(("T93P", "A114G")),
                               permutations(("F49L", "S118T")),
                               permutations(("F46T", "T52F")))
    }
    require({tuple(p["mutation_order"]) for p in paths["historically_favored_paths"]}
            == expected_historical <= accessible, "SWS1 historical paths changed")
    require(paths["published_counts"]["accessible"] == 1032,
            "Published SWS1 count was replaced by a recomputed value")
    print("sws1 assay: 128 genotypes; 8 historical paths; normalized table gives "
          "1076 accessible paths versus published 1032 (unresolved source discrepancy)")


def check_reference():
    offsets = json.loads((DATA / "numbering_offsets.json").read_text())
    models = {
        pdb: MMCIFParser(QUIET=True).get_structure(
            pdb, str(DATA / f"coordinates/{pdb}.cif")
        )[0]
        for pdb in ("2HHB", "1MBO", "7PRX", "1U19")
    }
    columns = {gene: reference_columns(gene) for gene in
               ("HBA1", "HBA2", "HBB", "MB", "NR3C1")}
    sws1_map = {row["source_position"]: row
                for row in offsets["sws1_benchmark"]["site_map"]}
    sws1_human = next(r for r in SeqIO.parse(
        DATA / "alignments/sws1/SWS1.codonified.fa", "fasta")
        if r.id == "Homo_sapiens")
    flat = []
    total_events = total_mappings = 0
    for study in STUDIES:
        document = json.loads((REFERENCE / f"{study}.json").read_text())
        require(document["study_id"] == study, f"Wrong study ID in {study}")
        require(document["schema_version"] == "1.0", f"Unknown schema in {study}")
        sources = {s["id"]: s for s in document["sources"]}
        for source in sources.values():
            if source.get("local_file"):
                path = ROOT / source["local_file"]
                require(path.is_file(), f"Missing source: {path}")
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                require(digest == source["local_sha256"], f"Changed source: {path}")
        events = {e["id"]: e for e in document["events"]}
        groups = {g["id"]: g for g in document["groups"]}
        require(len(events) == len(document["events"]), f"Duplicate event in {study}")
        require(len(groups) == len(document["groups"]), f"Duplicate group in {study}")
        for group in groups.values():
            for member in group["members"]:
                require(member in events, f"Unknown member {member} in {study}")
                require(events[member]["group"] == group["id"], f"Group mismatch: {member}")
        for constraint in document["order_constraints"]:
            require(constraint["locator"] and constraint["scope"], "Unsourced constraint")
            require(set(constraint["source_ids"]) <= sources.keys(), "Unknown constraint source")
            for member in constraint["prerequisites"] + constraint["dependent"]:
                unwrapped = re.sub(r"^(?:reverse|completion)\((.*)\)$", r"\1", member)
                require(unwrapped in events or unwrapped in groups,
                        f"Unknown constraint member: {study}/{member}")
        check_acyclic(document["order_constraints"])
        for event in events.values():
            total_events += 1
            require(event["locator"] and event["structural_region"], f"Incomplete event: {event['id']}")
            require(set(event["source_ids"]) <= sources.keys(), f"Unknown source: {event['id']}")
            require(event["group"] in groups, f"Unknown group: {event['id']}")
            for mapping in event["structure_mapping"]:
                position = mapping["author_position"]
                if event["kind"] == "deletion":
                    require(position is None and mapping["uniprot_position"] is None,
                            f"Deleted site assigned a modern residue: {event['id']}")
                if position is not None:
                    residue = models[mapping["pdb"]][mapping["chain"]][(" ", position, " ")]
                    require(seq1(residue.resname) == mapping["observed_residue"],
                            f"Residue mismatch: {event['id']}")
                    require(np.allclose(residue["CA"].coord, mapping["ca_xyz_angstrom"], atol=0.001),
                            f"Coordinate mismatch: {event['id']}")
                    chain = offsets["structures"][mapping["pdb"]]["chains"][mapping["chain"]]
                    canonical = position + chain["offset_to_uniprot"]
                    require(canonical == mapping["uniprot_position"], f"Offset mismatch: {event['id']}")
                    reference_position = canonical + (1 if mapping["pdb"] == "7PRX" else 0)
                    for gene, column in mapping["toga2_hg38_alignment_columns_1based"].items():
                        require(columns[gene][reference_position] == column,
                                f"Alignment column mismatch: {event['id']}/{gene}")
                    if study == "sws1":
                        crosswalk = sws1_map[event["position"]]
                        column = crosswalk["alignment_codon_column_1based"]
                        require(mapping["alignment_reference_columns_1based"] == {"SWS1": column},
                                f"SWS1 alignment column mismatch: {event['id']}")
                        require(mapping["human_cds_position_1based"] == crosswalk["human_cds_position_1based"],
                                f"SWS1 human position mismatch: {event['id']}")
                        require(str(sws1_human.seq[(column - 1) * 3:column * 3].translate()) == event["to"],
                                f"SWS1 derived human residue mismatch: {event['id']}")
                        require(position == crosswalk["template_author_position"],
                                f"SWS1 template mismatch: {event['id']}")
                    total_mappings += 1
                flat.append({
                    "case": document["case_number"], "study": study,
                    "event": event["id"], "label": event["label"], "kind": event["kind"],
                    "source_position": event["position"], "from": event["from"], "to": event["to"],
                    "interval": event["interval"], "group": event["group"],
                    "structural_region": event["structural_region"], "pdb": mapping["pdb"],
                    "chain": mapping["chain"], "author_position": position,
                    "modern_residue": mapping.get("observed_residue"),
                    "uniprot_position": mapping["uniprot_position"],
                    "hg38_codon_columns_1based": json.dumps(
                        mapping.get("toga2_hg38_alignment_columns_1based", {}), sort_keys=True),
                    "other_reference_codon_columns_1based": json.dumps(
                        mapping.get("alignment_reference_columns_1based", {}), sort_keys=True),
                    "other_reference_row": mapping.get("alignment_reference_row"),
                    "other_reference_residue_1based": mapping.get("human_cds_position_1based"),
                    "source_ids": ",".join(event["source_ids"]), "locator": event["locator"],
                    "epistasis_scoring": event.get("epistasis_scoring"),
                    "order_scoring": event.get("order_scoring"),
                })
        print(f"{study}: {len(events)} events, {len(groups)} groups, "
              f"{len(document['order_constraints'])} qualified constraints")
    check_sws1_panel()
    return flat, total_events, total_mappings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-tsv", action="store_true")
    args = parser.parse_args()
    rows, event_count, mapping_count = check_reference()
    path = REFERENCE / "sites.tsv"
    if args.write_tsv:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    else:
        with path.open(newline="") as handle:
            actual = list(csv.DictReader(handle, delimiter="\t"))
        expected = [{key: "" if value is None else str(value) for key, value in row.items()}
                    for row in rows]
        require(actual == expected, "sites.tsv is stale; regenerate with --write-tsv")
    print(f"PASS: {event_count} events, {mapping_count} coordinate mappings; "
          "source hashes and reference alignment columns checked. No model run performed.")


if __name__ == "__main__":
    main()
