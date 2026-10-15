"""Extract Yokoyama 2014 S1 Table without treating inferred spectra as assays.

Requires PyMuPDF (``import fitz``). Run from any directory. The outputs are
derived references, not additional experiments. The source's duplicated
genotype label and the disagreement between its table and path totals remain
explicit; the script never modifies phenotype labels to match a published sum.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import re

import fitz


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "epistaeon/data/experimental"
MUTATIONS = ["F46T", "F49L", "T52F", "F86L", "T93P", "A114G", "S118T"]
TABLE = DATA / "sws1_2014_S1_Table.pdf"


def extract_records():
    page = fitz.open(TABLE)[0]
    spans = [span for block in page.get_text("dict")["blocks"]
             if "lines" in block for line in block["lines"] for span in line["spans"]]
    rows = defaultdict(list)
    for word in page.get_text("words"):
        if 90 < word[1] < 680 and re.fullmatch(r"-?\d+", word[4]):
            rows[round(word[1], 1)].append(word)
    records = []
    for ordinal, (y, words) in enumerate(sorted(rows.items()), 1):
        for side in ("left", "right"):
            cells = sorted([w for w in words if (w[0] > 290) == (side == "right")])
            assert len(cells) == (8 if ordinal == 1 and side == "left" else 9)
            printed = "".join(w[4] for w in cells[:7])
            genotype = printed
            notes = []
            if side == "right" and ordinal == 53:
                # The published row repeats a fourth-order genotype. In its
                # fifth-order enumeration it occupies the sole missing label.
                assert printed == "0011011" and y == 570.8
                genotype = "0111011"
                notes.append("Probable genotype-label typo: printed0011011 is duplicated; "
                             "normalized to the sole missing genotype0111011 by exhaustive "
                             "combination enumeration and its position in the five-mutant block. "
                             "Not independently confirmed by an author correction.")
            lambda_word = cells[7]
            fonts = [s["font"] for s in spans
                     if s["bbox"][0] <= lambda_word[0] + 0.2
                     and s["bbox"][2] >= lambda_word[2] - 0.2
                     and abs(s["bbox"][1] - lambda_word[1]) < 1]
            assert len(fonts) == 1
            functional = "Bold" not in fonts[0]
            if fonts[0] == "Arial-BoldMT":
                assert genotype == "1010011"
                notes.append("Lambda is bold but not italic in the PDF; classified with "
                             "the source's inferred/nonfunctional values. This yields its "
                             "reported60 nonfunctional mutants. It is not a measured spectrum.")
            value = int(lambda_word[4])
            delta = int(cells[8][4]) if len(cells) == 9 else None
            if functional and delta is not None and value - 357 != delta:
                notes.append("Published delta-lambda differs from lambda-minus357; "
                             "both printed numbers are preserved and derived paths use lambda.")
            records.append({
                "genotype": genotype,
                "mask": sum(int(bit) << i for i, bit in enumerate(genotype)),
                "mutations": [m for m, bit in zip(MUTATIONS, genotype) if bit == "1"],
                "functional_pigment": functional,
                "lambda_max_nm": value if functional else None,
                "published_inferred_lambda_max_nm": None if functional else value,
                "lambda_status": "measured" if functional else "inferred_nonfunctional",
                "published_table_delta_lambda_nm": delta,
                "source": {"table": "S1 Table", "pdf_page": 1, "column": side,
                           "row_1based": ordinal, "pdf_y0": y,
                           "printed_genotype": printed, "lambda_font": fonts[0]},
                "curation_notes": notes,
            })
    assert len(records) == 128 and len({r["genotype"] for r in records}) == 128
    assert sum(not r["functional_pigment"] for r in records) == 60
    return sorted(records, key=lambda record: record["mask"])


def enumerate_paths(records):
    lookup = {r["genotype"]: r for r in records}
    historical = []
    for first, middle, last in product(permutations((4, 5)),
                                       permutations((1, 6)), permutations((0, 2))):
        historical.append(list(first) + [3] + list(middle) + list(last))
    accessible = []
    for order in permutations(range(7)):
        bits = ["0"] * 7
        genotypes = ["".join(bits)]
        for i in order:
            bits[i] = "1"
            genotypes.append("".join(bits))
        if not all(lookup[bits]["functional_pigment"] for bits in genotypes):
            continue
        wavelengths = [lookup[bits]["lambda_max_nm"] for bits in genotypes]
        differences = [b - a for a, b in zip(wavelengths, wavelengths[1:])]
        accessible.append({
            "mutation_order": [MUTATIONS[i] for i in order],
            "genotypes_including_start": genotypes,
            "lambda_max_nm": wavelengths,
            "step_delta_lambda_nm": differences,
            "every_absolute_step_below25_nm": all(abs(d) < 25 for d in differences),
            "historically_favored": list(order) in historical,
        })
    assert len(accessible) == 1076  # Table reproduction, not the published1032.
    assert sum(p["historically_favored"] for p in accessible) == 8
    # Explicitly test the disjunctive prerequisite within the assayed cube.
    assert all(r["genotype"][0] == "1" or r["genotype"][4] == "1"
               for r in records if r["functional_pigment"] and r["genotype"][2] == "1")
    return accessible


def main():
    records = extract_records()
    paths = enumerate_paths(records)
    source = {"doi": "10.1371/journal.pgen.1004884",
              "url": "https://ndownloader.figshare.com/files/1847638",
              "local_file": str(TABLE.relative_to(ROOT)),
              "sha256": sha256(TABLE.read_bytes()).hexdigest()}
    metadata = {
        "schema_version": "1.0", "study_id": "sws1", "source": source,
        "mutation_bit_order": MUTATIONS,
        "bit_convention": "Genotype strings run in mutation_bit_order; integer mask uses the first mutation as bit0. 0=ancestral,1=derived.",
        "assay_scope": "Chimeric AncBoreotheria S1 background with anole terminal segments; purified pigment absorption, not organismal fitness or dark-noise kinetics.",
        "published_mutant_counts": {"assayed": 127, "nonfunctional": 60},
        "curation_warning": "One probable duplicated genotype-label correction and one bold-only inferred value are explicit in individual records. Published versus table-recomputed path counts disagree; preserve the distinction.",
    }
    panel = {**metadata, "records_including_ancestral_control": records}
    path_reference = {
        **metadata,
        "enumeration": "Enumerate all7! orders, retaining paths for which every visited genotype is functional in normalized S1 Table. No inferred nonfunctional lambda value enters a path spectrum.",
        "published_counts": {"total_orders": 5040, "accessible": 1032, "blocked": 4008,
                             "every_absolute_step_below25_nm": 335,
                             "accessible_by_first_mutation": dict(zip(
                                 ["F46T", "F49L", "T52F", "F86L", "T93P", "A114G", "S118T"],
                                 [134, 74, 0, 252, 348, 102, 122]))},
        "recomputed_from_normalized_S1_Table": {
            "accessible": len(paths), "blocked": 5040 - len(paths),
            "every_absolute_step_below25_nm": sum(p["every_absolute_step_below25_nm"] for p in paths),
            "accessible_by_first_mutation": dict(sorted(Counter(
                p["mutation_order"][0] for p in paths).items()))},
        "historical_partial_order": [["T93P", "A114G"], ["F86L"],
                                     ["F49L", "S118T"], ["F46T", "T52F"]],
        "historical_order_evidence": "Parsimony reconstruction across15 mammals, main Results and Fig4/5; pair order is unresolved. This is not an experimentally obligatory four-stage order.",
        "historically_favored_paths": [p for p in paths if p["historically_favored"]],
        "all_accessible_paths_from_normalized_table": paths,
    }
    for suffix, document in (("genotypes", panel), ("paths", path_reference)):
        path = DATA / f"sws1_2014_{suffix}.json"
        path.write_text(json.dumps(document, indent=2) + "\n")
        print(f"Wrote {path.relative_to(ROOT)}")
    print("S1 Table:128 genotypes,60 nonfunctional;1076 accessible orders "
          "(paper reports1032),336 gradual paths (paper reports335);8 historical paths.")


if __name__ == "__main__":
    main()
