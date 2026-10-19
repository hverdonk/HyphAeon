"""Truth-free, standard-library preflight for the frozen benchmark profile.

This module deliberately imports no numerical/model code. The input itself is
the sole authority for which amino-acid changes enter the benchmark.
"""

from collections import Counter
import hashlib


SCORERS = (
    "potts-metropolis", "typicality-gate", "softmax-repaired", "dphi-product",
)
EXACT_UNIT_LIMIT = 22
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TO_AA = {
    a + b + c: _AA[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}


class BenchmarkInputError(ValueError):
    """Input does not meet the benchmark's fixed, exact-run profile."""


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def residue_of(codon):
    aa = CODON_TO_AA.get(codon.upper(), "-")
    return "-" if aa == "*" else aa


def read_codon_alignment(path):
    """Read strict, equal-length codon FASTA; retain every row in file order."""
    rows = {}
    name = None
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                tokens = line[1:].split()
                if not tokens:
                    raise BenchmarkInputError(f"empty FASTA header at line {line_number}")
                # Match aeon_core's FASTA identifier convention exactly.
                name = tokens[0].strip("'\"")
                if not name:
                    raise BenchmarkInputError(f"empty FASTA taxon at line {line_number}")
                if name in rows:
                    raise BenchmarkInputError(f"duplicate FASTA taxon {name!r}")
                rows[name] = ""
            else:
                if name is None:
                    raise BenchmarkInputError("benchmark input must be FASTA")
                seq = line.replace(" ", "").upper().replace("U", "T")
                if set(seq) - set("ACGTRYSWKMBDHVN?-"):
                    raise BenchmarkInputError(f"invalid nucleotide symbols in {name!r}")
                rows[name] += seq
    if len(rows) < 2:
        raise BenchmarkInputError("plurality inference needs at least two alignment taxa")
    lengths = {len(seq) for seq in rows.values()}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
        raise BenchmarkInputError("alignment rows must have the same nonzero length")
    if next(iter(lengths)) % 3:
        raise BenchmarkInputError("alignment length must be divisible by three")
    return {name: [seq[i:i + 3] for i in range(0, len(seq), 3)]
            for name, seq in rows.items()}


def _plurality(codons):
    # Keep background._plurality semantics, including first-seen codon ties.
    residues = Counter(residue_of(codon) for codon in codons)
    best = sorted(residues, key=lambda aa: (-residues[aa], aa == "-", aa))[0]
    codon = ("---" if best == "-" else
             Counter(c for c in codons if residue_of(c) == best).most_common(1)[0][0])
    return best, codon, residues[best] / len(codons)


def select_from_rows(rows, focal_taxon, max_genotypes=4096):
    """Select all resolved nonsynonymous changes, with no external site list.

    Every unselected column starts at the focal codon. Plurality codon ties use
    input-row order, just as the ordinary ``order`` command does.
    """
    if focal_taxon not in rows:
        raise BenchmarkInputError(f"focal taxon {focal_taxon!r} is not an exact FASTA row name")
    if len(rows) < 2:
        raise BenchmarkInputError("plurality inference needs at least two alignment taxa")
    if type(max_genotypes) is not int or max_genotypes < 1:
        raise BenchmarkInputError("max_genotypes must be a positive integer")
    n_codons = len(rows[focal_taxon])
    if not n_codons or any(len(row) != n_codons for row in rows.values()):
        raise BenchmarkInputError("alignment rows must have equal nonzero codon lengths")
    if any(len(codon) != 3 for row in rows.values() for codon in row):
        raise BenchmarkInputError("every alignment column must contain complete codons")
    others = [name for name in rows if name != focal_taxon]
    focal = [codon.upper() for codon in rows[focal_taxon]]
    start = list(focal)
    plurality_codons, columns, selected = [], [], []
    counts = Counter()
    residue_number = 0
    for col, derived_codon in enumerate(focal):
        nonfocal_codons = [rows[name][col].upper() for name in others]
        ancestral, ancestral_codon, support = _plurality(
            nonfocal_codons
        )
        derived = residue_of(derived_codon)
        plurality_codons.append(ancestral_codon)
        if derived != "-":
            residue_number += 1
        if ancestral_codon == derived_codon:
            kind = "unchanged"
        elif ancestral == "-" or derived == "-":
            # _plurality represents any unresolved winner as ---. An unknown
            # source codon must still be reported as unknown rather than a
            # reconstructed deletion simply because of that placeholder.
            unresolved_focal = derived == "-" and "-" not in derived_codon
            unresolved_source = ancestral == "-" and not any("-" in c for c in nonfocal_codons)
            kind = "unresolved" if unresolved_focal or unresolved_source else "indel"
        elif ancestral == derived:
            kind = "synonymous"
        else:
            kind = "nonsynonymous"
        counts[kind] += 1
        record = {
            "column0": col, "column1": col + 1,
            "residue_number": residue_number if derived != "-" else None,
            "ancestral": ancestral, "derived": derived,
            "ancestral_codon": ancestral_codon, "derived_codon": derived_codon,
            "plurality_support": support, "class": kind,
            "selected": kind == "nonsynonymous",
        }
        columns.append(record)
        if kind == "nonsynonymous":
            start[col] = ancestral_codon
            selected.append(dict(record, unit=len(selected),
                                 name=f"{ancestral}{residue_number}{derived}"))
    n_units = len(selected)
    if n_units < 1:
        raise BenchmarkInputError(
            "benchmark needs at least one resolved nonsynonymous unit; found none"
        )
    if n_units > EXACT_UNIT_LIMIT:
        raise BenchmarkInputError(
            f"exact benchmark supports at most {EXACT_UNIT_LIMIT} units; found {n_units}"
        )
    n_genotypes = 1 << n_units
    if n_genotypes > max_genotypes:
        raise BenchmarkInputError(
            f"exact benchmark needs {n_genotypes} genotypes for {n_units} units, "
            f"above max_genotypes {max_genotypes}; no sites are dropped"
        )
    full = list(start)
    for event in selected:
        full[event["column0"]] = event["derived_codon"]
    if full != focal:
        raise BenchmarkInputError("applying all selected units does not reproduce the focal sequence")
    return {
        "schema_version": "1.0", "focal_taxon": focal_taxon,
        "taxa": list(rows), "n_taxa": len(rows), "n_codons": n_codons,
        "event_policy": "resolved_nonsynonymous", "fixed_background": "focal",
        "plurality_residue_tie_break": "residue before gap, then alphabetical",
        "plurality_codon_tie_break": "first occurrence in input taxon order",
        "n_units": n_units, "n_genotypes": n_genotypes,
        "max_genotypes": max_genotypes, "selected_columns0": [s["column0"] for s in selected],
        "classification_counts": dict(sorted(counts.items())),
        "selected": selected, "columns": columns,
        "plurality_codons": plurality_codons, "starting_codons": start,
        "focal_codons": focal, "full_mask_reproduces_focal": full == focal,
    }


def select_benchmark_inputs(alignment, focal_taxon, max_genotypes=4096):
    """Public host/worker preflight, before any checkpoint is imported."""
    selection = select_from_rows(read_codon_alignment(alignment), focal_taxon, max_genotypes)
    selection["alignment_sha256"] = sha256_file(alignment)
    return selection
