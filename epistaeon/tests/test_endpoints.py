"""Endpoints come from alignment rows, and positions are numbered in one of them."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epistaeon.background import (  # noqa: E402
    GAP_CODON,
    GenotypeRows,
    codons_of,
    derive_endpoints,
    residue_of,
    resolve_positions,
    resolve_row,
)
from epistaeon.lattice import Lattice  # noqa: E402

# three columns: the second differs between anc and der, the third is a deletion
ANC = "AAA" + "AGC" + "TGC"          # K S C
DER = "AAA" + "CCC" + "---"          # K P -
OTHER = "AAA" + "AGC" + "TGC"        # K S C, same as the ancestor


def _rows(**seqs):
    return {name: codons_of(seq, 3) for name, seq in seqs.items()}


def test_substitutions_between_two_rows():
    rows = _rows(anc=ANC, der=DER, out=OTHER)
    ends = derive_endpoints(rows, "der", "anc")
    assert [s.name() for s in ends.substitutions] == ["S2P", "C3del"]
    assert ends.ancestral_protein == "KSC"
    assert ends.reference_name == "anc"
    assert ends.plurality_support == {}
    assert ends.class_counts() == {"nonsynonymous": 1, "indel": 1}


def test_a_synonymous_codon_change_is_a_substitution():
    # column 1 keeps lysine but changes codon; column 3 is nonsynonymous
    rows = _rows(anc="AAA" + "TGC", der="AAG" + "TGG")
    ends = derive_endpoints(rows, "der", "anc")
    assert [s.name() for s in ends.substitutions] == ["K1=AAA>AAG", "C2W"]
    assert ends.class_counts() == {"synonymous": 1, "nonsynonymous": 1}
    syn = ends.substitutions[0]
    assert syn.ancestral == syn.derived == "K"


def test_applying_every_unit_reproduces_the_derived_row():
    """The whole point of ordering codons: the trajectory ends at the row.

    With amino-acid units, the synonymous columns keep ancestral codons and the
    end of the trajectory stops short of the derived sequence -- which a
    nucleotide distance then measures as a gap that is not there.
    """
    anc, der = "AAA" + "AGC" + "TGC", "AAG" + "AGT" + "TGG"
    rows = _rows(anc=anc, der=der)
    ends = derive_endpoints(rows, "der", "anc")
    built = list(ends.ancestral_codons)
    for sub in ends.substitutions:
        built[sub.index] = ends.derived_codons[sub.index]
    assert "".join(built) == der
    assert "".join(ends.ancestral_codons) == anc


def test_a_column_neither_side_resolves_is_named_by_its_codons():
    rows = _rows(anc="NN-", der="---")
    ends = derive_endpoints(rows, "der", "anc")
    assert [s.name() for s in ends.substitutions] == ["c1:NN->---"]
    assert ends.class_counts() == {"unresolved": 1}


def test_deletion_is_an_indel_not_a_substitution_to_gap():
    rows = _rows(anc=ANC, der=DER)
    ends = derive_endpoints(rows, "der", "anc")
    deletion = ends.substitutions[-1]
    assert deletion.is_indel and deletion.derived == "-"
    assert ends.derived_codons[deletion.index] == GAP_CODON


def test_without_an_ancestor_states_come_from_column_plurality():
    # two rows carry S at column 2, one carries D; the derived row is excluded
    rows = _rows(der=DER, a="AAA" + "AGC" + "TGC", b="AAA" + "AGT" + "TGC",
                 c="AAA" + "GAC" + "TGC")
    ends = derive_endpoints(rows, "der")
    assert ends.reference_name == "der"
    # the deletion is named by its column: the derived row, which numbers
    # positions here, has no residue at it
    assert [s.name() for s in ends.substitutions] == ["S2P", "C@c3del"]
    # S in 2 of the 3 non-derived rows, and the codon is the commoner of AGC/AGT
    assert ends.plurality_support[1] == pytest.approx(2 / 3)
    assert residue_of(ends.ancestral_codons[1]) == "S"


def test_plurality_prefers_a_residue_over_a_gap_on_a_tie():
    rows = _rows(der="AAA", a="AGC", b="---")
    ends = derive_endpoints(rows, "der")
    assert residue_of(ends.ancestral_codons[0]) == "S"


def test_plurality_needs_another_row():
    with pytest.raises(SystemExit, match="at least one row"):
        derive_endpoints(_rows(der=DER), "der")


def test_columns_restrict_which_differences_count():
    rows = _rows(anc=ANC, der=DER)
    ends = derive_endpoints(rows, "der", "anc", columns=[1])
    assert [s.name() for s in ends.substitutions] == ["S2P"]


def test_positions_are_residues_of_the_reference_row():
    # a leading gap in the ancestor row shifts residue numbers off the columns
    rows = _rows(anc="---" + "AGC" + "TGC", der="AAA" + "CCC" + "TGC")
    ends = derive_endpoints(rows, "der", "anc")
    assert ends.residue_number == {1: 1, 2: 2}
    assert resolve_positions(["1"], ends) == [1]          # S at column 1
    assert resolve_positions(["c1"], ends) == [0]          # the gap column
    with pytest.raises(SystemExit, match="beyond the 2 residues"):
        resolve_positions(["3"], ends)


def test_insertion_at_a_column_the_reference_lacks():
    rows = _rows(anc="---", der="AAA")
    ends = derive_endpoints(rows, "der", "anc")
    assert [s.name() for s in ends.substitutions] == ["ins@c1K"]


def test_row_lookup_is_exact_then_loose():
    seqs = {"hg38": ANC, "HG38_alt": DER, "mm39": OTHER}
    assert resolve_row("hg38", seqs, "--x") == "hg38"
    assert resolve_row("mm", seqs, "--x") == "mm39"
    with pytest.raises(SystemExit, match="ambiguous"):
        resolve_row("hg", seqs, "--x")
    with pytest.raises(SystemExit, match="matches no sequence"):
        resolve_row("rn7", seqs, "--x")


def test_an_endpoint_need_not_be_an_analysed_taxon():
    # a row dropped by pruning or --max-species still supplies its codons
    seqs = {"kept": ANC, "dropped": DER}
    assert resolve_row("dropped", seqs, "--x") == "dropped"


def test_only_the_endpoint_taxa_are_read_from_their_own_rows():
    """An intermediate matching some other taxon is still a hypothetical genotype.

    Collapsing it onto that taxon would score a real sequence in place of the
    genotype on this branch, and would embed it among N taxa while its
    neighbours in the lattice are embedded among N + 1.
    """
    seqs = {"anc": "AAA" + "AGC", "der": "AAG" + "AGT", "bystander": "AAG" + "AGC"}
    rows = {name: codons_of(s, 2) for name, s in seqs.items()}
    ends = derive_endpoints(rows, "der", "anc")
    lat = Lattice(ends.ancestral_protein, ends.substitutions)

    class Ctx:                                  # only what real_row touches
        taxa = ["anc", "der", "bystander"]
        ancestor_name, derived_name = "anc", "der"
        derived_index = 1
    Ctx.seqs = {name: "".join(cs) for name, cs in rows.items()}

    g = GenotypeRows.__new__(GenotypeRows)
    g.ctx, g.ends, g.lattice, g.use_real_row = Ctx(), ends, lat, True
    g.host = Ctx.derived_index

    # only the host's own sequence leaves the alignment unmodified
    assert g.real_row(lat.full_mask) == 1
    # the ancestral genotype sits in the host's row, so the alignment differs
    # from the one the ancestor taxon belongs to
    assert g.real_row(0) is None
    # and a genotype matching some other taxon is still a hypothetical sequence
    assert g.nucleotides(1) == Ctx.seqs["bystander"]
    assert g.real_row(1) is None


