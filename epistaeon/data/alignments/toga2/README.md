# TOGA2 codon alignments (model input)

Downloaded 2026-09-18 from
`https://genome.senckenberg.de/download/TOGA2/MultiCodonAlignments/reference_human_hg38/v1/`
(`individualAlis/`, plus `speciesTree.nh`). `assemblies_and_species.tsv` is from
`https://genome.senckenberg.de/download/TOGA2/` and maps sequence names
(e.g. `HLmirAng2`) to species. Human-referenced, so **mammals only**.
Per the upstream README, only species with exactly one ortholog are included;
`N` marks masked (frameshift-affected) codons.

| File | Human reference transcript | Seqs | Columns (codons) | Human residues | Taxa after loader |
| --- | --- | --- | --- | --- | --- |
| `HBA1.codonified.fa` | ENST00000320868.9 | 149 | 174 | 142 | 130 |
| `HBA2.codonified.fa` | ENST00000251595.11 | 212 | 205 | 142 (+1 masked trailing codon) | 191 |
| `HBB.codonified.fa` | ENST00000647020.1 | 422 | 153 | 147 | 323 |
| `MB.codonified.fa` | NM_203377.1 | 863 | 186 | 154 | 665 |
| `NR3C1.codonified.fa` | NM_001364183.2 | 870 | 866 | **778** | 737 |
| `NR3C2.codonified.fa` | XM_047415708.1 | 867 | 1094 | 984 | 792 |
| `RHO.codonified.fa` | ENST00000296271.4 | 864 | 356 | 348 | 724 |

Every sequence name is a tip in `speciesTree.nh`; all seven load with
`aeon_core.dataset.load_alignment_and_tree(fa, "speciesTree.nh")`, which
prunes the tree to the alignment and collapses identical sequences.

## ⚠️ Position mapping — two steps, two traps

The model reports **alignment columns**, not residue numbers. To reach PDB
coordinates:

1. **Column → human residue.** Drop columns where `hg38` has a gap. Every
   alignment has such columns (HBB 6, RHO 8, MB 32, HBA1 32, HBA2 62,
   NR3C1 88), so column index ≠ residue number.
2. **Human residue → structure numbering.** Then apply
   `../../numbering_offsets.json`.

**NR3C1 is the GRγ isoform.** Its reference carries an extra Arg at position
452 (DNA-binding domain), so every human-reference position **after 451 is one
higher than UniProt P04150 / 7PRX numbering**: GR 637 is reference residue
638. Subtract 1 before applying the 7PRX offset. Fails silently if skipped.

MB and RHO: the human reference is the same length as the structure's
UniProt entry (sperm whale P02185, bovine P02699) with no indels, so human
residue *n* is structure-UniProt residue *n*; the sequences differ at 24 and
23 sites respectively.

## Known issues

Everything found while validating these files (2026-09-18). Items 1–2 change
results silently if ignored.

**Position mapping**

1. **Model positions are alignment columns, not residues.** Columns where
   `hg38` has a gap (insertions in other species): HBB 6, RHO 8, MB 32,
   HBA1 32, HBA2 62, NR3C1 88, NR3C2 110. Drop them to get human residue
   numbers before applying `../../numbering_offsets.json`.
2. **NR3C1 reference is the GRγ isoform** (NM_001364183.2, 778 aa). It has an
   extra Arg at 452, so human positions after 451 are one higher than UniProt
   P04150 / 7PRX. Subtract 1 there (GR 637 = reference 638).
3. **HBA2 reference ends in one masked codon** (translates to `X` at 143).
   Harmless: it falls after the last residue (142).
4. **MB and RHO references are human, the structures are not** (sperm whale
   1MBO, bovine 1U19). Neither has an indel relative to human, so the
   numbering schemes match: human residue *n* = 1U19 residue *n*, and = 1MBO
   residue *n* − 1 (the usual cleaved-Met offset). The sequences differ at 24
   (MB) and 23 (RHO) sites. Confirmed against human UniProt P02144 / P08100
   and inside the alignments (the sperm whale and cow rows have no indels
   relative to `hg38`); recorded under `human_equivalent` in
   `../../numbering_offsets.json` and reproduced by `scripts/numbering.py`.

**Taxon coverage**

5. **Mammals only.** All 3,847 sequences and all 882 species-tree tips are
   mammals; this is the human-referenced TOGA2 set. (The species table also
   lists birds, turtles and fish because it covers TOGA2's other references.)
6. **Deer mouse is missing from every globin alignment.** TOGA2 has three
   *Peromyscus maniculatus* assemblies, but none has an HBA1, HBA2 or HBB
   sequence here. Likely cause (inferred, not confirmed): TOGA2 keeps only
   species with exactly one ortholog, and deer mice carry multiple α- and
   β-globin copies. This removes Case 1's focal species.
7. **Andean waterfowl (*Merganetta armata*) is not in TOGA2 at all** — no
   assembly under any reference, not just this mammal set. Its sequences would
   have to come from elsewhere, e.g. the GenBank accessions behind Natarajan
   et al. 2015 (`../../experimental/`).
8. **Other white-paper taxa:** yak and snow leopard have HBB only (no HBA1 or
   HBA2). Sperm whale, southern elephant seal and Weddell seal are in MB.
   Naked mole-rat, star-nosed mole, *Myotis nattereri* and sperm whale are in
   RHO (with 75 cetaceans and 113 bats); *Rhinolophus unihastatus* has no TOGA2
   assembly, and no fish are present.
9. **Globins are thinly sampled.** HBA1 has 149, HBA2 212 and HBB 422 of the
   882 species, versus 863–870 for MB, NR3C1, NR3C2 and RHO — consistent with
   multi-copy species being excluded. The split between HBA1 and HBA2 is
   TOGA2's ortholog assignment; α-globin paralogs are hard to assign, so
   treat per-paralog species sets with caution.

**Data quality**

10. **HBA1 and HBA2 are sparse.** Share of codons that are gaps or masked
    (`N`): HBA1 34.5%, HBA2 44.0%, MB 19.5%, NR3C2 13.9%, NR3C1 13.5%,
    HBB 7.8%, RHO 5.0%. Expect weaker signal for Case 1.
11. **Internal stop codons:** grey seal (`HLhalGry1`) has 2 in NR3C1, a
    possible frameshift or pseudogene.
12. **The loader collapses identical sequences** before inference (HBA1 19,
    HBA2 21, HBB 99, MB 198, NR3C1 133, NR3C2 75, RHO 140). This doesn't affect
    site positions, but a named species can be merged into a duplicate.

**Provenance**

13. The upstream README points to `../assemblies_and_species.tsv`, which
    returns 404; the file is at the TOGA2 top level (URL above).
14. The white paper says the HBA/HBB loci span "115–349 vertebrate species";
    these alignments have 149 (HBA1), 212 (HBA2) and 422 (HBB) mammals.
