# Extant boreoeutherian SWS1 codon alignment

This is a **rebuilt** alignment of public nucleotide accessions cited by
[Yokoyama et al. (2014), *Epistatic Adaptive Evolution of Human Color Vision*](https://doi.org/10.1371/journal.pgen.1004884),
Materials and Methods, “Reconstruction of the ancestral SWS1 pigment.” It is
not the unavailable original author alignment. It contains extant sequences
only; no ancestral nucleotide sequence or ancestral codon has been invented.

`SWS1.codonified.fa` has **15 taxa × 350 codons (1,050 nucleotides)**. The human
row is `Homo_sapiens`. Protein homology was aligned with MAFFT G-INS-i, then
backtranslated by consuming the original CDS codons exactly once. All gaps are
whole `---` codons. Terminal stop codons are excluded from the alignment.

## Taxa and accession identity

The scope is the boreoeutherian union of the two Methods accession lists:
nine qualifying taxa from the 33-SWS1 reconstruction panel, plus six additional
primates. Elephant (Afrotheria), wallaby, dunnart, nonmammalian SWS1 sequences
and the RH1/RH2/SWS2 outgroups are excluded. This union contains ten primates
and five other boreoeutherians; it is not asserted to be the exact Figure 4
panel, which the text describes as nine primates and six other mammals.

| FASTA identifier | Paper accession | Extracted CDS record | Amino acids |
| --- | --- | --- | ---: |
| Homo_sapiens | M13295 | M13299.1, referencing M13295.1–M13298.1 | 348 |
| Pan_troglodytes | AF039433 | AF039435.1, referencing AF039433.1–AF039434.1 | 348 |
| Gorilla_gorilla | XM_004046176 | XM_004046176.1 | 348 |
| Pongo_abelii | XM_002818421 | XM_002818421.1 | 348 |
| Nomascus_leucogenys | XM_003261297 | XM_003261297.1 | 348 |
| Macaca_fascicularis | AF158977 | AF158977.1 | 349 |
| Saimiri_boliviensis | U53875 | U53875.1 | 349 |
| Callithrix_jacchus | L76201 | L76201.1, unannotated CDS fragment | 343 |
| Cephalopachus_bancanus | AB111463 | AB111463.1 | 346 |
| Daubentonia_madagascariensis | EF667285 | EF667285.1 | 348 |
| Bos_taurus | U92557 | U92557.1 | 349 |
| Mus_musculus | U49720 | U49720.1 | 346 |
| Rattus_norvegicus | U63972 | U63972.1 | 346 |
| Sciurus_carolinensis | DQ302163 | DQ302163.1 | 346 |
| Cavia_porcellus | AY552608 | AY552608.1 | 346 |

The paper's human and chimp accessions identify exon 1 records, rather than
full transcripts. Their companion records explicitly annotate remote CDS joins:
M13299.1 joins the original human five exons; AF039435.1 joins the original
chimp five exons. Those annotated joins were followed exactly. The three
RefSeq predictions use their historical **.1** versions (record dates 2012–2013),
rather than a current replacement transcript. Human's historical source
protein is 348 residues; current UniProt P03999 is 345 residues and must not
be silently substituted.

Two names differ between the paper and source records. The paper calls
U53875 *Saimiri sciureus*, but the record identifies *Saimiri boliviensis*;
the FASTA follows the accession's organism identity. *Tarsius bancanus* is
now named *Cephalopachus bancanus*. Both names are retained in `taxa.json`.

## Source limitations retained without imputation

- **L76201 (marmoset):** the 1,030-base record has no CDS feature. Frame 0
  starts with ATG, yields 343 complete codons without an internal stop, and
  matches the other SWS1 proteins along its length. Its last single T cannot
  define a codon and is excluded; the remaining C-terminal codons are missing
  and align as gaps. The nucleotide is retained in `raw/L76201.embl`. No
  completion or replacement sequence was supplied.
- **U63972 (rat):** the source CDS has one ambiguous nucleotide W, in codon 5
  (`GWC`, translating to X). The alignment retains W and X; no base is guessed.
- **U53875 (squirrel monkey):** its record title says “partial cds,” although
  the annotated CDS has an initiating ATG, terminal stop, and no internal stop
  or ambiguous base. This annotation discrepancy is retained in metadata.

The original paper uses **bovine rod-opsin homolog numbering** for the seven
substitutions (primary convention: Shi et al. 2001, Fig. 5). Source labels are
not sequential human SWS1 residues or alignment columns. The authoritative
seven-site crosswalk is in [../../numbering_offsets.json](../../numbering_offsets.json),
`sws1_benchmark.site_map`:

| Published event | Historical human CDS residue | Codon alignment column | 1U19 A homolog position |
| --- | ---: | ---: | ---: |
| F46T | 43 | 45 | 46 |
| F49L | 46 | 48 | 49 |
| T52F | 49 | 51 | 52 |
| F86L | 83 | 85 | 86 |
| T93P | 90 | 92 | 93 |
| A114G | 111 | 113 | 114 |
| S118T | 115 | 117 | 118 |

All seven sites are present as unambiguous codons in every row. The template
positions identify homologous locations, not experimentally measured SWS1
contacts. `human_coordinate_map.tsv` additionally maps all alignment columns
to the historical human sequence.

## Files and reproduction

- `raw/`: original EMBL/GenBank records, including companion exon records.
- `manifest.json`: exact source URLs, record versions, SHA256 input/output
  hashes, retrieval date, scope, software versions and aligner arguments.
- `taxa.json`: accession, organism, complete feature locations/qualifiers,
  length, trimming and ambiguity metadata for every row.
- `SWS1.source_cds.fa`: extracted nucleotide CDS, retaining annotated terminal
  stops where available; marmoset contains only its 343 complete codons.
- `SWS1.cds.fa` and `SWS1.protein.fa`: unaligned coding sequences and proteins,
  with terminal stops removed.
- `SWS1.protein.aligned.fa` and `SWS1.codonified.fa`: matched protein and codon
  alignments.
- `quality_report.json`: dimensions, translation, gap and completeness checks.
- `human_coordinate_map.tsv`: 1-based codon columns, nucleotide starts and
  ungapped historical human residue coordinates.

Requirements: Python with **Biopython 1.85**, and **MAFFT v7.526
(2024/Apr/26)**. The command is pinned to single-threaded G-INS-i:

```sh
mafft --quiet --amino --globalpair --maxiterate 1000 --thread 1 --inputorder SWS1.protein.fa
```

From the repository root, normal rebuilding and verification are offline:

```sh
python3 epistaeon/scripts/rebuild_sws1_alignment.py
python3 epistaeon/scripts/rebuild_sws1_alignment.py --check
```

`--download-missing` can recover absent raw records from their pinned public
URLs. Downloaded bytes must match their stored hashes; changed annotations
are rejected rather than silently replacing saved records. No phylogenetic
tree or ancestral reconstruction is needed to generate this alignment.

The integrity checks verify all 15 extant rows, original CDS translations,
absence of internal stop codons, exact codon backtranslation, whole-codon
gaps, expected source hashes and an exact offline rebuild of all products.
These checks do not run or evaluate the `epistaeon` model.

For selecting the seven focal units in an extant-only run with
`--focal-taxon Homo_sapiens`, use `--sites sites_human.txt`; `sites_columns.txt`
is the equivalent explicit `c<N>` codon-column list. These files do not supply
an ancestral row or impose an ancestral codon. Without `--ancestor`, the CLI
infers its own plurality states and must report that background rather than
calling it the published experimental ancestor. Keep all seven substitutions
separate: the historical stage sets denote unresolved chronology, not mutations
that occurred simultaneously.

All 15 nucleotide rows are unique, so duplicate pruning preserves 15 taxa;
human and chimp proteins are identical but their source codons differ. The
model tokenizer maps the retained ambiguous rat codon to its unknown token.
The source-preserving alignment itself does not require PyTorch; a full model
runtime does.
