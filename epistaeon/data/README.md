# Empirical data and coordinate inventory

> ### ⚠️ Numbering offsets — apply before any structural analysis
>
> `uniprot_position = pdb_author_position + offset`. **2HHB chains A/B/C/D and
> 1MBO chain A carry a +1 offset** (cleaved initiator Met); 1U19 is 1:1 to
> its own bovine rod-opsin sequence. SWS1 CDS/alignment positions require the
> explicit seven-site crosswalk, whose source labels use bovine opsin numbering. The
> ancestral receptors use local LBD numbering, `+531` to reach human GR
> NR3C1 through source position 210; downstream use the explicit crosswalk.
> Deleted S212 has no modern residue, and the 211/212 boundary is ambiguous.
> Machine-readable source of truth:
> [`numbering_offsets.json`](numbering_offsets.json). Details and rationale in
> [`../README.md`](../README.md). Getting this wrong fails silently.

These are published source files collected for reviewing
`white_paper_mutational_order_timing_epistasis.pdf`. Files are grouped by
empirical case, not by a claim that they constitute complete mutational-path
ground truth. No `axomeme_site_tensors.h5` was found under the parent directory
when this inventory was prepared (2026-09-15).

**Source review, 2026-10-15.** The verified substitution/site lists, partial
orders, structural locations and remaining validation gaps are now in
[`substitution_reference/README.md`](substitution_reference/README.md), with
per-case JSON and a flat [`sites.tsv`](substitution_reference/sites.tsv).
This extends the older `validation/studies/` lists after checking the original
papers. Case 2 uses the modern `1MBO` structure for residue-transition locations.
Case 4 now uses Yokoyama 2014 SWS1 cone opsin and a separately rebuilt
extant boreoeutherian codon alignment in [`alignments/sws1/`](alignments/sws1/).
The older structural audit below is not itself a source of evolutionary order.

## Experimental sources

| Case | Local files | Published source and scope |
| --- | --- | --- |
| Hemoglobin | `experimental/hemoglobin_2009_storz.pdf` | Storz et al., *PNAS* 2009, [doi:10.1073/pnas.0905224106](https://doi.org/10.1073/pnas.0905224106). High-altitude hemoglobin functional data; not a complete ancestral genotype panel. |
| Myoglobin | `experimental/myoglobin_2013_article.pdf`, `myoglobin_2013_supplement.pdf` | Mirceta et al., *Science* 2013, [doi:10.1126/science.1234192](https://doi.org/10.1126/science.1234192). Supplement includes species traits, myoglobin concentration/charge data, and ancestral inference; it does not supply a complete measured recombinant ancestral-mutation panel. |
| Steroid receptor | `experimental/steroid_receptor_2007_ortlund.pdf`, `steroid_receptor_2009_bridgham.pdf` | Ortlund et al., *Science* 2007, [doi:10.1126/science.1142819](https://doi.org/10.1126/science.1142819); Bridgham et al., *Nature* 2009, [doi:10.1038/nature08249](https://doi.org/10.1038/nature08249). Reconstructed receptor mutagenesis, ligand specificity, and epistasis. |
| Steroid receptor | `experimental/steroid_receptor_2011_plos.pdf`, `steroid_receptor_2011_ec50_data.doc`, `steroid_receptor_2011_sequence_accessions.doc`, `steroid_receptor_2011_ancestral_sequences.doc`, `steroid_receptor_2011_stability_data.doc` | Carroll et al., *PLOS Genetics* 2011, [doi:10.1371/journal.pgen.1002117](https://doi.org/10.1371/journal.pgen.1002117), [S1 EC50 data](https://doi.org/10.1371/journal.pgen.1002117.s003), [S2 accession list](https://doi.org/10.1371/journal.pgen.1002117.s004), [S3 ancestral sequences](https://doi.org/10.1371/journal.pgen.1002117.s005), [S5 stability calculations](https://doi.org/10.1371/journal.pgen.1002117.s007). The stability table is calculated, not measured thermal unfolding. |
| SWS1 cone opsin | `experimental/sws1_2014_yokoyama.pdf` and HTML; `sws1_2014_S1_Table.pdf`, `sws1_2014_S2_Table.pdf`, `sws1_2014_S3_Table.pdf`; `sws1_2014_S1_Fig.eps` through `sws1_2014_S4_Fig.eps` | Yokoyama et al., *PLOS Genetics* 2014, [doi:10.1371/journal.pgen.1004884](https://doi.org/10.1371/journal.pgen.1004884). All 127 mutant combinations of seven substitutions were assayed for pigment formation and absorption spectra. The paper reports 60 nonfunctional mutants and 1,032 accessible paths; phylogenetic/parsimony inference favors eight historical paths. No dark-noise or thermal-decay assay. |

The *Science* 2013 hemoglobin combinatorial experiment of Natarajan et al.,
[doi:10.1126/science.1236862](https://doi.org/10.1126/science.1236862),
assayed eight combinations of three allele blocks and additional single/double
variants. Its publication is available online, but its supporting genotype
data were not obtainable as a local file from the public endpoints tested.

## SWS1 alignment and ordering input

[`alignments/sws1/SWS1.codonified.fa`](alignments/sws1/SWS1.codonified.fa) is
a new codon-aware alignment of **extant boreoeutherian SWS1 sequences**, drawn
from the accession lists in Yokoyama 2014's Methods, including the six added
primates. Non-boreoeutherian mammals and the other vertebrates/outgroup opsins
are excluded. No ancestral sequence or anole-derived chimeric terminus is
inserted. This rebuild is not the authors' recovered alignment. Its accession
inventory and alignment provenance live beside the FASTA.

The seven events are **F46T, F49L, T52F, F86L, T93P, A114G and S118T**. The
historical partial order is
**`{T93P, A114G} -> F86L -> {F49L, S118T} -> {F46T, T52F}`**; each brace
is internally unordered, giving eight favored paths. This is distinct from
the reported 1,032 paths that avoid nonfunctional intermediates in the mutagenesis
panel. T52F-first is experimentally excluded on the assayed ancestral
background, but the entire historical partial order is not a biochemical
must-precede constraint. See [`substitution_reference/sws1.json`](substitution_reference/sws1.json).

**Accessibility counts require source reconciliation before exact scoring.**
The archived S1 table has a duplicated/missing genotype label and a
nonfunctionality-formatting anomaly. Enumeration from the normalized local
transcription gives **1,076** accessible paths, versus the paper's reported
**1,032** (and 336 versus 335 paths with every absolute step below 25 nm). The genotype/path artifacts retain this discrepancy; do not describe
the reported path count as independently reproduced or silently repair assay
labels to force agreement.

Paper site labels, historical human M13295 sequence positions, rebuilt
alignment columns and bovine 1U19 coordinates are separate. Load the SWS1
crosswalk from [`numbering_offsets.json`](numbering_offsets.json). Current
human UniProt P03999 has changed relative to M13295; a direct constant offset
is not a valid crosswalk. 1U19 is retained because the paper used it for
homology modeling, not because it is an experimentally determined SWS1
structure.

## Coordinate sources

All coordinate files are RCSB/wwPDB mmCIF downloads from
`https://files.rcsb.org/download/<PDB-ID>.cif`. Resolution values were checked
against each file's `_refine.ls_d_res_high` field.

| File | Structure | Resolution | Role |
| --- | --- | --- | --- |
| `coordinates/2HHB.cif` | Human deoxyhemoglobin | 1.74 Å | Section 4 accession |
| `coordinates/1MBO.cif` | Sperm-whale oxymyoglobin | 1.60 Å | Section 4 accession |
| `coordinates/1U19.cif` | Bovine rod opsin | 2.20 Å | Homology template used in the SWS1 paper; not a SWS1 structure |
| `coordinates/7PRX.cif` | Human GR LBD, wild type–velsecorat + PGC1α peptide | 2.20 Å | Active Case 3 structure |
| `alternate_coordinates/4LSJ.cif` | Human GR LBD–synthetic dibenzoxapine sulfonamide | 2.35 Å | Superseded: domain-swapped (see Case 3) |
| `ancestral_coordinates/2Q1H.cif` | Ancestral corticoid receptor–aldosterone | 1.90 Å | Backup |
| `ancestral_coordinates/2Q1V.cif` | AncCR–prednisone (entry titled "cortisol") | 1.95 Å | Backup |
| `ancestral_coordinates/2Q3Y.cif` | AncCR–DOC | 2.40 Å | Backup |
| `ancestral_coordinates/3RY9.cif` | AncGR1–DOC | 1.95 Å | Backup |
| `ancestral_coordinates/3GN8.cif` | AncGR2–dexamethasone | 2.50 Å | Backup |

The active set in `coordinates/` supplies structures for Cases 1–3 and the
bovine homolog template for Case 4, with `7PRX` in place of the receptor
entry. `4LSJ`, the likely intended §4 accession, is kept in
`alternate_coordinates/`; the Thornton ancestral series is kept in
`ancestral_coordinates/`. Both are backups.

These structures do not by themselves establish a tensor-to-coordinate
crosswalk. It requires the tensor's reference sequence and site-index
metadata, a sequence alignment to a chosen PDB polymer chain, and explicit
handling of alignment gaps and unresolved residues.

---

# Structure validation (2026-09-17)

Every coordinate file was re-parsed and checked against the residues and
quantities the white paper names. Scripts that reproduce every number below are
in `epistaeon/scripts/`. Distances are minimum heavy-atom unless stated.

## Resolution of the `4LS6` error

`4LS6` is *Bacillus subtilis* FabF and has no steroid receptor content. An RCSB
search for human glucocorticoid receptor entries returns **`4LSJ`**, one
character away, which is a genuine GR ligand-binding domain (2.35 Å,
*J. Med. Chem.* 2014) — so §4's accession is most likely a transcription slip.

The residue numbers §4 cites (Ser106, Leu111) are **AncCR local numbering**;
on a human GR structure they sit 531 positions higher. `4LSJ` itself turned out
to be domain-swapped, so the active Case 3 structure is the wild-type human GR
entry `7PRX` (see Case 3 below).

## Numbering crosswalk (established, not assumed)

Ancestral receptor entries number their LBD locally, `-2` to `247`. Human GR
entries use full-length NR3C1 numbering, `520`–`777`. Global pairwise alignment
of 3GN8 to 4P6X establishes the upstream offset (4P6X is no longer kept in the
repository; `scripts/_fetch.py` re-downloads it on demand so this stays
reproducible):

    human GR position = AncCR position + 531  (through source position 210)

    Anc  27 (Arg) == GR 558      Anc 106 (Pro) == GR 637
    Anc  36 (Gly) == GR 567      Anc 111 (Gln) == GR 642

Source positions 213 onward use +530. Deleted S212 has no modern homolog.
A generic alignment gaps Q211 instead of the source-identified Ser212, so the
211/212 boundary needs a curated, source-anchored alignment. Do not assign the
deleted serine to modern Asp742. The per-event reference records this explicitly.

## Coordinate-derived receptor ladder (not a historical event reference)

Raw alignment-based mismatch counts across the coordinate constructs:

| Transition | Substitutions | Contains |
| --- | --- | --- |
| AncCR (2Q1H) → AncGR1 (3RY9) | 38 (35 excl. 3 N-terminal tag residues) | Y27R, A36G |
| AncGR1 (3RY9) → AncGR2 (3GN8) | 39 | N26T, L29M, F98I, Q105L, **S106P**, **L111Q** |

These construct/reconstruction comparisons do not reproduce the original
historical event counts. Ortlund 2007 reports 25 substitutions in the early
interval and 36 substitutions plus one deletion in the later interval.
Carroll 2011 updates AncGR1 to AncGR1.1; keep its background distinct.

Residue identities confirm the epistatic ladder directly from coordinates:

| Position | AncCR 2Q1H | AncGR1 3RY9 | AncGR2 3GN8 |
| --- | --- | --- | --- |
| 27 | Tyr | Arg | Arg |
| 36 | **Ala** | Gly | Gly |
| 106 | Ser | Ser | **Pro** |
| 111 | Leu | Leu | **Gln** |

**`Thr36` in the white paper is wrong.** Position 36 is Ala in AncCR and Gly in
both derived receptors; Thr never occurs there. It is also 3.7 Å from
dexamethasone — pocket-lining, not a remote residue. Ortlund 2007 identifies
**S106P/L111Q** as group X, **N26T/Q105L** as permissive group Z, and
**L29M/F98I/S212del** as group Y; **Y27R** is an earlier permissive change.
S212del is a single-codon deletion, not S212D. `Thr36` is not a source target;
do not guess a replacement target from the white paper.

Measured pocket↔remote separations in 3GN8 (paper claims 15.2 Å):

    106 <-> 27   13.9 A Cb-Cb   11.6 A min-heavy
    111 <-> 27   14.5 A Cb-Cb    9.7 A min-heavy
    111 <-> 29    8.4 A Cb-Cb

Nothing reaches 15.2 Å, and none of these pairs exceeds the ">15 Å" the text
claims.

## §4 Table 2 claims, checked

| Claim | Measured | Verdict |
| --- | --- | --- |
| Hemoglobin α1β1, 55 contact pairs ≤5 Å | **55** | reproduces exactly |
| Hemoglobin α1β2, 37 contact pairs ≤5 Å | **37** | reproduces exactly |
| α-His87 ↔ β-site51, 18.4 Å | 18.4 Å for α1(A)↔β2(D); 31.3 Å for α1↔β1 | number is real but filed under the wrong interface row |
| Perutz α40Lys ↔ β146His | 2.5 Å (A↔D) | confirmed salt bridge |
| Steroid receptor remote network, 15.2 Å | 11.6–14.5 Å | not reproduced |

## Case-by-case sufficiency

**Case 1, hemoglobin — structure correct, site set in the paper is wrong.**
2HHB is clean (4 chains, 1–141/1–146, no gaps, 4 haems). But Storz 2009, our
local copy, defines the deer mouse system as **α sites 50, 57, 60, 64, 71** and
**β sites 62, 72, 128, 135**. These define isoforms, not exhaustive 2⁵/2⁴
mutant panels. Fig. 3 compares mixed-alpha hemolysates with the two beta
isoforms. The paper gives beta substitution polarity but no mutation order;
diagnostic sites are not individually demonstrated epistatic partners.
More seriously, the mechanism is
contradicted: α57 and α64 sit in the **haem pocket** (6.5–7.5 Å from haem),
15.3 Å and 18.9 Å from the nearest β1 residue. They cannot "disrupt" or
"tighten" the α1β1 packing face. The Table 1 pair α114↔β116 *is* a real α1β1
contact (2.7 Å) but neither residue is in the Storz site set.

**Case 2, myoglobin — modern scaffold accepted for transition locations.**
1MBO is *Physeter catodon* (the target taxon), 1.6 Å, 153 residues, no gaps,
haem + bound O₂, His93–Fe 2.1 Å. Use modern 1MBO to locate the substitutions
and protein regions. Ancestral coordinates are not required for this analysis.
Mirceta Fig. S2 supplies 199 labeled charge-changing events in 85 branch sets;
the reference maps them all onto 1MBO. Branch ancestry establishes a partial
order, with no within-branch order or measured epistatic mutant panel.

**Case 3, steroid receptor — 7PRX active; 4LSJ and the ancestral series as
backups.** 7PRX was chosen after checking all 28 human GR ligand-binding-domain
crystal structures against UniProt P04150. It is one of very few with the
wild-type sequence (most carry solubilizing mutations such as F602S/C638D),
2.20 Å, residues 529–776 fully modelled with no internal gaps, no incomplete
side chains, and a PGC1α coactivator peptide bound (active conformation). Its
contact calls agree closely with dexamethasone-bound 1M2Z and 4UDD.

- Numbering is full-length NR3C1 and matches UniProt exactly (offset 0). The
  white paper's residue names are AncCR numbering, so add 531: Tyr27 → 558,
  Leu29 → 560, Ala36 → 567, Ser106 → 637, Leu111 → 642. Human GR already
  carries the derived residues (Arg558, Pro637, Gln642).
- The ligand is velsecorat, a non-steroidal modulator. The contact rule uses
  only residue–residue distances; 7PRX and 1M2Z superimpose to 0.73 Å Cα RMSD.
  If a steroid-bound pocket is preferred, 1M2Z (dexamethasone, F602S only,
  2.50 Å) is the runner-up.

**Why not 4LSJ.** It is the human GR entry §4 most likely meant by `4LS6`, but
its crystal is **domain-swapped**: residues 526–551 of each chain pack onto the
neighbouring copy's core. Within-chain contacts there are wrong — GR 558
(Tyr27) has 1 contact partner in 4LSJ versus 6 in every other GR structure,
and 4LSJ calls about 80 fewer contacts overall. It also lacks residues 703–710
and carries F602Y and C638G.

The Thornton ancestral structures (`2Q1H`, `2Q1V`, `2Q3Y`, `3RY9`, `3GN8`) are
kept in `ancestral_coordinates/` as a backup. Note `2Q1V` is titled "in
complex with cortisol" but its deposited ligand is `PDN` = prednisone, and
`3GN8` numbering runs one lower than AncCR from about position 212 onward.

**Case 4, SWS1 — codon-aware extant input and source-qualified order.** The
rebuilt boreoeutherian alignment supplies the requested extant sequences;
the reference records seven source-labeled substitutions, the eight
phylogenetically favored paths and the separate experimental accessibility
criterion. It does not require ancestral sequences as comparative input.
An endpoint-based `order` run must separately declare how its endpoints match
the assayed background; column plurality is not the experimental ancestor.

1U19 provides a homolog template for locating mapped sites, following the
original study's modeling approach. Its residue identities and contacts
describe bovine rod opsin, not SWS1. Direct structural validation of the SWS1
mutants would require a declared SWS1 model or experimental structure; a
template contact is not itself evidence of epistasis or temporal order.

## Recommended coordinate set

| Case | Use |
| --- | --- |
| 1 | `2HHB` (only; re-derive the site set from Storz 2009) |
| 2 | `1MBO` |
| 3 | `7PRX` (backups: `4LSJ` in `alternate_coordinates/`, ancestral series in `ancestral_coordinates/`) |
| 4 | `1U19` as a homolog template only, using the SWS1 source-site crosswalk |


---

# Contact-mapping fitness of 2HHB / 1MBO / 1U19 / 7PRX

Validated against a narrower goal than the §4 accession check: given an
arbitrary epistatic site pair emitted by the model in a §3 case, can this
structure say reliably whether the two residues are in direct contact?
Reproduce with `epistaeon/scripts/numbering.py`, `mapfitness.py`, `bands.py`.
Case 3 uses `7PRX`; `4LSJ` and the ancestral series are backups. Case 4 uses
1U19 only as a mapped homolog template. Its baseline and power calculations
below describe all bovine template residues, not a recalibrated SWS1 site set.
`4P6X`, `1M2Z` and `4LS6` are not in the repository.

## Passes

**Numbering is exact and unambiguous.** Aligned to UniProt canonical sequences,
zero mismatches in the globin, myoglobin and bovine template chains:

| Chain | UniProt | Offset | Coverage |
| --- | --- | --- | --- |
| 2HHB A (α) | P69905 | UniProt = PDB **+1** | 141/141 modelled |
| 2HHB B (β) | P68871 | UniProt = PDB **+1** | 146/146 modelled |
| 1MBO A | P02185 | UniProt = PDB **+1** | 153/153 modelled |
| 1U19 A | P02699 | UniProt = PDB **+0** | 348/348 modelled |
| 7PRX A | P04150 | UniProt = PDB **+0** | 248/250 modelled (termini 528, 777 absent) |

The +1 offsets are initiator-Met cleavage in the mature protein. They must be
applied explicitly; an unconverted alignment index silently shifts every
globin and myoglobin contact call by one residue. 7PRX matches UniProt with
zero mismatches.

**Missing coordinates.** Every SEQRES position is modelled in 2HHB, 1MBO and
1U19. 7PRX lacks only its terminal residues 528 and 777, so GR 529–776 is
fully mappable. Zero
incomplete side chains in all 574 + 153 + 696 residues, so minimum heavy-atom
distances are trustworthy. (1MBO has 4 partial-occupancy/altloc residues;
1U19 carries an `ACE` N-terminal cap at position 0 that must be skipped when
iterating residues.)

**Class assignment is robust to crystallographic copy choice.**

| Comparison | median &#124;Δd&#124; | class disagreement | contact↔allosteric flips |
| --- | --- | --- | --- |
| 2HHB α copies A vs C | 0.10 Å | 0.73% | **0** of 9,316 |
| 2HHB β copies B vs D | 0.12 Å | 0.71% | **0** of 10,011 |
| 1U19 copies A vs B | 0.20 Å | 0.77% | **1** of 58,996 |

## Classification rule

The question is binary: of the epistatic pairs the model detects, how many are
in direct physical contact and how many are not?

**A pair is in direct contact if any heavy atom of one residue lies within
5.0 Å of any heavy atom of the other** (the white paper's own cutoff). Every
other pair is *not in direct contact*. There are no further classes.

Three details change the count and must be applied:

1. **Numbering offset.** Convert alignment positions with
   [`numbering_offsets.json`](numbering_offsets.json) before any lookup
   (+1 for 2HHB and 1MBO, 0 for 1U19 and 7PRX). Skipping it fails silently.
   For Case 4, apply the explicit SWS1 source-site/homolog crosswalk; the 1U19
   offset alone does not map SWS1. For Case 3, source residue names are AncCR numbering: apply the piecewise
   crosswalk and deletion exceptions; do not assign S212 a modern coordinate.
2. **Sequence separation.** Exclude pairs with |i − j| < 5. Chain neighbours
   are always within 5 Å, so they would count as contacts for free.
3. **Hemoglobin α–β pairs.** Count a pair as a contact if it touches at
   *either* interface of the tetramer (α1β1 or α1β2). The two interfaces share
   no contact pairs, so measuring on only one would miss half the true
   contacts. 2HHB is the deoxy (T) state; about half the α1β2 contacts differ
   in the oxy (R) state (`scripts/interface_states.py`).

**Baseline.** A contact fraction only means something against the rate among
all pairs (|i − j| ≥ 5):

| Structure | pairs in direct contact |
| --- | --- |
| 2HHB α | 2.3% |
| 2HHB β | 2.4% |
| 1MBO | 2.2% |
| 1U19 | 1.3% |
| 7PRX | 1.7% |

**Interpretation note.** Residues that couple through the haem or retinal
without touching each other are counted as not in direct contact. That is
correct under this rule, but some "not in contact" pairs are still
mechanistically close.

**Still open: cross-species transfer.** 1MBO is *Physeter catodon*, the Case 2
target taxon. But 2HHB is human while Case 1 targets *Peromyscus* and
other rodents, and 1U19 is bovine rod opsin while Case 4 concerns SWS1 cone
opsin. The template's own numbering check cannot validate SWS1 positions;
use the source-anchored homolog mapping. Structural-transfer uncertainty
remains even when the sequence mapping is explicit.

## Power requirement for §5

At the measured background contact rates, "Odds Ratio > 4.0, Fisher's exact
p < 10⁻⁵" requires a minimum number of detected pairs:

| Structure | background ≤5 Å rate | OR=4 implies | pairs needed |
| --- | --- | --- | --- |
| 2HHB α | 2.30% | 8.6% contacts | ≥ 181 |
| 2HHB β | 2.41% | 9.0% contacts | ≥ 173 |
| 1MBO | 2.22% | 8.3% contacts | ≥ 187 |
| 1U19 | 1.27% | 4.9% contacts | ≥ 317 |
| 7PRX | 1.65% | 6.3% contacts | ≥ 247 |

A detector emitting fewer pairs than this cannot meet the §5 standard on that
case regardless of accuracy.
