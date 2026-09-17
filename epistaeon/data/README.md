# Empirical data and coordinate inventory

These are published source files collected for reviewing
`white_paper_mutational_order_timing_epistasis.pdf`. Files are grouped by
empirical case, not by a claim that they constitute complete mutational-path
ground truth. No `axomeme_site_tensors.h5` was found under the parent directory
when this inventory was prepared (2026-09-15).

## Experimental sources

| Case | Local files | Published source and scope |
| --- | --- | --- |
| Hemoglobin | `experimental/hemoglobin_2009_storz.pdf` | Storz et al., *PNAS* 2009, [doi:10.1073/pnas.0905224106](https://doi.org/10.1073/pnas.0905224106). High-altitude hemoglobin functional data; not a complete ancestral genotype panel. |
| Hemoglobin | `experimental/hemoglobin_2015_andean_waterfowl.pdf`, `hemoglobin_2015_andean_waterfowl_oxygen_data.docx`, `hemoglobin_2015_andean_waterfowl_isoform_data.docx`, `hemoglobin_2015_andean_waterfowl_alignment.pdf` | Natarajan et al., *PLOS Genetics* 2015, [doi:10.1371/journal.pgen.1005681](https://doi.org/10.1371/journal.pgen.1005681), [S2 oxygen-affinity table](https://doi.org/10.1371/journal.pgen.1005681.s013), [S1 isoform table](https://doi.org/10.1371/journal.pgen.1005681.s012), [S2 alignment figure](https://doi.org/10.1371/journal.pgen.1005681.s002). The oxygen-affinity data use KCl and IHP conditions; they are not 2,3-DPG assays. |
| Myoglobin | `experimental/myoglobin_2013_article.pdf`, `myoglobin_2013_supplement.pdf` | Mirceta et al., *Science* 2013, [doi:10.1126/science.1234192](https://doi.org/10.1126/science.1234192). Supplement includes species traits, myoglobin concentration/charge data, and ancestral inference; it does not supply a complete measured recombinant ancestral-mutation panel. |
| Steroid receptor | `experimental/steroid_receptor_2007_ortlund.pdf`, `steroid_receptor_2009_bridgham.pdf` | Ortlund et al., *Science* 2007, [doi:10.1126/science.1142819](https://doi.org/10.1126/science.1142819); Bridgham et al., *Nature* 2009, [doi:10.1038/nature08249](https://doi.org/10.1038/nature08249). Reconstructed receptor mutagenesis, ligand specificity, and epistasis. |
| Steroid receptor | `experimental/steroid_receptor_2011_plos.pdf`, `steroid_receptor_2011_ec50_data.doc`, `steroid_receptor_2011_sequence_accessions.doc`, `steroid_receptor_2011_ancestral_sequences.doc`, `steroid_receptor_2011_stability_data.doc` | Carroll et al., *PLOS Genetics* 2011, [doi:10.1371/journal.pgen.1002117](https://doi.org/10.1371/journal.pgen.1002117), [S1 EC50 data](https://doi.org/10.1371/journal.pgen.1002117.s003), [S2 accession list](https://doi.org/10.1371/journal.pgen.1002117.s004), [S3 ancestral sequences](https://doi.org/10.1371/journal.pgen.1002117.s005), [S5 stability calculations](https://doi.org/10.1371/journal.pgen.1002117.s007). The stability table is calculated, not measured thermal unfolding. |
| Rhodopsin | `experimental/rhodopsin_2017_cetacean_supplement.pdf` | Dungan & Chang, *Proceedings of the Royal Society B* 2017, [doi:10.1098/rspb.2016.2743](https://doi.org/10.1098/rspb.2016.2743), [published supplement](https://doi.org/10.6084/m9.figshare.4653931.v1). Includes mutant spectral and retinal-release results; retinal release after photoactivation is not a direct dark-noise rate. The [Dryad raw archive](https://doi.org/10.5061/dryad.5k0s6) was located but returned HTTP 401 on download. |

The *Science* 2013 hemoglobin combinatorial experiment of Natarajan et al.,
[doi:10.1126/science.1236862](https://doi.org/10.1126/science.1236862),
assayed eight combinations of three allele blocks and additional single/double
variants. Its publication is available online, but its supporting genotype
data were not obtainable as a local file from the public endpoints tested.
Likewise, the 2008 ancestral-rhodopsin PNAS study of Yokoyama et al.,
[doi:10.1073/pnas.0802426105](https://doi.org/10.1073/pnas.0802426105),
was located online, but its publisher supplement could not be downloaded.

## Coordinate sources

All coordinate files are RCSB/wwPDB mmCIF downloads from
`https://files.rcsb.org/download/<PDB-ID>.cif`. Resolution values were checked
against each file's `_refine.ls_d_res_high` field.

| File | Structure | Resolution | Role |
| --- | --- | --- | --- |
| `coordinates/2HHB.cif` | Human deoxyhemoglobin | 1.74 Å | Section 4 accession |
| `coordinates/1MBO.cif` | Sperm-whale oxymyoglobin | 1.60 Å | Section 4 accession |
| `coordinates/1U19.cif` | Bovine rhodopsin | 2.20 Å | Section 4 accession |
| `coordinates/2Q1H.cif` | Ancestral corticoid receptor–aldosterone | 1.90 Å | Relevant receptor replacement |
| `coordinates/3RY9.cif` | AncGR1–DOC | 1.95 Å | Relevant receptor replacement |
| `coordinates/3GN8.cif` | AncGR2–dexamethasone | 2.50 Å | Relevant receptor replacement |
| `coordinates/2Q1V.cif` | AncCR–prednisone (entry titled "cortisol") | 1.95 Å | Added 2026-09-17, completes Ortlund 2007 set |
| `coordinates/2Q3Y.cif` | AncCR–DOC | 2.40 Å | Added 2026-09-17, completes Ortlund 2007 set |

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

`4LSJ` is nevertheless still the wrong structure for these analyses: its ligand
(`LSJ`) is a synthetic dibenzoxapine sulfonamide, not a corticosteroid, and the
residue numbers §4 cites (Ser106, Leu111) are **AncCR local numbering**, which
no human GR entry uses. The §4 receptor row therefore belongs to the ancestral
structures, not to any human GR entry.

## Numbering crosswalk (established, not assumed)

Ancestral receptor entries number their LBD locally, `-2` to `247`. Human GR
entries use full-length NR3C1 numbering, `520`–`777`. Global pairwise alignment
of 3GN8 to 4P6X gives a constant offset (4P6X is no longer kept in the
repository; `scripts/_fetch.py` re-downloads it on demand so this stays
reproducible):

    human GR position = AncCR position + 531

    Anc  27 (Arg) == GR 558      Anc 106 (Pro) == GR 637
    Anc  36 (Gly) == GR 567      Anc 111 (Gln) == GR 642

## The steroid receptor ladder is complete and correct

Alignment-based substitution counts across the resurrected series:

| Transition | Substitutions | Contains |
| --- | --- | --- |
| AncCR (2Q1H) → AncGR1 (3RY9) | 38 (35 excl. 3 N-terminal tag residues) | Y27R, A36G |
| AncGR1 (3RY9) → AncGR2 (3GN8) | 39 | N26T, L29M, F98I, Q105L, **S106P**, **L111Q** |

Residue identities confirm the epistatic ladder directly from coordinates:

| Position | AncCR 2Q1H | AncGR1 3RY9 | AncGR2 3GN8 |
| --- | --- | --- | --- |
| 27 | Tyr | Arg | Arg |
| 36 | **Ala** | Gly | Gly |
| 106 | Ser | Ser | **Pro** |
| 111 | Leu | Leu | **Gln** |

**`Thr36` in the white paper is wrong.** Position 36 is Ala in AncCR and Gly in
both derived receptors; Thr never occurs there. It is also 3.7 Å from
dexamethasone — pocket-lining, not a remote residue. Ortlund 2007 (our local
copy) names the permissive set as **Y27R, L29M, F98I, Q105L, N26T, S212D** with
**S106P/L111Q** as the switch; `Thr36` appears nowhere in that literature. The
intended second permissive residue is almost certainly **Leu29**.

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
| Rhodopsin Glu113 ↔ cytoplasmic face, 22.1 Å | 29.6 Å to Arg135; 19.7 Å from retinal | not reproduced |

## Case-by-case sufficiency

**Case 1, hemoglobin — structure correct, site set in the paper is wrong.**
2HHB is clean (4 chains, 1–141/1–146, no gaps, 4 haems). But Storz 2009, our
local copy, defines the deer mouse system as **α sites 50, 57, 60, 64, 71** and
**β sites 62, 72, 128, 135** — a 2⁵ and 2⁴ design, not the "2⁴ = 16 and
2⁸ = 256" the white paper claims. More seriously, the mechanism is
contradicted: α57 and α64 sit in the **haem pocket** (6.5–7.5 Å from haem),
15.3 Å and 18.9 Å from the nearest β1 residue. They cannot "disrupt" or
"tighten" the α1β1 packing face. The Table 1 pair α114↔β116 *is* a real α1β1
contact (2.7 Å) but neither residue is in the Storz site set.

**Case 2, myoglobin — correct scaffold, no ancestral coordinates exist.**
1MBO is *Physeter catodon* (the target taxon), 1.6 Å, 153 residues, no gaps,
haem + bound O₂, His93–Fe 2.1 Å. Mirceta 2013 crystallised no ancestors, so
1MBO can only serve as a geometric scaffold for charge mapping.

**Case 3, steroid receptor — now complete.** Added `2Q1V` and `2Q3Y` to finish
the Ortlund 2007 deposition set. Note `2Q1V` is titled "in complex with
cortisol" but its deposited ligand is `PDN` = prednisone (C21H26O5), *not*
cortisol (C21H30O5); no ancestral entry here binds cortisol itself. Modern
human GR entries (`4P6X` with `HCY`, `1M2Z`) were evaluated and **removed** —
the goal is to classify epistatic pairs on the resurrected backgrounds, and
modern structures use NR3C1 numbering that does not match the AncCR series.

**Case 4, rhodopsin — correct and sufficient as a scaffold.** 1U19 carries all
named residues in standard bovine numbering: Asp83, Glu113, Gly121, Glu122,
DRY motif Glu134-Arg135-Tyr136, Ala292, with retinal bound. Ala292 and Gly121
are 3.3 Å and 3.7 Å from retinal; they are 10.3 Å apart, so the
permissive/adaptive pair is a genuine second-shell interaction. The local
Dungan & Chang 2017 supplement assays sites 83, 292 and 299, consistent with
the paper's D83N/A292S.

## Recommended coordinate set

| Case | Use |
| --- | --- |
| 1 | `2HHB` (only; re-derive the site set from Storz 2009) |
| 2 | `1MBO` |
| 3 | `2Q1H`, `2Q1V`, `2Q3Y` (AncCR) → `3RY9` (AncGR1) → `3GN8` (AncGR2) |
| 4 | `1U19` |


---

# Contact-mapping fitness of 2HHB / 1MBO / 1U19 (2026-09-17)

Validated against a narrower goal than the §4 accession check: given an
arbitrary epistatic site pair emitted by the model in a §3 case, can this
structure assign it to *direct contact* or *allosteric network* reliably?
Reproduce with `epistaeon/scripts/numbering.py`, `mapfitness.py`, `bands.py`.
Case 3 now uses the Thornton ancestral structures, so `4P6X`, `1M2Z` and
`4LS6` have been removed from the repository.

## Passes

**Numbering is exact and unambiguous.** Aligned to UniProt canonical sequences,
zero mismatches in all four chains:

| Chain | UniProt | Offset | Coverage |
| --- | --- | --- | --- |
| 2HHB A (α) | P69905 | UniProt = PDB **+1** | 141/141 modelled |
| 2HHB B (β) | P68871 | UniProt = PDB **+1** | 146/146 modelled |
| 1MBO A | P02185 | UniProt = PDB **+1** | 153/153 modelled |
| 1U19 A | P02699 | UniProt = PDB **+0** | 348/348 modelled |

The +1 offsets are initiator-Met cleavage in the mature protein. They must be
applied explicitly; an unconverted alignment index silently shifts every
globin and myoglobin contact call by one residue.

**No missing coordinates.** Every SEQRES position is modelled in all three
structures — no unmodelled loops, so no detected pair is unmappable. Zero
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

## Blockers

**1. The two-band scheme does not cover most pairs.** §4 defines direct as
≤5 Å and allosteric as 10–25 Å. Measured over all pairs with |i−j| ≥ 5:

| Structure | direct ≤5 Å | *undefined* 5–10 Å | allosteric 10–25 Å | *undefined* >25 Å | total unclassified |
| --- | --- | --- | --- | --- | --- |
| 2HHB α | 2.3% | 15.6% | 68.8% | 13.2% | **28.9%** |
| 2HHB β | 2.4% | 15.3% | 65.7% | 16.6% | **31.9%** |
| 1MBO | 2.2% | 14.7% | 65.7% | 17.4% | **32.1%** |
| 1U19 | 1.3% | 7.4% | 44.2% | 47.1% | **54.5%** |

For rhodopsin the majority of pairs fall outside both bands. The scheme needs
a third intermediate class and an upper bound set by each protein's actual
diameter (38–73 Å here), not a fixed 25 Å.

**2. Hemoglobin α–β pairs cannot be classified at all without an interface
convention.** Of 20,586 α–β pairs, 92 are contacts at α1β1 and 92 at α1β2 —
but **zero are contacts at both**. The class of every inter-chain pair is
decided entirely by which interface you measure on. This must be fixed
explicitly (and reported) before any Case 1 inter-chain pair is scored.

**3. Cofactor-mediated coupling is invisible to a residue–residue criterion.**
Pairs where both residues line the cofactor (≤5 Å) yet sit >10 Å apart, and so
would be mis-called allosteric:

| Structure | cofactor | lining residues | mis-called pairs | share of lining pairs |
| --- | --- | --- | --- | --- |
| 2HHB α | haem | 23 | 103 | 46% |
| 1MBO | haem | 22 | 100 | 49% |
| 1U19 | retinal | 25 | 108 | 41% |

This is the most consequential gap: the phenotypes in Cases 1, 2 and 4 (P50,
λmax) are *all* cofactor-mediated, so the interactions the model most needs to
detect are exactly the ones this criterion misassigns. Contact definitions
should treat haem and retinal as bridging nodes.

**4. Cross-species transfer is unvalidated.** 1MBO is *Physeter catodon*, the
Case 2 target taxon. But 2HHB is human while Case 1 targets *Peromyscus* and
*Merganetta*, and 1U19 is bovine while Case 4 targets cetaceans and bats. The
numbering above is verified against the *human* and *bovine* references only;
an alignment from each target ortholog to the PDB chain, with explicit indel
handling, is still required.

## Power requirement for §5

At the measured background contact rates, "Odds Ratio > 4.0, Fisher's exact
p < 10⁻⁵" requires a minimum number of detected pairs:

| Structure | background ≤5 Å rate | OR=4 implies | pairs needed |
| --- | --- | --- | --- |
| 2HHB α | 2.30% | 8.6% contacts | ≥ 181 |
| 2HHB β | 2.41% | 9.0% contacts | ≥ 173 |
| 1MBO | 2.22% | 8.3% contacts | ≥ 187 |
| 1U19 | 1.27% | 4.9% contacts | ≥ 317 |

A detector emitting fewer pairs than this cannot meet the §5 standard on that
case regardless of accuracy.
