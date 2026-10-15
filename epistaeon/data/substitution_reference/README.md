# Source-verified substitution reference

Reviewed 2026-10-15 against the original studies. This reference extends and
corrects the older lists in `validation/studies/`. It is ground truth to compare
against model results, not a record of a model run.

| Case | Reference | Coverage and ordering evidence |
| --- | --- | --- |
| 1, hemoglobin | [hemoglobin.json](hemoglobin.json) | Storz 2009: five alpha diagnostic sites and four polarized beta substitutions. No measured order. Related Natarajan 2013 evidence is kept separately, with its remaining source-data gaps. |
| 2, myoglobin | [myoglobin.json](myoglobin.json) | Mirceta 2013 Fig. S2: all 199 labeled charge-changing substitutions in 85 branch groups. Published branch ancestry gives a partial order; changes on one branch remain unordered. Every site uses modern **1MBO A**. |
| 3, steroid receptor | [steroid_receptor.json](steroid_receptor.json) | Ortlund 2007 X/Y/Z and earlier Y27R; Bridgham 2009 W; Carroll 2011 sensitivity/buffering substitutions and controls. Seven qualified constraints distinguish assayed group paths, preferred paths, branch precedence and conjecture. |
| 4, SWS1 cone opsin | [sws1.json](sws1.json) | Yokoyama 2014: seven substitutions, a complete 127-mutant pigment-formation panel, a reported 1,032 accessible paths and eight phylogenetically favored paths. Historical order and experimental accessibility are distinct targets; local S1-table/path-count anomalies remain explicit. |

[sites.tsv](sites.tsv) is the flat, filterable residue-location list. The JSON
files hold the biological scope, sources, group memberships, constraints and
missing validation inputs. Source hashes identify the exact local files reviewed.
Mapped events record the structural chain, author position, observed template
residue and C-alpha coordinates, with their sequence/alignment crosswalks.
The original three cases use TOGA2 `hg38` codon columns. SWS1 uses the cited
human sequence M13295 and the rebuilt SWS1 alignment; bovine 1U19 is its
homology template, not an experimentally determined SWS1 structure. Source
site labels, sequence residues, alignment columns and template coordinates are
different coordinate systems.
S212del and the ambiguous Q211E boundary retain null coordinates with explicit
reasons. `order_scoring` marks eligibility only in the stated evidence context;
it never supplies a total order, and `order_evidence` preserves its qualification.

## How to read an ordered set

An array of events is **not** a total order. A group's members are unordered
unless a separately cited constraint says otherwise. Constraints apply only to
their stated background and resolution:

- `phylogenetic_branch_precedence`: an ancestral branch precedes its descendant;
  this does not demonstrate functional necessity or order within a branch.
- SWS1's four historical substitution sets summarize phylogenetic/parsimony
  inference. They do not claim that all other functional paths are impossible.
  Pigment-formation measurements define a separate accessibility test.
- `prerequisite_groups_before_group_completion`: tested prerequisites precede
  completion of a dependent set. This does not require every prerequisite member
  to precede every dependent member.
- `preferred_pair_order_from_functional_intermediates`: a better-supported
  functional path, not an absolute prohibition or directly observed history.
- `reverse_group_prerequisite`: applies to the assayed reverse path; do not
  silently convert it into forward chronology.
- `author_conjecture_supported_by_buffering_assay`: retain for interpretation;
  exclude from a score of experimentally established order.

`completion(Y)` refers to the whole Y set. `reverse(W)` refers to restoring the
ancestral W states, not introducing its forward substitutions. Do not compare
events on sister branches as if they belonged to one trajectory. Myoglobin
relations link the nearest annotated ancestor; transitive precedence follows
along the same lineage.

## Case-specific results

**Hemoglobin.** Storz Fig. 3 defines alpha positions 50/57/60/64/71 as
`PGAGS`, `HGAGS`, `HAGDG`, and beta positions 62/72/128/135 as `GGAA`, `ASSS`.
The beta outgroup polarization is **A62G, G72S, S128A, A135S**; its per-site
ancestral string is `AGSA`, neither sampled beta endpoint. Alpha direction is
not established here. The measured hemolysates contain mixed alpha isoforms.
Diagnostic variation is not evidence that every individual site is epistatic.

**Myoglobin.** The modern scaffold is sufficient for locating transitions, as
requested. For example, the cetacean-stem set
`{G1V, D53A, Q116H, N140K}` precedes the Physeteroidea `N12H` event on the
sperm-whale lineage. Those four stem events are internally unordered. Beaked
whales' separate `{N12H, N66H}` set and dolphin-lineage changes are on other
branches and must not be added to the sperm-whale trajectory. Net surface charge
uses the paper's site-specific pKa treatment at pH 6.5 (Table S2); it is not a
simple count of charged side chains. No mutant panel establishes site-pair
epistasis or obligatory within-branch orders in this study.

**Steroid receptor.** All positions below use original AncCR numbering:

| Set | Events | Supported constraint |
| --- | --- | --- |
| Earlier permissive | Y27R | Earlier reconstructed interval than X/Y/Z/W. |
| X | S106P, L111Q | Fig. 2A favors L111Q before S106P; low sensitivity on the opposite path does not establish an absolute impossibility. |
| Z | N26T, Q105L | Both participate in rescue of the tested XY background. |
| Y | L29M, F98I, S212del | X and Z precede full Y on the assayed functional group paths. No internal Y or X-versus-Z total order. |
| W | H84Q, Y91C, A107Y, G114Q, L197M | Full W is unfavorable before X remodeling. On the assayed whole-block reverse paths, restore Y before W, and W before X; partial W reversions can also partially rescue function. |
| 2011 earlier interval | V43A, R116H, C71S | C71S buffers the desensitizing pair. The authors **conjecture** C71S-first; they did not measure historical chronology. |

The 2011 updated AncGR1.1 background differs from 2007/2009 AncGR1.0, including
site 107. Do not concatenate their event labels into one endpoint trajectory.
The reference explicitly flags weaker single-mutant effects and negative or
excluded controls so they are not scored as demonstrated epistatic partners.

S212 is a **deletion**, with no homologous modern residue coordinate. A generic
sequence alignment gaps Q211 instead and pairs S212 with Asp742; the original
paper explicitly identifies the deleted serine. The reference records null
coordinates for S212 and flags the 211/212 boundary for source-anchored
curation. Candidate flanks 742/743 are context, not a confirmed deleted-residue
mapping. Other modern residues can also differ from the experimental endpoints,
for example N26T maps geometrically to modern Trp557.

**SWS1 cone opsin.** Yokoyama et al. (2014),
[doi:10.1371/journal.pgen.1004884](https://doi.org/10.1371/journal.pgen.1004884),
tested **F46T, F49L, T52F, F86L, T93P, A114G and S118T** on the reconstructed
AncBoreotheria S1 background. The seven-change mutant shifts the absorption
maximum from 357 to 411 nm. Sixty of 127 mutants failed to form functional
pigments; the paper reports 4,008 of 5,040 orders excluded and 1,032 accessible.
T52F alone is nonfunctional on that ancestral background, so it cannot be the
first step under the paper's pigment-formation criterion. This is not a
measurement of organismal fitness or dark-noise kinetics.

The phylogenetic/parsimony reconstruction narrows the historical candidates to
**`{T93P, A114G} -> F86L -> {F49L, S118T} -> {F46T, T52F}`**. Each brace is
internally unordered, yielding eight favored paths. The experiment does not
establish this partial order as an obligatory biochemical sequence. Keep the
paper-reported 1,032 accessible paths separate from these eight historically
favored paths. Interaction estimates for the nonfunctional mutants use
inferred spectral values; do not treat those values as measured absorption
peaks.

The local S1-table transcription exposes a duplicated/missing genotype label
and a nonfunctionality-formatting anomaly. Its normalized enumeration does
not yet reproduce the reported accessible-path count. Published counts and
locally enumerated counts remain separate in the source artifacts; resolve
that discrepancy before scoring exact agreement with the full accessibility
set. The historical four-set order is a separate source claim.

The input is an independently rebuilt codon-aware alignment of extant
boreoeutherian accessions listed in the paper, including its additional
primate sequences: [SWS1.codonified.fa](../alignments/sws1/SWS1.codonified.fa).
It contains no ancestral rows and is not the authors' recovered input
alignment. The experimental ancestral protein had anole-derived termini;
those segments must not be substituted into extant alignment rows. Source
site numbering must use the explicit SWS1 crosswalk rather than current human
UniProt numbering or a constant offset to 1U19.

## What is still needed to validate model output

1. **Matching model inputs.** Select the exact original branch, reconstruction
   version, taxa and experimental sequence backgrounds. Local receptor AA FASTAs
   and posterior probabilities already exist, but are not source-version-checked
   codon-aligned endpoint rows. The other original reconstructed branches also
   need explicit endpoints if an endpoint-based trajectory run is intended.
   The SWS1 extant alignment is sufficient as the requested comparative input;
   it does not supply the experimental ancestral background to the current
   endpoint-based `order` CLI. An extant-versus-extant contrast or plurality
   consensus is not the original ancestor. Historical codons and synonymous
   order are generally not supplied by these studies.
2. **Run-specific crosswalks.** Match each output column/unit to the source event,
   target ortholog, endpoint residue and active PDB/template. The included
   modern-reference columns are not CLI `--sites` positions on an ancestral
   row. Resolve the receptor deletion boundary and globin paralogs; keep the
   GR-gamma insertion, cleaved initiator methionines and SWS1 source/template
   crosswalk explicit.
3. **Experimental measurements and provenance.** Archive missing 2007/2009
   receptor supplements and the Storz supplement. Deer-mouse representatives
   are now in all three globin
   alignments, but they are not the whole segregating, phased experimental panel.
   Natarajan 2013 needs its verified supplement before its extra blocks and
   accessible paths can be scored. For SWS1, distinguish experimentally observed
   pigment formation from inferred spectral values for nonfunctional mutants;
   reconcile S1-table genotype/formatting anomalies and path counts before
   claiming exact recovery of the experimentally accessible set.
4. **A comparator with declared rules.** Join `order.json` unit members to these
   event IDs; evaluate supported pair/set-completion probabilities, preserve
   unknown orders, separate forward/reverse contexts and conjectures, handle
   deletion units and distinguish missing coverage from wrong predictions.
   Predeclare thresholds, uncertainty treatment, phenotype/assay condition and
   null baseline. There is no reference-aware automated comparator yet.
5. **Evidence the studies do not provide.** Storz has no mutation-order panel;
   Mirceta has no measured within-branch order or individual epistasis panel.
   SWS1's eight favored paths are historical inference, not the full measured
   accessibility set or a unique observed history. No additional scaffold or
   model result can supply missing ground-truth labels.

Use modern structures to check location and a separately declared structural
enrichment question. A 5-Angstrom contact is not itself evidence of epistasis or
temporal precedence. Myoglobin surface classification needs a specified exposure
rule. SWS1 contact claims based on 1U19 concern a mapped homolog template;
they are not direct measurements of the SWS1 mutant structures.

## Integrity check

From the repository root:

```bash
python3 epistaeon/scripts/validate_substitution_reference.py
```

This checks JSON links, source hashes, event/group references, noncyclic
myoglobin branch relations, residue identities and coordinates against the
active CIFs, numbering offsets, and the included reference-alignment columns.
It also regenerates the flat list with `--write-tsv`. It does not prove the
literature interpretation or test the model.

SWS1 machine-readable assay records are
[`sws1_2014_genotypes.json`](../experimental/sws1_2014_genotypes.json) and
[`sws1_2014_paths.json`](../experimental/sws1_2014_paths.json). The latter lists
all eight historical candidates and 1,076 paths from the explicitly normalized
S1 table; it does not claim to reproduce the paper's 1,032. A table-derived
rule suggesting **F46T or T93P before T52F** is provisional because it depends
on the unconfirmed genotype-label repair. The directly reported result that
**T52F cannot be first** does not depend on that repair.
