# Data facts for validation agents

Operational facts about the inputs in this repository, verified by the scripts
in `epistaeon/scripts/`. These are properties of the data, not claims about
biology. **The original published studies remain the only source of truth for
what is true about proteins and evolution.** Nothing here should be treated as
a finding to replicate.

## Numbering — the main silent-failure risk

Alignment positions, PDB positions and UniProt positions are three different
coordinate systems. Convert explicitly; a mistake produces a plausible-looking
wrong answer with no error.

**Model output is alignment columns, not residue numbers.** Drop the columns
where the human reference has a gap first: HBB 6, MB 32, HBA1 32,
HBA2 62, NR3C1 88, NR3C2 110. These counts apply only to the TOGA2 inputs.
The rebuilt SWS1 alignment has its own explicit source-site/reference-column
crosswalk; do not apply TOGA2 gap counts to it.

**Then convert to structure numbering** (`../data/numbering_offsets.json` is
the machine-readable source):

| Structure | Protein | Offset |
| --- | --- | --- |
| 2HHB chains A, C | alpha-globin | UniProt = PDB **+1** |
| 2HHB chains B, D | beta-globin | UniProt = PDB **+1** |
| 1MBO | myoglobin | UniProt = PDB **+1** |
| 1U19 | bovine rod-opsin homology template for SWS1 | **0** to its own bovine UniProt sequence only; use the SWS1 crosswalk for source sites |
| 7PRX | glucocorticoid receptor | **0** |
| ancestral receptors | source AncCR numbering | **+531** through210, **+530** from213; source211 ambiguous, deleted212 has no modern homolog |

Special cases:

- **NR3C1 alignment is the GR-gamma isoform.** An extra Arg at 452 means every
  human-reference position after 451 is one higher than UniProt / 7PRX.
  Subtract 1 before applying the 7PRX offset.
- **3GN8** numbering runs one lower than AncCR from about position 212, so the
  +531 rule does not hold there.
- **1U19** carries an `ACE` cap at author position 0; skip it when iterating.
- **1MBO is not human**, but has no indel relative to the human myoglobin
  reference: human myoglobin residue *n* is 1MBO residue *n* − 1.
- **1U19 is not SWS1.** It is the bovine rod-opsin template used for homology
  modeling in Yokoyama 2014. Its complete coordinates and own-sequence offset
  do not validate a direct SWS1 mapping. Historical human accession M13295 and
  current UniProt P03999 differ; use the source-anchored crosswalk.

## Contact classification

A pair is in direct contact if any heavy atom of one residue lies within 5.0 A
of any heavy atom of the other. Exclude pairs with |i − j| < 5, since chain
neighbours are always in contact.

Baseline contact rates, needed to judge whether any enrichment is meaningful:

| Structure | pairs in direct contact |
| --- | --- |
| 2HHB alpha | 2.3% |
| 2HHB beta | 2.4% |
| 1MBO | 2.2% |
| 1U19 | 1.3% |
| 7PRX | 1.7% |

Minimum detected pairs for an odds ratio of 4.0 at Fisher p < 1e-5: 2HHB alpha
181, 2HHB beta 173, 1MBO 187, 1U19 317, 7PRX 247. Below these counts the test
cannot reach significance regardless of accuracy. The 1U19 values describe
all template residues; they are not a calibrated SWS1 benchmark after homolog
mapping or restriction to seven focal sites.

## Species coverage

The TOGA2 alignments are human-referenced and **mammals only**; the shared
species tree has 882 tips. Birds, fish and reptiles are absent by construction.
The separate `data/alignments/sws1/` input is rebuilt from the paper-listed
extant boreoeutherian accessions, including additional primates. It contains
no outgroup opsins, non-boreoeutherian taxa or ancestral rows, and does not
claim to recover the authors' exact original alignment.

| Needed species | Status |
| --- | --- |
| *Peromyscus maniculatus* (deer mouse) | High-type representative now merged into HBA1, HBA2 and HBB as `HLperManSon3`; not the complete experimental allele panel |
| *Physeter macrocephalus* (sperm whale) | Present in MB |
| *Mirounga leonina*, *Leptonychotes weddellii* | Present in MB |
| Human, macaque, squirrel monkey, marmoset, cow, mouse, rat, gray squirrel, guinea pig | SWS1 accessions from the paper's original 33-sequence set |
| Chimpanzee, gorilla, orangutan, gibbon, tarsier, aye-aye | Additional SWS1 primate accessions from the paper's Methods |

Data quality: share of codons that are gaps or masked — HBA1 34.5%, HBA2 44.0%,
MB 19.5%, NR3C2 13.9%, NR3C1 13.5%, HBB 7.8%. Grey seal NR3C1 has two
internal stop codons. The loader collapses identical sequences before inference.

## Structures

All complete with no missing side chains unless noted. 2HHB, 1MBO and 1U19 have
every SEQRES position modelled; 7PRX lacks only its terminal residues 528 and
777. Contact calls are stable across crystallographic copies (median difference
0.10–0.20 A, essentially no class flips).

Two structural constraints that affect interpretation:

- **2HHB is the deoxy (T) state.** Roughly half the alpha1beta2 contacts differ
  in the oxy (R) state, so some alpha-beta pairs are state-dependent.
- **Haem and retinal are large.** Residues coupled through a cofactor without
  touching each other count as "not in direct contact" under the 5 A rule. That
  is correct under the rule, but worth noting when interpreting a negative.
- **SWS1 has a homology template, not a deposited focal structure.** Map each
  source site explicitly to 1U19 before reporting its location. Template
  residue identity need not match any SWS1 experimental background.

## SWS1 ordering evidence

Yokoyama et al. (2014), doi:10.1371/journal.pgen.1004884, assayed all 127 mutant
combinations of F46T, F49L, T52F, F86L, T93P, A114G and S118T. Sixty mutants
did not form functional pigments; the paper reports 1,032 accessible orders out
of 5,040. T52F alone fails on AncBoreotheria S1, so it cannot initiate a path
under that assay criterion. This does not measure organismal fitness.

The narrower historical inference is
`{T93P, A114G} -> F86L -> {F49L, S118T} -> {F46T, T52F}`: eight paths because
each of three pairs is internally unordered. This order comes from
phylogeny/parsimony, not a demonstration that all other functional paths are
impossible. No dark-noise or thermal-decay prerequisite is part of this case.
The local S1-table transcription has genotype/formatting anomalies and its
enumeration differs from the reported path count. Exact accessibility-set
scoring must wait for reconciliation; retain published and enumerated counts
separately.

## Model constraint

The published HyphAeon checkpoint scores each alignment column independently:
`window_size: 1`, no cross-site layers, and the batch dimension is the site.
Its only output head is an ordinal regressor predicting HyPhy-style statistics
(LRT, −log10 p, omega). There is no likelihood head and no energy function in
the released weights.

This matters when judging an ordering or contingency result: if such a result
appears, establish where the background dependence came from, because the
site-level scores alone cannot supply it.

## SWS1 input checks

The rebuilt input has 15 distinct nucleotide rows and 350 codon columns;
all seven focal sites are present and unambiguous in every row. Human and
chimp proteins are identical but their nucleotide sequences differ, so this
panel does not lose a taxon to nucleotide duplicate pruning. Marmoset has a
short missing C terminus and rat has one retained ambiguous codon near the
N terminus; neither affects the focal sites. The human-reference site list is
`data/alignments/sws1/sites_human.txt`; explicit codon columns are in
`sites_columns.txt`. Consult the dataset README before running without an
ancestral row: the model's plurality states are not the published ancestor.
