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
where the human reference has a gap first: HBB 6, RHO 8, MB 32, HBA1 32,
HBA2 62, NR3C1 88, NR3C2 110.

**Then convert to structure numbering** (`../data/numbering_offsets.json` is
the machine-readable source):

| Structure | Protein | Offset |
| --- | --- | --- |
| 2HHB chains A, C | alpha-globin | UniProt = PDB **+1** |
| 2HHB chains B, D | beta-globin | UniProt = PDB **+1** |
| 1MBO | myoglobin | UniProt = PDB **+1** |
| 1U19 | rhodopsin | **0** |
| 7PRX | glucocorticoid receptor | **0** |
| ancestral receptors | AncCR / AncGR1 / AncGR2 | local LBD numbering; **+531** to reach human GR |

Special cases:

- **NR3C1 alignment is the GR-gamma isoform.** An extra Arg at 452 means every
  human-reference position after 451 is one higher than UniProt / 7PRX.
  Subtract 1 before applying the 7PRX offset.
- **3GN8** numbering runs one lower than AncCR from about position 212, so the
  +531 rule does not hold there.
- **1U19** carries an `ACE` cap at author position 0; skip it when iterating.
- **1MBO and 1U19 are not human**, but neither has an indel relative to the
  human protein: human residue *n* is 1U19 residue *n* and 1MBO residue *n* − 1.

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
cannot reach significance regardless of accuracy.

## Species coverage

The alignments are TOGA2, human-referenced, and therefore **mammals only**: all
3,847 sequences and all 882 species-tree tips. Birds, fish and reptiles are
absent by construction.

| Needed species | Status |
| --- | --- |
| *Peromyscus maniculatus* (deer mouse) | In TOGA2 (3 assemblies) but **absent from HBA1, HBA2 and HBB** |
| *Physeter macrocephalus* (sperm whale) | Present in MB and RHO |
| *Mirounga leonina*, *Leptonychotes weddellii* | Present in MB |
| Naked mole-rat, star-nosed mole, *Myotis nattereri* | Present in RHO |
| *Rhinolophus unihastatus* | No TOGA2 assembly |

Data quality: share of codons that are gaps or masked — HBA1 34.5%, HBA2 44.0%,
MB 19.5%, NR3C2 13.9%, NR3C1 13.5%, HBB 7.8%, RHO 5.0%. Grey seal NR3C1 has two
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

## Model constraint

The published HyphAeon checkpoint scores each alignment column independently:
`window_size: 1`, no cross-site layers, and the batch dimension is the site.
Its only output head is an ordinal regressor predicting HyPhy-style statistics
(LRT, −log10 p, omega). There is no likelihood head and no energy function in
the released weights.

This matters when judging an ordering or contingency result: if such a result
appears, establish where the background dependence came from, because the
site-level scores alone cannot supply it.
