# Paper SWS1 protein reference for independent ASR evaluation

This directory contains **evaluation-only source truth** transcribed from
[Yokoyama et al. (2014), S4 Fig](https://doi.org/10.1371/journal.pgen.1004884.s004).
It is separate from extant alignment inputs and inferred ASR outputs. Do not
use its ancestral states to select inference parameters, models, trees or
root positions.

## What the original source provides

S4 displays three complete, aligned protein rows: **AncBoreotheria S1**,
**pigment g** and **human S1**. Each displayed row contains **347 residues**
without an alignment gap. AncBoreotheria is the reconstructed Boreoeutherian
SWS1 node, ancestral to Euarchontoglires and Laurasiatheria. It is not the
primate ancestor or the root of all mammals. Pigment g is the earlier ancestral
mammalian construct discussed in the paper, not the same inferred node or
sequence as AncBoreotheria.

The displayed constructs have **anole-derived terminal segments**. These are
red in S4; dots in the two lower rows mean identity with the ancestral row,
including the red standardized termini. Therefore the displayed full protein
must be described as an assayed/terminal-standardized construct, rather than a
fully natural ancestral sequence. The source does not provide ancestral codons.

The [main paper's Methods](https://doi.org/10.1371/journal.pgen.1004884) describe
the N and C termini as “sites 1–30” and “sites 313–348.” These source labels
are not sequential indices in the displayed protein. The actual S4 red-color
mask is displayed-protein positions **1–25 and 308–347**. The available natural
ancestral protein comparison is consequently restricted to the **282-residue
non-red core at positions 26–307**. The literal Methods bounds and this
distinction are preserved in `source_reference.json`.

## Transcription and independent checks

The ancestral row and every non-dot replacement in the other two rows were
transcribed directly from the original EPS rendered at 450 dpi. OCR was a
secondary visual aid. The saved rows were not derived from an ASR output or
from a human sequence with seven substitutions reversed.

Checks agree with the original paper:

- **25 differences** between AncBoreotheria and the S4-human row.
- **9 differences** between AncBoreotheria and pigment g.
- The entire **282-residue S4-human core** matches the historical
  accession-derived human protein exactly and uniquely at human CDS positions
  **28–309**.
- This bridge maps the core to the existing 350-codon benchmark alignment at
  codon columns **30–311**. No engineered terminus is used for matching.
- All seven ancestral focal states agree with the reviewed numbering
  crosswalk. All seven columns have unambiguous complete codons in all 15
  extant benchmark rows.

| Paper label | S4 protein index | Historical human CDS index | Benchmark codon column | Ancestor | Human |
| --- | ---: | ---: | ---: | --- | --- |
| F46T | 41 | 43 | 45 | F | T |
| F49L | 44 | 46 | 48 | F | L |
| T52F | 47 | 49 | 51 | T | F |
| F86L | 81 | 83 | 85 | F | L |
| T93P | 88 | 90 | 92 | T | P |
| A114G | 109 | 111 | 113 | A | G |
| S118T | 113 | 115 | 117 | S | T |

## Files

- `source_reference.json`: full transcriptions, color mask, original source
  hashes, node identity, coordinate bridge, focal states and caveats.
- `paper_S4.protein.fa`: three expanded 347-residue source protein rows.
- `paper_S4.dot_alignment.txt`: original block boundaries and dot notation.
- `paper_S4.core_positions.tsv`: source/human/alignment coordinates for every
  comparable residue, with all three protein states.
- `source_manifest.json`: SHA256 hashes for these evaluation reference files
  and their original sources.
- `verify_source_reference.py`: deterministic source-transcription checks;
  it performs no ancestral inference.

The original EPS and main HTML remain in `data/experimental/`. A fully natural
ancestral terminal sequence, the ancestral DNA, original author alignment and
ancestral posterior probabilities are unavailable in the reviewed source
package. Do not create supposed ancestral nucleotides by backtranslation.

## Rooting context

The benchmark's 15 taxa contain 14 Euarchontoglires and only one
Laurasiatherian, *Bos taurus*. An unrooted tree fitted under a reversible model
does not establish the biological Boreoeutherian root by itself.

The same paper explicitly lists three other mammalian SWS1 accessions suitable
for a **separate ASR-only outgroup augmentation**: elephant *Loxodonta africana*
**AY686753**, wallaby *Macropus eugenii* **AY286017**, and dunnart *Sminthopsis
crassicaudata* **AY442173**. Marsupials can provide an external placental root;
elephant provides a placental lineage outside Boreoeutheria. The target node
is the MRCA of the original boreoeutherian taxa, rather than the root of the
entire augmented tree. Preserve the original 15-taxon benchmark alignment.
