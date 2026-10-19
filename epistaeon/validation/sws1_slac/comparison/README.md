# SWS1 SLAC ancestral-state comparison

The frozen `Boreoeutheria` reconstruction matches **7/7 focal ancestral residues** and **278/282 available paper-core residues**.

Among 282 resolved core comparisons, there are 4 substitutions and 0 unresolved or unmapped positions. Resolved amino-acid identity is 98.58%.

| Substitution | Paper ancestor | SLAC ancestor | Reconstructed codon | Alignment column | Result |
|---|---|---|---|---:|---|
| F46T | F | F | TTT | 45 | match |
| F49L | F | F | TTT | 48 | match |
| T52F | T | T | ACA | 51 | match |
| F86L | F | F | TTC | 85 | match |
| T93P | T | T | ACT | 92 | match |
| A114G | A | A | GCC | 113 | match |
| S118T | S | S | TCT | 117 | match |

Coordinates were bridged by the unique identical S4-human core, independently of the reconstructed ancestor. Actual input gaps are retained. The figure's red terminal blocks are excluded; its sequential protein indices differ from the Methods prose bounds. [S4 Fig](https://ndownloader.figshare.com/files/1847637), [Yokoyama et al. 2014](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004884).

The paper provides a chimeric assay protein, rather than a complete natural ancestor or ancestral DNA. Consequently this comparison covers the 282-residue non-engineered core and the seven focal amino acids; it cannot assess true ancestral codons or natural termini. This is a separate SLAC reconstruction using rebuilt inputs, not an exact rerun of the original reconstruction. Agreement does not validate substitution order or epistasis.

Inference artifact checksums were verified before reading empirical truth and again before writing these outputs. No inference settings or predictions were modified. Full core discrepancies, unavailable states and coordinate mappings are in `core_residues.tsv`; hashes and detailed metrics are in `comparison.json`.
