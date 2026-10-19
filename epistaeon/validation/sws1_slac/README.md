# SWS1 ancestral reconstruction with HyPhy SLAC

The reconstructed **Boreoeutheria** node matches the paper at **all seven focal
ancestral residues**. Over the available non-engineered protein core it matches
**278/282 residues (98.58%)**. Four resolved amino-acid differences remain:

| Historical human CDS position | Original alignment codon column | Paper ancestor | SLAC ancestor |
| ---: | ---: | --- | --- |
| 54 | 56 | I | M |
| 76 | 78 | I | V |
| 94 | 96 | S | A |
| 295 | 297 | S | A |

All four differences equal the human state. These coordinates were derived
through the identical human core, rather than a constant numbering offset.
The focal ancestral pattern is **F, F, T, F, T, A, S** at paper positions
**46, 49, 52, 86, 93, 114, 118**, respectively. Consequently all five ancestral
states missed by the extant plurality are recovered in this test.

## Experiment

The original `data/alignments/sws1/SWS1.codonified.fa` has 15 extant taxa and
350 codon columns. All 15 rows and every original column were retained exactly.
Three SWS1 outgroups listed in the original paper were added only in this
separate experiment:

- Elephant, *Loxodonta africana*: AY686753.1.
- Wallaby, historical *Macropus eugenii*, current *Notamacropus eugenii*: AY286017.1.
- Dunnart, *Sminthopsis crassicaudata*: AY442173.1.

MAFFT **7.526** added translated outgroup sequences to the existing protein
alignment with `--add --keeplength --mapout --globalpair --thread 1`. Its residue
map was used to backtranslate the additions. One elephant N-terminal glycine
insertion was excluded to preserve the original alignment columns; its original
codon and position are recorded in `inputs/preparation.json`. Both marsupial
coding sequences were retained in full after terminal-stop removal.

The fixed species-scaffold tree groups primates and rodents as
Euarchontoglires, with cow as the sampled Laurasiatherian. Deep rodent
relationships were left as a polytomy. Branch lengths were fitted by HyPhy;
no topology search or truth-dependent tree selection was performed. Elephant
and the two marsupials locate the named **internal** Boreoeutheria node. The
degree-three computational root is on the marsupial side and is not the target
ancestor. Both input and SLAC output trees were checked to ensure that the
target subtends exactly the original 15 taxa.

HyPhy **2.5.96(MP)** ran SLAC **2.00**, using the Universal genetic code,
one CPU, all branches, `--samples 0`, `--pvalue 0.1`, and
`--kill-zero-lengths No`. SLAC fitted GTR and MG94xREV and exported joint
maximum-likelihood ancestral codon states. This is a point reconstruction;
per-residue posterior probabilities were not computed. The p-value threshold
controls SLAC's selection tests, not ancestral-state acceptance.

The fit completed successfully in about one second. Its raw JSON is unchanged.
An initial FASTA export encountered the retained ambiguous rat codon `GWC`:
SLAC represents its compatible amino-acid set as the single array cell `DV`.
Only export code was corrected, representing that unresolved cell as `X` and
preserving its raw states in `slac.json` and `sequence_export.json`. No model
refit was performed. The selected ancestor is fully resolved at all 350 codons.

## Comparison and limits

[Yokoyama et al. (2014), S4 Figure](https://doi.org/10.1371/journal.pgen.1004884.s004)
provides full protein rows for AncBoreotheria, pigment g and the human assay
construct. The displayed proteins have engineered anole-derived termini.
The comparison uses the **282-residue non-red core**, S4 sequential protein
positions 26–307, mapped through the exact human sequence into our codon
alignment. The original Methods terminal bounds and figure-color mask are
preserved separately in the [source reference](../sws1_slac_reference/README.md).

The paper supplies no ancestral nucleotide sequence or full natural ancestral
termini, so neither can be checked. Its reconstruction used PAML and a broader
vertebrate panel; this 18-taxon SLAC experiment is a comparison with the reported
protein, rather than an exact reproduction of the original analysis. Agreement
at ancestral states does not establish substitution order or epistasis. The
existing epistaeon benchmark and its plurality inference have not been changed.

## Artifacts and reproduction

- `inputs/`: added accession records and their hashes, 18-taxon alignments,
  MAFFT residue map, species scaffold and preparation checks.
- `inference_protocol.json`: command, versions, executable/source/library
  hashes and input checksums frozen before the fit.
- `results/slac.json`: original SLAC result, with all node states and model fits.
- `results/Boreoeutheria.codon.fa` and `.protein.fa`: the selected ancestor.
- `results/ancestors.codon.fa` and `.protein.fa`: SLAC-exported tips and internal
  nodes, preserving the original codon columns.
- `results/inference_script.py`: exact script that launched the fit.
- `results/postprocessing_script.py`: corrected export script snapshot.
- `results/inference_seal.json`: inference file checksums verified before and
  after the separate empirical comparison.
- [comparison/README.md](comparison/README.md): focal comparison, full
  per-residue TSV, numerical JSON and comparison hashes.

From the repository root, with the recorded HyPhy, MAFFT and Biopython 1.85
installed, create a new directory for an independent rerun:

```sh
sws1_trial=$(mktemp -d)
mkdir -p "$sws1_trial/inputs"
cp -R epistaeon/validation/sws1_slac/inputs/raw "$sws1_trial/inputs/"
cp epistaeon/validation/sws1_slac/inputs/outgroup_downloads.json "$sws1_trial/inputs/"
python3 epistaeon/scripts/run_sws1_slac_test.py --directory "$sws1_trial"
python3 epistaeon/scripts/compare_sws1_slac.py \
  --run-dir "$sws1_trial" \
  --alignment "$sws1_trial/inputs/SWS1.outgroups.codon.fa" \
  --ancestors-codon "$sws1_trial/results/ancestors.codon.fa" \
  --ancestors-protein "$sws1_trial/results/ancestors.protein.fa" \
  --selected-node "$sws1_trial/results/selected_node.json" \
  --inference-seal "$sws1_trial/results/inference_seal.json" \
  --paper-reference epistaeon/validation/sws1_slac_reference/source_reference.json \
  --substitution-reference epistaeon/data/substitution_reference/sws1.json \
  --output "$sws1_trial/comparison"
```

The inference script accepts no empirical sequence/reference argument. The
comparison runs only after inference files are sealed. Reproduction of the
source transcription and coordinate mapping can be checked independently:

```sh
python3 epistaeon/validation/sws1_slac_reference/verify_source_reference.py
python3 -m pytest epistaeon/tests/test_compare_sws1_slac.py -q
```

The local HyPhy package is Bioconda `hyphy-2.5.96-hadb5e62_0` for osx-arm64;
its package SHA256 is `838c592de1f9e57057034b2b6af7b4c445e2ea5d68b44c43b3a3da0fb0611d20`.
The saved protocol records the actual executable and batch-library hashes.
SLAC's node sequence fields are documented in the
[official HyPhy JSON reference](https://hyphy.org/resources/json-fields/#slac).
