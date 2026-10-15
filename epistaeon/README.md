# epistaeon

Infers **the order in which substitutions arrived** along one phylogenetic
branch, and whether some had to come before others. Given an ancestral and a
derived row of a codon alignment, it enumerates the orderings of the codon
substitutions between them, scores each intermediate genotype with the
hyphaeon model, and reports which orderings pass only through plausible
intermediates.

This directory also holds the empirical ground truth the engine is tested
against: the benchmark coordinate set, the published studies, and the scripts
that validate both. See [`data/README.md`](data/README.md) for the full audit
and [`AGENTS.md`](AGENTS.md) for how to run the engine.

The [source-verified substitution reference](data/substitution_reference/README.md)
lists expected events, partial orders, modern structural locations and the
remaining inputs needed for validation. It corrects the older study registries.

Case 4 uses **SWS1 cone opsin** from Yokoyama et al. (2014),
[doi:10.1371/journal.pgen.1004884](https://doi.org/10.1371/journal.pgen.1004884).
The [rebuilt codon-aware alignment](data/alignments/sws1/SWS1.codonified.fa)
contains extant boreoeutherian sequences from the paper's accession lists.
The [SWS1 reference](data/substitution_reference/sws1.json) keeps the eight
phylogenetically favored paths separate from the paper's reported 1,032
experimentally accessible orders. Table-transcription anomalies and a path-count
discrepancy remain explicit in the reference. No ancestral sequence is included
in this comparative input.

```bash
export PYTHONPATH=../aeon-core/src:../hyphaeon/src:src
python3 -m epistaeon.cli order --alignment <codon.fa> --use-tn93 \
  --focal-taxon <derived row> --ancestor <ancestral row> \
  --sites sites.txt --scorer typicality-gate \
  --weights data/model/model.safetensors --cpu -o out/
```

---

> # ⚠️ READ BEFORE RUNNING ANY STRUCTURAL ANALYSIS
>
> ## PDB author numbering is **not** alignment numbering
>
> Two of the four benchmark structures, 2HHB and 1MBO, are mature proteins
> whose initiator methionine was cleaved. Their PDB author numbering is therefore shifted by
> **one residue** from the UniProt-canonical numbering that alignment columns
> and site tensors use.
>
> ```
> uniprot_position = pdb_author_position + offset
> ```
>
> | Structure | Chain | Protein | UniProt | **Offset** |
> | --- | --- | --- | --- | --- |
> | `2HHB` | A, C | α-globin | P69905 | **+1** |
> | `2HHB` | B, D | β-globin | P68871 | **+1** |
> | `1MBO` | A | myoglobin | P02185 | **+1** |
> | `1U19` | A, B | bovine rod-opsin template for SWS1 | P02699 | **0** to its own sequence only |
> | `7PRX` | A | glucocorticoid receptor | P04150 | **0** |
>
> Skipping this shifts **every** hemoglobin and myoglobin contact call by one
> residue. It fails silently — no error, no gap, no mismatch — because the
> off-by-one lands on a real neighbouring residue. It is the single most
> likely way to get plausible-looking but wrong contact classifications.
>
> **Case 3 has a second trap.** The white paper names receptor residues in
> *ancestral* (AncCR) numbering — Tyr27, Ser106, Leu111. 7PRX uses human GR
> numbering, so add **+531** (Ser106 → GR 637, Leu111 → GR 642). The backup
> ancestral structures in `data/ancestral_coordinates/` use local LBD
> numbering (`-2` … `247`), and in 3GN8 the +531 rule breaks after about
> position 212. Deleted S212 itself has no modern residue coordinate; the
> source 211/212 boundary needs curated alignment.
>
> `1U19` carries an `ACE` acetyl cap at author position `0` that must be
> skipped when iterating residues. It is the homology template used by the
> SWS1 study, not a SWS1 structure. Apply the explicit source-site crosswalk
> before locating SWS1 substitutions; source labels are not direct current
> human UniProt positions or raw PDB numbers.
>
> **Alignment side, too.** The TOGA2 model inputs in
> `data/alignments/toga2/` report *alignment columns*, and the NR3C1 reference
> is the GRγ isoform (+1 Arg at 452), so human positions after 451 are one
> higher than 7PRX. See [`data/alignments/toga2/README.md`](data/alignments/toga2/README.md).
>
> **Do not hardcode these.** Load
> [`data/numbering_offsets.json`](data/numbering_offsets.json), the
> machine-readable source of truth. `python scripts/numbering.py` verifies the
> modern chain offsets; ancestral deletion exceptions require source curation.

---

## Scoring methods

Orderings are ranked by a product of **edge weights**: one weight per step,
`w(u | S)` for applying substitution unit `u` to genotype `S`. The inference
layer normalises over whole paths, so a weight need not be a per-step
probability; `-inf` blocks a step entirely. `--t-sel` sets the temperature `T`
(default 0.5), and `--scorer` picks the rule.

Each genotype is embedded by replacing the derived taxon's sequence in the
alignment, recomputing the TN93 distance matrix from that modified alignment,
and taking the host row's mean-over-sites hidden state, `z(S)`. Two scorers
also use `φ(S) = α‖z(S) − anchor‖`, the scalar divergence from a fixed anchor
(the ancestral genotype, or the model's `[ROOT]` token when no `--ancestor` is
given). **`α` defaults to 1, uncalibrated**, so φ is in raw latent units.

One constraint shapes all four. For any state function `f`, the increments
`Σ Δf` along a path telescope to `f(derived) − f(ancestor)` — the same for
every ordering. A scorer that is a plain exponential of such a difference
therefore ranks every ordering equally, and the engine will correctly report a
uniform distribution. Discrimination has to come from somewhere else:

| Scorer | Edge weight | Where discrimination comes from |
| --- | --- | --- |
| **`potts-metropolis`** (default) | `min(0, ΔE/T)` | the `min()`, and couplings making `ΔE` depend on the background |
| **`typicality-gate`** | `ΔV/T + V(S′)/T` | the absolute plausibility of the state entered |
| **`softmax-repaired`** | per-step log-softmax of `Δφ/T` | the normalisation changing from step to step |
| **`dphi-product`** | `log(Δφ/T)` | nowhere — a documented baseline |

### `potts-metropolis` (default)

A Potts energy over applied units, `E(S) = Σ_u f_u + Σ_{u<v} J_uv`, accepted
Metropolis-style. The site terms `f_u` are the model's own output: the change
in predicted selection evidence (LRT) at the unit's columns when its codons
switch from ancestral to derived. The couplings `J_uv` come from hyphaeon's
attribution co-selection network (CESI).

The `min(0, ·)` makes a deleterious step cost something wherever it falls in
the order, and `J` is what makes the cost depend on which substitutions are
already present — the epistasis.

Two things to know before trusting it:

- **CESI is a correlation statistic** over attribution patterns, not a fitted
  coupling constant. It is not learned site–site coupling.
- **If the network finds no coupling among your units**, `E` is additive, `ΔE`
  stops depending on the background, and every ordering ties. The run says so
  (`discrimination.degenerate_uniform`, `unit_pairs_with_coupling`). Check it
  every time.

### `typicality-gate`

Scores how much like a real protein each intermediate looks. Viability `V(S)`
is the negative Mahalanobis distance from `z(S)` to the cloud of extant taxon
embeddings, with a diagonally shrunk covariance.

The weight adds the change in viability to **the absolute viability of the
state being entered**, which is what breaks the telescoping: an ordering forced
through an atypical intermediate is penalised directly, not merely through a
difference. `--viability-floor` makes that a hard wall — intermediates below it
are blocked outright, expressing "selection does not traverse dead proteins".

### `softmax-repaired`

The white paper's rule: at each step, the chosen move's share of a softmax over
`Δφ/T` against the moves it competes with. Because the competing set shrinks as
the trajectory proceeds, the normalisation varies by step, which is what lets
it rank orderings at all.

**Caveat, as currently wired:** the competing set is only the units not yet
applied, so at the final step nothing competes and that step contributes
weight 1 regardless of its score — the defect the rule is supposed to fix. See
the open issue in [`AGENTS.md`](AGENTS.md).

### `dphi-product`

`log(Δφ/T)` per step, with non-positive increments blocked. **A baseline for
comparison, not a method to report.** Its temperature is inert — every ordering
has the same number of edges, so `T` scales every path identically — and what
it actually rewards is evenly sized steps, since the sum of `Δφ` is
path-independent and a product under a fixed sum is maximised when the terms
are equal.

## Layout

| Path | Contents |
| --- | --- |
| `data/coordinates/` | Active globin/receptor structures `2HHB`, `1MBO`, `7PRX`, plus SWS1 homolog template `1U19` |
| `data/alternate_coordinates/` | Backup: `4LSJ` (human GR; domain-swapped, superseded by 7PRX) |
| `data/ancestral_coordinates/` | Backup: Thornton ancestral receptor series (`2Q1H`, `2Q1V`, `2Q3Y`, `3RY9`, `3GN8`) |
| `data/experimental/` | Published source papers and supplements for the four cases |
| `data/uniprot/` | Canonical reference sequences, so the numbering check runs offline |
| `data/alignments/toga2/` | TOGA2 codon alignments + species tree for the globin/receptor cases |
| `data/alignments/sws1/` | Rebuilt codon-aware extant boreoeutherian SWS1 alignment and accession provenance |
| `data/numbering_offsets.json` | **Verified numbering crosswalk — load this, don't hardcode** |
| `data/substitution_reference/` | Source-verified events, ordered sets, structural locations and validation gaps |
| `data/README.md` | Full structure audit and contact-mapping fitness report |
| `data/model/` | Model checkpoint (gitignored; see `AGENTS.md` for the fetch) |
| `src/epistaeon/` | The engine: `lattice`, `scoring`, `trajectory`, `background`, `timing`, `cli` |
| `tests/` | 33 tests, incl. exact-DP-vs-brute-force and the inert-scorer regression |
| `validation/` | Replication harness: replicator, adversary, summarizer |
| `scripts/` | Reproduces every number in `data/README.md` |

## Scripts

Each runs standalone with no arguments (except `validate_cif.py`, which takes
mmCIF paths). They need `biopython`, `numpy` and `scipy`.

| Script | Purpose |
| --- | --- |
| `numbering.py` | Verifies the offsets in the box above against UniProt |
| `validate_cif.py` | Entity, chain, organism, ligand and gap inventory |
| `probe.py` | Residue identity at the positions the white paper names |
| `align_ladder.py` | AncCR → AncGR1 → AncGR2 substitution ladder; derives the +531 crosswalk |
| `contacts.py` | Reproduces the §4 Table 2 distance and contact-count claims |
| `storz.py` | Maps the real Storz 2009 deer-mouse sites onto 2HHB |
| `mapfitness.py` | Distance-class composition, copy stability, interface ambiguity, cofactor bridging |
| `bands.py` | Band coverage gaps and statistical power for the §5 benchmark |
| `interface_states.py` | Hemoglobin interface disjointness and T/R state dependence |
| `_fetch.py` | Resolves a structure locally, else caches a download from RCSB |

## Contact classification

Binary: a detected pair is in direct contact if any heavy atoms are within
5.0 Å, otherwise not. Apply the numbering offset first, exclude pairs with
|i − j| < 5, and count hemoglobin α–β pairs as contacts at either interface.
Baselines and details in [`data/README.md`](data/README.md).
