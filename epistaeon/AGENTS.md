# AGENTS.md

Orientation for agents working in this repository. The `epistaeon/` work is the
active project; `hyphaeon/`, `chronaeon/` and `aeon-core/` are the existing
packages it builds on.

---

## 1. The biology we are testing

When two homologous proteins differ by several substitutions along one
phylogenetic branch, the substitutions did not arrive simultaneously. **We are
trying to infer the order in which they occurred**, and specifically whether
some had to come before others.

Three archetypes matter:

- **Adaptive** — changes function directly. Often destabilising on the
  ancestral background.
- **Permissive** — functionally silent on its own, but stabilises the protein
  so that a later adaptive change is tolerated. Must come *first*.
- **Deleterious-in-isolation** — would complete a functional shift, but kills
  the protein unless the permissive changes are already present.

Order matters because a trajectory passing through a non-functional
intermediate is historically unlikely: selection does not traverse dead
proteins. So the question "which order?" is really "which orders pass only
through viable intermediates?"

---

## 2. Facts you must not get wrong

These have each caused a real error in this project. They fail silently.

**Numbering is three different coordinate systems.** Alignment columns, PDB
author numbering and UniProt numbering do not coincide. Load
`epistaeon/data/numbering_offsets.json`; never hardcode an offset.

- 2HHB and 1MBO: UniProt = PDB **+1** (cleaved initiator Met)
- 1U19, 7PRX: offset **0**
- AncCR → human GR: **+531** up to position 210, then **+530** — the S212
  deletion shifts it. A single constant is wrong.
- The TOGA2 **NR3C1 alignment is the GRγ isoform**: an extra Arg at 452 means
  positions after 451 are one higher than UniProt/7PRX.
- Model output is **alignment columns**, not residue numbers. Drop the columns
  where the human reference has a gap first.

**S212 is a deletion, not a substitution.** Ortlund's footnote 21: the deletion
replaces Ser at codon 212. It appears as a gap, so a site-based detector may
have no position to report. Do not score it as a miss.

**The white paper is not a source of truth.** `epistaeon/white_paper_*.pdf` was
AI-written and several claims are fabricated (`Thr36` does not exist; the
α57/α64 polarity is inverted; two cited distances do not reproduce). Do not
read it as ground truth, do not cite it, and do not test anything against it.
The original studies in `epistaeon/data/experimental/` are the only authority.

**The model scores each codon site independently.** `window_size: 1`, no
cross-site layers, and the batch dimension is the site. Therefore **no part of
this pipeline conditions a site on the rest of its own sequence**. If a result
appears to show background dependence, establish where it came from.

**Two consequences of that:**

1. A genotype's own residues cannot inform each other. What distinguishes two
   genotypes' forward passes is **where each sits among the other taxa**.

   **A genotype replaces the derived taxon's sequence**; it is never appended
   as an extra row. An extra row is an extra *taxon*, and the model reads taxa
   as evidence of a realized evolutionary process — so appending a genotype
   would feed the embedding of every site evidence for a sequence that may
   never have been viable. The alignment therefore keeps exactly its own taxa,
   and the lineage under study occupies the row it already has.

   Because that row's sequence changes, **the distance matrix is recomputed
   from the modified alignment** (TN93 for the host against every other taxon)
   and its MDS frame with it. Only the host's row and column can change, since
   every other pair is untouched and TN93 is pairwise; the patch is verified
   against a full recompute each run (`matrix_patch_max_error`, 0.0 observed).

   Each genotype costs one full forward pass, so exact mode over `n` units
   costs `2^n` — the binding cost of a run, and what `--max-genotypes` bounds.
   At the derived end the host row holds its own sequence, so the alignment,
   its matrix and its frame are the originals and the embedding is that
   taxon's own (verified identical). Measured on `HLmusEve1 → hg38`, 12 taxa:

   | what | relative ‖dz‖ |
   | --- | --- |
   | one substitution | 0.0030 – 0.0048 |
   | the recomputed frame, on a bystander taxon | 0.0001 |

   The frame moves ~30× less than the signal, so recomputing it per genotype
   is benign; `provenance.background_frame_shift` reports it per run. Check it
   against `per_unit_relative_delta` before trusting an ordering.
2. For any state function φ, `Σ Δφ` along a path telescopes to
   `φ(derived) − φ(ancestor)` — **identical for every ordering**. A scorer
   built on plain `exp(Δφ/T)` gives a uniform distribution and a before-matrix
   of 0.5 everywhere. Path discrimination must come from a non-linear
   combination rule (Metropolis, viability gate, bottleneck) or from varying
   normalisation.

---

## 3. Running the analyses

### Setup

```bash
# the packages are not installed; use the source tree
export PYTHONPATH=aeon-core/src:hyphaeon/src:epistaeon/src
```

Model weights are **gitignored** (7.7 MB). If `epistaeon/data/model/` is empty:

```bash
curl -sL -o epistaeon/data/model/model.safetensors \
  https://huggingface.co/datamonkey/hyphaeon/resolve/main/model.safetensors
```

### Infer substitution order

```bash
python3 -m epistaeon.cli order \
  --alignment epistaeon/data/alignments/toga2/NR3C1.codonified.fa \
  --use-tn93 \
  --focal-taxon AncGR2 \
  --ancestor AncGR1 \
  --sites sites.txt \
  --blocks blocks.json \
  --mode exact \
  --scorer typicality-gate \
  --weights epistaeon/data/model/model.safetensors \
  --cpu \
  --focal-unit 1 \
  --permissive-units 0,2 \
  -o out/
```

**Both endpoints are rows of the alignment**, named not supplied as FASTA.
`--focal-taxon` (aliases `--descendant`, `--derived`) names the derived row;
`--ancestor` names the ancestral one. Nothing is realigned, so no numbering
mapping is involved. `--max-species` **keeps both endpoints** whatever Faith's
PD would have chosen — without that, downsampling drops exactly the
well-sampled taxa (`hg38` among them) that make natural endpoints.

**`--focal-taxon` must be one of the analysed taxa**, because every genotype
replaces that row; a name that duplicate pruning collapsed resolves onto the
identical row that was kept. `--ancestor` need only be *in* the alignment: it
supplies codons, and if it is not an analysed taxon it simply is not part of
the background (`provenance.ancestor_in_background`).

**Omit `--ancestor`** and each ancestral state becomes the **plurality residue
of its alignment column** over the other rows, with φ anchored at the model's
`[ROOT]` embedding instead of an ancestral sequence. That is a consensus, not a
reconstruction: it can differ from the real ancestor at any column, so check
`provenance.plurality_support` before reading the result as history. To use a
reconstruction such as AncGR1, add it to the alignment as a row and name it.

**`order` requires `--use-tn93`**, and refuses a tree. Each genotype's
distances are recomputed from its sequence, and **a tree has no branch for a
sequence that is not in it** — that alone rules out the tree path, whatever the
scale.

Patristic and TN93 distances are both **substitutions per site**, so they are
the same kind of quantity, but they are not the same matrix here.
`speciesTree.nh` is one tree shipped once for all seven gene alignments, so its
branch lengths are a genome-wide scale, while TN93 measures the gene in front
of you — and NR3C1 is conserved. Measured over 40 taxa:

| | median | max |
| --- | --- | --- |
| patristic, from `speciesTree.nh` | 0.759 | 1.441 |
| TN93, from `NR3C1.codonified.fa` | 0.114 | 0.178 |

They correlate at r = 0.84 (Spearman 0.84), but **no single factor reconciles
them**: the best scalar fit is `tn93 ≈ 0.137 × patristic` at R² = 0.46, and the
ratio falls from 0.168 for the closest quartile of pairs to 0.123 for the
farthest, which is TN93 saturating. So a TN93 row patched into a patristic
matrix would read as far closer to every taxon than any real sequence is. The
other commands still take a tree.

**A substitution is a codon change, synonymous ones included**, so applying
every unit reproduces the derived row **exactly, nucleotide for nucleotide**
(verified both ways: TN93 0.0 from mask 0 to the ancestor row and from the full
mask to the derived row). Ordering only the amino-acid changes would leave the
end of the trajectory short of the derived sequence at every synonymous column
— 126 of them for `HLmusEve1 → hg38`, two thirds of the branch's nucleotide
distance — and TN93 measures nucleotides, so that gap would displace every
genotype. Synonymous units are named `K4=AAG>AAA`; `provenance` carries the
`substitution_classes` counts and each member's codons.

The consequence for `n`: a full mammalian branch is ~190 codon substitutions
(`HLmusEve1 → hg38`: 59 nonsynonymous, 126 synonymous, 9 indel), far past the
exact limit. **`--sites` or `--blocks` is not optional** on a real branch.

**Positions** in `--sites` and `--blocks` are 1-based residues of the ancestor
row, or of the derived row when there is no ancestor. Use `c<N>` for alignment
column `N` where that row is gapped — which is how a deletion's column is
named, since the row has no residue to number there.

### What `--blocks` does

`--sites` picks *which* substitutions to consider; `--blocks` groups them into
**units that move together**, and a unit is what gets ordered. The DP orders
`n` units, so grouping is what keeps `n` small enough to enumerate: the
receptor's X/Y/Z is 3 units covering 7 substitutions.

The file is JSON mapping a group name to positions, numbered exactly as
`--sites` is:

```json
{
  "X": [106, 111],
  "Y": [29, 98, 212],
  "Z": [26, 105]
}
```

- The group **name** is what appears in the MAP order, the before/after matrix
  and `--focal-unit`/`--permissive-units`, which take **unit indices** in the
  order the groups are read, singletons last.
- Any substitution not named in a group becomes **its own unit**, so nothing is
  silently dropped. With 194 substitutions and three groups covering 7 of them
  you get 3 + 187 = 190 units, not 3 — pair `--blocks` with `--sites` to drop
  the rest.
- A position named in a group that is not a substitution between these two
  endpoints is an error naming your own position, as is one that appears in two
  groups.
- Substitutions inside a unit are applied simultaneously and are never ordered
  against each other.

**Modes** (`--mode`, auto-selected by n):

| Mode | When | Notes |
| --- | --- | --- |
| `exact` | n ≤ 22 units **and** 2^n ≤ `--max-genotypes` (default 4096) | exact MAP, Z(M) and all marginals; 2^n forward passes, one TN93 call |
| `blocks` | large n, grouped | `--blocks` groups substitutions; the receptor's X/Y/Z is n=3 |
| `sample` | n beyond either limit | importance-sampled marginals with an ESS diagnostic; one TN93 call per genotype, still capped by `--max-genotypes` |

**Scorers** (`--scorer`): `potts-metropolis` (default), `typicality-gate`,
`softmax-repaired`, `dphi-product` (documented baseline only — its temperature
is inert and it rewards evenly sized steps).

### Reading the output

`out/order.json` carries the MAP order, its raw and **normalised** probability,
`Z(M)`, the before/after matrix `C`, the position matrix `R`, designated
marginals such as `P(Y after ≥ k of {X, Z})`, and a `provenance` block.

Check these before believing anything:

- **`discrimination.degenerate_uniform`** — if true, every ordering had equal
  weight and the scorer told you nothing. Compare `map_prob_normalised` against
  `uniform_baseline` (= 1/n!) every time.
- **`provenance.ancestral_states`** and **`plurality_support`** — whether the
  ancestral background is a real row or a per-column consensus, and how thin
  that consensus is. A support near 0.5 is a coin toss, not an ancestor.
- **`provenance.endpoints_in_background`** — false means that endpoint supplied
  codons without being one of the analysed taxa.
- **`provenance.forward_passes`** — one per distinct genotype; the run's cost.
  `genotypes_read_from_a_real_row` counts the genotypes that turned out to be
  the host's own sequence and needed no pass.
- **`provenance.background_frame_shift`** vs **`per_unit_relative_delta`** —
  the recomputed-frame wobble against the substitution signal. See §2.
- **`provenance.matrix_patch_max_error`** — 0.0 unless patching the host's row
  and column ever stops matching a full TN93 recompute.
- **`provenance.per_unit_relative_delta`** — how far one substitution moves the
  embedding, against float32 precision.
- **`relative_branch_position_monotone`** — φ can overshoot, so a "fraction"
  may fall outside [0, 1]. Read those as φ levels, not times.

### Tests

```bash
cd epistaeon && python3 -m pytest tests/ -q      # 30 tests
```

The suite includes exact-DP-versus-brute-force for n = 3–6 and an
**inert-scorer regression** — if you add a scorer, it must break
path-independence, and that test is what proves it.

### Validation harness

Separate from the engine: three agents replicate published studies, attack the
replication, and summarise. See `epistaeon/validation/README.md`.

```bash
/validate-epistaeon                                       # full flow
python3 epistaeon/validation/compact_reports.py --compact  # merge + validate
```

The harness is scoped to two axes only: **which epistatic sites** are found,
and **the inferred order** where a study measured it. `epistaeon`'s internal
machinery is not under test. `not_testable` is a first-class verdict and must
never be folded into a performance score.

---

## 4. Known open issues

- **The chimeric ancestral background is fixed.** Endpoints are now alignment
  rows, so the ancestral background is that whole row rather than the focal row
  with the substituted columns overwritten. Any result produced before that
  change (`_relattice`, since removed) is void.
- **No scorer yet recovers the measured receptor order** (Z and X before Y).
  `potts-metropolis` degenerates to uniform because the co-selection network
  found no coupling among those units. The scorer comparison has not been
  re-run since endpoints became alignment rows.
- **The ancestral reconstructions are not in the TOGA2 alignment.** AncGR1 and
  AncGR2 live in `epistaeon/data/sequences/ancestral/` as protein sequences,
  while endpoints must now be alignment rows, so the receptor analysis needs
  them codon-aligned into the alignment first. Column plurality is not a
  substitute for them.
- **Natarajan 2013** (*Science*, deer mouse) is the only measured hemoglobin
  epistasis and is not yet in the harness registry.
- **3RY9 is mislabelled.** It matches AncGR1.1 at 99.2% and AncGR1 at 92.0%, so
  it is the alternative reconstruction despite its PDB title.
- Stage 5 timing is relative only; absolute times need `--ne`/`--s` and raise
  when `2·Ne·s ≤ 1`, where the published formula returns a negative time.

## 5. Where things are

| Path | Contents |
| --- | --- |
| `epistaeon/src/epistaeon/` | the engine: `lattice`, `scoring`, `trajectory`, `background`, `timing`, `cli` |
| `epistaeon/IMPLEMENTATION_PLAN.md` | design rationale and the verification list |
| `epistaeon/data/README.md` | structure audit: what each PDB file is, and what reproduces |
| `epistaeon/data/alignments/toga2/README.md` | alignment provenance and 14 known issues |
| `epistaeon/data/numbering_offsets.json` | **the** numbering authority |
| `epistaeon/data/sequences/ancestral/` | AncCR, AncGR1, AncGR1.1, AncGR2 + posterior probabilities |
| `epistaeon/validation/` | the replication harness, registries and schemas |
| `epistaeon/scripts/` | standalone structure-validation scripts |
