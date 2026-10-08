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

### The one system with measured order

The resurrected corticoid receptors are the only case here where order was
*experimentally measured*, by building intermediate genotypes and assaying
them. Everything else in the repository has site-level ground truth only.

Receptor lineage: `AncCR → AncGR1 → AncGR2` (MR-like → still MR-like →
cortisol-specific). All mutagenesis was done **on the AncGR1 background**.

| Set | Substitutions | Role |
| --- | --- | --- |
| **X** | S106P, L111Q | adaptive, switches specificity |
| **Y** | L29M, F98I, S212del | completes the switch but destabilises |
| **Z** | N26T, Q105L | **permissive** for Y |
| Y27R | (earlier interval, AncCR→AncGR1) | permissive, predates the others |

Measured genotypes (Ortlund 2007, Fig. 3A):

```
AncGR1            functional     AncGR1+Y        NON-FUNCTIONAL
AncGR1+X          functional     AncGR1+X+Y      NON-FUNCTIONAL
AncGR1+Z          functional     AncGR1+X+Y+Z    functional, fully GR-like
```

So **Z and X must precede Y**, and Y27R precedes all of them. Bridgham 2009
adds the reverse direction: from AncGR2, reversing X gives a non-functional
receptor, while reversing the restrictive substitutions first restores nothing.

Numbering throughout is AncCR local LBD numbering, as both papers use.

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

**Two consequences of that, both verified numerically:**

1. Embeddings are *additive* across sites, so any genotype's embedding costs
   `1 + 2K` forward passes instead of one each (measured error ~1e-7).
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
  --tree epistaeon/data/alignments/toga2/speciesTree.nh \
  --ancestor anc.fa --descendant der.fa \
  --focal-taxon hg38 \
  --sites sites.txt --blocks blocks.json \
  --mode exact --scorer typicality-gate \
  --weights epistaeon/data/model/model.safetensors --cpu \
  --focal-unit 1 --permissive-units 0,2 \
  -o out/
```

`--ancestor` and `--descendant` must be **aligned to each other** (equal
length, gaps allowed). They may be a domain; the CLI maps them onto full-length
alignment columns through the focal row.

**Modes** (`--mode`, auto-selected by n):

| Mode | When | Notes |
| --- | --- | --- |
| `exact` | n ≤ 22 units | exact MAP, Z(M) and all marginals; 2^n memory |
| `blocks` | large n, grouped | `--blocks` groups substitutions; the receptor's X/Y/Z is n=3 |
| `sample` | n > 22 ungrouped | importance-sampled marginals with an ESS diagnostic |

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
- **`provenance.additivity_ok`** — gates the fast embedding path. The run aborts
  if it fails.
- **`provenance.per_unit_relative_delta`** — how far one substitution moves the
  embedding, against float32 precision.
- **`relative_branch_position_monotone`** — φ can overshoot, so a "fraction"
  may fall outside [0, 1]. Read those as φ levels, not times.

### Tests

```bash
cd epistaeon && python3 -m pytest tests/ -q      # 18 tests
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

- **The ancestral background is a chimera (bug).** `_relattice` in
  `epistaeon/src/epistaeon/cli.py` builds its ancestral string from the focal
  row and overwrites only the substituted columns, so the background stays
  ~73% human rather than AncGR1 (66 of 248 LBD positions differ). The whole
  claim under test is background dependence, so results are provisional until
  the full ancestral domain is written into the focal row.
- **No scorer yet recovers the measured receptor order** (Z and X before Y).
  `potts-metropolis` degenerates to uniform because the co-selection network
  found no coupling among those units. Treat the scorer comparison as
  provisional pending the fix above.
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
