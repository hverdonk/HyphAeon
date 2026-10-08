# Implement a trajectory-scoring engine (`epistaeon`) for ancestor→derived paths

## Context

The goal is to compare the likelihood of all mutational trajectories from an
ancestral sequence to a derived one, such that **less likely paths reflect poor
or inaccessible intermediate states**. The white paper's §2 is one proposal for
this, treated here as a starting idea rather than a specification, because two
of its components do not survive inspection:

- `Φ` is defined as an energy "derived from root-to-leaf attribution flow and
  cross-attention log-likelihoods", but the published checkpoint has a single
  ordinal-LRT head and no likelihood output.
- The encoder scores each codon column independently (`window_size: 1`, no
  cross-site layers), so **no per-site score can depend on the rest of the
  sequence**. Sequence context must come from pooled embeddings, explicit
  pairwise couplings, or an external model.

A third constraint is mathematical and binds every design: for any state
function φ, `Σ Δφ` along a path telescopes to `φ(derived) − φ(ancestor)`, which
is identical for every path. Verified numerically — summing `Δφ/T_sel` gives
spread 0.0000 across orderings, and unnormalised `Π exp(Δφ/T)` scores every
path the same. Discrimination must come from a rule that reads absolute state
levels (Metropolis, viability gate, bottleneck) or from varying normalisation.

The work therefore splits cleanly into a **scoring layer** (supplies a
background-dependent, non-telescoping `p(m|S)`) and an **inference layer**
(exact DP over the subset lattice). They are complements: the DP is provably
correct but *inert* if the scorer telescopes — with pure `exp(ΔE/T)` weights the
DP returns a uniform distribution over orderings and a `C` matrix of 0.5
everywhere.

## Design decisions (settled)

| Fork | Decision |
| --- | --- |
| φ anchor | Embedding of the **actual ancestral sequence** inserted as a taxon; `[ROOT]` token embedding as fallback |
| Stage 5 timing | Relative position along the branch only; absolute drift/sweep behind `--ne`/`--s`, raising rather than returning a negative time when `2·Ne·s ≤ 1` |
| Site selection | Auto-derive K from the ancestor/derived **codon** diff, synonymous changes included so the trajectory ends on the derived row exactly, **with block mode and an explicit site-subset flag also available** (a real branch is ~190 units, so one of them is required) |
| Scoring | Pluggable scorers behind one interface; default `potts-metropolis` |
| Trajectory inference | Exact subset DP; three selectable modes (below), auto-chosen by n with manual override |

## Package layout

New sibling package mirroring `chronaeon/` and `hyphaeon/`:

```
epistaeon/pyproject.toml           name="epistaeon", entry point epistaeon = "epistaeon.cli:main"
epistaeon/src/epistaeon/
    background.py   load alignment+model; insert a genotype into a focal taxon row
    lattice.py      derive substitutions; genotype = bitmask; closed-form embeddings; blocking
    scoring.py      scorer interface + four implementations
    trajectory.py   exact subset DP (Viterbi + sum-product + forward-backward); sampling fallback
    timing.py       relative branch position; guarded absolute conversion
    cli.py io.py
epistaeon/tests/
```

Existing `epistaeon/data/`, `scripts/` and `validation/` stay put.

## Reuse (do not reimplement)

- `aeon_core.dataset.load_alignment_and_tree` → `(c, a, d, z, invariable, taxa, L)`
- `aeon_core.inference.get_device`, `load_model`; `aeon_core.weights` for the HF checkpoint
- `aeon_core.splits.extract_cross_taxa_attentions_and_embeddings` → per-taxon embeddings
- `aeon_core.dataset.GENETIC_CODE`, `AA_MAP`, `CODON_TO_AA`; `hyphaeon.epistasis.CANONICAL_AA_TO_CODON`
- Substitution mechanics: the in-place pattern in `run_insilico_selection_dms`
  (`hyphaeon/src/hyphaeon/epistasis.py`), writing `c_tensor[site, focal_idx, 0]`
- Pairwise couplings: `compute_branch_coselection_network` (same file) gives CESI site-pair scores
- `aeon_core.io.write_json`, `write_csv`; numbering from `epistaeon/data/numbering_offsets.json`
  (never hardcode: AncCR→human GR is +531 to position 210, +530 from 212)

## Core mechanic — embeddings for any genotype

**Superseded.** The original plan exploited the encoder's lack of cross-site
mixing (`z = (1/L) Σ_s h_s`, each `h_s` depending only on column `s`) to get
every genotype's embedding from `1 + 2K` passes by arithmetic, gated on an
additivity test. That shortcut wrote the genotype into an existing taxon's row,
which gave it that taxon's phylogenetic position and left the ancestral
background a chimera of the focal row.

As implemented, each **intermediate** is instead its own row: appended to the
alignment, placed by TN93 distance to every taxon, with MDS coordinates
projected into the existing frame (Gower) so the real taxa and the `[ROOT]`
origin are identical across genotypes. Each intermediate's distance row is its
own, so nothing is shared between them and embeddings are **not** additive:

- one full forward pass per distinct intermediate, memoised by mask
- exact mode over `n` units costs `2^n` passes, bounded by `--max-genotypes`
- one TN93 call for the whole lattice in exact mode; one per genotype when
  sampling

The **endpoints are read from their own rows** rather than appended as copies
of themselves, so no genotype duplicates a real taxon. That leaves an endpoint
embedded among `N` taxa and an intermediate among `N + 1`, a frame difference
worth ~6x one substitution; it is measured per run as
`first_step_relative_delta` against `per_unit_relative_delta` and must be
checked before trusting an endpoint-adjacent step.

## Scoring layer (`scoring.py`)

Interface: `edge_weight(S, m) -> float > 0`, plus optional `state_score(S)`.
Weights need **not** be per-step-normalised — the DP normalises over paths.

1. **`potts-metropolis` (default).** `E(S) = Σ_s f_s(x_s) + Σ_{s<t} J_st(x_s,x_t)`;
   site terms `f` from per-site ΔLRT, couplings `J` from the co-selection
   network. Weight `min(1, exp(ΔE/T))`. Pairwise terms give real background
   dependence; deleterious steps suppress a path. Needs `J` sign/scale calibration.
2. **`typicality-gate`.** `V(S)` = Mahalanobis or k-NN distance of `z(S)` to the
   cloud of extant embeddings; weight combines `ΔV` with an absolute gate on
   `V(S∪{m})`, so atypical intermediates are penalised directly.
3. **`softmax-repaired`.** `softmax(Δφ/T_sel)` normalised over all
   single-residue moves at each step, times a viability gate.
4. **`dphi-product`.** Raw `Δφ/T` — baseline only; document that `T_sel` is
   inert (cancels as `T^-K`) and that it rewards evenly sized steps.

Zero weights are legal but handled explicitly: mask to `-inf` in log space, and
if `Z(M) = 0` report **"no accessible path"** rather than dividing by zero.

## Inference layer (`trajectory.py`)

Exact subset DP on the DAG whose nodes are subsets `S ⊆ M`, root `∅`, terminal
`M`, edges `S → S∪{m}` weighted `p(m|S)`. All in log space with log-sum-exp.

- **Viterbi / max-product**: `V(S) = max_{m∈S} [ V(S\{m}) + log p(m | S\{m}) ]`,
  with back-pointers; backtrack from `M` for the MAP order. `O(n·2ⁿ)`.
- **Forward sum-product**: `Z(S) = Σ_{m∈S} Z(S\{m}) · p(m | S\{m})`, giving
  `Z(M)` and hence `P(π) = P(π)/Z(M)` as a proper distribution over orderings.
- **Backward**: `B(S) = Σ_{m∉S} p(m|S) · B(S∪{m})`, `B(M) = 1`.
- **Edge marginals**: `Z(S)·p(m|S)·B(S∪{m}) / Z(M)`, from which:
  - `C_ij = P(m_i before m_j)` — sum over edges adding `j` with `i ∈ S`
  - `R_ik = P(m_i at step k)` — note `|S|` *is* the step index in this DAG, so
    position marginals need no extra state
  - `P(F after ≥k from permissive set P)` — sum over edges adding `F` where
    `|S ∩ P| ≥ k`. This maps directly onto the receptor ground truth.

**Three selectable modes** (`--mode`, auto-chosen by n, override allowed):

| Mode | Use | Cost |
| --- | --- | --- |
| `exact` | n ≤ 22 (default threshold) | `O(n·2ⁿ)`; 0.27 GB at n=24, 4.3 GB at n=28, infeasible at n=37 |
| `blocks` | group substitutions into ordered blocks and run the exact DP over blocks | trivial — the receptor's X/Y/Z/Y27R is n=4 |
| `sample` | n above threshold with no blocking | beam search for top paths + Monte-Carlo orderings, marginals with CIs |

`blocks` is the mode that matches what the experiments measured; `exact` on the
full auto-derived set is only reachable for shorter branches.

## CLI and outputs

```
epistaeon order --alignment <fa> [--tree <nwk> | --use-tn93] \
                --focal-taxon <row> [--ancestor <row>] \
                [--scorer potts-metropolis] [--mode auto] [--max-genotypes 4096] \
                [--blocks blocks.json] [--sites sites.txt] [--weights/--variant] -o out/
```

Both endpoints are **alignment rows**, named rather than supplied as FASTA
(`--focal-taxon` also accepts `--descendant`/`--derived`). With no `--ancestor`,
each ancestral state is its column's plurality residue and φ is anchored at
`[ROOT]`. Positions are 1-based residues of the ancestor row, or of the derived
row when there is none, with `c<N>` for a column that row is gapped at.

`order.json` must contain, at minimum:

1. MAP substitution (or block) order
2. its raw log-probability
3. its **normalised** probability `P(π_MAP | endpoint)`
4. the total path weight `Z(M)`
5. `C_ij` before/after matrix
6. `R_ik` position-probability matrix
7. designated marginals, e.g. `P(group Y after ≥k of {X, Z})`
8. `provenance`: model variant, scorer, mode, φ anchor, α, `T_sel`, distance
   source, both endpoint rows and whether each is in the background, how the
   ancestral states were obtained (row or plurality, with support), forward
   passes used, per-unit delta magnitudes vs float precision

`pairs.csv`: coupled sites, for the harness's site-identification targets.

## Verification

1. **DP vs brute force.** For n = 3–6, exact DP MAP, `Z(M)` and all path
   probabilities must match explicit enumeration of `n!` orderings. (Prototyped:
   agrees, and path probabilities sum to 1.000000.)
2. **Marginal sanity.** `R` row sums and column sums = 1; `C_ij + C_ji = 1` for
   `i ≠ j`.
3. **Inert-scorer regression.** Feeding pure `exp(ΔE/T)` weights must yield a
   uniform distribution and `C = 0.5` everywhere; the default scorer must *not*.
   This is the test that catches a telescoping scorer, the original design error.
4. **Genotype placement and endpoint exactness.** Re-projecting a point already
   in the MDS configuration returns its own coordinates; an endpoint resolves to
   its row whether or not it survived into the analysed taxa; and applying every
   unit reproduces the derived row nucleotide for nucleotide (TN93 0.0 at both
   ends of the branch, checked on two real ancestor rows).
5. **Zero-weight handling.** A blocked edge gives `-inf` not `NaN`; an all-blocked
   lattice reports "no accessible path" — and path weights too small to
   exponentiate must *not* be reported that way.
6. **Negative control.** Shuffle alignment columns within each taxon; inferred
   order must degrade to chance (`C → 0.5`).
7. **End-to-end on the receptor case.** `NR3C1.codonified.fa` + `speciesTree.nh`;
   ancestor AncGR1 (`steroid_receptor_2011_ancestral_sequences.doc`), descendant
   AncGR2 (`ancestral_coordinates/3GN8.cif`). Expected from
   `validation/studies/steroid_receptor.json`: Z (26, 105) and X (106, 111)
   before Y (29, 98, 212del); Y27R before all; AncGR1+X+Y non-functional while
   AncGR1+X and AncGR1+Z are functional. `S212del` is a gap, not a residue change.
8. **Scorer comparison.** Run all four scorers in `blocks` mode on the receptor
   case and report which recovers the measured order — the point of the
   pluggable interface is that the harness decides, not us.
9. **Harness**: `/validate-epistaeon`, then `compact_reports.py --compact`.

## Risks to report in the output, not paper over

- HyphAeon cannot condition a site on the rest of its own sequence. In
  `potts-metropolis` the epistasis comes from co-selection couplings; in
  `typicality-gate` it is geometric. Neither is learned site–site coupling.
- Per-site deltas scale as `1/L` (L = 866 for NR3C1), so single-substitution
  signal may sit near numerical noise; provenance reports the magnitudes.
- `P(π | endpoint)` is conditional on reaching the derived sequence by exactly
  these substitutions. It ranks orderings; calling a path "inaccessible" in
  absolute terms rests on the viability gate, which is a modelling choice.
- A normalised `P(π_MAP)` near `1/n!` means the model is not discriminating,
  not that all paths are equally viable. Report it next to the uniform baseline.
