# Blinded scorer benchmark

The benchmark separates prediction from empirical evaluation. A deterministic
container reconstructs an ancestor with HyPhy SLAC and runs inference with the
alignment, ASR tree/node, local checkpoint, checkpoint configuration and frozen
settings. Only after its outputs are sealed
does a separate evaluator read the empirical substitution reference. No
language-model agent participates in inference. The read-only evaluator agent
interprets completed results; it cannot run or modify the analysis.

## Run, seal, evaluate

Run from the repository root with Docker running. The first image build may
fetch the pinned base image, checksum-pinned HyPhy source and public package
wheels; reconstruction and inference have no network. The supplied checkpoint
must be local and use the safetensors format.
A nonempty `config.json` must be beside the checkpoint, including for a
checkpoint with embedded architecture metadata. Missing configuration must be
supplied before running rather than fetched inside the container.

```bash
export PYTHONPATH=aeon-core/src:hyphaeon/src:epistaeon/src
python3 -m epistaeon.cli benchmark-run \
  --alignment epistaeon/data/alignments/sws1/SWS1.codonified.fa \
  --focal-taxon Homo_sapiens \
  --asr-alignment epistaeon/validation/sws1_slac/inputs/SWS1.outgroups.codon.fa \
  --asr-tree epistaeon/validation/sws1_slac/inputs/species_scaffold.nwk \
  --asr-node Boreoeutheria \
  --weights epistaeon/data/model/model.safetensors \
  --output epistaeon/validation/runs/sws1-001
```

`benchmark-run` accepts no reference, manuscript, expected sites, blocks,
supplied ancestral sequence or expected order. The output directory must not
exist. It reconstructs once, freezes that result, and supplies identical
starting states to all four scoring methods. It records commands, logs, status, elapsed time,
input/source/runtime hashes, selected substitutions, starting codons, MAP order,
probability matrices and diagnostics. It writes a checksum seal before the
reference-aware stage can consume the results. A failed method is recorded as a
failure rather than quietly retried with a changed protocol.

The sealed run contains `preparation_protocol.json`, `protocol.json`,
`selection.json`, `asr/`, `asr_execution.json`, preparation logs, `inputs/`, the
allowlisted `source/` and `source_manifest.json`, copied runtime files,
`runtime_config.json`, `runtime_image.json`, `execution_manifest.json`,
`isolation_check.json`, `build.log`, `run_status.json` and `seal.json`.
Each `methods/<scorer>/` directory contains `status.json`, `stdout.log`,
`stderr.log` and, on successful inference, the engine's `order.json` and CSV
matrices. Per-method `selection.json` records the worker's independent input
preflight. Method commands, timing and exit codes are in `status.json`.
`asr/` contains the raw SLAC JSON, command/settings, executable and batch-library
hashes, input hashes, stdout/stderr, exit status and timing, fitted tree,
selected node and descendants, reconstructed codons and raw amino-acid states.
FASTA exports preserve ambiguous codons and use `X` for multi-state protein
cells. Reconstruction is stored outside scorer output mounts; its hash is
frozen in the final execution protocol before scoring starts.

`benchmark-run` exits `0` when every method completes and `1` when method
failures are recorded; input, setup or execution errors return `2`. Partial
outputs created after setup are sealed when worker termination is confirmed.
If interrupted-container cleanup cannot be confirmed, the run records
`cleanup_unconfirmed` and remains unsealed; the evaluator rejects it. A method failure can be
reported by the evaluator; a setup failure without the frozen protocol and
four method-status artifacts is not a completed benchmark.
An exact-budget failure returns `2` with status `preflight_failed`: it retains
and seals reconstruction and complete selection diagnostics, and records all
four methods as unrun. Such a run can be evaluated for its recorded failure,
but receives no historical-order accuracy.

After the sealed run completes:

```bash
python3 -m epistaeon.cli benchmark-evaluate \
  --run epistaeon/validation/runs/sws1-001 \
  --reference epistaeon/data/substitution_reference/sws1.json \
  --output epistaeon/validation/evaluations/sws1-001
```

The evaluator rejects changed, added, missing or symlinked run artifacts. It
reads the reference and prediction artifacts without loading a checkpoint or
rerunning inference. Its output is separate from the sealed run:

- `evaluation.json`: analysis protocol/ancestry, empirical constraints, coverage, per-method metrics and
  testability reasons.
- `comparison.csv` and `methods/<scorer>.json`: convenient method comparisons.
- `evaluation_manifest.json`: reference, run-seal, evaluator-code and numerical
  output hashes.

An invoking agent can then ask `epistaeon-evaluator` to interpret those outputs
and the original sources. That agent has read-only tools and returns a cited
Markdown assessment; the invoking orchestrator saves it as `assessment.md` in
the evaluation directory. The `/benchmark-epistaeon` command describes this
sequence. Do not use the existing manuscript-reading `epistaeon-replicator` for
this workflow.

A seal detects later changes; it is a checksum manifest, not a digital signature
or protection from a host administrator deliberately rewriting both files and
manifest. Preserve the original sealed run for audit. A changed experiment
requires a new run directory and a new frozen protocol.

## Fixed inference profile

New runs use protocol **2.0**, with HyPhy **2.5.96** SLAC joint maximum-likelihood
codons at the named internal node. `--asr-tree` and `--asr-node` are required.
`--asr-alignment` defaults to the scoring alignment; an augmented alignment can
add outgroups but must preserve every original scoring row and codon column
exactly. Its tree tips must match its taxa exactly. The selected node must be
internal, contain the focal taxon and retain the same descendants after
fitting. Missing, ambiguous, terminal or root-collapsed targets fail.

SLAC uses the Universal code, one CPU, seed `0`, all branches, zero reconstruction
samples, p-value threshold `0.1`, and preservation of zero-length internal
branches. It supplies a joint-ML point estimate, without residue posterior
probabilities; unresolved state sets remain in the reconstruction artifacts.
The inferred ancestor never becomes an additional scoring taxon. All original
scoring taxa are retained, independently of extra ASR outgroups.

The benchmark selects all resolved nonsynonymous SLAC-to-focal differences.
It excludes synonymous, indel and
unresolved differences, without using empirical sites or order. Excluded
columns retain the focal taxon's codons, so applying all selected substitutions
reproduces the complete focal nucleotide sequence exactly. For SWS1 these are
human codons. This benchmark policy does not change the existing `order`
command's default treatment of synonymous changes and endpoint backgrounds.

Settings are CPU, one numerical thread, seed `0`, all alignment taxa,
temperature `0.5`, alpha `1`, `[ROOT]`-token anchoring, no viability floor and
exact inference. The default maximum is 4,096 genotypes (`2^n` for `n` units).
There is also an unchanged hard limit of 22 units. The SWS1 inputs above select
**28 units**, requiring **268,435,456 genotypes**. This example therefore saves
and seals the reconstruction and selection, then stops before loading the
checkpoint. No automatic sampling, empirical filtering or site reduction occurs.

Native `order` also requires `--asr-tree` and `--asr-node` when `--ancestor` is
absent, and accepts `--asr-alignment` and an executable override `--hyphy`.
Otherwise it discovers the existing HyPhy 2.5.96 on PATH. Missing or incompatible
HyPhy stops execution. An explicit `--ancestor` alignment row bypasses SLAC and
conflicts with ASR arguments. Native substitution selection, grouping,
numbering and inference-mode policies remain unchanged. Its implicit anchor
remains `[ROOT]`; explicit ancestor rows retain their existing anchor default.

| Method | Interpretation |
| --- | --- |
| `potts-metropolis` | Existing Metropolis energy scorer; inspect coupling counts and uniformity. |
| `typicality-gate` | Existing extant-cloud typicality scorer with no hard viability floor. |
| `softmax-repaired` | Existing CLI semantics, including remaining-unit normalization and appended-row/projected-MDS placement. Its name does not imply those known limitations were repaired here. |
| `dphi-product` | Diagnostic baseline; its temperature does not discriminate order and nonpositive steps are blocked. |

The other three methods replace the focal row for each genotype. The benchmark
does not tune, repair or align these methods' semantics to the empirical answer.
Uniform predictions carry no identified MAP-order evidence, even if numerical
tie-breaking happens to display a historically compatible order. Blocked paths
and execution failures remain explicit diagnostics.

## Isolation and runtime provenance

`benchmark_runtime/runtime.json` records the platform, immutable official Python
base-image digest, resolver profile, CPU PyTorch wheel and TN93 implementation.
`requirements.lock` pins every transitive dependency to a version and the hash
of its selected Linux/amd64 CPython 3.11 wheel. Docker installs with
`--require-hashes --only-binary=:all:`. HyPhy is compiled from the upstream
2.5.96 commit and SHA-256-pinned tarball recorded in `runtime.json`, using a
fixed Debian package snapshot. Its batch libraries are installed in the same
image. The build disables OpenMP, MPI, CURL, BLAS and architecture-specific
vectorization; no network is available during SLAC execution.

The image is built from a staged context containing only the Dockerfile,
requirements lock, runtime metadata and allowlisted production Python modules
under `src/`. The runtime contains those production modules, pinned dependencies
and HyPhy with its libraries.
Repository documentation, empirical references, original papers, tests,
experimental data and agent histories are excluded. Input files have generic
container paths. At inference time the container has no network, a read-only
root filesystem, dropped capabilities, no privilege escalation and a temporary
`/tmp`; it receives only read-only `/inputs` and writable `/outputs` mounts.
The host's identity is used for output ownership. Each container is also limited
to one CPU, 4 GiB memory and 256 processes.

TN93 is pinned to the Python `tn93.tn93.TN93` implementation from package
`tn93==1.2.2`. The image rejects an external `tn93` executable so the engine's
binary-preference logic cannot silently choose a different implementation.
PyTorch is the official `2.5.1+cpu` CPython 3.11 Linux wheel. The base image is
`python:3.11.11-slim-bookworm` pinned to its linux/amd64 manifest rather than a
moving tag. An ARM host needs Docker's amd64 emulation support. Missing Docker
or runtime prerequisites stop execution; there is no unisolated fallback.

The preparation protocol freezes alignment, ASR alignment/tree, checkpoint,
configuration and source hashes before reconstruction, and records the built
runtime image ID. The scoring protocol additionally freezes the reconstructed
sequence artifact. `execution_manifest.json` also records the build
context digest, dependency-lock hash and Docker version. Repeatability concerns numerical artifacts within the
specified tolerances; elapsed times, timestamps and logs are audit metadata and
can differ. Isolation prevents runtime access to empirical answers. It does
not establish that training data were free of those answers, and it does not
prevent an operator from designing a later truth-informed experiment.

To intentionally update dependencies, revise `requirements.in`, resolve for
Python 3.11/linux amd64 using the resolver profile in `runtime.json`, download
all target wheels with hash verification, and retain only each selected wheel's
SHA-256 in the final lock. Validate dependency consistency and the container
smoke test, then commit the lock and runtime metadata as a new protocol. Never
regenerate the lock during a benchmark run.

## What evaluation means

Events match by **alignment codon column and ancestral-to-derived amino-acid
direction**, using the empirical crosswalk only in the evaluator. Coverage is
reported independently for target events and ordered pairs. Unmatched predicted
events are not automatically false positives. Empirically unresolved pairs are
excluded from order scoring.

For covered ordered pairs, the evaluator reports mean probability in the
empirical direction, Brier score and improvement over the uniform prediction
(`p=0.5`, Brier `0.25`), and MAP inversions when the run identifies a MAP order.
Uniform results receive probability-baseline scores and no credit from arbitrary
MAP ties. The projected MAP is compared with the eight historical candidates
only when all seven empirical substitutions are covered and MAP comparison is
meaningful.

For the SWS1 SLAC inputs, the complete selected candidate set represents **7/7
target events and all 18 ordered relationships**. Because the full 28-unit run
exceeds the budget, no method executes and evaluation reports
`not_testable: analysis_failed`, with no historical-order score. Candidate-set
coverage can be verified after sealing; it is not an inferred order.

Archived protocol **1.0** runs retain their original plurality semantics and
artifacts. Their SWS1 fixtures select 10 units and cover **2/7 target events and
0/18 ordered pairs**, yielding `not_testable: insufficient_order_coverage`.
The evaluator accepts both protocol versions and never relabels legacy runs as
SLAC. These acceptance values live only in tests and evaluator documentation;
they must never guide inference selection.

Yokoyama's eight favored paths summarize phylogenetic/parsimony history, not the
only biochemically accessible paths. Full 127-mutant assay-accessibility scoring
is deferred: the normalized supplementary table currently enumerates 1,076
accessible paths, whereas the paper reports 1,032. Preserve both figures and
the discrepancy rather than selecting the interpretation that favors a method.
See the original study, especially Figures 4–6 and supplementary Table S1:
[Yokoyama et al., 2014, PLOS Genetics](https://doi.org/10.1371/journal.pgen.1004884).

## Verification

Run the standard suite from the repository root:

```bash
PYTHONPATH=epistaeon/src:aeon-core/src:hyphaeon/src \
  python3 -m pytest epistaeon/tests -q
```

The real-checkpoint Docker integration is opt-in and requires the local model
and a running Docker daemon:

```bash
EPISTAEON_DOCKER_TEST=1 PYTHONPATH=epistaeon/src:aeon-core/src:hyphaeon/src \
  python3 -m pytest epistaeon/tests/test_benchmark_docker.py -q -s
```

It runs all four methods twice on a synthetic alignment with two selected
changes, verifies isolation and seals, and compares numerical predictions with
absolute tolerance `1e-7` and relative tolerance `1e-6`. It creates empirical
truth only after both analyses are sealed and verifies that evaluating changed
truth leaves every analysis artifact unchanged. This integration test does not
claim a full SWS1 inference run. A second Docker test reconstructs SWS1, verifies
the sealed 28-unit budget failure, and checks all seven empirical events only
after sealing. Native adapter tests also reconstruct the supplied SWS1 inputs.
Legacy fixtures continue to check protocol 1.0 coverage independently.
