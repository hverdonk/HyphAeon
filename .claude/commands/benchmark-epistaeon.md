---
description: Run a blinded, sealed epistaeon scorer benchmark, evaluate it separately, then request read-only source-grounded interpretation.
---

Read `epistaeon/validation/BLINDED_BENCHMARK.md`. Use the requested alignment,
focal taxon, local checkpoint, ASR tree, named internal ASR node, optional
column-identical ASR outgroup alignment, the required manuscript substitution
panel, and a new run directory. Default to the documented
SWS1/Homo_sapiens example when the user has not supplied a different case.

1. Run the deterministic `epistaeon benchmark-run` command. Do not invoke an
   analysis agent or the manuscript-reading `epistaeon-replicator`. The runner
   accepts an alignment-bound column-only `--substitution-panel`, but no truth
   reference, blocks, supplied ancestor sequence or expected order. For SWS1 use
   `epistaeon/data/alignments/sws1/substitution_panel.json`. Freeze this scope
   before inference; do not change it after inspecting model results. SLAC
   reconstruction happens inside its isolated runtime. Keep empirical references,
   papers, repository documentation and agent histories outside its container.
2. Require a valid seal, frozen protocol and all four method-status artifacts
   before evaluation. `completed_with_failures` and a sealed `preflight_failed`
   exact-budget rejection are evaluable: preserve failed or unrun methods
   explicitly without assigning order accuracy. If setup fails or required artifacts are missing, report
   the recorded failure. Do not fall back to an unisolated run, sampling or fewer
   sites. A failure is not permission to change settings against the known answer.
3. Only after sealing, run `epistaeon benchmark-evaluate --run <sealed-run>
   --reference <empirical-reference> --output <separate-evaluation-directory>`.
   Numerical evaluation must not load the checkpoint or rerun inference. The
   evaluation directory must be outside the sealed run.
4. Spawn `epistaeon-evaluator` with those three paths. Give it the deterministic
   evaluation artifacts and primary sources. It has read-only tools and returns
   Markdown; save its returned text as `assessment.md` under the evaluation
   directory. Never copy the assessment or empirical reference into the run.
5. Return the four-method comparison, coverage and testability status, links to
   run/evaluation/assessment artifacts, and material limitations. Preserve the
   distinction between a missing test and an incorrect prediction.
   State that protocol 3.0 tests ordering within manuscript-supplied sites,
   not site discovery; excluded columns retain focal codons.

The evaluator's knowledge must not flow back into the completed inference run.
If the user later requests a changed protocol, create a new run with its own
frozen settings and identify it as a separate analysis. Never overwrite the
original sealed predictions or quietly tune against the study.
