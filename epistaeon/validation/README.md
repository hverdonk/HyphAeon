# Validation harness for `epistaeon`

Tests whether `epistaeon` can reproduce the findings of four empirical cases,
on two axes: **which epistatic sites it identifies**, and **the order of
mutations, where a study measured it**. `epistaeon`'s internal machinery is not
under test — only the sites and order it reports.
Built before `epistaeon` exists, so the success criteria are fixed in advance
rather than written around whatever the implementation happens to produce.

**The original studies are the only source of truth.** The white paper is not
an input: agents are instructed not to read it, the registries contain no claim
drawn from it, and the validator rejects any report citing it as a source.
Its assertions are unreliable, so agreement or disagreement with it would carry
no information.

The [source-verified reference](../data/substitution_reference/README.md)
(2026-10-15) extends and corrects the original registries. It includes event
identities, structural locations, qualified ordered sets and explicit unknowns.
Candidate-site recall and structural enrichment are not evidence of recovering
individually demonstrated epistasis. Respect each target's `benchmark_scope`
and `score_as_individual_site_epistasis` fields; the existing `validates` field
is retained for report-schema compatibility, not as a claim of ground truth.

## Design

Three agents, defined in `.claude/agents/`:

| Agent | Runs | Job |
| --- | --- | --- |
| `epistaeon-replicator` | once per study | Reads the original study, pre-registers criteria, runs `epistaeon`, scores the result. Skeptical but fair. |
| `epistaeon-adversary` | once per study | Reads the same study **independently first**, then tries to overturn the replicator and expose flaws. |
| `epistaeon-summarizer` | once, at the end | Adjudicates between them across all four studies and writes the final assessment. |

Between the per-study agents and the summarizer, `compact_reports.py` validates
every report against its schema and merges them. That step is deterministic on
purpose: no judgement is delegated to code, and no agent's prose reaches the
summarizer unchecked.

```
studies/<id>.json ─► replicator ─► <id>.replication.json ─┐
                                                          ├─► compacted.json ─► summarizer ─► SUMMARY.md
original study ───► adversary  ─► <id>.rebuttal.json ─────┘
```

## Why it is built this way

- **Pre-registration.** The replicator must write its pass and fail conditions
  before running anything. The validator rejects a report that says otherwise,
  and the adversary is told to check for criteria that moved after the fact.
- **Independent reading.** The adversary reads the source before the
  replicator's report, so a misreading cannot propagate. Reading disputes are
  carried all the way into the summary.
- **Sourcing is enforced, not requested.** Every scored claim must trace to a
  page, table or figure in one of the study's `primary_sources`. A claim citing
  anything else is rejected by the validator, and the adversary treats an
  unsourced replication target as a `critical` target-drift challenge. Dropping
  a genuine published finding because `epistaeon` did badly on it counts the
  same way.
- **`not_testable` is a first-class verdict.** Several published findings
  cannot be tested against the available inputs: the full phased experimental
  deer-mouse allele panel is missing, and some phenotype panels cover only a
  handful of genotypes. SWS1 instead has all 127 combinations of seven
  substitutions; its eight historically favored paths remain distinct from
  the paper's reported 1,032 paths accessible under the pigment-formation
  assay. Scoring missing
  inputs as `epistaeon` failures would be
  wrong, and scoring them as successes would be worse. They are counted
  separately and the summarizer may not average them away.
- **False negatives are hunted too.** The adversary must verify reported
  failures, because a bad numbering offset looks exactly like a null result.

## Running it

```bash
# 1. per study (4x), in either order
#    Agent tool: subagent_type "epistaeon-replicator", prompt naming the study_id
#    then       subagent_type "epistaeon-adversary",  same study_id

# 2. merge and validate
python3 epistaeon/validation/compact_reports.py --compact

# 3. final assessment
#    Agent tool: subagent_type "epistaeon-summarizer"
```

Or run the whole flow with `/validate-epistaeon`.

Validate a single report while iterating:

```bash
python3 epistaeon/validation/compact_reports.py --validate reports/sws1.replication.json
```

## Layout

| Path | Contents |
| --- | --- |
| `studies/*.json` | Per-study registry: published findings to test, local sources, structures, known data limits, minimum pair counts for the odds-ratio benchmark |
| `../data/substitution_reference/` | Verified events and locations, partial-order semantics, source hashes and missing inputs |
| `DATA_FACTS.md` | Operational facts agents need — numbering offsets, contact baselines, species coverage, model constraints. Facts about the data, not findings to replicate |
| `schemas/*.json` | Report schemas for the two hand-offs |
| `reports/` | Agent output (gitignored except this structure) |
| `compact_reports.py` | Schema validation plus deterministic compaction |

## Study status before any run

From the structural and alignment audits already in `data/README.md` and
`data/alignments/toga2/README.md` and the rebuilt SWS1 input in
`data/alignments/sws1/README.md`:

| Study | Targets | Prospects |
| --- | --- |
| `steroid_receptor` | SR1–SR6 | Assayed intermediate constraints, preferred paths and reverse group paths; 2011 C71S-first remains author conjecture. Reconstruction versions must stay distinct. |
| `sws1` | SWS1 spectral tuning and path accessibility | Yokoyama 2014 tests all 127 seven-site mutant combinations. The paper reports 60 nonfunctional mutants and 1,032 accessible orders out of 5,040. The four historical substitution sets yield eight favored paths by phylogenetic inference, not eight uniquely viable paths. S1 transcription/path-count discrepancies must be reconciled before exact accessibility scoring. |
| `myoglobin` | MB1–MB2 | Modern1MBO supports surface/charge enrichment. FigS2 branch assignments support phylogenetic partial order, not measured within-branch order or epistasis. |
| `hemoglobin` | HB1–HB3 | Representative deer-mouse sequences are now present. Storz measures isoform/mixed-hemolysate effects, not individual-site epistasis or mutation order. |

Do not score unsupported site-pair or temporal labels. The original studies
provide different kinds of evidence: experimentally assayed intermediate paths,
phylogenetic branch precedence, contextual phenotype effects and unknown order.
The reference preserves those distinctions. A comparator that joins model units
to source events and scores conditional sets has not yet been implemented.

The SWS1 input contains only extant boreoeutherian sequences rebuilt from the
paper's accessions. It is not the recovered original alignment and contains no
ancestral sequences. Source sites 46, 49, 52, 86, 93, 114 and 118 require the
explicit crosswalk in `data/numbering_offsets.json`; they are not raw alignment
columns or automatically current human UniProt positions. Bovine 1U19 is only
the homolog template used in the original paper's structural modeling.
