# Validation harness for `epistaeon`

Tests whether `epistaeon` can reproduce the findings of the four published
studies the white paper rests on. Built before `epistaeon` exists, so that the
success criteria are fixed in advance rather than written around whatever the
implementation happens to produce.

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
- **The study decides what counts as a replication target, not the white paper.**
  Each registry splits its findings into `replication_targets` (sourced to the
  published study) and `white_paper_claims_not_in_source` (assertions the white
  paper makes with no published basis in the local materials). Only the former
  affect any verdict. The latter are recorded for information — several are
  known to be false, such as `Thr36`, a residue absent from every receptor
  structure — and the validator rejects any report that scores `epistaeon`
  against them. Testing a model against a claim no study made is a category
  error, not a test. The adversary treats target drift as a critical challenge.
- **`not_testable` is a first-class verdict.** Much of this white paper cannot
  be tested at all: the deer mouse is missing from the globin alignments, the
  Andean waterfowl is absent from TOGA2 entirely, no complete 2^K phenotype
  panel exists, and the trajectory-order ground truth the paper assumes was
  never measured by anyone. Scoring those as `epistaeon` failures would be
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
python3 epistaeon/validation/compact_reports.py --validate reports/rhodopsin.replication.json
```

## Layout

| Path | Contents |
| --- | --- |
| `studies/*.json` | Per-study registry: published findings to test, local sources, structures, known data limits, minimum pair counts for the odds-ratio benchmark |
| `schemas/*.json` | Report schemas for the two hand-offs |
| `reports/` | Agent output (gitignored except this structure) |
| `compact_reports.py` | Schema validation plus deterministic compaction |

## Study status before any run

From the structural and alignment audits already in `data/README.md` and
`data/alignments/toga2/README.md`:

| Study | Prospects |
| --- | --- |
| `steroid_receptor` | **Best case.** Real ancestral sequences, a real 16-genotype panel, full species coverage. |
| `rhodopsin` | **Best epistasis ground truth**: a published, statistically tested species-by-site interaction at site 83. |
| `myoglobin` | Species coverage good; ordering claims have no experimental ground truth. |
| `hemoglobin` | **Largely untestable.** The focal species is absent from the alignments and the waterfowl is absent from TOGA2. |

A `not_testable` verdict on `hemoglobin` is the expected honest outcome. The
harness is built to make that finding legible rather than to avoid it.
