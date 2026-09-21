# Validation harness for `epistaeon`

Tests whether `epistaeon` can reproduce the findings of four published studies.
Built before `epistaeon` exists, so the success criteria are fixed in advance
rather than written around whatever the implementation happens to produce.

**The original studies are the only source of truth.** The white paper is not
an input: agents are instructed not to read it, the registries contain no claim
drawn from it, and the validator rejects any report citing it as a source.
Its assertions are unreliable, so agreement or disagreement with it would carry
no information.

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
  cannot be tested against the available inputs: the deer mouse is missing from
  the globin alignments, the Andean waterfowl is absent from TOGA2 entirely,
  and no complete 2^K phenotype panel was ever published for any of these
  systems. Scoring those as `epistaeon` failures would be
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
| `DATA_FACTS.md` | Operational facts agents need — numbering offsets, contact baselines, species coverage, model constraints. Facts about the data, not findings to replicate |
| `schemas/*.json` | Report schemas for the two hand-offs |
| `reports/` | Agent output (gitignored except this structure) |
| `compact_reports.py` | Schema validation plus deterministic compaction |

## Study status before any run

From the structural and alignment audits already in `data/README.md` and
`data/alignments/toga2/README.md`:

| Study | Targets | Prospects |
| --- | --- |
| `steroid_receptor` | SR1–SR3 | **Best case.** Real ancestral sequences, a real 16-genotype panel, full species coverage. |
| `rhodopsin` | RH1–RH2 | **Best epistasis ground truth**: a published, statistically tested species-by-site interaction at site 83. |
| `myoglobin` | MB1–MB2 | Species coverage good; the focal taxon is also the structure's species. |
| `hemoglobin` | HB1–HB3 | **Largely untestable.** The focal species is absent from the alignments and the waterfowl is absent from TOGA2. |

Ten replication targets across the four studies. A `not_testable` verdict on
`hemoglobin` is the expected honest outcome; the harness is built to make that
finding legible rather than to avoid it.
