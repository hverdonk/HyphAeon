---
description: Run the full epistaeon validation harness — replicator and adversary per study, then the summarizer.
---

Run the validation harness in `epistaeon/validation/`. Read its README first.

The original published studies are the only source of truth. The white paper in `epistaeon/` is not an input and must not be read, cited, or tested against by any agent.

For each study in `epistaeon/validation/studies/` (hemoglobin, myoglobin, steroid_receptor, rhodopsin):

1. Spawn `epistaeon-replicator` with the study_id. Wait for `reports/<id>.replication.json`.
2. Spawn `epistaeon-adversary` with the same study_id, only after the replication report exists.

The four studies are independent, so run them in parallel; the two agents within a study are not — the adversary needs the replicator's report.

Then:

```bash
python3 epistaeon/validation/compact_reports.py --compact
```

Fix any validation failure by sending the agent back to correct its own report. Do not edit an agent's report yourself.

Finally spawn `epistaeon-summarizer` to write `reports/SUMMARY.md`, and relay its verdict table plus the overall assessment.

If `epistaeon` is not yet implemented, stop and say so rather than running the agents — they have nothing to test.
