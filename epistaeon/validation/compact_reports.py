#!/usr/bin/env python3
"""Validate agent reports against their schemas and compact them for the summarizer.

Two jobs, both deterministic so that no judgement is delegated to this layer:

  --validate FILE   check one report against its schema (agents call this)
  --compact         merge every study's replication + rebuttal into compacted.json

Compaction drops prose the summarizer does not need and keeps every verdict,
claim outcome, surviving challenge and reading dispute.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCHEMAS = HERE / "schemas"
REPORTS = HERE / "reports"
STUDIES = HERE / "studies"


def _schema_for(report: dict) -> Path:
    agent = report.get("agent")
    if agent == "replicator":
        return SCHEMAS / "replication.schema.json"
    if agent == "adversary":
        return SCHEMAS / "rebuttal.schema.json"
    raise SystemExit(f"report has unknown or missing 'agent' field: {agent!r}")


def validate(path: Path) -> dict:
    report = json.loads(path.read_text())
    schema = json.loads(_schema_for(report).read_text())
    try:
        import jsonschema
    except ImportError:
        _minimal_check(report, schema, path)
        print(f"[ok] {path.name} passed minimal checks (install jsonschema for full validation)")
        return report
    import jsonschema
    try:
        jsonschema.validate(report, schema)
    except jsonschema.ValidationError as exc:
        loc = "/".join(str(p) for p in exc.absolute_path) or "(root)"
        raise SystemExit(f"[FAIL] {path.name} at {loc}: {exc.message}")
    _content_checks(report, path)
    print(f"[ok] {path.name} valid")
    return report


def _minimal_check(report, schema, path):
    missing = [k for k in schema.get("required", []) if k not in report]
    if missing:
        raise SystemExit(f"[FAIL] {path.name} missing required fields: {missing}")
    _content_checks(report, path)


def _content_checks(report, path):
    """Checks the schema cannot express, aimed at the known failure modes."""
    problems = []
    if report.get("agent") == "replicator":
        crit = report.get("preregistered_criteria", {})
        if not crit.get("written_before_run"):
            problems.append("preregistered_criteria.written_before_run is false; criteria must precede the run")
        if not report.get("limitations"):
            problems.append("limitations is empty; state the weaknesses of your own analysis")
        for claim in report.get("claims", []):
            if claim.get("outcome") == "not_testable" and not claim.get("not_testable_reason"):
                problems.append(f"claim {claim.get('claim_id')} is not_testable without a reason")
            if not claim.get("evidence"):
                problems.append(f"claim {claim.get('claim_id')} has no evidence entries")
            if claim.get("claim_provenance") == "white_paper_only" and claim.get("scored_against_epistaeon"):
                problems.append(
                    f"claim {claim.get('claim_id')} is white_paper_only but marked as scored; "
                    "epistaeon may not be scored against claims no study made")
    if report.get("agent") == "adversary":
        for ch in report.get("challenges", []):
            if not ch.get("evidence"):
                problems.append(f"challenge {ch.get('challenge_id')} has no evidence entries")
        if report.get("verdict_on_replicator") != "overturned" and not report.get("credit_where_due"):
            problems.append("credit_where_due is empty; name what was sound or state that nothing was")
    if problems:
        raise SystemExit(f"[FAIL] {path.name}:\n  - " + "\n  - ".join(problems))


def compact() -> dict:
    studies = sorted(p.stem for p in STUDIES.glob("*.json"))
    out = {"studies": {}, "missing_reports": [], "totals": {}}
    counts = {"replicated": 0, "partially_replicated": 0, "not_replicated": 0, "not_testable": 0}
    testable = untestable = 0
    unsourced_total = 0

    for sid in studies:
        rep_p = REPORTS / f"{sid}.replication.json"
        reb_p = REPORTS / f"{sid}.rebuttal.json"
        if not rep_p.exists() or not reb_p.exists():
            out["missing_reports"].append(
                {"study_id": sid, "replication": rep_p.exists(), "rebuttal": reb_p.exists()}
            )
            continue
        rep, reb = validate(rep_p), validate(reb_p)

        claims, unsourced = [], []
        for c in rep.get("claims", []):
            prov = c.get("claim_provenance", "published_study")
            entry = {
                "claim_id": c["claim_id"],
                "published_finding": c["published_finding"],
                "claim_provenance": prov,
                "outcome": c["outcome"],
                "not_testable_reason": c.get("not_testable_reason"),
                "epistaeon_result": c["epistaeon_result"],
                "numbers": c.get("numbers", {}),
            }
            if prov == "white_paper_only":
                unsourced.append(entry)       # recorded, never scored
                continue
            claims.append(entry)
            if c["outcome"] == "not_testable":
                untestable += 1
            else:
                testable += 1

        surviving = [
            {k: ch[k] for k in ("challenge_id", "target_claim_id", "type", "severity", "argument", "survives")}
            for ch in reb.get("challenges", [])
            if ch.get("survives") in ("no", "weakened", "undetermined")
        ]
        counts[reb["revised_verdict"]] = counts.get(reb["revised_verdict"], 0) + 1

        unsourced_total += len(unsourced)
        out["studies"][sid] = {
            "replicator_verdict": rep["verdict"],
            "replicator_confidence": rep["confidence"],
            "adversary_verdict": reb["revised_verdict"],
            "adversary_confidence": reb["confidence"],
            "verdict_on_replicator": reb["verdict_on_replicator"],
            "verdicts_agree": rep["verdict"] == reb["revised_verdict"],
            "reading_disputes": reb["independent_reading"].get("disagreements_with_replicator_reading", []),
            "claims": claims,
            "white_paper_claims_recorded": unsourced,
            "surviving_challenges": surviving,
            "critical_unsurvived": [c for c in surviving if c["severity"] == "critical" and c["survives"] == "no"],
            "replicator_limitations": rep.get("limitations", []),
            "credit_where_due": reb.get("credit_where_due", []),
        }

    out["totals"] = {
        "studies_total": len(studies),
        "studies_reported": len(out["studies"]),
        "verdicts": counts,
        "claims_testable": testable,
        "claims_not_testable": untestable,
        "white_paper_claims_recorded_not_scored": unsourced_total,
    }
    REPORTS.mkdir(exist_ok=True)
    dest = REPORTS / "compacted.json"
    dest.write_text(json.dumps(out, indent=2))
    print(f"[ok] compacted {len(out['studies'])}/{len(studies)} studies -> {dest}")
    if out["missing_reports"]:
        print("[!] missing reports:", ", ".join(m["study_id"] for m in out["missing_reports"]))
    print(f"    testable claims: {testable} | not testable: {untestable} | "
          f"white-paper claims recorded but not scored: {unsourced_total}")
    print(f"    verdicts: {counts}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--validate", metavar="FILE", help="validate one report against its schema")
    g.add_argument("--compact", action="store_true", help="merge all reports into compacted.json")
    args = ap.parse_args()
    if args.validate:
        validate(Path(args.validate))
    else:
        compact()


if __name__ == "__main__":
    main()
