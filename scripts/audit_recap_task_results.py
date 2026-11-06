#!/usr/bin/env python3
"""Summarize RECAP-Iso runtime coverage and safety invariants for tasks."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any


def _connectivity(value: Any) -> str:
    return str(value or "").strip().upper().split("-")[0]


def _task_report(path: Path) -> dict[str, Any]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    matched = 0
    evaluated = 0
    formula_conditioned = 0
    ranking_changed = 0
    counterfactual_ranking_changed = 0
    statuses: Counter[str] = Counter()
    reasons: Counter[str] = Counter()
    context_sources: Counter[str] = Counter()
    context_reasons: Counter[str] = Counter()
    changed_examples: list[dict[str, Any]] = []
    violations: list[dict[str, Any]] = []

    for row in rows:
        match = row.get("match")
        if not isinstance(match, dict):
            continue
        matched += 1
        evidence = match.get("isopair_evidence") or {}
        status = str(evidence.get("status") or "missing")
        reason = str(evidence.get("reason") or "")
        statuses[status] += 1
        if reason:
            reasons[reason] += 1
        context = evidence.get("adduct_context") or {}
        if context:
            context_sources[str(context.get("evidence_source") or "none")] += 1
            context_reasons[str(context.get("reason") or "missing")] += 1
        is_evaluated = bool(evidence.get("evaluated"))
        is_formula_conditioned = bool(context.get("formula_conditioned"))
        changed = bool(evidence.get("ranking_changed"))
        counterfactual_changed = bool(
            evidence.get("counterfactual_ranking_changed")
        )
        evaluated += int(is_evaluated)
        formula_conditioned += int(is_evaluated and is_formula_conditioned)
        ranking_changed += int(changed)
        counterfactual_ranking_changed += int(counterfactual_changed)

        selected = _connectivity(evidence.get("selected_connectivity"))
        actual = _connectivity(match.get("inchikey"))
        serialized_confidence = str(
            match.get("compound_confidence")
            or match.get("confidence")
            or ""
        )
        row_violations = []
        if is_formula_conditioned and evidence.get(
            "high_precision_gate_applicable"
        ):
            row_violations.append("formula_conditioned_high_gate_enabled")
        if is_formula_conditioned and evidence.get("high_precision_support"):
            row_violations.append("formula_conditioned_high_support")
        if changed and selected and actual != selected:
            row_violations.append("selected_connectivity_not_serialized_top1")
        if changed and serialized_confidence.lower() == "high":
            row_violations.append("reranked_compound_remained_high")
        if row_violations:
            violations.append({
                "spectrum_id": row.get("spectrum_id"),
                "precursor_mz": row.get("precursor_mz"),
                "violations": row_violations,
            })
        if (changed or counterfactual_changed) and len(changed_examples) < 25:
            changed_examples.append({
                "spectrum_id": row.get("spectrum_id"),
                "precursor_mz": row.get("precursor_mz"),
                "retention_time": row.get("retention_time"),
                "formula": match.get("formula"),
                "resolved_adduct": context.get("resolved_adduct"),
                "provisional_connectivity": evidence.get(
                    "provisional_winner_connectivity"
                ),
                "selected_connectivity": evidence.get(
                    "selected_connectivity"
                ),
                "counterfactual_selected_connectivity": evidence.get(
                    "counterfactual_selected_connectivity"
                ),
                "runtime_ranking_changed": changed,
                "counterfactual_ranking_changed": counterfactual_changed,
                "selected_name": match.get("compound_name"),
                "compound_confidence": serialized_confidence or None,
            })

    return {
        "result_file": str(path.resolve()),
        "spectra": len(rows),
        "matched_spectra": matched,
        "evaluated": evaluated,
        "evaluated_pct_of_matches": (
            round(100.0 * evaluated / matched, 3) if matched else 0.0
        ),
        "formula_conditioned_evaluated": formula_conditioned,
        "ranking_changed": ranking_changed,
        "counterfactual_ranking_changed": counterfactual_ranking_changed,
        "statuses": dict(statuses.most_common()),
        "not_applicable_reasons": dict(reasons.most_common()),
        "adduct_context_sources": dict(context_sources.most_common()),
        "adduct_context_reasons": dict(context_reasons.most_common()),
        "safety_violation_count": len(violations),
        "safety_violations": violations[:25],
        "ranking_change_examples": changed_examples,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tasks-dir", type=Path, default=Path("output/tasks"))
    parser.add_argument("--task-ids", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    tasks = {}
    for task_id in args.task_ids:
        path = args.tasks_dir / f"{task_id}_results.json"
        if not path.exists():
            raise FileNotFoundError(path)
        tasks[task_id] = _task_report(path)
    report = {
        "schema_version": 1,
        "tasks": tasks,
        "totals": {
            "spectra": sum(task["spectra"] for task in tasks.values()),
            "matched_spectra": sum(
                task["matched_spectra"] for task in tasks.values()
            ),
            "evaluated": sum(task["evaluated"] for task in tasks.values()),
            "formula_conditioned_evaluated": sum(
                task["formula_conditioned_evaluated"]
                for task in tasks.values()
            ),
            "ranking_changed": sum(
                task["ranking_changed"] for task in tasks.values()
            ),
            "counterfactual_ranking_changed": sum(
                task["counterfactual_ranking_changed"]
                for task in tasks.values()
            ),
            "safety_violation_count": sum(
                task["safety_violation_count"] for task in tasks.values()
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
