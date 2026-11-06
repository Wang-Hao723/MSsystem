#!/usr/bin/env python3
"""Audit formula-conditioned RECAP-Iso coverage on locked cohorts.

This is a no-retuning audit.  It removes the explicit query adduct from each
previously locked same-formula retrieval case, reconstructs only the local
formula/mass-compatible adduct context, and reports the locked ranker's
counterfactual behavior.  MIX33 replay rejected production rank influence, so
the deployed policy is shadow-only and High precision gating is disabled for
every formula-conditioned case.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.recap_adduct import (  # noqa: E402
    canonical_adduct_key,
    resolve_isopair_adduct_context,
)
from utils.adduct import adduct_polarity  # noqa: E402


DEFAULT_OUTPUT = (
    PROJECT_ROOT / "eval" / "v9val" / "recap_formula_adduct_audit.json"
)
COHORTS = {
    "formula_disjoint_test": (
        PROJECT_ROOT / "models" / "isopair_metric_v2_noprec"
        / "test_predictions.json",
        PROJECT_ROOT / "models" / "isopair_ms_v1_full"
        / "pilot_spectra.json",
    ),
    "zero_overlap_confirmation": (
        PROJECT_ROOT / "models" / "isopair_confirmation_20260802"
        / "locked_predictions.json",
        PROJECT_ROOT / "models" / "isopair_confirmation_20260802"
        / "pilot_spectra.json",
    ),
}


def _rank(row: dict[str, Any], key: str) -> list[dict[str, Any]]:
    return sorted(
        row["candidates"],
        key=lambda item: (-float(item[key]), str(item["connectivity"])),
    )


def _ranking_facts(row: dict[str, Any], max_margin: float) -> dict[str, Any]:
    baseline = _rank(row, "baseline")
    fused = _rank(row, "fused")
    baseline_margin = float(
        baseline[0]["baseline"] - baseline[1]["baseline"]
    )
    runtime = fused[0] if baseline_margin <= max_margin else baseline[0]
    correct = str(row["correct_connectivity"])
    return {
        "baseline_correct": str(baseline[0]["connectivity"]) == correct,
        "runtime_correct": str(runtime["connectivity"]) == correct,
        "ranking_changed": (
            str(runtime["connectivity"]) != str(baseline[0]["connectivity"])
        ),
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    count = len(rows)
    baseline_correct = sum(bool(row["baseline_correct"]) for row in rows)
    runtime_correct = sum(bool(row["runtime_correct"]) for row in rows)
    return {
        "queries": count,
        "baseline_correct": baseline_correct,
        "runtime_correct": runtime_correct,
        "baseline_accuracy": baseline_correct / count if count else None,
        "runtime_accuracy": runtime_correct / count if count else None,
        "absolute_gain": (
            (runtime_correct - baseline_correct) / count if count else None
        ),
        "baseline_only_correct": sum(
            row["baseline_correct"] and not row["runtime_correct"]
            for row in rows
        ),
        "runtime_only_correct": sum(
            row["runtime_correct"] and not row["baseline_correct"]
            for row in rows
        ),
        "ranking_changes": sum(bool(row["ranking_changed"]) for row in rows),
    }


def _sum_metrics(reports: list[dict[str, Any]]) -> dict[str, Any]:
    queries = sum(int(report["queries"]) for report in reports)
    baseline_correct = sum(
        int(report["baseline_correct"]) for report in reports
    )
    runtime_correct = sum(
        int(report["runtime_correct"]) for report in reports
    )
    return {
        "queries": queries,
        "baseline_correct": baseline_correct,
        "runtime_correct": runtime_correct,
        "baseline_accuracy": baseline_correct / queries if queries else None,
        "runtime_accuracy": runtime_correct / queries if queries else None,
        "absolute_gain": (
            (runtime_correct - baseline_correct) / queries
            if queries else None
        ),
        "baseline_only_correct": sum(
            int(report["baseline_only_correct"]) for report in reports
        ),
        "runtime_only_correct": sum(
            int(report["runtime_only_correct"]) for report in reports
        ),
        "ranking_changes": sum(
            int(report["ranking_changes"]) for report in reports
        ),
    }


def _audit_cohort(
    prediction_path: Path,
    spectra_path: Path,
    *,
    max_margin: float,
) -> dict[str, Any]:
    predictions = json.loads(prediction_path.read_text(encoding="utf-8"))
    spectra = json.loads(spectra_path.read_text(encoding="utf-8"))
    reasons: Counter[str] = Counter()
    strata: Counter[str] = Counter()
    eligible_rows: list[dict[str, Any]] = []
    ineligible_rows: list[dict[str, Any]] = []
    neutral_invariance_failures = 0
    high_gate_violations = 0

    for row in predictions:
        query_index = int(row["query_index"])
        spectrum = spectra[query_index]
        stratum_adduct = str(row["adduct"])
        candidates = [
            {
                "inchikey": f"{candidate['connectivity']}-AAAAAAAAAA-A",
                "adduct": stratum_adduct,
            }
            for candidate in row["candidates"]
        ]
        kwargs = {
            "formula": row["formula"],
            "query_precursor_mz": spectrum["precursor_mz"],
            "query_adduct": "",
            "query_polarity": adduct_polarity(stratum_adduct),
        }
        context = resolve_isopair_adduct_context(
            **kwargs, candidates=candidates
        )
        context_without_candidates = resolve_isopair_adduct_context(
            **kwargs, candidates=[]
        )
        invariant_keys = (
            "status", "reason", "resolved_adduct", "adduct_key",
            "formula_conditioned", "ranking_eligible",
            "high_precision_gate_eligible",
        )
        if any(
            context.get(key) != context_without_candidates.get(key)
            for key in invariant_keys
        ):
            neutral_invariance_failures += 1
        if context.get("high_precision_gate_eligible"):
            high_gate_violations += 1

        facts = _ranking_facts(row, max_margin)
        resolved_matches_stratum = (
            context.get("adduct_key")
            == canonical_adduct_key(stratum_adduct)
        )
        eligible = bool(
            context.get("ranking_eligible")
            and resolved_matches_stratum
            and 2 <= len(row["candidates"]) <= 12
        )
        audit_row = {
            **facts,
            "formula": row["formula"],
            "stratum_adduct": stratum_adduct,
            "reason": context.get("reason"),
        }
        if eligible:
            eligible_rows.append(audit_row)
            strata[stratum_adduct] += 1
        else:
            if context.get("ranking_eligible") and not resolved_matches_stratum:
                reason = "resolved_adduct_does_not_match_locked_stratum"
            else:
                reason = str(context.get("reason") or "unknown")
            reasons[reason] += 1
            ineligible_rows.append(audit_row)

    all_rows = eligible_rows + ineligible_rows
    ineligible_baseline_rows = [
        {
            **row,
            "runtime_correct": row["baseline_correct"],
            "ranking_changed": False,
        }
        for row in ineligible_rows
    ]
    return {
        "queries": len(predictions),
        "formula_families": len({row["formula"] for row in predictions}),
        "formula_conditioned_ranking_eligible": len(eligible_rows),
        "formula_conditioned_ranking_coverage": (
            len(eligible_rows) / len(predictions) if predictions else 0.0
        ),
        "eligible_strata": dict(strata.most_common()),
        "ineligible_reasons": dict(reasons.most_common()),
        "candidate_label_neutrality_failures": neutral_invariance_failures,
        "formula_conditioned_high_gate_violations": high_gate_violations,
        "all_queries_if_explicit_adduct_available": _metrics(all_rows),
        "formula_conditioned_eligible_subset": _metrics(eligible_rows),
        "formula_conditioned_ineligible_subset_baseline_only": _metrics(
            ineligible_baseline_rows
        ),
        "formula_conditioned_rank_enabled_counterfactual_all_queries": _metrics(
            eligible_rows + ineligible_baseline_rows
        ),
        "deployed_shadow_policy_all_queries": _metrics(
            ineligible_baseline_rows + [
                {
                    **row,
                    "runtime_correct": row["baseline_correct"],
                    "ranking_changed": False,
                }
                for row in eligible_rows
            ]
        ),
    }


def _combined_report(cohorts: dict[str, dict[str, Any]]) -> dict[str, Any]:
    values = list(cohorts.values())
    eligible_strata: Counter[str] = Counter()
    ineligible_reasons: Counter[str] = Counter()
    for value in values:
        eligible_strata.update(value["eligible_strata"])
        ineligible_reasons.update(value["ineligible_reasons"])
    queries = sum(int(value["queries"]) for value in values)
    eligible = sum(
        int(value["formula_conditioned_ranking_eligible"])
        for value in values
    )
    return {
        "queries": queries,
        "formula_families": sum(
            int(value["formula_families"]) for value in values
        ),
        "formula_conditioned_ranking_eligible": eligible,
        "formula_conditioned_ranking_coverage": (
            eligible / queries if queries else 0.0
        ),
        "eligible_strata": dict(eligible_strata.most_common()),
        "ineligible_reasons": dict(ineligible_reasons.most_common()),
        "candidate_label_neutrality_failures": sum(
            int(value["candidate_label_neutrality_failures"])
            for value in values
        ),
        "formula_conditioned_high_gate_violations": sum(
            int(value["formula_conditioned_high_gate_violations"])
            for value in values
        ),
        "all_queries_if_explicit_adduct_available": _sum_metrics([
            value["all_queries_if_explicit_adduct_available"]
            for value in values
        ]),
        "formula_conditioned_eligible_subset": _sum_metrics([
            value["formula_conditioned_eligible_subset"]
            for value in values
        ]),
        "formula_conditioned_ineligible_subset_baseline_only": _sum_metrics([
            value["formula_conditioned_ineligible_subset_baseline_only"]
            for value in values
        ]),
        "formula_conditioned_rank_enabled_counterfactual_all_queries": _sum_metrics([
            value[
                "formula_conditioned_rank_enabled_counterfactual_all_queries"
            ]
            for value in values
        ]),
        "deployed_shadow_policy_all_queries": _sum_metrics([
            value["deployed_shadow_policy_all_queries"]
            for value in values
        ]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-baseline-margin", type=float, default=0.10)
    args = parser.parse_args()

    cohort_reports = {
        name: _audit_cohort(
            predictions,
            spectra,
            max_margin=args.max_baseline_margin,
        )
        for name, (predictions, spectra) in COHORTS.items()
    }
    report = {
        "schema_version": 1,
        "status": (
            "locked_no_retuning_formula_conditioned_adduct_shadow_audit"
        ),
        "identity_scope": "InChIKey first block (2-D connectivity)",
        "explicit_query_adduct_removed": True,
        "candidate_connectivity_neutrality_required": True,
        "formula_conditioned_high_precision_gate_allowed": False,
        "production_rank_influence_allowed": False,
        "production_policy": "shadow_only_after_mix33_validation_regressed",
        "ranking_max_baseline_margin": args.max_baseline_margin,
        "cohorts": cohort_reports,
        "combined_zero_overlap_cohorts": _combined_report(cohort_reports),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
