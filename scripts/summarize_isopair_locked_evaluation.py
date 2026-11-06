#!/usr/bin/env python3
"""Summarize the locked IsoPair-MS runtime policy without model retuning.

The script consumes predictions already produced by the formula-disjoint test
and the later zero-overlap confirmation cohort.  It does not fit a model or
select a threshold.  Formula-family bootstrap sampling is used because spectra
within one molecular formula are not independent observations.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
from scipy.stats import binomtest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ORIGINAL = (
    PROJECT_ROOT / "models" / "isopair_metric_v2_noprec"
    / "test_predictions.json"
)
DEFAULT_CONFIRMATION = (
    PROJECT_ROOT / "models" / "isopair_confirmation_20260802"
    / "locked_predictions.json"
)
DEFAULT_CONFIG = (
    PROJECT_ROOT / "models" / "isopair_metric_v2_noprec"
    / "runtime_config.json"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "eval" / "v9val" / "isopair_locked_runtime.json"


def _rank(row: dict[str, Any], key: str) -> list[dict[str, Any]]:
    return sorted(
        row["candidates"],
        key=lambda candidate: (-float(candidate[key]), candidate["connectivity"]),
    )


def _facts(row: dict[str, Any], max_margin: float) -> dict[str, Any]:
    baseline = _rank(row, "baseline")
    fused = _rank(row, "fused")
    baseline_margin = float(baseline[0]["baseline"] - baseline[1]["baseline"])
    fused_margin = float(fused[0]["fused"] - fused[1]["fused"])
    runtime = fused[0] if baseline_margin <= max_margin else baseline[0]
    correct = row["correct_connectivity"]
    maximum_wrong = max(
        float(candidate["baseline"])
        for candidate in row["candidates"]
        if candidate["connectivity"] != correct
    )
    return {
        "formula": row["formula"],
        "baseline_winner": baseline[0]["connectivity"],
        "fused_winner": fused[0]["connectivity"],
        "runtime_winner": runtime["connectivity"],
        "correct_connectivity": correct,
        "baseline_score": float(baseline[0]["baseline"]),
        "baseline_margin": baseline_margin,
        "fused_score": float(fused[0]["fused"]),
        "fused_margin": fused_margin,
        "maximum_wrong_baseline": maximum_wrong,
    }


def _correct(row: dict[str, Any], winner_key: str) -> bool:
    return row[winner_key] == row["correct_connectivity"]


def _wilson_lower(correct: int, total: int, z: float = 1.96) -> float:
    if total <= 0:
        return 0.0
    proportion = correct / total
    denominator = 1.0 + z * z / total
    centre = proportion + z * z / (2.0 * total)
    radius = z * np.sqrt(
        proportion * (1.0 - proportion) / total
        + z * z / (4.0 * total * total)
    )
    return float((centre - radius) / denominator)


def _paired_statistics(
    rows: Sequence[dict[str, Any]],
    *,
    winner_key: str,
    seed: int,
    iterations: int,
) -> dict[str, Any]:
    by_formula: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_formula[row["formula"]].append(row)
    formulas = sorted(by_formula)
    baseline = np.asarray(
        [_correct(row, "baseline_winner") for row in rows], dtype=bool
    )
    alternative = np.asarray(
        [_correct(row, winner_key) for row in rows], dtype=bool
    )
    rng = np.random.default_rng(seed)
    gains = np.empty(iterations, dtype=float)
    for index in range(iterations):
        chosen = rng.choice(formulas, len(formulas), replace=True)
        sample = [row for formula in chosen for row in by_formula[formula]]
        gains[index] = np.mean([
            _correct(row, winner_key) for row in sample
        ]) - np.mean([
            _correct(row, "baseline_winner") for row in sample
        ])
    lost = int(np.sum(baseline & ~alternative))
    gained = int(np.sum(~baseline & alternative))
    p_value = (
        float(binomtest(min(lost, gained), lost + gained, 0.5).pvalue)
        if lost + gained else 1.0
    )
    return {
        "queries": len(rows),
        "formula_families": len(formulas),
        "baseline_accuracy": float(np.mean(baseline)) if len(rows) else 0.0,
        "runtime_accuracy": float(np.mean(alternative)) if len(rows) else 0.0,
        "absolute_gain": float(np.mean(alternative) - np.mean(baseline))
        if len(rows) else 0.0,
        "family_cluster_bootstrap_95_ci": [
            float(value) for value in np.quantile(gains, [0.025, 0.975])
        ],
        "baseline_only_correct": lost,
        "runtime_only_correct": gained,
        "paired_exact_p": p_value,
    }


def _subset(
    rows: Sequence[dict[str, Any]],
    predicate: Callable[[dict[str, Any]], bool],
) -> list[dict[str, Any]]:
    return [row for row in rows if predicate(row)]


def _high_metrics(
    rows: Sequence[dict[str, Any]],
    *,
    score_threshold: float,
    margin_threshold: float,
) -> dict[str, Any]:
    baseline_calls = [
        row for row in rows
        if row["baseline_score"] >= 0.90 and row["baseline_margin"] > 0.05
    ]
    gated_calls = [
        row for row in baseline_calls
        if row["runtime_winner"] == row["baseline_winner"]
        and row["fused_winner"] == row["runtime_winner"]
        and row["fused_score"] >= score_threshold
        and row["fused_margin"] >= margin_threshold
    ]

    def metrics(selected: Sequence[dict[str, Any]], winner_key: str) -> dict[str, Any]:
        correct = sum(_correct(row, winner_key) for row in selected)
        return {
            "calls": len(selected),
            "correct": int(correct),
            "precision": correct / len(selected) if selected else 0.0,
            "wilson_lower_95": _wilson_lower(correct, len(selected)),
            "coverage": len(selected) / len(rows) if rows else 0.0,
        }

    baseline_report = metrics(baseline_calls, "baseline_winner")
    gated_report = metrics(gated_calls, "runtime_winner")
    return {
        "baseline_approximate_high": baseline_report,
        "isopair_precision_gated_existing_high": gated_report,
        "absolute_precision_gain": (
            gated_report["precision"] - baseline_report["precision"]
        ),
        "existing_high_call_retention": (
            gated_report["calls"] / baseline_report["calls"]
            if baseline_report["calls"] else 0.0
        ),
        "policy": (
            "IsoPair selective thresholds filter an otherwise baseline-High "
            "unchanged winner; they never promote a reranked near-tie"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, default=DEFAULT_ORIGINAL)
    parser.add_argument("--confirmation", type=Path, default=DEFAULT_CONFIRMATION)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-iterations", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=99173)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    policy = config["runtime_policy"]
    max_margin = float(policy["ranking_max_baseline_margin"])
    thresholds = config["selective_high_thresholds"]
    score_threshold = float(thresholds["score_threshold"])
    margin_threshold = float(thresholds["margin_threshold"])

    cohorts: dict[str, list[dict[str, Any]]] = {}
    for name, path in (
        ("formula_disjoint_test", args.original),
        ("zero_overlap_confirmation", args.confirmation),
    ):
        raw = json.loads(path.read_text(encoding="utf-8"))
        cohorts[name] = [_facts(row, max_margin) for row in raw]
    combined = [row for rows in cohorts.values() for row in rows]

    report: dict[str, Any] = {
        "schema_version": 1,
        "status": "locked_runtime_policy_no_retuning",
        "identity_scope": "InChIKey first block (2-D connectivity)",
        "stereochemistry_claimed": False,
        "ranking_max_baseline_margin": max_margin,
        "selective_high_thresholds": thresholds,
        "cohorts": {},
    }
    for cohort_index, (name, rows) in enumerate(cohorts.items()):
        ambiguous = _subset(rows, lambda row: row["baseline_margin"] <= 0.05)
        high_similarity = _subset(
            rows, lambda row: row["maximum_wrong_baseline"] >= 0.70
        )
        report["cohorts"][name] = {
            "all": _paired_statistics(
                rows,
                winner_key="runtime_winner",
                seed=args.seed + cohort_index,
                iterations=args.bootstrap_iterations,
            ),
            "baseline_ambiguous": _paired_statistics(
                ambiguous,
                winner_key="runtime_winner",
                seed=args.seed + 10 + cohort_index,
                iterations=args.bootstrap_iterations,
            ),
            "high_similarity_wrong_isomer": _paired_statistics(
                high_similarity,
                winner_key="runtime_winner",
                seed=args.seed + 20 + cohort_index,
                iterations=args.bootstrap_iterations,
            ),
            "high_precision_gate": _high_metrics(
                rows,
                score_threshold=score_threshold,
                margin_threshold=margin_threshold,
            ),
        }
    combined_ambiguous = _subset(
        combined, lambda row: row["baseline_margin"] <= 0.05
    )
    combined_high_similarity = _subset(
        combined, lambda row: row["maximum_wrong_baseline"] >= 0.70
    )
    report["combined_independent_tests"] = {
        "all": _paired_statistics(
            combined,
            winner_key="runtime_winner",
            seed=args.seed + 100,
            iterations=args.bootstrap_iterations,
        ),
        "baseline_ambiguous": _paired_statistics(
            combined_ambiguous,
            winner_key="runtime_winner",
            seed=args.seed + 110,
            iterations=args.bootstrap_iterations,
        ),
        "high_similarity_wrong_isomer": _paired_statistics(
            combined_high_similarity,
            winner_key="runtime_winner",
            seed=args.seed + 120,
            iterations=args.bootstrap_iterations,
        ),
        "high_precision_gate": _high_metrics(
            combined,
            score_threshold=score_threshold,
            margin_threshold=margin_threshold,
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
