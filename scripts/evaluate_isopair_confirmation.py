#!/usr/bin/env python3
"""Evaluate a locked IsoPair-MS artifact on fresh, excluded formula families."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

import lightgbm as lgb
import numpy as np
from scipy.stats import binomtest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = Path(__file__).resolve().parent
for path in (str(PROJECT_ROOT), str(SCRIPTS_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from train_isopair_metric_ranker import (
    _jsonable_metrics,
    _parallel_features,
    _rank_metrics,
    add_model_scores,
    apply_fusion,
    baseline_ambiguous_rows,
    baseline_high_metrics,
    build_retrieval_pairs,
    high_similarity_isomer_rows,
    load_selection,
    metric_subset_report,
    selective_metrics,
)


def _correct(row: dict[str, Any], key: str) -> bool:
    winner = sorted(
        row["candidates"],
        key=lambda candidate: (-float(candidate[key]), candidate["connectivity"]),
    )[0]
    return winner["connectivity"] == row["correct_connectivity"]


def paired_family_statistics(
    rows: Sequence[dict[str, Any]],
    *,
    seed: int = 99173,
    iterations: int = 10000,
) -> dict[str, Any]:
    formulas = sorted({row["formula"] for row in rows})
    by_formula = {
        formula: [row for row in rows if row["formula"] == formula]
        for formula in formulas
    }
    baseline = np.asarray([_correct(row, "baseline") for row in rows], dtype=bool)
    fused = np.asarray([_correct(row, "fused") for row in rows], dtype=bool)
    rng = np.random.default_rng(seed)
    differences = np.empty(iterations, dtype=float)
    for index in range(iterations):
        selected = rng.choice(formulas, len(formulas), replace=True)
        sample = [row for formula in selected for row in by_formula[formula]]
        differences[index] = (
            np.mean([_correct(row, "fused") for row in sample])
            - np.mean([_correct(row, "baseline") for row in sample])
        )
    lost = int(np.sum(baseline & ~fused))
    gained = int(np.sum(~baseline & fused))
    p_value = (
        float(binomtest(min(lost, gained), lost + gained, 0.5).pvalue)
        if lost + gained else 1.0
    )
    return {
        "n_queries": len(rows),
        "n_formula_families": len(formulas),
        "baseline_accuracy": float(np.mean(baseline)) if len(rows) else 0.0,
        "fused_accuracy": float(np.mean(fused)) if len(rows) else 0.0,
        "absolute_gain": float(np.mean(fused) - np.mean(baseline))
        if len(rows) else 0.0,
        "family_cluster_bootstrap_95_ci": [
            float(value) for value in np.quantile(differences, [0.025, 0.975])
        ],
        "baseline_only_correct": lost,
        "fused_only_correct": gained,
        "paired_exact_p": p_value,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--confirmation-cache",
        type=Path,
        default=PROJECT_ROOT / "models" / "isopair_confirmation_20260802",
    )
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=PROJECT_ROOT / "models" / "isopair_metric_v2_noprec",
    )
    parser.add_argument("--n-jobs", type=int, default=8)
    args = parser.parse_args()

    spectra, splits = load_selection(args.confirmation_cache)
    rows, pairs, owners = build_retrieval_pairs(spectra, splits, "test")
    features = _parallel_features(spectra, pairs, n_jobs=args.n_jobs)
    model = lgb.Booster(
        model_file=str(args.artifact_dir / "isopair_metric_v1.txt")
    )
    add_model_scores(rows, owners, model.predict(features))
    config = json.loads(
        (args.artifact_dir / "runtime_config.json").read_text(encoding="utf-8")
    )
    alpha = float(config["fusion_alpha_baseline"])
    apply_fusion(rows, alpha)
    thresholds = config["selective_high_thresholds"]
    ambiguous = baseline_ambiguous_rows(rows)
    high_similarity = high_similarity_isomer_rows(rows)
    report = {
        "schema_version": 1,
        "status": "locked_confirmation_no_retuning",
        "model_artifact": str(args.artifact_dir.resolve()),
        "confirmation_manifest": json.loads(
            (args.confirmation_cache / "manifest.json").read_text(encoding="utf-8")
        ),
        "fusion_alpha_baseline_locked": alpha,
        "fusion_alpha_metric_model_locked": 1.0 - alpha,
        "selective_high_thresholds_locked": thresholds,
        "all": metric_subset_report(rows),
        "baseline_ambiguous": metric_subset_report(ambiguous),
        "high_similarity_isomer": metric_subset_report(high_similarity),
        "baseline_high": baseline_high_metrics(rows),
        "selective_high": selective_metrics(rows, thresholds),
        "paired_family_statistics": {
            "all": paired_family_statistics(rows),
            "baseline_ambiguous": paired_family_statistics(ambiguous),
            "high_similarity_isomer": paired_family_statistics(high_similarity),
        },
    }
    (args.confirmation_cache / "locked_confirmation_evaluation.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.confirmation_cache / "locked_predictions.json").write_text(
        json.dumps(rows, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
