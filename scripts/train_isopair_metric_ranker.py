#!/usr/bin/env python3
"""Train the lightweight multi-metric IsoPair-MS candidate ranker.

Run ``train_isopair_ms.py`` first to create the formula-disjoint spectrum
selection cache.  This stage deliberately does not consume the DreaMS
embeddings: it learns from strict/broad product-ion, neutral-loss and matched
    intensity evidence that can be computed by the production V9 process itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import joblib
import lightgbm as lgb
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = Path(__file__).resolve().parent
for path in (str(PROJECT_ROOT), str(SCRIPTS_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from matching.isopair_features import FEATURE_NAMES, spectral_pair_features
from train_isopair_ms import (
    SpectrumRow,
    _cosine_score,
    _family_index,
    _jsonable_metrics,
    _prefer_cross_source_pair,
    _rank_metrics,
    _spectrum_signature,
    baseline_ambiguous_rows,
    baseline_high_metrics,
    enforce_formula_and_connectivity_disjoint_splits,
    high_similarity_isomer_rows,
    selective_metrics,
    split_overlap_audit,
    subset_report,
    tune_selective_high,
)


DEFAULT_CACHE = PROJECT_ROOT / "models" / "isopair_ms_v1_full"
DEFAULT_OUT = PROJECT_ROOT / "models" / "isopair_metric_v1"


def load_selection(
    cache_dir: Path,
    *,
    seed: int,
) -> tuple[list[SpectrumRow], list[str], dict[str, Any]]:
    payload = np.load(cache_dir / "pilot_embeddings.npz")
    splits = payload["splits"].astype(str).tolist()
    raw_rows = json.loads(
        (cache_dir / "pilot_spectra.json").read_text(encoding="utf-8")
    )
    spectra = [SpectrumRow(
        **{
            **row,
            "peaks": tuple(tuple(peak) for peak in row["peaks"]),
        }
    ) for row in raw_rows]
    if len(spectra) != len(splits):
        raise ValueError("Spectrum metadata and split cache have different lengths")
    splits = enforce_formula_and_connectivity_disjoint_splits(
        spectra, splits, seed=seed
    )
    audit = split_overlap_audit(spectra, splits)
    if not audit["zero_overlap"]:
        raise ValueError(f"IsoPair selection cache remains split-leaky: {audit}")
    return spectra, splits, audit


def _pair_feature(left: SpectrumRow, right: SpectrumRow) -> np.ndarray:
    return spectral_pair_features(
        left.peaks,
        right.peaks,
        query_precursor_mz=left.precursor_mz,
        reference_precursor_mz=right.precursor_mz,
        query_collision_energy=left.collision_energy,
        reference_collision_energy=right.collision_energy,
        query_collision_energy_unit=left.collision_energy_unit,
        reference_collision_energy_unit=right.collision_energy_unit,
        query_instrument_family=left.instrument_type,
        reference_instrument_family=right.instrument_type,
    )


def _parallel_features(
    spectra: Sequence[SpectrumRow],
    pairs: Sequence[tuple[int, int]],
    *,
    n_jobs: int,
) -> np.ndarray:
    started = time.time()
    rows = joblib.Parallel(n_jobs=n_jobs, prefer="threads", batch_size=128)(
        joblib.delayed(_pair_feature)(spectra[left], spectra[right])
        for left, right in pairs
    )
    print(
        f"features: {len(rows)} pairs in {time.time() - started:.1f}s",
        flush=True,
    )
    return np.asarray(rows, dtype=np.float32)


def make_pairs(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
    split: str,
    *,
    seed: int,
    positives_per_identity: int,
    hard_negatives: int,
) -> tuple[list[tuple[int, int]], np.ndarray, np.ndarray]:
    """Mine the most confusing same-formula identity by production cosine."""
    rng = random.Random(seed + {"train": 0, "validation": 1, "test": 2}[split])
    families = _family_index(spectra, splits, split)
    pairs: list[tuple[int, int]] = []
    labels: list[float] = []
    pair_groups: list[str] = []
    for (formula, adduct), identities in families.items():
        family_key = f"{formula}|{adduct}"
        items = [
            (identity, indices)
            for identity, indices in sorted(identities.items())
            if len(indices) >= 2
        ]
        if len(items) < 2:
            continue
        for identity, indices in items:
            for _ in range(positives_per_identity):
                anchor, positive = _prefer_cross_source_pair(indices, spectra, rng)
                anchor_signature = _spectrum_signature(spectra[anchor])
                # One representative per wrong connectivity prevents a heavily
                # replicated public source from monopolizing hard-negative mining.
                negative_pool = [
                    other_indices[0]
                    for other_identity, other_indices in items
                    if other_identity != identity
                    and _spectrum_signature(spectra[other_indices[0]])
                    != anchor_signature
                ]
                if not negative_pool:
                    continue
                similarities = [
                    _cosine_score(spectra[anchor], spectra[index])
                    for index in negative_pool
                ]
                order = np.argsort(-np.asarray(similarities), kind="stable")
                selected_negatives = [
                    negative_pool[int(position)]
                    for position in order[:hard_negatives]
                ]
                pairs.append((anchor, positive))
                labels.append(1.0)
                pair_groups.append(family_key)
                for negative in selected_negatives:
                    pairs.append((anchor, negative))
                    labels.append(0.0)
                    pair_groups.append(family_key)
    order = list(range(len(labels)))
    rng.shuffle(order)
    return (
        [pairs[index] for index in order],
        np.asarray([labels[index] for index in order], dtype=np.int8),
        np.asarray([pair_groups[index] for index in order], dtype=str),
    )


def build_retrieval_pairs(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
    split: str,
    *,
    max_queries_per_identity: int = 4,
) -> tuple[list[dict[str, Any]], list[tuple[int, int]], list[tuple[int, str]]]:
    """Flatten leave-one-spectrum-out retrieval comparisons for batch scoring."""
    families = _family_index(spectra, splits, split)
    rows: list[dict[str, Any]] = []
    pairs: list[tuple[int, int]] = []
    pair_owner: list[tuple[int, str]] = []
    for (formula, adduct), identities in sorted(families.items()):
        eligible = {
            identity: tuple(indices)
            for identity, indices in identities.items()
            if len(indices) >= 2
        }
        if len(eligible) < 2:
            continue
        for correct_identity, indices in sorted(eligible.items()):
            for query_index in indices[:max_queries_per_identity]:
                row_index = len(rows)
                rows.append({
                    "formula": formula,
                    "adduct": adduct,
                    "query_index": query_index,
                    "correct_connectivity": correct_identity,
                    "candidates": [],
                })
                for identity, references in sorted(eligible.items()):
                    available = [
                        index for index in references if index != query_index
                    ]
                    if not available:
                        continue
                    rows[row_index]["candidates"].append({
                        "connectivity": identity,
                        "baseline": max(
                            _cosine_score(spectra[query_index], spectra[index])
                            for index in available
                        ),
                        "metric_model": 0.0,
                    })
                    for reference in available:
                        pairs.append((query_index, reference))
                        pair_owner.append((row_index, identity))
    return rows, pairs, pair_owner


def add_model_scores(
    rows: list[dict[str, Any]],
    pair_owner: Sequence[tuple[int, str]],
    probabilities: np.ndarray,
) -> None:
    maxima: dict[tuple[int, str], float] = {}
    for owner, probability in zip(pair_owner, probabilities):
        maxima[owner] = max(maxima.get(owner, 0.0), float(probability))
    for row_index, row in enumerate(rows):
        for candidate in row["candidates"]:
            candidate["metric_model"] = maxima.get(
                (row_index, candidate["connectivity"]), 0.0
            )


def tune_fusion(rows: Sequence[dict[str, Any]]) -> float:
    best_key = (-1.0, -1.0, -1.0)
    best_alpha = 1.0
    for alpha in np.linspace(0.0, 1.0, 41):
        for row in rows:
            for candidate in row["candidates"]:
                candidate["fused"] = (
                    float(alpha) * float(candidate["baseline"])
                    + (1.0 - float(alpha)) * float(candidate["metric_model"])
                )
        metrics = _rank_metrics(rows, "fused")
        ambiguous = baseline_ambiguous_rows(rows)
        ambiguous_metrics = _rank_metrics(ambiguous, "fused")
        key = (
            float(metrics["top1_accuracy"]),
            float(ambiguous_metrics["top1_accuracy"]),
            float(metrics["mean_reciprocal_rank"]),
        )
        if key > best_key:
            best_key = key
            best_alpha = float(alpha)
    return best_alpha


def apply_fusion(rows: Sequence[dict[str, Any]], alpha: float) -> None:
    for row in rows:
        for candidate in row["candidates"]:
            candidate["fused"] = (
                alpha * float(candidate["baseline"])
                + (1.0 - alpha) * float(candidate["metric_model"])
            )
            # Compatibility keys for shared subset-report helpers.
            candidate["isopair"] = float(candidate["metric_model"])
            candidate["frozen_ssl"] = float(candidate["metric_model"])


def metric_subset_report(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    return {
        key: _jsonable_metrics(_rank_metrics(rows, key))
        for key in ("baseline", "metric_model", "fused")
    }


def _stable_fraction(text: str, seed: int) -> float:
    digest = hashlib.sha256(f"{seed}|{text}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") / float(2**64)


def split_validation_rows(
    rows: Sequence[dict[str, Any]],
    *,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Keep fusion tuning and selective-threshold calibration family-disjoint."""
    formulas = sorted({str(row.get("formula") or "") for row in rows})
    fusion_formulas = {
        formula for formula in formulas
        if _stable_fraction(f"validation-role|{formula}", seed) < 0.5
    }
    if not fusion_formulas or fusion_formulas == set(formulas):
        fusion_formulas = set(formulas[::2])
    fusion_rows = [
        row for row in rows if str(row.get("formula") or "") in fusion_formulas
    ]
    threshold_rows = [
        row for row in rows if str(row.get("formula") or "") not in fusion_formulas
    ]
    if not fusion_rows or not threshold_rows:
        raise ValueError(
            "Need at least two validation formula families to separate fusion "
            "and selective-threshold tuning"
        )
    return fusion_rows, threshold_rows


def group_early_stopping_indices(
    groups: np.ndarray,
    labels: np.ndarray,
    *,
    seed: int,
    validation_fraction: float = 0.15,
) -> tuple[np.ndarray, np.ndarray]:
    """Deterministically reserve complete train families for early stopping."""
    unique_groups = sorted(set(groups.astype(str).tolist()))
    validation_groups = {
        group for group in unique_groups
        if _stable_fraction(f"early-stop|{group}", seed) < validation_fraction
    }
    if not validation_groups or validation_groups == set(unique_groups):
        validation_groups = set(unique_groups[::max(2, int(1 / validation_fraction))])
    validation_mask = np.asarray([
        str(group) in validation_groups for group in groups
    ])
    train_index = np.where(~validation_mask)[0]
    validation_index = np.where(validation_mask)[0]
    if (
        not len(train_index)
        or not len(validation_index)
        or len(np.unique(labels[train_index])) < 2
        or len(np.unique(labels[validation_index])) < 2
    ):
        raise ValueError("Family-level early-stop split lacks both pair labels")
    return train_index, validation_index


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--seed", type=int, default=3407)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--positives-per-identity", type=int, default=2)
    parser.add_argument("--hard-negatives", type=int, default=3)
    parser.add_argument("--target-high-precision", type=float, default=0.95)
    parser.add_argument("--min-high-calls", type=int, default=50)
    parser.add_argument("--reuse-pair-cache", action="store_true")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    spectra, splits, selection_split_audit = load_selection(
        args.cache_dir, seed=args.seed
    )
    print(f"loaded {len(spectra)} formula-disjoint spectra")

    feature_cache = args.out_dir / "train_pair_features.npz"
    if args.reuse_pair_cache and feature_cache.exists():
        cache = np.load(feature_cache)
        if "groups" not in cache.files:
            raise ValueError(
                "Legacy pair cache lacks formula-family groups; rebuild without "
                "--reuse-pair-cache"
            )
        train_x = cache["x"].astype(np.float32)
        train_y = cache["y"].astype(np.int8)
        train_groups = cache["groups"].astype(str)
        print(f"reused {len(train_y)} training-pair features")
    else:
        train_pairs, train_y, train_groups = make_pairs(
            spectra,
            splits,
            "train",
            seed=args.seed,
            positives_per_identity=args.positives_per_identity,
            hard_negatives=args.hard_negatives,
        )
        print(
            f"training pairs={len(train_y)} positives={int(train_y.sum())} "
            f"negatives={int((train_y == 0).sum())}"
        )
        train_x = _parallel_features(
            spectra, train_pairs, n_jobs=args.n_jobs
        )
        np.savez_compressed(
            feature_cache,
            schema_version=np.asarray([2]),
            x=train_x,
            y=train_y,
            groups=train_groups,
        )

    model_parameters = dict(
        objective="binary",
        n_estimators=2000,
        learning_rate=0.025,
        num_leaves=31,
        max_depth=-1,
        min_child_samples=100,
        subsample=0.85,
        colsample_bytree=0.80,
        reg_alpha=0.20,
        reg_lambda=2.0,
        class_weight="balanced",
        random_state=args.seed,
        n_jobs=args.n_jobs,
        verbosity=-1,
    )
    early_train, early_validation = group_early_stopping_indices(
        train_groups, train_y, seed=args.seed
    )
    stopping_model = lgb.LGBMClassifier(**model_parameters)
    stopping_model.fit(
        train_x[early_train],
        train_y[early_train],
        feature_name=list(FEATURE_NAMES),
        eval_set=[(train_x[early_validation], train_y[early_validation])],
        eval_metric="binary_logloss",
        callbacks=[lgb.early_stopping(100, verbose=False)],
    )
    best_iteration = int(stopping_model.best_iteration_ or 1)
    refit_parameters = {**model_parameters, "n_estimators": best_iteration}
    model = lgb.LGBMClassifier(**refit_parameters)
    model.fit(train_x, train_y, feature_name=list(FEATURE_NAMES))
    train_probability = model.predict_proba(train_x)[:, 1]
    print(
        "pair train metrics:",
        f"AUROC={roc_auc_score(train_y, train_probability):.4f}",
        f"AUPRC={average_precision_score(train_y, train_probability):.4f}",
    )

    validation_rows, validation_pairs, validation_owner = build_retrieval_pairs(
        spectra, splits, "validation"
    )
    test_rows, test_pairs, test_owner = build_retrieval_pairs(
        spectra, splits, "test"
    )
    validation_x = _parallel_features(
        spectra, validation_pairs, n_jobs=args.n_jobs
    )
    test_x = _parallel_features(spectra, test_pairs, n_jobs=args.n_jobs)
    add_model_scores(
        validation_rows,
        validation_owner,
        model.predict_proba(validation_x)[:, 1],
    )
    add_model_scores(
        test_rows,
        test_owner,
        model.predict_proba(test_x)[:, 1],
    )

    fusion_tuning_rows, threshold_calibration_rows = split_validation_rows(
        validation_rows, seed=args.seed
    )
    alpha = tune_fusion(fusion_tuning_rows)
    apply_fusion(validation_rows, alpha)
    apply_fusion(test_rows, alpha)
    high_thresholds = tune_selective_high(
        threshold_calibration_rows,
        target_precision=args.target_high_precision,
        min_calls=args.min_high_calls,
    )
    validation_ambiguous = baseline_ambiguous_rows(validation_rows)
    test_ambiguous = baseline_ambiguous_rows(test_rows)
    validation_high_similarity = high_similarity_isomer_rows(validation_rows)
    test_high_similarity = high_similarity_isomer_rows(test_rows)

    report = {
        "schema_version": 2,
        "method": "IsoPair-MS strict/broad multi-metric gradient-boosted ranker",
        "created_at_unix": time.time(),
        "seed": args.seed,
        "split_policy": (
            "normalized-formula/connectivity-union-component disjoint; "
            "validation formula families separately assigned to fusion or "
            "selective-threshold tuning"
        ),
        "selection_split_audit": selection_split_audit,
        "identity_definition": "first 14 InChIKey characters (2-D connectivity)",
        "stereochemistry_claimed": False,
        "n_spectra": len(spectra),
        "n_training_pairs": len(train_y),
        "best_iteration_from_family_early_stopping": best_iteration,
        "validation_role_counts": {
            "fusion_tuning_queries": len(fusion_tuning_rows),
            "selective_threshold_queries": len(threshold_calibration_rows),
            "fusion_tuning_formulas": len({
                row["formula"] for row in fusion_tuning_rows
            }),
            "selective_threshold_formulas": len({
                row["formula"] for row in threshold_calibration_rows
            }),
        },
        "feature_names": list(FEATURE_NAMES),
        "fusion_alpha_baseline": alpha,
        "fusion_alpha_metric_model": 1.0 - alpha,
        "validation": {
            "all": metric_subset_report(validation_rows),
            "fusion_tuning_only": metric_subset_report(fusion_tuning_rows),
            "threshold_calibration_only": metric_subset_report(
                threshold_calibration_rows
            ),
            "baseline_ambiguous": metric_subset_report(validation_ambiguous),
            "high_similarity_isomer": metric_subset_report(
                validation_high_similarity
            ),
            "selective_high": selective_metrics(
                threshold_calibration_rows, high_thresholds
            ),
        },
        "test": {
            "all": metric_subset_report(test_rows),
            "baseline_ambiguous": metric_subset_report(test_ambiguous),
            "high_similarity_isomer": metric_subset_report(test_high_similarity),
            "baseline_high": baseline_high_metrics(test_rows),
            "selective_high": selective_metrics(test_rows, high_thresholds),
        },
        "selective_high_thresholds": high_thresholds,
        "limitations": [
            "The locked test is formula-disjoint but still derived from public-library records; an external same-method standard panel remains required for publication-grade external validation.",
            "The model targets constitutional/positional isomer connectivity and does not claim general stereoisomer resolution.",
            "Ranking influence must remain gated to true same-formula candidate groups and the locked applicability domain.",
            "Validation families used to tune fusion alpha are disjoint from those used to calibrate selective thresholds.",
            "Production ranking remains shadow-only until a new independent family cohort demonstrates a non-zero qualified gain.",
        ],
    }
    model.booster_.save_model(str(args.out_dir / "isopair_metric_v1.txt"))
    (args.out_dir / "locked_evaluation.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (args.out_dir / "runtime_config.json").write_text(
        json.dumps({
            "schema_version": 2,
            "model_name": "IsoPair-MS multi-metric ranker",
            "feature_names": list(FEATURE_NAMES),
            "fusion_alpha_baseline": alpha,
            "fusion_alpha_metric_model": 1.0 - alpha,
            "selective_high_thresholds": high_thresholds,
            "identity_definition": "InChIKey first block",
            "runtime_policy": {
                "ranking_influence_enabled": False,
                "ranking_deployment_state": "shadow_pending_independent_gain",
                "ranking_max_baseline_margin": 0.10,
                "ranking_scope": "same_formula_same_adduct_2d_connectivity",
                "high_precision_gate": "filter_existing_msms_high_only",
                "high_promotion_allowed": False,
                "reranked_near_tie_high_promotion_allowed": False,
                "stereochemistry_claimed": False,
            },
        }, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (args.out_dir / "validation_predictions.json").write_text(
        json.dumps(validation_rows, ensure_ascii=False),
        encoding="utf-8",
    )
    (args.out_dir / "test_predictions.json").write_text(
        json.dumps(test_rows, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
