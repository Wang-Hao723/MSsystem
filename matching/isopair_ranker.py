"""Runtime adapter for the formula-conditional IsoPair-MS ranker.

The model is evaluated inside the current production winner's same-formula,
same-adduct 2-D-connectivity group.  Ranking remains shadow-only because the
independent family confirmation interval crossed zero; the separately
qualified selective score may only veto an already-High call.  The module never
introduces a candidate, changes a production identity, or claims stereochemistry.
"""

from __future__ import annotations

import json
import math
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from config.settings import (
    ISOPAIR_RANK_INFLUENCE_ENABLED,
    ISOPAIR_ALLOWED_SOURCE_TIERS,
    ISOPAIR_CONFIG_PATH,
    ISOPAIR_MAX_BASELINE_MARGIN,
    ISOPAIR_MAX_CONNECTIVITIES,
    ISOPAIR_MIN_CONNECTIVITIES,
    ISOPAIR_MODEL_PATH,
)
from matching.honest_calibration import candidate_comparison_view
from matching.isopair_features import FEATURE_NAMES, spectral_pair_features
from matching.library_quality import normalize_inchikey
from matching.recap_adduct import (
    canonical_adduct_key,
    resolve_isopair_adduct_context,
)


def _connectivity(candidate: dict[str, Any]) -> str:
    key = normalize_inchikey(candidate.get("inchikey"))
    return key[:14] if key else ""


def _formula(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or ""))


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _row_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    comparison = candidate_comparison_view(candidate)
    ppm = candidate.get("mass_ppm")
    try:
        ppm_abs = abs(float(ppm)) if ppm is not None else 1e9
    except (TypeError, ValueError):
        ppm_abs = 1e9
    return (
        _safe_float(comparison.get("rank_score", comparison.get("final_score"))),
        _safe_float(comparison.get("cosine")),
        _safe_float(comparison.get("reverse_match")),
        -ppm_abs,
        str(candidate.get("lib_spectrum_id") or ""),
    )


@lru_cache(maxsize=4)
def _load_artifact(
    model_path_text: str,
    config_path_text: str,
) -> tuple[Any, dict[str, Any]]:
    import lightgbm as lgb

    model_path = Path(model_path_text)
    config_path = Path(config_path_text)
    if not model_path.exists() or not config_path.exists():
        raise FileNotFoundError("IsoPair-MS model or runtime config is missing")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    expected = list(config.get("feature_names") or [])
    if expected != list(FEATURE_NAMES):
        raise ValueError("IsoPair-MS feature schema does not match runtime code")
    policy = config.get("runtime_policy") or {}
    if (
        policy.get("high_precision_gate")
        != "filter_existing_msms_high_only"
        or policy.get("high_promotion_allowed") is not False
        or policy.get("stereochemistry_claimed") is not False
    ):
        raise ValueError("IsoPair-MS artifact lacks the locked safe runtime policy")
    model = lgb.Booster(model_file=str(model_path))
    if list(model.feature_name()) != list(FEATURE_NAMES):
        raise ValueError("IsoPair-MS model feature order is incompatible")
    return model, config


def artifact_available(
    model_path: Path = ISOPAIR_MODEL_PATH,
    config_path: Path = ISOPAIR_CONFIG_PATH,
) -> bool:
    try:
        _load_artifact(str(model_path), str(config_path))
        return True
    except Exception:
        return False


def _not_applicable(
    reason: str,
    *,
    formula: str = "",
    candidate_count: int = 0,
    adduct_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    evidence = {
        "schema_version": 2,
        "status": "not_applicable",
        "reason": reason,
        "formula": formula,
        "connectivity_candidate_count": int(candidate_count),
        "evaluated": False,
        "rank_influence_eligible": False,
        "runtime_influence": False,
        "ranking_changed": False,
        "high_precision_gate_applicable": False,
        "selective_threshold_pass": False,
        "high_precision_support": False,
        "high_eligible": False,
        "high_promotion_allowed": False,
        "identity_scope": "2d_connectivity_only",
        "stereochemistry_claimed": False,
    }
    if adduct_context is not None:
        evidence["adduct_context"] = adduct_context
    return evidence


def _assign_summary(
    rows: Iterable[dict[str, Any]],
    evidence: dict[str, Any],
) -> None:
    for row in rows:
        row["isopair_evidence"] = evidence


def apply_isopair_ranker(
    provisional_best: dict[str, Any] | None,
    results: list[dict[str, Any]],
    query_peaks: Iterable[Iterable[Any]],
    *,
    query_precursor_mz: Any,
    query_adduct: Any = None,
    query_adduct_reliable: bool | None = None,
    query_polarity: Any = None,
    query_charge: Any = None,
    query_collision_energy: Any = None,
    query_collision_energy_unit: str | None = None,
    query_instrument_family: str | None = None,
    model_path: Path = ISOPAIR_MODEL_PATH,
    config_path: Path = ISOPAIR_CONFIG_PATH,
    max_baseline_margin: float = ISOPAIR_MAX_BASELINE_MARGIN,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    """Return ``(selected_best, evidence)`` without mutating base scores."""
    if not provisional_best or not results:
        return provisional_best, _not_applicable("no_provisional_match")
    formula = _formula(provisional_best.get("formula"))
    if not formula:
        return provisional_best, _not_applicable("missing_formula")
    if _safe_float(query_precursor_mz) <= 0.0 or _safe_float(query_precursor_mz) > 1000.0:
        return provisional_best, _not_applicable(
            "query_precursor_outside_training_domain", formula=formula
        )

    adduct_context = resolve_isopair_adduct_context(
        formula=formula,
        query_precursor_mz=query_precursor_mz,
        query_adduct=query_adduct,
        query_adduct_reliable=query_adduct_reliable,
        query_polarity=query_polarity,
        query_charge=query_charge,
        candidates=results,
    )
    if not adduct_context.get("ranking_eligible"):
        return provisional_best, _not_applicable(
            str(adduct_context.get("reason") or "adduct_context_unresolved"),
            formula=formula,
            adduct_context=adduct_context,
        )
    query_adduct_key = str(adduct_context.get("adduct_key") or "")

    allowed_tiers = {str(value).lower() for value in ISOPAIR_ALLOWED_SOURCE_TIERS}
    group = [
        result for result in results
        if _formula(result.get("formula")) == formula
        and canonical_adduct_key(result.get("adduct")) == query_adduct_key
        and str(result.get("source_tier") or "").strip().lower() in allowed_tiers
        and _connectivity(result)
    ]
    rows_by_connectivity: dict[str, list[dict[str, Any]]] = {}
    for row in group:
        rows_by_connectivity.setdefault(_connectivity(row), []).append(row)
    candidate_count = len(rows_by_connectivity)
    if not ISOPAIR_MIN_CONNECTIVITIES <= candidate_count <= ISOPAIR_MAX_CONNECTIVITIES:
        evidence = _not_applicable(
            "connectivity_candidate_count_outside_training_domain",
            formula=formula,
            candidate_count=candidate_count,
            adduct_context=adduct_context,
        )
        _assign_summary(group, evidence)
        return provisional_best, evidence

    provisional_connectivity = _connectivity(provisional_best)
    if provisional_connectivity not in rows_by_connectivity:
        evidence = _not_applicable(
            "provisional_winner_outside_applicable_source_or_adduct_stratum",
            formula=formula,
            candidate_count=candidate_count,
            adduct_context=adduct_context,
        )
        _assign_summary(group, evidence)
        return provisional_best, evidence
    try:
        model, config = _load_artifact(str(model_path), str(config_path))
    except Exception as exc:
        evidence = _not_applicable(
            "model_unavailable_or_incompatible",
            formula=formula,
            candidate_count=candidate_count,
            adduct_context=adduct_context,
        )
        evidence["error"] = str(exc)[:160]
        _assign_summary(group, evidence)
        return provisional_best, evidence

    feature_rows: list[np.ndarray] = []
    scored_rows: list[dict[str, Any]] = []
    for connectivity in sorted(rows_by_connectivity):
        # Training used at most four reference spectra per connectivity.
        references = sorted(
            rows_by_connectivity[connectivity], key=_row_key, reverse=True
        )[:4]
        for reference in references:
            feature_rows.append(spectral_pair_features(
                query_peaks,
                reference.get("lib_peaks") or [],
                query_precursor_mz=query_precursor_mz,
                reference_precursor_mz=reference.get("lib_precursor_mz"),
                query_collision_energy=query_collision_energy,
                reference_collision_energy=reference.get("lib_collision_energy"),
                query_collision_energy_unit=query_collision_energy_unit,
                reference_collision_energy_unit=reference.get(
                    "lib_collision_energy_unit"
                ),
                query_instrument_family=query_instrument_family,
                reference_instrument_family=reference.get("lib_instrument_type"),
            ))
            scored_rows.append(reference)
    probabilities = model.predict(np.asarray(feature_rows, dtype=np.float32))
    probability_by_connectivity: dict[str, float] = {}
    for reference, probability in zip(scored_rows, probabilities):
        connectivity = _connectivity(reference)
        probability_by_connectivity[connectivity] = max(
            probability_by_connectivity.get(connectivity, 0.0),
            _safe_float(probability),
        )
    baseline_by_connectivity = {
        connectivity: max(
            _safe_float(candidate_comparison_view(row).get("cosine"))
            for row in rows
        )
        for connectivity, rows in rows_by_connectivity.items()
    }
    alpha = float(config["fusion_alpha_baseline"])
    fused_by_connectivity = {
        connectivity: (
            alpha * baseline_by_connectivity[connectivity]
            + (1.0 - alpha) * probability_by_connectivity[connectivity]
        )
        for connectivity in rows_by_connectivity
    }
    baseline_ranking = sorted(
        rows_by_connectivity,
        key=lambda connectivity: (
            -baseline_by_connectivity[connectivity], connectivity
        ),
    )
    fused_ranking = sorted(
        rows_by_connectivity,
        key=lambda connectivity: (
            -fused_by_connectivity[connectivity], connectivity
        ),
    )
    baseline_margin = (
        baseline_by_connectivity[baseline_ranking[0]]
        - baseline_by_connectivity[baseline_ranking[1]]
    )
    counterfactual_influence = baseline_margin <= float(max_baseline_margin)
    # Formula/mass compatibility is conditional on the retrieved formula, not
    # independent query evidence.  Locked public-library cohorts supported a
    # positive counterfactual, but the MIX33 standard replay reduced exact
    # Top-1 and replicate concordance.  Keep this branch fully evaluated in
    # shadow until independent prospective validation supports deployment.
    production_rank_influence_allowed = bool(
        ISOPAIR_RANK_INFLUENCE_ENABLED
        and not adduct_context.get("formula_conditioned")
    )
    influence = bool(
        counterfactual_influence and production_rank_influence_allowed
    )
    counterfactual_selected_connectivity = (
        fused_ranking[0]
        if counterfactual_influence else provisional_connectivity
    )
    selected_connectivity = (
        fused_ranking[0] if influence else provisional_connectivity
    )
    selected_rows = rows_by_connectivity[selected_connectivity]
    selected = (
        provisional_best
        if selected_connectivity == provisional_connectivity
        else max(selected_rows, key=_row_key)
    )
    ranking_changed = bool(
        influence and selected_connectivity != provisional_connectivity
    )
    counterfactual_ranking_changed = bool(
        counterfactual_influence
        and counterfactual_selected_connectivity != provisional_connectivity
    )
    fused_margin = (
        fused_by_connectivity[fused_ranking[0]]
        - fused_by_connectivity[fused_ranking[1]]
    )
    thresholds = config.get("selective_high_thresholds") or {}
    selective_threshold_pass = bool(
        selected_connectivity == fused_ranking[0]
        and fused_by_connectivity[fused_ranking[0]]
        >= float(thresholds.get("score_threshold", 1.0))
        and fused_margin >= float(thresholds.get("margin_threshold", 1.0))
    )
    counterfactual_selective_threshold_pass = bool(
        counterfactual_selected_connectivity == fused_ranking[0]
        and fused_by_connectivity[fused_ranking[0]]
        >= float(thresholds.get("score_threshold", 1.0))
        and fused_margin >= float(thresholds.get("margin_threshold", 1.0))
    )
    # The selective thresholds were locked over the complete formula-disjoint
    # validation cohort.  They are reliable as a precision filter for an
    # otherwise High, clear winner, but not as a stand-alone promotion rule in
    # the near-tie subset.  A connectivity changed by this module therefore
    # remains Medium/isomer-ambiguous until orthogonal evidence is available.
    high_precision_gate_applicable = bool(
        adduct_context.get("high_precision_gate_eligible")
    )
    high_precision_support = bool(
        high_precision_gate_applicable
        and selective_threshold_pass
        and not ranking_changed
    )
    candidate_scores = [
        {
            "connectivity": connectivity,
            "baseline_modified_cosine": round(
                baseline_by_connectivity[connectivity], 6
            ),
            "same_identity_probability": round(
                probability_by_connectivity[connectivity], 6
            ),
            "fused_score": round(fused_by_connectivity[connectivity], 6),
            "reference_spectra_scored": min(
                len(rows_by_connectivity[connectivity]), 4
            ),
        }
        for connectivity in fused_ranking
    ]
    if adduct_context.get("formula_conditioned"):
        status = (
            "formula_conditioned_shadow_near_tie"
            if counterfactual_influence
            else "formula_conditioned_shadow_clear_baseline_winner"
        )
    elif not ISOPAIR_RANK_INFLUENCE_ENABLED:
        status = (
            "explicit_adduct_shadow_near_tie"
            if counterfactual_influence
            else "explicit_adduct_shadow_clear_baseline_winner"
        )
    else:
        status = "applied" if influence else "clear_baseline_winner_no_rerank"
    evidence = {
        "schema_version": 2,
        "status": status,
        "model_name": str(config.get("model_name") or "IsoPair-MS"),
        "evaluated": True,
        "formula": formula,
        "adduct": str(adduct_context.get("resolved_adduct") or ""),
        "adduct_context": adduct_context,
        "connectivity_candidate_count": candidate_count,
        "baseline_winner_connectivity": baseline_ranking[0],
        "provisional_winner_connectivity": provisional_connectivity,
        "selected_connectivity": selected_connectivity,
        "baseline_margin": round(baseline_margin, 6),
        "maximum_influence_margin": float(max_baseline_margin),
        "fused_winner_score": round(
            fused_by_connectivity[fused_ranking[0]], 6
        ),
        "fused_winner_margin": round(fused_margin, 6),
        "fused_winner_connectivity": fused_ranking[0],
        "rank_influence_eligible": influence,
        "counterfactual_rank_influence_eligible": counterfactual_influence,
        "production_rank_influence_allowed": (
            production_rank_influence_allowed
        ),
        "rank_influence_feature_flag": bool(
            ISOPAIR_RANK_INFLUENCE_ENABLED
        ),
        "runtime_influence": influence,
        "ranking_changed": ranking_changed,
        "counterfactual_selected_connectivity": (
            counterfactual_selected_connectivity
        ),
        "counterfactual_ranking_changed": counterfactual_ranking_changed,
        "formula_conditioned_shadow": bool(
            adduct_context.get("formula_conditioned")
        ),
        "production_abstention_reason": (
            "formula_conditioned_adduct_not_independent_and_mix33_validation_regressed"
            if adduct_context.get("formula_conditioned")
            else (
                "independent_family_confirmation_gain_not_statistically_qualified"
                if not ISOPAIR_RANK_INFLUENCE_ENABLED else None
            )
        ),
        "high_precision_gate_applicable": high_precision_gate_applicable,
        "selective_threshold_pass": selective_threshold_pass,
        "counterfactual_selective_threshold_pass": (
            counterfactual_selective_threshold_pass
        ),
        "high_precision_support": high_precision_support,
        # Backward-compatible alias.  It is deliberately not a permission to
        # promote a Medium near-tie result to High.
        "high_eligible": high_precision_support,
        "high_promotion_allowed": False,
        "high_score_threshold": thresholds.get("score_threshold"),
        "high_margin_threshold": thresholds.get("margin_threshold"),
        "fusion_alpha_baseline": alpha,
        "fusion_alpha_metric_model": 1.0 - alpha,
        "candidate_scores": candidate_scores,
        "identity_scope": "2d_connectivity_only",
        "stereochemistry_claimed": False,
        "applicability_domain": (
            "same_formula_same_adduct_2_to_12_connectivities_"
            "experimental_reference_spectra"
        ),
    }
    for connectivity, rows in rows_by_connectivity.items():
        for row in rows:
            row["isopair_probability"] = round(
                probability_by_connectivity[connectivity], 6
            )
            row["isopair_fused_score"] = round(
                fused_by_connectivity[connectivity], 6
            )
            row["isopair_rank_applied"] = influence
            row["isopair_shadow_rank_applied"] = counterfactual_influence
            row["isopair_selected_connectivity"] = (
                connectivity == selected_connectivity
            )
            row["isopair_counterfactual_selected_connectivity"] = (
                connectivity == counterfactual_selected_connectivity
            )
            row["isopair_high_precision_support"] = bool(
                high_precision_support and connectivity == selected_connectivity
            )
            row["isopair_high_eligible"] = row[
                "isopair_high_precision_support"
            ]
            row["isopair_evidence"] = evidence
    return selected, evidence


ISOPAIR_RUNTIME_AVAILABLE = True
