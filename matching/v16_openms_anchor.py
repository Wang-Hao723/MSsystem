"""Fail-closed validation of OpenMS adduct groups for the V16 shadow graph.

The feature manifest is candidate-independent.  This module validates the
measured feature/adduct relationships only; it does not assign a compound
identity, modify a V15 score, or increment a compound count.
"""
from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = "mssystem.v16.openms-adduct-group-validation/1.0.0"
_REQUIRED_PIPELINE_STAGES = {
    "MassTraceDetection",
    "ElutionPeakDetection",
    "FeatureFindingMetabo",
    "MetaboliteFeatureDeconvolution",
}


def _text(value: Any) -> str:
    return str(value or "").strip()


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _adduct_polarity(value: Any) -> str:
    text = "".join(_text(value).split())
    if text.endswith("+"):
        return "POS"
    if text.endswith("-"):
        return "NEG"
    return ""


def _manifest_rejection(
    reasons: Sequence[str],
    *,
    task_polarity: str,
    manifest_sha256: str | None,
    feature_xml_sha256: str | None,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "rejected_manifest",
        "task_polarity": task_polarity,
        "manifest_sha256": manifest_sha256,
        "feature_xml_sha256": feature_xml_sha256,
        "manifest_rejection_reasons": sorted(set(reasons)),
        "eligible_group_count": 0,
        "eligible_feature_count": 0,
        "rejected_group_count": 0,
        "rejected_group_reason_counts": {},
        "evidence_by_feature_id": {},
        "candidate_identity_used": False,
        "production_runtime_influence": False,
    }


def validate_openms_adduct_group_manifest(
    manifest: Mapping[str, Any],
    *,
    task_polarity: str,
    manifest_sha256: str | None = None,
    feature_xml_sha256: str | None = None,
    feature_xml_path: str | Path | None = None,
    requirements: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return feature-bound shadow evidence plus a compact validation audit."""

    polarity = _text(task_polarity).upper()
    rules = dict(requirements or {})
    ppm = float(rules.get("neutral_mass_tolerance_ppm", 10.0))
    da_floor = float(rules.get("neutral_mass_tolerance_da_floor", 0.001))
    rt_limit = float(rules.get("maximum_member_rt_delta_seconds", 10.0))
    minimum_members = int(rules.get("minimum_member_feature_count", 2))
    minimum_adducts = int(rules.get("minimum_distinct_adduct_count", 2))

    rejection_reasons: list[str] = []
    if polarity not in {"POS", "NEG"}:
        rejection_reasons.append("task_polarity_not_exact_pos_or_neg")
    if _text(manifest.get("status")).lower() != "ok":
        rejection_reasons.append("feature_manifest_status_not_ok")
    if not _text(manifest.get("backend")).startswith("PyOpenMS "):
        rejection_reasons.append("backend_not_pyopenms")
    pipeline = {_text(value) for value in manifest.get("pipeline") or []}
    if not _REQUIRED_PIPELINE_STAGES.issubset(pipeline):
        rejection_reasons.append("required_openms_pipeline_stage_missing")
    if _text(manifest.get("adduct_deconvolution_status")).lower() != "ok":
        rejection_reasons.append("adduct_deconvolution_status_not_ok")
    parameters = manifest.get("parameters") or {}
    if not isinstance(parameters, Mapping):
        parameters = {}
        rejection_reasons.append("parameters_not_an_object")
    if parameters.get("group_adducts") is not True:
        rejection_reasons.append("group_adducts_not_true")
    if _text(parameters.get("ion_mode")).upper() != polarity:
        rejection_reasons.append("manifest_ion_mode_mismatch")
    if feature_xml_path is not None:
        declared_output = _text(
            manifest.get("output_features") or manifest.get("feature_xml")
        )
        if not declared_output:
            rejection_reasons.append("feature_xml_binding_missing")
        else:
            try:
                declared = Path(declared_output).expanduser().resolve()
                expected = Path(feature_xml_path).expanduser().resolve()
            except (OSError, RuntimeError, ValueError):
                rejection_reasons.append("feature_xml_binding_invalid")
            else:
                if declared != expected:
                    rejection_reasons.append("feature_xml_binding_mismatch")

    features = manifest.get("features") or []
    if not isinstance(features, list):
        features = []
        rejection_reasons.append("features_not_an_array")
    declared_feature_count = manifest.get("feature_count")
    try:
        count_matches = int(declared_feature_count) == len(features)
    except (TypeError, ValueError):
        count_matches = False
    if not count_matches:
        rejection_reasons.append("feature_count_mismatch")

    feature_ids = [
        _text(feature.get("feature_id"))
        for feature in features
        if isinstance(feature, Mapping)
    ]
    populated_ids = [value for value in feature_ids if value]
    if len(populated_ids) != len(set(populated_ids)):
        rejection_reasons.append("duplicate_feature_id")
    if rejection_reasons:
        return _manifest_rejection(
            rejection_reasons,
            task_polarity=polarity,
            manifest_sha256=manifest_sha256,
            feature_xml_sha256=feature_xml_sha256,
        )

    groups: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for feature in features:
        if not isinstance(feature, Mapping):
            continue
        group_id = _text(feature.get("adduct_group_id"))
        if group_id:
            groups[group_id].append(feature)

    evidence_by_feature_id: dict[str, dict[str, Any]] = {}
    rejected_reason_counts: Counter[str] = Counter()
    eligible_groups = 0
    for group_id, members in sorted(groups.items()):
        group_reasons: list[str] = []
        member_ids = [_text(member.get("feature_id")) for member in members]
        member_ids = sorted({value for value in member_ids if value})
        adducts = [_text(member.get("adduct")) for member in members]
        neutral_masses = [_finite(member.get("neutral_mass")) for member in members]
        retention_times = [_finite(member.get("rt_seconds")) for member in members]
        if len(member_ids) < minimum_members:
            group_reasons.append("insufficient_distinct_members")
        if any(not value for value in adducts):
            group_reasons.append("missing_explicit_adduct")
        if any(_adduct_polarity(value) != polarity for value in adducts):
            group_reasons.append("adduct_sign_mismatch")
        if len(set(adducts)) < minimum_adducts:
            group_reasons.append("insufficient_distinct_adducts")
        if any(value is None or value <= 0 for value in neutral_masses):
            group_reasons.append("invalid_neutral_mass")
        if any(value is None for value in retention_times):
            group_reasons.append("missing_retention_time")

        usable_masses = [float(value) for value in neutral_masses if value is not None]
        centre = statistics.median(usable_masses) if usable_masses else None
        mass_tolerance = (
            max(da_floor, abs(centre) * ppm * 1e-6)
            if centre is not None
            else None
        )
        if (
            mass_tolerance is not None
            and max(usable_masses) - min(usable_masses) > mass_tolerance
        ):
            group_reasons.append("neutral_mass_spread_exceeds_limit")
        usable_rts = [float(value) for value in retention_times if value is not None]
        if usable_rts and max(usable_rts) - min(usable_rts) > rt_limit:
            group_reasons.append("retention_time_spread_exceeds_limit")

        if group_reasons:
            rejected_reason_counts.update(set(group_reasons))
            continue
        eligible_groups += 1
        for member in members:
            feature_id = _text(member.get("feature_id"))
            evidence_by_feature_id[feature_id] = {
                "status": "validated",
                "network_id": group_id,
                "neutral_mass": round(float(centre), 10),
                "ion_identity": _text(member.get("adduct")),
                "runtime_influence": True,
                "production_runtime_influence": False,
                "identity_conflict": False,
                "neutral_mass_consistent": True,
                "member_feature_ids": member_ids,
                "member_feature_count": len(member_ids),
                "manifest_sha256": manifest_sha256,
                "feature_xml_sha256": feature_xml_sha256,
                "candidate_identity_used": False,
            }

    return {
        "schema_version": SCHEMA_VERSION,
        "status": "validated",
        "task_polarity": polarity,
        "manifest_sha256": manifest_sha256,
        "feature_xml_sha256": feature_xml_sha256,
        "manifest_rejection_reasons": [],
        "input_feature_count": len(features),
        "input_group_count": len(groups),
        "eligible_group_count": eligible_groups,
        "eligible_feature_count": len(evidence_by_feature_id),
        "rejected_group_count": len(groups) - eligible_groups,
        "rejected_group_reason_counts": dict(sorted(rejected_reason_counts.items())),
        "validation_limits": {
            "neutral_mass_tolerance_ppm": ppm,
            "neutral_mass_tolerance_da_floor": da_floor,
            "maximum_member_rt_delta_seconds": rt_limit,
            "minimum_member_feature_count": minimum_members,
            "minimum_distinct_adduct_count": minimum_adducts,
        },
        "evidence_by_feature_id": evidence_by_feature_id,
        "candidate_identity_used": False,
        "production_runtime_influence": False,
    }


def compact_openms_validation_audit(
    validation: Mapping[str, Any],
) -> dict[str, Any]:
    """Remove per-feature evidence before embedding an audit in a sidecar."""

    return {
        key: value
        for key, value in validation.items()
        if key != "evidence_by_feature_id"
    }
