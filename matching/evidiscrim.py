"""EviDiscrim-MS candidate-contrastive diagnostic evidence (shadow only).

This module asks a narrower question than another whole-spectrum reranker:
*which predicted ions actually differ between the unresolved candidates, how
uncertain are those differences, does the measured spectrum support either
direction, and at which instrument-native collision energy is the hardest pair
predicted to be most separable?*

All outputs are candidate-conditioned counterfactuals.  They never mutate the
production result and they are not posterior probabilities or independent
confirmation.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from matching.library_quality import compound_identity_key, normalize_formula_text
from schemas.evidiscrim import (
    CandidateContrastSupport,
    ContrastivePair,
    DiagnosticIon,
    EnergyPairSeparation,
    EnergySelection,
    EVIDISCRIM_MODULE_ID,
    EviDiscrimArtifact,
    EviDiscrimSpectrum,
)


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def results_sha256(results: Sequence[Mapping[str, Any]]) -> str:
    """Digest the complete legacy result without changing or normalizing it."""

    return hashlib.sha256(_canonical_bytes(results)).hexdigest()


def stream_select_evidiscrim_rows(
    path: Path,
    *,
    max_spectra: int,
    max_candidates_per_spectrum: int,
    supported_adducts: frozenset[str],
    policy: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], str, dict[str, Any]]:
    """Hash a large result array while retaining only bounded applicable rows.

    Streaming numbers use the same Python types as ``json.loads`` so this
    canonical digest matches :func:`results_sha256` without materializing the
    complete task in memory.
    """

    import ijson

    selected: list[dict[str, Any]] = []
    selected_indices: list[int] = []
    seen_pair_signatures: set[tuple[tuple[str, str], ...]] = set()
    duplicate_signatures_skipped = 0
    total_rows = 0
    digest = hashlib.sha256()
    digest.update(b"[")
    with path.open("rb") as handle:
        rows = ijson.items(handle, "item", use_float=True)
        for row_index, raw_row in enumerate(rows):
            if row_index:
                digest.update(b",")
            digest.update(_canonical_bytes(raw_row))
            total_rows += 1
            if len(selected) >= max(1, int(max_spectra)):
                continue
            if not isinstance(raw_row, Mapping):
                continue
            match = raw_row.get("match")
            if not isinstance(match, Mapping):
                continue
            if _number(raw_row.get("precursor_mz")) is None:
                continue
            if _number(raw_row.get("collision_energy")) is None:
                continue
            if not _observed_peaks(raw_row):
                continue
            candidates = _candidate_rows(
                match, max(2, int(max_candidates_per_spectrum))
            )
            candidates = [
                candidate for candidate in candidates
                if str(candidate.get("smiles") or "").strip()
                and str(candidate.get("adduct") or "").strip()
                in supported_adducts
            ]
            if len(candidates) < 2:
                continue
            structures = {
                str(candidate["_identity_key"]): _structure_profile(candidate)
                for candidate in candidates
            }
            eligible_pairs = [
                tuple(sorted((
                    str(left["_identity_key"]),
                    str(right["_identity_key"]),
                )))
                for left_index, left in enumerate(candidates)
                for right in candidates[left_index + 1:]
                if _pair_scope(
                    left,
                    right,
                    structures[str(left["_identity_key"])],
                    structures[str(right["_identity_key"])],
                    policy,
                )[0]
            ]
            if not eligible_pairs:
                continue
            pair_signature = tuple(sorted(eligible_pairs))
            bounded = policy.get("bounded_execution") or {}
            if (
                bounded.get("deduplicate_candidate_pair_signatures", True)
                and pair_signature in seen_pair_signatures
            ):
                duplicate_signatures_skipped += 1
                continue
            seen_pair_signatures.add(pair_signature)
            row = dict(raw_row)
            row["_evidiscrim_source_result_index"] = row_index
            selected.append(row)
            selected_indices.append(row_index)
    digest.update(b"]")
    return selected, digest.hexdigest(), {
        "mode": "streaming_bounded_selection_with_full_canonical_digest",
        "total_source_rows": total_rows,
        "selected_rows": len(selected),
        "selected_source_result_indices": selected_indices,
        "maximum_spectra": max(1, int(max_spectra)),
        "maximum_candidates_per_spectrum": max(
            2, int(max_candidates_per_spectrum)
        ),
        "selection_strategy": str(
            (policy.get("bounded_execution") or {}).get(
                "selection_strategy"
            ) or "source_order_unique_eligible_candidate_pair_signature"
        ),
        "duplicate_candidate_pair_signatures_skipped": (
            duplicate_signatures_skipped
        ),
    }


def load_evidiscrim_policy(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not payload.get("policy_id"):
        raise ValueError("EviDiscrim policy must be a versioned JSON object")
    if payload.get("module_state") != "shadow":
        raise ValueError("EviDiscrim policy must remain shadow")
    energy = payload.get("energy_selection") or {}
    if energy.get("cross_instrument_conversion_allowed"):
        raise ValueError("Cross-instrument collision-energy conversion is forbidden")
    if energy.get("universal_complete_fragmentation_assumed"):
        raise ValueError("A universal complete-fragmentation energy is forbidden")
    uncertainty = payload.get("prediction_uncertainty") or {}
    if "singleton_prediction_uncertainty_floor" in uncertainty:
        raise ValueError("Singleton uncertainty floors are forbidden in EviDiscrim v2")
    if int(uncertainty.get("minimum_independent_architecture_families", 0)) < 2:
        raise ValueError(
            "At least two independent architecture families are required"
        )
    return payload


def evidiscrim_artifact_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}_evidiscrim_shadow.json"


def read_evidiscrim_artifact(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("EviDiscrim artifact must be a JSON object")
    # Re-validate hard shadow invariants at every consumer boundary.
    return EviDiscrimArtifact.model_validate(payload).model_dump(mode="json")


def evidiscrim_artifact_summary(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {
            "requested": False,
            "applicable": False,
            "executed": False,
            "status": "not_started",
            "module_state": "shadow",
            "production_ranking_changed": False,
        }
    try:
        payload = read_evidiscrim_artifact(path)
    except (OSError, ValueError) as exc:
        return {
            "requested": True,
            "applicable": False,
            "executed": False,
            "status": "invalid_artifact",
            "reason": f"{type(exc).__name__}: {exc}"[:500],
            "module_state": "shadow",
            "production_ranking_changed": False,
        }
    summary = payload.get("summary") or {}
    return {
        "requested": True,
        "applicable": bool(summary.get("spectra_evaluated")),
        "executed": payload.get("status") == "completed_shadow",
        "status": payload.get("status"),
        "spectra_evaluated": int(summary.get("spectra_evaluated") or 0),
        "eligible_connectivity_pairs": int(
            summary.get("eligible_connectivity_pairs") or 0
        ),
        "resolved_counterfactual_pairs": int(
            summary.get("resolved_counterfactual_pairs") or 0
        ),
        "abstained_eligible_pairs": int(
            summary.get("abstained_eligible_pairs") or 0
        ),
        "ensemble_qualified_pairs": int(
            summary.get("ensemble_qualified_pairs") or 0
        ),
        "insufficient_architecture_pairs": int(
            summary.get("insufficient_architecture_pairs") or 0
        ),
        "same_connectivity_out_of_scope_pairs": int(
            summary.get("same_connectivity_out_of_scope_pairs") or 0
        ),
        "collision_energy_recommendations": int(
            summary.get("collision_energy_recommendations") or 0
        ),
        "instrument_domain_shifted_energy_audits": int(
            summary.get("instrument_domain_shifted_energy_audits") or 0
        ),
        "counterfactual_candidate_order_changed_spectra": int(
            summary.get("counterfactual_candidate_order_changed_spectra") or 0
        ),
        "validation_status": summary.get("validation_status"),
        "module_state": "shadow",
        "production_ranking_changed": False,
        "confidence_changed": False,
        "high_promotions": 0,
    }


def _number(value: Any) -> float | None:
    if isinstance(value, Mapping):
        value = value.get("value")
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _unit(value: Any) -> str:
    text = str(value or "").strip()
    return text or "instrument_native_unspecified"


def _candidate_rows(match: Mapping[str, Any], limit: int) -> list[dict[str, Any]]:
    rows: list[Any] = [match, *(match.get("top_candidates") or [])]
    passport = match.get("evidence_governance") or {}
    rows.extend(passport.get("ambiguity_set") or [])
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in rows:
        if not isinstance(raw, Mapping):
            continue
        candidate = dict(raw)
        identity = compound_identity_key(candidate)
        if not identity or identity in seen:
            continue
        seen.add(identity)
        candidate["_identity_key"] = identity
        candidate["_production_top1"] = len(output) == 0
        output.append(candidate)
        if len(output) >= max(2, int(limit)):
            break
    return output


def source_row_sha256(
    row: Mapping[str, Any], *, max_candidates_per_spectrum: int
) -> str:
    """Fingerprint only fields that determine pair evidence and row linkage."""

    match = row.get("match") if isinstance(row, Mapping) else None
    candidates = (
        _candidate_rows(match, max(2, int(max_candidates_per_spectrum)))
        if isinstance(match, Mapping) else []
    )
    payload = {
        "spectrum_id": str(row.get("spectrum_id") or ""),
        "precursor_mz": _number(row.get("precursor_mz")),
        "retention_time": _number(row.get("retention_time")),
        "collision_energy": _number(row.get("collision_energy")),
        "collision_energy_unit": _unit(row.get("collision_energy_unit")),
        "instrument_type": str(row.get("instrument_type") or ""),
        "instrument": str(row.get("instrument") or ""),
        "observed_peaks": [
            [round(mz, 6), round(intensity, 8)]
            for mz, intensity in _observed_peaks(row)
        ],
        "candidates": [{
            "identity_key": str(candidate["_identity_key"]),
            "smiles": str(candidate.get("smiles") or "").strip(),
            "formula": normalize_formula_text(
                candidate.get("formula") or candidate.get("molecular_formula")
            ),
            "adduct": str(candidate.get("adduct") or "").strip(),
        } for candidate in candidates],
    }
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _structure_profile(candidate: Mapping[str, Any]) -> dict[str, Any]:
    smiles = str(candidate.get("smiles") or "").strip()
    declared_formula = normalize_formula_text(
        candidate.get("formula") or candidate.get("molecular_formula")
    )
    if not smiles:
        return {
            "valid": False,
            "reason": "missing_candidate_smiles",
            "declared_formula": declared_formula,
        }
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None
    if mol is None:
        return {
            "valid": False,
            "reason": "invalid_candidate_smiles",
            "declared_formula": declared_formula,
        }
    try:
        isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        connectivity = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        derived_formula = normalize_formula_text(rdMolDescriptors.CalcMolFormula(mol))
    except Exception:
        return {
            "valid": False,
            "reason": "rdkit_structure_normalization_failed",
            "declared_formula": declared_formula,
        }
    if declared_formula and derived_formula and declared_formula != derived_formula:
        return {
            "valid": False,
            "reason": "declared_formula_structure_conflict",
            "declared_formula": declared_formula,
            "structure_formula": derived_formula,
            "canonical_isomeric_smiles": isomeric,
            "canonical_connectivity_smiles": connectivity,
        }
    return {
        "valid": True,
        "reason": "valid_rdkit_structure",
        "declared_formula": declared_formula,
        "structure_formula": derived_formula,
        "normalized_formula": declared_formula or derived_formula,
        "canonical_isomeric_smiles": isomeric,
        "canonical_connectivity_smiles": connectivity,
    }


def _pair_key(left_identity: str, right_identity: str) -> str:
    # Deliberately identical to acquisition.ambiguity_graph edge IDs so the
    # planning sidecar can be joined without using candidate names.
    payload = sorted([left_identity, right_identity])
    digest = hashlib.sha256(_canonical_bytes(payload)).hexdigest()[:16]
    return f"edge_{digest}"


def _pair_scope(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    left_structure: Mapping[str, Any],
    right_structure: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> tuple[bool, tuple[str, ...], str | None]:
    configured = policy.get("candidate_scope") or {}
    if not left_structure.get("valid"):
        return False, (), str(left_structure.get("reason") or "invalid_left_structure")
    if not right_structure.get("valid"):
        return False, (), str(right_structure.get("reason") or "invalid_right_structure")

    reasons: list[str] = []
    left_formula = str(left_structure.get("normalized_formula") or "")
    right_formula = str(right_structure.get("normalized_formula") or "")
    if configured.get("require_same_normalized_formula", True):
        if not left_formula or not right_formula:
            return False, (), "missing_normalized_formula"
        if left_formula != right_formula:
            return False, (), "different_normalized_formula"
        reasons.append("same_normalized_formula")

    left_adduct = str(left.get("adduct") or "").strip()
    right_adduct = str(right.get("adduct") or "").strip()
    if configured.get("require_same_explicit_adduct", True):
        if not left_adduct or not right_adduct:
            return False, tuple(reasons), "missing_explicit_adduct"
        if left_adduct != right_adduct:
            return False, tuple(reasons), "different_explicit_adduct"
        reasons.append("same_explicit_adduct")

    left_connectivity = str(left_structure.get("canonical_connectivity_smiles") or "")
    right_connectivity = str(right_structure.get("canonical_connectivity_smiles") or "")
    if configured.get("require_different_2d_connectivity", True):
        if left_connectivity == right_connectivity:
            left_iso = str(left_structure.get("canonical_isomeric_smiles") or "")
            right_iso = str(right_structure.get("canonical_isomeric_smiles") or "")
            reason = (
                "same_2d_connectivity_stereoisomer_out_of_scope"
                if left_iso != right_iso
                else "duplicate_same_2d_structure"
            )
            return False, tuple(reasons), reason
        reasons.append("different_2d_connectivity")
    return True, tuple(reasons), None


def _sanitize_peaks(peaks: Any) -> list[tuple[float, float]]:
    output: list[tuple[float, float]] = []
    if not isinstance(peaks, Sequence) or isinstance(peaks, (str, bytes)):
        return output
    for peak in peaks:
        if not isinstance(peak, Sequence) or isinstance(peak, (str, bytes)):
            continue
        if len(peak) < 2:
            continue
        mz = _number(peak[0])
        intensity = _number(peak[1])
        if mz is None or intensity is None or mz <= 0 or intensity <= 0:
            continue
        output.append((mz, intensity))
    if not output:
        return output
    maximum = max(intensity for _, intensity in output)
    if maximum <= 0:
        return []
    normalized = [(mz, intensity / maximum) for mz, intensity in output]
    return sorted(normalized, key=lambda item: (item[0], -item[1]))


def _prediction_members(record: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return prediction members without treating repeated runs as new models.

    ``prediction_members`` is the v2 contract: every member declares its
    provider and base architecture family.  The legacy
    ``prediction_replicates`` field remains readable, but all of those spectra
    inherit one architecture family and therefore cannot satisfy the
    independent-model requirement by themselves.
    """

    inherited_family = str(
        record.get("base_architecture_family")
        or record.get("architecture_family")
        or "legacy_unspecified_single_architecture"
    ).strip()
    inherited_provider = str(
        record.get("provider_id")
        or record.get("model_id")
        or "legacy_unspecified_provider"
    ).strip()

    def member(raw: Mapping[str, Any], peaks: Any, ordinal: int) -> dict[str, Any] | None:
        sanitized = _sanitize_peaks(peaks)
        if not sanitized:
            return None
        family = str(
            raw.get("base_architecture_family")
            or raw.get("architecture_family")
            or inherited_family
        ).strip() or "legacy_unspecified_single_architecture"
        provider = str(
            raw.get("provider_id")
            or raw.get("model_id")
            or inherited_provider
        ).strip() or "legacy_unspecified_provider"
        return {
            "peaks": sanitized,
            "architecture_family": family,
            "provider_id": provider,
            "member_id": str(
                raw.get("member_id")
                or raw.get("checkpoint_id")
                or f"{provider}::replicate-{ordinal}"
            ),
            "energy_transfer_qualified": bool(
                raw.get(
                    "energy_transfer_qualified",
                    record.get("energy_transfer_qualified", False),
                )
            ),
            "instrument_domain": dict(
                raw.get("instrument_domain")
                or record.get("instrument_domain")
                or {}
            ),
            "calibration_layer_id": str(
                raw.get("calibration_layer_id")
                or record.get("calibration_layer_id")
                or ""
            ),
        }

    raw_members = record.get("prediction_members")
    output: list[dict[str, Any]] = []
    if isinstance(raw_members, Sequence) and not isinstance(
        raw_members, (str, bytes)
    ):
        for ordinal, raw in enumerate(raw_members):
            if not isinstance(raw, Mapping):
                continue
            parsed = member(raw, raw.get("predicted_peaks") or [], ordinal)
            if parsed:
                output.append(parsed)
        if output:
            return output

    raw_replicates = record.get("prediction_replicates")
    if isinstance(raw_replicates, Sequence) and not isinstance(
        raw_replicates, (str, bytes)
    ):
        for ordinal, peaks in enumerate(raw_replicates):
            parsed = member(record, peaks, ordinal)
            if parsed:
                output.append(parsed)
        if output:
            return output
    parsed = member(record, record.get("predicted_peaks") or [], 0)
    return [parsed] if parsed else []


def _cluster_centers(
    spectra: Sequence[Sequence[tuple[float, float]]], tolerance: float
) -> list[float]:
    mz_values = sorted(mz for spectrum in spectra for mz, _ in spectrum)
    if not mz_values:
        return []
    clusters: list[list[float]] = [[mz_values[0]]]
    for mz in mz_values[1:]:
        center = sum(clusters[-1]) / len(clusters[-1])
        if abs(mz - center) <= tolerance:
            clusters[-1].append(mz)
        else:
            clusters.append([mz])
    return [sum(cluster) / len(cluster) for cluster in clusters]


def _vector_on_centers(
    peaks: Sequence[tuple[float, float]],
    centers: Sequence[float],
    tolerance: float,
) -> np.ndarray:
    vector = np.zeros(len(centers), dtype=float)
    for mz, intensity in peaks:
        if not centers:
            break
        distances = [abs(mz - center) for center in centers]
        index = int(np.argmin(distances))
        if distances[index] <= tolerance:
            vector[index] += intensity
    return np.clip(vector, 0.0, 1.0)


def _prediction_statistics(
    members: Sequence[Mapping[str, Any]],
    centers: Sequence[float],
    tolerance: float,
    policy: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    by_family: dict[str, list[np.ndarray]] = defaultdict(list)
    providers: set[str] = set()
    for item in members:
        peaks = item.get("peaks") or []
        family = str(
            item.get("architecture_family")
            or "legacy_unspecified_single_architecture"
        )
        by_family[family].append(_vector_on_centers(peaks, centers, tolerance))
        providers.add(str(item.get("provider_id") or "legacy_unspecified_provider"))
    if not by_family:
        zeros = np.zeros(len(centers), dtype=float)
        return zeros, zeros, {
            "qualified": False,
            "architecture_families": (),
            "providers": (),
            "between_architecture_uncertainty": 0.0,
            "within_architecture_uncertainty": 0.0,
            "uncertainty_applied": True,
        }
    family_means: list[np.ndarray] = []
    within_variances: list[np.ndarray] = []
    for family in sorted(by_family):
        matrix = np.vstack(by_family[family])
        family_means.append(np.mean(matrix, axis=0))
        if matrix.shape[0] > 1:
            within_variances.append(np.var(matrix, axis=0, ddof=1))
        else:
            within_variances.append(np.zeros(matrix.shape[1], dtype=float))
    family_matrix = np.vstack(family_means)
    mean = np.mean(family_matrix, axis=0)
    uncertainty = policy.get("prediction_uncertainty") or {}
    minimum_families = max(
        2,
        int(uncertainty.get("minimum_independent_architecture_families", 2)),
    )
    between = (
        np.std(family_matrix, axis=0, ddof=1)
        if family_matrix.shape[0] > 1
        else np.zeros(family_matrix.shape[1], dtype=float)
    )
    within_variance = np.mean(np.vstack(within_variances), axis=0)
    within = np.sqrt(np.clip(within_variance, 0.0, None))
    within_weight = max(
        0.0, float(uncertainty.get("within_architecture_variance_weight", 1.0))
    )
    apply_uncertainty = bool(uncertainty.get("apply_reliability_weight", True))
    if apply_uncertainty:
        sd = np.sqrt(
            between * between + within_weight * within_variance
        )
    else:
        sd = np.zeros(family_matrix.shape[1], dtype=float)
    families = tuple(sorted(by_family))
    return (
        np.clip(mean, 0.0, 1.0),
        np.clip(sd, 0.0, 1.0),
        {
            "qualified": len(families) >= minimum_families,
            "minimum_architecture_families": minimum_families,
            "architecture_families": families,
            "providers": tuple(sorted(providers)),
            "between_architecture_uncertainty": float(np.mean(between)),
            "within_architecture_uncertainty": float(np.mean(within)),
            "uncertainty_applied": apply_uncertainty,
        },
    )


def _jensen_shannon(left: np.ndarray, right: np.ndarray) -> float:
    left_sum = float(np.sum(left))
    right_sum = float(np.sum(right))
    if left_sum <= 0 or right_sum <= 0:
        return 0.0
    p = left / left_sum
    q = right / right_sum
    midpoint = 0.5 * (p + q)

    def divergence(source: np.ndarray) -> float:
        mask = source > 0
        return float(np.sum(source[mask] * np.log2(source[mask] / midpoint[mask])))

    return max(0.0, min(1.0, 0.5 * divergence(p) + 0.5 * divergence(q)))


def _diagnostic_statistics(
    left_members: Sequence[Mapping[str, Any]],
    right_members: Sequence[Mapping[str, Any]],
    observed: Sequence[tuple[float, float]],
    policy: Mapping[str, Any],
) -> tuple[
    list[DiagnosticIon],
    float,
    float,
    float,
    float,
    float,
    dict[str, Any],
]:
    alignment = policy.get("peak_alignment") or {}
    diagnostic = policy.get("diagnostic_evidence") or {}
    tolerance = float(alignment.get("mz_tolerance_da", 0.02))
    left_families = {
        str(item.get("architecture_family") or "") for item in left_members
    }
    right_families = {
        str(item.get("architecture_family") or "") for item in right_members
    }
    common_families = tuple(sorted((left_families & right_families) - {""}))
    left_common = [
        item for item in left_members
        if str(item.get("architecture_family") or "") in common_families
    ]
    right_common = [
        item for item in right_members
        if str(item.get("architecture_family") or "") in common_families
    ]
    centers = _cluster_centers(
        [
            *(item.get("peaks") or [] for item in left_common),
            *(item.get("peaks") or [] for item in right_common),
        ],
        tolerance,
    )
    base_audit = {
        "qualified": False,
        "common_architecture_families": common_families,
        "common_architecture_family_count": len(common_families),
        "left_architecture_families": tuple(sorted(left_families - {""})),
        "right_architecture_families": tuple(sorted(right_families - {""})),
        "left_providers": tuple(sorted({
            str(item.get("provider_id") or "legacy_unspecified_provider")
            for item in left_common
        })),
        "right_providers": tuple(sorted({
            str(item.get("provider_id") or "legacy_unspecified_provider")
            for item in right_common
        })),
        "instrument_domain_qualified": False,
        "between_architecture_uncertainty": 0.0,
        "within_architecture_uncertainty": 0.0,
    }
    if not centers:
        return [], 0.0, 0.0, 0.0, 0.0, 0.0, base_audit
    left_mean, left_sd, left_audit = _prediction_statistics(
        left_common, centers, tolerance, policy
    )
    right_mean, right_sd, right_audit = _prediction_statistics(
        right_common, centers, tolerance, policy
    )
    ensemble_qualified = bool(left_audit["qualified"] and right_audit["qualified"])
    instrument_domain_qualified = bool(left_common and right_common) and all(
        bool(item.get("energy_transfer_qualified", False))
        for item in [*left_common, *right_common]
    )
    ensemble_audit = {
        **base_audit,
        "qualified": ensemble_qualified,
        "instrument_domain_qualified": instrument_domain_qualified,
        "between_architecture_uncertainty": float(np.mean([
            left_audit["between_architecture_uncertainty"],
            right_audit["between_architecture_uncertainty"],
        ])),
        "within_architecture_uncertainty": float(np.mean([
            left_audit["within_architecture_uncertainty"],
            right_audit["within_architecture_uncertainty"],
        ])),
        "uncertainty_applied": bool(
            left_audit["uncertainty_applied"]
            and right_audit["uncertainty_applied"]
        ),
    }
    observed_vector = _vector_on_centers(observed, centers, tolerance)
    minimum_delta = float(
        diagnostic.get("minimum_absolute_prediction_delta", 0.12)
    )
    minimum_intensity = float(alignment.get("minimum_relative_intensity", 0.01))
    presence_threshold = float(diagnostic.get("observed_presence_threshold", 0.01))
    missing_weight = float(
        diagnostic.get("missing_observed_peak_penalty_weight", 0.35)
    )
    maximum_ions = max(1, int(alignment.get("maximum_diagnostic_ions_per_pair", 12)))
    ions: list[DiagnosticIon] = []
    for index, mz in enumerate(centers):
        delta = abs(float(left_mean[index] - right_mean[index]))
        if delta < minimum_delta:
            continue
        if max(float(left_mean[index]), float(right_mean[index])) < minimum_intensity:
            continue
        total_uncertainty = float(left_sd[index] + right_sd[index])
        reliability = delta / max(1e-12, delta + total_uncertainty)
        weight = max(0.0, min(1.0, delta * reliability))
        observed_intensity = float(observed_vector[index])
        present = observed_intensity >= presence_threshold
        absence_weight = 1.0 if present else missing_weight
        if left_mean[index] > right_mean[index]:
            direction = "left"
        elif right_mean[index] > left_mean[index]:
            direction = "right"
        else:
            direction = "indeterminate"
        ions.append(DiagnosticIon(
            mz=round(float(mz), 6),
            left_predicted_intensity=round(float(left_mean[index]), 6),
            right_predicted_intensity=round(float(right_mean[index]), 6),
            left_prediction_uncertainty=round(float(left_sd[index]), 6),
            right_prediction_uncertainty=round(float(right_sd[index]), 6),
            absolute_prediction_delta=round(delta, 6),
            reliability=round(reliability, 6),
            diagnostic_weight=round(weight, 6),
            observed_relative_intensity=round(observed_intensity, 6),
            observed_present=present,
            absence_penalty_weight=round(absence_weight, 6),
            predicted_direction=direction,
        ))
    ions.sort(key=lambda ion: (-ion.diagnostic_weight, ion.mz))
    ions = ions[:maximum_ions]
    if not ions:
        raw_js = _jensen_shannon(left_mean, right_mean)
        reliability_factor = (
            max(0.0, 1.0 - float(np.mean(0.5 * (left_sd + right_sd))))
            if ensemble_qualified else 0.0
        )
        return (
            [], 0.0, 0.0, 0.0, raw_js, reliability_factor, ensemble_audit
        )

    weighted_left_error = 0.0
    weighted_right_error = 0.0
    denominator = 0.0
    for ion in ions:
        effective_weight = ion.diagnostic_weight * ion.absence_penalty_weight
        weighted_left_error += effective_weight * abs(
            ion.observed_relative_intensity - ion.left_predicted_intensity
        )
        weighted_right_error += effective_weight * abs(
            ion.observed_relative_intensity - ion.right_predicted_intensity
        )
        denominator += effective_weight
    if denominator <= 0:
        left_fit = right_fit = 0.0
    else:
        left_fit = max(0.0, min(1.0, 1.0 - weighted_left_error / denominator))
        right_fit = max(0.0, min(1.0, 1.0 - weighted_right_error / denominator))
    total_weight = sum(ion.diagnostic_weight for ion in ions)
    raw_js = _jensen_shannon(left_mean, right_mean)
    reliability_factor = (
        max(0.0, 1.0 - float(np.mean(0.5 * (left_sd + right_sd))))
        if ensemble_qualified else 0.0
    )
    return (
        ions,
        left_fit,
        right_fit,
        total_weight,
        raw_js,
        reliability_factor,
        ensemble_audit,
    )


def _group_predictions(
    prediction_records: Sequence[Mapping[str, Any]],
) -> tuple[
    dict[tuple[int, str, float, str], list[dict[str, Any]]],
    dict[tuple[int, str, float, str], dict[str, Any]],
]:
    grouped: dict[tuple[int, str, float, str], list[dict[str, Any]]] = (
        defaultdict(list)
    )
    metadata: dict[tuple[int, str, float, str], dict[str, Any]] = {}
    for record in prediction_records:
        if not isinstance(record, Mapping):
            continue
        row_index = record.get("row_index")
        identity = str(record.get("identity_key") or "").strip()
        energy = _number(record.get("collision_energy"))
        unit = _unit(record.get("collision_energy_unit"))
        try:
            row_index = int(row_index)
        except (TypeError, ValueError):
            continue
        if row_index < 0 or not identity or energy is None or energy <= 0:
            continue
        members = _prediction_members(record)
        if not members:
            continue
        key = (row_index, identity, round(energy, 6), unit)
        grouped[key].extend(members)
        aggregate = metadata.setdefault(key, {"records": []})
        aggregate["records"].append(dict(record))
        aggregate["energy_transfer_qualified"] = all(
            bool(item.get("energy_transfer_qualified", False))
            for item in grouped[key]
        )
        aggregate["architecture_families"] = sorted({
            str(item.get("architecture_family") or "")
            for item in grouped[key]
            if str(item.get("architecture_family") or "")
        })
    return dict(grouped), metadata


def _common_energy_keys(
    row_index: int,
    left_identity: str,
    right_identity: str,
    grouped: Mapping[tuple[int, str, float, str], Any],
) -> list[tuple[float, str]]:
    left = {
        (energy, unit)
        for index, identity, energy, unit in grouped
        if index == row_index and identity == left_identity
    }
    right = {
        (energy, unit)
        for index, identity, energy, unit in grouped
        if index == row_index and identity == right_identity
    }
    return sorted(left & right, key=lambda item: (item[1], item[0]))


def _observed_peaks(row: Mapping[str, Any]) -> list[tuple[float, float]]:
    for key in ("raw_peaks", "query_peaks", "lcms_consensus_peaks"):
        peaks = _sanitize_peaks(row.get(key) or [])
        if peaks:
            return peaks
    return []


def _build_pair(
    *,
    row_index: int,
    row: Mapping[str, Any],
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    left_structure: Mapping[str, Any],
    right_structure: Mapping[str, Any],
    grouped: Mapping[tuple[int, str, float, str], list[dict[str, Any]]],
    policy: Mapping[str, Any],
) -> ContrastivePair:
    left_identity = str(left["_identity_key"])
    right_identity = str(right["_identity_key"])
    pair_key = _pair_key(left_identity, right_identity)
    eligible, scope_reasons, scope_abstention = _pair_scope(
        left, right, left_structure, right_structure, policy
    )
    common = _common_energy_keys(
        row_index, left_identity, right_identity, grouped
    )
    observed = _observed_peaks(row)
    observed_energy = _number(row.get("collision_energy"))
    observed_unit = _unit(row.get("collision_energy_unit"))
    uncertainty_method = str(
        (policy.get("prediction_uncertainty") or {}).get("method")
        or "independent_architecture_disagreement"
    )
    if not eligible:
        return ContrastivePair(
            pair_key=pair_key,
            left_identity_key=left_identity,
            right_identity_key=right_identity,
            left_compound_name=str(left.get("compound_name") or ""),
            right_compound_name=str(right.get("compound_name") or ""),
            eligible=False,
            status="ineligible",
            scope_reasons=scope_reasons,
            abstention_reason=scope_abstention,
            observed_collision_energy=observed_energy,
            collision_energy_unit=(observed_unit if observed_energy else None),
        )

    energy_separations: list[EnergyPairSeparation] = []
    observed_statistics: tuple[
        list[DiagnosticIon],
        float,
        float,
        float,
        float,
        float,
        dict[str, Any],
    ] | None = None
    energy_tolerance = float(
        (policy.get("peak_alignment") or {}).get("observed_energy_tolerance", 0.01)
    )
    for energy, unit in common:
        left_members = grouped[(row_index, left_identity, energy, unit)]
        right_members = grouped[(row_index, right_identity, energy, unit)]
        statistics = _diagnostic_statistics(
            left_members,
            right_members,
            observed if (
                observed_energy is not None
                and abs(energy - observed_energy) <= energy_tolerance
                and unit == observed_unit
            ) else (),
            policy,
        )
        ions, _, _, _, raw_js, reliability_factor, ensemble_audit = statistics
        adjusted_js = raw_js * reliability_factor
        energy_separations.append(EnergyPairSeparation(
            collision_energy=energy,
            collision_energy_unit=unit,
            raw_spectral_distribution_divergence=round(raw_js, 6),
            uncertainty_reliability_factor=round(reliability_factor, 6),
            spectral_distribution_divergence=round(adjusted_js, 6),
            diagnostic_ion_count=len(ions),
            diagnostic_ions=tuple(ion.mz for ion in ions),
            prediction_replicates_left=len(left_members),
            prediction_replicates_right=len(right_members),
            uncertainty_method=uncertainty_method,
            architecture_families_left=tuple(
                ensemble_audit["left_architecture_families"]
            ),
            architecture_families_right=tuple(
                ensemble_audit["right_architecture_families"]
            ),
            common_architecture_families=tuple(
                ensemble_audit["common_architecture_families"]
            ),
            common_architecture_family_count=int(
                ensemble_audit["common_architecture_family_count"]
            ),
            prediction_providers_left=tuple(ensemble_audit["left_providers"]),
            prediction_providers_right=tuple(ensemble_audit["right_providers"]),
            ensemble_qualified=bool(ensemble_audit["qualified"]),
            instrument_domain_qualified=bool(
                ensemble_audit["instrument_domain_qualified"]
            ),
            between_architecture_uncertainty=round(
                float(ensemble_audit["between_architecture_uncertainty"]), 6
            ),
            within_architecture_uncertainty=round(
                float(ensemble_audit["within_architecture_uncertainty"]), 6
            ),
        ))
        if (
            observed_energy is not None
            and abs(energy - observed_energy) <= energy_tolerance
            and unit == observed_unit
        ):
            observed_statistics = statistics

    base = dict(
        pair_key=pair_key,
        left_identity_key=left_identity,
        right_identity_key=right_identity,
        left_compound_name=str(left.get("compound_name") or ""),
        right_compound_name=str(right.get("compound_name") or ""),
        eligible=True,
        scope_reasons=scope_reasons,
        observed_collision_energy=observed_energy,
        collision_energy_unit=(observed_unit if observed_energy else None),
        energy_separations=tuple(energy_separations),
    )
    if not observed:
        return ContrastivePair(
            **base,
            status="abstained",
            abstention_reason="missing_observed_product_ions",
        )
    if observed_energy is None:
        return ContrastivePair(
            **base,
            status="abstained",
            abstention_reason="missing_observed_collision_energy",
        )
    if observed_statistics is None:
        return ContrastivePair(
            **base,
            status="abstained",
            abstention_reason="no_pair_prediction_at_observed_energy_and_unit",
        )
    (
        ions,
        left_fit,
        right_fit,
        total_weight,
        _,
        _,
        observed_ensemble_audit,
    ) = observed_statistics
    diagnostic = policy.get("diagnostic_evidence") or {}
    minimum_ions = int(diagnostic.get("minimum_diagnostic_ions", 2))
    minimum_weight = float(diagnostic.get("minimum_total_diagnostic_weight", 0.08))
    minimum_margin = float(
        diagnostic.get("minimum_absolute_pair_fit_margin", 0.05)
    )
    margin = left_fit - right_fit
    abstention_reason = None
    if not observed_ensemble_audit["qualified"]:
        abstention_reason = "insufficient_independent_architecture_families"
    elif len(ions) < minimum_ions:
        abstention_reason = "insufficient_candidate_diagnostic_ions"
    elif total_weight < minimum_weight:
        abstention_reason = "insufficient_uncertainty_adjusted_diagnostic_weight"
    elif abs(margin) < minimum_margin:
        abstention_reason = "diagnostic_fit_margin_below_shadow_threshold"
    winner = None
    status = "abstained"
    if abstention_reason is None:
        winner = left_identity if margin > 0 else right_identity
        status = "resolved_counterfactual"
    return ContrastivePair(
        **base,
        status=status,
        abstention_reason=abstention_reason,
        diagnostic_ions=tuple(ions),
        diagnostic_ion_count=len(ions),
        total_diagnostic_weight=round(total_weight, 6),
        left_diagnostic_fit=round(left_fit, 6),
        right_diagnostic_fit=round(right_fit, 6),
        signed_fit_margin_left_minus_right=round(margin, 6),
        counterfactual_winner_identity_key=winner,
        ensemble_qualified=bool(observed_ensemble_audit["qualified"]),
        common_architecture_family_count=int(
            observed_ensemble_audit["common_architecture_family_count"]
        ),
        common_architecture_families=tuple(
            observed_ensemble_audit["common_architecture_families"]
        ),
    )


def _candidate_support(
    candidates: Sequence[Mapping[str, Any]],
    pairs: Sequence[ContrastivePair],
) -> tuple[tuple[CandidateContrastSupport, ...], str | None, str]:
    support: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        identity = str(candidate["_identity_key"])
        support[identity] = {
            "identity_key": identity,
            "compound_name": str(candidate.get("compound_name") or ""),
            "production_top1": bool(candidate.get("_production_top1")),
            "eligible_pairs": 0,
            "resolved_pairs": 0,
            "pair_wins": 0,
            "pair_losses": 0,
            "pair_abstentions": 0,
            "net_diagnostic_support": 0.0,
            "counterfactual_rank": None,
        }
    for pair in pairs:
        if not pair.eligible:
            continue
        for identity in (pair.left_identity_key, pair.right_identity_key):
            support[identity]["eligible_pairs"] += 1
        if pair.status != "resolved_counterfactual":
            support[pair.left_identity_key]["pair_abstentions"] += 1
            support[pair.right_identity_key]["pair_abstentions"] += 1
            continue
        margin = float(pair.signed_fit_margin_left_minus_right or 0.0)
        support[pair.left_identity_key]["resolved_pairs"] += 1
        support[pair.right_identity_key]["resolved_pairs"] += 1
        support[pair.left_identity_key]["net_diagnostic_support"] += margin
        support[pair.right_identity_key]["net_diagnostic_support"] -= margin
        winner = str(pair.counterfactual_winner_identity_key or "")
        loser = (
            pair.right_identity_key
            if winner == pair.left_identity_key
            else pair.left_identity_key
        )
        support[winner]["pair_wins"] += 1
        support[loser]["pair_losses"] += 1

    eligible_pairs = [pair for pair in pairs if pair.eligible]
    resolved_pairs = [
        pair for pair in eligible_pairs if pair.status == "resolved_counterfactual"
    ]
    participating = [
        item for item in support.values() if item["eligible_pairs"] > 0
    ]
    complete = bool(
        eligible_pairs
        and len(resolved_pairs) == len(eligible_pairs)
        and participating
        and all(item["resolved_pairs"] > 0 for item in participating)
    )
    winner: str | None = None
    status: str
    if complete:
        ordered = sorted(
            participating,
            key=lambda item: (
                -int(item["pair_wins"]),
                -float(item["net_diagnostic_support"]),
                str(item["identity_key"]),
            ),
        )
        unique_lead = len(ordered) == 1 or (
            ordered[0]["pair_wins"], round(ordered[0]["net_diagnostic_support"], 9)
        ) != (
            ordered[1]["pair_wins"], round(ordered[1]["net_diagnostic_support"], 9)
        )
        if unique_lead:
            for rank, item in enumerate(ordered, start=1):
                item["counterfactual_rank"] = rank
            winner = str(ordered[0]["identity_key"])
            status = "resolved_counterfactual"
        else:
            status = "abstained_tied_pair_support"
    elif resolved_pairs:
        status = "partial_pair_support_no_counterfactual_top1"
    elif eligible_pairs:
        status = "abstained_no_resolved_candidate_pair"
    else:
        status = "not_applicable_no_eligible_candidate_pair"
    models = tuple(
        CandidateContrastSupport.model_validate({
            **item,
            "net_diagnostic_support": round(
                float(item["net_diagnostic_support"]), 6
            ),
        })
        for item in sorted(
            support.values(), key=lambda item: str(item["identity_key"])
        )
    )
    return models, winner, status


def _energy_selection(
    pairs: Sequence[ContrastivePair],
    prediction_metadata: Mapping[tuple[int, str, float, str], Mapping[str, Any]],
    row_index: int,
    policy: Mapping[str, Any],
) -> EnergySelection:
    eligible = [pair for pair in pairs if pair.eligible]
    energy_map: dict[tuple[float, str], dict[str, EnergyPairSeparation]] = defaultdict(dict)
    for pair in eligible:
        for separation in pair.energy_separations:
            energy_map[(
                separation.collision_energy,
                separation.collision_energy_unit,
            )][pair.pair_key] = separation
    per_energy: list[dict[str, Any]] = []
    for (energy, unit), by_pair in sorted(
        energy_map.items(), key=lambda item: (item[0][1], item[0][0])
    ):
        values = [
            float(item.spectral_distribution_divergence)
            for item in by_pair.values()
        ]
        diagnostic_ions = sorted({
            mz for item in by_pair.values() for mz in item.diagnostic_ions
        })[:20]
        per_energy.append({
            "collision_energy": energy,
            "collision_energy_unit": unit,
            "covered_pair_count": len(by_pair),
            "eligible_pair_count": len(eligible),
            "complete_pair_coverage": len(by_pair) == len(eligible),
            "independent_architecture_ensemble_qualified": bool(by_pair) and all(
                item.ensemble_qualified for item in by_pair.values()
            ),
            "instrument_domain_qualified": bool(by_pair) and all(
                item.instrument_domain_qualified for item in by_pair.values()
            ),
            "minimum_pair_separation": round(min(values), 6) if values else None,
            "mean_pair_separation": round(float(np.mean(values)), 6) if values else None,
            "diagnostic_ions": diagnostic_ions,
            "pairwise_separation": {
                pair_key: round(item.spectral_distribution_divergence, 6)
                for pair_key, item in sorted(by_pair.items())
            },
        })
    complete = [item for item in per_energy if item["complete_pair_coverage"]]
    energy_grid = tuple(sorted({
        float(item["collision_energy"]) for item in per_energy
    }))
    if not eligible:
        return EnergySelection(
            status="not_applicable_no_eligible_candidate_pair",
            eligible_pair_count=0,
            covered_pair_count=0,
            energy_grid=energy_grid,
            per_energy=tuple(per_energy),
        )
    if not complete:
        return EnergySelection(
            status="unavailable_incomplete_pair_prediction_coverage",
            eligible_pair_count=len(eligible),
            covered_pair_count=max(
                (int(item["covered_pair_count"]) for item in per_energy),
                default=0,
            ),
            energy_grid=energy_grid,
            per_energy=tuple(per_energy),
        )
    qualified_complete = [
        item for item in complete
        if item["independent_architecture_ensemble_qualified"]
    ]
    if not qualified_complete:
        return EnergySelection(
            status="unavailable_insufficient_independent_architecture_families",
            eligible_pair_count=len(eligible),
            covered_pair_count=len(eligible),
            energy_grid=energy_grid,
            per_energy=tuple(per_energy),
        )
    units = {str(item["collision_energy_unit"]) for item in qualified_complete}
    if len(units) != 1 or next(iter(units)).endswith("unspecified"):
        return EnergySelection(
            status="unavailable_unqualified_collision_energy_unit",
            eligible_pair_count=len(eligible),
            covered_pair_count=len(eligible),
            energy_grid=energy_grid,
            per_energy=tuple(per_energy),
        )
    chosen = sorted(
        qualified_complete,
        key=lambda item: (
            -float(item["minimum_pair_separation"] or 0.0),
            -float(item["mean_pair_separation"] or 0.0),
            float(item["collision_energy"]),
        ),
    )[0]
    if not chosen["instrument_domain_qualified"]:
        return EnergySelection(
            status="counterfactual_only_model_instrument_domain_shift",
            eligible_pair_count=len(eligible),
            covered_pair_count=len(eligible),
            energy_grid=energy_grid,
            per_energy=tuple(per_energy),
        )
    minimum_recommendable = float(
        (policy.get("energy_selection") or {}).get(
            "minimum_recommendable_worst_pair_separation", 0.10
        )
    )
    if float(chosen["minimum_pair_separation"] or 0.0) < minimum_recommendable:
        return EnergySelection(
            status=(
                "abstained_predicted_worst_pair_separation_below_threshold"
            ),
            eligible_pair_count=len(eligible),
            covered_pair_count=len(eligible),
            energy_grid=energy_grid,
            per_energy=tuple(per_energy),
        )
    return EnergySelection(
        status="recommended_shadow",
        recommended_collision_energy=float(chosen["collision_energy"]),
        collision_energy_unit=str(chosen["collision_energy_unit"]),
        minimum_pair_separation=float(chosen["minimum_pair_separation"]),
        eligible_pair_count=len(eligible),
        covered_pair_count=len(eligible),
        energy_grid=energy_grid,
        per_energy=tuple(per_energy),
    )


def build_evidiscrim_shadow(
    results: Sequence[Mapping[str, Any]],
    *,
    task_id: str,
    prediction_records: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
    model: Mapping[str, Any] | None = None,
    runtime: Mapping[str, Any] | None = None,
    skipped: Sequence[Mapping[str, Any]] = (),
    max_candidates_per_spectrum: int | None = None,
    source_results_sha256_override: str | None = None,
    source_results_immutability_verified: bool = False,
) -> dict[str, Any]:
    """Build a deterministic sidecar from predicted spectra and legacy rows."""

    if not task_id:
        raise ValueError("Task ID is required")
    if policy.get("module_state") != "shadow":
        raise ValueError("EviDiscrim policy must remain shadow")
    if source_results_sha256_override is not None:
        override = str(source_results_sha256_override).lower()
        if not re.fullmatch(r"[0-9a-f]{64}", override):
            raise ValueError("Source result digest override must be SHA-256 hex")
        if not source_results_immutability_verified:
            raise ValueError(
                "A source digest override requires an explicit immutability guard"
            )
        before = override
        source_scope = "bounded_rows_with_full_source_canonical_digest"
    else:
        before = results_sha256(results)
        source_scope = "provided_result_sequence"
    grouped, prediction_metadata = _group_predictions(prediction_records)
    policy_candidate_limit = int(
        (policy.get("candidate_scope") or {}).get(
            "maximum_candidates_per_spectrum", 5
        )
    )
    candidate_limit = (
        max(2, min(policy_candidate_limit, int(max_candidates_per_spectrum)))
        if max_candidates_per_spectrum is not None
        else policy_candidate_limit
    )
    spectra: list[EviDiscrimSpectrum] = []
    local_skipped = [dict(item) for item in skipped if isinstance(item, Mapping)]
    for local_row_index, row in enumerate(results):
        if not isinstance(row, Mapping):
            continue
        try:
            row_index = int(
                row.get("_evidiscrim_source_result_index", local_row_index)
            )
        except (TypeError, ValueError):
            row_index = local_row_index
        match = row.get("match")
        if not isinstance(match, Mapping):
            continue
        candidates = _candidate_rows(match, candidate_limit)
        if len(candidates) < 2:
            continue
        row_prediction_count = sum(
            index == row_index for index, _, _, _ in grouped
        )
        if not row_prediction_count:
            local_skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "no_candidate_forward_predictions",
            })
            continue
        structures = {
            str(candidate["_identity_key"]): _structure_profile(candidate)
            for candidate in candidates
        }
        sorted_candidates = sorted(
            candidates, key=lambda candidate: str(candidate["_identity_key"])
        )
        pairs: list[ContrastivePair] = []
        for left_index, left in enumerate(sorted_candidates):
            for right in sorted_candidates[left_index + 1:]:
                pairs.append(_build_pair(
                    row_index=row_index,
                    row=row,
                    left=left,
                    right=right,
                    left_structure=structures[str(left["_identity_key"])],
                    right_structure=structures[str(right["_identity_key"])],
                    grouped=grouped,
                    policy=policy,
                ))
        candidate_support, counterfactual_top1, spectrum_status = _candidate_support(
            candidates, pairs
        )
        production_top1 = str(candidates[0]["_identity_key"])
        eligible_count = sum(pair.eligible for pair in pairs)
        resolved_count = sum(
            pair.status == "resolved_counterfactual" for pair in pairs
        )
        abstained_count = sum(
            pair.eligible and pair.status == "abstained" for pair in pairs
        )
        ineligible_count = sum(not pair.eligible for pair in pairs)
        if spectrum_status == "resolved_counterfactual":
            spectrum_abstention = None
        else:
            spectrum_abstention = spectrum_status
        spectra.append(EviDiscrimSpectrum(
            row_index=row_index,
            source_row_sha256=source_row_sha256(
                row, max_candidates_per_spectrum=candidate_limit
            ),
            spectrum_id=str(row.get("spectrum_id") or row_index),
            query_precursor_mz=_number(row.get("precursor_mz")),
            observed_collision_energy=_number(row.get("collision_energy")),
            collision_energy_unit=(
                _unit(row.get("collision_energy_unit"))
                if _number(row.get("collision_energy")) is not None else None
            ),
            production_top1_identity_key=production_top1,
            candidate_count=len(candidates),
            eligible_pair_count=eligible_count,
            resolved_pair_count=resolved_count,
            abstained_pair_count=abstained_count,
            ineligible_pair_count=ineligible_count,
            pairwise_comparisons=tuple(pairs),
            candidate_support=candidate_support,
            counterfactual_top1_identity_key=counterfactual_top1,
            counterfactual_candidate_order_changed=bool(
                counterfactual_top1 and counterfactual_top1 != production_top1
            ),
            spectrum_status=spectrum_status,
            spectrum_abstention_reason=spectrum_abstention,
            energy_selection=_energy_selection(
                pairs, prediction_metadata, row_index, policy
            ),
        ))
    after = (
        before
        if source_results_sha256_override is not None
        else results_sha256(results)
    )
    all_pairs = [
        pair for spectrum in spectra for pair in spectrum.pairwise_comparisons
    ]
    summary = {
        "prediction_records_received": len(prediction_records),
        "maximum_candidates_per_spectrum": candidate_limit,
        "source_result_scope": source_scope,
        "spectra_evaluated": len(spectra),
        "candidate_pairs_total": len(all_pairs),
        "eligible_connectivity_pairs": sum(pair.eligible for pair in all_pairs),
        "resolved_counterfactual_pairs": sum(
            pair.status == "resolved_counterfactual" for pair in all_pairs
        ),
        "abstained_eligible_pairs": sum(
            pair.eligible and pair.status == "abstained" for pair in all_pairs
        ),
        "ensemble_qualified_pairs": sum(
            pair.eligible and pair.ensemble_qualified for pair in all_pairs
        ),
        "insufficient_architecture_pairs": sum(
            pair.abstention_reason
            == "insufficient_independent_architecture_families"
            for pair in all_pairs
        ),
        "ineligible_pairs": sum(not pair.eligible for pair in all_pairs),
        "same_connectivity_out_of_scope_pairs": sum(
            pair.abstention_reason
            in {
                "same_2d_connectivity_stereoisomer_out_of_scope",
                "duplicate_same_2d_structure",
            }
            for pair in all_pairs
        ),
        "counterfactual_candidate_order_changed_spectra": sum(
            spectrum.counterfactual_candidate_order_changed for spectrum in spectra
        ),
        "collision_energy_recommendations": sum(
            spectrum.energy_selection.status == "recommended_shadow"
            for spectrum in spectra
        ),
        "instrument_domain_shifted_energy_audits": sum(
            spectrum.energy_selection.status
            == "counterfactual_only_model_instrument_domain_shift"
            for spectrum in spectra
        ),
        "skipped_records": len(local_skipped),
        "production_ranking_changed": False,
        "confidence_changed": False,
        "high_promotions": 0,
        "validation_status": "unvalidated_shadow_no_accuracy_gain_claim",
        "uncertainty_gate": (
            "minimum_two_common_independent_architecture_families_no_singleton_floor"
        ),
    }
    status = "completed_shadow" if spectra else "not_applicable"
    artifact = EviDiscrimArtifact(
        task_id=task_id,
        policy_id=str(policy.get("policy_id") or "evidiscrim-unversioned"),
        status=status,
        source_results_sha256=before,
        source_results_sha256_after=after,
        source_results_mutated=before != after,
        model=dict(model or {
            "id": "unspecified_forward_spectrum_provider",
            "role": "candidate_structure_collision_energy_to_predicted_msms",
            "score_semantics": "candidate_conditioned_prediction_not_probability",
        }),
        runtime=dict(runtime or {}),
        spectra=tuple(spectra),
        skipped=tuple(local_skipped),
        summary=summary,
        interpretation=(
            "EviDiscrim-MS reports independent-architecture-disagreement-adjusted "
            "candidate-diagnostic "
            "counterfactuals and planning-only collision-energy separation. "
            "A singleton model cannot resolve a pair or recommend collision energy. "
            "Predicted spectra are candidate-conditioned simulations, not "
            "independent confirmation, and cannot alter production ranking, "
            "confidence, claim scope, or release."
        ),
    )
    return artifact.model_dump(mode="json")


__all__ = [
    "EVIDISCRIM_MODULE_ID",
    "build_evidiscrim_shadow",
    "evidiscrim_artifact_summary",
    "evidiscrim_artifact_path",
    "load_evidiscrim_policy",
    "read_evidiscrim_artifact",
    "results_sha256",
    "source_row_sha256",
    "stream_select_evidiscrim_rows",
]
