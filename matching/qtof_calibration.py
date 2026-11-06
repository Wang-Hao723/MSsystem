"""Identity-disjoint QTOF intensity calibration for forward spectra.

This adapter calibrates a base predictor to an exact QTOF acquisition domain.
It deliberately preserves the base architecture family and therefore never
counts as a second independent model in the EviDiscrim ensemble.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

import numpy as np
import sklearn
from sklearn.isotonic import IsotonicRegression

from matching.chemical_identity import rdkit_identity_profile
from matching.evidiscrim_spectrum import (
    aligned_predicted_observed_intensities,
    prepare_matchms_spectrum,
    processed_peak_pairs,
    spectrum_contract,
)


QTOF_CALIBRATION_SCHEMA_VERSION = "1.0.0"


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")).hexdigest()


def _text(value: Any) -> str:
    return " ".join(str(value or "").strip().lower().split())


def _energy(value: Any) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("collision energy must be finite and positive") from exc
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValueError("collision energy must be finite and positive")
    return round(parsed, 6)


def qtof_domain_key(record: Mapping[str, Any]) -> str:
    analyzer = _text(record.get("analyzer_type"))
    compact_analyzer = "".join(character for character in analyzer if character.isalnum())
    if not any(
        token in compact_analyzer for token in ("qtof", "quadrupoletimeofflight")
    ):
        raise ValueError("QTOF calibration requires an explicit QTOF analyzer")
    fields = [
        _text(record.get("instrument_vendor")),
        _text(record.get("instrument_model")),
        analyzer,
        _text(record.get("activation_type")),
        _text(record.get("polarity")),
        str(record.get("adduct") or "").strip(),
        f"{_energy(record.get('collision_energy')):g}",
        str(record.get("collision_energy_unit") or "").strip(),
    ]
    if any(not field for field in fields):
        raise ValueError("QTOF calibration domain fields must all be explicit")
    return "|".join(fields)


def fit_qtof_intensity_calibrator(
    examples: Sequence[Mapping[str, Any]],
    *,
    policy: Mapping[str, Any],
    blind_connectivity_keys: Sequence[str] = (),
) -> dict[str, Any]:
    """Fit exact-domain isotonic mappings on calibration-split standards only."""

    calibration_policy = policy.get("qtof_calibration") or {}
    processing = policy.get("spectrum_processing") or {}
    minimum_identities = int(
        calibration_policy.get("minimum_connectivity_identities_per_domain", 20)
    )
    minimum_peaks = int(
        calibration_policy.get("minimum_predicted_peaks_per_domain", 200)
    )
    blind = set(str(value) for value in blind_connectivity_keys)
    grouped: dict[str, dict[str, Any]] = defaultdict(
        lambda: {"x": [], "y": [], "identities": set(), "records": 0}
    )
    accepted_digest_rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    architecture_families: set[str] = set()
    providers: set[str] = set()
    for index, raw in enumerate(examples):
        if str(raw.get("split") or "") != "calibration":
            skipped.append({"index": index, "reason": "not_calibration_split"})
            continue
        profile = rdkit_identity_profile(
            raw.get("true_smiles"),
            declared_formula=raw.get("formula"),
            declared_inchikey=raw.get("true_inchikey"),
            require_single_component=True,
        )
        if not profile.get("valid"):
            skipped.append({
                "index": index,
                "reason": f"invalid_identity::{profile.get('reason')}",
            })
            continue
        connectivity = str(profile["connectivity_key"])
        if connectivity in blind:
            raise ValueError(
                "QTOF calibration identity overlaps the sealed blind test"
            )
        try:
            domain = qtof_domain_key(raw)
        except ValueError as exc:
            skipped.append({"index": index, "reason": str(exc)})
            continue
        observed = prepare_matchms_spectrum(
            raw.get("observed_peaks") or (),
            precursor_mz=float(raw.get("precursor_mz") or 0) or None,
            minimum_relative_intensity=float(
                processing.get("minimum_relative_intensity", 0.01)
            ),
            minimum_peak_count=int(processing.get("minimum_peak_count", 2)),
        )
        predicted = prepare_matchms_spectrum(
            raw.get("predicted_peaks") or (),
            precursor_mz=float(raw.get("precursor_mz") or 0) or None,
            minimum_relative_intensity=float(
                processing.get("minimum_relative_intensity", 0.01)
            ),
            minimum_peak_count=int(processing.get("minimum_peak_count", 2)),
        )
        if observed is None or predicted is None:
            skipped.append({
                "index": index,
                "reason": "insufficient_matchms_processed_peaks",
            })
            continue
        x, y = aligned_predicted_observed_intensities(
            predicted,
            observed,
            mz_tolerance_da=float(processing.get("mz_tolerance_da", 0.02)),
        )
        bucket = grouped[domain]
        bucket["x"].extend(float(value) for value in x)
        bucket["y"].extend(float(value) for value in y)
        bucket["identities"].add(connectivity)
        bucket["records"] += 1
        family = str(
            raw.get("base_architecture_family")
            or raw.get("architecture_family")
            or ""
        ).strip()
        provider = str(raw.get("provider_id") or raw.get("model_id") or "").strip()
        if not family or not provider:
            raise ValueError(
                "every calibration prediction must declare provider and base architecture"
            )
        architecture_families.add(family)
        providers.add(provider)
        accepted_digest_rows.append({
            "domain": domain,
            "connectivity_key": connectivity,
            "spectrum_id": str(raw.get("spectrum_id") or index),
            "provider_id": provider,
            "base_architecture_family": family,
            "observed_peaks": processed_peak_pairs(observed),
            "predicted_peaks": processed_peak_pairs(predicted),
        })
    if len(architecture_families) > 1:
        raise ValueError(
            "one calibration artifact may calibrate only one base architecture family"
        )
    curves: dict[str, Any] = {}
    rejected_domains: list[dict[str, Any]] = []
    for domain, bucket in sorted(grouped.items()):
        identity_count = len(bucket["identities"])
        peak_count = len(bucket["x"])
        if identity_count < minimum_identities or peak_count < minimum_peaks:
            rejected_domains.append({
                "domain_key": domain,
                "connectivity_identities": identity_count,
                "predicted_peaks": peak_count,
                "reason": "insufficient_domain_calibration_support",
            })
            continue
        x = np.asarray(bucket["x"], dtype=float)
        y = np.asarray(bucket["y"], dtype=float)
        model = IsotonicRegression(
            y_min=0.0,
            y_max=1.0,
            increasing=True,
            out_of_bounds="clip",
        ).fit(x, y)
        calibrated = model.predict(x)
        curves[domain] = {
            "x_thresholds": [round(float(value), 12) for value in model.X_thresholds_],
            "y_thresholds": [round(float(value), 12) for value in model.y_thresholds_],
            "connectivity_identities": identity_count,
            "records": int(bucket["records"]),
            "predicted_peaks": peak_count,
            "rmse_before": round(float(np.sqrt(np.mean((x - y) ** 2))), 10),
            "rmse_after": round(
                float(np.sqrt(np.mean((calibrated - y) ** 2))), 10
            ),
            "matched_peak_fraction": round(float(np.mean(y > 0)), 10),
        }
    base_family = next(iter(architecture_families), "")
    base_provider = next(iter(providers), "") if len(providers) == 1 else "multiple"
    training_digest = _canonical_sha256(accepted_digest_rows)
    calibrator_id = f"qtof-isotonic-{training_digest[:16]}"
    return {
        "schema_version": QTOF_CALIBRATION_SCHEMA_VERSION,
        "calibrator_id": calibrator_id,
        "status": "fitted" if curves else "unavailable_insufficient_domain_support",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "method": (
            "domain_specific_isotonic_predicted_to_observed_peak_intensity"
        ),
        "base_provider_id": base_provider,
        "base_architecture_family": base_family,
        "counts_as_independent_architecture": False,
        "architecture_family_after_calibration": base_family,
        "cross_instrument_transfer_allowed": False,
        "collision_energy_interpolation_allowed": False,
        "blind_identity_overlap_count": 0,
        "training_records_sha256": training_digest,
        "spectrum_contract": spectrum_contract(policy),
        "sklearn_version": sklearn.__version__,
        "curves": curves,
        "rejected_domains": rejected_domains,
        "skipped_records": skipped,
    }


def apply_qtof_intensity_calibrator(
    prediction_record: Mapping[str, Any],
    *,
    calibrator: Mapping[str, Any],
    instrument_context: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply one exact-domain curve without changing model-family identity."""

    if calibrator.get("status") != "fitted":
        raise ValueError("QTOF calibrator is not fitted")
    prediction_family = str(
        prediction_record.get("base_architecture_family")
        or prediction_record.get("architecture_family")
        or ""
    ).strip()
    calibrator_family = str(
        calibrator.get("base_architecture_family") or ""
    ).strip()
    if not prediction_family or prediction_family != calibrator_family:
        raise ValueError(
            "QTOF calibrator base architecture does not match the prediction"
        )
    domain = qtof_domain_key(instrument_context)
    curve = (calibrator.get("curves") or {}).get(domain)
    if not isinstance(curve, Mapping):
        raise ValueError("no exact QTOF calibration curve for this acquisition domain")
    x_thresholds = np.asarray(curve.get("x_thresholds") or (), dtype=float)
    y_thresholds = np.asarray(curve.get("y_thresholds") or (), dtype=float)
    if not len(x_thresholds) or len(x_thresholds) != len(y_thresholds):
        raise ValueError("invalid QTOF calibration curve")
    peaks = prediction_record.get("predicted_peaks") or ()
    mz_values: list[float] = []
    intensities: list[float] = []
    for raw in peaks:
        try:
            mz = float(raw[0])
            intensity = float(raw[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            mz_values.append(mz)
            intensities.append(intensity)
    if not intensities:
        raise ValueError("prediction has no usable peaks")
    maximum = max(intensities)
    normalized = np.asarray(intensities, dtype=float) / maximum
    calibrated = np.interp(
        normalized,
        x_thresholds,
        y_thresholds,
        left=float(y_thresholds[0]),
        right=float(y_thresholds[-1]),
    )
    calibrated_max = float(np.max(calibrated))
    if calibrated_max <= 0:
        raise ValueError("QTOF calibration suppressed every predicted peak")
    calibrated /= calibrated_max
    output = dict(prediction_record)
    output.update({
        "predicted_peaks": [
            [round(mz, 8), round(float(intensity), 10)]
            for mz, intensity in zip(mz_values, calibrated, strict=True)
            if intensity > 0
        ],
        "provider_id": (
            f"{prediction_record.get('provider_id') or 'unknown'}+"
            f"{calibrator['calibrator_id']}"
        ),
        "base_architecture_family": calibrator_family,
        "architecture_family": calibrator_family,
        "calibration_layer_id": calibrator["calibrator_id"],
        "qtof_calibration_domain_key": domain,
        "qtof_calibration_domain_qualified": True,
        "energy_transfer_qualified": True,
        "counts_as_independent_architecture": False,
        "calibration_residual_rmse": curve.get("rmse_after"),
    })
    return output


__all__ = [
    "QTOF_CALIBRATION_SCHEMA_VERSION",
    "apply_qtof_intensity_calibrator",
    "fit_qtof_intensity_calibrator",
    "qtof_domain_key",
]
