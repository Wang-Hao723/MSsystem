"""Conservative multi-energy support for coeluting authentic isomers.

This module answers a deliberately narrower question than a chromatographic
peak assignment: can one observed MS/MS entity be reconstructed by *both* of
two authentic-standard spectra at multiple collision energies?  A positive
answer supports joint presence of the pair, but never claims two separated
peaks, an exact Top-1 identity, or formal Level 1 identification.

The inference API accepts no truth labels or expected compound names.  Every
release gate is serialized so that a failed decision remains auditable.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import asdict, dataclass
from itertools import combinations
from typing import Any, Iterable, Mapping, Sequence

from matching.authentic_anchor import (
    _connectivity,
    greedy_cosine_similarity,
    validate_anchor_library,
)
from matching.chimeric_deconvolution import deconvolve_chimeric_spectrum


EVITWIN_COELUTION_SCHEMA_VERSION = "evitwin-coelution/1.0"
EVITWIN_COELUTION_ALGORITHM_ID = (
    "mssystem-v12-evitwin-multice-coelution-gate-1.0.0"
)


@dataclass(frozen=True)
class EviTwinCoelutionPolicy:
    """Frozen gates for dual-authentic-reference coelution support."""

    fragment_mz_tolerance_da: float = 0.02
    fragment_mz_tolerance_ppm: float = 20.0
    collision_energy_tolerance_ev: float = 1.0
    minimum_collision_energy_count: int = 2
    minimum_component_fraction: float = 0.10
    maximum_component_fraction_span: float = 0.20
    minimum_reconstruction_cosine: float = 0.90
    minimum_explained_energy: float = 0.75
    minimum_reconstruction_gain_over_single: float = 0.02
    maximum_component_similarity: float = 0.97
    minimum_specific_peaks_per_component_per_energy: int = 1
    minimum_aggregate_specific_peaks_per_component: int = 2
    minimum_reference_relative_intensity: float = 0.02
    l1_penalty: float = 0.02
    maximum_candidate_references: int = 10
    allow_exact_top1_change: bool = False
    allow_formal_level1_release: bool = False

    def validated(self) -> "EviTwinCoelutionPolicy":
        if not 0 < float(self.fragment_mz_tolerance_da) <= 0.2:
            raise ValueError("fragment_mz_tolerance_da must be in (0, 0.2]")
        if not 0 <= float(self.fragment_mz_tolerance_ppm) <= 100:
            raise ValueError("fragment_mz_tolerance_ppm must be in [0, 100]")
        if not 0 <= float(self.collision_energy_tolerance_ev) <= 10:
            raise ValueError("collision_energy_tolerance_ev must be in [0, 10]")
        if int(self.minimum_collision_energy_count) < 2:
            raise ValueError("at least two collision-energy conditions are required")
        if not 0 < float(self.minimum_component_fraction) < 0.5:
            raise ValueError("minimum_component_fraction must be in (0, 0.5)")
        if not 0 <= float(self.maximum_component_fraction_span) <= 1:
            raise ValueError("maximum_component_fraction_span must be in [0, 1]")
        for field_name in (
            "minimum_reconstruction_cosine",
            "minimum_explained_energy",
            "minimum_reconstruction_gain_over_single",
            "maximum_component_similarity",
            "minimum_reference_relative_intensity",
        ):
            value = float(getattr(self, field_name))
            if not 0 <= value <= 1:
                raise ValueError(f"{field_name} must be in [0, 1]")
        if int(self.minimum_specific_peaks_per_component_per_energy) < 1:
            raise ValueError("per-energy specific-peak minimum must be positive")
        if int(self.minimum_aggregate_specific_peaks_per_component) < 1:
            raise ValueError("aggregate specific-peak minimum must be positive")
        if not 0 < float(self.l1_penalty) <= 1:
            raise ValueError("l1_penalty must be in (0, 1]")
        if not 2 <= int(self.maximum_candidate_references) <= 50:
            raise ValueError("maximum_candidate_references must be in [2, 50]")
        if self.allow_exact_top1_change:
            raise ValueError("coelution support cannot change exact Top-1 identity")
        if self.allow_formal_level1_release:
            raise ValueError("coelution support cannot release formal Level 1")
        return self


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _collision_energy(row: Mapping[str, Any]) -> float | None:
    value = row.get("collision_energy_ev")
    if value is None:
        value = row.get("collision_energy")
    return _finite(value)


def _peaks(row: Mapping[str, Any]) -> Sequence[Sequence[Any]]:
    value = row.get("query_peaks") or row.get("peaks") or ()
    return value if isinstance(value, Sequence) else ()


def _normalised_formula(value: Any) -> str:
    return "".join(str(value or "").split()).upper()


def _reference_specific_masses(
    reference_peaks: Iterable[Sequence[Any]],
    competing_peaks: Iterable[Sequence[Any]],
    *,
    tolerance_da: float,
    minimum_relative_intensity: float,
) -> list[float]:
    prepared: list[tuple[float, float]] = []
    for peak in reference_peaks or ():
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            prepared.append((mz, intensity))
    other_masses = []
    for peak in competing_peaks or ():
        try:
            mz = float(peak[0])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(mz) and mz > 0:
            other_masses.append(mz)
    base = max((intensity for _, intensity in prepared), default=0.0)
    if base <= 0:
        return []
    return sorted(
        mz
        for mz, intensity in prepared
        if intensity / base >= minimum_relative_intensity
        and not any(abs(mz - other) <= tolerance_da for other in other_masses)
    )


def _matched_specific_masses(
    query_peaks: Iterable[Sequence[Any]],
    specific_masses: Iterable[float],
    *,
    tolerance_da: float,
) -> list[float]:
    query_masses = []
    for peak in query_peaks or ():
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if (
            math.isfinite(mz)
            and math.isfinite(intensity)
            and mz > 0
            and intensity > 0
        ):
            query_masses.append(mz)
    return [
        float(mass)
        for mass in specific_masses
        if any(abs(float(mass) - observed) <= tolerance_da for observed in query_masses)
    ]


def _best_reference_spectrum(
    entry: Mapping[str, Any],
    query_peaks: Sequence[Sequence[Any]],
    collision_energy: float,
    policy: EviTwinCoelutionPolicy,
) -> Mapping[str, Any] | None:
    eligible = []
    for spectrum in entry.get("spectra") or ():
        if not isinstance(spectrum, Mapping):
            continue
        reference_energy = _collision_energy(spectrum)
        reference_peaks = spectrum.get("peaks") or ()
        if (
            reference_energy is None
            or abs(reference_energy - collision_energy)
            > policy.collision_energy_tolerance_ev
            or len(reference_peaks) < 3
        ):
            continue
        similarity = greedy_cosine_similarity(
            query_peaks,
            reference_peaks,
            tolerance_da=policy.fragment_mz_tolerance_da,
        )
        eligible.append((similarity, -abs(reference_energy - collision_energy), spectrum))
    if not eligible:
        return None
    eligible.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return eligible[0][2]


def _component_fraction(result: Mapping[str, Any], reference_id: str) -> float:
    for row in result.get("components") or ():
        if str(row.get("id") or "") == reference_id:
            return float(row.get("fraction") or 0.0)
    return 0.0


def _evaluate_energy(
    query: Mapping[str, Any],
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    policy: EviTwinCoelutionPolicy,
) -> dict[str, Any] | None:
    energy = _collision_energy(query)
    query_peaks = _peaks(query)
    if energy is None or len(query_peaks) < 3:
        return None
    left_spectrum = _best_reference_spectrum(left, query_peaks, energy, policy)
    right_spectrum = _best_reference_spectrum(right, query_peaks, energy, policy)
    if left_spectrum is None or right_spectrum is None:
        return None
    left_id = str(left.get("reference_id") or "")
    right_id = str(right.get("reference_id") or "")
    left_peaks = left_spectrum.get("peaks") or ()
    right_peaks = right_spectrum.get("peaks") or ()
    components = [
        {"id": left_id, "peaks": left_peaks},
        {"id": right_id, "peaks": right_peaks},
    ]
    common_kwargs = {
        "penalty": policy.l1_penalty,
        "tolerance_da": policy.fragment_mz_tolerance_da,
        "tolerance_ppm": policy.fragment_mz_tolerance_ppm,
        "min_reconstruction_cosine": policy.minimum_reconstruction_cosine,
        "min_explained_energy": policy.minimum_explained_energy,
        "min_target_fraction": policy.minimum_component_fraction,
        "max_target_fraction_span": policy.maximum_component_fraction_span,
        "max_component_similarity": policy.maximum_component_similarity,
    }
    left_fit = deconvolve_chimeric_spectrum(
        query_peaks, components, left_id, **common_kwargs
    )
    right_fit = deconvolve_chimeric_spectrum(
        query_peaks, components, right_id, **common_kwargs
    )
    left_single = greedy_cosine_similarity(
        query_peaks,
        left_peaks,
        tolerance_da=policy.fragment_mz_tolerance_da,
    )
    right_single = greedy_cosine_similarity(
        query_peaks,
        right_peaks,
        tolerance_da=policy.fragment_mz_tolerance_da,
    )
    reconstruction = float(left_fit.get("reconstruction_cosine") or 0.0)
    gain = reconstruction - max(left_single, right_single)
    left_specific = _matched_specific_masses(
        query_peaks,
        _reference_specific_masses(
            left_peaks,
            right_peaks,
            tolerance_da=policy.fragment_mz_tolerance_da,
            minimum_relative_intensity=policy.minimum_reference_relative_intensity,
        ),
        tolerance_da=policy.fragment_mz_tolerance_da,
    )
    right_specific = _matched_specific_masses(
        query_peaks,
        _reference_specific_masses(
            right_peaks,
            left_peaks,
            tolerance_da=policy.fragment_mz_tolerance_da,
            minimum_relative_intensity=policy.minimum_reference_relative_intensity,
        ),
        tolerance_da=policy.fragment_mz_tolerance_da,
    )
    gates = {
        "left_sparse_fit": bool(left_fit.get("accepted")),
        "right_sparse_fit": bool(right_fit.get("accepted")),
        "reconstruction_gain": bool(
            gain >= policy.minimum_reconstruction_gain_over_single
        ),
        "left_specific_fragment": bool(
            len(left_specific)
            >= policy.minimum_specific_peaks_per_component_per_energy
        ),
        "right_specific_fragment": bool(
            len(right_specific)
            >= policy.minimum_specific_peaks_per_component_per_energy
        ),
    }
    return {
        "collision_energy_ev": round(energy, 3),
        "status": "supported" if all(gates.values()) else "rejected",
        "support_gates": gates,
        "left_single_similarity": round(left_single, 6),
        "right_single_similarity": round(right_single, 6),
        "reconstruction_cosine": round(reconstruction, 6),
        "explained_energy": left_fit.get("explained_energy"),
        "reconstruction_gain_over_best_single": round(gain, 6),
        "left_component_fraction": round(_component_fraction(left_fit, left_id), 6),
        "right_component_fraction": round(_component_fraction(left_fit, right_id), 6),
        "left_specific_evidence_mz": [round(value, 6) for value in left_specific],
        "right_specific_evidence_mz": [round(value, 6) for value in right_specific],
        "left_fit_status": left_fit.get("status"),
        "right_fit_status": right_fit.get("status"),
        "component_similarity": left_fit.get("target_max_component_similarity"),
    }


def _pair_result(
    queries: Sequence[Mapping[str, Any]],
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    policy: EviTwinCoelutionPolicy,
) -> dict[str, Any]:
    left_id = str(left.get("reference_id") or "")
    right_id = str(right.get("reference_id") or "")
    rows = [
        row
        for query in queries
        if (row := _evaluate_energy(query, left, right, policy)) is not None
    ]
    # Repeated DDA events at one CE are technical repeats, not independent
    # collision-energy evidence.  Retain only the strongest audit row per CE.
    best_by_energy: dict[float, dict[str, Any]] = {}
    for row in rows:
        energy = float(row["collision_energy_ev"])
        current = best_by_energy.get(energy)
        strength = (
            row["status"] == "supported",
            float(row["reconstruction_gain_over_best_single"]),
            float(row["reconstruction_cosine"]),
        )
        if current is None or strength > (
            current["status"] == "supported",
            float(current["reconstruction_gain_over_best_single"]),
            float(current["reconstruction_cosine"]),
        ):
            best_by_energy[energy] = row
    rows = [best_by_energy[key] for key in sorted(best_by_energy)]
    supported = [row for row in rows if row["status"] == "supported"]
    left_specific = {
        round(float(mass), 4)
        for row in supported
        for mass in row["left_specific_evidence_mz"]
    }
    right_specific = {
        round(float(mass), 4)
        for row in supported
        for mass in row["right_specific_evidence_mz"]
    }
    left_fractions = [float(row["left_component_fraction"]) for row in supported]
    right_fractions = [float(row["right_component_fraction"]) for row in supported]
    fraction_span = max(
        max(left_fractions, default=0.0) - min(left_fractions, default=0.0),
        max(right_fractions, default=0.0) - min(right_fractions, default=0.0),
    )
    aggregate_gates = {
        "multi_collision_energy_support": bool(
            len(supported) >= policy.minimum_collision_energy_count
        ),
        "left_aggregate_specific_fragments": bool(
            len(left_specific)
            >= policy.minimum_aggregate_specific_peaks_per_component
        ),
        "right_aggregate_specific_fragments": bool(
            len(right_specific)
            >= policy.minimum_aggregate_specific_peaks_per_component
        ),
        "cross_energy_fraction_stability": bool(
            fraction_span <= policy.maximum_component_fraction_span
        ),
    }
    status = "supported" if all(aggregate_gates.values()) else "abstained"
    return {
        "pair_id": "|".join(sorted((left_id, right_id))),
        "left_reference_id": left_id,
        "right_reference_id": right_id,
        "left_compound_code": left.get("compound_code"),
        "right_compound_code": right.get("compound_code"),
        "left_compound_name": left.get("compound_name"),
        "right_compound_name": right.get("compound_name"),
        "left_connectivity": _connectivity(left.get("inchikey")),
        "right_connectivity": _connectivity(right.get("inchikey")),
        "formula": left.get("formula"),
        "status": status,
        "supporting_collision_energy_count": len(supported),
        "supporting_collision_energies_ev": [
            row["collision_energy_ev"] for row in supported
        ],
        "aggregate_left_specific_peak_count": len(left_specific),
        "aggregate_right_specific_peak_count": len(right_specific),
        "maximum_component_fraction_span": round(fraction_span, 6),
        "median_reconstruction_cosine": round(
            statistics.median(
                [float(row["reconstruction_cosine"]) for row in supported]
            ),
            6,
        )
        if supported
        else 0.0,
        "median_reconstruction_gain_over_single": round(
            statistics.median(
                [
                    float(row["reconstruction_gain_over_best_single"])
                    for row in supported
                ]
            ),
            6,
        )
        if supported
        else 0.0,
        "aggregate_gates": aggregate_gates,
        "collision_energy_audit": rows,
    }


def evaluate_coeluting_reference_pairs(
    query_spectra: Sequence[Mapping[str, Any]],
    library: Mapping[str, Any],
    candidate_reference_ids: Sequence[str],
    *,
    policy: EviTwinCoelutionPolicy | None = None,
) -> dict[str, Any]:
    """Evaluate whether one entity contains a uniquely supported isomer pair."""

    validate_anchor_library(library)
    selected = (policy or EviTwinCoelutionPolicy()).validated()
    entries = {
        str(entry.get("reference_id") or ""): entry
        for entry in library.get("entries") or ()
        if isinstance(entry, Mapping)
    }
    candidate_ids = []
    for value in candidate_reference_ids:
        token = str(value or "")
        if token in entries and token not in candidate_ids:
            candidate_ids.append(token)
    truncated = len(candidate_ids) > selected.maximum_candidate_references
    candidate_ids = candidate_ids[: selected.maximum_candidate_references]
    pair_rows = []
    skipped_same_connectivity = 0
    skipped_formula = 0
    for left_id, right_id in combinations(candidate_ids, 2):
        left, right = entries[left_id], entries[right_id]
        left_connectivity = _connectivity(left.get("inchikey"))
        right_connectivity = _connectivity(right.get("inchikey"))
        if left_connectivity and left_connectivity == right_connectivity:
            skipped_same_connectivity += 1
            continue
        left_formula = _normalised_formula(left.get("formula"))
        right_formula = _normalised_formula(right.get("formula"))
        if not left_formula or left_formula != right_formula:
            skipped_formula += 1
            continue
        pair_rows.append(_pair_result(query_spectra, left, right, selected))
    pair_rows.sort(
        key=lambda row: (
            row["status"] != "supported",
            -int(row["supporting_collision_energy_count"]),
            -float(row["median_reconstruction_gain_over_single"]),
            -float(row["median_reconstruction_cosine"]),
            str(row["pair_id"]),
        )
    )
    supported = [row for row in pair_rows if row["status"] == "supported"]
    if len(supported) == 1:
        decision = {
            "status": "qualified_coeluting_pair_support",
            "reportable_scope": "coeluting_pair_presence_support",
            "reason_codes": ["unique_multi_ce_dual_authentic_fit_passed"],
            "supported_pair": supported[0],
        }
    elif len(supported) > 1:
        decision = {
            "status": "abstained",
            "reportable_scope": "none",
            "reason_codes": ["multiple_coeluting_pairs_remain_ambiguous"],
            "supported_pair": None,
        }
    else:
        decision = {
            "status": "abstained",
            "reportable_scope": "none",
            "reason_codes": [
                "no_pair_passed_multi_ce_dual_component_and_specific_fragment_gates"
            ],
            "supported_pair": None,
        }
    return {
        "schema_version": EVITWIN_COELUTION_SCHEMA_VERSION,
        "algorithm_id": EVITWIN_COELUTION_ALGORITHM_ID,
        "policy": asdict(selected),
        "candidate_reference_count": len(candidate_ids),
        "candidate_list_truncated": truncated,
        "evaluated_pair_count": len(pair_rows),
        "supported_pair_count": len(supported),
        "skipped_same_connectivity_pair_count": skipped_same_connectivity,
        "skipped_formula_mismatch_pair_count": skipped_formula,
        "decision": decision,
        "pair_audit": pair_rows,
        "invariants": {
            "truth_labels_used": False,
            "expected_compound_names_used": False,
            "exact_top1_changed": False,
            "separated_feature_assignment_claimed": False,
            "formal_level1_released": False,
            "coefficients_are_quantitative_abundances": False,
        },
    }


def validate_coelution_result(payload: Mapping[str, Any]) -> None:
    if payload.get("schema_version") != EVITWIN_COELUTION_SCHEMA_VERSION:
        raise ValueError("unsupported EviTwin coelution schema")
    if payload.get("algorithm_id") != EVITWIN_COELUTION_ALGORITHM_ID:
        raise ValueError("unexpected EviTwin coelution algorithm id")
    policy = EviTwinCoelutionPolicy(**dict(payload.get("policy") or {})).validated()
    if asdict(policy) != dict(payload.get("policy") or {}):
        raise ValueError("coelution policy must be complete and canonical")
    invariants = payload.get("invariants") or {}
    required_false = (
        "truth_labels_used",
        "expected_compound_names_used",
        "exact_top1_changed",
        "separated_feature_assignment_claimed",
        "formal_level1_released",
        "coefficients_are_quantitative_abundances",
    )
    if any(invariants.get(key) is not False for key in required_false):
        raise ValueError("coelution safety invariants are incomplete")
