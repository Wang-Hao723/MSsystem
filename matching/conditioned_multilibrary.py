"""Condition-aware aggregation of experimental library references for V15.

The legacy ranker retained only the highest-scoring reference spectrum for one
2-D identity.  V15 keeps every reference record, collapses peak-identical
mirrors, and can combine the legacy maximum with a robust median across
condition-compatible provenance strata.  An all-condition median is retained
for review only and can never alter ranking.  Collision energy is compared only
inside an explicit unit domain; no eV/NCE/percent/voltage conversion is inferred.

The returned utilities are ranking scores, not probabilities.  Compound names,
ground truth and biological annotations are deliberately absent from every
calculation.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from utils.collision_energy import parse_collision_energy
from utils.instrument_metadata import normalize_instrument_family


SCHEMA_VERSION = "mssystem-v15-conditioned-multilibrary/1.1"


def _text(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _unit_interval(value: Any) -> float:
    number = _finite(value)
    if number is None:
        raise ValueError("each reference requires a finite spectral score")
    return max(0.0, min(1.0, number))


def normalize_source_family(value: Any) -> str:
    """Normalize dataset labels without treating accessions as sources."""

    raw = _text(value).casefold()
    if not raw:
        return ""
    normalized = re.sub(r"[^a-z0-9]+", "-", raw).strip("-")
    normalized = re.sub(r"-(?:pos|neg|positive|negative)$", "", normalized)
    normalized = re.sub(r"-(?:v)?\d+(?:-\d+)*$", "", normalized)
    families = (
        ("spectraverse", "spectraverse"),
        ("massbank", "massbank"),
        ("bmdms", "bmdms"),
        ("vaniya", "vaniya-fiehn"),
        ("mona", "mona"),
        ("gnps", "gnps"),
        ("nist", "nist"),
    )
    return next(
        (family for marker, family in families if marker in normalized),
        normalized,
    )


def reference_content_fingerprint(peaks: Any) -> str | None:
    """Hash base-peak-normalized product ions for mirror detection."""

    normalized: list[tuple[float, float]] = []
    if not isinstance(peaks, Sequence) or isinstance(peaks, (str, bytes)):
        return None
    for peak in peaks:
        if not isinstance(peak, Sequence) or isinstance(peak, (str, bytes)):
            continue
        if len(peak) < 2:
            continue
        mz = _finite(peak[0])
        intensity = _finite(peak[1])
        if mz is None or intensity is None or mz <= 0.0 or intensity <= 0.0:
            continue
        normalized.append((mz, intensity))
    if not normalized:
        return None
    base_peak = max(intensity for _, intensity in normalized)
    if base_peak <= 0.0:
        return None
    canonical = sorted(
        (round(mz, 5), round(intensity / base_peak, 6))
        for mz, intensity in normalized
    )
    return hashlib.sha256(
        json.dumps(canonical, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _energy(value: Any, unit: Any) -> dict[str, Any]:
    return parse_collision_energy(value, unit_hint=unit)


def _normalized_model(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", _text(value).casefold())


def _meaningful_model(value: Any) -> str:
    model = _normalized_model(value)
    generic = {
        "", "qtof", "lcqtof", "lcesiqtof", "orbitrap", "qe",
        "triplequadrupole", "iontrap", "massspectrometer",
    }
    return "" if model in generic or len(model) < 5 else model


def _vendor(value: Any) -> str:
    text = _normalized_model(value)
    markers = (
        (("thermo", "qexactive", "exploris", "orbitrap"), "thermo"),
        (("agilent",), "agilent"),
        (("sciex", "absciex", "tripletof"), "sciex"),
        (("bruker", "impactii", "maxisii"), "bruker"),
        (("waters", "synapt", "xevo"), "waters"),
        (("shimadzu",), "shimadzu"),
        (("jeol",), "jeol"),
    )
    return next(
        (vendor for aliases, vendor in markers if any(alias in text for alias in aliases)),
        "",
    )


def _mode(value: Any) -> str:
    normalized = _text(value).casefold()
    return {
        "single": "fixed",
        "fixed": "fixed",
        "stepped": "stepped",
        "step": "stepped",
        "ramp": "ramp",
        "ramped": "ramp",
        "rolling": "rolling",
        "dynamic": "dynamic",
        "ces": "CES",
        "bbcid": "bbCID",
        "mse": "MSE",
    }.get(normalized, "")


def _activation(value: Any, analyzer_family: str) -> str:
    compact = _normalized_model(value)
    if not compact:
        return ""
    if "higherenergy" in compact or compact == "hcd" or compact.endswith("hcd"):
        return "HCD"
    if "collisioninduceddissociation" in compact or compact in {
        "cid", "beamcid", "beamtypecid", "beamcollisioninduceddissociation",
    }:
        return "beam_cid" if analyzer_family == "qtof" else "CID"
    return compact


def reference_condition(
    query: Mapping[str, Any],
    reference: Mapping[str, Any],
    *,
    near_delta_by_unit: Mapping[str, Any],
) -> dict[str, Any]:
    """Describe acquisition compatibility without cross-domain CE conversion."""

    query_energy = _energy(
        query.get("collision_energy"), query.get("collision_energy_unit")
    )
    reference_energy = _energy(
        reference.get("lib_collision_energy", reference.get("collision_energy")),
        reference.get(
            "lib_collision_energy_unit", reference.get("collision_energy_unit")
        ),
    )
    query_family = normalize_instrument_family(
        query.get("instrument_type"), query.get("instrument")
    )
    reference_family = normalize_instrument_family(
        reference.get("lib_instrument_type", reference.get("instrument_type")),
        reference.get("lib_instrument", reference.get("instrument")),
    )
    query_model = _meaningful_model(query.get("instrument"))
    reference_model = _meaningful_model(
        reference.get("lib_instrument", reference.get("instrument"))
    )
    same_model = bool(
        query_model
        and reference_model
        and (
            query_model == reference_model
            or query_model in reference_model
            or reference_model in query_model
        )
    )
    same_family = bool(
        query_family and reference_family and query_family == reference_family
    )
    query_vendor = _vendor(query.get("instrument"))
    reference_vendor = _vendor(
        reference.get("lib_instrument", reference.get("instrument"))
    )
    same_vendor_family = bool(
        same_family
        and query_vendor
        and reference_vendor
        and query_vendor == reference_vendor
    )
    query_domain = _text(query.get("instrument_domain_id"))
    reference_domain = _text(
        reference.get("lib_instrument_domain_id", reference.get("instrument_domain_id"))
    )
    same_domain = bool(query_domain and reference_domain and query_domain == reference_domain)
    if query_domain and reference_domain:
        platform_anchor = same_domain
        platform_anchor_scope = "exact_instrument_domain"
    elif query_model and reference_model:
        platform_anchor = same_model
        platform_anchor_scope = "instrument_model"
    elif query_vendor and reference_vendor:
        platform_anchor = same_vendor_family
        platform_anchor_scope = "vendor_analyzer_family"
    else:
        # Historical artifacts sometimes retained only QTOF/Orbitrap.  This is
        # an analyzer-family anchor, not a claim that vendor CE scales match.
        platform_anchor = same_family
        platform_anchor_scope = "analyzer_family_metadata_fallback"

    query_activation = _activation(
        query.get("activation_method", query.get("activation")), query_family
    )
    reference_activation = _activation(
        reference.get(
            "lib_activation_method", reference.get("activation_method")
        ),
        reference_family,
    )
    same_activation = bool(
        query_activation
        and reference_activation
        and query_activation == reference_activation
    )
    query_mode = _mode(query.get("collision_energy_mode"))
    reference_mode = _mode(
        reference.get(
            "lib_collision_energy_mode", reference.get("collision_energy_mode")
        )
    )
    same_mode = bool(query_mode and reference_mode and query_mode == reference_mode)
    query_region = _normalized_model(query.get("collision_cell_region"))
    reference_region = _normalized_model(
        reference.get(
            "lib_collision_cell_region", reference.get("collision_cell_region")
        )
    )
    region_required = query_vendor == "waters" or reference_vendor == "waters"
    region_compatible = bool(
        query_region and reference_region and query_region == reference_region
    ) if region_required else bool(
        not query_region or not reference_region or query_region == reference_region
    )
    same_unit = bool(
        query_energy.get("unit")
        and query_energy.get("unit") == reference_energy.get("unit")
    )
    query_value = _finite(query_energy.get("value"))
    reference_value = _finite(reference_energy.get("value"))
    delta = (
        abs(query_value - reference_value)
        if same_unit and query_value is not None and reference_value is not None
        else None
    )
    unit = str(query_energy.get("unit") or "")
    threshold = _finite(near_delta_by_unit.get(unit))
    near = bool(
        delta is not None
        and threshold is not None
        and delta <= threshold
        and same_mode
        and query_mode == "fixed"
    )

    conditioned_eligible = bool(
        platform_anchor
        and same_activation
        and region_compatible
        and near
    )
    if same_domain and conditioned_eligible:
        tier = "same_instrument_domain_near_ce"
    elif same_model and conditioned_eligible:
        tier = "same_model_near_ce"
    elif same_vendor_family and conditioned_eligible:
        tier = "same_vendor_analyzer_near_ce"
    elif conditioned_eligible:
        tier = "analyzer_family_metadata_fallback_near_ce"
    elif same_family and delta is not None and not same_mode:
        tier = "same_analyzer_family_ce_mode_not_comparable"
    elif same_family and not same_activation:
        tier = "same_analyzer_family_activation_not_comparable"
    elif same_family and delta is not None:
        tier = "same_analyzer_family_distant_or_domain_incompatible_ce"
    elif query_family and reference_family and delta is not None:
        tier = "cross_analyzer_same_ce_unit"
    elif same_family:
        tier = "same_analyzer_family_ce_not_comparable"
    elif query_family and reference_family:
        tier = "cross_analyzer_ce_not_comparable"
    else:
        tier = "instrument_or_ce_metadata_incomplete"

    return {
        "tier": tier,
        "query_instrument_family": query_family or "unknown",
        "reference_instrument_family": reference_family or "unknown",
        "same_instrument_model": same_model,
        "same_instrument_family": same_family,
        "same_instrument_domain": same_domain,
        "same_vendor_analyzer_family": same_vendor_family,
        "query_vendor": query_vendor or "unknown",
        "reference_vendor": reference_vendor or "unknown",
        "platform_anchor_eligible": platform_anchor,
        "platform_anchor_scope": platform_anchor_scope,
        "query_activation": query_activation or None,
        "reference_activation": reference_activation or None,
        "activation_compatible": same_activation,
        "query_collision_energy_mode": query_mode or None,
        "reference_collision_energy_mode": reference_mode or None,
        "collision_energy_mode_compatible": same_mode,
        "collision_cell_region_compatible": region_compatible,
        "query_collision_energy": query_value,
        "reference_collision_energy": reference_value,
        "collision_energy_unit": unit if same_unit else None,
        "collision_energy_same_unit": same_unit,
        "collision_energy_delta": round(delta, 6) if delta is not None else None,
        "near_collision_energy": near,
        "eligible_for_conditioned_consensus": conditioned_eligible,
        "cross_unit_conversion_applied": False,
        "composite_energy_reduced_to_midpoint": False,
        "cross_vendor_ce_equivalence_assumed": False,
    }


def _reference_record(
    raw: Mapping[str, Any],
    query: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    peaks = (
        raw.get("lib_peaks")
        or raw.get("_score_library_peaks")
        or raw.get("peaks")
        or []
    )
    declared_fingerprint = _text(raw.get("library_peak_fingerprint")).lower()
    fingerprint = (
        declared_fingerprint
        if re.fullmatch(r"[0-9a-f]{64}", declared_fingerprint)
        else reference_content_fingerprint(peaks)
    )
    source_family = normalize_source_family(raw.get("dataset_source"))
    return {
        "spectrum_id": _text(raw.get("lib_spectrum_id") or raw.get("spectrum_id")),
        "source_accession": _text(
            raw.get("lib_source_accession") or raw.get("source_accession")
        ),
        "dataset_source": _text(raw.get("dataset_source")),
        "source_family": source_family,
        "source_tier": _text(raw.get("source_tier")),
        "spectral_score": round(_unit_interval(raw.get("spectral_score")), 8),
        "score_components": dict(raw.get("score_components") or {}),
        "content_fingerprint": fingerprint,
        "content_fingerprint_available": bool(fingerprint),
        "condition": reference_condition(
            query,
            raw,
            near_delta_by_unit=policy.get("near_collision_energy_delta") or {},
        ),
    }


def _source_inventory(
    query: Mapping[str, Any],
    records: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> tuple[list[str], bool]:
    """Return searched source families and whether absence is interpretable.

    A source can be labelled ``no_candidate_reference`` only when the search
    artifact explicitly records that the source participated.  Configured
    display rows without that inventory are labelled ``search_state_unknown``;
    they are never reported as negative evidence.
    """

    match = query.get("match") if isinstance(query.get("match"), Mapping) else {}
    explicit = (
        query.get("searched_source_families")
        or query.get("searched_library_source_families")
        or match.get("searched_source_families")
        or match.get("searched_library_source_families")
        or []
    )
    if isinstance(explicit, str):
        explicit = [explicit]
    searched = {
        normalize_source_family(value)
        for value in explicit
        if normalize_source_family(value)
    }
    # Older artifacts did not persist a source inventory.  Observed sources are
    # still known to have participated, but absence from another source cannot
    # be interpreted.
    observed = {
        str(record.get("source_family") or "")
        for record in records
        if record.get("source_family")
    }
    inventory_complete = bool(searched)
    searched.update(observed)
    return sorted(searched), inventory_complete


def _library_evidence_matrix(
    query: Mapping[str, Any],
    records: Sequence[Mapping[str, Any]],
    profiles: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
    *,
    conditioned_qualified: bool,
) -> list[dict[str, Any]]:
    configured_raw = policy.get("report_source_families") or []
    if isinstance(configured_raw, str):
        configured_raw = [configured_raw]
    configured = {
        normalize_source_family(value)
        for value in configured_raw
        if normalize_source_family(value)
    }
    searched, inventory_complete = _source_inventory(query, records, policy)
    searched_set = set(searched)
    observed = {
        str(record.get("source_family") or "")
        for record in records
        if record.get("source_family")
    }
    universe = sorted(configured | searched_set | observed)
    profile_by_id = {
        str(profile.get("profile_id") or ""): profile for profile in profiles
    }
    output: list[dict[str, Any]] = []
    for source_family in universe:
        members = [
            record
            for record in records
            if record.get("source_family") == source_family
        ]
        compatible = [
            record
            for record in members
            if bool((record.get("condition") or {}).get(
                "eligible_for_conditioned_consensus"
            ))
        ]
        if members:
            search_status = (
                "matched_conditioned"
                if compatible
                else "matched_condition_incompatible_or_incomplete"
            )
        elif inventory_complete and source_family in searched_set:
            search_status = "no_candidate_reference"
        else:
            search_status = "search_state_unknown"
        best = min(
            members,
            key=lambda record: (
                -float(record.get("spectral_score") or 0.0),
                str(record.get("spectrum_id") or ""),
            ),
            default=None,
        )
        best_compatible = min(
            compatible,
            key=lambda record: (
                -float(record.get("spectral_score") or 0.0),
                str(record.get("spectrum_id") or ""),
            ),
            default=None,
        )
        profile_ids = sorted({
            str(record.get("provenance_profile_id") or "")
            for record in members
            if record.get("provenance_profile_id")
        })
        independent_profile_ids = [
            profile_id
            for profile_id in profile_ids
            if len(profile_by_id.get(profile_id, {}).get("source_families") or []) == 1
        ]
        output.append({
            "source_family": source_family,
            "search_status": search_status,
            "source_participation_explicit": bool(
                inventory_complete and source_family in searched_set
            ),
            "reference_record_count": len(members),
            "condition_compatible_record_count": len(compatible),
            "unique_content_profile_count": len(profile_ids),
            "independent_single_source_profile_count": len(independent_profile_ids),
            "best_score": (
                round(float(best["spectral_score"]), 8) if best else None
            ),
            "best_condition_compatible_score": (
                round(float(best_compatible["spectral_score"]), 8)
                if best_compatible else None
            ),
            "representative_spectrum_id": (
                best_compatible or best or {}
            ).get("spectrum_id"),
            "representative_source_accession": (
                best_compatible or best or {}
            ).get("source_accession"),
            "representative_condition": dict(
                (best_compatible or best or {}).get("condition") or {}
            ),
            "spectrum_ids": sorted({
                str(record.get("spectrum_id") or "")
                for record in members
                if record.get("spectrum_id")
            }),
            "source_accessions": sorted({
                str(record.get("source_accession") or "")
                for record in members
                if record.get("source_accession")
            }),
            "condition_tiers": sorted({
                str((record.get("condition") or {}).get("tier") or "")
                for record in members
                if (record.get("condition") or {}).get("tier")
            }),
            "rank_influence_eligible": bool(
                conditioned_qualified and compatible
            ),
            "absence_is_negative_evidence": False,
        })
    return output


def aggregate_references(
    query: Mapping[str, Any],
    references: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    """Aggregate all references belonging to one 2-D candidate identity."""

    if not references:
        raise ValueError("at least one reference record is required")
    best_fraction = float(policy.get("best_reference_fraction", 0.5))
    consensus_fraction = float(policy.get("consensus_fraction", 0.5))
    if not math.isclose(best_fraction + consensus_fraction, 1.0, abs_tol=1e-12):
        raise ValueError("V15 reference aggregation fractions must sum to one")
    if best_fraction < 0.0 or consensus_fraction < 0.0:
        raise ValueError("V15 reference aggregation fractions cannot be negative")
    minimum_strata = int(policy.get("minimum_provenance_strata", 2))
    minimum_profiles = int(policy.get("minimum_unique_content_profiles", 2))
    if minimum_strata < 2 or minimum_profiles < 2:
        raise ValueError("V15 corroboration requires at least two strata and profiles")

    records = [_reference_record(raw, query, policy) for raw in references]
    legacy_best = max(float(record["spectral_score"]) for record in records)

    by_fingerprint: dict[str, list[dict[str, Any]]] = defaultdict(list)
    missing_fingerprint = 0
    for record in records:
        fingerprint = record.get("content_fingerprint")
        if not fingerprint:
            missing_fingerprint += 1
            continue
        by_fingerprint[str(fingerprint)].append(record)

    profiles: list[dict[str, Any]] = []
    for fingerprint, members in sorted(by_fingerprint.items()):
        members.sort(
            key=lambda record: (
                -float(record["spectral_score"]),
                str(record.get("source_family") or ""),
                str(record.get("spectrum_id") or ""),
            )
        )
        representative = members[0]
        families = sorted(
            {str(record.get("source_family")) for record in members if record.get("source_family")}
        )
        profile = {
            "profile_id": fingerprint[:20],
            "content_fingerprint": fingerprint,
            "spectral_score": float(representative["spectral_score"]),
            "source_families": families,
            "mirrored_across_source_families": len(families) > 1,
            "mirrored_record_count": len(members),
            "representative_spectrum_id": representative.get("spectrum_id"),
            "representative_source_accession": representative.get("source_accession"),
            "condition": representative.get("condition") or {},
        }
        profiles.append(profile)
        for record in members:
            record["provenance_profile_id"] = profile["profile_id"]
            record["mirrored_record_count"] = len(members)
            record["mirrored_across_source_families"] = len(families) > 1

    # One source can deposit many spectra.  It contributes one maximum score.
    # Peak-identical profiles connect every source family that carries them.
    # All profiles from one connected source group then share one vote, so
    # copying several profiles between the same libraries cannot manufacture
    # one synthetic stratum per copied profile.
    source_parents: dict[str, str] = {}

    def source_root(source: str) -> str:
        source_parents.setdefault(source, source)
        while source_parents[source] != source:
            source_parents[source] = source_parents[source_parents[source]]
            source = source_parents[source]
        return source

    def connect_sources(left: str, right: str) -> None:
        left_root = source_root(left)
        right_root = source_root(right)
        if left_root == right_root:
            return
        first, second = sorted((left_root, right_root))
        source_parents[second] = first

    for profile in profiles:
        families = sorted(set(profile.get("source_families") or []))
        for family in families:
            source_root(family)
        for family in families[1:]:
            connect_sources(families[0], family)

    source_groups: dict[str, list[str]] = defaultdict(list)
    for family in sorted(source_parents):
        source_groups[source_root(family)].append(family)
    stratum_by_source = {
        family: (
            members[0]
            if len(members) == 1
            else f"source_group:{'+'.join(members)}"
        )
        for members in source_groups.values()
        for family in members
    }

    stratum_profiles: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for profile in profiles:
        families = profile["source_families"]
        if not families:
            continue
        stratum = stratum_by_source[families[0]]
        stratum_profiles[stratum].append(profile)

    def stratum_maxima(*, conditioned_only: bool) -> list[dict[str, Any]]:
        maxima = []
        for stratum, members in sorted(stratum_profiles.items()):
            eligible = [
                profile
                for profile in members
                if not conditioned_only
                or bool((profile.get("condition") or {}).get(
                    "eligible_for_conditioned_consensus"
                ))
            ]
            if not eligible:
                continue
            best = min(
                eligible,
                key=lambda profile: (
                    -float(profile["spectral_score"]), profile["profile_id"]
                ),
            )
            maxima.append(
                {
                    "provenance_stratum": stratum,
                    "profile_id": best["profile_id"],
                    "spectral_score": float(best["spectral_score"]),
                    "condition_tier": (best.get("condition") or {}).get("tier"),
                }
            )
        return maxima

    all_maxima = stratum_maxima(conditioned_only=False)
    conditioned_maxima = stratum_maxima(conditioned_only=True)
    unique_profiles = len(profiles)
    all_qualified = bool(
        unique_profiles >= minimum_profiles and len(all_maxima) >= minimum_strata
    )
    conditioned_qualified = bool(
        unique_profiles >= minimum_profiles
        and len(conditioned_maxima) >= minimum_strata
    )
    # Only condition-compatible corroboration can change ranking.  The
    # all-condition median remains useful for diagnosing library disagreement,
    # but it is never substituted when acquisition metadata are missing.
    selected_maxima = conditioned_maxima if conditioned_qualified else []
    selected_basis = (
        "compatible_instrument_domain_near_ce_stratum_median"
        if conditioned_qualified
        else "legacy_best_reference_no_condition_qualified_corroboration"
    )
    consensus_score = (
        float(statistics.median(row["spectral_score"] for row in selected_maxima))
        if selected_maxima
        else None
    )
    aggregate = (
        best_fraction * legacy_best + consensus_fraction * float(consensus_score)
        if consensus_score is not None
        else legacy_best
    )
    aggregate = max(0.0, min(1.0, aggregate))
    all_condition_consensus_score = (
        float(statistics.median(row["spectral_score"] for row in all_maxima))
        if all_qualified
        else None
    )
    library_matrix = _library_evidence_matrix(
        query,
        records,
        profiles,
        policy,
        conditioned_qualified=conditioned_qualified,
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "status": "scored" if selected_maxima else "legacy_score_retained",
        "legacy_best_reference_score": round(legacy_best, 8),
        "provenance_consensus_score": (
            round(consensus_score, 8) if consensus_score is not None else None
        ),
        "conditioned_multilibrary_score": round(aggregate, 8),
        "selected_consensus_basis": selected_basis,
        "rank_influence_basis": selected_basis,
        "rank_influence_allowed": bool(selected_maxima),
        "runtime_influence": not math.isclose(aggregate, legacy_best, abs_tol=1e-12),
        "unconditioned_corroboration_available": all_qualified,
        "display_only_unconditioned_consensus_score": (
            round(all_condition_consensus_score, 8)
            if all_condition_consensus_score is not None else None
        ),
        "unconditioned_consensus_affects_ranking": False,
        "score_is_probability": False,
        "reference_record_count": len(records),
        "unique_content_profile_count": unique_profiles,
        "mirror_duplicate_record_count": sum(
            max(0, int(profile["mirrored_record_count"]) - 1) for profile in profiles
        ),
        "missing_content_fingerprint_count": missing_fingerprint,
        "provenance_stratum_count": len(all_maxima),
        "conditioned_provenance_stratum_count": len(conditioned_maxima),
        "source_families": sorted(
            {
                family
                for profile in profiles
                for family in profile.get("source_families") or []
            }
        ),
        "selected_stratum_scores": selected_maxima,
        "display_only_all_condition_stratum_scores": all_maxima,
        "library_evidence_matrix": library_matrix,
        "provenance_profiles": profiles,
        "reference_records": records,
        "best_reference_fraction": best_fraction,
        "consensus_fraction": consensus_fraction,
        "minimum_provenance_strata": minimum_strata,
        "minimum_unique_content_profiles": minimum_profiles,
        "accessions_count_as_independent_evidence": False,
        "mirrored_sources_count_as_independent_evidence": False,
        "technical_replicates_count_as_independent_evidence": False,
        "cross_unit_collision_energy_conversion_applied": False,
        "compound_names_used_for_scoring": False,
    }


__all__ = [
    "SCHEMA_VERSION",
    "aggregate_references",
    "normalize_source_family",
    "reference_condition",
    "reference_content_fingerprint",
]
