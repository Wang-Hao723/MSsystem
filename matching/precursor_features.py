"""Deterministic MS2 precursor-event clustering and evidence accounting.

These clusters group repeated DDA / stepped-collision-energy spectra acquired
for the same precursor event.  They are deliberately named *precursor
features*, not quantified LC-MS features: robust chromatographic feature
finding and ion-identity grouping require the continuous MS1 signal.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

from config.settings import (
    PRECURSOR_FEATURE_MAX_RT_GAP_MIN,
    PRECURSOR_FEATURE_MAX_RT_SPAN_MIN,
    PRECURSOR_FEATURE_MZ_FLOOR_DA,
    PRECURSOR_FEATURE_MZ_PPM,
)
from matching.library_quality import compound_identity_key


def _number(value: Any) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _collision_energy(value: Any) -> float | None:
    if isinstance(value, dict):
        value = value.get("value")
    return _number(value)


def _reportable_scope_high(match: dict[str, Any]) -> bool:
    """Whether one match is High at its explicitly reportable identity scope."""
    scope = str(match.get("identification_scope") or "")
    if scope == "isomer_group":
        return str(match.get("isomer_group_confidence") or "") == "High"
    return str(
        match.get("two_dimensional_connectivity_confidence")
        or match.get("confidence")
        or ""
    ) == "High"


def _new_feature(index: int, row: dict, mz: float | None, rt: float | None) -> dict:
    return {
        "first_index": index,
        "polarity": str(row.get("ion_mode") or "").upper(),
        "mz": mz,
        "rt_start": rt,
        "rt_last": rt,
        "members": [(index, row)],
    }


def annotate_precursor_features(
    results: list[dict],
    *,
    mz_ppm: float = PRECURSOR_FEATURE_MZ_PPM,
    mz_floor_da: float = PRECURSOR_FEATURE_MZ_FLOOR_DA,
    max_rt_gap_min: float = PRECURSOR_FEATURE_MAX_RT_GAP_MIN,
    max_rt_span_min: float = PRECURSOR_FEATURE_MAX_RT_SPAN_MIN,
) -> dict:
    """Annotate rows in-place and return an honest precursor-feature summary.

    Clustering is answer-independent: only precursor m/z, retention time and
    polarity are used.  Candidate identity is summarized *after* the clusters
    are fixed and never changes membership or confidence.
    """
    if not results:
        return precursor_feature_summary(results)

    timed: list[tuple[float, float, int, dict]] = []
    untimed: list[tuple[int, dict, float | None, float | None]] = []
    for index, row in enumerate(results):
        if not isinstance(row, dict):
            continue
        mz = _number(row.get("precursor_mz"))
        rt = _number(row.get("retention_time"))
        if mz is None or mz <= 0 or rt is None:
            untimed.append((index, row, mz, rt))
        else:
            timed.append((rt, mz, index, row))

    timed.sort(key=lambda item: (item[0], item[1], item[2]))
    active: list[dict] = []
    features: list[dict] = []
    for rt, mz, index, row in timed:
        polarity = str(row.get("ion_mode") or "").upper()
        active = [
            feature for feature in active
            if feature["rt_last"] is not None
            and rt - feature["rt_last"] <= max_rt_gap_min
        ]
        candidates = []
        for feature in active:
            if feature["polarity"] != polarity:
                continue
            if rt - feature["rt_start"] > max_rt_span_min:
                continue
            tolerance = max(mz_floor_da, feature["mz"] * mz_ppm / 1e6)
            delta = abs(mz - feature["mz"])
            if delta <= tolerance:
                candidates.append((delta / tolerance, feature["first_index"], feature))
        if candidates:
            feature = min(candidates, key=lambda item: (item[0], item[1]))[2]
            n = len(feature["members"])
            feature["mz"] += (mz - feature["mz"]) / (n + 1)
            feature["rt_last"] = rt
            feature["members"].append((index, row))
        else:
            feature = _new_feature(index, row, mz, rt)
            active.append(feature)
            features.append(feature)

    # Missing RT/m/z cannot be safely joined to another acquisition event.
    for index, row, mz, rt in untimed:
        features.append(_new_feature(index, row, mz, rt))

    features.sort(key=lambda feature: feature["first_index"])
    for feature_number, feature in enumerate(features, start=1):
        feature_id = f"PF{feature_number:06d}"
        members = sorted(feature["members"], key=lambda item: item[0])
        rows = [row for _, row in members]
        rts = [value for value in (_number(row.get("retention_time")) for row in rows)
               if value is not None]
        energies = sorted({
            round(value, 4)
            for value in (_collision_energy(row.get("collision_energy")) for row in rows)
            if value is not None
        })
        matches = [row.get("match") for row in rows if isinstance(row.get("match"), dict)]
        identities = Counter(compound_identity_key(match) for match in matches)
        identities.pop("", None)
        dominant_identity, dominant_n = identities.most_common(1)[0] if identities else ("", 0)
        matched_n = len(matches)
        consensus_fraction = dominant_n / matched_n if matched_n else 0.0
        metadata = {
            "precursor_feature_id": feature_id,
            "precursor_feature_mz": round(feature["mz"], 6) if feature["mz"] is not None else None,
            "precursor_feature_rt_min": round(min(rts), 6) if rts else None,
            "precursor_feature_rt_max": round(max(rts), 6) if rts else None,
            "precursor_feature_rt_span_min": round(max(rts) - min(rts), 6) if rts else None,
            "precursor_feature_spectra": len(rows),
            "precursor_feature_collision_energies": energies,
            "precursor_feature_multi_energy": len(energies) >= 2,
            "precursor_feature_matched_spectra": matched_n,
            "precursor_feature_high_spectra": sum(
                _reportable_scope_high(match) for match in matches
            ),
            "precursor_feature_full_identity_high_spectra": sum(
                match.get("exact_structure_confidence", match.get("confidence"))
                == "High"
                for match in matches
            ),
            "precursor_feature_class_high_spectra": sum(
                match.get("confidence_class") == "High" for match in matches
            ),
            "precursor_feature_dominant_identity": dominant_identity,
            "precursor_feature_identity_consensus": round(consensus_fraction, 4),
            "precursor_feature_identity_consistent": (
                dominant_n >= 2 and consensus_fraction >= (2.0 / 3.0)
            ),
        }
        for row in rows:
            row.update(metadata)

    return precursor_feature_summary(results)


def precursor_feature_summary(results: list[dict]) -> dict:
    """Summarize already annotated rows without treating spectra as compounds."""
    lcms_groups: dict[str, list[dict]] = {}
    for row in results or []:
        if not isinstance(row, dict):
            continue
        feature_id = str(row.get("lcms_feature_id") or "")
        if feature_id:
            lcms_groups.setdefault(feature_id, []).append(row)
    has_lcms = bool(lcms_groups)

    # Prefer an observed chromatographic feature whenever OpenMS attached the
    # MS2 scan.  Only scans without such an attachment use the answer-neutral
    # m/z/RT/polarity precursor-event cluster.  The former implementation
    # advertised this hybrid definition but still counted every PF cluster,
    # so five repeated scans in F001_NEG_1 were silently counted as five extra
    # precursor events despite sharing measured LC-MS features.
    by_feature: dict[str, list[dict]] = {}
    fallback_feature_ids: set[str] = set()
    for index, row in enumerate(results or []):
        if not isinstance(row, dict):
            continue
        lcms_feature_id = str(row.get("lcms_feature_id") or "")
        if lcms_feature_id:
            feature_key = f"lcms:{lcms_feature_id}"
        else:
            precursor_feature_id = str(
                row.get("precursor_feature_id") or f"unassigned:{index}"
            )
            feature_key = f"fallback:{precursor_feature_id}"
            fallback_feature_ids.add(precursor_feature_id)
        by_feature.setdefault(feature_key, []).append(row)

    groups = list(by_feature.values())
    linked_rows = sum(len(group) for group in lcms_groups.values())
    return {
        "method": (
            "OpenMS continuous-MS1 chromatographic features with explicit "
            "MS2 precursor-event fallback"
            if has_lcms
            else "MS2 precursor-event clustering (m/z + RT + polarity)"
        ),
        "approximate": not has_lcms,
        "note": (
            "MS2 spectra linked to OpenMS LC-MS features use chromatographic "
            "feature identities; unlinked spectra remain explicit precursor-event "
            "fallback nodes."
            if has_lcms
            else
            "Repeated DDA/stepped-energy spectra are grouped; these are not "
            "continuous-MS1 quantified chromatographic features."
        ),
        "total_spectra": len(results or []),
        "precursor_features": len(groups),
        "lcms_feature_finding_available": has_lcms,
        "lcms_features_with_ms2": len(lcms_groups),
        "lcms_linked_ms2_spectra": linked_rows,
        "lcms_repeated_ms2_features": sum(
            len(group) >= 2 for group in lcms_groups.values()
        ),
        "fallback_precursor_features": len(fallback_feature_ids),
        "unlinked_ms2_spectra": max(
            0,
            len(results or []) - linked_rows,
        ),
        "repeated_precursor_features": sum(len(group) >= 2 for group in groups),
        "multi_energy_precursor_features": sum(
            len({
                round(energy, 4)
                for energy in (
                    _collision_energy(row.get("collision_energy"))
                    for row in group
                )
                if energy is not None
            }) >= 2
            for group in groups
        ),
        "matched_precursor_features": sum(
            any(isinstance(row.get("match"), dict) for row in group) for group in groups
        ),
        "high_precursor_features": sum(
            any(
                isinstance(row.get("match"), dict)
                and _reportable_scope_high(row["match"])
                for row in group
            )
            for group in groups
        ),
        "full_identity_high_precursor_features": sum(
            any(
                isinstance(row.get("match"), dict)
                and row["match"].get(
                    "exact_structure_confidence",
                    row["match"].get("confidence"),
                ) == "High"
                for row in group
            )
            for group in groups
        ),
        "class_high_precursor_features": sum(
            any(
                isinstance(row.get("match"), dict)
                and row["match"].get("confidence_class") == "High"
                for row in group
            )
            for group in groups
        ),
    }
