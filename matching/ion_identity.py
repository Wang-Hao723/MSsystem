"""Ion Identity Molecular Networking for one LC–MS(/MS) sample.

The module implements the evidence-bearing core of IIMN for the system's
already-detected MS2 precursor features:

* extract continuous-MS1 XICs from mzML;
* require co-elution and correlated chromatographic peak shapes;
* require a curated adduct/in-source-fragment neutral-mass relationship;
* assign one deterministic ion identity per feature and collapse connected
  ion forms to a representative feature.

No row is turned into an independent High-confidence compound solely by graph
propagation.  Linked annotations remain explicit ion-identity evidence.
"""
from __future__ import annotations

import bisect
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from matching.library_quality import compound_identity_key


PROTON = 1.007276466621
WATER = 18.01056468403
AMMONIA = 17.026549101
CARBON_DIOXIDE = 43.98982924


@dataclass(frozen=True)
class IonSpecies:
    label: str
    shift: float
    charge: int
    multimer: int = 1
    rank: int = 10
    relationship: str = "adduct"

    def neutral_mass(self, mz: float) -> float:
        return (float(mz) * abs(self.charge) - self.shift) / self.multimer


POSITIVE_IONS = (
    IonSpecies("[M+H]+", PROTON, 1, rank=1),
    IonSpecies("[M+Na]+", 22.989218, 1, rank=2),
    IonSpecies("[M+NH4]+", 18.033823, 1, rank=3),
    IonSpecies("[M+K]+", 38.963158, 1, rank=4),
    IonSpecies(
        "[M+H-H2O]+", PROTON - WATER, 1, rank=5,
        relationship="in_source_fragment",
    ),
    IonSpecies(
        "[M+H-2H2O]+", PROTON - 2 * WATER, 1, rank=6,
        relationship="in_source_fragment",
    ),
    IonSpecies(
        "[M+H-NH3]+", PROTON - AMMONIA, 1, rank=7,
        relationship="in_source_fragment",
    ),
)

NEGATIVE_IONS = (
    IonSpecies("[M-H]-", -PROTON, -1, rank=1),
    IonSpecies("[M+Cl]-", 34.969402, -1, rank=2),
    IonSpecies("[M+FA-H]-", 44.998201, -1, rank=3),
    IonSpecies("[M+Ac-H]-", 59.013851, -1, rank=4),
    IonSpecies(
        "[M-H-H2O]-", -PROTON - WATER, -1, rank=5,
        relationship="in_source_fragment",
    ),
    IonSpecies(
        "[M-H-CO2]-", -PROTON - CARBON_DIOXIDE, -1, rank=6,
        relationship="in_source_fragment",
    ),
)


def _number(value: Any) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _normalise_adduct(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


def _species_for_polarity(polarity: str) -> tuple[IonSpecies, ...]:
    return POSITIVE_IONS if str(polarity).upper() == "POS" else NEGATIVE_IONS


def infer_ion_pair(
    mz_left: float,
    mz_right: float,
    polarity: str,
    *,
    tolerance_da: float = 0.001,
    tolerance_ppm: float = 10.0,
) -> dict | None:
    """Return the best neutral-mass-consistent ion identity hypothesis."""
    if str(polarity).upper() not in {"POS", "NEG"}:
        return None
    candidates = []
    species = _species_for_polarity(polarity)
    for left in species:
        neutral_left = left.neutral_mass(mz_left)
        if neutral_left <= 0:
            continue
        for right in species:
            if left.label == right.label:
                continue
            # A source-fragment↔source-fragment mass difference is highly
            # non-specific without a detected intact molecular ion anchor.
            # Likewise, multimers/multiply charged ions require independent
            # charge/isotope evidence and are intentionally excluded from the
            # default species table.  This conservative boundary prevents the
            # dense false networks seen when every DDA precursor is treated as
            # an unconstrained ion hypothesis.
            if (
                left.relationship == "in_source_fragment"
                and right.relationship == "in_source_fragment"
            ):
                continue
            neutral_right = right.neutral_mass(mz_right)
            if neutral_right <= 0:
                continue
            neutral = (neutral_left + neutral_right) / 2.0
            tolerance = max(float(tolerance_da), neutral * float(tolerance_ppm) * 1e-6)
            error = abs(neutral_left - neutral_right)
            if error > tolerance:
                continue
            relationship = (
                "in_source_fragment"
                if "in_source_fragment" in {left.relationship, right.relationship}
                else "adduct"
            )
            # Numerical round-off can make two algebraically exact hypotheses
            # differ at ~1e-12.  Treat those as the same mass fit, then prefer
            # the chemically simpler/common ion pair by curated rank.
            candidates.append((
                round(error / tolerance, 6),
                left.rank + right.rank,
                left.label,
                right.label,
                {
                    "left_identity": left.label,
                    "right_identity": right.label,
                    "neutral_mass": neutral,
                    "neutral_mass_error_da": error,
                    "neutral_mass_error_ppm": error / neutral * 1e6,
                    "mass_tolerance_da": tolerance,
                    "mass_quality": max(0.0, 1.0 - error / tolerance),
                    "relationship": relationship,
                },
            ))
    if not candidates:
        return None
    result = min(candidates, key=lambda item: item[:4])[-1]
    return {
        key: round(value, 8) if isinstance(value, float) else value
        for key, value in result.items()
    }


def _scan_rt_minutes(spectrum: Mapping) -> float | None:
    scans = ((spectrum.get("scanList") or {}).get("scan") or [])
    scan = scans[0] if scans else {}
    raw = scan.get("scan start time")
    value = _number(raw)
    if value is None:
        return None
    unit_text = " ".join((
        str(getattr(raw, "unit_info", "") or ""),
        str(getattr(raw, "unit_name", "") or ""),
        str(scan),
    )).lower()
    if "second" in unit_text and "minute" not in unit_text:
        return value / 60.0
    return value


def _xic_intensity(
    mzs: np.ndarray,
    intensities: np.ndarray,
    target_mz: float,
    *,
    tolerance_da: float,
    tolerance_ppm: float,
) -> float:
    tolerance = max(
        float(tolerance_da), abs(target_mz) * float(tolerance_ppm) * 1e-6
    )
    start = int(np.searchsorted(mzs, target_mz - tolerance, side="left"))
    end = int(np.searchsorted(mzs, target_mz + tolerance, side="right"))
    if end <= start:
        return 0.0
    return float(np.max(intensities[start:end], initial=0.0))


def extract_targeted_ms1_traces(
    file_path: str | Path,
    targets: Sequence[Mapping],
    *,
    rt_half_window_min: float = 0.20,
    mz_tolerance_da: float = 0.001,
    mz_tolerance_ppm: float = 10.0,
) -> tuple[dict[str, dict], dict]:
    """Stream mzML once and extract XICs for supplied precursor features."""
    path = Path(file_path)
    if path.suffix.lower() not in {".mzml", ".mzxml", ".xml"}:
        return {}, {
            "status": "continuous_ms1_unavailable_for_format",
            "ms1_scan_count": 0,
        }
    try:
        from pyteomics import mzml
    except Exception as exc:
        return {}, {
            "status": "pyteomics_unavailable",
            "error": str(exc),
            "ms1_scan_count": 0,
        }

    valid_targets = []
    for target in targets:
        identifier = str(target.get("id") or "")
        target_mz = _number(target.get("mz"))
        target_rt = _number(target.get("rt"))
        if identifier and target_mz and target_mz > 0 and target_rt is not None:
            valid_targets.append((target_rt, identifier, target_mz))
    valid_targets.sort()
    centres = [item[0] for item in valid_targets]
    traces = {
        identifier: {"rt": [], "intensity": [], "target_mz": target_mz}
        for _, identifier, target_mz in valid_targets
    }
    ms1_scan_count = 0
    try:
        with mzml.read(str(path)) as reader:
            for spectrum in reader:
                if spectrum.get("ms level") != 1:
                    continue
                rt = _scan_rt_minutes(spectrum)
                mz_array = spectrum.get("m/z array")
                intensity_array = spectrum.get("intensity array")
                if rt is None or mz_array is None or intensity_array is None:
                    continue
                ms1_scan_count += 1
                mzs = np.asarray(mz_array, dtype=float)
                intensities = np.asarray(intensity_array, dtype=float)
                start = bisect.bisect_left(centres, rt - rt_half_window_min)
                end = bisect.bisect_right(centres, rt + rt_half_window_min)
                for _, identifier, target_mz in valid_targets[start:end]:
                    traces[identifier]["rt"].append(float(rt))
                    traces[identifier]["intensity"].append(_xic_intensity(
                        mzs,
                        intensities,
                        target_mz,
                        tolerance_da=mz_tolerance_da,
                        tolerance_ppm=mz_tolerance_ppm,
                    ))
    except Exception as exc:
        return {}, {
            "status": "ms1_trace_extraction_failed",
            "error": str(exc),
            "ms1_scan_count": ms1_scan_count,
        }

    usable = 0
    for trace in traces.values():
        intensities = np.asarray(trace["intensity"], dtype=float)
        trace["points"] = len(trace["rt"])
        trace["nonzero_points"] = int(np.count_nonzero(intensities > 0))
        trace["apex_intensity"] = float(np.max(intensities, initial=0.0))
        if trace["apex_intensity"] > 0 and trace["points"] >= 5:
            apex_index = int(np.argmax(intensities))
            trace["apex_rt"] = float(trace["rt"][apex_index])
            trace["status"] = "usable"
            usable += 1
        else:
            trace["apex_rt"] = None
            trace["status"] = "insufficient_trace"
    return traces, {
        "status": "ok" if ms1_scan_count else "no_ms1_scans",
        "ms1_scan_count": ms1_scan_count,
        "target_count": len(valid_targets),
        "usable_trace_count": usable,
        "backend": "pyteomics_streaming_xic",
    }


def chromatographic_shape_evidence(
    left: Mapping,
    right: Mapping,
    *,
    min_points: int = 5,
    min_overlap: float = 0.50,
    min_correlation: float = 0.85,
    max_apex_delta_min: float = 0.10,
) -> dict:
    """Evaluate IIMN-style peak-shape and co-elution evidence."""
    left_map = {
        round(float(rt), 8): float(value)
        for rt, value in zip(left.get("rt") or [], left.get("intensity") or [])
    }
    right_map = {
        round(float(rt), 8): float(value)
        for rt, value in zip(right.get("rt") or [], right.get("intensity") or [])
    }
    common = sorted(set(left_map).intersection(right_map))
    base = {"accepted": False, "points": len(common)}
    if len(common) < min_points:
        return {**base, "status": "insufficient_common_points"}
    left_values = np.asarray([left_map[rt] for rt in common], dtype=float)
    right_values = np.asarray([right_map[rt] for rt in common], dtype=float)
    either = np.count_nonzero((left_values > 0) | (right_values > 0))
    both = np.count_nonzero((left_values > 0) & (right_values > 0))
    overlap = float(both / either) if either else 0.0
    left_apex = int(np.argmax(left_values))
    right_apex = int(np.argmax(right_values))
    apex_delta = abs(common[left_apex] - common[right_apex])
    two_sided = (
        left_apex >= 2 and right_apex >= 2
        and len(common) - left_apex - 1 >= 2
        and len(common) - right_apex - 1 >= 2
    )
    if float(np.std(left_values)) <= 0 or float(np.std(right_values)) <= 0:
        correlation = 0.0
    else:
        correlation = float(np.corrcoef(left_values, right_values)[0, 1])
        if not np.isfinite(correlation):
            correlation = 0.0
    status = "accepted"
    if not two_sided:
        status = "apex_not_two_sided"
    elif apex_delta > max_apex_delta_min:
        status = "apex_rt_mismatch"
    elif overlap < min_overlap:
        status = "insufficient_overlap"
    elif correlation < min_correlation:
        status = "low_shape_correlation"
    return {
        **base,
        "accepted": status == "accepted",
        "status": status,
        "correlation": round(correlation, 6),
        "overlap_fraction": round(overlap, 6),
        "apex_rt_delta_min": round(apex_delta, 8),
        "two_points_each_side": two_sided,
    }


def _feature_rows(results: list[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for index, row in enumerate(results):
        if not isinstance(row, dict):
            continue
        identifier = str(
            row.get("lcms_feature_id")
            or row.get("precursor_feature_id")
            or f"ROW{index:06d}"
        )
        grouped[identifier].append(row)
    features = []
    for identifier, rows in grouped.items():
        mz_values = [_number(row.get("precursor_mz")) for row in rows]
        mz_values = [value for value in mz_values if value is not None and value > 0]
        rt_values = [_number(row.get("retention_time")) for row in rows]
        rt_values = sorted(value for value in rt_values if value is not None)
        polarities = [
            str(row.get("ion_mode") or "").upper() for row in rows
            if str(row.get("ion_mode") or "").upper() in {"POS", "NEG"}
        ]
        features.append({
            "id": identifier,
            "rows": rows,
            "mz": float(np.median(mz_values)) if mz_values else None,
            "rt": float(np.median(rt_values)) if rt_values else None,
            "polarity": Counter(polarities).most_common(1)[0][0] if polarities else "",
        })
    return sorted(features, key=lambda item: (
        item["rt"] if item["rt"] is not None else float("inf"), item["id"]
    ))


class _UnionFind:
    def __init__(self, identifiers: Sequence[str]):
        self.parent = {identifier: identifier for identifier in identifiers}

    def find(self, identifier: str) -> str:
        parent = self.parent[identifier]
        if parent != identifier:
            self.parent[identifier] = self.find(parent)
        return self.parent[identifier]

    def union(self, left: str, right: str) -> str:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return root_left
        keep, merge = sorted((root_left, root_right))
        self.parent[merge] = keep
        return keep


def _match_strength(feature: Mapping) -> tuple:
    matches = [
        row.get("match") for row in feature.get("rows") or []
        if isinstance(row.get("match"), dict)
    ]
    rank = {"High": 3, "Medium": 2, "Low": 1, "Unknown": 0}
    confidence = max((rank.get(str(match.get("confidence")), 0) for match in matches), default=0)
    score = max((float(match.get("final_score") or 0) for match in matches), default=0.0)
    trace_apex = float((feature.get("trace") or {}).get("apex_intensity") or 0.0)
    return confidence, score, trace_apex, str(feature.get("id"))


def _credible_identity_match(match: Mapping | None) -> bool:
    """Only sufficiently supported library hits may define a network identity."""
    if not isinstance(match, Mapping) or match.get("name_only_lead"):
        return False
    if str(match.get("confidence") or "") in {"High", "Medium"}:
        return True
    try:
        cosine = float(match.get("cosine") or 0.0)
        final = float(match.get("final_score") or 0.0)
        reverse = float(match.get("reverse_match") or 0.0)
        ppm = match.get("mass_ppm")
        ppm_ok = ppm is None or abs(float(ppm)) <= 50.0
    except (TypeError, ValueError):
        return False
    return cosine >= 0.80 and final >= 0.50 and reverse >= 0.55 and ppm_ok


def _best_network_match(features: Sequence[Mapping]) -> tuple[dict | None, str | None]:
    candidates = []
    for feature in features:
        for row in feature.get("rows") or []:
            match = row.get("match")
            if isinstance(match, dict) and _credible_identity_match(match):
                candidates.append((_match_strength(feature), match, feature["id"]))
    if not candidates:
        return None, None
    _, match, feature_id = max(candidates, key=lambda item: item[0])
    return match, feature_id


def _adduct_consistency(match: dict, identity: str) -> str:
    query_adduct = _normalise_adduct(match.get("query_adduct"))
    inferred = _normalise_adduct(identity)
    if query_adduct:
        return "consistent" if query_adduct == inferred else "conflict_with_query_metadata"
    if match.get("cross_adduct_match"):
        return "not_comparable_cross_adduct_library"
    library_adduct = _normalise_adduct(match.get("adduct"))
    if library_adduct:
        return "consistent" if library_adduct == inferred else "conflict_with_library_adduct"
    return "not_available"


def annotate_ion_identity_networks(
    results: list[dict],
    *,
    file_path: str | Path | None = None,
    traces: Mapping[str, Mapping] | None = None,
    max_rt_delta_min: float = 0.10,
    min_shape_correlation: float = 0.85,
    min_overlap_fraction: float = 0.50,
) -> dict:
    """Annotate serialized result rows in-place and return an IIMN summary."""
    features = _feature_rows(results)
    summary = {
        "method": "Ion Identity Molecular Networking (MS1 XIC + neutral mass)",
        "status": "not_run",
        "feature_count": len(features),
        "network_count": 0,
        "networked_features": 0,
        "collapsed_ion_forms": 0,
        "edge_count": 0,
        "adduct_edges": 0,
        "in_source_fragment_edges": 0,
        "propagated_annotations": 0,
        "identity_conflict_networks": 0,
        "adduct_conflicts": 0,
        "runtime_influence": False,
    }
    if not features:
        summary["status"] = "no_precursor_features"
        return summary

    trace_map: dict[str, Mapping]
    extraction = {"status": "provided_traces"}
    if traces is not None:
        trace_map = {str(key): value for key, value in traces.items()}
    elif file_path is not None:
        trace_map, extraction = extract_targeted_ms1_traces(file_path, features)
    else:
        trace_map = {}
        extraction = {"status": "continuous_ms1_not_supplied", "ms1_scan_count": 0}
    summary["trace_extraction"] = extraction
    summary["usable_trace_count"] = sum(
        len((trace_map.get(feature["id"]) or {}).get("rt") or []) >= 5
        and max(
            (trace_map.get(feature["id"]) or {}).get("intensity") or [0]
        ) > 0
        for feature in features
    )
    extraction_status = str(
        extraction.get("status") or "continuous_ms1_unavailable"
    )
    # ``extract_targeted_ms1_traces`` deliberately returns one empty
    # placeholder trace per target when an mzML contains no survey MS1 scans.
    # Testing only ``if not trace_map`` therefore mislabeled such files as
    # successfully evaluated (status=ok, zero networks).  Stop before graph
    # construction and expose the real applicability failure.
    if extraction_status not in {"ok", "provided_traces"}:
        summary["status"] = extraction_status
        for feature in features:
            for row in feature["rows"]:
                row["ion_identity_status"] = summary["status"]
        return summary
    if not trace_map or summary["usable_trace_count"] == 0:
        summary["status"] = "no_usable_ms1_traces"
        for feature in features:
            for row in feature["rows"]:
                row["ion_identity_status"] = summary["status"]
        return summary

    by_id = {feature["id"]: feature for feature in features}
    for feature in features:
        feature["trace"] = trace_map.get(feature["id"], {})

    candidate_edges = []
    for left_index, left in enumerate(features):
        if left["mz"] is None or left["rt"] is None or not left["polarity"]:
            continue
        for right in features[left_index + 1:]:
            if right["rt"] is None or right["rt"] - left["rt"] > max_rt_delta_min:
                break
            if right["polarity"] != left["polarity"] or right["mz"] is None:
                continue
            hypothesis = infer_ion_pair(left["mz"], right["mz"], left["polarity"])
            if hypothesis is None:
                continue
            shape = chromatographic_shape_evidence(
                left.get("trace") or {},
                right.get("trace") or {},
                min_overlap=min_overlap_fraction,
                min_correlation=min_shape_correlation,
                max_apex_delta_min=max_rt_delta_min,
            )
            if not shape["accepted"]:
                continue
            score = (
                0.55 * float(shape["correlation"])
                + 0.25 * float(shape["overlap_fraction"])
                + 0.20 * float(hypothesis["mass_quality"])
            )
            candidate_edges.append({
                "left": left["id"],
                "right": right["id"],
                **hypothesis,
                "shape": shape,
                "score": round(score, 6),
            })

    candidate_edges.sort(key=lambda edge: (
        -edge["score"], edge["left"], edge["right"],
        edge["left_identity"], edge["right_identity"],
    ))
    union = _UnionFind(list(by_id))
    assignments: dict[str, str] = {}
    component_neutral: dict[str, list[float]] = defaultdict(list)
    accepted_edges = []
    for edge in candidate_edges:
        left = edge["left"]
        right = edge["right"]
        if left in assignments and assignments[left] != edge["left_identity"]:
            continue
        if right in assignments and assignments[right] != edge["right_identity"]:
            continue
        root_left = union.find(left)
        root_right = union.find(right)
        known_neutral = component_neutral[root_left] + component_neutral[root_right]
        neutral = float(edge["neutral_mass"])
        if known_neutral:
            centre = float(np.median(known_neutral))
            tolerance = max(0.001, centre * 10.0e-6)
            if abs(neutral - centre) > tolerance:
                continue
        assignments[left] = edge["left_identity"]
        assignments[right] = edge["right_identity"]
        merged_values = [*known_neutral, neutral]
        new_root = union.union(left, right)
        component_neutral.pop(root_left, None)
        component_neutral.pop(root_right, None)
        component_neutral[new_root] = merged_values
        accepted_edges.append(edge)

    groups: dict[str, list[dict]] = defaultdict(list)
    for identifier in assignments:
        groups[union.find(identifier)].append(by_id[identifier])
    groups = {
        root: sorted(items, key=lambda item: item["id"])
        for root, items in groups.items() if len(items) >= 2
    }
    edge_by_feature: dict[str, list[dict]] = defaultdict(list)
    for edge in accepted_edges:
        if union.find(edge["left"]) not in groups:
            continue
        edge_by_feature[edge["left"]].append(edge)
        edge_by_feature[edge["right"]].append(edge)

    summary["status"] = "ok"
    summary["network_count"] = len(groups)
    summary["networked_features"] = sum(len(items) for items in groups.values())
    summary["collapsed_ion_forms"] = max(
        0, summary["networked_features"] - summary["network_count"]
    )
    summary["edge_count"] = sum(
        union.find(edge["left"]) in groups for edge in accepted_edges
    )
    summary["adduct_edges"] = sum(
        edge["relationship"] == "adduct" and union.find(edge["left"]) in groups
        for edge in accepted_edges
    )
    summary["in_source_fragment_edges"] = sum(
        edge["relationship"] == "in_source_fragment"
        and union.find(edge["left"]) in groups for edge in accepted_edges
    )
    summary["runtime_influence"] = bool(groups)

    ordered_groups = sorted(groups.values(), key=lambda items: items[0]["id"])
    for network_number, network_features in enumerate(ordered_groups, start=1):
        network_id = f"IIN{network_number:06d}"
        representative = max(network_features, key=_match_strength)
        representative_match, annotation_source_id = _best_network_match(network_features)
        identities = set()
        for feature in network_features:
            for row in feature["rows"]:
                match = row.get("match")
                if isinstance(match, dict) and _credible_identity_match(match):
                    identity = compound_identity_key(match)
                    if identity:
                        identities.add(identity)
        identity_conflict = len(identities) > 1
        if identity_conflict:
            summary["identity_conflict_networks"] += 1
        neutral_values = [
            float(edge["neutral_mass"])
            for feature in network_features for edge in edge_by_feature[feature["id"]]
        ]
        neutral_mass = float(np.median(neutral_values)) if neutral_values else None

        for feature in network_features:
            feature_id = feature["id"]
            identity = assignments.get(feature_id, "")
            feature_edges = edge_by_feature[feature_id]
            best_edge = max(feature_edges, key=lambda edge: edge["score"])
            evidence = {
                "status": "networked",
                "network_id": network_id,
                "ion_identity": identity,
                "neutral_mass": round(neutral_mass, 6) if neutral_mass is not None else None,
                "representative_feature_id": representative["id"],
                "is_representative": feature_id == representative["id"],
                "edge_count": len(feature_edges),
                "best_edge_score": best_edge["score"],
                "best_shape_correlation": best_edge["shape"]["correlation"],
                "best_overlap_fraction": best_edge["shape"]["overlap_fraction"],
                "identity_conflict": identity_conflict,
                "runtime_influence": True,
            }
            for row in feature["rows"]:
                row.update({
                    "ion_identity_status": "networked",
                    "ion_identity_network_id": network_id,
                    "ion_identity": identity,
                    "ion_identity_neutral_mass": evidence["neutral_mass"],
                    "ion_identity_representative_feature_id": representative["id"],
                    "ion_identity_is_representative": evidence["is_representative"],
                    "ion_identity_evidence": evidence,
                })
                match = row.get("match")
                if isinstance(match, dict):
                    consistency = _adduct_consistency(match, identity)
                    credible_match = _credible_identity_match(match)
                    match["ion_identity_evidence"] = {
                        **evidence,
                        "adduct_consistency": consistency,
                        "identity_defining_match": credible_match,
                    }
                    match["ion_identity_conflict"] = identity_conflict
                    match["ion_identity_adduct_conflict"] = bool(
                        credible_match and consistency.startswith("conflict")
                    )
                    if match["ion_identity_adduct_conflict"]:
                        summary["adduct_conflicts"] += 1
                elif representative_match is not None and not identity_conflict:
                    row["ion_identity_linked_annotation"] = {
                        "status": "linked_not_independent_identification",
                        "compound_name": representative_match.get("compound_name"),
                        "inchikey": representative_match.get("inchikey"),
                        "formula": representative_match.get("formula"),
                        "source_feature_id": annotation_source_id,
                        "network_id": network_id,
                        "ion_identity": identity,
                        "annotation_level": "Ion identity link; not an independent MS/MS High",
                    }
                    summary["propagated_annotations"] += 1

    for feature in features:
        if feature["id"] in assignments:
            continue
        trace = feature.get("trace") or {}
        status = "unlinked" if len(trace.get("rt") or []) >= 5 else "insufficient_ms1_trace"
        for row in feature["rows"]:
            row["ion_identity_status"] = status
    return summary
