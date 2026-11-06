"""Measured-MS1 ion networks, with explicit links back to real MS2 rows.

The input to the native IIMN builder contains one row per measured MS1
feature, not one per selected MS2 ion. Unfragmented partners may support an
ion-form hypothesis but never become MS/MS records or identified compounds.
This module returns evidence only: it does not mutate results, propagate
library identities, or change production ranking/counts.
"""
from __future__ import annotations

import copy
import math
import statistics
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence

from matching.ion_identity import (
    NEGATIVE_IONS, POSITIVE_IONS, annotate_ion_identity_networks,
)
from matching.v16_dual_polarity import semantic_sha256

ALGORITHM_ID = "mssystem-v16-measured-ms1-ion-network-1.0.0"
VARIANTS = ("linked_ms1_coordinates", "all_ms1_partners")
RESULT_FIELDS = (
    "precursor_mz", "retention_time", "retention_time_unit", "ion_mode",
    "polarity", "precursor_charge", "lcms_feature_id", "lcms_feature_mz",
    "lcms_feature_rt_minutes", "lcms_feature_link_status",
)


def _finite(value):
    if isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return value if math.isfinite(value) else None


def _polarity(value):
    return {"POSITIVE": "POS", "NEGATIVE": "NEG", "+": "POS", "-": "NEG"}.get(
        str(value or "").strip().upper(), str(value or "").strip().upper())


def prepare_ms1_network_input(
    task_id: str, polarity: str, features: Sequence[Mapping], rows: Sequence[Mapping],
) -> dict:
    """Project measured coordinates and validate existing, unique MS2 links.

    No matching by candidate mass/name is permitted. Malformed features and
    unusable row links are reported; duplicate feature IDs and mixed task
    polarity fail closed because their provenance cannot be resolved.
    """
    task_id, polarity = str(task_id or "").strip(), _polarity(polarity)
    if not task_id or polarity not in {"POS", "NEG"}:
        raise ValueError("An explicit task ID and POS/NEG polarity are required")
    nodes, rejected, seen = [], [], set()
    for index, feature in enumerate(features):
        if not isinstance(feature, Mapping):
            raise ValueError("Every feature must be a mapping")
        identifier = str(feature.get("feature_id") or "").strip()
        if not identifier or identifier in seen:
            raise ValueError("MS1 feature IDs must be nonempty and unique")
        seen.add(identifier)
        mz = _finite(feature.get("mz"))
        minutes, seconds = _finite(feature.get("rt_minutes")), _finite(feature.get("rt_seconds"))
        rt = minutes if minutes is not None else seconds / 60 if seconds is not None else None
        charge = _finite(feature.get("charge")) if feature.get("charge") is not None else 0.0
        reason = None
        if (mz is None or mz <= 0 or rt is None or rt < 0
                or (feature.get("rt_minutes") is not None and minutes is None)
                or (feature.get("rt_seconds") is not None and seconds is None)):
            reason = "invalid_measured_coordinates"
        elif minutes is not None and seconds is not None and abs(minutes - seconds / 60) > 1e-6:
            reason = "inconsistent_rt_units"
        elif charge is None or abs(charge) not in {0, 1}:
            reason = "unsupported_charge_magnitude"
        if reason:
            rejected.append({"feature_id": identifier, "feature_index": index, "reason": reason})
            continue
        nodes.append({"feature_id": identifier, "mz": mz, "rt_minutes": rt,
            "charge_magnitude": int(abs(charge)),
            "charge_evidence": "not_measured_single_charge_hypothesis" if charge == 0 else "measured_single_charge"})
    nodes.sort(key=lambda f: f["feature_id"])
    by_id = {node["feature_id"]: node for node in nodes}
    links = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError("Every result row must be a mapping; row indices cannot be dropped")
        if any(row.get(key) and _polarity(row[key]) != polarity for key in ("ion_mode", "polarity")):
            raise ValueError("Result polarity conflicts with the explicit task polarity")
        identifier = str(row.get("lcms_feature_id") or "").strip()
        node = by_id.get(identifier)
        mz, selected = _finite(row.get("lcms_feature_mz")), _finite(row.get("precursor_mz"))
        linked_rt = _finite(row.get("lcms_feature_rt_minutes"))
        charge = _finite(row.get("precursor_charge")) if row.get("precursor_charge") is not None else 0.0
        if row.get("lcms_feature_link_status") != "linked_openms_feature":
            status = "no_unique_ms1_link"
        elif node is None:
            status = "linked_feature_missing_or_ineligible"
        elif mz is None or abs(mz - node["mz"]) > 1e-6:
            status = "linked_feature_coordinate_mismatch"
        elif row.get("lcms_feature_rt_minutes") is not None and (
                linked_rt is None or abs(linked_rt - node["rt_minutes"]) > 1e-6):
            status = "linked_feature_rt_mismatch"
        elif selected is None or selected <= 0 or abs(selected - mz) > max(.02, mz * 50e-6):
            status = "selected_ion_outside_existing_link_guard"
        elif charge is None or abs(charge) not in {0, 1}:
            status = "unsupported_ms2_charge_magnitude"
        else:
            status = "verified_existing_ms1_link"
        links.append({"source_row_index": index, "feature_id": identifier or None, "status": status})
    linked = {link["feature_id"] for link in links if link["status"] == "verified_existing_ms1_link"}
    material = {"task_id": task_id, "polarity": polarity, "features": nodes, "row_bindings": links}
    return {**material, "algorithm_id": ALGORITHM_ID,
        "candidate_independent_projection_sha256": semantic_sha256(material),
        "rejected_features": rejected,
        "counts": {"input_ms1_features": len(features), "eligible_ms1_features": len(nodes),
            "linked_ms1_features": len(linked), "unfragmented_ms1_features": len(nodes) - len(linked),
            "raw_ms2_rows": len(rows), "row_binding_status_counts": dict(sorted(Counter(x["status"] for x in links).items())),
            "rejected_feature_reasons": dict(sorted(Counter(x["reason"] for x in rejected).items()))}}


def ms1_trace_targets(prepared: Mapping) -> list[dict]:
    """Use one trace per measured feature; both ablations share this extraction."""
    return [{"id": f["feature_id"], "mz": f["mz"], "rt": f["rt_minutes"]}
            for f in prepared["features"]]


def build_ms1_ion_network(prepared: Mapping, traces: Mapping, *, variant: str = "all_ms1_partners") -> dict:
    """Build candidate-free evidence without altering any input object."""
    if variant not in VARIANTS or prepared.get("algorithm_id") != ALGORITHM_ID:
        raise ValueError("Unsupported measured-MS1 network variant or input")
    material = {k: prepared[k] for k in ("task_id", "polarity", "features", "row_bindings")}
    if prepared.get("candidate_independent_projection_sha256") != semantic_sha256(material):
        raise ValueError("Measured-MS1 input projection changed after preparation")
    task_id, polarity = prepared["task_id"], prepared["polarity"]
    row_indices = defaultdict(list)
    for link in prepared["row_bindings"]:
        if link["status"] == "verified_existing_ms1_link":
            row_indices[link["feature_id"]].append(link["source_row_index"])
    selected = [f for f in prepared["features"]
                if variant == "all_ms1_partners" or f["feature_id"] in row_indices]
    # Native IIMN sees only measured coordinates: no library name, match, score,
    # formula, reference adduct, truth label or previously inferred group.
    network_rows = [{"lcms_feature_id": f["feature_id"], "precursor_mz": f["mz"],
        "retention_time": f["rt_minutes"], "ion_mode": polarity} for f in selected]
    native = annotate_ion_identity_networks(network_rows, traces=traces)
    if native.get("propagated_annotations") or native.get("identity_conflict_networks"):
        raise AssertionError("Candidate identities entered a measured-MS1 graph")
    by_id = {f["feature_id"]: f for f in selected}
    species = {ion.label: ion for ion in (POSITIVE_IONS if polarity == "POS" else NEGATIVE_IONS)}
    groups = defaultdict(list)
    for row in network_rows:
        if row.get("ion_identity_status") == "networked":
            groups[row["ion_identity_network_id"]].append(row)
    networks, anchors, rejected = [], [], []
    for _, members in sorted(groups.items()):
        members.sort(key=lambda r: r["lcms_feature_id"])
        ids = [r["lcms_feature_id"] for r in members]
        neutral_values = [species[r["ion_identity"]].neutral_mass(by_id[r["lcms_feature_id"]]["mz"])
                          for r in members]
        mass = statistics.median(neutral_values)
        spread = max(neutral_values) - min(neutral_values)
        if mass <= 0 or spread > max(.001, mass * 10e-6):
            rejected.append({"feature_ids": ids, "reason": "full_component_neutral_mass_spread",
                "neutral_mass": mass, "neutral_mass_spread_da": spread})
            continue
        actual_rows = sorted({index for fid in ids for index in row_indices[fid]})
        unfragmented = sorted(fid for fid in ids if not row_indices[fid])
        rt = statistics.median(by_id[fid]["rt_minutes"] for fid in ids)
        assignments = [{"feature_id": r["lcms_feature_id"], "ion_identity": r["ion_identity"],
            "measured_mz": by_id[r["lcms_feature_id"]]["mz"],
            "retention_time_min": by_id[r["lcms_feature_id"]]["rt_minutes"],
            "neutral_mass_hypothesis": neutral,
            "charge_evidence": by_id[r["lcms_feature_id"]]["charge_evidence"],
            "source_row_indices": list(row_indices[r["lcms_feature_id"]]),
            "best_shape_correlation": r["ion_identity_evidence"]["best_shape_correlation"],
            "best_overlap_fraction": r["ion_identity_evidence"]["best_overlap_fraction"]}
            for r, neutral in zip(members, neutral_values)]
        group_id = "V16MS1-" + semantic_sha256({"task_id": task_id, "polarity": polarity,
            "members": assignments})[:16]
        network = {"group_id": group_id, "feature_ids": ids, "source_row_indices": actual_rows,
            "ms1_only_feature_ids": unfragmented, "member_evidence": assignments,
            "neutral_mass": round(mass, 8), "neutral_mass_spread_da": round(spread, 10),
            "retention_time_min": round(rt, 8), "candidate_independent": True,
            "evidence_level": "ms1_ion_form_hypothesis_not_compound_identification"}
        networks.append(network)
        if not actual_rows:
            continue
        anchor_material = {"task_id": task_id, "polarity": polarity,
            "source": "reliable_iimn", "group_id": group_id}
        anchors.append({"anchor_id": "DPA-" + semantic_sha256(anchor_material)[:16],
            **anchor_material, "feature_ids": ids, "source_row_indices": actual_rows,
            "neutral_mass": round(mass, 8), "retention_time_min": round(rt, 8),
            "precursor_mz": round(statistics.median(by_id[fid]["mz"] for fid in ids), 8),
            "adducts": sorted({r["ion_identity"] for r in members}),
            "verified_member_feature_count": len(ids), "candidate_independent": True})
    networks.sort(key=lambda n: n["group_id"])
    anchors.sort(key=lambda a: a["anchor_id"])
    result = {"algorithm_id": ALGORITHM_ID, "variant": variant, "task_id": task_id,
        "polarity": polarity, "candidate_independent": True, "production_results_mutated": False,
        "runtime_effect": False, "compound_count_increment": 0,
        "input_projection_sha256": prepared["candidate_independent_projection_sha256"],
        "input_counts": copy.deepcopy(prepared["counts"]),
        "summary": {"status": native["status"], "tested_ms1_feature_count": len(selected),
            "network_count": len(networks), "ms2_bound_anchor_count": len(anchors),
            "ms1_only_network_count": sum(not n["source_row_indices"] for n in networks),
            "networks_with_ms1_only_partners": sum(bool(n["ms1_only_feature_ids"]) for n in networks),
            "rejected_component_count": len(rejected), "usable_trace_count": native.get("usable_trace_count", 0)},
        "networks": networks, "anchors": anchors, "rejected_components": rejected}
    result["semantic_sha256"] = semantic_sha256(result)
    return result
