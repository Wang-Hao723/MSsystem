"""Build the acquisition entity graph before any candidate-level aggregation.

Entity construction consumes only acquisition fields from result rows.  Match,
compound, score and confidence fields are intentionally inaccessible to the
graph builder.  Candidate consensus is computed afterwards as a separate shadow
view and cannot modify the production scan-level result.
"""

from __future__ import annotations

import hashlib
import json
from runtime.json_encoding import write_json_encoded
import math
import os
import statistics
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from matching.evidence_governance import build_feature_consensus
from pipeline.query_input_guards import charge_key, energy_condition
from schemas.entities import EntityEdge, EntityGraph, EntityLevel, EntityNode


ENTITY_BUNDLE_SCHEMA_VERSION = "1.0.1"
FALLBACK_EVENT_MZ_PPM = 15.0
FALLBACK_EVENT_MZ_FLOOR_DA = 0.01
FALLBACK_EVENT_RT_WINDOW_MIN = 0.03


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _text(value: Any) -> str | None:
    text = str(value or "").strip()
    return text or None


def _median(values: Iterable[Any]) -> float | None:
    numbers = [number for number in (_finite(value) for value in values) if number is not None]
    return statistics.median(numbers) if numbers else None


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _scan_signature(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        _text(row.get("spectrum_id")) or "",
        (_text(row.get("ion_mode")) or "").upper(),
        _finite(row.get("precursor_mz")),
        _finite(row.get("retention_time")),
        _finite(row.get("collision_energy")),
        _text(row.get("precursor_scan_ref")) or "",
        _text(row.get("precursor_feature_id")) or "",
        _text(row.get("lcms_feature_id")) or "",
        charge_key(row.get("precursor_charge", row.get("charge"))),
        json.dumps(energy_condition(row), sort_keys=True),
    )


def _scan_entries(task_id: str, rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted(
        (row for row in rows if isinstance(row, Mapping)),
        key=lambda row: tuple("" if value is None else str(value) for value in _scan_signature(row)),
    )
    occurrences: Counter[str] = Counter()
    entries = []
    for row in ordered:
        signature = _scan_signature(row)
        digest = hashlib.sha256(_canonical(signature)).hexdigest()[:16]
        occurrences[digest] += 1
        scan_entity_id = f"scan:{task_id}:{digest}:{occurrences[digest]}"
        entries.append(
            {
                "row": row,
                "scan_entity_id": scan_entity_id,
                "signature": signature,
                "polarity": (_text(row.get("ion_mode")) or "").upper(),
                "mz": _finite(row.get("precursor_mz")),
                "rt": _finite(row.get("retention_time")),
                "charge": charge_key(row.get("precursor_charge", row.get("charge"))),
            }
        )
    return entries


def _compatible_event(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    if left.get("polarity") != right.get("polarity"):
        return False
    if left.get("charge", "unknown") != right.get("charge", "unknown"):
        return False
    lmz, rmz = left.get("mz"), right.get("mz")
    lrt, rrt = left.get("rt"), right.get("rt")
    if None in {lmz, rmz, lrt, rrt}:
        return False
    mean_mz = (float(lmz) + float(rmz)) / 2.0
    mz_window = max(FALLBACK_EVENT_MZ_FLOOR_DA, mean_mz * FALLBACK_EVENT_MZ_PPM * 1e-6)
    return (
        abs(float(lmz) - float(rmz)) <= mz_window
        and abs(float(lrt) - float(rrt)) <= FALLBACK_EVENT_RT_WINDOW_MIN
    )


def _fallback_event_components(entries: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    if not entries:
        return []
    parent = list(range(len(entries)))
    members: dict[int, list[int]] = {index: [index] for index in range(len(entries))}

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    edges = []
    for left_index, left in enumerate(entries):
        for right_index in range(left_index + 1, len(entries)):
            right = entries[right_index]
            if _compatible_event(left, right):
                distance = (
                    abs(float(left["mz"]) - float(right["mz"])),
                    abs(float(left["rt"]) - float(right["rt"])),
                    left["scan_entity_id"],
                    right["scan_entity_id"],
                )
                edges.append((distance, left_index, right_index))
    edges.sort()
    for _, left_index, right_index in edges:
        left_root, right_root = find(left_index), find(right_index)
        if left_root == right_root:
            continue
        if any(
            not _compatible_event(entries[a], entries[b])
            for a in members[left_root]
            for b in members[right_root]
        ):
            continue
        keep, absorb = sorted((left_root, right_root))
        parent[absorb] = keep
        members[keep] = sorted(members[keep] + members[absorb])
        del members[absorb]
    components = [
        [entries[index] for index in indices]
        for _, indices in sorted(members.items())
    ]
    components.sort(key=lambda component: tuple(item["scan_entity_id"] for item in component))
    return components


def _event_groups(task_id: str, entries: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    fallback = []
    for entry in entries:
        row = entry["row"]
        stored = _text(row.get("precursor_feature_id"))
        reference = _text(row.get("precursor_scan_ref"))
        if stored:
            groups[f"event:{task_id}:stored:{stored}"].append(entry)
        elif reference:
            digest = hashlib.sha256(reference.encode("utf-8")).hexdigest()[:16]
            groups[f"event:{task_id}:reference:{digest}"].append(entry)
        else:
            fallback.append(entry)
    # A precursor_scan_ref identifies an MS1 scan, which can contain several
    # precursor ions. Its identifier alone cannot establish one MS2 event.
    resolved = {}
    for group_id, members in groups.items():
        if ':reference:' in group_id:
            components = _fallback_event_components(members)
        else:
            strata = defaultdict(list)
            for member in members:
                strata[(member['polarity'], member['charge'])].append(member)
            components = list(strata.values())
        if len(components) == 1:
            resolved[group_id] = members
        else:
            for component in components:
                signature = sorted(item['scan_entity_id'] for item in component)
                suffix = hashlib.sha256(_canonical(signature)).hexdigest()[:16]
                resolved[group_id + ':ion:' + suffix] = component
    groups = defaultdict(list, resolved)
    for component in _fallback_event_components(fallback):
        signature = sorted(item["scan_entity_id"] for item in component)
        digest = hashlib.sha256(_canonical(signature)).hexdigest()[:16]
        groups[f"event:{task_id}:fallback:{digest}"].extend(component)
    return dict(sorted(groups.items()))


def _node_payload(node: EntityNode) -> dict[str, Any]:
    return node.model_dump(mode="json", exclude_none=True)


def _edge_payload(edge: EntityEdge) -> dict[str, Any]:
    return edge.model_dump(mode="json", exclude_none=True)


def build_entity_graph(task_id: str, rows: Iterable[Mapping[str, Any]]) -> EntityGraph:
    entries = _scan_entries(task_id, rows)
    nodes: list[EntityNode] = []
    edges: list[EntityEdge] = []
    for entry in entries:
        row = entry["row"]
        nodes.append(
            EntityNode(
                entity_id=entry["scan_entity_id"],
                level=EntityLevel.MS2_SCAN,
                task_id=task_id,
                polarity=entry["polarity"] or None,
                representative_mz=entry["mz"],
                retention_time_min=entry["rt"],
                member_count=1,
                construction_source="measured_ms2_scan",
                attributes={
                    "spectrum_id": _text(row.get("spectrum_id")),
                    "collision_energy": _finite(row.get("collision_energy")),
                    "collision_energy_unit": _text(row.get("collision_energy_unit")),
                    "precursor_charge": entry['charge'] if isinstance(entry['charge'], int) else None,
                    "collision_energy_condition": energy_condition(row),
                    "precursor_scan_ref": _text(row.get("precursor_scan_ref")),
                },
            )
        )

    event_groups = _event_groups(task_id, entries)
    scan_to_event: dict[str, str] = {}
    event_to_entries: dict[str, list[dict[str, Any]]] = {}
    for event_id, members in event_groups.items():
        event_to_entries[event_id] = members
        source = (
            "stored_candidate_free_precursor_feature"
            if ":stored:" in event_id
            else "referenced_precursor_scan"
            if ":reference:" in event_id
            else "fallback_complete_link_mz_rt_cluster"
        )
        rt_values = [member["rt"] for member in members]
        nodes.append(
            EntityNode(
                entity_id=event_id,
                level=EntityLevel.PRECURSOR_EVENT,
                task_id=task_id,
                polarity=members[0]["polarity"] or None,
                representative_mz=_median(member["mz"] for member in members),
                retention_time_min=_median(rt_values),
                rt_start_min=min(value for value in rt_values if value is not None)
                if any(value is not None for value in rt_values)
                else None,
                rt_end_min=max(value for value in rt_values if value is not None)
                if any(value is not None for value in rt_values)
                else None,
                member_count=len(members),
                construction_source=source,
                attributes={
                    "collision_energy_conditions": list({
                        json.dumps(condition, sort_keys=True): condition
                        for member in members
                        if (condition := energy_condition(member['row'])) is not None
                    }.values()),
                    "collision_energies": sorted(
                        {
                            value
                            for value in (
                                _finite(member["row"].get("collision_energy"))
                                for member in members
                            )
                            if value is not None
                        }
                    ),
                    "pseudoreplicates_collapsed": max(0, len(members) - 1),
                },
            )
        )
        for member in members:
            scan_to_event[member["scan_entity_id"]] = event_id
            edges.append(
                EntityEdge(
                    source_id=member["scan_entity_id"],
                    target_id=event_id,
                    relation="belongs_to_precursor_event",
                    construction_basis=source,
                )
            )

    feature_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in entries:
        feature = _text(entry["row"].get("lcms_feature_id"))
        if feature:
            feature_groups[f"feature:{task_id}:{feature}"].append(entry)
    event_to_features: dict[str, set[str]] = defaultdict(set)
    feature_to_entries: dict[str, list[dict[str, Any]]] = {}
    for feature_id, members in sorted(feature_groups.items()):
        feature_to_entries[feature_id] = members
        rows_for_feature = [member["row"] for member in members]
        nodes.append(
            EntityNode(
                entity_id=feature_id,
                level=EntityLevel.CHROMATOGRAPHIC_FEATURE,
                task_id=task_id,
                polarity=members[0]["polarity"] or None,
                representative_mz=_median(
                    row.get("lcms_feature_mz") or row.get("precursor_mz")
                    for row in rows_for_feature
                ),
                retention_time_min=_median(
                    row.get("lcms_feature_rt_minutes") or row.get("retention_time")
                    for row in rows_for_feature
                ),
                rt_start_min=_median(row.get("lcms_feature_rt_start_minutes") for row in rows_for_feature),
                rt_end_min=_median(row.get("lcms_feature_rt_end_minutes") for row in rows_for_feature),
                member_count=len(members),
                construction_source="continuous_ms1_feature_link",
                attributes={
                    "charge": next(
                        (
                            row.get("lcms_feature_charge")
                            for row in rows_for_feature
                            if row.get("lcms_feature_charge") is not None
                        ),
                        None,
                    ),
                    "isotope_trace_count": max(
                        [int(row.get("lcms_isotope_trace_count") or 0) for row in rows_for_feature]
                        or [0]
                    ),
                },
            )
        )
        for member in members:
            event_id = scan_to_event[member["scan_entity_id"]]
            event_to_features[event_id].add(feature_id)
        for event_id in sorted(
            {scan_to_event[member["scan_entity_id"]] for member in members}
        ):
            edges.append(
                EntityEdge(
                    source_id=event_id,
                    target_id=feature_id,
                    relation="linked_to_chromatographic_feature",
                    construction_basis="measured_mz_rt_feature_link",
                )
            )

    ion_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in entries:
        row = entry["row"]
        group = _text(row.get("ion_identity_network_id") or row.get("lcms_adduct_group_id"))
        feature = _text(row.get("lcms_feature_id"))
        if group and feature:
            ion_groups[f"ion_group:{task_id}:{group}"].append(entry)
    neutral_nodes = 0
    for ion_id, members in sorted(ion_groups.items()):
        feature_ids = sorted(
            {f"feature:{task_id}:{member['row'].get('lcms_feature_id')}" for member in members}
        )
        rows_for_group = [member["row"] for member in members]
        nodes.append(
            EntityNode(
                entity_id=ion_id,
                level=EntityLevel.ION_IDENTITY_GROUP,
                task_id=task_id,
                polarity=members[0]["polarity"] or None,
                representative_mz=_median(member["mz"] for member in members),
                retention_time_min=_median(member["rt"] for member in members),
                member_count=len(feature_ids),
                construction_source="candidate_free_ion_identity_network",
                attributes={
                    "feature_count": len(feature_ids),
                    "adduct_labels": sorted(
                        {
                            str(row.get("lcms_adduct"))
                            for row in rows_for_group
                            if _text(row.get("lcms_adduct"))
                        }
                    ),
                },
            )
        )
        for feature_id in feature_ids:
            edges.append(
                EntityEdge(
                    source_id=feature_id,
                    target_id=ion_id,
                    relation="member_of_ion_identity_group",
                    construction_basis="ms1_coelution_isotope_adduct_relation",
                )
            )
        neutral_mass = _median(row.get("ion_identity_neutral_mass") or row.get("lcms_neutral_mass") for row in rows_for_group)
        conflict = any(
            str(row.get("ion_identity_status") or "").lower() in {"conflict", "identity_conflict"}
            for row in rows_for_group
        )
        if neutral_mass is not None and not conflict:
            neutral_id = ion_id.replace("ion_group:", "neutral:", 1)
            nodes.append(
                EntityNode(
                    entity_id=neutral_id,
                    level=EntityLevel.NEUTRAL_ANALYTE,
                    task_id=task_id,
                    polarity=None,
                    representative_mz=neutral_mass,
                    retention_time_min=_median(member["rt"] for member in members),
                    member_count=1,
                    construction_source="ion_identity_neutral_mass",
                    attributes={
                        "identity_claimed": False,
                        "neutral_mass_only": True,
                    },
                )
            )
            edges.append(
                EntityEdge(
                    source_id=ion_id,
                    target_id=neutral_id,
                    relation="supports_neutral_analyte_entity",
                    construction_basis="candidate_free_adduct_neutral_mass",
                )
            )
            neutral_nodes += 1

    nodes.sort(key=lambda node: (node.level.value, node.entity_id))
    edges.sort(key=lambda edge: (edge.source_id, edge.target_id, edge.relation))
    count_by_level = Counter(node.level.value for node in nodes)
    summary = {
        "entity_counts": dict(sorted(count_by_level.items())),
        "edge_count": len(edges),
        "scan_count": len(entries),
        "precursor_event_count": len(event_groups),
        "chromatographic_feature_count": len(feature_groups),
        "ion_identity_group_count": len(ion_groups),
        "neutral_analyte_entity_count": neutral_nodes,
        "events_without_ms1_feature": sum(
            event_id not in event_to_features for event_id in event_groups
        ),
        "events_with_multiple_ms1_features": sum(
            len(feature_ids) > 1 for feature_ids in event_to_features.values()
        ),
        "same_injection_extra_scans_collapsed": sum(
            max(0, len(members) - 1) for members in event_groups.values()
        ),
        "candidate_independent": True,
        "production_ranking_changed": False,
    }
    graph_payload = {
        "task_id": task_id,
        "nodes": [_node_payload(node) for node in nodes],
        "edges": [_edge_payload(edge) for edge in edges],
        "summary": summary,
    }
    graph_sha256 = hashlib.sha256(_canonical(graph_payload)).hexdigest()
    return EntityGraph(
        task_id=task_id,
        nodes=tuple(nodes),
        edges=tuple(edges),
        summary=summary,
        construction_inputs=(
            "spectrum_id",
            "ion_mode",
            "precursor_mz",
            "retention_time",
            "collision_energy",
            "precursor_scan_ref",
            "precursor_feature_id",
            "lcms_feature_id",
            "lcms_feature_mz",
            "lcms_feature_rt_minutes",
            "ion_identity_network_id",
            "lcms_adduct_group_id",
            "ion_identity_neutral_mass",
        ),
        candidate_independent=True,
        graph_sha256=graph_sha256,
    )


def entity_bundle_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}_entity_bundle.json"


def build_entity_bundle(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    # The graph and the downstream candidate consensus must observe the exact
    # same immutable row set.  Materialize once so callers can safely provide
    # the streaming legacy-artifact iterator for large QE tasks.
    materialized = [dict(row) for row in rows if isinstance(row, Mapping)]
    graph = build_entity_graph(task_id, materialized)
    consensus = build_feature_consensus(materialized, task_id)
    return {
        "schema_version": ENTITY_BUNDLE_SCHEMA_VERSION,
        "task_id": task_id,
        "entity_graph": graph.model_dump(mode="json", exclude_none=True),
        "feature_candidate_consensus": consensus,
        "summary": {
            **graph.summary,
            "feature_consensus_count": len(consensus),
            "feature_consensus_runtime_influence": False,
            "feature_consensus_ranking_changed": sum(
                bool(item.get("ranking_changed")) for item in consensus
            ),
            "candidate_aggregation_after_entity_construction": True,
            "expected_compounds_counted": False,
        },
        "invariants": {
            "entity_construction_uses_candidate_identity": False,
            "stepped_energy_scans_are_independent_votes": False,
            "feature_consensus_changes_production_top1": False,
            "neutral_entity_contains_library_identity": False,
        },
    }


def write_entity_bundle(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
    path: Path,
) -> dict[str, Any]:
    payload = build_entity_bundle(task_id, rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        write_json_encoded(payload, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    summary = dict(payload["summary"])
    summary.update(
        {
            "status": "available",
            "schema_version": payload["schema_version"],
            "graph_sha256": payload["entity_graph"]["graph_sha256"],
            "artifact": str(path.resolve()),
            "size_bytes": path.stat().st_size,
        }
    )
    return summary


def read_entity_bundle(
    path: Path,
    *,
    node_offset: int = 0,
    node_limit: int = 500,
    consensus_offset: int = 0,
    consensus_limit: int = 500,
) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    graph = dict(payload.get("entity_graph") or {})
    nodes = list(graph.pop("nodes", []) or [])
    edges = list(graph.pop("edges", []) or [])
    consensus = list(payload.get("feature_candidate_consensus") or [])
    return {
        "schema_version": payload.get("schema_version"),
        "task_id": payload.get("task_id"),
        "summary": payload.get("summary") or {},
        "invariants": payload.get("invariants") or {},
        "entity_graph": graph,
        "nodes": nodes[node_offset : node_offset + node_limit],
        "nodes_total": len(nodes),
        "node_offset": node_offset,
        "node_limit": node_limit,
        "nodes_has_more": node_offset + node_limit < len(nodes),
        "edges": edges,
        "feature_candidate_consensus": consensus[
            consensus_offset : consensus_offset + consensus_limit
        ],
        "feature_consensus_total": len(consensus),
        "consensus_offset": consensus_offset,
        "consensus_limit": consensus_limit,
        "consensus_has_more": consensus_offset + consensus_limit < len(consensus),
    }
