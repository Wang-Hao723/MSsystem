"""Read-only V16 cross-polarity evidence graph.

The graph joins one explicitly selected positive-mode task and one explicitly
selected negative-mode task.  Edges are constructed only from candidate-
independent neutral-mass evidence and retention time.  Candidate connectivity
and optional cross-polarity MS/MS similarity are attached *after* the graph is
frozen and therefore cannot create an edge, alter a winner, or increase the
reported compound count.

This module intentionally has no task-store, database, or web dependency.  It
accepts task metadata and already-materialized result rows and returns a
deterministic JSON-compatible sidecar payload without mutating either input.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


POLICY_SCHEMA_VERSION = "mssystem.v16.dual-polarity-shadow.policy/1.0.0"
SIDECAR_SCHEMA_VERSION = "mssystem.v16.dual-polarity-shadow.sidecar/1.0.0"
MODULE_ID = "v16_dual_polarity_evidence_graph"
ALGORITHM_ID = "mssystem-v16-dual-polarity-evidence-graph-1.0.0"
DEFAULT_POLICY_PATH = (
    Path(__file__).resolve().parent.parent
    / "config"
    / "v16_dual_polarity_policy.json"
)

_CONNECTIVITY = re.compile(r"^[A-Z]{14}$")
_FALSE_EXECUTION_FIELDS = (
    "runtime_effect",
    "production_top1_changed",
    "production_top5_changed",
    "production_score_changed",
    "production_decision_changed",
    "production_confidence_changed",
    "production_compound_count_changed",
    "automatic_activation_allowed",
)


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    """Return a deterministic SHA-256 for a JSON-compatible value."""

    return hashlib.sha256(_canonical(value)).hexdigest()


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _median(values: Sequence[Any]) -> float | None:
    numbers = [number for number in (_finite(value) for value in values) if number is not None]
    return float(statistics.median(numbers)) if numbers else None


def _normalise_polarity(value: Any) -> str:
    text = re.sub(r"[^A-Z]+", "", str(value or "").upper())
    if text in {"POS", "POSITIVE", "PLUS"}:
        return "POS"
    if text in {"NEG", "NEGATIVE", "MINUS"}:
        return "NEG"
    return ""


def _task_id(task: Mapping[str, Any]) -> str:
    return str(task.get("id") or task.get("task_id") or "").strip()


def _task_polarity(
    task: Mapping[str, Any], results: Sequence[Mapping[str, Any]]
) -> str:
    polarity = _normalise_polarity(task.get("ion_mode"))
    if not polarity:
        raise ValueError(
            "Each dual-polarity task must declare task['ion_mode'] as POS or NEG"
        )
    for row in results:
        if not isinstance(row, Mapping):
            continue
        row_polarities = {
            observed
            for key in ("ion_mode", "polarity")
            if (observed := _normalise_polarity(row.get(key)))
        }
        if any(observed != polarity for observed in row_polarities):
            raise ValueError("A result-row polarity conflicts with task['ion_mode']")
    return polarity


def _validate_policy(policy: Mapping[str, Any]) -> dict[str, Any]:
    selected = copy.deepcopy(dict(policy))
    if selected.get("schema_version") != POLICY_SCHEMA_VERSION:
        raise ValueError("V16 dual-polarity policy schema is not locked")
    if selected.get("module_id") != MODULE_ID:
        raise ValueError("V16 dual-polarity module id is not locked")
    if selected.get("algorithm_id") != ALGORITHM_ID:
        raise ValueError("V16 dual-polarity algorithm id is not locked")

    execution = selected.get("execution") or {}
    if execution.get("mode") != "shadow_only":
        raise ValueError("V16 dual-polarity execution must remain shadow-only")
    for field in _FALSE_EXECUTION_FIELDS:
        if execution.get(field) is not False:
            raise ValueError(f"V16 dual-polarity execution invariant changed: {field}")

    pairing = selected.get("task_pairing") or {}
    if pairing.get("explicit_task_ids_required") is not True:
        raise ValueError("V16 dual-polarity pairing must remain explicit")
    if pairing.get("filename_pairing_allowed") is not False:
        raise ValueError("V16 dual-polarity filename pairing is forbidden")
    if set(pairing.get("required_polarities") or []) != {"POS", "NEG"}:
        raise ValueError("V16 dual-polarity task polarities are not locked")

    anchors = selected.get("candidate_independent_anchor") or {}
    if anchors.get("unlinked_feature_neutral_mass_allowed") is not False:
        raise ValueError("Unlinked feature neutral mass cannot become a V16 anchor")
    if anchors.get("polarity_only_default_adduct_allowed") is not False:
        raise ValueError("Polarity-only default adduct inference is forbidden")
    for source in ("iimn_required", "openms_group_required"):
        requirements = anchors.get(source) or {}
        if requirements.get("explicit_adduct_required") is not True:
            raise ValueError("Reliable anchors require an explicit adduct")
        if requirements.get("adduct_sign_must_match_task_ion_mode") is not True:
            raise ValueError("Anchor adduct sign must match task ion mode")
        if _finite(requirements.get("minimum_member_feature_count")) != 2.0:
            raise ValueError("Reliable anchor groups require at least two features")
        if requirements.get("candidate_identity_conflict_affects_membership") is not False:
            raise ValueError("Candidate identity conflict cannot determine graph membership")
    openms = anchors.get("openms_group_required") or {}
    if _finite(openms.get("neutral_mass_tolerance_ppm")) != 10.0:
        raise ValueError("Validated OpenMS groups require a 10 ppm neutral-mass limit")
    if _finite(openms.get("neutral_mass_tolerance_da_floor")) != 0.001:
        raise ValueError("Validated OpenMS groups require a 0.001 Da mass floor")
    if _finite(openms.get("maximum_member_rt_delta_seconds")) != 10.0:
        raise ValueError("Validated OpenMS groups require a 10 second RT limit")
    if _finite(openms.get("minimum_distinct_adduct_count")) != 2.0:
        raise ValueError("Validated OpenMS groups require two distinct adducts")
    if openms.get("pyopenms_manifest_and_feature_xml_hash_binding_required") is not True:
        raise ValueError("Validated OpenMS groups require source-hash binding")
    if openms.get("deconvolution_status_required") != "ok":
        raise ValueError("Validated OpenMS groups require successful deconvolution")

    edge = selected.get("cross_polarity_edge") or {}
    if edge.get("candidate_identity_used_for_edge_construction") is not False:
        raise ValueError("Candidate identity cannot construct cross-polarity edges")
    if edge.get("spectrum_similarity_used_for_edge_construction") is not False:
        raise ValueError("Spectrum similarity cannot construct cross-polarity edges")
    if (_finite(edge.get("neutral_mass_tolerance_ppm")) or 0) <= 0:
        raise ValueError("Neutral-mass ppm tolerance must be positive")
    if (_finite(edge.get("neutral_mass_tolerance_da_floor")) or 0) <= 0:
        raise ValueError("Neutral-mass absolute tolerance must be positive")
    if _finite(edge.get("maximum_rt_delta_seconds")) != 20.0:
        raise ValueError("RT tolerance must remain fixed at 20 seconds")
    if edge.get("pairing_strategy") != "degree_one_unique_only":
        raise ValueError("V16 dual-polarity pairing strategy is not locked")
    if edge.get("one_to_many_action") != "ambiguous_no_merge":
        raise ValueError("One-to-many cross-polarity edges must remain ambiguous")

    support = selected.get("candidate_support") or {}
    for field in (
        "candidate_can_replace_winner",
        "candidate_can_receive_2d_connectivity_decision",
        "candidate_can_increase_compound_count",
        "same_name_is_connectivity_evidence",
        "same_formula_is_connectivity_evidence",
    ):
        if support.get(field) is not False:
            raise ValueError(f"V16 candidate-support invariant changed: {field}")
    if support.get("stable_2d_connectivity_required") is not True:
        raise ValueError("Stable 2-D connectivity must remain required")
    for field in (
        "candidate_attachment_occurs_after_graph_freeze",
        "potential_contaminants_excluded_from_display_support",
        "explicit_contaminant_or_exogenous_roles_excluded_from_display_support",
        "excluded_candidates_retained_for_quality_control",
    ):
        if support.get(field) is not True:
            raise ValueError(f"V16 candidate-support invariant changed: {field}")
    if support.get("flagship_binding_field") != "member_result_indices":
        raise ValueError("Flagship candidate support must bind by member_result_indices")

    auxiliary = selected.get("cross_polarity_ms2_similarity") or {}
    if auxiliary.get("state") != "auxiliary_only":
        raise ValueError("Cross-polarity MS/MS similarity must remain auxiliary")
    for field in (
        "can_create_candidate_independent_edge",
        "can_create_dual_connectivity_support",
        "can_change_ranking",
        "score_is_probability",
    ):
        if auxiliary.get(field) is not False:
            raise ValueError(f"V16 auxiliary-similarity invariant changed: {field}")
    return selected


def load_v16_dual_polarity_policy(path: str | Path | None = None) -> dict[str, Any]:
    selected = Path(path) if path is not None else DEFAULT_POLICY_PATH
    return _validate_policy(json.loads(selected.read_text(encoding="utf-8")))


def _task_binding_projection(
    task: Mapping[str, Any], *, task_id: str, polarity: str
) -> dict[str, Any]:
    workflow = task.get("untargeted_workflow") or {}
    if not isinstance(workflow, Mapping):
        workflow = {}
    return {
        "task_id": task_id,
        "status": str(task.get("status") or ""),
        "polarity": polarity,
        "result_run_id": str(task.get("result_run_id") or ""),
        "feature_xml": str(
            workflow.get("feature_xml") or workflow.get("output_features") or ""
        ),
    }


def _specific_text(value: Any) -> str:
    return str(value or "").strip()


def _row_rt_minutes(row: Mapping[str, Any]) -> float | None:
    feature_rt = _finite(row.get("lcms_feature_rt_minutes"))
    if feature_rt is not None:
        return feature_rt
    value = _finite(row.get("retention_time"))
    if value is None:
        return None
    unit = str(row.get("retention_time_unit") or "").strip().lower()
    if unit in {"s", "sec", "second", "seconds"}:
        return value / 60.0
    return value


def _normalise_formula(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "").strip())


def _stable_connectivity(candidate: Mapping[str, Any]) -> str:
    nested = candidate.get("candidate")
    values = []
    for source in (candidate, nested if isinstance(nested, Mapping) else {}):
        values.extend(
            source.get(key)
            for key in (
                "connectivity_identity",
                "candidate_connectivity",
                "identity_key",
                "inchikey",
            )
        )
    for raw in values:
        text = str(raw or "").strip().upper()
        for prefix in ("STD2D::", "STD2D:", "IK14::", "IK14:", "IK::"):
            if text.startswith(prefix):
                text = text[len(prefix) :]
                break
        first = text.split("-", 1)[0]
        if _CONNECTIVITY.fullmatch(first):
            return first
    return ""


def _candidate_material(candidate: Mapping[str, Any]) -> dict[str, Any]:
    nested = candidate.get("candidate")
    if isinstance(nested, Mapping):
        merged = dict(nested)
        merged.update(candidate)
        merged.pop("candidate", None)
        return merged
    return dict(candidate)


def _candidate_exclusion_reasons(material: Mapping[str, Any]) -> list[str]:
    reasons = []
    if material.get("potential_contaminant") is True:
        reasons.append("potential_contaminant=true")
    for field in ("compound_role", "candidate_role", "role"):
        value = _specific_text(material.get(field))
        normalised = re.sub(r"[^a-z]+", "_", value.lower()).strip("_")
        if value and ("contaminant" in normalised or "exogenous" in normalised):
            reasons.append(f"{field}={value}")
    return sorted(set(reasons))


def _candidate_provenance(
    candidate: Mapping[str, Any], *, origin: str, row_index: int, rank: int | None
) -> dict[str, Any]:
    material = _candidate_material(candidate)
    connectivity = _stable_connectivity(material)
    source = _specific_text(
        material.get("dataset_source")
        or material.get("candidate_source")
        or material.get("source")
    )
    accession = _specific_text(
        material.get("lib_source_accession")
        or material.get("source_accession")
        or material.get("lib_spectrum_id")
    )
    return {
        "connectivity": connectivity or None,
        "formula": _normalise_formula(
            material.get("formula") or material.get("molecular_formula")
        ),
        "display_name": _specific_text(
            material.get("compound_name") or material.get("name")
        ),
        "origin": origin,
        "row_index": int(row_index),
        "rank": int(rank) if rank is not None else None,
        "dataset_source": source,
        "source_accession": accession,
        "support_only": origin == "top40_support_only",
        "quality_control_exclusion_reasons": _candidate_exclusion_reasons(material),
    }


def _row_candidates(row: Mapping[str, Any], row_index: int) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    def add(candidate: Any, *, origin: str, rank: int | None) -> None:
        if not isinstance(candidate, Mapping):
            return
        item = _candidate_provenance(
            candidate, origin=origin, row_index=row_index, rank=rank
        )
        records.append(item)

    match = row.get("match")
    if isinstance(match, Mapping):
        add(match, origin="production_top1", rank=1)
        for rank, candidate in enumerate(match.get("top_candidates") or [], start=2):
            add(candidate, origin="native_candidate", rank=rank)
    winner = row.get("winner")
    if isinstance(winner, Mapping):
        add(winner, origin="production_top1", rank=1)
    for rank, candidate in enumerate(row.get("candidates") or [], start=1):
        add(candidate, origin="native_candidate", rank=rank)

    supports = []
    for container in (row, match if isinstance(match, Mapping) else {}):
        support = container.get("candidate_recall_support")
        if isinstance(support, Mapping):
            supports.extend(support.get("expanded_candidates") or [])
    for rank, candidate in enumerate(supports, start=1):
        add(candidate, origin="top40_support_only", rank=rank)
    return records


def _flagship_candidate_records_by_row(
    artifact: Mapping[str, Any],
) -> dict[int, list[dict[str, Any]]]:
    by_row: dict[int, list[dict[str, Any]]] = defaultdict(list)
    entities = artifact.get("entities") or []
    if isinstance(entities, (str, bytes)) or not isinstance(entities, Sequence):
        return by_row
    for entity in entities:
        if not isinstance(entity, Mapping):
            continue
        member_indices = sorted(
            {
                int(float(value))
                for value in entity.get("member_result_indices") or []
                if _finite(value) is not None and int(float(value)) >= 0
            }
        )
        if not member_indices:
            continue
        binding_index = member_indices[0]
        records: list[dict[str, Any]] = []

        def add(candidate: Any, *, origin: str, rank: int | None) -> None:
            if not isinstance(candidate, Mapping):
                return
            item = _candidate_provenance(
                candidate,
                origin=origin,
                row_index=binding_index,
                rank=rank,
            )
            item["flagship_entity_id"] = _specific_text(entity.get("entity_id"))
            item["member_result_indices"] = member_indices
            records.append(item)

        add(entity.get("winner"), origin="production_top1", rank=1)
        if entity.get("winner") is None:
            add(
                entity.get("leading_candidate_when_abstained"),
                origin="abstained_leading_candidate",
                rank=1,
            )
        for rank, candidate in enumerate(entity.get("candidates") or [], start=1):
            add(candidate, origin="native_candidate", rank=rank)
        recall = entity.get("candidate_recall_support")
        if isinstance(recall, Mapping):
            for rank, candidate in enumerate(
                recall.get("expanded_candidates") or [], start=1
            ):
                add(candidate, origin="top40_support_only", rank=rank)
        for member_index in member_indices:
            by_row[member_index].extend(copy.deepcopy(records))
    return by_row


def _candidate_records_for_anchor(
    anchor: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    *,
    flagship_by_row: Mapping[int, Sequence[Mapping[str, Any]]] | None,
) -> list[dict[str, Any]]:
    records = []
    for row_index in anchor.get("source_row_indices") or []:
        index = int(row_index)
        if flagship_by_row is not None:
            records.extend(copy.deepcopy(list(flagship_by_row.get(index) or [])))
        elif 0 <= index < len(rows) and isinstance(rows[index], Mapping):
            records.extend(_row_candidates(rows[index], index))
    # A flagship entity can bind the same record through several member rows.
    # Canonical de-duplication keeps the payload deterministic.
    unique: dict[bytes, dict[str, Any]] = {}
    for record in records:
        unique[_canonical(record)] = dict(record)
    return [unique[key] for key in sorted(unique)]


def _candidate_identity_conflict_observed(
    anchor: Mapping[str, Any], rows: Sequence[Mapping[str, Any]]
) -> bool:
    for raw_index in anchor.get("source_row_indices") or []:
        index = int(raw_index)
        if not (0 <= index < len(rows)) or not isinstance(rows[index], Mapping):
            continue
        row = rows[index]
        for field in ("ion_identity_evidence", "lcms_adduct_group_evidence"):
            evidence = row.get(field)
            if isinstance(evidence, Mapping) and evidence.get("identity_conflict") is True:
                return True
    return False


def _candidate_summary(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        connectivity = str(record.get("connectivity") or "")
        if connectivity and not record.get("quality_control_exclusion_reasons"):
            grouped[connectivity].append(record)
    output = []
    for connectivity in sorted(grouped):
        items = grouped[connectivity]
        origins = sorted({str(item.get("origin") or "") for item in items})
        formulas = sorted({str(item.get("formula") or "") for item in items if item.get("formula")})
        display_names = sorted(
            {
                str(item.get("display_name") or "")
                for item in items
                if item.get("display_name")
            }
        )
        provenance = sorted(
            {
                (
                    str(item.get("origin") or ""),
                    int(item.get("row_index") or 0),
                    int(item["rank"]) if item.get("rank") is not None else None,
                    str(item.get("dataset_source") or ""),
                    str(item.get("source_accession") or ""),
                    str(item.get("flagship_entity_id") or ""),
                    tuple(int(value) for value in item.get("member_result_indices") or []),
                )
                for item in items
            },
            key=lambda item: (
                item[0], item[1], item[2] if item[2] is not None else math.inf,
                item[3], item[4], item[5], item[6],
            ),
        )
        output.append(
            {
                "connectivity": connectivity,
                "formulas": formulas,
                "display_names": display_names,
                "display_names_used_for_structure_identity": False,
                "origins": origins,
                "has_native_candidate": any(
                    origin != "top40_support_only" for origin in origins
                ),
                "has_top40_support_only_candidate": (
                    "top40_support_only" in origins
                ),
                "provenance": [
                    {
                        "origin": item[0],
                        "row_index": item[1],
                        "rank": item[2],
                        "dataset_source": item[3],
                        "source_accession": item[4],
                        "flagship_entity_id": item[5] or None,
                        "member_result_indices": list(item[6]),
                    }
                    for item in provenance
                ],
            }
        )
    return output


def _quality_control_excluded_candidates(
    records: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    excluded = []
    seen = set()
    for record in records:
        reasons = tuple(record.get("quality_control_exclusion_reasons") or [])
        if not reasons:
            continue
        material = {
            "connectivity": record.get("connectivity"),
            "formula": str(record.get("formula") or ""),
            "display_name": str(record.get("display_name") or ""),
            "origin": str(record.get("origin") or ""),
            "row_index": int(record.get("row_index") or 0),
            "rank": (
                int(record["rank"]) if record.get("rank") is not None else None
            ),
            "dataset_source": str(record.get("dataset_source") or ""),
            "source_accession": str(record.get("source_accession") or ""),
            "exclusion_reasons": list(reasons),
            "excluded_from_display_candidate_support": True,
            "retained_for_quality_control_only": True,
            "display_name_used_for_structure_identity": False,
        }
        key = _canonical(material)
        if key not in seen:
            seen.add(key)
            excluded.append(material)
    return sorted(
        excluded,
        key=lambda item: (
            item["origin"],
            item["row_index"],
            item["rank"] if item["rank"] is not None else math.inf,
            item["connectivity"] or "",
            item["source_accession"],
        ),
    )


def _adduct_polarity(value: Any) -> str:
    text = re.sub(r"\s+", "", str(value or ""))
    if text.endswith("+"):
        return "POS"
    if text.endswith("-"):
        return "NEG"
    return ""


def _member_feature_ids(evidence: Mapping[str, Any]) -> list[str]:
    values = evidence.get("member_feature_ids") or []
    if isinstance(values, (str, bytes)):
        values = [values]
    if not isinstance(values, Sequence):
        return []
    return sorted({_specific_text(value) for value in values if _specific_text(value)})


def _row_anchor_evidence(
    row: Mapping[str, Any], task_polarity: str
) -> dict[str, Any] | None:
    iimn = row.get("ion_identity_evidence")
    if not isinstance(iimn, Mapping):
        iimn = {}
    iimn_group = _specific_text(
        row.get("ion_identity_network_id") or iimn.get("network_id")
    )
    iimn_neutral = _finite(
        row.get("ion_identity_neutral_mass")
        if row.get("ion_identity_neutral_mass") is not None
        else iimn.get("neutral_mass")
    )
    iimn_adduct = _specific_text(
        row.get("ion_identity") or iimn.get("ion_identity")
    )
    iimn_valid = bool(
        iimn_group
        and iimn_neutral is not None
        and iimn_neutral > 0
        and str(row.get("ion_identity_status") or "").lower() == "networked"
        and str(iimn.get("status") or "").lower() == "networked"
        and iimn.get("runtime_influence") is True
        and iimn_adduct
        and _adduct_polarity(iimn_adduct) == task_polarity
    )
    if iimn_valid:
        return {
            "source": "reliable_iimn",
            "group_id": iimn_group,
            "neutral_mass": iimn_neutral,
            "adduct": iimn_adduct,
            "evidence_member_feature_ids": _member_feature_ids(iimn),
            "evidence_member_feature_count": _finite(
                iimn.get("member_feature_count")
            ),
        }

    group_evidence = row.get("lcms_adduct_group_evidence")
    if not isinstance(group_evidence, Mapping):
        group_evidence = {}
    openms_group = _specific_text(row.get("lcms_adduct_group_id"))
    openms_neutral = _finite(row.get("lcms_neutral_mass"))
    openms_adduct = _specific_text(row.get("lcms_adduct"))
    openms_valid = bool(
        openms_group
        and openms_neutral is not None
        and openms_neutral > 0
        and str(group_evidence.get("status") or "").lower() == "validated"
        and group_evidence.get("runtime_influence") is True
        and group_evidence.get("neutral_mass_consistent") is True
        and openms_adduct
        and _adduct_polarity(openms_adduct) == task_polarity
    )
    if openms_valid:
        return {
            "source": "validated_openms_adduct_group",
            "group_id": openms_group,
            "neutral_mass": openms_neutral,
            "adduct": openms_adduct,
            "evidence_member_feature_ids": _member_feature_ids(group_evidence),
            "evidence_member_feature_count": _finite(
                group_evidence.get("member_feature_count")
            ),
        }
    return None


def _candidate_independent_group_conflicts(
    row: Mapping[str, Any],
) -> set[tuple[str, str]]:
    conflicts: set[tuple[str, str]] = set()
    group_evidence = row.get("lcms_adduct_group_evidence")
    if not isinstance(group_evidence, Mapping):
        group_evidence = {}
    openms_group = _specific_text(row.get("lcms_adduct_group_id"))
    if openms_group and group_evidence.get("neutral_mass_consistent") is False:
        conflicts.add(("validated_openms_adduct_group", openms_group))
    return conflicts


def _mass_tolerance(neutral_mass: float, policy: Mapping[str, Any]) -> float:
    edge = policy["cross_polarity_edge"]
    return max(
        float(edge["neutral_mass_tolerance_da_floor"]),
        abs(neutral_mass) * float(edge["neutral_mass_tolerance_ppm"]) * 1e-6,
    )


def _extract_anchors(
    task_id: str,
    polarity: str,
    results: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    conflicts: set[tuple[str, str]] = set()
    rejected_rows = 0
    for row_index, row in enumerate(results):
        if not isinstance(row, Mapping):
            rejected_rows += 1
            continue
        conflicts |= _candidate_independent_group_conflicts(row)
        evidence = _row_anchor_evidence(row, polarity)
        if evidence is None:
            rejected_rows += 1
            continue
        feature_id = _specific_text(
            row.get("lcms_feature_id")
            or row.get("precursor_feature_id")
            or row.get("spectrum_id")
            or f"row-{row_index}"
        )
        groups[(evidence["source"], evidence["group_id"])].append(
            {
                "row_index": row_index,
                "feature_id": feature_id,
                "neutral_mass": evidence["neutral_mass"],
                "retention_time_min": _row_rt_minutes(row),
                "precursor_mz": _finite(
                    row.get("lcms_feature_mz")
                    if row.get("lcms_feature_mz") is not None
                    else row.get("precursor_mz")
                ),
                "adduct": evidence["adduct"],
                "evidence_member_feature_ids": evidence[
                    "evidence_member_feature_ids"
                ],
                "evidence_member_feature_count": evidence[
                    "evidence_member_feature_count"
                ],
            }
        )

    anchors = []
    rejected_conflict_groups = 0
    rejected_inconsistent_groups = 0
    rejected_missing_rt_groups = 0
    rejected_insufficient_member_groups = 0
    for (source, group_id), entries in sorted(groups.items()):
        if (source, group_id) in conflicts:
            rejected_conflict_groups += 1
            continue
        neutral_values = [float(entry["neutral_mass"]) for entry in entries]
        centre = float(statistics.median(neutral_values))
        if max(neutral_values) - min(neutral_values) > _mass_tolerance(centre, policy):
            rejected_inconsistent_groups += 1
            continue
        retention_time = _median(
            [entry.get("retention_time_min") for entry in entries]
        )
        if retention_time is None:
            rejected_missing_rt_groups += 1
            continue
        feature_ids = sorted(
            {str(entry["feature_id"]) for entry in entries}
            | {
                str(feature_id)
                for entry in entries
                for feature_id in entry.get("evidence_member_feature_ids") or []
            }
        )
        declared_member_count = max(
            [
                int(count)
                for entry in entries
                if (count := _finite(entry.get("evidence_member_feature_count")))
                is not None
            ]
            or [0]
        )
        verified_member_count = max(len(feature_ids), declared_member_count)
        if verified_member_count < 2:
            rejected_insufficient_member_groups += 1
            continue
        anchor_material = {
            "task_id": task_id,
            "polarity": polarity,
            "source": source,
            "group_id": group_id,
        }
        anchors.append(
            {
                "anchor_id": f"DPA-{semantic_sha256(anchor_material)[:16]}",
                **anchor_material,
                "neutral_mass": round(centre, 8),
                "retention_time_min": round(retention_time, 8),
                "precursor_mz": (
                    round(_median([entry.get("precursor_mz") for entry in entries]), 8)
                    if _median([entry.get("precursor_mz") for entry in entries])
                    is not None
                    else None
                ),
                "feature_ids": feature_ids,
                "adducts": sorted(
                    {str(entry["adduct"]) for entry in entries if entry.get("adduct")}
                ),
                "source_row_indices": sorted(
                    {int(entry["row_index"]) for entry in entries}
                ),
                "verified_member_feature_count": verified_member_count,
                "candidate_independent": True,
            }
        )
    anchors.sort(key=lambda item: item["anchor_id"])
    return anchors, {
        "input_rows": len(results),
        "rows_without_reliable_anchor": rejected_rows,
        "rejected_conflict_groups": rejected_conflict_groups,
        "rejected_inconsistent_neutral_mass_groups": rejected_inconsistent_groups,
        "rejected_missing_rt_groups": rejected_missing_rt_groups,
        "rejected_insufficient_member_groups": rejected_insufficient_member_groups,
        "reliable_anchor_groups": len(anchors),
    }


def rt_alignment_semantic_sha256(alignment: Mapping[str, Any]) -> str:
    material = {
        key: value
        for key, value in alignment.items()
        if key not in {"semantic_sha256", "artifact_sha256"}
    }
    return semantic_sha256(material)


def _alignment_statuses_ok(alignment: Mapping[str, Any], task_ids: set[str]) -> bool:
    statuses = alignment.get("alignment_status")
    if statuses is None:
        return True
    if isinstance(statuses, Mapping):
        return all(
            str(statuses.get(task_id) or "").lower()
            in {"ok", "aligned", "reference"}
            for task_id in task_ids
        )
    if isinstance(statuses, Sequence) and not isinstance(statuses, (str, bytes)):
        values = []
        for item in statuses:
            if isinstance(item, Mapping):
                values.append(str(item.get("status") or "").lower())
            else:
                values.append(str(item or "").lower())
        return bool(values) and all(value in {"ok", "aligned", "reference"} for value in values)
    return False


def _apply_rt_alignment(
    anchors_by_polarity: dict[str, list[dict[str, Any]]],
    alignment: Mapping[str, Any] | None,
    *,
    source_hashes: Mapping[str, str],
    task_hashes: Mapping[str, str],
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    output = copy.deepcopy(anchors_by_polarity)
    if alignment is None:
        for anchors in output.values():
            for anchor in anchors:
                anchor["aligned_retention_time_min"] = anchor["retention_time_min"]
                anchor["rt_basis"] = "raw_task_rt"
        return output, {
            "status": "not_supplied_raw_rt_used",
            "semantic_sha256": None,
            "anchors_excluded_missing_mapping": 0,
        }

    task_ids = set(source_hashes)
    if str(alignment.get("status") or "").lower() != "ok":
        raise ValueError("RT alignment must have status=ok")
    if set(str(item) for item in alignment.get("task_ids") or []) != task_ids:
        raise ValueError("RT alignment is bound to different task IDs")
    if dict(alignment.get("source_results_sha256") or {}) != dict(source_hashes):
        raise ValueError("RT alignment is bound to different source results")
    if dict(alignment.get("source_task_metadata_sha256") or {}) != dict(task_hashes):
        raise ValueError("RT alignment is bound to different task metadata")
    if not _alignment_statuses_ok(alignment, task_ids):
        raise ValueError("RT alignment contains a failed or partial task alignment")
    observed_hash = rt_alignment_semantic_sha256(alignment)
    declared_hash = str(
        alignment.get("semantic_sha256") or alignment.get("artifact_sha256") or ""
    )
    if declared_hash and declared_hash != observed_hash:
        raise ValueError("RT alignment semantic hash does not match its content")
    mapping = alignment.get("aligned_rt_minutes")
    if not isinstance(mapping, Mapping):
        raise ValueError("RT alignment must provide aligned_rt_minutes by task")

    excluded = 0
    for polarity in ("POS", "NEG"):
        retained = []
        for anchor in output[polarity]:
            task_map = mapping.get(anchor["task_id"])
            if not isinstance(task_map, Mapping):
                excluded += 1
                continue
            values = []
            for key in (
                anchor["anchor_id"],
                anchor["group_id"],
                *anchor["feature_ids"],
            ):
                value = _finite(task_map.get(key))
                if value is not None:
                    values.append(value)
            if not values:
                excluded += 1
                continue
            anchor["aligned_retention_time_min"] = round(
                float(statistics.median(values)), 8
            )
            anchor["rt_basis"] = "hash_bound_rt_alignment"
            retained.append(anchor)
        output[polarity] = retained
    return output, {
        "status": "applied",
        "semantic_sha256": observed_hash,
        "anchors_excluded_missing_mapping": excluded,
    }


def _edge_candidate(
    positive: Mapping[str, Any],
    negative: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> dict[str, Any] | None:
    pos_mass = float(positive["neutral_mass"])
    neg_mass = float(negative["neutral_mass"])
    neutral = (pos_mass + neg_mass) / 2.0
    mass_error = abs(pos_mass - neg_mass)
    mass_tolerance = _mass_tolerance(neutral, policy)
    rt_delta = abs(
        float(positive["aligned_retention_time_min"])
        - float(negative["aligned_retention_time_min"])
    )
    rt_tolerance_seconds = float(
        policy["cross_polarity_edge"]["maximum_rt_delta_seconds"]
    )
    rt_tolerance = rt_tolerance_seconds / 60.0
    # Decimal thresholds in the policy are inclusive.  Absorb only binary
    # floating-point representation noise at the exact boundary; this does
    # not widen either scientific tolerance.
    if (
        mass_error - mass_tolerance > 1e-12
        or rt_delta - rt_tolerance > 1e-12
    ):
        return None
    distance = mass_error / mass_tolerance + rt_delta / rt_tolerance
    return {
        "positive_anchor_id": positive["anchor_id"],
        "negative_anchor_id": negative["anchor_id"],
        "distance": distance,
        "neutral_mass": neutral,
        "neutral_mass_error_da": mass_error,
        "neutral_mass_error_ppm": mass_error / neutral * 1e6,
        "neutral_mass_tolerance_da": mass_tolerance,
        "rt_delta_min": rt_delta,
        "rt_delta_seconds": rt_delta * 60.0,
        "rt_tolerance_min": rt_tolerance,
        "rt_tolerance_seconds": rt_tolerance_seconds,
    }


def _classify_edges(
    positive: Sequence[Mapping[str, Any]],
    negative: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    candidates = [
        edge
        for pos in positive
        for neg in negative
        if (edge := _edge_candidate(pos, neg, policy)) is not None
    ]
    positive_degree: dict[str, int] = defaultdict(int)
    negative_degree: dict[str, int] = defaultdict(int)
    for edge in candidates:
        positive_degree[str(edge["positive_anchor_id"])] += 1
        negative_degree[str(edge["negative_anchor_id"])] += 1

    compatible = []
    accepted = []
    ambiguous = []
    for raw in candidates:
        edge = dict(raw)
        pos_id = str(edge["positive_anchor_id"])
        neg_id = str(edge["negative_anchor_id"])
        edge["positive_anchor_degree"] = positive_degree[pos_id]
        edge["negative_anchor_degree"] = negative_degree[neg_id]
        if positive_degree[pos_id] == 1 and negative_degree[neg_id] == 1:
            edge["pairing_status"] = "accepted_unique"
            accepted.append(edge)
        else:
            edge["pairing_status"] = "ambiguous_no_merge"
            ambiguous.append(edge)
        compatible.append(edge)
    key = lambda item: (
        str(item["positive_anchor_id"]),
        str(item["negative_anchor_id"]),
    )
    compatible.sort(key=key)
    accepted.sort(key=key)
    ambiguous.sort(key=key)
    return compatible, accepted, ambiguous


def _candidate_support_for_entity(
    positive: Mapping[str, Any], negative: Mapping[str, Any]
) -> list[dict[str, Any]]:
    sides = {
        "POS": {
            item["connectivity"]: item
            for item in positive.get("candidate_support") or []
        },
        "NEG": {
            item["connectivity"]: item
            for item in negative.get("candidate_support") or []
        },
    }
    output = []
    for connectivity in sorted(set(sides["POS"]) | set(sides["NEG"])):
        pos = sides["POS"].get(connectivity)
        neg = sides["NEG"].get(connectivity)
        pos_formulas = set((pos or {}).get("formulas") or [])
        neg_formulas = set((neg or {}).get("formulas") or [])
        formula_conflict = bool(
            pos_formulas and neg_formulas and pos_formulas.isdisjoint(neg_formulas)
        )
        dual = bool(pos and neg and not formula_conflict)
        if formula_conflict:
            status = "cross_polarity_formula_conflict"
            scope = "conflicting_candidate_support_only"
        elif dual:
            status = "dual_polarity_support"
            scope = "dual_polarity_candidate_support_only"
        else:
            status = "single_polarity_support"
            scope = "candidate_hypothesis_support_only"
        polarity_provenance = {}
        for polarity, item in (("POS", pos), ("NEG", neg)):
            if item is not None:
                polarity_provenance[polarity] = copy.deepcopy(item)
        output.append(
            {
                "connectivity": connectivity,
                "status": status,
                "support_scope": scope,
                "dual_polarity_support": dual,
                "formula_conflict": formula_conflict,
                "polarities_observed": sorted(polarity_provenance),
                "polarity_provenance": polarity_provenance,
                "can_replace_winner": False,
                "can_receive_2d_connectivity_decision": False,
                "compound_count_increment": 0,
                "candidate_name_used_for_inference": False,
            }
        )
    return output


def _auxiliary_map(
    records: Sequence[Mapping[str, Any]] | None,
) -> dict[tuple[str, str], dict[str, Any]]:
    output = {}
    for record in records or []:
        if not isinstance(record, Mapping):
            continue
        pos_group = _specific_text(record.get("positive_group_id"))
        neg_group = _specific_text(record.get("negative_group_id"))
        if not pos_group or not neg_group:
            continue
        scores = {
            key: round(value, 8)
            for key in ("modified_cosine", "neutral_losses_cosine", "ms2deepscore")
            if (value := _finite(record.get(key))) is not None
        }
        output[(pos_group, neg_group)] = {
            "state": "auxiliary_only",
            "scores": scores,
            "used_for": ["edge_annotation", "analogue_hint"],
            "runtime_influence": False,
            "can_create_candidate_independent_edge": False,
            "can_create_dual_connectivity_support": False,
            "can_change_ranking": False,
            "score_is_probability": False,
        }
    return output


def _graph_anchor_projection(anchor: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "anchor_id": anchor["anchor_id"],
        "task_id": anchor["task_id"],
        "polarity": anchor["polarity"],
        "source": anchor["source"],
        "group_id": anchor["group_id"],
        "neutral_mass": anchor["neutral_mass"],
        "retention_time_min": anchor["retention_time_min"],
        "aligned_retention_time_min": anchor["aligned_retention_time_min"],
        "rt_basis": anchor["rt_basis"],
        "precursor_mz": anchor["precursor_mz"],
        "feature_ids": list(anchor["feature_ids"]),
        "source_row_indices": list(anchor["source_row_indices"]),
        "adducts": list(anchor["adducts"]),
        "verified_member_feature_count": anchor["verified_member_feature_count"],
        "candidate_independent": True,
    }


def _graph_edge_projection(edge: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "positive_anchor_id": edge["positive_anchor_id"],
        "negative_anchor_id": edge["negative_anchor_id"],
        "pairing_status": edge["pairing_status"],
        "positive_anchor_degree": int(edge["positive_anchor_degree"]),
        "negative_anchor_degree": int(edge["negative_anchor_degree"]),
        "neutral_mass": round(float(edge["neutral_mass"]), 8),
        "neutral_mass_error_da": round(float(edge["neutral_mass_error_da"]), 8),
        "neutral_mass_error_ppm": round(float(edge["neutral_mass_error_ppm"]), 6),
        "neutral_mass_tolerance_da": round(
            float(edge["neutral_mass_tolerance_da"]), 8
        ),
        "rt_delta_min": round(float(edge["rt_delta_min"]), 8),
        "rt_delta_seconds": round(float(edge["rt_delta_seconds"]), 6),
        "rt_tolerance_min": round(float(edge["rt_tolerance_min"]), 8),
        "rt_tolerance_seconds": round(float(edge["rt_tolerance_seconds"]), 6),
    }


def _frozen_candidate_independent_graph(
    anchors_by_polarity: Mapping[str, Sequence[Mapping[str, Any]]],
    compatible_edges: Sequence[Mapping[str, Any]],
    accepted_edges: Sequence[Mapping[str, Any]],
    ambiguous_edges: Sequence[Mapping[str, Any]],
    *,
    task_pair: Mapping[str, str],
    policy_sha256: str,
) -> dict[str, Any]:
    return {
        "graph_schema_version": "mssystem.v16.dual-polarity-candidate-independent-graph/1.0.0",
        "task_pair": {polarity: task_pair[polarity] for polarity in ("POS", "NEG")},
        "policy_sha256": policy_sha256,
        "anchors": {
            polarity: [
                _graph_anchor_projection(anchor)
                for anchor in anchors_by_polarity[polarity]
            ]
            for polarity in ("POS", "NEG")
        },
        "compatible_edges": [
            _graph_edge_projection(edge) for edge in compatible_edges
        ],
        "accepted_unique_edges": [
            _graph_edge_projection(edge) for edge in accepted_edges
        ],
        "ambiguous_no_merge_edges": [
            _graph_edge_projection(edge) for edge in ambiguous_edges
        ],
        "candidate_fields_present": False,
        "spectrum_similarity_present": False,
    }


def sidecar_semantic_sha256(payload: Mapping[str, Any]) -> str:
    material = copy.deepcopy(dict(payload))
    metadata = material.get("metadata")
    if isinstance(metadata, dict):
        metadata.pop("sidecar_semantic_sha256", None)
    return semantic_sha256(material)


def build_v16_dual_polarity_shadow(
    task_a: Mapping[str, Any],
    results_a: Sequence[Mapping[str, Any]],
    task_b: Mapping[str, Any],
    results_b: Sequence[Mapping[str, Any]],
    *,
    rt_alignment: Mapping[str, Any] | None = None,
    auxiliary_spectral_similarity: Sequence[Mapping[str, Any]] | None = None,
    candidate_support_by_task: Mapping[str, Mapping[str, Any]] | None = None,
    policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a deterministic dual-polarity shadow sidecar.

    ``task_a`` and ``task_b`` are explicit caller selections; filenames are
    never inspected.  Input ordering is normalized by polarity.
    """

    if not isinstance(task_a, Mapping) or not isinstance(task_b, Mapping):
        raise TypeError("Dual-polarity task metadata must be mappings")
    if isinstance(results_a, (str, bytes)) or isinstance(results_b, (str, bytes)):
        raise TypeError("Dual-polarity results must be row sequences")
    if candidate_support_by_task is not None and not isinstance(
        candidate_support_by_task, Mapping
    ):
        raise TypeError("candidate_support_by_task must map task IDs to artifacts")
    selected_policy = _validate_policy(policy) if policy is not None else load_v16_dual_polarity_policy()
    policy_hash = semantic_sha256(selected_policy)

    task_a_id, task_b_id = _task_id(task_a), _task_id(task_b)
    if not task_a_id or not task_b_id or task_a_id == task_b_id:
        raise ValueError("Dual-polarity analysis requires two distinct explicit task IDs")
    if str(task_a.get("status") or "") != "completed" or str(
        task_b.get("status") or ""
    ) != "completed":
        raise ValueError("Dual-polarity analysis requires two completed tasks")

    rows_a = list(results_a)
    rows_b = list(results_b)
    polarity_a = _task_polarity(task_a, rows_a)
    polarity_b = _task_polarity(task_b, rows_b)
    if {polarity_a, polarity_b} != {"POS", "NEG"}:
        raise ValueError("Dual-polarity analysis requires exactly one POS and one NEG task")

    supplied = {
        polarity_a: (task_a_id, task_a, rows_a),
        polarity_b: (task_b_id, task_b, rows_b),
    }
    source_hashes_before = {
        task_id: semantic_sha256(rows)
        for task_id, _, rows in supplied.values()
    }
    candidate_input = dict(candidate_support_by_task or {})
    candidate_support_hash_before = semantic_sha256(candidate_input)
    task_hashes = {
        task_id: semantic_sha256(
            _task_binding_projection(task, task_id=task_id, polarity=polarity)
        )
        for polarity, (task_id, task, _) in supplied.items()
    }

    anchors_by_polarity = {}
    extraction = {}
    for polarity in ("POS", "NEG"):
        task_id, _, rows = supplied[polarity]
        anchors, summary = _extract_anchors(
            task_id, polarity, rows, selected_policy
        )
        anchors_by_polarity[polarity] = anchors
        extraction[polarity] = summary

    aligned_anchors, alignment_summary = _apply_rt_alignment(
        anchors_by_polarity,
        rt_alignment,
        source_hashes=source_hashes_before,
        task_hashes=task_hashes,
    )
    compatible_edges, accepted_edges, ambiguous_edges = _classify_edges(
        aligned_anchors["POS"], aligned_anchors["NEG"], selected_policy
    )
    anchors_by_id = {
        anchor["anchor_id"]: anchor
        for anchors in aligned_anchors.values()
        for anchor in anchors
    }
    auxiliary = _auxiliary_map(auxiliary_spectral_similarity)
    task_pair = {
        "POS": supplied["POS"][0],
        "NEG": supplied["NEG"][0],
    }
    frozen_graph = _frozen_candidate_independent_graph(
        aligned_anchors,
        compatible_edges,
        accepted_edges,
        ambiguous_edges,
        task_pair=task_pair,
        policy_sha256=policy_hash,
    )
    frozen_graph_hash_before = semantic_sha256(frozen_graph)

    flagship_by_task: dict[
        str, Mapping[int, Sequence[Mapping[str, Any]]] | None
    ] = {}
    for polarity in ("POS", "NEG"):
        task_id = supplied[polarity][0]
        if candidate_support_by_task is None:
            flagship_by_task[task_id] = None
            continue
        artifact = candidate_input.get(task_id)
        if artifact is None:
            flagship_by_task[task_id] = {}
        elif not isinstance(artifact, Mapping):
            raise TypeError(
                f"Candidate-support artifact for task {task_id!r} must be a mapping"
            )
        else:
            flagship_by_task[task_id] = _flagship_candidate_records_by_row(artifact)

    entities = []
    for edge in accepted_edges:
        positive = anchors_by_id[edge["positive_anchor_id"]]
        negative = anchors_by_id[edge["negative_anchor_id"]]
        entity_material = {
            "positive_task_id": positive["task_id"],
            "positive_group_id": positive["group_id"],
            "negative_task_id": negative["task_id"],
            "negative_group_id": negative["group_id"],
        }
        entity_id = f"DPN-{semantic_sha256(entity_material)[:16]}"
        candidate_records = {}
        candidate_summaries = {}
        excluded_by_polarity = {}
        for side_polarity, anchor in (("POS", positive), ("NEG", negative)):
            task_id = str(anchor["task_id"])
            rows = supplied[side_polarity][2]
            records = _candidate_records_for_anchor(
                anchor,
                rows,
                flagship_by_row=flagship_by_task[task_id],
            )
            candidate_records[side_polarity] = records
            candidate_summaries[side_polarity] = _candidate_summary(records)
            excluded_by_polarity[side_polarity] = (
                _quality_control_excluded_candidates(records)
            )
        candidate_support = _candidate_support_for_entity(
            {"candidate_support": candidate_summaries["POS"]},
            {"candidate_support": candidate_summaries["NEG"]},
        )
        excluded_candidates = []
        for side_polarity in ("POS", "NEG"):
            for record in excluded_by_polarity[side_polarity]:
                excluded_candidates.append(
                    {"polarity": side_polarity, **copy.deepcopy(record)}
                )
        spectral = auxiliary.get((positive["group_id"], negative["group_id"]))
        precursor_records = []
        for anchor in (positive, negative):
            precursor_records.append(
                {
                    "task_id": anchor["task_id"],
                    "polarity": anchor["polarity"],
                    "anchor_id": anchor["anchor_id"],
                    "candidate_independent_source": anchor["source"],
                    "source_group_id": anchor["group_id"],
                    "feature_ids": list(anchor["feature_ids"]),
                    "source_row_indices": list(anchor["source_row_indices"]),
                    "precursor_mz": anchor["precursor_mz"],
                    "neutral_mass": anchor["neutral_mass"],
                    "retention_time_min": anchor["retention_time_min"],
                    "aligned_retention_time_min": anchor[
                        "aligned_retention_time_min"
                    ],
                    "rt_basis": anchor["rt_basis"],
                    "adducts": list(anchor["adducts"]),
                    "candidate_identity_conflict_observed": (
                        _candidate_identity_conflict_observed(
                            anchor,
                            supplied[str(anchor["polarity"])][2],
                        )
                    ),
                    "candidate_identity_conflict_used_for_graph_membership": False,
                    "identity_claimed": False,
                }
            )
        edge_payload = {
            "edge_type": "candidate_independent_cross_polarity",
            "construction_inputs": [
                "explicit_task_id",
                "task_polarity",
                "reliable_neutral_mass",
                "retention_time",
                "hash_bound_rt_alignment_when_supplied",
            ],
            "candidate_identity_used": False,
            "spectrum_similarity_used": False,
            "neutral_mass_error_da": round(edge["neutral_mass_error_da"], 8),
            "neutral_mass_error_ppm": round(edge["neutral_mass_error_ppm"], 6),
            "neutral_mass_tolerance_da": round(edge["neutral_mass_tolerance_da"], 8),
            "rt_delta_min": round(edge["rt_delta_min"], 8),
            "rt_delta_seconds": round(edge["rt_delta_seconds"], 6),
            "rt_tolerance_min": round(edge["rt_tolerance_min"], 8),
            "rt_tolerance_seconds": round(edge["rt_tolerance_seconds"], 6),
            "pairing_strategy": "degree_one_unique_only",
            "pairing_status": "accepted_unique",
            "positive_anchor_degree": int(edge["positive_anchor_degree"]),
            "negative_anchor_degree": int(edge["negative_anchor_degree"]),
        }
        entities.append(
            {
                "neutral_entity_id": entity_id,
                "neutral_mass": round(edge["neutral_mass"], 8),
                "identity_claimed": False,
                "precursor_event_evidence_record_count": 2,
                "precursor_event_evidence": precursor_records,
                "candidate_independent_edge": edge_payload,
                "candidate_support": candidate_support,
                "quality_control_excluded_candidates": excluded_candidates,
                "candidate_support_attachment": {
                    "state": "attached_after_candidate_independent_graph_freeze",
                    "source": (
                        "v15_flagship_member_result_indices"
                        if candidate_support_by_task is not None
                        else "legacy_result_rows"
                    ),
                    "binding_field": (
                        "member_result_indices"
                        if candidate_support_by_task is not None
                        else "source_row_indices"
                    ),
                    "candidate_fields_used_for_graph": False,
                },
                "auxiliary_cross_polarity_ms2_similarity": spectral,
                "compound_count_increment": 0,
                "production_ranking_changed": False,
            }
        )
    entities.sort(key=lambda item: item["neutral_entity_id"])

    frozen_graph_hash_after = semantic_sha256(
        _frozen_candidate_independent_graph(
            aligned_anchors,
            compatible_edges,
            accepted_edges,
            ambiguous_edges,
            task_pair=task_pair,
            policy_sha256=policy_hash,
        )
    )
    if frozen_graph_hash_after != frozen_graph_hash_before:
        raise RuntimeError(
            "Candidate attachment changed the frozen candidate-independent graph"
        )

    source_hashes_after = {
        task_id: semantic_sha256(rows)
        for task_id, _, rows in supplied.values()
    }
    if source_hashes_after != source_hashes_before:
        raise RuntimeError("Dual-polarity graph construction mutated source results")
    candidate_support_hash_after = semantic_sha256(candidate_input)
    if candidate_support_hash_after != candidate_support_hash_before:
        raise RuntimeError("Dual-polarity graph construction mutated candidate support")

    dual_connectivities = {
        item["connectivity"]
        for entity in entities
        for item in entity["candidate_support"]
        if item["dual_polarity_support"]
    }
    ambiguous_anchor_ids = {
        str(edge["positive_anchor_id"]) for edge in ambiguous_edges
    } | {str(edge["negative_anchor_id"]) for edge in ambiguous_edges}
    quality_control_excluded_candidate_count = sum(
        len(entity["quality_control_excluded_candidates"]) for entity in entities
    )
    not_evaluable_anchor_count = (
        sum(
            int(extraction[polarity][field])
            for polarity in ("POS", "NEG")
            for field in (
                "rejected_conflict_groups",
                "rejected_inconsistent_neutral_mass_groups",
                "rejected_missing_rt_groups",
                "rejected_insufficient_member_groups",
            )
        )
        + int(alignment_summary["anchors_excluded_missing_mapping"])
    )
    if not aligned_anchors["POS"] or not aligned_anchors["NEG"]:
        evaluation_status = "not_evaluable_missing_reliable_polarity_anchor"
    elif entities:
        evaluation_status = "evaluated_with_edges"
    else:
        evaluation_status = "evaluated_no_edge"

    payload = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "metadata": {
            "module_id": MODULE_ID,
            "algorithm_id": ALGORITHM_ID,
            "module_state": "shadow",
            "task_pair": {
                "POS": supplied["POS"][0],
                "NEG": supplied["NEG"][0],
                "explicit_selection_required": True,
                "filename_used_for_pairing": False,
            },
            "source_results_semantic_sha256": dict(sorted(source_hashes_before.items())),
            "source_results_semantic_sha256_after": dict(sorted(source_hashes_after.items())),
            "source_task_metadata_sha256": dict(sorted(task_hashes.items())),
            "source_results_unchanged": True,
            "policy_sha256": policy_hash,
            "rt_alignment": alignment_summary,
            "candidate_support_input_semantic_sha256": candidate_support_hash_before,
            "candidate_support_input_semantic_sha256_after": candidate_support_hash_after,
            "candidate_support_input_unchanged": True,
            "candidate_support_source": (
                "v15_flagship_member_result_indices"
                if candidate_support_by_task is not None
                else "legacy_result_rows"
            ),
            "candidate_independent_graph_sha256_before_candidate_attachment": (
                frozen_graph_hash_before
            ),
            "candidate_independent_graph_sha256_after_candidate_attachment": (
                frozen_graph_hash_after
            ),
            "auxiliary_spectral_similarity_input_sha256": (
                semantic_sha256(list(auxiliary_spectral_similarity or []))
            ),
        },
        "summary": {
            "status": evaluation_status,
            "reliable_anchor_counts": {
                polarity: len(aligned_anchors[polarity])
                for polarity in ("POS", "NEG")
            },
            "anchor_extraction": extraction,
            "candidate_independent_compatible_edge_count": len(compatible_edges),
            "candidate_independent_accepted_edge_count": len(accepted_edges),
            "ambiguous_edge_count": len(ambiguous_edges),
            "neutral_entity_count": len(entities),
            "accepted_neutral_entity_count": len(entities),
            "ambiguous_anchor_count": len(ambiguous_anchor_ids),
            "not_evaluable_anchor_count": not_evaluable_anchor_count,
            "input_rows_without_reliable_anchor_count": sum(
                int(extraction[polarity]["rows_without_reliable_anchor"])
                for polarity in ("POS", "NEG")
            ),
            "precursor_event_evidence_record_count": 2 * len(entities),
            "dual_supported_connectivity_count": len(dual_connectivities),
            "dual_polarity_supported_connectivity_count": len(
                dual_connectivities
            ),
            "quality_control_excluded_candidate_count": (
                quality_control_excluded_candidate_count
            ),
            "compound_count_increment": 0,
        },
        "candidate_independent_graph": {
            "semantic_sha256": frozen_graph_hash_before,
            **frozen_graph,
        },
        "ambiguity": {
            "status": "present" if ambiguous_edges else "none",
            "edge_count": len(ambiguous_edges),
            "anchor_count": len(ambiguous_anchor_ids),
            "edges": copy.deepcopy(frozen_graph["ambiguous_no_merge_edges"]),
            "automatic_merge_allowed": False,
        },
        "neutral_entities": entities,
        "invariants": {
            "ground_truth_used_at_inference": False,
            "candidate_name_used_at_inference": False,
            "sample_filename_used_for_pairing": False,
            "candidate_identity_used_for_edge_construction": False,
            "spectrum_similarity_used_for_edge_construction": False,
            "cross_polarity_ms2_similarity_is_auxiliary_only": True,
            "top40_single_polarity_candidates_are_support_only": True,
            "candidate_independent_graph_unchanged_after_candidate_attachment": (
                frozen_graph_hash_before == frozen_graph_hash_after
            ),
            "only_degree_one_edges_create_neutral_entities": all(
                edge["positive_anchor_degree"] == 1
                and edge["negative_anchor_degree"] == 1
                and edge["pairing_status"] == "accepted_unique"
                for edge in accepted_edges
            ),
            "one_to_many_edges_are_ambiguous_no_merge": all(
                edge["pairing_status"] == "ambiguous_no_merge"
                for edge in ambiguous_edges
            ),
            "explicit_adduct_required_and_sign_checked": True,
            "candidate_identity_conflict_used_for_graph_membership": False,
            "candidate_display_name_used_for_structure_identity": False,
            "contaminant_or_exogenous_candidates_excluded_from_display_support": True,
            "excluded_candidates_retained_for_quality_control": True,
            "production_top1_changed": False,
            "production_top5_changed": False,
            "production_score_changed": False,
            "production_decision_changed": False,
            "production_confidence_changed": False,
            "production_compound_count_changed": False,
            "one_neutral_entity_has_two_precursor_event_evidence_records": all(
                entity["precursor_event_evidence_record_count"] == 2
                and len(entity["precursor_event_evidence"]) == 2
                for entity in entities
            ),
        },
    }
    payload["metadata"]["sidecar_semantic_sha256"] = sidecar_semantic_sha256(payload)
    return payload


__all__ = [
    "ALGORITHM_ID",
    "DEFAULT_POLICY_PATH",
    "MODULE_ID",
    "POLICY_SCHEMA_VERSION",
    "SIDECAR_SCHEMA_VERSION",
    "build_v16_dual_polarity_shadow",
    "load_v16_dual_polarity_policy",
    "rt_alignment_semantic_sha256",
    "semantic_sha256",
    "sidecar_semantic_sha256",
]
