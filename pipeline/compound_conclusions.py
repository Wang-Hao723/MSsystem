"""Compact, compound-level conclusions for ordinary MSsystem users.

The scan result and the V11 precursor-entity decision are audit artifacts.  A
sample owner needs a different unit: one row per reportable compound identity
or unresolved molecular-formula isomer group.  This module builds that view
without changing candidate ranking, legacy confidence, or the final-claim
contract.

The view deliberately separates three questions that must not be collapsed:

* did a library spectrum support this identity;
* is the signal enriched in the sample rather than the paired blank; and
* is there enough orthogonal evidence to call the compound confirmed.

Classification, expert reasoning, biosynthetic reasoning and fragmentation
paths are retained as compact explanations.  They never enter the identity
ranking or the conclusion tier.
"""

from __future__ import annotations

import hashlib
import json
from runtime.json_encoding import write_json_encoded
import math
import os
import re
import uuid
from copy import deepcopy
from collections import Counter
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from matching.flagship_ranker import connectivity_identity
from pipeline.entity_confidence import ALGORITHM_ID as ENTITY_HIGH_ALGORITHM_ID, single_energy_high_evidence
from matching.biosource_context import biosource_context
from pipeline.compound_role_projection import (
    annotate_match_compound_roles, apply_report_compound_roles,
    plant_report_selection, report_role_catalog_summary,
)
from matching.compound_names import normalize_match_compound_names
from matching.fragmentation_rules import get_neutral_loss_by_mass
from matching.library_quality import formula_exact_mass, normalize_formula_text


from matching.v15_flagship_ranker import load_v15_flagship_policy
from pipeline.active_artifact_binding import active_implementation_sha256
from pipeline.derived_artifact_paths import (
    derived_artifact_path, binding_metadata, assert_derived_destination,
    metadata_binding_is_current, verified_primary_binding,
)


COMPOUND_CONCLUSION_SCHEMA_VERSION = "2.1.1"
COMPOUND_CONCLUSION_ALGORITHM_ID = "mssystem-v16-feature-result-3.2.0-bound-spectral-evidence"
COMPOUND_CONCLUSION_ARTIFACT_SUFFIX = "_v11_compounds.json"

_CONFIDENCE_RANK = {"High": 4, "Medium": 3, "Low": 2, "Unknown": 1}
_SAMPLE_SUPPORTED_CATEGORIES = {"sample_enriched", "not_detected"}
_BACKGROUND_CATEGORIES = {"blank_associated", "blank_dominant"}

_MATCH_DISPLAY_FIELDS = (
    "compound_name",
    "compound_name_cn",
    "source_compound_name",
    "formula",
    "smiles",
    "has_structure",
    "structure_quality",
    "inchikey",
    "reported_inchikey",
    "adduct",
    "candidate_adduct_hypothesis",
    "candidate_adduct_hypothesis_source",
    "mass_metadata_consistent",
    "precursor_metadata_audit",
    "dataset_source",
    "source_tier",
    "lib_spectrum_id",
    "lib_source_accession",
    "lib_collision_energy",
    "lib_collision_energy_raw",
    "lib_collision_energy_unit",
    "lib_collision_energy_mode",
    "ce_comparability",
    "cosine",
    "entropy",
    "reverse_match",
    "mass_ppm",
    "mass_score",
    "fragment_score",
    "intensity_score",
    "final_score",
    "rank_score",
    "confidence",
    "confidence_class",
    "exact_structure_confidence",
    "two_dimensional_connectivity_confidence",
    "isomer_group_confidence",
    "identification_scope",
    "honesty_tier",
    "evidence_tier",
    "annotation_level",
    "isomer_ambiguous",
    "isomer_n",
    "isomer_group_formula",
    "isomer_group_size",
    "isomer_ambiguity_type",
    "isomer_resolution_status",
    "candidate_margin",
    "precursor_high_evidence_gate_required",
    "precursor_high_evidence_gate_result",
    "precursor_high_evidence_gate_reason",
    "precursor_high_evidence_gate_reason_zh",
    "evidence_governance",
    "direct_matched_product_ions",
    "name_only_lead",
    "classification",
    "display_class",
    "classification_source",
    "display_classification_source",
    "classification_reliable",
    "classification_status",
    "classification_reliability_reason",
    "classification_confidence",
    "display_classification_confidence",
    "expert_score",
    "expert_score_source",
    "potential_contaminant",
    "report_role_projection",
    "role_annotation_scope",
    "role_annotation_source",
    "role_catalog_id",
    "role_catalog_version",
    "contaminant_note",
    "contaminant_category",
    "compound_role",
    "compound_role_note",
    "exclude_from_natural_product_count",
    "natural_product_candidate",
    "screening_domain",
    "pesticide_use",
    "chemical_family",
    "regulatory_status",
    "pathway_polarity_status",
    "pathway_adduct_status",
    "flagship_entity_score",
    "flagship_entity_rank",
    "primary_identity_source",
    "primary_entity_id",
    "primary_ranking_identity_key",
    "primary_decision_status",
    "primary_decision_scope",
    "primary_identity_resolution",
    "entity_high_evidence",
    "spectral_evidence_binding",
    "native_spectrum_confidence",
    "entity_candidate_confidence_source",
)

_ROW_DISPLAY_FIELDS = (
    "spectrum_id",
    "precursor_mz",
    "search_precursor_mz",
    "search_precursor_mz_source",
    "lcms_feature_used_for_search",
    "isotope_evidence_status",
    "retention_time",
    "retention_time_unit",
    "tic",
    "ion_mode",
    "precursor_charge",
    "collision_energy",
    "collision_energy_unit",
    "collision_energy_mode",
    "collision_energy_quantities",
    "activation_method",
    "acquisition_metadata_source",
    "effective_fragment_count",
    "fragmentation_status",
    "background_evidence",
    "precursor_feature_id",
    "lcms_feature_id",
    "ion_identity_network_id",
    "query_adduct",
    "query_adduct_reliable",
    "lcms_adduct",
    "lcms_neutral_mass",
    "lcms_adduct_group_id",
    "ion_identity_status",
    "ion_identity_neutral_mass",
    "blank_used",
    "blank_status",
    "blank_flagged",
    "blank_feature_matched",
    "blank_sample_ratio",
    "blank_evidence_category",
    "blank_evidence_role",
    "lcms_feature_inventory_only",
    "lcms_feature_quality",
    "lcms_feature_rt_start_minutes",
    "lcms_feature_rt_end_minutes",
)


def compound_conclusion_artifact_path(task_dir: Path, task_id: str) -> Path:
    policy = load_v15_flagship_policy()
    return derived_artifact_path(task_dir, task_id, COMPOUND_CONCLUSION_ARTIFACT_SUFFIX,
        policy=policy, actual_implementation_sha256=active_implementation_sha256(policy))


def _text(value: Any) -> str:
    return " ".join(str(value or "").split())


def _finite(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


_PRIMARY_CANDIDATE_FIELDS = (
    "compound_name",
    "compound_name_cn",
    "source_compound_name",
    "formula",
    "smiles",
    "inchikey",
    "reported_inchikey",
    "adduct",
    "has_structure",
    "structure_quality",
    "name_only_lead",
    "dataset_source",
    "source_tier",
    "lib_spectrum_id",
    "lib_source_accession",
    "classification",
    "display_class",
    "classification_source",
    "display_classification_source",
    "classification_reliable",
    "classification_status",
    "classification_reliability_reason",
    "classification_confidence",
    "display_classification_confidence",
    "feature_fragment_evidence",
    "potential_contaminant",
    "contaminant_category",
    "compound_role",
    "compound_role_note",
)


def _primary_payload(
    artifact: Mapping[str, Any] | Path | None,
) -> tuple[dict[str, Any] | None, Path | None]:
    if artifact is None:
        return None, None
    if isinstance(artifact, Mapping):
        return deepcopy(dict(artifact)), None
    path = Path(artifact)
    return json.loads(path.read_text(encoding="utf-8")), path


def _primary_entity_index(
    task_id: str,
    artifact: Mapping[str, Any] | Path | None,
) -> tuple[dict[int, dict[str, Any]], dict[str, Any], Path | None]:
    payload, path = _primary_payload(artifact)
    if payload is None:
        return {}, {}, path
    metadata = payload.get("metadata") or {}
    if str(metadata.get("task_id") or "") != str(task_id):
        raise ValueError("Primary V15 artifact belongs to another task")
    if metadata.get("production_identity_selector") is not True:
        raise ValueError("Primary V15 artifact is not a production identity selector")
    index: dict[int, dict[str, Any]] = {}
    for entity in payload.get("entities") or []:
        if not isinstance(entity, Mapping):
            continue
        compact = deepcopy(dict(entity))
        for result_index in entity.get("member_result_indices") or []:
            result_index = int(result_index)
            if result_index in index:
                raise ValueError(
                    "One result row is assigned to multiple V15 precursor entities"
                )
            index[result_index] = compact
    return index, dict(metadata), path


def _primary_candidate_rows(entity: Mapping[str, Any]) -> list[dict[str, Any]]:
    resolution = entity.get("identity_resolution") or {}
    competitive = set(
        str(value or "")
        for value in resolution.get("competitive_candidate_identity_keys") or []
    )
    rows = []
    for candidate in entity.get("candidates") or []:
        if not isinstance(candidate, Mapping):
            continue
        display = deepcopy(dict(candidate.get("candidate") or {}))
        display.update(
            {
                "flagship_identity_key": _text(candidate.get("identity_key")),
                "flagship_entity_rank": int(candidate.get("entity_rank") or 0),
                "flagship_entity_score": _finite(candidate.get("entity_score")),
                "entity_score": _finite(candidate.get("entity_score")),
                "rank_score": _finite(candidate.get("entity_score")),
                "final_score": _finite(candidate.get("entity_score")),
                "score_is_probability": False,
                "spectrally_competitive": _text(candidate.get("identity_key"))
                in competitive,
            }
        )
        rows.append(display)
    return rows


_SAVED_SPECTRAL_METRICS = (
    "cosine", "entropy", "reverse_match", "mass_ppm", "mass_score",
    "fragment_score", "intensity_score", "final_score", "rank_score",
)


def _selected_spectrum_evidence(row, selected, result_index=None):
    """Bind measured scores to one saved scan, full identity and reference ID.

    Entity scores cannot supply a spectral cosine or a precursor error. Names
    and 2-D connectivity alone also cannot transfer another reference's scores.
    """
    key = _text(selected.get("inchikey")).upper()
    reference_id = _text(selected.get("lib_spectrum_id"))
    scan = _text(row.get("spectrum_id"))
    binding = {
        "schema": "same-scan-selected-reference/1", "status": "unavailable",
        "spectrum_id": scan, "result_index": result_index,
        "candidate_inchikey": key, "lib_spectrum_id": reference_id,
        "metrics": {}, "saved_confidence": "Unknown",
    }
    if not scan or not key or not reference_id:
        binding["reason"] = "saved_scan_full_identity_or_reference_id_missing"
        return None, binding
    original = row.get("match") or {}
    candidates = [("match", original)] + [
        (f"match.top_candidates[{index}]", candidate)
        for index, candidate in enumerate(original.get("top_candidates") or [])
        if isinstance(candidate, Mapping)
    ]
    matched = []
    for source, candidate in candidates:
        if (_text(candidate.get("inchikey")).upper() != key
                or _text(candidate.get("lib_spectrum_id")) != reference_id):
            continue
        if candidate.get("spectrum_id") is not None and _text(candidate["spectrum_id"]) != scan:
            continue
        if (candidate.get("result_index") is not None and result_index is not None
                and candidate["result_index"] != result_index):
            continue
        metrics = {name: candidate[name] for name in _SAVED_SPECTRAL_METRICS
                   if _finite(candidate.get(name)) is not None}
        # A current scan's original top match is the authoritative saved
        # presentation. Its duplicate in top_candidates can have more digits.
        matched.append((source, candidate, metrics))
        if source == "match":
            break
    if not matched:
        binding["reason"] = "selected_full_identity_and_reference_not_saved_on_this_scan"
        return None, binding
    if len(matched) > 1:
        projections = {_canonical((metrics, candidate.get("confidence")))
                       for _, candidate, metrics in matched}
        if len(projections) != 1:
            binding["reason"] = "conflicting_saved_reference_metrics"
            return None, binding
    source, candidate, metrics = matched[0]
    binding.update(status="bound_saved_reference", source_field=source,
                   metrics=deepcopy(metrics),
                   saved_confidence=_text(candidate.get("confidence")) or "Unknown")
    return candidate, binding


def _candidate_only_saved_grade(binding):
    if binding.get("status") != "bound_saved_reference":
        return "Unknown"
    grade = binding.get("saved_confidence")
    return "Medium" if grade in {"High", "Medium"} else "Low" if grade == "Low" else "Unknown"


def _project_primary_identity(
    row: Mapping[str, Any], entity: Mapping[str, Any] | None,
    high_evidence: Mapping[str, Any] | None = None,
    *, result_index: int | None = None,
) -> dict[str, Any]:
    projected = deepcopy(dict(row))
    source_match = row.get("match")
    if not isinstance(source_match, Mapping) or not isinstance(entity, Mapping):
        return projected
    match = deepcopy(dict(source_match))
    selected = entity.get("winner") or entity.get("leading_candidate_when_abstained")
    if not isinstance(selected, Mapping):
        match.update(
            {
                "primary_identity_source": "v15_flagship_no_candidate",
                "primary_entity_id": _text(entity.get("entity_id")),
                "primary_ranking_identity_key": "",
                "primary_decision_status": _text(
                    (entity.get("decision") or {}).get("status")
                ),
                "primary_decision_scope": _text(
                    (entity.get("decision") or {}).get("scope")
                ),
                "primary_identity_resolution": deepcopy(
                    dict(entity.get("identity_resolution") or {})
                ),
                "confidence": "Unknown",
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": "Unknown",
                "identification_scope": "none",
                "isomer_ambiguous": False,
            }
        )
        projected["match"] = match
        return projected

    selected_display = deepcopy(dict(selected.get("candidate") or {}))
    saved_reference, spectral_binding = _selected_spectrum_evidence(
        row, selected_display, result_index,
    )
    # Never carry the original Top1 candidate's evidence into a different
    # entity-selected identity. The source row itself remains unchanged.
    match = deepcopy(dict(saved_reference)) if saved_reference is not None else {}
    for field in _PRIMARY_CANDIDATE_FIELDS:
        match.pop(field, None)
        if field in selected_display:
            match[field] = selected_display[field]
    for field in (
        "cosine",
        "entropy",
        "reverse_match",
        "mass_ppm",
        "mass_score",
        "fragment_score",
        "intensity_score",
    ):
        match.pop(field, None)
    for field in (
        "isomer_group_formula",
        "isomer_group_size",
        "isomer_group_evidence",
        "isomer_n",
        "isomer_ambiguity_type",
        "isomer_resolution_status",
    ):
        match.pop(field, None)

    decision = entity.get("decision") or {}
    resolution = deepcopy(dict(entity.get("identity_resolution") or {}))
    ranked = _primary_candidate_rows(entity)
    selected_identity = _text(selected.get("identity_key"))
    selected_row = next(
        (
            candidate
            for candidate in ranked
            if _text(candidate.get("flagship_identity_key")) == selected_identity
        ),
        deepcopy(selected_display),
    )
    match.update(selected_row)
    for field in _SAVED_SPECTRAL_METRICS:
        match.pop(field, None)
        if (field not in {"final_score", "rank_score"}
                and field in spectral_binding["metrics"]):
            match[field] = spectral_binding["metrics"][field]
    match.update(
        {
            "top_candidates": ranked,
            "flagship_entity_score": _finite(selected.get("entity_score")),
            "flagship_entity_rank": int(selected.get("entity_rank") or 1),
            "rank_score": _finite(selected.get("entity_score")),
            "final_score": _finite(selected.get("entity_score")),
            "primary_identity_source": "v15_flagship",
            "primary_entity_id": _text(entity.get("entity_id")),
            "primary_ranking_identity_key": selected_identity,
            "primary_decision_status": _text(decision.get("status")),
            "primary_decision_scope": _text(decision.get("scope")),
            "primary_identity_resolution": resolution,
            "score_is_probability": False,
            "spectral_evidence_binding": spectral_binding,
            "native_spectrum_confidence": spectral_binding["saved_confidence"],
        }
    )
    competitive_rows = [
        candidate
        for candidate in ranked
        if candidate.get("spectrally_competitive")
        and _text(candidate.get("flagship_identity_key")) != selected_identity
    ]
    match["isomer_set"] = competitive_rows
    match["same_formula_alternatives"] = [
        candidate
        for candidate in competitive_rows
        if normalize_formula_text(candidate.get("formula"))
        == normalize_formula_text(selected_row.get("formula"))
    ]
    match["same_formula_alternatives_n"] = len(
        match["same_formula_alternatives"]
    )

    upstream_status = _text(resolution.get("upstream_decision_status"))
    scope = _text(decision.get("scope"))
    if scope == "unresolved_precursor_event":
        match.update(
            {
                "confidence": "Unknown",
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": "Unknown",
                "identification_scope": "unresolved_precursor_event",
                "isomer_ambiguous": False,
            }
        )
    elif scope == "molecular_formula_isomer_group":
        group_confidence = ("High" if upstream_status == "selected"
                            else _candidate_only_saved_grade(spectral_binding))
        match.update(
            {
                "confidence": group_confidence,
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": group_confidence,
                "identification_scope": "isomer_group",
                "isomer_group_formula": resolution.get("reportable_formula")
                or selected_row.get("formula"),
                "isomer_group_size": max(2, len(competitive_rows) + 1),
                "isomer_ambiguous": True,
                "isomer_ambiguity_type": "v15_competitive_same_formula_connectivities",
            }
        )
    elif decision.get("status") == "selected" and scope == "2d_connectivity_candidate":
        match.update(
            {
                "confidence": "High",
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "High",
                "isomer_group_confidence": "Unknown",
                "identification_scope": "exact_2d_connectivity",
                "isomer_ambiguous": False,
            }
        )
    elif decision.get("status") == "selected_candidate_only":
        saved_grade = _candidate_only_saved_grade(spectral_binding)
        match.update(
            {
                "confidence": saved_grade,
                "confidence_class": saved_grade,
                "evidence_tier": "Strong Candidate" if saved_grade == "Medium" else "",
                "entity_candidate_confidence_source": (
                    "same_scan_saved_reference_capped_at_medium"
                    if spectral_binding["status"] == "bound_saved_reference"
                    else "unassessed_missing_saved_reference"
                ),
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": "Unknown",
                "identification_scope": scope or "candidate_only",
                "isomer_ambiguous": False,
            }
        )
    else:
        match.update(
            {
                "confidence": "Unknown",
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": "Unknown",
                "identification_scope": scope or "none",
                "isomer_ambiguous": False,
            }
        )
    if high_evidence:
        high_scope = high_evidence["scope"]
        match.update({
            "confidence": "High",
            "exact_structure_confidence": "Unknown",
            "two_dimensional_connectivity_confidence": "High" if high_scope == "2d_connectivity" else "Unknown",
            "isomer_group_confidence": "High" if high_scope == "isomer_group" else "Unknown",
            "identification_scope": "exact_2d_connectivity" if high_scope == "2d_connectivity" else "isomer_group",
            "entity_high_evidence": deepcopy(dict(high_evidence)),
        })
    projected["match"] = match
    return projected


def _blank_category(row: Mapping[str, Any]) -> str:
    category = _text(row.get("blank_evidence_category")).lower()
    if category:
        return category
    status = _text(row.get("blank_status")).lower()
    if "sample_enriched" in status:
        return "sample_enriched"
    if "blank_dominant" in status:
        return "blank_dominant"
    if "blank_associated" in status:
        return "blank_associated"
    if not bool(row.get("blank_used")):
        return "not_assessed"
    if not bool(row.get("blank_feature_matched")):
        return "not_detected"
    return "not_assessed"


def _sample_supported(row: Mapping[str, Any]) -> bool:
    category = _blank_category(row)
    if category in _SAMPLE_SUPPORTED_CATEGORIES:
        return True
    if category in _BACKGROUND_CATEGORIES:
        # An untyped blank selected only by ion polarity is not a verified
        # same-batch/same-method paired blank.  Its signal is valuable as a
        # background warning, but it is not strong enough to remove a sample
        # candidate from the primary result set.
        if _text(row.get("blank_evidence_role")).lower() == "unspecified":
            return True
        return False
    # A task without a paired blank cannot be treated as background-free, but
    # it must remain visible as a sample candidate.  The UI reports that the
    # background check was unavailable.
    return True


def _isomer_group_row(match: Mapping[str, Any]) -> bool:
    scope = _text(match.get("identification_scope")).lower()
    if scope == "exact_2d_connectivity":
        # Several full-InChIKey variants can share one 2-D connectivity.  They
        # remain unresolved above that layer, but they are not a formula-wide
        # constitutional/positional isomer group.
        return False
    return bool(match.get("isomer_ambiguous")) or scope in {
        "isomer_group",
        "molecular_formula_isomer_group",
        "stereoisomer_group",
    }


def _group_key(match: Mapping[str, Any]) -> tuple[str, str, str]:
    formula = normalize_formula_text(
        match.get("isomer_group_formula") or match.get("formula")
    )
    if _isomer_group_row(match) and formula:
        return f"isomer-formula:{formula}", "isomer_group", formula
    identity, source = connectivity_identity(match)
    # Library metadata occasionally reuses a structure string across entries
    # carrying incompatible neutral formulas.  A 2-D identity key must never
    # merge those records merely because their stored SMILES are identical.
    # Formula-scoping preserves legitimate multi-adduct deduplication while
    # making cross-formula metadata conflicts visible instead of silent.
    if formula:
        identity = f"{identity}|formula:{formula}"
        source = f"{source}+neutral_formula"
    return identity, "2d_connectivity", source


def _event_identity(match: Mapping[str, Any]) -> str:
    """Return the selected 2-D identity before any formula-group collapse."""

    identity, _ = connectivity_identity(match)
    return identity


def _analysis_feature_id(
    row: Mapping[str, Any], result_index: int | None = None
) -> str:
    """Return the identity-free per-run result unit for one observation.

    A measured continuous-MS1 feature is authoritative. The answer-neutral
    precursor-event cluster is an explicit fallback for formats or scans that
    could not be linked to a chromatographic feature. An unresolved row is
    kept separate rather than being merged by a candidate name.
    """

    lcms_feature_id = _text(row.get("lcms_feature_id"))
    if lcms_feature_id:
        return f"lcms:{lcms_feature_id}"
    precursor_feature_id = _text(row.get("precursor_feature_id"))
    if precursor_feature_id:
        return f"fallback:{precursor_feature_id}"
    spectrum_id = _text(row.get("spectrum_id"))
    if spectrum_id:
        return f"unresolved:{spectrum_id}"
    return f"unresolved-row:{int(result_index or 0)}"


def _materialize_feature_inventory(
    feature_inventory: Iterable[Mapping[str, Any]] | None,
) -> tuple[list[dict[str, Any]], set[str]]:
    """Translate the independent MS1 feature table into display-only rows.

    These rows make continuous-MS1 features without an associated DDA scan
    visible. They never enter library matching or identity selection. The
    upstream feature identifier must be unique because silently collapsing a
    malformed inventory would make the completeness check circular again.
    """

    if feature_inventory is None:
        return [], set()
    rows: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    for index, feature in enumerate(feature_inventory):
        if not isinstance(feature, Mapping):
            raise ValueError(f"Feature inventory row {index} is not an object")
        raw_id = _text(feature.get("feature_id"))
        if not raw_id:
            raise ValueError(f"Feature inventory row {index} has no feature_id")
        feature_id = f"lcms:{raw_id}"
        if feature_id in identifiers:
            raise ValueError(f"Duplicate feature inventory identifier: {raw_id}")
        identifiers.add(feature_id)
        rt_minutes = _finite(feature.get("rt_minutes"))
        rows.append({
            "lcms_feature_id": raw_id,
            "precursor_mz": _finite(feature.get("mz")),
            "retention_time": rt_minutes,
            "retention_time_unit": "minute" if rt_minutes is not None else None,
            "tic": _finite(feature.get("intensity")),
            "ion_mode": _text(feature.get("ion_mode")).upper() or None,
            "precursor_charge": feature.get("charge"),
            "lcms_adduct": feature.get("adduct"),
            "lcms_neutral_mass": _finite(feature.get("neutral_mass")),
            "lcms_adduct_group_id": feature.get("adduct_group_id"),
            "lcms_feature_quality": _finite(feature.get("quality")),
            "lcms_feature_rt_start_minutes": _finite(
                feature.get("rt_start_minutes")
            ),
            "lcms_feature_rt_end_minutes": _finite(
                feature.get("rt_end_minutes")
            ),
            "lcms_feature_inventory_only": True,
            "feature_inventory_index": index,
        })
    return rows, identifiers


def _best_event_scope_confidence(matches: Iterable[Mapping[str, Any]]) -> str:
    """Keep the strongest *single-scan* support at the event's honest scope.

    Stepped collision energies are repeated measurements of one precursor
    event, not independent votes.  Taking the best existing tier preserves the
    evidence from one useful spectrum without summing support across energies.
    """

    best_label = "Unknown"
    best_rank = _CONFIDENCE_RANK[best_label]
    for match in matches:
        for value in (
            match.get("isomer_group_confidence"),
            match.get("two_dimensional_connectivity_confidence"),
            match.get("confidence"),
        ):
            label = _text(value)
            rank = _CONFIDENCE_RANK.get(label, 0)
            if rank > best_rank:
                best_label = label
                best_rank = rank
    return best_label


def _precursor_event_plans(
    rows: Iterable[Mapping[str, Any]],
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Resolve scan-level candidates once per answer-independent precursor event.

    Each measured LC-MS feature, or explicit precursor-event fallback,
    receives exactly one grouping plan:

    * consistent identities retain their normal compound key;
    * conflicting structures with one known formula become one unresolved
      formula-isomer group; and
    * incompatible formula assignments become one unresolved formula-conflict
      event; and
    * structural disagreement with an incomplete formula assignment becomes
      one separately labelled unresolved event.

    Both unresolved-event variants are excluded from the compound count.  The
    distinction prevents a missing database formula from being reported as if
    two known molecular formulae had been shown to conflict.
    """

    event_rows: dict[str, list[Mapping[str, Any]]] = {}
    for result_index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            continue
        match = row.get("match")
        if not isinstance(match, Mapping):
            continue
        if _text(match.get("dataset_source")).lower() == "iceberg_fallback":
            continue
        event_id = _analysis_feature_id(row, result_index)
        if event_id:
            event_rows.setdefault(event_id, []).append(match)

    plans: dict[str, dict[str, Any]] = {}
    counts = Counter()
    for event_id, matches in event_rows.items():
        formulas = [
            normalize_formula_text(
                match.get("isomer_group_formula") or match.get("formula")
            )
            for match in matches
        ]
        known_formulas = sorted({formula for formula in formulas if formula})
        formula_complete = bool(formulas) and all(formulas)
        selected_identities = {_event_identity(match) for match in matches}
        selected_identities.discard("")
        group_keys = {_group_key(match)[0] for match in matches}
        disagrees = len(selected_identities) > 1 or len(group_keys) > 1
        formula_conflict = len(known_formulas) > 1
        primary_resolutions = [
            match.get("primary_identity_resolution") or {}
            for match in matches
            if isinstance(match.get("primary_identity_resolution"), Mapping)
        ]
        primary_statuses = {
            _text(resolution.get("status")) for resolution in primary_resolutions
        }
        primary_formulae = sorted(
            {
                normalize_formula_text(formula)
                for resolution in primary_resolutions
                for formula in resolution.get("candidate_formulae") or []
                if normalize_formula_text(formula)
            }
        )
        primary_incomplete = (
            "unresolved_incomplete_formula_event" in primary_statuses
        )
        primary_cross_formula = (
            "unresolved_cross_formula_event" in primary_statuses
        )
        primary_isomer_group = (
            "unresolved_same_formula_isomer_group" in primary_statuses
        )
        if primary_formulae:
            known_formulas = primary_formulae
        formula_complete = formula_complete and not primary_incomplete

        counts["matched_event_count"] += 1
        if len(matches) > 1:
            counts["multi_spectrum_event_count"] += 1

        if primary_cross_formula or formula_conflict:
            status = "cross_formula_conflict"
            plan = {
                "status": status,
                "group_key": f"unresolved-event:{event_id}",
                "identity_scope": "unresolved_precursor_event",
                "identity_source": "precursor_event_formula_conflict",
                "formula": "",
                "candidate_formulas": known_formulas,
                "scope_confidence": "Unknown",
            }
            counts["cross_formula_unresolved_event_count"] += 1
        elif (
            primary_isomer_group
            or (disagrees and formula_complete and len(known_formulas) == 1)
        ):
            formula = known_formulas[0]
            status = "same_formula_isomer_group"
            plan = {
                "status": status,
                "group_key": f"isomer-event:{event_id}|formula:{formula}",
                "identity_scope": "isomer_group",
                "identity_source": (
                    "v15_primary_competitive_same_formula_connectivities"
                    if primary_isomer_group
                    else "precursor_event_same_formula_conflict"
                ),
                "formula": formula,
                "candidate_formulas": known_formulas,
                "scope_confidence": _best_event_scope_confidence(matches),
            }
            counts["same_formula_unresolved_event_count"] += 1
        elif primary_incomplete or disagrees:
            status = "formula_assignment_incomplete"
            plan = {
                "status": status,
                "group_key": f"unresolved-event:{event_id}",
                "identity_scope": "unresolved_precursor_event",
                "identity_source": "precursor_event_formula_assignment_incomplete",
                "formula": "",
                "candidate_formulas": known_formulas,
                "scope_confidence": "Unknown",
            }
            counts["incomplete_formula_unresolved_event_count"] += 1
        else:
            status = "consistent_identity"
            plan = {
                "status": status,
                "group_key": None,
                "identity_scope": None,
                "identity_source": None,
                "formula": known_formulas[0] if len(known_formulas) == 1 else "",
                "candidate_formulas": known_formulas,
                "scope_confidence": _best_event_scope_confidence(matches),
            }
            counts["consistent_identity_event_count"] += 1
        plans[event_id] = plan

    summary = {
        "policy": "one_report_row_per_ion_feature",
        "matched_event_count": int(counts["matched_event_count"]),
        "multi_spectrum_event_count": int(counts["multi_spectrum_event_count"]),
        "consistent_identity_event_count": int(
            counts["consistent_identity_event_count"]
        ),
        "same_formula_unresolved_event_count": int(
            counts["same_formula_unresolved_event_count"]
        ),
        "cross_formula_unresolved_event_count": int(
            counts["cross_formula_unresolved_event_count"]
        ),
        "incomplete_formula_unresolved_event_count": int(
            counts["incomplete_formula_unresolved_event_count"]
        ),
    }
    return plans, summary


def _result_priority(row: Mapping[str, Any]) -> tuple[Any, ...]:
    match = row.get("match") if isinstance(row.get("match"), Mapping) else {}
    confidence = _CONFIDENCE_RANK.get(_text(match.get("confidence")), 0)
    isomer_confidence = _CONFIDENCE_RANK.get(
        _text(match.get("isomer_group_confidence")), 0
    )
    rank_score = _finite(match.get("rank_score") or match.get("final_score")) or 0.0
    cosine = _finite(match.get("cosine")) or 0.0
    ppm = _finite(match.get("mass_ppm"))
    return (
        -int(_sample_supported(row)),
        -max(confidence, isomer_confidence),
        -rank_score,
        -cosine,
        abs(ppm) if ppm is not None else 1e12,
        _text(row.get("spectrum_id")),
    )


def _compact_step(step: Mapping[str, Any]) -> dict[str, Any]:
    keys = (
        "type",
        "from_mz",
        "to_mz",
        "mass_difference",
        "neutral_loss",
        "neutral_loss_en",
        "formula",
        "mechanism",
        "assignment_level",
        "assignment_status",
        "mass_delta",
        "priority",
    )
    return {key: step.get(key) for key in keys if step.get(key) is not None}


def _specific_adduct(value: Any) -> str:
    text = _text(value)
    if text.lower() in {
        "", "nan", "none", "null", "n/a", "unknown", "pos", "positive",
        "+", "neg", "negative", "-",
    }:
        return ""
    return text


def _adduct_interpretation(
    match: Mapping[str, Any], row: Mapping[str, Any]
) -> dict[str, Any]:
    """Separate observed sample ion form from candidate-library metadata."""

    library = _specific_adduct(match.get("adduct"))
    candidate_hypothesis = _specific_adduct(
        match.get("candidate_adduct_hypothesis")
    )
    candidate_hypothesis_source = _text(
        match.get("candidate_adduct_hypothesis_source")
    )
    lcms = _specific_adduct(row.get("lcms_adduct"))
    query = _specific_adduct(row.get("query_adduct") or match.get("query_adduct"))
    if lcms:
        sample, source, independent, reliable = (
            lcms, "lcms_ion_identity", True, True
        )
    elif query and bool(
        row.get("query_adduct_reliable") or match.get("query_adduct_reliable")
    ):
        sample, source, independent, reliable = (
            query, "explicit_query_metadata", True, True
        )
    elif library:
        sample, source, independent, reliable = (
            library, "library_candidate_hypothesis", False, False
        )
    elif candidate_hypothesis:
        sample, source, independent, reliable = (
            candidate_hypothesis,
            "candidate_structure_formula_precursor_hypothesis",
            False,
            False,
        )
    else:
        sample, source, independent, reliable = (
            None, "unresolved", False, False
        )
    consistency = (
        "not_tested" if not sample or not library
        else "consistent" if re.sub(r"\s+", "", sample).upper()
        == re.sub(r"\s+", "", library).upper()
        else "inconsistent"
    )
    neutral_mass = _finite(
        row.get("lcms_neutral_mass") or row.get("ion_identity_neutral_mass")
    )
    return {
        "sample_ion_form": sample,
        "sample_ion_form_source": source,
        "sample_ion_form_reliable": reliable,
        "candidate_independent": independent,
        "neutral_mass": neutral_mass,
        "library_adduct": library or None,
        "candidate_ion_form_hypothesis": candidate_hypothesis or None,
        "candidate_ion_form_hypothesis_source": (
            candidate_hypothesis_source or None
        ),
        "sample_library_consistency": consistency,
        "status": (
            "assigned_from_sample_evidence" if independent
            else "candidate_hypothesis_only" if sample
            else "unresolved"
        ),
        "polarity_alone_used": False,
    }


def _validated_step(
    step: Mapping[str, Any], match: Mapping[str, Any], row: Mapping[str, Any]
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    """Revalidate legacy named losses under the current accurate-mass policy."""

    compact = _compact_step(step)
    formula = normalize_formula_text(step.get("formula"))
    observed = _finite(step.get("mass_difference"))
    if not formula or observed is None:
        return compact, None
    exact = formula_exact_mass(formula)
    error = abs(float(exact) - observed) if exact is not None else None
    adduct = _specific_adduct(row.get("lcms_adduct") or row.get("query_adduct") or match.get("adduct"))
    applicable = get_neutral_loss_by_mass(
        observed,
        tolerance=0.01,
        classification=_text(match.get("display_class") or match.get("classification")),
        strict_fg=True,
        polarity=_text(row.get("ion_mode")),
        candidate_smiles=_text(match.get("smiles")),
        query_adduct=adduct,
        query_collision_energy=row.get("collision_energy"),
    )
    supported = any(normalize_formula_text(item.formula) == formula for item in applicable)
    if exact is not None and error is not None and error <= 0.01 and supported:
        compact["assignment_status"] = "current_policy_supported"
        compact["mass_delta"] = round(error, 6)
        return compact, None
    rejected = {
        **compact,
        "assignment_status": "legacy_assignment_rejected_by_current_policy",
        "mass_delta": round(error, 6) if error is not None else None,
        "reason": (
            "accurate_mass_mismatch" if error is not None and error > 0.01
            else "candidate_structure_or_class_not_applicable"
        ),
    }
    return None, rejected


def _compact_product_structures(details: Mapping[str, Any]) -> list[dict[str, Any]]:
    analysis = details.get("structure_fragment_analysis")
    analysis = analysis if isinstance(analysis, Mapping) else {}
    values = analysis.get("plausible_fragments")
    values = values if isinstance(values, list) else []
    keys = (
        "observed_mz", "fragment_smiles", "fragment_mass", "expected_mz",
        "delta", "adduct", "assignment_type",
    )
    output = []
    for value in values:
        if not isinstance(value, Mapping) or not _text(value.get("fragment_smiles")):
            continue
        output.append({key: value.get(key) for key in keys if value.get(key) is not None})
        if len(output) >= 12:
            break
    return output


def _expert_dimensions(details: Mapping[str, Any]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    raw = details.get("dimension_assessments")
    raw = raw if isinstance(raw, Mapping) else {}
    dimensions: list[dict[str, Any]] = []
    counts = {
        "scored": 0,
        "no_evidence": 0,
        "not_applicable": 0,
        "unavailable": 0,
    }
    for key, value in raw.items():
        if not isinstance(value, Mapping):
            continue
        status = _text(value.get("status")) or "unavailable"
        normalized_status = status if status in counts else "unavailable"
        counts[normalized_status] += 1
        dimensions.append({
            "key": _text(value.get("key") or key),
            "label": _text(value.get("label")),
            "label_en": _text(value.get("label_en")),
            "status": normalized_status,
            "reason": _text(value.get("reason")),
            "score": _finite(value.get("score")),
            "weight": _finite(value.get("weight")),
            "evidence_count": int(value.get("evidence_count") or 0),
            "expected_count": int(value.get("expected_count") or 0),
        })
    return dimensions, counts


def _library_pathway_summary(
    match: Mapping[str, Any], row: Mapping[str, Any]
) -> dict[str, Any]:
    raw = match.get("pathway_json")
    try:
        pathway = json.loads(raw) if isinstance(raw, str) else raw
    except (TypeError, ValueError, json.JSONDecodeError):
        pathway = None
    if not isinstance(pathway, Mapping):
        return {"status": "unavailable", "candidate_conditioned": True}
    recurring_ions = [
        {
            key: item.get(key)
            for key in ("mz", "freq", "n_spectra", "rel_intensity", "annotation")
            if item.get(key) is not None
        }
        for item in (pathway.get("diagnostic_ions") or [])
        if isinstance(item, Mapping) and bool(item.get("recurring"))
    ][:20]
    losses = []
    for item in pathway.get("neutral_losses") or []:
        if not isinstance(item, Mapping) or not bool(item.get("recurring")):
            continue
        entry = {
            key: item.get(key)
            for key in (
                "loss_da", "freq", "n_spectra", "fragment_mz", "rel_intensity",
                "annotation", "nl_key",
            )
            if item.get(key) is not None
        }
        if item.get("annotation"):
            observed = _finite(item.get("loss_da"))
            matches = get_neutral_loss_by_mass(
                observed or -1.0,
                tolerance=0.01,
                classification=_text(match.get("display_class") or match.get("classification")),
                strict_fg=True,
                polarity=_text(row.get("ion_mode")),
                candidate_smiles=_text(match.get("smiles")),
                query_adduct=_specific_adduct(match.get("adduct")),
            )
            wanted = normalize_formula_text(
                next((candidate.formula for candidate in matches if candidate.name_en == item.get("annotation")), "")
            )
            entry["annotation_status"] = (
                "current_policy_supported" if wanted
                else "legacy_annotation_rejected_by_current_policy"
            )
        else:
            entry["annotation_status"] = "mass_difference_only"
        losses.append(entry)
        if len(losses) >= 20:
            break
    return {
        "status": "available" if recurring_ions or losses else "unavailable",
        "candidate_conditioned": True,
        "independent_evidence": False,
        "identity_rank_influence": False,
        "n_reference_spectra": int(pathway.get("n_spectra") or 0),
        "reference_adducts": list(pathway.get("adducts") or []),
        "reference_polarities": list(pathway.get("polarities") or []),
        "recurring_ions": recurring_ions,
        "recurring_losses": losses,
        "loss_chains": list(pathway.get("loss_chains") or [])[:8],
        "interpretation_boundary": "candidate_specific_library_summary_not_confirmed_mechanism",
    }


def _reasoning_strings(details: Mapping[str, Any]) -> list[str]:
    values: list[Any] = []
    v3 = details.get("v3") if isinstance(details.get("v3"), Mapping) else {}
    for candidate in (
        v3.get("reasoning"),
        details.get("chemical_reasoning"),
        details.get("reasoning"),
    ):
        if isinstance(candidate, str):
            values.append(candidate)
        elif isinstance(candidate, list):
            values.extend(candidate)
    output: list[str] = []
    for value in values:
        text = _text(value)[:240]
        if text and text not in output:
            output.append(text)
        if len(output) >= 8:
            break
    return output


def _auxiliary_interpretation(
    match: Mapping[str, Any], row: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    row = row if isinstance(row, Mapping) else {}
    details = match.get("expert_details")
    details = details if isinstance(details, Mapping) else {}
    reasoning = _reasoning_strings(details)
    biosource = [
        value for value in reasoning
        if "生源" in value or "biosynth" in value.lower()
    ][:4]
    neutral = details.get("neutral_loss_analysis")
    neutral = neutral if isinstance(neutral, Mapping) else {}
    path = neutral.get("best_pathway")
    if not isinstance(path, list) or not path:
        connected = neutral.get("connected_pathways")
        if isinstance(connected, list) and connected and isinstance(connected[0], list):
            path = connected[0]
        else:
            path = neutral.get("sequential_losses")
    if not isinstance(path, list):
        path = []
    steps: list[dict[str, Any]] = []
    rejected_steps: list[dict[str, Any]] = []
    for step in path:
        if not isinstance(step, Mapping):
            continue
        accepted, rejected = _validated_step(step, match, row)
        if accepted:
            steps.append(accepted)
        if rejected:
            rejected_steps.append(rejected)
        if len(steps) >= 8:
            break
    neutral_losses = neutral.get("neutral_losses")
    neutral_losses = neutral_losses if isinstance(neutral_losses, list) else []
    validated_neutral_losses: list[dict[str, Any]] = []
    for item in neutral_losses:
        if not isinstance(item, Mapping):
            continue
        accepted, rejected = _validated_step(item, match, row)
        if accepted:
            validated_neutral_losses.append(accepted)
        if rejected:
            rejected_steps.append(rejected)
    diagnostic = details.get("diagnostic_ion_analysis")
    diagnostic = diagnostic if isinstance(diagnostic, Mapping) else {}
    diagnostic_ions = diagnostic.get("matched_ions")
    diagnostic_ions = diagnostic_ions if isinstance(diagnostic_ions, list) else []
    expert_score = _finite(match.get("expert_score"))
    expert_source = _text(match.get("expert_score_source"))
    if not expert_source:
        fusion = details.get("score_fusion")
        if isinstance(fusion, Mapping):
            expert_source = _text(fusion.get("selected_source"))
    dimensions, dimension_counts = _expert_dimensions(details)
    product_structures = _compact_product_structures(details)
    pathway_summary = _library_pathway_summary(match, row)
    fragmentation_status = (
        "product_structure_candidates_available" if product_structures
        else "sample_neutral_loss_path_reconstructed" if steps
        else "mass_difference_evidence_only" if neutral_losses
        else "reference_library_pathway_only" if pathway_summary.get("status") == "available"
        else "unavailable"
    )
    total_dimensions = sum(dimension_counts.values())
    structured_dimensions_available = bool(dimensions)
    expert_status = (
        "unavailable" if not structured_dimensions_available
        else "partially_evaluated" if dimension_counts["unavailable"] > 0
        else "evaluated"
    )
    return {
        "classification": {
            "label": _text(match.get("display_class") or match.get("classification")),
            "source": _text(
                match.get("display_classification_source")
                or match.get("classification_source")
            ),
            "reliable": bool(match.get("classification_reliable")),
            "status": (
                "available"
                if _text(match.get("display_class") or match.get("classification"))
                and _text(match.get("display_class") or match.get("classification")).lower()
                not in {"unknown", "other", "unclassified"}
                else "unavailable"
            ),
        },
        "biosource": biosource_context(match, biosource),
        "expert": {
            "status": expert_status,
            "score": expert_score,
            "source": expert_source or "none",
            "score_semantics": "uncalibrated_secondary_summary_not_probability",
            "calibrated": False,
            "positive_support": any(
                item["status"] == "scored" and item["evidence_count"] > 0
                for item in dimensions
            ),
            "dimensions": dimensions,
            "coverage": {
                **dimension_counts,
                "total": total_dimensions,
            },
            "reasoning": reasoning[:4],
            "availability_reason": (
                None if structured_dimensions_available
                else "no_structured_dimension_assessments"
            ),
            "identity_rank_influence": False,
        },
        "adduct_assignment": _adduct_interpretation(match, row),
        "fragmentation": {
            "status": fragmentation_status,
            "steps": steps,
            "step_count": len(steps),
            "neutral_loss_count": len(neutral_losses),
            "neutral_loss_assignments": validated_neutral_losses[:24],
            "current_policy_assigned_neutral_loss_count": sum(
                bool(item.get("neutral_loss"))
                for item in validated_neutral_losses
            ),
            "diagnostic_ion_count": len(diagnostic_ions),
            "plausible_product_structures": product_structures,
            "product_structure_count": len(product_structures),
            "rejected_legacy_assignments": rejected_steps[:8],
            "candidate_specific_library_pathway": pathway_summary,
            "reason": (
                "reference_pathway_available_but_no_sample_pathway_reconstructed"
                if fragmentation_status == "reference_library_pathway_only"
                else None if fragmentation_status != "unavailable"
                else "no_supported_fragmentation_evidence_was_available"
            ),
            "assignment_semantics": (
                "candidate_structures_and_formulae_are_hypotheses_not_confirmed_products"
            ),
            "mechanism_confirmed": False,
            "identity_rank_influence": False,
        },
    }


def _compact_match(match: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: match.get(key)
        for key in _MATCH_DISPLAY_FIELDS
        if match.get(key) is not None
    }


def _retention_time_minutes(value: Any, unit: Any) -> float | None:
    parsed = _finite(value)
    if parsed is None:
        return None
    normalized_unit = _text(unit).lower().replace(" ", "")
    if normalized_unit in {"s", "sec", "second", "seconds"}:
        return parsed / 60.0
    if normalized_unit in {"ms", "msec", "millisecond", "milliseconds"}:
        return parsed / 60000.0
    return parsed


def _compact_row(row: Mapping[str, Any], result_index: int) -> dict[str, Any]:
    compact = {
        key: row.get(key)
        for key in _ROW_DISPLAY_FIELDS
        if row.get(key) is not None
    }
    retention_time = _retention_time_minutes(
        row.get("retention_time"), row.get("retention_time_unit")
    )
    if retention_time is not None:
        compact["retention_time"] = retention_time
        compact["retention_time_unit"] = "minute"
    compact["result_index"] = int(result_index)
    compact["blank_evidence_category"] = _blank_category(row)
    compact["sample_supported"] = _sample_supported(row)
    return compact


def _representative_payload(
    row: Mapping[str, Any], result_index: int
) -> dict[str, Any]:
    """Keep one display/audit payload per group instead of every full result row."""

    payload = _compact_row(row, result_index)
    match = row.get("match") if isinstance(row.get("match"), Mapping) else {}
    retained_match = {
        key: match.get(key)
        for key in (*_MATCH_DISPLAY_FIELDS, "expert_details")
        if match.get(key) is not None
    }
    payload["match"] = retained_match
    return payload


def _is_released(decision: Mapping[str, Any] | None) -> bool:
    if not isinstance(decision, Mapping):
        return False
    nested = decision.get("decision")
    if isinstance(nested, Mapping):
        return bool(nested.get("released"))
    return bool(decision.get("released"))


def _claim_by_spectrum(path: Path | None) -> dict[str, dict[str, Any]]:
    if path is None or not path.is_file():
        return {}
    try:
        import ijson

        with path.open("rb") as handle:
            return {
                _text(item.get("spectrum_id")): dict(item.get("decision") or {})
                for item in ijson.items(handle, "decisions.item", use_float=True)
                if _text(item.get("spectrum_id"))
            }
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return {
            _text(item.get("spectrum_id")): dict(item.get("decision") or {})
            for item in payload.get("decisions") or []
            if _text(item.get("spectrum_id"))
        }


def _candidate_names(match: Mapping[str, Any], formula: str) -> list[str]:
    values: list[str] = []

    def add(candidate: Mapping[str, Any]) -> None:
        candidate_formula = normalize_formula_text(candidate.get("formula"))
        if formula and candidate_formula and candidate_formula != formula:
            return
        name = _text(candidate.get("compound_name"))
        if name and name.lower() not in {"unknown", "unidentified"} and name not in values:
            values.append(name)

    add(match)
    for candidate in match.get("isomer_set") or []:
        if isinstance(candidate, Mapping):
            add(candidate)
    for candidate in match.get("same_formula_alternatives") or []:
        if isinstance(candidate, Mapping):
            add(candidate)
    for candidate in match.get("top_candidates") or []:
        if isinstance(candidate, Mapping):
            add(candidate)
        if len(values) >= 12:
            break
    return values[:12]


def _candidate_identity(candidate: Mapping[str, Any]) -> str:
    inchikey = _text(candidate.get("inchikey")).upper()
    if inchikey:
        return f"ik::{inchikey}"
    return "name::{}::{}".format(
        _text(candidate.get("compound_name")).casefold(),
        normalize_formula_text(candidate.get("formula")),
    )


def _candidate_member(
    candidate: Mapping[str, Any],
    *,
    rank: int,
    role: str,
    fallback: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a compact, ranked isomer member without implying resolution."""

    fallback = fallback if isinstance(fallback, Mapping) else {}

    def value(key: str) -> Any:
        current = candidate.get(key)
        if key == "rank_score":
            # Legacy ``isomer_set`` summaries wrote 0.0 when candidate-level
            # ranking evidence was not copied.  The retained Top-K row still
            # has the real score (or at minimum its final score), so export
            # that evidence instead of presenting the placeholder as measured.
            current_number = _finite(current)
            if current_number is None or (
                current_number == 0
                and (_finite(fallback.get("rank_score")) or _finite(
                    fallback.get("final_score")
                ))
            ):
                return (
                    fallback.get("rank_score")
                    if fallback.get("rank_score") not in (None, "")
                    else fallback.get("final_score")
                )
        return current if current not in (None, "") else fallback.get(key)

    return {
        key: item
        for key, item in {
            "rank": rank,
            "role": role,
            "is_ranking_lead": role == "ranking_lead",
            "compound_name": _text(value("compound_name")) or "Unnamed candidate",
            "compound_name_cn": _text(value("compound_name_cn")),
            # Keep the formula exactly as stored on the database candidate for
            # display and audit.  Formula normalisation is used only by the
            # grouping predicates.  Replacing a charged database formula such
            # as ``C22H23O11+`` with an empty normalised value would break the
            # atomic candidate identity (name/formula/structure) even though
            # the ranked candidate itself had not changed.
            "formula": _text(value("formula")),
            "inchikey": _text(value("inchikey")),
            "smiles": _text(value("smiles")),
            "adduct": _text(value("adduct")),
            "cosine": _finite(value("cosine")),
            "entropy": _finite(value("entropy")),
            "reverse_match": _finite(value("reverse_match")),
            "mass_ppm": _finite(value("mass_ppm")),
            "final_score": _finite(value("final_score")),
            "rank_score": _finite(value("rank_score")),
            "flagship_entity_rank": int(value("flagship_entity_rank") or 0),
            "flagship_entity_score": _finite(value("flagship_entity_score")),
            "confidence": _text(value("confidence")),
            "classification": _text(
                value("display_class") or value("classification")
            ),
            "dataset_source": _text(value("dataset_source")),
            "lib_spectrum_id": _text(value("lib_spectrum_id")),
            "relationship": (
                "ranking_lead" if role == "ranking_lead"
                else _text(value("relationship")) or "same_formula_competitor"
            ),
            # Every member admitted by ``ranked_candidate_members`` belongs to
            # the unresolved competitive set.  Older artifacts did not save
            # this boolean on Top-K rows, so carrying their missing value as
            # False would contradict the very reason the row is exported.
            "spectrally_competitive": True,
        }.items()
        if item not in (None, "")
    }


def ranked_candidate_members(
    match: Mapping[str, Any], formula: str | None = None
) -> list[dict[str, Any]]:
    """Expose the ranking lead plus every unresolved competitive isomer.

    This is the canonical presentation projection used by the interactive
    result page and every offline export.  Older result artifacts sometimes
    stored ``isomer_set`` as names only, or retained the competitors only in
    ``same_formula_alternatives``/``top_candidates``.  Those representations
    are normalised here so an export never silently loses a candidate that the
    saved analysis declared competitive.
    """

    formula = normalize_formula_text(
        (match.get("isomer_group_formula") or match.get("formula"))
        if formula is None
        else formula
    )

    top_candidates = [
        item for item in (match.get("top_candidates") or [])
        if isinstance(item, Mapping)
    ]
    top_by_identity = {
        _candidate_identity(item): item for item in top_candidates
    }
    top_by_name_formula = {
        (
            _text(item.get("compound_name")).casefold(),
            normalize_formula_text(item.get("formula")),
        ): item
        for item in top_candidates
        if _text(item.get("compound_name"))
    }
    members: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add(candidate: Mapping[str, Any], role: str) -> None:
        candidate_formula = normalize_formula_text(candidate.get("formula"))
        if formula and candidate_formula and candidate_formula != formula:
            return
        candidate_name = _text(candidate.get("compound_name"))
        if (
            not _text(candidate.get("inchikey"))
            and candidate_name
            and any(
                _text(member.get("compound_name")).casefold()
                == candidate_name.casefold()
                and normalize_formula_text(member.get("formula"))
                == candidate_formula
                for member in members
            )
        ):
            return
        identity = _candidate_identity(candidate)
        if identity in seen:
            return
        fallback = top_by_identity.get(identity) or top_by_name_formula.get((
            candidate_name.casefold(), candidate_formula,
        ))
        members.append(_candidate_member(
            candidate,
            rank=len(members) + 1,
            role=role,
            fallback=fallback,
        ))
        seen.add(identity)

    add(match, "ranking_lead")
    for candidate in match.get("isomer_set") or []:
        if isinstance(candidate, Mapping):
            add(candidate, "competitive_isomer")
        elif _text(candidate):
            add(
                {"compound_name": _text(candidate), "formula": formula},
                "competitive_isomer",
            )

    # Compatibility for saved artifacts that predate the explicit
    # ``isomer_set`` payload.  Only alternatives already marked competitive
    # are admitted; if the row itself declares unresolved ambiguity, missing
    # legacy flags are treated conservatively as competitive.
    unresolved = bool(match.get("isomer_ambiguous")) or _text(
        match.get("identification_scope")
    ).lower() in {
        "isomer_group",
        "molecular_formula_isomer_group",
        "stereoisomer_group",
        "unresolved_precursor_event",
    }
    for candidate in match.get("same_formula_alternatives") or []:
        if not isinstance(candidate, Mapping):
            continue
        if bool(candidate.get("spectrally_competitive")) or unresolved:
            add(candidate, "competitive_isomer")

    for candidate_name in (
        (match.get("isomer_group_evidence") or {}).get("candidate_names")
        if isinstance(match.get("isomer_group_evidence"), Mapping)
        else []
    ) or []:
        normalized_name = _text(candidate_name)
        if (
            unresolved
            and normalized_name
            and not any(
                _text(member.get("compound_name")).casefold()
                == normalized_name.casefold()
                for member in members
            )
        ):
            add(
                {"compound_name": normalized_name, "formula": formula},
                "competitive_isomer",
            )

    # Final legacy fallback: an unresolved row can only have been assigned
    # that scope from same-formula competitors.  Preserve those Top-K members
    # when neither of the newer explicit collections survived serialization.
    if unresolved and len(members) == 1:
        for candidate in top_candidates:
            if bool(candidate.get("spectrally_competitive")):
                add(candidate, "competitive_isomer")
    return members


# Private compatibility alias for older internal imports.
_ranked_candidate_members = ranked_candidate_members


def _peak_clusters(observations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    clusters: list[dict[str, Any]] = []
    ordered = sorted(
        observations,
        key=lambda item: (
            _finite(item.get("retention_time"))
            if _finite(item.get("retention_time")) is not None else math.inf,
            _finite(item.get("precursor_mz"))
            if _finite(item.get("precursor_mz")) is not None else math.inf,
            int(item.get("result_index") or 0),
        ),
    )
    for observation in ordered:
        rt = _finite(observation.get("retention_time"))
        mz = _finite(observation.get("precursor_mz"))
        target = None
        if rt is not None and mz is not None:
            for cluster in reversed(clusters):
                cluster_rt = _finite(cluster.get("retention_time"))
                cluster_mz = _finite(cluster.get("precursor_mz"))
                if cluster_rt is None or cluster_mz is None:
                    continue
                mz_tolerance = max(0.02, abs(mz) * 20e-6)
                if abs(rt - cluster_rt) <= 0.50 and abs(mz - cluster_mz) <= mz_tolerance:
                    target = cluster
                    break
        if target is None:
            clusters.append({
                **observation,
                "scan_count": 1,
                "result_indices": [int(observation["result_index"])],
                "spectrum_ids": [_text(observation.get("spectrum_id"))],
            })
            continue
        target["scan_count"] += 1
        target["result_indices"].append(int(observation["result_index"]))
        spectrum_id = _text(observation.get("spectrum_id"))
        if spectrum_id:
            target["spectrum_ids"].append(spectrum_id)
        if (_finite(observation.get("tic")) or 0.0) > (_finite(target.get("tic")) or 0.0):
            preserved = {
                "scan_count": target["scan_count"],
                "result_indices": target["result_indices"],
                "spectrum_ids": target["spectrum_ids"],
            }
            target.clear()
            target.update(observation)
            target.update(preserved)
    return clusters


def _conclusion(group: dict[str, Any]) -> dict[str, Any]:
    observations = group["observations"]
    sample_rows = [item for item in observations if item["sample_supported"]]
    considered = sample_rows if sample_rows else observations
    confirmed = any(item.get("released") for item in considered)
    connectivity_high = any(
        _text(item.get("two_dimensional_connectivity_confidence")) == "High"
        and not item.get("isomer_group")
        for item in considered
    )
    exact_full_identity_high = any(
        _text(item.get("exact_structure_confidence")) == "High"
        and not item.get("isomer_group")
        for item in considered
    )
    isomer_strength = max(
        (
            _CONFIDENCE_RANK.get(_text(item.get("isomer_group_confidence")), 0)
            for item in considered if item.get("isomer_group")
        ),
        default=0,
    )
    isomer_confidence = next(
        (
            label for label, rank in sorted(
                _CONFIDENCE_RANK.items(), key=lambda item: item[1], reverse=True
            )
            if rank == isomer_strength
        ),
        "Unknown",
    )
    scope_high = bool(
        connectivity_high
        or isomer_strength >= _CONFIDENCE_RANK["High"]
    )
    medium = any(
        _text(item.get("confidence")) in {"High", "Medium"}
        or _text(item.get("evidence_tier")) == "Strong Candidate"
        for item in considered
    )
    background_only = not sample_rows
    unverified_blank_association = any(
        item.get("blank_evidence_category") in _BACKGROUND_CATEGORIES
        and _text(item.get("blank_evidence_role")).lower() == "unspecified"
        for item in considered
    )
    unresolved_event = group["identity_scope"] == "unresolved_precursor_event"
    if unresolved_event:
        key = "unresolved_event"
        label = "未解析母离子事件"
        if "formula_assignment_incomplete" in set(
            group.get("precursor_event_resolution_statuses") or []
        ):
            note = (
                "同一母离子事件的候选结构不一致，且至少一个候选缺少可用分子式；"
                "保留候选供复核，但不将该事件计为独立化合物"
            )
        else:
            note = (
                "同一母离子事件的候选分子式互相冲突；保留候选供复核，"
                "但不将该事件计为独立化合物"
            )
    elif confirmed:
        key = "confirmed"
        label = "已确证存在"
        note = "已通过当前项目的正式结论条件"
    elif background_only:
        key = "background"
        label = "空白／背景相关"
        note = (
            "谱库有命中，但当前信号与任务所选空白证据相关，不计作样本化合物；"
            "若空白仅按极性自动分配，还需核验其批次与色谱方法可比性"
        )
    elif connectivity_high:
        key = "high"
        label = "高可信单一二维连接结构检出"
        note = (
            "样本信号在单一二维连接结构范围具有 High 级谱库支持；"
            + (
                "完整结构身份也达到 High，但尚未达到标准品确证"
                if exact_full_identity_high
                else "立体化学／完整结构定义层仍未解析，且尚未达到标准品确证"
            )
        )
    elif group["identity_scope"] == "isomer_group" and isomer_strength >= 3:
        key = "isomer"
        if isomer_strength >= _CONFIDENCE_RANK["High"]:
            label = "高可信异构体组检出"
            note = "High 级证据支持该分子式异构体组，但不能点名单一异构体"
        else:
            label = "中等证据异构体组"
            note = "中等证据支持该分子式异构体组，但不能点名单一异构体"
    elif medium:
        key = "candidate"
        label = "中等证据候选"
        note = "样本信号具有中等谱库支持；身份范围维持为候选，需补充正交证据"
    else:
        key = "insufficient"
        label = "证据不足"
        note = "仅保留为低强度谱库线索，不计入主要鉴定结果"
    blank_warning = None
    if unverified_blank_association and key != "background":
        blank_warning = (
            "同极性未定角色空白中亦观察到该信号，但批次／色谱方法配对未经核验；"
            "因此只作背景警示，不据此剔除样本候选"
        )
        note = f"{note}；{blank_warning}"
    return {
        "key": key,
        "label": label,
        "note": note,
        "confirmed_present": key == "confirmed",
        "counted_as_sample_compound": key in {"confirmed", "high", "isomer", "candidate"},
        "high_confidence_detected": key == "high",
        "high_confidence_2d_connectivity_detected": (
            connectivity_high and not unresolved_event
        ),
        "high_confidence_exact_structure_detected": (
            exact_full_identity_high and not unresolved_event
        ),
        "high_confidence_group_detected": scope_high and not unresolved_event,
        "isomer_group_only": key == "isomer",
        "isomer_group_confidence": (
            isomer_confidence if key == "isomer" else None
        ),
        "unverified_blank_association": unverified_blank_association,
        "blank_warning": blank_warning,
        "standard_confirmation_required": key != "confirmed",
    }


def build_compound_conclusions(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    claim_decisions_path: Path | None = None,
    source_artifact: Path | None = None,
    primary_identity_artifact: Mapping[str, Any] | Path | None = None,
    feature_inventory: Iterable[Mapping[str, Any]] | None = None,
    feature_inventory_source: Path | None = None,
) -> dict[str, Any]:
    """Build a deterministic, compact, compound-level result view."""

    source_rows_materialized = rows if isinstance(rows, list) else list(rows)
    primary_index, primary_metadata, primary_path = _primary_entity_index(
        task_id, primary_identity_artifact
    )
    high_evidence_by_entity = {}
    for entity in primary_index.values():
        entity_id = _text(entity.get("entity_id"))
        if entity_id not in high_evidence_by_entity:
            high_evidence_by_entity[entity_id] = single_energy_high_evidence(
                entity, source_rows_materialized
            )
    materialized_rows = [
        _project_primary_identity(
            row,
            primary_index.get(result_index),
            high_evidence_by_entity.get(_text((primary_index.get(result_index) or {}).get("entity_id"))),
            result_index=result_index,
        )
        if isinstance(row, Mapping)
        else row
        for result_index, row in enumerate(source_rows_materialized)
    ]
    feature_rows: dict[str, list[tuple[int, Mapping[str, Any]]]] = {}
    for result_index, row in enumerate(materialized_rows):
        if not isinstance(row, Mapping):
            continue
        feature_rows.setdefault(
            _analysis_feature_id(row, result_index), []
        ).append((result_index, row))
    inventory_rows, inventory_feature_ids = _materialize_feature_inventory(
        feature_inventory
    )
    result_lcms_feature_ids = {
        feature_id for feature_id in feature_rows if feature_id.startswith("lcms:")
    }
    for inventory_index, row in enumerate(inventory_rows):
        feature_id = _analysis_feature_id(row)
        if feature_id not in feature_rows:
            feature_rows[feature_id] = [
                (len(materialized_rows) + inventory_index, row)
            ]
    event_plans, event_aggregation = _precursor_event_plans(materialized_rows)
    claims = _claim_by_spectrum(claim_decisions_path)
    groups: dict[str, dict[str, Any]] = {}
    source_rows = 0
    matched_rows = 0
    source_digest = hashlib.sha256()
    legacy_high_identities: set[str] = set()

    for result_index, row in enumerate(materialized_rows):
        if not isinstance(row, Mapping):
            continue
        source_rows += 1
        source_match = row.get("match")
        if not isinstance(source_match, Mapping):
            continue
        # Old saved tasks may predate the current exact-name and compound-role
        # catalogs.  Re-project those presentation fields while building the
        # compound view so registry identifiers and reviewed external
        # surfactants do not leak into new exports.  The saved analysis and all
        # numeric scores remain untouched.
        match = deepcopy(dict(source_match))
        normalize_match_compound_names(match)
        annotate_match_compound_roles(match)
        if _text(source_match.get("primary_identity_source")) == "v15_flagship":
            # Name and structure fields are one atomic database candidate.  A
            # presentation normalizer may classify or flag that candidate, but
            # it cannot replace the V15 database name or splice legacy identity
            # fields back into the selected structure.
            for field in _PRIMARY_CANDIDATE_FIELDS:
                if field in source_match:
                    match[field] = deepcopy(source_match[field])
                else:
                    match.pop(field, None)
            for field in (
                "top_candidates",
                "isomer_set",
                "same_formula_alternatives",
                "same_formula_alternatives_n",
                "flagship_entity_score",
                "flagship_entity_rank",
                "rank_score",
                "final_score",
                "primary_identity_source",
                "primary_entity_id",
                "primary_ranking_identity_key",
                "primary_decision_status",
                "primary_decision_scope",
                "primary_identity_resolution",
                "entity_high_evidence",
                "confidence",
                "exact_structure_confidence",
                "two_dimensional_connectivity_confidence",
                "isomer_group_confidence",
                "identification_scope",
                "isomer_group_formula",
                "isomer_group_size",
                "isomer_ambiguous",
                "isomer_ambiguity_type",
            ):
                if field in source_match:
                    match[field] = deepcopy(source_match[field])
        # Roles describe the currently selected identity; a stale primary
        # artifact must not restore a pre-curation role over this projection.
        apply_report_compound_roles(match)
        normalized_row = dict(row)
        normalized_row["match"] = match
        if _text(match.get("dataset_source")).lower() == "iceberg_fallback":
            continue
        matched_rows += 1
        exact_identity, _ = connectivity_identity(match)
        if _text(match.get("confidence")) == "High":
            legacy_high_identities.add(exact_identity)
        candidate_identity_key, identity_scope, identity_source = _group_key(
            match
        )
        row_formula = normalize_formula_text(
            match.get("isomer_group_formula") or match.get("formula")
        )
        formula = row_formula
        event_id = _analysis_feature_id(normalized_row, result_index)
        event_plan = event_plans.get(event_id) if event_id else None
        event_status = (
            _text(event_plan.get("status"))
            if isinstance(event_plan, Mapping)
            else ""
        )
        if event_status in {
            "same_formula_isomer_group",
            "cross_formula_conflict",
            "formula_assignment_incomplete",
        }:
            candidate_identity_key = _text(event_plan.get("group_key"))
            identity_scope = _text(event_plan.get("identity_scope"))
            identity_source = _text(event_plan.get("identity_source"))
            formula = normalize_formula_text(event_plan.get("formula"))
        compact_row = _compact_row(normalized_row, result_index)
        decision = claims.get(_text(normalized_row.get("spectrum_id")))
        observation = {
            **compact_row,
            "confidence": _text(match.get("confidence")) or "Unknown",
            "exact_structure_confidence": _text(
                match.get("exact_structure_confidence")
            ) or _text(match.get("confidence")) or "Unknown",
            "two_dimensional_connectivity_confidence": _text(
                match.get("two_dimensional_connectivity_confidence")
            ) or _text(match.get("confidence")) or "Unknown",
            "isomer_group_confidence": _text(match.get("isomer_group_confidence")),
            "evidence_tier": _text(match.get("evidence_tier")),
            "identification_scope": _text(match.get("identification_scope")),
            "isomer_ambiguity_type": _text(match.get("isomer_ambiguity_type")),
            "isomer_group": _isomer_group_row(match),
            "released": _is_released(decision),
            "precursor_event_resolution_status": (
                event_status or "not_available"
            ),
        }
        if event_status == "same_formula_isomer_group":
            # Conflicting structures from the same event cannot retain a
            # scan-level exact/connectivity claim.  The strongest individual
            # spectrum may support the formula-isomer group, but repeated
            # collision energies are not added as independent votes.
            observation.update({
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": _text(
                    event_plan.get("scope_confidence")
                ) or "Unknown",
                "identification_scope": "isomer_group",
                "isomer_ambiguity_type": (
                    "precursor_event_cross_energy_structure_conflict"
                ),
                "isomer_group": True,
                "released": False,
            })
        elif event_status in {
            "cross_formula_conflict",
            "formula_assignment_incomplete",
        }:
            observation.update({
                "exact_structure_confidence": "Unknown",
                "two_dimensional_connectivity_confidence": "Unknown",
                "isomer_group_confidence": "Unknown",
                "identification_scope": "unresolved_precursor_event",
                "isomer_ambiguity_type": (
                    "precursor_event_cross_formula_conflict"
                    if event_status == "cross_formula_conflict"
                    else "precursor_event_formula_assignment_incomplete"
                ),
                "isomer_group": False,
                "released": False,
            })
        feature_group_key = f"ion-feature:{event_id}"
        group = groups.setdefault(feature_group_key, {
            "key": feature_group_key,
            "feature_id": event_id,
            "identity_key": candidate_identity_key,
            "identity_scope": identity_scope,
            "identity_source": identity_source,
            "formula": formula,
            "representative_row": None,
            "representative_priority": None,
            "observations": [],
            "candidate_names": [],
            "candidate_member_pool": {},
            "representative_member_order": [],
            "same_formula_alternative_count": 0,
            "biosynthetic_pathways": [],
            "expert_evaluated": False,
            "expert_positive": False,
            "fragmentation": None,
            "precursor_event_ids": set(),
            "precursor_event_resolution_statuses": set(),
            "candidate_formulas": set(),
        })
        source_precursor_event_id = _text(
            normalized_row.get("precursor_feature_id")
        )
        if source_precursor_event_id:
            group["precursor_event_ids"].add(source_precursor_event_id)
        elif event_id:
            group["precursor_event_ids"].add(event_id)
        if event_status:
            group["precursor_event_resolution_statuses"].add(event_status)
        candidate_formulas = (
            list(event_plan.get("candidate_formulas") or [])
            if isinstance(event_plan, Mapping)
            else []
        )
        if row_formula:
            candidate_formulas.append(row_formula)
        group["candidate_formulas"].update(
            normalize_formula_text(value)
            for value in candidate_formulas
            if normalize_formula_text(value)
        )
        representative_priority = (
            -int(bool(observation["released"])),
            *_result_priority(normalized_row),
        )
        if (
            group["representative_priority"] is None
            or representative_priority < group["representative_priority"]
        ):
            group["representative_row"] = _representative_payload(
                normalized_row, result_index
            )
            group["representative_priority"] = representative_priority
            group["same_formula_alternative_count"] = int(
                match.get("same_formula_alternatives_n")
                or len(match.get("same_formula_alternatives") or [])
            )
            member_formula = "" if identity_scope == (
                "unresolved_precursor_event"
            ) else formula
            group["representative_member_order"] = [
                _candidate_identity(member)
                for member in ranked_candidate_members(match, member_formula)
            ]
        member_formula = "" if identity_scope == (
            "unresolved_precursor_event"
        ) else formula
        row_candidate_members = ranked_candidate_members(match, member_formula)
        for member in row_candidate_members:
            identity = _candidate_identity(member)
            existing_member = group["candidate_member_pool"].get(identity)
            if (
                existing_member is None
                or (_finite(member.get("rank_score")) or -1.0)
                > (_finite(existing_member.get("rank_score")) or -1.0)
            ):
                group["candidate_member_pool"][identity] = member
        group["same_formula_alternative_count"] = max(
            int(group.get("same_formula_alternative_count") or 0),
            int(
                match.get("same_formula_alternatives_n")
                or len(match.get("same_formula_alternatives") or [])
            ),
        )
        group["observations"].append(observation)
        for name in _candidate_names(match, member_formula):
            if name not in group["candidate_names"]:
                group["candidate_names"].append(name)
        auxiliary = _auxiliary_interpretation(match, normalized_row)
        for pathway in auxiliary["biosource"]["pathways"]:
            if pathway not in group["biosynthetic_pathways"]:
                group["biosynthetic_pathways"].append(pathway)
        group["expert_evaluated"] = group["expert_evaluated"] or (
            auxiliary["expert"]["status"] in {"evaluated", "partially_evaluated"}
        )
        group["expert_positive"] = group["expert_positive"] or bool(
            auxiliary["expert"]["positive_support"]
        )
        current_fragmentation = group.get("fragmentation") or {}
        current_strength = (
            int(current_fragmentation.get("product_structure_count") or 0) * 100
            + int(current_fragmentation.get("step_count") or 0) * 10
            + int((current_fragmentation.get("candidate_specific_library_pathway") or {}).get("status") == "available")
        )
        new_strength = (
            int(auxiliary["fragmentation"].get("product_structure_count") or 0) * 100
            + int(auxiliary["fragmentation"].get("step_count") or 0) * 10
            + int((auxiliary["fragmentation"].get("candidate_specific_library_pathway") or {}).get("status") == "available")
        )
        if new_strength > current_strength:
            group["fragmentation"] = auxiliary["fragmentation"]
        source_digest.update(_canonical({
            "result_index": result_index,
            "spectrum_id": normalized_row.get("spectrum_id"),
            "group_key": feature_group_key,
            "library_identity_key": candidate_identity_key,
            "precursor_feature_id": event_id,
            "precursor_event_resolution_status": event_status,
            "confidence": match.get("confidence"),
            "isomer_group_confidence": match.get("isomer_group_confidence"),
            "rank_score": match.get("rank_score"),
            "blank_category": compact_row["blank_evidence_category"],
        }))
        source_digest.update(b"\n")

    compounds: list[dict[str, Any]] = []
    for key, group in groups.items():
        representative_row = group.pop("representative_row")
        group.pop("representative_priority", None)
        if not isinstance(representative_row, Mapping):
            continue
        representative_index = int(representative_row.get("result_index") or 0)
        representative_match = representative_row.get("match") or {}
        compact_match = _compact_match(representative_match)
        auxiliary = _auxiliary_interpretation(representative_match, representative_row)
        auxiliary["biosource"] = biosource_context(
            representative_match,
            group["biosynthetic_pathways"][:6],
        )
        auxiliary["expert"]["status"] = (
            auxiliary["expert"]["status"] if group["expert_evaluated"] else "unavailable"
        )
        auxiliary["expert"]["positive_support"] = group["expert_positive"]
        if group.get("fragmentation"):
            auxiliary["fragmentation"] = group["fragmentation"]
        flags = {
            "has_reasoning": bool(auxiliary["expert"].get("reasoning")),
            "has_fragmentation": bool(
                int(auxiliary["fragmentation"].get("step_count") or 0)
                or int(auxiliary["fragmentation"].get("product_structure_count") or 0)
            ),
            "has_fragmentation_evidence": bool(
                int(auxiliary["fragmentation"].get("step_count") or 0)
                or int(auxiliary["fragmentation"].get("product_structure_count") or 0)
                or int(auxiliary["fragmentation"].get("neutral_loss_count") or 0)
                or int(auxiliary["fragmentation"].get("diagnostic_ion_count") or 0)
            ),
            "has_library_pathway_reference": bool(
                (auxiliary["fragmentation"].get(
                    "candidate_specific_library_pathway"
                ) or {}).get("status") == "available"
            ),
            "has_neutral_loss": int(auxiliary["fragmentation"].get("neutral_loss_count") or 0) > 0,
            "fragmentation_step_count": int(auxiliary["fragmentation"].get("step_count") or 0),
            "neutral_loss_count": int(auxiliary["fragmentation"].get("neutral_loss_count") or 0),
            "diagnostic_ion_count": int(auxiliary["fragmentation"].get("diagnostic_ion_count") or 0),
            "product_structure_count": int(auxiliary["fragmentation"].get("product_structure_count") or 0),
            "sample_ion_form": (auxiliary.get("adduct_assignment") or {}).get("sample_ion_form"),
            "sample_ion_form_source": (auxiliary.get("adduct_assignment") or {}).get("sample_ion_form_source"),
            "biosynthetic_pathways": list(auxiliary["biosource"].get("pathways") or []),
            "pathway_status": str(
                (auxiliary["fragmentation"].get(
                    "candidate_specific_library_pathway"
                ) or {}).get("status") or "unavailable"
            ),
            "biosource_status": auxiliary["biosource"]["status"],
        }
        compact_match["display_evidence_flags"] = flags
        compact_match["auxiliary_interpretation"] = auxiliary
        compact_match["identity_key"] = group["identity_key"]
        claim = claims.get(_text(representative_row.get("spectrum_id")))
        if claim:
            compact_match["final_claim_decision"] = claim
        compact_rep = {
            key_name: value
            for key_name, value in representative_row.items()
            if key_name != "match"
        }
        compact_rep["match"] = compact_match
        conclusion = _conclusion(group)
        category_counts = Counter(
            item["blank_evidence_category"] for item in group["observations"]
        )
        ratios = [
            _finite(item.get("blank_sample_ratio"))
            for item in group["observations"]
        ]
        ratios = [value for value in ratios if value is not None]
        clusters = _peak_clusters([
            {
                key_name: item.get(key_name)
                for key_name in (
                    "result_index", "spectrum_id", "precursor_mz",
                    "retention_time", "retention_time_unit", "tic", "ion_mode",
                    "collision_energy", "collision_energy_unit",
                    "blank_evidence_category", "sample_supported",
                )
            }
            for item in group["observations"]
        ])
        # ``compound_name`` is the normalized database identity selected by the
        # matching layer.  Language aliases such as ``compound_name_cn`` remain
        # available as provenance metadata, but must never replace or decorate
        # the analytical result name.
        unresolved_event = group["identity_scope"] == (
            "unresolved_precursor_event"
        )
        display_name = _text(compact_match.get("compound_name")) or "Unnamed candidate"
        member_pool = dict(group.get("candidate_member_pool") or {})
        member_formula = "" if unresolved_event else group["formula"]
        representative_members = ranked_candidate_members(
            representative_match,
            member_formula,
        )
        # A connectivity can recur in several precursor entities while the
        # database records retain different full InChIKeys or structure
        # strings.  The report row is represented by one selected observation;
        # its displayed name/formula/InChIKey/SMILES must therefore come from
        # that same observation, even when another occurrence has a larger
        # score.  Other occurrences remain in ``occurrences`` and never
        # overwrite the representative candidate's atomic identity fields.
        for member in representative_members:
            member_pool[_candidate_identity(member)] = member
        representative_order = [
            _candidate_identity(member) for member in representative_members
            if _candidate_identity(member) in member_pool
        ]
        remaining_member_ids = sorted(
            (identity for identity in member_pool if identity not in representative_order),
            key=lambda identity: (
                int(member_pool[identity].get("flagship_entity_rank") or 10**9),
                -(_finite(member_pool[identity].get("rank_score")) or -1.0),
                _text(member_pool[identity].get("compound_name")).casefold(),
                identity,
            ),
        )
        candidate_members = []
        for member_index, identity in enumerate(
            representative_order + remaining_member_ids
        ):
            member = dict(member_pool[identity])
            member["rank"] = member_index + 1
            if member_index == 0:
                member["role"] = "ranking_lead"
            elif unresolved_event:
                member["role"] = "conflicting_candidate"
            else:
                member["role"] = "competitive_isomer"
            member["is_ranking_lead"] = member_index == 0
            if member_index == 0:
                member["relationship"] = "ranking_lead"
            elif unresolved_event:
                if not _text(member.get("formula")):
                    member["relationship"] = "candidate_missing_formula_assignment"
                elif "formula_assignment_incomplete" in set(
                    group.get("precursor_event_resolution_statuses") or []
                ):
                    member["relationship"] = "candidate_formula_assignment_available"
                else:
                    member["relationship"] = "conflicting_formula_candidate"
            apply_report_compound_roles(member)
            member["plant_report"] = plant_report_selection(member)
            candidate_members.append(member)
        lead_candidate = candidate_members[0] if candidate_members else _candidate_member(
            representative_match,
            rank=1,
            role="ranking_lead",
        )
        if group["identity_scope"] == "isomer_group":
            lead_name = _text(lead_candidate.get("compound_name"))
            # Group scope and member count already have dedicated structured
            # fields.  Concatenating localized prose into ``display_name`` made
            # the scientific value change when the interface language changed.
            display_name = lead_name or display_name
        elif unresolved_event:
            # No database candidate is promoted to the event identity when the
            # neutral-formula assignments conflict.  Exact source names remain
            # available in ``candidate_members`` for audit and manual review.
            display_name = "Unresolved precursor event"
        precursor_event_ids = sorted(group["precursor_event_ids"])
        event_resolution_statuses = sorted(
            group["precursor_event_resolution_statuses"]
        )
        candidate_formulas = sorted(group["candidate_formulas"])
        plant_report = plant_report_selection(compact_match)
        excluded_members = [
            member for member in candidate_members
            if member["plant_report"]["excluded_from_plant_candidate_table"]
        ]
        plant_report["excluded_candidate_member_count"] = len(excluded_members)
        plant_report["retained_candidate_member_count"] = len(candidate_members) - len(excluded_members)
        plant_report["mixed_candidate_roles"] = bool(
            excluded_members and len(excluded_members) < len(candidate_members)
        )
        # If an unresolved competitive group contains reviewed external
        # identities, retain the whole original group in review. Filtering its
        # members must not silently promote another winner or resolve isomers.
        if excluded_members:
            plant_report["excluded_from_plant_candidate_table"] = True
            if plant_report["mixed_candidate_roles"]:
                plant_report["review_category"] = "mixed_candidate_roles"
                plant_report["reason"] = "unresolved_group_contains_external_candidate"
        conclusion["counted_as_plant_candidate"] = bool(
            conclusion["counted_as_sample_compound"]
            and not plant_report["excluded_from_plant_candidate_table"]
        )
        compounds.append({
            "plant_report": plant_report,
            # Keep the chemical identity stable across runs/features. The
            # feature-scoped result has its own result_id below.
            "compound_id": "compound:" + hashlib.sha256(
                group["identity_key"].encode("utf-8")
            ).hexdigest()[:18],
            "result_id": "feature-result:" + hashlib.sha256(
                f"{task_id}|{group['feature_id']}".encode("utf-8")
            ).hexdigest()[:18],
            "run_id": task_id,
            "feature_id": group["feature_id"],
            "record_unit": "ion_feature",
            "identity_key": group["identity_key"],
            "identity_scope": group["identity_scope"],
            "identity_source": group["identity_source"],
            "record_type": (
                "unresolved_precursor_event"
                if unresolved_event
                else "compound_or_isomer_group"
            ),
            "display_name": display_name,
            "formula": group["formula"],
            "lead_candidate": lead_candidate,
            "lead_candidate_semantics": (
                "event_ranking_lead_not_reportable_compound_identity"
                if unresolved_event
                else (
                    "current_injection_ranking_lead_not_resolved_single_structure"
                    if group["identity_scope"] == "isomer_group"
                    else "current_injection_ranking_lead"
                )
            ),
            "candidate_members": candidate_members,
            "candidate_member_count": len(candidate_members),
            "other_same_formula_alternative_count": (
                0
                if unresolved_event
                else max(
                    0,
                    int(group.get("same_formula_alternative_count") or 0)
                    - max(0, len(candidate_members) - 1),
                )
            ),
            "candidate_names": group["candidate_names"][:12],
            "candidate_formulas": candidate_formulas,
            "unresolved_event_reason": (
                "incomplete_candidate_formula_assignment"
                if unresolved_event and "formula_assignment_incomplete" in set(
                    group.get("precursor_event_resolution_statuses") or []
                )
                else (
                    "conflicting_candidate_formulae" if unresolved_event else None
                )
            ),
            "precursor_event_ids": precursor_event_ids,
            "precursor_event_count": len(precursor_event_ids),
            "precursor_event_resolution_statuses": (
                event_resolution_statuses
            ),
            "conclusion": conclusion,
            "display": {
                "tier": (
                    "library_released"
                    if conclusion["key"] == "confirmed"
                    else "library_review"
                ),
                "primary_candidate_id": (
                    lead_candidate.get("library_identity_ref")
                    or lead_candidate.get("lib_spectrum_id")
                    or lead_candidate.get("inchikey")
                    or None
                ),
                "default_visible": True,
                "warnings": [],
            },
            "representative": compact_rep,
            "spectra_count": len(group["observations"]),
            "peak_cluster_count": len(clusters),
            "occurrences": clusters,
            "sample_evidence": {
                "status": (
                    "sample_supported"
                    if any(item["sample_supported"] for item in group["observations"])
                    else "background_related"
                ),
                "category_counts": dict(sorted(category_counts.items())),
                "sample_supported_spectra": sum(
                    bool(item["sample_supported"]) for item in group["observations"]
                ),
                "background_related_spectra": sum(
                    item["blank_evidence_category"] in _BACKGROUND_CATEGORIES
                    for item in group["observations"]
                ),
                "maximum_sample_to_blank_ratio": max(ratios) if ratios else None,
                "blank_evidence_available": any(
                    item["blank_evidence_category"] != "not_assessed"
                    for item in group["observations"]
                ),
            },
            "auxiliary_interpretation": auxiliary,
        })

    # Attach every observation to its feature result, including scans with no
    # library candidate. This preserves the raw evidence lineage while keeping
    # one primary row per measured feature.
    matched_feature_ids = {item["feature_id"] for item in compounds}
    for item in compounds:
        indexed_rows = feature_rows.get(item["feature_id"], [])
        all_observations = [
            _compact_row(row, result_index)
            for result_index, row in indexed_rows
        ]
        item["spectra_count"] = len(all_observations)
        item["occurrences"] = _peak_clusters(all_observations)
        item["peak_cluster_count"] = len(item["occurrences"])
        item["lineage"] = {
            "previous_scan_record_ids": [
                _text(row.get("spectrum_id"))
                for _, row in indexed_rows
                if _text(row.get("spectrum_id"))
            ],
            "candidate_names_used_for_feature_grouping": False,
        }

    for feature_id, indexed_rows in feature_rows.items():
        if feature_id in matched_feature_ids:
            continue
        ms2_spectrum_count = sum(
            bool(_text(row.get("spectrum_id"))) for _, row in indexed_rows
        )
        unresolved_reason = (
            "no_library_candidate"
            if ms2_spectrum_count
            else "no_associated_ms2"
        )
        ordered_rows = sorted(
            indexed_rows,
            key=lambda item: (
                -(_finite(item[1].get("tic")) or 0.0),
                item[0],
            ),
        )
        representative_index, representative_row = ordered_rows[0]
        compact_rep = _compact_row(representative_row, representative_index)
        observations = [
            _compact_row(row, result_index)
            for result_index, row in indexed_rows
        ]
        occurrences = _peak_clusters(observations)
        result_hash = hashlib.sha256(
            f"{task_id}|{feature_id}".encode("utf-8")
        ).hexdigest()[:18]
        compounds.append({
            "plant_report": {
                "excluded_from_plant_candidate_table": False,
                "review_category": "unannotated_feature",
                "reason": unresolved_reason,
            },
            "compound_id": f"compound:{result_hash}",
            "result_id": f"feature-result:{result_hash}",
            "run_id": task_id,
            "feature_id": feature_id,
            "record_unit": "ion_feature",
            "identity_key": f"unannotated-feature:{feature_id}",
            "identity_scope": "unannotated_feature",
            "identity_source": unresolved_reason,
            "record_type": "unannotated_ion_feature",
            "display_name": "Unannotated ion feature",
            "formula": "",
            "lead_candidate": {},
            "lead_candidate_semantics": unresolved_reason,
            "candidate_members": [],
            "candidate_member_count": 0,
            "other_same_formula_alternative_count": 0,
            "candidate_names": [],
            "candidate_formulas": [],
            "unresolved_event_reason": unresolved_reason,
            "precursor_event_ids": [feature_id],
            "precursor_event_count": 1,
            "precursor_event_resolution_statuses": [unresolved_reason],
            "conclusion": {
                "key": "unannotated",
                "label": "未注释特征",
                "note": (
                    "已检出离子特征，但没有关联的MS2谱"
                    if unresolved_reason == "no_associated_ms2"
                    else "已检出离子特征，但普通谱库检索没有可用结构候选"
                ),
                "confirmed_present": False,
                "counted_as_sample_compound": False,
                "counted_as_plant_candidate": False,
                "high_confidence_detected": False,
                "high_confidence_2d_connectivity_detected": False,
                "high_confidence_exact_structure_detected": False,
                "high_confidence_group_detected": False,
                "isomer_group_only": False,
                "isomer_group_confidence": None,
                "unverified_blank_association": False,
                "blank_warning": None,
                "standard_confirmation_required": True,
            },
            "display": {
                "tier": "unannotated",
                "primary_candidate_id": None,
                "default_visible": True,
                "warnings": [unresolved_reason],
            },
            "representative": compact_rep,
            "spectra_count": ms2_spectrum_count,
            "ms1_feature_observation_count": int(
                any(bool(row.get("lcms_feature_inventory_only")) for _, row in indexed_rows)
            ),
            "peak_cluster_count": len(occurrences),
            "occurrences": occurrences,
            "sample_evidence": {
                "status": "not_assessed",
                "category_counts": dict(Counter(
                    item.get("blank_evidence_category") or "not_assessed"
                    for item in observations
                )),
                "sample_supported_spectra": ms2_spectrum_count,
                "background_related_spectra": 0,
                "maximum_sample_to_blank_ratio": None,
                "blank_evidence_available": any(
                    item.get("blank_evidence_category") != "not_assessed"
                    for item in observations
                ),
            },
            "auxiliary_interpretation": {},
            "lineage": {
                "previous_scan_record_ids": [
                    _text(row.get("spectrum_id"))
                    for _, row in indexed_rows
                    if _text(row.get("spectrum_id"))
                ],
                "candidate_names_used_for_feature_grouping": False,
            },
        })

    status_priority = {
        "confirmed": 0,
        "high": 1,
        "isomer": 2,
        "candidate": 3,
        "background": 4,
        "unresolved_event": 5,
        "insufficient": 6,
        "unannotated": 7,
    }
    compounds.sort(key=lambda item: (
        status_priority.get(item["conclusion"]["key"], 99),
        int(bool(
            ((item.get("representative") or {}).get("match") or {}).get(
                "potential_contaminant"
            )
        )),
        -(_finite((item["representative"].get("match") or {}).get("rank_score")) or 0.0),
        int(
            (item["representative"].get("match") or {}).get(
                "flagship_entity_rank"
            )
            or 10**9
        ),
        _text(item.get("display_name")).lower(),
        item["identity_key"],
    ))
    for index, compound in enumerate(compounds):
        compound["index"] = index
    status_counts = Counter(item["conclusion"]["key"] for item in compounds)
    connectivity_high_count = int(status_counts.get("high", 0))
    exact_full_identity_high_count = sum(
        item["conclusion"].get("key") == "high"
        and bool(item["conclusion"].get("high_confidence_exact_structure_detected"))
        for item in compounds
    )
    isomer_group_count = int(status_counts.get("isomer", 0))
    isomer_group_high_count = sum(
        bool(item["conclusion"].get("high_confidence_group_detected"))
        and item["conclusion"].get("key") == "isomer"
        for item in compounds
    )
    background_high_count = sum(
        item["conclusion"].get("key") == "background"
        and bool(item["conclusion"].get("high_confidence_group_detected"))
        for item in compounds
    )
    spectral_high_scope_count = (
        connectivity_high_count + isomer_group_high_count + background_high_count
    )
    unverified_blank_warning_count = sum(
        bool(item["conclusion"].get("unverified_blank_association"))
        for item in compounds
    )
    unverified_blank_warning_high_count = sum(
        bool(item["conclusion"].get("unverified_blank_association"))
        and bool(item["conclusion"].get("high_confidence_group_detected"))
        and item["conclusion"].get("key") in {"high", "isomer"}
        for item in compounds
    )
    main_event_counts = Counter(
        _text(item.get("feature_id"))
        for item in compounds
        if bool(item["conclusion"].get("counted_as_sample_compound"))
        if _text(item.get("feature_id"))
    )
    main_event_violations = sorted(
        event_id for event_id, count in main_event_counts.items() if count > 1
    )
    atomic_identity_violations = []
    for item in compounds:
        representative_match = (
            (item.get("representative") or {}).get("match") or {}
        )
        lead_candidate = item.get("lead_candidate") or {}
        if _text(representative_match.get("primary_identity_source")) != (
            "v15_flagship"
        ):
            continue
        mismatched_fields = [
            field
            for field in ("compound_name", "formula", "inchikey", "smiles")
            if _text(representative_match.get(field))
            != _text(lead_candidate.get(field))
        ]
        if mismatched_fields:
            atomic_identity_violations.append(
                {
                    "compound_id": item.get("compound_id"),
                    "identity_key": item.get("identity_key"),
                    "primary_entity_id": representative_match.get(
                        "primary_entity_id"
                    ),
                    "fields": mismatched_fields,
                    "representative": {
                        field: representative_match.get(field)
                        for field in mismatched_fields
                    },
                    "lead_candidate": {
                        field: lead_candidate.get(field)
                        for field in mismatched_fields
                    },
                }
            )
    if atomic_identity_violations:
        raise RuntimeError(
            "Compound projection spliced fields from different primary candidates: "
            + json.dumps(
                atomic_identity_violations[:5],
                ensure_ascii=False,
                sort_keys=True,
            )
        )
    primary_stat = (
        primary_path.stat()
        if primary_path is not None and primary_path.is_file()
        else None
    )
    visible_lcms_feature_ids = {
        str(item.get("feature_id"))
        for item in compounds
        if str(item.get("feature_id") or "").startswith("lcms:")
    }
    visible_fallback_result_units = {
        str(item.get("feature_id"))
        for item in compounds
        if not str(item.get("feature_id") or "").startswith("lcms:")
    }
    missing_inventory_feature_ids = sorted(
        inventory_feature_ids - visible_lcms_feature_ids
    )
    unexpected_lcms_feature_ids = (
        sorted(visible_lcms_feature_ids - inventory_feature_ids)
        if feature_inventory is not None
        else []
    )
    inventory_tier_counts = Counter(
        str((item.get("display") or {}).get("tier") or "unknown")
        for item in compounds
        if str(item.get("feature_id") or "") in inventory_feature_ids
    )
    independent_partition_satisfied = (
        not missing_inventory_feature_ids
        and not unexpected_lcms_feature_ids
        and sum(inventory_tier_counts.values()) == len(inventory_feature_ids)
    ) if feature_inventory is not None else None
    inventory_stat = (
        feature_inventory_source.stat()
        if feature_inventory_source is not None
        and feature_inventory_source.is_file()
        else None
    )
    inventory_source_sha256 = (
        hashlib.sha256(feature_inventory_source.read_bytes()).hexdigest()
        if inventory_stat is not None
        else None
    )
    feature_inventory_digest = (
        hashlib.sha256(_canonical(inventory_rows)).hexdigest()
        if feature_inventory is not None
        else None
    )

    metadata = {
        "schema_version": COMPOUND_CONCLUSION_SCHEMA_VERSION,
        "algorithm_id": COMPOUND_CONCLUSION_ALGORITHM_ID,
        "task_id": task_id,
        "source_row_count": source_rows,
        "source_matched_row_count": matched_rows,
        "entity_high_evidence_algorithm_id": ENTITY_HIGH_ALGORITHM_ID,
        "single_energy_high_supported_entity_count": sum(bool(v) for v in high_evidence_by_entity.values()),
        "source_semantic_sha256": source_digest.hexdigest(),
        "report_role_catalog": report_role_catalog_summary(),
        "plant_candidate_count": sum(
            item["conclusion"]["counted_as_plant_candidate"] for item in compounds
        ),
        "external_review_count": sum(
            item["plant_report"]["excluded_from_plant_candidate_table"]
            for item in compounds
        ),
        "plant_candidate_selection": {
            "rule": "sample_supported_with_no_reviewed_external_candidate_in_competitive_group",
            "unclassified_origin_is_endogenous_evidence": False,
            "changes_spectral_confidence": False,
            "excluded_records_retained_in_all_and_external_review": True,
        },
        "unique_report_rows": len(compounds),
        "main_result_count": sum(
            item["conclusion"]["counted_as_sample_compound"] for item in compounds
        ),
        "confirmed_present_count": int(status_counts.get("confirmed", 0)),
        # A High result must always carry its reportable identity scope.  Keep
        # 2-D connectivity, full structure identity and unresolved isomer-group
        # counts separate so a zero in one layer cannot erase High evidence in
        # another layer.
        "high_confidence_detected_count": connectivity_high_count,
        "high_confidence_2d_connectivity_count": connectivity_high_count,
        "high_confidence_exact_structure_count": exact_full_identity_high_count,
        "high_confidence_isomer_group_count": isomer_group_high_count,
        "high_confidence_scope_group_count": (
            connectivity_high_count + isomer_group_high_count
        ),
        "background_high_confidence_scope_group_count": background_high_count,
        "spectral_high_confidence_scope_group_count": spectral_high_scope_count,
        "high_confidence_scope_group_definition": (
            "sample-supported High evidence groups at their honest reportable "
            "scope: single 2D connectivity or unresolved same-formula isomer group"
        ),
        "isomer_group_detected_count": isomer_group_count,
        "isomer_group_non_high_count": (
            isomer_group_count - isomer_group_high_count
        ),
        "candidate_count": int(status_counts.get("candidate", 0)),
        "unresolved_precursor_event_count": int(
            status_counts.get("unresolved_event", 0)
        ),
        "insufficient_count": int(status_counts.get("insufficient", 0)),
        "unannotated_feature_count": int(
            status_counts.get("unannotated", 0)
        ),
        "background_related_count": int(status_counts.get("background", 0)),
        "unverified_blank_warning_count": unverified_blank_warning_count,
        "unverified_blank_warning_high_count": (
            unverified_blank_warning_high_count
        ),
        "legacy_high_supported_identity_count": len(legacy_high_identities),
        "status_counts": dict(sorted(status_counts.items())),
        "precursor_event_aggregation": event_aggregation,
        "precursor_event_main_result_invariant": {
            "rule": "at_most_one_main_result_per_ion_feature",
            "satisfied": not main_event_violations,
            "violation_count": len(main_event_violations),
            "violating_event_ids": main_event_violations,
        },
        "feature_result_projection": {
            "record_unit": "ion_feature",
            "primary_key": ["run_id", "feature_id"],
            "candidate_names_used_for_grouping": False,
            "independent_inventory_available": feature_inventory is not None,
            "inventory_basis": (
                "openms_continuous_ms1_feature_table"
                if feature_inventory is not None
                else "unavailable_result_rows_are_not_an_independent_inventory"
            ),
            "n_detected_ion_features": (
                len(inventory_feature_ids) if feature_inventory is not None else None
            ),
            "n_visible_result_units": len(compounds),
            "n_visible_lcms_features": len(visible_lcms_feature_ids),
            "n_fallback_result_units": len(visible_fallback_result_units),
            "n_library_released_features": (
                int(inventory_tier_counts.get("library_released", 0))
                if feature_inventory is not None
                else int(status_counts.get("confirmed", 0))
            ),
            "n_library_review_features": (
                int(inventory_tier_counts.get("library_review", 0))
                if feature_inventory is not None
                else (
                    len(compounds)
                    - int(status_counts.get("confirmed", 0))
                    - int(status_counts.get("unannotated", 0))
                )
            ),
            "n_unannotated_features": (
                int(inventory_tier_counts.get("unannotated", 0))
                if feature_inventory is not None
                else int(status_counts.get("unannotated", 0))
            ),
            "inventory_tier_counts": dict(sorted(inventory_tier_counts.items())),
            "missing_inventory_feature_count": len(
                missing_inventory_feature_ids
            ),
            "missing_inventory_feature_ids": missing_inventory_feature_ids[:100],
            "unexpected_lcms_feature_count": len(
                unexpected_lcms_feature_ids
            ),
            "unexpected_lcms_feature_ids": unexpected_lcms_feature_ids[:100],
            "result_lcms_feature_count_before_inventory_projection": len(
                result_lcms_feature_ids
            ),
            "partition_satisfied": independent_partition_satisfied,
            "partition_requires_independent_inventory": True,
        },
        "primary_identity_projection": {
            "status": "applied" if primary_index else "legacy_fallback",
            "source": "v15_flagship" if primary_index else None,
            "algorithm_id": primary_metadata.get("algorithm_id"),
            "policy_sha256": primary_metadata.get("policy_sha256"),
            "entity_semantic_sha256": primary_metadata.get(
                "entity_semantic_sha256"
            ),
            "mapped_result_index_count": len(primary_index),
            "candidate_order_source": (
                "v15_entity_rank" if primary_index else "legacy_scan_rank_score"
            ),
            "candidate_names_used_for_resolution": False,
            "ground_truth_used_for_resolution": False,
        },
        "atomic_candidate_identity_invariant": {
            "rule": "name_formula_inchikey_smiles_from_one_v15_candidate",
            "satisfied": True,
            "violation_count": 0,
        },
        "identity_rank_influence": {
            "expert": False,
            "classification": False,
            "biosource": False,
            "fragmentation_explanation": False,
        },
        "deduplication_unit": (
            "2d_connectivity_or_unresolved_formula_isomer_group_or_"
            "unresolved_precursor_event"
        ),
        "confirmed_requires_formal_claim_release": True,
        "source_artifact": str(source_artifact.resolve()) if source_artifact else None,
        "source_artifact_size_bytes": source_artifact.stat().st_size if source_artifact and source_artifact.is_file() else None,
        "source_artifact_mtime_ns": source_artifact.stat().st_mtime_ns if source_artifact and source_artifact.is_file() else None,
        "feature_inventory_semantic_sha256": feature_inventory_digest,
        "feature_inventory_source": (
            str(feature_inventory_source.resolve())
            if feature_inventory_source is not None
            else None
        ),
        "feature_inventory_source_size_bytes": (
            inventory_stat.st_size if inventory_stat is not None else None
        ),
        "feature_inventory_source_mtime_ns": (
            inventory_stat.st_mtime_ns if inventory_stat is not None else None
        ),
        "feature_inventory_source_sha256": inventory_source_sha256,
        "primary_identity_artifact": (
            str(primary_path.resolve()) if primary_path is not None else None
        ),
        "primary_identity_artifact_size_bytes": (
            primary_stat.st_size if primary_stat is not None else None
        ),
        "primary_identity_artifact_mtime_ns": (
            primary_stat.st_mtime_ns if primary_stat is not None else None
        ),
    }
    return {
        "metadata": metadata,
        "compounds": compounds,
    }


def write_compound_conclusion_artifact(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
    output_path: Path,
    *,
    claim_decisions_path: Path | None = None,
    source_artifact: Path | None = None,
    primary_identity_artifact: Mapping[str, Any] | Path | None = None,
    feature_inventory: Iterable[Mapping[str, Any]] | None = None,
    feature_inventory_source: Path | None = None,
) -> dict[str, Any]:
    policy = load_v15_flagship_policy()
    implementation_digest = active_implementation_sha256(policy)
    assert_derived_destination(output_path, task_id, COMPOUND_CONCLUSION_ARTIFACT_SUFFIX,
        policy=policy, actual_implementation_sha256=implementation_digest)
    primary_binding = verified_primary_binding(primary_identity_artifact, policy,
        actual_implementation_sha256=implementation_digest)
    payload = build_compound_conclusions(
        task_id,
        rows,
        claim_decisions_path=claim_decisions_path,
        source_artifact=source_artifact,
        primary_identity_artifact=primary_identity_artifact,
        feature_inventory=feature_inventory,
        feature_inventory_source=feature_inventory_source,
    )
    if primary_binding != verified_primary_binding(primary_identity_artifact, policy,
            actual_implementation_sha256=active_implementation_sha256(policy)):
        raise ValueError("Primary artifact changed during compound materialization")
    payload["metadata"].update(primary_binding)
    payload["metadata"].update(binding_metadata(policy,
        actual_implementation_sha256=active_implementation_sha256(policy)))
    payload["metadata"]["generated_at"] = datetime.now(timezone.utc).isoformat()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(
        f".{output_path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    )
    with temporary.open("w", encoding="utf-8") as handle:
        write_json_encoded(payload, handle, ensure_ascii=False, separators=(",", ":"), default=str)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output_path)
    return {
        **payload["metadata"],
        "status": "available",
        "artifact": str(output_path.resolve()),
        "size_bytes": output_path.stat().st_size,
    }


def filter_compound_conclusions(
    payload: Mapping[str, Any], *, report_scope: str = "all",
) -> dict[str, Any]:
    """Select report rows without deleting evidence or changing identities.

    A mixed unresolved group stays intact in external review; biological-role
    filtering cannot choose a new winner or resolve its spectral ambiguity.
    """
    if report_scope not in {"all", "plant_candidates", "external_review"}:
        raise ValueError("Unknown compound report scope")
    compounds = list(payload.get("compounds") or [])
    if report_scope == "plant_candidates":
        compounds = [item for item in compounds if bool(
            (item.get("conclusion") or {}).get("counted_as_plant_candidate")
        )]
    elif report_scope == "external_review":
        compounds = [item for item in compounds if bool(
            (item.get("plant_report") or {}).get("excluded_from_plant_candidate_table")
        )]
    return {
        **payload,
        "metadata": {
            **dict(payload.get("metadata") or {}),
            "report_scope": report_scope,
            "selected_report_rows": len(compounds),
            "selection_changes_spectral_confidence": False,
        },
        "compounds": compounds,
    }


def read_compound_conclusion_artifact(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 2000,
    conclusion: str | None = None,
    report_scope: str = "all",
) -> dict[str, Any]:
    if offset < 0 or limit < 1:
        raise ValueError("offset must be >=0 and limit must be >=1")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload = filter_compound_conclusions(payload, report_scope=report_scope)
    compounds = list(payload.get("compounds") or [])
    wanted = _text(conclusion).lower()
    if wanted:
        compounds = [
            item for item in compounds
            if _text((item.get("conclusion") or {}).get("key")).lower() == wanted
        ]
    return {
        "metadata": payload.get("metadata") or {},
        "compounds": compounds[offset:offset + limit],
        "compounds_total": len(compounds),
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < len(compounds),
    }


def artifact_is_current(
    path: Path,
    source_artifact: Path,
    primary_identity_artifact: Path | None = None,
    feature_inventory_artifact: Path | None = None,
) -> bool:
    if not path.is_file() or not source_artifact.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        metadata = payload.get("metadata") or {}
        policy = load_v15_flagship_policy()
        if not metadata_binding_is_current(metadata, policy,
                actual_implementation_sha256=active_implementation_sha256(policy)):
            return False
        stat = source_artifact.stat()
        source_current = (
            metadata.get("algorithm_id") == COMPOUND_CONCLUSION_ALGORITHM_ID
            and metadata.get("report_role_catalog") == report_role_catalog_summary()
            and int(metadata.get("source_artifact_size_bytes") or -1) == stat.st_size
            and int(metadata.get("source_artifact_mtime_ns") or -1) == stat.st_mtime_ns
        )
        if not source_current:
            return False
        if feature_inventory_artifact is None:
            if metadata.get("feature_inventory_source") not in {None, ""}:
                return False
        else:
            if not feature_inventory_artifact.is_file():
                return False
            feature_stat = feature_inventory_artifact.stat()
            if (
                metadata.get("feature_inventory_source")
                != str(feature_inventory_artifact.resolve())
                or int(
                    metadata.get("feature_inventory_source_size_bytes") or -1
                ) != feature_stat.st_size
                or int(
                    metadata.get("feature_inventory_source_mtime_ns") or -1
                ) != feature_stat.st_mtime_ns
                or metadata.get("feature_inventory_source_sha256")
                != hashlib.sha256(
                    feature_inventory_artifact.read_bytes()
                ).hexdigest()
            ):
                return False
        if primary_identity_artifact is None:
            return metadata.get("primary_identity_artifact") in {None, ""}
        if not primary_identity_artifact.is_file():
            return False
        if "artifact_namespace" in policy:
            current_primary = verified_primary_binding(primary_identity_artifact, policy,
                actual_implementation_sha256=active_implementation_sha256(policy))
            if any(metadata.get(key) != value for key, value in current_primary.items()):
                return False
        primary_stat = primary_identity_artifact.stat()
        return (
            int(metadata.get("primary_identity_artifact_size_bytes") or -1)
            == primary_stat.st_size
            and int(metadata.get("primary_identity_artifact_mtime_ns") or -1)
            == primary_stat.st_mtime_ns
        )
    except (OSError, ValueError, TypeError):
        return False


__all__ = [
    "COMPOUND_CONCLUSION_ALGORITHM_ID",
    "COMPOUND_CONCLUSION_ARTIFACT_SUFFIX",
    "artifact_is_current",
    "build_compound_conclusions",
    "filter_compound_conclusions",
    "compound_conclusion_artifact_path",
    "read_compound_conclusion_artifact",
    "ranked_candidate_members",
    "write_compound_conclusion_artifact",
]
