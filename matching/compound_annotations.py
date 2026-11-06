"""Curated compound roles that must not be inferred from spectral similarity.

The mass-spectrum match answers *what ion resembles a library spectrum*.
It does not answer whether that compound is an endogenous/natural constituent
or an external contaminant.  This module keeps that orthogonal provenance role
as exact full-InChIKey curation and can also upgrade legacy cached results.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from .screening_catalog import screening_annotation_for


_ANNOTATIONS_PATH = Path(__file__).with_name("compound_annotations.json")
_SUPPLEMENTAL_ANNOTATIONS_PATH = Path(__file__).with_name(
    "compound_annotations_v12.json"
)

_ROLE_PROJECTION_FIELDS = (
    "compound_role",
    "potential_contaminant",
    "contaminant_category",
    "contaminant_note",
    "compound_role_note",
    "exclude_from_natural_product_count",
    "natural_product_candidate",
    "screening_domain",
    "pesticide_use",
    "chemical_family",
    "regulatory_status",
    "identity_source",
    "screening_catalog_id",
    "screening_catalog_released_at",
    "screening_identity_scope",
    "screening_identity_note",
    "role_catalog_id",
    "role_catalog_version",
    "role_annotation_scope",
    "role_annotation_source",
)


@lru_cache(maxsize=1)
def _runtime_rule_classifier():
    from mssystem.classification import RuleBasedClassifier

    return RuleBasedClassifier()


def project_reliable_structure_classification(
    item: dict[str, Any] | None,
) -> dict[str, Any] | None:
    """Project a reliable runtime structure class without rewriting raw metadata.

    Historical library rows may retain ``classification=Unknown`` even though
    the expert evaluator recorded a deterministic SMILES rule.  Older cached
    rows can also predate a newly added rule.  For a validated structure only,
    rerun the deterministic SMARTS layer when no usable class exists.  Reports
    consume ``display_class``; the raw library value remains untouched.
    """
    if not isinstance(item, dict):
        return item
    details = item.get("expert_details")
    provenance = (
        details.get("classification_provenance")
        if isinstance(details, dict) else None
    )
    effective = ""
    reason = ""
    confidence = None
    if isinstance(provenance, dict) and provenance.get("reliable") is True:
        effective = str(
            provenance.get("effective_classification")
            or provenance.get("analysis_class")
            or ""
        ).strip()
        reason = str(provenance.get("reliability_reason") or "structure_rule")
        confidence = provenance.get("input_confidence")

    if not effective or effective.lower() in {"unknown", "other", "none"}:
        current_labels = (
            str(item.get("display_class") or "").strip().lower(),
            str(item.get("classification") or "").strip().lower(),
        )
        labels_missing = all(
            label in {"", "unknown", "other", "unclassified", "none", "nan"}
            for label in current_labels
        )
        current_source = str(
            item.get("display_classification_source")
            or item.get("classification_source") or ""
        ).strip().lower()
        trusted_source = current_source.startswith((
            "curated", "manual", "smarts", "structure", "rule",
            "database_curated", "runtime_structure_rule",
        ))
        no_usable_class = bool(
            labels_missing
            or not trusted_source
        )
        structure_validated = (
            str(item.get("structure_quality") or "").strip().lower()
            == "validated"
            or item.get("structure_identity_verified") is True
        )
        smiles = str(item.get("smiles") or "").strip()
        if no_usable_class and structure_validated and smiles:
            try:
                prediction = _runtime_rule_classifier().predict_detailed(
                    smiles=smiles
                )
                if prediction.get("source") == "smarts":
                    effective = str(prediction.get("class") or "").strip()
                    reason = "cached_validated_smiles_structure_rule"
                    confidence = float(prediction.get("confidence") or 1.0)
            except Exception:
                effective = ""
        if not effective or effective.lower() in {"unknown", "other", "none"}:
            return item

    item["display_class"] = effective
    item["display_classification_source"] = f"runtime_structure_rule:{reason}"
    if reason.startswith("smiles_structure_rule") or "structure_rule" in reason:
        confidence = 1.0
    item["display_classification_confidence"] = confidence
    item["classification_reliable"] = True
    item["classification_status"] = "reliable"
    item["classification_reliability_reason"] = reason
    item["classification_projection_preserved_raw"] = True
    return item


@lru_cache(maxsize=1)
def load_compound_annotation_catalog() -> dict[str, Any]:
    try:
        payload = json.loads(_ANNOTATIONS_PATH.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return {
            "catalog_id": "unavailable",
            "catalog_version": "unavailable",
            "baseline_contract": {},
            "annotations": {},
        }
    annotations = payload.get("annotations", {}) if isinstance(payload, dict) else {}
    if not isinstance(annotations, dict):
        annotations = {}
    return {
        "catalog_id": str(payload.get("catalog_id") or "legacy-unversioned"),
        "catalog_version": str(payload.get("catalog_version") or "1"),
        "baseline_contract": dict(payload.get("baseline_contract") or {}),
        "annotations": {
            str(key).strip().upper(): dict(value)
            for key, value in annotations.items()
            if str(key).strip() and isinstance(value, dict)
        },
    }


def load_compound_annotations() -> dict[str, dict[str, Any]]:
    return load_compound_annotation_catalog()["annotations"]


@lru_cache(maxsize=1)
def load_supplemental_compound_annotation_catalog() -> dict[str, Any]:
    """Load post-baseline exact-identity additions without editing V9 history.

    ``compound_annotations.json`` is a frozen, digest-tested V9 baseline.  New
    exact identities belong in a separate versioned supplement so historical
    validation counts remain reproducible while current reports can still
    recognize newly reviewed external chemicals.
    """
    try:
        payload = json.loads(
            _SUPPLEMENTAL_ANNOTATIONS_PATH.read_text(encoding="utf-8")
        )
    except (OSError, TypeError, ValueError):
        return {
            "catalog_id": "unavailable",
            "catalog_version": "unavailable",
            "annotations": {},
        }
    annotations = payload.get("annotations", {}) if isinstance(payload, dict) else {}
    if not isinstance(annotations, dict):
        annotations = {}
    return {
        "catalog_id": str(payload.get("catalog_id") or "supplemental-unversioned"),
        "catalog_version": str(payload.get("catalog_version") or "1"),
        "annotations": {
            str(key).strip().upper(): dict(value)
            for key, value in annotations.items()
            if str(key).strip() and isinstance(value, dict)
        },
    }


def load_supplemental_compound_annotations() -> dict[str, dict[str, Any]]:
    return load_supplemental_compound_annotation_catalog()["annotations"]


def compound_annotation_catalog_summary() -> dict[str, Any]:
    """Return the release contract without exposing the full role catalog."""

    catalog = load_compound_annotation_catalog()
    annotations = catalog["annotations"]
    external = [
        value
        for value in annotations.values()
        if value.get("compound_role") == "external_contaminant"
    ]
    class_level = [
        value
        for value in external
        if value.get("role_annotation_scope") == "class_level_structure_rule"
    ]
    return {
        "catalog_id": catalog["catalog_id"],
        "catalog_version": catalog["catalog_version"],
        "annotation_count": len(annotations),
        "external_contaminant_count": len(external),
        "class_level_structure_rule_count": len(class_level),
        "baseline_contract": dict(catalog["baseline_contract"]),
        "supplemental_annotation_count": len(
            load_supplemental_compound_annotations()
        ),
        "supplemental_external_contaminant_count": sum(
            value.get("compound_role") == "external_contaminant"
            for value in load_supplemental_compound_annotations().values()
        ),
    }


def annotation_for(inchikey: Any) -> dict[str, Any]:
    """Return a defensive copy for an exact full-InChIKey annotation."""
    key = str(inchikey or "").strip().upper()
    value = load_compound_annotations().get(key)
    supplemental = load_supplemental_compound_annotations().get(key)
    screening = screening_annotation_for(key)
    merged: dict[str, Any] = {}
    if screening:
        merged.update(screening)
    if value:
        merged.update(value)
    if supplemental:
        merged.update(supplemental)
        catalog = load_supplemental_compound_annotation_catalog()
        merged.setdefault("role_catalog_id", catalog["catalog_id"])
        merged.setdefault("role_catalog_version", catalog["catalog_version"])
    return merged


def _clear_role_projection(item: dict[str, Any]) -> None:
    for field in _ROLE_PROJECTION_FIELDS:
        item.pop(field, None)


def annotate_compound_dict(item: dict[str, Any] | None) -> dict[str, Any] | None:
    """Attach role/contaminant fields to a match or Top-K item in place."""
    if not isinstance(item, dict):
        return item
    project_reliable_structure_classification(item)
    annotation = annotation_for(item.get("inchikey"))
    if not annotation:
        # A spectral-library hit does not establish biological origin.  Unknown
        # drugs, pesticides and industrial chemicals must not silently become
        # "natural-product candidates" merely because they are not in the
        # small contaminant curation list.
        previous_role = str(item.get("compound_role") or "").strip()
        if previous_role and previous_role != "unclassified":
            item.setdefault("source_compound_role", previous_role)
        _clear_role_projection(item)
        item["compound_role"] = "unclassified"
        item["potential_contaminant"] = False
        item["contaminant_category"] = ""
        item["contaminant_note"] = ""
        item["compound_role_note"] = ""
        item["exclude_from_natural_product_count"] = True
        item["natural_product_candidate"] = False
        return item

    # Catalog projections are authoritative for their version.  Clear fields
    # produced by an older catalog before applying the current annotation so a
    # cached task cannot retain a stale contaminant or pesticide role.
    _clear_role_projection(item)

    canonical_name = str(annotation.get("canonical_name") or "").strip()
    if canonical_name:
        current = str(item.get("compound_name") or "").strip()
        if current and current != canonical_name and not item.get("source_compound_name"):
            item["source_compound_name"] = current
        item["compound_name"] = canonical_name
    name_cn = str(annotation.get("name_cn") or "").strip()
    if name_cn:
        item["compound_name_cn"] = name_cn

    role = str(annotation.get("compound_role") or "unclassified")
    is_contaminant = role in {
        "external_contaminant",
        "curated_contaminant",
        "pesticide_screening_target",
        "environmental_contaminant",
        "industrial_contaminant",
        "transformation_product",
    }
    is_natural = role == "curated_natural"
    item["compound_role"] = role
    item["potential_contaminant"] = is_contaminant
    item["contaminant_category"] = (
        annotation.get("contaminant_category", "external_contaminant")
        if is_contaminant else ""
    )
    item["contaminant_note"] = annotation.get("note_cn", "") if is_contaminant else ""
    item["compound_role_note"] = annotation.get("note_cn", "")
    for field in (
        "screening_domain",
        "pesticide_use",
        "chemical_family",
        "regulatory_status",
        "identity_source",
        "screening_catalog_id",
        "screening_catalog_released_at",
        "screening_identity_scope",
        "screening_identity_note",
        "role_catalog_id",
        "role_catalog_version",
        "role_annotation_scope",
        "role_annotation_source",
    ):
        if annotation.get(field) not in (None, ""):
            item[field] = annotation[field]
    item["exclude_from_natural_product_count"] = bool(
        annotation.get("exclude_from_natural_product_count", not is_natural)
    )
    item["natural_product_candidate"] = bool(
        is_natural and not item["exclude_from_natural_product_count"]
    )
    override = str(annotation.get("classification_override") or "").strip()
    if override:
        existing_class = str(item.get("classification") or "").strip()
        existing_source = str(
            item.get("classification_source") or ""
        ).strip()
        if existing_class and existing_class != override:
            item.setdefault("source_classification", existing_class)
        if existing_source and existing_source != "curated_identity_annotation":
            item.setdefault("source_classification_source", existing_source)
        item["classification"] = override
        item["display_class"] = override
        item["classification_override_source"] = "curated_compound_role"
        item["classification_source"] = "curated_identity_annotation"
        item["display_classification_source"] = (
            "curated_identity_annotation"
        )
        item["classification_confidence"] = 1.0
        item["display_classification_confidence"] = 1.0
        item["classification_reliable"] = override.casefold() not in {
            "", "other", "unknown", "unclassified"
        }
        item["classification_status"] = (
            "reliable" if item["classification_reliable"]
            else "broad_or_unclassified"
        )
        item["classification_reliability_reason"] = (
            "exact_full_inchikey_curated_identity_annotation"
        )
    display_override = str(
        annotation.get("display_classification_override") or ""
    ).strip()
    if display_override and not override:
        item["display_class"] = display_override
        item["display_classification_source"] = "curated_identity_annotation"
        item["display_classification_confidence"] = 1.0
        item["classification_projection_preserved_raw"] = True
    return item


def annotate_match_compound_roles(match: dict[str, Any] | None) -> dict[str, Any] | None:
    """Annotate best match and nested candidate lists, including legacy JSON."""
    if not isinstance(match, dict):
        return match

    visited: set[int] = set()

    def _walk(item: dict[str, Any]) -> None:
        marker = id(item)
        if marker in visited:
            return
        visited.add(marker)
        annotate_compound_dict(item)
        for field in (
            "top_candidates",
            "isomer_set",
            "same_formula_alternatives",
        ):
            for child in item.get(field) or []:
                if isinstance(child, dict):
                    _walk(child)

    _walk(match)
    return match


def annotate_result_compound_roles(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    for row in results or []:
        if isinstance(row, dict):
            annotate_match_compound_roles(row.get("match"))
    return results
