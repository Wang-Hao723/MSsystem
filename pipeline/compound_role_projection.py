"""Report-only chemical-role curation, outside the immutable search cache.

Role review changes display selection, never spectral evidence or identity.
Unknown biological origin is not positive evidence of an endogenous compound.
"""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from matching.compound_annotations import (
    annotate_compound_dict as _base_annotate_compound_dict,
    annotate_match_compound_roles as _base_annotate_match_compound_roles,
    compound_annotation_catalog_summary as _base_catalog_summary,
)


_CATALOG_PATH = Path(__file__).resolve().parents[1] / "config/report_compound_roles.json"
_ROLE_FIELDS = (
    "compound_role", "potential_contaminant", "contaminant_category",
    "contaminant_note", "compound_role_note", "exclude_from_natural_product_count",
    "natural_product_candidate", "screening_domain", "pesticide_use",
    "chemical_family", "regulatory_status", "identity_source", "screening_catalog_id",
    "screening_catalog_released_at", "screening_identity_scope", "screening_identity_note",
    "role_catalog_id", "role_catalog_version", "role_annotation_scope", "role_annotation_source",
)
_CHILD_FIELDS = ("top_candidates", "isomer_set", "same_formula_alternatives")
_EXTERNAL_ROLES = frozenset({
    "external_contaminant", "curated_contaminant", "pesticide_screening_target",
    "environmental_contaminant", "industrial_contaminant", "transformation_product",
    "drug", "pharmaceutical", "animal_specific_metabolite",
})


@lru_cache(maxsize=4096)
def _structure_key(smiles: str) -> str:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(mol) if mol is not None else ""


@lru_cache(maxsize=4)
def _parse_catalog(raw: bytes) -> dict[str, Any]:
    payload = json.loads(raw)
    if payload.get("scope") != "candidate_role_for_report_selection_only":
        raise ValueError("Unsupported report compound-role scope")
    annotations = payload.get("annotations")
    if not isinstance(annotations, dict):
        raise ValueError("Report compound-role annotations must be a mapping")
    for key, annotation in annotations.items():
        if _structure_key(str(annotation.get("smiles") or "")) != key:
            raise ValueError("Report role catalog structure does not match its full InChIKey")
        if annotation.get("compound_role") not in _EXTERNAL_ROLES:
            raise ValueError("Report supplement currently permits external-role exclusions only")
        if not annotation.get("identity_source") or not annotation.get("role_annotation_source"):
            raise ValueError("Report role requires identity and role source citations")
    return payload


def report_role_catalog_summary() -> dict[str, Any]:
    raw = _CATALOG_PATH.read_bytes()
    payload = _parse_catalog(raw)
    return {
        "catalog_id": payload["catalog_id"],
        "catalog_version": payload["catalog_version"],
        "sha256": hashlib.sha256(raw).hexdigest(),
        "annotation_count": len(payload["annotations"]),
        "scope": payload["scope"],
        "changes_spectral_evidence": False,
    }


def compound_annotation_catalog_summary() -> dict[str, Any]:
    return {**_base_catalog_summary(), "report_supplement": report_role_catalog_summary()}


def _bound_report_prior_veto(item: dict[str, Any]) -> bool | None:
    """Read only a fully bound report projection's pre-existing veto.

    Older or partial report markers have no provenance for the original flag
    and cannot bypass the existing High veto. Changing the identity invalidates
    this marker rather than carrying an exemption to another candidate.
    """
    marker = item.get("report_role_projection") or {}
    if not isinstance(marker, dict):
        return None
    catalog = _parse_catalog(_CATALOG_PATH.read_bytes())
    key = str(item.get("inchikey") or "").strip().upper()
    smiles = str(item.get("smiles") or "").strip()
    if not key and smiles:
        key = _structure_key(smiles)
    if (
        marker.get("status") == "applied"
        and marker.get("projection_scope") == "report_selection_only"
        and marker.get("catalog_id") == catalog["catalog_id"]
        and marker.get("catalog_version") == catalog["catalog_version"]
        and marker.get("matched_inchikey") == key
        and key in catalog["annotations"]
        and (not smiles or _structure_key(smiles) == key)
        and isinstance(marker.get("pre_report_contaminant_veto"), bool)
    ):
        return marker["pre_report_contaminant_veto"]
    return None


def _prior_spectral_contaminant_veto(item: dict[str, Any]) -> bool:
    recorded = _bound_report_prior_veto(item)
    return bool(item.get("potential_contaminant")) if recorded is None else recorded


def spectral_contaminant_veto(item: dict[str, Any]) -> bool:
    """Keep an original veto; a report-only role cannot create a new one."""
    return _prior_spectral_contaminant_veto(item)


def _supplement(item: dict[str, Any], *, pre_report_contaminant_veto: bool) -> None:
    catalog = _parse_catalog(_CATALOG_PATH.read_bytes())
    key = str(item.get("inchikey") or "").strip().upper()
    reported_key = key
    smiles = str(item.get("smiles") or "").strip()
    # A structure is sufficient when an old row lacks a full key. Names alone
    # cannot bind the rule; nor can a matching name override a conflicting key.
    if not key and smiles:
        key = _structure_key(smiles)
    annotation = catalog["annotations"].get(key)
    item.pop("report_role_projection", None)
    if not annotation:
        if pre_report_contaminant_veto:
            item["potential_contaminant"] = True
        return
    if smiles and _structure_key(smiles) != key:
        if pre_report_contaminant_veto:
            item["potential_contaminant"] = True
        item["report_role_projection"] = {
            "status": "identity_conflict_not_applied",
            "catalog_id": catalog["catalog_id"],
            "matched_inchikey": key,
        }
        return
    for field in _ROLE_FIELDS:
        item.pop(field, None)
    note = annotation["note_cn"]
    item.update({
        "compound_role": annotation["compound_role"],
        "potential_contaminant": True,
        "contaminant_category": annotation["contaminant_category"],
        "contaminant_note": note,
        "compound_role_note": note,
        "exclude_from_natural_product_count": True,
        "natural_product_candidate": False,
        "identity_source": annotation["identity_source"],
        "role_catalog_id": catalog["catalog_id"],
        "role_catalog_version": catalog["catalog_version"],
        "role_annotation_scope": annotation["role_annotation_scope"],
        "role_annotation_source": annotation["role_annotation_source"],
        "report_role_projection": {
            "status": "applied",
            "catalog_id": catalog["catalog_id"],
            "catalog_version": catalog["catalog_version"],
            "matched_inchikey": key,
            "projection_scope": "report_selection_only",
            "pre_report_contaminant_veto": bool(pre_report_contaminant_veto),
            "identity_bound_by": (
                "smiles_and_full_inchikey" if smiles and reported_key
                else "verified_smiles" if smiles else "full_inchikey"
            ),
            "detection_confirmed_by_role": False,
            "spectral_evidence_changed": False,
        },
    })


def _walk(match: dict[str, Any]):
    stack = [match]
    visited: set[int] = set()
    while stack:
        item = stack.pop()
        if id(item) in visited:
            continue
        visited.add(id(item))
        yield item
        for field in _CHILD_FIELDS:
            stack.extend(child for child in item.get(field) or [] if isinstance(child, dict))


def apply_report_compound_roles(match: dict[str, Any] | None):
    """Refresh role fields only, after restoring an atomic ranked identity.

    The frozen annotator also normalizes names and classification. Its role
    projection is copied from a temporary dictionary, keeping those unrelated
    fields and all rank/confidence fields on the original candidate untouched.
    """
    if not isinstance(match, dict):
        return match
    for item in _walk(match):
        original_veto = _prior_spectral_contaminant_veto(item)
        temporary = dict(item)
        _base_annotate_compound_dict(temporary)
        for field in _ROLE_FIELDS:
            item.pop(field, None)
            if field in temporary:
                item[field] = temporary[field]
        _supplement(item, pre_report_contaminant_veto=(
            original_veto or bool(temporary.get("potential_contaminant"))
        ))
    return match


def annotate_match_compound_roles(match: dict[str, Any] | None):
    original_vetoes = {
        id(item): _prior_spectral_contaminant_veto(item)
        for item in _walk(match)
    } if isinstance(match, dict) else {}
    # The frozen annotator records an old non-curated role as source metadata.
    # Our previous report-only role is not an original source annotation.
    source_roles = {
        id(item): ("source_compound_role" in item, item.get("source_compound_role"))
        for item in _walk(match) if _bound_report_prior_veto(item) is not None
    } if isinstance(match, dict) else {}
    _base_annotate_match_compound_roles(match)
    if isinstance(match, dict):
        for item in _walk(match):
            if id(item) in source_roles:
                existed, original = source_roles[id(item)]
                if existed:
                    item["source_compound_role"] = original
                else:
                    item.pop("source_compound_role", None)
            _supplement(item, pre_report_contaminant_veto=(
                original_vetoes[id(item)] or bool(item.get("potential_contaminant"))
            ))
    return match


def annotate_result_compound_roles(results: list[dict[str, Any]]):
    for row in results or []:
        if isinstance(row, dict):
            annotate_match_compound_roles(row.get("match"))
    return results


def plant_report_selection(match: dict[str, Any] | None) -> dict[str, Any]:
    item = match if isinstance(match, dict) else {}
    external = bool(item.get("potential_contaminant") is True or item.get("compound_role") in _EXTERNAL_ROLES)
    conflict = (item.get("report_role_projection") or {}).get("status") == "identity_conflict_not_applied"
    return {
        "excluded_from_plant_candidate_table": bool(external or conflict),
        "review_category": "identity_conflict" if conflict else ("external_candidate" if external else None),
        "reason": "candidate_structure_key_conflict" if conflict else ("reviewed_external_role" if external else "no_reviewed_external_role"),
        # A natural-product role does not establish occurrence in this plant.
        "plant_origin_verified": False,
        "curated_natural_role": bool(item.get("compound_role") == "curated_natural"),
        "unclassified_origin_is_endogenous_evidence": False,
        "detection_confirmed_by_role": False,
        "role": item.get("compound_role") or "unclassified",
        "role_annotation_source": item.get("role_annotation_source") or None,
    }
