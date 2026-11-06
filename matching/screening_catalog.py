"""Versioned, identity-first catalog for contaminants and pesticide targets.

This module is deliberately separate from the spectral library.  A catalog
entry says that an exact chemical identity is relevant to a screening domain;
it does not say that a sample contains the compound, that an MRL is exceeded,
or that a library spectrum is fit for confirmation.  Runtime annotation is
activated only by a full InChIKey present in the versioned catalog.
"""

from __future__ import annotations

import json
import re
import sqlite3
from functools import lru_cache
from pathlib import Path
from typing import Any

from .library_quality import SEARCHABLE_SOURCE_TIERS


_CATALOG_PATH = Path(__file__).with_name("screening_catalog.json")
_FULL_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_SCREENING_ROLES = {
    "pesticide_screening_target",
    "environmental_contaminant",
    "industrial_contaminant",
    "transformation_product",
}


def _normalized_key(value: Any) -> str:
    key = str(value or "").strip().upper().replace("INCHIKEY=", "")
    return key if _FULL_INCHIKEY_RE.fullmatch(key) else ""


def validate_screening_catalog(
    payload: dict[str, Any], *, validate_structures: bool = False
) -> list[str]:
    """Return validation errors; an empty list means the catalog is valid."""
    errors: list[str] = []
    if not isinstance(payload, dict):
        return ["catalog root must be an object"]
    if payload.get("schema_version") != 1:
        errors.append("schema_version must equal 1")
    catalog_id = str(payload.get("catalog_id") or "").strip()
    if not catalog_id:
        errors.append("catalog_id is required")
    entries = payload.get("entries")
    if not isinstance(entries, list) or not entries:
        return [*errors, "entries must be a non-empty array"]

    seen_keys: dict[str, str] = {}
    for index, entry in enumerate(entries):
        path = f"entries[{index}]"
        if not isinstance(entry, dict):
            errors.append(f"{path} must be an object")
            continue
        name = str(entry.get("canonical_name") or "").strip()
        if not name:
            errors.append(f"{path}.canonical_name is required")
            name = path
        role = str(entry.get("compound_role") or "").strip()
        if role not in _SCREENING_ROLES:
            errors.append(f"{path}.compound_role is unsupported: {role!r}")
        if str(entry.get("regulatory_status") or "") != "not_encoded":
            errors.append(
                f"{path}.regulatory_status must be 'not_encoded'; "
                "matrix- and date-specific legal limits belong in a separate layer"
            )
        primary_key = _normalized_key(entry.get("full_inchikey"))
        if not primary_key:
            errors.append(f"{path}.full_inchikey is invalid")
        keys = [primary_key] if primary_key else []
        aliases = entry.get("alternate_full_inchikeys") or []
        if not isinstance(aliases, list):
            errors.append(f"{path}.alternate_full_inchikeys must be an array")
            aliases = []
        for alias in aliases:
            normalized = _normalized_key(alias)
            if not normalized:
                errors.append(f"{path}.alternate_full_inchikeys contains invalid key {alias!r}")
            else:
                keys.append(normalized)
        for key in keys:
            owner = seen_keys.get(key)
            if owner and owner != name:
                errors.append(f"full InChIKey {key} is duplicated by {owner!r} and {name!r}")
            seen_keys[key] = name
        if not str(entry.get("identity_source") or "").startswith("https://"):
            errors.append(f"{path}.identity_source must be an https URL")

        if validate_structures and primary_key and entry.get("smiles"):
            try:
                from rdkit import Chem

                mol = Chem.MolFromSmiles(str(entry["smiles"]))
                derived = Chem.MolToInchiKey(mol) if mol is not None else ""
                if derived != primary_key:
                    errors.append(
                        f"{path}.smiles derives {derived or 'no key'}, expected {primary_key}"
                    )
            except Exception as exc:  # pragma: no cover - environment-specific
                errors.append(f"{path}.smiles validation failed: {type(exc).__name__}: {exc}")
    return errors


@lru_cache(maxsize=1)
def load_screening_catalog() -> dict[str, Any]:
    """Load and minimally validate the bundled catalog."""
    try:
        payload = json.loads(_CATALOG_PATH.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return {"schema_version": 1, "catalog_id": "unavailable", "entries": []}
    errors = validate_screening_catalog(payload)
    if errors:
        return {
            "schema_version": 1,
            "catalog_id": "invalid",
            "entries": [],
            "validation_errors": errors,
        }
    return payload


@lru_cache(maxsize=1)
def _catalog_by_full_inchikey() -> dict[str, dict[str, Any]]:
    payload = load_screening_catalog()
    catalog_id = str(payload.get("catalog_id") or "")
    released_at = str(payload.get("released_at") or "")
    index: dict[str, dict[str, Any]] = {}
    for entry in payload.get("entries") or []:
        base = {
            "canonical_name": entry.get("canonical_name", ""),
            "name_cn": entry.get("name_cn", ""),
            "compound_role": entry.get("compound_role", ""),
            "contaminant_category": entry.get("contaminant_category", ""),
            "exclude_from_natural_product_count": True,
            "screening_domain": entry.get("screening_domain", ""),
            "pesticide_use": entry.get("pesticide_use", ""),
            "chemical_family": entry.get("chemical_family", ""),
            "regulatory_status": entry.get("regulatory_status", "not_encoded"),
            "identity_source": entry.get("identity_source", ""),
            "screening_catalog_id": catalog_id,
            "screening_catalog_released_at": released_at,
            "note_cn": entry.get("note_cn", ""),
        }
        primary = _normalized_key(entry.get("full_inchikey"))
        if primary:
            index[primary] = {**base, "screening_identity_scope": "full_inchikey"}
        for alias in entry.get("alternate_full_inchikeys") or []:
            key = _normalized_key(alias)
            if key:
                index[key] = {
                    **base,
                    "screening_identity_scope": "connectivity_alias_full_key",
                    "screening_identity_note": entry.get(
                        "alternate_key_note",
                        "Role annotation only; exact stereochemical identity remains unresolved.",
                    ),
                }
    return index


def screening_annotation_for(inchikey: Any) -> dict[str, Any]:
    """Return an exact-key screening-role annotation, never a name match."""
    key = _normalized_key(inchikey)
    annotation = _catalog_by_full_inchikey().get(key)
    return dict(annotation) if annotation else {}


def catalog_summary() -> dict[str, Any]:
    payload = load_screening_catalog()
    entries = payload.get("entries") or []
    domains: dict[str, int] = {}
    for entry in entries:
        domain = str(entry.get("screening_domain") or "unclassified")
        domains[domain] = domains.get(domain, 0) + 1
    return {
        "catalog_id": payload.get("catalog_id"),
        "released_at": payload.get("released_at"),
        "entries": len(entries),
        "indexed_full_inchikeys": len(_catalog_by_full_inchikey()),
        "domains": domains,
        "regulatory_claims_encoded": False,
    }


def audit_library_coverage(database_path: str | Path) -> dict[str, Any]:
    """Read-only audit of exact catalog identities and searchable spectra."""
    db_path = Path(database_path).expanduser().resolve()
    keys = sorted(_catalog_by_full_inchikey())
    if not keys:
        return {**catalog_summary(), "database_path": str(db_path), "identities": []}
    placeholders = ",".join("?" for _ in keys)
    searchable_tiers = sorted(SEARCHABLE_SOURCE_TIERS)
    tier_placeholders = ",".join("?" for _ in searchable_tiers)
    query = f"""
        SELECT upper(trim(c.inchikey)) AS inchikey,
               c.compound_name,
               c.structure_identity_verified,
               COUNT(DISTINCT s.id) AS spectrum_count,
               COUNT(DISTINCT CASE WHEN lower(coalesce(s.source_tier, ''))
                   IN ({tier_placeholders}) THEN s.id END)
                   AS searchable_spectrum_count,
               group_concat(DISTINCT s.source_tier) AS source_tiers,
               group_concat(DISTINCT s.ion_mode) AS ion_modes
        FROM Compounds c
        LEFT JOIN Mapping m ON m.compound_id = c.id
        LEFT JOIN Spectra s ON s.id = m.spectrum_id
        WHERE upper(trim(c.inchikey)) IN ({placeholders})
        GROUP BY upper(trim(c.inchikey)), c.id
        ORDER BY upper(trim(c.inchikey)), c.id
    """
    connection = sqlite3.connect(
        f"file:{db_path}?mode=ro&immutable=1", uri=True
    )
    connection.row_factory = sqlite3.Row
    try:
        rows = [
            dict(row)
            for row in connection.execute(query, [*searchable_tiers, *keys])
        ]
    finally:
        connection.close()

    by_key: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_key.setdefault(str(row["inchikey"]), []).append(row)
    identities = []
    for key in keys:
        records = by_key.get(key, [])
        annotation = _catalog_by_full_inchikey()[key]
        identities.append({
            "inchikey": key,
            "canonical_name": annotation.get("canonical_name"),
            "screening_identity_scope": annotation.get("screening_identity_scope"),
            "compound_records": len(records),
            "spectrum_count": sum(int(row.get("spectrum_count") or 0) for row in records),
            "searchable_spectrum_count": sum(
                int(row.get("searchable_spectrum_count") or 0) for row in records
            ),
            "structure_verified_records": sum(
                1 for row in records if int(row.get("structure_identity_verified") or 0) == 1
            ),
            "source_tiers": sorted({
                tier
                for row in records
                for tier in str(row.get("source_tiers") or "").split(",")
                if tier
            }),
            "ion_modes": sorted({
                mode
                for row in records
                for mode in str(row.get("ion_modes") or "").split(",")
                if mode
            }),
        })
    covered = [item for item in identities if item["searchable_spectrum_count"] > 0]
    return {
        **catalog_summary(),
        "database_path": str(db_path),
        "catalog_full_keys_with_searchable_spectra": len(covered),
        "catalog_full_key_coverage_fraction": round(len(covered) / len(keys), 6),
        "identities": identities,
    }
