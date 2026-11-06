"""Formula-conditioned structure enumeration from the read-only library.

The input is an msFIDDLE shadow artifact.  This module never edits the shared
spectral database and never feeds candidates back into production ranking.
"""

from __future__ import annotations

import hashlib
import json
import math
import sqlite3
from pathlib import Path
from typing import Any, Mapping

from matching.library_quality import compound_identity_key, normalize_formula_text


SCHEMA_VERSION = "1.0.0"
MODULE_ID = "formula_conditioned_library_enumeration"


def _json_sha256(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _base_artifact(formula_artifact: Mapping[str, Any], db_path: Path) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "task_id": str(formula_artifact.get("task_id") or ""),
        "module_id": MODULE_ID,
        "status": "unavailable",
        "deployment_state": "shadow",
        "runtime_influence": False,
        "ranking_changed": False,
        "confidence_changed": False,
        "claim_scope_changed": False,
        "high_promotions": 0,
        "source_formula_artifact_sha256": _json_sha256(formula_artifact),
        "library": {
            "path": str(db_path),
            "available": db_path.is_file(),
            "access_mode": "read_only",
        },
        "formula_catalog": {},
        "spectra": [],
        "summary": {},
        "interpretation": (
            "These are formula-conditioned structure hypotheses from the local "
            "catalog, not spectral matches. Catalog order is not a chemical "
            "identity ranking and cannot promote High."
        ),
    }


def _predicted_formulas(formula_artifact: Mapping[str, Any]) -> list[str]:
    formulas = {
        normalize_formula_text(candidate.get("formula"))
        for spectrum in formula_artifact.get("spectra") or []
        if isinstance(spectrum, Mapping)
        for candidate in spectrum.get("formula_candidates") or []
        if isinstance(candidate, Mapping)
    }
    formulas.discard("")
    return sorted(formulas)


def _load_formula_catalog(
    db_path: Path,
    formulas: list[str],
    *,
    max_candidates_per_formula: int,
) -> dict[str, dict[str, Any]]:
    catalog = {
        formula: {
            "formula": formula,
            "unique_structure_candidates": 0,
            "returned_candidates": [],
            "truncated": False,
        }
        for formula in formulas
    }
    connection = sqlite3.connect(
        f"file:{db_path.resolve()}?mode=ro&immutable=1", uri=True
    )
    connection.row_factory = sqlite3.Row
    try:
        columns = {
            row[1] for row in connection.execute("PRAGMA table_info(Compounds)")
        }
        verified = (
            "COALESCE(c.structure_identity_verified, 0)"
            if "structure_identity_verified" in columns
            else "0"
        )
        reported = (
            "COALESCE(c.reported_inchikey, '')"
            if "reported_inchikey" in columns
            else "''"
        )
        connection.execute(
            "CREATE TEMP TABLE WantedFormula (formula TEXT PRIMARY KEY)"
        )
        connection.executemany(
            "INSERT INTO WantedFormula(formula) VALUES (?)",
            ((formula,) for formula in formulas),
        )
        cursor = connection.execute(f"""
            SELECT c.id, c.compound_name, c.molecular_formula, c.smiles,
                   c.inchikey, c.exact_mass, c.classification,
                   {verified} AS structure_identity_verified,
                   {reported} AS reported_inchikey
            FROM Compounds AS c
            INNER JOIN WantedFormula AS w
              ON w.formula = c.molecular_formula
            WHERE COALESCE(TRIM(c.smiles), '') != ''
               OR COALESCE(TRIM(c.inchikey), '') != ''
            ORDER BY c.molecular_formula,
                     structure_identity_verified DESC,
                     CASE WHEN COALESCE(TRIM(c.inchikey), '') != '' THEN 0 ELSE 1 END,
                     c.id
        """)
        seen: dict[str, set[str]] = {formula: set() for formula in formulas}
        for row in cursor:
            formula = normalize_formula_text(row["molecular_formula"])
            if formula not in catalog:
                continue
            candidate = {
                "compound_name": str(row["compound_name"] or ""),
                "formula": formula,
                "smiles": str(row["smiles"] or ""),
                "inchikey": str(row["inchikey"] or ""),
                "reported_inchikey": str(row["reported_inchikey"] or ""),
                "exact_mass": _number(row["exact_mass"]),
                "classification": str(row["classification"] or ""),
                "structure_identity_verified": bool(
                    row["structure_identity_verified"]
                ),
            }
            identity = compound_identity_key(candidate)
            if not identity or identity in seen[formula]:
                continue
            seen[formula].add(identity)
            entry = catalog[formula]
            entry["unique_structure_candidates"] += 1
            if len(entry["returned_candidates"]) < max_candidates_per_formula:
                entry["returned_candidates"].append({
                    "identity_key": identity,
                    **candidate,
                    "catalog_order_semantics": (
                        "verified_structure_then_identifier_presence_not_spectral_rank"
                    ),
                    "runtime_influence": False,
                })
        for entry in catalog.values():
            entry["truncated"] = (
                entry["unique_structure_candidates"]
                > len(entry["returned_candidates"])
            )
    finally:
        connection.close()
    return catalog


def build_formula_conditioned_library_shadow(
    formula_artifact: Mapping[str, Any],
    *,
    db_path: Path,
    max_candidates_per_formula: int = 10,
) -> dict[str, Any]:
    artifact = _base_artifact(formula_artifact, db_path)
    if formula_artifact.get("status") != "completed_shadow":
        artifact["summary"] = {"reason": "formula_shadow_not_completed"}
        return artifact
    if not db_path.is_file():
        artifact["summary"] = {"reason": "library_database_unavailable"}
        return artifact
    formulas = _predicted_formulas(formula_artifact)
    if not formulas:
        artifact["status"] = "not_applicable"
        artifact["summary"] = {"reason": "no_refined_formula_candidates"}
        return artifact
    try:
        catalog = _load_formula_catalog(
            db_path,
            formulas,
            max_candidates_per_formula=max(1, int(max_candidates_per_formula)),
        )
    except sqlite3.Error as exc:
        artifact["status"] = "failed"
        artifact["summary"] = {
            "reason": "read_only_library_query_failed",
            "error": f"{type(exc).__name__}: {exc}"[:500],
        }
        return artifact
    artifact["formula_catalog"] = catalog
    unmatched_with_structure = 0
    matched_with_structure = 0
    for spectrum in formula_artifact.get("spectra") or []:
        if not isinstance(spectrum, Mapping):
            continue
        spectrum_formulas = [
            normalize_formula_text(candidate.get("formula"))
            for candidate in spectrum.get("formula_candidates") or []
            if isinstance(candidate, Mapping)
        ]
        spectrum_formulas = [formula for formula in spectrum_formulas if formula]
        formulas_with_candidates = [
            formula
            for formula in spectrum_formulas
            if catalog.get(formula, {}).get("unique_structure_candidates", 0) > 0
        ]
        unique_candidate_count = sum(
            int(catalog[formula]["unique_structure_candidates"])
            for formula in dict.fromkeys(formulas_with_candidates)
        )
        has_structures = unique_candidate_count > 0
        matched = bool(spectrum.get("matched"))
        matched_with_structure += int(matched and has_structures)
        unmatched_with_structure += int(not matched and has_structures)
        artifact["spectra"].append({
            "row_index": spectrum.get("row_index"),
            "spectrum_id": spectrum.get("spectrum_id"),
            "matched": matched,
            "predicted_formulas": spectrum_formulas,
            "formulas_with_catalog_structures": formulas_with_candidates,
            "unique_formula_conditioned_structure_candidates": unique_candidate_count,
            "runtime_influence": False,
        })
    artifact["status"] = "completed_shadow"
    artifact["summary"] = {
        "predicted_formulas_queried": len(formulas),
        "formulas_with_structure_candidates": sum(
            entry["unique_structure_candidates"] > 0 for entry in catalog.values()
        ),
        "unique_structure_candidates_across_formulas": sum(
            int(entry["unique_structure_candidates"]) for entry in catalog.values()
        ),
        "returned_structure_candidates": sum(
            len(entry["returned_candidates"]) for entry in catalog.values()
        ),
        "spectra_evaluated": len(artifact["spectra"]),
        "matched_spectra_with_structure_hypotheses": matched_with_structure,
        "unmatched_spectra_with_structure_hypotheses": unmatched_with_structure,
        "catalog_order_is_spectral_rank": False,
        "production_ranking_changed": False,
        "high_promotions": 0,
    }
    return artifact


__all__ = ["build_formula_conditioned_library_shadow"]
