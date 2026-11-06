#!/usr/bin/env python3
"""Repair eight verified target rows whose formula fields contradict identity.

The repair is deliberately narrow and idempotent.  Every row must have the
expected full InChIKey, a parseable SMILES that regenerates that InChIKey and
the authoritative formula below.  Reported precursor/adduct fields are never
changed; corrected active fields are written transactionally only with
``--apply``.  A complete before-image and receipt are produced for rollback.
Dry runs use a separate receipt and never create or replace the rollback image
or the applied receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.expected_compounds import _normalise_adduct
from matching.library_quality import formula_exact_mass, normalize_formula_text
from utils.adduct import ion_mz


DEFAULT_DATABASE = Path(
    "/Users/zizi/Desktop/Mssystem_shared_data/local_library.db"
)
DEFAULT_RECEIPT_DIR = ROOT / "docs" / "audits"
DEFAULT_BACKUP_DIR = Path(
    "/Users/zizi/Desktop/Mssystem/.v12_restore_points/"
    "20260818_pre_publishability_repair/targeted_db_rows"
)

REPAIR_EVENT_ID = "v12_verified_identity_formula_repair_20260818"
REPAIR_RECEIPT_ID = "mssystem-v12-verified-target-identity-conflicts-20260818"
APPLIED_RECEIPT_NAME = (
    "MSSYSTEM_V12_VERIFIED_TARGET_IDENTITY_REPAIR_20260818.json"
)
DRY_RUN_RECEIPT_NAME = (
    "MSSYSTEM_V12_VERIFIED_TARGET_IDENTITY_REPAIR_DRY_RUN_20260818.json"
)
BACKUP_NAME = "TargetCompounds_8_verified_identity_rows_before_20260818.json"

REPAIRS = {
    "F005R-38abd2f6992ded3d1580": {
        "name": "Methyl dihydrojasmonate",
        "inchikey": "KVWWIYGFBYDJQC-UHFFFAOYSA-N",
        "formula": "C13H22O3",
        "exact_mass": 226.15689456,
        "pubchem_cid": 102861,
    },
    "F005R-39e53b2b65aae785a539": {
        "name": "Lariciresinol",
        "inchikey": "MHXCIKYXNYCMHY-AUSJPIAWSA-N",
        "formula": "C20H24O6",
        "exact_mass": 360.15728848,
        "pubchem_cid": 332427,
    },
    "F005R-47a0d6418ed084070f70": {
        "name": "Forsythoside E",
        "inchikey": "QIMGUQIHCNDUKU-OJJLXHDHSA-N",
        "formula": "C20H30O12",
        "exact_mass": 462.17372639,
        "pubchem_cid": 69634125,
    },
    "F005R-48c0e80824bd65b20b29": {
        "name": "Forsythoside A",
        "inchikey": "DTOUWTJYUCZJQD-UJERWXFOSA-N",
        "formula": "C29H36O15",
        "exact_mass": 624.20542044,
        "pubchem_cid": 5281773,
        "canonical_adduct": "[M+Na]+",
    },
    "F005R-5266d2b0b4485f685b27": {
        "name": "Forsythoside A",
        "inchikey": "DTOUWTJYUCZJQD-UJERWXFOSA-N",
        "formula": "C29H36O15",
        "exact_mass": 624.20542044,
        "pubchem_cid": 5281773,
        "canonical_adduct": "[M+Na]+",
    },
    "F005R-bc3b71125305b7fbf58a": {
        "name": "Derhamnosyl suspensaside",
        "inchikey": "FQHPKISHHZQDFI-WZKZINOLSA-N",
        "formula": "C23H26O12",
        "exact_mass": 494.14242626,
        "pubchem_cid": 177826589,
    },
    "F005R-e0620ff1384fe0198ecf": {
        "name": "Eugenol",
        "inchikey": "RRAFCDWBNXTKKO-UHFFFAOYSA-N",
        "formula": "C10H12O2",
        "exact_mass": 164.083729621,
        "pubchem_cid": 3314,
    },
    "F006R-8ab8378779ebb1dd9513": {
        "name": "Hyperoside",
        "inchikey": "OVSQVDMCBVZWGM-DTGCRPNFSA-N",
        "formula": "C21H20O12",
        "exact_mass": 464.09547607,
        "pubchem_cid": 5281643,
    },
}


def _open(path: Path, read_only: bool) -> sqlite3.Connection:
    uri = f"file:{path.resolve()}?mode={'ro' if read_only else 'rw'}"
    connection = sqlite3.connect(uri, uri=True, timeout=30)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout=30000")
    return connection


def _rows(connection: sqlite3.Connection) -> list[dict[str, Any]]:
    marks = ",".join("?" for _ in REPAIRS)
    rows = connection.execute(
        f"SELECT * FROM TargetCompounds WHERE target_id IN ({marks}) "
        "ORDER BY target_id",
        tuple(REPAIRS),
    ).fetchall()
    if len(rows) != len(REPAIRS):
        found = {str(row["target_id"]) for row in rows}
        raise RuntimeError(f"Missing repair targets: {sorted(set(REPAIRS) - found)}")
    return [dict(row) for row in rows]


def _fingerprint(rows: list[dict[str, Any]]) -> str:
    canonical = json.dumps(
        rows, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _write_json(
    path: Path,
    payload: dict[str, Any],
    *,
    exclusive: bool,
) -> None:
    """Write JSON without silently replacing a protected artifact."""
    serialised = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    if exclusive:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(serialised)
        return
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(serialised, encoding="utf-8")
    temporary.replace(path)


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Invalid protected JSON artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise RuntimeError(f"Protected JSON artifact is not an object: {path}")
    return payload


def _has_repair_marker(row: dict[str, Any]) -> bool:
    try:
        payload = json.loads(str(row.get("quality_flags_json") or "{}"))
    except json.JSONDecodeError:
        return False
    if not isinstance(payload, dict):
        return False
    return any(
        isinstance(item, dict) and item.get("repair_id") == REPAIR_EVENT_ID
        for item in payload.get("v12_repairs", [])
    )


def _same_value(actual: Any, expected: Any) -> bool:
    if actual is None or expected is None:
        return actual is expected
    if isinstance(expected, float):
        try:
            return abs(float(actual) - expected) <= 1e-9
        except (TypeError, ValueError):
            return False
    return actual == expected


def _active_values_match(
    row: dict[str, Any],
    item: dict[str, Any],
) -> bool:
    return all((
        _same_value(row.get("molecular_formula"), item["after_formula"]),
        _same_value(row.get("exact_mass"), item["after_exact_mass"]),
        _same_value(row.get("precursor_mz"), item["after_precursor_mz"]),
        _same_value(row.get("adduct"), item["after_adduct"]),
    ))


def _validate_backup(
    backup_path: Path,
    *,
    expected_fingerprint: str | None = None,
) -> dict[str, Any]:
    payload = _read_json_object(backup_path)
    rows = payload.get("rows")
    if not isinstance(rows, list) or len(rows) != len(REPAIRS):
        raise RuntimeError(f"Protected rollback image has the wrong row count: {backup_path}")
    stored_fingerprint = str(payload.get("row_fingerprint_sha256") or "")
    calculated_fingerprint = _fingerprint(rows)
    if stored_fingerprint != calculated_fingerprint:
        raise RuntimeError(f"Protected rollback image fingerprint mismatch: {backup_path}")
    if expected_fingerprint and stored_fingerprint != expected_fingerprint:
        raise RuntimeError(
            "Existing rollback image belongs to a different database state: "
            f"{backup_path}"
        )
    target_ids = {str(row.get("target_id") or "") for row in rows}
    if target_ids != set(REPAIRS):
        raise RuntimeError(f"Protected rollback image target set mismatch: {backup_path}")
    return payload


def _load_existing_applied_receipt(
    receipt_path: Path,
    backup_path: Path,
    *,
    current_fingerprint: str,
) -> dict[str, Any]:
    if not receipt_path.is_file() or not backup_path.is_file():
        raise RuntimeError(
            "Rows contain the applied repair marker, but the protected applied "
            "receipt or rollback image is missing"
        )
    receipt = _read_json_object(receipt_path)
    if receipt.get("mode") != "applied":
        raise RuntimeError(f"Protected receipt is not an applied receipt: {receipt_path}")
    if receipt.get("repair_id") != REPAIR_RECEIPT_ID:
        raise RuntimeError(f"Protected receipt repair ID mismatch: {receipt_path}")
    if receipt.get("after_fingerprint_sha256") != current_fingerprint:
        raise RuntimeError(f"Applied receipt does not match the current rows: {receipt_path}")
    backup = _validate_backup(backup_path)
    if (
        receipt.get("before_fingerprint_sha256")
        != backup.get("row_fingerprint_sha256")
    ):
        raise RuntimeError("Applied receipt and rollback image disagree")
    return receipt


def _quality_payload(raw: Any, repair_record: dict[str, Any]) -> str:
    try:
        payload = json.loads(str(raw or "{}"))
    except json.JSONDecodeError:
        payload = {"legacy_unparsed_quality_flags": str(raw or "")}
    if isinstance(payload, list):
        payload = {"legacy_flags": payload}
    elif not isinstance(payload, dict):
        payload = {"legacy_quality_value": payload}
    repairs = payload.setdefault("v12_repairs", [])
    repair_id = repair_record["repair_id"]
    if not any(
        isinstance(item, dict) and item.get("repair_id") == repair_id
        for item in repairs
    ):
        repairs.append(repair_record)
    return json.dumps(payload, ensure_ascii=False, sort_keys=True)


def _validate_and_plan(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolDescriptors
    except ImportError as exc:
        raise RuntimeError("RDKit is required for guarded identity repair") from exc
    RDLogger.DisableLog("rdApp.*")

    plan = []
    for before in rows:
        target_id = str(before["target_id"])
        expected = REPAIRS[target_id]
        if not bool(before["structure_identity_verified"]):
            raise RuntimeError(f"{target_id}: structure_identity_verified is false")
        if str(before["inchikey"] or "").upper() != expected["inchikey"]:
            raise RuntimeError(f"{target_id}: unexpected InChIKey")
        mol = Chem.MolFromSmiles(str(before["smiles"] or ""))
        if mol is None:
            raise RuntimeError(f"{target_id}: invalid SMILES")
        generated_key = Chem.MolToInchiKey(mol)
        generated_formula = normalize_formula_text(
            rdMolDescriptors.CalcMolFormula(mol)
        )
        if generated_key != expected["inchikey"]:
            raise RuntimeError(
                f"{target_id}: SMILES InChIKey {generated_key} does not match expected"
            )
        if generated_formula != expected["formula"]:
            raise RuntimeError(
                f"{target_id}: SMILES formula {generated_formula} does not match expected"
            )
        calculated_mass = formula_exact_mass(expected["formula"])
        if calculated_mass is None or abs(calculated_mass - expected["exact_mass"]) > 2e-6:
            raise RuntimeError(f"{target_id}: exact-mass cross-check failed")

        canonical_adduct = str(expected.get("canonical_adduct") or "")
        if not canonical_adduct:
            canonical_adduct = _normalise_adduct(before["adduct"])
        active_precursor = (
            ion_mz(expected["exact_mass"], canonical_adduct)
            if canonical_adduct else None
        )
        source_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{expected['pubchem_cid']}"
        )
        repair_record = {
            "repair_id": REPAIR_EVENT_ID,
            "source": source_url,
            "basis": (
                "verified SMILES, full InChIKey, RDKit molecular formula, and "
                "PubChem identity record agree; legacy formula field conflicts"
            ),
            "before": {
                "molecular_formula": before["molecular_formula"],
                "exact_mass": before["exact_mass"],
                "precursor_mz": before["precursor_mz"],
                "adduct": before["adduct"],
            },
            "reported_fields_preserved": {
                "reported_precursor_mz": before["reported_precursor_mz"],
                "reported_adduct": before["reported_adduct"],
            },
        }
        provenance_note = (
            "v12 verified identity-field repair 2026-08-18; "
            f"PubChem CID {expected['pubchem_cid']}; original reported fields retained"
        )
        old_provenance = str(before["provenance"] or "").strip()
        provenance = old_provenance
        if provenance_note not in old_provenance:
            provenance = f"{old_provenance} | {provenance_note}".strip(" |")
        plan.append({
            "target_id": target_id,
            "compound_name": before["compound_name"],
            "before_formula": before["molecular_formula"],
            "after_formula": expected["formula"],
            "after_exact_mass": expected["exact_mass"],
            "after_adduct": canonical_adduct,
            "after_precursor_mz": active_precursor,
            "reported_precursor_mz_preserved": before["reported_precursor_mz"],
            "reported_adduct_preserved": before["reported_adduct"],
            "pubchem_cid": expected["pubchem_cid"],
            "pubchem_url": source_url,
            "quality_flags_json": _quality_payload(
                before["quality_flags_json"], repair_record
            ),
            "provenance": provenance,
        })
    return plan


def run(
    database: Path,
    receipt_dir: Path,
    backup_dir: Path,
    *,
    apply: bool,
) -> dict[str, Any]:
    if not database.is_file():
        raise FileNotFoundError(database)
    receipt_dir.mkdir(parents=True, exist_ok=True)
    if apply:
        backup_dir.mkdir(parents=True, exist_ok=True)
    backup_path = backup_dir / BACKUP_NAME
    applied_receipt_path = receipt_dir / APPLIED_RECEIPT_NAME
    dry_run_receipt_path = receipt_dir / DRY_RUN_RECEIPT_NAME

    with _open(database, read_only=not apply) as connection:
        before_rows = _rows(connection)
        before_fingerprint = _fingerprint(before_rows)
        plan = _validate_and_plan(before_rows)
        rows_by_id = {str(row["target_id"]): row for row in before_rows}
        markers = [
            _has_repair_marker(rows_by_id[item["target_id"]])
            for item in plan
        ]
        applied_states = [
            marker
            and _active_values_match(rows_by_id[item["target_id"]], item)
            for marker, item in zip(markers, plan)
        ]
        if any(markers) and not all(applied_states):
            raise RuntimeError(
                "Partial or inconsistent verified-identity repair state; refusing "
                "to write the database or protected artifacts"
            )

        already_applied = bool(applied_states) and all(applied_states)
        if apply and already_applied:
            return _load_existing_applied_receipt(
                applied_receipt_path,
                backup_path,
                current_fingerprint=before_fingerprint,
            )

        if apply:
            if applied_receipt_path.exists():
                raise RuntimeError(
                    "An applied receipt already exists but the database is not in "
                    "the corresponding applied state"
                )
            backup_payload = {
                "schema_version": "1.0.0",
                "database": str(database.resolve()),
                "table": "TargetCompounds",
                "row_count": len(before_rows),
                "row_fingerprint_sha256": before_fingerprint,
                "rows": before_rows,
            }
            if backup_path.exists():
                _validate_backup(
                    backup_path,
                    expected_fingerprint=before_fingerprint,
                )
            else:
                _write_json(backup_path, backup_payload, exclusive=True)

            connection.execute("BEGIN IMMEDIATE")
            try:
                for item in plan:
                    connection.execute(
                        """
                        UPDATE TargetCompounds
                        SET molecular_formula=?, exact_mass=?, precursor_mz=?,
                            adduct=?, quality_flags_json=?, provenance=?
                        WHERE target_id=?
                        """,
                        (
                            item["after_formula"], item["after_exact_mass"],
                            item["after_precursor_mz"], item["after_adduct"],
                            item["quality_flags_json"], item["provenance"],
                            item["target_id"],
                        ),
                    )
                connection.commit()
            except Exception:
                connection.rollback()
                raise
            after_rows = _rows(connection)
        else:
            after_rows = before_rows

    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    receipt = {
        "schema_version": "1.0.0",
        "repair_id": REPAIR_RECEIPT_ID,
        "generated_at_utc": generated,
        "mode": "applied" if apply else "dry_run",
        "database": str(database.resolve()),
        "transactional": True,
        "reported_fields_preserved": True,
        "before_fingerprint_sha256": before_fingerprint,
        "after_fingerprint_sha256": _fingerprint(after_rows),
        "backup_path": str(backup_path.resolve()) if apply else None,
        "protected_artifacts_modified": bool(apply),
        "rows": [
            {key: value for key, value in item.items()
             if key not in {"quality_flags_json", "provenance"}}
            for item in plan
        ],
    }
    receipt_path = applied_receipt_path if apply else dry_run_receipt_path
    _write_json(
        receipt_path,
        receipt,
        exclusive=apply,
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--receipt-dir", type=Path, default=DEFAULT_RECEIPT_DIR)
    parser.add_argument("--backup-dir", type=Path, default=DEFAULT_BACKUP_DIR)
    parser.add_argument(
        "--apply", action="store_true",
        help="Apply the guarded transaction; default is a read-only dry run.",
    )
    args = parser.parse_args()
    receipt = run(
        args.database, args.receipt_dir, args.backup_dir, apply=args.apply
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
