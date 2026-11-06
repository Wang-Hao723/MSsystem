from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from scripts import repair_verified_target_identity_conflicts as repair


def _fixture_database(path: Path) -> None:
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            CREATE TABLE TargetCompounds (
                target_id TEXT PRIMARY KEY,
                compound_name TEXT,
                molecular_formula TEXT,
                smiles TEXT,
                inchikey TEXT,
                exact_mass REAL,
                precursor_mz REAL,
                adduct TEXT,
                structure_identity_verified INTEGER,
                reported_precursor_mz REAL,
                reported_adduct TEXT,
                quality_flags_json TEXT,
                provenance TEXT
            )
            """
        )
        connection.execute(
            """
            INSERT INTO TargetCompounds VALUES (
                'T1', 'Ethanol', 'C2H4O', 'CCO',
                'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', NULL, 47.0, '[M+H]+',
                1, 47.0, '[M+H]+',
                '{"method_signature":"M-test","evidence_count":1}',
                'fixture source'
            )
            """
        )


def test_dry_run_cannot_replace_applied_receipt_or_backup(
    tmp_path: Path,
    monkeypatch,
) -> None:
    database = tmp_path / "library.sqlite"
    receipt_dir = tmp_path / "receipts"
    backup_dir = tmp_path / "backups"
    _fixture_database(database)

    exact_mass = repair.formula_exact_mass("C2H6O")
    assert exact_mass is not None
    monkeypatch.setattr(
        repair,
        "REPAIRS",
        {
            "T1": {
                "name": "Ethanol",
                "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "formula": "C2H6O",
                "exact_mass": exact_mass,
                "pubchem_cid": 702,
            }
        },
    )

    dry_receipt = repair.run(
        database, receipt_dir, backup_dir, apply=False
    )
    applied_path = receipt_dir / repair.APPLIED_RECEIPT_NAME
    dry_path = receipt_dir / repair.DRY_RUN_RECEIPT_NAME
    backup_path = backup_dir / repair.BACKUP_NAME
    assert dry_receipt["mode"] == "dry_run"
    assert dry_receipt["backup_path"] is None
    assert dry_receipt["protected_artifacts_modified"] is False
    assert dry_path.is_file()
    assert not applied_path.exists()
    assert not backup_dir.exists()

    applied_receipt = repair.run(
        database, receipt_dir, backup_dir, apply=True
    )
    assert applied_receipt["mode"] == "applied"
    assert applied_path.is_file()
    assert backup_path.is_file()
    backup = json.loads(backup_path.read_text(encoding="utf-8"))
    assert backup["rows"][0]["molecular_formula"] == "C2H4O"
    applied_bytes = applied_path.read_bytes()
    backup_bytes = backup_path.read_bytes()
    database_mtime = database.stat().st_mtime_ns

    second_dry_receipt = repair.run(
        database, receipt_dir, backup_dir, apply=False
    )
    assert second_dry_receipt["mode"] == "dry_run"
    assert applied_path.read_bytes() == applied_bytes
    assert backup_path.read_bytes() == backup_bytes

    repeated_apply = repair.run(
        database, receipt_dir, backup_dir, apply=True
    )
    assert repeated_apply == applied_receipt
    assert database.stat().st_mtime_ns == database_mtime
    assert applied_path.read_bytes() == applied_bytes
    assert backup_path.read_bytes() == backup_bytes


def test_applied_state_requires_matching_protected_artifacts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    database = tmp_path / "library.sqlite"
    receipt_dir = tmp_path / "receipts"
    backup_dir = tmp_path / "backups"
    _fixture_database(database)
    exact_mass = repair.formula_exact_mass("C2H6O")
    assert exact_mass is not None
    monkeypatch.setattr(
        repair,
        "REPAIRS",
        {
            "T1": {
                "name": "Ethanol",
                "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "formula": "C2H6O",
                "exact_mass": exact_mass,
                "pubchem_cid": 702,
            }
        },
    )
    repair.run(database, receipt_dir, backup_dir, apply=True)
    backup_path = backup_dir / repair.BACKUP_NAME
    backup = json.loads(backup_path.read_text(encoding="utf-8"))
    backup["rows"][0]["molecular_formula"] = "corrupted"
    backup_path.write_text(json.dumps(backup), encoding="utf-8")

    try:
        repair.run(database, receipt_dir, backup_dir, apply=True)
    except RuntimeError as exc:
        assert "fingerprint mismatch" in str(exc)
    else:
        raise AssertionError("corrupted protected backup was accepted")
