#!/usr/bin/env python3
"""Transactionally install a staged YDF target layer into an MSsystem DB.

This migration copies only ``TargetCompounds`` and ``DiagnosticIons``.  The
64,362 experimental spectra in the staged package were selected from the
existing shared library and therefore must not be copied back as duplicates.
The full evidence, provenance, exclusions, and target-spectrum links remain in
the staged sidecar database.

Dry-run is the default.  Applying directly to the canonical shared database
requires both ``--apply`` and ``--allow-canonical-shared``.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator


CANONICAL_SHARED_DB = Path("/Users/zizi/Desktop/Mssystem_shared_data/local_library.db")
EXPECTED_DATASET = "YDF_targeted_20260812"

TARGET_SCHEMA = """
CREATE TABLE IF NOT EXISTS TargetCompounds (
    id INTEGER PRIMARY KEY,
    target_id TEXT UNIQUE NOT NULL,
    source_file TEXT,
    source_content_sha256 TEXT,
    dataset_version TEXT,
    source_row_index INTEGER,
    dataset_source TEXT,
    compound_name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    inchi TEXT,
    inchikey TEXT,
    reported_inchikey TEXT,
    structure_identity_verified INTEGER DEFAULT 0,
    structure_identity_source TEXT,
    exact_mass REAL,
    precursor_mz REAL,
    adduct TEXT,
    charge INTEGER,
    reported_precursor_mz REAL,
    reported_adduct TEXT,
    reported_retention_time REAL,
    reported_retention_time_unit TEXT,
    reported_collision_energy_raw TEXT,
    instrument TEXT,
    target_searchability TEXT,
    plant_source TEXT,
    provenance TEXT,
    doi TEXT,
    dataset_accession TEXT,
    spectrum_accession TEXT,
    source_tier TEXT NOT NULL DEFAULT 'literature_presence_only',
    quality_flags_json TEXT,
    created_at TEXT
);
CREATE TABLE IF NOT EXISTS DiagnosticIons (
    id INTEGER PRIMARY KEY,
    target_id INTEGER NOT NULL,
    fragment_mz REAL NOT NULL,
    ordinal INTEGER,
    FOREIGN KEY(target_id) REFERENCES TargetCompounds(id)
);
CREATE INDEX IF NOT EXISTS ix_target_dataset ON TargetCompounds(dataset_source);
CREATE INDEX IF NOT EXISTS ix_target_precursor ON TargetCompounds(precursor_mz);
CREATE INDEX IF NOT EXISTS ix_target_name_formula ON TargetCompounds(compound_name, molecular_formula);
CREATE INDEX IF NOT EXISTS ix_diagnostic_target ON DiagnosticIons(target_id);
CREATE INDEX IF NOT EXISTS ix_diagnostic_mz ON DiagnosticIons(fragment_mz);
"""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_only(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row[1]) for row in connection.execute(f'PRAGMA table_info("{table}")')]


def ensure_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(TARGET_SCHEMA)


@contextmanager
def migration_lock(target_db: Path) -> Iterator[None]:
    lock_path = target_db.parent / ".ydf_target_migration.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def source_payload(source_db: Path) -> dict[str, Any]:
    with read_only(source_db) as source:
        dataset_rows = source.execute(
            "SELECT DISTINCT dataset_source FROM TargetCompounds"
        ).fetchall()
        datasets = {str(row[0] or "") for row in dataset_rows}
        if datasets != {EXPECTED_DATASET}:
            raise RuntimeError(f"unexpected staged target datasets: {sorted(datasets)}")
        columns = [column for column in table_columns(source, "TargetCompounds") if column != "id"]
        targets = [dict(row) for row in source.execute(
            f"SELECT id,{','.join(columns)} FROM TargetCompounds ORDER BY id"
        )]
        source_id_to_target = {int(row["id"]): row["target_id"] for row in targets}
        diagnostics = [
            {
                "target_id": source_id_to_target[int(row["target_id"])],
                "fragment_mz": float(row["fragment_mz"]),
                "ordinal": row["ordinal"],
            }
            for row in source.execute(
                "SELECT target_id,fragment_mz,ordinal FROM DiagnosticIons ORDER BY target_id,ordinal"
            )
        ]
        for row in targets:
            row.pop("id", None)
        return {
            "columns": columns,
            "targets": targets,
            "diagnostics": diagnostics,
            "source_sha256": sha256_file(source_db),
        }


def inspect_target(target_db: Path, target_ids: set[str]) -> dict[str, Any]:
    if not target_db.exists():
        return {
            "target_exists": False,
            "existing_dataset_targets": 0,
            "existing_dataset_diagnostic_ions": 0,
            "cross_dataset_target_id_conflicts": 0,
        }
    with read_only(target_db) as target:
        tables = {
            str(row[0]) for row in target.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        if "TargetCompounds" not in tables:
            return {
                "target_exists": True,
                "existing_dataset_targets": 0,
                "existing_dataset_diagnostic_ions": 0,
                "cross_dataset_target_id_conflicts": 0,
            }
        existing_ids = [
            int(row[0]) for row in target.execute(
                "SELECT id FROM TargetCompounds WHERE dataset_source=?",
                (EXPECTED_DATASET,),
            )
        ]
        diagnostics = 0
        if existing_ids and "DiagnosticIons" in tables:
            for start in range(0, len(existing_ids), 500):
                block = existing_ids[start:start + 500]
                marks = ",".join("?" for _ in block)
                diagnostics += int(target.execute(
                    f"SELECT COUNT(*) FROM DiagnosticIons WHERE target_id IN ({marks})",
                    block,
                ).fetchone()[0])
        conflicts = 0
        for start in range(0, len(target_ids), 500):
            block = sorted(target_ids)[start:start + 500]
            marks = ",".join("?" for _ in block)
            conflicts += int(target.execute(
                f"SELECT COUNT(*) FROM TargetCompounds WHERE target_id IN ({marks}) "
                "AND coalesce(dataset_source,'')<>?",
                [*block, EXPECTED_DATASET],
            ).fetchone()[0])
        return {
            "target_exists": True,
            "existing_dataset_targets": len(existing_ids),
            "existing_dataset_diagnostic_ions": diagnostics,
            "cross_dataset_target_id_conflicts": conflicts,
        }


def backup_database(source: Path, destination: Path) -> None:
    with read_only(source) as origin, sqlite3.connect(str(destination)) as backup:
        origin.backup(backup, pages=10000, sleep=0.05)


def migrate(
    source_db: Path,
    target_db: Path,
    *,
    apply: bool = False,
    allow_canonical_shared: bool = False,
    create_backup: bool = True,
) -> dict[str, Any]:
    source_db = source_db.expanduser().resolve()
    target_db = target_db.expanduser().resolve()
    if source_db == target_db:
        raise RuntimeError("source and target database must be different")
    if not source_db.is_file():
        raise FileNotFoundError(source_db)
    payload = source_payload(source_db)
    target_ids = {str(row["target_id"]) for row in payload["targets"]}
    before = inspect_target(target_db, target_ids)
    if before["cross_dataset_target_id_conflicts"]:
        raise RuntimeError(f"cross-dataset target_id conflicts: {before['cross_dataset_target_id_conflicts']}")
    receipt: dict[str, Any] = {
        "status": "dry_run" if not apply else "pending",
        "source_db": str(source_db),
        "source_sha256": payload["source_sha256"],
        "target_db": str(target_db),
        "dataset_source": EXPECTED_DATASET,
        "desired_targets": len(payload["targets"]),
        "desired_diagnostic_ions": len(payload["diagnostics"]),
        "spectra_copied": 0,
        "spectra_copy_policy": "not copied; staged experimental spectra originate from the existing shared library",
        "before": before,
        "backup": None,
    }
    if not apply:
        return receipt
    if target_db == CANONICAL_SHARED_DB.resolve() and not allow_canonical_shared:
        raise RuntimeError(
            "refusing canonical shared DB; add --allow-canonical-shared only after stopping writers and reviewing dry-run"
        )

    target_db.parent.mkdir(parents=True, exist_ok=True)
    with migration_lock(target_db):
        if target_db.exists() and create_backup:
            timestamp = time.strftime("%Y%m%dT%H%M%S")
            backup = target_db.with_name(f"{target_db.name}.pre_ydf_targeted_{timestamp}.bak")
            if backup.exists():
                raise FileExistsError(backup)
            backup_database(target_db, backup)
            receipt["backup"] = {
                "path": str(backup),
                "size_bytes": backup.stat().st_size,
                "sha256": sha256_file(backup),
            }

        connection = sqlite3.connect(str(target_db))
        try:
            connection.execute("PRAGMA busy_timeout=60000")
            connection.execute("PRAGMA foreign_keys=ON")
            ensure_schema(connection)
            target_columns = set(table_columns(connection, "TargetCompounds"))
            missing = [column for column in payload["columns"] if column not in target_columns]
            if missing:
                raise RuntimeError(f"target TargetCompounds lacks required columns: {missing}")
            connection.execute("BEGIN IMMEDIATE")
            old_ids = [
                int(row[0]) for row in connection.execute(
                    "SELECT id FROM TargetCompounds WHERE dataset_source=?",
                    (EXPECTED_DATASET,),
                )
            ]
            for start in range(0, len(old_ids), 500):
                block = old_ids[start:start + 500]
                marks = ",".join("?" for _ in block)
                connection.execute(
                    f"DELETE FROM DiagnosticIons WHERE target_id IN ({marks})", block
                )
            connection.execute(
                "DELETE FROM TargetCompounds WHERE dataset_source=?",
                (EXPECTED_DATASET,),
            )
            columns = payload["columns"]
            marks = ",".join("?" for _ in columns)
            quoted = ",".join(f'"{column}"' for column in columns)
            connection.executemany(
                f"INSERT INTO TargetCompounds({quoted}) VALUES ({marks})",
                [tuple(row[column] for column in columns) for row in payload["targets"]],
            )
            target_id_to_pk = {
                str(row[1]): int(row[0])
                for row in connection.execute(
                    "SELECT id,target_id FROM TargetCompounds WHERE dataset_source=?",
                    (EXPECTED_DATASET,),
                )
            }
            connection.executemany(
                "INSERT INTO DiagnosticIons(target_id,fragment_mz,ordinal) VALUES (?,?,?)",
                [
                    (target_id_to_pk[row["target_id"]], row["fragment_mz"], row["ordinal"])
                    for row in payload["diagnostics"]
                ],
            )
            installed_targets = int(connection.execute(
                "SELECT COUNT(*) FROM TargetCompounds WHERE dataset_source=?",
                (EXPECTED_DATASET,),
            ).fetchone()[0])
            installed_ions = int(connection.execute(
                """SELECT COUNT(*) FROM DiagnosticIons d JOIN TargetCompounds t
                ON t.id=d.target_id WHERE t.dataset_source=?""",
                (EXPECTED_DATASET,),
            ).fetchone()[0])
            if installed_targets != len(payload["targets"]) or installed_ions != len(payload["diagnostics"]):
                raise RuntimeError(
                    f"post-insert counts mismatch: {installed_targets}/{installed_ions}"
                )
            connection.execute("COMMIT")
        except Exception:
            try:
                connection.execute("ROLLBACK")
            except sqlite3.Error:
                pass
            raise
        finally:
            connection.close()

    after = inspect_target(target_db, target_ids)
    receipt.update({"status": "applied", "after": after})
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-db", type=Path, required=True)
    parser.add_argument("--target-db", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--allow-canonical-shared", action="store_true")
    parser.add_argument("--skip-backup", action="store_true")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    receipt = migrate(
        args.source_db,
        args.target_db,
        apply=args.apply,
        allow_canonical_shared=args.allow_canonical_shared,
        create_backup=not args.skip_backup,
    )
    if args.receipt:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
