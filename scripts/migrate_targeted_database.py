#!/usr/bin/env python3
"""Transactionally install a staged targeted MSsystem database into a shared library.

Generalization of migrate_ydf_targeted_database.py: the dataset identity comes
from ``--dataset-source`` (read from the staged DB, must be unique) instead of a
hard-coded constant, and ``--retire-dataset`` optionally removes superseded
dataset_source rows (e.g. old ``user_contrib:...`` imports) in the same
transaction after the new layer is installed.

Copies only TargetCompounds + DiagnosticIons.  Staged experimental spectra were
selected from the shared library and must not be copied back.  Dry-run is the
default; applying to the canonical shared DB needs --apply plus
--allow-canonical-shared.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

CANONICAL_SHARED_DB = Path("/Users/zizi/Desktop/Mssystem_shared_data/local_library.db")

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
def migration_lock(target_db: Path, dataset: str) -> Iterator[None]:
    safe = "".join(c if c.isalnum() else "_" for c in dataset)[:60]
    lock_path = target_db.parent / f".target_migration_{safe}.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def source_payload(source_db: Path, dataset: str | None) -> dict[str, Any]:
    with read_only(source_db) as source:
        dataset_rows = source.execute("SELECT DISTINCT dataset_source FROM TargetCompounds").fetchall()
        datasets = {str(row[0] or "") for row in dataset_rows}
        if dataset is None:
            if len(datasets) != 1:
                raise RuntimeError(f"staged DB holds multiple datasets {sorted(datasets)}; pass --dataset-source")
            dataset = next(iter(datasets))
        elif datasets != {dataset}:
            raise RuntimeError(f"unexpected staged target datasets: {sorted(datasets)} (wanted {dataset})")
        columns = [column for column in table_columns(source, "TargetCompounds") if column != "id"]
        targets = [dict(row) for row in source.execute(
            f"SELECT id,{','.join(columns)} FROM TargetCompounds ORDER BY id")]
        source_id_to_target = {int(row["id"]): row["target_id"] for row in targets}
        diagnostics = [
            {"target_id": source_id_to_target[int(row["target_id"])],
             "fragment_mz": float(row["fragment_mz"]), "ordinal": row["ordinal"]}
            for row in source.execute(
                "SELECT target_id,fragment_mz,ordinal FROM DiagnosticIons ORDER BY target_id,ordinal")
        ]
        for row in targets:
            row.pop("id", None)
        return {"dataset": dataset, "columns": columns, "targets": targets,
                "diagnostics": diagnostics, "source_sha256": sha256_file(source_db)}


def inspect_target(target_db: Path, dataset: str, target_ids: set[str]) -> dict[str, Any]:
    if not target_db.exists():
        return {"target_exists": False, "existing_dataset_targets": 0,
                "existing_dataset_diagnostic_ions": 0, "cross_dataset_target_id_conflicts": 0}
    with read_only(target_db) as target:
        tables = {str(r[0]) for r in target.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "TargetCompounds" not in tables:
            return {"target_exists": True, "existing_dataset_targets": 0,
                    "existing_dataset_diagnostic_ions": 0, "cross_dataset_target_id_conflicts": 0}
        existing_ids = [int(r[0]) for r in target.execute(
            "SELECT id FROM TargetCompounds WHERE dataset_source=?", (dataset,))]
        diagnostics = 0
        if existing_ids and "DiagnosticIons" in tables:
            for start in range(0, len(existing_ids), 500):
                block = existing_ids[start:start + 500]
                marks = ",".join("?" for _ in block)
                diagnostics += int(target.execute(
                    f"SELECT COUNT(*) FROM DiagnosticIons WHERE target_id IN ({marks})", block).fetchone()[0])
        conflicts = 0
        for start in range(0, len(target_ids), 500):
            block = sorted(target_ids)[start:start + 500]
            marks = ",".join("?" for _ in block)
            conflicts += int(target.execute(
                f"SELECT COUNT(*) FROM TargetCompounds WHERE target_id IN ({marks}) "
                "AND coalesce(dataset_source,'')<>?", [*block, dataset]).fetchone()[0])
        return {"target_exists": True, "existing_dataset_targets": len(existing_ids),
                "existing_dataset_diagnostic_ions": diagnostics,
                "cross_dataset_target_id_conflicts": conflicts}


def backup_database(source: Path, destination: Path) -> None:
    with read_only(source) as origin, sqlite3.connect(str(destination)) as backup:
        origin.backup(backup, pages=10000, sleep=0.05)


def delete_dataset(connection: sqlite3.Connection, dataset: str) -> tuple[int, int]:
    old_ids = [int(r[0]) for r in connection.execute(
        "SELECT id FROM TargetCompounds WHERE dataset_source=?", (dataset,))]
    for start in range(0, len(old_ids), 500):
        block = old_ids[start:start + 500]
        marks = ",".join("?" for _ in block)
        connection.execute(f"DELETE FROM DiagnosticIons WHERE target_id IN ({marks})", block)
    connection.execute("DELETE FROM TargetCompounds WHERE dataset_source=?", (dataset,))
    return len(old_ids), 0


def migrate(source_db: Path, target_db: Path, *, dataset: str | None, apply: bool,
            allow_canonical_shared: bool, create_backup: bool,
            retire_datasets: list[str], receipt_path: Path | None) -> dict[str, Any]:
    source_db = source_db.expanduser().resolve()
    target_db = target_db.expanduser().resolve()
    if source_db == target_db:
        raise RuntimeError("source and target database must be different")
    if not source_db.is_file():
        raise FileNotFoundError(source_db)
    payload = source_payload(source_db, dataset)
    dataset = payload["dataset"]
    target_ids = {str(row["target_id"]) for row in payload["targets"]}
    before = inspect_target(target_db, dataset, target_ids)
    retire_preview = {}
    if retire_datasets and target_db.exists():
        with read_only(target_db) as target:
            for ds in retire_datasets:
                retire_preview[ds] = int(target.execute(
                    "SELECT COUNT(*) FROM TargetCompounds WHERE dataset_source=?", (ds,)).fetchone()[0])
    if before["cross_dataset_target_id_conflicts"]:
        raise RuntimeError(f"cross-dataset target_id conflicts: {before['cross_dataset_target_id_conflicts']}")
    receipt: dict[str, Any] = {
        "status": "dry_run" if not apply else "pending",
        "source_db": str(source_db), "source_sha256": payload["source_sha256"],
        "target_db": str(target_db), "dataset_source": dataset,
        "desired_targets": len(payload["targets"]),
        "desired_diagnostic_ions": len(payload["diagnostics"]),
        "spectra_copied": 0,
        "spectra_copy_policy": "not copied; staged experimental spectra originate from the existing shared library",
        "before": before,
        "retire_datasets": retire_preview,
        "backup": None,
    }
    if not apply:
        if receipt_path:
            receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=1))
        return receipt
    if target_db == CANONICAL_SHARED_DB.resolve() and not allow_canonical_shared:
        raise RuntimeError("refusing canonical shared DB; add --allow-canonical-shared only after "
                           "stopping writers and reviewing dry-run")

    target_db.parent.mkdir(parents=True, exist_ok=True)
    with migration_lock(target_db, dataset):
        if target_db.exists() and create_backup:
            timestamp = time.strftime("%Y%m%dT%H%M%S")
            safe = "".join(c if c.isalnum() else "_" for c in dataset)[:40]
            backup = target_db.with_name(f"{target_db.name}.pre_{safe}_{timestamp}.bak")
            if backup.exists():
                raise FileExistsError(backup)
            backup_database(target_db, backup)
            receipt["backup"] = {"path": str(backup), "size_bytes": backup.stat().st_size,
                                 "sha256": sha256_file(backup)}
        connection = sqlite3.connect(str(target_db))
        try:
            connection.execute("PRAGMA busy_timeout=60000")
            connection.execute("PRAGMA foreign_keys=ON")
            ensure_schema(connection)
            target_columns = set(table_columns(connection, "TargetCompounds"))
            missing = [c for c in payload["columns"] if c not in target_columns]
            if missing:
                raise RuntimeError(f"target TargetCompounds lacks required columns: {missing}")
            connection.execute("BEGIN IMMEDIATE")
            retired = {}
            for ds in retire_datasets:
                removed, _ = delete_dataset(connection, ds)
                retired[ds] = removed
            receipt["retired_rows"] = retired
            delete_dataset(connection, dataset)  # idempotent replace of same dataset
            columns = payload["columns"]
            marks = ",".join("?" for _ in columns)
            quoted = ",".join(f'"{c}"' for c in columns)
            connection.executemany(
                f"INSERT INTO TargetCompounds({quoted}) VALUES ({marks})",
                [tuple(row[c] for c in columns) for row in payload["targets"]])
            target_id_to_pk = {str(r[1]): int(r[0]) for r in connection.execute(
                "SELECT id,target_id FROM TargetCompounds WHERE dataset_source=?", (dataset,))}
            connection.executemany(
                "INSERT INTO DiagnosticIons(target_id,fragment_mz,ordinal) VALUES (?,?,?)",
                [(target_id_to_pk[r["target_id"]], r["fragment_mz"], r["ordinal"]) for r in payload["diagnostics"]])
            installed = int(connection.execute(
                "SELECT COUNT(*) FROM TargetCompounds WHERE dataset_source=?", (dataset,)).fetchone()[0])
            installed_ions = int(connection.execute(
                """SELECT COUNT(*) FROM DiagnosticIons d JOIN TargetCompounds t ON t.id=d.target_id
                   WHERE t.dataset_source=?""", (dataset,)).fetchone()[0])
            if installed != len(payload["targets"]) or installed_ions != len(payload["diagnostics"]):
                raise RuntimeError(f"post-insert counts mismatch: {installed}/{installed_ions}")
            fk = connection.execute("PRAGMA foreign_key_check").fetchall()
            if fk:
                raise RuntimeError(f"foreign key violations: {fk[:5]}")
            connection.execute("COMMIT")
            receipt["status"] = "applied"
            receipt["installed_targets"] = installed
            receipt["installed_diagnostic_ions"] = installed_ions
        except Exception:
            try:
                connection.execute("ROLLBACK")
            except sqlite3.Error:
                pass
            raise
        finally:
            connection.close()
    if receipt_path:
        receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=1))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-db", required=True, type=Path)
    parser.add_argument("--target-db", required=True, type=Path)
    parser.add_argument("--dataset-source", default=None,
                        help="dataset identity inside the staged DB (default: single dataset present)")
    parser.add_argument("--apply", action="store_true", help="actually write (default dry-run)")
    parser.add_argument("--allow-canonical-shared", action="store_true")
    parser.add_argument("--retire-dataset", action="append", default=[],
                        help="remove superseded dataset_source rows in the same transaction")
    parser.add_argument("--no-backup", action="store_true")
    parser.add_argument("--receipt", type=Path, default=None)
    args = parser.parse_args()
    receipt = migrate(
        args.source_db, args.target_db, dataset=args.dataset_source, apply=args.apply,
        allow_canonical_shared=args.allow_canonical_shared, create_backup=not args.no_backup,
        retire_datasets=args.retire_dataset, receipt_path=args.receipt)
    print(json.dumps(receipt, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
