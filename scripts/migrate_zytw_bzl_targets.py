#!/usr/bin/env python3
"""Safely replace the zytw + BZL Expected Compounds dataset in V10.

The migration is deliberately narrow: it may write only ``TargetCompounds``
and ``DiagnosticIons``.  A JSON compensation ledger captures the exact rows
being replaced, and ``--rollback-ledger`` restores those rows transactionally.
By default the command is a read-only preview; writes require both ``--apply``
and ``MSSYSTEM_LIBRARY_WRITE_ENABLED=1`` while ports 8509/8510 are stopped.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.library_import import (  # noqa: E402
    _append_target_records_with_receipt,
    _assert_library_write_enabled,
    _target_identifier,
)
from matching.library_quality import assess_import_records  # noqa: E402


DATASET_SOURCE = "user_contrib:literature_template_zytw_bzl"
DEFAULT_DB = Path("/Users/zizi/Desktop/Mssystem_shared_data/local_library.db")
DEFAULT_PAYLOAD = Path(
    "/Users/zizi/Desktop/MsQCModel_v2/_local_only/mgf_msp_database/"
    "待新增数据库数据/zytw_bzl_target_payload_v4.json"
)
DEFAULT_LEDGER_DIR = Path(
    "/Users/zizi/Desktop/Mssystem_shared_data/migration_ledgers"
)
EXPECTED_SOURCE_ROWS = 346
EXPECTED_UNIQUE_TARGETS = 345
EXPECTED_STRUCTURES = 205
EXPECTED_FORMULA_ROWS = 345
EXPECTED_PERSISTED_FORMULA_TARGETS = 344
EXPECTED_PLANTS = {
    "Ligularia cymbulifera": {"source_rows": 118, "targets": 117, "structures": 87},
    "Scutellaria barbata": {"source_rows": 228, "targets": 228, "structures": 118},
}
SERVICE_PORTS = (8509, 8510)
_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _write_json_atomic(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _quote_identifier(value: str) -> str:
    if not _IDENTIFIER.fullmatch(value):
        raise ValueError(f"unsafe SQLite identifier: {value!r}")
    return f'"{value}"'


def _connect(path: Path, *, read_only: bool) -> sqlite3.Connection:
    if read_only:
        uri = path.resolve().as_uri() + "?mode=ro"
        connection = sqlite3.connect(uri, uri=True)
    else:
        connection = sqlite3.connect(str(path))
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout=60000")
    connection.execute("PRAGMA foreign_keys=ON")
    return connection


def _user_tables(connection: sqlite3.Connection) -> list[str]:
    return [
        str(row[0])
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%' ORDER BY name"
        )
    ]


def _protected_tables(connection: sqlite3.Connection) -> list[str]:
    writable = {"TargetCompounds", "DiagnosticIons"}
    return [table for table in _user_tables(connection) if table not in writable]


def _protected_snapshot(
    connection: sqlite3.Connection, tables: Sequence[str]
) -> dict[str, dict[str, Any]]:
    snapshot: dict[str, dict[str, Any]] = {}
    for table in tables:
        quoted = _quote_identifier(table)
        schema_row = connection.execute(
            "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", (table,)
        ).fetchone()
        schema = str(schema_row[0] or "") if schema_row else ""
        snapshot[table] = {
            "row_count": int(connection.execute(f"SELECT COUNT(*) FROM {quoted}").fetchone()[0]),
            "schema_sha256": _sha256_bytes(schema.encode("utf-8")),
        }
    return snapshot


def _quick_check(connection: sqlite3.Connection) -> str:
    row = connection.execute("PRAGMA quick_check(1)").fetchone()
    return str(row[0]) if row else "missing_result"


def _foreign_key_check(connection: sqlite3.Connection) -> list[list[Any]]:
    return [list(row) for row in connection.execute("PRAGMA foreign_key_check")]


def _dataset_summary(connection: sqlite3.Connection) -> dict[str, Any]:
    row = connection.execute(
        "SELECT COUNT(*) AS n, "
        "SUM(CASE WHEN COALESCE(TRIM(smiles),'')<>'' "
        "AND COALESCE(structure_identity_verified,0)=1 THEN 1 ELSE 0 END) AS structures, "
        "SUM(CASE WHEN COALESCE(TRIM(molecular_formula),'')<>'' THEN 1 ELSE 0 END) AS formulas "
        "FROM TargetCompounds WHERE dataset_source=?",
        (DATASET_SOURCE,),
    ).fetchone()
    plants = {
        str(item[0] or ""): {
            "targets": int(item[1]),
            "structures": int(item[2] or 0),
            "formulas": int(item[3] or 0),
        }
        for item in connection.execute(
            "SELECT plant_source, COUNT(*), "
            "SUM(CASE WHEN COALESCE(TRIM(smiles),'')<>'' "
            "AND COALESCE(structure_identity_verified,0)=1 THEN 1 ELSE 0 END), "
            "SUM(CASE WHEN COALESCE(TRIM(molecular_formula),'')<>'' THEN 1 ELSE 0 END) "
            "FROM TargetCompounds WHERE dataset_source=? GROUP BY plant_source "
            "ORDER BY plant_source",
            (DATASET_SOURCE,),
        )
    }
    searchability = {
        str(item[0] or ""): int(item[1])
        for item in connection.execute(
            "SELECT target_searchability, COUNT(*) FROM TargetCompounds "
            "WHERE dataset_source=? GROUP BY target_searchability",
            (DATASET_SOURCE,),
        )
    }
    tiers = {
        str(item[0] or ""): int(item[1])
        for item in connection.execute(
            "SELECT source_tier, COUNT(*) FROM TargetCompounds "
            "WHERE dataset_source=? GROUP BY source_tier",
            (DATASET_SOURCE,),
        )
    }
    versions = sorted(
        {str(item[0] or "") for item in connection.execute(
            "SELECT DISTINCT dataset_version FROM TargetCompounds WHERE dataset_source=?",
            (DATASET_SOURCE,),
        )}
    )
    target_ids = sorted(
        str(item[0]) for item in connection.execute(
            "SELECT target_id FROM TargetCompounds WHERE dataset_source=?",
            (DATASET_SOURCE,),
        )
    )
    source_indices = sorted(
        int(item[0]) for item in connection.execute(
            "SELECT source_row_index FROM TargetCompounds WHERE dataset_source=? "
            "AND source_row_index IS NOT NULL",
            (DATASET_SOURCE,),
        )
    )
    diagnostic_ions = int(
        connection.execute(
            "SELECT COUNT(*) FROM DiagnosticIons d JOIN TargetCompounds t "
            "ON t.id=d.target_id WHERE t.dataset_source=?",
            (DATASET_SOURCE,),
        ).fetchone()[0]
    )
    return {
        "targets": int(row[0] or 0),
        "structures": int(row[1] or 0),
        "formulas": int(row[2] or 0),
        "diagnostic_ions": diagnostic_ions,
        "plants": plants,
        "target_searchability": searchability,
        "source_tiers": tiers,
        "dataset_versions": versions,
        "target_ids_sha256": _sha256_bytes("\n".join(target_ids).encode("utf-8")),
        "source_row_indices": source_indices,
    }


def _database_summary(connection: sqlite3.Connection) -> dict[str, Any]:
    return {
        "TargetCompounds": int(
            connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0]
        ),
        "DiagnosticIons": int(
            connection.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0]
        ),
        "dataset": _dataset_summary(connection),
    }


def _compact_quality(quality: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in quality.items()
        if key not in {"audits", "rejected_record_samples"}
    }


def _load_and_prepare_payload(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    raw_bytes = path.read_bytes()
    payload_sha256 = _sha256_bytes(raw_bytes)
    payload = json.loads(raw_bytes)
    if payload.get("schema_version") != 1:
        raise ValueError("expected payload schema_version=1")
    if payload.get("dataset_source") != DATASET_SOURCE:
        raise ValueError("payload dataset_source does not match migration scope")
    records = payload.get("records")
    if not isinstance(records, list) or len(records) != EXPECTED_SOURCE_ROWS:
        raise ValueError(f"expected {EXPECTED_SOURCE_ROWS} payload rows")
    if payload.get("expected_record_count") != EXPECTED_SOURCE_ROWS:
        raise ValueError("payload expected_record_count is inconsistent")

    indices = [int(record.get("source_row_index")) for record in records]
    if indices != list(range(EXPECTED_SOURCE_ROWS)):
        raise ValueError("payload source_row_index must be contiguous and ordered")
    datasets = Counter(str(record.get("dataset") or "") for record in records)
    if datasets != Counter({"zytw": 118, "bzl": 228}):
        raise ValueError(f"unexpected source dataset split: {dict(datasets)}")
    plant_source_rows = Counter(str(record.get("plant_source") or "") for record in records)
    expected_source_split = Counter(
        {plant: expected["source_rows"] for plant, expected in EXPECTED_PLANTS.items()}
    )
    if plant_source_rows != expected_source_split:
        raise ValueError(f"unexpected plant split: {dict(plant_source_rows)}")
    if sum(bool(record.get("smiles")) for record in records) != EXPECTED_STRUCTURES:
        raise ValueError("payload structure count is not the audited 205")
    if sum(bool(record.get("formula")) for record in records) != EXPECTED_FORMULA_ROWS:
        raise ValueError("payload formula count is not the audited 345")
    for record in records:
        if record.get("dataset_source") != DATASET_SOURCE:
            raise ValueError("record escaped the single-dataset migration scope")
        if record.get("source_tier") != "literature_presence_only":
            raise ValueError("all rows must remain literature_presence_only")
        if record.get("peaks"):
            raise ValueError("literature presence rows may not invent MS2 peaks")

    quality_input: list[dict[str, Any]] = []
    source_by_index: dict[int, dict[str, Any]] = {}
    fields = (
        "dataset_source", "source_file", "source_row_index", "name", "formula",
        "smiles", "inchi", "exact_mass", "precursor_mz", "adduct", "charge",
        "plant_source", "source_tier", "intensity_origin", "peaks", "instrument",
        "provenance", "dataset_accession",
    )
    for source in records:
        source_index = int(source["source_row_index"])
        source_by_index[source_index] = source
        record = {field: source.get(field) for field in fields}
        record.update({
            "reported_inchikey": source.get("inchikey"),
            "source_content_sha256": payload_sha256,
            "doi": None,
            "spectrum_accession": None,
        })
        quality_input.append(record)

    sanitized, quality = assess_import_records(quality_input)
    if len(sanitized) != EXPECTED_SOURCE_ROWS:
        raise ValueError(f"V10 quality gate accepted {len(sanitized)}, expected 346")
    if quality.get("rejected") != 0:
        raise ValueError(f"V10 quality gate rejected rows: {quality.get('error_counts')}")
    record_types = quality.get("record_types") or {}
    if record_types != {
        "expected_target_input": 346,
        "experimental_spectrum_input": 0,
        "expected_targets_accepted": 346,
        "experimental_spectra_accepted": 0,
    }:
        raise ValueError(f"unexpected record routing: {record_types}")
    if quality.get("target_searchability_counts") != {
        "theoretical_ion_ready": 345,
        "identity_only_unsearchable": 1,
    }:
        raise ValueError("unexpected source-row searchability counts")

    # Keep the more specific resolver provenance after the V10 gate has
    # independently re-derived and verified every structural identity.
    for record in sanitized:
        source = source_by_index[int(record["source_row_index"])]
        if record.get("structure_identity_verified"):
            record["structure_identity_source"] = source.get(
                "structure_identity_source"
            ) or "smiles_derived"
        else:
            record["structure_identity_source"] = "none"

    structures = sum(
        bool(record.get("smiles")) and bool(record.get("structure_identity_verified"))
        for record in sanitized
    )
    formulas = sum(bool(record.get("formula")) for record in sanitized)
    if structures != EXPECTED_STRUCTURES or formulas != EXPECTED_FORMULA_ROWS:
        raise ValueError(
            f"quality gate changed audited coverage: structures={structures}, formulas={formulas}"
        )

    ids: dict[str, list[int]] = defaultdict(list)
    for record in sanitized:
        ids[_target_identifier(record)].append(int(record["source_row_index"]))
    duplicates = {key: value for key, value in ids.items() if len(value) > 1}
    if len(ids) != EXPECTED_UNIQUE_TARGETS or sorted(duplicates.values()) != [[18, 37]]:
        raise ValueError(
            f"unexpected identity collapse: unique={len(ids)}, duplicates={duplicates}"
        )
    expected_id_hash = _sha256_bytes("\n".join(sorted(ids)).encode("utf-8"))
    metadata = {
        "payload_path": str(path.resolve()),
        "payload_sha256": payload_sha256,
        "source_rows": len(records),
        "unique_targets": len(ids),
        "structures": structures,
        "formula_rows": formulas,
        "duplicate_source_rows": sorted(duplicates.values()),
        "expected_target_ids_sha256": expected_id_hash,
        "quality": _compact_quality(quality),
    }
    return sanitized, metadata


def _open_ports(ports: Iterable[int]) -> list[int]:
    active: list[int] = []
    for port in ports:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as client:
            client.settimeout(0.25)
            if client.connect_ex(("127.0.0.1", int(port))) == 0:
                active.append(int(port))
    return active


def _require_write_window() -> None:
    _assert_library_write_enabled()
    active = _open_ports(SERVICE_PORTS)
    if active:
        raise RuntimeError(
            "refusing database write while MSsystem service ports are active: "
            + ", ".join(map(str, active))
        )


def _dump_dataset(connection: sqlite3.Connection) -> dict[str, Any]:
    target_columns = [str(row[1]) for row in connection.execute(
        "PRAGMA table_info(TargetCompounds)"
    )]
    target_rows = [
        [row[column] for column in target_columns]
        for row in connection.execute(
            "SELECT * FROM TargetCompounds WHERE dataset_source=? ORDER BY id",
            (DATASET_SOURCE,),
        )
    ]
    target_ids = [int(row[target_columns.index("id")]) for row in target_rows]
    diagnostic_columns = [str(row[1]) for row in connection.execute(
        "PRAGMA table_info(DiagnosticIons)"
    )]
    diagnostic_rows: list[list[Any]] = []
    if target_ids:
        for offset in range(0, len(target_ids), 500):
            chunk = target_ids[offset:offset + 500]
            marks = ",".join("?" for _ in chunk)
            diagnostic_rows.extend(
                [row[column] for column in diagnostic_columns]
                for row in connection.execute(
                    f"SELECT * FROM DiagnosticIons WHERE target_id IN ({marks}) ORDER BY id",
                    chunk,
                )
            )
    dump = {
        "TargetCompounds": {"columns": target_columns, "rows": target_rows},
        "DiagnosticIons": {"columns": diagnostic_columns, "rows": diagnostic_rows},
    }
    dump["rows_sha256"] = _sha256_bytes(_json_bytes(dump))
    return dump


def _default_output_path(prefix: str, directory: Path) -> Path:
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S%fZ")
    return directory / f"{prefix}_{timestamp}.json"


def _assert_expected_after(
    after: dict[str, Any], before: dict[str, Any], metadata: dict[str, Any]
) -> None:
    dataset = after["dataset"]
    if dataset["targets"] != EXPECTED_UNIQUE_TARGETS:
        raise RuntimeError(f"post-migration target count mismatch: {dataset['targets']}")
    if dataset["structures"] != EXPECTED_STRUCTURES:
        raise RuntimeError(f"post-migration structure count mismatch: {dataset['structures']}")
    if dataset["formulas"] != EXPECTED_PERSISTED_FORMULA_TARGETS:
        raise RuntimeError(f"post-migration formula count mismatch: {dataset['formulas']}")
    if dataset["diagnostic_ions"] != 0:
        raise RuntimeError("presence-only dataset unexpectedly contains DiagnosticIons")
    if dataset["target_searchability"] != {
        "identity_only_unsearchable": 1,
        "theoretical_ion_ready": EXPECTED_PERSISTED_FORMULA_TARGETS,
    }:
        raise RuntimeError(
            f"post-migration searchability mismatch: {dataset['target_searchability']}"
        )
    if dataset["source_tiers"] != {"literature_presence_only": 345}:
        raise RuntimeError(f"post-migration source tiers mismatch: {dataset['source_tiers']}")
    if dataset["dataset_versions"] != [metadata["payload_sha256"]]:
        raise RuntimeError(f"post-migration version mismatch: {dataset['dataset_versions']}")
    if dataset["target_ids_sha256"] != metadata["expected_target_ids_sha256"]:
        raise RuntimeError("post-migration target identity set mismatch")
    for plant, expected in EXPECTED_PLANTS.items():
        observed = dataset["plants"].get(plant) or {}
        if observed.get("targets") != expected["targets"]:
            raise RuntimeError(f"{plant} target count mismatch: {observed}")
        if observed.get("structures") != expected["structures"]:
            raise RuntimeError(f"{plant} structure count mismatch: {observed}")
    if after["TargetCompounds"] != before["TargetCompounds"]:
        raise RuntimeError("global TargetCompounds count changed during replacement")
    if after["DiagnosticIons"] != before["DiagnosticIons"]:
        raise RuntimeError("global DiagnosticIons count changed during replacement")


def _apply_migration(
    records: list[dict[str, Any]],
    metadata: dict[str, Any],
    db_path: Path,
    ledger_dir: Path,
    receipt_path: Path | None,
) -> dict[str, Any]:
    _require_write_window()
    with closing(_connect(db_path, read_only=True)) as connection:
        before = _database_summary(connection)
        protected_tables = _protected_tables(connection)
        protected_before = _protected_snapshot(connection, protected_tables)
        quick_before = _quick_check(connection)
        foreign_keys_before = _foreign_key_check(connection)
        already_current = (
            before["dataset"]["target_ids_sha256"]
            == metadata["expected_target_ids_sha256"]
            and before["dataset"]["dataset_versions"] == [metadata["payload_sha256"]]
        )
        ledger_dump = None if already_current else _dump_dataset(connection)
        db_stat = db_path.stat()
    if quick_before != "ok" or foreign_keys_before:
        raise RuntimeError(
            f"database failed preflight: quick_check={quick_before}, "
            f"foreign_key_violations={foreign_keys_before[:5]}"
        )

    ledger_path: Path | None = None
    ledger_sha256: str | None = None
    if ledger_dump is not None:
        ledger_dir.mkdir(parents=True, exist_ok=True)
        ledger_path = _default_output_path(
            f"zytw_bzl_before_{metadata['payload_sha256'][:12]}", ledger_dir
        )
        ledger = {
            "schema_version": 1,
            "kind": "zytw_bzl_target_compensation_ledger",
            "created_at": _utc_now(),
            "db_path": str(db_path.resolve()),
            "db_file_size": db_stat.st_size,
            "db_file_mtime_ns": db_stat.st_mtime_ns,
            "dataset_source": DATASET_SOURCE,
            "expected_post_payload_sha256": metadata["payload_sha256"],
            "expected_post_target_ids_sha256": metadata[
                "expected_target_ids_sha256"
            ],
            "before": before,
            "protected_tables": protected_tables,
            "protected_snapshot": protected_before,
            "quick_check": quick_before,
            "foreign_key_check": foreign_keys_before,
            "dataset_dump": ledger_dump,
        }
        _write_json_atomic(ledger_path, ledger)
        ledger_sha256 = _sha256_bytes(ledger_path.read_bytes())

    library_receipt = _append_target_records_with_receipt(
        records,
        db_path,
        replace_datasets=True,
        protected_tables=protected_tables,
    )
    with closing(_connect(db_path, read_only=True)) as connection:
        after = _database_summary(connection)
        protected_after = _protected_snapshot(connection, protected_tables)
        quick_after = _quick_check(connection)
        foreign_keys_after = _foreign_key_check(connection)
    _assert_expected_after(after, before, metadata)
    if protected_after != protected_before:
        raise RuntimeError("a protected non-target table changed during migration")
    if quick_after != "ok" or foreign_keys_after:
        raise RuntimeError(
            f"database failed postflight: quick_check={quick_after}, "
            f"foreign_key_violations={foreign_keys_after[:5]}"
        )
    expected_no_op = already_current
    observed_no_op = (
        int(library_receipt.get("new_targets") or 0) == 0
        and not library_receipt.get("replaced_datasets")
    )
    if expected_no_op != observed_no_op:
        raise RuntimeError(
            f"idempotency result mismatch: expected_no_op={expected_no_op}, "
            f"receipt={library_receipt}"
        )
    result = {
        "schema_version": 1,
        "kind": "zytw_bzl_target_migration_receipt",
        "status": "no_op_already_current" if observed_no_op else "applied",
        "completed_at": _utc_now(),
        "db_path": str(db_path.resolve()),
        "dataset_source": DATASET_SOURCE,
        "payload": metadata,
        "before": before,
        "after": after,
        "library_receipt": library_receipt,
        "compensation_ledger": str(ledger_path) if ledger_path else None,
        "compensation_ledger_sha256": ledger_sha256,
        "protected_tables": protected_tables,
        "protected_snapshot_verified": protected_after == protected_before,
        "quick_check_before": quick_before,
        "quick_check_after": quick_after,
        "foreign_key_check_before": foreign_keys_before,
        "foreign_key_check_after": foreign_keys_after,
        "ports_verified_stopped": list(SERVICE_PORTS),
    }
    output = receipt_path or _default_output_path("zytw_bzl_migration_receipt", ledger_dir)
    result["receipt_path"] = str(output.resolve())
    _write_json_atomic(output, result)
    return result


def _install_authorizer(
    connection: sqlite3.Connection, protected_tables: Sequence[str]
) -> None:
    protected = {table.casefold() for table in protected_tables}
    protected_actions = {
        getattr(sqlite3, name)
        for name in (
            "SQLITE_INSERT", "SQLITE_UPDATE", "SQLITE_DELETE", "SQLITE_CREATE_INDEX",
            "SQLITE_CREATE_TABLE", "SQLITE_CREATE_TRIGGER", "SQLITE_CREATE_VIEW",
            "SQLITE_DROP_INDEX", "SQLITE_DROP_TABLE", "SQLITE_DROP_TRIGGER",
            "SQLITE_DROP_VIEW", "SQLITE_ALTER_TABLE", "SQLITE_REINDEX",
        )
        if hasattr(sqlite3, name)
    }

    def authorize(
        action: int,
        argument1: str | None,
        argument2: str | None,
        _database: str | None,
        _trigger: str | None,
    ) -> int:
        referenced = {
            str(argument1 or "").strip().casefold(),
            str(argument2 or "").strip().casefold(),
        }
        if action in protected_actions and protected.intersection(referenced):
            return sqlite3.SQLITE_DENY
        return sqlite3.SQLITE_OK

    connection.set_authorizer(authorize)


def _restore_rows(
    db_path: Path, ledger_path: Path, receipt_path: Path | None
) -> dict[str, Any]:
    _require_write_window()
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    if ledger.get("kind") != "zytw_bzl_target_compensation_ledger":
        raise ValueError("not a zytw/BZL compensation ledger")
    if ledger.get("dataset_source") != DATASET_SOURCE:
        raise ValueError("ledger dataset is outside this rollback scope")
    if Path(ledger.get("db_path", "")).resolve() != db_path.resolve():
        raise ValueError("ledger was created from a different database")
    payload_sha256 = str(ledger.get("expected_post_payload_sha256") or "")
    with closing(_connect(db_path, read_only=True)) as connection:
        before = _database_summary(connection)
        protected_tables = list(ledger.get("protected_tables") or [])
        protected_before = _protected_snapshot(connection, protected_tables)
        quick_before = _quick_check(connection)
        foreign_keys_before = _foreign_key_check(connection)
    if before["dataset"]["dataset_versions"] != [payload_sha256]:
        raise RuntimeError(
            "refusing rollback because the current dataset is not the exact "
            "version written by this ledger"
        )
    if before["dataset"]["target_ids_sha256"] != str(
        ledger.get("expected_post_target_ids_sha256") or ""
    ):
        raise RuntimeError(
            "refusing rollback because the current target identity set has changed"
        )
    if quick_before != "ok" or foreign_keys_before:
        raise RuntimeError("database failed rollback preflight")

    dump = ledger.get("dataset_dump") or {}
    target = dump.get("TargetCompounds") or {}
    diagnostic = dump.get("DiagnosticIons") or {}
    expected_dump_hash = str(dump.get("rows_sha256") or "")
    comparable_dump = {
        "TargetCompounds": target,
        "DiagnosticIons": diagnostic,
    }
    if _sha256_bytes(_json_bytes(comparable_dump)) != expected_dump_hash:
        raise ValueError("compensation ledger row checksum does not match")

    connection = _connect(db_path, read_only=False)
    committed = False
    try:
        _install_authorizer(connection, protected_tables)
        connection.execute("BEGIN IMMEDIATE")
        connection.execute(
            "DELETE FROM DiagnosticIons WHERE target_id IN "
            "(SELECT id FROM TargetCompounds WHERE dataset_source=?)",
            (DATASET_SOURCE,),
        )
        connection.execute(
            "DELETE FROM TargetCompounds WHERE dataset_source=?", (DATASET_SOURCE,)
        )
        target_columns = [str(value) for value in target.get("columns") or []]
        target_rows = list(target.get("rows") or [])
        if target_columns and target_rows:
            columns = ",".join(_quote_identifier(value) for value in target_columns)
            marks = ",".join("?" for _ in target_columns)
            connection.executemany(
                f"INSERT INTO TargetCompounds ({columns}) VALUES ({marks})", target_rows
            )
        diagnostic_columns = [str(value) for value in diagnostic.get("columns") or []]
        diagnostic_rows = list(diagnostic.get("rows") or [])
        if diagnostic_columns and diagnostic_rows:
            columns = ",".join(_quote_identifier(value) for value in diagnostic_columns)
            marks = ",".join("?" for _ in diagnostic_columns)
            connection.executemany(
                f"INSERT INTO DiagnosticIons ({columns}) VALUES ({marks})",
                diagnostic_rows,
            )
        connection.execute("COMMIT")
        committed = True
    finally:
        if not committed:
            try:
                connection.execute("ROLLBACK")
            except sqlite3.Error:
                pass
        connection.close()

    with closing(_connect(db_path, read_only=True)) as connection:
        after = _database_summary(connection)
        protected_after = _protected_snapshot(connection, protected_tables)
        quick_after = _quick_check(connection)
        foreign_keys_after = _foreign_key_check(connection)
    if after != ledger.get("before"):
        raise RuntimeError("rollback rows do not exactly reproduce the pre-migration summary")
    if protected_after != protected_before:
        raise RuntimeError("a protected table changed during rollback")
    if quick_after != "ok" or foreign_keys_after:
        raise RuntimeError("database failed rollback postflight")
    result = {
        "schema_version": 1,
        "kind": "zytw_bzl_target_rollback_receipt",
        "status": "rolled_back",
        "completed_at": _utc_now(),
        "db_path": str(db_path.resolve()),
        "dataset_source": DATASET_SOURCE,
        "ledger_path": str(ledger_path.resolve()),
        "ledger_sha256": _sha256_bytes(ledger_path.read_bytes()),
        "before": before,
        "after": after,
        "protected_snapshot_verified": protected_after == protected_before,
        "quick_check_before": quick_before,
        "quick_check_after": quick_after,
        "foreign_key_check_before": foreign_keys_before,
        "foreign_key_check_after": foreign_keys_after,
        "ports_verified_stopped": list(SERVICE_PORTS),
    }
    output = receipt_path or ledger_path.with_name(
        ledger_path.stem + "_rollback_receipt.json"
    )
    result["receipt_path"] = str(output.resolve())
    _write_json_atomic(output, result)
    return result


def _preview(records: list[dict[str, Any]], metadata: dict[str, Any], db_path: Path) -> dict[str, Any]:
    with closing(_connect(db_path, read_only=True)) as connection:
        current = _database_summary(connection)
        protected_tables = _protected_tables(connection)
    return {
        "kind": "zytw_bzl_target_migration_preview",
        "status": "read_only",
        "db_path": str(db_path.resolve()),
        "dataset_source": DATASET_SOURCE,
        "payload": metadata,
        "prepared_records": len(records),
        "current": current,
        "protected_tables": protected_tables,
        "open_service_ports": _open_ports(SERVICE_PORTS),
        "write_requested": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, default=DEFAULT_PAYLOAD)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--ledger-dir", type=Path, default=DEFAULT_LEDGER_DIR)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--rollback-ledger", type=Path)
    args = parser.parse_args()
    if args.apply and args.rollback_ledger:
        parser.error("--apply and --rollback-ledger are mutually exclusive")
    if not args.db.is_file():
        parser.error(f"database does not exist: {args.db}")

    if args.rollback_ledger:
        result = _restore_rows(args.db, args.rollback_ledger, args.receipt)
    else:
        if not args.payload.is_file():
            parser.error(f"payload does not exist: {args.payload}")
        records, metadata = _load_and_prepare_payload(args.payload)
        if args.apply:
            result = _apply_migration(
                records, metadata, args.db, args.ledger_dir, args.receipt
            )
        else:
            result = _preview(records, metadata, args.db)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
