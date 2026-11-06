#!/usr/bin/env python3
"""Atomically replace the legacy combined zytw/BZL target dataset.

The only writable tables are ``TargetCompounds`` and ``DiagnosticIons``.
Writes require ``--apply`` plus ``MSSYSTEM_LIBRARY_WRITE_ENABLED=1``.  For the
formal shared database, ports 8509 and 8510 must both be stopped.  Before the
transaction, an exact row-level compensation ledger is written; the same
command can restore it with ``--rollback-ledger``.
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from contextlib import closing
from pathlib import Path
from typing import Any, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.library_import import (  # noqa: E402
    _assert_library_write_enabled,
    _diagnostic_ions,
    _reported_full_inchikey,
    _target_identifier,
    _verified_inchikey,
)
from matching.library_quality import assess_import_records  # noqa: E402
from scripts.migrate_zytw_bzl_targets import (  # noqa: E402
    _connect,
    _default_output_path,
    _foreign_key_check,
    _install_authorizer,
    _json_bytes,
    _open_ports,
    _protected_snapshot,
    _protected_tables,
    _quick_check,
    _quote_identifier,
    _sha256_bytes,
    _utc_now,
    _write_json_atomic,
)


OLD_DATASET_SOURCE = "user_contrib:literature_template_zytw_bzl"
OLD_PAYLOAD_SHA256 = "3bba2a12f946731bea7cd893f96360e51fd4ec6edaa240702f7ad4eb4092f2aa"
DEFAULT_DB = Path("/Users/zizi/Desktop/Mssystem_shared_data/local_library.db")
DEFAULT_DATA_DIR = Path(
    "/Users/zizi/Desktop/MsQCModel_v2/_local_only/mgf_msp_database/"
    "待新增数据库数据"
)
DEFAULT_LEDGER_DIR = Path(
    "/Users/zizi/Desktop/Mssystem_shared_data/migration_ledgers"
)
SERVICE_PORTS = (8509, 8510)
SPECS: dict[str, dict[str, Any]] = {
    "zytw": {
        "dataset_source": "user_contrib:zytw",
        "payload": DEFAULT_DATA_DIR / "zytw_target_payload_v5.json",
        "source_rows": 118,
        "unique_targets": 117,
        "structures": 87,
        "formula_rows": 118,
        "persisted_formulas": 117,
        "identity_only": 0,
        "plant_source": "Ligularia cymbulifera",
        "duplicate_source_rows": [[18, 37]],
    },
    "bzl": {
        "dataset_source": "user_contrib:bzl",
        "payload": DEFAULT_DATA_DIR / "bzl_target_payload_v3.json",
        "source_rows": 228,
        "unique_targets": 228,
        "structures": 118,
        "formula_rows": 227,
        "persisted_formulas": 227,
        "identity_only": 1,
        "plant_source": "Scutellaria barbata",
        "duplicate_source_rows": [],
    },
}
SCOPED_SOURCES = (
    OLD_DATASET_SOURCE,
    SPECS["zytw"]["dataset_source"],
    SPECS["bzl"]["dataset_source"],
)


def _compact_quality(quality: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in quality.items()
        if key not in {"audits", "rejected_record_samples"}
    }


def _load_payload(
    dataset: str, path: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    spec = SPECS[dataset]
    raw_bytes = path.read_bytes()
    payload_sha256 = _sha256_bytes(raw_bytes)
    payload = json.loads(raw_bytes)
    records = payload.get("records")
    if payload.get("schema_version") != 1:
        raise ValueError(f"{dataset}: expected payload schema_version=1")
    if payload.get("dataset_source") != spec["dataset_source"]:
        raise ValueError(f"{dataset}: payload dataset_source mismatch")
    if payload.get("expected_record_count") != spec["source_rows"]:
        raise ValueError(f"{dataset}: expected_record_count mismatch")
    if not isinstance(records, list) or len(records) != spec["source_rows"]:
        raise ValueError(f"{dataset}: source row count mismatch")
    indices = [int(record.get("source_row_index")) for record in records]
    if indices != list(range(spec["source_rows"])):
        raise ValueError(f"{dataset}: source_row_index is not local and contiguous")
    if Counter(str(record.get("dataset") or "") for record in records) != Counter(
        {dataset: spec["source_rows"]}
    ):
        raise ValueError(f"{dataset}: payload contains the other contributor")
    if Counter(str(record.get("plant_source") or "") for record in records) != Counter(
        {spec["plant_source"]: spec["source_rows"]}
    ):
        raise ValueError(f"{dataset}: plant_source split mismatch")
    if sum(bool(record.get("smiles")) for record in records) != spec["structures"]:
        raise ValueError(f"{dataset}: structure coverage changed")
    if sum(bool(record.get("formula")) for record in records) != spec["formula_rows"]:
        raise ValueError(f"{dataset}: formula coverage changed")
    for record in records:
        if record.get("dataset_source") != spec["dataset_source"]:
            raise ValueError(f"{dataset}: record escaped its dataset scope")
        if record.get("source_tier") != "literature_presence_only":
            raise ValueError(f"{dataset}: non-presence source tier found")
        if record.get("peaks"):
            raise ValueError(f"{dataset}: literature target may not invent MS2 peaks")

    fields = (
        "dataset_source", "source_file", "source_row_index", "name", "formula",
        "smiles", "inchi", "exact_mass", "precursor_mz", "adduct", "charge",
        "plant_source", "source_tier", "intensity_origin", "peaks", "instrument",
        "provenance", "dataset_accession",
    )
    source_by_index: dict[int, dict[str, Any]] = {}
    quality_input: list[dict[str, Any]] = []
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
    if len(sanitized) != spec["source_rows"] or quality.get("rejected") != 0:
        raise ValueError(
            f"{dataset}: V10 quality gate changed row acceptance: "
            f"accepted={len(sanitized)}, errors={quality.get('error_counts')}"
        )
    expected_routing = {
        "expected_target_input": spec["source_rows"],
        "experimental_spectrum_input": 0,
        "expected_targets_accepted": spec["source_rows"],
        "experimental_spectra_accepted": 0,
    }
    if quality.get("record_types") != expected_routing:
        raise ValueError(f"{dataset}: unexpected record routing")
    expected_searchability = {"theoretical_ion_ready": spec["formula_rows"]}
    if spec["identity_only"]:
        expected_searchability["identity_only_unsearchable"] = spec["identity_only"]
    if quality.get("target_searchability_counts") != expected_searchability:
        raise ValueError(
            f"{dataset}: unexpected source-row searchability "
            f"{quality.get('target_searchability_counts')}"
        )

    for record in sanitized:
        source = source_by_index[int(record["source_row_index"])]
        record["structure_identity_source"] = (
            source.get("structure_identity_source")
            if record.get("structure_identity_verified")
            else "none"
        )

    grouped_ids: dict[str, list[int]] = defaultdict(list)
    persisted_by_id: dict[str, dict[str, Any]] = {}
    for record in sanitized:
        target_id = _target_identifier(record)
        grouped_ids[target_id].append(int(record["source_row_index"]))
        persisted_by_id.setdefault(target_id, record)
    duplicates = sorted(
        sorted(indices) for indices in grouped_ids.values() if len(indices) > 1
    )
    if len(persisted_by_id) != spec["unique_targets"]:
        raise ValueError(f"{dataset}: unique target count changed")
    if duplicates != spec["duplicate_source_rows"]:
        raise ValueError(f"{dataset}: duplicate collapse changed: {duplicates}")
    persisted = list(persisted_by_id.values())
    persisted_structures = sum(
        bool(record.get("smiles")) and bool(record.get("structure_identity_verified"))
        for record in persisted
    )
    persisted_formulas = sum(bool(record.get("formula")) for record in persisted)
    if persisted_structures != spec["structures"]:
        raise ValueError(f"{dataset}: persisted structure coverage changed")
    if persisted_formulas != spec["persisted_formulas"]:
        raise ValueError(f"{dataset}: persisted formula coverage changed")
    if any(_diagnostic_ions(record) for record in persisted):
        raise ValueError(f"{dataset}: diagnostic ions must remain empty")

    target_ids = sorted(persisted_by_id)
    metadata = {
        "dataset": dataset,
        "dataset_source": spec["dataset_source"],
        "payload_path": str(path.resolve()),
        "payload_sha256": payload_sha256,
        "source_rows": spec["source_rows"],
        "unique_targets": spec["unique_targets"],
        "structures": persisted_structures,
        "formula_rows": spec["formula_rows"],
        "persisted_formulas": persisted_formulas,
        "identity_only": spec["identity_only"],
        "duplicate_source_rows": duplicates,
        "expected_target_ids_sha256": _sha256_bytes(
            "\n".join(target_ids).encode("utf-8")
        ),
        "quality": _compact_quality(quality),
    }
    return persisted, metadata


def _source_summary(connection: sqlite3.Connection, source: str) -> dict[str, Any]:
    row = connection.execute(
        "SELECT COUNT(*), "
        "SUM(CASE WHEN COALESCE(TRIM(smiles),'')<>'' "
        "AND COALESCE(structure_identity_verified,0)=1 THEN 1 ELSE 0 END), "
        "SUM(CASE WHEN COALESCE(TRIM(molecular_formula),'')<>'' THEN 1 ELSE 0 END) "
        "FROM TargetCompounds WHERE dataset_source=?",
        (source,),
    ).fetchone()
    target_ids = sorted(
        str(item[0]) for item in connection.execute(
            "SELECT target_id FROM TargetCompounds WHERE dataset_source=?", (source,)
        )
    )
    plants = {
        str(item[0] or ""): int(item[1])
        for item in connection.execute(
            "SELECT plant_source, COUNT(*) FROM TargetCompounds "
            "WHERE dataset_source=? GROUP BY plant_source ORDER BY plant_source",
            (source,),
        )
    }
    searchability = {
        str(item[0] or ""): int(item[1])
        for item in connection.execute(
            "SELECT target_searchability, COUNT(*) FROM TargetCompounds "
            "WHERE dataset_source=? GROUP BY target_searchability",
            (source,),
        )
    }
    tiers = {
        str(item[0] or ""): int(item[1])
        for item in connection.execute(
            "SELECT source_tier, COUNT(*) FROM TargetCompounds "
            "WHERE dataset_source=? GROUP BY source_tier",
            (source,),
        )
    }
    versions = sorted({
        str(item[0] or "") for item in connection.execute(
            "SELECT DISTINCT dataset_version FROM TargetCompounds WHERE dataset_source=?",
            (source,),
        )
    })
    diagnostic_ions = int(connection.execute(
        "SELECT COUNT(*) FROM DiagnosticIons d JOIN TargetCompounds t "
        "ON t.id=d.target_id WHERE t.dataset_source=?", (source,)
    ).fetchone()[0])
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
    }


def _database_summary(connection: sqlite3.Connection) -> dict[str, Any]:
    return {
        "TargetCompounds": int(
            connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0]
        ),
        "DiagnosticIons": int(
            connection.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0]
        ),
        "sources": {
            source: _source_summary(connection, source) for source in SCOPED_SOURCES
        },
    }


def _is_current(summary: dict[str, Any], metadata: dict[str, dict[str, Any]]) -> bool:
    if summary["sources"][OLD_DATASET_SOURCE]["targets"] != 0:
        return False
    for dataset, expected in metadata.items():
        observed = summary["sources"][expected["dataset_source"]]
        if observed["targets"] != expected["unique_targets"]:
            return False
        if observed["structures"] != expected["structures"]:
            return False
        if observed["formulas"] != expected["persisted_formulas"]:
            return False
        if observed["diagnostic_ions"] != 0:
            return False
        if observed["dataset_versions"] != [expected["payload_sha256"]]:
            return False
        if observed["target_ids_sha256"] != expected["expected_target_ids_sha256"]:
            return False
    return True


def _assert_replaceable_before(summary: dict[str, Any]) -> None:
    old = summary["sources"][OLD_DATASET_SOURCE]
    if (
        old["targets"] != 345
        or old["structures"] != 205
        or old["formulas"] != 344
        or old["diagnostic_ions"] != 0
        or old["dataset_versions"] != [OLD_PAYLOAD_SHA256]
        or old["source_tiers"] != {"literature_presence_only": 345}
    ):
        raise RuntimeError(f"legacy combined dataset is not the audited version: {old}")
    for dataset in ("zytw", "bzl"):
        source = SPECS[dataset]["dataset_source"]
        if summary["sources"][source]["targets"] != 0:
            raise RuntimeError(f"refusing to overwrite unexpected pre-existing {source}")


def _assert_expected_after(
    after: dict[str, Any], before: dict[str, Any], metadata: dict[str, dict[str, Any]]
) -> None:
    if after["sources"][OLD_DATASET_SOURCE]["targets"] != 0:
        raise RuntimeError("legacy combined dataset still exists after split")
    for dataset, expected in metadata.items():
        observed = after["sources"][expected["dataset_source"]]
        required_searchability = {"theoretical_ion_ready": expected["persisted_formulas"]}
        if expected["identity_only"]:
            required_searchability["identity_only_unsearchable"] = expected["identity_only"]
        required = {
            "targets": expected["unique_targets"],
            "structures": expected["structures"],
            "formulas": expected["persisted_formulas"],
            "diagnostic_ions": 0,
            "plants": {SPECS[dataset]["plant_source"]: expected["unique_targets"]},
            "target_searchability": required_searchability,
            "source_tiers": {"literature_presence_only": expected["unique_targets"]},
            "dataset_versions": [expected["payload_sha256"]],
            "target_ids_sha256": expected["expected_target_ids_sha256"],
        }
        if observed != required:
            raise RuntimeError(
                f"post-split summary mismatch for {dataset}: {observed} != {required}"
            )
    if after["TargetCompounds"] != before["TargetCompounds"]:
        raise RuntimeError("global TargetCompounds count changed during split")
    if after["DiagnosticIons"] != before["DiagnosticIons"]:
        raise RuntimeError("global DiagnosticIons count changed during split")


def _dump_scoped(connection: sqlite3.Connection) -> dict[str, Any]:
    marks = ",".join("?" for _ in SCOPED_SOURCES)
    target_columns = [
        str(row[1]) for row in connection.execute("PRAGMA table_info(TargetCompounds)")
    ]
    target_rows = [
        [row[column] for column in target_columns]
        for row in connection.execute(
            f"SELECT * FROM TargetCompounds WHERE dataset_source IN ({marks}) ORDER BY id",
            SCOPED_SOURCES,
        )
    ]
    id_index = target_columns.index("id")
    target_ids = [int(row[id_index]) for row in target_rows]
    diagnostic_columns = [
        str(row[1]) for row in connection.execute("PRAGMA table_info(DiagnosticIons)")
    ]
    diagnostic_rows: list[list[Any]] = []
    for offset in range(0, len(target_ids), 500):
        chunk = target_ids[offset:offset + 500]
        if not chunk:
            continue
        chunk_marks = ",".join("?" for _ in chunk)
        diagnostic_rows.extend(
            [row[column] for column in diagnostic_columns]
            for row in connection.execute(
                f"SELECT * FROM DiagnosticIons WHERE target_id IN ({chunk_marks}) ORDER BY id",
                chunk,
            )
        )
    dump = {
        "TargetCompounds": {"columns": target_columns, "rows": target_rows},
        "DiagnosticIons": {"columns": diagnostic_columns, "rows": diagnostic_rows},
    }
    dump["rows_sha256"] = _sha256_bytes(_json_bytes(dump))
    return dump


def _require_write_window(db_path: Path) -> list[int]:
    _assert_library_write_enabled()
    active = _open_ports(SERVICE_PORTS)
    if db_path.resolve() == DEFAULT_DB.resolve() and active:
        raise RuntimeError(
            "refusing formal database write while MSsystem ports are active: "
            + ", ".join(map(str, active))
        )
    return active


def _delete_scoped(connection: sqlite3.Connection) -> None:
    marks = ",".join("?" for _ in SCOPED_SOURCES)
    connection.execute(
        "DELETE FROM DiagnosticIons WHERE target_id IN "
        f"(SELECT id FROM TargetCompounds WHERE dataset_source IN ({marks}))",
        SCOPED_SOURCES,
    )
    connection.execute(
        f"DELETE FROM TargetCompounds WHERE dataset_source IN ({marks})",
        SCOPED_SOURCES,
    )


def _insert_records(
    connection: sqlite3.Connection, records: Sequence[dict[str, Any]]
) -> None:
    sql = (
        "INSERT INTO TargetCompounds (target_id, source_file, source_content_sha256, "
        "dataset_version, source_row_index, dataset_source, compound_name, "
        "molecular_formula, smiles, inchi, inchikey, reported_inchikey, "
        "structure_identity_verified, structure_identity_source, exact_mass, "
        "precursor_mz, adduct, charge, reported_precursor_mz, reported_adduct, "
        "reported_retention_time, reported_retention_time_unit, "
        "reported_collision_energy_raw, instrument, target_searchability, "
        "plant_source, provenance, doi, dataset_accession, spectrum_accession, "
        "source_tier, quality_flags_json, created_at) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    )
    created_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    inserted: set[str] = set()
    for record in records:
        target_id = _target_identifier(record)
        if target_id in inserted:
            raise RuntimeError(f"duplicate persisted target ID: {target_id}")
        inserted.add(target_id)
        verified_key = _verified_inchikey(record)
        reported_key = _reported_full_inchikey(record)
        connection.execute(sql, (
            target_id,
            record.get("source_file"),
            record.get("source_content_sha256"),
            record.get("source_content_sha256"),
            record.get("source_row_index"),
            record.get("dataset_source"),
            record.get("name"),
            record.get("formula"),
            record.get("smiles"),
            record.get("inchi"),
            verified_key,
            reported_key,
            1 if verified_key else 0,
            record.get("structure_identity_source") or (
                "smiles_derived" if verified_key else (
                    "reported_unverified" if reported_key else "none"
                )
            ),
            record.get("exact_mass"),
            record.get("precursor_mz"),
            record.get("adduct"),
            record.get("charge"),
            record.get("reported_precursor_mz"),
            record.get("reported_adduct"),
            record.get("reported_retention_time"),
            record.get("reported_retention_time_unit"),
            record.get("reported_collision_energy_raw"),
            record.get("instrument"),
            record.get("target_searchability"),
            record.get("plant_source"),
            record.get("provenance"),
            record.get("doi"),
            record.get("dataset_accession"),
            record.get("spectrum_accession"),
            "literature_presence_only",
            json.dumps(record.get("quality_flags") or [], ensure_ascii=False),
            created_at,
        ))
        ions = _diagnostic_ions(record)
        if ions:
            raise RuntimeError("split payload unexpectedly contains DiagnosticIons")


def _apply_split(
    records: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
    db_path: Path,
    ledger_dir: Path,
    receipt_path: Path | None,
) -> dict[str, Any]:
    active_ports = _require_write_window(db_path)
    with closing(_connect(db_path, read_only=True)) as connection:
        before = _database_summary(connection)
        protected = _protected_tables(connection)
        protected_before = _protected_snapshot(connection, protected)
        quick_before = _quick_check(connection)
        foreign_before = _foreign_key_check(connection)
        already_current = _is_current(before, metadata)
        dump = None if already_current else _dump_scoped(connection)
        db_stat = db_path.stat()
    if quick_before != "ok" or foreign_before:
        raise RuntimeError("database failed split preflight")
    if not already_current:
        _assert_replaceable_before(before)

    ledger_path: Path | None = None
    ledger_sha256: str | None = None
    if dump is not None:
        ledger_dir.mkdir(parents=True, exist_ok=True)
        combined_version = _sha256_bytes(
            "\n".join(sorted(item["payload_sha256"] for item in metadata.values()))
            .encode("utf-8")
        )
        ledger_path = _default_output_path(
            f"zytw_bzl_split_before_{combined_version[:12]}", ledger_dir
        )
        ledger = {
            "schema_version": 1,
            "kind": "zytw_bzl_split_compensation_ledger",
            "created_at": _utc_now(),
            "db_path": str(db_path.resolve()),
            "db_file_size": db_stat.st_size,
            "db_file_mtime_ns": db_stat.st_mtime_ns,
            "scoped_sources": list(SCOPED_SOURCES),
            "expected_post": metadata,
            "before": before,
            "protected_tables": protected,
            "protected_snapshot": protected_before,
            "quick_check": quick_before,
            "foreign_key_check": foreign_before,
            "dataset_dump": dump,
        }
        _write_json_atomic(ledger_path, ledger)
        ledger_sha256 = _sha256_bytes(ledger_path.read_bytes())

    if not already_current:
        connection = _connect(db_path, read_only=False)
        committed = False
        try:
            _install_authorizer(connection, protected)
            connection.execute("BEGIN IMMEDIATE")
            _delete_scoped(connection)
            _insert_records(connection, records)
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
        protected_after = _protected_snapshot(connection, protected)
        quick_after = _quick_check(connection)
        foreign_after = _foreign_key_check(connection)
    _assert_expected_after(after, before, metadata)
    if protected_after != protected_before:
        raise RuntimeError("a protected non-target table changed during split")
    if quick_after != "ok" or foreign_after:
        raise RuntimeError("database failed split postflight")

    result = {
        "schema_version": 1,
        "kind": "zytw_bzl_split_migration_receipt",
        "status": "no_op_already_current" if already_current else "applied",
        "completed_at": _utc_now(),
        "db_path": str(db_path.resolve()),
        "scoped_sources": list(SCOPED_SOURCES),
        "payloads": metadata,
        "before": before,
        "after": after,
        "compensation_ledger": str(ledger_path) if ledger_path else None,
        "compensation_ledger_sha256": ledger_sha256,
        "protected_tables": protected,
        "protected_snapshot_verified": protected_after == protected_before,
        "quick_check_before": quick_before,
        "quick_check_after": quick_after,
        "foreign_key_check_before": foreign_before,
        "foreign_key_check_after": foreign_after,
        "active_service_ports_during_nonproduction_write": (
            active_ports if db_path.resolve() != DEFAULT_DB.resolve() else []
        ),
        "formal_service_ports_verified_stopped": (
            db_path.resolve() != DEFAULT_DB.resolve() or not active_ports
        ),
    }
    output = receipt_path or _default_output_path(
        "zytw_bzl_split_migration_receipt", ledger_dir
    )
    result["receipt_path"] = str(output.resolve())
    _write_json_atomic(output, result)
    return result


def _restore_rows(
    db_path: Path, ledger_path: Path, receipt_path: Path | None
) -> dict[str, Any]:
    active_ports = _require_write_window(db_path)
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    if ledger.get("kind") != "zytw_bzl_split_compensation_ledger":
        raise ValueError("not a zytw/BZL split compensation ledger")
    if tuple(ledger.get("scoped_sources") or ()) != SCOPED_SOURCES:
        raise ValueError("ledger scope does not match this migration")
    if Path(ledger.get("db_path", "")).resolve() != db_path.resolve():
        raise ValueError("ledger was created from a different database")
    expected_post = ledger.get("expected_post") or {}
    with closing(_connect(db_path, read_only=True)) as connection:
        before = _database_summary(connection)
        protected = list(ledger.get("protected_tables") or [])
        protected_before = _protected_snapshot(connection, protected)
        quick_before = _quick_check(connection)
        foreign_before = _foreign_key_check(connection)
    if not _is_current(before, expected_post):
        raise RuntimeError("refusing rollback because current split data changed")
    if protected_before != ledger.get("protected_snapshot"):
        raise RuntimeError("refusing rollback because a protected table changed")
    if quick_before != "ok" or foreign_before:
        raise RuntimeError("database failed rollback preflight")

    dump = ledger.get("dataset_dump") or {}
    target = dump.get("TargetCompounds") or {}
    diagnostic = dump.get("DiagnosticIons") or {}
    comparable = {"TargetCompounds": target, "DiagnosticIons": diagnostic}
    if _sha256_bytes(_json_bytes(comparable)) != str(dump.get("rows_sha256") or ""):
        raise ValueError("compensation ledger row checksum mismatch")

    connection = _connect(db_path, read_only=False)
    committed = False
    try:
        _install_authorizer(connection, protected)
        connection.execute("BEGIN IMMEDIATE")
        _delete_scoped(connection)
        target_columns = [str(value) for value in target.get("columns") or []]
        target_rows = list(target.get("rows") or [])
        if target_columns and target_rows:
            columns = ",".join(_quote_identifier(value) for value in target_columns)
            marks = ",".join("?" for _ in target_columns)
            connection.executemany(
                f"INSERT INTO TargetCompounds ({columns}) VALUES ({marks})",
                target_rows,
            )
        diagnostic_columns = [
            str(value) for value in diagnostic.get("columns") or []
        ]
        diagnostic_rows = list(diagnostic.get("rows") or [])
        if diagnostic_columns and diagnostic_rows:
            columns = ",".join(
                _quote_identifier(value) for value in diagnostic_columns
            )
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
        protected_after = _protected_snapshot(connection, protected)
        quick_after = _quick_check(connection)
        foreign_after = _foreign_key_check(connection)
    if after != ledger.get("before"):
        raise RuntimeError("rollback did not reproduce the exact pre-split summary")
    if protected_after != protected_before:
        raise RuntimeError("a protected table changed during rollback")
    if quick_after != "ok" or foreign_after:
        raise RuntimeError("database failed rollback postflight")

    result = {
        "schema_version": 1,
        "kind": "zytw_bzl_split_rollback_receipt",
        "status": "rolled_back",
        "completed_at": _utc_now(),
        "db_path": str(db_path.resolve()),
        "ledger_path": str(ledger_path.resolve()),
        "ledger_sha256": _sha256_bytes(ledger_path.read_bytes()),
        "before": before,
        "after": after,
        "protected_snapshot_verified": protected_after == protected_before,
        "quick_check_before": quick_before,
        "quick_check_after": quick_after,
        "foreign_key_check_before": foreign_before,
        "foreign_key_check_after": foreign_after,
        "active_service_ports_during_nonproduction_write": (
            active_ports if db_path.resolve() != DEFAULT_DB.resolve() else []
        ),
    }
    output = receipt_path or ledger_path.with_name(
        ledger_path.stem + "_rollback_receipt.json"
    )
    result["receipt_path"] = str(output.resolve())
    _write_json_atomic(output, result)
    return result


def _preview(
    records: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
    db_path: Path,
) -> dict[str, Any]:
    with closing(_connect(db_path, read_only=True)) as connection:
        current = _database_summary(connection)
        protected = _protected_tables(connection)
        quick = _quick_check(connection)
        foreign = _foreign_key_check(connection)
    return {
        "schema_version": 1,
        "kind": "zytw_bzl_split_migration_preview",
        "status": "read_only",
        "db_path": str(db_path.resolve()),
        "prepared_records": len(records),
        "payloads": metadata,
        "current": current,
        "already_current": _is_current(current, metadata),
        "protected_tables": protected,
        "quick_check": quick,
        "foreign_key_check": foreign,
        "open_service_ports": _open_ports(SERVICE_PORTS),
        "write_requested": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zytw-payload", type=Path, default=SPECS["zytw"]["payload"])
    parser.add_argument("--bzl-payload", type=Path, default=SPECS["bzl"]["payload"])
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
        payload_paths = {"zytw": args.zytw_payload, "bzl": args.bzl_payload}
        for path in payload_paths.values():
            if not path.is_file():
                parser.error(f"payload does not exist: {path}")
        prepared: dict[str, list[dict[str, Any]]] = {}
        metadata: dict[str, dict[str, Any]] = {}
        for dataset in ("zytw", "bzl"):
            prepared[dataset], metadata[dataset] = _load_payload(
                dataset, payload_paths[dataset]
            )
        records = prepared["zytw"] + prepared["bzl"]
        if len(records) != 345:
            raise RuntimeError(f"expected 345 persisted targets, found {len(records)}")
        if args.apply:
            result = _apply_split(
                records, metadata, args.db, args.ledger_dir, args.receipt
            )
        else:
            result = _preview(records, metadata, args.db)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
