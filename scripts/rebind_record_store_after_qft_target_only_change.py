#!/usr/bin/env python3
"""Rebind accepted read-only spectral assets after the QFT target-only insert.

The production SQLite file also contains the presence-only ``TargetCompounds``
table.  Adding the audited QFT panel changes the whole-file hash even though no
spectral record, compound, mapping, diagnostic ion, score input, or cached
normalized record changes.  This command verifies that narrow state, preserves
the original full-record and scoring acceptance, and updates only source-file
bindings that intentionally cover the whole SQLite file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sqlite3
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.library_record_store import (  # noqa: E402
    file_sha256,
    processing_identity,
    semantic_sha256,
    source_binding,
    validate_acceptance,
)
from runtime.manifest import sqlite_snapshot  # noqa: E402
from scripts.install_qft_negative_panel import (  # noqa: E402
    DATASET_SOURCE,
    DEFAULT_DB,
    DEFAULT_RECEIPT as PANEL_RECEIPT,
)


OLD_DATABASE_SHA256 = (
    "c44332bea84846a0dd7e7f9108af2dd5eaf180784395e246eb76a297b9a1a8d3"
)
OLD_DATABASE_SIZE_BYTES = 5_396_340_736
ACCEPTED_PROCESSING_SHA256 = (
    "469631f58584aab5381b0a56cb809ac7daa14a9b650e8f68315367bb204f7a44"
)
EXPECTED_TABLE_COUNTS = {
    "Spectra": 4_674_363,
    "Compounds": 1_436_463,
    "TargetCompounds": 38_107,
    "DiagnosticIons": 60_202,
}
RECORD_STORE = (
    PROJECT_ROOT
    / "output/library_record_store/full-20260915-formula-composition-v1"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "research/20260920_qft_negative_panel/"
    "TARGET_ONLY_CACHE_REBIND_RECEIPT.json"
)
CONFIG_PATHS = {
    "record_store": PROJECT_ROOT / "config/library_record_store.json",
    "metadata_registry": PROJECT_ROOT / "config/library_metadata_repairs.json",
    "sql_metadata": PROJECT_ROOT / "config/sql_metadata_repairs.json",
    "mass_registry": PROJECT_ROOT / "config/library_compound_mass_repairs.json",
    "mass_runtime": PROJECT_ROOT / "config/source_compound_mass_repairs.json",
}


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"expected a JSON object: {path}")
    return payload


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    handle, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _refresh_semantic(payload: dict[str, Any]) -> None:
    payload["semantic_sha256"] = semantic_sha256(
        {key: value for key, value in payload.items() if key != "semantic_sha256"}
    )


def _table_counts(connection: sqlite3.Connection) -> dict[str, int]:
    return {
        table: int(connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
        for table in (
            "Spectra",
            "Compounds",
            "Mapping",
            "TargetCompounds",
            "DiagnosticIons",
        )
    }


def _verify_qft_scope(db_path: Path) -> dict[str, Any]:
    connection = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        quick_check = str(connection.execute("PRAGMA quick_check").fetchone()[0])
        if quick_check != "ok":
            raise RuntimeError(f"SQLite quick_check failed: {quick_check}")
        counts = _table_counts(connection)
        for table, expected in EXPECTED_TABLE_COUNTS.items():
            if counts[table] != expected:
                raise RuntimeError(
                    f"unexpected {table} count: {counts[table]} != {expected}"
                )
        rows = connection.execute(
            """
            SELECT id,target_id,compound_name,inchikey,adduct,source_tier
            FROM TargetCompounds
            WHERE dataset_source=?
            ORDER BY id
            """,
            (DATASET_SOURCE,),
        ).fetchall()
        if len(rows) != 8:
            raise RuntimeError(f"QFT panel row count is {len(rows)}, expected 8")
        if len({str(row["inchikey"]) for row in rows}) != 8:
            raise RuntimeError("QFT panel full InChIKeys are not unique")
        if {str(row["adduct"]) for row in rows} != {"[M-H]-"}:
            raise RuntimeError("QFT panel contains an unexpected adduct")
        if {str(row["source_tier"]) for row in rows} != {
            "literature_presence_only"
        }:
            raise RuntimeError("QFT panel contains an unexpected source tier")
        diagnostics = int(
            connection.execute(
                """
                SELECT COUNT(*)
                FROM DiagnosticIons d
                JOIN TargetCompounds t ON t.id=d.target_id
                WHERE t.dataset_source=?
                """,
                (DATASET_SOURCE,),
            ).fetchone()[0]
        )
        if diagnostics != 0:
            raise RuntimeError("QFT install unexpectedly added diagnostic ions")
        return {
            "quick_check": quick_check,
            "table_row_counts": counts,
            "qft_target_row_ids": [int(row["id"]) for row in rows],
            "qft_target_ids": [str(row["target_id"]) for row in rows],
            "qft_compound_names": [str(row["compound_name"]) for row in rows],
            "qft_diagnostic_ion_rows": diagnostics,
        }
    finally:
        connection.close()


def _backup(paths: list[Path], destination: Path) -> dict[str, str]:
    destination.mkdir(parents=True, exist_ok=True)
    output: dict[str, str] = {}
    for path in paths:
        target = destination / path.name
        if target.exists():
            raise RuntimeError(f"backup already exists: {target}")
        shutil.copy2(path, target)
        output[str(path)] = str(target)
    return output


def rebind(db_path: Path, output_path: Path, *, apply: bool) -> dict[str, Any]:
    manifest_path = RECORD_STORE / "manifest.json"
    acceptance_path = RECORD_STORE / "acceptance.json"
    manifest = _read_json(manifest_path)
    acceptance = _read_json(acceptance_path)
    validate_acceptance(acceptance, manifest)

    old_binding = dict(manifest.get("source_db") or {})
    if old_binding.get("sha256") != OLD_DATABASE_SHA256:
        raise RuntimeError("record store is not bound to the expected pre-QFT database")
    if int(old_binding.get("size_bytes") or 0) != OLD_DATABASE_SIZE_BYTES:
        raise RuntimeError("record store pre-QFT database size is unexpected")
    identity = processing_identity()
    if identity != manifest.get("processing_identity"):
        raise RuntimeError("current normalization identity differs from the accepted cache")
    if identity.get("source_files_sha256") != ACCEPTED_PROCESSING_SHA256:
        raise RuntimeError("normalization source hash is not the accepted value")
    if not PANEL_RECEIPT.is_file():
        raise RuntimeError(f"missing QFT panel receipt: {PANEL_RECEIPT}")
    panel_receipt = _read_json(PANEL_RECEIPT)
    if (
        panel_receipt.get("status") not in {"installed", "already_installed_exact"}
        or int(panel_receipt.get("panel_row_count") or 0) != 8
        or panel_receipt.get("dataset_source") != DATASET_SOURCE
    ):
        raise RuntimeError("QFT panel receipt does not prove the expected insert")

    scope = _verify_qft_scope(db_path)
    current_binding = source_binding(db_path)
    if current_binding["sha256"] == OLD_DATABASE_SHA256:
        raise RuntimeError("database still has the pre-QFT hash; no rebind is needed")
    sqlite_state = sqlite_snapshot(db_path)
    current_snapshot_id = str(sqlite_state.get("snapshot_id") or "")
    if not current_snapshot_id:
        raise RuntimeError("could not compute the current SQLite snapshot id")

    metadata_registry = _read_json(CONFIG_PATHS["metadata_registry"])
    sql_metadata = _read_json(CONFIG_PATHS["sql_metadata"])
    mass_registry = _read_json(CONFIG_PATHS["mass_registry"])
    mass_runtime = _read_json(CONFIG_PATHS["mass_runtime"])
    record_store_config = _read_json(CONFIG_PATHS["record_store"])
    for name, payload in (
        ("metadata_registry", metadata_registry),
        ("sql_metadata", sql_metadata),
        ("mass_registry", mass_registry),
        ("mass_runtime", mass_runtime),
    ):
        if payload.get("source_database_sha256") != OLD_DATABASE_SHA256:
            raise RuntimeError(f"{name} is not bound to the pre-QFT database")

    timestamp = datetime.now(timezone.utc).isoformat()
    audit = {
        "schema_version": "mssystem-target-only-source-rebind/1.0",
        "rebound_at": timestamp,
        "dataset_source": DATASET_SOURCE,
        "inserted_target_rows": 8,
        "inserted_diagnostic_rows": 0,
        "spectra_rows_changed": 0,
        "compound_rows_changed": 0,
        "mapping_rows_changed": 0,
        "normalized_record_content_changed": False,
        "score_inputs_changed": False,
        "previous_database_sha256": OLD_DATABASE_SHA256,
        "current_database_sha256": current_binding["sha256"],
        "panel_receipt": str(PANEL_RECEIPT.resolve()),
        "panel_receipt_sha256": file_sha256(PANEL_RECEIPT),
    }

    new_manifest = dict(manifest)
    new_manifest["source_db"] = current_binding
    new_manifest["post_acceptance_target_only_source_rebind"] = audit
    _refresh_semantic(new_manifest)
    new_acceptance = dict(acceptance)
    new_acceptance["manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_acceptance["post_acceptance_target_only_source_rebind"] = audit
    _refresh_semantic(new_acceptance)

    metadata_registry["source_database_sha256"] = current_binding["sha256"]
    metadata_registry["post_acceptance_target_only_source_rebind"] = audit
    _refresh_semantic(metadata_registry)
    mass_registry["source_database_sha256"] = current_binding["sha256"]
    mass_registry["post_acceptance_target_only_source_rebind"] = audit
    _refresh_semantic(mass_registry)

    sql_metadata.update(
        {
            "source_database_sha256": current_binding["sha256"],
            "library_snapshot_id": current_snapshot_id,
            "registry_semantic_sha256": metadata_registry["semantic_sha256"],
            "cache_manifest_semantic_sha256": new_manifest[
                "semantic_sha256"
            ],
            "post_acceptance_target_only_source_rebind": audit,
        }
    )
    _refresh_semantic(sql_metadata)
    mass_runtime.update(
        {
            "source_database_sha256": current_binding["sha256"],
            "base_library_snapshot_id": current_snapshot_id,
            "registry_semantic_sha256": mass_registry["semantic_sha256"],
            "post_acceptance_target_only_source_rebind": audit,
        }
    )

    preview = {
        "schema_version": "mssystem-qft-target-only-cache-rebind/1.0",
        "status": "preview",
        "generated_at": timestamp,
        "database_path": str(db_path.resolve()),
        "previous_database_binding": old_binding,
        "current_database_binding": current_binding,
        "current_sqlite_snapshot_id": current_snapshot_id,
        "verified_scope": scope,
        "processing_identity": identity,
        "record_store_normalized_records_sha256": manifest[
            "normalized_records_sha256"
        ],
        "record_store_record_count": manifest["record_count"],
        "new_record_store_manifest_semantic_sha256": new_manifest[
            "semantic_sha256"
        ],
        "audit": audit,
        "backup_files": {},
    }
    if not apply:
        _write_json_atomic(output_path, preview)
        return preview

    backup_dir = output_path.parent / "pre_target_only_cache_rebind"
    paths = [
        manifest_path,
        acceptance_path,
        *CONFIG_PATHS.values(),
    ]
    preview["backup_files"] = _backup(paths, backup_dir)

    _write_json_atomic(manifest_path, new_manifest)
    _write_json_atomic(acceptance_path, new_acceptance)
    record_store_config[
        "accepted_cache_manifest_semantic_sha256"
    ] = new_manifest["semantic_sha256"]
    record_store_config["acceptance_receipt_file_sha256"] = file_sha256(
        acceptance_path
    )
    record_store_config["post_acceptance_target_only_source_rebind"] = audit
    _write_json_atomic(CONFIG_PATHS["record_store"], record_store_config)

    _write_json_atomic(CONFIG_PATHS["metadata_registry"], metadata_registry)
    sql_metadata["registry_semantic_sha256"] = metadata_registry[
        "semantic_sha256"
    ]
    _refresh_semantic(sql_metadata)
    _write_json_atomic(CONFIG_PATHS["sql_metadata"], sql_metadata)

    _write_json_atomic(CONFIG_PATHS["mass_registry"], mass_registry)
    mass_runtime["registry_sha256"] = file_sha256(
        CONFIG_PATHS["mass_registry"]
    )
    mass_runtime["registry_semantic_sha256"] = mass_registry[
        "semantic_sha256"
    ]
    _refresh_semantic(mass_runtime)
    _write_json_atomic(CONFIG_PATHS["mass_runtime"], mass_runtime)

    # Validate every coupled digest after publication.
    installed_manifest = _read_json(manifest_path)
    installed_acceptance = _read_json(acceptance_path)
    validate_acceptance(installed_acceptance, installed_manifest)
    installed_config = _read_json(CONFIG_PATHS["record_store"])
    if installed_config["accepted_cache_manifest_semantic_sha256"] != (
        installed_manifest["semantic_sha256"]
    ):
        raise RuntimeError("record store config manifest digest mismatch")
    if installed_config["acceptance_receipt_file_sha256"] != file_sha256(
        acceptance_path
    ):
        raise RuntimeError("record store config acceptance digest mismatch")

    preview["status"] = "applied"
    preview["applied_at"] = datetime.now(timezone.utc).isoformat()
    preview["installed_record_store_manifest_semantic_sha256"] = (
        installed_manifest["semantic_sha256"]
    )
    preview["installed_acceptance_semantic_sha256"] = installed_acceptance[
        "semantic_sha256"
    ]
    _write_json_atomic(output_path, preview)
    return preview


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = rebind(
        args.db.expanduser().resolve(),
        args.output.expanduser().resolve(),
        apply=args.apply,
    )
    print(json.dumps({
        key: result[key]
        for key in (
            "status",
            "previous_database_binding",
            "current_database_binding",
            "current_sqlite_snapshot_id",
            "record_store_record_count",
            "new_record_store_manifest_semantic_sha256",
        )
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
