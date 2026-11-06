#!/usr/bin/env python3
"""Rebind accepted spectral assets after reviewed Jiawei target-table inserts.

This command does not rebuild or edit normalized spectral records or NumPy
indices.  It accepts only an installer-receipt-bound ``TargetCompounds``
increment, or the separate source-specific TargetCompounds/DiagnosticIons
installer receipt. Each installer verifies its exact permitted increments,
unchanged schema, unrelated targets/diagnostics and immutable spectral tables.
This command then updates the whole-database bindings coupled to the
accepted record store.  Dry-run is the default.
"""

from __future__ import annotations

import argparse
import errno
import json
import os
import shutil
import sqlite3
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


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
from scripts.install_jiawei_sanwujiao_targets import (  # noqa: E402
    DEFAULT_CATALOG,
    DEFAULT_DB,
    DEFAULT_RECEIPT as INSTALL_RECEIPT,
    _database_scope_snapshot,
    _existing_dataset,
    _semantic_projection,
    _sha256_json,
    load_catalog,
)


DEFAULT_OUTPUT = (
    PROJECT_ROOT.parent
    / "research/20261020_jiawei_sanwujiao_targeted/target_only_rebind_receipt.json"
)
AUDIT_FIELD = "post_acceptance_jiawei_target_only_source_rebind"
UNCHANGED_TABLES = ("Spectra", "Compounds", "Mapping", "DiagnosticIons")


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"expected a JSON object: {path}")
    return payload


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
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


def _file_stat(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "st_dev": stat.st_dev,
        "st_ino": stat.st_ino,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
    }


def _binding_once(db_path: Path, output_path: Path, identity: str) -> tuple[dict, bool]:
    """Reuse an unchanged preview binding so a 5 GB database is hashed once."""

    current_stat = _file_stat(db_path)
    if output_path.is_file():
        try:
            prior = _read_json(output_path)
            if (
                prior.get("status") == "preview"
                and prior.get("rebind_identity") == identity
                and prior.get("binding_file_stat") == current_stat
                and isinstance(prior.get("current_database_binding"), dict)
            ):
                return dict(prior["current_database_binding"]), True
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            pass
    return source_binding(db_path), False


def _verify_install_scope(
    db_path: Path,
    catalog_path: Path,
    receipt: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    catalog_info, expected_rows = load_catalog(catalog_path)
    if receipt.get("status") != "installed":
        raise RuntimeError("installer receipt must prove a newly installed dataset")
    if int(receipt.get("inserted_target_rows") or 0) != len(expected_rows):
        raise RuntimeError("installer receipt target increment differs from catalog")
    for key in (
        "dataset_source", "dataset_version", "source_content_sha256",
        "target_semantic_sha256",
    ):
        expected = (
            _sha256_json(_semantic_projection(expected_rows))
            if key == "target_semantic_sha256" else catalog_info.get(key)
        )
        if receipt.get(key) != expected:
            raise RuntimeError(f"installer receipt {key} differs from catalog")
    pre = receipt.get("pre_database_state")
    if not isinstance(pre, Mapping):
        raise RuntimeError("installer receipt lacks the pre-database proof")

    connection = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        current = _database_scope_snapshot(
            connection, db_path, catalog_info["dataset_source"]
        )
        if current.get("schema_sha256") != pre.get("schema_sha256"):
            raise RuntimeError("SQLite schema changed during target installation")
        if current.get("schema_version") != pre.get("schema_version"):
            raise RuntimeError("SQLite schema version changed during target installation")
        before_counts = dict(pre.get("table_row_counts") or {})
        current_counts = dict(current.get("table_row_counts") or {})
        missing = [table for table in (*UNCHANGED_TABLES, "TargetCompounds") if table not in before_counts or table not in current_counts]
        if missing:
            raise RuntimeError("installer proof lacks table counts: " + ", ".join(missing))
        for table in UNCHANGED_TABLES:
            if current_counts[table] != before_counts[table]:
                raise RuntimeError(f"{table} row count changed during target installation")
        expected_total = before_counts["TargetCompounds"] + len(expected_rows)
        if current_counts["TargetCompounds"] != expected_total:
            raise RuntimeError("TargetCompounds increment is not exactly the catalog size")
        if current.get("non_dataset_target_rows") != pre.get("non_dataset_target_rows"):
            raise RuntimeError("non-Jiawei target row count changed")
        if current.get("non_dataset_targets_sha256") != pre.get("non_dataset_targets_sha256"):
            raise RuntimeError("non-Jiawei TargetCompounds content changed")
        installed = _existing_dataset(connection, catalog_info["dataset_source"])
        if _sha256_json(_semantic_projection(installed)) != _sha256_json(
            _semantic_projection(expected_rows)
        ):
            raise RuntimeError("installed Jiawei target rows differ from catalog")
        diagnostic_rows = int(connection.execute(
            "SELECT COUNT(*) FROM DiagnosticIons d JOIN TargetCompounds t "
            "ON t.id=d.target_id WHERE t.dataset_source=?",
            (catalog_info["dataset_source"],),
        ).fetchone()[0])
        if diagnostic_rows:
            raise RuntimeError("Jiawei target install added DiagnosticIons rows")
    finally:
        connection.close()
    return catalog_info, expected_rows, current


def _backup(paths: list[Path], destination: Path) -> dict[str, str]:
    destination.mkdir(parents=True, exist_ok=False)
    copied: dict[str, str] = {}
    for index, path in enumerate(paths):
        target = destination / f"{index:02d}_{path.name}"
        # Every bound file below is replaced atomically. Retain its old inode
        # without copying the 100+ MiB metadata registry into a second file.
        # Cross-filesystem backups fail before any binding is changed.
        try:
            os.link(path, target)
        except OSError as exc:
            # The accepted store can live on an external volume. Its small
            # manifest/acceptance files may be copied; large registries still
            # require hardlinks and never create a second large local copy.
            if exc.errno != errno.EXDEV or path.stat().st_size > 1024 * 1024:
                raise
            shutil.copy2(path, target)
        copied[str(path)] = str(target)
    return copied


def _default_config_paths(project_root: Path) -> dict[str, Path]:
    return {
        "record_store": project_root / "config/library_record_store.json",
        "metadata_registry": project_root / "config/library_metadata_repairs.json",
        "sql_metadata": project_root / "config/sql_metadata_repairs.json",
        "mass_registry": project_root / "config/library_compound_mass_repairs.json",
        "mass_runtime": project_root / "config/source_compound_mass_repairs.json",
    }


def rebind(
    db_path: Path,
    catalog_path: Path,
    install_receipt_path: Path,
    output_path: Path,
    *,
    apply: bool,
    project_root: Path = PROJECT_ROOT,
    config_paths: Mapping[str, Path] | None = None,
) -> dict[str, Any]:
    receipt = _read_json(install_receipt_path)
    consolidation = receipt.get("schema_version") == "mssystem-jiawei-target-consolidation/1.0"
    msms_install = receipt.get("schema_version") == "mssystem-jiawei-msms-target-install/1.0"
    if consolidation:
        from scripts.consolidate_jiawei_sanwujiao_targets import verify_merge_scope
        catalog_info, expected_rows, scope = verify_merge_scope(db_path, catalog_path, receipt)
    elif msms_install:
        from scripts.install_jiawei_sanwujiao_msms_evidence import verify_install_scope
        catalog_info, expected_rows, scope = verify_install_scope(
            db_path, catalog_path, receipt
        )
    else:
        catalog_info, expected_rows, scope = _verify_install_scope(
            db_path, catalog_path, receipt
        )
    configs = dict(config_paths or _default_config_paths(project_root))
    record_store_config = _read_json(configs["record_store"])
    record_store = Path(str(record_store_config.get("cache_dir") or "")).resolve()
    manifest_path = record_store / "manifest.json"
    acceptance_path = record_store / "acceptance.json"
    manifest = _read_json(manifest_path)
    acceptance = _read_json(acceptance_path)
    validate_acceptance(acceptance, manifest)
    if record_store_config.get("accepted_cache_manifest_semantic_sha256") != manifest.get("semantic_sha256"):
        raise RuntimeError("record-store config is not bound to the accepted manifest")
    if record_store_config.get("acceptance_receipt_file_sha256") != file_sha256(acceptance_path):
        raise RuntimeError("record-store config is not bound to the acceptance receipt")
    from runtime.library_cache_identity import processing_identity_compatibility
    processing_proof = processing_identity_compatibility(
        manifest.get("processing_identity"), processing_identity(),
        project_root=project_root,
    )
    if not processing_proof.get("compatible"):
        raise RuntimeError("normalization processing identity changed")
    old_binding = dict(manifest.get("source_db") or {})
    if str(old_binding.get("path") or "") != str(db_path.resolve()):
        raise RuntimeError("record store is bound to a different database")
    pre_stat = (receipt.get("pre_database_state") or {}).get("file_stat") or {}
    if int(old_binding.get("size_bytes") or -1) != int(pre_stat.get("size_bytes") or -2):
        raise RuntimeError("accepted record store size differs from installer pre-state")

    immutable_record_identity = {
        key: manifest.get(key)
        for key in (
            "record_count", "normalized_records_sha256", "record_alignment_sha256",
            "mass_counts", "files",
        )
    }
    rebind_identity = _sha256_json({
        "catalog_sha256": catalog_info["source_content_sha256"],
        "target_sha256": receipt["target_semantic_sha256"],
        "old_database_sha256": old_binding.get("sha256"),
    })
    current_binding, binding_reused = _binding_once(
        db_path, output_path, rebind_identity
    )
    if current_binding == old_binding:
        raise RuntimeError("database still has the pre-install binding")
    sqlite_state = sqlite_snapshot(db_path)
    snapshot_id = str(sqlite_state.get("snapshot_id") or "")
    if not snapshot_id:
        raise RuntimeError("current SQLite snapshot could not be computed")

    metadata_registry = _read_json(configs["metadata_registry"])
    sql_metadata = _read_json(configs["sql_metadata"])
    mass_registry = _read_json(configs["mass_registry"])
    mass_runtime = _read_json(configs["mass_runtime"])
    old_sha = str(old_binding.get("sha256") or "")
    for name, payload in (
        ("metadata_registry", metadata_registry),
        ("sql_metadata", sql_metadata),
        ("mass_registry", mass_registry),
        ("mass_runtime", mass_runtime),
    ):
        if payload.get("source_database_sha256") != old_sha:
            raise RuntimeError(f"{name} is not bound to the accepted pre-install database")
    if sql_metadata.get("processing_identity") != manifest.get("processing_identity"):
        raise RuntimeError("SQL metadata processing identity differs from record store")
    if sql_metadata.get("cache_manifest_semantic_sha256") != manifest.get("semantic_sha256"):
        raise RuntimeError("SQL metadata manifest binding differs from record store")

    timestamp = datetime.now(timezone.utc).isoformat()
    audit = {
        "schema_version": "mssystem-jiawei-target-consolidation-rebind/1.0" if consolidation else ("mssystem-jiawei-msms-target-source-rebind/1.0" if msms_install else "mssystem-jiawei-target-only-source-rebind/1.0"),
        "rebound_at": timestamp,
        "dataset_source": catalog_info["dataset_source"],
        "dataset_version": catalog_info["dataset_version"],
        "catalog_sha256": catalog_info["source_content_sha256"],
        "install_receipt": str(install_receipt_path.resolve()),
        "install_receipt_sha256": file_sha256(install_receipt_path),
        "inserted_target_rows": 0 if consolidation else len(expected_rows),
        "consolidated_target_rows": len(expected_rows) if consolidation else 0,
        "consolidated_source_datasets": receipt.get("source_datasets", []) if consolidation else [],
        "inserted_diagnostic_rows": int(receipt.get("inserted_diagnostic_rows") or 0) if msms_install else 0,
        "diagnostic_semantic_sha256": receipt.get("diagnostic_semantic_sha256") if msms_install else None,
        "spectra_rows_changed": 0,
        "compound_rows_changed": 0,
        "mapping_rows_changed": 0,
        "non_dataset_target_rows_changed": 0,
        "normalized_record_content_changed": False,
        "score_inputs_changed": False,
        "score_inputs_scope": "normalized_spectral_library_records",
        "expected_compounds_product_ion_rules_changed": bool(msms_install),
        "processing_identity_compatibility": processing_proof,
        "previous_database_sha256": old_sha,
        "current_database_sha256": current_binding["sha256"],
        "verified_scope": scope,
    }

    new_manifest = dict(manifest)
    new_manifest["source_db"] = current_binding
    new_manifest[AUDIT_FIELD] = audit
    _refresh_semantic(new_manifest)
    new_acceptance = dict(acceptance)
    new_acceptance["manifest_semantic_sha256"] = new_manifest["semantic_sha256"]
    new_acceptance[AUDIT_FIELD] = audit
    _refresh_semantic(new_acceptance)

    metadata_registry["source_database_sha256"] = current_binding["sha256"]
    metadata_registry[AUDIT_FIELD] = audit
    _refresh_semantic(metadata_registry)
    mass_registry["source_database_sha256"] = current_binding["sha256"]
    mass_registry[AUDIT_FIELD] = audit
    _refresh_semantic(mass_registry)
    sql_metadata.update({
        "source_database_sha256": current_binding["sha256"],
        "library_snapshot_id": snapshot_id,
        "registry_semantic_sha256": metadata_registry["semantic_sha256"],
        "cache_manifest_semantic_sha256": new_manifest["semantic_sha256"],
        AUDIT_FIELD: audit,
    })
    _refresh_semantic(sql_metadata)
    mass_runtime.update({
        "source_database_sha256": current_binding["sha256"],
        "base_library_snapshot_id": snapshot_id,
        "registry_semantic_sha256": mass_registry["semantic_sha256"],
        AUDIT_FIELD: audit,
    })

    result = {
        "schema_version": "mssystem-jiawei-target-only-cache-rebind/1.0",
        "status": "preview",
        "generated_at": timestamp,
        "rebind_identity": rebind_identity,
        "database_path": str(db_path.resolve()),
        "binding_file_stat": _file_stat(db_path),
        "binding_reused_from_preview": binding_reused,
        "previous_database_binding": old_binding,
        "current_database_binding": current_binding,
        "current_sqlite_snapshot_id": snapshot_id,
        "verified_scope": scope,
        "record_store_record_count": manifest.get("record_count"),
        "normalized_records_sha256": manifest.get("normalized_records_sha256"),
        "record_store_payload_identity": immutable_record_identity,
        "new_record_store_manifest_semantic_sha256": new_manifest["semantic_sha256"],
        "backup_files": {},
        "backup_method": "hardlinks_or_small_cross_volume_copies_before_atomic_replacement",
    }
    if not apply:
        _write_json_atomic(output_path, result)
        return result

    backup_dir = output_path.parent / "pre_jiawei_target_only_rebind"
    result["backup_files"] = _backup(
        [manifest_path, acceptance_path, *configs.values()], backup_dir
    )
    _write_json_atomic(manifest_path, new_manifest)
    _write_json_atomic(acceptance_path, new_acceptance)
    record_store_config["accepted_cache_manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    record_store_config["acceptance_receipt_file_sha256"] = file_sha256(
        acceptance_path
    )
    record_store_config[AUDIT_FIELD] = audit
    _write_json_atomic(configs["record_store"], record_store_config)
    _write_json_atomic(configs["metadata_registry"], metadata_registry)
    _refresh_semantic(sql_metadata)
    _write_json_atomic(configs["sql_metadata"], sql_metadata)
    _write_json_atomic(configs["mass_registry"], mass_registry)
    mass_runtime["registry_sha256"] = file_sha256(configs["mass_registry"])
    mass_runtime["registry_semantic_sha256"] = mass_registry["semantic_sha256"]
    _refresh_semantic(mass_runtime)
    _write_json_atomic(configs["mass_runtime"], mass_runtime)

    installed_manifest = _read_json(manifest_path)
    installed_acceptance = _read_json(acceptance_path)
    validate_acceptance(installed_acceptance, installed_manifest)
    installed_identity = {
        key: installed_manifest.get(key) for key in immutable_record_identity
    }
    if installed_identity != immutable_record_identity:
        raise RuntimeError("normalized records or index-array metadata changed")
    installed_config = _read_json(configs["record_store"])
    if installed_config.get("accepted_cache_manifest_semantic_sha256") != installed_manifest.get("semantic_sha256"):
        raise RuntimeError("record-store config manifest digest mismatch")
    if installed_config.get("acceptance_receipt_file_sha256") != file_sha256(acceptance_path):
        raise RuntimeError("record-store config acceptance digest mismatch")
    result["status"] = "applied"
    result["applied_at"] = datetime.now(timezone.utc).isoformat()
    result["record_store_payload_unchanged"] = True
    _write_json_atomic(output_path, result)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--install-receipt", type=Path, default=INSTALL_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = rebind(
        args.db.expanduser().resolve(),
        args.catalog.expanduser().resolve(),
        args.install_receipt.expanduser().resolve(),
        args.output.expanduser().resolve(),
        apply=args.apply,
    )
    print(json.dumps({
        key: result[key]
        for key in (
            "status", "current_sqlite_snapshot_id", "record_store_record_count",
            "normalized_records_sha256", "new_record_store_manifest_semantic_sha256",
            "binding_reused_from_preview",
        )
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
