#!/usr/bin/env python3
"""Rebind the accepted record cache after an Expected Compounds-only edit.

``matching.expected_compounds`` is included by the historical broad
``matching/**`` processing-identity rule even though it never normalizes a
spectral-library record.  This command proves that replacing only that file by
its repository version exactly reconstructs the accepted processing digest,
then updates the cache and SQL metadata identities without rebuilding or
changing any cached normalized record.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
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
    validate_acceptance,
)
from runtime.manifest import source_tree_snapshot  # noqa: E402


CHANGED_FILE = "matching/expected_compounds.py"
RECORD_STORE = (
    PROJECT_ROOT
    / "output/library_record_store/full-20260915-formula-composition-v1"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "research/20260920_qft_negative_panel/"
    "EXPECTED_WORKFLOW_SOURCE_REBIND_RECEIPT.json"
)
CONFIG_RECORD_STORE = PROJECT_ROOT / "config/library_record_store.json"
CONFIG_SQL_METADATA = PROJECT_ROOT / "config/sql_metadata_repairs.json"


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"expected JSON object: {path}")
    return payload


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _refresh_semantic(payload: dict[str, Any]) -> None:
    payload["semantic_sha256"] = semantic_sha256(
        {key: value for key, value in payload.items() if key != "semantic_sha256"}
    )


def _processing_files() -> list[dict[str, Any]]:
    return [
        dict(item)
        for item in source_tree_snapshot(PROJECT_ROOT)["files"]
        if item["path"].startswith(("matching/", "utils/", "mssystem/"))
        or item["path"]
        in {"config/settings.py", "config/retrieval_policy.json"}
    ]


def _accepted_identity_reconstruction(
    accepted_identity: dict[str, Any],
) -> dict[str, Any]:
    previous_bytes = subprocess.check_output(
        ["git", "show", f"HEAD:{CHANGED_FILE}"], cwd=PROJECT_ROOT
    )
    previous_sha = hashlib.sha256(previous_bytes).hexdigest()
    current_path = PROJECT_ROOT / CHANGED_FILE
    current_sha = file_sha256(current_path)
    files = _processing_files()
    matched = 0
    for item in files:
        if item["path"] != CHANGED_FILE:
            continue
        item["size_bytes"] = len(previous_bytes)
        item["sha256"] = previous_sha
        matched += 1
    if matched != 1:
        raise RuntimeError(f"could not isolate {CHANGED_FILE} in source snapshot")
    reconstructed = semantic_sha256(files)
    expected = str(accepted_identity.get("source_files_sha256") or "")
    if reconstructed != expected:
        raise RuntimeError(
            "Replacing only matching.expected_compounds does not reconstruct "
            f"the accepted processing identity: {reconstructed} != {expected}"
        )
    if current_sha == previous_sha:
        raise RuntimeError("Expected Compounds source has not changed")
    return {
        "changed_file": CHANGED_FILE,
        "repository_version_sha256": previous_sha,
        "repository_version_size_bytes": len(previous_bytes),
        "current_version_sha256": current_sha,
        "current_version_size_bytes": current_path.stat().st_size,
        "reconstructed_accepted_source_files_sha256": reconstructed,
        "other_processing_identity_files_changed": 0,
    }


def _backup(paths: list[Path], destination: Path) -> dict[str, str]:
    if destination.exists():
        raise RuntimeError(f"backup already exists: {destination}")
    destination.mkdir(parents=True)
    result = {}
    for path in paths:
        target = destination / path.name
        shutil.copy2(path, target)
        result[str(path)] = str(target)
    return result


def rebind(output_path: Path, *, apply: bool) -> dict[str, Any]:
    manifest_path = RECORD_STORE / "manifest.json"
    acceptance_path = RECORD_STORE / "acceptance.json"
    manifest = _read_json(manifest_path)
    acceptance = _read_json(acceptance_path)
    record_config = _read_json(CONFIG_RECORD_STORE)
    sql_config = _read_json(CONFIG_SQL_METADATA)
    validate_acceptance(acceptance, manifest)

    accepted_identity = dict(manifest.get("processing_identity") or {})
    proof = _accepted_identity_reconstruction(accepted_identity)
    current_identity = processing_identity()
    if current_identity == accepted_identity:
        raise RuntimeError("record cache already has the current identity")
    if {
        key: value
        for key, value in current_identity.items()
        if key != "source_files_sha256"
    } != {
        key: value
        for key, value in accepted_identity.items()
        if key != "source_files_sha256"
    }:
        raise RuntimeError("numeric runtime or library processing options changed")
    if sql_config.get("processing_identity") != accepted_identity:
        raise RuntimeError("SQL metadata is not bound to the accepted old identity")
    if sql_config.get("cache_manifest_semantic_sha256") != manifest.get(
        "semantic_sha256"
    ):
        raise RuntimeError("SQL metadata cache manifest binding differs")
    if record_config.get("accepted_cache_manifest_semantic_sha256") != (
        manifest.get("semantic_sha256")
    ):
        raise RuntimeError("record-store config manifest binding differs")
    if record_config.get("acceptance_receipt_file_sha256") != file_sha256(
        acceptance_path
    ):
        raise RuntimeError("record-store acceptance file binding differs")

    timestamp = datetime.now(timezone.utc).isoformat()
    audit = {
        "schema_version": "mssystem-non-normalization-source-rebind/1.0",
        "rebound_at": timestamp,
        **proof,
        "previous_processing_identity": accepted_identity,
        "current_processing_identity": current_identity,
        "normalized_record_content_changed": False,
        "normalized_records_sha256": manifest.get("normalized_records_sha256"),
        "record_alignment_sha256": manifest.get("record_alignment_sha256"),
        "record_count": manifest.get("record_count"),
        "candidate_scores_changed": False,
        "reason": (
            "matching.expected_compounds implements a separate targeted MS1 "
            "workflow and is not called while normalized library records are built"
        ),
    }
    new_manifest = dict(manifest)
    new_manifest["processing_identity"] = current_identity
    new_manifest["post_acceptance_non_normalization_source_rebind"] = audit
    _refresh_semantic(new_manifest)
    new_acceptance = dict(acceptance)
    new_acceptance["manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_acceptance["post_acceptance_non_normalization_source_rebind"] = audit
    _refresh_semantic(new_acceptance)
    new_record_config = dict(record_config)
    new_record_config["accepted_cache_manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_record_config["post_acceptance_non_normalization_source_rebind"] = audit
    new_sql_config = dict(sql_config)
    new_sql_config["processing_identity"] = current_identity
    new_sql_config["cache_manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_sql_config["post_acceptance_non_normalization_source_rebind"] = audit
    _refresh_semantic(new_sql_config)

    result = {
        "schema_version": "mssystem-expected-workflow-source-rebind/1.0",
        "status": "preview",
        "generated_at": timestamp,
        "record_store": str(RECORD_STORE),
        "proof": proof,
        "audit": audit,
        "previous_manifest_semantic_sha256": manifest["semantic_sha256"],
        "new_manifest_semantic_sha256": new_manifest["semantic_sha256"],
        "backup_files": {},
    }
    if not apply:
        _write_json_atomic(output_path, result)
        return result

    backup_dir = output_path.parent / "pre_expected_workflow_source_rebind"
    result["backup_files"] = _backup(
        [manifest_path, acceptance_path, CONFIG_RECORD_STORE, CONFIG_SQL_METADATA],
        backup_dir,
    )
    _write_json_atomic(manifest_path, new_manifest)
    _write_json_atomic(acceptance_path, new_acceptance)
    new_record_config["acceptance_receipt_file_sha256"] = file_sha256(
        acceptance_path
    )
    _write_json_atomic(CONFIG_RECORD_STORE, new_record_config)
    _write_json_atomic(CONFIG_SQL_METADATA, new_sql_config)

    installed_manifest = _read_json(manifest_path)
    installed_acceptance = _read_json(acceptance_path)
    validate_acceptance(installed_acceptance, installed_manifest)
    installed_record_config = _read_json(CONFIG_RECORD_STORE)
    installed_sql_config = _read_json(CONFIG_SQL_METADATA)
    if installed_manifest.get("processing_identity") != processing_identity():
        raise RuntimeError("installed record-store identity is not current")
    if installed_record_config.get(
        "accepted_cache_manifest_semantic_sha256"
    ) != installed_manifest.get("semantic_sha256"):
        raise RuntimeError("installed record-store config is inconsistent")
    if installed_record_config.get(
        "acceptance_receipt_file_sha256"
    ) != file_sha256(acceptance_path):
        raise RuntimeError("installed acceptance file hash is inconsistent")
    if installed_sql_config.get("processing_identity") != processing_identity():
        raise RuntimeError("installed SQL metadata identity is not current")
    if installed_sql_config.get(
        "cache_manifest_semantic_sha256"
    ) != installed_manifest.get("semantic_sha256"):
        raise RuntimeError("installed SQL metadata manifest binding is inconsistent")
    result["status"] = "applied"
    result["applied_at"] = datetime.now(timezone.utc).isoformat()
    result["installed_acceptance_semantic_sha256"] = installed_acceptance[
        "semantic_sha256"
    ]
    _write_json_atomic(output_path, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    result = rebind(args.output.resolve(), apply=args.apply)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
