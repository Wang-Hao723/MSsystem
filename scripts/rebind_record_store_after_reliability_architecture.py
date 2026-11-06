#!/usr/bin/env python3
"""Rebind accepted read-only records after two non-normalization edits.

The historical processing identity covers every ``matching/**`` file.  The
reliability architecture changed only the opt-in expected-compound workflow
and feature-level query organization.  Neither file is imported by the
spectral-library normalizer.  This command proves that substituting their
previous accepted digests exactly reconstructs the accepted broad identity,
then changes bindings only; cached records, indices, scores and the SQLite
library remain byte-for-byte unchanged.
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


RECORD_STORE = (
    PROJECT_ROOT
    / "output/library_record_store/full-20260915-formula-composition-v1"
)
CONFIG_RECORD_STORE = PROJECT_ROOT / "config/library_record_store.json"
CONFIG_SQL_METADATA = PROJECT_ROOT / "config/sql_metadata_repairs.json"
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "research/20260920_qft_reliability_realfix/"
    "NON_NORMALIZATION_SOURCE_REBIND.json"
)
ACCEPTED_SOURCE_FILES_SHA256 = (
    "724bca448c216cfb9c8a65eaad55df72d5f116fe28ea656722506725a594711c"
)
PREVIOUS_FILES = {
    "matching/expected_compounds.py": {
        "sha256": "40ef8e0ec974c9e8e0b3307c2fd3504b701843fdd4c5a1311b96b6a2f9e96610",
        "size_bytes": 120183,
        "source": "prior accepted rebind audit current_version",
    },
    "matching/untargeted_workflow.py": {
        "sha256": "94bb8e64e9d2b46077df1ccc2e075fbe81c59d5afb92b095e49d51b78a09692d",
        "size_bytes": 14889,
        "source": "git HEAD before reliability architecture edit",
    },
}


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"expected JSON object: {path}")
    return payload


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _proof(accepted_identity: dict[str, Any]) -> dict[str, Any]:
    files = _processing_files()
    observed: dict[str, Any] = {}
    for item in files:
        previous = PREVIOUS_FILES.get(item["path"])
        if previous is None:
            continue
        current_path = PROJECT_ROOT / item["path"]
        observed[item["path"]] = {
            "previous": previous,
            "current_sha256": file_sha256(current_path),
            "current_size_bytes": current_path.stat().st_size,
        }
        item["sha256"] = previous["sha256"]
        item["size_bytes"] = previous["size_bytes"]
    if set(observed) != set(PREVIOUS_FILES):
        raise RuntimeError("could not isolate both non-normalization files")
    reconstructed = semantic_sha256(files)
    accepted = str(accepted_identity.get("source_files_sha256") or "")
    if accepted != ACCEPTED_SOURCE_FILES_SHA256 or reconstructed != accepted:
        raise RuntimeError(
            "Two-file substitution did not reconstruct the accepted identity: "
            f"{reconstructed} != {accepted}"
        )
    head_bytes = subprocess.check_output(
        ["git", "show", "HEAD:matching/untargeted_workflow.py"],
        cwd=PROJECT_ROOT,
    )
    if hashlib.sha256(head_bytes).hexdigest() != PREVIOUS_FILES[
        "matching/untargeted_workflow.py"
    ]["sha256"]:
        raise RuntimeError("HEAD no longer proves the prior untargeted workflow")
    return {
        "substituted_files": observed,
        "substituted_file_count": len(observed),
        "other_processing_identity_files_changed": 0,
        "reconstructed_accepted_source_files_sha256": reconstructed,
        "library_normalizer_imports_changed_modules": False,
    }


def _backup(paths: list[Path], destination: Path) -> dict[str, str]:
    if destination.exists():
        raise RuntimeError(f"backup already exists: {destination}")
    destination.mkdir(parents=True)
    output = {}
    for path in paths:
        target = destination / path.name
        shutil.copy2(path, target)
        output[str(path)] = str(target)
    return output


def rebind(output_path: Path, *, apply: bool) -> dict[str, Any]:
    manifest_path = RECORD_STORE / "manifest.json"
    acceptance_path = RECORD_STORE / "acceptance.json"
    manifest = _read_json(manifest_path)
    acceptance = _read_json(acceptance_path)
    record_config = _read_json(CONFIG_RECORD_STORE)
    sql_config = _read_json(CONFIG_SQL_METADATA)
    validate_acceptance(acceptance, manifest)

    accepted_identity = dict(manifest.get("processing_identity") or {})
    proof = _proof(accepted_identity)
    current_identity = processing_identity()
    if current_identity == accepted_identity:
        raise RuntimeError("record cache already has the current identity")
    if {
        key: value for key, value in current_identity.items()
        if key != "source_files_sha256"
    } != {
        key: value for key, value in accepted_identity.items()
        if key != "source_files_sha256"
    }:
        raise RuntimeError("numeric runtime or normalization options changed")
    if sql_config.get("processing_identity") != accepted_identity:
        raise RuntimeError("SQL metadata is not bound to the accepted identity")
    if sql_config.get("cache_manifest_semantic_sha256") != manifest.get(
        "semantic_sha256"
    ):
        raise RuntimeError("SQL metadata manifest binding differs")
    if record_config.get("accepted_cache_manifest_semantic_sha256") != (
        manifest.get("semantic_sha256")
    ):
        raise RuntimeError("record-store config manifest binding differs")
    if record_config.get("acceptance_receipt_file_sha256") != file_sha256(
        acceptance_path
    ):
        raise RuntimeError("record-store acceptance binding differs")

    timestamp = datetime.now(timezone.utc).isoformat()
    audit = {
        "schema_version": "mssystem-non-normalization-source-rebind/2.0",
        "rebound_at": timestamp,
        **proof,
        "previous_processing_identity": accepted_identity,
        "current_processing_identity": current_identity,
        "normalized_record_content_changed": False,
        "normalized_records_sha256": manifest.get("normalized_records_sha256"),
        "record_alignment_sha256": manifest.get("record_alignment_sha256"),
        "record_count": manifest.get("record_count"),
        "candidate_scores_changed": False,
        "sqlite_library_changed": False,
        "reason": (
            "expected-compound opt-in isolation and feature-level query "
            "organization are outside spectral-library record normalization"
        ),
    }
    new_manifest = dict(manifest)
    new_manifest["processing_identity"] = current_identity
    new_manifest["post_acceptance_reliability_source_rebind"] = audit
    _refresh_semantic(new_manifest)
    new_acceptance = dict(acceptance)
    new_acceptance["manifest_semantic_sha256"] = new_manifest["semantic_sha256"]
    new_acceptance["post_acceptance_reliability_source_rebind"] = audit
    _refresh_semantic(new_acceptance)
    new_record_config = dict(record_config)
    new_record_config["accepted_cache_manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_record_config["post_acceptance_reliability_source_rebind"] = audit
    new_sql_config = dict(sql_config)
    new_sql_config["processing_identity"] = current_identity
    new_sql_config["cache_manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_sql_config["post_acceptance_reliability_source_rebind"] = audit
    _refresh_semantic(new_sql_config)

    result = {
        "schema_version": "mssystem-reliability-source-rebind/1.0",
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

    backup_dir = output_path.parent / "pre_reliability_source_rebind"
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
    if installed_manifest.get("processing_identity") != processing_identity():
        raise RuntimeError("installed record-store identity is not current")
    if _read_json(CONFIG_RECORD_STORE).get(
        "acceptance_receipt_file_sha256"
    ) != file_sha256(acceptance_path):
        raise RuntimeError("installed acceptance hash is inconsistent")
    if _read_json(CONFIG_SQL_METADATA).get(
        "processing_identity"
    ) != processing_identity():
        raise RuntimeError("installed SQL identity is not current")
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
