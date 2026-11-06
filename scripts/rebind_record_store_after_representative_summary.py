#!/usr/bin/env python3
"""Rebind the accepted record store after a reporting-only summary repair.

``matching.untargeted_workflow`` is covered by the historical broad source
identity even though the library normalizer does not import it.  This command
proves that the only covered-file change since the accepted generation is the
replacement of three misleading consensus counters with explicit
representative-spectrum counters.  It then updates bindings without changing
cached records, mass indices, candidate scores, or the SQLite library.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
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
CHANGED_PATH = "matching/untargeted_workflow.py"
PREVIOUS_FILE_SHA256 = (
    "f66ae7eb9c6f00f5cc99a39e4265e26eba68a9f08e9953c721d945a7a8b80d5a"
)
PREVIOUS_FILE_SIZE = 18_497
CURRENT_FILE_SHA256 = (
    "5bc165c2b6a9bc5cfc00ed5f951e5045c9d5f7bb86d54fd7d0b68a558d345998"
)
CURRENT_FILE_SIZE = 18_581
ACCEPTED_SOURCE_FILES_SHA256 = (
    "3fd33f519e886f79e99965903d79c5f1c8a38641ca003315da88b9a3554d6ee3"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "research/20260920_general_reliability_plan_review/"
    "REPRESENTATIVE_SUMMARY_SOURCE_REBIND.json"
)

CURRENT_BLOCK = """        \"consensus_spectra\": 0,
        \"actual_consensus_spectra\": 0,
        \"representative_spectra\": len(artifacts),
        \"multi_spectrum_consensus\": 0,
        \"single_spectrum_consensus\": 0,
        \"multi_spectrum_representative_groups\": repeated,
        \"single_spectrum_representatives\": len(artifacts) - repeated,
"""
PREVIOUS_BLOCK = """        \"consensus_spectra\": len(artifacts),
        \"actual_consensus_spectra\": 0,
        \"representative_spectra\": len(artifacts),
        \"multi_spectrum_consensus\": repeated,
        \"single_spectrum_consensus\": len(artifacts) - repeated,
"""


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
        or item["path"] in {"config/settings.py", "config/retrieval_policy.json"}
    ]


def _proof(manifest: dict[str, Any]) -> dict[str, Any]:
    changed = PROJECT_ROOT / CHANGED_PATH
    current_bytes = changed.read_bytes()
    if (
        len(current_bytes) != CURRENT_FILE_SIZE
        or hashlib.sha256(current_bytes).hexdigest() != CURRENT_FILE_SHA256
    ):
        raise RuntimeError("the reviewed representative-summary source changed")

    current_text = current_bytes.decode("utf-8")
    if current_text.count(CURRENT_BLOCK) != 1:
        raise RuntimeError("could not isolate the reviewed summary-only edit")
    previous_bytes = current_text.replace(
        CURRENT_BLOCK, PREVIOUS_BLOCK, 1
    ).encode("utf-8")
    if (
        len(previous_bytes) != PREVIOUS_FILE_SIZE
        or hashlib.sha256(previous_bytes).hexdigest() != PREVIOUS_FILE_SHA256
    ):
        raise RuntimeError("reverse patch did not reproduce the accepted source")

    prior_audit = manifest.get("post_acceptance_reliability_source_rebind") or {}
    prior_observation = (
        (prior_audit.get("substituted_files") or {}).get(CHANGED_PATH) or {}
    )
    if (
        prior_observation.get("current_sha256") != PREVIOUS_FILE_SHA256
        or int(prior_observation.get("current_size_bytes") or -1)
        != PREVIOUS_FILE_SIZE
    ):
        raise RuntimeError("accepted manifest lacks the prior source observation")

    files = _processing_files()
    located = 0
    for item in files:
        if item["path"] == CHANGED_PATH:
            item["sha256"] = PREVIOUS_FILE_SHA256
            item["size_bytes"] = PREVIOUS_FILE_SIZE
            located += 1
    if located != 1:
        raise RuntimeError("could not isolate the changed processing-identity file")
    reconstructed = semantic_sha256(files)
    accepted_identity = dict(manifest.get("processing_identity") or {})
    if (
        accepted_identity.get("source_files_sha256")
        != ACCEPTED_SOURCE_FILES_SHA256
        or reconstructed != ACCEPTED_SOURCE_FILES_SHA256
    ):
        raise RuntimeError("summary-only reverse patch did not recover accepted identity")

    normalization_entrypoints = [
        PROJECT_ROOT / "matching/library_record_store.py",
        PROJECT_ROOT / "matching/search.py",
        PROJECT_ROOT / "mssystem/data_loader.py",
    ]
    forbidden_imports = (
        "matching.untargeted_workflow",
        "from .untargeted_workflow",
        "import untargeted_workflow",
    )
    if any(
        marker in path.read_text(encoding="utf-8")
        for path in normalization_entrypoints
        for marker in forbidden_imports
    ):
        raise RuntimeError("library normalization now imports the changed module")

    return {
        "changed_file": CHANGED_PATH,
        "previous_file_sha256": PREVIOUS_FILE_SHA256,
        "previous_file_size_bytes": PREVIOUS_FILE_SIZE,
        "current_file_sha256": CURRENT_FILE_SHA256,
        "current_file_size_bytes": CURRENT_FILE_SIZE,
        "reverse_patch_reproduced_previous_file": True,
        "reconstructed_accepted_source_files_sha256": reconstructed,
        "other_processing_identity_files_changed": 0,
        "library_normalizer_imports_changed_module": False,
        "changed_fields": [
            "consensus_spectra",
            "multi_spectrum_consensus",
            "single_spectrum_consensus",
            "multi_spectrum_representative_groups",
            "single_spectrum_representatives",
        ],
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
    proof = _proof(manifest)
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
        raise RuntimeError("numeric runtime or normalization options changed")
    if sql_config.get("processing_identity") != accepted_identity:
        raise RuntimeError("SQL metadata does not match the accepted identity")
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
        "schema_version": "mssystem-non-normalization-source-rebind/3.0",
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
            "correct representative-spectrum summary labels without changing "
            "library normalization, retrieval, scoring, or cached records"
        ),
    }

    field = "post_acceptance_representative_summary_rebind"
    new_manifest = dict(manifest)
    new_manifest["processing_identity"] = current_identity
    new_manifest[field] = audit
    _refresh_semantic(new_manifest)
    new_acceptance = dict(acceptance)
    new_acceptance["manifest_semantic_sha256"] = new_manifest["semantic_sha256"]
    new_acceptance[field] = audit
    _refresh_semantic(new_acceptance)
    new_record_config = dict(record_config)
    new_record_config["accepted_cache_manifest_semantic_sha256"] = (
        new_manifest["semantic_sha256"]
    )
    new_record_config[field] = audit
    new_sql_config = dict(sql_config)
    new_sql_config["processing_identity"] = current_identity
    new_sql_config["cache_manifest_semantic_sha256"] = new_manifest[
        "semantic_sha256"
    ]
    new_sql_config[field] = audit
    _refresh_semantic(new_sql_config)

    result = {
        "schema_version": "mssystem-representative-summary-source-rebind/1.0",
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

    backup_dir = output_path.parent / "pre_representative_summary_rebind"
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
