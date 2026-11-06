"""Fail-closed compatibility for accepted spectral-library cache identities.

The record-store identity historically covered every module under ``matching``.
That conservative scope also captures new, isolated workflows which are not
reachable from spectral-library normalization.  This module recognizes only
reviewed additions whose exact bytes can be removed from the current snapshot
to reproduce the accepted source fingerprint.  Any other source or option
change remains incompatible.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any


_PROCESSING_PREFIXES = ("matching/", "utils/", "mssystem/")
_PROCESSING_FILES = {"config/settings.py", "config/retrieval_policy.json"}

# These modules were added after cache acceptance.  They consume feature-table
# inputs in an explicit MassCube replay workflow and are not imported by any
# other processing-identity source.  Exact hashes make this an allowlist of
# reviewed bytes, rather than a filename exclusion rule.
_ISOLATED_ADDITIONS = {
    "matching/masscube_peak_evidence.py": {
        "size_bytes": 16040,
        "sha256": "c8a7e5350b6dc10b4883551fc75b877fe4a3c413e9951521e85cb3703f18ff90",
        "import_token": b"masscube_peak_evidence",
    },
    "matching/masscube_peak_workflow.py": {
        "size_bytes": 10724,
        "sha256": "a4a774c7609ab9411ed132c1d8d304682b4f79b82d45eac8af9c3ae9397c6034",
        "import_token": b"masscube_peak_workflow",
    },
}


# Network preparation consumes already-normalized query results. This exact
# reviewed replacement changes neither library loading nor peak normalization.
# Both versions are bound; arbitrary edits still fail the accepted fingerprint.
_NON_NORMALIZING_REPLACEMENTS = {
    "matching/evidence_network.py": {
        "before": {"size_bytes": 45670,
                   "sha256": "2fde1efe4440f2619053adb96dafd39168182dd11a4a568c99d2c95cffe9492d"},
        "after": {"size_bytes": 48216,
                  "sha256": "85861050f03c36ad12dc328715c9f52e995a201d4ca804383b965025969cbfcd"},
    },
    "matching/molecular_network.py": {
        "before": {"size_bytes": 29814,
                   "sha256": "ba1ed046c89b41bc0768024feabc713e7be0ae60f6d8440866699932abc4d0b4"},
        "after": {"size_bytes": 33002,
                  "sha256": "fdd7017c39396562dbd37773c1eaa94d748d1174cbbe813bdb2385e179d6579c"},
    },
}


def _semantic_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _processing_files(project_root: Path) -> list[dict[str, Any]]:
    from runtime.manifest import source_tree_snapshot

    return [
        item
        for item in source_tree_snapshot(project_root)["files"]
        if item["path"].startswith(_PROCESSING_PREFIXES)
        or item["path"] in _PROCESSING_FILES
    ]


def processing_identity_compatibility(
    accepted: Any,
    current: Any,
    *,
    project_root: Path,
) -> dict[str, Any]:
    """Resolve one current identity against an immutable accepted identity.

    Exact equality always passes.  The sole non-exact case requires identical
    runtime/options, an exact current source snapshot, exact reviewed module
    bytes, no imports of those modules from retained processing sources, and a
    reconstructed fingerprint equal to the accepted fingerprint.
    """
    if accepted == current:
        return {"compatible": True, "mode": "exact"}
    if not isinstance(accepted, dict) or not isinstance(current, dict):
        return {"compatible": False, "reason": "identity_not_mapping"}
    accepted_source = accepted.get("source_files_sha256")
    current_source = current.get("source_files_sha256")
    if not isinstance(accepted_source, str) or not isinstance(current_source, str):
        return {"compatible": False, "reason": "source_fingerprint_missing"}
    accepted_options = {key: value for key, value in accepted.items() if key != "source_files_sha256"}
    current_options = {key: value for key, value in current.items() if key != "source_files_sha256"}
    if accepted_options != current_options:
        return {"compatible": False, "reason": "processing_options_changed"}

    files = _processing_files(Path(project_root).resolve())
    if _semantic_sha256(files) != current_source:
        return {"compatible": False, "reason": "current_source_fingerprint_not_reproducible"}
    by_path = {item["path"]: item for item in files}
    for relative, binding in _ISOLATED_ADDITIONS.items():
        record = by_path.get(relative)
        if record is None or record.get("size_bytes") != binding["size_bytes"] or record.get("sha256") != binding["sha256"]:
            return {"compatible": False, "reason": "reviewed_isolated_addition_changed", "path": relative}

    retained = [item for item in files if item["path"] not in _ISOLATED_ADDITIONS]
    tokens = tuple(binding["import_token"] for binding in _ISOLATED_ADDITIONS.values())
    root = Path(project_root).resolve()
    for item in retained:
        if not item["path"].endswith(".py"):
            continue
        try:
            content = (root / item["path"]).read_bytes()
        except OSError:
            return {"compatible": False, "reason": "retained_source_unreadable", "path": item["path"]}
        if any(token in content for token in tokens):
            return {"compatible": False, "reason": "isolated_addition_imported_by_processing_source", "path": item["path"]}
    reconstructed_records = []
    replacements = []
    for item in retained:
        replacement = _NON_NORMALIZING_REPLACEMENTS.get(item["path"])
        if replacement is not None and all(
            item.get(field) == replacement["after"][field]
            for field in ("size_bytes", "sha256")
        ):
            reconstructed_records.append({**item, **replacement["before"]})
            replacements.append({"path": item["path"], **replacement})
        else:
            reconstructed_records.append(item)
    reconstructed = _semantic_sha256(reconstructed_records)
    if reconstructed != accepted_source:
        return {
            "compatible": False,
            "reason": "unreviewed_processing_source_change",
            "reconstructed_source_files_sha256": reconstructed,
        }
    return {
        "compatible": True,
        "mode": "verified_isolated_non_normalizing_additions",
        "accepted_source_files_sha256": accepted_source,
        "current_source_files_sha256": current_source,
        "reconstructed_source_files_sha256": reconstructed,
        "excluded_additions": [
            {"path": path, "size_bytes": binding["size_bytes"], "sha256": binding["sha256"]}
            for path, binding in sorted(_ISOLATED_ADDITIONS.items())
        ],
        "normalization_or_scoring_inputs_changed": False,
        "reviewed_non_normalizing_replacements": replacements,
    }


def load_compatible_record_store(
    root: Path,
    *,
    source_db: Path,
    current_identity: dict[str, Any],
    project_root: Path,
    cache_bytes: int,
):
    """Load through the original strict loader after an exact compatibility proof."""
    from matching.library_record_store import load_record_store

    cache = Path(root).expanduser().resolve()
    try:
        manifest = json.loads((cache / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError) as exc:
        raise ValueError("Record cache manifest is unreadable") from exc
    accepted = manifest.get("processing_identity")
    resolution = processing_identity_compatibility(
        accepted, current_identity, project_root=project_root,
    )
    if resolution.get("compatible") is not True:
        raise ValueError(
            "Record cache processing identity mismatch: "
            + str(resolution.get("reason") or "unverified source change")
        )
    # The strict loader still checks the immutable manifest, acceptance receipt,
    # source database binding, and every payload hash against accepted_identity.
    return load_record_store(
        cache,
        source_db=source_db,
        identity=accepted,
        cache_bytes=cache_bytes,
    )


def load_runtime_compatible_record_store(
    db_path: Path,
    *,
    project_root: Path,
    cache_bytes: int | None = None,
):
    """Apply the configured runtime switch and load only a proven compatible cache."""
    from matching.library_record_store import DEFAULT_CACHE_BYTES, processing_identity

    root = Path(project_root).resolve()
    config_path = root / "config/library_record_store.json"
    config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.is_file() else {}
    override = os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "").strip()
    if override.lower() in {"off", "0", "false"}:
        return None
    if not override and config.get("enabled") is not True:
        return None
    bound_db = str(config.get("source_db") or "").strip()
    requested = Path(db_path).resolve()
    if not override and (not bound_db or Path(bound_db).resolve() != requested):
        return None
    cache = Path(override or config["cache_dir"]).expanduser().resolve()
    return load_compatible_record_store(
        cache,
        source_db=requested,
        current_identity=processing_identity(),
        project_root=root,
        cache_bytes=DEFAULT_CACHE_BYTES if cache_bytes is None else cache_bytes,
    )


__all__ = [
    "load_compatible_record_store",
    "load_runtime_compatible_record_store",
    "processing_identity_compatibility",
]
