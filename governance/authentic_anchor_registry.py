"""Fail-closed binding for the active authentic-standard MS/MS library.

The authentic-standard library is a local runtime artifact rather than source
code.  A registry entry is therefore not sufficient evidence that the lane is
deployable: the declared file must exist, validate as an AuthenticAnchor-MS
artifact, and match both its semantic content hash and (when declared) its raw
file hash.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Mapping

from matching.authentic_anchor import read_anchor_library


REGISTRY_SCHEMA_VERSION = "authentic-anchor-registry/1.0"
ACTIVE_REGISTRY_STATE = "active_authentic_reference_display_ranker"


def _sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def load_authentic_anchor_registry(path: str | Path) -> dict[str, Any]:
    registry_path = Path(path)
    payload = json.loads(registry_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("authentic-anchor registry must be a JSON object")
    if payload.get("schema_version") != REGISTRY_SCHEMA_VERSION:
        raise ValueError("unsupported authentic-anchor registry schema")
    if payload.get("state") != ACTIVE_REGISTRY_STATE:
        raise ValueError("authentic-anchor registry is not active")
    return payload


def resolve_active_library_path(
    registry: Mapping[str, Any],
    *,
    project_root: str | Path,
) -> Path:
    root = Path(project_root).resolve()
    declared_text = str(registry.get("active_library") or "").strip()
    if not declared_text:
        raise ValueError("active authentic-reference library path is missing")
    declared = Path(declared_text)
    resolved = declared.resolve() if declared.is_absolute() else (root / declared).resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError("active authentic-reference library is outside the project")
    if resolved == root:
        raise ValueError("active authentic-reference library cannot be the project root")
    return resolved


def audit_active_authentic_anchor_registry(
    registry_path: str | Path,
    *,
    project_root: str | Path,
) -> dict[str, Any]:
    """Return a JSON-safe, fail-closed deployment audit."""

    registry_file = Path(registry_path).resolve()
    report: dict[str, Any] = {
        "schema_version": "authentic-anchor-deployment-audit/1.0",
        "registry_path": str(registry_file),
        "registry_exists": registry_file.is_file(),
        "passed": False,
        "violations": [],
    }
    if not registry_file.is_file():
        report["violations"] = ["registry_missing"]
        return report
    try:
        registry = load_authentic_anchor_registry(registry_file)
        library_path = resolve_active_library_path(
            registry,
            project_root=project_root,
        )
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        report["violations"] = ["registry_invalid"]
        report["error"] = str(exc)
        return report

    report.update(
        {
            "registry_state": registry.get("state"),
            "library_path": str(library_path),
            "library_exists": library_path.is_file(),
            "expected_content_sha256": str(
                registry.get("active_library_content_sha256") or ""
            ),
            "expected_file_sha256": str(
                registry.get("active_library_file_sha256") or ""
            ),
            "instrument_family": registry.get("instrument_family"),
            "formal_level1_release_enabled": bool(
                registry.get("formal_level1_release_enabled", False)
            ),
        }
    )
    if not library_path.is_file():
        report["violations"] = ["active_library_missing"]
        return report

    try:
        actual_file_hash = _sha256_file(library_path)
        library = read_anchor_library(library_path)
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        report["violations"] = ["active_library_invalid"]
        report["error"] = str(exc)
        return report

    expected_content_hash = report["expected_content_sha256"]
    expected_file_hash = report["expected_file_sha256"]
    actual_content_hash = str(library.get("content_sha256") or "")
    violations: list[str] = []
    if not expected_content_hash:
        violations.append("registry_content_hash_missing")
    elif actual_content_hash != expected_content_hash:
        violations.append("active_library_content_hash_mismatch")
    if expected_file_hash and actual_file_hash != expected_file_hash:
        violations.append("active_library_file_hash_mismatch")
    entries = library.get("entries") or []
    if not isinstance(entries, list) or not entries:
        violations.append("active_library_empty")

    report.update(
        {
            "actual_content_sha256": actual_content_hash,
            "actual_file_sha256": actual_file_hash,
            "entry_count": len(entries) if isinstance(entries, list) else 0,
            "violations": violations,
            "passed": not violations,
        }
    )
    return report


def load_active_authentic_anchor_library(
    registry_path: str | Path,
    *,
    project_root: str | Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load the active library only after the deployment audit passes."""

    audit = audit_active_authentic_anchor_registry(
        registry_path,
        project_root=project_root,
    )
    if not audit.get("passed"):
        details = ",".join(audit.get("violations") or ()) or "unknown_failure"
        raise ValueError(f"active authentic-reference deployment failed: {details}")
    registry = load_authentic_anchor_registry(registry_path)
    library = read_anchor_library(Path(str(audit["library_path"])))
    return library, registry


def install_active_authentic_anchor_library(
    source_path: str | Path,
    registry_path: str | Path,
    *,
    project_root: str | Path,
) -> dict[str, Any]:
    """Atomically install one explicitly supplied, registry-bound artifact."""

    source = Path(source_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(f"authentic-reference source is missing: {source}")
    registry = load_authentic_anchor_registry(registry_path)
    destination = resolve_active_library_path(registry, project_root=project_root)
    source_file_hash_before = _sha256_file(source)
    source_library = read_anchor_library(source)
    expected_content_hash = str(registry.get("active_library_content_sha256") or "")
    expected_file_hash = str(registry.get("active_library_file_sha256") or "")
    if not expected_content_hash or source_library.get("content_sha256") != expected_content_hash:
        raise ValueError("source authentic-reference content hash does not match registry")
    if expected_file_hash and source_file_hash_before != expected_file_hash:
        raise ValueError("source authentic-reference file hash does not match registry")

    destination.parent.mkdir(parents=True, exist_ok=True)
    source_bytes = source.read_bytes()
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(source_bytes)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise

    source_file_hash_after = _sha256_file(source)
    if source_file_hash_after != source_file_hash_before:
        raise RuntimeError("source authentic-reference artifact changed during installation")
    audit = audit_active_authentic_anchor_registry(
        registry_path,
        project_root=project_root,
    )
    if not audit.get("passed"):
        raise RuntimeError(
            "installed authentic-reference artifact failed deployment audit: "
            + ",".join(audit.get("violations") or ())
        )
    return {
        "schema_version": "authentic-anchor-install-receipt/1.0",
        "source_path": str(source),
        "destination_path": str(destination),
        "source_file_sha256_before": source_file_hash_before,
        "source_file_sha256_after": source_file_hash_after,
        "source_mutated": False,
        "content_sha256": source_library.get("content_sha256"),
        "entry_count": len(source_library.get("entries") or []),
        "deployment_audit": audit,
        "installed": True,
    }


__all__ = [
    "ACTIVE_REGISTRY_STATE",
    "REGISTRY_SCHEMA_VERSION",
    "audit_active_authentic_anchor_registry",
    "install_active_authentic_anchor_library",
    "load_active_authentic_anchor_library",
    "load_authentic_anchor_registry",
    "resolve_active_library_path",
]
