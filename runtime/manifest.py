"""Build a machine-verifiable description of the running MSsystem V16 source state.

The V16 worktree contains intentionally preserved V15/V14/V12/V11/V10/V9-derived
uncommitted work.  A Git commit alone therefore cannot identify the executable
source.  The runtime manifest records both Git state and a content hash over the
actual production source tree, while large databases and model weights are
registered separately by immutable fingerprints.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import sqlite3
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MANIFEST_SCHEMA_VERSION = "1.0.0"
DEFAULT_RUNTIME_PACKAGES = (
    "fastapi",
    "uvicorn",
    "pydantic",
    "numpy",
    "scipy",
    "pandas",
    "scikit-learn",
    "rdkit",
    "matchms",
    "pyopenms",
    "ms-entropy",
    "lightgbm",
    "joblib",
    "beautifulsoup4",
    "lxml",
    "openpyxl",
    "reportlab",
    "ijson",
)
SOURCE_DIRECTORIES = (
    "acquisition",
    "calibration",
    "classifier",
    "config",
    "governance",
    "kb_mining",
    "kb_v3",
    "matching",
    "mssystem",
    "pipeline",
    "runtime",
    "schemas",
    "scripts",
    "utils",
    "web",
)
SOURCE_FILES = (
    "run.sh",
    "requirements.txt",
    "requirements.lock",
    "environment.yml",
)
SOURCE_SUFFIXES = {
    ".css",
    ".html",
    ".js",
    ".json",
    ".jsonl",
    ".py",
    ".sh",
    ".txt",
    ".yaml",
    ".yml",
}
SELF_REFERENTIAL_PATHS = {
    "config/golden_baseline.json",
    "config/runtime_manifest.json",
}
RUNTIME_IDENTITY_TABLES = (
    "Spectra",
    "Compounds",
    "TargetCompounds",
    "DiagnosticIons",
)
CONTENT_HASH_TABLES = ("TargetCompounds", "DiagnosticIons")


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_file(path: Path, *, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _iter_source_files(project_root: Path) -> Iterable[Path]:
    candidates: list[Path] = []
    for relative in SOURCE_FILES:
        path = project_root / relative
        if path.is_file():
            candidates.append(path)
    for relative in SOURCE_DIRECTORIES:
        root = project_root / relative
        if not root.is_dir():
            continue
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in SOURCE_SUFFIXES:
                continue
            relative_path = path.relative_to(project_root)
            if "__pycache__" in relative_path.parts or str(relative_path) in SELF_REFERENTIAL_PATHS:
                continue
            candidates.append(path)
    yield from sorted(set(candidates), key=lambda path: str(path.relative_to(project_root)))


def source_tree_snapshot(project_root: Path) -> dict[str, Any]:
    files: list[dict[str, Any]] = []
    tree_digest = hashlib.sha256()
    for path in _iter_source_files(project_root):
        relative = path.relative_to(project_root).as_posix()
        stat = path.stat()
        digest = sha256_file(path)
        record = {"path": relative, "size_bytes": stat.st_size, "sha256": digest}
        files.append(record)
        tree_digest.update(_canonical_json(record).encode("utf-8"))
        tree_digest.update(b"\n")
    return {
        "method": "sha256(canonical path,size,content-sha256 records)",
        "file_count": len(files),
        "source_tree_sha256": tree_digest.hexdigest(),
        "files": files,
    }


def source_tree_stat_signature(project_root: Path) -> str:
    """Cheaply identify source-tree changes before recomputing content hashes."""

    digest = hashlib.sha256()
    for path in _iter_source_files(project_root):
        relative = path.relative_to(project_root).as_posix()
        try:
            stat = path.stat()
        except OSError:
            continue
        digest.update(
            _canonical_json(
                {
                    "path": relative,
                    "size_bytes": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            ).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def runtime_manifest_dependency_signature(project_root: Path) -> str:
    """Identify source changes plus the DB stat bound by screening coverage."""

    digest = hashlib.sha256(source_tree_stat_signature(project_root).encode("ascii"))
    coverage_path = project_root / "config" / "screening_catalog_coverage.json"
    if coverage_path.is_file():
        try:
            payload = _load_json(coverage_path)
            fingerprint = payload.get("database_fingerprint") or {}
            database_path = Path(str(fingerprint.get("path") or ""))
            if database_path.is_file():
                stat = database_path.stat()
                digest.update(
                    _canonical_json(
                        {
                            "path": str(database_path.resolve()),
                            "size_bytes": stat.st_size,
                            "mtime_ns": stat.st_mtime_ns,
                        }
                    ).encode("utf-8")
                )
        except (OSError, TypeError, ValueError):
            digest.update(b"invalid-screening-coverage-artifact")
    authentic_registry_path = project_root / "config" / "authentic_anchor_registry.json"
    for path in (authentic_registry_path,):
        if path.is_file():
            stat = path.stat()
            digest.update(
                _canonical_json(
                    {
                        "path": str(path.resolve()),
                        "size_bytes": stat.st_size,
                        "mtime_ns": stat.st_mtime_ns,
                    }
                ).encode("utf-8")
            )
    if authentic_registry_path.is_file():
        try:
            from governance.authentic_anchor_registry import (
                load_authentic_anchor_registry,
                resolve_active_library_path,
            )

            registry = load_authentic_anchor_registry(authentic_registry_path)
            authentic_path = resolve_active_library_path(
                registry,
                project_root=project_root,
            )
            if authentic_path.is_file():
                stat = authentic_path.stat()
                digest.update(
                    _canonical_json(
                        {
                            "path": str(authentic_path),
                            "size_bytes": stat.st_size,
                            "mtime_ns": stat.st_mtime_ns,
                        }
                    ).encode("utf-8")
                )
            else:
                digest.update(b"active-authentic-reference-missing")
        except (OSError, TypeError, ValueError):
            digest.update(b"invalid-authentic-reference-registry")
    cache_config = project_root / "config" / "library_record_store.json"
    override = os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "").strip()
    digest.update(override.encode("utf-8"))
    try:
        config = _load_json(cache_config) if cache_config.is_file() else {}
        if (override.lower() not in {"off", "0", "false"}
                and (override or config.get("enabled") is True)):
            cache = Path(override or config["cache_dir"]).expanduser().resolve()
            for name in ("manifest.json", "acceptance.json"):
                artifact = cache / name
                if artifact.is_file():
                    stat = artifact.stat()
                    digest.update(_canonical_json({"path": str(artifact),
                        "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}).encode("utf-8"))
                else:
                    digest.update(f"missing-record-store-{name}".encode("utf-8"))
    except (OSError, TypeError, ValueError, KeyError):
        digest.update(b"invalid-record-store-binding")
    repair_registry = project_root / "config" / "library_metadata_repairs.json"
    if repair_registry.is_file():
        try:
            for source in _load_json(repair_registry)["source_files"]:
                stat = Path(source).stat()
                digest.update(_canonical_json({"path": source, "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns, "ctime_ns": stat.st_ctime_ns}).encode("utf-8"))
        except (OSError, TypeError, ValueError, KeyError):
            digest.update(b"invalid-metadata-repair-source")
    from runtime.verified_sql_retrieval import sql_metadata_dependency_signature
    digest.update(_canonical_json(sql_metadata_dependency_signature(project_root)).encode("utf-8"))
    return digest.hexdigest()


def _run_git(project_root: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=project_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def git_snapshot(project_root: Path) -> dict[str, Any]:
    try:
        commit = _run_git(project_root, "rev-parse", "HEAD")
        branch = _run_git(project_root, "branch", "--show-current")
        status_text = _run_git(
            project_root,
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        return {"available": False, "error": str(exc)}

    lines = [line for line in status_text.splitlines() if line]
    production_roots = tuple(f"{name}/" for name in SOURCE_DIRECTORIES)
    relevant = []
    for line in lines:
        path = line[3:].split(" -> ")[-1]
        if path in SOURCE_FILES or path.startswith(production_roots):
            relevant.append(line)
    status_digest = hashlib.sha256(status_text.encode("utf-8")).hexdigest()
    return {
        "available": True,
        "commit": commit,
        "branch": branch,
        "dirty": bool(lines),
        "dirty_entry_count": len(lines),
        "runtime_relevant_dirty_entry_count": len(relevant),
        "runtime_relevant_dirty_entries": relevant,
        "porcelain_sha256": status_digest,
        "identification_note": (
            "The Git commit is insufficient when dirty=true; source_tree_sha256 "
            "identifies the executable source content."
        ),
    }


def environment_snapshot(package_names: Sequence[str] = DEFAULT_RUNTIME_PACKAGES) -> dict[str, Any]:
    packages: dict[str, str | None] = {}
    for name in package_names:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "packages": packages,
    }


def _sampled_file_fingerprint(path: Path, *, sample_size: int = 1024 * 1024) -> dict[str, Any]:
    size = path.stat().st_size
    offsets = sorted({0, max(0, (size - sample_size) // 2), max(0, size - sample_size)})
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            block = handle.read(sample_size)
            digest.update(str(offset).encode("ascii"))
            digest.update(b"\0")
            digest.update(block)
    return {
        "method": "sha256(offset + first/middle/last up-to-1MiB blocks)",
        "sample_offsets": offsets,
        "sampled_sha256": digest.hexdigest(),
    }


def sqlite_snapshot(path: Path, *, count_table: str | None = "spectra") -> dict[str, Any]:
    if not path.exists():
        return {"path": str(path), "exists": False}
    stat = path.stat()
    record: dict[str, Any] = {
        "path": str(path.resolve()),
        "exists": True,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "fingerprint": _sampled_file_fingerprint(path),
    }
    try:
        uri = f"file:{path.resolve()}?mode=ro&immutable=1"
        with sqlite3.connect(uri, uri=True, timeout=10) as connection:
            geometry = {
                "page_count": int(connection.execute("PRAGMA page_count").fetchone()[0]),
                "page_size": int(connection.execute("PRAGMA page_size").fetchone()[0]),
                "schema_version": int(connection.execute("PRAGMA schema_version").fetchone()[0]),
                "user_version": int(connection.execute("PRAGMA user_version").fetchone()[0]),
            }
            record["sqlite"] = geometry
            tables = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            table_by_lower = {name.lower(): name for name in tables}
            if count_table:
                if count_table.lower() in table_by_lower:
                    actual_table = table_by_lower[count_table.lower()]
                    quoted = actual_table.replace('"', '""')
                    record["row_count"] = int(
                        connection.execute(f'SELECT count(*) FROM "{quoted}"').fetchone()[0]
                    )
            tracked_counts: dict[str, int] = {}
            tracked_content: dict[str, str] = {}
            for requested in RUNTIME_IDENTITY_TABLES:
                actual = table_by_lower.get(requested.lower())
                if not actual:
                    continue
                quoted_table = actual.replace('"', '""')
                tracked_counts[actual] = int(
                    connection.execute(
                        f'SELECT count(*) FROM "{quoted_table}"'
                    ).fetchone()[0]
                )
                if requested not in CONTENT_HASH_TABLES:
                    continue
                columns = [
                    str(row[1])
                    for row in connection.execute(
                        f'PRAGMA table_info("{quoted_table}")'
                    )
                ]
                digest = hashlib.sha256()
                if columns:
                    order_by = ", ".join(
                        f'"{column.replace(chr(34), chr(34) * 2)}"'
                        for column in columns
                    )
                    for row in connection.execute(
                        f'SELECT * FROM "{quoted_table}" ORDER BY {order_by}'
                    ):
                        digest.update(
                            _canonical_json(list(row)).encode("utf-8")
                        )
                        digest.update(b"\n")
                tracked_content[actual] = digest.hexdigest()
            record["tracked_table_row_counts"] = tracked_counts
            record["tracked_table_content_sha256"] = tracked_content
    except sqlite3.Error as exc:
        record["sqlite_error"] = str(exc)

    identity_payload = {
        "size_bytes": record["size_bytes"],
        "mtime_ns": record["mtime_ns"],
        "fingerprint": record["fingerprint"],
        "sqlite": record.get("sqlite"),
        "row_count": record.get("row_count"),
        "tracked_table_row_counts": record.get("tracked_table_row_counts"),
        "tracked_table_content_sha256": record.get(
            "tracked_table_content_sha256"
        ),
    }
    record["snapshot_id"] = f"sqlite-{hashlib.sha256(_canonical_json(identity_payload).encode()).hexdigest()[:24]}"
    record["snapshot_id_scope"] = (
        "Operational identity, not a full-file cryptographic checksum; combines SQLite geometry, "
        "file stat, tracked table counts/content and sampled file content."
    )
    return record


def library_snapshot_freshness(recorded: Mapping[str, Any]) -> dict[str, Any]:
    """Compare a frozen SQLite snapshot with the current read-only library."""

    path_text = str(recorded.get("path") or "").strip()
    expected = str(recorded.get("snapshot_id") or "").strip()
    if not path_text or not expected:
        return {
            "status": "unbound",
            "checked": False,
            "matches": True,
            "reason": "manifest_has_no_comparable_library_snapshot",
        }
    path = Path(path_text)
    current = sqlite_snapshot(path)
    matches = bool(current.get("exists")) and current.get("snapshot_id") == expected
    return {
        "status": "current" if matches else "stale",
        "checked": True,
        "matches": matches,
        "path": str(path.resolve()),
        "manifest_snapshot_id": expected,
        "current_snapshot_id": current.get("snapshot_id"),
        "manifest_table_row_counts": recorded.get("tracked_table_row_counts") or {},
        "current_table_row_counts": current.get("tracked_table_row_counts") or {},
        "manifest_table_content_sha256": recorded.get(
            "tracked_table_content_sha256"
        ) or {},
        "current_table_content_sha256": current.get(
            "tracked_table_content_sha256"
        ) or {},
    }


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def role_catalog_snapshot(project_root: Path) -> dict[str, Any]:
    """Fingerprint bundled role catalogs and expose their release cardinality."""

    role_path = project_root / "matching" / "compound_annotations.json"
    screening_path = project_root / "matching" / "screening_catalog.json"
    coverage_path = project_root / "config" / "screening_catalog_coverage.json"

    compound_roles: dict[str, Any] = {
        "path": str(role_path.resolve()),
        "exists": role_path.is_file(),
    }
    if role_path.is_file():
        payload = _load_json(role_path)
        annotations = payload.get("annotations") or {}
        role_counts: dict[str, int] = {}
        class_level = 0
        for annotation in annotations.values():
            if not isinstance(annotation, Mapping):
                continue
            role = str(annotation.get("compound_role") or "unclassified")
            role_counts[role] = role_counts.get(role, 0) + 1
            if annotation.get("role_annotation_scope") == "class_level_structure_rule":
                class_level += 1
        compound_roles.update(
            {
                "sha256": sha256_file(role_path),
                "catalog_id": payload.get("catalog_id"),
                "catalog_version": payload.get("catalog_version"),
                "annotation_count": len(annotations),
                "role_counts": role_counts,
                "external_contaminant_count": role_counts.get(
                    "external_contaminant", 0
                ),
                "class_level_structure_rule_count": class_level,
                "baseline_contract": payload.get("baseline_contract") or {},
            }
        )

    screening: dict[str, Any] = {
        "path": str(screening_path.resolve()),
        "exists": screening_path.is_file(),
    }
    if screening_path.is_file():
        payload = _load_json(screening_path)
        entries = payload.get("entries") or []
        indexed_keys: set[str] = set()
        role_counts: dict[str, int] = {}
        for entry in entries:
            if not isinstance(entry, Mapping):
                continue
            key = str(entry.get("full_inchikey") or "").strip().upper()
            if key:
                indexed_keys.add(key)
            indexed_keys.update(
                str(value or "").strip().upper()
                for value in entry.get("alternate_full_inchikeys") or []
                if str(value or "").strip()
            )
            role = str(entry.get("compound_role") or "unclassified")
            role_counts[role] = role_counts.get(role, 0) + 1
        screening.update(
            {
                "sha256": sha256_file(screening_path),
                "catalog_id": payload.get("catalog_id"),
                "released_at": payload.get("released_at"),
                "entry_count": len(entries),
                "indexed_full_inchikey_count": len(indexed_keys),
                "role_counts": role_counts,
                "regulatory_claims_encoded": False,
            }
        )
    coverage: dict[str, Any] = {
        "artifact_path": str(coverage_path.resolve()),
        "exists": coverage_path.is_file(),
        "status": "not_audited",
    }
    if coverage_path.is_file():
        try:
            payload = _load_json(coverage_path)
            fingerprint = payload.get("database_fingerprint") or {}
            database_path = Path(str(fingerprint.get("path") or ""))
            database_current = False
            if database_path.is_file():
                stat = database_path.stat()
                database_current = (
                    int(fingerprint.get("size_bytes") or -1) == stat.st_size
                    and int(fingerprint.get("mtime_ns") or -1) == stat.st_mtime_ns
                )
            catalog_current = payload.get("catalog_sha256") == screening.get("sha256")
            current = bool(database_current and catalog_current)
            coverage.update(
                {
                    "artifact_sha256": sha256_file(coverage_path),
                    "status": "current" if current else "stale",
                    "audit_generated_at": payload.get("audit_generated_at"),
                    "catalog_current": catalog_current,
                    "database_current": database_current,
                    "database_snapshot_id": payload.get("database_snapshot_id"),
                    "database_snapshot_binding": payload.get(
                        "database_snapshot_binding"
                    ),
                    "catalog_full_keys_with_searchable_spectra": payload.get(
                        "catalog_full_keys_with_searchable_spectra"
                    ),
                    "catalog_full_key_coverage_fraction": payload.get(
                        "catalog_full_key_coverage_fraction"
                    ),
                }
            )
        except (OSError, TypeError, ValueError) as exc:
            coverage.update({"status": "invalid", "error": str(exc)})
    screening["library_coverage"] = coverage
    return {"compound_roles": compound_roles, "screening": screening}


def _runtime_snapshot_id(manifest: Mapping[str, Any]) -> str:
    environment = manifest.get("environment") or {}
    library = manifest.get("library") or {}
    models = manifest.get("models") or {}
    source = manifest.get("source") or {}
    git = manifest.get("git") or {}
    production = library.get("production_full") or {}
    mapping = library.get("mapping") or {}
    identity_payload = {
        "git_commit": git.get("commit"),
        "source_tree_sha256": source.get("source_tree_sha256"),
        "python": (environment.get("python") or {}).get("version"),
        "packages": environment.get("packages") or {},
        "model_registry_sha256": models.get("registry_sha256"),
        "library_snapshot_id": production.get("snapshot_id"),
        "mapping_sha256": mapping.get("sha256"),
        "role_catalogs": manifest.get("role_catalogs") or {},
    }
    return "v16-" + hashlib.sha256(
        _canonical_json(identity_payload).encode("utf-8")
    ).hexdigest()[:24]


def model_snapshot(project_root: Path, registry_path: Path) -> dict[str, Any]:
    if not registry_path.exists():
        return {"registry": str(registry_path), "exists": False, "models": []}
    registry = _load_json(registry_path)
    models = []
    for entry in registry.get("models", []):
        resolved_artifacts = []
        for artifact in entry.get("artifacts", []):
            path = project_root / artifact["path"]
            resolved = dict(artifact)
            resolved.update({"exists": path.is_file()})
            if path.is_file():
                resolved.update({"size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
            resolved_artifacts.append(resolved)
        resolved_entry = dict(entry)
        resolved_entry["artifacts"] = resolved_artifacts
        models.append(resolved_entry)
    return {
        "registry": str(registry_path.resolve()),
        "registry_sha256": sha256_file(registry_path),
        "models": models,
    }


def fetch_health(url: str, *, timeout_seconds: float = 5.0) -> dict[str, Any]:
    try:
        with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
            payload = json.loads(response.read().decode("utf-8"))
        return {"url": url, "reachable": True, "payload": payload}
    except (OSError, ValueError, urllib.error.URLError) as exc:
        return {"url": url, "reachable": False, "error": str(exc)}


def fetch_api_schema(url: str, *, timeout_seconds: float = 5.0) -> dict[str, Any]:
    try:
        with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
            payload = json.loads(response.read().decode("utf-8"))
        paths = payload.get("paths", {}) if isinstance(payload, Mapping) else {}
        operations = 0
        for path_item in paths.values():
            if not isinstance(path_item, Mapping):
                continue
            operations += sum(
                1
                for method in path_item
                if str(method).lower()
                in {"get", "post", "put", "patch", "delete", "options", "head"}
            )
        canonical = _canonical_json(payload)
        return {
            "url": url,
            "reachable": True,
            "openapi_version": payload.get("openapi"),
            "api_title": (payload.get("info") or {}).get("title"),
            "api_version": (payload.get("info") or {}).get("version"),
            "path_count": len(paths),
            "operation_count": operations,
            "canonical_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        }
    except (OSError, ValueError, urllib.error.URLError) as exc:
        return {"url": url, "reachable": False, "error": str(exc)}


def build_runtime_manifest(
    project_root: Path,
    *,
    model_registry_path: Path | None = None,
    library_path: Path | None = None,
    mapping_path: Path | None = None,
    health_url: str | None = None,
    restore_archive: Path | None = None,
) -> dict[str, Any]:
    project_root = project_root.resolve()
    if model_registry_path is None:
        model_registry_path = project_root / "config" / "model_registry.json"
    if library_path is None or mapping_path is None:
        from config.settings import DB_PATH, LIBRARY_MAPPING_PATH

        library_path = library_path or Path(DB_PATH)
        mapping_path = mapping_path or Path(LIBRARY_MAPPING_PATH)

    mapping: dict[str, Any] = {"path": str(Path(mapping_path).resolve()), "exists": False}
    if Path(mapping_path).is_file():
        mapping.update(
            {
                "exists": True,
                "size_bytes": Path(mapping_path).stat().st_size,
                "sha256": sha256_file(Path(mapping_path)),
            }
        )

    source = source_tree_snapshot(project_root)
    from governance.authentic_anchor_registry import (
        audit_active_authentic_anchor_registry,
    )

    authentic_reference = audit_active_authentic_anchor_registry(
        project_root / "config" / "authentic_anchor_registry.json",
        project_root=project_root,
    )
    manifest: dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "project_root": str(project_root),
        "service_generation": "V16",
        "service_version": "16.0.0",
        "git": git_snapshot(project_root),
        "source": source,
        "environment": environment_snapshot(),
        "models": model_snapshot(project_root, model_registry_path),
        "role_catalogs": role_catalog_snapshot(project_root),
        "library": {
            "production_full": sqlite_snapshot(Path(library_path)),
            "mapping": mapping,
            "write_enabled": os.environ.get("MSSYSTEM_LIBRARY_WRITE_ENABLED", "0").strip().lower()
            in {"1", "true", "yes", "on"},
        },
        "authentic_reference": authentic_reference,
    }
    if health_url:
        manifest["live_service"] = fetch_health(health_url)
        if health_url.endswith("/api/health"):
            openapi_url = health_url[: -len("/api/health")] + "/openapi.json"
        else:
            openapi_url = health_url.rstrip("/") + "/openapi.json"
        manifest["api_schema"] = fetch_api_schema(openapi_url)
    if restore_archive is not None:
        archive = restore_archive.resolve()
        manifest["restore_archive"] = {
            "path": str(archive),
            "exists": archive.is_file(),
            "size_bytes": archive.stat().st_size if archive.is_file() else None,
            "sha256": sha256_file(archive) if archive.is_file() else None,
        }

    manifest["runtime_snapshot_id"] = _runtime_snapshot_id(manifest)
    return manifest


def refresh_runtime_manifest_source(
    project_root: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    """Refresh executable source identity without rescanning the shared library.

    Startup must identify the exact code and bundled catalogs it is about to
    execute.  The large production library fingerprint is preserved from the
    previously verified manifest and explicitly marked as not rescanned during
    this lightweight startup operation.
    """

    project_root = project_root.resolve()
    manifest_path = manifest_path.resolve()
    try:
        existing = _load_json(manifest_path)
    except (OSError, TypeError, ValueError):
        existing = {}
    manifest = dict(existing) if isinstance(existing, Mapping) else {}
    manifest.update(
        {
            "schema_version": MANIFEST_SCHEMA_VERSION,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "project_root": str(project_root),
            "service_generation": "V16",
            "service_version": "16.0.0",
            "git": git_snapshot(project_root),
            "source": source_tree_snapshot(project_root),
            "environment": environment_snapshot(),
            "role_catalogs": role_catalog_snapshot(project_root),
            "refresh_scope": "source_environment_and_bundled_role_catalogs",
        }
    )
    from governance.authentic_anchor_registry import (
        audit_active_authentic_anchor_registry,
    )

    manifest["authentic_reference"] = audit_active_authentic_anchor_registry(
        project_root / "config" / "authentic_anchor_registry.json",
        project_root=project_root,
    )
    # V16 may be created by a copy-on-write clone.  Rebind all lightweight,
    # checkout-local dependencies so a refreshed manifest never continues to
    # advertise the parent V10 checkout as the executable source of truth.
    model_registry_path = project_root / "config" / "model_registry.json"
    manifest["models"] = model_snapshot(project_root, model_registry_path)
    golden_path = project_root / "config" / "golden_baseline.json"
    golden_record: dict[str, Any] = {
        "path": str(golden_path.resolve()),
        "exists": golden_path.is_file(),
    }
    if golden_path.is_file():
        try:
            golden_payload = _load_json(golden_path)
            golden_record.update(
                {
                    "sha256": sha256_file(golden_path),
                    "schema_version": golden_payload.get("schema_version"),
                    "task_count": len(golden_payload.get("tasks") or {}),
                }
            )
        except (OSError, TypeError, ValueError) as exc:
            golden_record.update({"status": "invalid", "error": str(exc)})
    manifest["golden_baseline"] = golden_record
    if isinstance(manifest.get("golden_baseline_refresh"), Mapping):
        manifest["golden_baseline_refresh"] = {
            **dict(manifest["golden_baseline_refresh"]),
            "reason": (
                "V16 output/ does not contain the frozen V9-derived golden task artifacts"
            ),
        }
    if isinstance(manifest.get("screening_coverage_refresh"), Mapping):
        manifest["screening_coverage_refresh"] = {
            **dict(manifest["screening_coverage_refresh"]),
            "coverage_artifact": str(
                (project_root / "config" / "screening_catalog_coverage.json").resolve()
            ),
        }
    library = dict(manifest.get("library") or {})
    library["startup_verification"] = (
        "preserved_previous_read_only_snapshot; shared library not rescanned"
    )
    manifest["library"] = library
    manifest.pop("live_service", None)
    manifest.pop("api_schema", None)
    manifest["runtime_snapshot_id"] = _runtime_snapshot_id(manifest)
    write_json_atomic(manifest_path, manifest)
    return manifest


def library_record_store_compatibility(project_root: Path) -> dict[str, Any]:
    """Check accepted cache/code compatibility without loading spectral records.

    Source freshness alone does not prove the next library load can succeed.
    Full record checksums are still verified by the actual library loader.
    """
    from matching.library_record_store import (
        processing_identity, semantic_sha256, validate_acceptance,
    )
    from runtime.library_cache_identity import processing_identity_compatibility

    path = project_root / "config" / "library_record_store.json"
    override = os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "").strip()
    if override.lower() in {"off", "0", "false"}:
        return {"status": "disabled", "compatible": True, "required": False}
    try:
        config = _load_json(path) if path.is_file() else {}
        if not override and config.get("enabled") is not True:
            return {"status": "disabled", "compatible": True, "required": False}
        cache = Path(override or config["cache_dir"]).expanduser().resolve()
        manifest = _load_json(cache / "manifest.json")
        accepted_identity = manifest.get("processing_identity")
        current_identity = processing_identity()
        identity_resolution = processing_identity_compatibility(
            accepted_identity, current_identity, project_root=project_root,
        )
        if identity_resolution.get("compatible") is not True:
            return {
                "status": "processing_identity_mismatch", "compatible": False,
                "required": True, "cache_dir": str(cache),
                "accepted_processing_identity": accepted_identity,
                "current_processing_identity": current_identity,
                "processing_identity_resolution": identity_resolution,
            }
        body = {k: v for k, v in manifest.items() if k != "semantic_sha256"}
        if manifest.get("semantic_sha256") != semantic_sha256(body):
            raise ValueError("Cache manifest checksum mismatch")
        acceptance = _load_json(cache / "acceptance.json")
        validate_acceptance(acceptance, manifest)
        return {
            "status": "compatible", "compatible": True, "required": True,
            "cache_dir": str(cache), "record_count": manifest.get("record_count"),
            "processing_identity": current_identity,
            "accepted_processing_identity": accepted_identity,
            "processing_identity_resolution": identity_resolution,
            "full_record_integrity_verified_by": "library_load",
        }
    except (OSError, TypeError, ValueError, KeyError) as exc:
        return {"status": "missing_or_invalid", "compatible": False,
                "required": True, "error": str(exc)}


def runtime_manifest_freshness(
    project_root: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    """Compare the manifest with current executable source and role catalogs."""

    try:
        manifest = _load_json(manifest_path)
    except (OSError, TypeError, ValueError) as exc:
        return {"status": "missing_or_invalid", "fresh": False, "error": str(exc)}
    current_source = source_tree_snapshot(project_root)
    current_catalogs = role_catalog_snapshot(project_root)
    recorded_source = manifest.get("source") or {}
    recorded_catalogs = manifest.get("role_catalogs") or {}
    recorded_library = ((manifest.get("library") or {}).get("production_full") or {})
    catalog_matches = {
        name: bool((recorded_catalogs.get(name) or {}).get("sha256"))
        and (recorded_catalogs.get(name) or {}).get("sha256")
        == (current_catalogs.get(name) or {}).get("sha256")
        for name in ("compound_roles", "screening")
    }
    source_matches = recorded_source.get("source_tree_sha256") == current_source.get(
        "source_tree_sha256"
    )
    library_status = library_snapshot_freshness(recorded_library)
    library_matches = bool(library_status.get("matches"))
    fresh = bool(
        source_matches and all(catalog_matches.values()) and library_matches
    )
    from governance.authentic_anchor_registry import (
        audit_active_authentic_anchor_registry,
    )

    authentic_reference = audit_active_authentic_anchor_registry(
        project_root / "config" / "authentic_anchor_registry.json",
        project_root=project_root,
    )
    screening_coverage_current = (
        ((current_catalogs.get("screening") or {}).get("library_coverage") or {}).get(
            "status"
        )
        == "current"
    )
    record_store = library_record_store_compatibility(project_root)
    from runtime.library_metadata_repairs import metadata_repair_preflight
    metadata_repairs = metadata_repair_preflight(project_root)
    from runtime.verified_sql_retrieval import sql_metadata_repair_preflight
    sql_metadata_repairs = sql_metadata_repair_preflight(project_root)
    from runtime.source_compound_mass_repairs import source_compound_mass_preflight
    compound_mass_repairs = source_compound_mass_preflight(project_root)
    release_ready = bool(
        fresh and screening_coverage_current and authentic_reference.get("passed")
        and record_store.get("compatible")
        and metadata_repairs.get("compatible")
        and sql_metadata_repairs.get("compatible")
        and compound_mass_repairs.get("compatible")
    )
    return {
        "status": "fresh" if fresh else "stale",
        "fresh": fresh,
        "release_ready": release_ready,
        "screening_coverage_current": screening_coverage_current,
        "runtime_snapshot_id": manifest.get("runtime_snapshot_id"),
        "generated_at": manifest.get("generated_at"),
        "manifest_source_tree_sha256": recorded_source.get("source_tree_sha256"),
        "current_source_tree_sha256": current_source.get("source_tree_sha256"),
        "manifest_source_file_count": recorded_source.get("file_count"),
        "current_source_file_count": current_source.get("file_count"),
        "catalog_matches": catalog_matches,
        "library_matches": library_matches,
        "library": library_status,
        "library_record_store": record_store,
        "library_metadata_repairs": metadata_repairs,
        "sql_metadata_repairs": sql_metadata_repairs,
        "source_compound_mass_repairs": compound_mass_repairs,
        "authentic_reference": authentic_reference,
        "role_catalogs": current_catalogs,
    }


def write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
