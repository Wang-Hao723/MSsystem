"""Content-addressed storage and recovery state for vendor RAW conversion.

This module deliberately has no dependency on the web server or on the
converter implementation.  A caller can therefore add caching around a
conversion without importing ProteoWizard or PyOpenMS into the server process.

Cache entries are accepted only when all of the following still agree:

* the canonical cache-key material and its SHA-256 key;
* the cache and conversion validation status;
* the cached mzML size and SHA-256 digest.

Materialized task files never share an inode with the cache payload.  APFS
``clonefile`` is preferred because it provides copy-on-write independence.  A
hard link is used only as a short-lived source pin and is copied to a new inode
before the task path is published.
"""

from __future__ import annotations

import contextlib
import ctypes
import ctypes.util
import errno
import fcntl
import hashlib
import hmac
import json
import os
import re
import stat
import sys
import time
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterator, Mapping


CACHE_KEY_SCHEMA_VERSION = "mssystem-vendor-conversion-cache-key-v1"
CACHE_MANIFEST_SCHEMA_VERSION = "mssystem-vendor-conversion-cache-manifest-v1"
LEDGER_SCHEMA_VERSION = "mssystem-vendor-conversion-ledger-v1"

PEAK_PICKING_STRATEGIES = frozenset({"preserve", "vendor", "openms"})
LEDGER_STATES = frozenset({"queued", "converting", "validated", "registered", "failed"})
DEFAULT_LEDGER_LEASE_SECONDS = 5 * 60.0

_CACHE_KEY_PATTERN = re.compile(r"[0-9a-f]{64}")
_CACHE_PAYLOAD_NAME = "converted.mzML"
_CACHE_MANIFEST_NAME = "manifest.json"


class VendorConversionCacheError(RuntimeError):
    """Base error for cache or ledger operations."""


class SourceMaterialChangedError(VendorConversionCacheError):
    """The source no longer agrees with the signatures in ``raw_report``."""


class CacheEntryInvalidError(VendorConversionCacheError):
    """A cache entry or conversion manifest failed integrity validation."""


class LedgerStateError(VendorConversionCacheError):
    """A ledger is corrupt or a requested state transition is invalid."""


@dataclass(frozen=True)
class ConversionCacheKey:
    """Canonical key material and its content-addressed SHA-256 digest."""

    digest: str
    material: dict[str, Any]


@dataclass(frozen=True)
class CacheHit:
    """One fully verified cache entry."""

    key: str
    payload_path: Path
    manifest_path: Path
    manifest: dict[str, Any]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _json_clone(value: Any) -> Any:
    """Return a detached JSON value and reject non-persistable metadata."""

    return json.loads(_canonical_json(value))


def _sha256_file(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        try:
            os.fsync(descriptor)
        except OSError as exc:
            if exc.errno not in {errno.EINVAL, errno.ENOTSUP}:
                raise
    finally:
        os.close(descriptor)


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Durably publish one JSON document by replacing it on the same volume."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    )
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        _fsync_directory(destination.parent)
    finally:
        with contextlib.suppress(FileNotFoundError):
            temporary.unlink()


@contextlib.contextmanager
def _file_lock(path: Path, *, shared: bool = False) -> Iterator[None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as handle:
        operation = fcntl.LOCK_SH if shared else fcntl.LOCK_EX
        fcntl.flock(handle.fileno(), operation)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _normalized_signature(signature: Mapping[str, Any]) -> dict[str, Any]:
    fields = (
        "size_bytes",
        "mtime_ns",
        "sha256",
        "sha1_for_mzml_provenance",
    )
    normalized = {field: signature[field] for field in fields if field in signature}
    sha256 = str(normalized.get("sha256") or "")
    if not _CACHE_KEY_PATTERN.fullmatch(sha256):
        raise ValueError("vendor source signature must include a SHA-256 digest")
    for field in ("size_bytes", "mtime_ns"):
        value = normalized.get(field)
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"vendor source signature {field} must be non-negative")
    return normalized


def _normalized_relative_path(value: Any) -> str:
    if value is None:
        raise ValueError("relative vendor source path must not be null")
    normalized = str(value).replace("\\", "/")
    relative = PurePosixPath(normalized)
    if (
        not normalized
        or relative.is_absolute()
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise ValueError(f"invalid relative vendor source path: {value!r}")
    return relative.as_posix()


def _critical_signatures(raw_report: Mapping[str, Any]) -> dict[str, Any]:
    required = raw_report.get("required_file_signatures") or {}
    if not isinstance(required, Mapping):
        raise ValueError("raw_report.required_file_signatures must be a mapping")
    normalized_required: dict[str, Any] = {}
    for relative, signature in sorted(required.items(), key=lambda item: str(item[0])):
        if not isinstance(signature, Mapping):
            raise ValueError("each required vendor source signature must be a mapping")
        normalized_required[_normalized_relative_path(relative)] = (
            _normalized_signature(signature)
        )
    source_signature = raw_report.get("source_file_signature") or {}
    if not isinstance(source_signature, Mapping):
        raise ValueError("raw_report.source_file_signature must be a mapping")
    return {
        "required_files": normalized_required,
        "source_file": (
            _normalized_signature(source_signature) if source_signature else {}
        ),
    }


def _stat_record(path: Path, relative: str) -> dict[str, Any]:
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode):
        raise SourceMaterialChangedError(
            "Internal symlinks are not cacheable vendor source material: " f"{relative}"
        )
    if not stat.S_ISREG(before.st_mode):
        raise SourceMaterialChangedError(
            f"Non-regular vendor source entry is not cacheable: {relative}"
        )
    after = path.lstat()
    if (
        before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or before.st_ino != after.st_ino
    ):
        raise SourceMaterialChangedError(
            f"Vendor source changed while inventory was collected: {relative}"
        )
    return {
        "path": relative,
        "size_bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
    }


def _directory_material(root: Path) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for directory, directory_names, file_names in os.walk(
        root, topdown=True, followlinks=False
    ):
        directory_path = Path(directory)
        for name in tuple(directory_names):
            candidate = directory_path / name
            if candidate.is_symlink():
                relative = candidate.relative_to(root).as_posix()
                raise SourceMaterialChangedError(
                    "Internal directory symlinks are not cacheable vendor source "
                    f"material: {relative}"
                )
        for name in file_names:
            candidate = directory_path / name
            relative = candidate.relative_to(root).as_posix()
            records.append(_stat_record(candidate, relative))
    records.sort(key=lambda record: record["path"])
    return {
        "kind": "directory-size-mtime-inventory",
        "file_count": len(records),
        "total_bytes": sum(record["size_bytes"] for record in records),
        "inventory": records,
        "inventory_sha256": hashlib.sha256(
            _canonical_json(records).encode("utf-8")
        ).hexdigest(),
    }


def _file_material(path: Path) -> dict[str, Any]:
    record = _stat_record(path, path.name)
    return {
        "kind": "file-size-mtime-inventory",
        "file_count": 1,
        "total_bytes": record["size_bytes"],
        "inventory": [record],
        "inventory_sha256": hashlib.sha256(
            _canonical_json([record]).encode("utf-8")
        ).hexdigest(),
    }


def _validated_persisted_directory_material(value: Mapping[str, Any]) -> dict[str, Any]:
    persisted = _json_clone(value)
    inventory = persisted.get("inventory")
    if not isinstance(inventory, list):
        raise SourceMaterialChangedError(
            "persisted directory_material.inventory must be a list"
        )
    normalized_inventory: list[dict[str, Any]] = []
    seen_paths: set[str] = set()
    for record in inventory:
        if not isinstance(record, Mapping):
            raise SourceMaterialChangedError(
                "persisted directory material contains a malformed record"
            )
        try:
            relative = _normalized_relative_path(record.get("path"))
        except ValueError as exc:
            raise SourceMaterialChangedError(str(exc)) from exc
        size = record.get("size_bytes")
        mtime = record.get("mtime_ns")
        if (
            relative in seen_paths
            or not isinstance(size, int)
            or size < 0
            or not isinstance(mtime, int)
            or mtime < 0
        ):
            raise SourceMaterialChangedError(
                "persisted directory material contains invalid size, mtime, or path"
            )
        seen_paths.add(relative)
        normalized_inventory.append(
            {"path": relative, "size_bytes": size, "mtime_ns": mtime}
        )
    normalized_inventory.sort(key=lambda record: record["path"])
    observed_count = len(normalized_inventory)
    observed_total = sum(record["size_bytes"] for record in normalized_inventory)
    observed_digest = hashlib.sha256(
        _canonical_json(normalized_inventory).encode("utf-8")
    ).hexdigest()
    if (
        persisted.get("file_count") != observed_count
        or persisted.get("total_bytes") != observed_total
        or persisted.get("inventory_sha256") != observed_digest
    ):
        raise SourceMaterialChangedError(
            "persisted directory material summary does not match its inventory"
        )
    return {
        "kind": str(persisted.get("kind") or "directory-size-mtime-inventory"),
        "file_count": observed_count,
        "total_bytes": observed_total,
        "inventory": normalized_inventory,
        "inventory_sha256": observed_digest,
    }


def _check_critical_signature_stats(
    source: Path,
    raw_report: Mapping[str, Any],
) -> None:
    required = raw_report.get("required_file_signatures") or {}
    if source.is_dir():
        candidates = (
            (source / _normalized_relative_path(relative), signature)
            for relative, signature in required.items()
        )
    else:
        source_signature = raw_report.get("source_file_signature") or {}
        candidates = ((source, source_signature),) if source_signature else ()

    for path, signature in candidates:
        if not path.is_file() or not isinstance(signature, Mapping):
            raise SourceMaterialChangedError(
                f"Signed vendor source file is missing: {path.name}"
            )
        observed = path.stat()
        expected_size = signature.get("size_bytes")
        expected_mtime = signature.get("mtime_ns")
        if expected_size is not None and int(expected_size) != observed.st_size:
            raise SourceMaterialChangedError(
                f"Signed vendor source size changed: {path.name}"
            )
        if expected_mtime is not None and int(expected_mtime) != observed.st_mtime_ns:
            raise SourceMaterialChangedError(
                f"Signed vendor source mtime changed: {path.name}"
            )
        inspection_stat = signature.get("inspection_stat")
        if isinstance(inspection_stat, Mapping):
            expected_identity = (
                inspection_stat.get("device"),
                inspection_stat.get("inode"),
                inspection_stat.get("ctime_ns"),
            )
            if not all(isinstance(value, int) for value in expected_identity):
                raise SourceMaterialChangedError(
                    f"Signed vendor source inspection identity is invalid: {path.name}"
                )
            observed_identity = (
                observed.st_dev,
                observed.st_ino,
                observed.st_ctime_ns,
            )
            if observed_identity != expected_identity:
                raise SourceMaterialChangedError(
                    f"Signed vendor source identity changed: {path.name}"
                )
            continue
        expected_sha256 = str(signature.get("sha256") or "")
        observed_sha256 = _sha256_file(path)
        after_hash = path.stat()
        if (
            observed.st_dev != after_hash.st_dev
            or observed.st_ino != after_hash.st_ino
            or observed.st_size != after_hash.st_size
            or observed.st_mtime_ns != after_hash.st_mtime_ns
        ):
            raise SourceMaterialChangedError(
                f"Signed vendor source changed while hashing: {path.name}"
            )
        if not hmac.compare_digest(observed_sha256, expected_sha256):
            raise SourceMaterialChangedError(
                f"Signed vendor source content changed: {path.name}"
            )


def collect_raw_source_material(raw_report: Mapping[str, Any]) -> dict[str, Any]:
    """Collect path-free cache material for an inspected vendor source.

    Agilent directories bind every relative filename, size and nanosecond mtime,
    in addition to the cryptographic signatures already present in
    ``raw_report``.  If the source path is no longer available, callers may
    persist that result under ``raw_report['directory_material']`` and reuse it.
    """

    if not isinstance(raw_report, Mapping):
        raise TypeError("raw_report must be a mapping")
    source_path_value = str(raw_report.get("source_path") or "").strip()
    source_path = Path(source_path_value).expanduser() if source_path_value else None
    if source_path is not None and source_path.exists():
        resolved = source_path.resolve()
        _check_critical_signature_stats(resolved, raw_report)
        if resolved.is_dir():
            directory = _directory_material(resolved)
        elif resolved.is_file():
            directory = _file_material(resolved)
        else:
            raise SourceMaterialChangedError(
                "Vendor source is neither a regular file nor a directory"
            )
    else:
        persisted = raw_report.get("directory_material")
        if not isinstance(persisted, Mapping):
            raise SourceMaterialChangedError(
                "Vendor source is unavailable and raw_report has no persisted "
                "directory_material"
            )
        directory = _validated_persisted_directory_material(persisted)

    source_kind = str(raw_report.get("source_kind") or "").strip()
    raw_name = str(raw_report.get("raw_name") or "").strip()
    if not source_kind or not raw_name:
        raise ValueError("raw_report source_kind and raw_name must be non-empty")
    critical_signatures = _critical_signatures(raw_report)
    if source_kind == "agilent_d" and not critical_signatures["required_files"]:
        raise ValueError("Agilent raw_report must include required file signatures")
    if source_kind == "thermo_raw" and not critical_signatures["source_file"]:
        raise ValueError("Thermo raw_report must include a source file signature")
    acquisition = raw_report.get("acquisition") or {}
    if not isinstance(acquisition, Mapping):
        raise ValueError("raw_report.acquisition must be a mapping")
    polarities = sorted(
        {
            str(value).strip().upper()
            for value in (acquisition.get("polarities") or [])
            if str(value).strip()
        }
    )
    validation_expectations = {
        "ms1_expected": acquisition.get("ms1_expected"),
        "ms2_expected": acquisition.get("ms2_expected"),
        "polarities": polarities,
        "isolation_width_fallback_da": acquisition.get("isolation_width_da"),
    }

    return {
        "source_kind": source_kind,
        # Vendor readers can embed this name in mzML provenance, so it is
        # semantic input even though the absolute source path is deliberately
        # excluded.
        "raw_name": raw_name,
        "critical_signatures": critical_signatures,
        "directory_material": directory,
        "validation_expectations": validation_expectations,
    }


def _normalized_converter_identity(
    converter_identity: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(converter_identity, Mapping) or not converter_identity:
        raise ValueError("converter_identity must be a non-empty mapping")
    normalized_converter = _json_clone(converter_identity)
    backend = str(normalized_converter.get("backend") or "").strip().lower()
    if not backend:
        raise ValueError("converter_identity.backend must be non-empty")
    if (
        backend == "docker"
        and not str(normalized_converter.get("image_identity") or "").strip()
    ):
        raise ValueError(
            "Docker converter_identity.image_identity must be an immutable image ID"
        )
    if (
        backend == "local"
        and not str(normalized_converter.get("version") or "").strip()
    ):
        raise ValueError("Local converter_identity.version must be non-empty")
    return normalized_converter


def _normalized_conversion_policy(
    resolved_peak_picking: str,
    conversion_policy: Mapping[str, Any] | None,
) -> tuple[str, dict[str, Any]]:
    peak_picking = str(resolved_peak_picking or "").strip().lower()
    if peak_picking not in PEAK_PICKING_STRATEGIES:
        raise ValueError("resolved_peak_picking must be preserve, vendor, or openms")
    default_policy = {
        "output_format": "mzML",
        "precision_bits": 64,
        "zlib": True,
        "ms_level_filter": None,
        "peak_picking_strategy_resolved": peak_picking,
        "peak_picking_filter": (
            "peakPicking true 1-" if peak_picking == "vendor" else None
        ),
        "peak_picking_filter_order": (
            "first_filter" if peak_picking == "vendor" else None
        ),
        "downstream_openms_peak_picking": peak_picking == "openms",
        "preserves_ms1_and_ms2": True,
    }
    if conversion_policy:
        normalized_policy = _json_clone(conversion_policy)
        policy_peak_picking = normalized_policy.get("peak_picking_strategy_resolved")
        if (
            policy_peak_picking is not None
            and str(policy_peak_picking).strip().lower() != peak_picking
        ):
            raise ValueError(
                "conversion_policy peak-picking strategy disagrees with the "
                "resolved strategy"
            )
        # The requested spelling and audited command/output path are not
        # converter semantics.  Binding either would make ``auto -> vendor``
        # miss an explicit ``vendor`` entry or bind a task UUID into the key.
        normalized_policy.pop("peak_picking_strategy_requested", None)
        normalized_policy.pop("command_audit_without_source_path", None)
        # Reusing a process-local preflight changes latency, not mzML bytes.
        normalized_policy.pop("converter_preflight_reused", None)
        default_policy.update(normalized_policy)
    return peak_picking, default_policy


def build_conversion_cache_key(
    raw_report: Mapping[str, Any],
    *,
    converter_identity: Mapping[str, Any],
    resolved_peak_picking: str,
    validation_rules_version: str,
    conversion_policy: Mapping[str, Any] | None = None,
) -> ConversionCacheKey:
    """Build a cache key from source, converter, policy and validation identity."""

    normalized_converter = _normalized_converter_identity(converter_identity)
    peak_picking, normalized_policy = _normalized_conversion_policy(
        resolved_peak_picking, conversion_policy
    )
    rules_version = str(validation_rules_version or "").strip()
    if not rules_version:
        raise ValueError("validation_rules_version must be non-empty")
    material = {
        "schema_version": CACHE_KEY_SCHEMA_VERSION,
        "source": collect_raw_source_material(raw_report),
        "converter_identity": normalized_converter,
        "resolved_peak_picking": peak_picking,
        "conversion_policy": normalized_policy,
        "validation_rules_version": rules_version,
    }
    canonical = _canonical_json(material).encode("utf-8")
    return ConversionCacheKey(
        digest=hashlib.sha256(canonical).hexdigest(),
        material=material,
    )


def _validated_cache_key(cache_key: ConversionCacheKey) -> str:
    if not isinstance(cache_key, ConversionCacheKey):
        raise TypeError("cache_key must be a ConversionCacheKey")
    if not _CACHE_KEY_PATTERN.fullmatch(cache_key.digest):
        raise ValueError("cache key must be a lowercase SHA-256 digest")
    expected = hashlib.sha256(
        _canonical_json(cache_key.material).encode("utf-8")
    ).hexdigest()
    if not hmac.compare_digest(expected, cache_key.digest):
        raise ValueError("cache key does not match its canonical material")
    return cache_key.digest


def _copy_and_hash(source: Path, destination: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    before = source.stat()
    with source.open("rb") as input_handle, destination.open("xb") as output_handle:
        while chunk := input_handle.read(4 * 1024 * 1024):
            output_handle.write(chunk)
            digest.update(chunk)
            size += len(chunk)
        output_handle.flush()
        os.fsync(output_handle.fileno())
    after = source.stat()
    if (
        before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or before.st_ino != after.st_ino
    ):
        with contextlib.suppress(FileNotFoundError):
            destination.unlink()
        raise CacheEntryInvalidError("mzML changed while it was copied")
    return size, digest.hexdigest()


def _try_apfs_clone(source: Path, destination: Path) -> bool:
    """Use Darwin clonefile when available; return False for safe fallbacks."""

    if sys.platform != "darwin" or destination.exists():
        return False
    library_name = ctypes.util.find_library("c") or "libc.dylib"
    try:
        libc = ctypes.CDLL(library_name, use_errno=True)
        clonefile = libc.clonefile
    except (AttributeError, OSError):
        return False
    clonefile.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
    clonefile.restype = ctypes.c_int
    result = clonefile(os.fsencode(source), os.fsencode(destination), 0)
    if result == 0:
        return True
    with contextlib.suppress(FileNotFoundError):
        destination.unlink()
    return False


class VendorConversionCache:
    """Verified content-addressed cache for converted mzML files."""

    def __init__(self, root: str | Path):
        self.root = Path(root).expanduser().resolve(strict=False)

    def entry_directory(self, key: str) -> Path:
        if not _CACHE_KEY_PATTERN.fullmatch(str(key)):
            raise ValueError("cache key must be a lowercase SHA-256 digest")
        return self.root / key[:2] / key

    def manifest_path(self, key: str) -> Path:
        return self.entry_directory(key) / _CACHE_MANIFEST_NAME

    def payload_path(self, key: str) -> Path:
        return self.entry_directory(key) / _CACHE_PAYLOAD_NAME

    def _lock_path(self, key: str) -> Path:
        return self.root / key[:2] / f".{key}.lock"

    def _validate_entry_unlocked(self, key: str) -> CacheHit:
        manifest_path = self.manifest_path(key)
        payload_path = self.payload_path(key)
        if manifest_path.is_symlink() or payload_path.is_symlink():
            raise CacheEntryInvalidError("cache manifest or payload is a symlink")
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CacheEntryInvalidError("cache manifest is unreadable") from exc
        if not isinstance(manifest, dict):
            raise CacheEntryInvalidError("cache manifest must be a JSON object")
        if manifest.get("schema_version") != CACHE_MANIFEST_SCHEMA_VERSION:
            raise CacheEntryInvalidError("cache manifest schema is unsupported")
        if manifest.get("status") != "completed":
            raise CacheEntryInvalidError("cache entry status is not completed")
        if manifest.get("output_validation_status") != "valid":
            raise CacheEntryInvalidError("cached mzML validation status is not valid")
        if manifest.get("cache_key") != key:
            raise CacheEntryInvalidError("cache manifest key does not match its path")
        key_material = manifest.get("key_material")
        if not isinstance(key_material, dict):
            raise CacheEntryInvalidError("cache key material is missing")
        material_digest = hashlib.sha256(
            _canonical_json(key_material).encode("utf-8")
        ).hexdigest()
        if not hmac.compare_digest(material_digest, key):
            raise CacheEntryInvalidError("cache key material digest does not match")

        payload = manifest.get("payload")
        if not isinstance(payload, dict):
            raise CacheEntryInvalidError("cache payload metadata is missing")
        if payload.get("filename") != _CACHE_PAYLOAD_NAME:
            raise CacheEntryInvalidError("cache payload filename is invalid")
        if not payload_path.is_file():
            raise CacheEntryInvalidError("cached mzML is missing")
        expected_size = payload.get("size_bytes")
        expected_hash = str(payload.get("sha256") or "")
        if not isinstance(expected_size, int) or expected_size < 0:
            raise CacheEntryInvalidError("cached mzML size is invalid")
        if not _CACHE_KEY_PATTERN.fullmatch(expected_hash):
            raise CacheEntryInvalidError("cached mzML SHA-256 is invalid")

        conversion_audit = manifest.get("conversion_audit")
        if not isinstance(conversion_audit, dict):
            raise CacheEntryInvalidError("cache conversion audit is missing")
        if conversion_audit.get("status") != "completed":
            raise CacheEntryInvalidError(
                "cache conversion audit status is not completed"
            )
        if conversion_audit.get("output_status") != "valid":
            raise CacheEntryInvalidError(
                "cache conversion audit output status is not valid"
            )
        audit_output = conversion_audit.get("output")
        if not isinstance(audit_output, dict):
            raise CacheEntryInvalidError("cache conversion audit output is missing")
        if audit_output.get("status") != "valid":
            raise CacheEntryInvalidError(
                "cache conversion audit output validation status is not valid"
            )
        audit_size = audit_output.get("size_bytes")
        if not isinstance(audit_size, int) or audit_size < 0:
            raise CacheEntryInvalidError(
                "cache conversion audit output size is invalid"
            )
        if audit_size != expected_size:
            raise CacheEntryInvalidError(
                "cache conversion audit output size does not match the payload"
            )
        audit_hash = str(audit_output.get("sha256") or "")
        if not _CACHE_KEY_PATTERN.fullmatch(audit_hash):
            raise CacheEntryInvalidError(
                "cache conversion audit output SHA-256 is invalid"
            )
        if not hmac.compare_digest(audit_hash, expected_hash):
            raise CacheEntryInvalidError(
                "cache conversion audit output SHA-256 does not match the payload"
            )

        if payload_path.stat().st_size != expected_size:
            raise CacheEntryInvalidError("cached mzML size no longer matches")
        observed_hash = _sha256_file(payload_path)
        if not hmac.compare_digest(observed_hash, expected_hash):
            raise CacheEntryInvalidError("cached mzML SHA-256 no longer matches")
        return CacheHit(
            key=key,
            payload_path=payload_path,
            manifest_path=manifest_path,
            manifest=manifest,
        )

    def lookup(self, cache_key: ConversionCacheKey) -> CacheHit | None:
        """Return a hit only after manifest, status, size and hash validation."""

        key = _validated_cache_key(cache_key)
        entry = self.entry_directory(key)
        if not entry.is_dir() or entry.is_symlink():
            return None
        with _file_lock(self._lock_path(key), shared=True):
            try:
                return self._validate_entry_unlocked(key)
            except CacheEntryInvalidError:
                return None

    @staticmethod
    def _validate_conversion_manifest(
        mzml_path: Path,
        conversion_manifest: Mapping[str, Any],
    ) -> tuple[int, str]:
        if conversion_manifest.get("status") != "completed":
            raise CacheEntryInvalidError("conversion manifest status must be completed")
        output = conversion_manifest.get("output")
        if not isinstance(output, Mapping) or output.get("status") != "valid":
            raise CacheEntryInvalidError(
                "conversion output validation status must be valid"
            )
        if not mzml_path.is_file() or mzml_path.is_symlink():
            raise CacheEntryInvalidError("converted mzML must be a regular file")
        expected_size = output.get("size_bytes")
        expected_hash = str(output.get("sha256") or "")
        if not isinstance(expected_size, int) or expected_size < 0:
            raise CacheEntryInvalidError("conversion output size_bytes must be present")
        if not _CACHE_KEY_PATTERN.fullmatch(expected_hash):
            raise CacheEntryInvalidError("conversion output sha256 must be present")
        return expected_size, expected_hash

    @staticmethod
    def _validate_manifest_binding(
        cache_key: ConversionCacheKey,
        conversion_manifest: Mapping[str, Any],
        validation_rules_version: str,
    ) -> None:
        key_material = cache_key.material
        if key_material.get("schema_version") != CACHE_KEY_SCHEMA_VERSION:
            raise CacheEntryInvalidError("cache key schema is unsupported")
        source_report = conversion_manifest.get("source")
        if not isinstance(source_report, Mapping):
            raise CacheEntryInvalidError(
                "conversion manifest source report must be present"
            )
        observed_source = collect_raw_source_material(source_report)
        if _canonical_json(observed_source) != _canonical_json(
            key_material.get("source")
        ):
            raise CacheEntryInvalidError(
                "conversion manifest source does not match the cache key"
            )

        converter = conversion_manifest.get("converter")
        if not isinstance(converter, Mapping):
            raise CacheEntryInvalidError(
                "conversion manifest converter identity must be present"
            )
        observed_converter = _normalized_converter_identity(converter)
        if _canonical_json(observed_converter) != _canonical_json(
            key_material.get("converter_identity")
        ):
            raise CacheEntryInvalidError(
                "conversion manifest converter does not match the cache key"
            )

        conversion_policy = conversion_manifest.get("conversion_policy")
        if not isinstance(conversion_policy, Mapping):
            raise CacheEntryInvalidError("conversion manifest policy must be present")
        resolved = str(conversion_policy.get("peak_picking_strategy_resolved") or "")
        _, observed_policy = _normalized_conversion_policy(resolved, conversion_policy)
        if _canonical_json(observed_policy) != _canonical_json(
            key_material.get("conversion_policy")
        ):
            raise CacheEntryInvalidError(
                "conversion manifest policy does not match the cache key"
            )
        if resolved.strip().lower() != key_material.get("resolved_peak_picking"):
            raise CacheEntryInvalidError(
                "conversion peak-picking result does not match the cache key"
            )

        rules_version = str(validation_rules_version or "").strip()
        if not rules_version:
            raise CacheEntryInvalidError(
                "validation_rules_version must be provided when storing"
            )
        if rules_version != key_material.get("validation_rules_version"):
            raise CacheEntryInvalidError(
                "validation rules version does not match the cache key"
            )

    def store(
        self,
        cache_key: ConversionCacheKey,
        mzml_path: str | Path,
        conversion_manifest: Mapping[str, Any],
        *,
        validation_rules_version: str,
    ) -> CacheHit:
        """Atomically store one already validated converter output."""

        key = _validated_cache_key(cache_key)
        source = Path(mzml_path).expanduser().resolve()
        expected_size, expected_hash = self._validate_conversion_manifest(
            source, conversion_manifest
        )
        self._validate_manifest_binding(
            cache_key, conversion_manifest, validation_rules_version
        )
        entry = self.entry_directory(key)
        entry.parent.mkdir(parents=True, exist_ok=True)
        with _file_lock(self._lock_path(key)):
            if entry.exists() and (not entry.is_dir() or entry.is_symlink()):
                raise CacheEntryInvalidError("cache entry path is not a directory")
            entry.mkdir(parents=True, exist_ok=True)
            try:
                return self._validate_entry_unlocked(key)
            except CacheEntryInvalidError:
                pass

            temporary_payload = entry / (
                f".{_CACHE_PAYLOAD_NAME}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
            )
            try:
                observed_size, observed_hash = _copy_and_hash(source, temporary_payload)
                if observed_size != expected_size:
                    raise CacheEntryInvalidError(
                        "converted mzML size disagrees with its validation manifest"
                    )
                if not hmac.compare_digest(observed_hash, expected_hash):
                    raise CacheEntryInvalidError(
                        "converted mzML SHA-256 disagrees with its validation manifest"
                    )
                os.chmod(temporary_payload, 0o444)
                with temporary_payload.open("rb") as handle:
                    os.fsync(handle.fileno())
                os.replace(temporary_payload, self.payload_path(key))
                _fsync_directory(entry)
            finally:
                with contextlib.suppress(FileNotFoundError):
                    temporary_payload.unlink()

            cache_manifest = {
                "schema_version": CACHE_MANIFEST_SCHEMA_VERSION,
                "status": "completed",
                "cache_key": key,
                "key_material": _json_clone(cache_key.material),
                "created_at": _utc_now(),
                "output_validation_status": "valid",
                "payload": {
                    "filename": _CACHE_PAYLOAD_NAME,
                    "size_bytes": observed_size,
                    "sha256": observed_hash,
                },
                "conversion_audit": {
                    "status": conversion_manifest.get("status"),
                    "output_status": (conversion_manifest.get("output") or {}).get(
                        "status"
                    ),
                    "converter": _json_clone(
                        conversion_manifest.get("converter") or {}
                    ),
                    "conversion_policy": _json_clone(
                        conversion_manifest.get("conversion_policy") or {}
                    ),
                    "output": {
                        field: value
                        for field, value in _json_clone(
                            conversion_manifest.get("output") or {}
                        ).items()
                        if field != "path"
                    },
                },
            }
            _atomic_write_json(self.manifest_path(key), cache_manifest)
            return self._validate_entry_unlocked(key)

    def materialize(
        self,
        cache_key: ConversionCacheKey,
        destination: str | Path,
        *,
        overwrite: bool = False,
    ) -> dict[str, Any]:
        """Publish a verified, independently mutable task mzML from the cache."""

        key = _validated_cache_key(cache_key)
        requested_target = Path(destination).expanduser()
        if requested_target.suffix.lower() != ".mzml":
            raise ValueError("materialized destination must end in .mzML")
        requested_target.parent.mkdir(parents=True, exist_ok=True)
        # Resolve the containing directory so a deliberate managed-volume
        # symlink remains supported, but never resolve the leaf.  Resolving a
        # dangling leaf symlink would publish outside the caller's directory.
        target = requested_target.parent.resolve() / requested_target.name
        if target.is_symlink():
            raise CacheEntryInvalidError(
                "materialized destination must not be a symbolic link"
            )
        temporary = target.parent / (
            f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
        )
        source_link = target.parent / (
            f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.cachelink"
        )
        if target.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite existing output: {target}")

        with _file_lock(self._lock_path(key), shared=True):
            hit = self._validate_entry_unlocked(key)
            payload = hit.manifest["payload"]
            method = ""
            try:
                if _try_apfs_clone(hit.payload_path, temporary):
                    method = "apfs_clone"
                else:
                    try:
                        os.link(
                            hit.payload_path,
                            source_link,
                            follow_symlinks=False,
                        )
                    except OSError:
                        _copy_and_hash(hit.payload_path, temporary)
                        method = "copy"
                    else:
                        # The hard link pins the verified cache inode only while
                        # copying.  It is never published as the task output.
                        _copy_and_hash(source_link, temporary)
                        method = "hardlink_then_copy"
                    finally:
                        with contextlib.suppress(FileNotFoundError):
                            source_link.unlink()

                temporary_stat = temporary.stat()
                payload_stat = hit.payload_path.stat()
                if (temporary_stat.st_dev, temporary_stat.st_ino) == (
                    payload_stat.st_dev,
                    payload_stat.st_ino,
                ):
                    raise CacheEntryInvalidError(
                        "materialized task output still shares the cache inode"
                    )
                os.chmod(
                    temporary,
                    stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH,
                )
                with temporary.open("rb") as handle:
                    os.fsync(handle.fileno())
                observed_size = temporary.stat().st_size
                observed_hash = _sha256_file(temporary)
                if observed_size != payload["size_bytes"]:
                    raise CacheEntryInvalidError(
                        "materialized mzML size does not match the cache"
                    )
                if not hmac.compare_digest(observed_hash, payload["sha256"]):
                    raise CacheEntryInvalidError(
                        "materialized mzML SHA-256 does not match the cache"
                    )

                if overwrite:
                    os.replace(temporary, target)
                else:
                    # link(temp, target) is an atomic no-clobber publish.  Once
                    # temp is unlinked, target remains the sole name for its
                    # task-owned inode.
                    os.link(temporary, target, follow_symlinks=False)
                    temporary.unlink()
                _fsync_directory(target.parent)
            finally:
                with contextlib.suppress(FileNotFoundError):
                    temporary.unlink()
                with contextlib.suppress(FileNotFoundError):
                    source_link.unlink()

        return {
            "status": "materialized",
            "cache_key": key,
            "cache_hit": True,
            "method": method,
            "path": str(target),
            "size_bytes": payload["size_bytes"],
            "sha256": payload["sha256"],
        }


class ConversionLedger:
    """Atomic per-item conversion state with restart recovery helpers."""

    _allowed_transitions = {
        "queued": frozenset({"converting", "failed"}),
        "converting": frozenset({"validated", "failed"}),
        "validated": frozenset({"registered", "failed"}),
        "registered": frozenset(),
        "failed": frozenset({"queued"}),
    }

    def __init__(
        self,
        path: str | Path,
        *,
        job_id: str | None = None,
        clock: Callable[[], str] = _utc_now,
        worker_id: str | None = None,
        lease_seconds: float = DEFAULT_LEDGER_LEASE_SECONDS,
        time_source: Callable[[], float] = time.time,
    ):
        self.path = Path(path).expanduser().resolve(strict=False)
        self.lock_path = self.path.with_name(f".{self.path.name}.lock")
        self.job_id = str(job_id or self.path.stem).strip()
        if not self.job_id:
            raise ValueError("job_id must be non-empty")
        self._clock = clock
        self.worker_id = str(
            worker_id or f"pid:{os.getpid()}:{uuid.uuid4().hex}"
        ).strip()
        if not self.worker_id:
            raise ValueError("worker_id must be non-empty")
        self.lease_seconds = float(lease_seconds)
        if self.lease_seconds <= 0:
            raise ValueError("lease_seconds must be positive")
        self._time_source = time_source

    def _empty(self) -> dict[str, Any]:
        now = self._clock()
        return {
            "schema_version": LEDGER_SCHEMA_VERSION,
            "job_id": self.job_id,
            "created_at": now,
            "updated_at": now,
            "items": {},
        }

    def _read_unlocked(self) -> dict[str, Any]:
        if not self.path.exists():
            return self._empty()
        try:
            document = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise LedgerStateError("conversion ledger is unreadable") from exc
        if not isinstance(document, dict):
            raise LedgerStateError("conversion ledger must be a JSON object")
        if document.get("schema_version") != LEDGER_SCHEMA_VERSION:
            raise LedgerStateError("conversion ledger schema is unsupported")
        if str(document.get("job_id") or "") != self.job_id:
            raise LedgerStateError("conversion ledger job_id does not match")
        items = document.get("items")
        if not isinstance(items, dict):
            raise LedgerStateError("conversion ledger items must be an object")
        for item_id, item in items.items():
            if not isinstance(item_id, str) or not isinstance(item, dict):
                raise LedgerStateError("conversion ledger item is malformed")
            if item.get("item_id") != item_id:
                raise LedgerStateError(
                    f"conversion ledger item {item_id!r} has a mismatched item_id"
                )
            if item.get("state") not in LEDGER_STATES:
                raise LedgerStateError(
                    f"conversion ledger item {item_id!r} has an invalid state"
                )
            if not isinstance(item.get("attempt"), int) or item["attempt"] < 0:
                raise LedgerStateError(
                    f"conversion ledger item {item_id!r} has an invalid attempt"
                )
            if not isinstance(item.get("history"), list):
                raise LedgerStateError(
                    f"conversion ledger item {item_id!r} has invalid history"
                )
            if item.get("state") in {"converting", "validated"}:
                if not str(item.get("claim_token") or ""):
                    raise LedgerStateError(
                        f"conversion ledger item {item_id!r} has no claim token"
                    )
                if not str(item.get("owner_id") or ""):
                    raise LedgerStateError(
                        f"conversion ledger item {item_id!r} has no owner"
                    )
            if item.get("state") == "converting" and not isinstance(
                item.get("lease_expires_at"), (int, float)
            ):
                raise LedgerStateError(
                    f"conversion ledger item {item_id!r} has no lease expiry"
                )
            if item.get("state") == "validated" and not str(
                item.get("registration_token") or ""
            ):
                raise LedgerStateError(
                    f"conversion ledger item {item_id!r} has no registration token"
                )
        return document

    def _write_unlocked(self, document: dict[str, Any]) -> None:
        document["updated_at"] = self._clock()
        _atomic_write_json(self.path, document)

    def snapshot(self) -> dict[str, Any]:
        with _file_lock(self.lock_path, shared=True):
            return _json_clone(self._read_unlocked())

    def get(self, item_id: str) -> dict[str, Any] | None:
        item = self.snapshot()["items"].get(str(item_id))
        return _json_clone(item) if item is not None else None

    def queue(
        self,
        item_id: str,
        *,
        metadata: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Create a queued item, or explicitly retry an existing failed item."""

        identifier = str(item_id).strip()
        if not identifier:
            raise ValueError("item_id must be non-empty")
        with _file_lock(self.lock_path):
            document = self._read_unlocked()
            existing = document["items"].get(identifier)
            now = self._clock()
            if existing is None:
                item = {
                    "item_id": identifier,
                    "state": "queued",
                    "attempt": 0,
                    "metadata": _json_clone(metadata or {}),
                    "updated_at": now,
                    "history": [{"state": "queued", "at": now}],
                }
                document["items"][identifier] = item
            elif existing.get("state") == "failed":
                item = existing
                item["state"] = "queued"
                item["updated_at"] = now
                item.pop("error", None)
                item.pop("details", None)
                item.pop("claim_token", None)
                item.pop("owner_id", None)
                item.pop("lease_expires_at", None)
                item.pop("registration_token", None)
                item.pop("recovered_from", None)
                if metadata is not None:
                    item["metadata"] = _json_clone(metadata)
                item.setdefault("history", []).append(
                    {"state": "queued", "at": now, "reason": "explicit_retry"}
                )
            elif existing.get("state") == "queued":
                item = existing
                if metadata is not None:
                    item["metadata"] = _json_clone(metadata)
                    item["updated_at"] = now
            else:
                raise LedgerStateError(
                    f"cannot queue {identifier!r} from {existing.get('state')!r}"
                )
            self._write_unlocked(document)
            return _json_clone(item)

    def claim(
        self,
        item_id: str,
        *,
        metadata: Mapping[str, Any],
        details: Mapping[str, Any] | None = None,
        owner_id: str | None = None,
        lease_seconds: float | None = None,
    ) -> dict[str, Any]:
        """Atomically create/retry a queue item and claim it once.

        Combining these operations prevents two HTTP requests from both
        observing an unclaimed token between separate ``queue`` and
        ``mark_converting`` writes.
        """

        identifier = str(item_id).strip()
        if not identifier:
            raise ValueError("item_id must be non-empty")
        normalized_metadata = _json_clone(metadata)
        duration = float(self.lease_seconds if lease_seconds is None else lease_seconds)
        if duration <= 0:
            raise ValueError("lease_seconds must be positive")
        with _file_lock(self.lock_path):
            document = self._read_unlocked()
            item = document["items"].get(identifier)
            now = self._clock()
            if item is None:
                item = {
                    "item_id": identifier,
                    "state": "queued",
                    "attempt": 0,
                    "metadata": normalized_metadata,
                    "updated_at": now,
                    "history": [{"state": "queued", "at": now}],
                }
                document["items"][identifier] = item
            elif item.get("state") == "failed":
                item["state"] = "queued"
                item["updated_at"] = now
                item["metadata"] = normalized_metadata
                item.pop("error", None)
                item.pop("details", None)
                item.pop("claim_token", None)
                item.pop("owner_id", None)
                item.pop("lease_expires_at", None)
                item.pop("registration_token", None)
                item.pop("recovered_from", None)
                item.setdefault("history", []).append(
                    {"state": "queued", "at": now, "reason": "explicit_retry"}
                )
            elif item.get("state") == "queued":
                if _canonical_json(item.get("metadata") or {}) != _canonical_json(
                    normalized_metadata
                ):
                    raise LedgerStateError(
                        f"queued item {identifier!r} metadata does not match"
                    )
            else:
                raise LedgerStateError(
                    f"cannot claim {identifier!r} from {item.get('state')!r}"
                )

            item["state"] = "converting"
            item["attempt"] = int(item.get("attempt") or 0) + 1
            item["claim_token"] = uuid.uuid4().hex
            item["owner_id"] = str(owner_id or self.worker_id)
            item["lease_expires_at"] = self._time_source() + duration
            item["updated_at"] = now
            item.pop("registration_token", None)
            item.pop("recovered_from", None)
            if details is not None:
                item["details"] = _json_clone(details)
            else:
                item.pop("details", None)
            item.setdefault("history", []).append(
                {
                    "state": "converting",
                    "at": now,
                    **(
                        {"details": _json_clone(details)}
                        if details is not None
                        else {}
                    ),
                }
            )
            self._write_unlocked(document)
            return _json_clone(item)

    def _transition(
        self,
        item_id: str,
        state: str,
        *,
        details: Mapping[str, Any] | None = None,
        error: str | None = None,
        claim_token: str | None = None,
        registration_token: str | None = None,
        owner_id: str | None = None,
        lease_seconds: float | None = None,
        allow_unclaimed_queued_failure: bool = False,
    ) -> dict[str, Any]:
        if state not in LEDGER_STATES:
            raise ValueError(f"unknown conversion state: {state}")
        identifier = str(item_id).strip()
        with _file_lock(self.lock_path):
            document = self._read_unlocked()
            item = document["items"].get(identifier)
            if item is None:
                raise LedgerStateError(f"unknown conversion item: {identifier!r}")
            previous = item.get("state")
            if state not in self._allowed_transitions.get(previous, frozenset()):
                raise LedgerStateError(
                    f"invalid conversion transition {previous!r} -> {state!r}"
                )
            if previous == "queued":
                if claim_token is not None:
                    raise LedgerStateError(
                        f"stale claim token supplied for queued item {identifier!r}"
                    )
                if state == "failed" and not allow_unclaimed_queued_failure:
                    raise LedgerStateError(
                        f"queued item {identifier!r} requires mark_queued_failed"
                    )
            if previous in {"converting", "validated"}:
                supplied_claim = str(claim_token or "")
                active_claim = str(item.get("claim_token") or "")
                if not supplied_claim or not hmac.compare_digest(
                    supplied_claim, active_claim
                ):
                    raise LedgerStateError(
                        f"stale or missing claim token for {identifier!r}"
                    )
            if state == "registered":
                supplied_registration = str(registration_token or "")
                expected_registration = str(item.get("registration_token") or "")
                if not supplied_registration or not hmac.compare_digest(
                    supplied_registration, expected_registration
                ):
                    raise LedgerStateError(
                        f"stale or missing registration token for {identifier!r}"
                    )
            now = self._clock()
            item["state"] = state
            item["updated_at"] = now
            if state == "converting":
                duration = float(
                    self.lease_seconds if lease_seconds is None else lease_seconds
                )
                if duration <= 0:
                    raise ValueError("lease_seconds must be positive")
                item["attempt"] = int(item.get("attempt") or 0) + 1
                item["claim_token"] = uuid.uuid4().hex
                item["owner_id"] = str(owner_id or self.worker_id)
                item["lease_expires_at"] = self._time_source() + duration
                item.pop("registration_token", None)
                item.pop("recovered_from", None)
            elif state == "validated":
                item["registration_token"] = uuid.uuid4().hex
                item.pop("lease_expires_at", None)
            elif state in {"registered", "failed"} and previous in {
                "converting",
                "validated",
            }:
                item["last_claim_token"] = item.pop("claim_token", None)
                item["last_owner_id"] = item.pop("owner_id", None)
                item.pop("lease_expires_at", None)
            if details is not None:
                item["details"] = _json_clone(details)
            else:
                item.pop("details", None)
            if error is not None:
                item["error"] = str(error)
            elif state != "failed":
                item.pop("error", None)
            event: dict[str, Any] = {"state": state, "at": now}
            if details is not None:
                event["details"] = _json_clone(details)
            if error is not None:
                event["error"] = str(error)
            item.setdefault("history", []).append(event)
            self._write_unlocked(document)
            return _json_clone(item)

    def mark_converting(
        self,
        item_id: str,
        *,
        details: Mapping[str, Any] | None = None,
        owner_id: str | None = None,
        lease_seconds: float | None = None,
    ) -> dict[str, Any]:
        return self._transition(
            item_id,
            "converting",
            details=details,
            owner_id=owner_id,
            lease_seconds=lease_seconds,
        )

    def mark_validated(
        self,
        item_id: str,
        claim_token: str,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        return self._transition(
            item_id, "validated", details=details, claim_token=claim_token
        )

    def mark_registered(
        self,
        item_id: str,
        claim_token: str,
        registration_token: str,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        return self._transition(
            item_id,
            "registered",
            details=details,
            claim_token=claim_token,
            registration_token=registration_token,
        )

    def mark_failed(
        self,
        item_id: str,
        error: str,
        claim_token: str,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        return self._transition(
            item_id,
            "failed",
            details=details,
            error=str(error),
            claim_token=claim_token,
        )

    def mark_queued_failed(
        self,
        item_id: str,
        error: str,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Fail an item before a worker claim, such as selection validation."""

        return self._transition(
            item_id,
            "failed",
            details=details,
            error=str(error),
            allow_unclaimed_queued_failure=True,
        )

    def renew_lease(
        self,
        item_id: str,
        claim_token: str,
        *,
        lease_seconds: float | None = None,
    ) -> dict[str, Any]:
        """Extend one active claim without changing its attempt token."""

        identifier = str(item_id).strip()
        duration = float(self.lease_seconds if lease_seconds is None else lease_seconds)
        if duration <= 0:
            raise ValueError("lease_seconds must be positive")
        with _file_lock(self.lock_path):
            document = self._read_unlocked()
            item = document["items"].get(identifier)
            if item is None or item.get("state") != "converting":
                raise LedgerStateError(
                    f"cannot renew a non-converting item: {identifier!r}"
                )
            active_claim = str(item.get("claim_token") or "")
            if not claim_token or not hmac.compare_digest(
                str(claim_token), active_claim
            ):
                raise LedgerStateError(
                    f"stale or missing claim token for {identifier!r}"
                )
            item["lease_expires_at"] = self._time_source() + duration
            item["updated_at"] = self._clock()
            self._write_unlocked(document)
            return _json_clone(item)

    def recover_interrupted(self, *, force: bool = False) -> list[str]:
        """Requeue expired claims, or all claims after a verified process restart."""

        recovered: list[str] = []
        with _file_lock(self.lock_path):
            document = self._read_unlocked()
            for identifier, item in sorted(document["items"].items()):
                if item.get("state") != "converting":
                    continue
                lease_expires_at = float(item.get("lease_expires_at") or 0.0)
                if not force and lease_expires_at > self._time_source():
                    continue
                now = self._clock()
                item["state"] = "queued"
                item["updated_at"] = now
                item["recovered_from"] = "converting"
                item.pop("details", None)
                item["last_claim_token"] = item.pop("claim_token", None)
                item["last_owner_id"] = item.pop("owner_id", None)
                item.pop("lease_expires_at", None)
                item.pop("registration_token", None)
                item.setdefault("history", []).append(
                    {
                        "state": "queued",
                        "at": now,
                        "reason": (
                            "forced_process_recovery" if force else "lease_expired"
                        ),
                    }
                )
                recovered.append(identifier)
            if recovered:
                self._write_unlocked(document)
        return recovered

    def recovery_plan(self) -> dict[str, list[str]]:
        """Group persisted item IDs by state without mutating the ledger."""

        items = self.snapshot()["items"]
        return {
            state: sorted(
                identifier
                for identifier, item in items.items()
                if item.get("state") == state
            )
            for state in sorted(LEDGER_STATES)
        }


__all__ = [
    "CACHE_KEY_SCHEMA_VERSION",
    "CACHE_MANIFEST_SCHEMA_VERSION",
    "LEDGER_SCHEMA_VERSION",
    "LEDGER_STATES",
    "CacheEntryInvalidError",
    "CacheHit",
    "ConversionCacheKey",
    "ConversionLedger",
    "LedgerStateError",
    "SourceMaterialChangedError",
    "VendorConversionCache",
    "VendorConversionCacheError",
    "build_conversion_cache_key",
    "collect_raw_source_material",
]
