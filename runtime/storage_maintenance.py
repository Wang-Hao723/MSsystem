"""Conservative storage-retention engine for MSsystem runtime artifacts.

The engine keeps recoverable tasks and reproducible export scratch data in
separate policies.  It never scans raw-data, task-result, library, model, log,
or active-cache roots.  Deletion is limited to recycle entries supplied by the
server's validated inventory callback and direct children of ``export_tmp``.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import stat
import threading
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


POLICY_SCHEMA_VERSION = 1
RECYCLE_RETENTION_CHOICES = frozenset({7, 30, 90})
CACHE_RETENTION_CHOICES = frozenset({1, 7, 30})
ACTIVE_TASK_STATES = frozenset({
    "pending", "uploading", "uploaded", "converting", "queued", "running",
    "finalizing", "cancelling", "cancel_requested", "stopping",
})
TERMINAL_TASK_STATES = frozenset({"completed", "failed", "cancelled"})
PREVIEW_TTL_SECONDS = 120
OVERVIEW_CACHE_SECONDS = 60
MAX_RECEIPT_BYTES = 1024 * 1024
_IDENTIFIER = re.compile(r"[A-Za-z0-9_-]{1,80}\Z")
_RECEIPT_NAME = re.compile(
    r"storage_maintenance_(recycle|cache)_\d{8}T\d{6}Z-[0-9a-f]{8}\.json\Z"
)


class StorageMaintenanceError(RuntimeError):
    def __init__(self, detail: str, *, status_code: int = 409):
        super().__init__(detail)
        self.detail = detail
        self.status_code = status_code


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _iso(value: datetime) -> str:
    return value.astimezone(timezone.utc).isoformat()


def _parse_time(value: Any) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    token = value.strip().replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(token)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def default_policy() -> dict[str, Any]:
    return {
        "schema_version": POLICY_SCHEMA_VERSION,
        "recycle": {
            "automatic_cleanup_enabled": False,
            "retention_days": 30,
            "last_enabled_at": None,
        },
        "cache": {
            "automatic_cleanup_enabled": False,
            "retention_days": 7,
            "last_enabled_at": None,
        },
        "last_cleanup": None,
        "last_auto_cleanup_date": None,
    }


def _validate_policy(payload: Any) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise ValueError("storage policy must be an object")
    allowed_root = {
        "schema_version", "recycle", "cache", "last_cleanup",
        "last_auto_cleanup_date",
    }
    if set(payload) - allowed_root:
        raise ValueError("storage policy contains unsupported fields")
    if type(payload.get("schema_version")) is not int or payload["schema_version"] != 1:
        raise ValueError("unsupported storage policy schema version")
    result = default_policy()
    for name, choices in (
        ("recycle", RECYCLE_RETENTION_CHOICES),
        ("cache", CACHE_RETENTION_CHOICES),
    ):
        section = payload.get(name)
        if not isinstance(section, dict):
            raise ValueError(f"storage policy {name} must be an object")
        if set(section) != {
            "automatic_cleanup_enabled", "retention_days", "last_enabled_at",
        }:
            raise ValueError(f"storage policy {name} fields are invalid")
        enabled = section["automatic_cleanup_enabled"]
        days = section["retention_days"]
        enabled_at = section["last_enabled_at"]
        if type(enabled) is not bool:
            raise ValueError(f"storage policy {name} enabled must be boolean")
        if type(days) is not int or days not in choices:
            raise ValueError(f"storage policy {name} retention is invalid")
        if enabled_at is not None and _parse_time(enabled_at) is None:
            raise ValueError(f"storage policy {name} enabled time is invalid")
        result[name] = {
            "automatic_cleanup_enabled": enabled,
            "retention_days": days,
            "last_enabled_at": enabled_at,
        }
    last_cleanup = payload.get("last_cleanup")
    if last_cleanup is not None and not isinstance(last_cleanup, dict):
        raise ValueError("last_cleanup must be null or an object")
    auto_date = payload.get("last_auto_cleanup_date")
    if auto_date is not None and (
        not isinstance(auto_date, str)
        or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", auto_date)
    ):
        raise ValueError("last_auto_cleanup_date is invalid")
    result["last_cleanup"] = last_cleanup
    result["last_auto_cleanup_date"] = auto_date
    return result


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(
        f"{path.name}.{os.getpid()}.{threading.get_ident()}.tmp"
    )
    try:
        with open(partial, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        partial.replace(path)
    finally:
        partial.unlink(missing_ok=True)


def _canonical_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()


def _default_open_file_detector() -> set[Path] | None:
    """Return files opened by this service process, or None on uncertainty."""

    try:
        import psutil

        return {Path(item.path).resolve(strict=False) for item in psutil.Process().open_files()}
    except Exception:
        return None


@dataclass(frozen=True)
class _CacheItem:
    name: str
    path: Path
    size_bytes: int
    file_count: int
    mtime: datetime
    expires_at: datetime
    expired: bool
    fingerprint: str
    device: int
    inode: int
    safe: bool
    reason: str | None

    def public(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "size_bytes": self.size_bytes,
            "file_count": self.file_count,
            "mtime": _iso(self.mtime),
            "expires_at": _iso(self.expires_at),
            "expired": self.expired,
            "safe": self.safe,
            "blocked_reason": self.reason,
            "fingerprint": self.fingerprint,
        }


class StorageMaintenanceManager:
    """State-bound preview and cleanup manager.

    Callback contract:
    ``tasks_snapshot()`` returns a mapping of task IDs to task dictionaries (or
    a sequence of task dictionaries). ``recycle_inventory_unlocked()`` must be
    called only while ``recycle_lock`` is held and returns the existing server
    inventory schema. ``cache_lock`` must expose non-blocking ``acquire`` and
    ``release`` and must be the same lock held by report/export generation.
    """

    def __init__(
        self,
        *,
        output_dir: Path,
        tasks_snapshot: Callable[[], Any],
        recycle_inventory_unlocked: Callable[[], Mapping[str, Any]],
        recycle_lock: Any,
        cache_lock: Any = None,
        append_audit: Callable[[dict[str, Any]], None] | None = None,
        clock: Callable[[], datetime] = _utc_now,
        open_file_detector: Callable[[], set[Path] | None] = _default_open_file_detector,
    ):
        self.output_dir = Path(output_dir)
        self.runtime_dir = self.output_dir / "runtime"
        self.policy_path = self.runtime_dir / "storage_maintenance_policy.json"
        self.recycle_root = self.output_dir / "task_recycle_bin"
        self.export_tmp_root = self.output_dir / "export_tmp"
        self.tasks_snapshot = tasks_snapshot
        self.recycle_inventory_unlocked = recycle_inventory_unlocked
        self.recycle_lock = recycle_lock
        self.cache_lock = cache_lock
        self.append_audit = append_audit
        self.clock = clock
        self.open_file_detector = open_file_detector
        self._preview_lock = threading.Lock()
        self._previews: dict[str, dict[str, Any]] = {}
        self._policy_lock = threading.RLock()
        self._overview_lock = threading.Lock()
        self._overview_cached_at = 0.0
        self._overview_cached: dict[str, Any] | None = None

    def load_policy(self) -> dict[str, Any]:
        with self._policy_lock:
            if not self.policy_path.exists():
                return default_policy()
            if self.policy_path.is_symlink() or not self.policy_path.is_file():
                raise StorageMaintenanceError("存储策略文件类型异常，已拒绝清理")
            try:
                payload = json.loads(self.policy_path.read_text(encoding="utf-8"))
                return _validate_policy(payload)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                raise StorageMaintenanceError(f"存储策略无效：{exc}") from exc

    def save_policy(
        self,
        *,
        recycle_enabled: bool,
        recycle_retention_days: int,
        cache_enabled: bool,
        cache_retention_days: int,
        confirm_auto_delete: bool,
    ) -> dict[str, Any]:
        if type(recycle_enabled) is not bool or type(cache_enabled) is not bool:
            raise StorageMaintenanceError("自动清理开关必须是布尔值", status_code=400)
        if type(recycle_retention_days) is not int or recycle_retention_days not in RECYCLE_RETENTION_CHOICES:
            raise StorageMaintenanceError("回收站保留天数只能为7、30或90天", status_code=400)
        if type(cache_retention_days) is not int or cache_retention_days not in CACHE_RETENTION_CHOICES:
            raise StorageMaintenanceError("导出临时缓存保留天数只能为1、7或30天", status_code=400)
        with self._policy_lock:
            current = self.load_policy()
            newly_enabled = (
                (recycle_enabled and not current["recycle"]["automatic_cleanup_enabled"])
                or (cache_enabled and not current["cache"]["automatic_cleanup_enabled"])
            )
            shortened_while_enabled = (
                recycle_enabled
                and current["recycle"]["automatic_cleanup_enabled"]
                and recycle_retention_days < current["recycle"]["retention_days"]
            ) or (
                cache_enabled
                and current["cache"]["automatic_cleanup_enabled"]
                and cache_retention_days < current["cache"]["retention_days"]
            )
            if (newly_enabled or shortened_while_enabled) and confirm_auto_delete is not True:
                raise StorageMaintenanceError(
                    "启用自动删除或缩短保留期必须显式确认", status_code=400
                )
            now = _iso(self.clock())
            updated = default_policy()
            updated["last_cleanup"] = current.get("last_cleanup")
            updated["last_auto_cleanup_date"] = current.get("last_auto_cleanup_date")
            for name, enabled, days in (
                ("recycle", recycle_enabled, recycle_retention_days),
                ("cache", cache_enabled, cache_retention_days),
            ):
                was_enabled = current[name]["automatic_cleanup_enabled"]
                updated[name] = {
                    "automatic_cleanup_enabled": enabled,
                    "retention_days": days,
                    "last_enabled_at": (
                        now if enabled and not was_enabled
                        else current[name].get("last_enabled_at")
                    ),
                }
            _atomic_json(self.policy_path, _validate_policy(updated))
        self.invalidate_overview()
        return updated

    def invalidate_overview(self) -> None:
        with self._overview_lock:
            self._overview_cached = None
            self._overview_cached_at = 0.0

    @staticmethod
    def _history_item(raw: Any, requested: Sequence[Any], *, failed: bool) -> dict[str, Any] | None:
        if not isinstance(raw, Mapping):
            return None
        task_id = str(raw.get("task_id") or "").strip()[:80]
        name = str(raw.get("name") or "").strip()[:240]
        if task_id and not name:
            deletion_id = str(raw.get("deletion_id") or "")
            for item in requested:
                if not isinstance(item, Mapping):
                    continue
                if (
                    str(item.get("task_id") or "") == task_id
                    and str(item.get("deletion_id") or "") == deletion_id
                ):
                    name = str(item.get("filename") or "").strip()[:240]
                    break
        if not name:
            name = task_id
        if not name:
            return None
        result: dict[str, Any] = {"name": name}
        if task_id:
            result["task_id"] = task_id
        if not failed:
            size = raw.get("size_bytes")
            if type(size) is int and size >= 0:
                result["size_bytes"] = size
        else:
            result["error"] = str(raw.get("error") or "unknown_error")[:500]
        return result

    def cleanup_history(self, limit: int = 10) -> dict[str, Any]:
        if type(limit) is not int or not 1 <= limit <= 20:
            raise StorageMaintenanceError("历史记录数量必须为1到20", status_code=400)
        runtime = self.runtime_dir
        if not runtime.exists():
            return {"records": [], "rejected_count": 0}
        if runtime.is_symlink() or not runtime.is_dir():
            raise StorageMaintenanceError("运行时目录类型异常，已拒绝读取清理历史")
        records: list[tuple[datetime, dict[str, Any]]] = []
        rejected = 0
        try:
            candidates = list(runtime.iterdir())
        except OSError as exc:
            raise StorageMaintenanceError(f"无法读取清理历史：{exc}") from exc
        for path in candidates:
            match = _RECEIPT_NAME.fullmatch(path.name)
            if not match:
                continue
            fd = None
            try:
                before = path.lstat()
                if path.is_symlink() or not stat.S_ISREG(before.st_mode):
                    raise OSError("receipt is not a plain file")
                fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
                opened = os.fstat(fd)
                if (
                    not stat.S_ISREG(opened.st_mode)
                    or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)
                    or opened.st_size > MAX_RECEIPT_BYTES
                ):
                    raise OSError("receipt identity or size rejected")
                chunks = []
                remaining = MAX_RECEIPT_BYTES + 1
                while remaining:
                    chunk = os.read(fd, min(65536, remaining))
                    if not chunk:
                        break
                    chunks.append(chunk)
                    remaining -= len(chunk)
                data = b"".join(chunks)
                if len(data) > MAX_RECEIPT_BYTES:
                    raise OSError("receipt too large")
                payload = json.loads(data.decode("utf-8"))
                if not isinstance(payload, Mapping) or payload.get("schema_version") != 1:
                    raise ValueError("unsupported receipt")
                target = match.group(1)
                if payload.get("target") != target:
                    raise ValueError("receipt target mismatch")
                finished = _parse_time(payload.get("finished_at"))
                if finished is None:
                    raise ValueError("unfinished receipt")
                requested = payload.get("requested_items") or []
                removed = [
                    item for raw in (payload.get("removed") or [])
                    if (item := self._history_item(raw, requested, failed=False)) is not None
                ]
                failures = [
                    item for raw in (payload.get("failures") or [])
                    if (item := self._history_item(raw, requested, failed=True)) is not None
                ]
                record = {
                    "target": target,
                    "finished_at": _iso(finished),
                    "removed_count": len(removed),
                    "removed_bytes": sum(int(item.get("size_bytes") or 0) for item in removed),
                    "failure_count": len(failures),
                    "removed_items": removed,
                    "failed_items": failures,
                }
                records.append((finished, record))
            except (OSError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
                rejected += 1
            finally:
                if fd is not None:
                    os.close(fd)
        records.sort(key=lambda item: item[0], reverse=True)
        return {"records": [item[1] for item in records[:limit]], "rejected_count": rejected}

    def _task_rows(self) -> list[Mapping[str, Any]]:
        snapshot = self.tasks_snapshot()
        if isinstance(snapshot, Mapping):
            values = snapshot.values()
        elif isinstance(snapshot, Sequence) and not isinstance(snapshot, (str, bytes)):
            values = snapshot
        else:
            raise StorageMaintenanceError("任务状态回调返回值无效")
        return [row for row in values if isinstance(row, Mapping)]

    def active_tasks(self) -> dict[str, Any]:
        counts: dict[str, int] = {}
        for task in self._task_rows():
            status = str(task.get("status") or "").strip().casefold()
            if status in TERMINAL_TASK_STATES:
                continue
            # Fail closed for launcher/server states added after this module:
            # only an explicitly terminal task permits export-cache deletion.
            bucket = status if status in ACTIVE_TASK_STATES else f"unknown:{status or 'missing'}"
            counts[bucket] = counts.get(bucket, 0) + 1
        return {"count": sum(counts.values()), "by_status": dict(sorted(counts.items()))}

    def _retention_anchor(
        self, base_time: datetime, policy_section: Mapping[str, Any]
    ) -> datetime:
        enabled_at = _parse_time(policy_section.get("last_enabled_at"))
        return max(base_time, enabled_at) if enabled_at else base_time

    def _recycle_state(
        self, policy: Mapping[str, Any], *, now: datetime
    ) -> dict[str, Any]:
        with self.recycle_lock:
            inventory = dict(self.recycle_inventory_unlocked())
            entries = []
            snapshot = self.tasks_snapshot()
            if isinstance(snapshot, Mapping):
                current_ids = {
                    str(task.get("id") or key)
                    for key, task in snapshot.items()
                    if isinstance(task, Mapping)
                }
            elif isinstance(snapshot, Sequence) and not isinstance(
                snapshot, (str, bytes)
            ):
                current_ids = {
                    str(task.get("id") or "")
                    for task in snapshot
                    if isinstance(task, Mapping) and task.get("id")
                }
            else:
                raise StorageMaintenanceError("任务状态回调返回值无效")
            section = policy["recycle"]
            for raw in inventory.get("active_entries") or []:
                if not isinstance(raw, Mapping):
                    continue
                deletion_id = str(raw.get("deletion_id") or "")
                task_id = str(raw.get("task_id") or "")
                reasons = list(raw.get("purge_block_reasons") or [])
                safe = bool(raw.get("safe_to_purge"))
                if not _IDENTIFIER.fullmatch(deletion_id) or not _IDENTIFIER.fullmatch(task_id):
                    safe = False
                    reasons.append("invalid_recycle_identifier")
                task_root = self.recycle_root / deletion_id / task_id
                try:
                    recycle_resolved = self.recycle_root.resolve(strict=False)
                    generation_root = self.recycle_root / deletion_id
                    stat = task_root.lstat()
                    task_root.resolve(strict=False).relative_to(recycle_resolved)
                    if (
                        self.recycle_root.is_symlink()
                        or generation_root.is_symlink()
                        or task_root.is_symlink()
                        or not task_root.is_dir()
                    ):
                        raise OSError("recycle root is not a plain directory")
                    base = _parse_time(raw.get("deleted_at")) or datetime.fromtimestamp(
                        stat.st_mtime, timezone.utc
                    )
                    device, inode = int(stat.st_dev), int(stat.st_ino)
                except (OSError, RuntimeError, ValueError):
                    safe = False
                    reasons.append("recycle_path_unavailable")
                    base = _parse_time(raw.get("deleted_at")) or now
                    device, inode = -1, -1
                if task_id in current_ids:
                    safe = False
                    reasons.append("task_id_exists_in_current_registry")
                expires = self._retention_anchor(base, section) + timedelta(
                    days=int(section["retention_days"])
                )
                entry = {
                    "deletion_id": deletion_id,
                    "task_id": task_id,
                    "filename": str(raw.get("filename") or task_id),
                    "deleted_at": raw.get("deleted_at"),
                    "expires_at": _iso(expires),
                    "expired": now >= expires,
                    "size_bytes": int(raw.get("size_bytes") or 0),
                    "file_count": int(raw.get("file_count") or 0),
                    "tree_signature": raw.get("tree_signature"),
                    "manifest_sha256": raw.get("manifest_sha256"),
                    "safe": safe,
                    "blocked_reasons": sorted(set(reasons)),
                    "device": device,
                    "inode": inode,
                }
                entry["fingerprint"] = _canonical_hash(entry)
                entries.append(entry)
            eligible = [item for item in entries if item["expired"] and item["safe"]]
            blocked = [item for item in entries if item["expired"] and not item["safe"]]
            return {
                "inventory_signature": inventory.get("signature"),
                "structural_errors": list(inventory.get("structural_errors") or []),
                "unsafe_active_count": len(inventory.get("unsafe_active_entries") or []),
                "entries": entries,
                "eligible": eligible,
                "blocked": blocked,
                "total_count": len(entries),
                "total_bytes": sum(item["size_bytes"] for item in entries),
                "expired_count": len(eligible),
                "expired_bytes": sum(item["size_bytes"] for item in eligible),
            }

    def _cache_tree(self, entry: Path, section: Mapping[str, Any], now: datetime) -> _CacheItem:
        rows: list[tuple[Any, ...]] = []
        size = 0
        files = 0
        latest_mtime = datetime.fromtimestamp(0, timezone.utc)
        safe = True
        reason: str | None = None
        try:
            root = self.export_tmp_root.resolve(strict=False)
            entry_stat = entry.lstat()
            if entry.is_symlink():
                modified = datetime.fromtimestamp(entry_stat.st_mtime, timezone.utc)
                expires = self._retention_anchor(modified, section) + timedelta(
                    days=int(section["retention_days"])
                )
                return _CacheItem(
                    name=entry.name,
                    path=entry,
                    size_bytes=0,
                    file_count=0,
                    mtime=modified,
                    expires_at=expires,
                    expired=now >= expires,
                    fingerprint=_canonical_hash([
                        (
                            ".",
                            "link",
                            entry_stat.st_dev,
                            entry_stat.st_ino,
                            entry_stat.st_size,
                            entry_stat.st_mtime_ns,
                        )
                    ]),
                    device=int(entry_stat.st_dev),
                    inode=int(entry_stat.st_ino),
                    safe=False,
                    reason="symbolic_link_present",
                )
            entry.resolve(strict=False).relative_to(root)
            stack = [entry]
            while stack:
                current = stack.pop()
                stat = current.lstat()
                relative = str(current.relative_to(entry)) or "."
                modified = datetime.fromtimestamp(stat.st_mtime, timezone.utc)
                latest_mtime = max(latest_mtime, modified)
                if current.is_symlink():
                    safe, reason = False, "symbolic_link_present"
                    rows.append((relative, "link", stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns))
                    continue
                if current.is_file():
                    size += int(stat.st_size)
                    files += 1
                    kind = "file"
                elif current.is_dir():
                    kind = "dir"
                    with os.scandir(current) as iterator:
                        stack.extend(Path(item.path) for item in iterator)
                else:
                    safe, reason = False, "unsupported_filesystem_entry"
                    kind = "other"
                rows.append((relative, kind, stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns))
            root_stat = entry.lstat()
        except (OSError, RuntimeError, ValueError):
            safe, reason = False, "cache_path_unavailable_or_outside_allowlist"
            root_stat = type("MissingStat", (), {"st_dev": -1, "st_ino": -1})()
            if latest_mtime.year == 1970:
                latest_mtime = now
        expires = self._retention_anchor(latest_mtime, section) + timedelta(
            days=int(section["retention_days"])
        )
        return _CacheItem(
            name=entry.name,
            path=entry,
            size_bytes=size,
            file_count=files,
            mtime=latest_mtime,
            expires_at=expires,
            expired=now >= expires,
            fingerprint=_canonical_hash(rows),
            device=int(root_stat.st_dev),
            inode=int(root_stat.st_ino),
            safe=safe,
            reason=reason,
        )

    def _cache_state(self, policy: Mapping[str, Any], *, now: datetime) -> dict[str, Any]:
        entries: list[_CacheItem] = []
        structural_errors: list[str] = []
        try:
            if self.export_tmp_root.exists() and self.export_tmp_root.is_symlink():
                structural_errors.append("export_tmp_root_is_symbolic_link")
            elif self.export_tmp_root.is_dir():
                entries = [
                    self._cache_tree(path, policy["cache"], now)
                    for path in sorted(self.export_tmp_root.iterdir(), key=lambda item: item.name)
                ]
        except OSError:
            structural_errors.append("export_tmp_inventory_failed")
        eligible = [item for item in entries if item.expired and item.safe]
        blocked = [item for item in entries if item.expired and not item.safe]
        return {
            "entries": entries,
            "eligible": eligible,
            "blocked": blocked,
            "structural_errors": structural_errors,
            "total_count": len(entries),
            "total_bytes": sum(item.size_bytes for item in entries),
            "expired_count": len(eligible),
            "expired_bytes": sum(item.size_bytes for item in eligible),
        }

    def overview(self, *, force: bool = False) -> dict[str, Any]:
        monotonic = time.monotonic()
        with self._overview_lock:
            if (
                not force and self._overview_cached is not None
                and monotonic - self._overview_cached_at < OVERVIEW_CACHE_SECONDS
            ):
                payload = json.loads(json.dumps(self._overview_cached))
                payload["active_tasks"] = self.active_tasks()
                return payload
        policy = self.load_policy()
        now = self.clock()
        recycle = self._recycle_state(policy, now=now)
        cache = self._cache_state(policy, now=now)
        active = self.active_tasks()
        payload = {
            "schema_version": 1,
            "policy": policy,
            "recycle": {
                key: recycle[key] for key in (
                    "total_count", "total_bytes", "expired_count", "expired_bytes",
                    "unsafe_active_count", "structural_errors",
                )
            },
            "cache": {
                "allowlist_root": "export_tmp",
                "total_count": cache["total_count"],
                "total_bytes": cache["total_bytes"],
                "expired_count": cache["expired_count"],
                "expired_bytes": cache["expired_bytes"],
                "blocked_count": len(cache["blocked"]),
                "structural_errors": cache["structural_errors"],
            },
            "active_tasks": active,
            "last_cleanup": policy.get("last_cleanup"),
        }
        with self._overview_lock:
            self._overview_cached = payload
            self._overview_cached_at = monotonic
        return json.loads(json.dumps(payload))

    @contextmanager
    def _cache_guard(self):
        lock = self.cache_lock
        if lock is None or not hasattr(lock, "acquire") or not hasattr(lock, "release"):
            yield False
            return
        try:
            acquired = bool(lock.acquire(blocking=False))
        except TypeError:
            acquired = bool(lock.acquire(False))
        if not acquired:
            yield False
            return
        try:
            yield True
        finally:
            lock.release()

    def _cache_use_block(self, candidates: Sequence[_CacheItem]) -> list[str]:
        opened = self.open_file_detector()
        if opened is None:
            return ["open_file_detection_unavailable"]
        candidate_roots = [item.path.resolve(strict=False) for item in candidates]
        for opened_path in opened:
            resolved = Path(opened_path).resolve(strict=False)
            for root in candidate_roots:
                try:
                    resolved.relative_to(root)
                    return ["export_tmp_contains_open_file"]
                except ValueError:
                    continue
        return []

    def preview(self, target: str) -> dict[str, Any]:
        if target not in {"recycle", "cache"}:
            raise StorageMaintenanceError("未知清理目标", status_code=400)
        policy = self.load_policy()
        now = self.clock()
        blocked: list[str] = []
        if target == "recycle":
            state = self._recycle_state(policy, now=now)
            if state["structural_errors"]:
                blocked.extend(f"recycle:{item}" for item in state["structural_errors"])
            blocked.extend(
                f"recycle:{item['task_id']}:{','.join(item['blocked_reasons'])}"
                for item in state["blocked"]
            )
            candidates = list(state["eligible"])
            summary = [
                {
                    "deletion_id": item["deletion_id"],
                    "task_id": item["task_id"],
                    "filename": item["filename"],
                    "size_bytes": item["size_bytes"],
                    "expires_at": item["expires_at"],
                }
                for item in candidates
            ]
            fingerprints = [item["fingerprint"] for item in candidates]
            state_signature = state["inventory_signature"]
        else:
            active = self.active_tasks()
            if active["count"]:
                blocked.append("active_tasks_present")
            with self._cache_guard() as locked:
                if not locked:
                    blocked.append("cache_lock_unavailable_or_busy")
                    state = self._cache_state(policy, now=now)
                else:
                    state = self._cache_state(policy, now=now)
                    blocked.extend(self._cache_use_block(state["eligible"]))
            blocked.extend(f"cache:{item.name}:{item.reason}" for item in state["blocked"])
            blocked.extend(f"cache:{item}" for item in state["structural_errors"])
            candidates = list(state["eligible"])
            summary = [item.public() for item in candidates]
            fingerprints = [
                _canonical_hash({
                    "name": item.name, "fingerprint": item.fingerprint,
                    "device": item.device, "inode": item.inode,
                })
                for item in candidates
            ]
            state_signature = _canonical_hash(fingerprints)
        content_fingerprint = _canonical_hash({
            "target": target,
            "policy": policy[target],
            "state_signature": state_signature,
            "items": fingerprints,
        })
        token: str | None = None
        if not blocked:
            token = uuid.uuid4().hex
            with self._preview_lock:
                current = time.time()
                self._previews = {
                    key: value for key, value in self._previews.items()
                    if value["expires_at_epoch"] >= current
                }
                self._previews[token] = {
                    "target": target,
                    "content_fingerprint": content_fingerprint,
                    "state_signature": state_signature,
                    "items": summary,
                    "item_fingerprints": fingerprints,
                    "expires_at_epoch": current + PREVIEW_TTL_SECONDS,
                }
        return {
            "target": target,
            "token": token,
            "expires_in_seconds": PREVIEW_TTL_SECONDS if token else 0,
            "item_count": len(summary),
            "size_bytes": sum(int(item.get("size_bytes") or 0) for item in summary),
            "items": summary,
            "blocked_reasons": blocked,
            "policy": policy[target],
            "content_fingerprint": content_fingerprint,
        }

    def _consume_preview(self, token: str) -> dict[str, Any]:
        with self._preview_lock:
            preview = self._previews.pop(token, None)
        if preview is None:
            raise StorageMaintenanceError("清理预览已失效，请重新预览")
        if preview["expires_at_epoch"] < time.time():
            raise StorageMaintenanceError("清理预览已过期，请重新预览")
        return preview

    def clean(self, token: str, *, confirmed: bool, automatic: bool = False) -> dict[str, Any]:
        if confirmed is not True:
            raise StorageMaintenanceError("必须显式确认清理", status_code=400)
        preview = self._consume_preview(token)
        target = preview["target"]
        if target == "recycle":
            result = self._clean_recycle(preview, automatic=automatic)
        else:
            result = self._clean_cache(preview, automatic=automatic)
        self.invalidate_overview()
        return result

    def _receipt(self, target: str, preview: Mapping[str, Any], automatic: bool) -> tuple[Path, dict[str, Any]]:
        cleanup_id = f"{self.clock().strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
        path = self.runtime_dir / f"storage_maintenance_{target}_{cleanup_id}.json"
        receipt = {
            "schema_version": 1,
            "cleanup_id": cleanup_id,
            "target": target,
            "automatic": automatic,
            "started_at": _iso(self.clock()),
            "content_fingerprint": preview["content_fingerprint"],
            "requested_items": preview["items"],
            "removed": [],
            "failures": [],
        }
        _atomic_json(path, receipt)
        return path, receipt

    def _finish_receipt(self, path: Path, receipt: dict[str, Any]) -> dict[str, Any]:
        receipt.update({
            "finished_at": _iso(self.clock()),
            "removed_count": len(receipt["removed"]),
            "removed_bytes": sum(int(item.get("size_bytes") or 0) for item in receipt["removed"]),
            "failure_count": len(receipt["failures"]),
        })
        _atomic_json(path, receipt)
        with self._policy_lock:
            policy = self.load_policy()
            policy["last_cleanup"] = {
                "target": receipt["target"],
                "at": receipt["finished_at"],
                "removed_count": receipt["removed_count"],
                "removed_bytes": receipt["removed_bytes"],
                "failure_count": receipt["failure_count"],
                "automatic": receipt["automatic"],
                "receipt": path.name,
            }
            _atomic_json(self.policy_path, _validate_policy(policy))
        audit = {
            "event": "storage_maintenance_cleanup",
            **policy["last_cleanup"],
            "cleanup_id": receipt["cleanup_id"],
        }
        if self.append_audit:
            try:
                self.append_audit(audit)
            except Exception as exc:
                receipt["audit_append_error"] = str(exc)
                _atomic_json(path, receipt)
        return {
            "cleanup_id": receipt["cleanup_id"],
            "target": receipt["target"],
            "removed_count": receipt["removed_count"],
            "removed_bytes": receipt["removed_bytes"],
            "failure_count": receipt["failure_count"],
            "failures": receipt["failures"],
            "receipt_path": str(path),
        }

    def _clean_recycle(self, preview: Mapping[str, Any], *, automatic: bool) -> dict[str, Any]:
        with self.recycle_lock:
            policy = self.load_policy()
            state = self._recycle_state(policy, now=self.clock())
            if state["structural_errors"]:
                raise StorageMaintenanceError("回收站结构检查失败；没有删除任何项目")
            current = {
                (item["deletion_id"], item["task_id"]): item
                for item in state["eligible"]
            }
            targets = []
            fingerprints = []
            for requested in preview["items"]:
                key = (requested["deletion_id"], requested["task_id"])
                item = current.get(key)
                if not item:
                    raise StorageMaintenanceError("回收站内容或保留期已变化，请重新预览")
                task_root = self.recycle_root / key[0] / key[1]
                stat = task_root.lstat()
                if task_root.is_symlink() or not task_root.is_dir() or (
                    int(stat.st_dev), int(stat.st_ino)
                ) != (item["device"], item["inode"]):
                    raise StorageMaintenanceError("回收站目录身份已变化，请重新预览")
                targets.append((item, task_root))
                fingerprints.append(item["fingerprint"])
            current_fingerprint = _canonical_hash({
                "target": "recycle", "policy": policy["recycle"],
                "state_signature": state["inventory_signature"],
                "items": fingerprints,
            })
            if current_fingerprint != preview["content_fingerprint"]:
                raise StorageMaintenanceError("回收站内容已变化，请重新预览")
            path, receipt = self._receipt("recycle", preview, automatic)
            for item, task_root in targets:
                try:
                    shutil.rmtree(task_root)
                    receipt["removed"].append({
                        "deletion_id": item["deletion_id"],
                        "task_id": item["task_id"],
                        "size_bytes": item["size_bytes"],
                    })
                    try:
                        task_root.parent.rmdir()
                    except OSError:
                        pass
                except Exception as exc:
                    receipt["failures"].append({
                        "deletion_id": item["deletion_id"],
                        "task_id": item["task_id"], "error": str(exc),
                    })
                _atomic_json(path, receipt)
            return self._finish_receipt(path, receipt)

    def _clean_cache(self, preview: Mapping[str, Any], *, automatic: bool) -> dict[str, Any]:
        if self.active_tasks()["count"]:
            raise StorageMaintenanceError("存在上传、排队或运行任务；已拒绝清理导出临时缓存")
        with self._cache_guard() as locked:
            if not locked:
                raise StorageMaintenanceError("导出锁正在使用或不可用；已拒绝清理缓存")
            policy = self.load_policy()
            state = self._cache_state(policy, now=self.clock())
            if state["structural_errors"]:
                raise StorageMaintenanceError("export_tmp结构检查失败；没有删除任何项目")
            current = {item.name: item for item in state["eligible"]}
            targets: list[_CacheItem] = []
            fingerprints = []
            for requested in preview["items"]:
                item = current.get(str(requested.get("name") or ""))
                if item is None:
                    raise StorageMaintenanceError("缓存内容或保留期已变化，请重新预览")
                stat = item.path.lstat()
                if item.path.is_symlink() or (
                    int(stat.st_dev), int(stat.st_ino)
                ) != (item.device, item.inode):
                    raise StorageMaintenanceError("缓存目录身份已变化，请重新预览")
                targets.append(item)
                fingerprints.append(_canonical_hash({
                    "name": item.name, "fingerprint": item.fingerprint,
                    "device": item.device, "inode": item.inode,
                }))
            use_blocks = self._cache_use_block(targets)
            if use_blocks:
                raise StorageMaintenanceError("导出临时缓存占用状态无法安全确认；没有删除任何项目")
            current_fingerprint = _canonical_hash({
                "target": "cache", "policy": policy["cache"],
                "state_signature": _canonical_hash(fingerprints),
                "items": fingerprints,
            })
            if current_fingerprint != preview["content_fingerprint"]:
                raise StorageMaintenanceError("缓存内容已变化，请重新预览")
            path, receipt = self._receipt("cache", preview, automatic)
            for item in targets:
                try:
                    if item.path.is_dir():
                        shutil.rmtree(item.path)
                    else:
                        item.path.unlink()
                    receipt["removed"].append({
                        "name": item.name, "size_bytes": item.size_bytes,
                    })
                except Exception as exc:
                    receipt["failures"].append({"name": item.name, "error": str(exc)})
                _atomic_json(path, receipt)
            return self._finish_receipt(path, receipt)

    def run_auto_cleanup(self) -> dict[str, Any]:
        policy = self.load_policy()
        enabled = [
            name for name in ("recycle", "cache")
            if policy[name]["automatic_cleanup_enabled"]
        ]
        if not enabled:
            return {"status": "disabled", "results": []}
        active = self.active_tasks()
        if active["count"]:
            return {"status": "skipped_active_tasks", "active_tasks": active, "results": []}
        today = self.clock().date().isoformat()
        if policy.get("last_auto_cleanup_date") == today:
            return {"status": "already_checked_today", "results": []}
        results = []
        for target in enabled:
            preview = self.preview(target)
            if preview["blocked_reasons"]:
                results.append({
                    "target": target, "status": "blocked",
                    "blocked_reasons": preview["blocked_reasons"],
                })
            elif preview["token"]:
                results.append(self.clean(preview["token"], confirmed=True, automatic=True))
        with self._policy_lock:
            current = self.load_policy()
            current["last_auto_cleanup_date"] = today
            _atomic_json(self.policy_path, _validate_policy(current))
        return {"status": "checked", "results": results}


__all__ = [
    "ACTIVE_TASK_STATES",
    "CACHE_RETENTION_CHOICES",
    "POLICY_SCHEMA_VERSION",
    "RECYCLE_RETENTION_CHOICES",
    "StorageMaintenanceError",
    "StorageMaintenanceManager",
    "TERMINAL_TASK_STATES",
    "default_policy",
]
