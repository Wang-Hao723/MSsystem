"""Read-only presentation of the frozen XBT completed batch.

This module deliberately has no dependency on ``web.server`` or analysis code.
It only reads the frozen delivery receipts and lazily materialises their nine
verified exports in an external-volume cache.
"""

from __future__ import annotations

import asyncio
import csv
import hashlib
import io
import json
import lzma
import os
import queue
import shutil
import subprocess
import threading
import zipfile
from pathlib import Path
from typing import Any, Iterator
from urllib.parse import quote

from fastapi import HTTPException, Query
from fastapi.responses import FileResponse, StreamingResponse


COMPLETED_BATCH_ID = "xbt-20261021-ms1-msms"
_TITLE = "宣痹汤 / XBT · 149份 MS1＋MS/MS＋配对空白（2026-10-21）"
_REMOTE_EXPORT_DIR = Path("/data/Mssystem/xbt_real_batch_20261021_ms2/exports")
_CACHE_FLOOR_BYTES = 50 * 1024**3
_OPERATION_CAP_BYTES = 192 * 1024**2
_COPY_TIMEOUT_SECONDS = 180
_ALLOWED_CACHE_PREFIXES = (Path("/Volumes"),)

_EXPORT_FILES = {
    "file_manifest": "all_149_files_manifest.csv.xz",
    "xbt_csv": "all_xbt_named_candidates.csv.xz",
    "xbt_jsonl": "all_xbt_named_candidates.jsonl.xz",
    "qc": "all_candidate_qc.jsonl.xz",
    "ms1": "all_134620_ms1_features.jsonl.xz",
    "msms_table": "all_msms_precursor_events.csv.xz",
    "ms1_table": "all_134620_ms1_features.csv.xz",
    "qc_table": "all_candidate_qc.csv.xz",
    "summary": "FULL_EXPORT_SUMMARY.json",
}
_TABLE_KEYS = {"file_manifest", "xbt_csv", "msms_table", "ms1_table", "qc_table"}
_EXPORT_LABELS = {
    "file_manifest": "149份文件清单",
    "xbt_csv": "XBT联合附加筛选（CSV，含未解析事件）",
    "xbt_jsonl": "XBT联合附加筛选完整记录（JSONL）",
    "qc": "候选质量控制完整记录（JSONL）",
    "ms1": "MS1特征完整记录（JSONL）",
    "msms_table": "MS/MS前体事件表",
    "ms1_table": "MS1特征表",
    "qc_table": "候选质量控制表",
    "summary": "完整导出摘要",
}
_LOCAL_FILES = {
    "report": "FINAL_RESULTS_OVERVIEW.html",
    "index": "FINAL_DELIVERY_INDEX.json",
    "packaging": "FINAL_PACKAGING_RECEIPT.json",
    "verification": "FINAL_DELIVERY_VERIFICATION.json",
    "native_state": "STATE.json",
    "conclusions_state": "CONCLUSIONS_STATE.json",
}
_locks_guard = threading.Lock()
_cache_locks: dict[str, threading.Lock] = {}


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"冻结交付文件无法读取：{path.name}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"冻结交付文件格式错误：{path.name}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _disposition(filename: str, *, inline: bool = False) -> str:
    mode = "inline" if inline else "attachment"
    fallback = "".join(c if ord(c) < 128 and c not in '\"\\' else "_" for c in filename)
    return f"{mode}; filename=\"{fallback}\"; filename*=UTF-8''{quote(filename)}"


def _cache_lock(path: Path) -> threading.Lock:
    key = str(path)
    with _locks_guard:
        return _cache_locks.setdefault(key, threading.Lock())


def _validate_cache_volume(cache_root: Path) -> None:
    resolved = cache_root.resolve()
    matched = next(
        (prefix for prefix in _ALLOWED_CACHE_PREFIXES if resolved == prefix or prefix in resolved.parents),
        None,
    )
    if matched is None:
        raise RuntimeError("完成批次缓存必须位于已挂载的 /Volumes 外置卷，禁止写入系统盘")
    if matched == Path("/Volumes"):
        relative = resolved.relative_to(matched)
        if not relative.parts:
            raise RuntimeError("完成批次缓存必须位于具体的已挂载外置卷")
        volume_root = matched / relative.parts[0]
        if not volume_root.is_dir() or not os.path.ismount(volume_root):
            raise RuntimeError("完成批次外置缓存卷未挂载，禁止回退到系统盘")


class CompletedBatchService:
    def __init__(self, batch_root: Path, cache_root: Path, remote_host: str) -> None:
        self.batch_root = Path(batch_root).resolve()
        self.cache_root = Path(cache_root).resolve()
        self.remote_host = remote_host
        self.paths = {key: self.batch_root / name for key, name in _LOCAL_FILES.items()}
        self.index = _json(self.paths["index"])
        self.packaging = _json(self.paths["packaging"])
        self.verification = _json(self.paths["verification"])
        self._validate_local_binding()
        self.exports = self._load_exports()

    def _validate_local_binding(self) -> None:
        if self.index.get("status") != "complete":
            raise RuntimeError("冻结批次尚未完成")
        if self.packaging.get("status") != "completed" or self.verification.get("status") != "passed":
            raise RuntimeError("冻结批次交付校验尚未通过")
        expected_index = self.packaging.get("index_sha256")
        expected_report = self.packaging.get("html_sha256")
        if expected_index != _sha256(self.paths["index"]):
            raise RuntimeError("冻结交付索引 SHA-256 不匹配")
        if expected_report != _sha256(self.paths["report"]):
            raise RuntimeError("冻结结果报告 SHA-256 不匹配")

    def _load_exports(self) -> dict[str, dict[str, Any]]:
        records = self.index.get("remote_exports")
        if not isinstance(records, list) or len(records) != len(_EXPORT_FILES):
            raise RuntimeError("冻结索引必须包含9个导出文件")
        found: dict[str, dict[str, Any]] = {}
        for raw in records:
            if not isinstance(raw, dict):
                raise RuntimeError("冻结导出记录格式错误")
            name = raw.get("name")
            if name not in _EXPORT_FILES or name in found:
                raise RuntimeError("冻结索引含未知或重复导出名")
            expected_path = _REMOTE_EXPORT_DIR / _EXPORT_FILES[name]
            if raw.get("path") != str(expected_path):
                raise RuntimeError(f"冻结导出路径不符合白名单：{name}")
            size = raw.get("size")
            sha = raw.get("sha256")
            if not isinstance(size, int) or size < 0 or not isinstance(sha, str) or len(sha) != 64:
                raise RuntimeError(f"冻结导出校验信息无效：{name}")
            found[name] = {**raw, "filename": _EXPORT_FILES[name]}
        if set(found) != set(_EXPORT_FILES):
            raise RuntimeError("冻结索引导出集合不完整")
        return found

    def _summary(self) -> dict[str, Any]:
        totals = dict(self.index.get("totals") or {})
        return {
            "id": COMPLETED_BATCH_ID,
            "title": _TITLE,
            "status": "completed",
            "completed_samples": totals.get("files", 149),
            "total_samples": 149,
            "paired_blank_files": 102,
            "paired_groups": self.index.get("paired_groups", 50),
            "totals": totals,
            "report_url": f"/api/completed-batches/{COMPLETED_BATCH_ID}/report",
            "detail_url": f"/api/completed-batches/{COMPLETED_BATCH_ID}",
            "download_all_url": f"/api/completed-batches/{COMPLETED_BATCH_ID}/download-all",
        }

    def list_payload(self) -> dict[str, Any]:
        return {"batches": [self._summary()]}

    def detail_payload(self) -> dict[str, Any]:
        native = _json(self.paths["native_state"])
        conclusions = _json(self.paths["conclusions_state"])
        native_samples = native.get("samples") or {}
        conclusion_samples = conclusions.get("samples") or {}
        if (
            native.get("state") != "completed"
            or conclusions.get("state") != "completed"
            or len(native_samples) != 149
            or set(native_samples) != set(conclusion_samples)
        ):
            raise RuntimeError("冻结批次的149份原生结果与结论状态不完整或不一致")
        samples = []
        for sample_key in sorted(conclusion_samples):
            entry = conclusion_samples[sample_key]
            summary = dict(entry.get("summary") or {})
            source = dict(summary.get("source_sample") or {})
            native_entry = dict(native_samples.get(sample_key) or {})
            samples.append({
                "sample_key": sample_key,
                "task_id": hashlib.sha256((COMPLETED_BATCH_ID + ":" + sample_key).encode()).hexdigest()[:12],
                "native_conclusion_records": summary.get("compound_result_records", 0),
                "sample_id": source.get("sample_id"),
                "polarity": source.get("polarity"),
                "replicate": source.get("replicate"),
                "searched_ms2_spectra": summary.get("parsed_searched_ms2_spectra", native_entry.get("spectra_completed", 0)),
                "matched_ms2_spectra": summary.get("matched_ms2_spectra", native_entry.get("matched_spectra", 0)),
                "ms1_features": summary.get("ms1_feature_records", native_entry.get("ms1_features", 0)),
                "xbt_named_candidate_rows": summary.get("xbt_named_candidate_rows", 0),
                "blank_flagged_ms2": summary.get("blank_flagged_ms2", 0),
                "blank_clean_view_excluded_ms2": summary.get("blank_clean_view_excluded_ms2", 0),
                "ambiguity_group_rows": summary.get("ambiguity_group_rows", 0),
                "background_qc_rows": summary.get("background_qc_rows", 0),
                "external_qc_rows": summary.get("external_review_rows", 0),
                "name_qc_rows": summary.get("database_title_or_unnamed_qc_rows", 0),
                "unlinked_qc_rows": summary.get("xbt_link_review_qc_rows", 0),
                "ms1_evidence_counts": summary.get("ms1_evidence_level_counts", {}),
            })
        export_rows = self.index.get("table_row_counts") or {}
        exports = []
        for key, record in self.exports.items():
            is_csv = record["filename"].endswith(".csv.xz")
            exports.append({
                "name": key,
                "filename": record["filename"],
                "path": record["path"],
                "label": _EXPORT_LABELS[key],
                "size": record["size"],
                "sha256": record["sha256"],
                "rows": export_rows.get(record["filename"]),
                "download_url": f"/api/completed-batches/{COMPLETED_BATCH_ID}/exports/{key}",
                "csv_url": (
                    f"/api/completed-batches/{COMPLETED_BATCH_ID}/exports/{key}?decompress=true"
                    if is_csv else None
                ),
            })
        return {
            **self._summary(),
            "ms1_evidence_counts": self.index.get("ms1_evidence_counts", {}),
            "declared_vs_searched_ms2": self.index.get("declared_vs_searched_ms2", {}),
            "msclassify": self.index.get("msclassify", {}),
            "xbt_dataset": self.index.get("xbt_dataset"),
            "exports": exports,
            "samples": samples,
            "index_url": f"/api/completed-batches/{COMPLETED_BATCH_ID}/index",
            "verification_url": f"/api/completed-batches/{COMPLETED_BATCH_ID}/verification",
            "counting_unit": self.index.get("counting_unit"),
            "models_written_by_workflow": False,
            "no_analysis_rerun": True,
        }

    def _cached_valid(self, record: dict[str, Any], target: Path) -> bool:
        return target.is_file() and target.stat().st_size == record["size"] and _sha256(target) == record["sha256"]

    def _copy_one(self, name: str) -> Path:
        record = self.exports[name]
        target = self.cache_root / record["filename"]
        with _cache_lock(target):
            _validate_cache_volume(self.cache_root)
            if self._cached_valid(record, target):
                return target
            if target.exists():
                target.unlink()
            self.cache_root.mkdir(parents=True, exist_ok=True)
            if shutil.disk_usage(self.cache_root).free - record["size"] < _CACHE_FLOOR_BYTES:
                raise RuntimeError("外置缓存卷剩余空间不足，复制后必须保留至少50 GiB")
            partial = target.with_name(f".{target.name}.part.{os.getpid()}.{threading.get_ident()}")
            try:
                result = subprocess.run(
                    ["scp", "-q", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12", "--", f"{self.remote_host}:{record['path']}", str(partial)],
                    check=False,
                    capture_output=True,
                    timeout=_COPY_TIMEOUT_SECONDS,
                )
                if result.returncode != 0:
                    raise RuntimeError("远端已完成结果暂时无法读取，请确认 msqc-server 连接后重试")
                if not self._cached_valid(record, partial):
                    raise RuntimeError(f"下载文件校验失败：{name} 的大小或 SHA-256 不匹配")
                os.replace(partial, target)
                return target
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError("读取远端已完成结果超时，请稍后重试") from exc
            finally:
                partial.unlink(missing_ok=True)

    def ensure_exports(self, names: list[str]) -> dict[str, Path]:
        missing_bytes = 0
        for name in names:
            record = self.exports[name]
            target = self.cache_root / record["filename"]
            if not self._cached_valid(record, target):
                missing_bytes += record["size"]
        if missing_bytes > _OPERATION_CAP_BYTES:
            raise RuntimeError("本次缓存新增量超过192 MiB限制")
        return {name: self._copy_one(name) for name in names}

    def table_payload(self, name: str, offset: int, limit: int, sample_key: str | None) -> dict[str, Any]:
        source = self.ensure_exports([name])[name]
        selected: list[dict[str, str]] = []
        total = self._expected_table_total(name, sample_key)
        if total is None:
            raise RuntimeError(f"冻结表缺少登记行数：{name}")
        matched = 0
        with lzma.open(source, "rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = list(reader.fieldnames or [])
            if sample_key is not None and "sample_key" not in columns:
                raise ValueError("该表不含 sample_key 字段，不能按样本筛选")
            for row in reader:
                if sample_key is not None and row.get("sample_key") != sample_key:
                    continue
                if offset <= matched < offset + limit:
                    selected.append(dict(row))
                matched += 1
                if len(selected) == limit:
                    break
        return {
            "columns": columns,
            "rows": selected,
            "total": total,
            "offset": offset,
            "limit": limit,
            "has_more": offset + len(selected) < total,
            "table_key": name,
            "sample_key": sample_key,
        }

    def _expected_table_total(self, name: str, sample_key: str | None) -> int | None:
        if sample_key is None:
            filename = self.exports[name]["filename"]
            value = (self.index.get("table_row_counts") or {}).get(filename)
            return int(value) if isinstance(value, int) else None
        conclusions = (_json(self.paths["conclusions_state"]).get("samples") or {}).get(sample_key)
        if conclusions is None:
            return 0
        summary = conclusions.get("summary") or {}
        if name == "file_manifest":
            return 1
        if name == "xbt_csv":
            return int(summary.get("xbt_named_candidate_rows", 0))
        if name == "msms_table":
            return int(summary.get("parsed_searched_ms2_spectra", 0))
        if name == "ms1_table":
            return int(summary.get("ms1_feature_records", 0))
        if name == "qc_table":
            return sum(int(summary.get(key, 0)) for key in (
                "xbt_link_review_qc_rows", "external_review_rows", "background_qc_rows", "database_title_or_unnamed_qc_rows"
            ))
        return None


class _ZipCancelled(Exception):
    pass


class _QueueSink(io.RawIOBase):
    def __init__(self, output: queue.Queue[Any], cancelled: threading.Event) -> None:
        self.output = output
        self.cancelled = cancelled
        self.position = 0

    def writable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return False

    def tell(self) -> int:
        return self.position

    def write(self, data: bytes | bytearray) -> int:
        chunk = bytes(data)
        self.position += len(chunk)
        while chunk and not self.cancelled.is_set():
            try:
                self.output.put(chunk, timeout=0.1)
                break
            except queue.Full:
                continue
        if self.cancelled.is_set():
            raise _ZipCancelled
        return len(chunk)


def _zip_stream(members: list[tuple[Path, str]]) -> Iterator[bytes]:
    output: queue.Queue[Any] = queue.Queue(maxsize=8)
    cancelled = threading.Event()
    done = object()

    def put_control(item: Any) -> None:
        while not cancelled.is_set():
            try:
                output.put(item, timeout=0.1)
                return
            except queue.Full:
                continue

    def produce() -> None:
        try:
            sink = _QueueSink(output, cancelled)
            with zipfile.ZipFile(sink, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as archive:
                for path, arcname in members:
                    if cancelled.is_set():
                        raise _ZipCancelled
                    info = zipfile.ZipInfo(arcname)
                    info.compress_type = zipfile.ZIP_STORED
                    info.flag_bits |= 0x800
                    with archive.open(info, "w", force_zip64=True) as destination, path.open("rb") as source:
                        shutil.copyfileobj(source, destination, length=1024 * 1024)
        except _ZipCancelled:
            return
        except BaseException as exc:  # propagated by the response iterator
            put_control(exc)
        finally:
            put_control(done)

    producer = threading.Thread(target=produce, name="completed-batch-zip", daemon=True)
    producer.start()
    try:
        while True:
            item = output.get()
            if item is done:
                return
            if isinstance(item, BaseException):
                raise item
            yield item
    finally:
        cancelled.set()
        while True:
            try:
                output.get_nowait()
            except queue.Empty:
                break
        producer.join(timeout=2.0)


def register_completed_batch_routes(
    app: Any,
    *,
    batch_root: Path,
    cache_root: Path,
    remote_host: str = "msqc-server",
) -> CompletedBatchService | None:
    """Register read-only completed-batch endpoints and return their service."""
    batch_root = Path(batch_root).resolve()
    required_paths = [batch_root / name for name in _LOCAL_FILES.values()]
    present = [path.exists() for path in required_paths]
    service: CompletedBatchService | None = None
    service_error: str | None = None
    if any(present):
        try:
            service = CompletedBatchService(batch_root, cache_root, remote_host)
        except (OSError, RuntimeError) as exc:
            service_error = str(exc) or "冻结交付文件无法读取"

    def require_service(*, list_route: bool = False) -> CompletedBatchService | None:
        if service is not None:
            return service
        if service_error is not None:
            raise HTTPException(status_code=502, detail=f"已完成批次交付无效：{service_error}")
        if list_route:
            return None
        raise HTTPException(status_code=404, detail="已完成批次尚未安装或不可用")

    def require_batch(batch_id: str) -> None:
        if batch_id != COMPLETED_BATCH_ID:
            raise HTTPException(status_code=404, detail="已完成批次不存在")

    def require_export(name: str, *, table: bool = False) -> None:
        current = require_service()
        assert current is not None
        allowed = _TABLE_KEYS if table else set(current.exports)
        if name not in allowed:
            raise HTTPException(status_code=404, detail="冻结导出文件不存在")

    async def run_blocking(function: Any, *args: Any) -> Any:
        try:
            return await asyncio.to_thread(function, *args)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        except RuntimeError as exc:
            raise HTTPException(status_code=502, detail=str(exc)) from exc

    @app.get("/api/completed-batches")
    async def completed_batches() -> dict[str, Any]:
        current = require_service(list_route=True)
        return current.list_payload() if current is not None else {"batches": []}

    @app.get("/api/completed-batches/{batch_id}")
    async def completed_batch_detail(batch_id: str) -> dict[str, Any]:
        require_batch(batch_id)
        current = require_service()
        assert current is not None
        return await run_blocking(current.detail_payload)

    @app.get("/api/completed-batches/{batch_id}/report")
    async def completed_batch_report(batch_id: str, download: bool = False) -> FileResponse:
        require_batch(batch_id)
        current = require_service()
        assert current is not None
        return FileResponse(
            current.paths["report"],
            media_type="text/html; charset=utf-8",
            headers={"Content-Disposition": _disposition("宣痹汤_XBT_完整结果.html", inline=not download), "Cache-Control": "no-store"},
        )

    @app.get("/api/completed-batches/{batch_id}/index")
    @app.get("/api/completed-batches/{batch_id}/FINAL_DELIVERY_INDEX.json")
    async def completed_batch_index(batch_id: str) -> FileResponse:
        require_batch(batch_id)
        current = require_service()
        assert current is not None
        return FileResponse(current.paths["index"], media_type="application/json", headers={"Cache-Control": "no-store"})

    @app.get("/api/completed-batches/{batch_id}/verification")
    async def completed_batch_verification(batch_id: str) -> FileResponse:
        require_batch(batch_id)
        current = require_service()
        assert current is not None
        return FileResponse(current.paths["verification"], media_type="application/json", headers={"Cache-Control": "no-store"})

    @app.get("/api/completed-batches/{batch_id}/exports/{export_name}")
    async def completed_batch_export(batch_id: str, export_name: str, decompress: bool = False):
        require_batch(batch_id)
        require_export(export_name)
        current = require_service()
        assert current is not None
        record = current.exports[export_name]
        if decompress:
            if not record["filename"].endswith(".csv.xz"):
                raise HTTPException(status_code=422, detail="仅 CSV.XZ 导出支持解压为 CSV")
        path = (await run_blocking(current.ensure_exports, [export_name]))[export_name]
        if decompress:

            def csv_chunks() -> Iterator[bytes]:
                yield b"\xef\xbb\xbf"
                with lzma.open(path, "rb") as source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        yield chunk

            filename = record["filename"][:-3]
            return StreamingResponse(csv_chunks(), media_type="text/csv; charset=utf-8", headers={"Content-Disposition": _disposition(filename), "Cache-Control": "no-store"})
        return FileResponse(path, media_type="application/x-xz" if path.suffix == ".xz" else "application/json", filename=record["filename"], headers={"Cache-Control": "no-store"})

    @app.get("/api/completed-batches/{batch_id}/tables/{export_name}")
    async def completed_batch_table(
        batch_id: str,
        export_name: str,
        offset: int = Query(0, ge=0),
        limit: int = Query(50, ge=1, le=200),
        sample_key: str | None = None,
    ) -> dict[str, Any]:
        require_batch(batch_id)
        require_export(export_name, table=True)
        current = require_service()
        assert current is not None
        return await run_blocking(current.table_payload, export_name, offset, limit, sample_key)

    @app.get("/api/completed-batches/{batch_id}/download-all")
    async def completed_batch_download_all(batch_id: str):
        require_batch(batch_id)
        current = require_service()
        assert current is not None
        cached = await run_blocking(current.ensure_exports, list(_EXPORT_FILES))
        members = [(cached[name], current.exports[name]["filename"]) for name in _EXPORT_FILES]
        members.extend([
            (current.paths["report"], "宣痹汤_XBT_完整结果.html"),
            (current.paths["index"], "FINAL_DELIVERY_INDEX.json"),
            (current.paths["packaging"], "FINAL_PACKAGING_RECEIPT.json"),
            (current.paths["verification"], "FINAL_DELIVERY_VERIFICATION.json"),
        ])
        return StreamingResponse(
            _zip_stream(members),
            media_type="application/zip",
            headers={"Content-Disposition": _disposition("宣痹汤_XBT_149份完整冻结交付.zip"), "Cache-Control": "no-store"},
        )

    return service


__all__ = ["COMPLETED_BATCH_ID", "CompletedBatchService", "register_completed_batch_routes"]
