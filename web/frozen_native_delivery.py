"""Read-only access to the verified 149-sample native-conclusion delivery."""

from __future__ import annotations

import asyncio
import csv
import hashlib
import json
import lzma
import os
import shutil
import subprocess
import threading
from pathlib import Path
from typing import Any

from fastapi import HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse

from web.completed_batches import COMPLETED_BATCH_ID, _disposition, _zip_stream


REMOTE_ROOT = Path("/data/Mssystem/xbt_real_batch_20261021_ms2/delivery_access_audit_20261022")
REMOTE_FILES = {
    "native_compounds": REMOTE_ROOT / "all_native_conclusions.csv.xz",
    "report": REMOTE_ROOT / "FULL_NATIVE_RESULTS_149.html",
    "index": REMOTE_ROOT / "FULL_NATIVE_DELIVERY_INDEX.json",
}
SCHEMA = "mssystem.frozen-native-delivery/1.0"
_CACHE_FLOOR = 50 * 1024**3
_MAX_COPY = 320 * 1024**2


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _semantic_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
    ).encode()).hexdigest()


class FrozenNativeDelivery:
    def __init__(self, bridge):
        self.bridge = bridge
        self.cache_root = Path(bridge.service.cache_root).resolve() / "native_delivery"
        self._index = None
        self._index_raw = None
        self._index_sha256 = None
        self._lock = threading.Lock()

    def _ssh_bytes(self, path: Path, *, timeout=180) -> bytes:
        if path not in REMOTE_FILES.values():
            raise RuntimeError("冻结原生交付路径不在白名单")
        completed = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
             self.bridge.service.remote_host, f"cat -- {path}"],
            capture_output=True, timeout=timeout,
        )
        if completed.returncode:
            raise RuntimeError("冻结原生交付尚未生成或远端读取失败")
        return completed.stdout

    def _validate_index(self, value: dict) -> dict:
        if value.get("schema_version") != SCHEMA or value.get("status") != "completed":
            raise RuntimeError("冻结原生交付索引状态或版本无效")
        semantic = value.get("semantic_sha256")
        if semantic != _semantic_sha256({k: v for k, v in value.items() if k != "semantic_sha256"}):
            raise RuntimeError("冻结原生交付索引语义SHA不匹配")
        original = self.bridge.service.paths["index"]
        if value.get("original_delivery_index_sha256") != _sha256(original):
            raise RuntimeError("冻结原生交付未绑定当前原始交付索引")
        samples = value.get("samples")
        if not isinstance(samples, list) or len(samples) != 149:
            raise RuntimeError("冻结原生交付索引必须包含149份样本")
        expected = {entry["sample_key"]: (task_id, entry)
                    for task_id, entry in self.bridge.entries.items()}
        if set(row.get("sample_key") for row in samples) != set(expected):
            raise RuntimeError("冻结原生交付样本集合与任务登记不一致")
        seen_tasks = set()
        for row in samples:
            task_id, entry = expected[row["sample_key"]]
            seen_tasks.add(row.get("task_id"))
            if (row.get("task_id") != task_id
                    or row.get("source_sha256") != entry["source_sha256"]
                    or row.get("native_result") != entry["results"]
                    or row.get("conclusion_result") != entry["conclusions"]):
                raise RuntimeError("冻结原生交付样本输入绑定不一致")
        if seen_tasks != set(self.bridge.tasks):
            raise RuntimeError("冻结原生交付任务ID集合不完整")
        artifacts = value.get("artifacts") or {}
        for key in ("native_compounds", "report"):
            record = artifacts.get(key) or {}
            if (record.get("path") != str(REMOTE_FILES[key])
                    or not isinstance(record.get("size"), int)
                    or record["size"] < 1 or record["size"] > _MAX_COPY
                    or not isinstance(record.get("sha256"), str)
                    or len(record["sha256"]) != 64):
                raise RuntimeError("冻结原生交付产物绑定无效：" + key)
        if value.get("sample_count") != 149 or value.get("counting_unit") != "record_occurrences":
            raise RuntimeError("冻结原生交付计数单位或样本数无效")
        if value.get("standard_confirmed_count") != 0:
            raise RuntimeError("冻结原生交付标准品计数与冻结批次不一致")
        return value

    def index(self, *, refresh=False) -> dict:
        with self._lock:
            if self._index is not None and not refresh:
                return self._index
            try:
                raw = self._ssh_bytes(REMOTE_FILES["index"])
                value = json.loads(raw.decode("utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                raise RuntimeError("冻结原生交付索引无法解析") from exc
            self._index = self._validate_index(value)
            self._index_raw = raw
            self._index_sha256 = hashlib.sha256(raw).hexdigest()
            return self._index

    def detail(self) -> dict:
        index = self.index()
        artifacts = index["artifacts"]
        base = f"/api/completed-batches/{COMPLETED_BATCH_ID}"
        return {
            "id": COMPLETED_BATCH_ID,
            "status": "completed",
            "sample_count": 149,
            "native_conclusion_rows": index["native_conclusion_rows"],
            "record_type_counts": index["record_type_counts"],
            "conclusion_status_counts": index["conclusion_status_counts"],
            "counting_unit": "record_occurrences",
            "standard_confirmed_count": 0,
            "original_delivery_index_sha256": index["original_delivery_index_sha256"],
            "artifacts": artifacts,
            "report_url": base + "/native-report",
            "index_url": base + "/native-index",
            "native_compounds_xz_url": base + "/native-exports/native_compounds",
            "native_compounds_csv_url": base + "/native-exports/native_compounds?decompress=true",
            "native_table_url": base + "/native-table",
            "download_all_url": base + "/native-download-all",
            "no_analysis_rerun": True,
        }

    def index_artifact(self) -> Path:
        self.index()
        assert self._index_raw is not None and self._index_sha256 is not None
        target = self.cache_root / REMOTE_FILES["index"].name
        with self._lock:
            if (target.is_file() and target.stat().st_size == len(self._index_raw)
                    and _sha256(target) == self._index_sha256):
                return target
            resolved = target.resolve()
            try:
                volume = Path("/Volumes") / resolved.relative_to("/Volumes").parts[0]
            except (ValueError, IndexError) as exc:
                raise RuntimeError("冻结原生交付索引缓存必须位于/Volumes外置卷") from exc
            if not volume.is_mount() or shutil.disk_usage(volume).free - len(self._index_raw) < _CACHE_FLOOR:
                raise RuntimeError("冻结原生交付索引缓存卷未挂载或空间不足")
            self.cache_root.mkdir(parents=True, exist_ok=True)
            partial = target.with_name(target.name + f".partial.{os.getpid()}")
            try:
                partial.write_bytes(self._index_raw)
                if _sha256(partial) != self._index_sha256:
                    raise RuntimeError("冻结原生交付索引缓存SHA校验失败")
                os.replace(partial, target)
            finally:
                partial.unlink(missing_ok=True)
            return target

    def _cache_file(self, key: str) -> Path:
        if key not in {"native_compounds", "report"}:
            raise RuntimeError("未知冻结原生交付产物")
        record = self.index()["artifacts"][key]
        target = self.cache_root / Path(record["path"]).name
        with self._lock:
            if target.is_file() and target.stat().st_size == record["size"] and _sha256(target) == record["sha256"]:
                return target
            volume = Path("/Volumes") / target.resolve().relative_to("/Volumes").parts[0]
            if not volume.is_mount():
                raise RuntimeError("冻结原生交付缓存卷未挂载")
            if record["size"] > _MAX_COPY or shutil.disk_usage(volume).free - record["size"] < _CACHE_FLOOR:
                raise RuntimeError("冻结原生交付缓存空间不足或文件超过320MiB限制")
            self.cache_root.mkdir(parents=True, exist_ok=True)
            partial = target.with_name(target.name + f".partial.{os.getpid()}")
            try:
                with partial.open("wb") as handle:
                    completed = subprocess.run(
                        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
                         self.bridge.service.remote_host, f"cat -- {REMOTE_FILES[key]}"],
                        stdout=handle, stderr=subprocess.PIPE, timeout=240,
                    )
                if completed.returncode or partial.stat().st_size != record["size"] or _sha256(partial) != record["sha256"]:
                    raise RuntimeError("冻结原生交付缓存SHA/size校验失败")
                os.replace(partial, target)
            finally:
                partial.unlink(missing_ok=True)
            return target

    def artifact(self, key: str) -> Path:
        return self._cache_file(key)

    def table(self, *, offset: int, limit: int, sample_key: str | None) -> dict:
        if offset < 0 or limit < 1 or limit > 200:
            raise ValueError("冻结原生结论分页参数无效；limit最大200")
        index = self.index()
        if sample_key is None:
            total = int(index["native_conclusion_rows"])
        else:
            matching = [row for row in index["samples"] if row["sample_key"] == sample_key]
            if len(matching) != 1:
                raise ValueError("未知冻结样本键")
            total = int(matching[0]["native_conclusion_rows"])
        path = self.artifact("native_compounds")
        page = []
        columns = []
        selected_index = 0
        with lzma.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = list(reader.fieldnames or [])
            for row in reader:
                if sample_key is not None and row.get("sample_key") != sample_key:
                    continue
                if selected_index >= offset and len(page) < limit:
                    page.append(row)
                selected_index += 1
                if len(page) >= limit:
                    break
        return {
            "columns": columns, "rows": page, "total": total,
            "offset": offset, "limit": limit, "has_more": offset + len(page) < total,
            "table_key": "native_compounds", "sample_key": sample_key,
            "counting_unit": "record_occurrences", "no_analysis_rerun": True,
        }


def register_frozen_native_delivery_routes(app, bridge):
    """Register lazy, read-only routes for the verified full-native delivery."""
    if bridge is None:
        return None
    delivery = FrozenNativeDelivery(bridge)
    base = f"/api/completed-batches/{COMPLETED_BATCH_ID}"

    def fail(exc):
        return HTTPException(status_code=502, detail=str(exc))

    @app.get(base + "/native-delivery")
    async def native_delivery_detail():
        try:
            return await asyncio.to_thread(delivery.detail)
        except (RuntimeError, OSError, KeyError, ValueError) as exc:
            raise fail(exc)

    @app.get(base + "/native-index")
    async def native_delivery_index():
        try:
            return JSONResponse(await asyncio.to_thread(delivery.index))
        except (RuntimeError, OSError, KeyError, ValueError) as exc:
            raise fail(exc)

    @app.get(base + "/native-report")
    async def native_delivery_report():
        try:
            path = await asyncio.to_thread(delivery.artifact, "report")
            return HTMLResponse(path.read_text(encoding="utf-8"))
        except (RuntimeError, OSError, KeyError, ValueError) as exc:
            raise fail(exc)

    @app.get(base + "/native-exports/native_compounds")
    async def native_compounds_export(decompress: bool = False):
        try:
            path = await asyncio.to_thread(delivery.artifact, "native_compounds")
            if not decompress:
                def compressed():
                    with path.open("rb") as handle:
                        for block in iter(lambda: handle.read(1024 * 1024), b""):
                            yield block
                return StreamingResponse(compressed(), media_type="application/x-xz", headers={
                    "Content-Disposition": _disposition(path.name)})
            def csv_stream():
                with lzma.open(path, "rb") as handle:
                    for block in iter(lambda: handle.read(1024 * 1024), b""):
                        yield block
            return StreamingResponse(csv_stream(), media_type="text/csv; charset=utf-8", headers={
                "Content-Disposition": _disposition("all_native_conclusions.csv")})
        except (RuntimeError, OSError, KeyError, ValueError, lzma.LZMAError) as exc:
            raise fail(exc)

    @app.get(base + "/native-table")
    async def native_compounds_table(offset: int = 0, limit: int = 50, sample_key: str | None = None):
        try:
            return await asyncio.to_thread(
                delivery.table, offset=offset, limit=limit, sample_key=sample_key,
            )
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc))
        except (RuntimeError, OSError, KeyError, lzma.LZMAError) as exc:
            raise fail(exc)

    @app.get(base + "/native-download-all")
    async def native_delivery_download_all():
        try:
            legacy = await asyncio.to_thread(
                bridge.service.ensure_exports, list(bridge.service.exports),
            )
            native_csv, native_report, native_index = await asyncio.gather(
                asyncio.to_thread(delivery.artifact, "native_compounds"),
                asyncio.to_thread(delivery.artifact, "report"),
                asyncio.to_thread(delivery.index_artifact),
            )
            members = [(legacy[name], "legacy/" + bridge.service.exports[name]["filename"])
                       for name in bridge.service.exports]
            members.extend([
                (bridge.service.paths["report"], "legacy/既有XBT附加筛选摘要.html"),
                (bridge.service.paths["index"], "legacy/FINAL_DELIVERY_INDEX.json"),
                (bridge.service.paths["packaging"], "legacy/FINAL_PACKAGING_RECEIPT.json"),
                (bridge.service.paths["verification"], "legacy/FINAL_DELIVERY_VERIFICATION.json"),
                (native_csv, "all_native_conclusions.csv.xz"),
                (native_report, "FULL_NATIVE_RESULTS_149.html"),
                (native_index, "FULL_NATIVE_DELIVERY_INDEX.json"),
            ])
            return StreamingResponse(
                _zip_stream(members), media_type="application/zip",
                headers={"Content-Disposition": _disposition("宣痹汤_149份完整原生与附加筛选交付.zip")},
            )
        except (RuntimeError, OSError, KeyError, ValueError) as exc:
            raise fail(exc)

    return delivery


__all__ = ["FrozenNativeDelivery", "register_frozen_native_delivery_routes"]
