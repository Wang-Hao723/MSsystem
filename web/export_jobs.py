"""Persist batch archives independently of the browser's HTTP connection."""
from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor
import hashlib
import inspect
import shutil
import json
from pathlib import Path
import re
import resource
import sys
import threading
import time
import uuid
import zipfile

from web.export_identity import file_sha256

_context = threading.local()


def has_progress_callback():
    return bool(getattr(_context, "progress", None))


def report_progress(**fields):
    callback = getattr(_context, "progress", None)
    if callback:
        callback(**fields)


def rebind_request(request):
    """Replace a running job's provisional request with its final identity."""
    callback = getattr(_context, "rebind_request", None)
    if not callback:
        raise RuntimeError("Export request rebinding is only available inside a job")
    return callback(request)


def canonical_sql_signature(signature, root):
    """Bind the same dependency set before and after SQL initialization."""
    root = Path(root)
    config_path = root / "config/sql_metadata_repairs.json"
    if not config_path.is_file() or len(signature) != 2:
        return signature
    from runtime.verified_sql_retrieval import _path, _signature
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not config.get("enabled"):
        return signature
    registry = json.loads(_path(root, config["registry_path"]).read_text(encoding="utf-8"))
    paths = [Path(item[0]) for item in signature[1]]
    paths.extend(Path(name) for name in registry["source_files"])
    database = _path(root, config["source_database_path"])
    paths.extend(Path(str(database) + suffix) for suffix in ("-wal", "-journal"))
    return signature[0], _signature(list(dict.fromkeys(paths)))


class ExportJobs:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock = threading.RLock()
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="export-job")
        self.active = {}
        self.completed = {}
        self.owned = set()
        self._restore_completed_index()

    @staticmethod
    def _request_key(request):
        canonical = {key: value for key, value in request.items()
                     if not str(key).startswith("_")}
        canonical["task_ids"] = sorted(set(request["task_ids"]))
        encoded = json.dumps(
            canonical, sort_keys=True, separators=(",", ":"), default=str,
        ).encode()
        return hashlib.sha256(encoded).hexdigest()

    def _path(self, job_id):
        if not re.fullmatch(r"[a-f0-9]{32}", job_id):
            raise KeyError("Invalid export identifier")
        return self.root / f"{job_id}.json"

    def _write(self, state):
        path = self._path(state["id"])
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(state, ensure_ascii=False), encoding="utf-8")
        temporary.replace(path)

    def _artifact_path(self, job_id):
        return self.root / f"{job_id}.zip"

    @staticmethod
    def _resources():
        usage = resource.getrusage(resource.RUSAGE_SELF)
        return {"process_cpu_seconds": time.process_time(),
                "max_rss_kib": int(usage.ru_maxrss // 1024
                                   if sys.platform == "darwin"
                                   else usage.ru_maxrss)}

    def _validate_completed(self, state, *, verify_container=True):
        if state.get("status") != "completed" or not state.get("artifact_sha256"):
            return False
        artifact = self._artifact_path(state["id"])
        try:
            if (not artifact.is_file()
                    or artifact.stat().st_size != int(state.get("size_bytes") or -1)
                    or file_sha256(artifact) != state["artifact_sha256"]):
                return False
            if verify_container:
                media = str(state.get("media_type") or "")
                native = state.get("request", {}).get("mode") == "single"
                if (not native or media in {
                    "application/zip",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                }):
                    with zipfile.ZipFile(artifact) as archive:
                        if archive.testzip():
                            return False
                        if (native and media.endswith("sheet")
                                and "xl/workbook.xml" not in archive.namelist()):
                            return False
                elif media == "application/pdf":
                    with artifact.open("rb") as source:
                        if source.read(5) != b"%PDF-":
                            return False
        except (OSError, ValueError, KeyError, zipfile.BadZipFile):
            return False
        return True

    def _adopt_legacy_download_hash(self, state):
        """Validate an old file for download without making it cache-reusable."""
        if state.get("status") != "completed":
            return False
        artifact = self._artifact_path(state["id"])
        try:
            if not artifact.is_file():
                return False
            declared_size = state.get("size_bytes")
            if declared_size is not None and artifact.stat().st_size != int(declared_size):
                return False
            media = str(state.get("media_type") or "application/zip")
            native = state.get("request", {}).get("mode") == "single"
            if (not native or media in {
                "application/zip",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            }):
                with zipfile.ZipFile(artifact) as archive:
                    if archive.testzip():
                        return False
                    if (native and media.endswith("sheet")
                            and "xl/workbook.xml" not in archive.namelist()):
                        return False
            elif media == "application/pdf":
                with artifact.open("rb") as source:
                    if source.read(5) != b"%PDF-":
                        return False
            state.update(
                artifact_sha256=file_sha256(artifact),
                size_bytes=artifact.stat().st_size,
                cache_reusable=False,
                legacy_hash_adopted_at=time.time(),
            )
            self._write(state)
            return True
        except (OSError, ValueError, KeyError, zipfile.BadZipFile):
            return False

    def _restore_completed_index(self):
        """Recover only SHA-verified completed artifacts after a restart."""
        for path in self.root.glob("*.json"):
            try:
                state = json.loads(path.read_text(encoding="utf-8"))
                key = self._request_key(state["request"])
            except (OSError, ValueError, TypeError, KeyError):
                continue
            if state.get("cache_reusable", True) and self._validate_completed(state):
                previous_id = self.completed.get(key)
                if previous_id is None:
                    self.completed[key] = state["id"]
                    continue
                try:
                    previous = json.loads(self._path(previous_id).read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    previous = {}
                if float(state.get("updated_at") or 0) > float(previous.get("updated_at") or 0):
                    self.completed[key] = state["id"]

    def get(self, job_id):
        with self.lock:
            try:
                state = json.loads(self._path(job_id).read_text(encoding="utf-8"))
            except FileNotFoundError:
                raise KeyError("Export not found") from None
            if state["status"] in {"queued", "running"} and job_id not in self.owned:
                state.update(status="failed", error="服务重启中断了未完成的导出，请重新生成；原分析结果仍保留。")
                self._write(state)
            return state

    def submit(self, request, runner):
        key = self._request_key(request)
        with self.lock:
            if key in self.active:
                return self.get(self.active[key])
            completed_id = self.completed.get(key)
            if completed_id is not None:
                validation_started = time.perf_counter()
                try:
                    completed_state = self.get(completed_id)
                except KeyError:
                    completed_state = None
                if completed_state is not None and self._validate_completed(completed_state):
                    return {**completed_state, "cache_hit": True,
                            "cache_validation_wall_seconds":
                                time.perf_counter() - validation_started}
                self.completed.pop(key, None)
            job_id = uuid.uuid4().hex
            state = dict(id=job_id, status="queued", completed=0,
                         total=len(request["task_ids"]), current_task=None,
                         created_at=time.time(), request=request, cache_hit=False)
            self._write(state)
            self.owned.add(job_id)
            self.active[key] = job_id
            self.executor.submit(self._run, key, state, runner)
            return dict(state)

    def _run(self, key, state, runner):
        active_keys = {key}
        final_key = [key]

        def update(**fields):
            with self.lock:
                phase_metric = fields.get("export_profile")
                if isinstance(phase_metric, dict):
                    history = list(state.get("phase_metrics") or [])
                    history.append(dict(phase_metric))
                    fields["phase_metrics"] = history
                state.update(fields, updated_at=time.time())
                self._write(state)

        def rebind(request):
            new_key = self._request_key(request)
            with self.lock:
                owner = self.active.get(new_key)
                if owner not in {None, state["id"]}:
                    raise RuntimeError("An identical final export is already active")
                previous_signature = state.get("request", {}).get("source_signature")
                state["request"] = dict(request)
                state["initial_source_signature"] = previous_signature
                state["final_source_signature"] = request.get("source_signature")
                state["updated_at"] = time.time()
                self._write(state)
                self.active[new_key] = state["id"]
                active_keys.add(new_key)
                final_key[0] = new_key
            return dict(state)
        temporary = self.root / f"{state['id']}.zip.part"
        run_started = time.perf_counter()
        start_resources = self._resources()
        try:
            update(status="running", phase="preparing")
            _context.progress = update
            _context.rebind_request = rebind
            response = runner()
            if inspect.isawaitable(response):
                response = asyncio.run(response)
            if getattr(response, "status_code", 200) >= 400:
                raise ValueError("Export renderer returned an error")
            render_finished = time.perf_counter()
            update(phase="saving", render_wall_seconds=render_finished - run_started)

            async def save():
                with temporary.open("wb") as output:
                    if hasattr(response, "body_iterator"):
                        async for chunk in response.body_iterator:
                            output.write(chunk.encode("utf-8") if isinstance(chunk, str) else chunk)
                    elif hasattr(response, "path"):
                        with open(response.path, "rb") as source:
                            shutil.copyfileobj(source, output)
                    else:
                        output.write(response.body)
                if response.background:
                    await response.background()

            asyncio.run(save())
            persisted_at = time.perf_counter()
            if not temporary.stat().st_size:
                raise ValueError("导出文件为空")
            native = state["request"].get("mode") == "single"
            media = response.headers.get("content-type", "application/zip").split(";")[0]
            if not native or media in {"application/zip", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}:
                with zipfile.ZipFile(temporary) as archive:
                    bad = archive.testzip()
                    if bad:
                        raise ValueError("Invalid ZIP member: " + bad)
                    if native and media.endswith("sheet") and "xl/workbook.xml" not in archive.namelist():
                        raise ValueError("Invalid Excel workbook")
            if media == "application/pdf":
                with temporary.open("rb") as source:
                    if source.read(5) != b"%PDF-":
                        raise ValueError("PDF renderer did not return a PDF")
            artifact = self._artifact_path(state["id"])
            temporary.replace(artifact)
            artifact_sha256 = file_sha256(artifact)
            end_resources = self._resources()
            update(status="completed", phase="ready", headers=dict(response.headers),
                   media_type=media, completed=state["total"],
                   exported=state.get("exported", state["total"]),
                   size_bytes=artifact.stat().st_size,
                   artifact_sha256=artifact_sha256, cache_reusable=True,
                   persist_wall_seconds=time.perf_counter() - persisted_at,
                   total_wall_seconds=time.perf_counter() - run_started,
                   process_cpu_seconds=(end_resources["process_cpu_seconds"]
                                        - start_resources["process_cpu_seconds"]),
                   max_rss_kib=end_resources["max_rss_kib"])
            with self.lock:
                self.completed[final_key[0]] = state["id"]
        except Exception as exc:
            temporary.unlink(missing_ok=True)
            update(status="failed", phase="failed", error=str(getattr(exc, "detail", exc)))
        finally:
            _context.progress = None
            _context.rebind_request = None
            with self.lock:
                for active_key in active_keys:
                    if self.active.get(active_key) == state["id"]:
                        self.active.pop(active_key, None)

    def archive_path(self, job_id):
        state = self.get(job_id)
        if state["status"] != "completed":
            raise ValueError("Export is not ready")
        if not state.get("artifact_sha256"):
            if not self._adopt_legacy_download_hash(state):
                raise ValueError("Historical export artifact failed validation")
            state = self.get(job_id)
        if not self._validate_completed(state):
            raise ValueError("Export artifact failed its saved SHA256 or format validation")
        return self._artifact_path(job_id)
