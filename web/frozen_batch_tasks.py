"""Expose completed, externally executed V16 tasks through the native task UI.

Only task metadata is registered locally. Results remain in their frozen
remote artifacts; the reader verifies their original bindings on every fetch.
No saved-result migration, scorer, sidecar builder or analysis queue is called.
"""

from __future__ import annotations

import asyncio
import csv
import hashlib
import io
import json
import lzma
import shlex
import subprocess
import threading
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path

from fastapi import HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, Response, StreamingResponse

from web.completed_batches import COMPLETED_BATCH_ID, CompletedBatchService, _json, _sha256

_NATIVE_PROTOCOL = "1d09e21e37dad4c856a70c6b09a8ccecf611aa5cadc40567f8acd82767fbe779"
_CACHE_BYTES = 32 * 1024**2
_READERS = threading.BoundedSemaphore(2)


def _frozen_run_configuration(service: CompletedBatchService) -> tuple[dict, dict]:
    """Resolve the saved target and blank inputs without inventing defaults."""
    manifest_path = getattr(
        service, "input_manifest_path",
        service.batch_root.parent / "20261021_new_real_samples_preflight" / "INPUT_MANIFEST.json",
    )
    pairing_path = getattr(
        service, "blank_pairing_path",
        service.batch_root.parent / "20261021_new_real_samples_preflight" / "BLANK_PAIRING.json",
    )
    try:
        manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
        pairings = json.loads(Path(pairing_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError("冻结输入清单或空白配对记录无法读取") from exc
    if not isinstance(manifest, list) or not isinstance(pairings, list):
        raise RuntimeError("冻结输入清单或空白配对记录格式错误")
    inputs = {str(row.get("path") or ""): row for row in manifest}
    if len(inputs) != len(manifest) or "" in inputs:
        raise RuntimeError("冻结输入清单路径缺失或重复")
    by_sample = {}
    for pairing in pairings:
        blanks = list(pairing.get("blanks") or [])
        for source_path in pairing.get("samples") or []:
            if source_path in by_sample:
                raise RuntimeError("冻结样品存在重复空白配对记录")
            by_sample[source_path] = {**pairing, "blanks": blanks}
    target = {
        "selection_source": "user_explicit",
        "selected_datasets": [service.index["xbt_dataset"]],
        "selection_changes_saved_results": False,
    }
    return {"inputs": inputs, "by_sample": by_sample}, target


class FrozenBatchTasks:
    def __init__(self, service: CompletedBatchService):
        self.service = service
        binding = _json(service.batch_root / "V16_INTEGRATION_STATE.json")["baseline_frozen_artifact_sha256"]
        for name in ("STATE.json", "CONCLUSIONS_STATE.json", "CONCLUSIONS_TERMINAL_RECEIPT.json"):
            if _sha256(service.batch_root / name) != binding[name]:
                raise RuntimeError("冻结任务登记文件SHA不匹配：" + name)
        native = _json(service.paths["native_state"])
        conclusions = _json(service.paths["conclusions_state"])
        service.detail_payload()  # completion and matching 149-file sets
        self.entries = {}
        self.tasks = {}
        self._legacy_tasks = {}
        self._legacy_registration_sha256 = {}
        self._lock = threading.Lock()
        self._cache = OrderedDict()
        self._cache_size = 0
        self._read_locks = {}
        blank_sources, target_selection = _frozen_run_configuration(service)
        for key, conclusion in sorted(conclusions["samples"].items()):
            row = native["samples"][key]
            if row["status"] != "completed" or conclusion["status"] != "completed":
                raise RuntimeError("Only completed frozen samples may enter the task registry")
            if (conclusion["native_remote_sha256"] != row["verified_remote_sha256"]
                    or conclusion["native_remote_path"] != row["verified_remote_path"]):
                raise RuntimeError("Frozen native/conclusion registry binding mismatch")
            task_id = hashlib.sha256((COMPLETED_BATCH_ID + ":" + key).encode()).hexdigest()[:12]
            summary = conclusion["summary"]
            total = int(summary["parsed_searched_ms2_spectra"])
            matched = int(summary["matched_ms2_spectra"])
            source = summary["source_sample"]
            source_path = row["source_path"]
            source_input = blank_sources["inputs"].get(source_path)
            pairing = blank_sources["by_sample"].get(source_path)
            if not source_input or source_input.get("role") != "sample" or not pairing:
                raise RuntimeError("冻结样品缺少可核验的输入或空白配对记录")
            blank_rows = [blank_sources["inputs"].get(path) for path in pairing["blanks"]]
            if any(not item or item.get("role") != "blank" for item in blank_rows):
                raise RuntimeError("冻结空白配对引用了未登记的空白文件")
            blank_names = [Path(item["path"]).name for item in blank_rows]
            blank_policy = str(pairing.get("blank_filter_policy") or "")
            if blank_policy != "soft_evidence":
                raise RuntimeError("冻结空白配对不是已登记的soft_evidence比较操作")
            method_comparable = pairing.get("method_and_batch_comparability_user_confirmed") is True
            hard_exclusion_authorized = pairing.get("hard_background_exclusion_confirmed") is True
            task = {
                "id": task_id, "filename": Path(row["source_path"]).name,
                "file_path": row["source_path"], "source_sha256": row["source_sha256"],
                "status": "completed", "progress": 100, "results": [],
                "upload_time": row["started_at"], "start_time": row["started_at"],
                "end_time": conclusion["finished_at"], "state_version": 1,
                "run_attempt": 0, "active_run_id": None,
                "result_run_id": "frozen:" + key, "error": None,
                "total_spectra": total, "matched_spectra": matched,
                "match_rate": round(matched / total * 100, 1) if total else 0,
                "ion_mode": source["polarity"],
                "ms1_feature_count": summary["ms1_feature_records"],
                "unique_compounds": None,
                "compound_count_semantics": "saved_XBT_candidate_records_per_file_not_unique_chemicals",
                "params": {"db_mode": "full", "ion_mode": source["polarity"],
                           "enable_expected_compounds": True,
                           "expected_compound_datasets": [service.index["xbt_dataset"]],
                           "expected_compound_dataset_selection_source": "user_explicit",
                           "blank_filter_policy": blank_policy},
                "expected_compounds_summary": {"requested": True, "status": "ok",
                    "detected_feature_count": summary["ms1_feature_records"],
                    "evidence_level_counts": summary["ms1_evidence_level_counts"],
                    "standard_confirmed_count": 0,
                    "target_dataset_selection": deepcopy(target_selection)},
                "blank_summary": {"requested": True, "applicable": True,
                    "filter_policy": blank_policy,
                    "method_comparability_verified": method_comparable,
                    "background_exclusion_authorized": hard_exclusion_authorized,
                    "blank_operation": "comparison_only",
                    "blank_assignment_mode": "frozen_explicit_pairing",
                    "assigned_blank_replicate_count": len(blank_rows),
                    "assigned_blank_filenames": blank_names,
                    "assigned_blank_polarities": sorted({
                        str(item.get("observed_mode") or item.get("filename_mode") or "").upper()
                        for item in blank_rows if item
                    } - {""}),
                    "skipped_blank_replicate_count": 0,
                    "skipped_blank_filenames": [],
                    "blank_roles": {"background_blank": len(blank_rows), "unspecified": 0},
                    "flagged_spectra": summary["blank_flagged_ms2"],
                    "clean_view_excluded_spectra": summary["blank_clean_view_excluded_ms2"],
                    "excluded_spectra": 0},
                "message": "V16原生批处理已完成；MS1、MS/MS及配对空白结果只读接入",
                "completed_batch": {"id": COMPLETED_BATCH_ID, "sample_key": key,
                    "origin": "external_v16_native_batch", "read_only": True,
                    "registered_after_execution": True,
                    "report_url": service._summary()["report_url"],
                    "download_all_url": service._summary()["download_all_url"],
                    "summary": summary},
            }
            self.tasks[task_id] = task
            self.entries[task_id] = {"sample_key": key, "source_sha256": row["source_sha256"],
                "native_protocol_sha256": _NATIVE_PROTOCOL,
                **{kind: {"path": entry["verified_remote_path"],
                           "size": entry["verified_remote_size"],
                           "sha256": entry["verified_remote_sha256"]}
                   for kind, entry in (("results", row), ("conclusions", conclusion))}}
            legacy_task = deepcopy(task)
            legacy_task["params"].pop("expected_compound_dataset_selection_source", None)
            legacy_task["expected_compounds_summary"].pop("target_dataset_selection", None)
            for field in (
                "blank_operation", "blank_assignment_mode",
                "assigned_blank_replicate_count", "assigned_blank_filenames",
                "assigned_blank_polarities", "skipped_blank_replicate_count",
                "skipped_blank_filenames", "blank_roles", "excluded_spectra",
            ):
                legacy_task["blank_summary"].pop(field, None)
            legacy_task["blank_summary"]["method_comparability_verified"] = True
            legacy_task["blank_summary"]["background_exclusion_authorized"] = True
            legacy_sha256 = hashlib.sha256(
                json.dumps({"task": legacy_task, "artifacts": self.entries[task_id]},
                           sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()
            ).hexdigest()
            self._legacy_tasks[task_id] = legacy_task
            self._legacy_registration_sha256[task_id] = legacy_sha256
            task["completed_batch"]["registration_sha256"] = hashlib.sha256(
                json.dumps({"task": task, "artifacts": self.entries[task_id]},
                           sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()
            ).hexdigest()

    def install(self, tasks):
        # Removing a task card does not mutate its frozen remote artifacts.
        # Keep that user action persistent across service restarts.
        from config.settings import RUNTIME_ROOT
        deleted = set()
        audit = RUNTIME_ROOT / "runtime" / "task_deletion_audit.jsonl"
        if audit.is_file():
            for line in audit.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("event") in {"single_task_recycled", "batch_task_recycled"}:
                    deleted.update(row.get("task_ids") or [])
                    deleted.update(row.get("deleted_task_ids") or [])
        added = 0
        for task_id, task in self.tasks.items():
            existing = tasks.get(task_id)
            if existing is None and task_id in deleted:
                continue
            if existing is not None:
                existing_sha256 = (existing.get("completed_batch") or {}).get("registration_sha256")
                if existing_sha256 == self._legacy_registration_sha256[task_id]:
                    legacy = self._legacy_tasks[task_id]
                    for field, value in legacy.items():
                        existing_value = existing.get(field)
                        if field == "completed_batch":
                            existing_value = dict(existing_value or {})
                            existing_value.pop("registration_sha256", None)
                        if field != "results" and existing_value != value:
                            raise RuntimeError("Frozen registered task metadata changed: " + field)
                    transient_results = existing.get("results", [])
                    existing.update(deepcopy(task))
                    existing["results"] = transient_results
                    continue
                if existing_sha256 != task["completed_batch"]["registration_sha256"]:
                    raise RuntimeError("Frozen task ID collides with an existing task")
                # A task card may gain transient display fields; immutable identity
                # and completion facts must still agree with the frozen registration.
                for field, value in task.items():
                    if field != "results" and existing.get(field) != value:
                        raise RuntimeError("Frozen registered task metadata changed: " + field)
                continue
            tasks[task_id] = deepcopy(task)
            added += 1
        return added

    def read(self, task_id, operation, **options):
        request = {**self.entries[task_id], "operation": operation, **options}
        cache_key = json.dumps([task_id, operation, options], sort_keys=True)
        with self._lock:
            lock = self._read_locks.setdefault(cache_key, threading.Lock())
        with lock:
            with self._lock:
                if cache_key in self._cache:
                    raw = self._cache.pop(cache_key)
                    self._cache[cache_key] = raw
                    return json.loads(raw)
            program = Path(__file__).with_name("frozen_batch_reader.py").read_text()
            try:
                with _READERS:
                    completed = subprocess.run(
                        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
                         self.service.remote_host, "python3 -c " + shlex.quote(program)],
                        input=json.dumps(request).encode(), capture_output=True, timeout=180,
                    )
            except (OSError, subprocess.TimeoutExpired) as exc:
                raise RuntimeError("读取远端已完成结果失败，请稍后重试；不会重新分析") from exc
            if completed.returncode:
                # Tracebacks are operational evidence, not browser instructions.
                error = completed.stderr.decode(errors="replace")[-700:]
                raise RuntimeError("冻结结果读取或SHA绑定校验失败：" + error)
            raw = completed.stdout
            if len(raw) > 64 * 1024**2:
                raise RuntimeError("冻结结果投影超过64 MiB读取限制")
            value = json.loads(raw)
            if len(raw) <= _CACHE_BYTES:
                with self._lock:
                    while self._cache and self._cache_size + len(raw) > _CACHE_BYTES:
                        _, old = self._cache.popitem(last=False)
                        self._cache_size -= len(old)
                    self._cache[cache_key] = raw
                    self._cache_size += len(raw)
            return value

    def read_export(self, task_id, operation):
        """Return a verified saved export payload without creating local files."""
        request = {**self.entries[task_id], "operation": operation}
        program = Path(__file__).with_name("frozen_batch_reader.py").read_text()
        try:
            with _READERS:
                completed = subprocess.run(
                    ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
                     self.service.remote_host, "python3 -c " + shlex.quote(program)],
                    input=json.dumps(request).encode(), capture_output=True, timeout=180,
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise RuntimeError("读取远端已完成结果失败，请稍后重试；不会重新分析") from exc
        if completed.returncode:
            error = completed.stderr.decode(errors="replace")[-700:]
            raise RuntimeError("冻结结果读取或SHA绑定校验失败：" + error)
        return completed.stdout

    @staticmethod
    def export_response(payload, *, filename, media_type):
        def chunks():
            view = memoryview(payload)
            for offset in range(0, len(view), 256 * 1024):
                yield view[offset:offset + 256 * 1024]
        return StreamingResponse(chunks(), media_type=media_type, headers={
            "Content-Disposition": f'attachment; filename="{filename}"'})

    def csv_response(self, task_id, key, *, allowed_spectra=None, isomer_only=False):
        sample_key = self.entries[task_id]["sample_key"]
        path = self.service.ensure_exports([key])[key]
        def stream():
            yield b"\xef\xbb\xbf"
            with lzma.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                buffer = io.StringIO()
                writer = csv.DictWriter(buffer, fieldnames=reader.fieldnames)
                writer.writeheader()
                yield buffer.getvalue().encode()
                buffer.seek(0); buffer.truncate(0)
                for row in reader:
                    if row.get("sample_key") == sample_key:
                        if allowed_spectra is not None and row.get("spectrum_id") not in allowed_spectra:
                            continue
                        if isomer_only and str(row.get("isomer_limited", "")).lower() not in {"true", "1"}:
                            continue
                        writer.writerow(row)
                        yield buffer.getvalue().encode()
                        buffer.seek(0); buffer.truncate(0)
        return StreamingResponse(stream(), media_type="text/csv; charset=utf-8", headers={
            "Content-Disposition": f'attachment; filename="{task_id}_{key}.csv"'})


def install_frozen_task_access(app, bridge, tasks, row_fields, match_fields):
    """Route imported tasks to saved artifacts before native migration handlers."""
    if bridge is None:
        return

    @app.middleware("http")
    async def frozen_task_access(request, call_next):
        parts = request.url.path.strip("/").split("/")
        task_id = parts[2] if len(parts) >= 3 and parts[:2] == ["api", "tasks"] else None
        if task_id not in bridge.tasks:
            if request.method == "POST" and request.url.path.startswith("/api/export/jobs/"):
                selected = await request.json()
                if isinstance(selected, list) and any(key in bridge.tasks for key in selected):
                    return JSONResponse({"detail": "本批完整报告与全部CSV/JSONL已生成，请使用宣痹汤完整交付包",
                        "download_url": bridge.service._summary()["download_all_url"]}, status_code=409)
            if request.method == "POST" and request.url.path == "/api/tasks/batch/analyze":
                payload = await request.json()
                if any(key in bridge.tasks for key in payload.get("task_ids", [])):
                    return JSONResponse({"detail": "宣痹汤已完成任务只读保留，请从所选任务中移除后操作其他任务"}, status_code=409)
            return await call_next(request)
        suffix = "/".join(parts[3:])
        registry = tasks() if callable(tasks) else tasks
        if task_id not in registry:
            return JSONResponse({"detail": "Task not found"}, status_code=404)
        if request.method == "DELETE" and not suffix:
            return await call_next(request)
        if request.method != "GET":
            return JSONResponse({"detail": "宣痹汤任务为已完成结果的只读接入；分析、重试及结果改写均已禁用"}, status_code=409)
        query = request.query_params
        try:
            def integer(name, default, maximum=5000):
                value = int(query.get(name, default))
                if value < (1 if name == "limit" else 0) or value > maximum:
                    raise ValueError("分页参数超出允许范围")
                return value
            def choice(name, default, allowed):
                value = query.get(name, default)
                if value not in allowed:
                    raise ValueError("未知筛选参数：" + name)
                return value
            offset = integer("offset", 0, 1_000_000)
            limit = integer("limit", 500, 5000)
            async def read(operation, **options):
                return await asyncio.to_thread(bridge.read, task_id, operation, **options)
            if suffix in {"", "status"}:
                registry = tasks() if callable(tasks) else tasks
                task = {key: value for key, value in registry[task_id].items() if key != "results"}
                return JSONResponse(task)
            if suffix == "results":
                if query.get("projection", "full") not in {"summary", "full"}:
                    raise ValueError("未知结果投影")
                data = await read("summary", row_fields=sorted(row_fields), match_fields=sorted(match_fields))
                rows = data["results"]
                blank = choice("blank_view", "all", {"all", "clean", "flagged"})
                if blank == "clean": rows = [row for row in rows if not row.get("blank_clean_view_excluded")]
                elif blank == "flagged": rows = [row for row in rows if row.get("blank_flagged")]
                state = choice("match_state", "all", {"all", "matched", "unmatched"})
                if query.get("matched_only") == "true": state = "matched"
                counts = {"total": len(rows), "matched": sum(bool(row.get("match")) for row in rows),
                          "unmatched": sum(not row.get("match") for row in rows),
                          "blank_flagged": sum(bool(row.get("blank_flagged")) for row in rows)}
                if state == "matched": rows = [row for row in rows if row.get("match")]
                elif state == "unmatched": rows = [row for row in rows if not row.get("match")]
                if query.get("min_cosine"):
                    rows = [row for row in rows if row.get("match") and float(row["match"].get("cosine", 0)) >= float(query["min_cosine"])]
                governance = choice("governance_view", "all", {"all", "preeligible", "blocked", "isomer_group"})
                if governance != "all":
                    def include(row):
                        match = row.get("match") or {}
                        evidence = match.get("evidence_governance") or {}
                        if governance == "preeligible": return bool(evidence.get("release_eligible_before_replication"))
                        if governance == "blocked": return bool(match) and not evidence.get("release_eligible_before_replication")
                        return evidence.get("claim_status") == "isomer_group_only"
                    rows = [row for row in rows if include(row)]
                page = rows[offset:offset + limit]
                # Full saved rows are fetched individually, never regenerated.
                if query.get("projection", "full") == "full":
                    if len(page) > 500:
                        return JSONResponse({"detail": "完整谱图分页最多500条；列表可使用 projection=summary"}, status_code=422)
                    page = (await read("items", indices=[row["result_index"] for row in page]))["results"]
                return JSONResponse({"results": page, "total": len(rows), "offset": offset, "limit": limit,
                    "has_more": offset + limit < len(rows), "view_counts": counts,
                    "projection": query.get("projection", "full"), "saved_results_unchanged": True})
            if suffix.startswith("results/item/"):
                data = await read("item", index=int(parts[-1]))
                return JSONResponse({**data, "task_id": task_id, "projection": "full_single_item",
                                     "saved_legacy_result_unchanged": True})
            if suffix.startswith("compounds/item/"):
                result_view = choice("result_view", "native", {"native", "xbt_linked"})
                locator = suffix[len("compounds/item/"):]
                options = ({"index": int(locator)} if locator.isdigit()
                           else {"compound_id": locator})
                data = await read("compound_item", result_view=result_view, **options)
                return JSONResponse({**data, "task_id": task_id,
                    "saved_compound_conclusion_unchanged": True,
                    "no_analysis_rerun": True})
            if suffix == "compounds":
                scope = choice("report_scope", "all", {"all", "plant_candidates", "external_review"})
                result_view = choice("result_view", "native", {"native", "xbt_linked"})
                if limit > 500:
                    raise ValueError("冻结化合物结论分页最多500条")
                data = await read("compound_page", report_scope=scope, result_view=result_view,
                                  conclusion=query.get("conclusion"), offset=offset, limit=limit)
                return JSONResponse({**data, "compounds_total": data["filtered_total"], "task_id": task_id,
                    "governance_blocked_candidates": [], "governance_blocked_has_more": False,
                    "result_unit": ("saved_native_compound_conclusion_records"
                                    if result_view == "native" else "saved_XBT_linked_subset_records"),
                    "no_analysis_rerun": True})
            if suffix == "flagship":
                data = await read("flagship", **({"entity_id": query["entity_id"]} if query.get("entity_id") else {}))
                rows = data["entities"]
                if query.get("entity_id"): rows = [row for row in rows if row.get("entity_id") == query["entity_id"]]
                if query.get("decision_status"): rows = [row for row in rows if row.get("decision", {}).get("status") == query["decision_status"]]
                return JSONResponse({**data, "entities": rows[offset:offset + limit], "entities_total": len(rows),
                    "offset": offset, "limit": limit, "has_more": offset + limit < len(rows),
                    "task_id": task_id, "legacy_scan_results_unchanged": True})
            if suffix == "flagship/structure":
                if not query.get("entity_id"):
                    raise ValueError("缺少entity_id")
                data = await read("flagship", entity_id=query["entity_id"])
                if not data["entities"]:
                    return JSONResponse({"detail": "冻结实体不存在"}, status_code=404)
                entity = data["entities"][0]
                selected = entity.get("winner") or entity.get("leading_candidate_when_abstained") or {}
                from rdkit import Chem
                from rdkit.Chem.Draw import rdMolDraw2D
                molecule = Chem.MolFromSmiles(str((selected.get("candidate") or {}).get("smiles") or ""))
                if molecule is None:
                    return JSONResponse({"detail": "冻结候选没有可绘制结构"}, status_code=404)
                drawing = rdMolDraw2D.MolDraw2DSVG(420, 300)
                drawing.DrawMolecule(molecule); drawing.FinishDrawing()
                return Response(drawing.GetDrawingText(), media_type="image/svg+xml")
            if suffix in {"expected-compounds", "expected-compounds/groups"}:
                data = await read("expected")
                features = data["features"]
                if suffix.endswith("/groups"):
                    from web.expected_compounds_export import expected_group_rows
                    groups = expected_group_rows(data, language="zh")
                    supported = [row for row in groups if int(row.get("supported_feature_count") or 0) > 0]
                    view = choice("evidence_view", "all", {"all", "supported", "mass_only"})
                    selected = supported if view == "supported" else ([row for row in groups if row not in supported] if view == "mass_only" else groups)
                    return JSONResponse({"groups": selected[offset:offset + limit], "group_total": len(groups),
                        "filtered_group_total": len(selected), "supported_group_total": len(supported),
                        "mass_only_group_total": len(groups) - len(supported), "summary": data["summary"],
                        "status": "ok", "offset": offset, "limit": limit,
                        "has_more": offset + limit < len(selected), "ordinary_library_match_counted": False})
                if query.get("supported_only") == "true": features = [row for row in features if row.get("cd_pattern_accepted")]
                if query.get("evidence_level"): features = [row for row in features if row.get("evidence_level") == query["evidence_level"]]
                return JSONResponse({**data, "features": features[offset:offset + limit], "filtered_total": len(features),
                    "offset": offset, "limit": limit, "has_more": offset + limit < len(features),
                    "ordinary_library_match_counted": False})
            if suffix == "download/csv":
                layer = choice("result_layer", "ms2", {"ms2", "targeted", "targeted_features",
                              "plant_candidates", "isomer_candidates"})
                blank_view = choice("blank_view", "all", {"all", "clean", "flagged"})
                spectra = None
                if blank_view != "all":
                    if layer != "ms2":
                        raise ValueError("该冻结表按原导出保留全部记录；空白筛选仅用于MS/MS事件表")
                    saved = await read("summary", row_fields=["spectrum_id", "blank_clean_view_excluded", "blank_flagged"], match_fields=[])
                    spectra = {row["spectrum_id"] for row in saved["results"]
                               if (not row.get("blank_clean_view_excluded") if blank_view == "clean" else row.get("blank_flagged"))}
                key = {"targeted": "ms1_table", "targeted_features": "ms1_table",
                       "plant_candidates": "xbt_csv", "isomer_candidates": "xbt_csv",
                       "external_review": "qc_table"}.get(layer, "msms_table")
                return await asyncio.to_thread(bridge.csv_response, task_id, key,
                                               allowed_spectra=spectra, isomer_only=layer == "isomer_candidates")
            if suffix.startswith("frozen-tables/"):
                key = suffix.split("/")[1]
                if key == "native_compounds":
                    payload = await asyncio.to_thread(bridge.read_export, task_id, "native_compounds_csv")
                    return bridge.export_response(payload, filename=f"{task_id}_native_compounds.csv",
                                                  media_type="text/csv; charset=utf-8")
                if key not in {"xbt_csv", "msms_table", "ms1_table", "qc_table"}:
                    raise ValueError("未知冻结结果表")
                return await asyncio.to_thread(bridge.csv_response, task_id, key)
            if suffix == "download/native-conclusions.json.xz":
                payload = await asyncio.to_thread(bridge.read_export, task_id, "native_conclusions_xz")
                return bridge.export_response(payload, filename=f"{task_id}_native_conclusions.json.xz",
                                              media_type="application/x-xz")
            if suffix in {"report", "download/html"}:
                from web.frozen_batch_report import render_frozen_sample_report
                html = await asyncio.to_thread(render_frozen_sample_report, bridge, task_id)
                headers = ({"Content-Disposition": f'attachment; filename="{task_id}_full_report.html"'}
                           if suffix == "download/html" else None)
                return HTMLResponse(html, headers=headers)
            if suffix == "replicate-consistency":
                return JSONResponse({"status": "not_computed_in_frozen_batch", "compounds": []})
            if suffix == "acquisition-plans":
                return JSONResponse({"plans": [], "metadata": {"status": "not_run_in_frozen_batch"}})
            if suffix == "identification-records":
                return JSONResponse({"records": [], "metadata": {"status": "not_saved_in_frozen_batch"}})
            return JSONResponse({"detail": "本批未保存该独立模块产物。已完成MS1、MS/MS和配对空白结果可在当前结果页与完整批次报告读取",
                                 "completed_batch": bridge.tasks[task_id]["completed_batch"]}, status_code=404)
        except (ValueError, IndexError) as exc:
            return JSONResponse({"detail": str(exc)}, status_code=422)
        except (RuntimeError, OSError, KeyError) as exc:
            return JSONResponse({"detail": str(exc)}, status_code=502)
