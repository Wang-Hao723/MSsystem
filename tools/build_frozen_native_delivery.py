#!/usr/bin/env python3
"""Build the frozen 149-sample native-conclusion delivery on msqc-server."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import html
import io
import json
import lzma
import os
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


SCHEMA = "mssystem.frozen-native-delivery/1.0"
REMOTE_OUTPUT = Path("/data/Mssystem/xbt_real_batch_20261021_ms2/delivery_access_audit_20261022")
BUDGET = 320 * 1024**2
FLOOR = 25 * 1024**3


def _sha_bytes(value):
    return hashlib.sha256(value).hexdigest()


def _sha_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _semantic(value):
    return _sha_bytes(json.dumps(value, ensure_ascii=False, sort_keys=True,
                                  separators=(",", ":")).encode())


def _atomic_bytes(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    with temporary.open("wb") as handle:
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _html_report(samples, totals):
    rows = []
    for row in samples:
        task_url = f"/api/tasks/{row['task_id']}/report"
        rows.append("<tr>" + "".join(
            f"<td>{html.escape(str(value))}</td>" for value in (
                row["sample_id"], row["polarity"], row["native_conclusion_rows"],
                row["matched_ms2_spectra"], row["ms1_features"],
                row["xbt_linked_rows"], row["external_review_rows"],
            )
        ) + f'<td><a href="{task_url}">逐样本完整报告</a></td></tr>')
    body = "".join(rows)
    return f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>宣痹汤149份完整原生结论交付</title>
<style>body{{font:15px/1.6 system-ui,sans-serif;margin:30px;color:#24323c}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #d8e0e3;padding:7px;text-align:left}}th{{background:#eef4f4;position:sticky;top:0}}.scroll{{overflow:auto;max-height:75vh}}code{{overflow-wrap:anywhere}}</style></head><body>
<h1>宣痹汤149份完整原生结论交付</h1>
<p>本页汇总冻结工作流保存的谱库候选、未解析事件、背景与未注释记录，共 {totals['native_conclusion_rows']} 条记录出现次数；标准品确证 0。来源未知、外源或身份冲突记录保留在QC证据中，不列作已核验植物内源成分。</p>
<p>已检索MS/MS记录 {totals['searched_ms2_spectra']}，产生谱库候选的MS/MS记录 {totals['matched_ms2_spectra']}，独立MS1特征 {totals['ms1_features']}，XBT关联附加筛选记录 {totals['xbt_linked_rows']}。XBT数量是附加筛选记录数，不是命名化合物数。</p>
<div class="scroll"><table><thead><tr><th>样本</th><th>极性</th><th>原生结论记录</th><th>匹配MS/MS</th><th>MS1特征</th><th>XBT附加筛选</th><th>外源/冲突QC</th><th>链接</th></tr></thead><tbody>{body}</tbody></table></div>
<p>计数单位为跨样本记录出现次数。该交付仅格式化已验证冻结文件，没有重新检索、评分或分析。</p></body></html>'''.encode("utf-8")


def _build_delivery(payload, output_dir, reader_source, *, workers=2):
    output_dir = Path(output_dir)
    if workers < 1 or workers > 2:
        raise ValueError("workers must be 1 or 2")
    if len(payload.get("samples") or []) != 149:
        raise ValueError("delivery requires exactly 149 registered samples")
    if len({row["sample_key"] for row in payload["samples"]}) != 149:
        raise ValueError("delivery sample keys are not unique")
    source_total = sum(row["results"]["size"] + row["conclusions"]["size"]
                       for row in payload["samples"])
    if source_total >= 4 * 1024**3:
        raise ValueError("frozen source batch exceeds 4 GiB delivery boundary")
    output_dir.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(output_dir).free - BUDGET < FLOOR:
        raise OSError("remote delivery volume would fall below 25 GiB floor")
    started = datetime.now(timezone.utc).isoformat()
    log_path = output_dir / "BUILD_LOG.jsonl"
    with log_path.open("a", encoding="utf-8") as log:
        log.write(json.dumps({"event": "started", "at": started, "pid": os.getpid(),
                              "workers": workers, "source_compressed_bytes": source_total}) + "\n")

    namespace = {"__name__": "frozen_delivery_reader"}
    exec(reader_source, namespace)
    def one(row):
        request = {key: row[key] for key in (
            "sample_key", "source_sha256", "native_protocol_sha256", "results", "conclusions"
        )}
        document = namespace["read_document"](request, "conclusions")
        native = document["compound_conclusions"]
        compounds = native.get("compounds")
        if not isinstance(compounds, list):
            raise ValueError(row["sample_key"] + ": native compounds missing")
        expected_rows = row["summary"].get("compound_result_records")
        if expected_rows is not None and int(expected_rows) != len(compounds):
            raise ValueError(row["sample_key"] + ": native row count differs from frozen summary")
        # Reuse the already verified sidecar document.  Calling ``project``
        # here would decompress the same potentially large file a second time.
        csv_bytes = namespace["_compound_csv"](document, request)
        conclusion_counts = Counter((item.get("conclusion") or {}).get("key") or "unknown"
                                    for item in compounds)
        record_counts = Counter(item.get("record_type") or "unknown" for item in compounds)
        xbt_rows = len(((document.get("report_views") or {}).get("xbt_named_candidates") or {}).get("compounds") or [])
        external = sum(bool((item.get("plant_report") or {}).get(
            "excluded_from_plant_candidate_table")) for item in compounds)
        return row, csv_bytes, dict(record_counts), dict(conclusion_counts), xbt_rows, external

    ordered = sorted(payload["samples"], key=lambda row: row["sample_key"])
    csv_path = output_dir / "all_native_conclusions.csv.xz"
    temporary_csv = csv_path.with_name(csv_path.name + f".partial.{os.getpid()}")
    samples = []
    all_record_types = Counter()
    all_conclusions = Counter()
    totals = Counter()
    expected_header = None
    try:
        with lzma.open(temporary_csv, "wb", preset=6) as output:
          with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            iterator = iter(ordered)
            pending = []
            for _ in range(workers):
                try:
                    pending.append(pool.submit(one, next(iterator)))
                except StopIteration:
                    break
            while pending:
                future = pending.pop(0)
                row, csv_bytes, record_counts, conclusion_counts, xbt_rows, external = future.result()
                text = csv_bytes.decode("utf-8-sig")
                header, body = text.split("\n", 1)
                if expected_header is None:
                    expected_header = header
                    output.write(b"\xef\xbb\xbf" + (header + "\n").encode())
                elif header != expected_header:
                    raise ValueError("native CSV headers differ across samples")
                output.write(body.encode())
                native_rows = sum(record_counts.values())
                all_record_types.update(record_counts)
                all_conclusions.update(conclusion_counts)
                totals.update({
                    "native_conclusion_rows": native_rows,
                    "searched_ms2_spectra": int(row["summary"]["parsed_searched_ms2_spectra"]),
                    "matched_ms2_spectra": int(row["summary"]["matched_ms2_spectra"]),
                    "ms1_features": int(row["summary"]["ms1_feature_records"]),
                    "xbt_linked_rows": xbt_rows,
                    "external_review_rows": external,
                })
                samples.append({
                    "sample_key": row["sample_key"], "task_id": row["task_id"],
                    "sample_id": row["summary"]["source_sample"].get("sample_id"),
                    "polarity": row["summary"]["source_sample"].get("polarity"),
                    "source_sha256": row["source_sha256"],
                    "native_result": row["results"], "conclusion_result": row["conclusions"],
                    "native_conclusion_rows": native_rows,
                    "record_type_counts": record_counts, "conclusion_status_counts": conclusion_counts,
                    "searched_ms2_spectra": int(row["summary"]["parsed_searched_ms2_spectra"]),
                    "matched_ms2_spectra": int(row["summary"]["matched_ms2_spectra"]),
                    "ms1_features": int(row["summary"]["ms1_feature_records"]),
                    "xbt_linked_rows": xbt_rows, "external_review_rows": external,
                    "standard_confirmed_count": 0,
                })
                try:
                    pending.append(pool.submit(one, next(iterator)))
                except StopIteration:
                    pass
        os.replace(temporary_csv, csv_path)
    finally:
        temporary_csv.unlink(missing_ok=True)
    html_bytes = _html_report(samples, totals)
    html_path = output_dir / "FULL_NATIVE_RESULTS_149.html"
    _atomic_bytes(html_path, html_bytes)
    artifacts = {
        "native_compounds": {"path": str(csv_path), "size": csv_path.stat().st_size,
                             "sha256": _sha_file(csv_path), "rows": totals["native_conclusion_rows"]},
        "report": {"path": str(html_path), "size": html_path.stat().st_size,
                   "sha256": _sha_file(html_path)},
    }
    index = {
        "schema_version": SCHEMA, "status": "completed",
        "created_at": datetime.now(timezone.utc).isoformat(), "builder_pid": os.getpid(),
        "workers": workers, "sample_count": 149,
        "native_conclusion_rows": totals["native_conclusion_rows"],
        "record_type_counts": dict(all_record_types),
        "conclusion_status_counts": dict(all_conclusions),
        "counting_unit": "record_occurrences", "standard_confirmed_count": 0,
        "source_compressed_bytes": source_total,
        "original_delivery_index_path": payload["original_delivery_index_path"],
        "original_delivery_index_sha256": payload["original_delivery_index_sha256"],
        "task_ids": [row["task_id"] for row in samples], "sample_keys": [row["sample_key"] for row in samples],
        "totals": dict(totals), "samples": samples, "artifacts": artifacts,
        "analysis_rerun": False, "source_files_written": False,
    }
    index["semantic_sha256"] = _semantic(index)
    index_path = output_dir / "FULL_NATIVE_DELIVERY_INDEX.json"
    _atomic_bytes(index_path, json.dumps(index, ensure_ascii=False, indent=2, sort_keys=True).encode())
    receipt = {
        "status": "completed", "finished_at": datetime.now(timezone.utc).isoformat(),
        "pid": os.getpid(), "sample_count": 149,
        "native_conclusion_rows": totals["native_conclusion_rows"],
        "index_path": str(index_path), "index_sha256": _sha_file(index_path),
        "artifacts": artifacts, "analysis_rerun": False,
    }
    _atomic_bytes(output_dir / "BUILD_RECEIPT.json",
                  json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True).encode())
    with log_path.open("a", encoding="utf-8") as log:
        log.write(json.dumps({"event": "completed", "at": receipt["finished_at"],
                              "pid": os.getpid(), "rows": totals["native_conclusion_rows"]}) + "\n")
    return index


def _payload(project_root):
    import importlib.util
    import types
    root = Path(project_root).resolve()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    # Standalone delivery must not execute web/__init__.py, which imports the
    # live server and its heavyweight runtime.  Load only the two read-only
    # modules required to construct the frozen registry.
    web_package = types.ModuleType("web")
    web_package.__path__ = [str(root / "web")]
    sys.modules["web"] = web_package
    def load(name):
        spec = importlib.util.spec_from_file_location("web." + name, root / "web" / (name + ".py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules["web." + name] = module
        spec.loader.exec_module(module)
        return module
    completed = load("completed_batches")
    frozen = load("frozen_batch_tasks")
    CompletedBatchService, _sha256 = completed.CompletedBatchService, completed._sha256
    FrozenBatchTasks = frozen.FrozenBatchTasks
    batch_root = root.parent / "research/20261021_xbt_ms2_integrated"
    cache_root = Path("/Volumes/LaCie/Mssystem_V16/completed_batch_exports")
    service = CompletedBatchService(batch_root, cache_root, "msqc-server")
    bridge = FrozenBatchTasks(service)
    samples = []
    for task_id, entry in sorted(bridge.entries.items(), key=lambda item: item[1]["sample_key"]):
        samples.append({**entry, "task_id": task_id,
                        "summary": bridge.tasks[task_id]["completed_batch"]["summary"]})
    return {
        "original_delivery_index_path": str(service.paths["index"]),
        "original_delivery_index_sha256": _sha256(service.paths["index"]),
        "samples": samples,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--host", default="msqc-server")
    parser.add_argument("--workers", type=int, default=2, choices=(1, 2))
    parser.add_argument("--output", default=str(REMOTE_OUTPUT))
    args = parser.parse_args(argv)
    if Path(args.output) != REMOTE_OUTPUT:
        raise SystemExit("remote output path is fixed")
    payload = _payload(args.project_root)
    reader_source = (Path(args.project_root) / "web/frozen_batch_reader.py").read_text()
    source = Path(__file__).read_text()
    wrapper = (
        "import json\nns={'__name__':'remote_delivery_builder'}\n"
        + "exec(" + repr(source) + ",ns)\n"
        + "result=ns['_build_delivery'](" + repr(payload) + "," + repr(args.output)
        + "," + repr(reader_source) + ",workers=" + str(args.workers) + ")\n"
        + "print(json.dumps({'status':'completed','sample_count':result['sample_count'],"
          "'rows':result['native_conclusion_rows'],'semantic_sha256':result['semantic_sha256']}))\n"
    )
    completed = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12", args.host, "python3", "-"],
        input=wrapper.encode(), timeout=3600,
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
