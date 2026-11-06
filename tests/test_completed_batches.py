import csv
import hashlib
import io
import json
import lzma
import subprocess
import sys
import threading
import time
import zipfile
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from web import completed_batches as module


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _csv_xz(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return lzma.compress(output.getvalue().encode())


@pytest.fixture()
def completed_fixture(tmp_path, monkeypatch):
    batch = tmp_path / "batch"
    cache = tmp_path / "external-cache"
    batch.mkdir()
    cache.mkdir()
    sample_keys = [f"sample-{number:03d}" for number in range(149)]
    common = [{"sample_key": sample_keys[0], "value": "first"}, {"sample_key": sample_keys[0], "value": "second"}, {"sample_key": sample_keys[1], "value": "third"}]
    table_rows = {
        "file_manifest": [{"sample_key": key, "source": f"{key}.mzML"} for key in sample_keys],
        "xbt_csv": [{**row, "candidate_name": f"candidate-{index}"} for index, row in enumerate(common)],
        "msms_table": [{**row, "matched": str(index == 0).lower()} for index, row in enumerate(common)],
        "ms1_table": [{**row, "ms2_association": "" if index == 1 else f"ms2-{index}"} for index, row in enumerate(common)],
        "qc_table": [{**row, "qc_view": "unlinked" if index < 2 else "external"} for index, row in enumerate(common)],
    }
    payloads = {
        key: _csv_xz(rows) for key, rows in table_rows.items()
    }
    payloads.update({
        "xbt_jsonl": lzma.compress(b"{\"sample_key\":\"sample-000\"}\n"),
        "qc": lzma.compress(b"{\"sample_key\":\"sample-000\"}\n"),
        "ms1": lzma.compress(b"{\"sample_key\":\"sample-000\"}\n"),
        "summary": b'{"status":"completed"}\n',
    })
    exports = []
    table_counts = {}
    for name, filename in module._EXPORT_FILES.items():
        data = payloads[name]
        exports.append({
            "name": name,
            "path": str(module._REMOTE_EXPORT_DIR / filename),
            "sha256": _sha(data),
            "size": len(data),
        })
        if name in table_rows:
            table_counts[filename] = len(table_rows[name])

    totals = {
        "files": 149,
        "declared_ms2_scans": 3,
        "searched_ms2_spectra": 3,
        "matched_ms2_spectra": 1,
        "ms1_features": 3,
        "xbt_named_candidate_rows": 3,
        "ambiguity_group_rows": 1,
        "standard_confirmed": 0,
        "unlinked_qc_rows": 2,
        "external_qc_rows": 1,
        "background_qc_rows": 0,
        "name_qc_rows": 0,
    }
    index = {
        "status": "complete",
        "totals": totals,
        "paired_groups": 50,
        "counting_unit": "跨文件记录出现次数，不是去重化合物数",
        "ms1_evidence_counts": {"MS1_EXACT_MASS_CANDIDATE": 2, "MS1_ISOTOPE_SUPPORTED": 1},
        "declared_vs_searched_ms2": {"declared": 3, "searched": 3, "native_parser_omitted_records": 0},
        "msclassify": {"status": "未运行"},
        "xbt_dataset": "XBT_targeted_20260817",
        "models_written_by_workflow": False,
        "remote_exports": exports,
        "table_row_counts": table_counts,
    }
    index_bytes = json.dumps(index, ensure_ascii=False).encode()
    report_bytes = b"<!doctype html><title>frozen report</title>"
    (batch / "FINAL_DELIVERY_INDEX.json").write_bytes(index_bytes)
    (batch / "FINAL_RESULTS_OVERVIEW.html").write_bytes(report_bytes)
    packaging = {"status": "completed", "index_sha256": _sha(index_bytes), "html_sha256": _sha(report_bytes)}
    (batch / "FINAL_PACKAGING_RECEIPT.json").write_text(json.dumps(packaging), encoding="utf-8")
    (batch / "FINAL_DELIVERY_VERIFICATION.json").write_text(json.dumps({"status": "passed"}), encoding="utf-8")

    native_samples = {}
    conclusion_samples = {}
    for index_number, key in enumerate(sample_keys):
        is_first = index_number == 0
        is_second = index_number == 1
        native_samples[key] = {"spectra_completed": 2 if is_first else (1 if is_second else 0), "matched_spectra": int(is_first), "ms1_features": 2 if is_first else (1 if is_second else 0)}
        conclusion_samples[key] = {"summary": {
            "source_sample": {"sample_key": key, "sample_id": f"S{index_number:03d}", "polarity": "POS" if index_number % 2 == 0 else "NEG", "replicate": 1},
            "parsed_searched_ms2_spectra": native_samples[key]["spectra_completed"],
            "matched_ms2_spectra": native_samples[key]["matched_spectra"],
            "ms1_feature_records": native_samples[key]["ms1_features"],
            "xbt_named_candidate_rows": 2 if is_first else (1 if is_second else 0),
            "blank_flagged_ms2": 1 if is_first else 0,
            "blank_clean_view_excluded_ms2": 1 if is_first else 0,
            "ambiguity_group_rows": 1 if is_first else 0,
            "xbt_link_review_qc_rows": 2 if is_first else 0,
            "external_review_rows": 1 if is_second else 0,
            "background_qc_rows": 0,
            "database_title_or_unnamed_qc_rows": 0,
            "ms1_evidence_level_counts": {"MS1_EXACT_MASS_CANDIDATE": native_samples[key]["ms1_features"]},
        }}
    (batch / "STATE.json").write_text(json.dumps({"state": "completed", "samples": native_samples}), encoding="utf-8")
    (batch / "CONCLUSIONS_STATE.json").write_text(json.dumps({"state": "completed", "samples": conclusion_samples}), encoding="utf-8")

    monkeypatch.setattr(module, "_ALLOWED_CACHE_PREFIXES", (tmp_path.resolve(),))
    app = FastAPI()
    service = module.register_completed_batch_routes(app, batch_root=batch, cache_root=cache, remote_host="test-host")
    return {
        "batch": batch,
        "cache": cache,
        "payloads": payloads,
        "service": service,
        "client": TestClient(app),
        "sample_keys": sample_keys,
        "report_bytes": report_bytes,
        "index_bytes": index_bytes,
    }


def _populate_cache(fixture, names=None):
    names = names or module._EXPORT_FILES.keys()
    for name in names:
        filename = module._EXPORT_FILES[name]
        (fixture["cache"] / filename).write_bytes(fixture["payloads"][name])


def test_list_and_detail_are_read_only_and_expose_149_samples(completed_fixture, monkeypatch):
    calls = []
    modules_before = set(sys.modules)
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: calls.append(args))
    client = completed_fixture["client"]
    listed = client.get("/api/completed-batches")
    assert listed.status_code == 200
    batch = listed.json()["batches"][0]
    assert batch["id"] == module.COMPLETED_BATCH_ID
    assert batch["completed_samples"] == batch["total_samples"] == 149
    assert batch["paired_blank_files"] == 102

    detail = client.get(f"/api/completed-batches/{module.COMPLETED_BATCH_ID}")
    assert detail.status_code == 200
    payload = detail.json()
    assert len(payload["samples"]) == 149
    assert payload["samples"][0]["searched_ms2_spectra"] == 2
    assert payload["models_written_by_workflow"] is False
    assert payload["no_analysis_rerun"] is True
    assert len(payload["exports"]) == 9
    assert calls == []
    newly_imported = set(sys.modules) - modules_before
    assert not any(name.startswith(("matching", "retriever")) for name in newly_imported)


def test_table_pagination_sample_filter_and_ms1_only_are_lossless(completed_fixture):
    _populate_cache(completed_fixture, ["ms1_table", "xbt_csv"])
    client = completed_fixture["client"]
    batch = module.COMPLETED_BATCH_ID
    sample = completed_fixture["sample_keys"][0]
    page = client.get(f"/api/completed-batches/{batch}/tables/ms1_table", params={"sample_key": sample, "offset": 1, "limit": 1})
    assert page.status_code == 200
    body = page.json()
    assert body["total"] == 2
    assert body["rows"] == [{"sample_key": sample, "value": "second", "ms2_association": ""}]
    assert body["has_more"] is False

    export = client.get(f"/api/completed-batches/{batch}/exports/xbt_csv", params={"decompress": "true"})
    assert export.status_code == 200
    assert export.content.startswith(b"\xef\xbb\xbf")
    assert export.content[3:] == lzma.decompress(completed_fixture["payloads"]["xbt_csv"])


def test_unknown_batch_export_and_path_are_rejected(completed_fixture):
    client = completed_fixture["client"]
    batch = module.COMPLETED_BATCH_ID
    assert client.get("/api/completed-batches/not-this-batch").status_code == 404
    assert client.get(f"/api/completed-batches/{batch}/exports/../../STATE.json").status_code in {404, 422}
    assert client.get(f"/api/completed-batches/{batch}/exports/not-an-export").status_code == 404
    assert client.get(f"/api/completed-batches/{batch}/tables/xbt_jsonl").status_code == 404


def test_wrong_hash_fails_and_never_publishes_partial(completed_fixture, monkeypatch):
    def fake_run(args, **kwargs):
        Path(args[-1]).write_bytes(b"corrupt")
        return subprocess.CompletedProcess(args, 0, b"", b"")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(module, "_CACHE_FLOOR_BYTES", 0)
    response = completed_fixture["client"].get(
        f"/api/completed-batches/{module.COMPLETED_BATCH_ID}/exports/xbt_csv"
    )
    assert response.status_code == 502
    assert "SHA-256" in response.json()["detail"]
    assert not (completed_fixture["cache"] / module._EXPORT_FILES["xbt_csv"]).exists()
    assert not list(completed_fixture["cache"].glob("*.part.*"))


def test_unmounted_external_volume_never_falls_back_to_system_disk(completed_fixture, monkeypatch):
    forbidden = Path("/Volumes/Mssystem-test-volume-must-not-exist/completed-batches")
    assert not forbidden.parent.exists()
    monkeypatch.setattr(module, "_ALLOWED_CACHE_PREFIXES", (Path("/Volumes"),))
    completed_fixture["service"].cache_root = forbidden
    with pytest.raises(RuntimeError, match="未挂载"):
        completed_fixture["service"].ensure_exports(["xbt_csv"])
    assert not forbidden.exists()


def test_streamed_zip_contains_all_nine_exports_and_frozen_receipts(completed_fixture):
    _populate_cache(completed_fixture)
    response = completed_fixture["client"].get(
        f"/api/completed-batches/{module.COMPLETED_BATCH_ID}/download-all"
    )
    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        names = archive.namelist()
        assert len(names) == len(set(names)) == 13
        for key, filename in module._EXPORT_FILES.items():
            assert archive.read(filename) == completed_fixture["payloads"][key]
        assert archive.read("宣痹汤_XBT_完整结果.html") == completed_fixture["report_bytes"]
        assert archive.read("FINAL_DELIVERY_INDEX.json") == completed_fixture["index_bytes"]
        csv_rows = list(csv.DictReader(io.StringIO(lzma.decompress(archive.read(module._EXPORT_FILES["file_manifest"])).decode())))
        assert len(csv_rows) == 149


def test_report_index_alias_and_validation_limits(completed_fixture):
    client = completed_fixture["client"]
    batch = module.COMPLETED_BATCH_ID
    report = client.get(f"/api/completed-batches/{batch}/report")
    assert report.content == completed_fixture["report_bytes"]
    assert "inline" in report.headers["content-disposition"]
    assert client.get(f"/api/completed-batches/{batch}/index").content == completed_fixture["index_bytes"]
    assert client.get(f"/api/completed-batches/{batch}/FINAL_DELIVERY_INDEX.json").content == completed_fixture["index_bytes"]
    assert client.get(f"/api/completed-batches/{batch}/tables/xbt_csv?limit=201").status_code == 422


def test_zip_stream_cancellation_stops_producer_and_keeps_queue_bounded(tmp_path, monkeypatch):
    source = tmp_path / "small.bin"
    source.write_bytes(b"source")
    queues = []
    queue_class = module.queue.Queue

    def queue_factory(*args, **kwargs):
        instance = queue_class(*args, **kwargs)
        queues.append(instance)
        return instance

    def sustained_copy(_source, destination, length):
        for _ in range(10_000):
            destination.write(b"x" * 1024)

    monkeypatch.setattr(module.queue, "Queue", queue_factory)
    monkeypatch.setattr(module.shutil, "copyfileobj", sustained_copy)
    stream = module._zip_stream([(source, "small.bin")])
    assert next(stream)
    stream.close()
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline and any(
        thread.name == "completed-batch-zip" and thread.is_alive()
        for thread in threading.enumerate()
    ):
        time.sleep(0.01)
    assert not any(
        thread.name == "completed-batch-zip" and thread.is_alive()
        for thread in threading.enumerate()
    )
    assert len(queues) == 1
    assert queues[0].maxsize == 8
    assert queues[0].qsize() <= 8


def test_missing_batch_does_not_block_registration_and_routes_are_empty_or_404(tmp_path):
    app = FastAPI()
    service = module.register_completed_batch_routes(
        app,
        batch_root=tmp_path / "absent-batch",
        cache_root=tmp_path / "cache",
        remote_host="unused",
    )
    assert service is None
    client = TestClient(app)
    assert client.get("/api/completed-batches").json() == {"batches": []}
    response = client.get(f"/api/completed-batches/{module.COMPLETED_BATCH_ID}")
    assert response.status_code == 404
    assert "尚未安装" in response.json()["detail"]


def test_corrupt_present_batch_does_not_block_registration_and_routes_return_502(tmp_path):
    batch = tmp_path / "corrupt-batch"
    batch.mkdir()
    (batch / "FINAL_DELIVERY_INDEX.json").write_text("not-json", encoding="utf-8")
    app = FastAPI()
    service = module.register_completed_batch_routes(
        app,
        batch_root=batch,
        cache_root=tmp_path / "cache",
        remote_host="unused",
    )
    assert service is None
    client = TestClient(app)
    for path in (
        "/api/completed-batches",
        f"/api/completed-batches/{module.COMPLETED_BATCH_ID}",
    ):
        response = client.get(path)
        assert response.status_code == 502
        assert "交付无效" in response.json()["detail"]
