import io
from pathlib import Path
import threading
import zipfile

from fastapi import HTTPException
from fastapi.responses import StreamingResponse
import pytest

from web.export_jobs import ExportJobs, report_progress
from web.export_identity import build_export_identity


def response():
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as archive:
        archive.writestr("result.csv", "name,score\nexample,0.9\n")
    return StreamingResponse(iter([output.getvalue()]), headers={"X-MSsystem-Exported-Tasks": "1"})


def test_connection_independent_archive_and_duplicate_submission(tmp_path):
    jobs = ExportJobs(tmp_path)
    started, release = threading.Event(), threading.Event()

    def runner():
        report_progress(completed=0, current_task="sample")
        started.set()
        assert release.wait(5)
        report_progress(completed=1)
        return response()

    request = {"task_ids": ["sample"], "format": "csv_ms2"}
    job = jobs.submit(request, runner)
    assert started.wait(5)
    assert jobs.submit(request, runner)["id"] == job["id"]
    assert jobs.get(job["id"])["current_task"] == "sample"
    release.set()
    jobs.executor.shutdown(wait=True)
    assert jobs.get(job["id"])["status"] == "completed"
    assert jobs.get(job["id"])["completed"] == 1
    with zipfile.ZipFile(jobs.archive_path(job["id"])) as archive:
        assert archive.read("result.csv").startswith(b"name,score")
    # A new service instance can download a completed file without rerunning.
    restarted = ExportJobs(tmp_path)
    assert restarted.archive_path(job["id"]).read_bytes() == jobs.archive_path(job["id"]).read_bytes()
    restarted.executor.shutdown()


def test_failure_retains_real_cause_and_is_not_downloadable(tmp_path):
    jobs = ExportJobs(tmp_path)

    def runner():
        raise HTTPException(status_code=503, detail="SQL evidence invalid")

    job = jobs.submit({"task_ids": ["sample"]}, runner)
    jobs.executor.shutdown(wait=True)
    assert jobs.get(job["id"])["error"] == "SQL evidence invalid"
    with pytest.raises(ValueError):
        jobs.archive_path(job["id"])


def test_interrupted_job_and_path_traversal(tmp_path):
    jobs = ExportJobs(tmp_path)
    job_id = "a" * 32
    jobs._write({"id": job_id, "status": "running"})
    assert jobs.get(job_id)["status"] == "failed"
    with pytest.raises(KeyError):
        jobs.get("../../anything")
    jobs.executor.shutdown()


def test_completed_request_is_reused_without_running_again(tmp_path):
    jobs = ExportJobs(tmp_path)
    calls = []
    request = {"task_ids": ["sample"], "format": "csv_ms2",
               "source_signature": "same-content"}

    def runner():
        calls.append(True)
        return response()

    first = jobs.submit(request, runner)
    jobs.executor.shutdown(wait=True)
    reused = jobs.submit(request, runner)
    assert reused["id"] == first["id"]
    assert reused["cache_hit"] is True
    assert reused["cache_validation_wall_seconds"] >= 0
    assert calls == [True]


def test_content_identity_changes_with_same_stat_and_with_derived_artifact(tmp_path):
    raw = tmp_path / "sample_results.json"
    derived = tmp_path / "sample_compounds.json"
    raw.write_bytes(b"AAAA")
    derived.write_bytes(b"1111")
    raw_stat = raw.stat()
    derived_stat = derived.stat()
    first = build_export_identity(
        {"raw": raw, "derived": derived}, bindings={"language": "zh"}
    )["sha256"]
    raw.write_bytes(b"BBBB")
    raw.touch()
    import os
    os.utime(raw, ns=(raw_stat.st_atime_ns, raw_stat.st_mtime_ns))
    second = build_export_identity(
        {"raw": raw, "derived": derived}, bindings={"language": "zh"}
    )["sha256"]
    assert second != first
    derived.write_bytes(b"2222")
    os.utime(derived, ns=(derived_stat.st_atime_ns, derived_stat.st_mtime_ns))
    third = build_export_identity(
        {"raw": raw, "derived": derived}, bindings={"language": "zh"}
    )["sha256"]
    assert third != second


def test_restart_reuses_only_sha_valid_completed_file(tmp_path):
    request = {"task_ids": ["sample"], "format": "csv_ms2",
               "source_signature": "same-content"}
    jobs = ExportJobs(tmp_path)
    first = jobs.submit(request, response)
    jobs.executor.shutdown(wait=True)

    restarted = ExportJobs(tmp_path)
    assert restarted.submit(request, lambda: pytest.fail("must reuse"))["id"] == first["id"]
    restarted.executor.shutdown(wait=True)

    jobs.archive_path(first["id"]).write_bytes(b"broken")
    repaired = ExportJobs(tmp_path)
    replacement = repaired.submit(request, response)
    repaired.executor.shutdown(wait=True)
    assert replacement["id"] != first["id"]
    assert repaired.get(replacement["id"])["status"] == "completed"


def test_legacy_completed_file_remains_downloadable_but_not_cache_reusable(tmp_path):
    job_id = "d" * 32
    request = {"task_ids": ["sample"], "format": "csv_ms2"}
    jobs = ExportJobs(tmp_path)
    artifact = tmp_path / f"{job_id}.zip"
    with zipfile.ZipFile(artifact, "w") as archive:
        archive.writestr("result.csv", "name,score\nexample,0.9\n")
    jobs._write({
        "id": job_id, "status": "completed", "request": request,
        "size_bytes": artifact.stat().st_size, "media_type": "application/zip",
    })
    jobs.executor.shutdown(wait=True)

    restarted = ExportJobs(tmp_path)
    assert restarted.archive_path(job_id) == artifact
    adopted = restarted.get(job_id)
    assert adopted["artifact_sha256"]
    assert adopted["cache_reusable"] is False
    replacement = restarted.submit(request, response)
    restarted.executor.shutdown(wait=True)
    assert replacement["id"] != job_id


def test_bad_zip_in_completed_index_does_not_crash_restart(tmp_path):
    job_id = "e" * 32
    request = {"task_ids": ["sample"], "format": "csv_ms2"}
    artifact = tmp_path / f"{job_id}.zip"
    artifact.write_bytes(b"not-a-zip")
    state = {
        "id": job_id, "status": "completed", "request": request,
        "size_bytes": artifact.stat().st_size,
        "artifact_sha256": __import__("hashlib").sha256(artifact.read_bytes()).hexdigest(),
        "media_type": "application/zip", "cache_reusable": True,
    }
    path = tmp_path / f"{job_id}.json"
    path.write_text(__import__("json").dumps(state), encoding="utf-8")
    restarted = ExportJobs(tmp_path)
    assert restarted.completed == {}
    restarted.executor.shutdown(wait=True)


def test_part_and_failed_files_are_never_reused_after_restart(tmp_path):
    request = {"task_ids": ["sample"], "format": "csv_ms2"}
    failed_id = "b" * 32
    part_id = "c" * 32
    jobs = ExportJobs(tmp_path)
    jobs._write({"id": failed_id, "status": "failed", "request": request})
    jobs._write({"id": part_id, "status": "running", "request": request})
    (tmp_path / f"{part_id}.zip.part").write_bytes(b"partial")
    jobs.executor.shutdown(wait=True)

    restarted = ExportJobs(tmp_path)
    created = restarted.submit(request, response)
    restarted.executor.shutdown(wait=True)
    assert created["id"] not in {failed_id, part_id}
    assert restarted.get(created["id"])["status"] == "completed"


def test_excel_post_returns_before_slow_refresh_and_rebinds_active_identity(
    monkeypatch, tmp_path
):
    import time
    from fastapi.responses import Response
    from fastapi.testclient import TestClient
    import web.server as server

    jobs = ExportJobs(tmp_path / "jobs")
    refresh_started = threading.Event()
    release_refresh = threading.Event()
    render_started = threading.Event()
    release_render = threading.Event()
    identity_generation = {"value": "before"}
    task_id = "excel-refresh-fixture"
    task = {"id": task_id, "filename": "fixture.mzML", "status": "completed"}

    def identity(_task_ids, _parameters=None):
        return {"sha256": identity_generation["value"],
                "identity_wall_seconds": 0.001}

    def refresh(_task_id, _task):
        refresh_started.set()
        assert release_refresh.wait(5)
        identity_generation["value"] = "after"
        return {}

    def workbook(_task_id, *, evidence_prevalidated=False):
        assert evidence_prevalidated is True
        render_started.set()
        assert release_render.wait(5)
        output = io.BytesIO()
        with zipfile.ZipFile(output, "w") as archive:
            archive.writestr("xl/workbook.xml", "<workbook/>")
        return Response(
            output.getvalue(),
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

    monkeypatch.setattr(server, "_export_jobs", jobs)
    monkeypatch.setattr(server, "_tasks", {task_id: task})
    monkeypatch.setattr(server, "_export_source_identity", identity)
    monkeypatch.setattr(server, "_compound_scope_summary_for_export", refresh)
    monkeypatch.setattr(server, "_download_excel_impl", workbook)
    monkeypatch.setattr(server, "_report_export_lock", threading.RLock())

    client = TestClient(server.app)
    submitted_at = time.perf_counter()
    first = client.post(f"/api/tasks/{task_id}/export/jobs/excel")
    submit_seconds = time.perf_counter() - submitted_at
    assert first.status_code == 202
    assert submit_seconds < 1.0
    first_id = first.json()["id"]
    assert refresh_started.wait(2)
    assert jobs.get(first_id)["status"] == "running"

    release_refresh.set()
    assert render_started.wait(2)
    second = client.post(f"/api/tasks/{task_id}/export/jobs/excel")
    assert second.status_code == 202
    assert second.json()["id"] == first_id

    release_render.set()
    jobs.executor.shutdown(wait=True)
    completed = jobs.get(first_id)
    assert completed["status"] == "completed"
    assert completed["initial_source_signature"] == "before"
    assert completed["final_source_signature"] == "after"
    assert completed["request"]["_preparation_wall_seconds"] >= 0

    third = client.post(f"/api/tasks/{task_id}/export/jobs/excel")
    assert third.json()["id"] == first_id
    assert third.json()["cache_hit"] is True


def test_invalid_archive_is_never_published(tmp_path):
    jobs = ExportJobs(tmp_path)
    job = jobs.submit({"task_ids": ["sample"]},
                      lambda: StreamingResponse(iter([b"not a zip"])))
    jobs.executor.shutdown(wait=True)
    assert jobs.get(job["id"])["status"] == "failed"
    assert not (tmp_path / f"{job['id']}.zip").exists()


def test_invalid_sql_fails_before_loading_any_sample(monkeypatch):
    import web.server as server
    import runtime.verified_sql_retrieval as sql
    monkeypatch.setattr(sql, "sql_metadata_repair_status", lambda: {
        "status": "invalid_binding", "error": "changed dependency"})
    monkeypatch.setattr(server, "_load_task_results_data",
                        lambda *_: pytest.fail("A failed preflight must not load results"))
    with pytest.raises(HTTPException) as error:
        server.batch_export("html_detailed", ["sample"], "zh")
    assert error.value.status_code == 503
    assert "changed dependency" in error.value.detail


def test_historical_download_is_preserved_but_cannot_be_reused(monkeypatch, tmp_path):
    import web.server as server
    jobs = ExportJobs(tmp_path)
    monkeypatch.setattr(server, "_export_jobs", jobs)
    monkeypatch.setattr(server, "_export_source_signature", lambda *_: "generation-two")
    job = jobs.submit({"task_ids": ["sample"], "source_signature": "generation-one"}, response)
    jobs.executor.shutdown(wait=True)
    state = server.get_export_job(job["id"])
    assert state["status"] == "completed"
    assert state["stale"] is True and state["can_reuse"] is False
    download = server.download_export_job(job["id"])
    assert download.status_code == 200
    assert download.headers["X-MSsystem-Export-Stale"] == "true"
    assert Path(download.path).read_bytes() == jobs.archive_path(job["id"]).read_bytes()
    js = Path(server.PROJECT_ROOT, "web/static/js/export.js").read_text()
    assert "job?.status === 'failed' || job?.stale" in js


def test_export_identity_ignores_manifest_timestamp_but_tracks_source(monkeypatch, tmp_path):
    import json
    import web.server as server
    import runtime.verified_sql_retrieval as sql
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(sql, "sql_metadata_dependency_signature", lambda: ("verified",))
    (tmp_path / "config").mkdir()
    manifest = tmp_path / "config/runtime_manifest.json"
    source = tmp_path / "sample_results.json"
    source.write_text("[]")
    manifest.write_text(json.dumps({"runtime_snapshot_id": "same", "generated_at": "before"}))
    before = server._export_source_signature(["sample"])
    manifest.write_text(json.dumps({"runtime_snapshot_id": "same", "generated_at": "after"}))
    assert server._export_source_signature(["sample"]) == before
    source.write_text('[{"changed":true}]')
    assert server._export_source_signature(["sample"]) != before


def test_busy_export_reports_queue_and_legacy_request_fails_fast(monkeypatch, tmp_path):
    import time
    import web.server as server
    jobs = ExportJobs(tmp_path)
    lock = threading.RLock()
    monkeypatch.setattr(server, "_report_export_lock", lock)
    guarded = server._serialized_report_export(response)
    failures = []
    with lock:
        def legacy():
            try:
                guarded()
            except HTTPException as exc:
                failures.append(exc.status_code)
        legacy_thread = threading.Thread(target=legacy)
        legacy_thread.start()
        legacy_thread.join(timeout=2)
        assert not legacy_thread.is_alive()
        assert failures == [409]
        job = jobs.submit({"task_ids": ["sample"]}, guarded)
        deadline = time.monotonic() + 2
        while time.monotonic() < deadline:
            state = jobs.get(job["id"])
            if state.get("phase") == "waiting_export_slot":
                break
            time.sleep(.01)
        assert state["status"] == "queued"
        assert state["phase"] == "waiting_export_slot"
    jobs.executor.shutdown(wait=True)
    assert jobs.get(job["id"])["status"] == "completed"


def test_uninitialized_invalid_library_stops_before_sample_loop(monkeypatch, tmp_path):
    import web.server as server
    import runtime.verified_sql_retrieval as sql
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    (tmp_path / "sample_results.json").write_text("[]")
    monkeypatch.setattr(sql, "sql_metadata_repair_status", lambda: {"status": "not_initialized"})
    checks = []
    def preflight():
        checks.append(True)
        return {"compatible": False, "error": "cache generation mismatch"}
    monkeypatch.setattr(sql, "sql_metadata_repair_preflight", preflight)
    monkeypatch.setattr(server, "_load_task_results_data", lambda *_: pytest.fail("Must stop before loading"))
    with pytest.raises(HTTPException) as error:
        server.batch_export("csv", ["sample"], "zh")
    assert checks == [True]
    assert error.value.status_code == 503

@pytest.mark.parametrize('kind', ['csv', 'html', 'pdf', 'file', 'async'])
def test_native_artifact_response_types(tmp_path, kind):
    from fastapi.responses import Response, FileResponse
    jobs = ExportJobs(tmp_path / 'jobs')
    payload = b'%PDF-1.4\n%%EOF' if kind == 'pdf' else b'name,score\nexample,0.9\n'
    media = {'pdf': 'application/pdf', 'html': 'text/html'}.get(kind, 'text/csv')
    headers = {'Content-Disposition': 'attachment; filename="test.csv"'}
    def runner():
        if kind == 'file':
            source = tmp_path / 'source.csv'; source.write_bytes(payload)
            return FileResponse(source, media_type=media, headers=headers)
        return Response(payload, media_type=media, headers=headers)
    async def async_runner():
        return runner()
    job = jobs.submit({'task_ids': ['sample'], 'format': kind, 'mode': 'single'}, async_runner if kind == 'async' else runner)
    jobs.executor.shutdown(wait=True)
    state = jobs.get(job['id'])
    assert state['status'] == 'completed', state
    assert state['completed'] == 1
    assert jobs.archive_path(job['id']).read_bytes() == payload
    assert state['media_type'] == media


def test_false_pdf_is_not_published(tmp_path):
    from fastapi.responses import Response
    jobs = ExportJobs(tmp_path)
    job = jobs.submit({'task_ids': ['sample'], 'format': 'pdf', 'mode': 'single'},
                      lambda: Response(b'<html>Error</html>', media_type='application/pdf'))
    jobs.executor.shutdown(wait=True)
    assert jobs.get(job['id'])['status'] == 'failed'


def test_sql_signature_identical_before_and_after_verification(tmp_path):
    import json
    from runtime.verified_sql_retrieval import _signature
    from web.export_jobs import canonical_sql_signature
    config = tmp_path / 'config'; config.mkdir()
    registry = config / 'registry.json'
    source = tmp_path / 'source.msp'; source.write_text('unchanged')
    database = tmp_path / 'library.db'; database.write_bytes(b'database')
    registry.write_text(json.dumps({'source_files': {str(source): {}}}))
    policy = config / 'sql_metadata_repairs.json'
    policy.write_text(json.dumps({'enabled': True, 'registry_path': str(registry), 'source_database_path': str(database)}))
    cold = ('', _signature([policy, registry, database]))
    warm = ('', _signature([policy, registry, database, source, Path(str(database)+'-wal'), Path(str(database)+'-journal')]))
    before = canonical_sql_signature(cold, tmp_path)
    assert before == canonical_sql_signature(warm, tmp_path)
    source.write_text('changed source')
    assert canonical_sql_signature(cold, tmp_path) != before


def test_png_download_remains_png_without_native_rdkit_cairo(monkeypatch):
    import base64
    from PIL import Image
    import web.structure as structure
    pytest.importorskip('cairosvg')
    def unavailable(*args, **kwargs):
        raise RuntimeError('RDKit was built without Cairo')
    monkeypatch.setattr(structure.Draw, 'MolToImage', unavailable)
    data = base64.b64decode(structure.smiles_to_image_base64('CCO', size=(160, 120)))
    assert data.startswith(b'\x89PNG\r\n\x1a\n')
    with Image.open(io.BytesIO(data)) as image:
        assert image.size == (160, 120)
        lo, hi = image.convert('L').getextrema()
        assert lo < hi
