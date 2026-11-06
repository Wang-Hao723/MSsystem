"""Batch task execution and ion-mode safety regression tests."""

import hashlib
import json
import time

import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient

from web.server import AnalysisParams, _resolve_spectrum_polarities, app


def test_blank_evidence_defaults_to_soft_identify_first_policy():
    params = AnalysisParams(
        blank_files=["1234abcd_blank.mzML"],
        blank_roles={"1234abcd_blank.mzML": "method_blank"},
    )

    assert params.blank_filter_policy == "soft_evidence"
    assert params.blank_clean_max_ratio == pytest.approx(1.0)
    assert params.blank_clean_min_replicates == 2
    assert params.blank_roles == {
        "1234abcd_blank.mzML": "method_blank",
    }
    assert params.blank_method_comparability_confirmed is False
    assert params.blank_background_exclusion_confirmed is False


def test_hard_blank_exclusion_requires_double_confirmation_and_typed_role():
    from pydantic import ValidationError

    base = {
        "blank_files": ["1234abcd_blank.mzML"],
        "blank_roles": {"1234abcd_blank.mzML": "method_blank"},
        "blank_filter_policy": "hard_legacy",
    }
    with pytest.raises(ValidationError, match="explicit confirmation"):
        AnalysisParams(**base)

    confirmed = AnalysisParams(
        **base,
        blank_method_comparability_confirmed=True,
        blank_background_exclusion_confirmed=True,
    )
    assert confirmed.blank_filter_policy == "hard_legacy"

    generic_background = AnalysisParams(
        blank_files=["1234abcd_blank.mzML"],
        blank_roles={"1234abcd_blank.mzML": "background_blank"},
        blank_filter_policy="hard_legacy",
        blank_method_comparability_confirmed=True,
        blank_background_exclusion_confirmed=True,
    )
    assert generic_background.blank_roles == {
        "1234abcd_blank.mzML": "background_blank",
    }

    with pytest.raises(ValidationError, match="accepts only explicitly typed"):
        AnalysisParams(
            blank_files=["1234abcd_blank.mzML"],
            blank_roles={"1234abcd_blank.mzML": "unspecified"},
            blank_filter_policy="hard_legacy",
            blank_method_comparability_confirmed=True,
            blank_background_exclusion_confirmed=True,
        )


def test_browser_blank_workflow_does_not_silently_force_zero_exclusion():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    html = (root / "web/static/index.html").read_text(encoding="utf-8")
    script = (root / "web/static/js/upload.js").read_text(encoding="utf-8")

    assert '<option value="hard_legacy" selected>' in html
    assert 'id="blankMethodComparable"' in html
    assert 'id="blankExclusionConfirmed"' in html
    assert "Object.assign(params, _blankAnalysisParameters())" in script
    assert "function analysisParamsForRun(" in script
    assert "params.blank_filter_policy = 'soft_evidence'" not in script
    assert "blank_method_comparability_confirmed: hard" in script
    assert "blank_background_exclusion_confirmed: hard" in script


def test_auto_polarity_fills_unknown_scans_from_unambiguous_file_consensus():
    spectra = [{"polarity": "NEG"}, {"polarity": None}, {"polarity": "NEG"}]

    detected, counts = _resolve_spectrum_polarities(spectra, "auto")

    assert detected == "NEG"
    assert counts == {"POS": 0, "NEG": 3, "unknown": 0}
    assert spectra[1]["polarity"] == "NEG"
    assert spectra[1]["polarity_source"] == "file_consensus"


def test_auto_polarity_rejects_fully_unknown_input():
    with pytest.raises(ValueError, match="could not be identified"):
        _resolve_spectrum_polarities([{"polarity": None}], "auto")


def test_auto_polarity_rejects_ambiguous_mixed_file():
    spectra = [{"polarity": "POS"}, {"polarity": "NEG"}, {"polarity": None}]
    with pytest.raises(ValueError, match="mixed-polarity"):
        _resolve_spectrum_polarities(spectra, "auto")


def test_explicit_polarity_overrides_every_scan():
    spectra = [{"polarity": "POS"}, {"polarity": None}]

    detected, counts = _resolve_spectrum_polarities(spectra, "NEG")

    assert detected == "NEG"
    assert counts == {"NEG": 2, "unknown": 0}
    assert all(s["polarity"] == "NEG" for s in spectra)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cosine_threshold": -0.1},
        {"final_score_threshold": 1.1},
        {"mz_tolerance_da": 0},
        {"mz_tolerance_ppm": -1},
        {"ms2_feature_link_tolerance_ppm": 0},
        {"top_n_peaks": 10000},
        {"blank_mode": "silent_skip"},
        {"blank_intensity_mode": "median"},
        {"db_mode": "legacy"},
        {"ion_mode": "mixed"},
    ],
)
def test_analysis_parameter_ranges_are_validated(kwargs):
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        AnalysisParams(**kwargs)


def test_batch_analyze_route_precedes_dynamic_task_route():
    paths = [route.path for route in app.routes]
    assert "/api/tasks/batch/analyze" in paths
    assert paths.index("/api/tasks/batch/analyze") < paths.index("/api/tasks/{task_id}/analyze")


def test_web_exposes_batch_run_and_ion_mode_controls():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    html = (root / "web/static/index.html").read_text(encoding="utf-8")
    tasks_js = (root / "web/static/js/tasks.js").read_text(encoding="utf-8")
    upload_js = (root / "web/static/js/upload.js").read_text(encoding="utf-8")

    assert "<h1>MSsystem</h1>" in html
    assert 'id="batchIonMode"' in html
    assert 'id="btnBatchRun"' in html
    assert "'/tasks/batch/analyze'" in tasks_js
    assert "Upload.analysisParamsForRun('batchIonMode')" in tasks_js
    assert "params," in tasks_js
    assert "blankRunLabel(params)" in tasks_js
    assert 'id="fileInput"' in html and " multiple" in html
    assert 'id="folderInput"' in html and "webkitdirectory" in html
    assert "['fileInput', 'folderInput']" in upload_js
    assert "task_ids: uploaded.map" in upload_js
    assert 'id="runtimeInstanceBadge"' in html
    assert 'id="ms2FeatureLinkTolerancePpm"' in html
    assert "ms2_feature_link_tolerance_ppm: numberValue(" in upload_js
    assert "service_generation" in upload_js
    assert "const expectedGeneration = 'V16'" in upload_js


def _patch_isolated_upload_state(monkeypatch, tmp_path):
    """Keep upload integration tests away from the user's real task index."""
    import web.server as server

    upload_dir = tmp_path / "uploads"
    task_dir = tmp_path / "tasks"
    upload_dir.mkdir()
    task_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASK_DIR", task_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks_index.json")
    monkeypatch.setattr(server, "_tasks", {})
    return server, upload_dir


def test_multi_file_upload_to_batch_queue_real_api_chain(monkeypatch, tmp_path):
    """Mixed-format multipart files become tasks and one batch call queues all."""
    server, upload_dir = _patch_isolated_upload_state(monkeypatch, tmp_path)
    queued = []

    async def fake_start_analysis(task_id, params):
        queued.append((task_id, params.db_mode, params.ion_mode))
        return {"task_id": task_id, "status": "started"}

    monkeypatch.setattr(server, "start_analysis", fake_start_analysis)
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/upload",
            files=[
                ("files", ("sample_A.mgf", b"BEGIN IONS\nPEPMASS=101\n50 1\nEND IONS\n", "text/plain")),
                ("files", ("sample_B.mzML", b"<mzML></mzML>", "application/xml")),
                ("files", ("sample_C.html", "<html>前体 m/z=101; MS2 碎片</html>".encode(), "text/html")),
            ],
        )
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["count"] == 3
        assert body["failed_count"] == 0
        assert [item["filename"] for item in body["uploaded"]] == [
            "sample_A.mgf", "sample_B.mzML", "sample_C.html",
        ]
        assert len(list(upload_dir.iterdir())) == 3

        task_ids = [item["task_id"] for item in body["uploaded"]]
        batch = client.post(
            "/api/tasks/batch/analyze",
            json={
                "task_ids": task_ids,
                "params": {"db_mode": "full", "ion_mode": "auto"},
            },
        )
        assert batch.status_code == 200, batch.text
        assert batch.json()["queued"] == task_ids
        assert batch.json()["queued_count"] == 3
        assert queued == [
            (task_ids[0], "full", "auto"),
            (task_ids[1], "full", "auto"),
            (task_ids[2], "full", "auto"),
        ]
    finally:
        client.close()


def test_batch_upload_isolates_bad_file_and_sanitizes_filename(monkeypatch, tmp_path):
    """One bad item does not discard its valid neighbours or escape upload dir."""
    server, upload_dir = _patch_isolated_upload_state(monkeypatch, tmp_path)
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/upload",
            files=[
                ("files", ("../../valid.mgf", b"BEGIN IONS\nEND IONS\n", "text/plain")),
                ("files", ("ordinary.html", b"<html>not a spectrum report</html>", "text/html")),
            ],
        )
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["count"] == 1
        assert body["failed_count"] == 1
        assert body["uploaded"][0]["filename"] == "valid.mgf"
        assert body["failed"][0]["filename"] == "ordinary.html"
        stored = list(upload_dir.iterdir())
        assert len(stored) == 1
        assert stored[0].parent == upload_dir
        assert not any(path.name.startswith(".") for path in upload_dir.iterdir())
    finally:
        client.close()


def test_batch_upload_enforces_aggregate_byte_limit(monkeypatch, tmp_path):
    server, upload_dir = _patch_isolated_upload_state(monkeypatch, tmp_path)
    monkeypatch.setattr(server, "MAX_BATCH_UPLOAD_SIZE_BYTES", 50)
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/upload",
            files=[
                ("files", ("first.mgf", b"x" * 30, "text/plain")),
                ("files", ("second.mgf", b"y" * 30, "text/plain")),
            ],
        )
        assert response.status_code == 200, response.text
        assert response.json()["count"] == 1
        assert response.json()["failed_count"] == 1
        assert "Batch exceeds" in response.json()["failed"][0]["error"]
        assert len(list(upload_dir.iterdir())) == 1
    finally:
        client.close()


def test_batch_queue_returns_structured_failures_without_forcing_reupload(monkeypatch, tmp_path):
    """A queue rejection stays inspectable instead of turning upload into a retry loop."""
    server, _ = _patch_isolated_upload_state(monkeypatch, tmp_path)
    server._tasks.update({
        "a": {"id": "a", "status": "pending"},
        "b": {"id": "b", "status": "pending"},
    })

    async def reject_start(task_id, params):
        raise HTTPException(status_code=409, detail=f"cannot queue {task_id}")

    monkeypatch.setattr(server, "start_analysis", reject_start)
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/tasks/batch/analyze",
            json={"task_ids": ["a", "b"], "params": {"db_mode": "full"}},
        )
        assert response.status_code == 200, response.text
        assert response.json()["queued_count"] == 0
        assert [item["id"] for item in response.json()["failed"]] == ["a", "b"]
    finally:
        client.close()


def test_batch_queue_isolates_unexpected_per_task_error(monkeypatch, tmp_path):
    server, _ = _patch_isolated_upload_state(monkeypatch, tmp_path)
    server._tasks.update({
        "a": {"id": "a", "status": "pending"},
        "b": {"id": "b", "status": "pending"},
    })

    async def mixed_start(task_id, params):
        if task_id == "a":
            raise RuntimeError("unexpected queue error")
        return {"task_id": task_id, "status": "queued"}

    monkeypatch.setattr(server, "start_analysis", mixed_start)
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/tasks/batch/analyze",
            json={"task_ids": ["a", "b"], "params": {"db_mode": "full"}},
        )
        assert response.status_code == 200, response.text
        assert response.json()["queued"] == ["b"]
        assert response.json()["failed"][0]["id"] == "a"
    finally:
        client.close()


def test_start_analysis_claims_task_before_executor_submit(monkeypatch, tmp_path):
    """Two immediate starts cannot submit the same task twice."""
    server, upload_dir = _patch_isolated_upload_state(monkeypatch, tmp_path)
    sample = upload_dir / "sample.mgf"
    sample.write_text("BEGIN IONS\nPEPMASS=101\n50 1\nEND IONS\n")
    server._tasks["one"] = {
        "id": "one", "filename": "sample.mgf", "file_path": str(sample),
        "status": "pending", "progress": 0, "message": "", "results": [],
    }
    submitted = []

    class FakeExecutor:
        def submit(self, fn, task_id):
            submitted.append((fn, task_id))

    monkeypatch.setattr(server, "_lib_loaded", True)
    monkeypatch.setattr(server, "_lib_name", "Full Library")
    monkeypatch.setattr(server, "_executor_full", FakeExecutor())
    client = TestClient(server.app)
    try:
        first = client.post(
            "/api/tasks/one/analyze",
            json={"db_mode": "full", "ion_mode": "NEG"},
        )
        second = client.post(
            "/api/tasks/one/analyze",
            json={"db_mode": "full", "ion_mode": "NEG"},
        )
        assert first.status_code == 200, first.text
        assert first.json()["status"] == "queued"
        assert second.status_code == 409, second.text
        assert len(submitted) == 1
        assert submitted[0][1] == "one"
        assert server._tasks["one"]["message"] == (
            "Queued for Full-library serial analysis..."
        )
    finally:
        client.close()


def test_long_stage_heartbeat_reports_elapsed_time_without_fake_progress(
    monkeypatch,
):
    import web.server as server

    events = []
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            "heartbeat": {
                "id": "heartbeat",
                "status": "running",
                "progress": 3,
                "message": "",
            }
        },
    )
    monkeypatch.setattr(server, "_save_tasks", lambda **_kwargs: True)
    monkeypatch.setattr(
        server,
        "_sse_broadcast",
        lambda event, payload: events.append((event, payload)),
    )

    stop, thread = server._start_task_stage_heartbeat(
        "heartbeat",
        stage="openms_feature_detection",
        message="Detecting MS1 chromatographic features.",
        timeout_seconds=300,
        interval_seconds=0.01,
    )
    time.sleep(0.04)
    stop.set()
    thread.join(timeout=1.0)

    task = server._tasks["heartbeat"]
    assert task["progress"] == 3
    assert "Still active" in task["message"]
    assert "bounded at 300s" in task["message"]
    assert task["stage_heartbeat"]["stage"] == "openms_feature_detection"
    assert events
    assert all(payload["progress"] == 3 for _, payload in events)


def test_invalid_blank_upload_is_rejected_and_deleted(monkeypatch, tmp_path):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/blank/upload",
            files={"file": ("../../bad.mgf", b"not an mgf", "text/plain")},
        )
        assert response.status_code == 400, response.text
        assert list(blank_dir.iterdir()) == []
    finally:
        client.close()


def test_valid_blank_upload_returns_context_peak_count(monkeypatch, tmp_path):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(
        server, "BLANK_UPLOADS_FILE", tmp_path / "blank_uploads_index.json"
    )
    mgf = (
        b"BEGIN IONS\nTITLE=blank\nPEPMASS=200\nRTINSECONDS=60\n"
        b"50 100\n75 50\n100 25\nEND IONS\n"
    )
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/blank/upload",
            files={"file": ("../../blank.mgf", mgf, "text/plain")},
        )
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["blank_filename"] == "blank.mgf"
        assert body["n_peaks"] == 3
        assert body["n_ms1_scans"] == 0
        assert body["n_ms2_contexts"] == 1
        assert body["polarities"] == []
        assert body["blank_role"] == "background_blank"
        assert body["content_sha256"] == hashlib.sha256(mgf).hexdigest()
        assert "/" not in body["blank_path"]
        assert server._resolve_blank_file_token(body["blank_path"]).is_file()
        assert len(list(blank_dir.iterdir())) == 1
    finally:
        client.close()


def test_byte_identical_managed_blanks_are_one_analytical_replicate(
    monkeypatch, tmp_path,
):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    blank_index = tmp_path / "blank_uploads_index.json"
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(server, "BLANK_UPLOADS_FILE", blank_index)
    mgf = (
        b"BEGIN IONS\nTITLE=blank\nPEPMASS=200\nRTINSECONDS=60\n"
        b"50 100\n75 50\n100 25\nEND IONS\n"
    )

    client = TestClient(server.app)
    try:
        first_response = client.post(
            "/api/blank/upload",
            files={"file": ("blank-a.mgf", mgf, "text/plain")},
        )
        second_response = client.post(
            "/api/blank/upload",
            files={"file": ("blank-b.mgf", mgf, "text/plain")},
        )
        assert first_response.status_code == 200, first_response.text
        assert second_response.status_code == 200, second_response.text
        first = first_response.json()
        second = second_response.json()
        assert first["content_sha256"] == second["content_sha256"]

        unique, ignored = server._deduplicate_managed_blank_tokens([
            first["blank_path"], second["blank_path"],
        ])
        assert unique == [first["blank_path"]]
        assert ignored == [{
            "blank_path": second["blank_path"],
            "duplicate_of": first["blank_path"],
            "content_sha256": first["content_sha256"],
        }]
    finally:
        client.close()


def test_uploaded_blank_can_be_deleted_by_managed_token(monkeypatch, tmp_path):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    token = "1234abcd_blank.mgf"
    blank_path = blank_dir / token
    blank_path.write_text("BEGIN IONS\nEND IONS\n")
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(
        server, "BLANK_UPLOADS_FILE", tmp_path / "blank_uploads_index.json"
    )
    monkeypatch.setattr(server, "_tasks", {})

    client = TestClient(server.app)
    try:
        response = client.delete(f"/api/blanks/{token}")
        assert response.status_code == 200, response.text
        assert response.json()["deleted"] == token
        assert not blank_path.exists()
    finally:
        client.close()


def test_uploaded_blank_role_is_persisted_and_listed(monkeypatch, tmp_path):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    token = "1234abcd_blank.mgf"
    (blank_dir / token).write_text("BEGIN IONS\nEND IONS\n")
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(
        server, "BLANK_UPLOADS_FILE", tmp_path / "blank_uploads_index.json"
    )
    monkeypatch.setattr(server, "_tasks", {})

    client = TestClient(server.app)
    try:
        response = client.post(
            f"/api/blanks/{token}/role",
            json={"blank_role": "carryover_blank"},
        )
        assert response.status_code == 200, response.text
        assert response.json()["blank_role"] == "carryover_blank"

        listed = client.get("/api/blanks")
        assert listed.status_code == 200, listed.text
        assert listed.json()["blanks"][0]["blank_role"] == "carryover_blank"

        payload = json.loads(
            (tmp_path / "blank_uploads_index.json").read_text()
        )
        assert payload["schema_version"] == 2
        assert payload["uploads"][token]["blank_role"] == "carryover_blank"
    finally:
        client.close()


def test_uploaded_blank_cannot_be_deleted_while_active_task_uses_it(
    monkeypatch, tmp_path,
):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    token = "1234abcd_blank.mgf"
    blank_path = blank_dir / token
    blank_path.write_text("BEGIN IONS\nEND IONS\n")
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(
        server, "BLANK_UPLOADS_FILE", tmp_path / "blank_uploads_index.json"
    )
    monkeypatch.setattr(server, "_tasks", {
        "active": {
            "id": "active",
            "status": "queued",
            "params": {"blank_files": [token]},
        }
    })

    client = TestClient(server.app)
    try:
        response = client.delete(f"/api/blanks/{token}")
        assert response.status_code == 409, response.text
        assert blank_path.exists()
    finally:
        client.close()


def test_blank_uploads_remain_editable_after_page_refresh(monkeypatch, tmp_path):
    """Legacy/current server files must be listable without browser memory."""
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    token = "1234abcd_NEG_blank.mgf"
    blank_path = blank_dir / token
    blank_path.write_text("BEGIN IONS\nEND IONS\n")
    blank_index = tmp_path / "blank_uploads_index.json"
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(server, "BLANK_UPLOADS_FILE", blank_index)
    monkeypatch.setattr(server, "_tasks", {
        "active": {
            "id": "active",
            "status": "queued",
            "params": {"blank_files": [token]},
        }
    })

    client = TestClient(server.app)
    try:
        response = client.get("/api/blanks")
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["count"] == 1
        assert body["blanks"][0]["blank_path"] == token
        assert body["blanks"][0]["blank_filename"] == "NEG_blank.mgf"
        assert body["blanks"][0]["in_use"] is True
        assert body["blanks"][0]["active_task_ids"] == ["active"]
        assert blank_index.is_file()

        server._tasks = {}
        deleted = client.delete(f"/api/blanks/{token}")
        assert deleted.status_code == 200, deleted.text
        assert not blank_path.exists()
        assert client.get("/api/blanks").json()["count"] == 0
    finally:
        client.close()


def test_restart_recovers_queued_and_running_but_preserves_pending_upload(
    monkeypatch, tmp_path,
):
    import json
    import web.server as server

    queued_file = tmp_path / "queued.mzML"
    running_file = tmp_path / "running.mzML"
    finalizing_file = tmp_path / "finalizing.mzML"
    for path in (queued_file, running_file, finalizing_file):
        path.write_text("<mzML/>", encoding="utf-8")
    tasks_file = tmp_path / "tasks_index.json"
    payload = {
        "queued": {
            "id": "queued", "status": "queued",
            "file_path": str(queued_file), "params": {"db_mode": "full"},
        },
        "running": {
            "id": "running", "status": "running",
            "file_path": str(running_file), "params": {"db_mode": "full"},
        },
        "finalizing": {
            "id": "finalizing", "status": "finalizing",
            "file_path": str(finalizing_file), "params": {"db_mode": "small"},
        },
        "pending": {
            "id": "pending", "status": "pending", "file_path": "",
        },
    }
    tasks_file.write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(server, "TASKS_FILE", tasks_file)
    monkeypatch.setattr(server, "_tasks", {})

    server._load_tasks()

    assert server._tasks["queued"]["status"] == "pending"
    assert server._tasks["running"]["status"] == "pending"
    assert server._tasks["finalizing"]["status"] == "pending"
    assert server._tasks["pending"]["status"] == "pending"
    assert server._tasks["queued"]["error"] is None
    assert "自动恢复排队" in server._tasks["queued"]["message"]
    assert server._restart_recovery_task_ids == {
        "queued", "running", "finalizing",
    }


def test_restart_marks_active_task_failed_only_when_managed_input_is_missing(
    monkeypatch, tmp_path,
):
    import json
    import web.server as server

    tasks_file = tmp_path / "tasks_index.json"
    tasks_file.write_text(json.dumps({
        "missing": {
            "id": "missing",
            "status": "running",
            "file_path": str(tmp_path / "missing.mzML"),
        },
    }), encoding="utf-8")
    monkeypatch.setattr(server, "TASKS_FILE", tasks_file)
    monkeypatch.setattr(server, "_tasks", {})

    server._load_tasks()

    assert server._tasks["missing"]["status"] == "failed"
    assert "原始上传文件缺失" in server._tasks["missing"]["error"]


def test_task_registry_refuses_unconfirmed_empty_overwrite_and_keeps_backup(
    monkeypatch, tmp_path,
):
    import json
    import web.server as server

    tasks_file = tmp_path / "tasks_index.json"
    tasks_file.write_text(json.dumps({
        "kept": {"id": "kept", "status": "pending", "file_path": ""},
    }), encoding="utf-8")
    monkeypatch.setattr(server, "TASKS_FILE", tasks_file)
    monkeypatch.setattr(server, "_tasks", {})

    assert server._save_tasks() is False
    assert "kept" in json.loads(tasks_file.read_text(encoding="utf-8"))
    assert server._save_tasks(allow_empty=True) is True
    assert json.loads(tasks_file.read_text(encoding="utf-8")) == {}
    backup = tmp_path / "tasks_index.backup.json"
    assert "kept" in json.loads(backup.read_text(encoding="utf-8"))


def test_corrupt_primary_task_registry_recovers_from_last_good_backup(
    monkeypatch, tmp_path,
):
    import json
    import web.server as server

    tasks_file = tmp_path / "tasks_index.json"
    tasks_file.write_text("{broken", encoding="utf-8")
    backup = tmp_path / "tasks_index.backup.json"
    backup.write_text(json.dumps({
        "restored": {
            "id": "restored", "status": "pending", "file_path": "",
        },
    }), encoding="utf-8")
    monkeypatch.setattr(server, "TASKS_FILE", tasks_file)
    monkeypatch.setattr(server, "_tasks", {})

    server._load_tasks()

    assert "restored" in server._tasks
    assert json.loads(tasks_file.read_text(encoding="utf-8"))["restored"][
        "status"
    ] == "pending"


def test_loading_requested_library_automatically_requeues_restart_recovery(
    monkeypatch, tmp_path,
):
    import asyncio
    import web.server as server

    sample = tmp_path / "resume.mzML"
    sample.write_text("<mzML/>", encoding="utf-8")
    task = {
        "id": "resume",
        "status": "pending",
        "file_path": str(sample),
        "params": {"db_mode": "full", "ion_mode": "NEG"},
        "restart_recovery": {"status": "awaiting_library"},
    }
    monkeypatch.setattr(server, "_tasks", {"resume": task})
    monkeypatch.setattr(server, "_restart_recovery_task_ids", {"resume"})
    monkeypatch.setattr(server, "_save_tasks", lambda **_kwargs: True)
    queued = []

    async def fake_start(task_id, params):
        queued.append((task_id, params.db_mode, params.ion_mode))
        server._tasks[task_id]["status"] = "queued"

    monkeypatch.setattr(server, "start_analysis", fake_start)

    result = asyncio.run(server._resume_recovered_tasks("full"))

    assert result == {"resumed": ["resume"], "deferred": [], "failed": []}
    assert queued == [("resume", "full", "NEG")]
    assert server._tasks["resume"]["restart_recovery"]["status"] == "resumed"
    assert server._restart_recovery_task_ids == set()


def test_analysis_fails_closed_before_search_when_blank_polarity_is_opposite(
    monkeypatch, tmp_path,
):
    import matching.blank_subtraction as blank_module
    import web.server as server

    sample_path = tmp_path / "sample.bin"
    sample_path.write_bytes(b"placeholder")
    blank_path = tmp_path / "blank.mzML"
    blank_path.write_bytes(b"placeholder")
    neg_profile = blank_module.create_blank_profile_from_list([(100.0, 100.0)])
    neg_profile.available_polarities = frozenset({"NEG"})
    neg_profile.ms1_scan_count = 1

    monkeypatch.setattr(server, "_tasks", {
        "polarity-mismatch": {
            "id": "polarity-mismatch",
            "filename": sample_path.name,
            "file_path": str(sample_path),
            "status": "queued",
            "params": {"ion_mode": "auto", "blank_files": ["blank-token"]},
            "results": [],
        }
    })
    monkeypatch.setattr(
        server,
        "_parse_spectrum_file",
        lambda *_args, **_kwargs: [{
            "spectrum_id": "scan=1",
            "precursor_mz": 279.159,
            "retention_time": 10.0,
            "retention_time_unit": "minutes",
            "polarity": "POS",
            "peaks": [[149.0, 1000.0], [205.0, 250.0]],
        }],
    )
    monkeypatch.setattr(server, "_resolve_blank_file_token", lambda _token: blank_path)
    monkeypatch.setattr(blank_module, "load_blank_file", lambda _path: neg_profile)
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_sse_broadcast", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        server,
        "search",
        lambda *_args, **_kwargs: pytest.fail("search must not run after polarity mismatch"),
    )

    server._run_analysis("polarity-mismatch")

    task = server._tasks["polarity-mismatch"]
    assert task["status"] == "failed"
    assert task["blank_summary"]["applicable"] is False
    assert task["blank_summary"]["sample_polarities"] == ["POS"]
    assert task["blank_summary"]["blank_polarities"] == ["NEG"]
    assert "缺少同极性空白" in task["error"]


def test_blank_filtered_spectra_are_exported_in_separate_audit_sections(
    monkeypatch, tmp_path,
):
    import io

    from openpyxl import load_workbook

    server, _ = _patch_isolated_upload_state(monkeypatch, tmp_path)
    server._tasks["blank-audit"] = {
        "id": "blank-audit",
        "filename": "sample.mzML",
        "status": "completed",
        "total_spectra": 1,
        "matched_spectra": 0,
        "match_rate": 0.0,
        "results": [{
            "spectrum_id": "scan=7",
            "precursor_mz": 500.0,
            "retention_time": 5.0,
            "peaks_count": 3,
            "match": None,
            "blank_excluded": True,
            "blank_status": "excluded_below_ratio",
            "blank_evidence_source": "blank_ms1",
            "blank_sample_intensity": 200.0,
            "blank_intensity": 100.0,
            "blank_sample_ratio": 2.0,
            "blank_mz_error_ppm": 1.2,
            "blank_rt_delta_minutes": 0.05,
        }],
    }
    client = TestClient(server.app)
    try:
        csv_response = client.get(
            "/api/tasks/blank-audit/download/csv?format_type=report"
        )
        assert csv_response.status_code == 200
        assert "空白特征过滤排除" in csv_response.text
        assert "blank_ms1" in csv_response.text
        assert "2.0000" in csv_response.text

        excel_response = client.get("/api/tasks/blank-audit/download/excel")
        assert excel_response.status_code == 200
        workbook = load_workbook(io.BytesIO(excel_response.content), read_only=True)
        assert "08_空白过滤" in workbook.sheetnames
        summary_values = {
            workbook["01_结果导航"].cell(row, 1).value:
            workbook["01_结果导航"].cell(row, 2).value
            for row in range(2, workbook["01_结果导航"].max_row + 1)
        }
        assert summary_values["空白过滤排除谱图"] == 1
    finally:
        client.close()


def test_soft_blank_views_keep_full_identification_and_export_clean_sheet(
    monkeypatch, tmp_path,
):
    import io

    from openpyxl import load_workbook

    server, _ = _patch_isolated_upload_state(monkeypatch, tmp_path)

    def matched_row(scan_id, name, inchikey, *, clean_excluded):
        return {
            "spectrum_id": scan_id,
            "precursor_mz": 500.0,
            "retention_time": 5.0,
            "peaks_count": 3,
            "match": {
                "compound_name": name,
                "inchikey": inchikey,
                "smiles": "CCO",
                "formula": "C2H6O",
                "cosine": 0.90,
                "final_score": 0.70,
                "rank_score": 0.70,
                "confidence": "High",
                "confidence_class": "High",
                "source_tier": "experimental",
                "dataset_source": "test",
            },
            "blank_feature_matched": True,
            "blank_flagged": clean_excluded,
            "blank_clean_view_excluded": clean_excluded,
            "blank_excluded": False,
            "blank_status": (
                "retained_blank_dominant"
                if clean_excluded else "retained_sample_enriched"
            ),
            "blank_evidence_category": (
                "blank_dominant" if clean_excluded else "sample_enriched"
            ),
            "blank_evidence_role": "method_blank",
            "blank_sample_ratio": 0.8 if clean_excluded else 4.0,
        }

    server._tasks["soft-blank"] = {
        "id": "soft-blank",
        "filename": "sample.mzML",
        "status": "completed",
        "total_spectra": 2,
        "matched_spectra": 2,
        "match_rate": 100.0,
        "results": [
            matched_row(
                "scan=1", "Blank-dominant target",
                "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                clean_excluded=True,
            ),
            matched_row(
                "scan=2", "Sample-enriched target",
                "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                clean_excluded=False,
            ),
        ],
    }
    client = TestClient(server.app)
    try:
        all_view = client.get(
            "/api/tasks/soft-blank/results?limit=10&blank_view=all"
        )
        clean_view = client.get(
            "/api/tasks/soft-blank/results?limit=10&blank_view=clean"
        )
        flagged_view = client.get(
            "/api/tasks/soft-blank/results?limit=10&blank_view=flagged"
        )
        assert all_view.json()["total"] == 2
        assert clean_view.json()["total"] == 1
        assert clean_view.json()["results"][0]["match"]["compound_name"] == (
            "Sample-enriched target"
        )
        assert flagged_view.json()["total"] == 1
        assert flagged_view.json()["results"][0]["match"]["compound_name"] == (
            "Blank-dominant target"
        )

        full_csv = client.get(
            "/api/tasks/soft-blank/download/csv?blank_view=all"
        )
        clean_csv = client.get(
            "/api/tasks/soft-blank/download/csv?blank_view=clean"
        )
        assert "Blank-dominant target" in full_csv.text
        assert "Blank-dominant target" not in clean_csv.text.split(
            "=== 空白证据标记"
        )[0]
        assert "Sample-enriched target" in clean_csv.text

        excel_response = client.get("/api/tasks/soft-blank/download/excel")
        workbook = load_workbook(
            io.BytesIO(excel_response.content), read_only=True
        )
        assert "03_MS2谱库鉴定" in workbook.sheetnames
        assert "07_MS2清洁视图" in workbook.sheetnames
        assert "09_空白证据" in workbook.sheetnames
        summary_values = {
            workbook["01_结果导航"].cell(row, 1).value:
            workbook["01_结果导航"].cell(row, 2).value
            for row in range(2, workbook["01_结果导航"].max_row + 1)
        }
        assert summary_values["空白清洁视图 MS2 结果"] == 1
    finally:
        client.close()


def test_library_invalidation_clears_web_and_search_caches(monkeypatch):
    import importlib
    import web.server as server

    search_module = importlib.import_module("matching.search")

    monkeypatch.setattr(server, "_lib_records", [{"old": True}])
    monkeypatch.setattr(server, "_mz_index", {1: [0]})
    monkeypatch.setattr(server, "_lib_loaded", True)
    monkeypatch.setattr(server, "_lib_name", "Full (old)")
    monkeypatch.setattr(server, "_lib_spectra_count", 1)
    monkeypatch.setattr(search_module, "_lib_records", [{"old": True}])
    monkeypatch.setattr(search_module, "_mz_index", {1: [0]})
    monkeypatch.setattr(search_module, "_nm_index", {1: [0]})
    monkeypatch.setattr(search_module, "_lib_db_path", __file__)

    server._invalidate_library_state()

    assert server._lib_loaded is False
    assert server._lib_records == []
    assert server._mz_index == {}
    assert server._lib_name == ""
    assert server._lib_spectra_count == 0
    assert search_module._lib_records is None
    assert search_module._mz_index is None
    assert search_module._nm_index is None
    assert search_module._lib_db_path is None


def test_blank_token_cannot_escape_managed_directory(monkeypatch, tmp_path):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    outside = tmp_path / "outside.mgf"
    outside.write_text("secret")
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)

    with pytest.raises(ValueError, match="outside"):
        server._resolve_blank_file_token(str(outside))


def test_active_task_cannot_be_deleted(monkeypatch, tmp_path):
    server, upload_dir = _patch_isolated_upload_state(monkeypatch, tmp_path)
    sample = upload_dir / "running.mgf"
    sample.write_text("BEGIN IONS\nEND IONS\n")
    server._tasks["running"] = {
        "id": "running", "filename": sample.name,
        "file_path": str(sample), "status": "running",
    }
    client = TestClient(server.app)
    try:
        response = client.delete("/api/tasks/running")
        assert response.status_code == 409, response.text
        assert "running" in server._tasks
        assert sample.exists()
    finally:
        client.close()


def test_library_preview_uses_atomic_temp_files_and_cleans_them(monkeypatch, tmp_path):
    import matching.library_import as library_import
    import web.server as server

    raw_dir = tmp_path / "user_contrib" / "raw"
    monkeypatch.setattr(server, "USER_CONTRIB_RAW_DIR", raw_dir)
    monkeypatch.setattr(
        library_import,
        "preview_import",
        lambda paths: {
            "quality": {"accepted": 1, "rejected": 0},
            "spectral_records": 1,
            "literature_target_records": 0,
            "provenance_held_records": 0,
            "importable_records": 1,
            "received": [path.name for path in paths],
        },
    )
    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/library/import/preview",
            files={"files": ("../../library.csv", b"name,mz\nA,1\n", "text/csv")},
        )
        assert response.status_code == 200, response.text
        assert response.json()["importable_records"] == 1
        assert raw_dir.exists()
        assert list(raw_dir.iterdir()) == []
    finally:
        client.close()


def test_corrupt_saved_results_are_not_reported_as_zero_hits(monkeypatch, tmp_path):
    server, _ = _patch_isolated_upload_state(monkeypatch, tmp_path)
    server._tasks["broken"] = {
        "id": "broken", "status": "completed", "results": [],
    }
    (server.TASK_DIR / "broken_results.json").write_text("{not-json")

    with pytest.raises(HTTPException) as exc_info:
        server._get_task_results("broken")
    assert exc_info.value.status_code == 500
    assert "unreadable" in str(exc_info.value.detail)


def test_v16_health_identity_and_single_endpoint_are_explicit():
    from pathlib import Path

    client = TestClient(app)
    try:
        response = client.get("/api/health")
        assert response.status_code == 200
        assert response.json()["service_generation"] == "V16"
        assert response.json()["default_port"] == 8516
        assert Path(response.json()["project_root"]).resolve() == (
            Path(__file__).resolve().parents[1]
        )
        assert response.json()["capabilities"]["blank_evidence_v2"] is True
        assert response.json()["capabilities"]["blank_filter_default"] == (
            "soft_evidence"
        )
        assert response.json()["storage"]["library_write_enabled"] is False
        system = client.get("/api/system/info")
        assert system.status_code == 200
        assert system.json()["ui_and_api_same_origin"] is True
        assert system.json()["network_engine"]["singletons_are_families"] is False
        page = client.get("/")
        assert page.status_code == 200
        assert page.headers["cache-control"] == "no-store, max-age=0"
        assert "MSsystem · I.D.E.A platform" in page.text
        assert "MSsystem V7" not in page.text
        assert 'id="folderInput"' in page.text
        assert 'id="vendorRawInput"' in page.text
        assert "/static/js/upload.js?v=16.0" in page.text
        script = client.get("/static/js/upload.js")
        assert script.status_code == 200
        assert "refreshRuntimeIdentity" in script.text
        assert "convertVendorFolder" in script.text
        assert "'/tasks/batch/analyze'" in script.text
        assert "refreshBlankUploads" in script.text
        assert "handleIonModeChange" in script.text
        server_source = (
            Path(__file__).resolve().parents[1] / "web/server.py"
        ).read_text(encoding="utf-8")
        assert "DEFAULT_PORT = 8516" in server_source
        run_source = (
            Path(__file__).resolve().parents[1] / "run.sh"
        ).read_text(encoding="utf-8")
        assert "kill -9" not in run_source
        assert "V16 will not kill it" in run_source
    finally:
        client.close()


def test_persisted_task_upload_paths_point_to_current_v16_tree():
    """Migration is local; missing historical uploads are explicit and terminal."""
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    from web.server import _repair_migrated_upload_path

    index = root / "output/tasks_index.json"
    if not index.exists():
        pytest.skip("No persisted local tasks")
    tasks = json.loads(index.read_text(encoding="utf-8"))
    if not tasks:
        pytest.skip("Persisted task index is currently empty")
    for task in tasks.values():
        path = Path(_repair_migrated_upload_path(task["file_path"]))
        assert str(path).startswith(str(root / "output/tmp_uploads"))
        available = path.is_file()
        state = {
            "status": "available" if available else "missing",
            "available": available,
            "managed_upload_path": path.parent == root / "output/tmp_uploads",
            "reanalysis_available": available,
        }
        assert state["managed_upload_path"] is True
        if not path.is_file():
            assert task["status"] in {"completed", "failed", "cancelled"}
            assert state == {
                "status": "missing",
                "available": False,
                "managed_upload_path": True,
                "reanalysis_available": False,
            }
