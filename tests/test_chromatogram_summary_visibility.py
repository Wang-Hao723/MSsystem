"""TIC/BPC summary remains visible and preserves acquisition scope."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OVERVIEW = (ROOT / "web/static/js/sample-overview.js").read_text(encoding="utf-8")
REPORT = (ROOT / "web/static/js/report-generator.js").read_text(encoding="utf-8")
I18N = (ROOT / "web/static/js/i18n-results-catalog.js").read_text(encoding="utf-8")


def test_chromatogram_summary_is_outside_collapsed_details_for_saved_results():
    render = OVERVIEW[OVERVIEW.index("function renderOverview(bundle)"):
                      OVERVIEW.index("function matchesAnalysisFilter(group)")]
    assignment = render[render.index("area.innerHTML = `"):]
    summary = assignment.index("${chromatogramHtml}")
    compounds = assignment.index('<section class="card sample-compound-card">')
    disclosure = assignment.index('<details class="vr-disclosure">')

    assert summary < compounds < disclosure
    assert render.count("${chromatogramHtml}") == 1


def test_summary_labels_ms2_only_trace_without_calling_it_ms1():
    assert "当前保存结果为 MS2-only 采集事件" in REPORT
    assert "TIC 仅按 MS2 谱图强度近似展示，不属于 MS1 色谱" in REPORT
    assert "系统不会用 MS2 基峰冒充真实 BPC" in REPORT
    assert "const hasMs1Trace = ms1TraceCount >= 2" in REPORT
    assert "The TIC is an MS2-spectrum intensity approximation, not an MS1 chromatogram" in I18N


def test_tic_api_returns_saved_ms1_and_bpc_payload_without_reanalysis(
    tmp_path, monkeypatch,
):
    from web import server

    task_id = "registered-tic-payload"
    payload = {
        "chromatogram_schema_version": 3,
        "ms1": [
            {"rt": 1.0, "tic": 100.0, "bpc": 80.0, "base_peak_mz": 101.1},
            {"rt": 1.2, "tic": 250.0, "bpc": 190.0, "base_peak_mz": 202.2},
        ],
        "ms2": [{"rt": 1.1, "tic": 30.0, "scan_id": 7, "precursor_mz": 202.2}],
        "rt_unit": "minute",
        "tic_source": "mzml_ms1_survey_scans",
        "bpc_available": True,
        "bpc_source": "mzml_ms1_survey_scans",
        "bpc_backfill_attempted": True,
    }
    (tmp_path / f"{task_id}_tic_overview.json").write_text(
        json.dumps(payload), encoding="utf-8",
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setitem(server._tasks, task_id, {
        "id": task_id, "filename": "saved-result.mzML", "status": "completed",
    })

    result = asyncio.run(server.get_task_tic_overview(task_id))

    assert result == payload
    assert result["ms1"][1]["bpc"] == 190.0
    assert result["tic_source"] == "mzml_ms1_survey_scans"


def test_tic_api_marks_legacy_mgf_fallback_as_ms2_only(tmp_path, monkeypatch):
    from web import server

    task_id = "registered-ms2-only"
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setitem(server._tasks, task_id, {
        "id": task_id, "filename": "saved-result.mgf", "status": "completed",
    })
    monkeypatch.setattr(server, "_get_task_results", lambda _task_id: [{
        "retention_time": 60.0,
        "tic": 1200.0,
        "bpc": 700.0,
        "base_peak_mz": 149.02,
        "spectrum_id": 4,
        "precursor_mz": 249.15,
    }])

    result = asyncio.run(server.get_task_tic_overview(task_id))

    assert result["ms1"] == []
    assert result["ms1_spectra_available"] is False
    assert result["tic_source"] == "ms2_spectrum_metadata_approximation"
    assert result["bpc_available"] is False
    assert result["bpc_scope"] == "unavailable_without_continuous_ms1"
