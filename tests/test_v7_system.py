import json
from pathlib import Path

from fastapi.testclient import TestClient


ROOT = Path(__file__).resolve().parents[1]


def test_public_negative_source_registry_is_exposed():
    from web.server import app

    with TestClient(app) as client:
        response = client.get("/api/library/source-registry")
        legacy_alias = client.get("/api/library/public-sources")
        system = client.get("/api/system/info")
    assert response.status_code == 200
    assert legacy_alias.status_code == 200
    payload = response.json()
    snapshot = payload["database_snapshot"]
    assert snapshot["positive_only_compounds"] > 0
    assert snapshot["negative_only_compounds"] > 0
    assert (
        snapshot["positive_only_compounds"]
        + snapshot["negative_only_compounds"]
        + snapshot["both_polarities"]
        <= snapshot["mapped_compounds"]
    )
    assert payload["coverage_summary"]["positive_compound_coverage_pct"] > 70
    assert payload["coverage_summary"]["negative_compound_coverage_pct"] > 40
    assert payload["coverage_summary"]["dual_polarity_coverage_pct"] > 20
    source_ids = {source["id"] for source in payload["sources"]}
    assert "riken_plasma_authentic_negative" in source_ids
    assert "gnps_nih_natural_products_round2_negative" in source_ids
    assert payload["ingestion_policy"]["experimental_only_default"] is True
    assert system.status_code == 200
    assert system.json()["storage"]["free_gib"] > 0


def test_agilent_input_audit_requires_ms_levels_and_paired_hashes(monkeypatch):
    from web import server

    task = {
        "id": "agilent-v7",
        "filename": "sample.mzML",
        "source_raw_filename": "sample.d",
        "vendor_conversion": {
            "status": "completed",
            "source": {
                "acquisition": {
                    "instrument": "Q-TOF",
                    "acquisition_mode": "AutoMS2",
                    "polarities": ["NEG"],
                    "ms1_expected": True,
                    "ms2_expected": True,
                    "isolation_width_da": 4.0,
                    "collision_energies": [30.0, 40.0],
                }
            },
            "output": {
                "status": "valid",
                "sha256": "abc",
                "actual_spectrum_count": 12,
                "mass_spectrum_count": 12,
                "auxiliary_spectrum_count": 0,
                "ms1_count": 5,
                "nonempty_ms1_count": 5,
                "ms2_count": 7,
                "nonempty_ms2_count": 7,
                "precursor_count": 7,
                "polarities": {"NEG": 12, "POS": 0},
                "source_checksum_matches": {
                    "MSScan.bin": True,
                    "MSPeak.bin": True,
                    "AcqMethod.xml": True,
                },
                "warnings": [],
                "errors": [],
            },
            "converter": {"backend": "test"},
        },
    }
    monkeypatch.setattr(server, "_load_tasks", lambda: None)
    monkeypatch.setattr(server, "_tasks", {"agilent-v7": task})

    with TestClient(server.app) as client:
        verified = client.get("/api/tasks/agilent-v7/input-audit")
        assert verified.status_code == 200
        body = verified.json()
        assert body["verdict"] == "verified"
        assert body["preserves_expected_ms_levels"] is True
        assert body["paired_raw_provenance"] is True

        server._tasks["agilent-v7"]["vendor_conversion"]["output"]["ms1_count"] = 0
        failed = client.get("/api/tasks/agilent-v7/input-audit")
        assert failed.json()["verdict"] == "failed_closed"


def test_v7_frontend_uses_polarity_neutral_workbench_and_library_center():
    html = (ROOT / "web" / "static" / "index.html").read_text(encoding="utf-8")
    assert 'data-page="dashboard"' in html
    assert 'id="page-dashboard"' in html
    assert 'data-page="library"' in html
    assert 'id="page-library"' in html
    assert 'data-page="data-sources"' not in html
    assert "负离子数据源" not in html
    assert "/static/js/dashboard.js?v=15.1" in html
    assert "/static/js/results.js?v=16.3-reliability-20260920" in html
    assert "/static/js/library-center.js?v=7.2" in html
    assert "/static/css/workspace-v7.css?v=16.0" in html
    assert (
        'LC–MS/MS 数据处理与<span class="dashboard-title-term">化合物鉴定</span>'
        in html
    )
    assert "POS 与 NEG 使用同一分析流程和判定标准" in html
    assert "从原始数据到可审计的分子证据" not in html
    assert "Evidence Workbench" not in html
    registry = json.loads(
        (ROOT / "config" / "public_source_registry.json").read_text(
            encoding="utf-8"
        )
    )
    assert registry["audited_at"] == "2026-07-29"


def test_run_script_has_safe_port_ownership_contract():
    source = (ROOT / "run.sh").read_text(encoding="utf-8")
    assert 'PORT="${MSSYSTEM_PORT:-8516}"' in source
    assert "lsof -tiTCP" in source
    assert "pid_is_ours" in source
    assert "kill -9" not in source
    assert "pip install" not in source
