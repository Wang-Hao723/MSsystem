from __future__ import annotations

import json
import re
from dataclasses import asdict
from pathlib import Path

from fastapi.testclient import TestClient

from matching.authentic_anchor import AuthenticAnchorPolicy, finalize_anchor_library
from web import server


ROOT = Path(__file__).resolve().parents[1]
OVERVIEW = (ROOT / "web/static/js/sample-overview.js").read_text(encoding="utf-8")
REPORT = (ROOT / "web/static/js/report-generator.js").read_text(encoding="utf-8")
EVIDENCE_REVIEW = (ROOT / "web/static/js/evidence-review.js").read_text(
    encoding="utf-8"
)
INDEX = (ROOT / "web/static/index.html").read_text(encoding="utf-8")


def _library() -> dict:
    policy = AuthenticAnchorPolicy().validated()
    return finalize_anchor_library(
        {
            "library_id": "evitwin-api-test",
            "runtime_mode": "shadow_only",
            "construction_policy": {
                "explicit_task_id_only": True,
                "candidate_generation_allowed": False,
                "authentic_library_shadow_search_allowed": True,
            },
            "release_policy": asdict(policy),
            "entries": [
                {
                    "reference_id": "RIGHT:NEG",
                    "compound_code": "RIGHT",
                    "compound_name": "Right authentic candidate",
                    "formula": "C10H10O4",
                    "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    "neutral_mass": 194.05791,
                    "polarity": "NEG",
                    "instrument_family": "agilent-qtof-6546",
                    "method_id": "c18-test-method",
                    "authentic_standard": True,
                    "source_task_id": "independent-standard-task",
                    "reference_qc": {"selected_adduct": "[M-H]-"},
                    "spectra": [
                        {
                            "collision_energy_ev": 10.0,
                            "adduct": "[M-H]-",
                            "retention_time_min": 3.0,
                            "peaks": [[80.0, 100.0], [100.0, 60.0], [120.0, 40.0]],
                        },
                        {
                            "collision_energy_ev": 20.0,
                            "adduct": "[M-H]-",
                            "retention_time_min": 3.0,
                            "peaks": [[80.0, 80.0], [100.0, 100.0], [140.0, 30.0]],
                        },
                    ],
                }
            ],
        }
    )


def test_evitwin_api_materializes_active_ranker_from_v15_without_touching_primary(
    tmp_path, monkeypatch
):
    task_id = "evitwin-api"
    library = _library()
    library_path = tmp_path / "authentic_library.json"
    library_path.write_text(json.dumps(library), encoding="utf-8")
    registry_path = tmp_path / "authentic_anchor_registry.json"
    registry_path.write_text(
        json.dumps(
            {
                "schema_version": "authentic-anchor-registry/1.0",
                "active_library": library_path.name,
                "active_library_content_sha256": library["content_sha256"],
                "state": "active_authentic_reference_display_ranker",
                "instrument_family": "agilent-qtof-6546",
            }
        ),
        encoding="utf-8",
    )
    rows = [
        {
            "collision_energy": 10.0,
            "query_adduct": "[M-H]-",
            "query_adduct_reliable": True,
            "rt_method_id": "c18-test-method",
            "query_peaks": [[80.0, 100.0], [100.0, 60.0], [120.0, 40.0]],
        },
        {
            "collision_energy": 20.0,
            "query_adduct": "[M-H]-",
            "query_adduct_reliable": True,
            "rt_method_id": "c18-test-method",
            "query_peaks": [[80.0, 80.0], [100.0, 100.0], [140.0, 30.0]],
        },
    ]
    v15 = {
        "metadata": {
            "module_id": "v15_conditioned_multilibrary_flagship",
            "algorithm_id": "mssystem-v15-test",
        },
        "entities": [
            {
                "entity_id": "entity-1",
                "representative_precursor_mz": 193.0506335,
                "representative_retention_time": 4.2,
                "member_result_indices": [0, 1],
                "identity_resolution": {
                    "status": "unresolved_same_formula_isomer_group",
                    "reportable_formula": "C10H10O4",
                },
                "candidates": [
                    {
                        "identity_key": "std2d:CCCCCCCCCCCCCC",
                        "vetoed": False,
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                        },
                    },
                    {
                        "identity_key": "std2d:AAAAAAAAAAAAAA",
                        "vetoed": False,
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                        },
                    },
                ],
                "winner": {
                    "entity_score": 0.8,
                    "candidate": {
                        "compound_name": "Wrong public winner",
                        "formula": "C10H10O4",
                        "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                    },
                },
            }
        ]
    }
    (tmp_path / f"{task_id}_results.json").write_text(
        json.dumps(rows), encoding="utf-8"
    )
    v15_path = server._flagship_artifact_path(tmp_path, task_id)
    original_v15 = json.dumps(v15, sort_keys=True)
    v15_path.write_text(json.dumps(v15), encoding="utf-8")
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(server, "EVITWIN_REGISTRY_FILE", registry_path)
    # This test isolates the read-only EviTwin projection.  Production now
    # rejects stale V15 identity artifacts before building any downstream
    # view; the synthetic fixture intentionally omits the full production
    # policy and therefore declares its primary sidecar current explicitly.
    monkeypatch.setattr(
        server,
        "_flagship_artifact_is_current",
        lambda path: path == v15_path,
    )
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "filename": "family_NEG.mzML",
                "ion_mode": "NEG",
                "params": {"rt_method_id": "c18-test-method"},
                "vendor_conversion": {
                    "source": {
                        "acquisition": {
                            "instrument_models": ["Agilent Q-TOF"]
                        }
                    }
                },
            }
        },
    )

    response = TestClient(server.app).get(f"/api/tasks/{task_id}/evitwin-ms")

    assert response.status_code == 200
    payload = response.json()
    assert payload["candidate_display_order_active"] is True
    assert payload["primary_v15_artifact_unchanged"] is True
    assert payload["formal_level1_release_increment"] == 0
    assert payload["compounds"][0]["evitwin_winner"]["compound_code"] == "RIGHT"
    assert payload["metadata"]["source_primary_is_v15"] is True
    assert json.dumps(json.loads(v15_path.read_text()), sort_keys=True) == original_v15


def test_evitwin_is_visible_first_and_is_included_in_the_overall_report():
    assert "/evitwin-ms" in OVERVIEW
    assert "EviTwin-MS 标品谱联合排序" in OVERVIEW
    assert "真实标品谱优先候选" in REPORT
    assert "sample-report-authentic-ranking" in REPORT
    assert "样品 RT" in REPORT
    assert "confidence-guide.js?v=14.0.4" in INDEX
    assert re.search(r'<script\s+src="/static/js/sample-overview\.js\?v=[^"\s]+"', INDEX)
    assert re.search(r'<script\s+src="/static/js/report-generator\.js\?v=[^"\s]+"', INDEX)
    assert "组内竞争" in OVERVIEW
    assert "组内竞争" in REPORT
    assert "evidence-review.js?v=16.0.2" in INDEX
    assert "Top40 扩展候选支持" in EVIDENCE_REVIEW
    assert "逐 CE 查询谱质量审计" in EVIDENCE_REVIEW
    assert "项目 target–decoy" in EVIDENCE_REVIEW
    assert "q-value ${calibration.q_value == null ? '未生成'" in EVIDENCE_REVIEW
