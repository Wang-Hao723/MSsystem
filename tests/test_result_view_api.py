from fastapi.testclient import TestClient

from web import server


def _task(rows):
    return {
        "id": "task-view",
        "filename": "sample.mzML",
        "status": "completed",
        "progress": 100,
        "upload_time": "2026-08-03T10:30:00Z",
        "total_spectra": len(rows),
        "matched_spectra": sum(bool(row.get("match")) for row in rows),
        "results": rows,
        "presentation_calibration_version": server.PRESENTATION_CALIBRATION_VERSION,
    }


def test_result_list_is_compact_and_detail_is_lazy(monkeypatch):
    rows = [
        {
            "spectrum_id": "scan-1",
            "precursor_mz": 101.1,
            "query_peaks": [[50.0, 100.0], [75.0, 20.0]],
            "display_class": "alkaloid",
            "dl_confidence": 0.7,
            "match": None,
        },
        {
            "spectrum_id": "scan-2",
            "precursor_mz": 202.2,
            "query_peaks": [[80.0, 100.0]],
            "raw_peaks": [[80.0, 100.0]],
            "match": {
                "compound_name": "Example",
                "formula": "C10H12O2",
                "cosine": 0.91,
                "final_score": 0.87,
                "rank_score": 0.87,
                "candidate_adduct_hypothesis": "[M-H]-",
                "candidate_adduct_hypothesis_source": (
                    "structure_formula_library_precursor"
                ),
                "precursor_metadata_audit": {
                    "status": "consistent_via_structure_formula_adduct_hypothesis"
                },
                "mass_metadata_consistent": False,
                "top_candidates": [
                    {"compound_name": "Example", "score_trace": [1] * 1000}
                ],
                "expert_details": {
                    "chemical_reasoning": "large explanation",
                    "v3": {
                        "reasoning": [
                            "分类: phenolic acid",
                            "生源: 莽草酸途径 → 苯丙素",
                        ]
                    },
                    "diagnostic_ion_analysis": {
                        "matched_ions": [{"observed_mz": 151.0}]
                    },
                    "neutral_loss_analysis": {
                        "neutral_losses": [{"mz": 18.0}],
                        "best_pathway": [{"from_mz": 180.0, "to_mz": 162.0}],
                    },
                },
                "evidence_governance": {
                    "policy_id": "evigate",
                    "claim_status": "provisional_needs_independent_injection",
                    "release_blockers": ["needs_replication"],
                    "nodes": [{"payload": "large" * 1000}],
                },
            },
        },
    ]
    monkeypatch.setattr(server, "_tasks", {"task-view": _task(rows)})
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_prepare_results_for_display", lambda value: value)
    monkeypatch.setattr(
        server, "_attach_claim_decisions", lambda value, path: {"status": "absent"}
    )
    monkeypatch.setattr(server, "_refresh_completed_task_presentation", lambda *args: None)

    client = TestClient(server.app)

    metadata = client.get("/api/tasks/task-view")
    assert metadata.status_code == 200
    assert "results" not in metadata.json()

    response = client.get(
        "/api/tasks/task-view/results?projection=summary&limit=10"
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["view_counts"] == {
        "total": 2,
        "matched": 1,
        "unmatched": 1,
        "blank_flagged": 0,
        "blank_excluded": 0,
        "blank_clean_view_excluded": 0,
    }
    assert [row["result_index"] for row in payload["results"]] == [0, 1]
    assert "query_peaks" not in payload["results"][0]
    compact_match = payload["results"][1]["match"]
    assert "top_candidates" not in compact_match
    assert "expert_details" not in compact_match
    assert compact_match["display_evidence_flags"]["has_neutral_loss"] is True
    assert compact_match["display_evidence_flags"]["biosynthetic_pathways"] == [
        "生源: 莽草酸途径 → 苯丙素"
    ]
    assert compact_match["display_evidence_flags"]["neutral_loss_count"] == 1
    assert compact_match["display_evidence_flags"]["fragmentation_step_count"] == 1
    assert compact_match["display_evidence_flags"]["diagnostic_ion_count"] == 1
    assert compact_match["candidate_adduct_hypothesis"] == "[M-H]-"
    assert compact_match["mass_metadata_consistent"] is False
    assert compact_match["display_evidence_flags"][
        "candidate_ion_form_hypothesis"
    ] == "[M-H]-"
    assert compact_match["display_evidence_flags"]["sample_ion_form"] == ""
    assert "nodes" not in compact_match["evidence_governance"]

    unmatched = client.get(
        "/api/tasks/task-view/results?projection=summary&match_state=unmatched"
    ).json()
    assert unmatched["total"] == 1
    assert unmatched["results"][0]["result_index"] == 0

    detail = client.get("/api/tasks/task-view/results/item/1")
    assert detail.status_code == 200
    assert detail.json()["result"]["query_peaks"] == [[80.0, 100.0]]
    assert detail.json()["result"]["match"]["top_candidates"]
