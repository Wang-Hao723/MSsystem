from fastapi.testclient import TestClient

from governance.module_effects import build_module_effect_summary
from web import server


def _example_task():
    return {
        "id": "effect-task",
        "filename": "sample.mzML",
        "status": "completed",
        "progress": 100,
        "total_spectra": 100,
        "matched_spectra": 20,
        "unique_compounds": 8,
        "compound_confidence_summary": {
            "high_supported_unique_compounds": 3,
        },
        "identification_funnel": {"unmatched_with_nonzero_query_peaks": 79},
        "module_availability": {
            "ms1_scan_count": 50,
            "ms1_feature_detection": {
                "requested": True,
                "applicable": True,
                "executed_with_ms2_links": True,
                "lcms_features_with_ms2": 12,
                "linked_ms2_spectra": 30,
            },
            "iimn": {"requested": True, "applicable": True, "executed": True},
            "chimeric_deconvolution": {
                "requested": True,
                "applicable": True,
                "executed": True,
            },
        },
        "deconvolution_impact": {
            "accepted_score_increase_spectra": 2,
            "top1_changed": False,
            "delta": {
                "high_supported_unique_compounds": 1,
                "representative_high_unique_compounds": 1,
                "high_spectra": 2,
            },
        },
        "ion_identity_networks": 2,
        "ion_identity_networked_features": 5,
        "ion_identity_linked_annotations": 1,
        "expected_compounds_summary": {
            "requested": True,
            "status": "ok",
            "detected_feature_count": 10,
            "supported_features": 4,
            "mass_only_features": 6,
            "standard_spectral_library_match_increment": 0,
        },
        "identification_contract_summary": {
            "record_count": 20,
            "production_ranking_changed_records": 0,
            "module_state_counts": {"collision_energy_profile": {"shadow": 20}},
        },
        "claim_contract_summary": {
            "record_count": 20,
            "released_count": 0,
            "calibrated_count": 0,
            "production_ranking_changed_records": 0,
        },
        "evograph_fragmentation_summary": {
            "status": "completed",
            "spectra_processed": 20,
            "truncated_by_task_limit": True,
            "formula_context_available": 8,
            "nodes": 120,
            "selected_edges": 70,
            "heuristically_supported_edges": 45,
            "ambiguous_edges": 25,
            "unassigned_edges": 30,
            "cycle_count": 0,
            "high_promotions": 0,
        },
        "evidiscrim_summary": {
            "requested": True,
            "applicable": True,
            "executed": True,
            "spectra_evaluated": 3,
            "eligible_connectivity_pairs": 18,
            "resolved_counterfactual_pairs": 12,
            "abstained_eligible_pairs": 6,
            "same_connectivity_out_of_scope_pairs": 0,
            "collision_energy_recommendations": 1,
            "instrument_domain_shifted_energy_audits": 0,
            "counterfactual_candidate_order_changed_spectra": 0,
        },
    }


def test_module_effects_distinguish_nonranking_outputs():
    payload = build_module_effect_summary(_example_task())
    effects = {effect["module_id"]: effect for effect in payload["effects"]}

    assert payload["core_ranking_changed"] is False
    assert effects["chimeric_deconvolution"]["metrics"][
        "high_supported_identity_delta"
    ] == 1
    assert effects["chimeric_deconvolution"]["ranking_effect"] == (
        "selected_identity_preserved"
    )
    assert effects["expected_compounds"]["metrics"][
        "ordinary_ms2_match_increment"
    ] == 0
    assert effects["ms1_feature_detection"]["executed"] is True
    assert effects["evograph_fragmentation"]["executed"] is True
    assert effects["evograph_fragmentation"]["ranking_effect"] == "forbidden"
    assert effects["evograph_fragmentation"]["metrics"]["high_promotions"] == 0
    assert "达到任务影子上限" in effects["evograph_fragmentation"]["headline"]
    assert effects["evidiscrim_candidate_contrast"]["executed"] is True
    assert effects["evidiscrim_candidate_contrast"]["metrics"][
        "resolved_counterfactual_pairs"
    ] == 12
    assert effects["evidiscrim_candidate_contrast"]["ranking_effect"] == "forbidden"


def test_module_effect_api(monkeypatch):
    monkeypatch.setattr(
        server, "_ensure_task_presentation_current", lambda task_id: _example_task()
    )

    response = TestClient(server.app).get("/api/tasks/effect-task/module-effects")

    assert response.status_code == 200
    assert response.json()["task_id"] == "effect-task"
    assert len(response.json()["effects"]) >= 8
