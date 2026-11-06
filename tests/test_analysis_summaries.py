from web import server


def _row(scan_id, confidence, rank_score):
    return {
        "spectrum_id": scan_id,
        "precursor_mz": 301.0073,
        "precursor_feature_id": f"PF{scan_id}",
        "precursor_isolation_evidence": {
            "status": "referenced_ms1_not_in_file",
            "purity": None,
        },
        "match": {
            "compound_name": "Example compound",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "smiles": "CCO",
            "formula": "C2H6O",
            "confidence": confidence,
            "rank_score": rank_score,
            "final_score": rank_score,
            "cosine": rank_score,
            "chimeric_deconvolution": {
                "status": "isolation_window_unavailable",
                "isolation_status": "referenced_ms1_not_in_file",
                "runtime_influence": False,
            },
        },
    }


def test_high_supported_and_representative_high_are_separate_counts():
    rows = [
        _row("000001", "High", 0.72),
        _row("000002", "Medium", 0.91),
    ]

    summary = server._compound_confidence_summary(rows)

    assert summary["unique_library_candidates"] == 1
    assert summary["high_reportable_scope_unique_compounds"] == 0
    assert summary["representative_high_reportable_scope_unique_compounds"] == 0
    assert summary["full_identity_high_supported_unique_compounds"] == 1
    assert summary["representative_full_identity_high_unique_compounds"] == 0
    assert summary["full_identity_high_spectra"] == 1
    assert summary["high_supported_unique_compounds"] == 1
    assert summary["representative_high_unique_compounds"] == 0
    assert summary["supported_but_representative_not_high"] == 1
    assert summary["ranking_changed"] is False


def test_reportable_high_scope_includes_unresolved_isomer_group_without_full_identity_high():
    row = _row("000003", "Medium", 0.91)
    row["match"].update({
        "identification_scope": "isomer_group",
        "isomer_group_confidence": "High",
        "two_dimensional_connectivity_confidence": "Medium",
    })

    summary = server._compound_confidence_summary([row])

    assert summary["high_reportable_scope_unique_compounds"] == 1
    assert summary["representative_isomer_group_high_unique_compounds"] == 1
    assert summary["full_identity_high_supported_unique_compounds"] == 0
    assert summary["representative_full_identity_high_unique_compounds"] == 0
    assert summary["legacy_unqualified_high_fields_scope"] == "full_structure_identity"


def test_deconvolution_impact_uses_frozen_raw_scores_without_reranking():
    row = _row("000001", "High", 0.72)
    row["match"].update({
        "cosine": 0.93,
        "raw_cosine": 0.86,
        "final_score": 0.62,
        "raw_final_score": 0.48,
        "rank_score": 0.48,
        "raw_rank_score": 0.48,
        "reverse_match": 0.90,
        "raw_reverse_match": 0.90,
        "classification": "flavonoid",
        "class_agreement": 1.0,
        "source_tier": "experimental",
        "chimeric_deconvolution": {
            "status": "accepted_score_applied",
            "accepted": True,
            "runtime_influence": True,
            "score_contribution": 0.14,
            "isolation_status": "mixed",
        },
    })

    impact = server._deconvolution_impact_summary([row])

    assert impact["baseline"]["high_supported_unique_compounds"] == 0
    assert impact["current"]["high_supported_unique_compounds"] == 1
    assert impact["delta"]["high_supported_unique_compounds"] == 1
    assert impact["accepted_score_increase_spectra"] == 1
    assert impact["top1_changed"] is False
    assert impact["candidate_margin_changed"] is False
    assert impact["raw_rank_score_changed"] is False
    assert impact["invariant_violations"] == {
        "top1_identity": 0,
        "candidate_margin": 0,
        "raw_rank_score": 0,
    }


def test_deconvolution_impact_detects_saved_invariant_violations():
    row = _row("000001", "High", 0.72)
    row["match"].update({
        "cosine": 0.93,
        "raw_cosine": 0.86,
        "final_score": 0.62,
        "raw_final_score": 0.48,
        "rank_score": 0.60,
        "raw_rank_score": 0.48,
        "candidate_margin": 0.15,
        "classification": "flavonoid",
        "class_agreement": 1.0,
        "source_tier": "experimental",
        "chimeric_deconvolution": {
            "status": "accepted_score_applied",
            "accepted": True,
            "runtime_influence": True,
            "score_contribution": 0.14,
            "isolation_status": "mixed",
            "selected_identity_before_deconvolution": "ik::DIFFERENT",
            "top1_identity_preserved": False,
            "deconvolved_candidate_margin_used": True,
        },
    })

    impact = server._deconvolution_impact_summary([row])

    assert impact["top1_changed"] is True
    assert impact["candidate_margin_changed"] is True
    assert impact["raw_rank_score_changed"] is True
    assert impact["invariant_violations"] == {
        "top1_identity": 1,
        "candidate_margin": 1,
        "raw_rank_score": 1,
    }


def test_identification_funnel_separates_scans_events_and_identities():
    high = _row("000001", "High", 0.95)
    high["precursor_feature_id"] = "PF000001"
    low = _row("000002", "Unknown", 0.65)
    low["precursor_feature_id"] = "PF000002"
    low["match"]["inchikey"] = "CCCCCCCCCCCCCC-DDDDDDDDDD-E"
    low["match"]["cosine"] = 0.65
    low["match"]["final_score"] = 0.40
    unmatched = {
        "spectrum_id": "000003",
        "precursor_feature_id": "PF000003",
        "query_peaks": [[100.0, 1.0], [120.0, 0.5], [140.0, 0.2]],
    }

    funnel = server._identification_funnel_summary([high, low, unmatched])

    assert funnel["ms2_scans"] == 3
    assert funnel["ms2_precursor_events_approximate"] == 3
    assert funnel["matched_ms2_scans"] == 2
    assert funnel["unique_library_candidates"] == 2
    assert funnel["unmatched_ms2_scans"] == 1
    assert funnel["unmatched_with_nonzero_query_peaks"] == 1
    assert funnel["non_high_primary_blockers"]["cosine_below_0_70"] == 1
    assert funnel["units_are_not_interchangeable"] is True


def test_identification_funnel_excludes_blank_removed_rows_from_search_units():
    high = _row("000001", "High", 0.95)
    excluded = _row("000002", "High", 0.96)
    excluded["blank_excluded"] = True

    funnel = server._identification_funnel_summary([high, excluded])

    assert funnel["ms2_scans"] == 2
    assert funnel["blank_excluded_ms2_scans"] == 1
    assert funnel["ms2_precursor_events_approximate"] == 1
    assert funnel["matched_ms2_scans"] == 1
    assert funnel["matched_ms2_precursor_events"] == 1
    assert funnel["unique_library_candidates"] == 1


def test_ms2_only_input_marks_deconvolution_and_iimn_unavailable():
    rows = [_row("000001", "High", 0.72)]
    legacy_iimn_summary = {
        "status": "ok",
        "feature_count": 1,
        "network_count": 0,
        "usable_trace_count": 0,
        "trace_extraction": {
            "status": "no_ms1_scans",
            "ms1_scan_count": 0,
            "usable_trace_count": 0,
        },
    }

    availability = server._analysis_module_availability(
        {
            "enable_chimeric_deconvolution": True,
            "enable_ion_identity_networking": True,
        },
        legacy_iimn_summary,
        rows,
    )

    assert availability["continuous_ms1_available"] is False
    assert availability["ms1_scan_count"] == 0
    assert availability["iimn"]["status"] == "unavailable_no_ms1_scans"
    assert availability["iimn"]["executed"] is False
    assert availability["chimeric_deconvolution"]["status"] == (
        "unavailable_no_ms1_isolation_context"
    )
    assert availability["chimeric_deconvolution"]["executed"] is False
    assert [item["code"] for item in availability["warnings"]] == [
        "continuous_ms1_missing"
    ]
    assert "仅填写隔离窗宽度不能替代缺失的 MS1" in (
        availability["warnings"][0]["action"]
    )


def test_legacy_iimn_ok_status_is_repaired_from_trace_extraction():
    effective = server._effective_ion_identity_summary({
        "status": "ok",
        "usable_trace_count": 0,
        "trace_extraction": {
            "status": "no_ms1_scans",
            "ms1_scan_count": 0,
        },
    })

    assert effective["status"] == "no_ms1_scans"
