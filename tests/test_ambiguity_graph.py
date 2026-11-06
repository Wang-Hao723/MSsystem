from acquisition.ambiguity_graph import build_acquisition_needs


def _row(selected_name="Candidate A"):
    return {
        "spectrum_id": "scan-10",
        "precursor_mz": 301.1234,
        "retention_time": 4.2,
        "ion_mode": "NEG",
        "lcms_feature_id": "feature-7",
        "match": {
            "compound_name": selected_name,
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "formula": "C15H10O7",
            "rank_score": 0.81,
            "confidence": "Medium",
            "top_candidates": [
                {
                    "compound_name": "Candidate B",
                    "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                    "formula": "C15H10O7",
                    "rank_score": 0.78,
                }
            ],
            "chimeric_deconvolution": {"isolation_status": "mixed"},
        },
    }


def test_ambiguity_graph_builds_same_formula_candidate_pair():
    needs = build_acquisition_needs("task-a", [_row()])

    assert len(needs) == 1
    assert len(needs[0].candidates) == 2
    assert len(needs[0].ambiguity_edges) == 1
    assert "same_formula_connectivity_competition" in (
        needs[0].ambiguity_edges[0].reasons
    )
    assert "narrow_isolation_confirmation" in needs[0].missing_evidence
    assert needs[0].candidate_independent_linkage_fields == (
        "polarity", "precursor_mz", "retention_time_min"
    )


def test_entity_id_does_not_depend_on_selected_candidate_name():
    first = build_acquisition_needs("task-a", [_row("Candidate A")])[0]
    second = build_acquisition_needs("task-a", [_row("Renamed Candidate")])[0]

    assert first.entity_id == second.entity_id

