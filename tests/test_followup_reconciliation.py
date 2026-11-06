from acquisition.reconciliation import link_followup_observations


def _plan():
    return {
        "plan_id": "plan-1",
        "task_id": "baseline",
        "need": {
            "precursor_mz": 300.0,
            "retention_time_min": 4.0,
            "polarity": "NEG",
            "candidates": [
                {"identity_key": "ik::A"},
                {"identity_key": "ik::B"},
            ],
        },
    }


def test_followup_linkage_uses_observation_coordinates_before_candidates():
    followup = [{
        "precursor_mz": 300.001,
        "retention_time": 4.05,
        "ion_mode": "NEG",
        "match": {"compound_name": "A", "inchikey": "A"},
    }]

    delta = link_followup_observations(
        _plan(), "followup", followup, mz_tolerance_ppm=10, rt_tolerance_min=0.1
    )

    assert delta.linked_result_indices == (0,)
    assert delta.candidate_fields_used_for_linkage == ()
    assert delta.linkage_fields_used == (
        "polarity", "precursor_mz", "retention_time_min"
    )
    assert delta.production_ranking_changed is False
    assert delta.automatic_claim_release is False
