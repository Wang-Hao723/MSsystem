from eval.evaluate_masscube_mix33_identification import mass_pool, representative, top_candidates


def test_all_arms_share_the_existing_mass_cluster():
    rows = [{"precursor_mz": mz} for mz in (300., 300.005, 300.02, 300.6)]
    assert mass_pool(rows, 300.) == rows[:2]
    assert mass_pool(rows, None) == []


def test_quality_and_isf_selection_use_measurements_without_compound_identity():
    rows = [
        {"spectrum_id": "1", "feature_id": "a", "precursor_mz": 300., "tic": 1000, "truth": True},
        {"spectrum_id": "2", "feature_id": "b", "precursor_mz": 300., "tic": 500, "truth": False},
    ]
    features = {
        "a": {"quality": {"status": "ok", "gaussian_similarity": .7, "noise_score": .3, "asymmetry_factor": 2.}, "candidate_isf_parent_ids": ["parent"]},
        "b": {"quality": {"status": "ok", "gaussian_similarity": .99, "noise_score": 0., "asymmetry_factor": 1.}, "candidate_isf_parent_ids": []},
    }
    assert representative(rows, features, 300., "baseline")["spectrum_id"] == "1"
    for arm in ("isf", "quality", "combined"):
        assert representative(rows, features, 300., arm)["spectrum_id"] == "2"
    rows[0]["truth"], rows[1]["truth"] = False, True
    assert representative(rows, features, 300., "combined")["spectrum_id"] == "2"


def test_missing_quality_or_all_fragment_flags_preserve_output():
    rows = [{"spectrum_id": "1", "feature_id": "a", "precursor_mz": 300., "tic": 1000}]
    assert representative(rows, {}, 300., "quality") == rows[0]
    assert representative(rows, {"a": {"candidate_isf_parent_ids": ["p"]}}, 300., "combined") == rows[0]
    assert representative([], {}, 300., "combined") is None


def test_top5_preserves_top1_and_deduplicates_connectivity():
    first = {"inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "compound_name": "a"}
    other_stereo = {"inchikey": "AAAAAAAAAAAAAA-CCCCCCCCCC-D", "compound_name": "other stereochemistry"}
    remaining = [{"inchikey": chr(66+i)*14 + "-BBBBBBBBBB-C"} for i in range(5)]
    actual = top_candidates(dict(first, top_candidates=[other_stereo] + remaining))
    assert len(actual) == 5
    assert actual[0]["inchikey"] == first["inchikey"]
    assert other_stereo["inchikey"] not in [item["inchikey"] for item in actual]
