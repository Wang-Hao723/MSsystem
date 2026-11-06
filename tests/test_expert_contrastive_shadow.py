from matching.expert_contrastive import evaluate_contrastive_shadow


def test_contrastive_evidence_is_explicitly_shadow_only():
    evidence = evaluate_contrastive_shadow(
        [[151.003, 1000.0], [121.029, 600.0], [107.013, 300.0]],
        301.071,
        "flavonoid",
        query_polarity="NEG",
        candidate_adduct="[M-H]-",
    )

    assert evidence["status"] == "shadow_only_uncalibrated"
    assert evidence["runtime_influence"] is False
    assert evidence["score_contribution"] == 0.0
    assert evidence["rule_rows_used"] > 3_000_000
    assert evidence["limitations"]


def test_contrastive_evidence_rejects_unknown_class():
    evidence = evaluate_contrastive_shadow(
        [[100.0, 100.0]],
        200.0,
        "MassBank_NIST",
        query_polarity="POS",
    )

    assert evidence["status"] == "not_applicable"
    assert evidence["runtime_influence"] is False
    assert evidence["matched_ions"] == []
