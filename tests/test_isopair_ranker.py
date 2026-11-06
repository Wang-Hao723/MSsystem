"""Runtime safety and ranking tests for IsoPair-MS."""

import numpy as np

import matching.isopair_ranker as isopair
from matching.library_quality import formula_exact_mass
from utils.adduct import ion_mz


IK_A1 = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_A2 = "AAAAAAAAAAAAAA-CCCCCCCCCC-D"
IK_B = "DDDDDDDDDDDDDD-EEEEEEEEEE-F"


class _FakeModel:
    def predict(self, matrix):
        matrix = np.asarray(matrix)
        return np.where(matrix[:, 0] >= 60.0, 0.99, 0.10)


def _artifact():
    return _FakeModel(), {
        "model_name": "IsoPair-MS test artifact",
        "fusion_alpha_baseline": 0.275,
        "selective_high_thresholds": {
            "score_threshold": 0.75,
            "margin_threshold": 0.20,
        },
    }


def _candidate(
    inchikey,
    cosine,
    peak_mz,
    *,
    formula="C10H12O2",
    adduct="[M+H]+",
):
    return {
        "compound_name": inchikey[:4],
        "inchikey": inchikey,
        "formula": formula,
        "adduct": adduct,
        "source_tier": "experimental",
        "cosine": cosine,
        "final_score": cosine,
        "rank_score": cosine,
        "reverse_match": 0.9,
        "mass_ppm": 2.0,
        "lib_precursor_mz": 101.0,
        "lib_peaks": [(peak_mz, 100.0), (30.0, 20.0), (40.0, 10.0)],
    }


def _install_fake_artifact(monkeypatch):
    monkeypatch.setattr(isopair, "_load_artifact", lambda *_args: _artifact())
    # Most tests exercise the legacy ranking mechanics deliberately.  The
    # separate fail-closed test below verifies the production default.
    monkeypatch.setattr(isopair, "ISOPAIR_RANK_INFLUENCE_ENABLED", True)
    monkeypatch.setattr(
        isopair,
        "spectral_pair_features",
        lambda _query, reference, **_kwargs: np.asarray(
            [float(reference[0][0])], dtype=np.float32
        ),
    )


def _run(best, rows):
    return isopair.apply_isopair_ranker(
        best,
        rows,
        [(20.0, 100.0), (40.0, 50.0), (60.0, 20.0)],
        query_precursor_mz=101.0,
        query_adduct="[M+H]+",
    )


def test_near_tie_can_change_only_the_same_formula_connectivity(monkeypatch):
    _install_fake_artifact(monkeypatch)
    first = _candidate(IK_A1, 0.90, 50.0)
    second = _candidate(IK_B, 0.86, 60.0)

    selected, evidence = _run(first, [first, second])

    assert selected is second
    assert evidence["ranking_changed"] is True
    assert evidence["selected_connectivity"] == IK_B[:14]
    assert evidence["identity_scope"] == "2d_connectivity_only"
    assert evidence["stereochemistry_claimed"] is False
    assert evidence["high_precision_support"] is False
    assert evidence["high_promotion_allowed"] is False


def test_production_default_keeps_explicit_adduct_ranking_in_shadow(monkeypatch):
    _install_fake_artifact(monkeypatch)
    monkeypatch.setattr(isopair, "ISOPAIR_RANK_INFLUENCE_ENABLED", False)
    first = _candidate(IK_A1, 0.90, 50.0)
    second = _candidate(IK_B, 0.86, 60.0)

    selected, evidence = _run(first, [first, second])

    assert selected is first
    assert evidence["status"] == "explicit_adduct_shadow_near_tie"
    assert evidence["ranking_changed"] is False
    assert evidence["counterfactual_ranking_changed"] is True
    assert evidence["runtime_influence"] is False
    assert evidence["production_rank_influence_allowed"] is False
    assert evidence["production_abstention_reason"] == (
        "independent_family_confirmation_gain_not_statistically_qualified"
    )


def test_clear_baseline_winner_is_never_reranked(monkeypatch):
    _install_fake_artifact(monkeypatch)
    first = _candidate(IK_A1, 0.95, 50.0)
    second = _candidate(IK_B, 0.70, 60.0)

    selected, evidence = _run(first, [first, second])

    assert selected is first
    assert evidence["status"] == "clear_baseline_winner_no_rerank"
    assert evidence["runtime_influence"] is False
    assert evidence["ranking_changed"] is False


def test_selective_score_supports_but_never_creates_high(monkeypatch):
    _install_fake_artifact(monkeypatch)
    first = _candidate(IK_A1, 0.95, 60.0)
    second = _candidate(IK_B, 0.75, 50.0)

    selected, evidence = _run(first, [first, second])

    assert selected is first
    assert evidence["high_precision_support"] is True
    assert first["isopair_high_precision_support"] is True
    assert evidence["high_promotion_allowed"] is False


def test_mass_inconsistent_missing_query_adduct_fails_closed(monkeypatch):
    _install_fake_artifact(monkeypatch)
    first = _candidate(IK_A1, 0.95, 60.0)
    second = _candidate(IK_B, 0.75, 50.0)

    selected, evidence = isopair.apply_isopair_ranker(
        first,
        [first, second],
        [(20.0, 1.0), (40.0, 0.5), (60.0, 0.2)],
        query_precursor_mz=101.0,
        query_adduct="",
        query_polarity="POS",
    )

    assert selected is first
    assert evidence["status"] == "not_applicable"
    assert evidence["reason"] == "no_formula_mass_consistent_adduct"


def test_formula_conditioned_adduct_is_scored_in_shadow_without_reranking(
    monkeypatch,
):
    _install_fake_artifact(monkeypatch)
    first = _candidate(IK_A1, 0.90, 50.0)
    second = _candidate(IK_B, 0.86, 60.0)
    precursor = ion_mz(formula_exact_mass("C10H12O2"), "[M+H]+")

    selected, evidence = isopair.apply_isopair_ranker(
        first,
        [first, second],
        [(20.0, 1.0), (40.0, 0.5), (60.0, 0.2)],
        query_precursor_mz=precursor,
        query_adduct="",
        query_polarity="POS",
    )

    assert selected is first
    assert evidence["evaluated"] is True
    assert evidence["status"] == "formula_conditioned_shadow_near_tie"
    assert evidence["ranking_changed"] is False
    assert evidence["counterfactual_ranking_changed"] is True
    assert evidence["counterfactual_selected_connectivity"] == IK_B[:14]
    assert evidence["production_rank_influence_allowed"] is False
    assert evidence["runtime_influence"] is False
    assert evidence["adduct_context"]["formula_conditioned"] is True
    assert evidence["high_precision_gate_applicable"] is False
    assert evidence["selective_threshold_pass"] is False
    assert evidence["counterfactual_selective_threshold_pass"] is True
    assert evidence["high_precision_support"] is False
    assert first["isopair_high_precision_support"] is False
    assert first["isopair_selected_connectivity"] is True
    assert second["isopair_counterfactual_selected_connectivity"] is True


def test_same_connectivity_full_key_variants_are_not_ranked(monkeypatch):
    _install_fake_artifact(monkeypatch)
    first = _candidate(IK_A1, 0.95, 60.0)
    stereo_variant = _candidate(IK_A2, 0.94, 50.0)
    other_connectivity = _candidate(IK_B, 0.70, 50.0)

    selected, evidence = _run(
        first, [first, stereo_variant, other_connectivity]
    )

    assert selected is first
    assert evidence["connectivity_candidate_count"] == 2
    assert first["isopair_selected_connectivity"] is True
    assert stereo_variant["isopair_selected_connectivity"] is True
    assert evidence["stereochemistry_claimed"] is False


def test_candidate_count_above_training_domain_fails_closed(monkeypatch):
    _install_fake_artifact(monkeypatch)
    rows = []
    for index in range(13):
        letter = chr(ord("A") + index)
        key = f"{letter * 14}-BBBBBBBBBB-C"
        rows.append(_candidate(key, 0.95 - index * 0.01, 60.0))

    selected, evidence = _run(rows[0], rows)

    assert selected is rows[0]
    assert evidence["status"] == "not_applicable"
    assert evidence["reason"] == (
        "connectivity_candidate_count_outside_training_domain"
    )
