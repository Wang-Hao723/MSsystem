"""Regression tests for High-confidence candidate ordering and ambiguity."""

import asyncio

import matching.honest_calibration as honest_calibration
from matching.honest_calibration import annotate
from matching.search import _select_best_result, _top_candidates_with_best
from eval.eval_reproducibility import _candidate_list


IK_BEST = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_ALT = "CCCCCCCCCCCCCC-DDDDDDDDDD-E"


def _candidate(
    name,
    inchikey,
    *,
    formula="C10H12O2",
    cosine=0.95,
    final=0.70,
    rank=None,
    classification="flavonoid",
):
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "formula": formula,
        "smiles": "CCOc1ccccc1",
        "has_structure": True,
        "cosine": cosine,
        "final_score": final,
        "rank_score": final if rank is None else rank,
        "reverse_match": 0.90,
        "mass_ppm": 5.0,
        "classification": classification,
        "source_tier": "experimental",
    }


def test_weak_same_formula_alternative_is_not_forced_isomer_ambiguity():
    best = _candidate("best", IK_BEST, cosine=0.95, final=0.70)
    weak = _candidate("weak isomer", IK_ALT, cosine=0.60, final=0.45)
    best["top_candidates"] = [best.copy(), weak]

    result = annotate(best)

    assert result["same_formula_alternatives_n"] == 1
    assert result["noncompetitive_same_formula_n"] == 1
    assert result["isomer_ambiguous"] is False
    assert result["compound_confidence"] == "High"


def test_close_same_formula_alternative_remains_ambiguous():
    best = _candidate("best", IK_BEST, cosine=0.95, final=0.70)
    close = _candidate("close isomer", IK_ALT, cosine=0.935, final=0.69)
    best["top_candidates"] = [best.copy(), close]

    result = annotate(best)

    assert result["isomer_ambiguous"] is True
    assert result["isomer_n"] == 1
    assert result["compound_confidence"] == "Medium"


def test_cross_library_isomer_with_moderate_cosine_gap_remains_ambiguous():
    """Atractylenolide-II/Costunolide-like score gaps are not exact-ID proof."""
    best = _candidate("best", IK_BEST, cosine=0.9384, final=0.8601)
    plausible = _candidate("plausible isomer", IK_ALT, cosine=0.8890, final=0.8200)
    best["top_candidates"] = [best.copy(), plausible]

    result = annotate(best)

    assert result["isomer_ambiguous"] is True
    assert result["compound_confidence"] == "Medium"


def test_formula_is_derived_from_smiles_for_isomer_guard():
    best = _candidate("best", IK_BEST, formula="", cosine=0.95, final=0.70)
    same_formula = _candidate(
        "same-formula structure", IK_ALT, formula="", cosine=0.93, final=0.69
    )
    best["top_candidates"] = [best.copy(), same_formula]

    result = annotate(best)

    assert result["same_formula_alternatives_n"] == 1
    assert result["isomer_ambiguous"] is True


def test_ionized_and_neutral_formulae_share_the_isomer_guard():
    best = _candidate(
        "deprotonated record",
        IK_BEST,
        formula="C21H19O12-",
        cosine=0.95,
        final=0.70,
    )
    best.update({
        "smiles": "O=c1cc(-c2ccc(O)c([O-])c2)oc2cc(O)cc(O)c12",
        "adduct": "[M-H]-",
    })
    neutral_isomer = _candidate(
        "neutral-formula isomer",
        IK_ALT,
        formula="C21H20O12",
        cosine=0.94,
        final=0.69,
    )
    neutral_isomer["adduct"] = "[M-H]-"
    best["top_candidates"] = [best.copy(), neutral_isomer]

    result = annotate(best)

    assert result["same_formula_alternatives_n"] == 1
    assert result["isomer_ambiguous"] is True
    assert result["compound_confidence"] == "Medium"


def test_precursor_structure_adduct_conflict_cannot_be_compound_high():
    best = _candidate("best", IK_BEST, cosine=0.96, final=0.80)
    best["mass_metadata_consistent"] = False
    best["top_candidates"] = [best.copy()]

    result = annotate(best)

    assert result["compound_confidence"] == "Medium"
    assert result["class_confidence"] == "High"


def test_sparse_three_peak_match_without_adjacent_ms1_precursor_is_not_high():
    best = _candidate("sparse false-high guard", IK_BEST, cosine=0.96, final=0.80)
    best.update({
        "query_search_precursor_mz_source": "ms2_selected_ion_mz_fallback",
        "query_adduct": "",
        "query_adduct_reliable": False,
        "isotope_evidence_status": "precursor_not_observed_in_ms1",
        "multimetric_shadow": {
            "direct_cosine_hungarian": {
                "status": "available",
                "score": 0.96,
                "matched_peaks": 3,
            }
        },
    })
    best["top_candidates"] = [best.copy()]

    result = annotate(best)

    assert result["compound_confidence"] == "Medium"
    assert result["isomer_group_confidence"] == "Medium"
    assert result["precursor_high_evidence_gate_required"] is True
    assert result["precursor_high_evidence_gate_result"] == (
        "vetoed_to_medium_sparse_ms2_without_precursor_support"
    )
    assert "仅有 3 个直接匹配的产物离子" in result[
        "precursor_high_evidence_gate_reason_zh"
    ]
    assert result["direct_matched_product_ions"] == 3


def test_four_peak_match_passes_sparse_precursor_minimum_support_gate():
    best = _candidate("four product ions", IK_BEST, cosine=0.96, final=0.80)
    best.update({
        "query_search_precursor_mz_source": "ms2_selected_ion_mz_fallback",
        "query_adduct": "",
        "query_adduct_reliable": False,
        "isotope_evidence_status": "precursor_not_observed_in_ms1",
        "multimetric_shadow": {
            "direct_cosine_hungarian": {
                "status": "available",
                "score": 0.96,
                "matched_peaks": 4,
            }
        },
    })
    best["top_candidates"] = [best.copy()]

    result = annotate(best)

    assert result["compound_confidence"] == "High"
    assert result["precursor_high_evidence_gate_required"] is False
    assert result["precursor_high_evidence_gate_result"] == (
        "passed_minimum_direct_product_ion_support"
    )


def test_margin_uses_distinct_competitor_when_selected_best_was_not_topk_zero():
    best = _candidate("selected", IK_BEST, cosine=0.85, final=0.56, rank=0.62)
    far = _candidate(
        "far competitor", IK_ALT, formula="C11H14O2",
        cosine=0.60, final=0.40, rank=0.63,
    )
    # Reproduces an old cache after isomer tie-break: another entry was stored
    # first and the selected best appeared at tc[1].
    best["top_candidates"] = [far, best.copy()]

    result = annotate(best)

    # A large margin is useful descriptive evidence, but the historic 443-query
    # sweep did not calibrate a production High tier.  It no longer promotes a
    # base Medium result.
    assert result["compound_confidence"] == "Medium"
    assert result["class_agreement"] == 1.0


def test_untrusted_ml_class_votes_cannot_unlock_margin_promotion(monkeypatch):
    monkeypatch.setattr(honest_calibration, "USE_MARGIN_UPGRADE", True)
    best = _candidate("selected", IK_BEST, cosine=0.85, final=0.56)
    best["classification_source"] = "ml"
    competitor = _candidate(
        "far", IK_ALT, formula="C11H14O2", cosine=0.60, final=0.40
    )
    competitor["classification_source"] = "ml_research_unvalidated"
    best["top_candidates"] = [best.copy(), competitor]

    result = annotate(best)

    assert result["class_consensus"] == "other"
    assert result["class_agreement"] == 0.0
    assert result["compound_confidence"] == "Medium"


def test_even_trusted_margin_signal_is_counterfactual_only(monkeypatch):
    monkeypatch.setattr(honest_calibration, "USE_MARGIN_UPGRADE", True)
    best = _candidate("selected", IK_BEST, cosine=0.85, final=0.56)
    best["classification_source"] = "curated"
    competitor = _candidate(
        "far", IK_ALT, formula="C11H14O2", cosine=0.60, final=0.40
    )
    competitor["classification_source"] = "curated"
    best["top_candidates"] = [best.copy(), competitor]

    result = annotate(best)

    assert result["margin_upgrade_counterfactual_eligible"] is True
    assert result["margin_upgrade_runtime_influence"] is False
    assert result["compound_confidence"] == "Medium"


def test_shadow_pathway_cannot_promote_medium_to_high(monkeypatch):
    monkeypatch.setattr(honest_calibration, "USE_PATHWAY_CONFIDENCE", True)
    best = _candidate("selected", IK_BEST, cosine=0.85, final=0.56)
    best["expert_details"] = {
        "pathway_align": {
            "score": 0.99,
            "candidate_conditioned": True,
            "independent_evidence": False,
            "eligible_for_fusion": False,
            "runtime_influence": False,
        }
    }
    best["top_candidates"] = [best.copy()]

    result = annotate(best)

    assert result["compound_confidence"] == "Medium"
    assert result["confidence_semantics"] == (
        "rule_based_operational_evidence_tier"
    )
    assert result["confidence_is_probability"] is False
    assert result["confidence_is_project_qvalue"] is False


def test_legacy_cached_pathway_influence_fails_closed_until_reanalysis():
    best = _candidate("legacy", IK_BEST, cosine=0.96, final=0.80)
    best["expert_details"] = {
        "expert_total_score": 0.10,
        "pathway_align": {"score": 0.90},
        "score_fusion": {"method": "max_of_independent_expert_sources"},
    }
    best["top_candidates"] = [best.copy()]

    result = annotate(best)

    assert result["legacy_unvalidated_expert_influence_detected"] is True
    assert result["compound_confidence"] == "Medium"


def test_deconvolved_cosine_cannot_expand_candidate_margin():
    best = _candidate("selected", IK_BEST, cosine=0.89, final=0.60, rank=0.56)
    best.update({
        "raw_cosine": 0.80,
        "raw_final_score": 0.50,
        "raw_rank_score": 0.56,
        "chimeric_deconvolution": {
            "status": "accepted_score_applied",
            "accepted": True,
            "runtime_influence": True,
            "isolation_status": "mixed",
        },
    })
    competitor = _candidate(
        "competitor",
        IK_ALT,
        formula="C11H14O2",
        cosine=0.75,
        final=0.48,
        rank=0.50,
    )
    best["top_candidates"] = [best.copy(), competitor]

    result = annotate(best)

    # The deconvolved view would create a 0.14 margin and trigger the enabled
    # margin-upgrade path.  Candidate comparison must retain the raw 0.05 gap.
    assert result["candidate_comparison_source"] == "raw_measured_spectrum"
    assert result["candidate_comparison_cosine"] == 0.80
    assert result["candidate_margin"] == 0.05
    assert result["compound_confidence"] == "Medium"


def test_deconvolution_cannot_hide_raw_same_formula_competitor():
    best = _candidate("selected", IK_BEST, cosine=0.94, final=0.72, rank=0.70)
    best.update({
        "raw_cosine": 0.86,
        "raw_final_score": 0.62,
        "raw_rank_score": 0.70,
        "chimeric_deconvolution": {
            "status": "accepted_score_applied",
            "accepted": True,
            "runtime_influence": True,
            "isolation_status": "mixed",
        },
    })
    same_formula = _candidate(
        "raw competitive isomer",
        IK_ALT,
        cosine=0.84,
        final=0.60,
        rank=0.60,
    )
    best["top_candidates"] = [best.copy(), same_formula]

    result = annotate(best)

    assert result["candidate_margin"] == 0.02
    assert result["isomer_ambiguous"] is True
    assert result["isomer_set"][0]["cosine"] == 0.84
    assert result["compound_confidence"] == "Medium"


def test_class_consensus_does_not_double_vote_best():
    best = _candidate("best", IK_BEST, classification="flavonoid")
    other = _candidate(
        "other", IK_ALT, formula="C11H14O2", classification="alkaloid",
        cosine=0.60, final=0.40,
    )
    best["top_candidates"] = [best.copy(), other]

    result = annotate(best)

    assert result["class_consensus"] == "flavonoid"
    assert result["class_agreement"] == 0.5


def test_public_top_candidates_always_begin_with_selected_best():
    selected = _candidate("selected", IK_BEST, rank=0.899)
    higher_ranked = _candidate("rank winner", IK_ALT, rank=0.900)

    ordered = _top_candidates_with_best(selected, [higher_ranked, selected])

    assert ordered[0]["inchikey"] == IK_BEST
    assert ordered[1]["inchikey"] == IK_ALT


def test_public_top_candidates_includes_selected_best_missing_from_cache():
    selected = _candidate("selected", IK_BEST, cosine=0.92, rank=0.61)
    selected.update({
        "raw_cosine": 0.81,
        "raw_final_score": 0.53,
        "raw_rank_score": 0.61,
    })
    competitor = _candidate("cached only", IK_ALT, cosine=0.85, rank=0.60)

    ordered = _top_candidates_with_best(selected, [competitor])

    assert ordered[0]["inchikey"] == IK_BEST
    assert ordered[0]["cosine"] == 0.81
    assert ordered[1]["inchikey"] == IK_ALT


def test_any_post_deconvolution_reselection_still_uses_raw_tie_breaks():
    deconvolved = _candidate(
        "deconvolved", IK_BEST, cosine=0.99, final=0.82, rank=0.60
    )
    deconvolved.update({
        "raw_cosine": 0.80,
        "raw_final_score": 0.58,
        "raw_rank_score": 0.60,
        "raw_reverse_match": 0.80,
    })
    raw_winner = _candidate(
        "raw winner", IK_ALT, formula="C11H14O2",
        cosine=0.85, final=0.60, rank=0.60,
    )

    selected = _select_best_result([deconvolved, raw_winner])

    assert selected["inchikey"] == IK_ALT


def test_selected_best_is_not_lost_when_raw_rank_is_below_topk_limit():
    selected = _candidate("selected", IK_BEST, rank=0.50)
    candidates = [
        _candidate(
            f"candidate-{index}",
            f"{chr(67 + index) * 14}-{chr(68 + index) * 10}-{chr(69 + index)}",
            formula=f"C{11 + index}H14O2",
            rank=0.90 - index * 0.01,
        )
        for index in range(5)
    ]

    ordered = _top_candidates_with_best(selected, [*candidates, selected], limit=5)

    assert len(ordered) == 5
    assert ordered[0]["inchikey"] == IK_BEST


def test_historical_evaluator_uses_production_match_as_top1():
    selected = _candidate("selected", IK_BEST, rank=0.899)
    stale_first = _candidate("stale first", IK_ALT, rank=0.900)
    selected["top_candidates"] = [stale_first, selected.copy()]

    candidates = _candidate_list(selected)

    assert candidates[0]["inchikey"] == IK_BEST
    assert candidates[1]["inchikey"] == IK_ALT


def test_cached_display_reorders_stale_topk_before_honesty_calibration():
    from web.server import _prepare_results_for_display

    selected = _candidate("selected", IK_BEST, rank=0.899)
    stale_first = _candidate("stale first", IK_ALT, rank=0.900)
    selected["top_candidates"] = [stale_first, selected.copy()]
    rows = [{
        "spectrum_id": "1",
        "precursor_mz": 300.0,
        "retention_time": 1.0,
        "ion_mode": "POS",
        "collision_energy": 20.0,
        "match": selected,
    }]

    prepared = _prepare_results_for_display(rows)

    assert prepared[0]["match"]["inchikey"] == IK_BEST
    assert prepared[0]["match"]["top_candidates"][0]["inchikey"] == IK_BEST
    assert prepared[0]["match"]["top_candidates"][1]["inchikey"] == IK_ALT


def test_cached_display_upgrades_ao4703_name_role_and_wrong_ml_class():
    from web.server import _prepare_results_for_display

    selected = _candidate(
        "AO4703",
        "VMZVBRIIHDRYGK-UHFFFAOYSA-N",
        formula="C17H29NO",
        classification="lignan",
    )
    selected["smiles"] = "CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)CN(C)C"
    selected["display_class"] = "lignan"
    selected["top_candidates"] = [selected.copy()]
    rows = [{
        "spectrum_id": "ao4703",
        "precursor_mz": 264.2322,
        "retention_time": 1.0,
        "ion_mode": "POS",
        "collision_energy": 20.0,
        "match": selected,
    }]

    match = _prepare_results_for_display(rows)[0]["match"]

    assert match["compound_name"] == (
        "2,6-Di-tert-butyl-4-(dimethylaminomethyl)phenol (AO4703)"
    )
    assert match["source_compound_name"] == "AO4703"
    assert match["compound_name_cn"].endswith("（AO4703）")
    assert match["classification"] == "other"
    assert match["display_class"] == "other"
    assert match["potential_contaminant"] is True
    assert match["contaminant_category"] == "synthetic_phenolic_antioxidant"
    assert "外源污染物" in match["honest_headline"]
    assert "lignan" not in match["honest_headline"].lower()


def test_cached_display_keeps_topk_selected_score_on_raw_spectrum_view():
    from web.server import _prepare_results_for_display

    selected = _candidate("selected", IK_BEST, cosine=0.96, final=0.76, rank=0.58)
    selected.update({
        "raw_cosine": 0.82,
        "raw_final_score": 0.54,
        "raw_rank_score": 0.58,
        "raw_reverse_match": 0.81,
        "chimeric_deconvolution": {
            "status": "accepted_score_applied",
            "accepted": True,
            "runtime_influence": True,
            "isolation_status": "mixed",
        },
    })
    competitor = _candidate(
        "competitor", IK_ALT, formula="C11H14O2",
        cosine=0.79, final=0.52, rank=0.55,
    )
    selected["top_candidates"] = [selected.copy(), competitor]
    rows = [{
        "spectrum_id": "deconvolved",
        "precursor_mz": 300.0,
        "retention_time": 1.0,
        "ion_mode": "POS",
        "collision_energy": 20.0,
        "match": selected,
    }]

    prepared = _prepare_results_for_display(rows)
    match = prepared[0]["match"]

    assert match["cosine"] == 0.96
    assert match["top_candidates"][0]["cosine"] == 0.82
    assert (
        match["top_candidates"][0]["candidate_comparison_source"]
        == "raw_measured_spectrum"
    )
    assert match["candidate_margin"] == 0.03


def test_completed_status_refreshes_stale_confidence_summary(monkeypatch):
    from web import server

    selected = _candidate(
        "selected", IK_BEST, formula="C8H10O", cosine=0.97, final=0.80
    )
    selected["top_candidates"] = [selected.copy()]
    selected["chimeric_deconvolution"] = {
        "status": "poor_reconstruction",
        "accepted": False,
        "runtime_influence": False,
        "isolation_status": "severe",
    }
    task_id = "stale-summary"
    task = {
        "id": task_id,
        "filename": "sample.mzML",
        "status": "completed",
        "progress": 100,
        "message": "old",
        "total_spectra": 1,
        "matched_spectra": 1,
        "match_rate": 100.0,
        "compound_confidence_summary": {
            "high_supported_unique_compounds": 0,
            "representative_high_unique_compounds": 0,
            "high_spectra": 0,
        },
        "params": {},
        "ion_identity_summary": {},
        "results": [{
            "spectrum_id": "1",
            "precursor_mz": 300.0,
            "retention_time": 1.0,
            "ion_mode": "POS",
            "collision_energy": 20.0,
            "match": selected,
        }],
    }
    monkeypatch.setattr(server, "_tasks", {task_id: task})
    monkeypatch.setattr(server, "_save_tasks", lambda: None)

    status = asyncio.run(server.get_task_status(task_id))

    summary = status["compound_confidence_summary"]
    assert summary["high_supported_unique_compounds"] == 1
    assert summary["representative_high_unique_compounds"] == 1
    assert summary["high_spectra"] == 1
    assert status["high_precursor_features"] == 1
    assert (
        task["presentation_calibration_version"]
        == server.PRESENTATION_CALIBRATION_VERSION
    )
