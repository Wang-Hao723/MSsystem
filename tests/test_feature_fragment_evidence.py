import copy
import json
import math
from pathlib import Path

import pytest

from matching.feature_fragments.catalog import load_catalog
from matching.feature_fragments.evidence import (
    conditional_log_likelihood_ratio,
    evaluate_feature_fragment_evidence,
    mass_window,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = (
    ROOT
    / "tests"
    / "fixtures"
    / "feature_fragments"
    / "massbank_antwerp_phospholipids_v1.json"
)


def _fixture():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _record(accession_suffix):
    return next(
        record
        for record in _fixture()["spectra"]
        if record["spectrum_id"].endswith(accession_suffix)
    )


def _context(record, chemical_class=None):
    acquisition = record["acquisition"]
    return {
        "chemical_class": chemical_class or record["compound"]["chemical_class"],
        "classification_reliable": True,
        "classification_source": "ChemOnt_from_candidate_structure",
        "polarity": acquisition["polarity"],
        "adduct": acquisition["adduct"],
        "charge": acquisition["charge"],
        "instrument_domain_id": acquisition["instrument_domain_id"],
        "activation": acquisition["activation"],
        "collision_energy": acquisition["collision_energy"],
        "quality": record["quality"],
    }


def _peaks(record):
    return [[peak[0], peak[1]] for peak in record["peaks"]]


def test_pc_184_is_positive_class_support_but_never_identity_confirmation():
    record = _record("N105906_EF88")
    result = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=_context(record),
        catalog=load_catalog(),
    )

    informative = [
        item for item in result["rule_evidence"] if item["status"] == "present"
    ]
    assert len(informative) == 1
    assert informative[0]["rule_id"] == "FF-PC-184-POS-H-AG6560-20EV"
    assert informative[0]["log_likelihood_ratio"] > 0
    assert result["combined_log_likelihood_ratio"] > 0
    assert result["evidence_stage"] == "msms_rule_support"
    assert result["runtime_influence"] is False
    assert result["identity_rank_influence"] is False
    assert result["score_contribution"] == 0.0
    assert result["exact_identity_confirmed"] is False
    assert result["standard_confirmed"] is False
    assert result["isomer_discriminating"] is False
    standard = result["feature_fragment_score_standard"]
    assert standard["primary_unit"] == "natural_log_likelihood_ratio_nats"
    assert standard["count_alone_has_fixed_points"] is False
    assert result["present_independent_fragment_count"] == 1
    assert result["combined_score_nats"] == pytest.approx(
        result["combined_log_likelihood_ratio"]
    )
    assert result["combined_score_bits"] == pytest.approx(
        result["combined_score_nats"] / math.log(2.0)
    )
    step = result["cumulative_evidence_by_independent_count"][0]
    assert step["independent_evidence_count"] == 1
    assert step["selected_rule_id"] == "FF-PC-184-POS-H-AG6560-20EV"
    assert step["cumulative_score_nats"] == pytest.approx(
        result["combined_score_nats"]
    )


def test_nce_and_ev_are_not_compared_and_composite_ce_does_not_match_fixed_rule():
    record = _record("N105906_EF88")
    nce_context = _context(record)
    nce_context["collision_energy"] = {"value": 20, "unit": "NCE", "mode": "fixed"}
    nce = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=nce_context,
        catalog=load_catalog(),
    )
    assert nce["combined_log_likelihood_ratio"] == 0.0
    assert any(
        item["reason"] == "collision_energy_unit_mismatch"
        for item in nce["rule_evidence"]
    )

    ramp_context = _context(record)
    ramp_context["collision_energy"] = {
        "value": None,
        "unit": "eV",
        "mode": "ramp",
        "minimum": 10,
        "maximum": 30,
    }
    ramp = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=ramp_context,
        catalog=load_catalog(),
    )
    assert ramp["combined_log_likelihood_ratio"] == 0.0
    assert any(
        item["reason"] == "collision_energy_mode_mismatch"
        for item in ramp["rule_evidence"]
    )


def test_absence_is_uninformative_until_every_quality_gate_passes():
    record = _record("N105906_EF88")
    peaks_without_184 = [
        peak for peak in _peaks(record) if abs(peak[0] - 184.0733) > 0.02
    ]
    context = _context(record)
    context["quality"] = {**record["quality"], "precursor_isolated": False}
    gated = evaluate_feature_fragment_evidence(
        peaks=peaks_without_184,
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=context,
        catalog=load_catalog(),
    )
    assert gated["combined_log_likelihood_ratio"] == 0.0
    assert any(
        item["status"] == "uninformative_absence"
        and "precursor_isolated" in item["missing_quality_fields"]
        for item in gated["rule_evidence"]
    )

    eligible = evaluate_feature_fragment_evidence(
        peaks=peaks_without_184,
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=_context(record),
        catalog=load_catalog(),
    )
    assert eligible["combined_log_likelihood_ratio"] < 0
    assert any(
        item["status"] == "eligible_absence" for item in eligible["rule_evidence"]
    )


def test_pe_141_response_preserves_the_observed_40_ev_target_absence():
    record = _record("N106006_FB57")
    result = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=_context(record),
        catalog=load_catalog(),
    )
    informative = [
        item
        for item in result["rule_evidence"]
        if item["status"] in {"present", "eligible_absence"}
    ]
    assert len(informative) == 1
    assert informative[0]["rule_id"] == "FF-PE-NL141-POS-H-AG6560-40EV"
    assert informative[0]["status"] == "eligible_absence"
    assert informative[0]["log_likelihood_ratio"] < 0


def test_circular_or_unreliable_classification_cannot_unlock_rules():
    record = _record("N105906_EF88")
    circular = _context(record)
    circular["classification_source"] = "feature_fragments"
    result = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=circular,
        catalog=load_catalog(),
    )
    assert result["status"] == "not_applicable"
    assert result["reason"] == "circular_classification_source"


def test_likelihood_ratio_and_mass_window_use_counts_and_rule_specific_tolerance():
    lr = conditional_log_likelihood_ratio(
        target_present=6,
        target_total=6,
        background_present=0,
        background_total=4,
        observed_present=True,
    )
    assert lr["likelihood_ratio"] > 1
    absence = conditional_log_likelihood_ratio(
        target_present=6,
        target_total=6,
        background_present=0,
        background_total=4,
        observed_present=False,
    )
    assert absence["likelihood_ratio"] < 1
    assert mass_window(1000.0, {"absolute_da": 0.005, "ppm": 20}) == pytest.approx(0.02)
    assert mass_window(100.0, {"absolute_da": 0.005, "ppm": 20}) == pytest.approx(0.005)


def test_correlated_rules_contribute_only_one_fragmentation_group_term():
    record = _record("N105906_EF88")
    catalog = load_catalog()
    applicable = next(
        rule
        for rule in catalog["rules"]
        if rule["rule_id"] == "FF-PC-184-POS-H-AG6560-20EV"
    )
    duplicate = copy.deepcopy(applicable)
    duplicate["rule_id"] = "TEST-CORRELATED-DUPLICATE"
    duplicate["occurrence_statistics"] = {
        "present": 4,
        "total": 6,
        "frequency": 4 / 6,
        "counting_unit": "independent_2d_identity",
    }
    duplicate_catalog = {**catalog, "rules": [applicable, duplicate]}
    result = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=_context(record),
        catalog=duplicate_catalog,
    )
    present = [item for item in result["rule_evidence"] if item["status"] == "present"]
    assert len(present) == 2
    assert result["independent_group_representatives"] == [
        "phosphocholine_headgroup_184"
    ]
    assert result["combined_log_likelihood_ratio"] == max(
        item["log_likelihood_ratio"] for item in present
    )


def test_two_independent_fragments_sum_their_rule_specific_scores():
    record = _record("N105906_EF88")
    catalog = load_catalog()
    applicable = next(
        rule
        for rule in catalog["rules"]
        if rule["rule_id"] == "FF-PC-184-POS-H-AG6560-20EV"
    )
    second = copy.deepcopy(applicable)
    second["rule_id"] = "TEST-INDEPENDENT-SECOND"
    second["independent_fragmentation_group"] = "test_independent_second"
    second["occurrence_statistics"] = {
        "present": 4,
        "total": 6,
        "frequency": 4 / 6,
        "counting_unit": "independent_2d_identity",
    }
    result = evaluate_feature_fragment_evidence(
        peaks=_peaks(record),
        precursor_mz=record["acquisition"]["precursor_mz"],
        context=_context(record),
        catalog={**catalog, "rules": [applicable, second]},
    )
    present = [item for item in result["rule_evidence"] if item["status"] == "present"]

    assert result["present_independent_fragment_count"] == 2
    assert len(result["cumulative_evidence_by_independent_count"]) == 2
    assert result["combined_score_nats"] == pytest.approx(
        sum(item["score_nats"] for item in present)
    )
    assert result["cumulative_evidence_by_independent_count"][-1][
        "cumulative_score_nats"
    ] == pytest.approx(result["combined_score_nats"])
