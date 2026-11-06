import copy

import pytest

from matching.ion_identity import PROTON
from matching.v16_ms1_ion_network import (
    VARIANTS, build_ms1_ion_network, ms1_trace_targets, prepare_ms1_network_input,
)


def inputs(*, linked_partner=False, polarity="NEG"):
    shifts = (-PROTON, 44.998201) if polarity == "NEG" else (PROTON, 22.989218)
    features = [{"feature_id": fid, "mz": 300 + shift, "rt_minutes": 1.0,
                 "rt_seconds": 60., "charge": -1 if polarity == "NEG" else 1}
                for fid, shift in zip(("F1", "F2"), shifts)]
    rows = [{"precursor_mz": f["mz"] - .003, "retention_time": .99,
             "ion_mode": polarity, "lcms_feature_id": f["feature_id"],
             "lcms_feature_mz": f["mz"], "lcms_feature_rt_minutes": 1.,
             "lcms_feature_link_status": "linked_openms_feature"}
            for f in (features if linked_partner else features[:1])]
    trace = {"rt": [.94, .96, .98, 1., 1.02, 1.04, 1.06],
             "intensity": [1., 2., 4., 8., 4., 2., 1.], "apex_intensity": 8.}
    return features, rows, {f["feature_id"]: copy.deepcopy(trace) for f in features}


def build(features, rows, traces, *, variant="all_ms1_partners", polarity="NEG"):
    prepared = prepare_ms1_network_input("task", polarity, features, rows)
    return build_ms1_ion_network(prepared, traces, variant=variant)


@pytest.mark.parametrize("polarity", ["POS", "NEG"])
def test_unfragmented_partner_supplies_ms1_evidence_not_ms2_identification(polarity):
    f, r, t = inputs(polarity=polarity)
    result = build(f, r, t, polarity=polarity)
    assert len(result["anchors"]) == 1
    anchor = result["anchors"][0]
    assert anchor["neutral_mass"] == pytest.approx(300.)
    assert anchor["source_row_indices"] == [0]
    assert anchor["feature_ids"] == ["F1", "F2"]
    assert anchor["verified_member_feature_count"] == 2
    assert result["networks"][0]["ms1_only_feature_ids"] == ["F2"]
    assert result["compound_count_increment"] == 0
    assert result["production_results_mutated"] is result["runtime_effect"] is False


def test_registered_ablation_excludes_unfragmented_partner():
    f, r, t = inputs()
    assert build(f, r, t, variant=VARIANTS[0])["anchors"] == []
    assert len(build(f, r, t, variant=VARIANTS[1])["anchors"]) == 1


def test_mass_coordinate_is_measured_ms1_not_selected_ion():
    f, r, t = inputs(linked_partner=True)
    r[0]["precursor_mz"] -= .009
    r[1]["precursor_mz"] += .008
    assert build(f, r, t, variant=VARIANTS[0])["anchors"][0]["neutral_mass"] == 300.


def test_candidate_poisoning_and_input_mutation_are_impossible():
    f, r, t = inputs(linked_partner=True)
    original = copy.deepcopy((f, r, t))
    baseline = build(f, r, t)
    assert (f, r, t) == original
    for row in r:
        row.update({"match": {"confidence": "High", "compound_name": "POISON",
            "final_score": 1., "inchikey": "POISON"}, "truth_connectivity": "POISON",
            "ion_identity_neutral_mass": 999., "ion_identity_network_id": "POISON"})
    for feature in f:
        feature.update({"adduct": "POISON", "neutral_mass": 999., "adduct_group_id": "POISON"})
    assert build(f, r, t) == baseline


def test_raw_row_indices_are_not_feature_or_compacted_indices():
    f, r, t = inputs()
    r = [{"ion_mode": "NEG"}] * 4 + r + [copy.deepcopy(r[0])]
    result = build(list(reversed(f)), r, t)
    assert result["anchors"][0]["source_row_indices"] == [4, 5]
    assert result["networks"][0]["member_evidence"][1]["source_row_indices"] == []


def test_ms1_only_network_is_retained_but_not_promoted_to_ms2_anchor():
    f, _, t = inputs()
    result = build(f, [], t)
    assert len(result["networks"]) == 1
    assert not result["anchors"]
    assert result["summary"]["ms1_only_network_count"] == 1


@pytest.mark.parametrize("field,value,reason", [
    ("lcms_feature_link_status", "ambiguous_openms_feature", "no_unique_ms1_link"),
    ("lcms_feature_id", "unknown", "linked_feature_missing_or_ineligible"),
    ("lcms_feature_mz", 302., "linked_feature_coordinate_mismatch"),
    ("lcms_feature_rt_minutes", 2., "linked_feature_rt_mismatch"),
    ("precursor_mz", 100., "selected_ion_outside_existing_link_guard"),
    ("precursor_charge", 2, "unsupported_ms2_charge_magnitude"),
    ("precursor_charge", 1.5, "unsupported_ms2_charge_magnitude"),
])
def test_invalid_ms2_links_never_become_provenance(field, value, reason):
    f, r, t = inputs()
    r[0][field] = value
    prepared = prepare_ms1_network_input("task", "NEG", f, r)
    assert prepared["row_bindings"][0]["status"] == reason
    assert build_ms1_ion_network(prepared, t)["anchors"] == []


@pytest.mark.parametrize("field,value,reason", [
    ("mz", float("nan"), "invalid_measured_coordinates"),
    ("mz", -1., "invalid_measured_coordinates"),
    ("mz", True, "invalid_measured_coordinates"),
    ("rt_minutes", -1., "invalid_measured_coordinates"),
    ("rt_minutes", float("nan"), "invalid_measured_coordinates"),
    ("rt_seconds", float("inf"), "invalid_measured_coordinates"),
    ("rt_seconds", 100., "inconsistent_rt_units"),
    ("charge", 2, "unsupported_charge_magnitude"),
    ("charge", float("inf"), "unsupported_charge_magnitude"),
])
def test_invalid_ms1_features_are_reported(field, value, reason):
    f, r, t = inputs()
    f[1][field] = value
    prepared = prepare_ms1_network_input("task", "NEG", f, r)
    assert prepared["rejected_features"][0]["reason"] == reason
    assert not build_ms1_ion_network(prepared, t)["anchors"]


def test_unknown_charge_remains_an_explicit_hypothesis():
    f, r, t = inputs()
    f[1]["charge"] = 0
    evidence = build(f, r, t)["networks"][0]["member_evidence"][1]
    assert evidence["charge_evidence"] == "not_measured_single_charge_hypothesis"


def test_seconds_conversion_and_deterministic_feature_order():
    f, r, t = inputs()
    first = build(f, r, t)
    for feature in f:
        feature.pop("rt_minutes")
    assert build(list(reversed(f)), r, t) == first
    assert ms1_trace_targets(prepare_ms1_network_input("task", "NEG", f, r))[0]["rt"] == 1.


@pytest.mark.parametrize("mutation", ["duplicate", "empty_id", "polarity", "missing_task", "bad_variant", "tamper"])
def test_ambiguous_or_changed_input_fails_closed(mutation):
    f, r, t = inputs()
    if mutation == "duplicate":
        f[1]["feature_id"] = f[0]["feature_id"]
    elif mutation == "empty_id":
        f[0]["feature_id"] = ""
    elif mutation == "polarity":
        r[0]["ion_mode"] = "POS"
    with pytest.raises(ValueError):
        prepared = prepare_ms1_network_input("" if mutation == "missing_task" else "task", "NEG", f, r)
        if mutation == "tamper":
            prepared["features"][0]["mz"] += 1.
        build_ms1_ion_network(prepared, t, variant="invalid" if mutation == "bad_variant" else VARIANTS[1])


@pytest.mark.parametrize("kind", ["absent", "flat", "anti_correlated", "one_sided"])
def test_ms1_shape_gates_are_not_relaxed(kind):
    f, r, t = inputs()
    if kind == "absent":
        t = {}
    elif kind == "flat":
        t["F2"]["intensity"] = [1.] * 7
    elif kind == "anti_correlated":
        t["F2"]["intensity"] = [8., 4., 2., 1., 2., 4., 8.]
    else:
        t["F2"]["intensity"] = [8., 4., 2., 1., 0., 0., 0.]
    assert not build(f, r, t)["anchors"]


def test_mass_and_rt_gates_are_not_relaxed():
    f, r, t = inputs()
    f[1]["mz"] += .02
    assert not build(f, r, t)["anchors"]
    f, r, t = inputs()
    f[1]["rt_minutes"] = 1.2
    f[1]["rt_seconds"] = 72.
    assert not build(f, r, t)["anchors"]


def test_full_component_neutral_spread_is_checked(monkeypatch):
    f, r, t = inputs()
    f[1]["mz"] += .01
    def fake_native(rows, **kwargs):
        for row, adduct in zip(rows, ("[M-H]-", "[M+FA-H]-")):
            row.update({"ion_identity_status": "networked", "ion_identity_network_id": "IIN1",
                "ion_identity": adduct})
        return {"status": "ok"}
    monkeypatch.setattr("matching.v16_ms1_ion_network.annotate_ion_identity_networks", fake_native)
    result = build(f, r, t)
    assert result["anchors"] == []
    assert result["rejected_components"][0]["reason"] == "full_component_neutral_mass_spread"


def test_none_task_id_is_not_a_literal_task_named_none():
    f, r, _ = inputs()
    with pytest.raises(ValueError, match="explicit task"):
        prepare_ms1_network_input(None, "NEG", f, r)
