from copy import deepcopy

from eval.v16_fragment_discrimination_20260907.projection_repair import (
    normalize_projection, reverse_row, select_with_frozen_baseline,
)


def fixture():
    fallback = {"smiles":"CCO", "inchikey":"LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "lib_spectrum_id":"measured_A",
        "candidate_comparison_source":"measured_spectrum", "cosine":.9}
    explicit = {"smiles":"CCO", "inchikey":"LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "lib_spectrum_id":"measured_B",
        "cosine":.8, "reference_evidence_records":[{"lib_spectrum_id":"measured_B", "cosine":.8}]}
    return {"match":{**fallback, "top_candidates":[deepcopy(fallback), explicit]}}


def ref_ids(row):
    match = row["match"]
    return [tuple(r["lib_spectrum_id"] for r in c.get("reference_evidence_records", [])) for c in [match]+match["top_candidates"]]


def test_compact_precedence_is_independent_of_alias_order():
    a, count = normalize_projection(fixture(), {}, {}, "compact_only")
    b, other = normalize_projection(reverse_row(fixture()), {}, {}, "compact_only")
    assert ref_ids(a) == ref_ids(b) == [("measured_B",)] * 3
    assert count == other == 0


def test_distinct_resolved_measured_reference_is_not_silently_discarded():
    row, count = normalize_projection(fixture(), {}, {"measured_A":"fp_a", "measured_B":"fp_b"}, "measured_union")
    assert count == 1
    assert all(set(ids)=={"measured_A", "measured_B"} for ids in ref_ids(row))


def test_compact_reference_wins_over_fallback_for_the_same_record_id():
    row = fixture()
    row["match"]["lib_spectrum_id"] = "measured_B"
    row["match"]["top_candidates"][0]["lib_spectrum_id"] = "measured_B"
    fixed, count = normalize_projection(row, {}, {"measured_B":"fp_b"}, "measured_union")
    assert count == 0
    assert ref_ids(fixed) == [("measured_B",)] * 3
    assert fixed["match"]["reference_evidence_records"][0]["cosine"] == .8


def test_unresolved_or_predicted_reference_does_not_enter_measured_union():
    row = fixture()
    for c in [row["match"]]+row["match"]["top_candidates"]:
        if not c.get("reference_evidence_records"):
            c["candidate_comparison_source"] = "predicted_spectrum"
    fixed, count = normalize_projection(row, {}, {"measured_A":"fp_a"}, "measured_union")
    assert count == 0
    assert ref_ids(fixed) == [("measured_B",)] * 3
    unresolved, count = normalize_projection(fixture(), {}, {}, "measured_union")
    assert count == 0 and ref_ids(unresolved) == [("measured_B",)] * 3


def test_no_compact_evidence_keeps_existing_fallback():
    row = fixture()
    row["match"]["top_candidates"].pop()
    fixed, count = normalize_projection(row, {}, {"measured_A":"fp_a"}, "measured_union")
    assert count == 0 and fixed == row


def test_input_and_rank_metadata_remain_unchanged():
    row = fixture()
    before = deepcopy(row)
    normalize_projection(row, {}, {"measured_A":"fp_a", "measured_B":"fp_b"}, "measured_union")
    assert row == before


def test_cross_occurrence_decision_is_protected_by_saved_native_rank():
    candidates = [{"connectivity":"A", "rank":1, "origin":"native"}, {"connectivity":"B", "rank":2, "origin":"native"}]
    assert select_with_frozen_baseline("B", candidates, "A") == "B"
    assert select_with_frozen_baseline("A", candidates, "B") == "B"
    assert select_with_frozen_baseline("A", list(reversed(candidates)), "B") == "B"


def test_background_identity_cannot_be_promoted_into_native_pool():
    candidates = [{"connectivity":"A", "rank":1, "origin":"native"}, {"connectivity":"C", "rank":1, "origin":"expanded"}]
    assert select_with_frozen_baseline("A", candidates, "C") == "A"
