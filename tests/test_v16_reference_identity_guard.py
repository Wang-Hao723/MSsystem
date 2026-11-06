"""Offline library-identity checks; no service, database, or user data writes."""
from copy import deepcopy

from eval.v16_fragment_discrimination_20260907.identity_guard import (
    clean_candidate, conflicting_contents, key, native_first,
)


def panel():
    return {"candidates":[
        {"connectivity":"A", "references":[{"fingerprint":"shared"}, {"fingerprint":"a"}]},
        {"connectivity":"B", "references":[{"fingerprint":"shared"}, {"fingerprint":"b"}]},
    ]}


def test_different_structures_sharing_content_are_flagged_symmetrically():
    assert conflicting_contents(panel(), {"A":"A", "B":"B"}) == {"shared":["A", "B"]}


def test_duplicate_content_with_one_tautomer_identity_is_not_a_label_conflict():
    assert conflicting_contents(panel(), {"A":"A", "B":"A"}) == {}


def test_conflict_detection_is_reference_count_and_order_invariant():
    e = panel()
    expected = conflicting_contents(e, {"A":"A", "B":"B"})
    e["candidates"].reverse()
    for c in e["candidates"]:
        c["references"] *= 7
        c["references"].reverse()
    assert conflicting_contents(e, {"A":"A", "B":"B"}) == expected


def test_both_conflicting_reference_labels_lose_their_spectral_evidence():
    fingerprints = {"ref_A":"shared", "ref_B":"shared"}
    for sid in fingerprints:
        assert clean_candidate({"lib_spectrum_id":sid, "cosine":1.0}, {"shared":[]}, fingerprints, {}) is None


def test_exhausted_compact_evidence_cannot_reuse_a_saved_score():
    c = {"lib_spectrum_id":"old", "cosine":1., "reference_evidence_records":[{"lib_spectrum_id":"bad", "cosine":1.}]}
    assert clean_candidate(c, {"shared":[]}, {"bad":"shared"}, {}) is None


def test_unresolved_fingerprint_does_not_prove_a_conflict():
    c = {"lib_spectrum_id":"unresolved", "cosine":.8}
    assert clean_candidate(c, {"shared":[]}, {}, {}) == c


def test_surviving_reference_replaces_all_old_score_fields_without_mutation():
    c = {"lib_spectrum_id":"bad", "raw_cosine":1., "raw_entropy":1., "smiles":"CCO",
        "reference_evidence_records":[{"lib_spectrum_id":"good", "cosine":.4}, {"lib_spectrum_id":"bad", "cosine":1.}]}
    original = deepcopy(c)
    cleaned = clean_candidate(c, {"shared":[]}, {"bad":"shared", "good":"unique"}, {"scan_score":{"cosine":1.}})
    assert c == original
    assert cleaned["lib_spectrum_id"] == "good"
    assert cleaned["cosine"] == .4
    assert "raw_cosine" not in cleaned and "raw_entropy" not in cleaned
    assert len(cleaned["reference_evidence_records"]) == 1


def test_native_abstention_is_not_dropped_from_top1_analysis():
    assert native_first({"winner":None, "candidates":[{"identity_key":"std2d:A"}]}) == "A"
    assert native_first({"winner":None, "candidates":[]}) == ""


def test_production_cyperone_tautomer_transform_is_not_new_structure_evidence():
    a = {"smiles":"C=C(C)C1C=C2C(C)C(=O)CCC2(C)CC1"}
    b = {"smiles":"C=C(C)C1CCC2(C)CCC(=O)C(C)=C2C1"}
    assert key(a) == key(b)


def test_ferulate_and_acetyl_mandelate_are_not_collapsed():
    a = {"smiles":"COc1cc(C=CC(=O)O)ccc1O"}
    b = {"smiles":"CC(=O)OC(C(=O)O)c1ccccc1"}
    assert key(a) != key(b)
