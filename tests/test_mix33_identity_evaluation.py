from eval.eval_reproducibility import evaluate_standard, targeted_feature_audit
import pytest

from eval.entropy_rerank_cv import _net_at
from eval.entropy_rerank_standards import _truth_rank
from eval.standard_names import (
    ground_truth_match,
    identifier_match,
    load_gt,
    require_complete_ground_truth,
    validate_ground_truth,
)


def _standard():
    return {
        "zh_name": "桃叶珊瑚苷",
        "en_name": "Aucubin",
        "mz": "391.1234",
        "mode": "NEG",
        "formula": "C15H22O9",
        "inchikey": "RJWJHRPNHPHBRN-FKVJWERZSA-N",
        "cas": "479-98-1",
        "pubchem_cid": "91458",
    }


def _record(scan, rt, tic, match=None):
    return {
        "spectrum_id": scan,
        "precursor_mz": 391.12341,
        "retention_time": rt,
        "tic": tic,
        "match": match,
    }


def test_aucubin_gt_carries_authoritative_identifiers():
    rows = load_gt(__import__("pathlib").Path("eval/33components_gt.csv"))
    aucubin = next(row for row in rows if row["en_name"] == "Aucubin")
    assert aucubin["inchikey"] == "RJWJHRPNHPHBRN-FKVJWERZSA-N"
    assert aucubin["cas"] == "479-98-1"
    assert aucubin["pubchem_cid"] == "91458"


def test_identifier_match_distinguishes_full_key_from_connectivity():
    std = _standard()
    assert identifier_match(std, {"inchikey": "RJWJHRPNHPHBRN-FKVJWERZSA-N"}) == "inchikey_exact"
    assert identifier_match(std, {"inchikey": "RJWJHRPNHPHBRN-UHFFFAOYSA-N"}) == "inchikey_connectivity"
    assert identifier_match(std, {"inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N"}) == ""


def test_deployment_gt_gate_reports_current_mix33_as_incomplete():
    rows = load_gt(__import__("pathlib").Path("eval/33components_gt.csv"))
    report = validate_ground_truth(rows, "connectivity")
    assert report["total"] == 33
    assert report["eligible"] == 1
    assert report["complete"] is False
    with pytest.raises(ValueError, match="1/33 eligible"):
        require_complete_ground_truth(rows, "connectivity")


def test_loose_truth_requires_formula_and_cannot_accept_name_only():
    std = {"en_name": "Aucubin", "formula": "C15H22O9", "inchikey": ""}
    assert not ground_truth_match(
        std, {"compound_name": "Aucubin", "formula": ""}, "loose"
    )
    assert ground_truth_match(
        std, {"compound_name": "Aucubin", "formula": "C15H22O9"}, "loose"
    )


def test_production_truth_rank_uses_stored_top_candidate_order():
    std = _standard()
    wrong_rank_high = {
        "compound_name": "wrong", "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        "rank_score": 0.99,
    }
    true_stored_top1 = {
        "compound_name": "Aucubin",
        "inchikey": "RJWJHRPNHPHBRN-FKVJWERZSA-N",
        "rank_score": 0.90,
    }
    candidates = [true_stored_top1, wrong_rank_high]
    assert _truth_rank(candidates, "production", std, "connectivity") == 0
    assert _truth_rank(candidates, "rank_score", std, "connectivity") == 1


def test_cv_counts_repeated_scan_flips_as_one_compound_vote():
    records = [
        ("Apigenin", 0.03, False, True),
        ("Apigenin", 0.04, False, True),
        ("Apigenin", 0.05, False, True),
        ("Monotropein", 0.03, True, False),
    ]
    assert _net_at(records, 0.02) == (1, 1, 0)


def test_targeted_audit_does_not_let_global_highest_tic_hide_later_rt_feature():
    generic_aucubin = {
        "compound_name": "CHEMBL1873725",
        "inchikey": "RJWJHRPNHPHBRN-UHFFFAOYSA-N",
        "formula": "",
        "final_score": 0.4108,
        "cosine": 0.3737,
    }
    records = [
        _record("early-high-tic", 1.8665, 10000, None),
        _record("early-second", 1.8722, 8000, None),
        _record("aucubin-feature", 2.1735, 7500, {"top_candidates": [generic_aucubin]}),
        _record("aucubin-feature-second", 2.1791, 5000, None),
    ]

    untargeted = evaluate_standard(_standard(), records)
    assert untargeted["spectrum_id"] == "early-high-tic"
    assert not untargeted["hit_name_top5"]

    audit = untargeted["targeted_feature_audit"]
    assert audit["n_precursor_features"] == 2
    assert audit["matched"] is True
    assert audit["match_kind"] == "inchikey_connectivity"
    assert audit["spectrum_id"] == "aucubin-feature"


def test_formula_only_candidate_remains_isomer_level():
    records = [_record("formula", 2.0, 100, {
        "top_candidates": [{"compound_name": "other isomer", "formula": "C15H22O9"}]
    })]
    audit = targeted_feature_audit(_standard(), records)
    assert audit["matched"] is True
    assert audit["match_kind"] == "formula_isomer"
