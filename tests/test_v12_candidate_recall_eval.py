from __future__ import annotations

import pytest

from eval.evaluate_v12_candidate_recall_mix33 import (
    _anchor_audit,
    _mode_formula_source_audit,
    exact_precursor_tolerance,
    formula_anchor_audit,
)


ATRACTYLENOLIDE_III = {
    "en_name": "Atractylenolide III",
    "formula": "C15H20O3",
    "mz": 247.1334,
    "mode": "POS",
    "inchikey": "FBMORZZOJSDNRQ-GLQYFDAESA-N",
    "cas": "",
    "pubchem_cid": "155948",
}


def test_exact_precursor_tolerance_is_bounded_da_or_ppm():
    assert exact_precursor_tolerance(100.0) == pytest.approx(0.02)
    assert exact_precursor_tolerance(3000.0) == pytest.approx(0.03)


def test_formula_audit_detects_atractylenolide_polarity_anchor_conflict():
    audit = formula_anchor_audit(ATRACTYLENOLIDE_III)

    assert audit["closest_common_adduct"]["adduct"] == "[M-H]-"
    assert audit["closest_common_adduct"]["polarity"] == "NEG"
    assert audit["closest_common_adduct_within_exact_tolerance"] is True
    assert audit["closest_adduct_matches_declared_mode"] is False
    assert audit["closest_declared_mode_adduct"]["adduct"] == "[M+H]+"
    assert audit["audit_only_not_used_for_retrieval"] is True

    anchor = _anchor_audit(
        ATRACTYLENOLIDE_III,
        {
            "entity_id": "wrong-anchor",
            "representative_precursor_mz": 247.16766866,
        },
    )
    assert anchor["status"] == "invalid"
    assert (
        anchor["reason"]
        == "ground_truth_mz_and_declared_polarity_inconsistent_with_formula"
    )
    assert anchor["blind_window_widening_used"] is False


def test_formula_consistent_sitosterol_is_classed_as_missing_measured_anchor():
    standard = {
        "en_name": "beta-Sitosterol",
        "formula": "C29H50O",
        "mz": 415.3934,
        "mode": "POS",
        "inchikey": "KZJWDPNRJALLNS-VJSFXXLFSA-N",
    }
    anchor = _anchor_audit(
        standard,
        {
            "entity_id": "off-target-entity",
            "representative_precursor_mz": 415.34344482,
        },
    )

    assert anchor["status"] == "invalid"
    assert anchor["reason"] == "measured_entity_anchor_missing_at_ground_truth_mz"
    assert anchor["formula_audit"]["closest_common_adduct"]["adduct"] == "[M+H]+"


def test_mode_formula_source_audit_reports_all_matching_scans_post_hoc():
    anchor = _anchor_audit(
        ATRACTYLENOLIDE_III,
        {
            "entity_id": "wrong-anchor",
            "representative_precursor_mz": 247.16766866,
        },
    )
    rows = [
        {
            "spectrum_id": "ce35",
            "precursor_mz": 249.1490,
            "collision_energy": 35,
            "match": {"inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C"},
        },
        {
            "spectrum_id": "ce45",
            "precursor_mz": 249.1495,
            "collision_energy": 45,
            "match": {"inchikey": ATRACTYLENOLIDE_III["inchikey"]},
        },
    ]

    audit = _mode_formula_source_audit(ATRACTYLENOLIDE_III, rows, anchor)

    assert audit is not None
    assert audit["matching_source_scan_count"] == 2
    assert audit["any_source_scan_top1_matches_truth"] is True
    assert [item["spectrum_id"] for item in audit["matching_source_scans"]] == [
        "ce35",
        "ce45",
    ]
    assert audit["audit_only_not_used_for_retrieval"] is True

