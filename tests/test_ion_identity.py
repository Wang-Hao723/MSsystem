import numpy as np

import matching.ion_identity as ion_identity_module
from matching.ion_identity import (
    PROTON,
    WATER,
    annotate_ion_identity_networks,
    chromatographic_shape_evidence,
    infer_ion_pair,
)


def _trace(scale=1.0, shift=0):
    rt = [round(5.00 + i * 0.01, 2) for i in range(9)]
    values = np.roll(np.asarray([1, 5, 15, 45, 100, 45, 15, 5, 1], dtype=float), shift)
    return {"rt": rt, "intensity": (values * scale).tolist()}


def _row(feature_id, mz, match=None):
    return {
        "precursor_feature_id": feature_id,
        "precursor_mz": mz,
        "retention_time": 5.04,
        "ion_mode": "POS",
        "match": match,
    }


def _match(name="Example", inchikey="AAAAAAAAAAAAAA-BBBBBBBBBB-C"):
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "formula": "C10H20O2",
        "final_score": 0.82,
        "confidence": "High",
        "adduct": "[M+H]+",
        "cross_adduct_match": False,
    }


def test_infer_adduct_and_in_source_fragment_pairs():
    neutral = 300.0
    sodium = infer_ion_pair(neutral + PROTON, neutral + 22.989218, "POS")
    assert sodium is not None
    assert sodium["left_identity"] == "[M+H]+"
    assert sodium["right_identity"] == "[M+Na]+"
    assert sodium["relationship"] == "adduct"

    water_loss = infer_ion_pair(
        neutral + PROTON, neutral + PROTON - WATER, "POS"
    )
    assert water_loss is not None
    assert water_loss["right_identity"] == "[M+H-H2O]+"
    assert water_loss["relationship"] == "in_source_fragment"


def test_peak_shape_requires_correlation_overlap_and_two_sided_apex():
    accepted = chromatographic_shape_evidence(_trace(), _trace(scale=3.0))
    assert accepted["accepted"] is True
    assert accepted["correlation"] == 1.0
    rejected = chromatographic_shape_evidence(_trace(), _trace(shift=3))
    assert rejected["accepted"] is False


def test_network_collapses_adduct_and_source_fragment_and_links_unmatched_rows():
    neutral = 300.0
    rows = [
        _row("PF000001", neutral + PROTON, _match()),
        _row("PF000002", neutral + 22.989218),
        _row("PF000003", neutral + PROTON - WATER),
    ]
    summary = annotate_ion_identity_networks(
        rows,
        traces={
            "PF000001": _trace(),
            "PF000002": _trace(scale=0.5),
            "PF000003": _trace(scale=0.2),
        },
    )
    assert summary["network_count"] == 1
    assert summary["networked_features"] == 3
    assert summary["adduct_edges"] >= 1
    assert summary["in_source_fragment_edges"] >= 1
    assert summary["propagated_annotations"] == 2
    assert {row["ion_identity_network_id"] for row in rows} == {"IIN000001"}
    assert rows[1]["ion_identity_linked_annotation"]["status"] == (
        "linked_not_independent_identification"
    )
    assert rows[1]["match"] is None


def test_mass_relation_without_correlated_shape_is_not_networked():
    neutral = 300.0
    rows = [
        _row("PF000001", neutral + PROTON, _match()),
        _row("PF000002", neutral + 22.989218),
    ]
    bad_shape = {"rt": _trace()["rt"], "intensity": [100, 1, 80, 2, 60, 3, 40, 4, 20]}
    summary = annotate_ion_identity_networks(
        rows,
        traces={"PF000001": _trace(), "PF000002": bad_shape},
    )
    assert summary["network_count"] == 0
    assert rows[0]["ion_identity_status"] == "unlinked"


def test_no_ms1_placeholder_traces_are_not_reported_as_success(monkeypatch):
    rows = [_row("PF000001", 301.007276, _match())]

    def no_ms1_traces(_file_path, _features):
        return {
            "PF000001": {
                "rt": [],
                "intensity": [],
                "status": "insufficient_trace",
            },
        }, {
            "status": "no_ms1_scans",
            "ms1_scan_count": 0,
            "target_count": 1,
            "usable_trace_count": 0,
        }

    monkeypatch.setattr(
        ion_identity_module, "extract_targeted_ms1_traces", no_ms1_traces
    )
    summary = annotate_ion_identity_networks(
        rows, file_path="ms2_only.mzML"
    )

    assert summary["status"] == "no_ms1_scans"
    assert summary["usable_trace_count"] == 0
    assert summary["network_count"] == 0
    assert rows[0]["ion_identity_status"] == "no_ms1_scans"


def test_conflicting_compound_identifications_are_flagged_not_propagated():
    neutral = 300.0
    rows = [
        _row("PF000001", neutral + PROTON, _match()),
        _row(
            "PF000002",
            neutral + 22.989218,
            _match("Different", "CCCCCCCCCCCCCC-DDDDDDDDDD-E"),
        ),
    ]
    summary = annotate_ion_identity_networks(
        rows,
        traces={"PF000001": _trace(), "PF000002": _trace(scale=2)},
    )
    assert summary["identity_conflict_networks"] == 1
    assert rows[0]["match"]["ion_identity_conflict"] is True
    assert rows[1]["match"]["ion_identity_conflict"] is True
    assert summary["propagated_annotations"] == 0
