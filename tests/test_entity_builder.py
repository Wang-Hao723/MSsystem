from __future__ import annotations

import json

import pytest

from pipeline.entity_builder import (
    build_entity_bundle,
    build_entity_graph,
    read_entity_bundle,
    write_entity_bundle,
)


def _row(scan, *, mz=301.0712, rt=5.0, ce=20.0, feature="F1", ion_group="IG1"):
    return {
        "spectrum_id": str(scan),
        "precursor_mz": mz,
        "retention_time": rt,
        "ion_mode": "NEG",
        "collision_energy": ce,
        "precursor_feature_id": "PF1",
        "lcms_feature_id": feature,
        "lcms_feature_mz": 301.0712,
        "lcms_feature_rt_minutes": 5.0,
        "lcms_feature_rt_start_minutes": 4.9,
        "lcms_feature_rt_end_minutes": 5.1,
        "ion_identity_network_id": ion_group,
        "ion_identity_neutral_mass": 302.0785,
        "lcms_adduct": "[M-H]-",
        "match": {
            "compound_name": "Candidate A",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "formula": "C15H10O7",
            "cosine": 0.9,
            "final_score": 0.7,
            "rank_score": 0.7,
            "source_tier": "experimental",
            "top_candidates": [],
        },
    }


def test_entity_graph_is_independent_of_candidate_identity_and_score():
    rows = [_row(1, ce=20.0), _row(2, ce=40.0), _row(3, ce=60.0)]
    before = build_entity_graph("task", rows)
    changed = json.loads(json.dumps(rows))
    for row in changed:
        row["match"]["compound_name"] = "Completely different"
        row["match"]["inchikey"] = "BBBBBBBBBBBBBB-UHFFFAOYSA-N"
        row["match"]["formula"] = "C99H99"
        row["match"]["rank_score"] = 0.01
    after = build_entity_graph("task", changed)

    assert before.graph_sha256 == after.graph_sha256
    assert before.model_dump() == after.model_dump()


def test_scan_event_feature_ion_and_neutral_hierarchy_is_explicit():
    graph = build_entity_graph(
        "task", [_row(1, ce=20.0), _row(2, ce=40.0), _row(3, ce=60.0)]
    )
    counts = graph.summary["entity_counts"]

    assert counts["ms2_scan"] == 3
    assert counts["precursor_event"] == 1
    assert counts["chromatographic_feature"] == 1
    assert counts["ion_identity_group"] == 1
    assert counts["neutral_analyte"] == 1
    assert graph.summary["same_injection_extra_scans_collapsed"] == 2
    assert all(edge.candidate_independent for edge in graph.edges)
    neutral = next(node for node in graph.nodes if node.level.value == "neutral_analyte")
    assert neutral.attributes["identity_claimed"] is False


def test_fallback_complete_link_is_order_invariant_and_does_not_chain():
    def fallback(scan, mz, rt):
        row = _row(scan, mz=mz, rt=rt, feature=None, ion_group=None)
        row["precursor_feature_id"] = None
        return row

    rows = [
        fallback("a", 100.0000, 1.00),
        fallback("b", 100.0005, 1.02),
        fallback("c", 100.0010, 1.04),
    ]
    forward = build_entity_graph("task", rows)
    reverse = build_entity_graph("task", list(reversed(rows)))

    assert forward.graph_sha256 == reverse.graph_sha256
    assert forward.summary["precursor_event_count"] == 2


def test_entity_bundle_keeps_candidate_consensus_shadow_only(tmp_path):
    rows = [_row(1), _row(2, ce=40.0)]
    path = tmp_path / "bundle.json"
    summary = write_entity_bundle("task", rows, path)
    payload = read_entity_bundle(path, node_limit=2, consensus_limit=1)

    assert summary["feature_consensus_runtime_influence"] is False
    assert summary["feature_consensus_ranking_changed"] == 0
    assert payload["invariants"]["feature_consensus_changes_production_top1"] is False
    assert payload["nodes_total"] == 6
    assert payload["feature_consensus_total"] == 1


def test_real_mix33_entity_graph_has_no_runtime_ranking_change():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    path = root / "output" / "tasks" / "4515e871f701_results.json"
    if not path.exists():
        pytest.skip("V9 golden MIX33 result is unavailable")
    rows = json.loads(path.read_text(encoding="utf-8"))
    bundle = build_entity_bundle("4515e871f701", rows)

    assert bundle["summary"]["scan_count"] == 4189
    assert bundle["summary"]["precursor_event_count"] == 1397
    assert bundle["summary"]["feature_consensus_ranking_changed"] == 0
    assert bundle["invariants"]["entity_construction_uses_candidate_identity"] is False
