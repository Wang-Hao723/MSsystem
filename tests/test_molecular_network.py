"""Regression tests for the V7 evidence-aware molecular-network module."""

from pathlib import Path

import numpy as np
import scipy.sparse as sp

from matching.molecular_network import (
    build_network,
    component_size_filter,
    connected_components,
    pairwise_scores,
    to_edges_csv,
    to_nodes_csv,
    topk_filter,
)


ROOT = Path(__file__).resolve().parent.parent


def test_precursor_window_uses_sorted_mass_axis():
    """Unsorted scan order must not break precursor-window filtering."""
    spectra = sp.csr_matrix(np.ones((3, 3), dtype=np.float32))
    spectra = sp.diags(1.0 / np.sqrt(np.asarray(spectra.multiply(spectra).sum(axis=1)).ravel())) @ spectra
    occupied = sp.csr_matrix(np.ones((3, 3), dtype=np.float32))
    precursor_mz = np.array([500.0, 100.0, 510.0])

    edges = pairwise_scores(
        spectra.tocsr(), occupied, precursor_mz,
        cos_thresh=0.9, min_peaks=3, precursor_window=20.0,
    )

    assert {(i, j) for i, j, _, _ in edges} == {(0, 2)}


def test_public_edge_ids_survive_low_information_singletons():
    """Low-information spectra stay visible without changing public edge IDs."""
    shared_peaks = [[75.0, 40.0], [100.0, 100.0], [125.0, 70.0], [150.0, 55.0]]
    results = [
        {"precursor_mz": 400.0, "query_peaks": [[100.0, 1.0]], "match": None},
        {"precursor_mz": 500.0, "query_peaks": shared_peaks, "match": None},
        {"precursor_mz": 505.0, "query_peaks": shared_peaks, "match": None},
    ]

    network = build_network(
        results, cos_thresh=0.8, min_peaks=3, topk=10,
        mz_tol=0.1, precursor_window=20.0,
    )

    assert [node["idx"] for node in network["nodes"]] == [0, 1, 2]
    assert [(edge[0], edge[1]) for edge in network["edges"]] == [(1, 2)]
    assert [node["degree"] for node in network["nodes"]] == [0, 1, 1]
    assert (
        network["nodes"][0]["network_eligibility"]
        == "singleton_insufficient_fragments"
    )


def test_network_prefers_full_raw_peaks_over_top15_search_view():
    shared_raw = [[50.0 + i, 100.0 - i] for i in range(8)]
    results = [
        {
            "precursor_mz": 400.0 + i,
            "query_peaks": [[100.0 + i, 100.0], [150.0 + i, 50.0], [200.0 + i, 25.0]],
            "raw_peaks": shared_raw,
            "match": None,
        }
        for i in range(2)
    ]

    network = build_network(
        results, cos_thresh=0.9, min_peaks=6, topk=10,
        mz_tol=0.1, precursor_window=20.0,
    )

    assert len(network["edges"]) == 1
    assert all(node["peaks_count"] == 8 for node in network["nodes"])


def test_singletons_are_not_reported_as_molecular_families():
    results = [
        {
            "precursor_mz": 300.0 + index * 100,
            "raw_peaks": [
                [50.0 + index * 7 + peak, 100.0 - peak]
                for peak in range(6)
            ],
            "match": None,
        }
        for index in range(3)
    ]
    network = build_network(
        results,
        cos_thresh=0.95,
        min_peaks=6,
        precursor_window=20.0,
    )

    assert network["stats"]["n_molecular_families"] == 0
    assert network["stats"]["singleton_count"] == 3
    assert network["families"] == []
    assert all(node["family_id"] is None for node in network["nodes"])


def test_repeated_precursor_event_becomes_one_consensus_node():
    shared = [[75.0 + index * 10, 100.0 - index * 5] for index in range(7)]
    results = [
        {
            "precursor_feature_id": "PF000001",
            "spectrum_id": "scan-1",
            "precursor_mz": 500.0,
            "retention_time": 2.0,
            "raw_peaks": shared,
            "match": None,
        },
        {
            "precursor_feature_id": "PF000001",
            "spectrum_id": "scan-2",
            "precursor_mz": 500.0004,
            "retention_time": 2.01,
            "raw_peaks": shared,
            "match": None,
        },
        {
            "precursor_feature_id": "PF000002",
            "spectrum_id": "scan-3",
            "precursor_mz": 505.0,
            "retention_time": 2.2,
            "raw_peaks": shared,
            "match": None,
        },
    ]

    network = build_network(results, cos_thresh=0.8, min_peaks=4)

    assert network["stats"]["n_nodes"] == 2
    assert network["stats"]["consensus_nodes"] == 1
    assert network["stats"]["collapsed_spectra"] == 1
    assert network["nodes"][0]["spectrum_count"] == 2
    assert network["nodes"][0]["member_result_indices"] == [0, 1]


def test_precursor_shifted_analog_is_found_by_exact_modified_cosine():
    precursor_a = 400.0
    precursor_b = 414.0
    peaks_a = [[60.0, 100], [80.0, 90], [100.0, 80], [120.0, 70], [140.0, 60], [160.0, 50]]
    peaks_b = [[60.0, 100], [80.0, 90], [100.0, 80], [134.0, 70], [154.0, 60], [174.0, 50]]
    network = build_network(
        [
            {"precursor_mz": precursor_a, "raw_peaks": peaks_a, "match": None},
            {"precursor_mz": precursor_b, "raw_peaks": peaks_b, "match": None},
        ],
        cos_thresh=0.80,
        min_peaks=6,
        mz_tol=0.05,
    )

    assert len(network["edge_details"]) == 1
    edge = network["edge_details"][0]
    assert edge["edge_type"] == "analog"
    assert edge["modified_cosine"] >= 0.80
    assert edge["direct_cosine"] < edge["modified_cosine"]


def test_validated_iimn_is_a_separate_network_evidence_layer():
    results = [
        {
            "precursor_feature_id": f"PF00000{index + 1}",
            "precursor_mz": 300.0 + index * 22.0,
            "raw_peaks": [[60.0 + index * 20 + peak * 5, 100 - peak] for peak in range(6)],
            "ion_identity_status": "networked",
            "ion_identity_network_id": "IIN000001",
            "ion_identity": "[M-H]-" if index == 0 else "[M+Cl]-",
            "match": None,
        }
        for index in range(2)
    ]
    network = build_network(
        results,
        cos_thresh=0.99,
        min_peaks=6,
        precursor_window=1.0,
    )

    assert network["stats"]["ion_identity_edges"] == 1
    assert network["stats"]["n_molecular_families"] == 1
    assert network["edge_details"][0]["edge_type"] == "ion_identity"


def test_node_and_edge_csv_exports_are_joinable():
    shared = [[75.0, 40.0], [100.0, 100.0], [125.0, 70.0], [150.0, 55.0],
              [175.0, 35.0], [200.0, 25.0]]
    results = [
        {"spectrum_id": "scan-1", "precursor_mz": 500.0, "raw_peaks": shared, "match": None},
        {"spectrum_id": "scan-2", "precursor_mz": 505.0, "raw_peaks": shared, "match": None},
    ]
    network = build_network(results, cos_thresh=0.8, min_peaks=3, max_component_size=0)

    nodes_csv = to_nodes_csv(network)
    edges_csv = to_edges_csv(network)

    assert "id,spectrum_id,precursor_mz" in nodes_csv
    assert "0,scan-1,500.0" in nodes_csv
    assert "source,target,cosine,matched_peaks,family_id" in edges_csv
    assert "0,1,1.0,6,0" in edges_csv


def test_topk_filter_requires_mutual_neighbours():
    edges = [
        (0, 1, 0.90, 8),
        (0, 2, 0.80, 8),
        (1, 2, 0.70, 8),
    ]
    assert topk_filter(edges, n_nodes=3, topk=1) == [(0, 1, 0.90, 8)]


def test_component_size_filter_breaks_hairballs():
    edges = [
        (0, 1, 0.99, 9),
        (1, 2, 0.98, 9),
        (2, 3, 0.97, 9),
        (3, 4, 0.96, 9),
    ]
    filtered = component_size_filter(edges, n_nodes=5, max_component_size=2)
    _, families = connected_components(5, filtered)
    assert max(map(len, families.values())) <= 2
    assert (0, 1, 0.99, 9) in filtered


def test_network_page_and_api_are_mounted():
    html = (ROOT / "web/static/index.html").read_text(encoding="utf-8")
    assert 'data-page="network"' in html
    assert 'id="page-network"' in html
    assert '/static/js/network.js?v=20261022-core-evidence-repair-v1' in html

    from web.server import app
    paths = {route.path for route in app.routes}
    assert "/api/tasks/{task_id}/network" in paths
    assert "/api/tasks/{task_id}/network/node/{node_id}" in paths
    assert "/api/tasks/{task_id}/network/graphml" in paths
    assert "/api/tasks/{task_id}/network/nodes.csv" in paths
    assert "/api/tasks/{task_id}/network/edges.csv" in paths
