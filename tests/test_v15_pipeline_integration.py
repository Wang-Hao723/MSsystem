from __future__ import annotations

import json
import copy
import pytest
import matching.v15_flagship_ranker as policy_loader
import pipeline.v15_flagship as flagship_pipeline

from matching.v15_candidate_recall import (
    load_v15_candidate_recall_policy,
    semantic_sha256,
)
from matching.flagship_ranker import flagship_policy_sha256
from pipeline.v15_flagship import (
    _entity_semantic_sha256,
    artifact_is_current,
    upgrade_v15_flagship_identity_scope_artifact,
    write_v15_flagship_artifact,
    v15_flagship_artifact_path,
)


@pytest.fixture
def legacy_scope_upgrade_policy(monkeypatch):
    """Exercise the retained pre-namespace identity-only migration contract."""
    policy = copy.deepcopy(policy_loader.load_v15_flagship_policy())
    policy.pop('artifact_namespace', None)
    for module in (policy_loader, flagship_pipeline):
        monkeypatch.setattr(module, 'load_v15_flagship_policy', lambda: copy.deepcopy(policy))
    return policy


def _candidate() -> dict:
    return {
        "compound_name": "Primary candidate",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "cosine": 0.95,
        "entropy": 0.95,
        "reverse_match": 0.95,
        "fragment_score": 0.95,
        "intensity_score": 0.95,
        "mass_ppm": 0.0,
        "final_score": 0.95,
        "rank_score": 0.95,
        "dataset_source": "synthetic-test",
        "source_tier": "experimental",
    }


def _rows() -> list[dict]:
    rows = []
    for ce in (20.0, 30.0):
        candidate = _candidate()
        candidate["top_candidates"] = [_candidate()]
        rows.append(
            {
                "spectrum_id": f"scan-{ce:g}",
                "precursor_feature_id": "feature-1",
                "precursor_mz": 195.0,
                "search_precursor_mz": 195.0,
                "ion_mode": "POS",
                "collision_energy": ce,
                "collision_energy_unit": "eV",
                "raw_peaks": [[80.0, 10.0], [100.0, 100.0], [140.0, 30.0]],
                "match": candidate,
            }
        )
    return rows


class _RecallRetriever:
    library_snapshot_id = "synthetic-library-snapshot"

    def __init__(self) -> None:
        self.policy = load_v15_candidate_recall_policy()
        self.policy_sha256 = semantic_sha256(self.policy)
        self.calls = 0

    def retrieve_entity(
        self,
        entity_id,
        scans,
        *,
        existing_candidates,
        detection_method_id,
    ):
        self.calls += 1
        return {
            "schema_version": "2.0.0",
            "module_id": "v12_evirank_candidate_recall",
            "algorithm_id": self.policy["algorithm_id"],
            "entity_id": entity_id,
            "detection_method_id": detection_method_id,
            "library_snapshot_id": self.library_snapshot_id,
            "policy_sha256": self.policy_sha256,
            "existing_connectivity_count": len(existing_candidates),
            "target_total_connectivities": 40,
            "expanded_candidate_count": 1,
            "expanded_candidates": [
                {
                    "identity_key": "std2d:CCCCCCCCCCCCCC",
                    "connectivity_identity": "CCCCCCCCCCCCCC",
                    "candidate": {
                        "compound_name": "Expanded hypothesis",
                        "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDD-E",
                    },
                    "retrieval_rank": 1,
                    "retrieval_rrf_utility": 0.031,
                    "distinct_condition_support": 2,
                    "spectral_view_support": ["library_match", "sensitive"],
                }
            ],
            "source_paired_consensus": {"status": "not_applicable"},
            "scan_audit": [],
        }


def test_v15_artifact_materializes_top40_support_without_second_library_search(
    tmp_path,
):
    retriever = _RecallRetriever()
    destination = v15_flagship_artifact_path(tmp_path, "task")

    summary = write_v15_flagship_artifact(
        "task",
        _rows(),
        destination,
        ce_evidence_retriever=retriever,
        candidate_recall_retriever=retriever,
        detection_method_id="qtof-dda",
    )
    payload = json.loads(destination.read_text(encoding="utf-8"))
    entity = payload["entities"][0]

    assert summary["candidate_recall_support_status"] == "evaluated"
    assert retriever.calls == 1
    assert entity["winner"]["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    assert entity["candidate_recall_support"]["expanded_candidates"][0][
        "can_replace_winner"
    ] is False
    assert payload["invariants"][
        "candidate_recall_support_affects_primary_ranking"
    ] is False
    assert payload["metadata"]["decision_assurance_status"] == "available"
    assert entity["decision_assurance"]["candidate_recall"][
        "expanded_candidate_count"
    ] == 1
    calibration = entity["decision_assurance"][
        "project_target_decoy_calibration"
    ]
    assert calibration["status"] == "unavailable"
    assert calibration["q_value"] is None


def test_compatible_v15_artifact_upgrades_scope_without_reranking(tmp_path, legacy_scope_upgrade_policy):
    destination = tmp_path / "upgrade_v15_flagship.json"
    write_v15_flagship_artifact("upgrade", _rows(), destination)
    payload = json.loads(destination.read_text(encoding="utf-8"))
    ranking_before = [
        [
            (
                candidate.get("identity_key"),
                candidate.get("entity_rank"),
                candidate.get("entity_score"),
            )
            for candidate in entity.get("candidates") or []
        ]
        for entity in payload["entities"]
    ]

    payload["policy"].pop("identity_resolution")
    payload["policy"]["algorithm_id"] = (
        "mssystem-v15-conditioned-multilibrary-1.2.0"
    )
    for entity in payload["entities"]:
        resolution = entity.pop("identity_resolution")
        decision = entity["decision"]
        decision["status"] = resolution["upstream_decision_status"]
        decision["scope"] = resolution["upstream_decision_scope"]
        for key in (
            "reportable_identity_key",
            "reportable_formula",
            "ranking_lead_identity_key",
            "ranking_lead_is_reportable_single_connectivity",
            "identity_resolution_module_id",
            "identity_resolution_algorithm_id",
        ):
            decision.pop(key, None)
        selected = entity.get("winner") or entity.get(
            "leading_candidate_when_abstained"
        )
        if isinstance(selected, dict):
            selected.pop("reporting_semantics", None)
    metadata = payload["metadata"]
    for key in list(metadata):
        if key.startswith("identity_resolution") or key.startswith(
            "unresolved_"
        ):
            metadata.pop(key)
    metadata["algorithm_id"] = payload["policy"]["algorithm_id"]
    metadata["policy_sha256"] = flagship_policy_sha256(payload["policy"])
    metadata["entity_semantic_sha256"] = _entity_semantic_sha256(payload)
    destination.write_text(json.dumps(payload), encoding="utf-8")

    assert artifact_is_current(destination) is False
    summary = upgrade_v15_flagship_identity_scope_artifact(
        destination,
        task_id="upgrade",
        runtime_snapshot_id="runtime-upgrade-test",
    )
    upgraded = json.loads(destination.read_text(encoding="utf-8"))
    ranking_after = [
        [
            (
                candidate.get("identity_key"),
                candidate.get("entity_rank"),
                candidate.get("entity_score"),
            )
            for candidate in entity.get("candidates") or []
        ]
        for entity in upgraded["entities"]
    ]

    assert artifact_is_current(destination) is True
    assert summary["upgrade_status"] == "upgraded_without_library_retrieval"
    assert ranking_after == ranking_before
    assert upgraded["metadata"]["identity_resolution_upgrade_changes_winner"] is False
    assert upgraded["metadata"]["runtime_snapshot_id"] == "runtime-upgrade-test"
    assert all(entity.get("identity_resolution") for entity in upgraded["entities"])


def test_v15_scope_upgrade_rejects_unverified_entity_digest(tmp_path, legacy_scope_upgrade_policy):
    destination = tmp_path / "corrupt_v15_flagship.json"
    write_v15_flagship_artifact("corrupt", _rows(), destination)
    payload = json.loads(destination.read_text(encoding="utf-8"))
    payload["policy"].pop("identity_resolution")
    payload["policy"]["algorithm_id"] = (
        "mssystem-v15-conditioned-multilibrary-1.2.0"
    )
    payload["metadata"]["algorithm_id"] = payload["policy"]["algorithm_id"]
    payload["metadata"]["policy_sha256"] = flagship_policy_sha256(
        payload["policy"]
    )
    payload["metadata"]["identity_resolution_status"] = None
    payload["metadata"]["entity_semantic_sha256"] = "0" * 64
    destination.write_text(json.dumps(payload), encoding="utf-8")

    try:
        upgrade_v15_flagship_identity_scope_artifact(
            destination,
            task_id="corrupt",
        )
    except ValueError as exc:
        assert "entity semantic digest" in str(exc)
    else:
        raise AssertionError("A corrupt V15 artifact must not be upgraded")
