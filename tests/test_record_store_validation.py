import copy
import json
import os

import pytest

from eval.v16_evidence_validation import record_store_validation as validation


def test_scoring_projection_removes_only_storage_method_flag():
    score = {"compound_name": "Example", "final_score": .91, "confidence": "High",
        "top_k": [{"compound_name": "Other", "score": .90},
                  {"compound_name": "Third", "score": .89}],
        "retrieval_audit": {"mmap_mass_index_used": True, "candidates_scored": 3}}
    original = copy.deepcopy(score)
    changed = validation.score_projection(score)
    assert score == original
    assert changed == {**score, "retrieval_audit": {"candidates_scored": 3}}
    other = copy.deepcopy(score)
    other["final_score"] += .001
    assert validation.score_projection(other) != changed
    other = copy.deepcopy(score)
    other["top_k"].reverse()
    assert validation.score_projection(other) != changed
    assert validation.score_projection(None) is None


def test_query_argument_projection_is_complete_and_unambiguous():
    assert len(validation.QUERY_ARGS) == len(validation.QUERY_FIELDS)
    assert len(set(validation.QUERY_ARGS)) == len(validation.QUERY_ARGS)
    assert "match" not in validation.QUERY_FIELDS
    assert "truth_connectivity" not in validation.QUERY_FIELDS


def test_per_query_baseline_records_do_not_force_a_full_preload():
    class Backing:
        def __init__(self):
            self.read = []
        def __len__(self):
            return 2000000
        def read_uncached(self, i):
            self.read.append(i)
            return {"index": i}
    backing = Backing()
    selected = validation.QueryLocalRecords(backing)
    assert len(selected) == 2000000
    assert selected[1500000] == {"index": 1500000}
    assert selected[1500000] == {"index": 1500000}
    assert backing.read == [1500000]


def test_protocol_distinguishes_partial_memory_from_full_startup():
    protocol = json.loads(validation.PROTOCOL.read_text())
    assert protocol["identity_accuracy_claim_allowed"] is False
    assert protocol["candidate_order"]["allowed_mismatches"] == 0
    assert protocol["resources"]["partial_memory_comparison_must_be_labelled_partial"] is True
    assert protocol["resources"]["os_cache_eviction_or_other_process_termination_allowed"] is False


def acceptance_fixture(tmp_path):
    manifest = {"semantic_sha256": "cache-hash", "record_count": 1234,
                "normalized_records_sha256": "record-hash"}
    payload = {"schema_version": "mssystem.v16.library-record-store-validation/1.0.0",
        "manifest_semantic_sha256": manifest["semantic_sha256"],
        "protocol_sha256": validation.file_sha256(validation.PROTOCOL),
        "roundtrip": {"verified_record_count": 1234, "full_record_roundtrip_sha256": "record-hash"},
        "retrieval": {"candidate_order_cases": 10384, "candidate_order_mismatches": 0},
        "scoring": {"paired_scoring_queries": 138, "score_mismatches": 0},
        "shared_database_unchanged": True}
    payload["semantic_sha256"] = validation.semantic_sha256(payload)
    output = tmp_path / "validation.json"
    output.write_text(json.dumps(payload))
    return manifest, payload, output


def test_acceptance_publication_supports_filesystems_without_hardlinks(tmp_path, monkeypatch):
    manifest, payload, output = acceptance_fixture(tmp_path)
    def no_hardlinks(*args, **kwargs):
        raise OSError(45, "Operation not supported")
    monkeypatch.setattr(os, "link", no_hardlinks)
    assert validation.publish_acceptance(payload, output, tmp_path, manifest)
    saved = (tmp_path / "acceptance.json").read_bytes()
    assert validation.publish_acceptance(payload, output, tmp_path, manifest) is False
    assert (tmp_path / "acceptance.json").read_bytes() == saved


@pytest.mark.parametrize("change", ["scores", "record_count", "protocol", "receipt_hash"])
def test_acceptance_publication_does_not_reinterpret_failed_validation(tmp_path, change):
    manifest, payload, output = acceptance_fixture(tmp_path)
    if change == "scores":
        payload["scoring"]["score_mismatches"] = 1
    elif change == "record_count":
        payload["roundtrip"]["verified_record_count"] -= 1
    elif change == "protocol":
        payload["protocol_sha256"] = "wrong"
    if change != "receipt_hash":
        payload["semantic_sha256"] = validation.semantic_sha256(
            {key: value for key, value in payload.items() if key != "semantic_sha256"})
    else:
        payload["semantic_sha256"] = "wrong"
    output.write_text(json.dumps(payload))
    with pytest.raises(ValueError):
        validation.publish_acceptance(payload, output, tmp_path, manifest)
    assert not (tmp_path / "acceptance.json").exists()
