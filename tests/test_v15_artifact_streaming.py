from __future__ import annotations

import json
from pathlib import Path

from matching.flagship_ranker import flagship_policy_sha256
from matching.v15_flagship_ranker import load_v15_flagship_policy
from pipeline import v15_flagship
from pipeline.flagship import read_flagship_artifact


def test_v15_current_check_streams_metadata_without_full_json_load(
    tmp_path, monkeypatch
):
    policy = load_v15_flagship_policy()
    artifact = v15_flagship.v15_flagship_artifact_path(tmp_path, 'large')
    v15_flagship.write_v15_flagship_artifact('large', [], artifact)
    verified_metadata = json.loads(artifact.read_text())['metadata']
    artifact.write_text(
        json.dumps(
            {
                "metadata": {
                    **verified_metadata,
                    "algorithm_id": policy["algorithm_id"],
                    "policy_sha256": flagship_policy_sha256(policy),
                    "identity_resolution_status": "available",
                    "entity_count": 2000,
                },
                "invariants": {"legacy_result_mutated": False},
                "entities": [
                    {"entity_id": f"E{index}", "padding": "x" * 1000}
                    for index in range(2000)
                ],
            },
            separators=(",", ":"),
        ),
        encoding="utf-8",
    )
    v15_flagship._read_artifact_metadata_cached.cache_clear()

    original_read_text = Path.read_text

    def guarded_read_text(path, *args, **kwargs):
        if path.resolve() == artifact.resolve():
            raise AssertionError("artifact_is_current loaded the complete artifact")
        return original_read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", guarded_read_text)

    assert v15_flagship.artifact_is_current(artifact) is True


def test_unfiltered_flagship_page_uses_metadata_total(tmp_path):
    artifact = tmp_path / "v11_flagship.json"
    artifact.write_text(
        json.dumps(
            {
                "metadata": {"entity_count": 100},
                "invariants": {},
                "entities": [
                    {"entity_id": f"E{index}"} for index in range(100)
                ],
            }
        ),
        encoding="utf-8",
    )

    page = read_flagship_artifact(artifact, offset=10, limit=2)

    assert [row["entity_id"] for row in page["entities"]] == ["E10", "E11"]
    assert page["entities_total"] == 100
    assert page["has_more"] is True
