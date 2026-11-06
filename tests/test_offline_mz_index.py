from __future__ import annotations

import json

import pytest

from matching.offline_mz_index import (
    build_library_mass_indices,
    load_library_mass_indices,
    record_alignment_sha256,
)
from scripts.build_library_mmap_index import _snapshot_id
from matching.search import _neutral_candidates, _query_candidates


def _bucket_index(records, bucket_size=0.1, field="precursor_mz"):
    output = {}
    for index, record in enumerate(records):
        mass = record.get(field)
        if mass is None:
            continue
        output.setdefault(int(mass / bucket_size), []).append(index)
    for indices in output.values():
        indices.sort(key=lambda index: records[index][field])
    return output


def test_mmap_precursor_index_is_candidate_order_equivalent(tmp_path):
    records = [
        {
            "spectrum_id": f"s{index}",
            "precursor_mz": mass,
            "neutral_mass": mass - 1.0073,
        }
        for index, mass in enumerate(
            [301.000, 300.991, 301.009, 100.0, 301.021, 300.980, 500.0]
        )
    ]
    manifest = build_library_mass_indices(
        records,
        tmp_path,
        library_snapshot_id="snapshot-test",
        complete_library=True,
    )
    mmap_mz, mmap_neutral, loaded = load_library_mass_indices(
        tmp_path,
        expected_library_snapshot_id="snapshot-test",
        expected_record_count=len(records),
        expected_record_alignment_sha256=record_alignment_sha256(records),
    )
    bucket_mz = _bucket_index(records)
    bucket_neutral = _bucket_index(records, field="neutral_mass")

    for query in (100.0, 300.99, 301.0, 301.01, 500.0):
        assert _query_candidates(
            query, records, mmap_mz, tol_da=0.02, tol_ppm=10
        ) == _query_candidates(
            query, records, bucket_mz, tol_da=0.02, tol_ppm=10
        )
    for query in (299.98, 299.99, 300.0):
        assert _neutral_candidates(
            query, records, mmap_neutral, tol_da=0.02, tol_ppm=10
        ) == _neutral_candidates(
            query, records, bucket_neutral, tol_da=0.02, tol_ppm=10
        )
    assert manifest == loaded
    assert mmap_mz.mapped_bytes > 0


def test_mmap_index_refuses_alignment_or_partial_cache(tmp_path):
    records = [{"spectrum_id": "s1", "precursor_mz": 100.0}]
    build_library_mass_indices(
        records,
        tmp_path,
        library_snapshot_id="snapshot-test",
        complete_library=False,
    )
    with pytest.raises(ValueError, match="Partial benchmark"):
        load_library_mass_indices(tmp_path)
    with pytest.raises(ValueError, match="alignment"):
        load_library_mass_indices(
            tmp_path,
            expected_record_alignment_sha256="bad",
            require_complete_library=False,
        )


def test_retrieval_policy_requires_equivalence_before_cutover():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    policy = json.loads(
        (root / "config" / "retrieval_policy.json").read_text(encoding="utf-8")
    )
    gate = policy["offline_mmap_cutover"]
    assert gate["default_enabled"] is False
    assert gate["requires_record_alignment_hash_match"] is True
    assert gate["requires_golden_topk_equivalence"] is True


def test_full_index_builder_uses_production_library_snapshot_id():
    assert _snapshot_id().startswith("sqlite-")
