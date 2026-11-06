from __future__ import annotations

import pytest

from eval.evaluate_v12_source_paired_mix33 import (
    _raw_v11_join_audit,
    paired_exact_two_sided_p,
)


def test_raw_replay_joins_by_unique_spectrum_id_not_legacy_array_index():
    rows = [
        {"spectrum_id": "scan-b", "peaks": [[10.0, 1.0]]},
        {"spectrum_id": "scan-a", "peaks": [[20.0, 1.0]]},
    ]
    artifact = {
        "metadata": {"source_row_count": 2},
        "entities": [
            {
                "member_result_indices": [0, 1],
                "member_spectrum_ids": ["scan-a", "scan-b"],
            }
        ],
    }

    lookup, audit = _raw_v11_join_audit(rows, artifact)

    assert set(lookup) == {"scan-a", "scan-b"}
    assert lookup["scan-a"]["peaks"] == [[20.0, 1.0]]
    assert audit["join_key"] == "unique_spectrum_id"
    assert audit["array_index_used_for_raw_join"] is False
    assert audit["legacy_array_index_alignment_complete"] is False
    assert all(audit["checks"].values())


def test_raw_replay_refuses_nonunique_spectrum_ids():
    artifact = {
        "metadata": {"source_row_count": 2},
        "entities": [
            {
                "member_result_indices": [0, 1],
                "member_spectrum_ids": ["scan-a", "scan-b"],
            }
        ],
    }
    with pytest.raises(ValueError, match="unique and non-empty"):
        _raw_v11_join_audit(
            [{"spectrum_id": "scan-a"}, {"spectrum_id": "scan-a"}],
            artifact,
        )


def test_paired_exact_two_sided_p_uses_all_discordant_pairs():
    assert paired_exact_two_sided_p(3, 0) == pytest.approx(0.25)
    assert paired_exact_two_sided_p(2, 1) == pytest.approx(1.0)
    assert paired_exact_two_sided_p(0, 0) is None
    with pytest.raises(ValueError, match="cannot be negative"):
        paired_exact_two_sided_p(-1, 0)
