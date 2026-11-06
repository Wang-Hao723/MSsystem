from __future__ import annotations

import json
from pathlib import Path

import pytest

from matching.feature_fragment_benchmark.corpus import (
    RecordReject,
    iter_json_array,
    stable_compound_name,
)
from matching.feature_fragment_benchmark.mirror import (
    cluster_mirror_records,
    exact_peak_hash,
)


def test_iter_json_array_streams_ten_thousand_records_in_small_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "ten_thousand.json"
    with path.open("w", encoding="utf-8") as handle:
        handle.write("[")
        for index in range(10_000):
            if index:
                handle.write(",")
            json.dump(
                {"index": index, "payload": f"record-{index:05d}"},
                handle,
                separators=(",", ":"),
            )
        handle.write("]")

    def whole_file_read_forbidden(*args, **kwargs):  # noqa: ANN002, ANN003
        raise AssertionError("iter_json_array must not call Path.read_text")

    monkeypatch.setattr(Path, "read_text", whole_file_read_forbidden)
    count = 0
    checksum = 0
    for row in iter_json_array(path, chunk_size=17):
        count += 1
        checksum += row["index"]

    assert count == 10_000
    assert checksum == 49_995_000


def test_exact_peak_hash_and_mirror_cluster_ignore_order_and_uniform_scaling(
    benchmark_record_factory,
) -> None:
    ordered = [[50.0, 10.0], [100.0, 5.0], [200.0, 1.0]]
    reordered_scaled = [[200.0, 7.0], [50.0, 70.0], [100.0, 35.0]]
    assert exact_peak_hash(ordered) == exact_peak_hash(reordered_scaled)

    left = benchmark_record_factory("left", peaks=ordered)
    right = benchmark_record_factory("right", peaks=reordered_scaled)
    clustered = cluster_mirror_records([left, right])

    assert clustered["exact_duplicate_edge_count"] == 1
    assert clustered["component_count"] == 1
    assert clustered["mapping"]["left"] == clustered["mapping"]["right"]


@pytest.mark.parametrize(
    "composite_name",
    [
        "NCGC001234 C10H12O2 3,4-dihydroxybenzoic acid",
        "CCMSLIB000000001 C6H12O6 D-glucose",
        "HMDB0000123 C9H8O4 acetylsalicylic acid",
        "MSBNK-RIKEN-001 C8H10N4O2 caffeine",
    ],
)
def test_database_accession_composite_cannot_be_selected_as_compound_name(
    composite_name: str,
) -> None:
    with pytest.raises(RecordReject, match="stable_compound_name_missing"):
        stable_compound_name([composite_name])

    assert stable_compound_name([composite_name, "Caffeine"]) == ("Caffeine", 1)
