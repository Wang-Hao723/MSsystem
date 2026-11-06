"""Peak hashing and conservative mirror-spectrum clustering.

Exact hashes are insensitive to input peak order and uniform intensity
scaling.  Near-mirror links are emitted only when acquisition metadata are
compatible and both cosine similarity and matched-peak coverage pass fixed
thresholds.  The thresholds are engineering defaults, not empirical claims;
the workflow labels such links as candidates until a labelled duplicate set
has been calibrated independently.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from typing import Any


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def exact_peak_hash(
    peaks: Sequence[Sequence[float]], *, mz_digits: int = 5, intensity_digits: int = 6
) -> str:
    """Return an order- and scale-invariant SHA-256 peak hash."""
    cleaned = [
        (float(peak[0]), float(peak[1]))
        for peak in peaks
        if len(peak) >= 2
        and math.isfinite(float(peak[0]))
        and math.isfinite(float(peak[1]))
        and float(peak[0]) > 0
        and float(peak[1]) > 0
    ]
    if not cleaned:
        raise ValueError("a peak hash requires at least one positive finite peak")
    base = max(intensity for _, intensity in cleaned)
    normalized = [
        [round(mz, mz_digits), round(intensity / base, intensity_digits)]
        for mz, intensity in sorted(cleaned)
    ]
    return hashlib.sha256(_canonical_bytes(normalized)).hexdigest()


def sparse_peak_signature(
    peaks: Sequence[Sequence[float]], *, top_n: int = 8, mz_bin: float = 0.05
) -> tuple[int, ...]:
    """Return a deterministic blocking signature from the strongest peaks."""
    ranked = sorted(
        ((float(p[0]), float(p[1])) for p in peaks if len(p) >= 2),
        key=lambda item: (-item[1], item[0]),
    )[:top_n]
    return tuple(sorted({int(round(mz / mz_bin)) for mz, _ in ranked}))


def greedy_peak_cosine(
    left: Sequence[Sequence[float]],
    right: Sequence[Sequence[float]],
    *,
    mz_tolerance: float = 0.01,
) -> tuple[float, float, float]:
    """Compute a high-resolution greedy cosine and bilateral peak coverage.

    Intensities are L2-normalized.  Candidate matches within ``mz_tolerance``
    are selected by decreasing intensity product, matching the central
    behaviour of greedy cosine implementations used for library spectra.
    """
    a = [(float(p[0]), float(p[1])) for p in left if len(p) >= 2 and p[1] > 0]
    b = [(float(p[0]), float(p[1])) for p in right if len(p) >= 2 and p[1] > 0]
    if not a or not b:
        return 0.0, 0.0, 0.0
    norm_a = math.sqrt(sum(intensity * intensity for _, intensity in a))
    norm_b = math.sqrt(sum(intensity * intensity for _, intensity in b))
    candidates: list[tuple[float, int, int]] = []
    for i, (mz_a, int_a) in enumerate(a):
        for j, (mz_b, int_b) in enumerate(b):
            delta = mz_b - mz_a
            if abs(delta) <= mz_tolerance:
                candidates.append((int_a * int_b, i, j))
    used_a: set[int] = set()
    used_b: set[int] = set()
    dot = 0.0
    for product, i, j in sorted(candidates, reverse=True):
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        dot += product
    return dot / (norm_a * norm_b), len(used_a) / len(a), len(used_b) / len(b)


class UnionFind:
    """Small deterministic union-find implementation."""

    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        while parent != self.parent[parent]:
            self.parent[parent] = self.parent[self.parent[parent]]
            parent = self.parent[parent]
        self.parent[value] = parent
        return parent

    def union(self, left: str, right: str) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        keep, merge = sorted((root_left, root_right))
        self.parent[merge] = keep


def cluster_mirror_records(
    records: Sequence[Mapping[str, Any]],
    *,
    allow_near_candidates: bool = False,
    mz_tolerance: float = 0.01,
    cosine_threshold: float = 0.995,
    coverage_threshold: float = 0.90,
) -> dict[str, Any]:
    """Cluster exact mirrors and optionally conservative near-mirror candidates."""
    ids = [str(record["spectrum_id"]) for record in records]
    union = UnionFind(ids)
    exact_buckets: dict[str, list[str]] = defaultdict(list)
    record_by_id = {str(record["spectrum_id"]): record for record in records}
    for record in records:
        exact_buckets[str(record["independence"]["exact_peak_hash"])].append(
            str(record["spectrum_id"])
        )
    exact_edges = 0
    for members in exact_buckets.values():
        anchor = members[0]
        for member in members[1:]:
            union.union(anchor, member)
            exact_edges += 1

    candidate_edges: list[dict[str, Any]] = []
    if allow_near_candidates:
        blocks: dict[tuple[Any, ...], list[str]] = defaultdict(list)
        for record in records:
            acquisition = record["acquisition"]
            compound = record["compound"]
            key = (
                compound["connectivity"],
                acquisition["polarity"],
                acquisition["adduct"],
                acquisition["instrument_domain_id"],
                acquisition["collision_energy"]["unit"],
                acquisition["collision_energy"]["mode"],
                sparse_peak_signature(record["peaks"]),
            )
            blocks[key].append(str(record["spectrum_id"]))
        for members in blocks.values():
            if len(members) < 2:
                continue
            for i, left_id in enumerate(members):
                for right_id in members[i + 1 :]:
                    left = record_by_id[left_id]
                    right = record_by_id[right_id]
                    cosine, coverage_left, coverage_right = greedy_peak_cosine(
                        left["peaks"], right["peaks"], mz_tolerance=mz_tolerance
                    )
                    if (
                        cosine >= cosine_threshold
                        and min(coverage_left, coverage_right) >= coverage_threshold
                    ):
                        union.union(left_id, right_id)
                        candidate_edges.append(
                            {
                                "left": left_id,
                                "right": right_id,
                                "cosine": cosine,
                                "coverage_left": coverage_left,
                                "coverage_right": coverage_right,
                                "status": "uncalibrated_near_mirror_candidate",
                            }
                        )

    components: dict[str, list[str]] = defaultdict(list)
    for spectrum_id in ids:
        components[union.find(spectrum_id)].append(spectrum_id)
    mapping: dict[str, str] = {}
    for members in components.values():
        digest = hashlib.sha256("\n".join(sorted(members)).encode("utf-8")).hexdigest()
        cluster_id = f"mirror:{digest[:20]}"
        for member in members:
            mapping[member] = cluster_id
    return {
        "mapping": mapping,
        "component_count": len(components),
        "exact_duplicate_edge_count": exact_edges,
        "near_candidate_edges": candidate_edges,
        "near_candidate_status": (
            "uncalibrated_candidates_included"
            if allow_near_candidates
            else "not_evaluated_exact_hash_only"
        ),
    }
