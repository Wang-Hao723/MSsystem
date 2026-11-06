"""Content-validated NumPy mmap indices for precursor and neutral masses.

The cache stores only mass-to-record-position mappings.  It is deliberately
separate from candidate scoring, so enabling it cannot change scores or Top-K
ordering when the record-alignment digest matches.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np


MMAP_INDEX_SCHEMA_VERSION = "1.0.0"


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0 else None


def _record_precursor_mass(record: dict) -> float | None:
    return _finite(record.get("search_precursor_mz")) or _finite(
        record.get("precursor_mz")
    )


def record_alignment_sha256(records: Sequence[dict]) -> str:
    digest = hashlib.sha256()
    for index, record in enumerate(records):
        payload = (
            index,
            _record_precursor_mass(record),
            _finite(record.get("neutral_mass")),
            str(record.get("spectrum_id") or ""),
        )
        digest.update(
            json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode(
                "utf-8"
            )
        )
        digest.update(b"\n")
    return digest.hexdigest()


class MMapMassIndex:
    """Read-only globally sorted mass index backed by `.npy` memory maps."""

    def __init__(self, masses_path: Path, record_indices_path: Path):
        self.masses_path = Path(masses_path)
        self.record_indices_path = Path(record_indices_path)
        self._masses = np.load(self.masses_path, mmap_mode="r", allow_pickle=False)
        self._record_indices = np.load(
            self.record_indices_path, mmap_mode="r", allow_pickle=False
        )
        if self._masses.ndim != 1 or self._record_indices.ndim != 1:
            raise ValueError("Mmap mass-index arrays must be one-dimensional")
        if len(self._masses) != len(self._record_indices):
            raise ValueError("Mmap mass and record-index arrays are misaligned")
        if len(self._masses) > 1 and np.any(self._masses[1:] < self._masses[:-1]):
            raise ValueError("Mmap mass array is not sorted")

    def __len__(self) -> int:
        return int(len(self._masses))

    def range_query(self, center: float, tolerance: float) -> list[int]:
        center = float(center)
        tolerance = float(tolerance)
        if not math.isfinite(center) or not math.isfinite(tolerance) or tolerance < 0:
            return []
        left = int(np.searchsorted(self._masses, center - tolerance, side="left"))
        right = int(np.searchsorted(self._masses, center + tolerance, side="right"))
        return [int(value) for value in self._record_indices[left:right]]

    @property
    def mapped_bytes(self) -> int:
        return self.masses_path.stat().st_size + self.record_indices_path.stat().st_size


def _write_npy_atomic(path: Path, array: np.ndarray) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("wb") as handle:
        np.save(handle, array, allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _sorted_mass_arrays(
    records: Sequence[dict],
    field: str,
) -> tuple[np.ndarray, np.ndarray]:
    masses = []
    indices = []
    for index, record in enumerate(records):
        mass = (
            _record_precursor_mass(record)
            if field == "precursor"
            else _finite(record.get("neutral_mass"))
        )
        if mass is None:
            continue
        masses.append(mass)
        indices.append(index)
    mass_array = np.asarray(masses, dtype=np.float64)
    index_dtype = np.uint32 if len(records) <= np.iinfo(np.uint32).max else np.uint64
    index_array = np.asarray(indices, dtype=index_dtype)
    if len(mass_array):
        order = np.lexsort((index_array, mass_array))
        mass_array = mass_array[order]
        index_array = index_array[order]
    return mass_array, index_array


def build_library_mass_indices(
    records: Sequence[dict],
    output_dir: Path,
    *,
    library_snapshot_id: str,
    source_db: Path | None = None,
    complete_library: bool = True,
) -> dict:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {}
    counts = {}
    for field in ("precursor", "neutral"):
        masses, record_indices = _sorted_mass_arrays(records, field)
        masses_path = output_dir / f"{field}_masses.npy"
        indices_path = output_dir / f"{field}_record_indices.npy"
        _write_npy_atomic(masses_path, masses)
        _write_npy_atomic(indices_path, record_indices)
        files[field] = {
            "masses": masses_path.name,
            "record_indices": indices_path.name,
        }
        counts[field] = len(masses)
    source_db = Path(source_db).resolve() if source_db else None
    manifest = {
        "schema_version": MMAP_INDEX_SCHEMA_VERSION,
        "library_snapshot_id": library_snapshot_id,
        "record_count": len(records),
        "record_alignment_sha256": record_alignment_sha256(records),
        "complete_library": bool(complete_library),
        "source_db": str(source_db) if source_db else None,
        "source_db_size": source_db.stat().st_size if source_db and source_db.exists() else None,
        "source_db_mtime_ns": source_db.stat().st_mtime_ns if source_db and source_db.exists() else None,
        "counts": counts,
        "files": files,
        "ordering": "ascending_mass_then_record_index",
        "score_or_rank_fields_stored": False,
    }
    manifest_path = output_dir / "manifest.json"
    temporary = manifest_path.with_name(f".{manifest_path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, manifest_path)
    return manifest


def load_library_mass_indices(
    output_dir: Path,
    *,
    expected_library_snapshot_id: str | None = None,
    expected_record_count: int | None = None,
    expected_record_alignment_sha256: str | None = None,
    require_complete_library: bool = True,
) -> tuple[MMapMassIndex, MMapMassIndex, dict]:
    output_dir = Path(output_dir)
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != MMAP_INDEX_SCHEMA_VERSION:
        raise ValueError("Unsupported mmap index schema")
    if (
        expected_library_snapshot_id is not None
        and manifest.get("library_snapshot_id") != expected_library_snapshot_id
    ):
        raise ValueError("Mmap index library snapshot mismatch")
    if (
        expected_record_count is not None
        and int(manifest.get("record_count") or -1) != expected_record_count
    ):
        raise ValueError("Mmap index record count mismatch")
    if (
        expected_record_alignment_sha256 is not None
        and manifest.get("record_alignment_sha256")
        != expected_record_alignment_sha256
    ):
        raise ValueError("Mmap index record alignment mismatch")
    if require_complete_library and not manifest.get("complete_library"):
        raise ValueError("Partial benchmark index cannot be used in production")
    precursor = MMapMassIndex(
        output_dir / manifest["files"]["precursor"]["masses"],
        output_dir / manifest["files"]["precursor"]["record_indices"],
    )
    neutral = MMapMassIndex(
        output_dir / manifest["files"]["neutral"]["masses"],
        output_dir / manifest["files"]["neutral"]["record_indices"],
    )
    if len(precursor) != int(manifest["counts"]["precursor"]):
        raise ValueError("Mmap precursor index count mismatch")
    if len(neutral) != int(manifest["counts"]["neutral"]):
        raise ValueError("Mmap neutral index count mismatch")
    return precursor, neutral, manifest
