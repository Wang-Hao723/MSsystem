"""Read-only, content-bound, on-demand normalized spectral-library records.

JSON records preserve Python numeric values and every normalized library field.
Only offsets and two small mass indices are NumPy arrays; the spectrum records
are decoded on demand into a byte-budgeted LRU. No pickle, database writes, or
candidate/scoring changes are involved. Building never activates a cache.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import mmap
import operator
import os
import sys
import tempfile
import threading
from array import array
from collections import OrderedDict
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

from matching.offline_mz_index import MMapMassIndex, _finite, _record_precursor_mass

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "mssystem.library-record-store/1.0.0"
CONFIG_PATH = ROOT / "config/library_record_store.json"
FILE_NAMES = ("records.jsonl", "record_offsets.npy", "precursor_masses.npy",
    "precursor_record_indices.npy", "neutral_masses.npy", "neutral_record_indices.npy")
DEFAULT_CACHE_BYTES = 128 * 1024 * 1024


def canonical_bytes(value) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                      allow_nan=True).encode("utf-8")


def semantic_sha256(value) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def source_binding(path: Path) -> dict:
    path = Path(path).resolve()
    before = path.stat()
    wal = Path(str(path) + "-wal")
    journal = Path(str(path) + "-journal")
    if any(p.exists() and p.stat().st_size for p in (wal, journal)):
        raise RuntimeError("Record cache cannot bind a database with an active WAL/journal")
    digest = file_sha256(path)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if (any(getattr(before, key) != getattr(after, key) for key in fields)
            or any(p.exists() and p.stat().st_size for p in (wal, journal))):
        raise RuntimeError("Source database changed while computing its content binding")
    return {"path": str(path), "size_bytes": after.st_size, "sha256": digest}


def processing_identity() -> dict:
    """Bind processing code, curated names/roles, options and numeric runtime.

    Web/UI/evaluation files are excluded because they do not normalize library
    records. Matching, utilities, classification and their JSON assets are not.
    """
    from runtime.manifest import source_tree_snapshot
    search = importlib.import_module("matching.search")
    files = [item for item in source_tree_snapshot(ROOT)["files"]
        if item["path"].startswith(("matching/", "utils/", "mssystem/"))
        or item["path"] in {"config/settings.py", "config/retrieval_policy.json"}]
    return {"source_files_sha256": semantic_sha256(files),
        "python": list(sys.version_info[:2]), "numpy": np.__version__,
        "load_name_only_leads": bool(search.LOAD_NAME_ONLY_LEADS),
        "library_top_n_peaks": search.TOP_N_PEAKS, "library_relative_threshold": search.REL_THRESHOLD}


def _deep_size(value, seen=None) -> int:
    seen = set() if seen is None else seen
    identity = id(value)
    if identity in seen:
        return 0
    seen.add(identity)
    size = sys.getsizeof(value)
    if isinstance(value, dict):
        size += sum(_deep_size(k, seen) + _deep_size(v, seen) for k, v in value.items())
    elif isinstance(value, (list, tuple)):
        size += sum(_deep_size(v, seen) for v in value)
    return size


def _write_json_new(path: Path, payload: dict) -> None:
    with path.open("xb") as handle:
        handle.write(canonical_bytes(payload) + b"\n")
        handle.flush()
        os.fsync(handle.fileno())


def build_record_store(records: Iterable[dict], output_dir: Path, *, source_db: Path,
                       identity: dict, complete_library: bool = True, progress=None) -> dict:
    """Stream a normalized library into a new immutable cache directory.

    An interrupted build leaves its uniquely named staging directory for
    inspection. Existing cache directories and shared databases are untouched.
    """
    from matching.conditioned_multilibrary import normalize_source_family
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Record cache already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    source = source_binding(source_db)
    stage = Path(tempfile.mkdtemp(prefix=".record-store-build-", dir=output_dir.parent))
    offsets = array("Q", [0])
    masses = {name: array("d") for name in ("precursor", "neutral")}
    indices = {name: array("Q") for name in masses}
    record_digest, alignment = hashlib.sha256(), hashlib.sha256()
    families = set()
    count = 0
    with (stage / "records.jsonl").open("xb") as handle:
        for record in records:
            if not isinstance(record, dict):
                raise TypeError("Normalized library records must be dictionaries")
            encoded = canonical_bytes(record) + b"\n"
            handle.write(encoded)
            record_digest.update(encoded)
            offsets.append(offsets[-1] + len(encoded))
            precursor, neutral = _record_precursor_mass(record), _finite(record.get("neutral_mass"))
            alignment.update(json.dumps((count, precursor, neutral, str(record.get("spectrum_id") or "")),
                separators=(",", ":"), ensure_ascii=False).encode("utf-8") + b"\n")
            for field, value in (("precursor", precursor), ("neutral", neutral)):
                if value is not None:
                    masses[field].append(value)
                    indices[field].append(count)
            family = normalize_source_family(record.get("dataset_source"))
            if family:
                families.add(family)
            count += 1
            if progress is not None and count % 100000 == 0:
                progress(count)
        handle.flush()
        os.fsync(handle.fileno())
    if not count:
        raise ValueError(f"An empty record cache cannot be published; staging: {stage}")
    np.save(stage / "record_offsets.npy", np.asarray(offsets, dtype=np.uint64), allow_pickle=False)
    counts = {}
    for field in masses:
        mz = np.asarray(masses[field], dtype=np.float64)
        idx = np.asarray(indices[field], dtype=np.uint32 if count < 2**32 else np.uint64)
        order = np.lexsort((idx, mz))
        np.save(stage / f"{field}_masses.npy", mz[order], allow_pickle=False)
        np.save(stage / f"{field}_record_indices.npy", idx[order], allow_pickle=False)
        counts[field] = len(mz)
    if source_binding(source_db) != source:
        raise RuntimeError(f"Database changed while building cache; staging: {stage}")
    files = {name: {"sha256": file_sha256(stage / name), "size_bytes": (stage / name).stat().st_size}
             for name in FILE_NAMES}
    manifest = {"schema_version": SCHEMA, "complete_library": complete_library is True,
        "source_db": source, "processing_identity": identity, "record_count": count,
        "normalized_records_sha256": record_digest.hexdigest(),
        "record_alignment_sha256": alignment.hexdigest(), "mass_counts": counts,
        "source_families": sorted(families), "files": files,
        "record_format": "canonical UTF-8 JSON with LF; offsets index original runtime order",
        "score_fields_changed": False, "production_activation_allowed_without_acceptance": False}
    manifest["semantic_sha256"] = semantic_sha256(manifest)
    _write_json_new(stage / "manifest.json", manifest)
    # rename is atomic on the same volume. The destination must not already
    # exist; no previous cache generation is removed or replaced.
    if output_dir.exists():
        raise FileExistsError(f"Concurrent cache publication: {output_dir}; staging: {stage}")
    stage.rename(output_dir)
    return manifest


class MMapRecordStore(Sequence):
    """A list-compatible immutable backing store with a bounded decoded cache."""

    def __init__(self, root: Path, manifest: dict, *, cache_bytes: int = DEFAULT_CACHE_BYTES):
        if isinstance(cache_bytes, bool) or cache_bytes < 0:
            raise ValueError("Decoded cache byte budget must be nonnegative")
        self.root, self.manifest = Path(root), manifest
        self.source_families = tuple(manifest["source_families"])
        self._offsets = np.load(self.root / "record_offsets.npy", mmap_mode="r", allow_pickle=False)
        expected = int(manifest["record_count"])
        if (self._offsets.dtype != np.dtype("uint64") or self._offsets.ndim != 1
                or len(self._offsets) != expected + 1 or self._offsets[0] != 0
                or np.any(self._offsets[1:] <= self._offsets[:-1])
                or self._offsets[-1] != (self.root / "records.jsonl").stat().st_size):
            raise ValueError("Record offset array is invalid or misaligned")
        self._handle = (self.root / "records.jsonl").open("rb")
        self._mapped = mmap.mmap(self._handle.fileno(), 0, access=mmap.ACCESS_READ)
        self._cache = OrderedDict()
        self._cache_bytes = self._hits = self._misses = 0
        self._budget = int(cache_bytes)
        self._lock = threading.RLock()

    def __len__(self):
        return int(self.manifest["record_count"])

    def read_uncached(self, index: int) -> dict:
        index = operator.index(index)
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        left, right = int(self._offsets[index]), int(self._offsets[index + 1])
        record = json.loads(self._mapped[left:right])
        if not isinstance(record, dict):
            raise ValueError("Cached record is not a dictionary")
        return record

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        index = operator.index(index)
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        with self._lock:
            if index in self._cache:
                self._hits += 1
                self._cache.move_to_end(index)
                return self._cache[index][0]
        record = self.read_uncached(index)
        size = _deep_size(record) + sys.getsizeof(index) + 128
        with self._lock:
            self._misses += 1
            if index in self._cache:
                self._cache.move_to_end(index)
                return self._cache[index][0]
            if size <= self._budget:
                while self._cache and self._cache_bytes + size > self._budget:
                    _, (_, removed_size) = self._cache.popitem(last=False)
                    self._cache_bytes -= removed_size
                self._cache[index] = (record, size)
                self._cache_bytes += size
        return record

    def cache_info(self) -> dict:
        with self._lock:
            return {"decoded_records": len(self._cache), "estimated_decoded_bytes": self._cache_bytes,
                "decoded_byte_budget": self._budget, "hits": self._hits, "misses": self._misses,
                "record_file_bytes": self.manifest["files"]["records.jsonl"]["size_bytes"]}

    def close(self):
        if hasattr(self, "_lock"):
            with self._lock:
                self._cache.clear()
                self._cache_bytes = 0
        if hasattr(self, "_mapped"):
            self._mapped.close()
        if hasattr(self, "_handle"):
            self._handle.close()

    def __del__(self):
        self.close()


def validate_acceptance(payload: dict, manifest: dict) -> None:
    if (payload.get("schema_version") != "mssystem.library-record-store-acceptance/1.0.0"
            or payload.get("semantic_sha256") != semantic_sha256({k: v for k, v in payload.items() if k != "semantic_sha256"})
            or payload.get("manifest_semantic_sha256") != manifest["semantic_sha256"]
            or payload.get("full_record_roundtrip_sha256") != manifest["normalized_records_sha256"]
            or payload.get("verified_record_count") != manifest["record_count"]
            or payload.get("candidate_order_cases", 0) < 10000
            or payload.get("candidate_order_mismatches") != 0
            or payload.get("paired_scoring_queries", 0) < 100
            or payload.get("score_mismatches") != 0
            or payload.get("shared_database_unchanged") is not True):
        raise ValueError("Record cache lacks complete record/retrieval/scoring acceptance")


def load_record_store(root: Path, *, source_db: Path, identity: dict,
                      require_acceptance: bool = True, cache_bytes: int = DEFAULT_CACHE_BYTES):
    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    digest = semantic_sha256({k: v for k, v in manifest.items() if k != "semantic_sha256"})
    if (manifest.get("schema_version") != SCHEMA or manifest.get("semantic_sha256") != digest
            or manifest.get("complete_library") is not True
            or manifest.get("processing_identity") != identity
            or set(manifest.get("files") or {}) != set(FILE_NAMES)):
        raise ValueError("Record cache manifest, completeness or processing identity mismatch")
    if require_acceptance:
        validate_acceptance(json.loads((root / "acceptance.json").read_text(encoding="utf-8")), manifest)
    if manifest["source_db"] != source_binding(source_db):
        raise ValueError("Record cache source database content mismatch")
    for name in FILE_NAMES:
        path = root / name
        expected = manifest["files"][name]
        if (path.is_symlink() or path.stat().st_size != expected["size_bytes"]
                or file_sha256(path) != expected["sha256"]):
            raise ValueError(f"Record cache file hash mismatch: {name}")
    indices = []
    for field in ("precursor", "neutral"):
        index = MMapMassIndex(root / f"{field}_masses.npy", root / f"{field}_record_indices.npy")
        if (index._masses.dtype != np.dtype("float64") or index._record_indices.dtype.kind != "u"
                or len(index) != manifest["mass_counts"][field]
                or np.any(~np.isfinite(index._masses)) or np.any(index._masses <= 0)
                or (len(index) and int(np.max(index._record_indices)) >= manifest["record_count"])):
            raise ValueError("Record cache mass arrays are invalid")
        indices.append(index)
    return MMapRecordStore(root, manifest, cache_bytes=cache_bytes), *indices


def try_load_runtime_record_store(db_path: Path):
    """Load only an explicitly enabled, independently accepted cache generation."""
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8")) if CONFIG_PATH.is_file() else {}
    override = os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "").strip()
    if override.lower() in {"off", "0", "false"}:
        return None
    if not override and config.get("enabled") is not True:
        return None
    bound_db = str(config.get("source_db") or "").strip()
    if not override and (not bound_db or Path(bound_db).resolve() != Path(db_path).resolve()):
        return None
    root = Path(override or config["cache_dir"]).expanduser().resolve()
    return load_record_store(root, source_db=db_path, identity=processing_identity())
