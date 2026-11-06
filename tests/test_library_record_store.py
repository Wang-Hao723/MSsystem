import copy
import ast
import hashlib
import importlib
import json
import math
import sqlite3
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

from matching import library_record_store as store
from matching.offline_mz_index import record_alignment_sha256
from matching.search import _neutral_candidates, _query_candidates


def records():
    return [{"spectrum_id": str(i), "precursor_mz": mz, "neutral_mass": mz - 1.,
        "dataset_source": "GNPS", "compound_name": "测试 + original name",
        "peaks": [[20., 3], [35., 12.5], [45., 4]],
        "clean_peaks": [[20., 3.], [35., 12.5], [45., 4.]],
        "adduct": "[M+H]+", "quality_flags": ["keep"], "name_only_lead": False,
        "collision_energy_max": None}
        for i, mz in enumerate([301., 299.99, 301., 100., 300.98, 301.02])]


def fixture(tmp_path, *, complete=True):
    source = tmp_path / "test.db"
    source.write_bytes(b"read-only source placeholder")
    root = tmp_path / "cache"
    original = records()
    manifest = store.build_record_store(iter(original), root, source_db=source,
        identity={"test": 1}, complete_library=complete)
    return source, root, manifest, original


def open_fixture(source, root, **kwargs):
    return store.load_record_store(root, source_db=source, identity={"test": 1},
                                   require_acceptance=False, **kwargs)


def buckets(values, field):
    output = {}
    for i, row in enumerate(values):
        output.setdefault(int(row[field] / .1), []).append(i)
    for indices in output.values():
        indices.sort(key=lambda i: values[i][field])
    return output


def test_roundtrip_preserves_every_field_order_and_alignment(tmp_path):
    source, root, manifest, original = fixture(tmp_path)
    loaded, mz, neutral = open_fixture(source, root)
    assert list(loaded) == original
    assert loaded[-1] == original[-1]
    assert loaded[1:5:2] == original[1:5:2]
    assert manifest["record_alignment_sha256"] == record_alignment_sha256(original)
    assert loaded.source_families
    digest = hashlib.sha256()
    for i in range(len(loaded)):
        digest.update(store.canonical_bytes(loaded.read_uncached(i)) + b"\n")
    assert digest.hexdigest() == manifest["normalized_records_sha256"]
    loaded.close()


def test_cached_and_bucket_candidate_order_equal_including_ties(tmp_path):
    source, root, _, original = fixture(tmp_path)
    loaded, mz, neutral = open_fixture(source, root)
    for value in [100., 301., 300.98, 299.99, 499.]:
        assert _query_candidates(value, original, buckets(original, "precursor_mz")) == _query_candidates(value, loaded, mz)
        assert _neutral_candidates(value, original, buckets(original, "neutral_mass")) == _neutral_candidates(value, loaded, neutral)
    loaded.close()


def test_cache_uses_a_byte_budget_and_does_not_mutate_backing_records(tmp_path):
    source, root, _, original = fixture(tmp_path)
    size = store._deep_size(original[0]) + 300
    loaded, _, _ = open_fixture(source, root, cache_bytes=size * 2)
    for repeat in range(3):
        for i in range(len(loaded)):
            assert loaded[i] == original[i]
            assert loaded.cache_info()["estimated_decoded_bytes"] <= size * 2
    assert loaded.cache_info()["decoded_records"] < len(loaded)
    loaded[0]["compound_name"] = "in-memory mutation"
    assert loaded.read_uncached(0)["compound_name"] == original[0]["compound_name"]
    assert store.file_sha256(root / "records.jsonl") == loaded.manifest["files"]["records.jsonl"]["sha256"]
    loaded.close()


def test_zero_budget_and_parallel_reads(tmp_path):
    source, root, _, original = fixture(tmp_path)
    loaded, _, _ = open_fixture(source, root, cache_bytes=0)
    with ThreadPoolExecutor(max_workers=4) as pool:
        assert list(pool.map(lambda i: loaded[i % 6], range(120))) == [original[i % 6] for i in range(120)]
    assert loaded.cache_info()["decoded_records"] == 0
    assert loaded.cache_info()["estimated_decoded_bytes"] == 0
    with pytest.raises(IndexError):
        loaded[6]
    with pytest.raises(IndexError):
        loaded[-7]
    with pytest.raises(TypeError):
        loaded[1.2]
    loaded.close()


@pytest.mark.parametrize("filename", store.FILE_NAMES)
def test_corrupt_payload_or_index_fails_before_loading(tmp_path, filename):
    source, root, _, _ = fixture(tmp_path)
    with (root / filename).open("ab") as handle:
        handle.write(b"bad")
    with pytest.raises(ValueError, match="hash mismatch"):
        open_fixture(source, root)


def test_missing_readiness_partial_wrong_source_or_code_are_rejected(tmp_path):
    source, root, _, _ = fixture(tmp_path, complete=False)
    with pytest.raises(ValueError, match="completeness"):
        open_fixture(source, root)


def test_no_automatic_activation_and_mismatched_identity(tmp_path):
    source, root, _, _ = fixture(tmp_path)
    with pytest.raises(FileNotFoundError):
        store.load_record_store(root, source_db=source, identity={"test": 1})
    with pytest.raises(ValueError, match="identity"):
        store.load_record_store(root, source_db=source, identity={"test": 2}, require_acceptance=False)
    source.write_bytes(b"changed source")
    with pytest.raises(ValueError, match="source database content"):
        open_fixture(source, root)


@pytest.mark.parametrize("suffix", ["-wal", "-journal"])
def test_active_sqlite_sidecars_are_not_ignored(tmp_path, suffix):
    source, root, _, _ = fixture(tmp_path)
    Path(str(source) + suffix).write_bytes(b"uncheckpointed content")
    with pytest.raises(RuntimeError, match="WAL/journal"):
        open_fixture(source, root)


def test_rebuilding_never_overwrites_previous_cache(tmp_path):
    source, root, manifest, original = fixture(tmp_path)
    before = store.file_sha256(root / "manifest.json")
    with pytest.raises(FileExistsError):
        store.build_record_store(iter(original), root, source_db=source, identity={"test": 1})
    assert store.file_sha256(root / "manifest.json") == before


def test_failed_build_never_publishes_a_complete_manifest(tmp_path):
    source = tmp_path / "source.db"
    source.write_bytes(b"source")
    root = tmp_path / "broken"
    def broken():
        yield records()[0]
        raise RuntimeError("interrupted")
    with pytest.raises(RuntimeError, match="interrupted"):
        store.build_record_store(broken(), root, source_db=source, identity={"test": 1})
    assert not root.exists()


def valid_acceptance(manifest):
    payload = {"schema_version": "mssystem.library-record-store-acceptance/1.0.0",
        "manifest_semantic_sha256": manifest["semantic_sha256"],
        "full_record_roundtrip_sha256": manifest["normalized_records_sha256"],
        "verified_record_count": manifest["record_count"], "candidate_order_cases": 10000,
        "candidate_order_mismatches": 0, "paired_scoring_queries": 138, "score_mismatches": 0,
        "shared_database_unchanged": True}
    payload["semantic_sha256"] = store.semantic_sha256(payload)
    return payload


@pytest.mark.parametrize("field,value", [("score_mismatches", 1), ("candidate_order_mismatches", 1),
    ("paired_scoring_queries", 99), ("candidate_order_cases", 9999), ("verified_record_count", 1),
    ("shared_database_unchanged", False), ("manifest_semantic_sha256", "wrong")])
def test_no_cutover_with_incomplete_or_failed_acceptance(tmp_path, field, value):
    _, _, manifest, _ = fixture(tmp_path)
    payload = valid_acceptance(manifest)
    store.validate_acceptance(payload, manifest)
    payload[field] = value
    payload["semantic_sha256"] = store.semantic_sha256(
        {key: item for key, item in payload.items() if key != "semantic_sha256"}
    )
    with pytest.raises(ValueError, match="acceptance"):
        store.validate_acceptance(payload, manifest)


def test_runtime_config_and_explicit_off_switch(tmp_path, monkeypatch):
    source, root, manifest, original = fixture(tmp_path)
    (root / "acceptance.json").write_text(json.dumps(valid_acceptance(manifest)))
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"enabled": True, "source_db": str(source), "cache_dir": str(root)}))
    monkeypatch.setattr(store, "CONFIG_PATH", config)
    monkeypatch.setattr(store, "processing_identity", lambda: {"test": 1})
    monkeypatch.delenv("MSSYSTEM_LIBRARY_RECORD_STORE", raising=False)
    loaded, _, _ = store.try_load_runtime_record_store(source)
    assert loaded[0] == original[0]
    loaded.close()
    assert store.try_load_runtime_record_store(tmp_path / "other.db") is None
    monkeypatch.setenv("MSSYSTEM_LIBRARY_RECORD_STORE", "off")
    assert store.try_load_runtime_record_store(source) is None


def make_database(path):
    with sqlite3.connect(path) as conn:
        conn.executescript('''
            CREATE TABLE Spectra (id INTEGER PRIMARY KEY, spectrum_id TEXT, precursor_mz REAL,
                peaks_json TEXT, adduct TEXT, retention_time REAL, charge INTEGER, source_tier TEXT);
            CREATE TABLE Compounds (id INTEGER PRIMARY KEY, compound_name TEXT, smiles TEXT,
                inchikey TEXT, molecular_formula TEXT, exact_mass REAL, classification TEXT);
            CREATE TABLE Mapping (spectrum_id INTEGER, compound_id INTEGER);
        ''')
        for i, (mass, exact, adduct, tier, peaks) in enumerate([
            (101., 100., "[M+H]+", "experimental", [[20., 3], [35., 12.5], [45., 4]]),
            (301., None, "[M+H]+", "experimental", [[20., 3], [35., 12.5], [45., 4]]),
            (301., None, "Positive", "experimental", [[20., 3], [35., 12.5], [45., 4]]),
            (301., 300., "[M+H]+", "predicted", [[20., 3], [35., 12.5], [45., 4]]),
            (301., 300., "[M+H]+", "experimental", [[20., 3]]),
        ], 1):
            conn.execute("INSERT INTO Spectra VALUES (?,?,?,?,?,?,?,?)", (i, str(i), mass, json.dumps(peaks), adduct, 1., 1, tier))
            conn.execute("INSERT INTO Compounds VALUES (?,?,?,?,?,?,?)", (i, "Example", "CCO", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "C2H6O", exact, "Alcohol"))
            conn.execute("INSERT INTO Mapping VALUES (?,?)", (i, i))


def test_streaming_and_legacy_loaders_keep_all_records_and_gates(tmp_path, monkeypatch):
    search = importlib.import_module("matching.search")
    path = tmp_path / "library.db"
    make_database(path)
    monkeypatch.setenv("MSSYSTEM_LIBRARY_RECORD_STORE", "off")
    for name, value in (("_lib_records", None), ("_lib_db_path", None), ("_mz_index", None), ("_nm_index", None), ("_lib_source_families", ())):
        monkeypatch.setattr(search, name, value)
    before = store.file_sha256(path)
    streamed = list(search.iter_normalized_library_records(path))
    loaded, _ = search.load_library(path)
    assert len(streamed) == 3
    assert streamed == loaded
    assert streamed[0]["neutral_mass"] == 100.
    assert streamed[1]["neutral_mass_source"] == "library_adduct"
    assert "neutral_mass" not in streamed[2]
    assert store.file_sha256(path) == before


def test_native_loader_uses_accepted_cache_and_keeps_small_database_switch(tmp_path, monkeypatch):
    search = importlib.import_module("matching.search")
    full, small = tmp_path / "full.db", tmp_path / "small.db"
    make_database(full)
    make_database(small)
    original = list(search.iter_normalized_library_records(full))
    root = tmp_path / "full-cache"
    manifest = store.build_record_store(original, root, source_db=full, identity={"test": 1})
    (root / "acceptance.json").write_text(json.dumps(valid_acceptance(manifest)))
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"enabled": True, "source_db": str(full), "cache_dir": str(root)}))
    monkeypatch.setattr(store, "CONFIG_PATH", config)
    monkeypatch.setattr(store, "processing_identity", lambda: {"test": 1})
    monkeypatch.delenv("MSSYSTEM_LIBRARY_RECORD_STORE", raising=False)
    for name, value in (("_lib_records", None), ("_lib_db_path", None), ("_mz_index", None),
                        ("_nm_index", None), ("_lib_source_families", ())):
        monkeypatch.setattr(search, name, value)

    def reject_preload(*args, **kwargs):
        raise AssertionError("An enabled full cache must not call the legacy preloader")

    with monkeypatch.context() as cutover:
        cutover.setattr(search, "iter_normalized_library_records", reject_preload)
        loaded, index = search.load_library(full)
        assert isinstance(loaded, store.MMapRecordStore)
        assert isinstance(index, store.MMapMassIndex)
        assert isinstance(search._nm_index, store.MMapMassIndex)
        assert len(loaded) == len(original)
        assert store.canonical_bytes(loaded[0]) == store.canonical_bytes(original[0])
        assert search._lib_source_families == loaded.source_families
        again, again_index = search.load_library(full)
        assert again is loaded and again_index is index
    switched, small_index = search.load_library(small)
    assert isinstance(switched, list) and isinstance(small_index, dict)
    assert search._lib_db_path == small.resolve()
    assert switched == list(search.iter_normalized_library_records(small))
    loaded.close()


def test_native_loader_does_not_silently_preload_when_enabled_cache_is_missing(tmp_path, monkeypatch):
    search = importlib.import_module("matching.search")
    full = tmp_path / "full.db"
    make_database(full)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"enabled": True, "source_db": str(full),
                                  "cache_dir": str(tmp_path / "missing-cache")}))
    monkeypatch.setattr(store, "CONFIG_PATH", config)
    monkeypatch.setattr(store, "processing_identity", lambda: {"test": 1})
    monkeypatch.delenv("MSSYSTEM_LIBRARY_RECORD_STORE", raising=False)
    monkeypatch.setattr(search, "_lib_records", None)
    monkeypatch.setattr(search, "_lib_db_path", None)
    with pytest.raises(FileNotFoundError):
        search.load_library(full)
    assert search._lib_records is None


def test_legacy_nonfinite_raw_peaks_are_preserved_not_silently_rewritten(tmp_path):
    source = tmp_path / "source.db"
    source.write_bytes(b"source")
    original = records()
    original[0]["peaks"].append([float("nan"), float("inf")])
    root = tmp_path / "nan-peaks"
    store.build_record_store(original, root, source_db=source, identity={"test": 1})
    loaded, _, _ = open_fixture(source, root)
    assert math.isnan(loaded[0]["peaks"][-1][0])
    assert loaded[0]["peaks"][-1][1] == float("inf")
    assert loaded[0]["clean_peaks"] == original[0]["clean_peaks"]
    assert store.canonical_bytes(loaded[0]) == store.canonical_bytes(original[0])
    loaded.close()


def test_processing_identity_resolves_the_module_not_matching_search_function():
    identity = store.processing_identity()
    assert len(identity["source_files_sha256"]) == 64
    assert identity["library_top_n_peaks"] > 0


def test_per_record_filtering_ast_matches_frozen_pre_refactor_loader():
    # SHA-256 of the 22 per-row statements before append in commit 67e0fe7.
    # Only neutral-mass preparation and yielding replace the final append.
    search = importlib.import_module("matching.search")
    tree = ast.parse(Path(search.__file__).read_text())
    func = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "iter_normalized_library_records")
    loop = next(x for x in ast.walk(func) if isinstance(x, ast.For) and isinstance(x.iter, ast.Name) and x.iter.id == "rows")
    source = ast.dump(ast.Module(body=loop.body[:-4], type_ignores=[]), include_attributes=False).encode()
    assert hashlib.sha256(source).hexdigest() == "10994dd3c75608961fb6259e91bf5a2e4edf9a2047239320b89c6909a9b27982"
