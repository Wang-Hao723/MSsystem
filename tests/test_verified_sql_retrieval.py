"""Real SQLite boundary checks; no production library or sample tasks are used."""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import sqlite3

import pytest
from rdkit import Chem

from matching.evirank_retrieval import SQLiteCandidateRecallRetriever
from runtime import verified_sql_retrieval as M


class Backing(list):
    def read_uncached(self, index):
        return self[index]


def seal(payload):
    payload["semantic_sha256"] = M.digest({k: v for k, v in payload.items() if k != "semantic_sha256"})
    return payload


@pytest.fixture
def bound(tmp_path):
    source = tmp_path / "source.mgf"
    source.write_text("CHARGE=1-\nION_MODE=NEGATIVE\n")
    db = tmp_path / "library.db"
    peaks = [[20.0, 100.0], [30.0, 80.0], [40.0, 60.0]]
    key = Chem.MolToInchiKey(Chem.MolFromSmiles("CC(=O)O"))
    with sqlite3.connect(db) as connection:
        connection.executescript('''
            CREATE TABLE Compounds(id INTEGER PRIMARY KEY,compound_name TEXT,molecular_formula TEXT,
                smiles TEXT,inchikey TEXT,exact_mass REAL);
            CREATE TABLE Spectra(id INTEGER PRIMARY KEY,spectrum_id TEXT,source_file TEXT,
                precursor_mz REAL,peaks_json TEXT,adduct TEXT,charge INTEGER,ion_mode TEXT,
                source_accession TEXT,source_tier TEXT,dataset_source TEXT,quality_flags_json TEXT);
            CREATE TABLE Mapping(id INTEGER PRIMARY KEY,compound_id INTEGER,spectrum_id INTEGER);
            CREATE INDEX ix_s_pmz ON Spectra(precursor_mz);
        ''')
        for pk in (1, 2):
            connection.execute("INSERT INTO Compounds VALUES(?,?,?,?,?,?)",
                (pk, "acetate candidate", "C2H4O2", "CC(=O)O", key, 60.0211))
            connection.execute("INSERT INTO Spectra VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                (pk, f"spectrum-{pk}", str(source), 59.01, json.dumps(peaks), None, 1,
                 None if pk == 1 else "POSITIVE", None, "experimental", "fixture", "[]"))
            connection.execute("INSERT INTO Mapping VALUES(?,?,?)", (pk, pk, pk))
    backing = Backing([{"spectrum_id": "spectrum-1", "precursor_mz": 59.01,
        "charge": 1, "adduct": None, "peaks": peaks, "smiles": "CC(=O)O", "inchikey": key}])
    identity = {"source_files_sha256": "fixture-version", "python": [3, 13]}
    dbhash = hashlib.sha256(db.read_bytes()).hexdigest()
    backing.manifest = {"semantic_sha256": "fixture-cache", "record_count": 1,
        "normalized_records_sha256": "fixture-records", "source_db": {"sha256": dbhash},
        "processing_identity": identity}
    repair = {"record_index": 0, "spectrum_id": "spectrum-1", "original_record_sha256": M.digest(backing[0]),
        "updates": {"charge": -1, "source_ion_mode": "NEGATIVE"}, "source_polarity": "NEG",
        "source_binding": {"source_file": str(source), "block_sha256": "source-block",
            "structure_key_equal": True, "raw_peaks_equal": True,
            "existing_importer_peak_transform_equal": False, "verified_adduct_not_activated": "[M-H]-"}}
    registry = seal({"schema_version": "mssystem.source-metadata-repairs/1.0.0", "registry_id": "fixture",
        "base_normalized_records_sha256": "fixture-records", "base_record_count": 1,
        "source_database_sha256": dbhash, "rows": [repair],
        "source_files": {str(source): {"sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                                      "size_bytes": source.stat().st_size}}})
    binding = {"database_row": {"id": 1, "spectrum_id": "spectrum-1", "source_file": str(source),
        "precursor_mz": 59.01, "charge": 1, "adduct": None, "ion_mode": None,
        "source_accession": None, "inchikey": key, "structure_key": key.split("-")[0],
        "peak_fingerprint": M.peak_fingerprint(peaks)}, "record_index": 0,
        "updates": repair["updates"], "source_polarity": "NEG",
        "original_record_sha256": repair["original_record_sha256"], "source_block_sha256": "source-block"}
    (tmp_path / "config").mkdir()
    bridge = tmp_path / "bindings.jsonl"
    bridge.write_text(json.dumps(binding) + "\n")
    (tmp_path / "config/library_metadata_repairs.json").write_text(json.dumps(registry))
    (tmp_path / "config/library_record_store.json").write_text(json.dumps({"enabled": True, "cache_dir": str(tmp_path / "cache")}))
    config = seal({"schema_version": M.SCHEMA, "enabled": True, "database_writes_allowed": False,
        "source_database_path": str(db), "source_database_sha256": dbhash,
        "library_snapshot_id": "fixture-snapshot", "registry_path": "config/library_metadata_repairs.json",
        "registry_semantic_sha256": registry["semantic_sha256"], "bindings_path": str(bridge),
        "bindings_sha256": hashlib.sha256(bridge.read_bytes()).hexdigest(),
        "bindings_size_bytes": bridge.stat().st_size, "registered_reference_count": 1,
        "cache_manifest_semantic_sha256": "fixture-cache", "base_normalized_records_sha256": "fixture-records",
        "processing_identity": identity})
    (tmp_path / "config/sql_metadata_repairs.json").write_text(json.dumps(config))
    raw = SQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot")
    return tmp_path, db, backing, registry, [binding], config, raw


def state(bound, **kwargs):
    root, db, backing, registry, bindings, config, raw = bound
    return M.VerifiedSQLBindings(config, bindings, registry, backing, raw,
                                 project_root=root, dependencies=[db], **kwargs)


def test_restoration_before_polarity_filter_preserves_original_algorithm(bound):
    root, db, backing, registry, bindings, config, raw = bound
    before_db = db.read_bytes()
    original, _ = raw._fetch_window(59.01, "POS")
    fixed = M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=state(bound))
    positive, _ = fixed._fetch_window(59.01, "POS")
    negative, audit = fixed._fetch_window(59.01, "NEG")
    assert {x["spectrum_id"] for x in original} == {"spectrum-1", "spectrum-2"}
    assert {x["spectrum_id"] for x in positive} == {"spectrum-2"}
    assert {x["spectrum_id"] for x in negative} == {"spectrum-1"}
    row = negative[0]
    assert row["charge"] == -1 and row["ion_mode"] == "NEGATIVE"
    assert row["adduct"] is None  # A source-verified but deferred adduct stays inactive.
    assert row["peaks"] == backing[0]["peaks"] and row["precursor_mz"] == 59.01
    assert row["source_metadata_repair"]["record_index"] == 0
    assert audit["sql_metadata_repairs"]["complete_sql_and_original_cache_verified_count"] == 1
    assert fixed.library_snapshot_id == M.verified_sql_snapshot_id("fixture-snapshot", root)
    assert fixed.library_snapshot_id != raw.library_snapshot_id
    assert db.read_bytes() == before_db
    cached, cache_audit = fixed._fetch_window(59.01, "NEG")
    assert cached is negative and cache_audit["cache_hit"] is True
    assert fixed.metadata_bindings.rows_applied_count == 2


@pytest.mark.parametrize("field,value", [
    ("spectrum_id", "other"), ("source_file", "other"), ("precursor_mz", 60.0),
    ("charge", 2), ("adduct", "[M+H]+"), ("ion_mode", "POS"),
    ("source_accession", "changed"), ("inchikey", "WRONG"),
    ("peaks_json", "[[1,2]]"), ("smiles", "CCO"),
])
def test_complete_verification_rejects_changed_sql_original(bound, field, value):
    _, db, *_ = bound
    table = "Compounds" if field in {"inchikey", "smiles"} else "Spectra"
    with sqlite3.connect(db) as connection:
        connection.execute(f"UPDATE {table} SET {field}=? WHERE id=1", (value,))
    with pytest.raises(ValueError, match="SQL metadata original"):
        state(bound)


@pytest.mark.parametrize("damage", ["registry_hash", "extra_adduct", "original_cache", "membership", "cache_version", "snapshot"])
def test_reject_wrong_registry_generation_membership_and_cache(bound, damage):
    _, _, backing, registry, bindings, config, raw = bound
    if damage == "registry_hash": registry["registry_id"] = "different"
    if damage == "extra_adduct": bindings[0]["updates"] = {**bindings[0]["updates"], "adduct": "[M-H]-"}
    if damage == "original_cache": backing[0]["precursor_mz"] = 61.0
    if damage == "membership": bindings.append(deepcopy(bindings[0]))
    if damage == "cache_version": backing.manifest["processing_identity"] = {"source_files_sha256": "new"}
    if damage == "snapshot": config["library_snapshot_id"] = "different"
    with pytest.raises(ValueError): state(bound)


def test_sql_cached_windows_reject_database_change(bound):
    _, db, *_ = bound
    fixed = M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=state(bound))
    fixed._fetch_window(59.01, "NEG")
    with sqlite3.connect(db) as connection:
        connection.execute("UPDATE Spectra SET charge=2 WHERE id=1")
    with pytest.raises(ValueError, match="dependencies changed"):
        fixed._fetch_window(59.01, "NEG")


def test_new_verified_retrievers_reuse_only_the_checked_original_plan(bound, monkeypatch):
    _, db, *_, raw = bound
    verified = state(bound)
    expected = deepcopy(raw._sqlite_query_strategy)
    monkeypatch.setattr(SQLiteCandidateRecallRetriever, "_resolve_sqlite_query_strategy",
                        lambda *_: pytest.fail("verified generation must not repeat full Mapping scans"))
    first = M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=verified)
    second = M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=verified)
    assert first._sqlite_query_strategy == second._sqlite_query_strategy == expected
    assert first._sqlite_query_strategy is not second._sqlite_query_strategy
    records, audit = first._fetch_window(59.01, "NEG")
    assert [record["spectrum_id"] for record in records] == ["spectrum-1"]
    assert audit["sql_metadata_repairs"]["complete_sql_and_original_cache_verified_count"] == 1
    assert first._window_cache is not second._window_cache


def test_new_retriever_rejects_changed_database_before_plan_reuse(bound):
    _, db, *_ = bound
    verified = state(bound)
    with sqlite3.connect(db) as connection:
        connection.execute("UPDATE Spectra SET charge=2 WHERE id=1")
    with pytest.raises(ValueError, match="dependencies changed"):
        M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=verified)


def test_local_connection_non_candidate_sql_is_unchanged_and_readonly(bound):
    _, db, *_, raw = bound
    fixed = M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=state(bound))
    a, b = raw._connect(), fixed._connect()
    try:
        for sql in ("PRAGMA table_info(Spectra)", "SELECT count(*) FROM Spectra"):
            assert [tuple(x) for x in a.execute(sql)] == [tuple(x) for x in b.execute(sql)]
        with pytest.raises(sqlite3.OperationalError):
            b.execute("UPDATE Spectra SET charge=2")
        with pytest.raises(ValueError, match="Unrecognized SQL"):
            b.execute("SELECT s.peaks_json FROM Spectra s")
    finally:
        a.close(); b.close()


def test_per_candidate_digest_rejects_changed_row_after_full_check(bound):
    _, _, _, _, bindings, _, raw = bound
    verified = state(bound)
    connection = raw._connect()
    try:
        row = dict(connection.execute('''SELECT s.id AS database_spectrum_id,s.spectrum_id,s.source_file,
            s.precursor_mz,s.peaks_json,s.charge,s.adduct,s.ion_mode,s.source_accession,c.smiles,c.inchikey
            FROM Spectra s JOIN Mapping m ON m.id=s.id JOIN Compounds c ON c.id=m.compound_id WHERE s.id=1''').fetchone())
    finally:
        connection.close()
    row["peaks_json"] = "[[2,3]]"
    with pytest.raises(ValueError, match="original candidate changed"):
        verified.apply(row)


@pytest.mark.parametrize("damage", ["config", "bridge", "source", "processing"])
def test_preflight_rejects_hash_and_processing_mismatch_before_loading_cache(bound, damage):
    root, *_ = bound
    path = root / "config/sql_metadata_repairs.json"
    if damage == "config":
        payload = json.loads(path.read_text()); payload["registered_reference_count"] = 2; path.write_text(json.dumps(payload))
    if damage == "bridge": (root / "bindings.jsonl").write_text("changed")
    if damage == "source": (root / "source.mgf").write_text("changed")
    # The fixture version intentionally differs from actual production code.
    result = M.sql_metadata_repair_preflight(root)
    assert result["compatible"] is False
    if damage == "processing": assert "processing version" in result["error"]


def test_sql_preflight_uses_shared_compatible_loader_without_rebinding_cache(bound, monkeypatch):
    root, db, backing, _, _, config, _ = bound
    current = {**config["processing_identity"], "source_files_sha256": "current-with-reviewed-additions"}
    calls = {}

    import matching.library_record_store as record_store
    import runtime.library_cache_identity as identity_resolver
    import runtime.manifest as runtime_manifest

    monkeypatch.setattr(record_store, "processing_identity", lambda: current)

    def resolve(accepted, observed, *, project_root):
        calls["resolve"] = (accepted, observed, project_root)
        return {"compatible": True, "mode": "verified_isolated_non_normalizing_additions"}

    def load(cache, *, source_db, current_identity, project_root, cache_bytes):
        calls["load"] = (cache, source_db, current_identity, project_root, cache_bytes)
        backing.close = lambda: None
        return backing, object(), object()

    monkeypatch.setattr(identity_resolver, "processing_identity_compatibility", resolve)
    monkeypatch.setattr(identity_resolver, "load_compatible_record_store", load)
    monkeypatch.setattr(runtime_manifest, "sqlite_snapshot", lambda path: {"snapshot_id": "fixture-snapshot"})
    result = M.sql_metadata_repair_preflight(root)
    assert result["compatible"] is True
    assert calls["resolve"] == (config["processing_identity"], current, root.resolve())
    assert calls["load"] == (
        (root / "cache").resolve(), db.resolve(), current, root.resolve(), 0,
    )
    assert backing.manifest["processing_identity"] == config["processing_identity"]


def test_sql_preflight_rejects_unproved_identity_before_cache_load(bound, monkeypatch):
    root, *_ = bound
    imported_loader = {"called": False}
    import runtime.library_cache_identity as identity_resolver

    monkeypatch.setattr(
        identity_resolver,
        "processing_identity_compatibility",
        lambda accepted, current, *, project_root: {
            "compatible": False,
            "reason": "unreviewed_processing_source_change",
        },
    )

    def should_not_load(*args, **kwargs):
        imported_loader["called"] = True
        raise AssertionError("cache loader must not run")

    monkeypatch.setattr(identity_resolver, "load_compatible_record_store", should_not_load)
    result = M.sql_metadata_repair_preflight(root)
    assert result == {
        "compatible": False,
        "status": "invalid_binding",
        "error": "SQL metadata cache processing version mismatch: unreviewed_processing_source_change",
    }
    assert imported_loader["called"] is False


def test_nonproduction_and_unconfigured_fixtures_keep_base_retriever(bound, tmp_path):
    root, db, *_ = bound
    other = tmp_path / "other-root"; other.mkdir()
    unconfigured = M.create_verified_sql_retriever(db, library_snapshot_id="fixture-snapshot", project_root=other)
    assert type(unconfigured) is SQLiteCandidateRecallRetriever
    another_db = tmp_path / "another.db"; another_db.write_bytes(db.read_bytes())
    different_db = M.create_verified_sql_retriever(another_db, library_snapshot_id="other-db", project_root=root)
    assert type(different_db) is SQLiteCandidateRecallRetriever
    assert M.sql_metadata_repair_preflight(other) == {"compatible": True, "status": "not_configured"}


def test_state_is_shared_and_status_counts_actual_instances_only(bound, monkeypatch):
    root, db, *_ = bound
    verified = state(bound)
    monkeypatch.setitem(M._STATES, str(root.resolve()), verified)
    assert M.sql_metadata_repair_status(project_root=root)["verified_retrievers_created"] == 0
    a = M.create_verified_sql_retriever(db, library_snapshot_id="fixture-snapshot", project_root=root)
    b = M.create_verified_sql_retriever(db, library_snapshot_id="fixture-snapshot", project_root=root)
    assert a.metadata_bindings is b.metadata_bindings is verified
    assert M.sql_metadata_repair_status(project_root=root)["verified_retrievers_created"] == 2
    assert M.sql_metadata_repair_status(project_root=root)["candidate_row_applications"] == 0
    assert M.sql_metadata_repair_preflight(root)["complete_sql_and_original_cache_verified_count"] == 1
    assert db in M.sql_metadata_dependency_paths(root)
    assert "metadata:" in M.verified_sql_snapshot_id("fixture-snapshot", root)
    with pytest.raises(ValueError, match="different library snapshot"):
        M.verified_sql_snapshot_id("old", root)


def test_deleted_configuration_cannot_silently_restore_uncorrected_sql(bound, monkeypatch):
    root, db, backing, registry, bindings, config, raw = bound
    verified = M.VerifiedSQLBindings(config, bindings, registry, backing, raw,
        project_root=root, dependencies=M.sql_metadata_dependency_paths(root))
    monkeypatch.setitem(M._STATES, str(root.resolve()), verified)
    (root / "config/sql_metadata_repairs.json").unlink()
    assert M.sql_metadata_repair_preflight(root)["compatible"] is False
    with pytest.raises(ValueError, match="dependencies changed"):
        M.create_verified_sql_retriever(db, library_snapshot_id="fixture-snapshot", project_root=root)
    with pytest.raises(ValueError, match="dependencies changed"):
        M.verified_sql_snapshot_id("fixture-snapshot", root)


def test_cold_production_snapshot_rejects_missing_sql_configuration(bound, monkeypatch):
    root, *_ = bound
    monkeypatch.setattr(M, "ROOT", root)
    (root / "config/sql_metadata_repairs.json").unlink()
    with pytest.raises(ValueError, match="configuration is missing"):
        M.verified_sql_snapshot_id("fixture-snapshot", root)


def test_source_and_bridge_changes_invalidate_shared_state(bound, monkeypatch):
    root, db, backing, registry, bindings, config, raw = bound
    paths = M.sql_metadata_dependency_paths(root) + [root / "source.mgf"]
    verified = M.VerifiedSQLBindings(config, bindings, registry, backing, raw,
        project_root=root, dependencies=paths)
    monkeypatch.setitem(M._STATES, str(root.resolve()), verified)
    before = M.sql_metadata_dependency_signature(root)
    (root / "bindings.jsonl").write_text("modified")
    assert M.sql_metadata_dependency_signature(root) != before
    assert M.sql_metadata_repair_preflight(root)["compatible"] is False
    assert M.sql_metadata_repair_status(project_root=root)["status"] == "invalid_binding"


def test_dense_window_budget_still_fails_closed(bound):
    _, db, *_ = bound
    fixed = M.VerifiedSQLiteCandidateRecallRetriever(db, library_snapshot_id="fixture-snapshot", bindings=state(bound))
    fixed.policy["library_gate"]["maximum_sql_rows_per_precursor_window"] = 1
    records, audit = fixed._fetch_window(59.01, "NEG")
    assert records == []
    assert audit["status"] == "dense_window_budget_exceeded_fail_closed"
    assert audit["sql_metadata_repairs"]["candidate_row_applications"] == 0


def _dependency_probe(module, paths):
    state = object.__new__(module.VerifiedSQLBindings)
    state.mass_bindings = None  # This fixture exercises only dependency freshness.
    state.dependencies = tuple(paths)
    state.signature = module._signature(paths)
    state.cache_override = os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE', '')
    return state


@pytest.mark.parametrize('suffix', ['-wal', '-journal'])
def test_empty_sidecar_permission_only_change_preserves_generation(tmp_path, suffix):
    p = tmp_path/('library.db'+suffix)
    p.write_bytes(b'')
    new = _dependency_probe(M, [p])
    before = p.stat()
    p.chmod(0o600 if before.st_mode & 0o777 != 0o600 else 0o640)
    assert p.stat().st_ctime_ns != before.st_ctime_ns
    assert p.stat().st_mtime_ns == before.st_mtime_ns and p.stat().st_size == 0
    new.assert_fresh()


@pytest.mark.parametrize('kind', ['empty_wal_grows', 'nonempty_same_length_restore_mtime',
                                  'main_empty_ctime', 'empty_wal_inode', 'empty_wal_mtime', 'deleted'])
def test_meaningful_dependency_changes_still_rejected(tmp_path, kind):
    p = tmp_path/('library.db' if kind == 'main_empty_ctime' else 'library.db-wal')
    p.write_bytes(b'old' if kind == 'nonempty_same_length_restore_mtime' else b'')
    before = p.stat(); state = _dependency_probe(M, [p])
    if kind == 'empty_wal_grows':
        p.write_bytes(b'database pages')
    elif kind == 'nonempty_same_length_restore_mtime':
        p.write_bytes(b'new'); os.utime(p, ns=(before.st_atime_ns, before.st_mtime_ns))
    elif kind == 'main_empty_ctime':
        p.chmod(0o600)
    elif kind == 'empty_wal_inode':
        replacement = tmp_path/'replacement'; replacement.write_bytes(b'')
        os.utime(replacement, ns=(before.st_atime_ns, before.st_mtime_ns)); replacement.replace(p)
    elif kind == 'empty_wal_mtime':
        os.utime(p, ns=(before.st_atime_ns, before.st_mtime_ns+1_000_000))
    else:
        p.unlink()
    with pytest.raises(ValueError, match='dependencies changed'):
        state.assert_fresh()


def test_real_wal_commit_invalidates_without_checkpoint(tmp_path):
    p = tmp_path/'library.db'
    with sqlite3.connect(p) as connection:
        assert connection.execute('PRAGMA journal_mode=WAL').fetchone()[0] == 'wal'
        connection.execute('CREATE TABLE test(value)'); connection.commit()
        connection.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        wal = Path(str(p)+'-wal'); assert wal.stat().st_size == 0
        state = _dependency_probe(M, [p, wal]); before = p.read_bytes()
        connection.execute('INSERT INTO test VALUES (1)'); connection.commit()
        assert p.read_bytes() == before and wal.stat().st_size > 0
        with pytest.raises(ValueError, match='dependencies changed'):
            state.assert_fresh()


def test_real_main_database_commit_invalidates(tmp_path):
    p = tmp_path/'library.db'
    with sqlite3.connect(p) as connection:
        connection.execute('CREATE TABLE test(value)'); connection.commit()
        state = _dependency_probe(M, [p])
        connection.execute('INSERT INTO test VALUES (1)'); connection.commit()
        with pytest.raises(ValueError, match='dependencies changed'):
            state.assert_fresh()
