"""Source-bound metadata at the read-only SQL candidate retrieval boundary.

The original retriever still owns its SELECT, window budget, ordering, polarity
filter and scoring. A local connection adapter adds two provenance columns and
restores only fields already active in the normalized-reference registry.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sqlite3
import threading
from typing import Any

from matching.evirank_retrieval import SQLiteCandidateRecallRetriever
from runtime.library_metadata_repairs import digest, validate_registry
from runtime.source_file_verification import verify_registered_sources

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "mssystem.sql-metadata-repairs/1.0.0"
_LOCK = threading.RLock()
_STATES: dict[str, "VerifiedSQLBindings"] = {}
_CREATED = 0
_APPLIED = 0
_SELECT_PREFIX = "SELECT s.precursor_mz,s.peaks_json,s.spectrum_id,s.adduct,"
_ROW_FIELDS = (
    "database_spectrum_id", "spectrum_id", "source_file", "precursor_mz",
    "charge", "adduct", "ion_mode", "source_accession", "inchikey",
    "smiles", "peaks_json",
)


def _path(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def _stat(path: Path) -> tuple:
    try:
        stat = path.stat()
        # An empty SQLite sidecar has no database pages. Permission-only
        # ctime changes cannot change its contents; still bind its inode,
        # size and mtime, and keep the full signature for all nonempty files.
        ctime = (None if stat.st_size == 0 and str(path).endswith(("-wal", "-journal"))
                 else stat.st_ctime_ns)
        return (str(path), stat.st_dev, stat.st_ino, stat.st_size,
                stat.st_mtime_ns, ctime)
    except FileNotFoundError:
        return (str(path), "missing")


def _signature(paths) -> tuple:
    return tuple(_stat(Path(path)) for path in paths)


def _hash(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _read_config(root: Path) -> dict | None:
    path = root / "config/sql_metadata_repairs.json"
    if not path.exists():
        return None
    config = json.loads(path.read_text(encoding="utf-8"))
    if (config.get("schema_version") != SCHEMA
            or config.get("semantic_sha256") != digest({k: v for k, v in config.items()
                                                       if k != "semantic_sha256"})
            or config.get("enabled") is not True
            or config.get("database_writes_allowed") is not False):
        raise ValueError("Invalid SQL metadata repair configuration")
    return config


def sql_metadata_dependency_paths(project_root=ROOT) -> list[Path]:
    """Paths whose stat changes invalidate a cached readiness/repair result."""
    root = Path(project_root).resolve()
    with _LOCK:
        previous = _STATES.get(str(root))
        if previous is not None:
            # Retain old dependencies after a config change: readiness rejects
            # the changed generation until restart and full verification.
            return list(previous.dependencies)
    from runtime.source_compound_mass_repairs import source_compound_mass_dependency_paths
    paths = source_compound_mass_dependency_paths(root) + [root / "config/sql_metadata_repairs.json",
             root / "config/library_record_store.json",
             root / "config/library_metadata_repairs.json"]
    try:
        config = _read_config(root)
        if config:
            paths.extend([_path(root, config["bindings_path"]),
                          _path(root, config["registry_path"]),
                          _path(root, config["source_database_path"])])
            store = json.loads((root / "config/library_record_store.json").read_text())
            cache = _path(root, os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "").strip()
                          or store["cache_dir"])
            from matching.library_record_store import FILE_NAMES
            paths.extend(cache / name for name in (*FILE_NAMES, "manifest.json", "acceptance.json"))
    except (OSError, ValueError, TypeError, KeyError):
        pass  # The preflight reports the precise invalid configuration.
    return list(dict.fromkeys(paths))


def sql_metadata_dependency_signature(project_root=ROOT) -> tuple:
    return (os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", ""),
            _signature(sql_metadata_dependency_paths(project_root)))


def peak_fingerprint(peaks) -> str:
    return hashlib.sha256(json.dumps(sorted([[float(a), float(b)] for a, b in peaks]),
                                    separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _row_digest(row: dict) -> str:
    return digest({field: row.get(field) for field in _ROW_FIELDS})


def _verify_original_sql_row(raw: dict, binding: dict) -> None:
    expected = binding["database_row"]
    for field in _ROW_FIELDS[:9]:
        expected_field = "id" if field == "database_spectrum_id" else field
        if raw.get(field) != expected.get(expected_field):
            raise ValueError("SQL metadata original row mismatch: " + field)
    if peak_fingerprint(json.loads(raw["peaks_json"])) != expected["peak_fingerprint"]:
        raise ValueError("SQL metadata original peak fingerprint mismatch")
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(raw.get("smiles") or "")
    key = Chem.MolToInchiKey(molecule).split("-")[0] if molecule else None
    if not key or key != expected["structure_key"]:
        raise ValueError("SQL metadata original SMILES structure mismatch")


class VerifiedSQLBindings:
    """One completely checked generation, shared by retriever instances."""

    def __init__(self, config: dict, bindings: list[dict], registry: dict,
                 backing, raw_retriever: SQLiteCandidateRecallRetriever, *,
                 dependencies=(), project_root=ROOT, progress=None):
        from matching.library_polarity import classify_spectrum_polarity
        self.config = dict(config)
        self.project_root = Path(project_root).resolve()
        self.dependencies = tuple(dependencies)
        self.signature = _signature(self.dependencies)
        self.cache_override = os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "")
        self.db_path = raw_retriever.db_path
        self.base_snapshot_id = raw_retriever.library_snapshot_id
        self.library_snapshot_id = (self.base_snapshot_id + ":metadata:"
            + config["registry_semantic_sha256"] + ":sql:" + config["bindings_sha256"])
        self.mass_bindings = None
        self.rows: dict[int, dict] = {}
        self.instance_count = 0
        self.fields_applied_count = 0
        self.rows_applied_count = 0
        registry = validate_registry(registry)
        if registry["semantic_sha256"] != config["registry_semantic_sha256"]:
            raise ValueError("SQL and normalized metadata registry generations differ")
        if raw_retriever.library_snapshot_id != config["library_snapshot_id"]:
            raise ValueError("SQL repair database snapshot mismatch")
        manifest = backing.manifest
        if (manifest["normalized_records_sha256"] != config["base_normalized_records_sha256"]
                or manifest["normalized_records_sha256"] != registry["base_normalized_records_sha256"]
                or manifest["record_count"] != registry["base_record_count"]
                or manifest["source_db"]["sha256"] != config["source_database_sha256"]
                or manifest["source_db"]["sha256"] != registry["source_database_sha256"]
                or manifest["semantic_sha256"] != config["cache_manifest_semantic_sha256"]
                or manifest["processing_identity"] != config["processing_identity"]):
            raise ValueError("SQL repair normalized cache binding mismatch")
        by_index = {row["record_index"]: row for row in registry["rows"]}
        if len(bindings) != len(by_index) or len(bindings) != config["registered_reference_count"]:
            raise ValueError("SQL repair registry membership mismatch")
        by_pk = {}
        seen_indices = set()
        for binding in bindings:
            index = binding["record_index"]
            expected = by_index.get(index)
            primary_key = binding["database_row"]["id"]
            if (not isinstance(primary_key, int) or isinstance(primary_key, bool)
                    or primary_key < 1 or primary_key in by_pk or index in seen_indices
                    or expected is None
                    or binding["database_row"]["spectrum_id"] != expected["spectrum_id"]
                    or binding["database_row"]["source_file"] != expected["source_binding"]["source_file"]
                    or binding["updates"] != expected["updates"]
                    or binding["source_polarity"] != expected["source_polarity"]
                    or binding["original_record_sha256"] != expected["original_record_sha256"]
                    or binding["source_block_sha256"] != expected["source_binding"]["block_sha256"]):
                raise ValueError("SQL repair row differs from active source-bound registry")
            by_pk[primary_key] = binding
            seen_indices.add(index)
        if seen_indices != set(by_index):
            raise ValueError("SQL repair registry coverage mismatch")
        # Keep the original verified Mapping join. Remove only the precursor
        # index hint for these release checks targeting exact primary keys.
        strategy = raw_retriever._sqlite_query_strategy
        # This exact plan already passed the original provider's complete
        # Mapping coverage/alignment checks against this database generation.
        # Keep its original precursor index hint for candidate retrieval; only
        # the release-check SELECT below removes that hint for exact PK reads.
        self.verified_query_strategy = dict(strategy)
        from_clause = strategy["from_clause"]
        if strategy.get("precursor_index"):
            from_clause = from_clause.replace(" INDEXED BY " +
                raw_retriever._quoted_identifier(strategy["precursor_index"]), "")
        self.query_strategy = {k: v for k, v in strategy.items() if k != "from_clause"}
        connection = raw_retriever._connect()
        keys = sorted(by_pk)
        try:
            for offset in range(0, len(keys), 400):
                batch = keys[offset:offset + 400]
                placeholders = ",".join("?" for _ in batch)
                raw_rows = connection.execute(
                    "SELECT s.id AS database_spectrum_id,s.spectrum_id,s.source_file,"
                    "s.precursor_mz,s.peaks_json,s.charge,s.adduct,s.ion_mode,s.source_accession,"
                    "m.spectrum_id AS mapping_spectrum_id,c.smiles,c.inchikey FROM "
                    + from_clause + " WHERE s.id IN (" + placeholders + ")", batch).fetchall()
                if len(raw_rows) != len(batch):
                    raise ValueError("SQL repair primary-key join is not one-to-one")
                seen_batch = set()
                for row in raw_rows:
                    raw = dict(row)
                    primary_key = raw["database_spectrum_id"]
                    if primary_key not in batch or primary_key in seen_batch or raw["mapping_spectrum_id"] != primary_key:
                        raise ValueError("SQL repair mapping identity mismatch")
                    seen_batch.add(primary_key)
                    binding = by_pk[primary_key]
                    _verify_original_sql_row(raw, binding)
                    original = backing.read_uncached(binding["record_index"])
                    if (digest(original) != binding["original_record_sha256"]
                            or original["spectrum_id"] != raw["spectrum_id"]
                            or peak_fingerprint(original["peaks"]) != binding["database_row"]["peak_fingerprint"]):
                        raise ValueError("SQL repair original normalized reference mismatch")
                    fixed = dict(raw)
                    for field, value in binding["updates"].items():
                        fixed["ion_mode" if field == "source_ion_mode" else field] = value
                    polarity = classify_spectrum_polarity(fixed.get("adduct"), fixed.get("charge"), fixed.get("ion_mode"))
                    if polarity != binding["source_polarity"]:
                        raise ValueError("SQL repaired polarity differs from normalized reference")
                    self.rows[primary_key] = {
                        "original_sql_row_sha256": _row_digest(raw),
                        "record_index": binding["record_index"],
                        "source_block_sha256": binding["source_block_sha256"],
                        "updates": dict(binding["updates"]),
                    }
                if progress:
                    progress(min(offset + len(batch), len(keys)), len(keys))
        finally:
            connection.close()
        if len(self.rows) != config["registered_reference_count"]:
            raise ValueError("Incomplete SQL original-record verification")
        from runtime.source_compound_mass_repairs import get_source_compound_mass_state, source_compound_mass_snapshot_id
        self.mass_bindings = get_source_compound_mass_state(
            self.project_root, backing=backing, raw_retriever=raw_retriever)
        self.library_snapshot_id = source_compound_mass_snapshot_id(self.library_snapshot_id, self.project_root)
        self.assert_fresh()

    def assert_fresh(self):
        if self.mass_bindings is not None:
            self.mass_bindings.assert_fresh()
        if (self.cache_override != os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "")
                or self.signature != _signature(self.dependencies)):
            raise ValueError("Verified SQL metadata dependencies changed; restart and reverify")

    def apply(self, row) -> dict:
        """Check each standalone or streamed row before exposing its result."""
        self.assert_fresh()
        result = self._apply_checked_row(row)
        self.assert_fresh()
        return result

    def _apply_checked_row(self, row) -> dict:
        # Internal to an entry/exit generation-checked cursor batch or apply().
        # The per-repair original-row digest remains mandatory.
        global _APPLIED
        raw = dict(row)
        binding = self.rows.get(raw["database_spectrum_id"])
        if binding is None:
            return self.mass_bindings._apply_sql_checked_row(raw) if self.mass_bindings is not None else raw
        if _row_digest(raw) != binding["original_sql_row_sha256"]:
            raise ValueError("SQL original candidate changed after complete verification")
        for field, value in binding["updates"].items():
            raw["ion_mode" if field == "source_ion_mode" else field] = value
        raw["source_metadata_repair"] = {
            "registry_semantic_sha256": self.config["registry_semantic_sha256"],
            "bindings_sha256": self.config["bindings_sha256"],
            "record_index": binding["record_index"],
            "source_block_sha256": binding["source_block_sha256"],
            "fields_restored": sorted(binding["updates"]),
            "source_database_written": False,
        }
        with _LOCK:
            self.fields_applied_count += len(binding["updates"])
            self.rows_applied_count += 1
            _APPLIED += 1
        return self.mass_bindings._apply_sql_checked_row(raw) if self.mass_bindings is not None else raw

    def status(self):
        return {
            "status": "applied_to_retriever" if self.instance_count else "bindings_verified_no_retriever_created",
            "registered_reference_count": len(self.rows),
            "complete_sql_and_original_cache_verified_count": len(self.rows),
            "verified_retrievers_created": self.instance_count,
            "candidate_row_applications": self.rows_applied_count,
            "metadata_field_applications": self.fields_applied_count,
            "registry_semantic_sha256": self.config["registry_semantic_sha256"],
            "bindings_sha256": self.config["bindings_sha256"],
            "base_library_snapshot_id": self.base_snapshot_id,
            "library_snapshot_id": self.library_snapshot_id,
            "source_database_written": False,
            "peaks_structures_precursor_masses_changed": False,
            "source_compound_mass_repairs": (self.mass_bindings.status() if self.mass_bindings is not None else {"status": "not_configured"}),
        }


class _BoundCursor:
    def __init__(self, cursor, bindings):
        self._cursor, self._bindings = cursor, bindings

    def fetchall(self):
        self._bindings.assert_fresh()
        result = [self._bindings._apply_checked_row(row) for row in self._cursor.fetchall()]
        self._bindings.assert_fresh()
        return result

    def fetchone(self):
        self._bindings.assert_fresh()
        row = self._cursor.fetchone()
        result = None if row is None else self._bindings._apply_checked_row(row)
        self._bindings.assert_fresh()
        return result

    def fetchmany(self, size=None):
        self._bindings.assert_fresh()
        rows = self._cursor.fetchmany() if size is None else self._cursor.fetchmany(size)
        result = [self._bindings._apply_checked_row(row) for row in rows]
        self._bindings.assert_fresh()
        return result

    def __iter__(self):
        for row in self._cursor:
            yield self._bindings.apply(row)

    def __getattr__(self, name):
        return getattr(self._cursor, name)


class _BoundConnection:
    def __init__(self, connection, bindings):
        self._connection, self._bindings = connection, bindings

    def execute(self, sql, parameters=()):
        if sql.startswith(_SELECT_PREFIX):
            self._bindings.assert_fresh()
            # A changed base SELECT must be reviewed explicitly, not guessed.
            if ("c.smiles,c.inchikey,c.exact_mass FROM " not in sql
                    or "ORDER BY s.precursor_mz,s.spectrum_id,c.id" not in sql):
                raise ValueError("SQL candidate projection changed; adapter review required")
            sql = "SELECT s.id AS database_spectrum_id,s.source_file," + sql[len("SELECT "):]
            return _BoundCursor(self._connection.execute(sql, parameters), self._bindings)
        if sql.lstrip().upper().startswith("SELECT") and "s.peaks_json" in sql:
            raise ValueError("Unrecognized SQL candidate projection; adapter review required")
        return self._connection.execute(sql, parameters)

    def __getattr__(self, name):
        return getattr(self._connection, name)


class VerifiedSQLiteCandidateRecallRetriever(SQLiteCandidateRecallRetriever):
    def __init__(self, db_path, *, library_snapshot_id, bindings, policy=None):
        global _CREATED
        self.metadata_bindings = bindings
        bindings.assert_fresh()
        if (Path(db_path).resolve() != bindings.db_path
                or library_snapshot_id != bindings.base_snapshot_id):
            raise ValueError("Verified SQL retriever requested a different database snapshot")
        super().__init__(db_path, library_snapshot_id=bindings.library_snapshot_id, policy=policy)
        with _LOCK:
            bindings.instance_count += 1
            _CREATED += 1

    def _resolve_sqlite_query_strategy(self):
        # A new sample needs a fresh window cache, not three more complete
        # Mapping scans. Reuse only the fully verified binding's original plan,
        # with the same database/config/source freshness guard on both sides.
        self.metadata_bindings.assert_fresh()
        strategy = dict(self.metadata_bindings.verified_query_strategy)
        self.metadata_bindings.assert_fresh()
        return strategy

    def _connect(self):
        self.metadata_bindings.assert_fresh()
        return _BoundConnection(super()._connect(), self.metadata_bindings)

    def _fetch_window(self, precursor_mz, query_polarity):
        # Check even cache hits; a stale generation cannot keep returning rows.
        self.metadata_bindings.assert_fresh()
        records, audit = super()._fetch_window(precursor_mz, query_polarity)
        result = records, {**audit, "sql_metadata_repairs": self.metadata_bindings.status()}
        # Also cover normalization and cached-window return, after cursor use.
        self.metadata_bindings.assert_fresh()
        return result


def _load_state(root: Path, *, progress=None) -> VerifiedSQLBindings | None:
    with _LOCK:
        previous = _STATES.get(str(root))
        if previous is not None:
            previous.assert_fresh()
            return previous
        config = _read_config(root)
        if config is None:
            if root == ROOT and (root / "config/library_metadata_repairs.json").is_file():
                raise ValueError("Production SQL metadata configuration is missing")
            return None
        registry_path = _path(root, config["registry_path"])
        registry = validate_registry(json.loads(registry_path.read_text(encoding="utf-8")))
        if registry["semantic_sha256"] != config["registry_semantic_sha256"]:
            raise ValueError("Configured SQL and normalized registry generations differ")
        bridge_path = _path(root, config["bindings_path"])
        db_path = _path(root, config["source_database_path"])
        dependencies = sql_metadata_dependency_paths(root)
        dependencies.extend(Path(name) for name in registry["source_files"])
        dependencies.extend(Path(str(db_path) + suffix) for suffix in ("-wal", "-journal"))
        dependencies = list(dict.fromkeys(dependencies))
        before = _signature(dependencies)
        if (bridge_path.stat().st_size != config["bindings_size_bytes"]
                or _hash(bridge_path) != config["bindings_sha256"]):
            raise ValueError("SQL metadata bridge file hash mismatch")
        verify_registered_sources(registry["source_files"])
        with bridge_path.open(encoding="utf-8") as handle:
            bindings = [json.loads(line) for line in handle]
        from matching.library_record_store import processing_identity
        from runtime.library_cache_identity import (
            load_compatible_record_store,
            processing_identity_compatibility,
        )
        store_config = json.loads((root / "config/library_record_store.json").read_text())
        override = os.environ.get("MSSYSTEM_LIBRARY_RECORD_STORE", "").strip()
        if override.lower() in {"off", "0", "false"} or (not override and store_config.get("enabled") is not True):
            raise ValueError("SQL metadata repairs require their accepted normalized cache")
        cache_path = _path(root, override or store_config["cache_dir"])
        current_identity = processing_identity()
        identity_resolution = processing_identity_compatibility(
            config["processing_identity"], current_identity, project_root=root,
        )
        if identity_resolution.get("compatible") is not True:
            raise ValueError(
                "SQL metadata cache processing version mismatch: "
                + str(identity_resolution.get("reason") or "unverified source change")
            )
        # Resolve only the reviewed, byte-exact isolated additions.  The shared
        # loader then invokes the original strict loader with the cache's
        # accepted identity, retaining every payload, acceptance, source-DB,
        # and manifest check.  VerifiedSQLBindings below separately requires
        # that accepted identity to equal this SQL generation's binding.
        backing, precursor_index, neutral_index = load_compatible_record_store(
            cache_path, source_db=db_path, current_identity=current_identity,
            project_root=root, cache_bytes=0)
        try:
            from runtime.manifest import sqlite_snapshot
            if sqlite_snapshot(db_path)["snapshot_id"] != config["library_snapshot_id"]:
                raise ValueError("SQL metadata operational database snapshot mismatch")
            raw_retriever = SQLiteCandidateRecallRetriever(db_path, library_snapshot_id=config["library_snapshot_id"])
            state = VerifiedSQLBindings(config, bindings, registry, backing, raw_retriever,
                                        dependencies=dependencies, project_root=root, progress=progress)
        finally:
            backing.close()
            del precursor_index, neutral_index
        if before != _signature(dependencies):
            raise ValueError("SQL metadata dependencies changed during verification")
        _STATES[str(root)] = state
        return state


def create_verified_sql_retriever(db_path, *, library_snapshot_id, policy=None, project_root=ROOT):
    """Use repairs only for the precisely configured production database."""
    root = Path(project_root).resolve()
    previous = _STATES.get(str(root))
    if previous is not None:
        previous.assert_fresh()
    config = _read_config(root)
    if config is None and root == ROOT:
        from config.settings import DB_PATH
        if Path(db_path).resolve() == Path(DB_PATH).resolve() and (root / "config/library_metadata_repairs.json").is_file():
            raise ValueError("Production SQL metadata configuration is missing")
    if config is None or Path(db_path).resolve() != _path(root, config["source_database_path"]):
        return SQLiteCandidateRecallRetriever(db_path, library_snapshot_id=library_snapshot_id, policy=policy)
    state = _load_state(root)
    return VerifiedSQLiteCandidateRecallRetriever(db_path, library_snapshot_id=library_snapshot_id,
                                                  policy=policy, bindings=state)


def verified_sql_snapshot_id(base_id, project_root=ROOT) -> str:
    """Return the configured generation identity without creating a retriever.

    This is an artifact compatibility key, not a claim of database verification;
    preflight/status separately report that actual verification state.
    """
    root = Path(project_root).resolve()
    state = _STATES.get(str(root))
    if state is not None:
        state.assert_fresh()
    config = _read_config(root)
    if config is None:
        if root == ROOT and (root / "config/library_metadata_repairs.json").is_file():
            raise ValueError("Production SQL metadata configuration is missing")
        return str(base_id)
    if str(base_id) != config["library_snapshot_id"]:
        raise ValueError("SQL metadata requested a different library snapshot")
    from runtime.source_compound_mass_repairs import source_compound_mass_snapshot_id
    return source_compound_mass_snapshot_id(
        str(base_id) + ":metadata:" + config["registry_semantic_sha256"]
        + ":sql:" + config["bindings_sha256"], root)


def sql_metadata_repair_preflight(project_root=ROOT, *, progress=None) -> dict[str, Any]:
    try:
        state = _load_state(Path(project_root).resolve(), progress=progress)
        if state is None:
            return {"compatible": True, "status": "not_configured"}
        return {"compatible": True, **state.status()}
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, sqlite3.Error) as exc:
        return {"compatible": False, "status": "invalid_binding", "error": str(exc)}


def sql_metadata_repair_status(retriever=None, *, project_root=ROOT) -> dict[str, Any]:
    if retriever is not None and isinstance(retriever, VerifiedSQLiteCandidateRecallRetriever):
        state = retriever.metadata_bindings
    else:
        state = _STATES.get(str(Path(project_root).resolve()))
    if state is None:
        return {"status": "not_initialized", "verified_retrievers_created": 0,
                "candidate_row_applications": 0}
    try:
        state.assert_fresh()
    except ValueError as exc:
        return {"status": "invalid_binding", "error": str(exc),
                "verified_retrievers_created": state.instance_count}
    return state.status()
