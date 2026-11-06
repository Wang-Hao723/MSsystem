from __future__ import annotations

import json
import sqlite3

from web import server


def _counts(seed: int) -> dict:
    return {
        "spectra": seed,
        "physical_spectra": seed,
        "source_searchable_spectra": seed - 1,
        "runtime_searchable_spectra": seed - 2,
        "excluded_by_source": 1,
        "excluded_by_runtime_structure_gate": 1,
        "runtime_projection_exact": True,
        "compounds": seed // 2,
        "source_tiers": {"experimental": seed - 1},
    }


def test_current_count_snapshot_avoids_opening_sqlite(tmp_path, monkeypatch):
    small = tmp_path / "small.db"
    full = tmp_path / "full.db"
    small.write_bytes(b"small-library-identity")
    full.write_bytes(b"full-library-identity")
    snapshot = tmp_path / "library_counts_snapshot.json"
    snapshot.write_text(
        json.dumps(
            {
                "schema_version": "mssystem.library-counts-snapshot/1.0.0",
                "generated_at": "2026-08-31T00:00:00Z",
                "load_name_only_leads": bool(server.LOAD_NAME_ONLY_LEADS),
                "libraries": {
                    "small": {
                        "path": str(small.resolve()),
                        "size_bytes": small.stat().st_size,
                        "mtime_ns": small.stat().st_mtime_ns,
                        "counts": _counts(100),
                    },
                    "full": {
                        "path": str(full.resolve()),
                        "size_bytes": full.stat().st_size,
                        "mtime_ns": full.stat().st_mtime_ns,
                        "counts": _counts(200),
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    def forbidden_connect(*_args, **_kwargs):
        raise AssertionError("current library counts snapshot opened SQLite")

    monkeypatch.setattr(server, "DB_PATH_SMALL", small)
    monkeypatch.setattr(server, "DB_PATH", full)
    monkeypatch.setattr(server, "_LIBRARY_COUNTS_SNAPSHOT_PATH", snapshot)
    monkeypatch.setattr(sqlite3, "connect", forbidden_connect)
    monkeypatch.setattr(server, "_lib_loaded", False)
    server._library_count_cache.clear()

    result = server._library_counts_sync_unlocked()

    assert result["small"]["runtime_searchable_spectra"] == 98
    assert result["full"]["runtime_searchable_spectra"] == 198
    assert result["small"]["count_snapshot_status"] == "current"
    assert result["full"]["count_snapshot_status"] == "current"


def test_changed_library_stat_invalidates_snapshot(tmp_path, monkeypatch):
    library = tmp_path / "small.db"
    library.write_bytes(b"before")
    snapshot = tmp_path / "library_counts_snapshot.json"
    snapshot.write_text(
        json.dumps(
            {
                "schema_version": "mssystem.library-counts-snapshot/1.0.0",
                "generated_at": "2026-08-31T00:00:00Z",
                "load_name_only_leads": bool(server.LOAD_NAME_ONLY_LEADS),
                "libraries": {
                    "small": {
                        "path": str(library.resolve()),
                        "size_bytes": library.stat().st_size,
                        "mtime_ns": library.stat().st_mtime_ns,
                        "counts": _counts(100),
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(server, "_LIBRARY_COUNTS_SNAPSHOT_PATH", snapshot)
    assert server._library_counts_snapshot_entry("small", library) is not None

    library.write_bytes(b"after-stat-change")

    assert server._library_counts_snapshot_entry("small", library) is None
