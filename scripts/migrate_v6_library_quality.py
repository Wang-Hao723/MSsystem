#!/usr/bin/env python3
"""Apply the reversible V6 library provenance/contaminant migration.

This migration deliberately does not delete spectra.  It:

* quarantines the old ``upload_template_literature`` pseudo-spectra, whose four
  reported ions have no experimental intensities;
* quarantines one malformed AO1680 precursor record while preserving the other
  experimental spectra;
* corrects oxidized/unoxidized Irgafos 168 names and class metadata;
* adds provenance/quality columns and indexes used by the V6 loader;
* labels pre-provenance Small-library rows ``legacy_unresolved`` so they remain
  usable only under the corresponding confidence cap.

Before any semantic write, an APFS copy-on-write clone is created and checked.
The script is idempotent through ``LibraryQualityMigrations``.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import time
from pathlib import Path


MIGRATION_ID = "v6_20260715_library_quality_roles_targets"
OXIDIZED_IK = "AZSKHRTUXHLAHS-UHFFFAOYSA-N"
IRGAFOS_IK = "JKIJEFPNVSHHEI-UHFFFAOYSA-N"
BAD_AO_SPECTRUM_ID = "MassBank_NIST_2399986"
LITERATURE_DATASET = "upload_template_literature"


def _connect_readonly(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)


def _columns(conn: sqlite3.Connection, table: str) -> set[str]:
    return {str(row[1]) for row in conn.execute(f"PRAGMA table_info({table})")}


def _quick_check(path: Path) -> None:
    conn = _connect_readonly(path)
    try:
        result = conn.execute("PRAGMA quick_check").fetchone()
    finally:
        conn.close()
    if not result or result[0] != "ok":
        raise RuntimeError(f"SQLite quick_check failed for {path}: {result}")


def _table_row_counts(path: Path) -> dict[str, int]:
    conn = _connect_readonly(path)
    try:
        tables = [
            str(row[0]) for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' "
                "AND name NOT LIKE 'sqlite_%' ORDER BY name"
            )
        ]
        counts = {}
        for table in tables:
            quoted = table.replace('"', '""')
            counts[table] = int(
                conn.execute(f'SELECT COUNT(*) FROM "{quoted}"').fetchone()[0]
            )
        return counts
    finally:
        conn.close()


def _foreign_key_errors(path: Path) -> list[tuple]:
    conn = _connect_readonly(path)
    try:
        return list(conn.execute("PRAGMA foreign_key_check"))
    finally:
        conn.close()


def _already_applied(path: Path) -> bool:
    conn = _connect_readonly(path)
    try:
        tables = {
            row[0] for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        if "LibraryQualityMigrations" not in tables:
            return False
        return conn.execute(
            "SELECT 1 FROM LibraryQualityMigrations WHERE migration_id=?",
            (MIGRATION_ID,),
        ).fetchone() is not None
    finally:
        conn.close()


def _clone_backup(path: Path) -> Path:
    """Create one online-consistent SQLite backup and audit its contents."""
    stamp = time.strftime("%Y%m%dT%H%M%S")
    backup = path.with_name(f"{path.name}.pre_v6_quality_{stamp}.bak")
    counter = 1
    while backup.exists():
        backup = path.with_name(
            f"{path.name}.pre_v6_quality_{stamp}_{counter}.bak"
        )
        counter += 1
    source_counts = _table_row_counts(path)
    source = _connect_readonly(path)
    destination = sqlite3.connect(str(backup))
    try:
        source.execute("PRAGMA busy_timeout=60000")
        destination.execute("PRAGMA busy_timeout=60000")
        source.backup(destination)
        destination.commit()
    except Exception:
        destination.close()
        source.close()
        backup.unlink(missing_ok=True)
        raise
    finally:
        try:
            destination.close()
        except Exception:
            pass
        try:
            source.close()
        except Exception:
            pass

    _quick_check(backup)
    fk_errors = _foreign_key_errors(backup)
    if fk_errors:
        backup.unlink(missing_ok=True)
        raise RuntimeError(f"SQLite backup foreign_key_check failed: {fk_errors[:5]}")
    backup_counts = _table_row_counts(backup)
    if backup_counts != source_counts:
        backup.unlink(missing_ok=True)
        raise RuntimeError(
            f"SQLite backup table/row counts differ: {source_counts} != {backup_counts}"
        )
    return backup


def _expected_literature_count(path: Path) -> int:
    return 123 if "small" in path.stem.lower() else 180


def audit(path: Path) -> dict[str, object]:
    conn = _connect_readonly(path)
    try:
        spectra_columns = _columns(conn, "Spectra")
        source_expr = "source_tier" if "source_tier" in spectra_columns else "NULL"
        return {
            "database": str(path),
            "spectra": int(conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0]),
            "literature_rows": int(conn.execute(
                "SELECT COUNT(*) FROM Spectra WHERE dataset_source=?",
                (LITERATURE_DATASET,),
            ).fetchone()[0]),
            "bad_ao_rows": int(conn.execute(
                "SELECT COUNT(*) FROM Spectra WHERE spectrum_id=?",
                (BAD_AO_SPECTRUM_ID,),
            ).fetchone()[0]),
            "oxidized_compounds": int(conn.execute(
                "SELECT COUNT(*) FROM Compounds WHERE inchikey=?", (OXIDIZED_IK,)
            ).fetchone()[0]),
            "irgafos_compounds": int(conn.execute(
                "SELECT COUNT(*) FROM Compounds WHERE inchikey=?", (IRGAFOS_IK,)
            ).fetchone()[0]),
            "source_tiers": dict(conn.execute(
                f"SELECT coalesce({source_expr}, '(missing)'), COUNT(*) "
                "FROM Spectra GROUP BY 1"
            ).fetchall()),
        }
    finally:
        conn.close()


def _add_column(conn: sqlite3.Connection, table: str, definition: str) -> None:
    name = definition.split()[0]
    if name not in _columns(conn, table):
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {definition}")


def _ensure_performance_indexes(path: Path) -> bool:
    """Add reconstructible indexes without changing semantic row data."""
    conn = sqlite3.connect(str(path))
    conn.execute("PRAGMA busy_timeout=60000")
    try:
        existed = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='index' "
            "AND name='ix_c_name_formula_norm'"
        ).fetchone() is not None
        conn.execute(
            "CREATE INDEX IF NOT EXISTS ix_c_name_formula_norm "
            "ON Compounds(lower(trim(compound_name)), molecular_formula)"
        )
        conn.commit()
        return not existed
    finally:
        conn.close()


def _update_compound_role(
    conn: sqlite3.Connection, inchikey: str, canonical_name: str
) -> int:
    columns = _columns(conn, "Compounds")
    assignments: list[str] = []
    values: list[object] = []
    for column, value in (
        ("compound_name", canonical_name),
        ("classification", "other"),
        ("compound_class", "other"),
        ("classification_v3", "other"),
        ("classification_source", "curated_compound_role"),
        ("classification_confidence", 1.0),
    ):
        if column in columns:
            assignments.append(f"{column}=?")
            values.append(value)
    values.append(inchikey)
    cursor = conn.execute(
        f"UPDATE Compounds SET {', '.join(assignments)} WHERE inchikey=?", values
    )
    return int(cursor.rowcount)


def migrate(path: Path, dry_run: bool = False) -> dict[str, object]:
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    _quick_check(path)
    before = audit(path)
    expected_literature = _expected_literature_count(path)
    if before["literature_rows"] != expected_literature:
        raise RuntimeError(
            f"Unexpected {LITERATURE_DATASET} count in {path.name}: "
            f"{before['literature_rows']} (expected {expected_literature})"
        )
    if before["oxidized_compounds"] != 1 or before["irgafos_compounds"] != 1:
        raise RuntimeError(f"Unexpected Irgafos compound cardinality: {before}")
    if before["bad_ao_rows"] not in (0, 1):
        raise RuntimeError(f"Unexpected malformed AO spectrum cardinality: {before}")
    if dry_run:
        return {"dry_run": True, "before": before}
    if _already_applied(path):
        index_added = _ensure_performance_indexes(path)
        _quick_check(path)
        return {
            "already_applied": True,
            "performance_index_added": index_added,
            "before": before,
        }

    backup = _clone_backup(path)
    conn = sqlite3.connect(str(path))
    conn.execute("PRAGMA busy_timeout=60000")
    try:
        conn.execute("BEGIN IMMEDIATE")
        _add_column(conn, "Spectra", "source_tier TEXT")
        _add_column(conn, "Spectra", "intensity_origin TEXT")
        _add_column(conn, "Spectra", "quality_flags_json TEXT")

        # Only the Small DB lacks historic tier metadata.  Do not invent
        # experimental provenance; keep it usable under a Medium confidence cap.
        conn.execute(
            "UPDATE Spectra SET source_tier='legacy_unresolved' "
            "WHERE source_tier IS NULL OR trim(source_tier)=''"
        )
        conn.execute(
            "UPDATE Spectra SET intensity_origin='unknown' "
            "WHERE intensity_origin IS NULL OR trim(intensity_origin)=''"
        )
        conn.execute(
            "UPDATE Spectra SET quality_flags_json='[]' "
            "WHERE quality_flags_json IS NULL OR trim(quality_flags_json)=''"
        )

        ao_count = _update_compound_role(
            conn, OXIDIZED_IK,
            "Tris(2,4-di-tert-butylphenyl) phosphate (oxidized Irgafos 168)",
        )
        irgafos_count = _update_compound_role(conn, IRGAFOS_IK, "Irgafos 168")
        if ao_count != 1 or irgafos_count != 1:
            raise RuntimeError(
                f"Compound update cardinality changed: AO={ao_count}, Irgafos={irgafos_count}"
            )

        literature_flag = json.dumps(
            [
                "literature_presence_only",
                "fragment_intensities_missing",
                "excluded_from_ordinary_spectral_search",
            ],
            ensure_ascii=False,
        )
        lit_cursor = conn.execute(
            "UPDATE Spectra SET source_tier='quarantined', intensity_origin='default', "
            "quality_flags_json=? WHERE dataset_source=?",
            (literature_flag, LITERATURE_DATASET),
        )
        if int(lit_cursor.rowcount) != expected_literature:
            raise RuntimeError(
                f"Literature quarantine touched {lit_cursor.rowcount}; expected {expected_literature}"
            )

        bad_ao_flag = json.dumps(
            [
                "precursor_formula_adduct_mismatch",
                "excluded_from_ordinary_spectral_search",
            ],
            ensure_ascii=False,
        )
        bad_cursor = conn.execute(
            "UPDATE Spectra SET source_tier='quarantined', quality_flags_json=? "
            "WHERE spectrum_id=?",
            (bad_ao_flag, BAD_AO_SPECTRUM_ID),
        )
        if int(bad_cursor.rowcount) != int(before["bad_ao_rows"]):
            raise RuntimeError("Malformed AO quarantine cardinality changed")

        conn.execute(
            "CREATE INDEX IF NOT EXISTS ix_s_source_tier ON Spectra(source_tier)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS ix_s_dataset_source ON Spectra(dataset_source)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS ix_c_name_formula_norm "
            "ON Compounds(lower(trim(compound_name)), molecular_formula)"
        )
        conn.execute(
            "CREATE TABLE IF NOT EXISTS LibraryQualityMigrations ("
            "migration_id TEXT PRIMARY KEY, applied_at TEXT NOT NULL, "
            "backup_path TEXT NOT NULL, details_json TEXT NOT NULL)"
        )
        details = {
            "literature_quarantined": expected_literature,
            "bad_ao_quarantined": int(before["bad_ao_rows"]),
            "oxidized_irgafos_compounds_corrected": ao_count,
            "irgafos_compounds_corrected": irgafos_count,
        }
        conn.execute(
            "INSERT INTO LibraryQualityMigrations VALUES (?,?,?,?)",
            (
                MIGRATION_ID,
                time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                str(backup),
                json.dumps(details, ensure_ascii=False, sort_keys=True),
            ),
        )
        conn.execute("COMMIT")
    except Exception:
        try:
            conn.execute("ROLLBACK")
        except Exception:
            pass
        conn.close()
        raise
    conn.close()
    _quick_check(path)
    after = audit(path)
    return {"backup": str(backup), "before": before, "after": after}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db", action="append", default=None,
        help="Database path (repeatable); defaults to Full and Small",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    project = Path(__file__).resolve().parent.parent
    paths = [Path(value) for value in (args.db or [
        project / "data/local_library.db",
        project / "data/local_library_small.db",
    ])]
    for path in paths:
        result = migrate(path, dry_run=args.dry_run)
        print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
