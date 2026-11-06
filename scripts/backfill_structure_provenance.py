#!/usr/bin/env python3
"""Backfill structure-derived InChIKeys and explicit verification provenance.

``structure_identity_verified`` means that the stored structural identifier was
derived reproducibly from a parseable SMILES/InChI.  It does *not* claim that an
external authority authenticated the submitted compound identity.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sqlite3
import sys
from collections import Counter
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.settings import DB_PATH, LIBRARY_WRITE_ENABLED  # noqa: E402
from matching.library_quality import inchi_metadata, smiles_metadata  # noqa: E402
from rdkit import rdBase  # noqa: E402


FULL_KEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
MIGRATION_ID = "v7_20260729_structure_provenance_backfill_v1"


def _text(value: Any) -> str:
    return str(value or "").strip()


def _key(value: Any) -> str:
    text = _text(value).upper()
    return text if FULL_KEY.fullmatch(text) else ""


@lru_cache(maxsize=250_000)
def _derive(source_kind: str, structure: str) -> tuple[str, str]:
    metadata = (
        smiles_metadata(structure)
        if source_kind == "smiles"
        else inchi_metadata(structure)
    )
    return str(metadata.get("status") or "unverifiable"), _key(
        metadata.get("inchikey")
    )


def _connection(path: Path, writable: bool) -> sqlite3.Connection:
    if writable:
        connection = sqlite3.connect(str(path), timeout=120)
    else:
        connection = sqlite3.connect(
            f"file:{path.resolve()}?mode=ro", uri=True, timeout=120
        )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout=120000")
    return connection


def _new_values(row: sqlite3.Row) -> tuple[dict[str, Any], str]:
    smiles = _text(row["smiles"])
    inchi = _text(row["inchi"])
    source_kind = "smiles" if smiles else "inchi" if inchi else "none"
    structure = smiles or inchi
    old_key = _key(row["inchikey"])
    old_reported = _key(row["reported_inchikey"])
    if source_kind == "none":
        values = {
            "inchikey": row["inchikey"],
            "reported_inchikey": row["reported_inchikey"],
            "structure_identity_verified": 0,
            "structure_identity_source": "no_parseable_structure",
        }
        return values, "no_structure"

    status, derived_key = _derive(source_kind, structure)
    if status != "valid" or not derived_key:
        values = {
            "inchikey": row["inchikey"],
            "reported_inchikey": row["reported_inchikey"],
            "structure_identity_verified": 0,
            "structure_identity_source": (
                f"{source_kind}_{status}:RDKit-{rdBase.rdkitVersion}"
            ),
        }
        return values, f"{source_kind}_{status}"

    if not old_key:
        relationship = "derived_filled_missing_key"
    elif old_key == derived_key:
        relationship = "reported_full_key_exact"
    elif old_key[:14] == derived_key[:14]:
        relationship = "reported_connectivity_only_match"
    else:
        relationship = "reported_connectivity_mismatch"
    reported = old_reported or old_key or None
    values = {
        "inchikey": derived_key,
        "reported_inchikey": reported,
        "structure_identity_verified": 1,
        "structure_identity_source": (
            f"{source_kind}_derived:RDKit-{rdBase.rdkitVersion};"
            f"{relationship};not_external_authentication"
        ),
    }
    return values, relationship


def iter_updates(
    connection: sqlite3.Connection,
    *,
    limit: int = 0,
) -> Iterator[tuple[dict[str, Any], str]]:
    cursor = connection.execute(
        """
        SELECT id, compound_name, smiles, inchi, inchikey, reported_inchikey,
               structure_identity_verified, structure_identity_source
        FROM Compounds ORDER BY id
        """
    )
    emitted = 0
    for row in cursor:
        new, relationship = _new_values(row)
        old = {
            "inchikey": row["inchikey"],
            "reported_inchikey": row["reported_inchikey"],
            "structure_identity_verified": int(
                row["structure_identity_verified"] or 0
            ),
            "structure_identity_source": row["structure_identity_source"],
        }
        if all(old[key] == new[key] for key in new):
            yield {
                "changed": False,
                "compound_id": int(row["id"]),
                "relationship": relationship,
            }, relationship
            continue
        update = {
            "changed": True,
            "compound_id": int(row["id"]),
            "compound_name": str(row["compound_name"] or ""),
            "relationship": relationship,
            "old": old,
            "new": new,
        }
        yield update, relationship
        emitted += 1
        if limit and emitted >= limit:
            return


def default_ledger(path: Path) -> Path:
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    return (
        path.parent
        / "migration_ledgers"
        / f"structure_provenance_{path.stem}_{stamp}.jsonl.gz"
    )


def write_ledger(
    path: Path,
    db_path: Path,
    updates: Iterator[tuple[dict[str, Any], str]],
) -> tuple[int, Counter]:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter = Counter()
    changed = 0
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        header = {
            "type": "structure_provenance_backfill",
            "schema_version": 1,
            "migration_id": MIGRATION_ID,
            "database": str(db_path.resolve()),
            "created_at": datetime.now().astimezone().isoformat(),
            "rdkit_version": rdBase.rdkitVersion,
            "verification_semantics": (
                "structure-derived internal consistency, not external identity "
                "authentication"
            ),
        }
        handle.write(json.dumps({"header": header}, ensure_ascii=False) + "\n")
        for update, relationship in updates:
            counts[relationship] += 1
            counts["rows_scanned"] += 1
            if not update.get("changed"):
                counts["unchanged"] += 1
                continue
            handle.write(json.dumps(update, ensure_ascii=False) + "\n")
            changed += 1
        counts["eligible_updates"] = changed
    return changed, counts


def apply_ledger(connection: sqlite3.Connection, path: Path) -> int:
    changed = 0
    connection.execute("BEGIN IMMEDIATE")
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            header = json.loads(next(handle)).get("header", {})
            if header.get("type") != "structure_provenance_backfill":
                raise ValueError("Not a structure-provenance ledger")
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                old = row["old"]
                new = row["new"]
                cursor = connection.execute(
                    """
                    UPDATE Compounds
                    SET inchikey=?, reported_inchikey=?,
                        structure_identity_verified=?,
                        structure_identity_source=?
                    WHERE id=?
                      AND COALESCE(inchikey,'')=COALESCE(?,'')
                      AND COALESCE(reported_inchikey,'')=COALESCE(?,'')
                      AND COALESCE(structure_identity_verified,0)=?
                      AND COALESCE(structure_identity_source,'')=COALESCE(?,'')
                    """,
                    (
                        new["inchikey"],
                        new["reported_inchikey"],
                        new["structure_identity_verified"],
                        new["structure_identity_source"],
                        row["compound_id"],
                        old["inchikey"],
                        old["reported_inchikey"],
                        old["structure_identity_verified"],
                        old["structure_identity_source"],
                    ),
                )
                if cursor.rowcount != 1:
                    raise RuntimeError(
                        "Concurrent-change guard failed for compound "
                        f"{row['compound_id']}"
                    )
                changed += 1
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS LibraryQualityMigrations (
                migration_id TEXT PRIMARY KEY, applied_at TEXT NOT NULL,
                backup_path TEXT NOT NULL, details_json TEXT NOT NULL
            )
            """
        )
        connection.execute(
            "INSERT OR REPLACE INTO LibraryQualityMigrations VALUES (?,?,?,?)",
            (
                MIGRATION_ID,
                datetime.now().astimezone().isoformat(),
                str(path.resolve()),
                json.dumps(
                    {
                        "changed": changed,
                        "rdkit_version": rdBase.rdkitVersion,
                        "semantics": "structure-derived; not external authentication",
                    },
                    sort_keys=True,
                ),
            ),
        )
        connection.commit()
        return changed
    except Exception:
        connection.rollback()
        raise


def rollback_ledger(connection: sqlite3.Connection, path: Path) -> int:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = json.loads(next(handle)).get("header", {})
        if header.get("type") != "structure_provenance_backfill":
            raise ValueError("Not a structure-provenance ledger")
        rows = [json.loads(line) for line in handle if line.strip()]
    changed = 0
    connection.execute("BEGIN IMMEDIATE")
    try:
        for row in rows:
            old, new = row["old"], row["new"]
            cursor = connection.execute(
                """
                UPDATE Compounds
                SET inchikey=?, reported_inchikey=?,
                    structure_identity_verified=?,
                    structure_identity_source=?
                WHERE id=?
                  AND COALESCE(inchikey,'')=COALESCE(?,'')
                  AND COALESCE(reported_inchikey,'')=COALESCE(?,'')
                  AND COALESCE(structure_identity_verified,0)=?
                  AND COALESCE(structure_identity_source,'')=COALESCE(?,'')
                """,
                (
                    old["inchikey"],
                    old["reported_inchikey"],
                    old["structure_identity_verified"],
                    old["structure_identity_source"],
                    row["compound_id"],
                    new["inchikey"],
                    new["reported_inchikey"],
                    new["structure_identity_verified"],
                    new["structure_identity_source"],
                ),
            )
            changed += max(0, int(cursor.rowcount))
        connection.commit()
        return changed
    except Exception:
        connection.rollback()
        raise


def summary(connection: sqlite3.Connection) -> dict[str, int]:
    row = connection.execute(
        """
        SELECT COUNT(*),
               SUM(COALESCE(structure_identity_verified,0)=1),
               SUM(structure_identity_source IS NOT NULL
                   AND trim(structure_identity_source)!=''),
               SUM(reported_inchikey IS NOT NULL
                   AND trim(reported_inchikey)!=''),
               SUM(inchikey IS NULL OR trim(inchikey)='')
        FROM Compounds
        """
    ).fetchone()
    return {
        "total_compounds": int(row[0] or 0),
        "structure_identity_verified": int(row[1] or 0),
        "with_structure_identity_source": int(row[2] or 0),
        "with_reported_inchikey": int(row[3] or 0),
        "missing_inchikey": int(row[4] or 0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, default=DB_PATH)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--rollback", type=Path)
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    db_path = args.db.expanduser().resolve()
    if (args.apply or args.rollback) and not LIBRARY_WRITE_ENABLED:
        raise SystemExit("Library writes are disabled by MSSYSTEM_LIBRARY_WRITE_ENABLED")
    if args.apply and args.rollback:
        raise SystemExit("Choose either --apply or --rollback")

    if args.rollback:
        connection = _connection(db_path, True)
        changed = rollback_ledger(connection, args.rollback.expanduser().resolve())
        result = {"rolled_back": changed, "after": summary(connection)}
        connection.close()
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return

    before_connection = _connection(db_path, args.apply)
    before = summary(before_connection)
    if args.apply:
        ledger = (args.ledger or default_ledger(db_path)).expanduser().resolve()
        changed_expected, counts = write_ledger(
            ledger,
            db_path,
            iter_updates(before_connection, limit=max(0, args.limit)),
        )
        changed = apply_ledger(before_connection, ledger)
        if changed != changed_expected:
            raise RuntimeError(
                f"Expected {changed_expected} updates but applied {changed}"
            )
        check = before_connection.execute("PRAGMA quick_check").fetchone()
        result = {
            "mode": "apply",
            "database": str(db_path),
            "before": before,
            "audit": dict(sorted(counts.items())),
            "changed": changed,
            "rollback_ledger": str(ledger),
            "quick_check": check[0] if check else None,
            "after": summary(before_connection),
        }
    else:
        counts: Counter = Counter()
        examples = []
        for update, relationship in iter_updates(
            before_connection, limit=max(0, args.limit)
        ):
            counts[relationship] += 1
            counts["rows_scanned"] += 1
            if update.get("changed"):
                counts["eligible_updates"] += 1
                if len(examples) < 20:
                    examples.append(update)
        result = {
            "mode": "dry_run",
            "database": str(db_path),
            "before": before,
            "audit": dict(sorted(counts.items())),
            "examples": examples,
        }
    before_connection.close()
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
