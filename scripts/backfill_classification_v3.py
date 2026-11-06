#!/usr/bin/env python3
"""Audit and safely backfill missing ``Compounds.classification_v3`` values.

This migration is intentionally precision-first:

* existing v3 labels are never overwritten;
* placeholder/invalid structures are skipped;
* broad SMARTS classes are written only when an independently stored raw class
  normalizes to the same class;
* structure-only writes are limited to high-specificity SMARTS classes;
* ``--apply`` writes a compressed rollback ledger before one atomic transaction.

The default action is a read-only dry run.  This tool does not use the
MS2Classifier to fill the database because its softmax output is not calibrated
as compound-level structural ground truth.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.settings import DB_PATH, DB_PATH_SMALL, LIBRARY_WRITE_ENABLED  # noqa: E402
from mssystem.classification import (  # noqa: E402
    RuleBasedClassifier,
    normalize_filename_class,
)
from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

PLACEHOLDERS = {"", "n/a", "na", "none", "null", "-"}
STRUCTURE_ONLY_CLASSES = {
    "flavonoid",
    "coumarin",
    "xanthone",
    "stilbene",
    "quinone",
    "phenolic_acid",
    "nucleoside",
    "sugar",
    "glycoside",
}


@dataclass(frozen=True)
class ClassificationUpdate:
    compound_id: int
    compound_name: str
    old_classification_v3: Optional[str]
    old_classification_source: Optional[str]
    old_classification_confidence: Optional[float]
    new_classification_v3: str
    new_classification_source: str
    new_classification_confidence: float
    matched_smarts: str
    raw_classification: str


def _connect(db_path: Path, writable: bool) -> sqlite3.Connection:
    if writable:
        connection = sqlite3.connect(str(db_path), timeout=60)
    else:
        connection = sqlite3.connect(
            f"file:{db_path.resolve()}?mode=ro", uri=True, timeout=60
        )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout=60000")
    return connection


def _valid_structure_text(value: object) -> bool:
    return (
        isinstance(value, str)
        and value.strip().lower() not in PLACEHOLDERS
    )


def classify_missing_rows(
    connection: sqlite3.Connection,
    *,
    limit: int = 0,
) -> tuple[list[ClassificationUpdate], Counter]:
    """Return eligible missing-only updates and transparent audit counters."""
    classifier = RuleBasedClassifier()
    sql = """
        SELECT id, compound_name, smiles, classification,
               classification_v3, classification_source,
               classification_confidence
        FROM Compounds
        WHERE classification_v3 IS NULL OR trim(classification_v3) = ''
        ORDER BY id
    """
    rows: Iterable[sqlite3.Row] = connection.execute(sql)
    updates: list[ClassificationUpdate] = []
    counters: Counter = Counter()
    for row in rows:
        counters["rows_missing_v3"] += 1
        smiles = row["smiles"]
        if not _valid_structure_text(smiles):
            counters["skipped_no_usable_structure"] += 1
            continue
        if Chem.MolFromSmiles(str(smiles)) is None:
            counters["skipped_invalid_smiles"] += 1
            continue
        detail = classifier.predict_detailed(smiles=str(smiles))
        predicted = str(detail.get("class") or "Unknown")
        if predicted == "Unknown":
            counters["skipped_no_smarts_class"] += 1
            continue

        raw = str(row["classification"] or "").strip()
        raw_normalized = normalize_filename_class(raw)
        if raw_normalized == predicted:
            source = "smarts_raw_agreement_v6"
            confidence = 0.99
            counters["eligible_smarts_raw_agreement"] += 1
        elif predicted in STRUCTURE_ONLY_CLASSES:
            source = "smarts_high_specific_v6"
            confidence = 0.98
            counters["eligible_high_specific_structure"] += 1
        else:
            counters[f"skipped_structure_only_broad:{predicted}"] += 1
            continue

        updates.append(
            ClassificationUpdate(
                compound_id=int(row["id"]),
                compound_name=str(row["compound_name"] or ""),
                old_classification_v3=row["classification_v3"],
                old_classification_source=row["classification_source"],
                old_classification_confidence=row["classification_confidence"],
                new_classification_v3=predicted,
                new_classification_source=source,
                new_classification_confidence=confidence,
                matched_smarts=str(detail.get("matched_smarts") or ""),
                raw_classification=raw,
            )
        )
        counters[f"new_class:{predicted}"] += 1
        if limit and len(updates) >= limit:
            counters["limited"] = 1
            break
    counters["eligible_updates"] = len(updates)
    return updates, counters


def _default_ledger_path(db_path: Path) -> Path:
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    directory = db_path.parent / "migration_ledgers"
    return directory / f"classification_v3_backfill_{db_path.stem}_{stamp}.jsonl.gz"


def write_ledger(path: Path, db_path: Path, updates: list[ClassificationUpdate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "type": "classification_v3_backfill",
        "schema_version": 1,
        "database": str(db_path.resolve()),
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "rows": len(updates),
    }
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps({"header": header}, ensure_ascii=False) + "\n")
        for update in updates:
            handle.write(json.dumps(asdict(update), ensure_ascii=False) + "\n")


def apply_updates(
    connection: sqlite3.Connection,
    updates: list[ClassificationUpdate],
) -> int:
    """Apply missing-only updates atomically and return the changed row count."""
    before = connection.total_changes
    connection.execute("BEGIN IMMEDIATE")
    try:
        connection.executemany(
            """
            UPDATE Compounds
            SET classification_v3 = ?,
                classification_source = ?,
                classification_confidence = ?
            WHERE id = ?
              AND (classification_v3 IS NULL OR trim(classification_v3) = '')
            """,
            [
                (
                    update.new_classification_v3,
                    update.new_classification_source,
                    update.new_classification_confidence,
                    update.compound_id,
                )
                for update in updates
            ],
        )
        changed = connection.total_changes - before
        if changed != len(updates):
            raise RuntimeError(
                f"Atomic safety check failed: expected {len(updates)} writes, got {changed}"
            )
        connection.commit()
        return changed
    except Exception:
        connection.rollback()
        raise


def rollback_ledger(connection: sqlite3.Connection, ledger_path: Path) -> int:
    """Restore only rows still carrying this migration's exact new values."""
    updates = []
    with gzip.open(ledger_path, "rt", encoding="utf-8") as handle:
        header_line = json.loads(next(handle))
        if header_line.get("header", {}).get("type") != "classification_v3_backfill":
            raise ValueError("Not a classification_v3 backfill ledger")
        for line in handle:
            if line.strip():
                updates.append(json.loads(line))
    before = connection.total_changes
    connection.execute("BEGIN IMMEDIATE")
    try:
        for row in updates:
            connection.execute(
                """
                UPDATE Compounds
                SET classification_v3 = ?,
                    classification_source = ?,
                    classification_confidence = ?
                WHERE id = ?
                  AND classification_v3 = ?
                  AND classification_source = ?
                """,
                (
                    row.get("old_classification_v3"),
                    row.get("old_classification_source"),
                    row.get("old_classification_confidence"),
                    row["compound_id"],
                    row["new_classification_v3"],
                    row["new_classification_source"],
                ),
            )
        changed = connection.total_changes - before
        connection.commit()
        return changed
    except Exception:
        connection.rollback()
        raise


def coverage(connection: sqlite3.Connection) -> dict:
    row = connection.execute(
        """
        SELECT COUNT(*) AS total,
               SUM(CASE WHEN classification_v3 IS NOT NULL
                         AND trim(classification_v3) != '' THEN 1 ELSE 0 END) AS with_v3,
               SUM(CASE WHEN smiles IS NOT NULL AND trim(smiles) != ''
                         AND lower(trim(smiles)) NOT IN
                             ('n/a','na','none','null','-')
                        THEN 1 ELSE 0 END) AS usable_structure_text
        FROM Compounds
        """
    ).fetchone()
    total = int(row["total"] or 0)
    with_v3 = int(row["with_v3"] or 0)
    return {
        "total": total,
        "with_classification_v3": with_v3,
        "missing_classification_v3": total - with_v3,
        "coverage_percent": round(100.0 * with_v3 / total, 3) if total else 0.0,
        "usable_structure_text": int(row["usable_structure_text"] or 0),
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
    if not db_path.exists():
        raise SystemExit(f"Database not found: {db_path}")
    if (args.apply or args.rollback) and not LIBRARY_WRITE_ENABLED:
        raise SystemExit("Library writes are disabled by MSSYSTEM_LIBRARY_WRITE_ENABLED")
    if args.apply and args.rollback:
        raise SystemExit("Choose either --apply or --rollback")

    if args.rollback:
        connection = _connect(db_path, writable=True)
        changed = rollback_ledger(connection, args.rollback.expanduser().resolve())
        print(json.dumps({"rolled_back": changed, "coverage": coverage(connection)}, indent=2))
        connection.close()
        return

    connection = _connect(db_path, writable=args.apply)
    before = coverage(connection)
    updates, counters = classify_missing_rows(connection, limit=max(0, args.limit))
    result = {
        "database": str(db_path),
        "mode": "apply" if args.apply else "dry_run",
        "before": before,
        "audit": dict(sorted(counters.items())),
    }
    if args.apply:
        ledger = (args.ledger or _default_ledger_path(db_path)).expanduser().resolve()
        write_ledger(ledger, db_path, updates)
        changed = apply_updates(connection, updates)
        result.update(
            {
                "changed": changed,
                "rollback_ledger": str(ledger),
                "after": coverage(connection),
            }
        )
    else:
        projected = before["with_classification_v3"] + len(updates)
        result["projected_coverage_percent"] = round(
            100.0 * projected / before["total"], 3
        )
    connection.close()
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
