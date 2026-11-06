#!/usr/bin/env python3
"""Read-only audit of production compound-name resolution quality.

The audit never edits the spectral library.  It reports how every distinct
``(compound_name, full InChIKey)`` pair will be presented by the runtime,
including unresolved registry identifiers, compact development/library codes,
verified aliases and source-title cleanups.
"""

from __future__ import annotations

import argparse
import collections
import importlib.util
import json
import re
import sqlite3
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
_NAMES_SPEC = importlib.util.spec_from_file_location(
    "compound_names_audit", ROOT / "matching" / "compound_names.py"
)
if _NAMES_SPEC is None or _NAMES_SPEC.loader is None:
    raise RuntimeError("Cannot load matching/compound_names.py")
_NAMES_MODULE = importlib.util.module_from_spec(_NAMES_SPEC)
_NAMES_SPEC.loader.exec_module(_NAMES_MODULE)

resolve_compound_name = _NAMES_MODULE.resolve_compound_name
_CODE_PREFIX = re.compile(r"^[A-Z]+")


def _connect_readonly(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)


def audit_names(database: Path, *, example_limit: int = 30) -> dict[str, Any]:
    status_rows: collections.Counter[str] = collections.Counter()
    status_pairs: collections.Counter[str] = collections.Counter()
    code_prefix_rows: collections.Counter[str] = collections.Counter()
    code_examples: list[dict[str, Any]] = []
    ao_source_rows: list[dict[str, Any]] = []
    total_rows = 0
    total_pairs = 0

    with _connect_readonly(database) as conn:
        for raw_name, inchikey, n_rows in conn.execute(
            """
            SELECT compound_name, inchikey, COUNT(*)
            FROM Compounds
            GROUP BY compound_name, inchikey
            """
        ):
            row_count = int(n_rows)
            total_rows += row_count
            total_pairs += 1
            display_name, status = resolve_compound_name(raw_name, inchikey)
            status_rows[status] += row_count
            status_pairs[status] += 1

            if status == "code_identifier_unresolved":
                source_name = str(raw_name or "").strip()
                prefix_match = _CODE_PREFIX.match(source_name)
                prefix = prefix_match.group(0) if prefix_match else "<none>"
                code_prefix_rows[prefix] += row_count
                if len(code_examples) < example_limit:
                    code_examples.append(
                        {
                            "source_name": source_name,
                            "inchikey": str(inchikey or ""),
                            "display_name": display_name,
                            "compound_rows": row_count,
                        }
                    )

        for row in conn.execute(
            """
            SELECT compound_name, inchikey, molecular_formula, smiles,
                   classification_v3, classification_source
            FROM Compounds
            WHERE compound_name GLOB 'AO[0-9]*'
            ORDER BY compound_name
            """
        ):
            source_name, inchikey, formula, smiles, class_v3, class_source = row
            display_name, status = resolve_compound_name(source_name, inchikey)
            ao_source_rows.append(
                {
                    "source_name": source_name,
                    "display_name": display_name,
                    "resolution_status": status,
                    "inchikey": inchikey,
                    "formula": formula,
                    "has_structure": bool(str(smiles or "").strip()),
                    "stored_classification_v3": class_v3,
                    "stored_classification_source": class_source,
                }
            )

    return {
        "database": str(database.resolve()),
        "read_only": True,
        "compound_rows": total_rows,
        "distinct_name_inchikey_pairs": total_pairs,
        "resolution_status_rows": dict(status_rows.most_common()),
        "resolution_status_pairs": dict(status_pairs.most_common()),
        "unresolved_compact_code_rows": status_rows.get(
            "code_identifier_unresolved", 0
        ),
        "unresolved_compact_code_prefixes": dict(code_prefix_rows.most_common(50)),
        "unresolved_compact_code_examples": code_examples,
        "ao_source_records": ao_source_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=ROOT.parents[1] / "Mssystem_shared_data" / "local_library.db",
    )
    parser.add_argument("--examples", type=int, default=30)
    args = parser.parse_args()
    if not args.db.is_file():
        parser.error(f"database does not exist: {args.db}")
    payload = audit_names(args.db, example_limit=max(0, args.examples))
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
