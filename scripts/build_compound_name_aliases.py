#!/usr/bin/env python3
"""Build conservative identifier display aliases from verified sources.

Only exact *full* InChIKey matches are accepted.  This intentionally avoids
connectivity-layer and fuzzy-name propagation because stereoisomers can share
the first InChIKey block while requiring different scientific names.
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from matching.compound_names import is_identifier_name, is_registry_labeled_name  # noqa: E402

MISSING = {"", "nan", "none", "null", "na", "n/a", "n.a."}
ADDUCT_SUFFIX = re.compile(
    r"\s+(?:\[M[^\]]+\][+-]?|M(?:[+-][A-Za-z0-9]+)+(?:[+-][A-Za-z0-9]+)*)\s*$",
    re.IGNORECASE,
)
IDENTIFIER_NAME = re.compile(
    r"^(?:NCGC|CHEMBL|SCHEMBL|CHEBI|HMDB|CID\s*\d|DB\d|MLS\d|DTXSID|NPASS)",
    re.IGNORECASE,
)
IDENTIFIER_ONLY = re.compile(r"^[A-Z]{2,}[A-Z0-9:_-]*\d[A-Z0-9:_-]*$")


def clean_name(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if text.lower() in MISSING:
        return ""
    return ADDUCT_SUFFIX.sub("", text).strip(" _!")


def is_readable_alias(value: object) -> bool:
    text = clean_name(value)
    return bool(
        text
        and len(text) <= 120
        and not IDENTIFIER_NAME.match(text)
        and not is_identifier_name(text)
        and not text.isdigit()
        and not (text == text.upper() and IDENTIFIER_ONLY.match(text))
        and not text.lower().startswith("suspect related to")
    )


def alias_score(value: str) -> tuple[int, int, int, str]:
    """Prefer concise natural/product names over long systematic descriptions."""
    return (
        int(len(value) <= 60),
        int(any(ch.isspace() for ch in value)),
        -len(value),
        value.casefold(),
    )


def connect_readonly(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)


def build_aliases(
    current_db: Path, legacy_db: Path, curated_path: Path | None = None
) -> dict[str, str]:
    with connect_readonly(current_db) as current:
        raw_keys = {
            str(inchikey).strip().upper()
            for row in current.execute(
                "SELECT inchikey, compound_name FROM Compounds WHERE inchikey IS NOT NULL"
            )
            for inchikey, compound_name in [row]
            if inchikey and is_registry_labeled_name(compound_name)
        }

    aliases: dict[str, str] = {}
    with connect_readonly(legacy_db) as legacy:
        for inchikey, name in legacy.execute(
            "SELECT inchikey, compound_name FROM Compounds WHERE inchikey IS NOT NULL"
        ):
            key = str(inchikey or "").strip().upper()
            if key not in raw_keys or not is_readable_alias(name):
                continue
            cleaned = clean_name(name)
            previous = aliases.get(key)
            if previous is None or alias_score(cleaned) > alias_score(previous):
                aliases[key] = cleaned

    if curated_path and curated_path.exists():
        payload = json.loads(curated_path.read_text(encoding="utf-8"))
        curated = payload.get("aliases", payload) if isinstance(payload, dict) else {}
        for inchikey, name in (curated.items() if isinstance(curated, dict) else []):
            key = str(inchikey or "").strip().upper()
            if key in raw_keys and is_readable_alias(name):
                aliases[key] = clean_name(name)
    return dict(sorted(aliases.items()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-db", type=Path, required=True)
    parser.add_argument("--legacy-db", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--curated-aliases", type=Path,
        default=ROOT / "matching" / "compound_name_aliases_curated.json",
    )
    args = parser.parse_args()

    aliases = build_aliases(args.current_db, args.legacy_db, args.curated_aliases)
    payload = {
        "schema_version": 1,
        "matching_rule": "exact full InChIKey only",
        "sources": [args.legacy_db.name, args.curated_aliases.name],
        "aliases": aliases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp.replace(args.output)
    print(f"wrote {len(aliases)} exact aliases to {args.output}")


if __name__ == "__main__":
    main()
