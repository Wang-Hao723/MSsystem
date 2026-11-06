#!/usr/bin/env python3
"""Build exact-full-InChIKey display aliases for unresolved library names.

The production ``compound_name`` column is source provenance and is never
overwritten.  This script reads it in SQLite read-only mode, asks PubChem only
for records whose source label is missing/identifier-like/unsafe, verifies that
PubChem returned the exact requested *full* InChIKey, and writes a separate
runtime alias file.  Curated aliases continue to take precedence.

PubChem requests are sent in small POST batches below the published request
rate.  A durable cache makes the operation resumable and auditable.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import requests


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.compound_names import (  # noqa: E402
    alias_quality,
    load_compound_name_aliases,
    resolve_compound_name,
)


FULL_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
FAKE_LOOKUP_KEY = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
PUG_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/"
    "Title,IUPACName,MolecularFormula,MonoisotopicMass,InChI,InChIKey/JSON"
)
USER_AGENT = "Mssystem-v7-compound-name-audit/2026.07"
DEFAULT_CACHE = (
    ROOT.parents[1]
    / "Mssystem_shared_data"
    / "migration_ledgers"
    / "compound_name_pubchem_cache_20260731.json"
)
DEFAULT_OUTPUT = ROOT / "matching" / "compound_name_aliases_pubchem.json"


def _clean(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "null", "na", "n/a"} else text


def _full_key(value: Any) -> str:
    text = _clean(value).upper()
    return text if FULL_INCHIKEY_RE.fullmatch(text) else ""


def _chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _connect_readonly(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)


def _source_label_needs_alias(raw_name: Any) -> tuple[bool, str, str]:
    """Classify a source label without allowing its real key to hit an alias."""

    display, status = resolve_compound_name(raw_name, FAKE_LOOKUP_KEY)
    needs_alias = bool(
        status in {
            "missing",
            "identifier_unresolved",
            "code_identifier_unresolved",
            "computational_candidate",
        }
        or status.startswith("unsafe_")
        or display.startswith("Unresolved compound (")
        or display == "Unresolved library compound"
    )
    return needs_alias, display, status


def _local_alias_score(name: str) -> tuple[int, int, int, str]:
    """Prefer concise source synonyms when PubChem has no usable title."""

    return (
        int(len(name) <= 80),
        int(len(name) <= 40),
        -len(name),
        name.casefold(),
    )


def collect_targets(
    database: Path,
    existing_aliases: dict[str, str] | None = None,
) -> tuple[set[str], dict[str, str], dict[str, int]]:
    """Return verified unresolved keys and safe same-key local fallbacks."""

    existing_aliases = existing_aliases or {}
    target_keys: set[str] = set()
    local_candidates: dict[str, set[str]] = defaultdict(set)
    status_counts: dict[str, int] = defaultdict(int)

    with _connect_readonly(database) as connection:
        columns = {
            row[1]
            for row in connection.execute("PRAGMA table_info(Compounds)")
        }
        verified_expression = (
            "COALESCE(structure_identity_verified, 0)"
            if "structure_identity_verified" in columns
            else "1"
        )
        for raw_name, raw_key, verified, row_count in connection.execute(
            f"""
            SELECT compound_name, inchikey, {verified_expression}, COUNT(*)
            FROM Compounds
            WHERE inchikey IS NOT NULL
            GROUP BY compound_name, inchikey, {verified_expression}
            """
        ):
            key = _full_key(raw_key)
            if not key or not bool(verified):
                continue
            needs_alias, display, status = _source_label_needs_alias(raw_name)
            if needs_alias:
                status_counts[status] += int(row_count)
                if key not in existing_aliases:
                    target_keys.add(key)
                continue
            safe, _reason = alias_quality(display)
            if safe:
                local_candidates[key].add(display)

    selected_local = {
        key: max(names, key=_local_alias_score)
        for key, names in local_candidates.items()
        if names
    }
    return target_keys, selected_local, dict(status_counts)


def _request_with_retry(
    session: requests.Session,
    keys: list[str],
    *,
    max_attempts: int = 6,
) -> requests.Response:
    delay = 1.0
    last_response: requests.Response | None = None
    for attempt in range(1, max_attempts + 1):
        response = session.post(
            PUG_URL,
            data={"inchikey": ",".join(keys)},
            timeout=120,
        )
        last_response = response
        if response.status_code not in {429, 500, 502, 503, 504}:
            return response
        if attempt < max_attempts:
            time.sleep(delay)
            delay = min(delay * 2.0, 30.0)
    assert last_response is not None
    return last_response


def fetch_keys(
    session: requests.Session,
    keys: list[str],
    request_log: list[dict[str, Any]],
    *,
    request_interval: float,
) -> dict[str, dict[str, Any] | None]:
    """Fetch a key batch, bisecting 404s and verifying returned full keys."""

    if not keys:
        return {}
    started = time.monotonic()
    response = _request_with_retry(session, keys)
    request_log.append(
        {
            "keys": keys,
            "key_count": len(keys),
            "status_code": response.status_code,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "throttling_control": response.headers.get(
                "X-Throttling-Control", ""
            ),
        }
    )
    if response.status_code == 404 and len(keys) > 1:
        middle = len(keys) // 2
        time.sleep(request_interval)
        left = fetch_keys(
            session,
            keys[:middle],
            request_log,
            request_interval=request_interval,
        )
        time.sleep(request_interval)
        left.update(
            fetch_keys(
                session,
                keys[middle:],
                request_log,
                request_interval=request_interval,
            )
        )
        return left
    if response.status_code == 404:
        return {keys[0]: None}
    response.raise_for_status()

    result: dict[str, dict[str, Any] | None] = {key: None for key in keys}
    properties = (
        response.json().get("PropertyTable", {}).get("Properties", [])
    )
    for item in properties:
        returned_key = _full_key(item.get("InChIKey"))
        if returned_key in result:
            result[returned_key] = item
    return result


def _safe_pubchem_name(
    key: str, properties: dict[str, Any] | None
) -> tuple[str, str] | None:
    if not properties or _full_key(properties.get("InChIKey")) != key:
        return None
    for field in ("Title", "IUPACName"):
        candidate = _clean(properties.get(field))
        safe, _reason = alias_quality(candidate)
        if safe:
            return candidate, f"PubChem_{field}"
    return None


def build_alias_payload(
    *,
    target_keys: set[str],
    local_candidates: dict[str, str],
    properties_by_key: dict[str, dict[str, Any] | None],
    database: Path,
    status_counts: dict[str, int],
) -> dict[str, Any]:
    aliases: dict[str, str] = {}
    provenance: dict[str, str] = {}
    for key in sorted(target_keys):
        pubchem = _safe_pubchem_name(key, properties_by_key.get(key))
        if pubchem is not None:
            aliases[key], provenance[key] = pubchem
        elif key in local_candidates:
            aliases[key] = local_candidates[key]
            provenance[key] = "local_exact_full_inchikey_synonym"

    source_counts: dict[str, int] = defaultdict(int)
    for source in provenance.values():
        source_counts[source] += 1
    return {
        "schema_version": 1,
        "matching_rule": "exact full InChIKey only",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "database": str(database.resolve()),
        "source": "PubChem PUG REST plus safe local exact-full-InChIKey synonym fallback",
        "safety_rules": [
            "only structure_identity_verified=1 source records are eligible",
            "PubChem response full InChIKey must exactly equal requested full InChIKey",
            "identifier/noise/candidate metadata is rejected by alias_quality",
            "raw source compound_name is retained separately for provenance",
        ],
        "target_full_inchikeys": len(target_keys),
        "resolved_aliases": len(aliases),
        "unresolved_after_enrichment": len(target_keys) - len(aliases),
        "alias_source_counts": dict(sorted(source_counts.items())),
        "source_status_rows_targeted": dict(
            sorted(status_counts.items(), key=lambda item: (-item[1], item[0]))
        ),
        "aliases": aliases,
        "provenance": provenance,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=ROOT.parents[1] / "Mssystem_shared_data" / "local_library.db",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--request-interval", type=float, default=0.35)
    parser.add_argument(
        "--max-new-keys",
        type=int,
        default=0,
        help="Testing/resume aid; 0 fetches every currently missing key.",
    )
    args = parser.parse_args()
    if not args.db.is_file():
        parser.error(f"database does not exist: {args.db}")
    if args.batch_size < 1 or args.batch_size > 100:
        parser.error("--batch-size must be between 1 and 100")
    if args.request_interval < 0.2:
        parser.error("--request-interval must be >= 0.2 seconds")

    existing_aliases = load_compound_name_aliases()
    target_keys, local_candidates, status_counts = collect_targets(
        args.db, existing_aliases
    )
    cache_payload: dict[str, Any] = {}
    if args.cache.exists():
        try:
            loaded = json.loads(args.cache.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                cache_payload = loaded
        except (OSError, ValueError):
            cache_payload = {}
    properties_by_key: dict[str, dict[str, Any] | None] = dict(
        cache_payload.get("properties_by_inchikey") or {}
    )
    missing = sorted(target_keys - properties_by_key.keys())
    if args.max_new_keys > 0:
        missing = missing[: args.max_new_keys]

    print(
        json.dumps(
            {
                "verified_unresolved_full_inchikeys": len(target_keys),
                "safe_local_fallbacks": sum(
                    key in local_candidates for key in target_keys
                ),
                "cache_entries": len(properties_by_key),
                "keys_to_fetch": len(missing),
            },
            ensure_ascii=False,
        ),
        flush=True,
    )

    session = requests.Session()
    session.headers.update(
        {"User-Agent": USER_AGENT, "Accept": "application/json"}
    )
    request_log: list[dict[str, Any]] = []
    for batch_number, batch in enumerate(
        _chunks(missing, args.batch_size), start=1
    ):
        fetched = fetch_keys(
            session,
            batch,
            request_log,
            request_interval=args.request_interval,
        )
        # A 200 response can omit an individual identifier. Retry only those
        # omissions singly so a valid but unusual record is not cached as a miss.
        omitted = [key for key, value in fetched.items() if value is None]
        for key in omitted:
            time.sleep(args.request_interval)
            fetched[key] = fetch_keys(
                session,
                [key],
                request_log,
                request_interval=args.request_interval,
            ).get(key)
        properties_by_key.update(fetched)
        _atomic_json(
            args.cache,
            {
                "schema_version": 1,
                "updated_at": datetime.now(timezone.utc).isoformat(),
                "api_url": PUG_URL,
                "properties_by_inchikey": properties_by_key,
                "request_count": int(cache_payload.get("request_count") or 0)
                + len(request_log),
                "latest_request_log": request_log[-100:],
            },
        )
        resolved = sum(
            bool(properties_by_key.get(key)) for key in target_keys
        )
        print(
            f"pubchem_batch={batch_number} fetched={len(batch)} "
            f"cached={len(properties_by_key)} resolved={resolved}",
            flush=True,
        )
        time.sleep(args.request_interval)

    payload = build_alias_payload(
        target_keys=target_keys,
        local_candidates=local_candidates,
        properties_by_key=properties_by_key,
        database=args.db,
        status_counts=status_counts,
    )
    _atomic_json(args.output, payload)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "target_full_inchikeys": payload["target_full_inchikeys"],
                "resolved_aliases": payload["resolved_aliases"],
                "unresolved_after_enrichment": payload[
                    "unresolved_after_enrichment"
                ],
                "alias_source_counts": payload["alias_source_counts"],
            },
            ensure_ascii=False,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
