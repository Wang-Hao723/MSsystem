#!/usr/bin/env python3
"""Validate the screening catalog and audit exact-key spectral coverage."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.screening_catalog import (
    audit_library_coverage,
    load_screening_catalog,
    validate_screening_catalog,
)
from runtime.manifest import sha256_file, write_json_atomic


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the identity-first contaminant/pesticide catalog and "
            "inspect a library database in read-only immutable mode."
        )
    )
    parser.add_argument(
        "--database",
        type=Path,
        help="Optional local_library.db path for an exact full-InChIKey coverage audit.",
    )
    parser.add_argument(
        "--validate-structures",
        action="store_true",
        help="Recompute full InChIKeys from catalog SMILES with RDKit.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional JSON artifact path for a reproducible release audit.",
    )
    args = parser.parse_args()

    payload = load_screening_catalog()
    errors = validate_screening_catalog(
        payload,
        validate_structures=args.validate_structures,
    )
    if errors:
        print(json.dumps({"status": "invalid", "errors": errors}, ensure_ascii=False, indent=2))
        return 1
    result = (
        audit_library_coverage(args.database)
        if args.database
        else {
            "status": "valid",
            "catalog_id": payload.get("catalog_id"),
            "entries": len(payload.get("entries") or []),
            "regulatory_claims_encoded": False,
        }
    )
    result["status"] = "valid"
    result["audit_generated_at"] = datetime.now(timezone.utc).isoformat()
    catalog_path = PROJECT_ROOT / "matching" / "screening_catalog.json"
    result["catalog_sha256"] = sha256_file(catalog_path)
    if args.database:
        database = args.database.expanduser().resolve()
        stat = database.stat()
        result["database_fingerprint"] = {
            "path": str(database),
            "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
        manifest_path = PROJECT_ROOT / "config" / "runtime_manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            production = (manifest.get("library") or {}).get("production_full") or {}
        except (OSError, TypeError, ValueError):
            production = {}
        snapshot_matches = (
            str(production.get("path") or "") == str(database)
            and int(production.get("size_bytes") or -1) == stat.st_size
            and int(production.get("mtime_ns") or -1) == stat.st_mtime_ns
        )
        result["database_snapshot_id"] = (
            production.get("snapshot_id") if snapshot_matches else None
        )
        result["database_snapshot_binding"] = (
            "matched_runtime_manifest" if snapshot_matches else "unbound"
        )
    if args.output is not None:
        write_json_atomic(args.output, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
