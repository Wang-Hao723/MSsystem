#!/usr/bin/env python3
"""Quarantine generated aliases that cannot safely name a whole compound."""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.compound_names import alias_quality  # noqa: E402


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def audit_aliases(path: Path) -> tuple[dict[str, str], list[dict], Counter]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    aliases = payload.get("aliases", payload)
    if not isinstance(aliases, dict):
        raise ValueError("Alias file does not contain an aliases object")
    safe: dict[str, str] = {}
    quarantine = []
    counts: Counter = Counter()
    for inchikey, alias in sorted(aliases.items()):
        accepted, reason = alias_quality(alias)
        counts[reason] += 1
        if accepted:
            safe[str(inchikey)] = str(alias)
        else:
            quarantine.append(
                {
                    "inchikey": str(inchikey),
                    "alias": str(alias),
                    "reason": reason,
                }
            )
    return safe, quarantine, counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--aliases",
        type=Path,
        default=ROOT / "matching" / "compound_name_aliases.json",
    )
    parser.add_argument(
        "--quarantine",
        type=Path,
        default=ROOT / "matching" / "compound_name_aliases_quarantine.json",
    )
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    safe, quarantine, counts = audit_aliases(args.aliases)
    result = {
        "mode": "apply" if args.apply else "dry_run",
        "input_aliases": len(safe) + len(quarantine),
        "safe_aliases": len(safe),
        "quarantined_aliases": len(quarantine),
        "reason_counts": dict(sorted(counts.items())),
    }
    if args.apply:
        source_payload = json.loads(args.aliases.read_text(encoding="utf-8"))
        cleaned_payload = {
            key: value
            for key, value in source_payload.items()
            if key != "aliases"
        }
        cleaned_payload.update(
            {
                "schema_version": 2,
                "alias_quality_rule": (
                    "exact full InChIKey plus whole-compound-name safety gate"
                ),
                "audited_at": datetime.now(timezone.utc).isoformat(),
                "quarantined_aliases": len(quarantine),
                "aliases": safe,
            }
        )
        quarantine_payload = {
            "schema_version": 1,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "source": str(args.aliases.resolve()),
            "policy": (
                "Quarantined aliases remain provenance only and are never "
                "presented as canonical compound names."
            ),
            "reason_counts": dict(sorted(counts.items())),
            "aliases": quarantine,
        }
        _atomic_json(args.aliases, cleaned_payload)
        _atomic_json(args.quarantine, quarantine_payload)
        result["outputs"] = [
            str(args.aliases.resolve()),
            str(args.quarantine.resolve()),
        ]
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
