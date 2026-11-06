#!/usr/bin/env python3
"""Install the registry-bound V14 authentic-standard MS/MS artifact."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from governance.authentic_anchor_registry import (  # noqa: E402
    install_active_authentic_anchor_library,
)


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        required=True,
        help="Explicit source artifact; directories are never scanned.",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=PROJECT_ROOT / "config" / "authentic_anchor_registry.json",
    )
    parser.add_argument(
        "--receipt",
        type=Path,
        default=(
            PROJECT_ROOT
            / "docs"
            / "audits"
            / "V14_AUTHENTIC_ANCHOR_INSTALL_RECEIPT_20260825.json"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    receipt = install_active_authentic_anchor_library(
        args.source,
        args.registry,
        project_root=PROJECT_ROOT,
    )
    receipt["installed_at"] = datetime.now(timezone.utc).isoformat()
    receipt["registry_path"] = str(args.registry.resolve())
    _write_json_atomic(args.receipt.resolve(), receipt)
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
