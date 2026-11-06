#!/usr/bin/env python3
"""Apply the fail-closed sealed-test export gate to V2 rules."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.export import export_empirically_validated_rules


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sealed-test-rules", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    tested = json.loads(args.sealed_test_rules.read_text(encoding="utf-8"))
    result = export_empirically_validated_rules(tested)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
