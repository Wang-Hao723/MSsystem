#!/usr/bin/env python3
"""Audit an explicitly supplied validation truth/event pair before scoring."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.validation_anchor import (  # noqa: E402
    audit_validation_anchor_batch,
    require_metric_eligible_anchors,
)


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _truth_records(payload: Any) -> list[Mapping[str, Any]]:
    records = payload.get("records") if isinstance(payload, Mapping) else payload
    if not isinstance(records, list):
        raise ValueError("truth JSON must be a list or an object with records")
    return records


def _events_by_record(payload: Any) -> dict[str, list[Mapping[str, Any]]]:
    if isinstance(payload, Mapping):
        source = payload.get("events_by_record", payload)
        if not isinstance(source, Mapping):
            raise ValueError("events_by_record must be an object")
        output = {}
        for key, events in source.items():
            if not isinstance(events, list):
                raise ValueError(f"events for {key} must be a list")
            output[str(key)] = events
        return output
    if not isinstance(payload, list):
        raise ValueError("event JSON must be a list or keyed object")
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for index, event in enumerate(payload):
        if not isinstance(event, Mapping):
            raise ValueError(f"event {index} is not an object")
        record_id = str(event.get("record_id") or "").strip()
        if not record_id:
            raise ValueError(f"event {index} has no record_id")
        grouped[record_id].append(
            {key: value for key, value in event.items() if key != "record_id"}
        )
    return dict(grouped)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--require-all-valid",
        action="store_true",
        help="Return a non-zero status while any anchor needs manual review.",
    )
    args = parser.parse_args()
    truth_path = args.truth.resolve()
    event_path = args.events.resolve()
    result = audit_validation_anchor_batch(
        _truth_records(_load(truth_path)),
        _events_by_record(_load(event_path)),
    )
    result["inputs"] = {
        "truth": {"path": str(truth_path), "sha256": _sha256(truth_path)},
        "events": {"path": str(event_path), "sha256": _sha256(event_path)},
    }
    result["inference_called"] = False
    result["retrieval_called"] = False
    result["ranking_called"] = False
    destination = args.output.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, destination)
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print(f"result={destination}")
    if args.require_all_valid:
        require_metric_eligible_anchors(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
