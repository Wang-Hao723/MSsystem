#!/usr/bin/env python3
"""Create a fail-closed project calibration artifact from paired competitions."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from calibration.models import CalibrationProject, CompetitionPair
from calibration.project import calibrate_project, load_calibration_method


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--method",
        type=Path,
        default=PROJECT_ROOT / "config" / "calibration_method.json",
    )
    parser.add_argument("--require-calibrated", action="store_true")
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    project = CalibrationProject.model_validate(payload.get("project"))
    pairs = tuple(
        CompetitionPair.model_validate(item) for item in payload.get("pairs") or []
    )
    artifact = calibrate_project(project, load_calibration_method(args.method), pairs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(
            artifact.model_dump(mode="json"),
            handle,
            ensure_ascii=False,
            indent=2,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, args.output)
    print(
        json.dumps(
            {
                "status": artifact.status.value,
                "project_id": artifact.project.project_id,
                "competitions": artifact.readiness.competition_count,
                "strata": len(artifact.strata),
                "assignments": len(artifact.assignments),
                "unassigned_target_winners": artifact.target_winners_unassigned,
                "reason_codes": artifact.readiness.reason_codes,
            },
            ensure_ascii=False,
        )
    )
    if args.require_calibrated and artifact.status.value != "calibrated":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
