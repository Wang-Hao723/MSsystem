#!/usr/bin/env python3
"""Run the EviDiscrim-MS flagship module for one committed V10 result."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.evidiscrim import (
    evidiscrim_artifact_path,
    load_evidiscrim_policy,
    stream_select_evidiscrim_rows,
)
from matching.openworld_shadow import (
    FIORA_ALLOWED_ADDUCTS,
    FioraRuntime,
    build_fiora_evidiscrim_shadow,
    write_shadow_artifact,
)
from schemas.evidiscrim import EviDiscrimArtifact


DEFAULT_FIORA_EXECUTABLE = Path("/opt/anaconda3/envs/fiora/bin/fiora-predict")
DEFAULT_FIORA_MODEL = Path(
    "/Users/zizi/Desktop/fiora/fiora/resources/models/fiora_OS_v1.0.0.pt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate candidate-contrastive diagnostic ions, prediction "
            "uncertainty, and planning-only collision-energy separation as a "
            "shadow sidecar. Production ranking is never modified."
        )
    )
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--policy",
        type=Path,
        default=ROOT / "config" / "evidiscrim_policy.json",
    )
    parser.add_argument(
        "--collision-energies",
        type=float,
        nargs="*",
        help=(
            "Instrument-native CE grid. Omit to use the versioned policy; an "
            "unqualified CE unit never receives an expanded numeric grid."
        ),
    )
    parser.add_argument(
        "--fiora-executable", type=Path, default=DEFAULT_FIORA_EXECUTABLE
    )
    parser.add_argument("--fiora-model", type=Path, default=DEFAULT_FIORA_MODEL)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--max-candidates", type=int, default=5)
    parser.add_argument("--max-spectra", type=int, default=25)
    parser.add_argument(
        "--stream-threshold-mb",
        type=float,
        default=128.0,
        help="Stream large result arrays above this size instead of loading them whole.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_path = args.results or (
        ROOT / "output" / "tasks" / f"{args.task_id}_results.json"
    )
    output_path = args.output or evidiscrim_artifact_path(
        ROOT / "output" / "tasks", args.task_id
    )
    policy = load_evidiscrim_policy(args.policy)
    initial_stat = results_path.stat()
    stat_token = (initial_stat.st_size, initial_stat.st_mtime_ns)
    streaming = initial_stat.st_size >= max(
        0.0, float(args.stream_threshold_mb)
    ) * 1024 * 1024
    source_digest_override = None
    selection_audit = None
    if streaming:
        results, source_digest_override, selection_audit = (
            stream_select_evidiscrim_rows(
                results_path,
                max_spectra=args.max_spectra,
                max_candidates_per_spectrum=args.max_candidates,
                supported_adducts=FIORA_ALLOWED_ADDUCTS,
                policy=policy,
            )
        )
    else:
        results = json.loads(results_path.read_text(encoding="utf-8"))
        if not isinstance(results, list):
            raise SystemExit("Results artifact must contain a JSON array")
    selected_stat = results_path.stat()
    if (selected_stat.st_size, selected_stat.st_mtime_ns) != stat_token:
        raise SystemExit("Source result changed during EviDiscrim input selection")
    payload = build_fiora_evidiscrim_shadow(
        results,
        task_id=args.task_id,
        runtime=FioraRuntime(
            executable=args.fiora_executable,
            model_path=args.fiora_model,
            timeout_seconds=args.timeout,
            device=args.device,
        ),
        policy=policy,
        collision_energy_grid=args.collision_energies,
        max_candidates_per_spectrum=args.max_candidates,
        max_spectra=args.max_spectra,
        source_results_sha256_override=source_digest_override,
        source_results_immutability_verified=streaming,
    )
    final_stat = results_path.stat()
    if (final_stat.st_size, final_stat.st_mtime_ns) != stat_token:
        raise SystemExit("Source result changed during EviDiscrim model execution")
    if selection_audit is not None:
        payload["runtime"]["input_selection"] = selection_audit
        payload["runtime"]["source_results_path"] = str(results_path.resolve())
        payload = EviDiscrimArtifact.model_validate(payload).model_dump(mode="json")
    write_shadow_artifact(payload, output_path)
    print(json.dumps({
        "status": payload["status"],
        "module_id": payload["module_id"],
        "output": str(output_path.resolve()),
        "summary": payload["summary"],
        "production_results_modified": False,
    }, ensure_ascii=False, indent=2))
    return 0 if payload["status"] in {"completed_shadow", "not_applicable"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
