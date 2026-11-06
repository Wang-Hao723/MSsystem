#!/usr/bin/env python3
"""Freeze, calibrate, and evaluate the EviDiscrim-MS blind benchmark."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from eval.evidiscrim_benchmark import (  # noqa: E402
    evaluate_locked_blind_predictions,
    freeze_identity_disjoint_blind_set,
    load_authentic_standard_manifest,
    load_benchmark_policy,
    lock_prediction_bundle,
)
from matching.qtof_calibration import (  # noqa: E402
    apply_qtof_intensity_calibrator,
    fit_qtof_intensity_calibrator,
)


DEFAULT_POLICY = ROOT / "config" / "evidiscrim_benchmark_policy.json"


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_rows(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".jsonl":
        return [
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    payload = _read_json(path)
    if isinstance(payload, dict):
        if "records" in payload:
            payload = payload["records"]
        elif "observations" in payload:
            payload = payload["observations"]
    if not isinstance(payload, list):
        raise ValueError(f"{path} must contain a JSON list")
    return [dict(item) for item in payload if isinstance(item, dict)]


def _atomic_json(path: Path, payload: Any, *, mode: int = 0o644) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _require_writable_targets(paths: list[Path], *, force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        joined = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to overwrite without --force: {joined}")


def freeze(args: argparse.Namespace) -> dict[str, Any]:
    policy = load_benchmark_policy(args.policy)
    records = load_authentic_standard_manifest(args.manifest)
    bundle = freeze_identity_disjoint_blind_set(records, policy=policy)
    targets = {
        "train": args.output_dir / "train_authentic_standards.json",
        "calibration": args.output_dir / "calibration_authentic_standards.json",
        "blind": args.output_dir / "blind_test_public.json",
        "labels": args.output_dir / "blind_test_labels.sealed.json",
        "audit": args.output_dir / "identity_disjoint_split_audit.json",
        "receipt": args.output_dir / "freeze_receipt.json",
    }
    _require_writable_targets(list(targets.values()), force=args.force)
    _atomic_json(targets["train"], {"records": bundle["train_records"]})
    _atomic_json(
        targets["calibration"], {"records": bundle["calibration_records"]}
    )
    _atomic_json(targets["blind"], bundle["blind_test_public"])
    _atomic_json(targets["labels"], bundle["sealed_labels"], mode=0o600)
    _atomic_json(targets["audit"], bundle["leakage_audit"])
    receipt = {
        key: value for key, value in bundle.items()
        if key not in {
            "train_records", "calibration_records", "blind_test_public",
            "sealed_labels",
        }
    }
    receipt["artifacts"] = {key: str(path) for key, path in targets.items()}
    _atomic_json(targets["receipt"], receipt)
    return receipt


def calibrate_qtof(args: argparse.Namespace) -> dict[str, Any]:
    policy = load_benchmark_policy(args.policy)
    examples = _read_rows(args.examples)
    blind_keys: list[str] = []
    if args.sealed_labels:
        labels = _read_json(args.sealed_labels)
        blind_keys = sorted({
            str(connectivity_key)
            for item in labels.get("labels") or ()
            for connectivity_key in (
                item.get("candidate_connectivity_by_id") or {}
            ).values()
            if str(connectivity_key)
        })
    artifact = fit_qtof_intensity_calibrator(
        examples,
        policy=policy,
        blind_connectivity_keys=blind_keys,
    )
    _require_writable_targets([args.output], force=args.force)
    _atomic_json(args.output, artifact)
    return artifact


def apply_qtof(args: argparse.Namespace) -> dict[str, Any]:
    calibrator = _read_json(args.calibrator)
    predictions = _read_rows(args.predictions)
    calibrated = [
        apply_qtof_intensity_calibrator(
            prediction,
            calibrator=calibrator,
            instrument_context=prediction,
        )
        for prediction in predictions
    ]
    artifact = {
        "schema_version": "1.0.0",
        "status": "calibrated",
        "calibrator_id": calibrator.get("calibrator_id"),
        "counts_as_independent_architecture": False,
        "base_architecture_family": calibrator.get(
            "base_architecture_family"
        ),
        "records": calibrated,
    }
    _require_writable_targets([args.output], force=args.force)
    _atomic_json(args.output, artifact)
    return artifact


def lock_predictions(args: argparse.Namespace) -> dict[str, Any]:
    public_blind = _read_json(args.public_blind)
    predictions = _read_rows(args.predictions)
    bundle = lock_prediction_bundle(public_blind, predictions)
    _require_writable_targets([args.output], force=args.force)
    _atomic_json(args.output, bundle)
    return bundle


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    policy = load_benchmark_policy(args.policy)
    report = evaluate_locked_blind_predictions(
        public_blind_artifact=_read_json(args.public_blind),
        sealed_labels=_read_json(args.sealed_labels),
        predictions=_read_json(args.predictions),
        project_fdr_observations=_read_rows(args.project_fdr),
        split_leakage_audit=_read_json(args.split_audit),
        policy=policy,
    )
    _require_writable_targets([args.output], force=args.force)
    _atomic_json(args.output, report)
    return report


def parser() -> argparse.ArgumentParser:
    main = argparse.ArgumentParser(description=__doc__)
    main.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    commands = main.add_subparsers(dest="command", required=True)

    freeze_parser = commands.add_parser("freeze")
    freeze_parser.add_argument("--manifest", type=Path, required=True)
    freeze_parser.add_argument("--output-dir", type=Path, required=True)
    freeze_parser.add_argument("--force", action="store_true")
    freeze_parser.set_defaults(handler=freeze)

    calibration_parser = commands.add_parser("calibrate-qtof")
    calibration_parser.add_argument("--examples", type=Path, required=True)
    calibration_parser.add_argument("--sealed-labels", type=Path)
    calibration_parser.add_argument("--output", type=Path, required=True)
    calibration_parser.add_argument("--force", action="store_true")
    calibration_parser.set_defaults(handler=calibrate_qtof)

    apply_calibration_parser = commands.add_parser("apply-qtof")
    apply_calibration_parser.add_argument(
        "--predictions", type=Path, required=True
    )
    apply_calibration_parser.add_argument(
        "--calibrator", type=Path, required=True
    )
    apply_calibration_parser.add_argument("--output", type=Path, required=True)
    apply_calibration_parser.add_argument("--force", action="store_true")
    apply_calibration_parser.set_defaults(handler=apply_qtof)

    lock_parser = commands.add_parser("lock-predictions")
    lock_parser.add_argument("--public-blind", type=Path, required=True)
    lock_parser.add_argument("--predictions", type=Path, required=True)
    lock_parser.add_argument("--output", type=Path, required=True)
    lock_parser.add_argument("--force", action="store_true")
    lock_parser.set_defaults(handler=lock_predictions)

    evaluate_parser = commands.add_parser("evaluate")
    evaluate_parser.add_argument("--public-blind", type=Path, required=True)
    evaluate_parser.add_argument("--sealed-labels", type=Path, required=True)
    evaluate_parser.add_argument("--predictions", type=Path, required=True)
    evaluate_parser.add_argument("--project-fdr", type=Path, required=True)
    evaluate_parser.add_argument("--split-audit", type=Path, required=True)
    evaluate_parser.add_argument("--output", type=Path, required=True)
    evaluate_parser.add_argument("--force", action="store_true")
    evaluate_parser.set_defaults(handler=evaluate)
    return main


def main() -> int:
    args = parser().parse_args()
    result = args.handler(args)
    print(json.dumps({
        "status": result.get("status"),
        "policy_id": result.get("policy_id"),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
