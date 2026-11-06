#!/usr/bin/env python3
"""Isolated msFIDDLE batch worker.

This process is intentionally launched with the dedicated Python 3.10
environment.  The V10 web/runtime environment never imports the model's
PyTorch dependency directly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one msFIDDLE batch")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--rescore-model", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--instrument-type", default="qtof")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--top-k", type=int, default=5)
    return parser.parse_args()


def _load_spectra(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    spectra = payload.get("spectra") if isinstance(payload, dict) else None
    if not isinstance(spectra, list):
        raise ValueError("Input must contain a 'spectra' JSON array")
    return spectra


def main() -> int:
    args = parse_args()
    from msfiddle import MsFiddlePredictor

    predictor = MsFiddlePredictor(
        instrument_type=args.instrument_type,
        device=args.device,
        no_cuda=args.device == "cpu",
        config_path=str(args.config) if args.config else None,
        resume_path=str(args.model),
        rescore_resume_path=str(args.rescore_model),
        download_models=False,
        batch_size=max(1, int(args.batch_size)),
        verbose=False,
    )
    predictions = predictor.predict_batch(
        _load_spectra(args.input),
        top_k=max(1, int(args.top_k)),
        batch_size=max(1, int(args.batch_size)),
    )
    args.output.write_text(
        json.dumps({"predictions": predictions}, ensure_ascii=False),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
