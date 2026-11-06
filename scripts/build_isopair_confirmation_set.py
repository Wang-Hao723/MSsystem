#!/usr/bin/env python3
"""Build a fresh same-formula confirmation set excluding prior IsoPair families."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = Path(__file__).resolve().parent
for path in (str(PROJECT_ROOT), str(SCRIPTS_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from train_isopair_ms import (
    DEFAULT_DB,
    _connect_readonly,
    fetch_spectra,
    filter_usable_families,
    select_family_strata,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument(
        "--exclude-cache",
        type=Path,
        default=PROJECT_ROOT / "models" / "isopair_ms_v1_full",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PROJECT_ROOT / "models" / "isopair_confirmation_20260802",
    )
    parser.add_argument("--families", type=int, default=500)
    parser.add_argument("--seed", type=int, default=99173)
    parser.add_argument("--spectra-per-identity", type=int, default=4)
    args = parser.parse_args()

    old_rows = json.loads(
        (args.exclude_cache / "pilot_spectra.json").read_text(encoding="utf-8")
    )
    excluded_formulas = {str(row["formula"]) for row in old_rows}
    connection = _connect_readonly(args.db)
    # Oversample every deterministic fold, then collapse all remaining families
    # into a confirmation-only test set.  No parameter or threshold is tuned on it.
    oversampled = select_family_strata(
        connection,
        seed=args.seed,
        train_families=args.families,
        validation_families=args.families,
        test_families=args.families,
        min_identities=2,
        max_identities=12,
    )
    candidates = sorted({
        (formula, adduct)
        for strata in oversampled.values()
        for formula, adduct in strata
        if formula not in excluded_formulas
    })
    candidates.sort(
        key=lambda item: (
            __import__("hashlib").sha256(
                f"{args.seed}|confirmation|{item[0]}|{item[1]}".encode()
            ).hexdigest(),
            item,
        )
    )
    chosen = candidates[:args.families]
    if len(chosen) < args.families:
        raise RuntimeError(
            f"Only {len(chosen)} non-overlapping families were available"
        )
    selected = {"train": [], "validation": [], "test": chosen}
    spectra, splits = fetch_spectra(
        connection,
        selected,
        max_spectra_per_identity=args.spectra_per_identity,
    )
    connection.close()
    spectra, splits = filter_usable_families(spectra, splits)
    retained_formulas = {row.formula for row in spectra}
    overlap = retained_formulas & excluded_formulas
    if overlap:
        raise RuntimeError(f"Confirmation family leakage detected: {sorted(overlap)[:3]}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "pilot_spectra.json").write_text(
        json.dumps([asdict(row) for row in spectra], ensure_ascii=False),
        encoding="utf-8",
    )
    np.savez_compressed(
        args.out_dir / "pilot_embeddings.npz",
        embeddings=np.empty((len(spectra), 0), dtype=np.float16),
        splits=np.asarray(splits),
    )
    manifest = {
        "schema_version": 1,
        "purpose": "locked external-family confirmation only",
        "seed": args.seed,
        "requested_families": args.families,
        "retained_families": len(retained_formulas),
        "spectra": len(spectra),
        "excluded_prior_families": len(excluded_formulas),
        "overlap_with_prior_families": 0,
        "split": "test only",
    }
    (args.out_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
