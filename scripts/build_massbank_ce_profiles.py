#!/usr/bin/env python3
"""Backward-compatible entry point for the multi-source CE profile builder."""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_ce_profiles import (  # noqa: F401
    build_database,
    iter_mgf_records,
    iter_msp_records,
    iter_source_records,
    main,
    resolve_input_files,
)


if __name__ == "__main__":
    raise SystemExit(main())
