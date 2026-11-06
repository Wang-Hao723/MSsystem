"""Load the V15 feature-fragment catalogs from explicit paths."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .schema import validate_catalog


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_DIR = ROOT / "config" / "feature_fragments"


def _read(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_catalog(config_dir: Path | str | None = None) -> dict[str, Any]:
    """Load and validate sources, instrument domains, and rules.

    Only these three named files are read.  No data directory discovery or
    recursive scan occurs.
    """
    directory = Path(config_dir) if config_dir is not None else DEFAULT_CONFIG_DIR
    sources = _read(directory / "source_catalog.json")
    domains = _read(directory / "instrument_domains.json")
    rules = _read(directory / "rules_v1.json")
    catalog = {
        "schema_version": rules["schema_version"],
        "sources": sources["sources"],
        "instrument_domains": domains["instrument_domains"],
        "rules": rules["rules"],
    }
    validate_catalog(catalog)
    return catalog
