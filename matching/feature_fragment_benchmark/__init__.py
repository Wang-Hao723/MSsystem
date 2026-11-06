"""Public-spectrum benchmark for auditable feature-fragment evidence.

This package is deliberately separate from the production rankers.  It emits
shadow evidence and validation receipts; it cannot promote an MS1 candidate,
an MS/MS library candidate, an isomer group, or a standard confirmation.
"""

from .corpus import (
    CorpusBuildError,
    build_massbank_corpus,
    iter_json_array,
    iter_jsonl,
)
from .export import export_empirically_validated_rules
from .split import build_isolated_split, validate_split_isolation
from .statistics import calibrate_rules, discover_rules, test_frozen_rules

__all__ = [
    "CorpusBuildError",
    "build_isolated_split",
    "build_massbank_corpus",
    "calibrate_rules",
    "discover_rules",
    "export_empirically_validated_rules",
    "iter_json_array",
    "iter_jsonl",
    "test_frozen_rules",
    "validate_split_isolation",
]
