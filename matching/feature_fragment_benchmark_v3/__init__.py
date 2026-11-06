"""Fail-closed V3 benchmark for instrument-conditioned fragment evidence.

This package is research-only.  It does not import or mutate the runtime
feature-fragment catalog and cannot contribute to search ranking.
"""

from .parsing import parse_collision_energy_program, parse_scan_range
from .provenance import (
    MetadataValidationError,
    make_acquisition_replication_group,
    normalize_public_spectrum,
)
from .statistics import estimate_conditional_lr, fit_ce_response

__all__ = [
    "MetadataValidationError",
    "estimate_conditional_lr",
    "fit_ce_response",
    "make_acquisition_replication_group",
    "normalize_public_spectrum",
    "parse_collision_energy_program",
    "parse_scan_range",
]
