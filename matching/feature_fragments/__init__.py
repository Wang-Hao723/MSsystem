"""Auditable feature-fragment research and fail-closed runtime controls.

The research evaluator remains importable for deterministic tests and future
benchmark development.  Search first enforces the sealed runtime export gate;
the current gate exports zero rules, so the research catalog is not evaluated
for production results.  Feature-fragment output never confirms an exact
identity, changes a candidate rank, or releases an authentic-standard claim.
"""

from .catalog import load_catalog
from .evidence import (
    evaluate_feature_fragment_evidence,
    feature_fragment_score_standard,
)
from .instrument_context import (
    normalize_activation_method,
    normalize_collision_energy_mode,
    resolve_instrument_domain,
)
from .runtime_gate import (
    FeatureFragmentRuntimeGateError,
    load_runtime_gate,
    runtime_gate_audit,
    validate_runtime_gate,
)
from .schema import CatalogValidationError, validate_catalog

__all__ = [
    "CatalogValidationError",
    "FeatureFragmentRuntimeGateError",
    "evaluate_feature_fragment_evidence",
    "feature_fragment_score_standard",
    "load_catalog",
    "load_runtime_gate",
    "normalize_activation_method",
    "normalize_collision_energy_mode",
    "resolve_instrument_domain",
    "runtime_gate_audit",
    "validate_catalog",
    "validate_runtime_gate",
]
