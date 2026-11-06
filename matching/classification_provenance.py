"""Lightweight trust rules for candidate classification provenance.

Production ranking and confidence code must never treat a class inferred from
the same spectrum, or from an unvalidated structure model, as independent
chemical evidence.  Keep the policy in this dependency-light module so search,
confidence annotation, mining scripts, and tests all enforce the same rule.
"""

from __future__ import annotations

from typing import Any, Mapping


TRUSTED_CLASSIFICATION_SOURCE_PREFIXES = (
    "curated",
    "manual",
    "smarts",
    "structure",
    "rule",
    "database_curated",
)

UNTRUSTED_CLASSIFICATION_SOURCES = {
    "",
    "unknown",
    "other",
    "ml",
    "msclass",
    "msclass_spectrum",
    "fragment",
    "fragment_inference",
    "spectrum",
    "existing",
}

_NON_CLASSES = {"", "unknown", "other", "unclassified", "none", "nan"}


def classification_source_trust(
    source: Any,
    *,
    allow_legacy_missing: bool = False,
) -> tuple[bool, str]:
    """Return whether a stored class source may unlock class-specific rules."""
    if source is None:
        if allow_legacy_missing:
            return True, "legacy_caller_without_source_metadata"
        return False, "missing_source_metadata"
    normalized = str(source or "").strip().lower()
    if normalized in UNTRUSTED_CLASSIFICATION_SOURCES:
        return False, f"untrusted_source:{normalized or 'missing'}"
    if normalized.startswith(("ml_", "model_", "msclass_", "fragment_")):
        return False, f"untrusted_source:{normalized}"
    if normalized.startswith(TRUSTED_CLASSIFICATION_SOURCE_PREFIXES):
        return True, f"trusted_source:{normalized}"
    return False, f"unreviewed_source:{normalized}"


def trusted_analysis_class(
    candidate: Mapping[str, Any],
    *,
    allow_legacy_missing: bool = False,
) -> tuple[str, bool, str]:
    """Resolve the only class allowed to influence scoring or confidence.

    A class independently confirmed by the structure-aware expert layer wins
    even when the stored source was ML.  Otherwise the stored class is accepted
    only when its provenance is explicitly trusted.
    """
    details = candidate.get("expert_details")
    if isinstance(details, Mapping):
        provenance = details.get("classification_provenance")
        analysis_class = str(details.get("analysis_class") or "").strip()
        if isinstance(provenance, Mapping):
            analysis_class = str(
                provenance.get("analysis_class") or analysis_class
            ).strip()
            if (
                bool(provenance.get("reliable"))
                and analysis_class.lower() not in _NON_CLASSES
            ):
                return (
                    analysis_class,
                    True,
                    str(
                        provenance.get("reliability_reason")
                        or "independently_confirmed_analysis_class"
                    ),
                )

    raw_class = str(candidate.get("classification") or "").strip()
    trusted, reason = classification_source_trust(
        candidate.get("classification_source"),
        allow_legacy_missing=allow_legacy_missing,
    )
    if trusted and raw_class.lower() not in _NON_CLASSES:
        return raw_class, True, reason

    # Deterministic SMILES rules are independently reproducible from the
    # candidate structure and remain valid for historical caches whose expert
    # payload predates a newly added rule.  Import lazily so this policy module
    # stays cheap for callers that have no structures or no RDKit runtime.
    smiles = str(candidate.get("smiles") or "").strip()
    if smiles:
        try:
            from mssystem.classification import RuleBasedClassifier

            structure_class = str(
                RuleBasedClassifier().predict(smiles=smiles) or ""
            ).strip()
        except Exception:
            structure_class = ""
        if structure_class.lower() not in _NON_CLASSES:
            return (
                structure_class,
                True,
                "smiles_structure_rule"
                if raw_class.lower() in _NON_CLASSES
                else "smiles_structure_rule_overrode_untrusted_class",
            )
    return "Unknown", False, reason


__all__ = [
    "TRUSTED_CLASSIFICATION_SOURCE_PREFIXES",
    "UNTRUSTED_CLASSIFICATION_SOURCES",
    "classification_source_trust",
    "trusted_analysis_class",
]
