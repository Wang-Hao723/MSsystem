"""Pure helpers for describing same-formula candidate relationships.

The first InChIKey block represents molecular connectivity.  A different full
key with the same first block is therefore not the same problem as a
same-formula candidate with different connectivity.  Keeping those cases
separate prevents ordinary achiral LC-MS/MS evidence from overclaiming
stereochemical resolution.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Iterable

from matching.library_quality import (
    compound_identity_key,
    normalize_inchikey,
    resolve_structure_formula,
)

_UNSPECIFIED_STEREO_ISOTOPE_BLOCK = "UHFFFAOYSA"


def connectivity_block(item: dict[str, Any] | None) -> str:
    if not isinstance(item, dict):
        return ""
    key = normalize_inchikey(item.get("inchikey"))
    return key.split("-", 1)[0] if key else ""


def inchikey_layers(
    item: dict[str, Any] | None,
) -> tuple[str, str, str]:
    """Return connectivity, stereo/isotope and protonation key layers.

    The second InChIKey block is deliberately described conservatively.  It can
    encode more than stereochemistry, so a differing block is not automatically
    labelled as an enantiomer/diastereomer.
    """
    if not isinstance(item, dict):
        return "", "", ""
    key = normalize_inchikey(item.get("inchikey"))
    parts = key.split("-") if key else []
    if len(parts) != 3:
        return "", "", ""
    return parts[0], parts[1], parts[2]


def same_connectivity_variant_subtype(
    left: dict[str, Any] | None,
    right: dict[str, Any] | None,
) -> str:
    """Describe why two full keys differ despite shared connectivity."""
    left_connectivity, left_definition, left_protonation = inchikey_layers(left)
    right_connectivity, right_definition, right_protonation = inchikey_layers(
        right
    )
    if (
        not left_connectivity
        or left_connectivity != right_connectivity
        or compound_identity_key(left) == compound_identity_key(right)
    ):
        return ""
    if (
        left_definition
        and left_definition == right_definition
        and left_protonation
        and left_protonation != right_protonation
    ):
        return "protonation_or_charge_layer_variant"
    if left_definition and right_definition and left_definition != right_definition:
        if _UNSPECIFIED_STEREO_ISOTOPE_BLOCK in {
            left_definition,
            right_definition,
        }:
            return "structure_definition_specificity_variant"
        return "stereochemical_isotopic_or_tautomeric_variant"
    return "same_connectivity_variant_unresolved"


def structure_definition_is_unspecified(
    item: dict[str, Any] | None,
) -> bool:
    """Whether the full key carries the common unspecified stereo/isotope block."""
    _, definition, _ = inchikey_layers(item)
    return definition == _UNSPECIFIED_STEREO_ISOTOPE_BLOCK


def isomer_relationship(
    left: dict[str, Any] | None,
    right: dict[str, Any] | None,
) -> str:
    """Classify a candidate pair without inventing stereochemical certainty."""
    if not isinstance(left, dict) or not isinstance(right, dict):
        return "unknown_relationship"
    if compound_identity_key(left) == compound_identity_key(right):
        return "same_identity"
    left_formula = resolve_structure_formula(
        str(left.get("smiles") or ""),
        str(left.get("formula") or ""),
        str(left.get("adduct") or ""),
    )
    right_formula = resolve_structure_formula(
        str(right.get("smiles") or ""),
        str(right.get("formula") or ""),
        str(right.get("adduct") or ""),
    )
    if not left_formula or left_formula != right_formula:
        return "not_same_formula"
    left_block = connectivity_block(left)
    right_block = connectivity_block(right)
    if left_block and right_block and left_block == right_block:
        return "same_connectivity_variant"
    if left_block and right_block and left_block != right_block:
        return "constitutional_or_positional_isomer"
    return "same_formula_identity_unresolved"


def summarize_isomer_relationships(
    selected: dict[str, Any],
    alternatives: Iterable[dict[str, Any]],
) -> dict[str, Any]:
    alternatives = [
        candidate for candidate in alternatives if isinstance(candidate, dict)
    ]
    relationships = [
        isomer_relationship(selected, candidate)
        for candidate in alternatives
    ]
    relationships = [
        value for value in relationships
        if value not in {"same_identity", "not_same_formula"}
    ]
    counts = Counter(relationships)
    variant_subtypes = Counter(
        subtype
        for candidate in alternatives
        for subtype in [
            same_connectivity_variant_subtype(selected, candidate)
        ]
        if subtype
    )
    if not relationships:
        group_type = "none"
    elif set(relationships) == {"same_connectivity_variant"}:
        group_type = "same_connectivity_variant_group"
    elif set(relationships) == {"constitutional_or_positional_isomer"}:
        group_type = "constitutional_or_positional_isomer_group"
    elif len(set(relationships)) > 1:
        group_type = "mixed_isomer_group"
    else:
        group_type = "same_formula_identity_unresolved_group"
    return {
        "ambiguity_type": group_type,
        "relationship_counts": dict(sorted(counts.items())),
        "same_connectivity_variant_subtype_counts": dict(
            sorted(variant_subtypes.items())
        ),
        "msms_exact_resolution_capable": (
            bool(relationships)
            and "same_connectivity_variant" not in relationships
            and "same_formula_identity_unresolved" not in relationships
        ),
    }
