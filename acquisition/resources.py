"""Fail-closed resource declarations for active-evidence planning."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any

from schemas.acquisition import AcquisitionNeed, AcquisitionResourceDeclaration


RESOURCE_GATED_ACTIONS = frozenset({
    "authentic_standard_rt",
    "ion_mobility_ccs",
})

_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def normalize_resource_identity_key(value: str) -> str:
    """Normalize full InChIKeys to the ambiguity graph's canonical key."""

    cleaned = str(value or "").strip()
    if not cleaned:
        raise ValueError("resource identity key cannot be blank")
    if cleaned.casefold().startswith("ik::"):
        suffix = cleaned[4:].strip().upper()
        if not _INCHIKEY.fullmatch(suffix):
            raise ValueError(f"invalid InChIKey resource identity: {cleaned}")
        return f"ik::{suffix}"
    upper = cleaned.upper()
    if _INCHIKEY.fullmatch(upper):
        return f"ik::{upper}"
    if cleaned.startswith(("sm::", "nf::")):
        return cleaned
    raise ValueError(
        "resource identities must be full InChIKeys or canonical ik::/sm::/nf:: keys"
    )


def _digest(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def resolve_acquisition_resources(
    context: Mapping[str, Any] | None,
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    """Resolve one task's explicit resources; absence always means unavailable."""

    raw = (context or {}).get("acquisition_resources") or {}
    declaration = AcquisitionResourceDeclaration.model_validate(raw)
    resource_policy = policy.get("resource_policy") or {}
    standard_policy = resource_policy.get("authentic_standard_rt") or {}
    mobility_policy = resource_policy.get("ion_mobility_ccs") or {}

    standard_keys = tuple(sorted({
        normalize_resource_identity_key(value)
        for value in declaration.authentic_standard_identity_keys
    }))
    ccs_keys = tuple(sorted({
        normalize_resource_identity_key(value)
        for value in declaration.ccs_reference_identity_keys
    }))
    standard_cap = max(0, int(standard_policy.get("maximum_action_budget", 200)))
    mobility_cap = max(0, int(mobility_policy.get("maximum_action_budget", 200)))
    standard_budget = min(
        declaration.authentic_standard_action_budget, standard_cap
    )
    mobility_budget = min(declaration.ion_mobility_action_budget, mobility_cap)
    if not standard_keys:
        standard_budget = 0
    if not declaration.ion_mobility_available or not ccs_keys:
        mobility_budget = 0

    canonical = {
        "declaration_id": declaration.declaration_id,
        "authentic_standard_identity_keys": standard_keys,
        "authentic_standard_action_budget": standard_budget,
        "ion_mobility_available": declaration.ion_mobility_available,
        "ccs_reference_identity_keys": ccs_keys,
        "ion_mobility_action_budget": mobility_budget,
        "default_fail_closed": True,
    }
    resource_sha256 = _digest(canonical)
    explicit = bool(
        standard_budget
        or mobility_budget
        or standard_keys
        or ccs_keys
        or declaration.ion_mobility_available
    )
    return {
        **canonical,
        "mode": "explicit_declaration" if explicit else "sample_only_fail_closed",
        "resource_context_sha256": resource_sha256,
        "effective_declaration_id": (
            declaration.declaration_id or f"resource-{resource_sha256[:12]}"
        ),
    }


def reference_keys_for_need(
    need: AcquisitionNeed,
    resources: Mapping[str, Any],
    action_kind: str,
) -> tuple[str, ...]:
    if action_kind == "authentic_standard_rt":
        available = set(resources.get("authentic_standard_identity_keys") or ())
    elif action_kind == "ion_mobility_ccs":
        available = set(resources.get("ccs_reference_identity_keys") or ())
    else:
        return ()
    return tuple(sorted(
        candidate.identity_key
        for candidate in need.candidates
        if candidate.identity_key in available
    ))


def allocate_resource_actions(
    needs: Sequence[AcquisitionNeed],
    resources: Mapping[str, Any],
) -> dict[str, frozenset[str]]:
    """Allocate global per-task budgets without inventing unavailable actions."""

    allocations: dict[str, set[str]] = {need.need_id: set() for need in needs}
    specifications = (
        ("authentic_standard_rt", "authentic_standard_action_budget"),
        ("ion_mobility_ccs", "ion_mobility_action_budget"),
    )
    for action_kind, budget_key in specifications:
        remaining = max(0, int(resources.get(budget_key) or 0))
        if not remaining:
            continue
        eligible: list[tuple[int, int, AcquisitionNeed]] = []
        for order, need in enumerate(needs):
            keys = set(reference_keys_for_need(need, resources, action_kind))
            if not keys:
                continue
            covered_edges = sum(
                edge.left.identity_key in keys or edge.right.identity_key in keys
                for edge in need.ambiguity_edges
            )
            if covered_edges:
                eligible.append((-covered_edges, order, need))
        for _, _, need in sorted(eligible)[:remaining]:
            allocations[need.need_id].add(action_kind)
    return {
        need_id: frozenset(action_kinds)
        for need_id, action_kinds in allocations.items()
    }


def public_resource_summary(
    resources: Mapping[str, Any],
    *,
    standard_actions_planned: int,
    mobility_actions_planned: int,
) -> dict[str, Any]:
    """Return audit metadata without dumping a laboratory inventory."""

    omitted = []
    if not standard_actions_planned:
        omitted.append("authentic_standard_rt")
    if not mobility_actions_planned:
        omitted.append("ion_mobility_ccs")
    return {
        "mode": resources.get("mode"),
        "default_fail_closed": True,
        "declaration_id": resources.get("declaration_id"),
        "resource_context_sha256": resources.get("resource_context_sha256"),
        "authentic_standard_inventory_count": len(
            resources.get("authentic_standard_identity_keys") or ()
        ),
        "authentic_standard_action_budget": int(
            resources.get("authentic_standard_action_budget") or 0
        ),
        "authentic_standard_actions_planned": standard_actions_planned,
        "ion_mobility_available": bool(resources.get("ion_mobility_available")),
        "ccs_reference_count": len(resources.get("ccs_reference_identity_keys") or ()),
        "ion_mobility_action_budget": int(
            resources.get("ion_mobility_action_budget") or 0
        ),
        "ion_mobility_actions_planned": mobility_actions_planned,
        "omitted_action_kinds": omitted,
        "sample_only_primary_workflow": not (
            standard_actions_planned or mobility_actions_planned
        ),
    }
