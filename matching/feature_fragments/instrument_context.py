"""Fail-closed acquisition-domain resolution for feature-fragment rules.

The resolver deliberately requires an exact instrument model plus the
reported activation, collision-energy unit, and acquisition mode.  Analyzer
families such as ``QTOF`` or ``Orbitrap`` are not sufficient because collision
cells and vendor energy scales are not interchangeable.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any


_MODE_ALIASES = {
    "single": "fixed",
    "fixed": "fixed",
    "stepped": "stepped",
    "step": "stepped",
    "ramp": "ramp",
    "ramped": "ramp",
    "rolling": "rolling",
    "dynamic": "dynamic",
    "ces": "CES",
    "bbcid": "bbCID",
    "mse": "MSE",
}


def _compact(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def normalize_activation_method(value: Any) -> str:
    """Normalize only activation labels with unambiguous semantics."""
    text = str(value or "").strip().lower()
    compact = _compact(text)
    if not compact:
        return ""
    if (
        "higherenergycollisionaldissociation" in compact
        or "higherenergycollisiondissociation" in compact
        or compact == "hcd"
        or compact.endswith("hcd")
    ):
        return "HCD"
    if "beamtypecollisioninduceddissociation" in compact or compact in {
        "beamcid",
        "beamcollisioninduceddissociation",
    }:
        return "beam_cid"
    if "collisioninduceddissociation" in compact or compact == "cid":
        # Plain CID is retained as unresolved.  It can be mapped to beam CID
        # only after an exact QTOF model selects a unique catalog domain.
        return "CID"
    return ""


def normalize_collision_energy_mode(value: Any) -> str:
    return _MODE_ALIASES.get(str(value or "").strip().lower(), "")


def _model_matches(instrument_text: str, model: str) -> bool:
    instrument = _compact(instrument_text)
    expected = _compact(model)
    return bool(instrument and expected and expected in instrument)


def resolve_instrument_domain(
    *,
    catalog: Mapping[str, Any],
    instrument: Any = None,
    instrument_type: Any = None,
    activation: Any = None,
    collision_energy_unit: Any = None,
    collision_energy_mode: Any = None,
    collision_cell_region: Any = None,
    explicit_domain_id: Any = None,
) -> dict[str, Any]:
    """Resolve one exact catalog domain without converting CE quantities."""
    reported_activation = normalize_activation_method(activation)
    reported_mode = normalize_collision_energy_mode(collision_energy_mode)
    unit = str(collision_energy_unit or "").strip()
    region = str(collision_cell_region or "").strip()
    explicit = str(explicit_domain_id or "").strip()
    instrument_text = " ".join(
        part for part in (str(instrument or ""), str(instrument_type or "")) if part
    )

    base = {
        "status": "unavailable",
        "instrument_domain_id": "",
        "activation": reported_activation,
        "collision_energy_unit": unit,
        "collision_energy_mode": reported_mode,
        "cross_unit_conversion_applied": False,
        "reported_instrument": str(instrument or ""),
        "reported_instrument_type": str(instrument_type or ""),
        "reported_activation": str(activation or ""),
        "reported_collision_cell_region": region,
    }
    if not unit:
        return {**base, "reason": "collision_energy_unit_missing"}
    if not reported_mode:
        return {**base, "reason": "collision_energy_mode_missing_or_unsupported"}
    if not reported_activation:
        return {**base, "reason": "activation_method_missing_or_unsupported"}

    domains = list(catalog.get("instrument_domains") or [])
    if explicit:
        domains = [domain for domain in domains if domain.get("domain_id") == explicit]
        if not domains:
            return {**base, "reason": "explicit_instrument_domain_unknown"}
    elif not instrument_text.strip():
        return {**base, "reason": "exact_instrument_model_missing"}

    candidates: list[tuple[dict[str, Any], str]] = []
    for domain in domains:
        matched_model = ""
        if explicit:
            matched_model = "explicit_domain_id"
        else:
            matched_model = next(
                (
                    str(model)
                    for model in domain.get("models") or []
                    if _model_matches(instrument_text, str(model))
                ),
                "",
            )
            if not matched_model:
                continue
        domain_activation = str(domain.get("activation") or "")
        activation_matches = reported_activation == domain_activation
        if (
            reported_activation == "CID"
            and domain_activation == "beam_cid"
            and str(domain.get("analyzer") or "").upper() == "QTOF"
        ):
            activation_matches = True
        if not activation_matches:
            continue
        if str(domain.get("ce_unit") or "") != unit:
            continue
        if str(domain.get("ce_mode") or "") != reported_mode:
            continue
        domain_region = str(domain.get("collision_cell_region") or "")
        if region and _compact(region) != _compact(domain_region):
            continue
        candidates.append((domain, matched_model))

    if not candidates:
        return {
            **base,
            "reason": "no_exact_model_activation_ce_domain_match",
            "candidate_domain_ids": [],
        }
    if len(candidates) > 1:
        return {
            **base,
            "reason": "instrument_domain_ambiguous",
            "candidate_domain_ids": sorted(
                str(domain.get("domain_id") or "") for domain, _ in candidates
            ),
        }

    domain, matched_model = candidates[0]
    resolved_activation = str(domain.get("activation") or "")
    return {
        **base,
        "status": "resolved",
        "reason": "exact_model_activation_unit_mode_domain_match",
        "instrument_domain_id": str(domain.get("domain_id") or ""),
        "activation": resolved_activation,
        "matched_model": matched_model,
        "collision_cell_region": str(domain.get("collision_cell_region") or ""),
        "catalog_vendor": str(domain.get("vendor") or ""),
        "catalog_family": str(domain.get("family") or ""),
        "catalog_analyzer": str(domain.get("analyzer") or ""),
        "candidate_domain_ids": [str(domain.get("domain_id") or "")],
    }


__all__ = [
    "normalize_activation_method",
    "normalize_collision_energy_mode",
    "resolve_instrument_domain",
]
