"""Extract planning-only candidate observables without promoting predictions."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from acquisition.ambiguity_graph import candidate_identity_key


EVIDISCRIM_CONTEXT_KEY = "__evidiscrim__"


def _mz_values(value: Any, *, peak_pair: bool = False) -> set[float]:
    """Read mass lists, nested [mass, intensity] peaks, or mass dictionaries.

    A flat list denotes multiple masses; only a nested sequence denotes a
    peak pair.  This preserves diagnostic-ion/loss lists without counting
    spectral intensities as additional masses.
    """
    output: set[float] = set()
    if isinstance(value, Mapping):
        for key in ("mz", "fragment_mz", "product_mz", "mass"):
            try:
                parsed = float(value.get(key))
            except (TypeError, ValueError):
                continue
            if math.isfinite(parsed) and parsed > 0:
                output.add(round(parsed, 2))
                break
    elif isinstance(value, (list, tuple)):
        if peak_pair and value and not isinstance(value[0], (Mapping, list, tuple)):
            output.update(_mz_values(value[0]))
        else:
            for item in value:
                output.update(_mz_values(item, peak_pair=isinstance(item, (list, tuple))))
    else:
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            return output
        if math.isfinite(parsed) and parsed > 0:
            output.add(round(parsed, 2))
    return output


def candidate_observables(candidate: Mapping[str, Any]) -> dict[str, Any]:
    fragments: set[float] = set()
    neutral_losses: set[float] = set()
    for key in (
        "diagnostic_ions", "matched_fragments", "library_peaks", "lib_peaks",
        "peaks", "fragment_mz", "product_ions",
    ):
        fragments.update(_mz_values(candidate.get(key)))
    detail = candidate.get("isomer_fragment_detail") or {}
    if isinstance(detail, Mapping):
        for key in ("candidate_specific_product_ions",):
            fragments.update(_mz_values(detail.get(key)))
        for key in ("candidate_specific_neutral_losses",):
            neutral_losses.update(_mz_values(detail.get(key)))
    for key in ("neutral_losses", "diagnostic_losses"):
        neutral_losses.update(_mz_values(candidate.get(key)))
    return {
        "identity_key": candidate_identity_key(candidate),
        "diagnostic_fragments": sorted(fragments)[:30],
        "neutral_losses": sorted(neutral_losses)[:20],
        "formula": str(candidate.get("formula") or ""),
        "adduct": str(candidate.get("adduct") or ""),
        "collision_energy": (
            candidate.get("collision_energy")
            or candidate.get("lib_collision_energy")
        ),
        "source": "experimental_library_or_candidate_metadata",
        "confirmatory": False,
    }


def build_observable_catalog(
    row: Mapping[str, Any],
    *,
    evidiscrim_spectrum: Mapping[str, Any] | None = None,
) -> dict[str, dict[str, Any]]:
    match = row.get("match") or {}
    candidates = [match, *(
        item for item in match.get("top_candidates") or []
        if isinstance(item, Mapping)
    )]
    catalog: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        if not isinstance(candidate, Mapping):
            continue
        observable = candidate_observables(candidate)
        catalog.setdefault(observable["identity_key"], observable)
    if isinstance(evidiscrim_spectrum, Mapping):
        pairwise: dict[str, dict[str, Any]] = {}
        for raw_pair in evidiscrim_spectrum.get("pairwise_comparisons") or []:
            if not isinstance(raw_pair, Mapping):
                continue
            pair_key = str(raw_pair.get("pair_key") or "")
            if not pair_key or not raw_pair.get("eligible"):
                continue
            pairwise[pair_key] = {
                "pair_key": pair_key,
                "left_identity_key": raw_pair.get("left_identity_key"),
                "right_identity_key": raw_pair.get("right_identity_key"),
                "energy_separations": list(raw_pair.get("energy_separations") or []),
                "observed_diagnostic_ions": [
                    item.get("mz")
                    for item in raw_pair.get("diagnostic_ions") or []
                    if isinstance(item, Mapping) and item.get("mz") is not None
                ],
                "score_semantics": raw_pair.get("score_semantics"),
            }
        energy_selection = evidiscrim_spectrum.get("energy_selection") or {}
        catalog[EVIDISCRIM_CONTEXT_KEY] = {
            "module_id": "evidiscrim_candidate_contrast",
            "module_state": "shadow",
            "runtime_influence": False,
            "candidate_conditioned": True,
            "planning_only": True,
            "energy_selection_status": energy_selection.get("status"),
            "collision_energy_unit": energy_selection.get(
                "collision_energy_unit"
            ),
            "pairs": pairwise,
        }
    return catalog


def evidiscrim_pair_energy_separation(
    pair_key: str,
    proposal: Mapping[str, Any],
    catalog: Mapping[str, Mapping[str, Any]],
) -> tuple[float, tuple[float, ...]] | None:
    """Return a planning-only CE separation when every shadow gate is valid."""

    context = catalog.get(EVIDISCRIM_CONTEXT_KEY) or {}
    if not isinstance(context, Mapping):
        return None
    if (
        context.get("module_state") != "shadow"
        or context.get("runtime_influence") is not False
        or context.get("planning_only") is not True
        or context.get("energy_selection_status") != "recommended_shadow"
    ):
        return None
    pair = (context.get("pairs") or {}).get(pair_key) or {}
    if not isinstance(pair, Mapping):
        return None
    proposed_energies = tuple(
        float(value) for value in proposal.get("collision_energies") or ()
    )
    proposed_unit = str(proposal.get("collision_energy_unit") or "")
    matched: list[Mapping[str, Any]] = []
    for raw in pair.get("energy_separations") or []:
        if not isinstance(raw, Mapping):
            continue
        try:
            energy = float(raw.get("collision_energy"))
            separation = float(raw.get("spectral_distribution_divergence"))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(energy) or not math.isfinite(separation):
            continue
        if str(raw.get("collision_energy_unit") or "") != proposed_unit:
            continue
        if (
            raw.get("ensemble_qualified") is not True
            or raw.get("instrument_domain_qualified") is not True
        ):
            continue
        if any(abs(energy - proposed) <= 0.01 for proposed in proposed_energies):
            if 0 <= separation <= 1:
                matched.append(raw)
    if not matched:
        return None
    best = sorted(
        matched,
        key=lambda item: (
            -float(item.get("spectral_distribution_divergence") or 0.0),
            float(item.get("collision_energy") or 0.0),
        ),
    )[0]
    ions = tuple(sorted({
        round(float(value), 4)
        for value in best.get("diagnostic_ions") or ()
        if math.isfinite(float(value)) and float(value) > 0
    }))[:12]
    return float(best["spectral_distribution_divergence"]), ions


def pair_fragment_separation(
    left_key: str,
    right_key: str,
    catalog: Mapping[str, Mapping[str, Any]],
) -> tuple[float, tuple[float, ...]]:
    left = _mz_values((catalog.get(left_key) or {}).get("diagnostic_fragments") or [])
    right = _mz_values((catalog.get(right_key) or {}).get("diagnostic_fragments") or [])
    if not left or not right:
        # Missing observables do not establish absence of a peak.  Withhold
        # separability credit and diagnostic ions unless both sides have data.
        return 0.0, ()
    union = left | right
    symmetric = left ^ right
    score = len(symmetric) / max(1, len(union))
    return min(1.0, max(0.05, score)), tuple(sorted(symmetric)[:12])
