"""EvoGraph-F: an additive, uncalibrated fragmentation-graph prototype.

The module deliberately lives outside the production search path.  It builds a
sparse, acyclic peak-to-peak neutral-loss graph from one immutable result row.
Formula assignments are only attempted for explicit singly charged adducts
whose charge-retention transform is defined below.  Normalized edge weights are
descriptive heuristics, not posterior probabilities or release evidence.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from matching.fragmentation_rules import COMMON_NEUTRAL_LOSSES
from matching.library_quality import (
    MONOISOTOPIC_MASS,
    formula_exact_mass,
    normalize_formula_text,
)
from utils.adduct import NEUTRAL_SHIFTS


SCHEMA_VERSION = "1.0.0"
MODULE_ID = "evograph_fragmentation"
PROTON_MASS = 1.007276466621
ELECTRON_MASS = 0.000548579909
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")

_CANONICAL_ADDUCT_ALIASES = {
    "[M+H]+": {"[M+H]+", "M+H"},
    "[M-H]-": {"[M-H]-", "M-H"},
    "[M]+": {"[M]+", "M+"},
    "[M]-": {"[M]-", "M-"},
    "[M+Na]+": {"[M+Na]+", "M+Na"},
    "[M+K]+": {"[M+K]+", "M+K"},
    "[M+NH4]+": {"[M+NH4]+", "M+NH4"},
    "[M+Cl]-": {"[M+Cl]-", "M+Cl"},
    "[M+FA-H]-": {"[M+FA-H]-", "[M+HCOO]-"},
    "[M+Hac-H]-": {"[M+Hac-H]-", "[M+CH3COO]-"},
}


def _clean_key(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


_EXACT_ADDUCT_LOOKUP = {
    _clean_key(alias): canonical
    for canonical, aliases in _CANONICAL_ADDUCT_ALIASES.items()
    for alias in aliases
    if alias in NEUTRAL_SHIFTS
}


def default_policy_path() -> Path:
    return Path(__file__).resolve().parent.parent / "config" / "evograph_policy.json"


def load_evograph_policy(path: Path | None = None) -> dict[str, Any]:
    payload = json.loads((path or default_policy_path()).read_text(encoding="utf-8"))
    if payload.get("module_id") != MODULE_ID:
        raise ValueError("EvoGraph policy module_id does not match runtime module")
    influence = payload.get("runtime_influence") or {}
    if any(bool(value) for value in influence.values()):
        raise ValueError("EvoGraph shadow policy cannot declare production influence")
    return payload


def canonical_explicit_adduct(value: Any) -> str | None:
    """Canonicalize only an exact table entry; never guess from polarity."""

    return _EXACT_ADDUCT_LOOKUP.get(_clean_key(value))


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _search_precursor_mz(row: Mapping[str, Any]) -> float | None:
    """Return the analytical mass used by the HRMS search.

    ``precursor_mz`` remains the mzML selected-ion/isolation target for raw
    traceability.  When continuous MS1 feature linking supplied a measured
    centroid, every downstream mass-consistency calculation must use that same
    search mass axis.
    """

    search_mz = _finite(row.get("search_precursor_mz"))
    if search_mz is not None and search_mz > 0:
        return search_mz
    return _finite(row.get("precursor_mz"))


def _charge_magnitude(value: Any) -> int | None:
    if value is None or value == "":
        return None
    match = re.search(r"\d+", str(value))
    if not match:
        return 1 if str(value).strip() in {"+", "-"} else None
    magnitude = int(match.group())
    return magnitude if magnitude > 0 else None


def _formula_counts(formula: Any) -> dict[str, int] | None:
    normalized = normalize_formula_text(formula)
    if not normalized:
        return None
    counts: dict[str, int] = {}
    position = 0
    for match in _FORMULA_TOKEN.finditer(normalized):
        if match.start() != position:
            return None
        element, raw_count = match.groups()
        if element not in MONOISOTOPIC_MASS:
            return None
        counts[element] = counts.get(element, 0) + (
            int(raw_count) if raw_count else 1
        )
        position = match.end()
    return counts if position == len(normalized) and counts else None


def _counts_formula(counts: Mapping[str, int]) -> str:
    positive = {str(k): int(v) for k, v in counts.items() if int(v) > 0}
    ordered: list[tuple[str, int]] = []
    if "C" in positive:
        ordered.append(("C", positive.pop("C")))
        if "H" in positive:
            ordered.append(("H", positive.pop("H")))
    ordered.extend((element, positive[element]) for element in sorted(positive))
    return "".join(
        f"{element}{count if count != 1 else ''}" for element, count in ordered
    )


def _mass_from_counts(counts: Mapping[str, int]) -> float:
    return sum(MONOISOTOPIC_MASS[element] * count for element, count in counts.items())


def _dbe(counts: Mapping[str, int]) -> float | None:
    if not counts.get("C"):
        return None
    halogens = sum(counts.get(element, 0) for element in ("F", "Cl", "Br", "I"))
    return 1.0 + counts.get("C", 0) + 0.5 * counts.get("N", 0) - 0.5 * (
        counts.get("H", 0) + halogens
    )


def _formula_valence_plausible(
    counts: Mapping[str, int], *, minimum_dbe: float = -0.5
) -> bool:
    """Conservatively reject impossible neutral-composition formulae.

    EvoGraph formulae are candidate neutral compositions under a declared
    charge-retention transform, not confirmed product structures.  The old
    subset-only check admitted losses such as CH27O and C2H17N.  A carbon-
    containing neutral composition must at least satisfy the ordinary DBE
    hydrogen bound; inorganic small losses (H2O, NH3, CO2, ...) remain valid.
    """

    if any(int(value) < 0 for value in counts.values()):
        return False
    if counts.get("C", 0):
        value = _dbe(counts)
        if value is None or value < minimum_dbe:
            return False
    return True


def _ion_to_neutral_mass(mz: float, adduct: str) -> float | None:
    if adduct == "[M+H]+":
        return mz - PROTON_MASS
    if adduct == "[M-H]-":
        return mz + PROTON_MASS
    if adduct == "[M]+":
        return mz + ELECTRON_MASS
    if adduct == "[M]-":
        return mz - ELECTRON_MASS
    return None


def _neutral_to_ion_mass(neutral_mass: float, adduct: str) -> float | None:
    if adduct == "[M+H]+":
        return neutral_mass + PROTON_MASS
    if adduct == "[M-H]-":
        return neutral_mass - PROTON_MASS
    if adduct == "[M]+":
        return neutral_mass - ELECTRON_MASS
    if adduct == "[M]-":
        return neutral_mass + ELECTRON_MASS
    return None


def _mass_tolerance(mass: float, policy: Mapping[str, Any], prefix: str) -> float:
    ppm = float(policy.get(f"{prefix}_mass_tolerance_ppm") or 0.0)
    absolute = float(policy.get(f"{prefix}_mass_tolerance_da") or 0.0)
    return max(absolute, abs(mass) * ppm * 1e-6)


def _enumerate_subformulas(
    parent_counts: Mapping[str, int],
    target_mass: float,
    tolerance: float,
    *,
    max_candidates: int,
    max_states: int,
    minimum_dbe: float,
    adduct: str,
) -> tuple[list[dict[str, Any]], bool, int]:
    """Enumerate bounded subformulae, deriving hydrogen at the leaf."""

    if target_mass <= 0:
        return [], False, 0
    hydrogen_max = int(parent_counts.get("H", 0))
    heavy = [
        element for element in parent_counts
        if element != "H" and parent_counts[element] > 0
    ]
    # Small-count and high-mass dimensions first make interval pruning sharper.
    heavy.sort(key=lambda element: (parent_counts[element], -MONOISOTOPIC_MASS[element], element))
    remaining_max = [0.0] * (len(heavy) + 1)
    for index in range(len(heavy) - 1, -1, -1):
        element = heavy[index]
        remaining_max[index] = remaining_max[index + 1] + (
            parent_counts[element] * MONOISOTOPIC_MASS[element]
        )

    visited = 0
    truncated = False
    found: list[dict[str, Any]] = []
    current: dict[str, int] = {}
    hydrogen_mass = MONOISOTOPIC_MASS["H"]

    def visit(index: int, partial_mass: float) -> None:
        nonlocal visited, truncated
        if visited >= max_states:
            truncated = True
            return
        visited += 1
        if index == len(heavy):
            residual = target_mass - partial_mass
            approximate_h = int(round(residual / hydrogen_mass))
            for hydrogen in range(max(0, approximate_h - 1), min(hydrogen_max, approximate_h + 1) + 1):
                if adduct == "[M-H]-" and hydrogen < 1:
                    continue
                counts = dict(current)
                if hydrogen:
                    counts["H"] = hydrogen
                exact_mass = partial_mass + hydrogen * hydrogen_mass
                error = exact_mass - target_mass
                if abs(error) > tolerance:
                    continue
                dbe = _dbe(counts)
                if dbe is not None and dbe < minimum_dbe:
                    continue
                formula = _counts_formula(counts)
                if not formula:
                    continue
                found.append({
                    "formula": formula,
                    "exact_neutral_mass": round(exact_mass, 8),
                    "target_neutral_mass": round(target_mass, 8),
                    "mass_error_da": round(error, 8),
                    "mass_error_ppm": round(error / target_mass * 1e6, 3),
                    "dbe": round(dbe, 3) if dbe is not None else None,
                })
            return

        element = heavy[index]
        element_mass = MONOISOTOPIC_MASS[element]
        maximum = int(parent_counts[element])
        maximum_after = remaining_max[index + 1] + hydrogen_max * hydrogen_mass
        lower = max(0, math.ceil((target_mass - tolerance - partial_mass - maximum_after) / element_mass))
        upper = min(maximum, math.floor((target_mass + tolerance - partial_mass) / element_mass))
        if lower > upper:
            return
        for count in range(lower, upper + 1):
            if count:
                current[element] = count
            else:
                current.pop(element, None)
            visit(index + 1, partial_mass + count * element_mass)
            if truncated:
                break
        current.pop(element, None)

    visit(0, 0.0)
    found.sort(key=lambda item: (abs(item["mass_error_da"]), item["formula"]))
    deduplicated: list[dict[str, Any]] = []
    seen: set[str] = set()
    for item in found:
        if item["formula"] in seen:
            continue
        seen.add(item["formula"])
        sigma = max(tolerance / 2.0, 1e-9)
        item["heuristic_weight"] = round(
            math.exp(-0.5 * (float(item["mass_error_da"]) / sigma) ** 2), 6
        )
        item["calibrated"] = False
        deduplicated.append(item)
        if len(deduplicated) >= max_candidates:
            break
    return deduplicated, truncated, visited


def _known_loss_catalogue() -> list[dict[str, Any]]:
    catalogue: dict[tuple[str, int], dict[str, Any]] = {}
    for rule_id, rule in COMMON_NEUTRAL_LOSSES.items():
        formula = normalize_formula_text(rule.formula)
        exact = formula_exact_mass(formula)
        if not formula or exact is None or exact <= 0:
            continue
        key = (formula, round(float(exact) * 100000))
        candidate = {
            "rule_id": str(rule_id),
            "formula": formula,
            "exact_mass": float(exact),
            "name": str(rule.name_en or rule.name or rule_id),
            "priority": int(rule.priority),
            "polarity": str(rule.polarity or "both"),
            "role": "soft_prior_only",
        }
        previous = catalogue.get(key)
        if previous is None or candidate["priority"] > previous["priority"]:
            catalogue[key] = candidate
    return sorted(catalogue.values(), key=lambda item: (-item["priority"], item["formula"]))


_KNOWN_LOSSES = _known_loss_catalogue()


def _match_known_loss(
    mass: float,
    tolerance: float,
    *,
    formula: str | None = None,
    ion_mode: str = "",
) -> dict[str, Any] | None:
    polarity = "positive" if str(ion_mode).upper().startswith("POS") else (
        "negative" if str(ion_mode).upper().startswith("NEG") else ""
    )
    matches = []
    for loss in _KNOWN_LOSSES:
        if formula and loss["formula"] != formula:
            continue
        error = float(loss["exact_mass"]) - mass
        if abs(error) > tolerance:
            continue
        polarity_ok = loss["polarity"] in {"", "both", polarity}
        matches.append((not polarity_ok, abs(error), -loss["priority"], loss, error))
    if not matches:
        return None
    _, _, _, loss, error = min(matches, key=lambda item: item[:3])
    return {
        **loss,
        "mass_error_da": round(error, 8),
        "calibrated": False,
    }


def _select_context(row: Mapping[str, Any], policy: Mapping[str, Any]) -> dict[str, Any]:
    match = row.get("match") if isinstance(row.get("match"), Mapping) else {}
    formula_value = row.get("parent_formula") or match.get("formula")
    formula = normalize_formula_text(formula_value)
    formula_source = (
        "candidate_independent_input" if row.get("parent_formula") else
        "selected_library_candidate" if formula else "unavailable"
    )
    formula_independent = bool(row.get("parent_formula") and row.get("parent_formula_independent"))

    context_candidates = []
    lcms_adduct = canonical_explicit_adduct(row.get("lcms_adduct"))
    lcms_neutral_mass = _finite(row.get("lcms_neutral_mass"))
    if lcms_adduct and lcms_neutral_mass is not None:
        context_candidates.append((lcms_adduct, "lcms_ion_identity", True, lcms_neutral_mass))
    query_adduct = canonical_explicit_adduct(row.get("query_adduct"))
    if query_adduct and bool(row.get("query_adduct_reliable")):
        context_candidates.append((query_adduct, "explicit_query_metadata", True, None))
    match_adduct = canonical_explicit_adduct(match.get("adduct"))
    if match_adduct:
        context_candidates.append((match_adduct, "selected_library_candidate", False, None))

    adduct, adduct_source, adduct_independent, neutral_reference = (
        context_candidates[0] if context_candidates else (None, "unavailable", False, None)
    )
    precursor_mz = _search_precursor_mz(row)
    charge = _charge_magnitude(row.get("precursor_charge"))
    if charge is None and adduct:
        charge = 1
    supported_adducts = set(
        (policy.get("formula_policy") or {}).get("supported_charge_retention_adducts") or []
    )
    reasons = []
    formula_status = "available"
    if not formula:
        formula_status = "unavailable"
        reasons.append("parent_formula_unavailable")
    if not adduct:
        formula_status = "unavailable"
        reasons.append("explicit_adduct_unavailable")
    elif adduct not in supported_adducts:
        formula_status = "ambiguous"
        reasons.append("charge_retention_transform_not_defined_for_adduct")
    if charge != int((policy.get("formula_policy") or {}).get("supported_charge_magnitude") or 1):
        formula_status = "unavailable"
        reasons.append("unsupported_precursor_charge")
    ion_mode = str(row.get("ion_mode") or "").strip().upper()
    adduct_mode = (
        "POS" if adduct and str(adduct).endswith("+") else
        "NEG" if adduct and str(adduct).endswith("-") else ""
    )
    if ion_mode.startswith(("POS", "NEG")) and adduct_mode and not ion_mode.startswith(adduct_mode):
        formula_status = "contradicted"
        reasons.append("ion_mode_adduct_polarity_inconsistent")

    formula_mass = formula_exact_mass(formula) if formula else None
    consistency = "not_tested"
    consistency_error_ppm = None
    if formula_mass is not None and precursor_mz is not None and adduct in supported_adducts:
        expected = _neutral_to_ion_mass(float(formula_mass), str(adduct))
        if expected is not None and expected > 0:
            error = precursor_mz - expected
            consistency_error_ppm = error / expected * 1e6
            tolerance = _mass_tolerance(expected, policy.get("formula_policy") or {}, "parent")
            consistency = "consistent" if abs(error) <= tolerance else "contradicted"
            if consistency == "contradicted":
                formula_status = "contradicted"
                reasons.append("parent_formula_precursor_mass_inconsistent")
    if formula_mass is not None and neutral_reference is not None:
        error = neutral_reference - float(formula_mass)
        tolerance = _mass_tolerance(float(formula_mass), policy.get("formula_policy") or {}, "parent")
        if abs(error) > tolerance:
            formula_status = "contradicted"
            consistency = "contradicted"
            reasons.append("lcms_neutral_mass_formula_inconsistent")

    return {
        "formula": formula or None,
        "formula_source": formula_source,
        "formula_candidate_independent": formula_independent,
        "adduct": adduct,
        "adduct_source": adduct_source,
        "adduct_candidate_independent": adduct_independent,
        "adduct_polarity": adduct_mode or None,
        "neutral_mass_reference": neutral_reference,
        "precursor_charge_magnitude": charge,
        "formula_assignment_status": formula_status,
        "mass_consistency": consistency,
        "mass_consistency_error_ppm": (
            round(consistency_error_ppm, 3) if consistency_error_ppm is not None else None
        ),
        "reasons": sorted(set(reasons)),
    }


def _clean_fragment_peaks(row: Mapping[str, Any], policy: Mapping[str, Any]) -> list[dict[str, float]]:
    input_policy = policy.get("input_policy") or {}
    raw = row.get("raw_peaks") or []
    processed = row.get("query_peaks") or []
    peaks = (
        raw or processed
        if input_policy.get("prefer_raw_peaks", True)
        else processed or raw
    )
    precursor = _search_precursor_mz(row)
    valid: list[tuple[float, float]] = []
    for peak in peaks:
        if not isinstance(peak, (list, tuple)) or len(peak) < 2:
            continue
        mz, intensity = _finite(peak[0]), _finite(peak[1])
        if mz is None or intensity is None or mz <= 0 or intensity <= 0:
            continue
        if precursor is not None and mz >= precursor - float(input_policy.get("minimum_loss_da") or 1.0):
            continue
        valid.append((mz, intensity))
    if not valid:
        return []
    base = max(intensity for _, intensity in valid)
    threshold = float(input_policy.get("minimum_relative_intensity") or 0.0)
    filtered = [peak for peak in valid if peak[1] / base >= threshold]
    filtered.sort(key=lambda peak: (-peak[1], peak[0]))
    filtered = filtered[: int(input_policy.get("max_peaks") or 24)]
    filtered.sort(key=lambda peak: -peak[0])
    return [
        {"mz": round(mz, 8), "intensity": round(intensity, 6), "relative_intensity": round(intensity / base, 6)}
        for mz, intensity in filtered
    ]


def _is_subset(child: Mapping[str, int], parent: Mapping[str, int]) -> bool:
    return all(int(count) <= int(parent.get(element, 0)) for element, count in child.items())


def _loss_counts(parent: Mapping[str, int], child: Mapping[str, int]) -> dict[str, int]:
    return {
        element: int(parent.get(element, 0)) - int(child.get(element, 0))
        for element in parent
        if int(parent.get(element, 0)) - int(child.get(element, 0)) > 0
    }


def _candidate_edge(
    parent: Mapping[str, Any],
    child: Mapping[str, Any],
    *,
    ion_mode: str,
    known_tolerance: float,
    mass_tolerance_ppm: float,
    minimum_loss_dbe: float,
) -> dict[str, Any] | None:
    observed_loss = float(parent["mz"]) - float(child["mz"])
    if observed_loss <= 0:
        return None
    best: dict[str, Any] | None = None
    for parent_formula in parent.get("formula_candidates") or []:
        parent_counts = _formula_counts(parent_formula.get("formula"))
        if not parent_counts:
            continue
        for child_formula in child.get("formula_candidates") or []:
            child_counts = _formula_counts(child_formula.get("formula"))
            if not child_counts or not _is_subset(child_counts, parent_counts):
                continue
            loss_counts = _loss_counts(parent_counts, child_counts)
            if not _formula_valence_plausible(
                loss_counts, minimum_dbe=minimum_loss_dbe
            ):
                continue
            loss_formula = _counts_formula(loss_counts)
            if not loss_formula:
                continue
            theoretical_loss = _mass_from_counts(loss_counts)
            error = theoretical_loss - observed_loss
            tolerance = max(known_tolerance, observed_loss * mass_tolerance_ppm * 1e-6)
            if abs(error) > tolerance:
                continue
            known = _match_known_loss(
                theoretical_loss,
                known_tolerance,
                formula=loss_formula,
                ion_mode=ion_mode,
            )
            sigma = max(tolerance / 2.0, 1e-9)
            fit = math.exp(-0.5 * (error / sigma) ** 2)
            raw_score = fit * float(parent_formula.get("heuristic_weight") or 1.0) * float(
                child_formula.get("heuristic_weight") or 1.0
            )
            if known:
                raw_score *= 1.0 + min(0.25, max(0, known["priority"]) / 50.0)
            candidate = {
                "parent_id": parent["node_id"],
                "child_id": child["node_id"],
                "kind": "formula_conserving",
                "observed_loss_da": round(observed_loss, 8),
                "loss_formula": loss_formula,
                "theoretical_loss_da": round(theoretical_loss, 8),
                "mass_error_da": round(error, 8),
                "mass_error_ppm": round(error / theoretical_loss * 1e6, 3),
                "parent_formula": parent_formula["formula"],
                "child_formula": child_formula["formula"],
                "common_loss_prior": known,
                "raw_heuristic_weight": round(raw_score, 8),
                "formula_conservation": True,
                "formula_semantics": "candidate_neutral_compositions_not_confirmed_structures",
            }
            if best is None or candidate["raw_heuristic_weight"] > best["raw_heuristic_weight"]:
                best = candidate
    if best is not None:
        return best

    known = _match_known_loss(observed_loss, known_tolerance, ion_mode=ion_mode)
    if not known:
        return None
    sigma = max(known_tolerance / 2.0, 1e-9)
    fit = math.exp(-0.5 * (float(known["mass_error_da"]) / sigma) ** 2)
    return {
        "parent_id": parent["node_id"],
        "child_id": child["node_id"],
        "kind": "mass_only_common_loss",
        "observed_loss_da": round(observed_loss, 8),
        "loss_formula": known["formula"],
        "theoretical_loss_da": round(float(known["exact_mass"]), 8),
        "mass_error_da": known["mass_error_da"],
        "mass_error_ppm": round(float(known["mass_error_da"]) / float(known["exact_mass"]) * 1e6, 3),
        "parent_formula": None,
        "child_formula": None,
        "common_loss_prior": known,
        "raw_heuristic_weight": round(0.55 * fit * (1.0 + min(0.2, known["priority"] / 50.0)), 8),
        "formula_conservation": None,
        "formula_semantics": "mass_only_common_loss_not_structure_assignment",
    }


def _cycle_count(nodes: list[Mapping[str, Any]], selected_edges: list[Mapping[str, Any]]) -> int:
    parents = {str(edge["child_id"]): str(edge["parent_id"]) for edge in selected_edges}
    cycles = 0
    for node in nodes:
        cursor = str(node["node_id"])
        seen: set[str] = set()
        while cursor in parents:
            if cursor in seen:
                cycles += 1
                break
            seen.add(cursor)
            cursor = parents[cursor]
    return cycles


def build_evograph_graph(
    row: Mapping[str, Any],
    *,
    result_index: int = 0,
    policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one deterministic shadow graph without mutating ``row``."""

    policy = dict(policy or load_evograph_policy())
    formula_policy = policy.get("formula_policy") or {}
    graph_policy = policy.get("graph_policy") or {}
    context = _select_context(row, policy)
    precursor_mz = _search_precursor_mz(row)
    fragments = _clean_fragment_peaks(row, policy)
    match = row.get("match") if isinstance(row.get("match"), Mapping) else {}

    base = {
        "schema_version": SCHEMA_VERSION,
        "module_id": MODULE_ID,
        "module_state": "shadow",
        "result_index": int(result_index),
        "spectrum_id": str(row.get("spectrum_id") or result_index),
        # Backward-compatible analytical precursor field.  The two explicit
        # fields below prevent the selected-ion target and measured MS1 search
        # centroid from being conflated in reports or audits.
        "precursor_mz": precursor_mz,
        "selected_ion_mz": _finite(row.get("precursor_mz")),
        "search_precursor_mz": precursor_mz,
        "search_precursor_mz_source": str(
            row.get("search_precursor_mz_source") or
            "ms2_selected_ion_mz_fallback"
        ),
        "lcms_feature_used_for_search": bool(
            row.get("lcms_feature_used_for_search")
        ),
        "ion_mode": str(row.get("ion_mode") or ""),
        "collision_energy": row.get("collision_energy"),
        "collision_energy_unit": str(row.get("collision_energy_unit") or ""),
        "instrument": str(row.get("instrument_type") or row.get("instrument") or ""),
        "selected_candidate": {
            "formula": context.get("formula"),
            "adduct": str(match.get("adduct") or ""),
            "source": context.get("formula_source"),
            "candidate_independent": context.get("formula_candidate_independent", False),
        },
        "context": context,
        "calibration": {
            **dict(policy.get("calibration") or {}),
            "log_bayes_factor": None,
            "expert_evidence_state": "unavailable",
        },
        "runtime_influence": dict(policy.get("runtime_influence") or {}),
        "production_ranking_changed": False,
        "legacy_confidence_changed": False,
        "high_promotions": 0,
    }
    if precursor_mz is None or precursor_mz <= 0:
        return {
            **base,
            "status": "unavailable",
            "reason": "invalid_precursor_mz",
            "nodes": [],
            "edges": [],
            "graph_audit": {"cycle_count": 0, "selected_edge_count": 0},
        }

    parent_candidates: list[dict[str, Any]] = []
    parent_counts = _formula_counts(context.get("formula"))
    formula_enabled = bool(
        parent_counts and context.get("formula_assignment_status") == "available"
    )
    if formula_enabled:
        exact = _mass_from_counts(parent_counts)
        parent_candidates = [{
            "formula": _counts_formula(parent_counts),
            "exact_neutral_mass": round(exact, 8),
            "target_neutral_mass": round(exact, 8),
            "mass_error_da": 0.0,
            "mass_error_ppm": 0.0,
            "dbe": round(_dbe(parent_counts), 3) if _dbe(parent_counts) is not None else None,
            "heuristic_weight": 1.0,
            "calibrated": False,
        }]

    nodes: list[dict[str, Any]] = [{
        "node_id": "p0",
        "node_type": "precursor",
        "mz": round(precursor_mz, 8),
        "intensity": _finite(row.get("precursor_intensity")),
        "relative_intensity": None,
        "formula_candidates": parent_candidates,
        "formula_enumeration_truncated": False,
        "formula_enumeration_states": 0,
    }]
    adduct = context.get("adduct")
    for index, peak in enumerate(fragments, start=1):
        candidates: list[dict[str, Any]] = []
        truncated = False
        states = 0
        if formula_enabled and adduct:
            neutral_target = _ion_to_neutral_mass(float(peak["mz"]), str(adduct))
            if neutral_target is not None:
                tolerance = _mass_tolerance(neutral_target, formula_policy, "fragment")
                candidates, truncated, states = _enumerate_subformulas(
                    parent_counts,
                    neutral_target,
                    tolerance,
                    max_candidates=int(formula_policy.get("max_candidates_per_peak") or 5),
                    max_states=int(formula_policy.get("max_enumeration_states_per_peak") or 15000),
                    minimum_dbe=float(formula_policy.get("minimum_dbe") or -0.5),
                    adduct=str(adduct),
                )
        nodes.append({
            "node_id": f"f{index}",
            "node_type": "fragment",
            **peak,
            "formula_candidates": candidates,
            "formula_enumeration_truncated": truncated,
            "formula_enumeration_states": states,
        })

    retained_edges: list[dict[str, Any]] = []
    selected_edges: list[dict[str, Any]] = []
    edge_status_counts = {"heuristically_supported": 0, "ambiguous": 0, "unassigned": 0}
    max_parents = int(graph_policy.get("max_parents_per_fragment") or 8)
    max_incoming = int(graph_policy.get("max_incoming_edges_retained") or 3)
    known_tolerance = float(graph_policy.get("known_loss_tolerance_da") or 0.015)
    edge_mass_tolerance_ppm = float(
        graph_policy.get("edge_mass_tolerance_ppm")
        or formula_policy.get("fragment_mass_tolerance_ppm")
        or 10.0
    )
    minimum_loss_dbe = float(
        graph_policy.get("minimum_loss_dbe")
        if graph_policy.get("minimum_loss_dbe") is not None
        else -0.5
    )
    null_weight = float(graph_policy.get("null_edge_weight") or 0.35)
    support_threshold = float(graph_policy.get("heuristic_support_threshold") or 0.55)
    ambiguity_margin = float(graph_policy.get("ambiguity_margin") or 0.15)

    for child_index in range(1, len(nodes)):
        child = nodes[child_index]
        higher = nodes[:child_index]
        root = higher[0]
        nonroot = sorted(
            higher[1:],
            key=lambda item: (float(item["mz"]) - float(child["mz"]), -float(item.get("relative_intensity") or 0.0)),
        )
        parents = [root, *nonroot[: max(0, max_parents - 1)]]
        candidates = [
            edge for edge in (
                _candidate_edge(
                    parent,
                    child,
                    ion_mode=str(row.get("ion_mode") or ""),
                    known_tolerance=known_tolerance,
                    mass_tolerance_ppm=edge_mass_tolerance_ppm,
                    minimum_loss_dbe=minimum_loss_dbe,
                )
                for parent in parents
            )
            if edge is not None
        ]
        candidates.sort(key=lambda edge: (-float(edge["raw_heuristic_weight"]), edge["parent_id"]))
        total_weight = null_weight + sum(float(edge["raw_heuristic_weight"]) for edge in candidates)
        for edge in candidates:
            edge["heuristic_marginal"] = round(float(edge["raw_heuristic_weight"]) / total_weight, 6)
            edge["calibrated"] = False
            edge["posterior_semantics"] = "normalized_local_heuristic_weight_not_probability"

        if not candidates:
            edge_status_counts["unassigned"] += 1
            child["incoming_edge_status"] = "unassigned"
            child["null_edge_heuristic_marginal"] = 1.0
            continue
        best = candidates[0]
        second = float(candidates[1]["heuristic_marginal"]) if len(candidates) > 1 else 0.0
        margin = float(best["heuristic_marginal"]) - second
        if float(best["heuristic_marginal"]) >= support_threshold and margin >= ambiguity_margin:
            status = "heuristically_supported"
        else:
            status = "ambiguous"
        edge_status_counts[status] += 1
        child["incoming_edge_status"] = status
        child["null_edge_heuristic_marginal"] = round(null_weight / total_weight, 6)
        best["selected"] = True
        best["selection_status"] = status
        selected_edges.append(best)
        for edge in candidates[:max_incoming]:
            edge.setdefault("selected", False)
            edge.setdefault("selection_status", "alternative")
            retained_edges.append(edge)

    cycle_count = _cycle_count(nodes, selected_edges)
    direction_violations = sum(
        float(next(node["mz"] for node in nodes if node["node_id"] == edge["parent_id"]))
        <= float(next(node["mz"] for node in nodes if node["node_id"] == edge["child_id"]))
        for edge in selected_edges
    )
    conservation_violations = sum(
        edge.get("formula_conservation") is False for edge in selected_edges
    )
    return {
        **base,
        "status": "completed" if fragments else "unavailable",
        "reason": None if fragments else "no_eligible_fragment_peaks",
        "graph_semantics": "sparse_charge_retention_neutral_loss_hypothesis",
        "formula_graph_available": formula_enabled,
        "nodes": nodes,
        "edges": retained_edges,
        "edge_status_counts": edge_status_counts,
        "graph_audit": {
            "node_count": len(nodes),
            "fragment_node_count": max(0, len(nodes) - 1),
            "candidate_edge_count": len(retained_edges),
            "selected_edge_count": len(selected_edges),
            "cycle_count": cycle_count,
            "direction_violations": direction_violations,
            "formula_conservation_violations": conservation_violations,
            "formula_enumeration_truncated_nodes": sum(
                bool(node.get("formula_enumeration_truncated")) for node in nodes
            ),
        },
    }


__all__ = [
    "MODULE_ID",
    "SCHEMA_VERSION",
    "build_evograph_graph",
    "canonical_explicit_adduct",
    "default_policy_path",
    "load_evograph_policy",
]
