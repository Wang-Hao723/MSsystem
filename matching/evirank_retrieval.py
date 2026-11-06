"""V12 EviRecall-MV: bounded, measured-spectrum candidate expansion.

The retriever is deliberately separate from V11 production ranking.  It reads
an explicitly versioned spectral-library snapshot in SQLite read-only mode,
uses only measured precursor/polarity/MS2/condition metadata, and returns a
shadow candidate-recall bundle.  Compound names are display metadata only;
ground truth, standards lists and expert/class/pathway outputs are rejected.

Two symmetrically processed matchms ModifiedCosine views are fused with
collision-energy-conditioned reciprocal-rank fusion (RRF).  In addition, a
strict source-paired consensus compares same-formula connectivities within the
same provenance stratum, collision-energy condition and spectral view.  Extra
copies from one stratum are median-collapsed and never create extra source
votes.  Provenance labels whose selected candidate-pair references have the
same normalized peak-content profile are collapsed across sources.  Missing
content fingerprints or other evidence causes abstention.  These utilities are not
probabilities and cannot alter V11 production output.
"""

from __future__ import annotations

import hashlib
import json
import math
import sqlite3
import statistics
from copy import deepcopy
from collections import Counter, OrderedDict, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any
from urllib.parse import quote

from rdkit import rdBase as _rdBase

from matching.chemical_identity import rdkit_identity_profile
from matching.library_polarity import classify_spectrum_polarity
from matching.library_quality import (
    normalize_formula_text,
    normalize_inchikey,
    source_tier_is_searchable,
)
from matching.conditioned_multilibrary import reference_condition
from matching.scores import (
    batch_modified_cosine_scores,
    batch_neutral_loss_cosine_scores,
)
from utils.spectrum import process_peaks, resolve_peak_processing_config


EVIRECALL_SCHEMA_VERSION = "2.0.0"
MODULE_ID = "v12_evirank_candidate_recall"
FORBIDDEN_REQUEST_FIELDS = {
    "compound_name",
    "expected_compound",
    "expected_compounds",
    "ground_truth",
    "ground_truth_identity",
    "standard_id",
    "standard_name",
    "expert_score",
    "classification",
    "biosynthetic_pathway",
    "candidate_pathway_alignment",
    "language_model_output",
}


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _text(value: Any) -> str:
    return str(value or "").strip()


def _polarity(value: Any) -> str | None:
    text = _text(value).upper()
    if text in {"POS", "POSITIVE", "+"}:
        return "POS"
    if text in {"NEG", "NEGATIVE", "-"}:
        return "NEG"
    return None


def default_evirecall_policy_path() -> Path:
    return Path(__file__).resolve().parent.parent / "config" / "v12_candidate_recall_policy.json"


def _validate_evirecall_policy(policy: Mapping[str, Any]) -> dict[str, Any]:
    payload = deepcopy(dict(policy))
    if payload.get("deployment_state") != "shadow":
        raise ValueError("V12 EviRecall-MV must remain shadow")
    if payload.get("display_name") != "V12 EviRecall-MV":
        raise ValueError("Unexpected EviRecall-MV display name")
    candidate_pool = payload.get("candidate_pool") or {}
    for field in (
        "target_total_connectivities",
        "maximum_additions_per_entity",
        "maximum_identities_per_scan_view",
        "minimum_library_fragment_peaks",
    ):
        if int(candidate_pool.get(field) or 0) <= 0:
            raise ValueError(f"candidate_pool.{field} must be positive")
    mass = payload.get("mass_retrieval") or {}
    if float(mass.get("direct_precursor_da") or 0.0) <= 0.0:
        raise ValueError("A positive direct precursor Da tolerance is required")
    if float(mass.get("direct_precursor_ppm") or 0.0) <= 0.0:
        raise ValueError("A positive direct precursor ppm tolerance is required")
    if mass.get("blind_window_widening_allowed") is not False:
        raise ValueError("Blind precursor-window widening must remain forbidden")
    if mass.get("blind_adduct_enumeration_allowed") is not False:
        raise ValueError("Blind adduct enumeration must remain forbidden")
    if mass.get("candidate_fitted_mass_correction_allowed") is not False:
        raise ValueError("Candidate-fitted mass correction must remain forbidden")
    views = payload.get("spectral_views") or []
    names = [_text(view.get("name")) for view in views]
    if not names or any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("Spectral views require unique non-empty names")
    for view in views:
        if view.get("similarity") not in {
            "matchms_modified_cosine_greedy",
            "matchms_neutral_losses_cosine",
        }:
            raise ValueError("EviRecall-MV view similarity is not supported")
        if float(view.get("fragment_tolerance_da") or 0.0) <= 0.0:
            raise ValueError("Each spectral view needs a positive fragment tolerance")
        resolve_peak_processing_config(view.get("peak_processing_preset"))
    fusion = payload.get("fusion") or {}
    if float(fusion.get("rrf_constant") or 0.0) <= 0.0:
        raise ValueError("A positive RRF constant is required")
    if fusion.get("score_is_probability") is not False:
        raise ValueError("EviRecall-MV RRF utility cannot be a probability")
    consensus = payload.get("source_paired_consensus") or {}
    if consensus.get("module_id") != "v12_source_paired_ce_consensus":
        raise ValueError("The source-paired consensus module id is not locked")
    if consensus.get("algorithm_id") != "mssystem-v12-source-paired-ce-consensus-0.1.0":
        raise ValueError("The source-paired consensus algorithm id is not locked")
    if consensus.get("deployment_state") != "shadow":
        raise ValueError("Source-paired consensus must remain shadow")
    if consensus.get("candidate_scope") != (
        "v11_same_formula_connectivities_else_bounded_evirecall_bridge"
    ):
        raise ValueError("The source-paired candidate scope is not locked")
    bounded_positive = (
        "maximum_v11_top1_margin",
        "minimum_same_formula_connectivities",
        "maximum_same_formula_connectivities",
        "maximum_expanded_same_formula_challengers",
        "minimum_common_provenance_strata",
        "minimum_unique_provenance_evidence_profiles",
        "minimum_distinct_collision_energy_conditions",
        "minimum_spectral_views",
        "minimum_paired_observations",
        "minimum_median_score_delta",
        "minimum_observation_win_fraction",
        "minimum_provenance_stratum_win_fraction",
        "minimum_unique_provenance_profile_win_fraction",
    )
    if any(float(consensus.get(field) or 0.0) <= 0.0 for field in bounded_positive):
        raise ValueError("Every source-paired consensus threshold must be positive")
    if float(consensus["maximum_v11_top1_margin"]) > 0.10:
        raise ValueError("The source-paired V11 margin gate cannot exceed 0.10")
    if int(consensus["minimum_same_formula_connectivities"]) < 2:
        raise ValueError("Source-paired consensus requires at least two connectivities")
    if int(consensus["maximum_same_formula_connectivities"]) < int(
        consensus["minimum_same_formula_connectivities"]
    ):
        raise ValueError("The source-paired connectivity bounds are inconsistent")
    if int(consensus["maximum_expanded_same_formula_challengers"]) >= int(
        consensus["maximum_same_formula_connectivities"]
    ):
        raise ValueError("The expanded challenger bridge exceeds the connectivity bound")
    for field in (
        "minimum_observation_win_fraction",
        "minimum_provenance_stratum_win_fraction",
        "minimum_unique_provenance_profile_win_fraction",
    ):
        if float(consensus[field]) > 1.0:
            raise ValueError(f"source_paired_consensus.{field} cannot exceed one")
    if consensus.get("technical_repeats_within_condition_view_stratum") != "median_score_delta":
        raise ValueError("Technical repeats must be median-collapsed")
    if consensus.get("reference_records_within_provenance_stratum") != "maximum_similarity_per_candidate":
        raise ValueError("Reference copies within one provenance stratum must be collapsed")
    if consensus.get("reference_content_fingerprint") != (
        "normalized_raw_peaks_base_peak_scaled_mz5_intensity6_sha256"
    ):
        raise ValueError("The source-paired reference-content fingerprint is not locked")
    if consensus.get("mirrored_provenance_strata") != (
        "collapse_identical_candidate_pair_reference_content_profiles"
    ):
        raise ValueError("Mirrored provenance strata must be content-collapsed")
    if consensus.get("unique_profile_delta") != (
        "median_of_member_provenance_stratum_medians"
    ):
        raise ValueError("The unique provenance-profile delta is not locked")
    if consensus.get("provenance_strata_are_independent_sources") is not False:
        raise ValueError("Provenance strata may not be labelled independent sources")
    if consensus.get("score_is_probability") is not False:
        raise ValueError("Source-paired score deltas cannot be probabilities")
    if consensus.get("production_effect") is not False:
        raise ValueError("Source-paired consensus may not affect production")
    library = payload.get("library_gate") or {}
    if library.get("require_stable_2d_connectivity_identity") is not True:
        raise ValueError("Expanded candidates require a stable connectivity identity")
    if library.get("same_polarity_or_unknown_library_polarity") is not True:
        raise ValueError("The same-polarity library gate must remain enabled")
    if int(library.get("maximum_sql_rows_per_precursor_window") or 0) <= 0:
        raise ValueError("A positive SQL row budget is required")
    if int(library.get("maximum_precursor_windows_in_memory") or 0) <= 0:
        raise ValueError("A positive bounded precursor-window cache size is required")
    tiers = [_text(value).lower() for value in library.get("searchable_source_tiers") or []]
    if not tiers or any(not source_tier_is_searchable(value) for value in tiers):
        raise ValueError("Every configured searchable source tier must pass the runtime gate")
    forbidden = set(payload.get("forbidden_inputs") or [])
    if not {
        "compound_name_for_retrieval_or_scoring",
        "ground_truth_identity",
        "standard_id",
        "expected_compound_list",
    } <= forbidden:
        raise ValueError("EviRecall-MV leakage barriers are incomplete")
    invariants = payload.get("invariants") or {}
    if invariants.get("candidate_names_are_display_only") is not True:
        raise ValueError("Candidate names must remain display-only")
    if invariants.get("production_ranking_changes") is not False:
        raise ValueError("EviRecall-MV may not change production ranking")
    if invariants.get("expanded_candidates_enter_shadow_candidate_union") is not True:
        raise ValueError("EviRecall candidates must enter the audited shadow union")
    if invariants.get("expanded_candidates_enter_residual_model_scoring") is not False:
        raise ValueError("Expanded candidates may not enter residual-model scoring")
    if invariants.get("candidate_recall_and_top1_are_reported_separately") is not True:
        raise ValueError("Candidate recall and Top-1 must be reported separately")
    return payload


def load_evirecall_policy(path: Path | None = None) -> dict[str, Any]:
    selected = path or default_evirecall_policy_path()
    payload = json.loads(selected.read_text(encoding="utf-8"))
    return _validate_evirecall_policy(payload)


def evirecall_policy_sha256(policy: Mapping[str, Any]) -> str:
    return _sha256(policy)


def _forbidden_request_paths(value: Any, prefix: str = "") -> list[str]:
    paths: list[str] = []
    if isinstance(value, Mapping):
        for raw_key, child in value.items():
            key = str(raw_key)
            path = f"{prefix}.{key}" if prefix else key
            if key in FORBIDDEN_REQUEST_FIELDS:
                paths.append(path)
            paths.extend(_forbidden_request_paths(child, path))
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for index, child in enumerate(value):
            paths.extend(_forbidden_request_paths(child, f"{prefix}[{index}]"))
    return paths


def _stable_connectivity(record: Mapping[str, Any]) -> tuple[str | None, str]:
    key = normalize_inchikey(record.get("inchikey"))
    if key:
        return key.split("-", 1)[0], "reported_full_inchikey_connectivity"
    smiles = _text(record.get("smiles"))
    if not smiles or smiles.lower() in {"n/a", "na", "nan", "none", "null", "-"}:
        return None, "missing_stable_structure_identity"
    blocker = _rdBase.BlockLogs()
    try:
        profile = rdkit_identity_profile(
            smiles,
            declared_formula=record.get("molecular_formula"),
        )
    finally:
        del blocker
    if not profile.get("valid"):
        return None, str(profile.get("reason") or "invalid_structure_identity")
    return _text(profile.get("connectivity_key")) or None, "rdkit_2d_connectivity"


def _condition_id(scan: Mapping[str, Any]) -> str:
    value = _finite(scan.get("collision_energy"))
    unit = _text(scan.get("collision_energy_unit")).lower() or "unspecified"
    mode = _text(scan.get("collision_energy_mode")).lower() or "unspecified"
    if value is None:
        return f"ce:unknown|{unit}|{mode}"
    return f"ce:{value:.6g}|{unit}|{mode}"


def _existing_connectivities(candidates: Sequence[Mapping[str, Any]]) -> set[str]:
    identities: set[str] = set()
    for candidate in candidates:
        display = candidate.get("candidate") or candidate
        key = normalize_inchikey(display.get("inchikey"))
        if key:
            identities.add(key.split("-", 1)[0])
            continue
        identity_key = _text(candidate.get("identity_key"))
        if identity_key.startswith("std2d:"):
            identities.add(identity_key.split(":", 1)[1])
        elif identity_key.startswith("ik14:"):
            identities.add(identity_key.split(":", 1)[1])
    return identities


def _candidate_connectivity(candidate: Mapping[str, Any]) -> str | None:
    direct = _text(candidate.get("connectivity_identity"))
    if direct:
        return direct
    identity_key = _text(candidate.get("identity_key"))
    if identity_key.startswith(("std2d:", "ik14:")):
        return identity_key.split(":", 1)[1] or None
    display = candidate.get("candidate") or candidate
    key = normalize_inchikey(display.get("inchikey"))
    return key.split("-", 1)[0] if key else None


def _candidate_identity_key(candidate: Mapping[str, Any], connectivity: str) -> str:
    identity_key = _text(candidate.get("identity_key"))
    if identity_key.startswith(("std2d:", "ik14:")):
        return identity_key
    return f"std2d:{connectivity}"


def _candidate_formula(candidate: Mapping[str, Any]) -> str:
    display = candidate.get("candidate") or candidate
    return normalize_formula_text(display.get("formula") or display.get("molecular_formula"))


def _provenance_stratum(value: Any) -> str | None:
    """Canonicalize one declared dataset source without claiming independence."""

    text = " ".join(_text(value).split()).casefold()
    return text or None


def _reference_content_fingerprint(record: Mapping[str, Any]) -> str | None:
    """Hash normalized reference peaks so mirrored source labels can collapse.

    The fingerprint deliberately excludes names, identifiers and provenance.
    Intensities are normalized to the base peak and peak order is canonicalized,
    making copies with different absolute intensity scales compare identically.
    """

    normalized: list[tuple[float, float]] = []
    for peak in record.get("peaks") or []:
        if not isinstance(peak, Sequence) or isinstance(peak, (str, bytes)):
            continue
        if len(peak) < 2:
            continue
        mz = _finite(peak[0])
        intensity = _finite(peak[1])
        if mz is None or intensity is None or mz <= 0.0 or intensity <= 0.0:
            continue
        normalized.append((mz, intensity))
    if not normalized:
        return None
    base_peak = max(intensity for _, intensity in normalized)
    if base_peak <= 0.0:
        return None
    canonical_peaks = sorted(
        (round(mz, 5), round(intensity / base_peak, 6))
        for mz, intensity in normalized
    )
    return _sha256(canonical_peaks)


def _source_paired_consensus(
    existing_candidates: Sequence[Mapping[str, Any]],
    scan_results: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
    *,
    expanded_candidate_connectivities: Sequence[str] = (),
) -> dict[str, Any]:
    """Compare close same-formula challengers in matched provenance strata.

    The routine never sees names, standards lists or ground truth.  Each raw
    delta is matched on query condition, spectral view and declared dataset
    provenance.  Technical repeats are median-collapsed before gate counts are
    computed, and multiple reference copies within one provenance stratum have
    already been reduced to one maximum-similarity score in ``retrieve_scan``.
    """

    config = policy["source_paired_consensus"]
    thresholds = {
        "maximum_v11_top1_margin": float(config["maximum_v11_top1_margin"]),
        "minimum_same_formula_connectivities": int(
            config["minimum_same_formula_connectivities"]
        ),
        "maximum_same_formula_connectivities": int(
            config["maximum_same_formula_connectivities"]
        ),
        "maximum_expanded_same_formula_challengers": int(
            config["maximum_expanded_same_formula_challengers"]
        ),
        "minimum_common_provenance_strata": int(
            config["minimum_common_provenance_strata"]
        ),
        "minimum_unique_provenance_evidence_profiles": int(
            config["minimum_unique_provenance_evidence_profiles"]
        ),
        "minimum_distinct_collision_energy_conditions": int(
            config["minimum_distinct_collision_energy_conditions"]
        ),
        "minimum_spectral_views": int(config["minimum_spectral_views"]),
        "minimum_paired_observations": int(config["minimum_paired_observations"]),
        "minimum_median_score_delta": float(config["minimum_median_score_delta"]),
        "minimum_observation_win_fraction": float(
            config["minimum_observation_win_fraction"]
        ),
        "minimum_provenance_stratum_win_fraction": float(
            config["minimum_provenance_stratum_win_fraction"]
        ),
        "minimum_unique_provenance_profile_win_fraction": float(
            config["minimum_unique_provenance_profile_win_fraction"]
        ),
    }

    base_candidates = [
        candidate
        for candidate in existing_candidates
        if isinstance(candidate, Mapping)
        and not candidate.get("vetoed")
        and _candidate_connectivity(candidate)
    ]
    baseline = base_candidates[0] if base_candidates else None
    baseline_connectivity = (
        _candidate_connectivity(baseline) if baseline is not None else None
    )
    baseline_identity = (
        _candidate_identity_key(baseline, baseline_connectivity)
        if baseline is not None and baseline_connectivity
        else None
    )
    result: dict[str, Any] = {
        "schema_version": "1.0.0",
        "module_id": config["module_id"],
        "algorithm_id": config["algorithm_id"],
        "deployment_state": "shadow",
        "status": "not_applicable",
        "reason": "no_v11_baseline_candidate",
        "evaluated": False,
        "baseline_identity_key": baseline_identity,
        "baseline_connectivity_identity": baseline_connectivity,
        "baseline_margin": None,
        "formula": None,
        "same_formula_connectivity_count": 0,
        "native_same_formula_connectivity_count": 0,
        "expanded_same_formula_challenger_count": 0,
        "candidate_scope": config["candidate_scope"],
        "thresholds": thresholds,
        "comparisons": [],
        "selected_identity_key": baseline_identity,
        "selected_connectivity_identity": baseline_connectivity,
        "ranking_changed_from_v11": False,
        "technical_repeats_collapsed_by_median": True,
        "reference_copies_collapsed_within_provenance_stratum": True,
        "mirrored_reference_content_collapsed_across_provenance_strata": True,
        "provenance_strata_are_independent_sources": False,
        "ground_truth_used_at_inference": False,
        "candidate_names_affect_scoring": False,
        "score_is_probability": False,
        "production_effect": False,
    }
    if baseline is None or baseline_connectivity is None:
        return result
    if len(base_candidates) < 2:
        result["reason"] = "v11_runner_up_unavailable"
        return result
    baseline_score = _finite(baseline.get("entity_score"))
    runner_score = _finite(base_candidates[1].get("entity_score"))
    if baseline_score is None or runner_score is None or baseline_score < runner_score:
        result["reason"] = "invalid_v11_top1_margin"
        return result
    margin = baseline_score - runner_score
    formula = _candidate_formula(baseline)
    result["baseline_margin"] = round(margin, 8)
    result["formula"] = formula or None
    if not formula:
        result["reason"] = "v11_top1_formula_unavailable"
        return result

    candidates_by_connectivity: dict[str, Mapping[str, Any]] = {}
    identity_keys: dict[str, str] = {}
    vetoed_connectivities = {
        connectivity
        for candidate in existing_candidates
        if isinstance(candidate, Mapping)
        and candidate.get("vetoed")
        and (connectivity := _candidate_connectivity(candidate))
    }
    for candidate in base_candidates:
        connectivity = _candidate_connectivity(candidate)
        if not connectivity:
            continue
        candidates_by_connectivity.setdefault(connectivity, candidate)
        identity_keys.setdefault(
            connectivity, _candidate_identity_key(candidate, connectivity)
        )
    native_connectivities = set(candidates_by_connectivity)
    native_same_formula = {
        connectivity
        for connectivity, candidate in candidates_by_connectivity.items()
        if _candidate_formula(candidate) == formula
    }
    result["native_same_formula_connectivity_count"] = len(native_same_formula)
    if len(native_same_formula) >= 2:
        allowed_expanded_connectivities: set[str] = set()
        candidate_scope_applied = "v11_native_same_formula_only"
    else:
        ranked_same_formula_expansions: list[str] = []
        for raw_connectivity in expanded_candidate_connectivities:
            connectivity = _text(raw_connectivity)
            if not connectivity or connectivity in ranked_same_formula_expansions:
                continue
            candidate = next(
                (
                    (scan_result.get("candidates") or {}).get(connectivity)
                    for scan_result in scan_results
                    if isinstance(
                        (scan_result.get("candidates") or {}).get(connectivity),
                        Mapping,
                    )
                ),
                None,
            )
            if isinstance(candidate, Mapping) and _candidate_formula(candidate) == formula:
                ranked_same_formula_expansions.append(connectivity)
        allowed_expanded_connectivities = set(
            ranked_same_formula_expansions[
                : thresholds["maximum_expanded_same_formula_challengers"]
            ]
        )
        candidate_scope_applied = "v11_baseline_plus_bounded_evirecall_bridge"
    result["candidate_scope_applied"] = candidate_scope_applied
    result["expanded_bridge_connectivities_considered"] = sorted(
        allowed_expanded_connectivities
    )
    for scan_result in scan_results:
        for connectivity, candidate in (scan_result.get("candidates") or {}).items():
            if not isinstance(candidate, Mapping):
                continue
            connectivity_text = _text(connectivity)
            if (
                not connectivity_text
                or connectivity_text in vetoed_connectivities
                or (
                    connectivity_text not in native_connectivities
                    and connectivity_text not in allowed_expanded_connectivities
                )
            ):
                continue
            candidates_by_connectivity.setdefault(connectivity_text, candidate)
            identity_keys.setdefault(connectivity_text, f"std2d:{connectivity_text}")

    same_formula = {
        connectivity: candidate
        for connectivity, candidate in candidates_by_connectivity.items()
        if _candidate_formula(candidate) == formula
    }
    result["expanded_same_formula_challenger_count"] = len(
        set(same_formula) - native_connectivities
    )
    result["same_formula_connectivity_count"] = len(same_formula)
    minimum_connectivities = thresholds["minimum_same_formula_connectivities"]
    maximum_connectivities = thresholds["maximum_same_formula_connectivities"]
    if len(same_formula) < minimum_connectivities:
        result["reason"] = "insufficient_same_formula_connectivities"
        return result
    if len(same_formula) > maximum_connectivities:
        result["reason"] = "same_formula_candidate_set_exceeds_bound"
        return result
    if margin > thresholds["maximum_v11_top1_margin"]:
        result.update({
            "status": "baseline_retained",
            "reason": "v11_margin_above_gate",
        })
        return result

    comparisons: list[dict[str, Any]] = []
    for challenger_connectivity in sorted(same_formula):
        if challenger_connectivity == baseline_connectivity:
            continue
        collapsed_groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        platform_groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        conditioned_groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        reference_pairs: dict[
            tuple[str, str, str], set[tuple[str, str]]
        ] = defaultdict(set)
        conditioned_reference_pairs: dict[
            tuple[str, str, str], set[tuple[str, str]]
        ] = defaultdict(set)
        platform_reference_pairs: dict[
            tuple[str, str, str], set[tuple[str, str]]
        ] = defaultdict(set)
        condition_tier_counts: Counter[str] = Counter()
        for scan_result in scan_results:
            condition = _text(scan_result.get("condition_id"))
            scan_candidates = scan_result.get("candidates") or {}
            baseline_scan = scan_candidates.get(baseline_connectivity) or {}
            challenger_scan = scan_candidates.get(challenger_connectivity) or {}
            baseline_views = baseline_scan.get("source_view_scores") or {}
            challenger_views = challenger_scan.get("source_view_scores") or {}
            baseline_fingerprint_views = (
                baseline_scan.get("source_reference_content_fingerprints") or {}
            )
            challenger_fingerprint_views = (
                challenger_scan.get("source_reference_content_fingerprints") or {}
            )
            baseline_reference_views = (
                baseline_scan.get("source_reference_records") or {}
            )
            challenger_reference_views = (
                challenger_scan.get("source_reference_records") or {}
            )
            query_condition = scan_result.get("query_condition") or {}
            for view_name in sorted(set(baseline_views) & set(challenger_views)):
                baseline_sources = baseline_views.get(view_name) or {}
                challenger_sources = challenger_views.get(view_name) or {}
                baseline_fingerprints = (
                    baseline_fingerprint_views.get(view_name) or {}
                )
                challenger_fingerprints = (
                    challenger_fingerprint_views.get(view_name) or {}
                )
                baseline_references = baseline_reference_views.get(view_name) or {}
                challenger_references = challenger_reference_views.get(view_name) or {}
                for stratum in sorted(set(baseline_sources) & set(challenger_sources)):
                    baseline_value = _finite(baseline_sources.get(stratum))
                    challenger_value = _finite(challenger_sources.get(stratum))
                    baseline_fingerprint = _text(
                        baseline_fingerprints.get(stratum)
                    )
                    challenger_fingerprint = _text(
                        challenger_fingerprints.get(stratum)
                    )
                    if (
                        baseline_value is None
                        or challenger_value is None
                        or not baseline_fingerprint
                        or not challenger_fingerprint
                    ):
                        continue
                    group_key = (condition, view_name, stratum)
                    collapsed_groups[group_key].append(
                        challenger_value - baseline_value
                    )
                    reference_pairs[group_key].add(
                        (baseline_fingerprint, challenger_fingerprint)
                    )
                    baseline_condition = reference_condition(
                        query_condition,
                        baseline_references.get(stratum) or {},
                        near_delta_by_unit={"eV": 5.0, "NCE": 5.0},
                    )
                    challenger_condition = reference_condition(
                        query_condition,
                        challenger_references.get(stratum) or {},
                        near_delta_by_unit={"eV": 5.0, "NCE": 5.0},
                    )
                    pair_tier = (
                        f"{baseline_condition.get('tier') or 'unknown'}|"
                        f"{challenger_condition.get('tier') or 'unknown'}"
                    )
                    condition_tier_counts[pair_tier] += 1
                    if (
                        baseline_condition.get("platform_anchor_eligible")
                        and challenger_condition.get("platform_anchor_eligible")
                    ):
                        platform_groups[group_key].append(
                            challenger_value - baseline_value
                        )
                        platform_reference_pairs[group_key].add(
                            (baseline_fingerprint, challenger_fingerprint)
                        )
                    if (
                        baseline_condition.get("eligible_for_conditioned_consensus")
                        and challenger_condition.get("eligible_for_conditioned_consensus")
                    ):
                        conditioned_groups[group_key].append(
                            challenger_value - baseline_value
                        )
                        conditioned_reference_pairs[group_key].add(
                            (baseline_fingerprint, challenger_fingerprint)
                        )

        observations = [
            (
                condition,
                view_name,
                stratum,
                float(statistics.median(values)),
                len(values),
                tuple(sorted(reference_pairs[(condition, view_name, stratum)])),
            )
            for (condition, view_name, stratum), values in sorted(collapsed_groups.items())
        ]
        deltas = [item[3] for item in observations]
        strata = sorted({item[2] for item in observations})
        conditions = sorted({item[0] for item in observations})
        views = sorted({item[1] for item in observations})
        stratum_deltas: dict[str, list[float]] = defaultdict(list)
        stratum_reference_pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
        for _, _, stratum, delta, _, pairs in observations:
            stratum_deltas[stratum].append(delta)
            stratum_reference_pairs[stratum].update(pairs)
        stratum_medians = {
            stratum: float(statistics.median(values))
            for stratum, values in stratum_deltas.items()
        }
        profile_strata: dict[str, list[str]] = defaultdict(list)
        for stratum, pairs in stratum_reference_pairs.items():
            if pairs:
                profile_strata[_sha256(sorted(pairs))].append(stratum)
        profile_medians = {
            profile: float(
                statistics.median(stratum_medians[stratum] for stratum in members)
            )
            for profile, members in profile_strata.items()
        }
        conditioned_observations = [
            (
                condition,
                view_name,
                stratum,
                float(statistics.median(values)),
                len(values),
                tuple(
                    sorted(
                        conditioned_reference_pairs[(condition, view_name, stratum)]
                    )
                ),
            )
            for (condition, view_name, stratum), values in sorted(
                conditioned_groups.items()
            )
        ]
        platform_observations = [
            (
                condition,
                view_name,
                stratum,
                float(statistics.median(values)),
                len(values),
                tuple(
                    sorted(platform_reference_pairs[(condition, view_name, stratum)])
                ),
            )
            for (condition, view_name, stratum), values in sorted(
                platform_groups.items()
            )
        ]
        platform_deltas = [item[3] for item in platform_observations]
        platform_strata = sorted({item[2] for item in platform_observations})
        platform_conditions = sorted({item[0] for item in platform_observations})
        platform_views = sorted({item[1] for item in platform_observations})
        platform_stratum_deltas: dict[str, list[float]] = defaultdict(list)
        platform_stratum_pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
        for _, _, stratum, delta, _, pairs in platform_observations:
            platform_stratum_deltas[stratum].append(delta)
            platform_stratum_pairs[stratum].update(pairs)
        platform_stratum_medians = {
            stratum: float(statistics.median(values))
            for stratum, values in platform_stratum_deltas.items()
        }
        platform_profile_strata: dict[str, list[str]] = defaultdict(list)
        for stratum, pairs in platform_stratum_pairs.items():
            if pairs:
                platform_profile_strata[_sha256(sorted(pairs))].append(stratum)
        platform_profile_medians = {
            profile: float(
                statistics.median(
                    platform_stratum_medians[stratum] for stratum in members
                )
            )
            for profile, members in platform_profile_strata.items()
        }
        platform_median_delta = (
            float(statistics.median(platform_deltas)) if platform_deltas else None
        )
        platform_observation_wins = sum(delta > 0.0 for delta in platform_deltas)
        platform_stratum_wins = sum(
            value > 0.0 for value in platform_stratum_medians.values()
        )
        platform_profile_wins = sum(
            value > 0.0 for value in platform_profile_medians.values()
        )
        conditioned_deltas = [item[3] for item in conditioned_observations]
        conditioned_strata = sorted({item[2] for item in conditioned_observations})
        conditioned_conditions = sorted(
            {item[0] for item in conditioned_observations}
        )
        conditioned_views = sorted({item[1] for item in conditioned_observations})
        conditioned_stratum_deltas: dict[str, list[float]] = defaultdict(list)
        conditioned_stratum_pairs: dict[
            str, set[tuple[str, str]]
        ] = defaultdict(set)
        for _, _, stratum, delta, _, pairs in conditioned_observations:
            conditioned_stratum_deltas[stratum].append(delta)
            conditioned_stratum_pairs[stratum].update(pairs)
        conditioned_stratum_medians = {
            stratum: float(statistics.median(values))
            for stratum, values in conditioned_stratum_deltas.items()
        }
        conditioned_profile_strata: dict[str, list[str]] = defaultdict(list)
        for stratum, pairs in conditioned_stratum_pairs.items():
            if pairs:
                conditioned_profile_strata[_sha256(sorted(pairs))].append(stratum)
        conditioned_profile_medians = {
            profile: float(
                statistics.median(
                    conditioned_stratum_medians[stratum] for stratum in members
                )
            )
            for profile, members in conditioned_profile_strata.items()
        }
        conditioned_median_delta = (
            float(statistics.median(conditioned_deltas))
            if conditioned_deltas
            else None
        )
        conditioned_observation_wins = sum(
            delta > 0.0 for delta in conditioned_deltas
        )
        conditioned_stratum_wins = sum(
            value > 0.0 for value in conditioned_stratum_medians.values()
        )
        conditioned_profile_wins = sum(
            value > 0.0 for value in conditioned_profile_medians.values()
        )
        median_delta = float(statistics.median(deltas)) if deltas else None
        observation_wins = sum(delta > 0.0 for delta in deltas)
        observation_win_fraction = observation_wins / len(deltas) if deltas else 0.0
        stratum_wins = sum(value > 0.0 for value in stratum_medians.values())
        stratum_win_fraction = stratum_wins / len(stratum_medians) if stratum_medians else 0.0
        profile_wins = sum(value > 0.0 for value in profile_medians.values())
        profile_win_fraction = (
            profile_wins / len(profile_medians) if profile_medians else 0.0
        )
        failed_gates: list[str] = []
        if len(strata) < thresholds["minimum_common_provenance_strata"]:
            failed_gates.append("insufficient_common_provenance_strata")
        if len(profile_medians) < thresholds[
            "minimum_unique_provenance_evidence_profiles"
        ]:
            failed_gates.append("insufficient_unique_provenance_evidence_profiles")
        if len(conditions) < thresholds["minimum_distinct_collision_energy_conditions"]:
            failed_gates.append("insufficient_distinct_collision_energy_conditions")
        if len(views) < thresholds["minimum_spectral_views"]:
            failed_gates.append("insufficient_spectral_views")
        if len(deltas) < thresholds["minimum_paired_observations"]:
            failed_gates.append("insufficient_paired_observations")
        if median_delta is None or median_delta < thresholds["minimum_median_score_delta"]:
            failed_gates.append("median_score_delta_below_gate")
        if observation_win_fraction < thresholds["minimum_observation_win_fraction"]:
            failed_gates.append("observation_win_fraction_below_gate")
        if stratum_win_fraction < thresholds["minimum_provenance_stratum_win_fraction"]:
            failed_gates.append("provenance_stratum_win_fraction_below_gate")
        if profile_win_fraction < thresholds[
            "minimum_unique_provenance_profile_win_fraction"
        ]:
            failed_gates.append("unique_provenance_profile_win_fraction_below_gate")
        comparisons.append({
            "challenger_identity_key": identity_keys[challenger_connectivity],
            "challenger_connectivity_identity": challenger_connectivity,
            "paired_observation_count": len(deltas),
            "technical_repeat_observation_count": sum(item[4] for item in observations),
            "common_provenance_stratum_count": len(strata),
            "common_provenance_strata": strata,
            "unique_provenance_evidence_profile_count": len(profile_medians),
            "unique_provenance_evidence_profile_sizes": sorted(
                len(members) for members in profile_strata.values()
            ),
            "distinct_collision_energy_condition_count": len(conditions),
            "spectral_view_count": len(views),
            "spectral_views": views,
            "median_score_delta": round(median_delta, 8) if median_delta is not None else None,
            "observation_win_count": observation_wins,
            "observation_win_fraction": round(observation_win_fraction, 8),
            "provenance_stratum_win_count": stratum_wins,
            "provenance_stratum_win_fraction": round(stratum_win_fraction, 8),
            "unique_provenance_profile_win_count": profile_wins,
            "unique_provenance_profile_win_fraction": round(
                profile_win_fraction, 8
            ),
            "reference_condition_audit": {
                "definition": (
                    "both_candidate_references_same_query_analyzer_family_and_"
                    "within_5_units_in_same_explicit_eV_or_NCE_domain"
                ),
                "paired_observation_count": len(conditioned_deltas),
                "common_provenance_stratum_count": len(conditioned_strata),
                "distinct_query_collision_energy_condition_count": len(
                    conditioned_conditions
                ),
                "spectral_view_count": len(conditioned_views),
                "unique_provenance_evidence_profile_count": len(
                    conditioned_profile_medians
                ),
                "median_score_delta": (
                    round(conditioned_median_delta, 8)
                    if conditioned_median_delta is not None
                    else None
                ),
                "observation_win_count": conditioned_observation_wins,
                "observation_win_fraction": round(
                    conditioned_observation_wins / len(conditioned_deltas), 8
                ) if conditioned_deltas else 0.0,
                "provenance_stratum_win_count": conditioned_stratum_wins,
                "provenance_stratum_win_fraction": round(
                    conditioned_stratum_wins / len(conditioned_stratum_medians), 8
                ) if conditioned_stratum_medians else 0.0,
                "unique_provenance_profile_win_count": conditioned_profile_wins,
                "unique_provenance_profile_win_fraction": round(
                    conditioned_profile_wins / len(conditioned_profile_medians), 8
                ) if conditioned_profile_medians else 0.0,
                "condition_pair_tier_counts": dict(
                    sorted(condition_tier_counts.items())
                ),
                "cross_unit_conversion_applied": False,
            },
            "reference_platform_audit": {
                "definition": (
                    "both_candidate_references_match_the_same_available_"
                    "instrument_domain_model_vendor_analyzer_or_explicitly_"
                    "labelled_analyzer_family_fallback;_"
                    "missing_reference_CE_remains_missing_and_is_not_treated_as_"
                    "an_exact_CE_match"
                ),
                "paired_observation_count": len(platform_deltas),
                "common_provenance_stratum_count": len(platform_strata),
                "distinct_query_collision_energy_condition_count": len(
                    platform_conditions
                ),
                "spectral_view_count": len(platform_views),
                "unique_provenance_evidence_profile_count": len(
                    platform_profile_medians
                ),
                "median_score_delta": (
                    round(platform_median_delta, 8)
                    if platform_median_delta is not None
                    else None
                ),
                "observation_win_count": platform_observation_wins,
                "observation_win_fraction": round(
                    platform_observation_wins / len(platform_deltas), 8
                ) if platform_deltas else 0.0,
                "provenance_stratum_win_count": platform_stratum_wins,
                "provenance_stratum_win_fraction": round(
                    platform_stratum_wins / len(platform_stratum_medians), 8
                ) if platform_stratum_medians else 0.0,
                "unique_provenance_profile_win_count": platform_profile_wins,
                "unique_provenance_profile_win_fraction": round(
                    platform_profile_wins / len(platform_profile_medians), 8
                ) if platform_profile_medians else 0.0,
                "reference_collision_energy_required_for_this_count": False,
                "reference_collision_energy_missing_treated_as_match": False,
                "cross_instrument_conversion_applied": False,
                "cross_vendor_ce_equivalence_assumed": False,
            },
            "failed_gates": failed_gates,
            "qualified": not failed_gates,
            "score_is_probability": False,
        })

    result["comparisons"] = comparisons
    result["evaluated"] = bool(comparisons)
    qualified = [record for record in comparisons if record["qualified"]]
    if not qualified:
        result.update({
            "status": "baseline_retained",
            "reason": (
                "no_qualified_source_paired_challenger"
                if comparisons
                else "no_same_formula_challenger_with_paired_evidence"
            ),
        })
        return result
    selected = min(
        qualified,
        key=lambda record: (
            -float(record["median_score_delta"]),
            -float(record["observation_win_fraction"]),
            -float(record["unique_provenance_profile_win_fraction"]),
            -float(record["provenance_stratum_win_fraction"]),
            record["challenger_identity_key"],
        ),
    )
    result.update({
        "status": "qualified_switch",
        "reason": "qualified_source_paired_same_formula_challenger",
        "selected_identity_key": selected["challenger_identity_key"],
        "selected_connectivity_identity": selected["challenger_connectivity_identity"],
        "ranking_changed_from_v11": True,
    })
    return result


def build_candidate_recall_bundle(
    v11_payload: Mapping[str, Any],
    source_rows: Sequence[Mapping[str, Any]],
    *,
    retriever: "SQLiteCandidateRecallRetriever",
    detection_method_id: str = "unspecified",
) -> dict[str, Any]:
    """Build an all-entity shadow expansion bundle bound to V11 and raw rows.

    The V11 member-index partition must cover every raw result exactly once,
    and every indexed spectrum id must match.  This prevents a caller from
    silently pairing an entity artifact with another acquisition.  Only the
    measured scan projection created by :meth:`retrieve_entity` reaches the
    spectral retriever; legacy matches and names remain outside its inputs.
    """

    metadata = v11_payload.get("metadata")
    entities = v11_payload.get("entities")
    if not isinstance(metadata, Mapping):
        raise ValueError("V11 metadata is required for candidate recall")
    if not isinstance(entities, Sequence) or isinstance(entities, (str, bytes)):
        raise ValueError("V11 entities are required for candidate recall")
    if len(source_rows) != int(metadata.get("source_row_count") or -1):
        raise ValueError("Raw result row count does not match the V11 artifact")
    source_hash = _text(metadata.get("source_semantic_sha256"))
    if len(source_hash) != 64:
        raise ValueError("V11 source semantic hash is missing")
    method_id = _text(detection_method_id) or "unspecified"
    seen_indices: set[int] = set()
    output_entities: list[dict[str, Any]] = []
    for entity in entities:
        if not isinstance(entity, Mapping):
            raise ValueError("Every V11 entity must be a JSON object")
        entity_id = _text(entity.get("entity_id"))
        if not entity_id:
            raise ValueError("Every V11 entity needs an entity_id")
        member_indices = [
            int(index) for index in (entity.get("member_result_indices") or [])
        ]
        member_ids = [
            _text(value) for value in (entity.get("member_spectrum_ids") or [])
        ]
        if len(member_indices) != len(member_ids):
            raise ValueError(f"V11 member index/id length mismatch for {entity_id}")
        scans = []
        for index, expected_id in zip(member_indices, member_ids):
            if index < 0 or index >= len(source_rows):
                raise ValueError(f"V11 member index is out of range for {entity_id}")
            if index in seen_indices:
                raise ValueError(f"Raw result index {index} belongs to multiple V11 entities")
            row = source_rows[index]
            actual_id = _text(row.get("spectrum_id"))
            if actual_id != expected_id:
                raise ValueError(
                    f"Raw spectrum id mismatch at index {index}: "
                    f"{actual_id!r} != {expected_id!r}"
                )
            seen_indices.add(index)
            scans.append(row)
        output_entities.append(
            retriever.retrieve_entity(
                entity_id,
                scans,
                existing_candidates=[
                    candidate
                    for candidate in (entity.get("candidates") or [])
                    if isinstance(candidate, Mapping)
                ],
                detection_method_id=method_id,
            )
        )
    if seen_indices != set(range(len(source_rows))):
        raise ValueError("V11 entity membership does not partition every raw result row")
    return {
        "metadata": {
            "schema_version": EVIRECALL_SCHEMA_VERSION,
            "module_id": MODULE_ID,
            "algorithm_id": retriever.policy["algorithm_id"],
            "algorithm_name": retriever.policy["display_name"],
            "deployment_state": "shadow",
            "task_id": metadata.get("task_id"),
            "detection_method_id": method_id,
            "base_algorithm_id": metadata.get("algorithm_id"),
            "base_policy_sha256": metadata.get("policy_sha256"),
            "base_source_semantic_sha256": source_hash,
            "base_artifact_sha256": _sha256(v11_payload),
            "library_snapshot_id": retriever.library_snapshot_id,
            "policy_sha256": retriever.policy_sha256,
            "source_row_count": len(source_rows),
            "entity_count": len(output_entities),
            "entities_with_expansion": sum(
                bool(entity.get("expanded_candidate_count"))
                for entity in output_entities
            ),
            "expanded_candidate_count": sum(
                int(entity.get("expanded_candidate_count") or 0)
                for entity in output_entities
            ),
            "entities_with_qualified_source_paired_switch": sum(
                (entity.get("source_paired_consensus") or {}).get("status")
                == "qualified_switch"
                for entity in output_entities
            ),
            "source_membership_partition_validated": True,
            "score_is_probability": False,
            "production_effect": False,
        },
        "policy": retriever.policy,
        "invariants": {
            "shadow_only": True,
            "ground_truth_used_at_inference": False,
            "candidate_names_affect_retrieval_or_scoring": False,
            "expanded_candidates_enter_shadow_candidate_union": True,
            "expanded_candidates_enter_residual_model_scoring": False,
            "source_paired_consensus_ground_truth_used_at_inference": False,
            "provenance_strata_are_independent_sources": False,
            "mirrored_reference_content_profiles_create_extra_votes": False,
            "same_connectivity_different_detection_methods_preexcluded": False,
            "existing_candidate_exclusion_scope": "within_one_entity_and_detection_method_only",
            "technical_repeats_increase_condition_count": False,
            "blind_window_widening": False,
            "blind_adduct_enumeration": False,
            "score_is_probability": False,
            "production_effect": False,
        },
        "entities": output_entities,
    }


class SQLiteCandidateRecallRetriever:
    """Read-only, bounded spectral candidate retriever for one library snapshot."""

    def __init__(
        self,
        db_path: Path,
        *,
        library_snapshot_id: str,
        policy: Mapping[str, Any] | None = None,
    ) -> None:
        self.db_path = Path(db_path).resolve()
        if not self.db_path.is_file():
            raise FileNotFoundError(self.db_path)
        self.library_snapshot_id = _text(library_snapshot_id)
        if not self.library_snapshot_id:
            raise ValueError("A frozen library_snapshot_id is required")
        self.policy = (
            load_evirecall_policy()
            if policy is None
            else _validate_evirecall_policy(policy)
        )
        self.policy_sha256 = evirecall_policy_sha256(self.policy)
        self._window_cache: OrderedDict[
            tuple[float, float, str | None],
            tuple[list[dict[str, Any]], dict[str, Any]],
        ] = OrderedDict()
        self._sqlite_query_strategy = self._resolve_sqlite_query_strategy()

    @staticmethod
    def _quoted_identifier(value: str) -> str:
        return '"' + str(value).replace('"', '""') + '"'

    @classmethod
    def _table_indexes(
        cls,
        connection: sqlite3.Connection,
        table: str,
    ) -> dict[str, tuple[str, ...]]:
        indexes: dict[str, tuple[str, ...]] = {}
        quoted_table = cls._quoted_identifier(table)
        for row in connection.execute(f"PRAGMA index_list({quoted_table})"):
            index_name = str(row[1])
            quoted_index = cls._quoted_identifier(index_name)
            indexes[index_name] = tuple(
                str(info[2])
                for info in connection.execute(
                    f"PRAGMA index_info({quoted_index})"
                )
                if info[2] is not None
            )
        return indexes

    def _resolve_sqlite_query_strategy(self) -> dict[str, Any]:
        """Choose an exact bounded join plan without mutating the library.

        Some frozen production snapshots predate the explicit
        ``Mapping(spectrum_id)`` index.  SQLite otherwise starts the candidate
        recall query from Mapping and scans millions of rows for every narrow
        precursor window.  Those snapshots were imported with the verified
        invariant ``Mapping.id == Mapping.spectrum_id``.  When the invariant
        and bidirectional row coverage both hold, joining through the Mapping
        primary key is exactly equivalent and keeps the database read-only.
        A non-aligned snapshot must provide a real spectrum_id index; it fails
        closed instead of falling back to an unbounded full-table scan.
        """

        connection = self._connect()
        try:
            spectra_indexes = self._table_indexes(connection, "Spectra")
            mapping_indexes = self._table_indexes(connection, "Mapping")
            precursor_index = next(
                (
                    name
                    for name, columns in spectra_indexes.items()
                    if columns and columns[0] == "precursor_mz"
                ),
                None,
            )
            mapping_spectrum_index = next(
                (
                    name
                    for name, columns in mapping_indexes.items()
                    if columns and columns[0] == "spectrum_id"
                ),
                None,
            )

            if mapping_spectrum_index:
                mapping_join = "m.spectrum_id=s.id"
                mapping_strategy = "indexed_mapping_spectrum_id"
            else:
                mapping_misaligned = connection.execute(
                    "SELECT 1 FROM Mapping AS m "
                    "WHERE m.spectrum_id IS NULL OR m.id != m.spectrum_id "
                    "LIMIT 1"
                ).fetchone()
                spectrum_without_mapping = connection.execute(
                    "SELECT 1 FROM Spectra AS s "
                    "LEFT JOIN Mapping AS m ON m.id=s.id "
                    "WHERE m.id IS NULL LIMIT 1"
                ).fetchone()
                mapping_without_spectrum = connection.execute(
                    "SELECT 1 FROM Mapping AS m "
                    "LEFT JOIN Spectra AS s ON s.id=m.id "
                    "WHERE s.id IS NULL LIMIT 1"
                ).fetchone()
                if (
                    mapping_misaligned is not None
                    or spectrum_without_mapping is not None
                    or mapping_without_spectrum is not None
                ):
                    raise RuntimeError(
                        "Candidate recall requires an index whose leading "
                        "column is Mapping.spectrum_id. The read-only library "
                        "does not provide that index and its Mapping primary "
                        "key is not an exact verified spectrum_id surrogate."
                    )
                mapping_join = "m.id=s.id"
                mapping_strategy = "verified_mapping_primary_key_alignment"

            spectra_hint = (
                " INDEXED BY " + self._quoted_identifier(precursor_index)
                if precursor_index
                else ""
            )
            from_clause = (
                f"Spectra AS s{spectra_hint} "
                f"JOIN Mapping AS m ON {mapping_join} "
                "JOIN Compounds AS c ON c.id=m.compound_id"
            )
            return {
                "from_clause": from_clause,
                "precursor_index": precursor_index,
                "mapping_spectrum_index": mapping_spectrum_index,
                "mapping_join_strategy": mapping_strategy,
                "library_read_only": True,
            }
        finally:
            connection.close()

    def _cache_window(
        self,
        key: tuple[float, float, str | None],
        records: list[dict[str, Any]],
        audit: dict[str, Any],
    ) -> None:
        self._window_cache[key] = (records, audit)
        self._window_cache.move_to_end(key)
        maximum = int(
            self.policy["library_gate"]["maximum_precursor_windows_in_memory"]
        )
        while len(self._window_cache) > maximum:
            self._window_cache.popitem(last=False)

    def _connect(self) -> sqlite3.Connection:
        uri = f"file:{quote(str(self.db_path), safe='/')}?mode=ro"
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only=ON")
        return connection

    def _fetch_window(
        self,
        precursor_mz: float,
        query_polarity: str | None,
    ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
        mass = self.policy["mass_retrieval"]
        tolerance = max(
            float(mass["direct_precursor_da"]),
            precursor_mz * float(mass["direct_precursor_ppm"]) / 1e6,
        )
        low, high = precursor_mz - tolerance, precursor_mz + tolerance
        cache_key = (round(low, 8), round(high, 8), query_polarity)
        cached = self._window_cache.get(cache_key)
        if cached is not None:
            self._window_cache.move_to_end(cache_key)
            records, audit = cached
            return records, {**audit, "cache_hit": True}

        tiers = tuple(self.policy["library_gate"]["searchable_source_tiers"])
        placeholders = ",".join("?" for _ in tiers)
        where = (
            "s.precursor_mz BETWEEN ? AND ? "
            f"AND lower(coalesce(s.source_tier,'')) IN ({placeholders})"
        )
        maximum_rows = int(
            self.policy["library_gate"]["maximum_sql_rows_per_precursor_window"]
        )
        connection = self._connect()
        try:
            from_clause = str(self._sqlite_query_strategy["from_clause"])
            spectra_columns = {
                str(row[1])
                for row in connection.execute("PRAGMA table_info(Spectra)")
            }
            row_count = int(connection.execute(
                "SELECT count(*) FROM " + from_clause + " WHERE " + where,
                (low, high, *tiers),
            ).fetchone()[0])
            if row_count > maximum_rows:
                audit = {
                    "status": "dense_window_budget_exceeded_fail_closed",
                    "precursor_mz": precursor_mz,
                    "effective_tolerance_da": tolerance,
                    "sql_rows": row_count,
                    "maximum_sql_rows": maximum_rows,
                    "sqlite_query_strategy": {
                        key: value
                        for key, value in self._sqlite_query_strategy.items()
                        if key != "from_clause"
                    },
                    "cache_hit": False,
                }
                self._cache_window(cache_key, [], audit)
                return [], audit
            condition_fields = (
                "collision_energy",
                "collision_energy_raw",
                "collision_energy_unit",
                "collision_energy_mode",
                "collision_energy_min",
                "collision_energy_max",
                "instrument_type",
                "instrument",
                "source_accession",
                "activation_method",
                "instrument_domain_id",
                "collision_cell_region",
            )
            condition_select = ",".join(
                (
                    f"s.{field} AS {field}"
                    if field in spectra_columns
                    else f"NULL AS {field}"
                )
                for field in condition_fields
            )
            rows = connection.execute(
                "SELECT s.precursor_mz,s.peaks_json,s.spectrum_id,s.adduct,"
                "s.charge,s.ion_mode,s.source_tier,s.dataset_source,"
                "s.quality_flags_json," + condition_select + ","
                "c.compound_name,c.molecular_formula,"
                "c.smiles,c.inchikey,c.exact_mass "
                "FROM " + from_clause + " WHERE " + where + " "
                "ORDER BY s.precursor_mz,s.spectrum_id,c.id",
                (low, high, *tiers),
            ).fetchall()
        finally:
            connection.close()

        records: list[dict[str, Any]] = []
        excluded_polarity = 0
        excluded_identity = 0
        excluded_quality = 0
        for row in rows:
            record = dict(row)
            tier = _text(record.get("source_tier")).lower()
            if not source_tier_is_searchable(tier):
                continue
            library_polarity = classify_spectrum_polarity(
                record.get("adduct"), record.get("charge"), record.get("ion_mode")
            )
            if query_polarity and library_polarity not in {None, query_polarity}:
                excluded_polarity += 1
                continue
            try:
                flags = json.loads(_text(record.get("quality_flags_json")) or "[]")
            except (TypeError, ValueError):
                flags = []
            if (
                "excluded_from_ordinary_spectral_search" in flags
                or "flat_fragment_intensities_reverse_only" in flags
            ):
                excluded_quality += 1
                continue
            connectivity, identity_source = _stable_connectivity(record)
            if not connectivity:
                excluded_identity += 1
                continue
            try:
                peaks = json.loads(_text(record.get("peaks_json")) or "[]")
            except (TypeError, ValueError):
                peaks = []
            record.update({
                "connectivity_identity": connectivity,
                "identity_source": identity_source,
                "library_polarity": library_polarity,
                "quality_flags": flags,
                "peaks": peaks,
                "processed_views": {},
            })
            records.append(record)
        audit = {
            "status": "available",
            "precursor_mz": precursor_mz,
            "effective_tolerance_da": tolerance,
            "sql_rows": row_count,
            "eligible_structure_rows": len(records),
            "excluded_wrong_polarity": excluded_polarity,
            "excluded_unstable_identity": excluded_identity,
            "excluded_quality": excluded_quality,
            "blind_window_widening": False,
            "blind_adduct_enumeration": False,
            "sqlite_query_strategy": {
                key: value
                for key, value in self._sqlite_query_strategy.items()
                if key != "from_clause"
            },
            "cache_hit": False,
        }
        self._cache_window(cache_key, records, audit)
        return records, audit

    def retrieve_scan(self, query: Mapping[str, Any]) -> dict[str, Any]:
        leaked = _forbidden_request_paths(query)
        if leaked:
            raise ValueError(
                "Leakage-prone fields are forbidden in an EviRecall request: "
                + ", ".join(sorted(set(leaked)))
            )
        precursor = _finite(query.get("precursor_mz"))
        if precursor is None or precursor <= 0:
            raise ValueError("A positive measured precursor_mz is required")
        query_polarity = _polarity(query.get("polarity") or query.get("ion_mode"))
        records, window_audit = self._fetch_window(precursor, query_polarity)
        if not records:
            return {
                "status": window_audit["status"],
                "query_id": _text(query.get("spectrum_id") or query.get("query_id")),
                "condition_id": _condition_id(query),
                "rankings": {},
                "candidates": {},
                "window_audit": window_audit,
                "score_is_probability": False,
            }

        raw_peaks = query.get("peaks") or query.get("raw_peaks") or query.get("query_peaks") or []
        minimum_peaks = int(
            self.policy["candidate_pool"]["minimum_library_fragment_peaks"]
        )
        candidates: dict[str, dict[str, Any]] = {}
        rankings: dict[str, list[str]] = {}
        for view in self.policy["spectral_views"]:
            view_name = _text(view["name"])
            config = resolve_peak_processing_config(view["peak_processing_preset"])
            query_view, query_audit = process_peaks(
                raw_peaks,
                precursor,
                config=config,
                precursor_charge=query.get("charge"),
            )
            if len(query_view) < minimum_peaks:
                rankings[view_name] = []
                continue
            usable_records: list[dict[str, Any]] = []
            library_views: list[list[list[float]]] = []
            library_precursors: list[float] = []
            for record in records:
                processed = record["processed_views"].get(view_name)
                if processed is None:
                    processed, _ = process_peaks(
                        record.get("peaks") or [],
                        record.get("precursor_mz"),
                        config=config,
                        precursor_charge=record.get("charge"),
                    )
                    record["processed_views"][view_name] = processed
                if len(processed) < minimum_peaks:
                    continue
                usable_records.append(record)
                library_views.append(processed)
                library_precursors.append(float(record["precursor_mz"]))
            score_function = (
                batch_neutral_loss_cosine_scores
                if view["similarity"] == "matchms_neutral_losses_cosine"
                else batch_modified_cosine_scores
            )
            scores = score_function(
                query_view,
                library_views,
                precursor,
                library_precursors,
                mz_tol=float(view["fragment_tolerance_da"]),
            )
            best: dict[str, tuple[tuple[Any, ...], dict[str, Any], float]] = {}
            best_by_source: dict[
                tuple[str, str], tuple[tuple[Any, ...], dict[str, Any], float]
            ] = {}
            source_record_counts: dict[tuple[str, str], int] = defaultdict(int)
            for record, raw_score in zip(usable_records, scores):
                score = float(raw_score)
                identity = record["connectivity_identity"]
                mass_error = abs(float(record["precursor_mz"]) - precursor)
                order = (-score, mass_error, _text(record.get("spectrum_id")))
                previous = best.get(identity)
                if previous is None or order < previous[0]:
                    best[identity] = (order, record, score)
                stratum = _provenance_stratum(record.get("dataset_source"))
                if stratum:
                    source_key = (identity, stratum)
                    source_record_counts[source_key] += 1
                    previous_source = best_by_source.get(source_key)
                    if previous_source is None or order < previous_source[0]:
                        best_by_source[source_key] = (order, record, score)
            ordered = sorted(
                best,
                key=lambda identity: (
                    best[identity][0],
                    identity,
                ),
            )[: int(self.policy["candidate_pool"]["maximum_identities_per_scan_view"])]
            rankings[view_name] = ordered
            for rank, identity in enumerate(ordered, 1):
                _, record, score = best[identity]
                candidate = candidates.setdefault(identity, {
                    "connectivity_identity": identity,
                    "identity_source": record["identity_source"],
                    "candidate": {
                        "compound_name": _text(record.get("compound_name")),
                        "formula": _text(record.get("molecular_formula")),
                        "smiles": _text(record.get("smiles")),
                        "inchikey": normalize_inchikey(record.get("inchikey")),
                    },
                    "view_scores": {},
                    "view_ranks": {},
                    "reference_records": {},
                    "source_view_scores": {},
                    "source_reference_record_counts": {},
                    "source_reference_content_fingerprints": {},
                    "source_reference_spectrum_ids": {},
                })
                candidate["view_scores"][view_name] = round(score, 8)
                candidate["view_ranks"][view_name] = rank
                candidate["reference_records"][view_name] = {
                    "spectrum_id": _text(record.get("spectrum_id")),
                    "dataset_source": _text(record.get("dataset_source")),
                    "source_tier": _text(record.get("source_tier")),
                    "library_precursor_mz": float(record["precursor_mz"]),
                    "mass_error_da": round(
                        float(record["precursor_mz"]) - precursor, 8
                    ),
                    "library_polarity": record.get("library_polarity"),
                    "collision_energy": record.get("collision_energy"),
                    "collision_energy_raw": _text(
                        record.get("collision_energy_raw")
                    ),
                    "collision_energy_unit": _text(
                        record.get("collision_energy_unit")
                    ),
                    "collision_energy_mode": _text(
                        record.get("collision_energy_mode")
                    ),
                    "collision_energy_min": record.get("collision_energy_min"),
                    "collision_energy_max": record.get("collision_energy_max"),
                    "instrument_type": _text(record.get("instrument_type")),
                    "instrument": _text(record.get("instrument")),
                    "source_ms_level": _text(record.get("source_ms_level")),
                    "database_spectrum_id": record.get("database_spectrum_id"),
                    "source_accession": _text(record.get("source_accession")),
                    "content_fingerprint": _reference_content_fingerprint(record),
                }
                source_rows = sorted(
                    (
                        stratum,
                        source_record[2],
                        source_record_counts[(identity, stratum)],
                        _reference_content_fingerprint(source_record[1]),
                        _text(source_record[1].get("spectrum_id")),
                        source_record[1],
                    )
                    for (source_identity, stratum), source_record in best_by_source.items()
                    if source_identity == identity
                )
                candidate["source_view_scores"][view_name] = {
                    stratum: round(score_value, 8)
                    for stratum, score_value, _, _, _, _ in source_rows
                }
                candidate["source_reference_record_counts"][view_name] = {
                    stratum: count
                    for stratum, _, count, _, _, _ in source_rows
                }
                candidate["source_reference_content_fingerprints"][view_name] = {
                    stratum: fingerprint
                    for stratum, _, _, fingerprint, _, _ in source_rows
                    if fingerprint
                }
                candidate["source_reference_spectrum_ids"][view_name] = {
                    stratum: spectrum_id
                    for stratum, _, _, _, spectrum_id, _ in source_rows
                    if spectrum_id
                }
                candidate.setdefault("source_reference_records", {})[
                    view_name
                ] = {
                    stratum: {
                        "spectrum_id": spectrum_id,
                        "source_accession": _text(source_record.get("source_accession")),
                        "dataset_source": _text(source_record.get("dataset_source")),
                        "source_tier": _text(source_record.get("source_tier")),
                        "score": round(score_value, 8),
                        "reference_record_count_in_source": count,
                        "content_fingerprint": fingerprint,
                        "collision_energy": source_record.get("collision_energy"),
                        "collision_energy_raw": _text(
                            source_record.get("collision_energy_raw")
                        ),
                        "collision_energy_unit": _text(
                            source_record.get("collision_energy_unit")
                        ),
                        "collision_energy_mode": _text(
                            source_record.get("collision_energy_mode")
                        ),
                        "collision_energy_min": source_record.get(
                            "collision_energy_min"
                        ),
                        "collision_energy_max": source_record.get(
                            "collision_energy_max"
                        ),
                        "instrument_type": _text(
                            source_record.get("instrument_type")
                        ),
                        "instrument": _text(source_record.get("instrument")),
                        "source_ms_level": _text(source_record.get("source_ms_level")),
                        "database_spectrum_id": source_record.get("database_spectrum_id"),
                        "activation_method": _text(
                            source_record.get("activation_method")
                        ),
                        "instrument_domain_id": _text(
                            source_record.get("instrument_domain_id")
                        ),
                        "collision_cell_region": _text(
                            source_record.get("collision_cell_region")
                        ),
                    }
                    for (
                        stratum,
                        score_value,
                        count,
                        fingerprint,
                        spectrum_id,
                        source_record,
                    ) in source_rows
                }
                candidate.setdefault("query_peak_audits", {})[view_name] = query_audit
        return {
            "status": "available" if any(rankings.values()) else "query_low_information",
            "query_id": _text(query.get("spectrum_id") or query.get("query_id")),
            "condition_id": _condition_id(query),
            "query_precursor_mz": precursor,
            "query_polarity": query_polarity,
            "query_condition": {
                "collision_energy": query.get("collision_energy"),
                "collision_energy_unit": _text(
                    query.get("collision_energy_unit")
                ),
                "collision_energy_mode": _text(
                    query.get("collision_energy_mode")
                ),
                "instrument_type": _text(query.get("instrument_type")),
                "instrument": _text(query.get("instrument")),
                "activation_method": _text(
                    query.get("activation_method") or query.get("activation")
                ),
                "instrument_domain_id": _text(query.get("instrument_domain_id")),
                "collision_cell_region": _text(query.get("collision_cell_region")),
            },
            "rankings": rankings,
            "candidates": candidates,
            "window_audit": window_audit,
            "score_is_probability": False,
        }

    def retrieve_entity(
        self,
        entity_id: str,
        scans: Sequence[Mapping[str, Any]],
        *,
        existing_candidates: Sequence[Mapping[str, Any]] = (),
        detection_method_id: str = "unspecified",
    ) -> dict[str, Any]:
        scan_results = []
        for scan in scans:
            projection = {
                "spectrum_id": scan.get("spectrum_id"),
                "precursor_mz": scan.get("precursor_mz"),
                "polarity": scan.get("ion_mode") or scan.get("polarity"),
                "charge": scan.get("charge") or scan.get("precursor_charge"),
                "collision_energy": scan.get("collision_energy"),
                "collision_energy_unit": scan.get("collision_energy_unit"),
                "collision_energy_mode": scan.get("collision_energy_mode"),
                "instrument_type": scan.get("instrument_type"),
                "instrument": scan.get("instrument"),
                "activation_method": scan.get("activation_method"),
                "instrument_domain_id": scan.get("instrument_domain_id"),
                "collision_cell_region": scan.get("collision_cell_region"),
                "peaks": scan.get("raw_peaks") or scan.get("query_peaks") or scan.get("peaks") or [],
            }
            scan_results.append(self.retrieve_scan(projection))

        rrf_constant = float(self.policy["fusion"]["rrf_constant"])
        contributions: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        candidate_records: dict[str, dict[str, Any]] = {}
        best_display_order: dict[str, tuple[float, str]] = {}
        for result in scan_results:
            condition = result["condition_id"]
            for view_name, ranking in result.get("rankings", {}).items():
                for rank, identity in enumerate(ranking, 1):
                    contributions[(condition, view_name, identity)].append(
                        1.0 / (rrf_constant + rank)
                    )
                    record = result["candidates"][identity]
                    score = float(record["view_scores"].get(view_name) or 0.0)
                    order = (-score, _text((record.get("candidate") or {}).get("inchikey")))
                    if identity not in best_display_order or order < best_display_order[identity]:
                        best_display_order[identity] = order
                        candidate_records[identity] = record

        fused: dict[str, dict[str, Any]] = {}
        for (condition, view_name, identity), values in contributions.items():
            record = fused.setdefault(identity, {
                "connectivity_identity": identity,
                "condition_view_contributions": [],
                "conditions": set(),
                "views": set(),
                "rrf_score": 0.0,
            })
            collapsed = float(statistics.median(values))
            record["condition_view_contributions"].append({
                "condition_id": condition,
                "view": view_name,
                "technical_repeat_count": len(values),
                "collapsed_rrf_contribution": round(collapsed, 10),
            })
            record["conditions"].add(condition)
            record["views"].add(view_name)
            record["rrf_score"] += collapsed

        existing = _existing_connectivities(existing_candidates)
        ranked = sorted(
            (record for identity, record in fused.items() if identity not in existing),
            key=lambda record: (
                -float(record["rrf_score"]),
                -len(record["conditions"]),
                -len(record["views"]),
                record["connectivity_identity"],
            ),
        )
        target = int(self.policy["candidate_pool"]["target_total_connectivities"])
        maximum_additions = int(
            self.policy["candidate_pool"]["maximum_additions_per_entity"]
        )
        slots = max(0, min(maximum_additions, target - len(existing)))
        additions = []
        for rank, fused_record in enumerate(ranked[:slots], 1):
            identity = fused_record["connectivity_identity"]
            source = candidate_records[identity]
            additions.append({
                "identity_key": f"std2d:{identity}",
                "connectivity_identity": identity,
                "identity_scope": "two_dimensional_connectivity",
                "identity_source": source.get("identity_source"),
                "candidate": source.get("candidate") or {},
                "candidate_origin": "evirecall_expansion",
                "enters_shadow_candidate_union": True,
                "enters_residual_model_scoring": False,
                "source_paired_selection_qualified": False,
                "retrieval_rank": rank,
                "retrieval_rrf_utility": round(float(fused_record["rrf_score"]), 10),
                "score_is_probability": False,
                "production_effect": False,
                "distinct_condition_support": len(fused_record["conditions"]),
                "spectral_view_support": sorted(fused_record["views"]),
                "condition_view_contributions": sorted(
                    fused_record["condition_view_contributions"],
                    key=lambda item: (item["condition_id"], item["view"]),
                ),
                "best_reference_records": source.get("reference_records") or {},
            })
        source_paired_consensus = _source_paired_consensus(
            existing_candidates,
            scan_results,
            self.policy,
            expanded_candidate_connectivities=[
                item["connectivity_identity"] for item in additions
            ],
        )
        selected_connectivity = _text(
            source_paired_consensus.get("selected_connectivity_identity")
        )
        qualified_switch = source_paired_consensus.get("status") == "qualified_switch"
        for candidate in additions:
            candidate["source_paired_selection_qualified"] = bool(
                qualified_switch
                and candidate["connectivity_identity"] == selected_connectivity
            )
        return {
            "schema_version": EVIRECALL_SCHEMA_VERSION,
            "module_id": MODULE_ID,
            "algorithm_id": self.policy["algorithm_id"],
            "algorithm_name": self.policy["display_name"],
            "deployment_state": "shadow",
            "entity_id": _text(entity_id),
            "detection_method_id": _text(detection_method_id) or "unspecified",
            "library_snapshot_id": self.library_snapshot_id,
            "policy_sha256": self.policy_sha256,
            "existing_connectivity_count": len(existing),
            "target_total_connectivities": target,
            "expanded_candidate_count": len(additions),
            "expanded_candidates": additions,
            "source_paired_consensus": source_paired_consensus,
            "scan_count": len(scans),
            "distinct_condition_count": len({result["condition_id"] for result in scan_results}),
            "technical_repeats_increase_condition_count": False,
            "mirrored_reference_content_profiles_create_extra_votes": False,
            "candidate_names_affect_retrieval_or_scoring": False,
            "expanded_candidates_enter_shadow_candidate_union": True,
            "expanded_candidates_enter_residual_model_scoring": False,
            "blind_window_widening": False,
            "blind_adduct_enumeration": False,
            "score_is_probability": False,
            "production_effect": False,
            "scan_audit": [
                {
                    "query_id": result["query_id"],
                    "condition_id": result["condition_id"],
                    "status": result["status"],
                    "window_audit": result["window_audit"],
                    "ranked_identity_counts": {
                        view: len(ranking)
                        for view, ranking in result.get("rankings", {}).items()
                    },
                }
                for result in scan_results
            ],
        }


__all__ = [
    "EVIRECALL_SCHEMA_VERSION",
    "FORBIDDEN_REQUEST_FIELDS",
    "MODULE_ID",
    "SQLiteCandidateRecallRetriever",
    "build_candidate_recall_bundle",
    "default_evirecall_policy_path",
    "evirecall_policy_sha256",
    "load_evirecall_policy",
]
