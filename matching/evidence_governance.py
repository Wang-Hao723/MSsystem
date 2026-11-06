"""EviGate-MS: applicability-aware evidence passports and replicate governance.

This module is deliberately *downstream* of spectral retrieval.  It never changes
the production Top-1 identity, candidate order, or any stored similarity score.
Instead it answers a different question: which scientific claim is justified by
the evidence that is actually applicable, sufficiently independent, and free of
obvious circular reuse?

Two units are kept separate:

* :func:`build_evidence_passport` audits one reported spectrum/candidate.
* :func:`build_replicate_consensus` collapses stepped-energy scans inside one
  precursor event, then seeks agreement across separately submitted injections.

The second step prevents three collision-energy scans from being counted as
three independent confirmations.  Candidate aggregation is a shadow/read-only
view and is exposed as an ambiguity set; it is not a replacement ranker and it
does not claim conformal coverage until an external calibration pack exists.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from typing import Any, Iterable

from matching.library_quality import compound_identity_key


SCHEMA_VERSION = "1.0.0"
POLICY_ID = "evigate-ms-locked-v1-2026-08-02"

# Locked before evaluation.  These are evidence-release thresholds, not search
# thresholds; changing them requires a new policy identifier and a new ablation.
MIN_SPECTRAL_COSINE = 0.80
MIN_FINAL_SCORE = 0.50
MIN_SCOPED_EXPERT_SCORE = 0.40
MAX_SUPPORTING_MASS_PPM = 50.0
COHORT_MZ_PPM = 15.0
COHORT_RT_WINDOW_MIN = 0.15
MIN_INDEPENDENT_INJECTIONS = 2
MIN_IDENTITY_AGREEMENT = 2.0 / 3.0
RRF_OFFSET = 20.0
MAX_AMBIGUITY_SET = 5

_EXPERIMENTAL_SOURCE_TIERS = {
    "experimental",
    "curated_experimental",
    "local_experimental",
}
_SEVERE_QUALITY_FLAGS = {
    "precursor_adduct_conflict_manual_review_required",
    "invalid_structure",
    "formula_structure_conflict",
}


def _number(value: Any) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _median(values: Iterable[Any]) -> float | None:
    numbers = [number for number in (_number(value) for value in values)
               if number is not None]
    return statistics.median(numbers) if numbers else None


def _identity(candidate: dict) -> str:
    try:
        return compound_identity_key(candidate)
    except Exception:
        return ""


def _connectivity_identity(candidate: dict) -> str:
    """Connectivity-level key for replicate voting.

    Ordinary achiral LC-MS/MS cannot use a differing stereochemical InChIKey
    block as evidence that two injections disagree.  The first 14 characters
    preserve constitutional connectivity while keeping stereochemical/full
    identity variants visible in the evidence passport.
    """
    inchikey = str(candidate.get("inchikey") or "").strip().upper()
    first_block = inchikey.split("-")[0]
    if len(first_block) == 14 and first_block.isalpha():
        return f"ik14::{first_block}"
    return _identity(candidate)


def _candidate_summary(candidate: dict, identity: str | None = None) -> dict:
    return {
        "identity_key": identity if identity is not None else _identity(candidate),
        "compound_name": str(candidate.get("compound_name") or ""),
        "formula": str(candidate.get("formula") or ""),
        "inchikey": str(candidate.get("inchikey") or ""),
        "cosine": _number(candidate.get("cosine")),
        "final_score": _number(candidate.get("final_score")),
        "expert_score": _number(candidate.get("expert_score")),
        "source_tier": str(candidate.get("source_tier") or "unknown"),
    }


def _node(
    key: str,
    label: str,
    *,
    status: str,
    direction: str,
    independence_group: str,
    independent_measurement: bool,
    used_for: list[str],
    reason: str,
    value: Any = None,
    applicability: str = "applicable",
    provenance: dict | None = None,
    warnings: list[str] | None = None,
) -> dict:
    return {
        "key": key,
        "label": label,
        "status": status,
        "direction": direction,
        "applicability": applicability,
        "reason": reason,
        "value": value,
        "provenance": provenance or {},
        "independence_group": independence_group,
        "independent_measurement": independent_measurement,
        "used_for": used_for,
        "warnings": warnings or [],
    }


def _scoped_expert(candidate: dict) -> dict:
    details = candidate.get("expert_details") or {}
    dimensions = details.get("dimension_assessments") or {}
    scored = []
    rejected = []
    for key, assessment in dimensions.items():
        if not isinstance(assessment, dict):
            continue
        status = str(assessment.get("status") or "")
        if status == "scored":
            scored.append(str(key))
        elif status in {"not_applicable", "unavailable"}:
            rejected.append(str(key))
    score = _number(candidate.get("expert_score")) or 0.0
    applicable = bool(scored)
    passes = applicable and score >= MIN_SCOPED_EXPERT_SCORE
    return {
        "score": round(score, 4),
        "applicable": applicable,
        "passes_locked_gate": passes,
        "scored_dimensions": sorted(scored),
        "non_scored_dimensions": sorted(rejected),
        "threshold": MIN_SCOPED_EXPERT_SCORE,
        "selected_source": str(
            ((details.get("score_fusion") or {}).get("selected_source"))
            or candidate.get("expert_score_source")
            or "none"
        ),
    }


def _candidate_release_audit(candidate: dict, row: dict | None = None) -> dict:
    """Apply the locked, scope-aware release gate to one candidate.

    Passing this gate means "eligible for replicate governance", not confirmed
    identity.  Independent-injection agreement or same-method RT evidence is
    still required before EviGate releases a connectivity-level candidate.
    """
    row = row or {}
    cosine = _number(candidate.get("cosine")) or 0.0
    final_score = _number(candidate.get("final_score")) or 0.0
    ppm = _number(candidate.get("mass_ppm"))
    source_tier = str(candidate.get("source_tier") or "unknown").lower()
    smiles = str(candidate.get("smiles") or "").strip().lower()
    quality_flags = {str(flag) for flag in candidate.get("quality_flags") or []}
    expert = _scoped_expert(candidate)
    blockers: list[str] = []
    if cosine < MIN_SPECTRAL_COSINE:
        blockers.append("cosine_below_locked_release_threshold")
    if final_score < MIN_FINAL_SCORE:
        blockers.append("final_score_below_locked_release_threshold")
    if not expert["passes_locked_gate"]:
        blockers.append(
            "expert_not_applicable"
            if not expert["applicable"]
            else "scoped_expert_support_below_threshold"
        )
    if candidate.get("has_structure") is False or smiles in {
        "", "nan", "none", "null", "n/a"
    }:
        blockers.append("no_validated_structure")
    if source_tier not in _EXPERIMENTAL_SOURCE_TIERS:
        blockers.append("library_source_not_experimental")
    if candidate.get("mass_metadata_consistent") is False:
        blockers.append("precursor_formula_adduct_conflict")
    if ppm is not None and abs(ppm) > MAX_SUPPORTING_MASS_PPM:
        blockers.append("mass_error_above_locked_release_limit")
    if quality_flags & _SEVERE_QUALITY_FLAGS:
        blockers.append("severe_library_quality_flag")
    # A soft flag from an unspecified/carryover blank is context, not a hard
    # analyte veto.  Only the role-aware clean-view decision represents a
    # directional background contradiction at this layer.
    if bool(row.get("blank_clean_view_excluded")):
        blockers.append("blank_background_conflict")
    if bool(candidate.get("potential_contaminant")):
        blockers.append("potential_external_contaminant")
    return {
        "eligible": not blockers,
        "blockers": blockers,
        "expert": expert,
        "thresholds": {
            "cosine": MIN_SPECTRAL_COSINE,
            "final_score": MIN_FINAL_SCORE,
            "scoped_expert": MIN_SCOPED_EXPERT_SCORE,
            "max_abs_mass_ppm": MAX_SUPPORTING_MASS_PPM,
        },
    }


def _ambiguity_set(match: dict) -> list[dict]:
    """Build a transparent, non-calibrated set around the immutable Top-1."""
    candidates: list[dict] = []
    seen: set[str] = set()
    best_score = _number(match.get("rank_score"))
    if best_score is None:
        best_score = _number(match.get("final_score")) or 0.0
    best_formula = str(match.get("formula") or "")
    pool = [match, *(match.get("top_candidates") or [])]
    for rank, candidate in enumerate(pool, start=1):
        if not isinstance(candidate, dict):
            continue
        identity = _identity(candidate) or f"rank:{rank}"
        if identity in seen:
            continue
        seen.add(identity)
        score = _number(candidate.get("rank_score"))
        if score is None:
            score = _number(candidate.get("final_score")) or 0.0
        same_formula = bool(
            best_formula and str(candidate.get("formula") or "") == best_formula
        )
        competitive = rank == 1 or same_formula or best_score - score <= 0.05
        if not competitive:
            continue
        item = _candidate_summary(candidate, identity)
        item.update({
            "source_rank": rank,
            "rank_score": round(score, 4),
            "delta_from_top": round(max(0.0, best_score - score), 4),
            "same_formula_as_top1": same_formula,
        })
        candidates.append(item)
        if len(candidates) >= MAX_AMBIGUITY_SET:
            break
    return candidates


def build_evidence_passport(
    row: dict,
    *,
    sample_role: str = "sample",
    cohort_evidence: dict | None = None,
) -> dict | None:
    """Return a provenance-rich, scope-aware passport for one reported match."""
    match = row.get("match") if isinstance(row, dict) else None
    if not isinstance(match, dict):
        return None

    nodes: list[dict] = []
    cosine = _number(match.get("cosine")) or 0.0
    entropy = _number(match.get("entropy"))
    reverse = _number(match.get("reverse_match"))
    nodes.append(_node(
        "spectral_similarity",
        "Library MS/MS similarity",
        status="supporting" if cosine >= MIN_SPECTRAL_COSINE else "insufficient",
        direction="supports_candidate" if cosine >= MIN_SPECTRAL_COSINE else "weak",
        independence_group="query_ms2_vs_library_spectrum",
        independent_measurement=False,
        used_for=["ranking", "confidence"],
        reason=(
            "The correlated cosine/entropy/reverse metrics are counted as one "
            "spectral evidence family, not three confirmations."
        ),
        value={
            "cosine": round(cosine, 4),
            "entropy": round(entropy, 4) if entropy is not None else None,
            "reverse_match": round(reverse, 4) if reverse is not None else None,
        },
        provenance={
            "query_spectrum_id": str(row.get("spectrum_id") or ""),
            "library_spectrum_id": str(match.get("lib_spectrum_id") or ""),
            "library_accession": str(match.get("lib_source_accession") or ""),
        },
        warnings=(
            ["collision_energy_or_instrument_metadata_incomplete"]
            if not (
                row.get("collision_energy") is not None
                and match.get("lib_collision_energy") is not None
                and (row.get("instrument_type") or row.get("instrument"))
                and (match.get("lib_instrument_type") or match.get("lib_instrument"))
            ) else []
        ),
    ))

    ppm = _number(match.get("mass_ppm"))
    mass_conflict = match.get("mass_metadata_consistent") is False
    if mass_conflict or (ppm is not None and abs(ppm) > MAX_SUPPORTING_MASS_PPM):
        mass_status, mass_direction = "contradictory", "contradicts_candidate"
        mass_reason = "Formula/adduct metadata conflict or excessive precursor error."
    elif ppm is None:
        mass_status, mass_direction = "unavailable", "neutral"
        mass_reason = "No auditable precursor-mass error was stored."
    else:
        mass_status, mass_direction = "supporting", "supports_formula"
        mass_reason = "Precursor mass is compatible with the stored formula/adduct basis."
    nodes.append(_node(
        "precursor_mass",
        "Precursor mass and adduct",
        status=mass_status,
        direction=mass_direction,
        independence_group="precursor_measurement",
        independent_measurement=False,
        used_for=["retrieval", "confidence"],
        reason=mass_reason,
        value={
            "mass_ppm": ppm,
            "basis": str(match.get("mass_match_basis") or ""),
            "metadata_consistent": match.get("mass_metadata_consistent"),
        },
    ))

    source_tier = str(match.get("source_tier") or "unknown").lower()
    source_ok = source_tier in _EXPERIMENTAL_SOURCE_TIERS
    nodes.append(_node(
        "library_provenance",
        "Library source provenance",
        status="supporting" if source_ok else "scope_limited",
        direction="supports_library_candidate" if source_ok else "limits_claim",
        independence_group="library_record_provenance",
        independent_measurement=False,
        used_for=["confidence", "claim_scope"],
        reason=(
            "Experimental spectrum provenance supports a Level-2 library candidate."
            if source_ok else
            "Predicted, presence-only, or unknown provenance cannot release an "
            "experimental structural claim."
        ),
        value={
            "source_tier": source_tier,
            "dataset_source": str(match.get("dataset_source") or ""),
            "intensity_origin": str(match.get("intensity_origin") or "unknown"),
            "quality_flags": list(match.get("quality_flags") or []),
        },
    ))

    expert = _scoped_expert(match)
    if not expert["applicable"]:
        expert_status, expert_direction, expert_applicability = (
            "not_applicable", "neutral", "not_applicable"
        )
    elif expert["passes_locked_gate"]:
        expert_status, expert_direction, expert_applicability = (
            "supporting", "supports_candidate", "applicable"
        )
    else:
        expert_status, expert_direction, expert_applicability = (
            "insufficient", "weak", "applicable"
        )
    nodes.append(_node(
        "scoped_expert_reasoning",
        "Scope-checked expert reasoning",
        status=expert_status,
        direction=expert_direction,
        independence_group="query_ms2_derived_reasoning",
        independent_measurement=False,
        used_for=["ranking", "release_gate", "explanation"],
        reason=(
            "Only dimensions explicitly marked scored are eligible; unavailable "
            "or out-of-scope rules contribute no pseudo-evidence."
        ),
        value=expert,
        applicability=expert_applicability,
        warnings=[
            "expert_score_already_contributes_to_final_score_do_not_count_as_independent"
        ],
    ))

    feature_fragments = match.get("feature_fragment_evidence") or {}
    ff_status = str(feature_fragments.get("status") or "unavailable")
    ff_groups = list(
        feature_fragments.get("independent_group_representatives") or []
    )
    ff_log_lr = _number(
        feature_fragments.get("combined_log_likelihood_ratio")
    ) or 0.0
    if ff_status == "available" and ff_groups:
        ff_applicability = "applicable"
        if ff_log_lr > 0:
            ff_node_status, ff_direction = "supporting", "supports_class"
        elif ff_log_lr < 0:
            ff_node_status, ff_direction = "contradictory", "contradicts_class"
        else:
            ff_node_status, ff_direction = "neutral", "neutral"
        ff_reason = (
            "Conditional class-level log-likelihood ratio was calculated from "
            "target/background identity counts in an exact acquisition domain."
        )
    elif ff_status == "not_applicable":
        ff_node_status, ff_direction = "not_applicable", "neutral"
        ff_applicability = "not_applicable"
        ff_reason = str(
            feature_fragments.get("reason")
            or "No independently established class matched a validated rule."
        )
    else:
        ff_node_status, ff_direction = "unavailable", "neutral"
        ff_applicability = "unavailable"
        ff_reason = str(
            feature_fragments.get("reason")
            or "Exact instrument, activation, collision-energy metadata, or rule data were unavailable."
        )
    informative_rules = [
        item for item in feature_fragments.get("rule_evidence") or []
        if isinstance(item, dict)
        and item.get("status") in {"present", "eligible_absence"}
    ]
    nodes.append(_node(
        "feature_fragment_rule_support",
        "Condition-specific feature-fragment evidence",
        status=ff_node_status,
        direction=ff_direction,
        independence_group="query_ms2_rule_evidence",
        independent_measurement=False,
        used_for=["explanation"],
        reason=ff_reason,
        value={
            "module": str(feature_fragments.get("module") or "v15_feature_fragments"),
            "evidence_stage": str(
                feature_fragments.get("evidence_stage") or "msms_rule_support"
            ),
            "combined_log_likelihood_ratio": ff_log_lr,
            "likelihood_ratio": _number(feature_fragments.get("likelihood_ratio")),
            "informative_rule_count": len(informative_rules),
            "independent_fragmentation_groups": ff_groups,
            "instrument_domain_id": str(
                ((feature_fragments.get("acquisition_domain") or {}).get(
                    "instrument_domain_id"
                )) or ""
            ),
            "runtime_export_gate": dict(
                feature_fragments.get("runtime_export_gate") or {}
            ),
            "identity_rank_influence": False,
            "exact_identity_confirmed": False,
            "standard_confirmed": False,
        },
        applicability=ff_applicability,
        warnings=[
            "derived_from_same_query_ms2",
            "does_not_confirm_identity",
            "does_not_affect_ranking_or_release",
        ],
    ))

    if match.get("isomer_ambiguous"):
        isomer_status, isomer_direction = "contradictory", "limits_exact_structure"
        isomer_reason = (
            "Competitive same-formula candidates prevent a single exact-structure claim."
        )
    else:
        isomer_status, isomer_direction = "supporting", "supports_current_scope"
        isomer_reason = "No competitive same-formula alternative passed the stored guard."
    nodes.append(_node(
        "candidate_competition",
        "Candidate and isomer competition",
        status=isomer_status,
        direction=isomer_direction,
        independence_group="same_query_candidate_pool",
        independent_measurement=False,
        used_for=["confidence", "claim_scope"],
        reason=isomer_reason,
        value={
            "candidate_margin": _number(match.get("candidate_margin")),
            "isomer_ambiguous": bool(match.get("isomer_ambiguous")),
            "isomer_group_size": int(match.get("isomer_group_size") or 0),
            "ambiguity_type": str(match.get("isomer_ambiguity_type") or "none"),
        },
    ))

    blank_used = bool(row.get("blank_used")) or str(row.get("blank_status") or "") not in {
        "", "not_requested", "not_run"
    }
    if sample_role == "blank":
        blank_status, blank_direction, blank_app = (
            "role_veto", "background_context_only", "applicable"
        )
        blank_reason = "The query is declared as a blank; analyte-presence claims are refused."
    elif not blank_used:
        blank_status, blank_direction, blank_app = (
            "not_run", "neutral", "unavailable"
        )
        blank_reason = "No independent blank comparison was requested; absence is not evidence."
    elif bool(row.get("blank_clean_view_excluded")):
        blank_status, blank_direction, blank_app = (
            "contradictory", "contradicts_sample_specificity", "applicable"
        )
        blank_reason = "The feature recurs in the blank without sufficient sample enrichment."
    elif bool(row.get("blank_flagged")):
        blank_status, blank_direction, blank_app = (
            "context_warning", "limits_sample_specificity", "applicable"
        )
        blank_reason = (
            "The feature is blank-associated, but the declared blank role/policy "
            "does not justify a hard analyte veto."
        )
    elif str(row.get("blank_status") or "") == "no_blank_feature":
        blank_status, blank_direction, blank_app = (
            "supporting", "supports_sample_specificity", "applicable"
        )
        blank_reason = "Applicable blank replicates did not contain a matching feature."
    else:
        blank_status, blank_direction, blank_app = (
            "unevaluable", "neutral", "applicable"
        )
        blank_reason = "Blank evidence was requested but did not support a directional claim."
    nodes.append(_node(
        "blank_contrast",
        "Independent blank contrast",
        status=blank_status,
        direction=blank_direction,
        independence_group="blank_acquisition",
        independent_measurement=blank_app == "applicable",
        used_for=["sample_specificity", "claim_scope"],
        reason=blank_reason,
        value={
            "blank_status": str(row.get("blank_status") or "not_requested"),
            "blank_flagged": bool(row.get("blank_flagged")),
            "sample_to_blank_ratio": _number(row.get("blank_sample_ratio")),
            "replicates_applicable": int(row.get("blank_replicates_applicable") or 0),
        },
        applicability=blank_app,
    ))

    rt = match.get("rt_evidence") or {}
    rt_eligible = bool(match.get("rt_confirmed") or rt.get("eligible_for_exact_identity"))
    if rt_eligible:
        rt_status, rt_direction, rt_app = "supporting", "supports_exact_identity", "applicable"
        rt_reason = "A same-method local standard provides orthogonal RT evidence."
    elif rt:
        rt_status, rt_direction, rt_app = "insufficient", "neutral", "applicable"
        rt_reason = "RT metadata exists but is not eligible for exact-identity confirmation."
    else:
        rt_status, rt_direction, rt_app = "not_run", "neutral", "unavailable"
        rt_reason = "No same-method local RT standard was available."
    nodes.append(_node(
        "local_rt_standard",
        "Same-method local RT standard",
        status=rt_status,
        direction=rt_direction,
        independence_group="local_standard_rt",
        independent_measurement=rt_eligible,
        used_for=["confidence", "isomer_resolution"],
        reason=rt_reason,
        value={
            "status": str(rt.get("status") or "not_run"),
            "delta_seconds": _number(rt.get("rt_delta_sec")),
            "method_status": str(rt.get("method_status") or "unverified"),
        },
        applicability=rt_app,
    ))

    feature_spectra = int(row.get("precursor_feature_spectra") or 1)
    feature_consistency = bool(row.get("precursor_feature_identity_consistent"))
    nodes.append(_node(
        "within_injection_stability",
        "Within-injection stepped-energy stability",
        status=(
            "supporting" if feature_consistency
            else "conflicting" if feature_spectra >= 2
            else "single_observation"
        ),
        direction=(
            "supports_repeatability" if feature_consistency
            else "limits_repeatability" if feature_spectra >= 2
            else "neutral"
        ),
        independence_group=f"precursor_event:{row.get('precursor_feature_id') or row.get('spectrum_id')}",
        independent_measurement=False,
        used_for=["stability_audit"],
        reason=(
            "Stepped collision energies share one injection/precursor event and "
            "are collapsed to one vote; they are not independent replicates."
        ),
        value={
            "spectra": feature_spectra,
            "collision_energies": list(row.get("precursor_feature_collision_energies") or []),
            "identity_consensus": _number(row.get("precursor_feature_identity_consensus")),
        },
        warnings=["pseudoreplicate_guard_applied"],
    ))

    cohort_evidence = cohort_evidence or {}
    if cohort_evidence:
        nodes.append(_node(
            "independent_injection_consensus",
            "Independent-injection consensus",
            status=str(cohort_evidence.get("status") or "unevaluable"),
            direction=str(cohort_evidence.get("direction") or "neutral"),
            independence_group="declared_injection_groups",
            independent_measurement=True,
            used_for=["release_gate", "stability_audit"],
            reason=str(cohort_evidence.get("reason") or ""),
            value=cohort_evidence,
        ))

    audit = _candidate_release_audit(match, row)
    ambiguity = _ambiguity_set(match)
    hard_blockers = list(audit["blockers"])
    if sample_role == "blank":
        hard_blockers.append("declared_blank_query")
    if match.get("isomer_ambiguous"):
        claim_status = "isomer_group_only"
        claim_scope = "molecular_formula_isomer_group"
    elif match.get("rt_confirmed") and not hard_blockers:
        claim_status = "orthogonally_supported_exact_candidate"
        claim_scope = "exact_structure_candidate"
    elif audit["eligible"]:
        claim_status = "provisional_needs_independent_injection"
        claim_scope = "connectivity_level_library_candidate"
    elif str(match.get("confidence") or "") in {"High", "Medium"}:
        claim_status = "reported_but_governance_blocked"
        claim_scope = "spectral_candidate_only"
    else:
        claim_status = "insufficient_for_release"
        claim_scope = "putative_or_class_level"

    if sample_role == "blank":
        claim_status = "blank_background_context_only"
        claim_scope = "qc_background_identity"
    elif bool(row.get("blank_clean_view_excluded")):
        claim_status = "blank_conflicted_candidate"
        claim_scope = "background_or_carryover_candidate"

    applicable_nodes = [node for node in nodes if node["applicability"] == "applicable"]
    independent_groups = {
        node["independence_group"] for node in applicable_nodes
        if node["independent_measurement"] and node["status"] == "supporting"
    }
    warnings = sorted({
        warning for node in nodes for warning in node.get("warnings") or []
    })
    return {
        "schema_version": SCHEMA_VERSION,
        "policy_id": POLICY_ID,
        "ranking_changed": False,
        "confidence_field_overwritten": False,
        "claim_status": claim_status,
        "claim_scope": claim_scope,
        "release_eligible_before_replication": audit["eligible"],
        "release_blockers": sorted(set(hard_blockers)),
        "applicable_evidence_nodes": len(applicable_nodes),
        "independent_support_groups": len(independent_groups),
        "independence_groups": sorted(independent_groups),
        "warnings": warnings,
        "ambiguity_set": ambiguity,
        "ambiguity_set_semantics": (
            "score-near and same-formula candidates; no conformal coverage guarantee"
        ),
        "coverage_guarantee": False,
        "counterfactuals": {
            "spectral_only_reports": str(match.get("confidence") or "") in {"High", "Medium"},
            "expert_unscoped_reports": (
                str(match.get("confidence") or "") in {"High", "Medium"}
                and (_number(match.get("expert_score")) or 0.0) >= MIN_SCOPED_EXPERT_SCORE
            ),
            "expert_scoped_reports": (
                str(match.get("confidence") or "") in {"High", "Medium"}
                and expert["passes_locked_gate"]
            ),
            "full_governed_release_eligible": audit["eligible"],
        },
        "nodes": nodes,
    }


def annotate_evidence_passports(
    results: list[dict], *, sample_role: str = "sample"
) -> list[dict]:
    """Attach a deterministic passport to every matched result in-place."""
    for row in results or []:
        if not isinstance(row, dict) or not isinstance(row.get("match"), dict):
            continue
        row["match"]["evidence_governance"] = build_evidence_passport(
            row, sample_role=sample_role
        )
    return results


def summarize_evidence_governance(results: list[dict]) -> dict:
    passports = []
    identities: dict[str, dict] = {}
    for row in results or []:
        match = row.get("match") if isinstance(row, dict) else None
        if not isinstance(match, dict):
            continue
        passport = match.get("evidence_governance")
        if not isinstance(passport, dict):
            passport = build_evidence_passport(row)
        if not isinstance(passport, dict):
            continue
        passports.append(passport)
        identity = _identity(match)
        if identity:
            identities.setdefault(identity, passport)
    claim_counts = Counter(str(item.get("claim_status") or "unknown") for item in passports)
    identity_claim_counts = Counter(
        str(item.get("claim_status") or "unknown") for item in identities.values()
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "policy_id": POLICY_ID,
        "matched_spectra_audited": len(passports),
        "unique_identities_audited": len(identities),
        "claim_status_spectra": dict(sorted(claim_counts.items())),
        "claim_status_unique_identities": dict(sorted(identity_claim_counts.items())),
        "release_eligible_before_replication_spectra": sum(
            bool(item.get("release_eligible_before_replication")) for item in passports
        ),
        "release_eligible_before_replication_unique_identities": sum(
            bool(item.get("release_eligible_before_replication")) for item in identities.values()
        ),
        "pseudoreplicate_guarded_spectra": sum(
            "pseudoreplicate_guard_applied" in (item.get("warnings") or [])
            for item in passports
        ),
        "ranking_changed": False,
        "release_requires": (
            "same-method RT confirmation or agreement across at least two "
            "declared independent injection groups"
        ),
    }


def _feature_group_key(row: dict, index: int) -> str:
    lcms = str(row.get("lcms_feature_id") or "")
    if lcms:
        return f"lcms:{lcms}"
    precursor = str(row.get("precursor_feature_id") or "")
    if precursor:
        return f"precursor:{precursor}"
    return f"spectrum:{row.get('spectrum_id') or index}"


def build_feature_consensus(results: list[dict], task_id: str) -> list[dict]:
    """Collapse one precursor event to one shadow vote and an ambiguity set."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for index, row in enumerate(results or []):
        if (
            isinstance(row, dict)
            and not bool(row.get("blank_excluded"))
            and isinstance(row.get("match"), dict)
        ):
            groups[_feature_group_key(row, index)].append(row)

    features: list[dict] = []
    for feature_id, rows in groups.items():
        states: dict[str, dict] = {}
        n_spectra = len(rows)
        for row in rows:
            match = row["match"]
            per_spectrum_seen: set[str] = set()
            for rank, candidate in enumerate(
                [match, *(match.get("top_candidates") or [])], start=1
            ):
                if not isinstance(candidate, dict):
                    continue
                identity = _connectivity_identity(candidate)
                if not identity or identity in per_spectrum_seen:
                    continue
                per_spectrum_seen.add(identity)
                score = _number(candidate.get("rank_score"))
                if score is None:
                    score = _number(candidate.get("final_score")) or 0.0
                state = states.setdefault(identity, {
                    "identity_key": identity,
                    "top1_votes": 0,
                    "ranks": [],
                    "scores": [],
                    "observations": [],
                    "full_identity_keys": set(),
                })
                state["top1_votes"] += int(rank == 1)
                state["ranks"].append(rank)
                state["scores"].append(score)
                state["observations"].append((candidate, row))
                full_identity = _identity(candidate)
                if full_identity:
                    state["full_identity_keys"].add(full_identity)

        ranked: list[dict] = []
        for state in states.values():
            observed = len(state["ranks"])
            representative, representative_row = max(
                state["observations"],
                key=lambda pair: (
                    _number(pair[0].get("rank_score"))
                    if _number(pair[0].get("rank_score")) is not None
                    else _number(pair[0].get("final_score")) or 0.0
                ),
            )
            rrf = sum(1.0 / (RRF_OFFSET + rank) for rank in state["ranks"]) / n_spectra
            mean_with_missing = sum(state["scores"]) / n_spectra
            audits = [
                _candidate_release_audit(candidate, row)
                for candidate, row in state["observations"]
            ]
            item = _candidate_summary(representative, state["identity_key"])
            item.update({
                "top1_votes": state["top1_votes"],
                "observed_spectra": observed,
                "feature_spectra": n_spectra,
                "coverage_fraction": round(observed / n_spectra, 4),
                "rrf20": round(rrf, 6),
                "mean_rank_score_with_missing_zero": round(mean_with_missing, 4),
                "eligible_observations": sum(audit["eligible"] for audit in audits),
                "eligible_for_replicate_gate": any(audit["eligible"] for audit in audits),
                "release_blockers": sorted({
                    blocker for audit in audits for blocker in audit["blockers"]
                }),
                "isomer_ambiguous": any(
                    bool(candidate.get("isomer_ambiguous"))
                    for candidate, _ in state["observations"]
                ),
                "connectivity_variant_count": len(state["full_identity_keys"]),
            })
            ranked.append(item)

        ranked.sort(key=lambda item: (
            -int(item["top1_votes"]),
            -float(item["coverage_fraction"]),
            -float(item["rrf20"]),
            -float(item["mean_rank_score_with_missing_zero"]),
            str(item["identity_key"]),
        ))
        lead = ranked[0] if ranked else None
        ambiguity = []
        if lead:
            for item in ranked:
                same_formula = bool(
                    lead.get("formula")
                    and item.get("formula") == lead.get("formula")
                )
                if (
                    item is lead
                    or item["top1_votes"] > 0
                    or same_formula
                    or (
                        lead["mean_rank_score_with_missing_zero"]
                        - item["mean_rank_score_with_missing_zero"] <= 0.05
                    )
                ):
                    ambiguity.append(item)
                if len(ambiguity) >= MAX_AMBIGUITY_SET:
                    break
        rt = _median(row.get("retention_time") for row in rows)
        mz = _median(row.get("precursor_mz") for row in rows)
        features.append({
            "task_id": task_id,
            "feature_id": feature_id,
            "feature_basis": "continuous_ms1" if feature_id.startswith("lcms:") else "ms2_precursor_event",
            "polarity": str(rows[0].get("ion_mode") or "").upper(),
            "precursor_mz": round(mz, 6) if mz is not None else None,
            "retention_time_min": round(rt, 6) if rt is not None else None,
            "n_spectra": n_spectra,
            "collision_energies": sorted({
                round(value, 4)
                for value in (_number(row.get("collision_energy")) for row in rows)
                if value is not None
            }),
            "consensus_lead": lead,
            "ambiguity_set": ambiguity,
            "candidate_count": len(ranked),
            "ranking_changed": False,
            "shadow_consensus_only": True,
            "pseudoreplicates_collapsed": max(0, n_spectra - 1),
        })
    features.sort(key=lambda item: (
        item["polarity"],
        float(item["retention_time_min"] or math.inf),
        float(item["precursor_mz"] or math.inf),
        item["feature_id"],
    ))
    return features


def _feature_distance(left: dict, right: dict, mz_ppm: float, rt_window: float):
    if left.get("polarity") != right.get("polarity"):
        return None
    lmz, rmz = _number(left.get("precursor_mz")), _number(right.get("precursor_mz"))
    lrt, rrt = _number(left.get("retention_time_min")), _number(right.get("retention_time_min"))
    if None in {lmz, rmz, lrt, rrt} or not lmz or not rmz:
        return None
    # Use the pair mean as the denominator.  A left-hand denominator makes the
    # same pair infinitesimally direction dependent and is inappropriate for a
    # cohort relation that must be invariant to submission order.
    ppm = abs(lmz - rmz) / ((lmz + rmz) / 2.0) * 1e6
    rt_delta = abs(lrt - rrt)
    if ppm > mz_ppm or rt_delta > rt_window:
        return None
    return ppm / mz_ppm + rt_delta / rt_window, ppm, rt_delta


def _feature_alignment_key(feature: dict) -> tuple:
    """Stable key containing no candidate identity or score information."""
    return (
        str(feature.get("polarity") or ""),
        float(_number(feature.get("retention_time_min")) or math.inf),
        float(_number(feature.get("precursor_mz")) or math.inf),
        str(feature.get("task_id") or ""),
        str(feature.get("feature_id") or ""),
    )


def _cluster_features(
    features_by_task: dict[str, list[dict]], *, mz_ppm: float, rt_window: float
) -> list[list[dict]]:
    """Order-invariant complete-link alignment with one feature per task.

    Candidate edges are formed exclusively from polarity, precursor m/z and RT;
    compound identities and match scores are deliberately absent.  A globally
    sorted constrained-union pass then joins only components whose task sets are
    disjoint and whose every cross-pair remains inside the locked window.  This
    prevents both task-order drift and single-linkage chaining.
    """
    features = sorted(
        (
            feature
            for task_id in sorted(features_by_task)
            for feature in features_by_task[task_id]
        ),
        key=_feature_alignment_key,
    )
    if not features:
        return []

    # A narrow m/z sweep avoids the quadratic all-pairs cost for larger cohorts.
    # The exact relation is still evaluated by _feature_distance before an edge
    # enters the graph.
    edges: list[tuple[float, float, float, int, int]] = []
    by_polarity: dict[str, list[tuple[float, int]]] = defaultdict(list)
    for index, feature in enumerate(features):
        mz = _number(feature.get("precursor_mz"))
        rt = _number(feature.get("retention_time_min"))
        if mz is not None and mz > 0 and rt is not None:
            by_polarity[str(feature.get("polarity") or "")].append((mz, index))
    for polarity in sorted(by_polarity):
        ordered = sorted(by_polarity[polarity])
        for position, (left_mz, left_index) in enumerate(ordered):
            left = features[left_index]
            for right_mz, right_index in ordered[position + 1:]:
                # Symmetric ppm distance grows monotonically as right_mz grows.
                ppm_lower_bound = (
                    (right_mz - left_mz) / ((left_mz + right_mz) / 2.0) * 1e6
                )
                if ppm_lower_bound > mz_ppm:
                    break
                right = features[right_index]
                if left.get("task_id") == right.get("task_id"):
                    continue
                distance = _feature_distance(left, right, mz_ppm, rt_window)
                if distance is not None:
                    edges.append((*distance, left_index, right_index))
    edges.sort()

    parent = list(range(len(features)))
    component_members: dict[int, list[int]] = {
        index: [index] for index in range(len(features))
    }
    component_tasks: dict[int, set[str]] = {
        index: {str(feature.get("task_id") or "")}
        for index, feature in enumerate(features)
    }

    def find(index: int) -> int:
        root = index
        while parent[root] != root:
            root = parent[root]
        while parent[index] != index:
            next_index = parent[index]
            parent[index] = root
            index = next_index
        return root

    for _, _, _, left_index, right_index in edges:
        left_root, right_root = find(left_index), find(right_index)
        if left_root == right_root:
            continue
        if component_tasks[left_root] & component_tasks[right_root]:
            continue
        # Complete-link compatibility prevents a chain A~B~C from joining when
        # A and C themselves fall outside the fixed m/z/RT window.
        if any(
            _feature_distance(
                features[left_member], features[right_member], mz_ppm, rt_window
            ) is None
            for left_member in component_members[left_root]
            for right_member in component_members[right_root]
        ):
            continue
        keep, absorb = sorted((left_root, right_root))
        parent[absorb] = keep
        component_members[keep] = sorted(
            component_members[keep] + component_members[absorb]
        )
        component_tasks[keep] |= component_tasks[absorb]
        del component_members[absorb]
        del component_tasks[absorb]

    clusters = [
        sorted((features[index] for index in member_indices), key=_feature_alignment_key)
        for _, member_indices in sorted(component_members.items())
    ]
    clusters.sort(key=lambda cluster: (
        str(cluster[0].get("polarity") or ""),
        float(_median(item.get("retention_time_min") for item in cluster) or math.inf),
        float(_median(item.get("precursor_mz") for item in cluster) or math.inf),
        tuple(_feature_alignment_key(item) for item in cluster),
    ))
    return clusters


def build_replicate_consensus(
    task_results: dict[str, list[dict]],
    *,
    task_roles: dict[str, str] | None = None,
    independence_groups: dict[str, str] | None = None,
    mz_ppm: float = COHORT_MZ_PPM,
    rt_window_min: float = COHORT_RT_WINDOW_MIN,
    min_independent_injections: int = MIN_INDEPENDENT_INJECTIONS,
    min_identity_agreement: float = MIN_IDENTITY_AGREEMENT,
) -> dict:
    """Build injection-level consensus without merging spectra or reranking Top-1."""
    task_roles = task_roles or {}
    independence_groups = independence_groups or {}
    task_ids = sorted(task_results)
    features_by_task = {
        task_id: build_feature_consensus(task_results[task_id], task_id)
        for task_id in task_ids
    }
    clusters = _cluster_features(
        features_by_task, mz_ppm=mz_ppm, rt_window=rt_window_min
    )
    all_blank_features = [
        feature
        for task_id, features in features_by_task.items()
        if task_roles.get(task_id, "sample") == "blank"
        for feature in features
    ]

    output_clusters = []
    released_identities: set[str] = set()
    released_scope_units: set[str] = set()
    counterfactual_identities = {
        "spectral_only": set(),
        "expert_unscoped": set(),
        "expert_scoped": set(),
        "full_governed": set(),
    }
    for task_id in task_ids:
        results = task_results[task_id]
        for row in results or []:
            match = row.get("match") if isinstance(row, dict) else None
            if not isinstance(match, dict):
                continue
            identity = _connectivity_identity(match)
            if not identity:
                continue
            confidence_ok = str(match.get("confidence") or "") in {"High", "Medium"}
            if confidence_ok:
                counterfactual_identities["spectral_only"].add(identity)
            if confidence_ok and (_number(match.get("expert_score")) or 0.0) >= MIN_SCOPED_EXPERT_SCORE:
                counterfactual_identities["expert_unscoped"].add(identity)
            if confidence_ok and _scoped_expert(match)["passes_locked_gate"]:
                counterfactual_identities["expert_scoped"].add(identity)

    for cluster_index, members in enumerate(clusters, start=1):
        sample_members = [
            member for member in members
            if task_roles.get(member["task_id"], "sample") != "blank"
        ]
        blank_members = [
            member for member in members
            if task_roles.get(member["task_id"], "sample") == "blank"
        ]
        votes: dict[str, list[dict]] = defaultdict(list)
        for member in sample_members:
            lead = member.get("consensus_lead")
            if isinstance(lead, dict) and lead.get("identity_key"):
                votes[str(lead["identity_key"])].append(member)
        if votes:
            dominant_identity, dominant_members = max(
                votes.items(),
                key=lambda item: (
                    len({
                        independence_groups.get(member["task_id"], member["task_id"])
                        for member in item[1]
                    }),
                    len(item[1]),
                    sum(
                        float((member.get("consensus_lead") or {}).get("rrf20") or 0.0)
                        for member in item[1]
                    ),
                    item[0],
                ),
            )
        else:
            dominant_identity, dominant_members = "", []
        independent_groups = {
            independence_groups.get(member["task_id"], member["task_id"])
            for member in dominant_members
        }
        evaluable_groups = {
            independence_groups.get(member["task_id"], member["task_id"])
            for member in sample_members
        }
        agreement = (
            len(independent_groups) / len(evaluable_groups)
            if evaluable_groups else 0.0
        )
        qualified_groups = {
            independence_groups.get(member["task_id"], member["task_id"])
            for member in dominant_members
            if bool((member.get("consensus_lead") or {}).get(
                "eligible_for_replicate_gate"
            ))
        }
        # Blank alignment is searched against all declared blank features, not
        # only the one-to-one transport cluster.  A dense DDA region may contain
        # several nearly coincident precursor events; letting a greedy sample-
        # sample pairing consume the nearest blank event would otherwise hide a
        # real background recurrence.  Proximity is still fixed without using
        # identity; identity is inspected only after the eligible neighbourhood
        # has been formed.
        sample_centre = {
            "polarity": members[0].get("polarity") if members else "",
            "precursor_mz": _median(
                member.get("precursor_mz") for member in sample_members
            ),
            "retention_time_min": _median(
                member.get("retention_time_min") for member in sample_members
            ),
        }
        nearby_blank_members = [
            member for member in all_blank_features
            if _feature_distance(
                sample_centre, member, mz_ppm, rt_window_min
            ) is not None
        ] if sample_members else []
        blank_conflict = any(
            any(
                str(candidate.get("identity_key") or "") == dominant_identity
                for candidate in member.get("ambiguity_set") or []
            )
            for member in nearby_blank_members
        ) if dominant_identity else False
        ambiguous = any(
            bool((member.get("consensus_lead") or {}).get("isomer_ambiguous"))
            for member in dominant_members
        )

        if not dominant_identity:
            action = "no_candidate"
        elif blank_conflict:
            action = "abstain_blank_recurrent"
        elif len(independent_groups) < min_independent_injections:
            action = "needs_more_independent_injections"
        elif agreement < min_identity_agreement:
            action = "abstain_identity_conflict"
        elif len(qualified_groups) < min_independent_injections:
            action = "abstain_scope_or_quality_gate"
        elif ambiguous:
            action = "release_isomer_group_candidate"
            released_scope_units.add(dominant_identity)
            counterfactual_identities["full_governed"].add(dominant_identity)
        else:
            action = "release_replicate_stable_candidate"
            released_identities.add(dominant_identity)
            released_scope_units.add(dominant_identity)
            counterfactual_identities["full_governed"].add(dominant_identity)

        lead = None
        if dominant_members:
            lead = max(
                (member["consensus_lead"] for member in dominant_members),
                key=lambda candidate: (
                    int(candidate.get("eligible_for_replicate_gate", False)),
                    int(candidate.get("top1_votes") or 0),
                    float(candidate.get("rrf20") or 0.0),
                ),
            )
        output_clusters.append({
            "cluster_id": f"EGC{cluster_index:06d}",
            "polarity": str(members[0].get("polarity") or "") if members else "",
            "precursor_mz": _median(member.get("precursor_mz") for member in members),
            "retention_time_min": _median(
                member.get("retention_time_min") for member in members
            ),
            "member_tasks": [member["task_id"] for member in members],
            "sample_member_tasks": [member["task_id"] for member in sample_members],
            "blank_member_tasks": sorted({
                member["task_id"] for member in nearby_blank_members
            }),
            "dominant_identity": dominant_identity,
            "dominant_candidate": lead,
            "independent_injection_groups": len(independent_groups),
            "evaluable_injection_groups": len(evaluable_groups),
            "qualified_injection_groups": len(qualified_groups),
            "identity_agreement": round(agreement, 4),
            "blank_conflict": blank_conflict,
            "action": action,
            "ranking_changed": False,
            "member_features": [
                {
                    "task_id": member["task_id"],
                    "feature_id": member["feature_id"],
                    "n_spectra": member["n_spectra"],
                    "pseudoreplicates_collapsed": member["pseudoreplicates_collapsed"],
                    "consensus_lead": member.get("consensus_lead"),
                    "ambiguity_set": member.get("ambiguity_set"),
                }
                for member in members
            ],
        })

    action_counts = Counter(cluster["action"] for cluster in output_clusters)
    return {
        "schema_version": SCHEMA_VERSION,
        "policy_id": POLICY_ID,
        "method": "feature_collapsed_independent_injection_stability_selection",
        "task_ids": task_ids,
        "task_roles": {
            task_id: task_roles.get(task_id, "sample") for task_id in task_ids
        },
        "independence_groups": {
            task_id: independence_groups.get(task_id, task_id) for task_id in task_ids
        },
        "independence_scope": "submitted_injection_groups_not_biological_replicates",
        "parameters": {
            "mz_ppm": mz_ppm,
            "rt_window_min": rt_window_min,
            "min_independent_injections": min_independent_injections,
            "min_identity_agreement": min_identity_agreement,
            "scoped_expert_threshold": MIN_SCOPED_EXPERT_SCORE,
        },
        "summary": {
            "input_tasks": len(task_ids),
            "input_features": sum(len(items) for items in features_by_task.values()),
            "aligned_feature_clusters": len(output_clusters),
            "action_counts": dict(sorted(action_counts.items())),
            "released_unique_identities": len(released_identities),
            "released_scope_unique_units": len(released_scope_units),
            "counterfactual_unique_identities": {
                key: len(value) for key, value in counterfactual_identities.items()
            },
            "pseudoreplicate_scans_collapsed": sum(
                int(feature.get("pseudoreplicates_collapsed") or 0)
                for features in features_by_task.values() for feature in features
            ),
            "ranking_changed": False,
        },
        "clusters": output_clusters,
    }
