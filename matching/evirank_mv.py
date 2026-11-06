"""V12-Flagship / EviRank-MV shadow ranker.

EviRank-MV is a bounded residual challenger layered on an immutable
V11-Flagship entity artifact.  The default runtime has no frozen residual
model.  EviRecall candidates enter an audited shadow candidate union, but they
remain ineligible for residual-model scoring and cannot become Top-1 unless a
locked source-paired, multi-collision-energy consensus explicitly qualifies
the same-formula switch.  A residual model may be applied to native V11
candidates only when its training contract proves connectivity-disjoint
splits, retention of same-identity cross-method cases, and non-independence of
technical repeats.

Scores produced here are ranking utilities, never probabilities.  Missing or
inapplicable views are neutral and cannot penalize a candidate.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from matching.library_quality import normalize_formula_text, normalize_inchikey


EVIRANK_SCHEMA_VERSION = "2.0.0"
MODULE_ID = "v12_flagship_evirank_mv"
EXPECTED_EVIRECALL_POLICY_SHA256 = (
    "0bd5ce8b9a7c03771287841c491fcc346162d97097bca0ea06fbf1f2fcb2fa6a"
)

NATIVE_FEATURES = (
    "v11_entity_score",
    "ce_support_mean",
    "ce_support_median",
    "ce_coverage",
    "ce_top1_fraction",
    "ce_support_floor",
    "ce_support_range",
    "ce_support_stability",
    "reference_source_diversity",
    "chemistry_consistency",
    "candidate_structure_available",
)
OPTIONAL_FEATURES = (
    "forward_spectrum_ensemble",
    "fragment_graph_fit",
    "retention_time_fit",
    "ccs_fit",
    "isotope_fit",
    "ion_identity_fit",
)
FEATURE_ORDER = (*NATIVE_FEATURES, *OPTIONAL_FEATURES)


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


def _clip(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, float(value)))


def _candidate_connectivity(candidate: Mapping[str, Any]) -> str | None:
    direct = str(candidate.get("connectivity_identity") or "").strip()
    if direct:
        return direct
    identity_key = str(candidate.get("identity_key") or "").strip()
    if identity_key.startswith(("std2d:", "ik14:")):
        return identity_key.split(":", 1)[1] or None
    display = candidate.get("candidate") or candidate
    key = normalize_inchikey(display.get("inchikey"))
    return key.split("-", 1)[0] if key else None


def _candidate_identity_key(candidate: Mapping[str, Any]) -> str | None:
    identity_key = str(candidate.get("identity_key") or "").strip()
    if identity_key.startswith(("std2d:", "ik14:")):
        return identity_key
    connectivity = _candidate_connectivity(candidate)
    return f"std2d:{connectivity}" if connectivity else None


def default_evirank_policy_path() -> Path:
    return Path(__file__).resolve().parent.parent / "config" / "v12_flagship_policy.json"


def _validate_evirank_policy(policy: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(policy))
    if payload.get("deployment_state") != "shadow":
        raise ValueError("V12 EviRank-MV must remain shadow until external promotion")
    if payload.get("display_name") != "V12-Flagship / EviRank-MV":
        raise ValueError("Unexpected V12 display name")
    residual = payload.get("residual_model") or {}
    maximum = _finite(residual.get("maximum_absolute_delta"))
    if maximum is None or not 0.0 <= maximum <= 0.10:
        raise ValueError("V12 residual bound must be between 0 and 0.10")
    method = payload.get("method_case_contract") or {}
    if not method.get("same_identity_different_methods_are_retained"):
        raise ValueError("V12 must retain the same identity across detection methods")
    if method.get("technical_repeats_increase_independent_identity_n"):
        raise ValueError("Technical repeats cannot increase independent identity n")
    candidate_recall = payload.get("candidate_recall") or {}
    if candidate_recall.get("module_id") != "v12_evirank_candidate_recall":
        raise ValueError("V12 candidate-recall module id is not locked")
    if candidate_recall.get("algorithm_id") != "mssystem-v12-evirecall-mv-0.2.0":
        raise ValueError("V12 candidate-recall algorithm id is not locked")
    if candidate_recall.get("policy_sha256") != EXPECTED_EVIRECALL_POLICY_SHA256:
        raise ValueError("V12 candidate-recall frozen policy hash is not locked")
    if candidate_recall.get("deployment_state") != "shadow":
        raise ValueError("V12 candidate recall must remain shadow")
    if candidate_recall.get("expanded_candidates_enter_shadow_candidate_union") is not True:
        raise ValueError("Expanded candidates must enter the audited shadow candidate union")
    if candidate_recall.get("expanded_candidates_enter_residual_model_scoring") is not False:
        raise ValueError("Expanded candidates may not enter residual-model scoring")
    if candidate_recall.get(
        "expanded_candidate_selection_requires_qualified_source_paired_consensus"
    ) is not True:
        raise ValueError("Expanded candidate selection requires qualified source-paired consensus")
    if candidate_recall.get("source_paired_consensus_module_id") != "v12_source_paired_ce_consensus":
        raise ValueError("Source-paired consensus module id is not locked")
    if candidate_recall.get("source_paired_consensus_algorithm_id") != (
        "mssystem-v12-source-paired-ce-consensus-0.1.0"
    ):
        raise ValueError("Source-paired consensus algorithm id is not locked")
    if candidate_recall.get("production_effect") is not False:
        raise ValueError("V12 candidate recall may not affect production")
    forbidden = set(payload.get("forbidden_ranking_inputs") or [])
    required = {
        "compound_name",
        "ground_truth_identity",
        "standard_id",
        "expert_score",
        "classification",
        "biosynthetic_pathway",
        "candidate_pathway_alignment",
        "language_model_output",
        "legacy_confidence_label",
    }
    if not required <= forbidden:
        raise ValueError("V12 policy does not forbid every leakage-prone rank input")
    return payload


def load_evirank_policy(path: Path | None = None) -> dict[str, Any]:
    policy_path = path or default_evirank_policy_path()
    payload = json.loads(policy_path.read_text(encoding="utf-8"))
    return _validate_evirank_policy(payload)


def evirank_policy_sha256(policy: Mapping[str, Any]) -> str:
    return _sha256(policy)


def _validate_v11_payload(
    payload: Mapping[str, Any], policy: Mapping[str, Any]
) -> tuple[Mapping[str, Any], str]:
    metadata = payload.get("metadata")
    embedded_policy = payload.get("policy")
    entities = payload.get("entities")
    if not isinstance(metadata, Mapping):
        raise ValueError("V11 artifact metadata is missing")
    if not isinstance(embedded_policy, Mapping):
        raise ValueError("V11 artifact must embed its frozen policy")
    if not isinstance(entities, Sequence) or isinstance(entities, (str, bytes)):
        raise ValueError("V11 artifact entities are missing")
    expected = str((policy.get("base_algorithm") or {}).get("algorithm_id") or "")
    actual = str(metadata.get("algorithm_id") or "")
    if actual != expected:
        raise ValueError(f"V11 base algorithm mismatch: {actual or 'missing'}")
    embedded_hash = _sha256(embedded_policy)
    declared_hash = str(metadata.get("policy_sha256") or "")
    if embedded_hash != declared_hash:
        raise ValueError("V11 embedded policy hash does not match metadata")
    source_hash = str(metadata.get("source_semantic_sha256") or "")
    if len(source_hash) != 64:
        raise ValueError("V11 source semantic hash is missing")
    return metadata, source_hash


def _validate_model(
    model: Mapping[str, Any] | None, policy: Mapping[str, Any]
) -> tuple[dict[str, Any] | None, str | None]:
    if model is None:
        return None, None
    artifact = dict(model)
    if artifact.get("state") != "qualified_shadow":
        raise ValueError("EviRank residual model must declare qualified_shadow")
    if artifact.get("algorithm_id") != policy.get("algorithm_id"):
        raise ValueError("Residual model algorithm_id mismatch")
    if tuple(artifact.get("feature_order") or ()) != FEATURE_ORDER:
        raise ValueError("Residual model feature order does not match V12 contract")
    training = artifact.get("training_contract") or {}
    required_true = (
        "connectivity_identity_disjoint",
        "same_identity_different_methods_retained",
        "technical_repeats_nonindependent",
        "ground_truth_unavailable_at_inference",
    )
    if not all(training.get(key) is True for key in required_true):
        raise ValueError("Residual model training contract is not leakage-safe")
    maximum = _finite(artifact.get("maximum_absolute_delta"))
    policy_maximum = float((policy.get("residual_model") or {})["maximum_absolute_delta"])
    if maximum is None or maximum < 0.0 or maximum > policy_maximum:
        raise ValueError("Residual model exceeds the V12 delta bound")
    temperature = _finite(artifact.get("temperature"))
    if temperature is None or temperature <= 0.0:
        raise ValueError("Residual model temperature must be positive")
    for field in ("feature_mean", "feature_scale", "coefficients"):
        values = artifact.get(field)
        if not isinstance(values, Mapping) or set(values) != set(FEATURE_ORDER):
            raise ValueError(f"Residual model {field} is incomplete")
        if any(_finite(values[key]) is None for key in FEATURE_ORDER):
            raise ValueError(f"Residual model {field} contains non-finite values")
    if any(float(artifact["feature_scale"][key]) <= 0.0 for key in FEATURE_ORDER):
        raise ValueError("Residual model feature scales must be positive")
    return artifact, _sha256(artifact)


def _validate_candidate_recall_bundle(
    bundle: Mapping[str, Any] | None,
    *,
    v11_payload: Mapping[str, Any],
    source_hash: str,
    detection_method_id: str,
    entity_ids: set[str],
    policy: Mapping[str, Any],
) -> tuple[dict[str, Mapping[str, Any]], str | None, Mapping[str, Any] | None]:
    if bundle is None:
        return {}, None, None
    metadata = bundle.get("metadata")
    embedded_policy = bundle.get("policy")
    invariants = bundle.get("invariants")
    entities = bundle.get("entities")
    if not isinstance(metadata, Mapping):
        raise ValueError("Candidate-recall bundle metadata is missing")
    if not isinstance(embedded_policy, Mapping):
        raise ValueError("Candidate-recall bundle must embed its frozen policy")
    if not isinstance(invariants, Mapping):
        raise ValueError("Candidate-recall bundle invariants are missing")
    if not isinstance(entities, Sequence) or isinstance(entities, (str, bytes)):
        raise ValueError("Candidate-recall bundle entities are missing")
    recall_contract = policy.get("candidate_recall") or {}
    if metadata.get("module_id") != recall_contract.get("module_id"):
        raise ValueError("Candidate-recall bundle module id mismatch")
    if metadata.get("algorithm_id") != recall_contract.get("algorithm_id"):
        raise ValueError("Candidate-recall bundle algorithm id mismatch")
    if metadata.get("deployment_state") != "shadow":
        raise ValueError("Candidate-recall bundle must remain shadow")
    if metadata.get("production_effect") is not False:
        raise ValueError("Candidate-recall bundle declares a production effect")
    if metadata.get("base_source_semantic_sha256") != source_hash:
        raise ValueError("Candidate-recall bundle is bound to another V11 source")
    if metadata.get("base_artifact_sha256") != _sha256(v11_payload):
        raise ValueError("Candidate-recall bundle is bound to another V11 artifact")
    if str(metadata.get("detection_method_id") or "") != detection_method_id:
        raise ValueError("Candidate-recall bundle detection-method binding mismatch")
    if _sha256(embedded_policy) != str(metadata.get("policy_sha256") or ""):
        raise ValueError("Candidate-recall embedded policy hash mismatch")
    if str(metadata.get("policy_sha256") or "") != str(
        recall_contract.get("policy_sha256") or ""
    ):
        raise ValueError("Candidate-recall bundle does not use the frozen policy hash")
    if not str(metadata.get("library_snapshot_id") or "").strip():
        raise ValueError("Candidate-recall library snapshot binding is missing")
    if metadata.get("source_membership_partition_validated") is not True:
        raise ValueError("Candidate-recall source membership was not validated")
    required_invariants = {
        "shadow_only": True,
        "ground_truth_used_at_inference": False,
        "candidate_names_affect_retrieval_or_scoring": False,
        "expanded_candidates_enter_shadow_candidate_union": True,
        "expanded_candidates_enter_residual_model_scoring": False,
        "source_paired_consensus_ground_truth_used_at_inference": False,
        "provenance_strata_are_independent_sources": False,
        "mirrored_reference_content_profiles_create_extra_votes": False,
        "same_connectivity_different_detection_methods_preexcluded": False,
        "technical_repeats_increase_condition_count": False,
        "blind_window_widening": False,
        "blind_adduct_enumeration": False,
        "score_is_probability": False,
        "production_effect": False,
    }
    if any(invariants.get(key) is not value for key, value in required_invariants.items()):
        raise ValueError("Candidate-recall bundle invariants are incomplete")
    embedded_consensus_policy = embedded_policy.get("source_paired_consensus") or {}
    if embedded_consensus_policy.get("module_id") != recall_contract.get(
        "source_paired_consensus_module_id"
    ):
        raise ValueError("Candidate-recall source-paired module contract mismatch")
    if embedded_consensus_policy.get("algorithm_id") != recall_contract.get(
        "source_paired_consensus_algorithm_id"
    ):
        raise ValueError("Candidate-recall source-paired algorithm contract mismatch")
    v11_entities = {
        str(entity.get("entity_id") or ""): entity
        for entity in (v11_payload.get("entities") or [])
        if isinstance(entity, Mapping)
    }
    indexed: dict[str, Mapping[str, Any]] = {}
    for record in entities:
        if not isinstance(record, Mapping):
            raise ValueError("Every candidate-recall entity must be a JSON object")
        entity_id = str(record.get("entity_id") or "")
        if not entity_id or entity_id in indexed:
            raise ValueError("Candidate-recall entity ids must be unique and non-empty")
        if str(record.get("detection_method_id") or "") != detection_method_id:
            raise ValueError("Candidate-recall entity method binding mismatch")
        if record.get("production_effect") is not False:
            raise ValueError("Expanded candidate entity declares a production effect")
        if record.get("expanded_candidates_enter_shadow_candidate_union") is not True:
            raise ValueError("Expanded candidate entity is missing the shadow-union contract")
        if record.get("expanded_candidates_enter_residual_model_scoring") is not False:
            raise ValueError("Expanded candidate entity declares residual-model eligibility")
        identities: set[str] = set()
        for candidate in record.get("expanded_candidates") or []:
            identity = str(candidate.get("identity_key") or "")
            if not identity.startswith("std2d:") or identity in identities:
                raise ValueError("Expanded candidates need unique stable 2-D identities")
            if candidate.get("score_is_probability") is not False:
                raise ValueError("Candidate-recall utility cannot be a probability")
            if candidate.get("production_effect") is not False:
                raise ValueError("Expanded candidate declares a production effect")
            if candidate.get("candidate_origin") != "evirecall_expansion":
                raise ValueError("Expanded candidate origin is missing")
            if candidate.get("enters_shadow_candidate_union") is not True:
                raise ValueError("Expanded candidate does not enter the shadow union")
            if candidate.get("enters_residual_model_scoring") is not False:
                raise ValueError("Expanded candidate may not enter residual-model scoring")
            identities.add(identity)
        consensus = record.get("source_paired_consensus")
        if not isinstance(consensus, Mapping):
            raise ValueError("Candidate-recall entity source-paired consensus is missing")
        if consensus.get("module_id") != recall_contract.get(
            "source_paired_consensus_module_id"
        ):
            raise ValueError("Source-paired consensus module id mismatch")
        if consensus.get("algorithm_id") != recall_contract.get(
            "source_paired_consensus_algorithm_id"
        ):
            raise ValueError("Source-paired consensus algorithm id mismatch")
        if consensus.get("deployment_state") != "shadow":
            raise ValueError("Source-paired consensus must remain shadow")
        required_consensus_false = (
            "provenance_strata_are_independent_sources",
            "ground_truth_used_at_inference",
            "candidate_names_affect_scoring",
            "score_is_probability",
            "production_effect",
        )
        if any(consensus.get(key) is not False for key in required_consensus_false):
            raise ValueError("Source-paired consensus invariants are incomplete")
        required_consensus_true = (
            "technical_repeats_collapsed_by_median",
            "reference_copies_collapsed_within_provenance_stratum",
            "mirrored_reference_content_collapsed_across_provenance_strata",
        )
        if any(consensus.get(key) is not True for key in required_consensus_true):
            raise ValueError("Source-paired consensus collapse invariants are incomplete")
        status = str(consensus.get("status") or "")
        if status not in {"not_applicable", "baseline_retained", "qualified_switch"}:
            raise ValueError("Invalid source-paired consensus status")
        base_entity = v11_entities.get(entity_id) or {}
        base_candidates = [
            candidate
            for candidate in (base_entity.get("candidates") or [])
            if isinstance(candidate, Mapping)
        ]
        eligible_base_candidates = [
            candidate
            for candidate in base_candidates
            if not candidate.get("vetoed") and _candidate_identity_key(candidate)
        ]
        base_identities = {
            identity
            for candidate in base_candidates
            if (identity := _candidate_identity_key(candidate))
        }
        union_identities = base_identities | identities
        baseline_identity = str(consensus.get("baseline_identity_key") or "")
        selected_identity = str(consensus.get("selected_identity_key") or "")
        if baseline_identity and baseline_identity not in base_identities:
            raise ValueError("Source-paired baseline is not a V11 candidate")
        expected_baseline_identity = (
            _candidate_identity_key(eligible_base_candidates[0])
            if eligible_base_candidates
            else None
        )
        if baseline_identity and baseline_identity != expected_baseline_identity:
            raise ValueError("Source-paired baseline is not the V11 Top-1 candidate")
        reported_margin = _finite(consensus.get("baseline_margin"))
        if reported_margin is not None:
            if len(eligible_base_candidates) < 2:
                raise ValueError("Source-paired V11 margin has no runner-up")
            expected_margin = float(
                eligible_base_candidates[0].get("entity_score") or 0.0
            ) - float(eligible_base_candidates[1].get("entity_score") or 0.0)
            if expected_margin < 0.0 or not math.isclose(
                reported_margin, expected_margin, abs_tol=1e-8
            ):
                raise ValueError("Source-paired V11 margin does not match the base artifact")
        if selected_identity and selected_identity not in union_identities:
            raise ValueError("Source-paired selection is outside the shadow candidate union")
        comparisons = consensus.get("comparisons") or []
        if not isinstance(comparisons, Sequence) or isinstance(comparisons, (str, bytes)):
            raise ValueError("Source-paired comparisons must be a JSON array")
        thresholds = consensus.get("thresholds") or {}
        expected_thresholds = {
            "maximum_v11_top1_margin": float(
                embedded_consensus_policy["maximum_v11_top1_margin"]
            ),
            "minimum_same_formula_connectivities": int(
                embedded_consensus_policy["minimum_same_formula_connectivities"]
            ),
            "maximum_same_formula_connectivities": int(
                embedded_consensus_policy["maximum_same_formula_connectivities"]
            ),
            "maximum_expanded_same_formula_challengers": int(
                embedded_consensus_policy[
                    "maximum_expanded_same_formula_challengers"
                ]
            ),
            "minimum_common_provenance_strata": int(
                embedded_consensus_policy["minimum_common_provenance_strata"]
            ),
            "minimum_unique_provenance_evidence_profiles": int(
                embedded_consensus_policy[
                    "minimum_unique_provenance_evidence_profiles"
                ]
            ),
            "minimum_distinct_collision_energy_conditions": int(
                embedded_consensus_policy[
                    "minimum_distinct_collision_energy_conditions"
                ]
            ),
            "minimum_spectral_views": int(
                embedded_consensus_policy["minimum_spectral_views"]
            ),
            "minimum_paired_observations": int(
                embedded_consensus_policy["minimum_paired_observations"]
            ),
            "minimum_median_score_delta": float(
                embedded_consensus_policy["minimum_median_score_delta"]
            ),
            "minimum_observation_win_fraction": float(
                embedded_consensus_policy["minimum_observation_win_fraction"]
            ),
            "minimum_provenance_stratum_win_fraction": float(
                embedded_consensus_policy[
                    "minimum_provenance_stratum_win_fraction"
                ]
            ),
            "minimum_unique_provenance_profile_win_fraction": float(
                embedded_consensus_policy[
                    "minimum_unique_provenance_profile_win_fraction"
                ]
            ),
        }
        if thresholds != expected_thresholds:
            raise ValueError("Source-paired consensus thresholds differ from frozen policy")
        truly_qualified: set[str] = set()
        for comparison in comparisons:
            if not isinstance(comparison, Mapping):
                raise ValueError("Source-paired comparison must be a JSON object")
            challenger = str(comparison.get("challenger_identity_key") or "")
            if challenger not in union_identities:
                raise ValueError("Source-paired challenger is outside the candidate union")
            metrics_pass = (
                int(comparison.get("common_provenance_stratum_count") or 0)
                >= int(thresholds.get("minimum_common_provenance_strata") or 0)
                and int(
                    comparison.get("unique_provenance_evidence_profile_count") or 0
                )
                >= int(
                    thresholds.get("minimum_unique_provenance_evidence_profiles")
                    or 0
                )
                and int(comparison.get("distinct_collision_energy_condition_count") or 0)
                >= int(thresholds.get("minimum_distinct_collision_energy_conditions") or 0)
                and int(comparison.get("spectral_view_count") or 0)
                >= int(thresholds.get("minimum_spectral_views") or 0)
                and int(comparison.get("paired_observation_count") or 0)
                >= int(thresholds.get("minimum_paired_observations") or 0)
                and (_finite(comparison.get("median_score_delta")) or -math.inf)
                >= float(thresholds.get("minimum_median_score_delta") or math.inf)
                and (_finite(comparison.get("observation_win_fraction")) or 0.0)
                >= float(thresholds.get("minimum_observation_win_fraction") or math.inf)
                and (_finite(comparison.get("provenance_stratum_win_fraction")) or 0.0)
                >= float(
                    thresholds.get("minimum_provenance_stratum_win_fraction") or math.inf
                )
                and (
                    _finite(
                        comparison.get(
                            "unique_provenance_profile_win_fraction"
                        )
                    )
                    or 0.0
                )
                >= float(
                    thresholds.get(
                        "minimum_unique_provenance_profile_win_fraction"
                    )
                    or math.inf
                )
            )
            if comparison.get("qualified") is True:
                if not metrics_pass or comparison.get("failed_gates"):
                    raise ValueError("Source-paired qualified comparison fails its gates")
                truly_qualified.add(challenger)
        if status == "qualified_switch":
            baseline_margin = _finite(consensus.get("baseline_margin"))
            same_formula_count = int(
                consensus.get("same_formula_connectivity_count") or 0
            )
            if (
                baseline_margin is None
                or baseline_margin > expected_thresholds["maximum_v11_top1_margin"]
                or same_formula_count
                < expected_thresholds["minimum_same_formula_connectivities"]
                or same_formula_count
                > expected_thresholds["maximum_same_formula_connectivities"]
            ):
                raise ValueError("Qualified source-paired switch fails its entity gates")
            if selected_identity == baseline_identity or selected_identity not in truly_qualified:
                raise ValueError("Qualified source-paired switch lacks a qualified challenger")
            if consensus.get("ranking_changed_from_v11") is not True:
                raise ValueError("Qualified source-paired switch must declare a rank change")
        elif consensus.get("ranking_changed_from_v11") is not False:
            raise ValueError("Unqualified source-paired consensus cannot change ranking")
        indexed[entity_id] = record
    if set(indexed) != entity_ids:
        raise ValueError("Candidate-recall bundle must cover exactly the V11 entities")
    if int(metadata.get("entity_count") or -1) != len(indexed):
        raise ValueError("Candidate-recall entity count metadata mismatch")
    return indexed, _sha256(bundle), metadata


def _optional_view_record(
    bundle: Mapping[str, Any] | None,
    *,
    entity_id: str,
    identity_key: str,
    view_name: str,
) -> Mapping[str, Any] | None:
    if not isinstance(bundle, Mapping):
        return None
    entity = (bundle.get("entities") or {}).get(entity_id)
    if not isinstance(entity, Mapping):
        return None
    candidate = (entity.get("candidates") or {}).get(identity_key)
    if not isinstance(candidate, Mapping):
        return None
    record = candidate.get(view_name)
    return record if isinstance(record, Mapping) else None


def _optional_view(
    bundle: Mapping[str, Any] | None,
    *,
    entity_id: str,
    identity_key: str,
    view_name: str,
    source_hash: str,
    policy: Mapping[str, Any],
) -> tuple[float | None, dict[str, Any]]:
    record = _optional_view_record(
        bundle,
        entity_id=entity_id,
        identity_key=identity_key,
        view_name=view_name,
    )
    audit = {"view": view_name, "available": record is not None, "eligible": False}
    if record is None:
        audit["reason"] = "missing_view_is_neutral"
        return None, audit
    if str(record.get("source_semantic_sha256") or "") != source_hash:
        audit["reason"] = "source_binding_mismatch"
        return None, audit
    if record.get("eligible") is not True:
        audit["reason"] = str(record.get("reason") or "view_not_eligible")
        return None, audit
    contract = ((policy.get("multi_view") or {}).get("optional_view_contracts") or {}).get(view_name) or {}
    if contract.get("candidate_conditioned") and record.get("candidate_conditioned") is not True:
        audit["reason"] = "candidate_conditioning_semantics_missing"
        return None, audit
    if contract.get("require_instrument_domain_match") and record.get("instrument_domain_match") is not True:
        audit["reason"] = "instrument_domain_mismatch"
        return None, audit
    minimum_architectures = int(contract.get("minimum_independent_architecture_families") or 0)
    if int(record.get("independent_architecture_families") or 0) < minimum_architectures:
        audit["reason"] = "insufficient_independent_architecture_families"
        return None, audit
    score = _finite(record.get("score"))
    if score is None or not 0.0 <= score <= 1.0:
        audit["reason"] = "invalid_view_score"
        return None, audit
    audit.update({"eligible": True, "reason": "qualified", "score": round(score, 6)})
    return score, audit


def candidate_feature_vector(
    candidate: Mapping[str, Any],
    *,
    entity_id: str,
    source_hash: str,
    view_bundle: Mapping[str, Any] | None,
    policy: Mapping[str, Any],
) -> tuple[dict[str, float | None], list[dict[str, Any]], set[str]]:
    components = candidate.get("score_components") or {}
    ce_records = candidate.get("collision_energy_evidence") or []
    ce_supports = [
        value
        for record in ce_records
        if isinstance(record, Mapping)
        and (value := _finite(record.get("support"))) is not None
    ]
    stability = 1.0 - min(1.0, statistics.pstdev(ce_supports)) if len(ce_supports) > 1 else 1.0
    chemistry_status = str((candidate.get("chemistry_gate") or {}).get("status") or "")
    display = candidate.get("candidate") or {}
    features: dict[str, float | None] = {
        "v11_entity_score": _finite(candidate.get("entity_score")),
        "ce_support_mean": _finite(components.get("mean_collision_energy_support")),
        "ce_support_median": _finite(components.get("median_collision_energy_support")),
        "ce_coverage": _finite(components.get("collision_energy_coverage")),
        "ce_top1_fraction": _finite(components.get("collision_energy_top1_fraction")),
        "ce_support_floor": min(ce_supports) if ce_supports else None,
        "ce_support_range": max(ce_supports) - min(ce_supports) if ce_supports else None,
        "ce_support_stability": stability if ce_supports else None,
        "reference_source_diversity": _finite(components.get("reference_source_diversity")),
        "chemistry_consistency": 1.0 if chemistry_status == "consistent" else 0.0 if chemistry_status == "contradicted" else None,
        "candidate_structure_available": 1.0 if str(display.get("smiles") or "").strip() else 0.0,
    }
    families = {"measured_spectral"}
    if ce_supports:
        families.add("collision_energy_consistency")
    if int(candidate.get("reference_source_family_count") or 0) > 0:
        families.add("reference_provenance")
    if chemistry_status in {"consistent", "contradicted"}:
        families.add("candidate_independent_chemistry")
    audits: list[dict[str, Any]] = []
    identity_key = str(candidate.get("identity_key") or "")
    for view_name in OPTIONAL_FEATURES:
        value, audit = _optional_view(
            view_bundle,
            entity_id=entity_id,
            identity_key=identity_key,
            view_name=view_name,
            source_hash=source_hash,
            policy=policy,
        )
        features[view_name] = value
        audits.append(audit)
        if value is not None:
            families.add(view_name)
    return features, audits, families


def _model_utility(features: Mapping[str, float | None], model: Mapping[str, Any]) -> float:
    utility = float(model.get("intercept") or 0.0)
    means = model["feature_mean"]
    scales = model["feature_scale"]
    coefficients = model["coefficients"]
    for feature in FEATURE_ORDER:
        value = features.get(feature)
        # Missing means exactly the training mean, giving a zero standardized
        # contribution.  It therefore cannot become negative evidence.
        standardized = 0.0 if value is None else (
            float(value) - float(means[feature])
        ) / float(scales[feature])
        utility += float(coefficients[feature]) * standardized
    if not math.isfinite(utility):
        raise ValueError("Residual model produced a non-finite utility")
    return utility


def _shadow_decision(
    entity: Mapping[str, Any],
    candidates: Sequence[Mapping[str, Any]],
    policy: Mapping[str, Any],
    source_paired_consensus: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    valid = [
        candidate
        for candidate in candidates
        if not candidate.get("vetoed")
        and candidate.get("shadow_selection_eligible", True) is not False
    ]
    decision = policy.get("decision") or {}
    minimum = float(decision.get("minimum_entity_score") or 0.50)
    minimum_margin = float(decision.get("minimum_connectivity_margin") or 0.04)
    minimum_ce = int(decision.get("minimum_distinct_collision_energies_for_connectivity") or 2)
    minimum_supported_ce = int(decision.get("minimum_winner_collision_energy_support") or 2)
    isomer_margin = float(decision.get("isomer_group_margin") or 0.05)
    if (
        isinstance(source_paired_consensus, Mapping)
        and source_paired_consensus.get("status") == "qualified_switch"
    ):
        selected_identity = str(source_paired_consensus.get("selected_identity_key") or "")
        comparison = next(
            (
                record
                for record in (source_paired_consensus.get("comparisons") or [])
                if isinstance(record, Mapping)
                and str(record.get("challenger_identity_key") or "") == selected_identity
            ),
            {},
        )
        return {
            "status": "selected_candidate_only",
            "scope": "source_paired_same_formula_connectivity_candidate",
            "selection_basis": "qualified_source_paired_ce_consensus",
            "winner_margin": None,
            "v11_top1_margin": source_paired_consensus.get("baseline_margin"),
            "source_paired_median_score_delta": comparison.get("median_score_delta"),
            "source_paired_observation_win_fraction": comparison.get(
                "observation_win_fraction"
            ),
            "source_paired_provenance_stratum_win_fraction": comparison.get(
                "provenance_stratum_win_fraction"
            ),
            "source_paired_unique_provenance_profile_count": comparison.get(
                "unique_provenance_evidence_profile_count"
            ),
            "source_paired_unique_provenance_profile_win_fraction": comparison.get(
                "unique_provenance_profile_win_fraction"
            ),
            "minimum_score": minimum,
            "minimum_margin": minimum_margin,
            "blockers": ["same_formula_connectivity_requires_orthogonal_confirmation"],
            "stereochemistry_claimed": False,
            "score_is_probability": False,
            "production_effect": False,
        }
    leader = valid[0] if valid else None
    runner = valid[1] if len(valid) > 1 else None
    leader_score = float((leader or {}).get("evirank_mv_score") or 0.0)
    margin = leader_score - float((runner or {}).get("evirank_mv_score") or 0.0) if runner else leader_score
    if leader is None:
        status = "unmatched" if not candidates else "abstained"
        scope = "none"
        blockers = (
            ["no_library_candidate"]
            if not candidates
            else ["no_qualified_shadow_selection_candidate"]
        )
    elif leader_score < minimum:
        status, scope, blockers = "abstained", "none", ["entity_score_below_minimum"]
    else:
        formula = normalize_formula_text((leader.get("candidate") or {}).get("formula"))
        same_formula = [
            candidate for candidate in valid[1:]
            if formula
            and normalize_formula_text((candidate.get("candidate") or {}).get("formula")) == formula
            and leader_score - float(candidate.get("evirank_mv_score") or 0.0) <= isomer_margin
        ]
        total_ce = int(entity.get("distinct_collision_energy_bins") or 0)
        supported_ce = int(leader.get("collision_energy_bins_supported") or 0)
        if same_formula:
            status, scope, blockers = (
                "selected_isomer_group",
                "molecular_formula_isomer_group",
                ["competitive_same_formula_connectivity"],
            )
        elif total_ce >= minimum_ce and supported_ce >= minimum_supported_ce and margin >= minimum_margin:
            status, scope, blockers = "selected", "2d_connectivity_candidate", []
        elif total_ce < minimum_ce:
            status, scope, blockers = (
                "selected_candidate_only",
                "single_collision_energy_candidate",
                ["insufficient_distinct_collision_energies"],
            )
        else:
            status, scope, blockers = (
                "selected_candidate_only",
                "ambiguous_connectivity_candidate",
                ["insufficient_entity_margin_or_energy_coverage"],
            )
    return {
        "status": status,
        "scope": scope,
        "winner_margin": round(margin, 6),
        "minimum_score": minimum,
        "minimum_margin": minimum_margin,
        "blockers": blockers,
        "stereochemistry_claimed": False,
        "score_is_probability": False,
        "production_effect": False,
    }


def _rank_entity(
    entity: Mapping[str, Any],
    *,
    source_hash: str,
    view_bundle: Mapping[str, Any] | None,
    candidate_recall: Mapping[str, Any] | None,
    model: Mapping[str, Any] | None,
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    entity_id = str(entity.get("entity_id") or "")
    base_candidates = [
        copy.deepcopy(candidate)
        for candidate in (entity.get("candidates") or [])
        if isinstance(candidate, Mapping)
    ]
    prepared: list[dict[str, Any]] = []
    utilities: list[float] = []
    minimum_families = int(
        ((policy.get("multi_view") or {}).get("minimum_eligible_view_families_for_residual") or 2)
    )
    for base_rank, candidate in enumerate(base_candidates, start=1):
        features, view_audit, families = candidate_feature_vector(
            candidate,
            entity_id=entity_id,
            source_hash=source_hash,
            view_bundle=view_bundle,
            policy=policy,
        )
        utility = _model_utility(features, model) if model is not None else 0.0
        eligible = bool(model is not None and len(families) >= minimum_families and not candidate.get("vetoed"))
        prepared.append({
            **candidate,
            "candidate_origin": "v11",
            "enters_shadow_candidate_union": True,
            "v11_entity_rank": base_rank,
            "v11_entity_score": round(float(candidate.get("entity_score") or 0.0), 6),
            "evirank_mv_features": {
                key: (round(value, 6) if value is not None else None)
                for key, value in features.items()
            },
            "evirank_mv_view_audit": view_audit,
            "eligible_view_families": sorted(families),
            "residual_eligible": eligible,
            "shadow_selection_eligible": not bool(candidate.get("vetoed")),
            "_utility": utility,
        })
        if eligible:
            utilities.append(utility)
    center = statistics.median(utilities) if utilities else 0.0
    maximum = float((model or {}).get("maximum_absolute_delta") or 0.0)
    temperature = float((model or {}).get("temperature") or 1.0)
    for candidate in prepared:
        delta = 0.0
        if candidate["residual_eligible"]:
            delta = maximum * math.tanh((float(candidate["_utility"]) - center) / temperature)
        baseline = float(candidate["v11_entity_score"])
        candidate["evirank_mv_residual"] = round(delta, 6)
        candidate["evirank_mv_score"] = round(_clip(baseline + delta), 6)
        candidate["score_is_probability"] = False
        candidate["production_effect"] = False
        candidate.pop("_utility", None)
    prepared.sort(
        key=lambda candidate: (
            bool(candidate.get("vetoed")),
            -float(candidate.get("evirank_mv_score") or 0.0),
            int(candidate.get("v11_entity_rank") or 999999),
            str(candidate.get("identity_key") or ""),
        )
    )
    for rank, candidate in enumerate(prepared, start=1):
        candidate["evirank_mv_rank"] = rank

    v11_winner = entity.get("winner") or entity.get("leading_candidate_when_abstained")
    v11_identity = str((v11_winner or {}).get("identity_key") or "") or None
    recall_record = copy.deepcopy(candidate_recall) if candidate_recall is not None else None
    expanded_count = int((recall_record or {}).get("expanded_candidate_count") or 0)
    source_paired_consensus = copy.deepcopy(
        (recall_record or {}).get("source_paired_consensus")
    )
    qualified_switch = bool(
        isinstance(source_paired_consensus, Mapping)
        and source_paired_consensus.get("status") == "qualified_switch"
    )
    selected_identity = (
        str(source_paired_consensus.get("selected_identity_key") or "")
        if isinstance(source_paired_consensus, Mapping)
        else ""
    )

    base_identity_keys = {
        str(candidate.get("identity_key") or "") for candidate in prepared
    }
    expanded: list[dict[str, Any]] = []
    for candidate in (recall_record or {}).get("expanded_candidates") or []:
        if not isinstance(candidate, Mapping):
            continue
        identity_key = str(candidate.get("identity_key") or "")
        if identity_key in base_identity_keys:
            continue
        selection_eligible = bool(qualified_switch and identity_key == selected_identity)
        expanded.append({
            **copy.deepcopy(dict(candidate)),
            "candidate_origin": "evirecall_expansion",
            "enters_shadow_candidate_union": True,
            "enters_residual_model_scoring": False,
            "v11_entity_rank": None,
            "v11_entity_score": None,
            "evirank_mv_features": None,
            "evirank_mv_view_audit": [],
            "eligible_view_families": [],
            "residual_eligible": False,
            "evirank_mv_residual": None,
            "evirank_mv_score": None,
            "evirank_mv_rank": None,
            "shadow_selection_eligible": selection_eligible,
            "source_paired_selection_qualified": selection_eligible,
            "score_is_probability": False,
            "production_effect": False,
        })
    expanded.sort(
        key=lambda candidate: (
            int(candidate.get("retrieval_rank") or 999999),
            str(candidate.get("identity_key") or ""),
        )
    )

    for candidate in prepared:
        candidate["source_paired_selection_qualified"] = bool(
            qualified_switch
            and str(candidate.get("identity_key") or "") == selected_identity
        )
    candidate_union = [*prepared, *expanded]
    if qualified_switch:
        selected = next(
            (
                candidate
                for candidate in candidate_union
                if str(candidate.get("identity_key") or "") == selected_identity
                and candidate.get("shadow_selection_eligible") is not False
            ),
            None,
        )
        if selected is None:
            raise ValueError("Qualified source-paired selection is absent from the shadow union")
        candidate_union = [
            selected,
            *(
                candidate
                for candidate in candidate_union
                if candidate is not selected
            ),
        ]
    for rank, candidate in enumerate(candidate_union, start=1):
        candidate["shadow_candidate_union_rank"] = rank
    winner = next(
        (
            candidate
            for candidate in candidate_union
            if not candidate.get("vetoed")
            and candidate.get("shadow_selection_eligible", True) is not False
        ),
        None,
    )
    winner_identity = str((winner or {}).get("identity_key") or "") or None
    return {
        "entity_id": entity_id,
        "entity_source": entity.get("entity_source"),
        "entity_source_value": entity.get("entity_source_value"),
        "representative_precursor_mz": entity.get("representative_precursor_mz"),
        "representative_retention_time": entity.get("representative_retention_time"),
        "member_scan_count": entity.get("member_scan_count"),
        "member_result_indices": entity.get("member_result_indices") or [],
        "distinct_collision_energy_bins": entity.get("distinct_collision_energy_bins"),
        "v11_decision": copy.deepcopy(entity.get("decision") or {}),
        "v11_winner_identity_key": v11_identity,
        "decision": _shadow_decision(
            entity,
            candidate_union,
            policy,
            source_paired_consensus=source_paired_consensus,
        ),
        "winner": winner,
        "candidates": candidate_union,
        "candidate_count": len(candidate_union),
        "v11_candidate_count": len(prepared),
        "ranking_changed_from_v11": bool(winner_identity and v11_identity and winner_identity != v11_identity),
        "source_paired_consensus": source_paired_consensus,
        "source_paired_switch_applied": qualified_switch,
        "candidate_recall_expanded": expanded_count > 0,
        "candidate_recall_expanded_candidate_count": expanded_count,
        "candidate_recall_expansion": recall_record,
        "expanded_candidates_entered_shadow_candidate_union": bool(expanded),
        "expanded_candidates_entered_residual_model_scoring": False,
        "expanded_candidates_entered_evirank_scoring": False,
        "production_v11_artifact_mutated": False,
        "production_effect": False,
    }


def build_evirank_mv_result(
    v11_payload: Mapping[str, Any],
    *,
    detection_method_id: str = "unspecified",
    view_bundle: Mapping[str, Any] | None = None,
    candidate_recall_bundle: Mapping[str, Any] | None = None,
    residual_model: Mapping[str, Any] | None = None,
    policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one deterministic V12 shadow result from an immutable V11 artifact."""

    selected_policy = (
        load_evirank_policy()
        if policy is None
        else _validate_evirank_policy(policy)
    )
    v11_metadata, source_hash = _validate_v11_payload(v11_payload, selected_policy)
    model, model_hash = _validate_model(residual_model, selected_policy)
    method_id = str(detection_method_id or "").strip() or "unspecified"
    source_entities = [
        entity
        for entity in (v11_payload.get("entities") or [])
        if isinstance(entity, Mapping)
    ]
    entity_ids = {str(entity.get("entity_id") or "") for entity in source_entities}
    if "" in entity_ids or len(entity_ids) != len(source_entities):
        raise ValueError("V11 entity ids must be unique and non-empty")
    recall_index, recall_hash, recall_metadata = _validate_candidate_recall_bundle(
        candidate_recall_bundle,
        v11_payload=v11_payload,
        source_hash=source_hash,
        detection_method_id=method_id,
        entity_ids=entity_ids,
        policy=selected_policy,
    )
    if isinstance(view_bundle, Mapping):
        declared_source = str((view_bundle.get("metadata") or {}).get("source_semantic_sha256") or "")
        if declared_source and declared_source != source_hash:
            raise ValueError("Multi-view bundle is bound to another V11 source artifact")
    entities = [
        _rank_entity(
            entity,
            source_hash=source_hash,
            view_bundle=view_bundle,
            candidate_recall=recall_index.get(str(entity.get("entity_id") or "")),
            model=model,
            policy=selected_policy,
        )
        for entity in source_entities
    ]
    decisions = Counter(str((entity.get("decision") or {}).get("status") or "unknown") for entity in entities)
    return {
        "metadata": {
            "schema_version": EVIRANK_SCHEMA_VERSION,
            "module_id": MODULE_ID,
            "algorithm_id": selected_policy["algorithm_id"],
            "algorithm_name": selected_policy["display_name"],
            "deployment_state": "shadow",
            "task_id": v11_metadata.get("task_id"),
            "detection_method_id": method_id,
            "method_context_available": method_id != "unspecified",
            "base_algorithm_id": v11_metadata.get("algorithm_id"),
            "base_policy_sha256": v11_metadata.get("policy_sha256"),
            "base_source_semantic_sha256": source_hash,
            "base_entity_semantic_sha256": v11_metadata.get("entity_semantic_sha256"),
            "base_artifact_sha256": _sha256(v11_payload),
            "policy_sha256": evirank_policy_sha256(selected_policy),
            "view_bundle_sha256": _sha256(view_bundle) if view_bundle is not None else None,
            "candidate_recall_bundle_sha256": recall_hash,
            "candidate_recall_algorithm_id": (recall_metadata or {}).get("algorithm_id"),
            "candidate_recall_policy_sha256": (recall_metadata or {}).get("policy_sha256"),
            "candidate_recall_library_snapshot_id": (recall_metadata or {}).get("library_snapshot_id"),
            "candidate_recall_state": (
                "attached_shadow_union_source_paired_gated"
                if recall_metadata is not None
                else "unavailable_no_candidate_recall_bundle"
            ),
            "entities_with_candidate_recall_expansion": sum(
                bool(entity.get("candidate_recall_expanded")) for entity in entities
            ),
            "entities_with_expanded_candidates_in_shadow_union": sum(
                bool(entity.get("expanded_candidates_entered_shadow_candidate_union"))
                for entity in entities
            ),
            "entities_with_qualified_source_paired_switch": sum(
                bool(entity.get("source_paired_switch_applied")) for entity in entities
            ),
            "source_paired_consensus_algorithm_id": (
                (selected_policy.get("candidate_recall") or {}).get(
                    "source_paired_consensus_algorithm_id"
                )
            ),
            "expanded_candidates_entered_shadow_candidate_union": any(
                bool(entity.get("expanded_candidates_entered_shadow_candidate_union"))
                for entity in entities
            ),
            "expanded_candidates_entered_residual_model_scoring": False,
            "expanded_candidates_entered_evirank_scoring": False,
            "residual_model_id": (model or {}).get("model_id"),
            "residual_model_sha256": model_hash,
            "residual_model_state": (model or {}).get("state") or "unavailable_no_frozen_model",
            "entity_count": len(entities),
            "entities_reranked_from_v11": sum(bool(entity.get("ranking_changed_from_v11")) for entity in entities),
            "decision_status_counts": dict(sorted(decisions.items())),
            "score_is_probability": False,
            "production_identity_selector": False,
            "production_ranking_changed": False,
            "legacy_result_mutated": False,
            "v11_artifact_mutated": False,
        },
        "policy": selected_policy,
        "invariants": {
            "shadow_only": True,
            "v11_production_result_immutable": True,
            "missing_views_are_negative_evidence": False,
            "score_is_probability": False,
            "same_identity_different_methods_retained": True,
            "technical_repeats_are_independent_identity_votes": False,
            "method_cases_are_independent_identity_votes": False,
            "ground_truth_available_at_inference": False,
            "compound_name_affects_ranking": False,
            "expert_score_affects_ranking": False,
            "classification_affects_ranking": False,
            "biosynthetic_pathway_affects_ranking": False,
            "candidate_pathway_affects_ranking": False,
            "candidate_recall_names_affect_ranking": False,
            "candidate_recall_expanded_candidates_enter_shadow_candidate_union": True,
            "candidate_recall_expanded_candidates_enter_residual_model_scoring": False,
            "expanded_candidate_top1_requires_qualified_source_paired_consensus": True,
            "source_paired_consensus_ground_truth_available_at_inference": False,
            "source_paired_provenance_strata_are_independent_sources": False,
            "same_connectivity_different_detection_methods_preexcluded": False,
            "stereochemistry_claimed": False,
        },
        "entities": entities,
    }


def build_method_case_index(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Index validation/training records without deleting method variants.

    Connectivity identity is ground truth here and is strictly a cohort-design
    field; it is never passed to :func:`build_evirank_mv_result`.
    """

    seen_records: set[str] = set()
    cases: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for raw in records:
        record = dict(raw)
        record_id = str(record.get("record_id") or "").strip()
        identity = str(record.get("connectivity_identity") or "").strip()
        method = str(record.get("detection_method_id") or "").strip()
        role = str(record.get("role") or "primary").strip()
        if not record_id or not identity or not method:
            raise ValueError("Every method record needs record_id, connectivity_identity and detection_method_id")
        if record_id in seen_records:
            raise ValueError(f"Duplicate method record_id: {record_id}")
        if role not in {"primary", "technical_repeat"}:
            raise ValueError(f"Unsupported method record role: {role}")
        seen_records.add(record_id)
        record["role"] = role
        cases[(identity, method)].append(record)
    method_cases = []
    for (identity, method), members in sorted(cases.items()):
        primaries = [member for member in members if member["role"] == "primary"]
        if len(primaries) != 1:
            raise ValueError(
                "Each connectivity identity x detection method case requires exactly one primary"
            )
        method_cases.append({
            "method_case_id": hashlib.sha256(f"{identity}|{method}".encode("utf-8")).hexdigest()[:20],
            "connectivity_identity": identity,
            "detection_method_id": method,
            "primary_record_id": primaries[0]["record_id"],
            "technical_repeat_record_ids": sorted(
                member["record_id"] for member in members if member["role"] == "technical_repeat"
            ),
            "all_record_ids": sorted(member["record_id"] for member in members),
            "independent_identity_vote_count": 1,
            "method_case_vote_count": 1,
        })
    identities = sorted({identity for identity, _ in cases})
    return {
        "records_retained": len(records),
        "connectivity_identity_count": len(identities),
        "method_case_count": len(method_cases),
        "technical_repeat_count": sum(
            len(case["technical_repeat_record_ids"]) for case in method_cases
        ),
        "same_identity_different_methods_retained": True,
        "technical_repeats_increase_independent_identity_n": False,
        "method_cases_increase_independent_identity_n": False,
        "method_cases": method_cases,
    }


__all__ = [
    "EVIRANK_SCHEMA_VERSION",
    "FEATURE_ORDER",
    "MODULE_ID",
    "build_evirank_mv_result",
    "build_method_case_index",
    "candidate_feature_vector",
    "default_evirank_policy_path",
    "evirank_policy_sha256",
    "load_evirank_policy",
]
