"""Atomic artifact contract for the V15 production entity decisions."""

from __future__ import annotations

import hashlib
import json
from runtime.json_encoding import write_json_encoded
import os
import uuid
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

from matching.v15_flagship_ranker import (
    load_v15_flagship_policy,
)
from pipeline.spectral_tie_ranker import (
    build_v15_flagship_result,
    artifact_implementation_is_current,
)
from pipeline.spectral_tie_scope import apply_spectral_tie_scope_guard
from matching.flagship_ranker import flagship_policy_sha256
from matching.evirank_retrieval import build_candidate_recall_bundle
from matching.v15_candidate_recall import attach_v15_candidate_recall
from matching.v15_decision_assurance import attach_v15_decision_assurance
from matching.v15_identity_resolution import apply_v15_identity_resolution
from matching.v15_ce_consensus import (
    apply_v15_ce_evidence,
    build_v15_ce_evidence_bundle,
)
from pipeline.flagship import read_flagship_artifact
from pipeline.v16_candidate_missingness import attach_replicate_missingness_audit


from pipeline.active_artifact_binding import active_implementation_sha256
from pipeline.derived_artifact_paths import (
    derived_artifact_path, binding_metadata, assert_derived_destination,
    metadata_binding_is_current,
)


V15_FLAGSHIP_ARTIFACT_SUFFIX = "_v15_flagship.json"
V16_FLAGSHIP_ARTIFACT_SUFFIX = "_v16_flagship.json"


def v15_flagship_artifact_path(task_dir: Path, task_id: str) -> Path:
    # Preserve saved V15 rankings. The active V16 rule materializes a separate
    # derived artifact; it never silently upgrades the old winner in place.
    policy = load_v15_flagship_policy()
    suffix = (V16_FLAGSHIP_ARTIFACT_SUFFIX
              if str(policy.get("algorithm_id") or "").startswith("mssystem-v16-")
              else V15_FLAGSHIP_ARTIFACT_SUFFIX)
    return derived_artifact_path(task_dir, task_id, suffix, policy=policy,
        actual_implementation_sha256=active_implementation_sha256(policy))


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _entity_semantic_sha256(payload: Mapping[str, Any]) -> str:
    digest = hashlib.sha256()
    for entity in payload.get("entities") or []:
        digest.update(_canonical(entity))
        digest.update(b"\n")
    return digest.hexdigest()


def _write_payload_atomic(payload: Mapping[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(
        f".{output_path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    )
    with temporary.open("w", encoding="utf-8") as handle:
        write_json_encoded(
            payload,
            handle,
            ensure_ascii=False,
            separators=(",", ":"),
            default=str,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output_path)


@lru_cache(maxsize=256)
def _read_artifact_metadata_cached(
    path_text: str, size_bytes: int, mtime_ns: int
) -> dict[str, Any]:
    """Read only metadata; stat arguments make the cache generation-safe."""

    path = Path(path_text)
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return dict(payload.get("metadata") or {})
    with path.open("rb") as handle:
        return dict(next(ijson.items(handle, "metadata", use_float=True), {}) or {})


def artifact_is_current(path: Path) -> bool:
    """Return whether an artifact uses the active V15 policy and resolver."""

    if not path.is_file():
        return False
    try:
        stat = path.stat()
        metadata = _read_artifact_metadata_cached(
            str(path.resolve()), int(stat.st_size), int(stat.st_mtime_ns)
        )
        policy = load_v15_flagship_policy()
        return bool(
            metadata.get("algorithm_id") == policy.get("algorithm_id")
            and metadata.get("policy_sha256") == flagship_policy_sha256(policy)
            and metadata.get("identity_resolution_status") == "available"
            and artifact_implementation_is_current(metadata, policy)
            and metadata_binding_is_current(metadata, policy,
                actual_implementation_sha256=active_implementation_sha256(policy))
        )
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return False


def upgrade_v15_flagship_identity_scope_artifact(
    path: Path,
    *,
    task_id: str | None = None,
    runtime_snapshot_id: str | None = None,
) -> dict[str, Any]:
    """Upgrade a compatible V15 artifact without repeating library search.

    The 1.4.x line changes only the identity-scope resolver.  A saved artifact
    may be upgraded in place when every pre-existing ranking policy field is
    byte-semantically equal to the active policy and its stored entity digest
    still verifies.  Candidate order, scores and winner selection are then
    guarded again by :func:`apply_v15_identity_resolution`.

    Any older artifact with another ranking policy, a corrupt digest, or an
    unexpected task identity raises and must be rebuilt from immutable scan
    evidence through the full production path.
    """

    payload = json.loads(path.read_text(encoding="utf-8"))
    metadata = payload.get("metadata") or {}
    if metadata.get("production_identity_selector") is not True:
        raise ValueError("Saved artifact is not a production V15 identity selector")
    if str(metadata.get("module_id") or "") != (
        "v15_conditioned_multilibrary_flagship"
    ):
        raise ValueError("Saved artifact is not a conditioned V15 flagship result")
    saved_task_id = str(metadata.get("task_id") or "")
    if task_id is not None and saved_task_id != str(task_id):
        raise ValueError("Saved V15 artifact belongs to another task")
    active_policy = load_v15_flagship_policy()
    if "artifact_namespace" in active_policy:
        if not metadata_binding_is_current(metadata, active_policy,
                actual_implementation_sha256=active_implementation_sha256(active_policy)):
            raise ValueError("Versioned artifact requires a new bound materialization")
        suffix = (V16_FLAGSHIP_ARTIFACT_SUFFIX if str(active_policy.get("algorithm_id") or "").startswith("mssystem-v16-") else V15_FLAGSHIP_ARTIFACT_SUFFIX)
        assert_derived_destination(path, saved_task_id, suffix, policy=active_policy,
            actual_implementation_sha256=active_implementation_sha256(active_policy))
    saved_policy = payload.get("policy") or {}
    if (saved_policy.get("spectral_tie_handling") or {}).get("enabled") is True:
        if not artifact_implementation_is_current(metadata, saved_policy):
            raise ValueError("Saved spectral tie implementation changed; rebuild from scan evidence")
    if bool(
        metadata.get("algorithm_id") == active_policy.get("algorithm_id")
        and metadata.get("policy_sha256")
        == flagship_policy_sha256(active_policy)
        and metadata.get("identity_resolution_status") == "available"
        and artifact_implementation_is_current(metadata, active_policy)
    ):
        return {
            **metadata,
            "status": "available",
            "upgrade_status": "already_current",
            "artifact": str(path.resolve()),
            "size_bytes": path.stat().st_size,
        }

    saved_policy = payload.get("policy")
    if not isinstance(saved_policy, Mapping):
        raise ValueError("Saved V15 artifact has no embedded policy")
    if metadata.get("algorithm_id") != saved_policy.get("algorithm_id"):
        raise ValueError("Saved V15 artifact algorithm and policy disagree")
    if metadata.get("policy_sha256") != flagship_policy_sha256(saved_policy):
        raise ValueError("Saved V15 policy digest does not verify")

    saved_ranking_policy = dict(saved_policy)
    active_ranking_policy = dict(active_policy)
    saved_ranking_policy.pop("algorithm_id", None)
    active_ranking_policy.pop("algorithm_id", None)
    saved_ranking_policy.pop("identity_resolution", None)
    active_ranking_policy.pop("identity_resolution", None)
    if _canonical(saved_ranking_policy) != _canonical(active_ranking_policy):
        raise ValueError(
            "Saved V15 ranking policy differs from the active ranking policy"
        )
    saved_resolution_policy = saved_policy.get("identity_resolution") or {}
    saved_resolution_algorithm = str(
        saved_resolution_policy.get("algorithm_id") or ""
    )
    if saved_resolution_policy.get("enabled") and saved_resolution_algorithm not in {
        "mssystem-v15-identity-scope-resolution-1.0.0",
        "mssystem-v15-identity-scope-resolution-1.1.0",
    }:
        raise ValueError(
            "Saved artifact has an unsupported identity-resolution generation"
        )

    saved_entity_digest = str(metadata.get("entity_semantic_sha256") or "")
    if not saved_entity_digest:
        raise ValueError("Saved V15 artifact has no entity semantic digest")
    if _entity_semantic_sha256(payload) != saved_entity_digest:
        raise ValueError("Saved V15 entity semantic digest does not verify")

    source_algorithm_id = str(metadata.get("algorithm_id") or "")
    source_policy_sha256 = str(metadata.get("policy_sha256") or "")
    if saved_resolution_policy.get("enabled"):
        resolver_blockers = {
            "competitive_connectivity_formula_assignment_incomplete",
            "competitive_connectivities_have_conflicting_formulae",
            "competitive_same_formula_connectivities",
        }
        resolver_decision_fields = {
            "reportable_identity_key",
            "reportable_formula",
            "ranking_lead_identity_key",
            "ranking_lead_is_reportable_single_connectivity",
            "identity_resolution_module_id",
            "identity_resolution_algorithm_id",
        }
        for entity in payload.get("entities") or []:
            if not isinstance(entity, dict):
                continue
            resolution = entity.pop("identity_resolution", None)
            if not isinstance(resolution, Mapping):
                raise ValueError(
                    "Saved resolved entity has no verifiable upstream decision"
                )
            if str(resolution.get("algorithm_id") or "") != (
                saved_resolution_algorithm
            ):
                raise ValueError(
                    "Saved entity resolution generation differs from its policy"
                )
            decision = dict(entity.get("decision") or {})
            decision["status"] = resolution.get("upstream_decision_status")
            decision["scope"] = resolution.get("upstream_decision_scope")
            decision["blockers"] = [
                blocker
                for blocker in decision.get("blockers") or []
                if blocker not in resolver_blockers
            ]
            for field in resolver_decision_fields:
                decision.pop(field, None)
            entity["decision"] = decision
            selected = entity.get("winner") or entity.get(
                "leading_candidate_when_abstained"
            )
            if isinstance(selected, dict):
                selected.pop("reporting_semantics", None)
    payload["policy"] = active_policy
    payload = apply_v15_identity_resolution(payload)
    payload = apply_spectral_tie_scope_guard(payload)
    payload = attach_v15_decision_assurance(payload)
    metadata = payload.setdefault("metadata", {})
    metadata.update(
        {
            "algorithm_id": active_policy["algorithm_id"],
            "policy_sha256": flagship_policy_sha256(active_policy),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "runtime_snapshot_id": (
                runtime_snapshot_id
                if runtime_snapshot_id is not None
                else metadata.get("runtime_snapshot_id")
            ),
            "identity_resolution_upgrade_status": "verified_policy_extension",
            "identity_resolution_upgrade_mode": (
                "existing_v15_shortlist_without_library_retrieval"
            ),
            "identity_resolution_upgrade_source_algorithm_id": (
                source_algorithm_id
            ),
            "identity_resolution_upgrade_source_policy_sha256": (
                source_policy_sha256
            ),
            "identity_resolution_upgrade_source_entity_semantic_sha256": (
                saved_entity_digest
            ),
            "identity_resolution_upgrade_replaced_algorithm_id": (
                saved_resolution_algorithm or None
            ),
            "identity_resolution_upgrade_changes_candidate_order": False,
            "identity_resolution_upgrade_changes_candidate_score": False,
            "identity_resolution_upgrade_changes_winner": False,
        }
    )
    metadata["entity_semantic_sha256"] = _entity_semantic_sha256(payload)
    _write_payload_atomic(payload, path)
    return {
        **metadata,
        "status": "available",
        "upgrade_status": "upgraded_without_library_retrieval",
        "artifact": str(path.resolve()),
        "size_bytes": path.stat().st_size,
    }


def write_v15_flagship_artifact(
    task_id: str,
    legacy_results: Iterable[Mapping[str, Any]],
    output_path: Path,
    *,
    runtime_snapshot_id: str | None = None,
    ce_evidence_retriever: Any | None = None,
    candidate_recall_retriever: Any | None = None,
    detection_method_id: str = "unspecified",
) -> dict[str, Any]:
    """Build V15 without modifying scan results or the V11 audit artifact."""

    policy = load_v15_flagship_policy()
    implementation_digest = active_implementation_sha256(policy)
    suffix = (V16_FLAGSHIP_ARTIFACT_SUFFIX if str(policy.get("algorithm_id") or "").startswith("mssystem-v16-") else V15_FLAGSHIP_ARTIFACT_SUFFIX)
    assert_derived_destination(output_path, task_id, suffix, policy=policy,
        actual_implementation_sha256=implementation_digest)
    source_rows = [dict(row) for row in legacy_results if isinstance(row, Mapping)]
    base_payload = build_v15_flagship_result(task_id, source_rows, policy=policy)
    base_payload = attach_replicate_missingness_audit(base_payload, source_rows)
    payload = base_payload
    candidate_recall_bundle = None
    if candidate_recall_retriever is not None:
        candidate_recall_bundle = build_candidate_recall_bundle(
            base_payload,
            source_rows,
            retriever=candidate_recall_retriever,
            detection_method_id=detection_method_id,
        )
    if ce_evidence_retriever is not None:
        ce_evidence = build_v15_ce_evidence_bundle(
            base_payload,
            source_rows,
            retriever=ce_evidence_retriever,
            candidate_recall_bundle=candidate_recall_bundle,
            detection_method_id=detection_method_id,
        )
        payload = apply_v15_ce_evidence(base_payload, ce_evidence)
        payload["metadata"]["source_paired_ce_evidence_status"] = "evaluated"
    else:
        payload["metadata"].update(
            {
                "source_paired_ce_evidence_status": "unavailable",
                "source_paired_ce_gate_evaluated_entities": 0,
                "entities_reranked_by_source_paired_ce_gate": 0,
                "source_paired_ce_gate_production_effect": False,
            }
        )
    if candidate_recall_bundle is not None:
        payload = attach_v15_candidate_recall(
            payload,
            candidate_recall_bundle,
            base_payload=base_payload,
        )
    else:
        payload["metadata"].update(
            {
                "candidate_recall_support_status": "unavailable",
                "entities_with_candidate_recall_expansion": 0,
                "expanded_candidate_hypothesis_count": 0,
                "candidate_recall_support_affects_primary_ranking": False,
            }
        )
    payload = apply_v15_identity_resolution(payload)
    payload = apply_spectral_tie_scope_guard(payload)
    payload = attach_v15_decision_assurance(payload)
    metadata = payload["metadata"]
    metadata.update(binding_metadata(policy,
        actual_implementation_sha256=active_implementation_sha256(policy)))
    metadata.update(
        {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "runtime_snapshot_id": runtime_snapshot_id,
            "entity_semantic_sha256": _entity_semantic_sha256(payload),
            "artifact_role": (
                "primary_v16_multilibrary_with_missingness_audit"
                if str(metadata.get("algorithm_id") or "").startswith("mssystem-v16-")
                else "primary_v15_conditioned_multilibrary_and_source_paired_ce_identification"
            ),
            "legacy_scan_result_role": "candidate_generation_and_raw_evidence",
            "v11_v14_score_lane_role": "embedded_read_only_audit_baseline",
        }
    )
    _write_payload_atomic(payload, output_path)
    return {
        **metadata,
        "status": "available",
        "artifact": str(output_path.resolve()),
        "size_bytes": output_path.stat().st_size,
    }


__all__ = [
    "V15_FLAGSHIP_ARTIFACT_SUFFIX",
    "artifact_is_current",
    "read_flagship_artifact",
    "upgrade_v15_flagship_identity_scope_artifact",
    "v15_flagship_artifact_path",
    "write_v15_flagship_artifact",
]
