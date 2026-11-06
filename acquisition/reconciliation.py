"""Candidate-independent linkage of follow-up observations to acquisition plans."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

from acquisition.ambiguity_graph import candidate_identity_key
from schemas.acquisition import EvidenceDelta


def _polarity(value: Any) -> str:
    text = str(value or "").strip().upper()
    if text.startswith("POS") or text == "+":
        return "POS"
    if text.startswith("NEG") or text == "-":
        return "NEG"
    return text


def link_followup_observations(
    plan: Mapping[str, Any],
    followup_task_id: str,
    followup_results: Sequence[Mapping[str, Any]],
    *,
    mz_tolerance_ppm: float = 10.0,
    rt_tolerance_min: float = 0.20,
) -> EvidenceDelta:
    need = plan.get("need") or {}
    target_mz = float(need.get("precursor_mz"))
    target_rt = need.get("retention_time_min")
    target_rt = float(target_rt) if target_rt is not None else None
    target_polarity = _polarity(need.get("polarity"))
    linked_indices: list[int] = []
    # Linkage is completed before reading candidate identities.  Only physical
    # observation coordinates participate in this loop.
    for index, row in enumerate(followup_results):
        try:
            observed_mz = float(row.get("precursor_mz"))
        except (TypeError, ValueError):
            continue
        ppm = abs(observed_mz - target_mz) / target_mz * 1e6
        if ppm > mz_tolerance_ppm:
            continue
        observed_polarity = _polarity(row.get("ion_mode") or row.get("polarity"))
        if target_polarity and observed_polarity and target_polarity != observed_polarity:
            continue
        if target_rt is not None:
            try:
                observed_rt = float(row.get("retention_time"))
            except (TypeError, ValueError):
                continue
            if abs(observed_rt - target_rt) > rt_tolerance_min:
                continue
        linked_indices.append(index)

    baseline_keys = tuple(sorted({
        str(candidate.get("identity_key") or "")
        for candidate in need.get("candidates") or []
        if str(candidate.get("identity_key") or "")
    }))
    followup_keys: set[str] = set()
    for index in linked_indices:
        match = followup_results[index].get("match")
        if not isinstance(match, Mapping):
            continue
        followup_keys.add(candidate_identity_key(match))
        for candidate in match.get("top_candidates") or []:
            if isinstance(candidate, Mapping):
                followup_keys.add(candidate_identity_key(candidate))
    overlap = tuple(sorted(set(baseline_keys).intersection(followup_keys)))
    if not linked_indices:
        change = "no_candidate_independent_followup_observation"
    elif len(overlap) == 1:
        change = "candidate_set_narrowed_to_one_observed_identity"
    elif overlap:
        change = "ambiguity_persists_after_followup"
    elif followup_keys:
        change = "followup_candidate_set_conflicts_with_baseline"
    else:
        change = "followup_observed_without_library_candidate"
    digest = hashlib.sha256(json.dumps(
        [plan.get("plan_id"), followup_task_id, linked_indices],
        separators=(",", ":"),
    ).encode("utf-8")).hexdigest()[:16]
    return EvidenceDelta(
        delta_id=f"delta_{digest}",
        plan_id=str(plan.get("plan_id") or ""),
        baseline_task_id=str(plan.get("task_id") or ""),
        followup_task_id=followup_task_id,
        linked_result_indices=tuple(linked_indices),
        linkage_fields_used=("polarity", "precursor_mz", "retention_time_min"),
        candidate_fields_used_for_linkage=(),
        baseline_candidate_keys=baseline_keys,
        followup_candidate_keys=tuple(sorted(followup_keys)),
        resolved_candidate_keys=overlap if len(overlap) == 1 else (),
        evidence_change=change,
        observed_at=datetime.now(timezone.utc).isoformat(),
    )
