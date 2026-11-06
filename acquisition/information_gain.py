"""Probability-free maximin utility for active evidence acquisition actions."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from typing import Any

from acquisition.observable_catalog import (
    evidiscrim_pair_energy_separation,
    pair_fragment_separation,
)
from schemas.acquisition import AcquisitionAction, AcquisitionNeed


def _clip(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def _action_id(need_id: str, proposal: Mapping[str, Any]) -> str:
    payload = json.dumps(
        [need_id, proposal.get("action_kind"), proposal.get("collision_energies"),
         proposal.get("collision_energy_design"), proposal.get("isolation_width_da"),
         proposal.get("isolation_width_profile_id"),
         proposal.get("max_injection_time_ms"),
         proposal.get("resource_identity_keys"),
         proposal.get("resource_declaration_id")],
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return "action_" + hashlib.sha256(payload).hexdigest()[:16]


def _pair_separation(
    need: AcquisitionNeed,
    proposal: Mapping[str, Any],
    catalog: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, float], tuple[float, ...], bool]:
    kind = str(proposal.get("action_kind") or "")
    resource_identity_keys = set(proposal.get("resource_identity_keys") or ())
    scores: dict[str, float] = {}
    diagnostic: set[float] = set()
    evidiscrim_used = False
    for edge in need.ambiguity_edges:
        fragment_score, ions = pair_fragment_separation(
            edge.left.identity_key, edge.right.identity_key, catalog
        )
        if kind == "authentic_standard_rt":
            diagnostic.update(ions)
            covered = sum(
                identity_key in resource_identity_keys
                for identity_key in (
                    edge.left.identity_key, edge.right.identity_key
                )
            )
            separation = 0.95 if covered == 2 else 0.78 if covered == 1 else 0.0
        elif kind == "ion_mobility_ccs":
            diagnostic.update(ions)
            covered = sum(
                identity_key in resource_identity_keys
                for identity_key in (
                    edge.left.identity_key, edge.right.identity_key
                )
            )
            if covered == 2:
                separation = (
                    0.62
                    if "same_formula_connectivity_competition" in edge.reasons
                    else 0.45
                )
            elif covered == 1:
                separation = 0.28
            else:
                separation = 0.0
        elif kind == "targeted_ms2_multi_ce":
            predicted = evidiscrim_pair_energy_separation(
                edge.edge_id, proposal, catalog
            )
            if predicted is not None:
                separation, predicted_ions = predicted
                diagnostic.update(predicted_ions)
                evidiscrim_used = True
            else:
                diagnostic.update(ions)
                separation = min(0.92, fragment_score + 0.18)
        elif kind == "narrow_isolation_ms2":
            diagnostic.update(ions)
            separation = min(0.88, fragment_score + 0.12)
        elif kind == "high_sensitivity_targeted_ms2":
            diagnostic.update(ions)
            separation = min(0.82, fragment_score + 0.08)
        elif kind == "opposite_polarity_targeted_ms2":
            diagnostic.update(ions)
            separation = min(0.75, fragment_score + 0.05)
        else:
            diagnostic.update(ions)
            separation = fragment_score
        scores[edge.edge_id] = round(_clip(separation), 4)
    return scores, tuple(sorted(diagnostic)[:12]), evidiscrim_used


def score_action(
    need: AcquisitionNeed,
    row: Mapping[str, Any],
    proposal: Mapping[str, Any],
    catalog: Mapping[str, Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> AcquisitionAction:
    """Score an action by worst-pair separability, not a posterior probability."""

    kind = str(proposal.get("action_kind") or "")
    pairwise, diagnostic_ions, evidiscrim_used = _pair_separation(
        need, proposal, catalog
    )
    minimum_separation = min(pairwise.values(), default=0.0)
    defaults = policy.get("utility_defaults") or {}
    peaks_count = int(row.get("peaks_count") or len(row.get("query_peaks") or []))
    measurement_success = _clip(
        float(defaults.get("measurement_success_base", 0.72))
        + min(peaks_count, 20) / 100.0
    )
    if kind == "high_sensitivity_targeted_ms2":
        measurement_success = _clip(measurement_success + 0.12)
    if kind == "authentic_standard_rt":
        measurement_success = 0.90
    transferability = {
        "targeted_ms2_multi_ce": 0.72,
        "narrow_isolation_ms2": 0.78,
        "high_sensitivity_targeted_ms2": 0.70,
        "authentic_standard_rt": 0.96,
        "ion_mobility_ccs": 0.70,
        "opposite_polarity_targeted_ms2": 0.48,
    }.get(kind, 0.55)
    energy_design = proposal.get("collision_energy_design") or {}
    numeric_schedule = bool(energy_design.get("numeric_schedule_available", True))
    energy_evidence_status = str(energy_design.get("evidence_status") or "")
    if kind in {
        "targeted_ms2_multi_ce",
        "narrow_isolation_ms2",
        "opposite_polarity_targeted_ms2",
    } and not numeric_schedule:
        transferability = max(0.20, transferability - 0.18)
    if str(proposal.get("collision_energy_unit") or "").endswith("unspecified"):
        transferability = max(0.25, transferability - 0.12)
    if kind == "narrow_isolation_ms2" and proposal.get("isolation_width_da") is None:
        transferability = max(0.20, transferability - 0.20)
    isolation_status = str(
        ((row.get("match") or {}).get("chimeric_deconvolution") or {}).get(
            "isolation_status"
        ) or ""
    )
    purity = {"clean": 0.95, "mixed": 0.62, "severe": 0.38}.get(
        isolation_status, 0.72
    )
    if kind == "narrow_isolation_ms2":
        purity = max(purity, 0.88)
    if kind in {"authentic_standard_rt", "ion_mobility_ccs"}:
        purity = 0.95
    costs = policy.get("duty_cycle_costs") or {}
    risks = policy.get("risk_penalties") or {}
    duty_cost = _clip(float(costs.get(kind, 0.30)))
    risk = _clip(float(risks.get(kind, 0.20)))
    weights = policy.get("utility_weights") or {}
    cost_weight = float(weights.get("duty_cycle_cost", 0.18))
    risk_weight = float(weights.get("risk", 0.12))
    raw_utility = (
        minimum_separation * measurement_success * transferability * purity
        - cost_weight * duty_cost
        - risk_weight * risk
    )
    utility = _clip(raw_utility)
    uncertainty = 0.06
    if not diagnostic_ions and kind not in {"authentic_standard_rt", "ion_mobility_ccs"}:
        uncertainty += 0.09
    if energy_evidence_status in {
        "single_energy_only_no_response_curve",
        "no_observed_energy_response_curve",
    }:
        uncertainty += 0.04
    if kind in {
        "targeted_ms2_multi_ce",
        "narrow_isolation_ms2",
        "opposite_polarity_targeted_ms2",
    } and not numeric_schedule:
        uncertainty += 0.08
    if str(proposal.get("collision_energy_unit") or "").endswith("unspecified"):
        uncertainty += 0.05
    if kind == "narrow_isolation_ms2" and proposal.get("isolation_width_da") is None:
        uncertainty += 0.08
    if kind in {"authentic_standard_rt", "ion_mobility_ccs"}:
        covered_candidates = set(proposal.get("resource_identity_keys") or ())
        if len(covered_candidates) < len(need.candidates):
            uncertainty += 0.06
    interval = (_clip(utility - uncertainty), _clip(utility + uncertainty))
    rationale = tuple(proposal.get("rationale") or ()) + (
        "效用是有界优先级，不是身份概率、置信区间或 conformal coverage。",
        f"最难区分候选对的预期可分辨度为 {minimum_separation:.2f}。",
    )
    if evidiscrim_used:
        rationale += (
            "EviDiscrim-MS 的不确定性校正候选对谱分布差异仅用于碰撞能量规划；预测谱不是独立确认，也未改变当前排序。",
        )
    polarity = str(proposal.get("polarity_override") or need.polarity or "").strip() or None
    return AcquisitionAction(
        action_id=_action_id(need.need_id, proposal),
        action_kind=kind,
        label=str(proposal.get("label") or kind),
        precursor_mz=need.precursor_mz,
        retention_time_min=need.retention_time_min,
        rt_window_min=float(policy.get("default_rt_window_min", 0.15)),
        polarity=polarity,
        collision_energies=tuple(
            float(value) for value in proposal.get("collision_energies") or ()
        ),
        collision_energy_unit=proposal.get("collision_energy_unit"),
        collision_energy_design=proposal.get("collision_energy_design"),
        isolation_width_da=proposal.get("isolation_width_da"),
        isolation_width_profile_id=proposal.get("isolation_width_profile_id"),
        isolation_width_evidence_status=proposal.get(
            "isolation_width_evidence_status"
        ),
        max_injection_time_ms=proposal.get("max_injection_time_ms"),
        diagnostic_ions=diagnostic_ions,
        pairwise_separation=pairwise,
        minimum_pair_separation=round(minimum_separation, 4),
        measurement_success_factor=round(measurement_success, 4),
        transferability_factor=round(transferability, 4),
        purity_factor=round(purity, 4),
        duty_cycle_cost=round(duty_cost, 4),
        risk_penalty=round(risk, 4),
        utility=round(utility, 4),
        utility_interval=(round(interval[0], 4), round(interval[1], 4)),
        rationale=rationale,
        feasibility=str(proposal.get("feasibility") or "manual_review_required"),
        required_resources=tuple(proposal.get("required_resources") or ()),
        resource_status=str(proposal.get("resource_status") or "not_required"),
        resource_identity_keys=tuple(
            proposal.get("resource_identity_keys") or ()
        ),
        resource_declaration_id=proposal.get("resource_declaration_id"),
    )
