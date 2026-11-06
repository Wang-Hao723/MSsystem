"""Generate feasible follow-up actions from missing evidence and instrument context."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from acquisition.energy_design import (
    build_collision_energy_recommendation,
    classify_instrument_domain,
)
from acquisition.resources import (
    reference_keys_for_need,
    resolve_acquisition_resources,
)
from schemas.acquisition import AcquisitionNeed


def _number(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed


def _isolation_width_recommendation(
    row: Mapping[str, Any],
    policy: Mapping[str, Any],
    context: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Keep isolation widths inside a verified instrument-native domain."""

    context = context or {}
    vendor = context.get("vendor_conversion") or {}
    source = vendor.get("source") or {}
    acquisition = source.get("acquisition") or {}
    params = context.get("params") or {}
    observed = None
    for value in (
        row.get("isolation_window_width_da"),
        params.get("isolation_window_width_da"),
        context.get("recommended_isolation_window_width_da"),
        acquisition.get("isolation_width_da"),
    ):
        parsed = _number(value)
        if parsed is not None and parsed > 0:
            observed = parsed
            break
    domain, _ = classify_instrument_domain(row, context)
    profiles = policy.get("isolation_width_design") or {}
    if domain == "agilent_qtof":
        profile = profiles.get("agilent_qtof") or {}
        routine = observed or _number(profile.get("medium_isolation_width_da"))
        narrow = _number(profile.get("narrow_isolation_width_da"))
        return {
            "routine_width_da": routine,
            "narrow_width_da": narrow,
            "profile_id": str(
                profile.get("profile_id") or "agilent_qtof_native_isolation"
            ),
            "evidence_status": (
                "observed_method_width_plus_vendor_native_narrow"
                if observed is not None
                else "vendor_native_medium_and_narrow"
            ),
            "rationale": (
                "Agilent Q-TOF uses the instrument-native Narrow setting (~1.3 amu), not a 0.4 Da Orbitrap-style value.",
                "The routine width preserves the observed method setting when available.",
            ),
        }
    return {
        "routine_width_da": observed,
        "narrow_width_da": None,
        "profile_id": "method_owner_defined_native_isolation",
        "evidence_status": (
            "observed_method_width_only" if observed is not None
            else "instrument_domain_unqualified_no_numeric_width"
        ),
        "rationale": (
            "No numeric narrow isolation width is emitted outside a qualified instrument domain.",
        ),
    }


def candidate_actions(
    need: AcquisitionNeed,
    row: Mapping[str, Any],
    policy: Mapping[str, Any],
    *,
    context: Mapping[str, Any] | None = None,
    enabled_resource_actions: frozenset[str] | None = None,
    resolved_resources: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Enumerate feasible advisory actions; unavailable resources are omitted."""

    defaults = policy.get("action_defaults") or {}
    resources = dict(
        resolved_resources or resolve_acquisition_resources(context, policy)
    )
    isolation = _isolation_width_recommendation(row, policy, context)
    if enabled_resource_actions is None:
        enabled_resource_actions = frozenset(
            action_kind
            for action_kind, budget_key in (
                ("authentic_standard_rt", "authentic_standard_action_budget"),
                ("ion_mobility_ccs", "ion_mobility_action_budget"),
            )
            if int(resources.get(budget_key) or 0) > 0
        )
    energy = build_collision_energy_recommendation(
        row, policy, context=context
    )
    energies = tuple(energy["initial_energies"])
    ce_unit = str(energy["collision_energy_unit"])
    energy_design = energy["design"]
    energy_rationale = tuple(energy["rationale"])
    proposals: list[dict[str, Any]] = [{
        "action_kind": "targeted_ms2_multi_ce",
        "label": "能量分辨靶向 MS/MS（序贯侦察）",
        "collision_energies": energies,
        "collision_energy_unit": ce_unit,
        "collision_energy_design": energy_design,
        "isolation_width_da": isolation["routine_width_da"],
        "isolation_width_profile_id": isolation["profile_id"],
        "isolation_width_evidence_status": isolation["evidence_status"],
        "max_injection_time_ms": float(defaults.get("routine_injection_time_ms", 100.0)),
        "feasibility": "instrument_method_review_required",
        "rationale": (
            "先测量前体存活、产物离子丰富度和跨能量新颖度，再决定是否补低、中或高能量。",
        ) + energy_rationale + tuple(isolation["rationale"]),
    }]
    missing = set(need.missing_evidence)
    deconvolution = (row.get("match") or {}).get("chimeric_deconvolution") or {}
    if "narrow_isolation_confirmation" in missing or str(
        deconvolution.get("isolation_status") or ""
    ) in {"mixed", "severe"}:
        proposals.append({
            "action_kind": "narrow_isolation_ms2",
            "label": "窄隔离窗能量分辨复采 MS/MS",
            "collision_energies": energies,
            "collision_energy_unit": ce_unit,
            "collision_energy_design": energy_design,
            "isolation_width_da": isolation["narrow_width_da"],
            "isolation_width_profile_id": isolation["profile_id"],
            "isolation_width_evidence_status": isolation["evidence_status"],
            "max_injection_time_ms": float(defaults.get("narrow_injection_time_ms", 150.0)),
            "feasibility": "verify_precursor_density_and_instrument_limit",
            "rationale": (
                "当前隔离窗存在混合/严重共隔离，先降低嵌合风险。",
                "窄窗可能损失离子流，需与注入时间联动。",
            ) + tuple(isolation["rationale"]) + energy_rationale,
        })
    peaks_count = int(_number(row.get("peaks_count")) or len(row.get("query_peaks") or []))
    if peaks_count < int(defaults.get("weak_fragment_peak_count", 8)):
        proposals.append({
            "action_kind": "high_sensitivity_targeted_ms2",
            "label": "高灵敏长注入靶向 MS/MS（能量侦察后）",
            "collision_energies": (),
            "collision_energy_unit": ce_unit,
            "isolation_width_da": isolation["routine_width_da"],
            "isolation_width_profile_id": isolation["profile_id"],
            "isolation_width_evidence_status": isolation["evidence_status"],
            "max_injection_time_ms": float(defaults.get("long_injection_time_ms", 250.0)),
            "feasibility": "run_after_energy_scout_then_check_peak_width_and_cycle_time",
            "rationale": (
                "当前碎片峰较少；先完成能量侦察，再在最低且仍有诊断信息的能量延长注入。",
                "长注入动作不预填固定能量，避免把低信号问题与欠碎裂问题混为一谈。",
            ),
        })
    declaration_id = str(resources.get("effective_declaration_id") or "") or None
    standard_keys = reference_keys_for_need(
        need, resources, "authentic_standard_rt"
    )
    if (
        "authentic_standard_rt" in enabled_resource_actions
        and standard_keys
    ):
        proposals.append({
            "action_kind": "authentic_standard_rt",
            "label": "同方法标品 RT 共洗脱",
            "collision_energies": (),
            "collision_energy_unit": None,
            "isolation_width_da": None,
            "max_injection_time_ms": None,
            "feasibility": "explicit_authentic_standard_inventory_confirmed",
            "required_resources": ("authentic_standard",),
            "resource_status": "explicitly_available",
            "resource_identity_keys": standard_keys,
            "resource_declaration_id": declaration_id,
            "rationale": (
                "仅对资源清单中明确登记的候选标品安排同方法保留时间/共洗脱验证。",
                "未登记的候选不假定有标品；预测 RT 不能替代实测共洗脱。",
            ),
        })
    ccs_keys = reference_keys_for_need(need, resources, "ion_mobility_ccs")
    if (
        "ion_mobility_ccs" in enabled_resource_actions
        and bool(resources.get("ion_mobility_available"))
        and ccs_keys
    ):
        proposals.append({
            "action_kind": "ion_mobility_ccs",
            "label": "离子淌度/CCS 正交测量",
            "collision_energies": (),
            "collision_energy_unit": None,
            "isolation_width_da": None,
            "max_injection_time_ms": None,
            "feasibility": "explicit_ion_mobility_and_ccs_reference_confirmed",
            "required_resources": ("ion_mobility_instrument", "ccs_reference"),
            "resource_status": "explicitly_available",
            "resource_identity_keys": ccs_keys,
            "resource_declaration_id": declaration_id,
            "rationale": (
                "仪器能力与候选特异 CCS 参比均已显式登记；按实际覆盖候选对评价。",
                "没有候选参比的 CCS 测量不能被包装成身份确认。",
            ),
        })
    if bool(policy.get("allow_opposite_polarity", True)):
        opposite = "NEG" if str(need.polarity or "").upper().startswith("POS") else "POS"
        opposite_energy = build_collision_energy_recommendation(
            row,
            policy,
            context=context,
            ignore_observed_energy=True,
        )
        proposals.append({
            "action_kind": "opposite_polarity_targeted_ms2",
            "label": f"{opposite} 反极性能量分辨靶向 MS/MS",
            "collision_energies": tuple(opposite_energy["initial_energies"]),
            "collision_energy_unit": opposite_energy["collision_energy_unit"],
            "collision_energy_design": opposite_energy["design"],
            "isolation_width_da": isolation["routine_width_da"],
            "isolation_width_profile_id": isolation["profile_id"],
            "isolation_width_evidence_status": isolation["evidence_status"],
            "max_injection_time_ms": float(defaults.get("routine_injection_time_ms", 100.0)),
            "polarity_override": opposite,
            "feasibility": "requires_opposite_polarity_detectability",
            "rationale": (
                "候选在另一极性可能产生互补碎裂，但不可把未检出当成身份反证。",
                "反极性必须重新测量自己的能量响应，不能沿用当前极性的最佳能量。",
            ) + tuple(isolation["rationale"]) + tuple(opposite_energy["rationale"]),
        })
    return proposals
