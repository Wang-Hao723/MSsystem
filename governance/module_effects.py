"""Explain the observable effect of each analytical module for one task.

The effect inspector is deliberately descriptive.  It consumes persisted task
summaries and never reranks candidates, recalibrates confidence, or promotes a
claim.  Its purpose is to make non-ranking improvements visible instead of
letting users infer that an unchanged Top-1 table means a module did not run.
"""

from __future__ import annotations

from typing import Any


SCHEMA_VERSION = "1.0.0"


def _number(value: Any, default: int | float = 0) -> int | float:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)):
        return value
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return int(parsed) if parsed.is_integer() else parsed


def _effect(
    module_id: str,
    display_name: str,
    *,
    state: str,
    requested: bool,
    applicable: bool,
    executed: bool,
    headline: str,
    metrics: dict[str, Any],
    ranking_effect: str,
    interpretation: str,
) -> dict[str, Any]:
    return {
        "module_id": module_id,
        "display_name": display_name,
        "state": state,
        "requested": requested,
        "applicable": applicable,
        "executed": executed,
        "headline": headline,
        "metrics": metrics,
        "ranking_effect": ranking_effect,
        "interpretation": interpretation,
    }


def build_module_effect_summary(task: dict[str, Any]) -> dict[str, Any]:
    """Build a browser-ready, audit-only module effect summary."""

    availability = task.get("module_availability") or {}
    confidence = task.get("compound_confidence_summary") or {}
    funnel = task.get("identification_funnel") or {}
    deconvolution = task.get("deconvolution_impact") or {}
    ms1 = availability.get("ms1_feature_detection") or {}
    iimn = availability.get("iimn") or {}
    expected = task.get("expected_compounds_summary") or {}
    identification = task.get("identification_contract_summary") or {}
    claims = task.get("claim_contract_summary") or {}
    active_evidence = task.get("active_evidence_summary") or {}
    evograph = task.get("evograph_fragmentation_summary") or {}
    evidiscrim = task.get("evidiscrim_summary") or {}
    shadow_states = identification.get("module_state_counts") or {}

    top1_changed = bool(deconvolution.get("top1_changed"))
    production_ranking_changed = int(
        _number(identification.get("production_ranking_changed_records"), 0)
    )
    high_delta = (deconvolution.get("delta") or {}).get(
        "high_supported_unique_compounds", 0
    )
    effects = [
        _effect(
            "core_spectral_ranker",
            "核心实验谱库检索",
            state="production",
            requested=True,
            applicable=int(_number(task.get("total_spectra"), 0)) > 0,
            executed=str(task.get("status") or "") == "completed",
            headline=(
                f"{int(_number(task.get('matched_spectra'), 0))}/"
                f"{int(_number(task.get('total_spectra'), 0))} 张 MS2 有候选，"
                f"{int(_number(task.get('unique_compounds'), 0))} 个去重身份"
            ),
            metrics={
                "matched_ms2": int(_number(task.get("matched_spectra"), 0)),
                "total_ms2": int(_number(task.get("total_spectra"), 0)),
                "unique_library_candidates": int(
                    _number(task.get("unique_compounds"), 0)
                ),
                "high_supported_identities": int(
                    _number(confidence.get("high_supported_unique_compounds"), 0)
                ),
                "unmatched_with_fragments": int(
                    _number(funnel.get("unmatched_with_nonzero_query_peaks"), 0)
                ),
            },
            ranking_effect="production_ranker",
            interpretation="这是普通 MS2 候选表的来源；其他证据层不会自动并入该计数。",
        ),
        _effect(
            "ms1_feature_detection",
            "连续 MS1 特征检测与 MS2 归属",
            state="operational",
            requested=bool(ms1.get("requested")),
            applicable=bool(ms1.get("applicable")),
            executed=bool(ms1.get("executed_with_ms2_links")),
            headline=(
                f"{int(_number(ms1.get('lcms_features_with_ms2'), 0))} 个 LC–MS 特征"
                f"关联 {int(_number(ms1.get('linked_ms2_spectra'), 0))} 张 MS2"
            ),
            metrics={
                "ms1_scans": int(_number(availability.get("ms1_scan_count"), 0)),
                "lcms_features_with_ms2": int(
                    _number(ms1.get("lcms_features_with_ms2"), 0)
                ),
                "linked_ms2_spectra": int(
                    _number(ms1.get("linked_ms2_spectra"), 0)
                ),
            },
            ranking_effect="none",
            interpretation="增加色谱特征、峰界与共识谱上下文，不等同于新增化合物身份。",
        ),
        _effect(
            "chimeric_deconvolution",
            "嵌合谱受限去卷积",
            state="qualified",
            requested=bool((availability.get("chimeric_deconvolution") or {}).get("requested")),
            applicable=bool((availability.get("chimeric_deconvolution") or {}).get("applicable")),
            executed=bool((availability.get("chimeric_deconvolution") or {}).get("executed")),
            headline=(
                f"{int(_number(deconvolution.get('accepted_score_increase_spectra'), 0))} 张通过守卫；"
                f"High 支持身份净增 {int(_number(high_delta, 0)):+d}"
            ),
            metrics={
                "accepted_spectra": int(
                    _number(deconvolution.get("accepted_score_increase_spectra"), 0)
                ),
                "high_supported_identity_delta": int(_number(high_delta, 0)),
                "representative_high_delta": int(_number(
                    (deconvolution.get("delta") or {}).get(
                        "representative_high_unique_compounds"
                    ), 0
                )),
                "high_spectra_delta": int(_number(
                    (deconvolution.get("delta") or {}).get("high_spectra"), 0
                )),
                "top1_changed": top1_changed,
            },
            ranking_effect=("invariant_violation" if top1_changed else "selected_identity_preserved"),
            interpretation="只允许在已选身份内增加受守卫分数；Top-1 不变是设计约束，不是未生效。",
        ),
        _effect(
            "ion_identity_networking",
            "IIMN 离子身份网络",
            state="operational_audit",
            requested=bool(iimn.get("requested")),
            applicable=bool(iimn.get("applicable")),
            executed=bool(iimn.get("executed")),
            headline=(
                f"{int(_number(task.get('ion_identity_networks'), 0))} 个网络，"
                f"{int(_number(task.get('ion_identity_linked_annotations'), 0))} 个链接注释"
            ),
            metrics={
                "networks": int(_number(task.get("ion_identity_networks"), 0)),
                "networked_features": int(
                    _number(task.get("ion_identity_networked_features"), 0)
                ),
                "linked_annotations": int(
                    _number(task.get("ion_identity_linked_annotations"), 0)
                ),
            },
            ranking_effect="none",
            interpretation="合并共洗脱加合物/源内碎片上下文；链接离子不冒充独立 High。",
        ),
        _effect(
            "expected_compounds",
            "Expected Compounds 连续 MS1 靶向层",
            state="independent_evidence_layer",
            requested=bool(expected.get("requested")),
            applicable=str(expected.get("status") or "") == "ok",
            executed=str(expected.get("status") or "") == "ok",
            headline=(
                f"{int(_number(expected.get('detected_feature_count'), 0))} 个原始 XIC 峰，"
                f"{int(_number(expected.get('supported_features'), 0))} 个有正交支持"
            ),
            metrics={
                "detected_features": int(
                    _number(expected.get("detected_feature_count"), 0)
                ),
                "supported_features": int(
                    _number(expected.get("supported_features"), 0)
                ),
                "mass_only_features": int(
                    _number(expected.get("mass_only_features"), 0)
                ),
                "ordinary_ms2_match_increment": int(_number(
                    expected.get("standard_spectral_library_match_increment"), 0
                )),
            },
            ranking_effect="separate_from_ordinary_ms2",
            interpretation="它回答文献靶标在 MS1 中是否有色谱证据，明确不改普通 MS2 命中数。",
        ),
        _effect(
            "evograph_fragmentation",
            "EvoGraph-F 结构化中性丢失图",
            state="shadow",
            requested=str(evograph.get("status") or "") not in {"", "not_started"},
            applicable=int(_number(task.get("total_spectra"), 0)) > 0,
            executed=str(evograph.get("status") or "") == "completed",
            headline=(
                f"{int(_number(evograph.get('spectra_processed'), 0))} 张谱构图，"
                f"{int(_number(evograph.get('heuristically_supported_edges'), 0))} 条启发式支持边"
                + ("（达到任务影子上限）" if evograph.get("truncated_by_task_limit") else "")
            ),
            metrics={
                "spectra_processed": int(
                    _number(evograph.get("spectra_processed"), 0)
                ),
                "truncated_by_task_limit": bool(
                    evograph.get("truncated_by_task_limit")
                ),
                "formula_context_available": int(
                    _number(evograph.get("formula_context_available"), 0)
                ),
                "nodes": int(_number(evograph.get("nodes"), 0)),
                "selected_edges": int(
                    _number(evograph.get("selected_edges"), 0)
                ),
                "heuristically_supported_edges": int(
                    _number(evograph.get("heuristically_supported_edges"), 0)
                ),
                "ambiguous_edges": int(
                    _number(evograph.get("ambiguous_edges"), 0)
                ),
                "unassigned_edges": int(
                    _number(evograph.get("unassigned_edges"), 0)
                ),
                "cycle_count": int(_number(evograph.get("cycle_count"), 0)),
                "high_promotions": int(
                    _number(evograph.get("high_promotions"), 0)
                ),
            },
            ranking_effect="forbidden",
            interpretation=(
                "边权是未校准的局部启发式量，不是概率；外部标准校准前不生成 "
                "Bayes factor、不改 Top-1/置信度，也不触发 High。"
            ),
        ),
        _effect(
            "evidiscrim_candidate_contrast",
            "EviDiscrim-MS 候选对比诊断证据",
            state="shadow",
            requested=bool(evidiscrim.get("requested")),
            applicable=bool(evidiscrim.get("applicable")),
            executed=bool(evidiscrim.get("executed")),
            headline=(
                f"{int(_number(evidiscrim.get('resolved_counterfactual_pairs'), 0))}/"
                f"{int(_number(evidiscrim.get('eligible_connectivity_pairs'), 0))} 个连接体候选对形成可判别反事实，"
                f"{int(_number(evidiscrim.get('collision_energy_recommendations'), 0))} 个 CE 规划建议"
            ),
            metrics={
                "spectra_evaluated": int(
                    _number(evidiscrim.get("spectra_evaluated"), 0)
                ),
                "eligible_connectivity_pairs": int(
                    _number(evidiscrim.get("eligible_connectivity_pairs"), 0)
                ),
                "resolved_counterfactual_pairs": int(
                    _number(evidiscrim.get("resolved_counterfactual_pairs"), 0)
                ),
                "abstained_eligible_pairs": int(
                    _number(evidiscrim.get("abstained_eligible_pairs"), 0)
                ),
                "same_connectivity_out_of_scope_pairs": int(_number(
                    evidiscrim.get("same_connectivity_out_of_scope_pairs"), 0
                )),
                "collision_energy_recommendations": int(_number(
                    evidiscrim.get("collision_energy_recommendations"), 0
                )),
                "instrument_domain_shifted_energy_audits": int(_number(
                    evidiscrim.get(
                        "instrument_domain_shifted_energy_audits"
                    ), 0
                )),
                "counterfactual_order_changed_spectra": int(_number(
                    evidiscrim.get(
                        "counterfactual_candidate_order_changed_spectra"
                    ), 0
                )),
            },
            ranking_effect="forbidden",
            interpretation=(
                "只比较同分子式、不同二维连接体候选的预测诊断峰；模型方差会降权，"
                "同连接体立体异构体和域外 CE 会弃权。它不改当前 Top-1 或置信度。"
            ),
        ),
        _effect(
            "shadow_evidence_modules",
            "影子证据与反事实模块",
            state="shadow",
            requested=True,
            applicable=bool(shadow_states),
            executed=bool(shadow_states),
            headline=f"{len(shadow_states)} 类模块留下逐记录状态",
            metrics={
                "module_state_counts": shadow_states,
                "identification_records": int(
                    _number(identification.get("record_count"), 0)
                ),
                "production_ranking_changed_records": production_ranking_changed,
            },
            ranking_effect="forbidden",
            interpretation="CE、对比专家、RECAP 等仅显示/审计；影子态本来就必须保持排序不变。",
        ),
        _effect(
            "evigate_v2_claim_compiler",
            "EviGate 2 最终声称编译",
            state="qualified",
            requested=True,
            applicable=bool(claims.get("record_count")),
            executed=bool(claims.get("record_count")),
            headline=(
                f"编译 {int(_number(claims.get('record_count'), 0))} 条，"
                f"释放 {int(_number(claims.get('released_count'), 0))} 条"
            ),
            metrics={
                "compiled": int(_number(claims.get("record_count"), 0)),
                "released": int(_number(claims.get("released_count"), 0)),
                "calibrated": int(_number(claims.get("calibrated_count"), 0)),
                "production_ranking_changed_records": int(
                    _number(claims.get("production_ranking_changed_records"), 0)
                ),
            },
            ranking_effect="none",
            interpretation="最终声称出口可以保守阻断/弃权，但不回写旧 Top-1 或旧置信字段。",
        ),
        _effect(
            "eviacquire_active_evidence",
            "EviAcquire-MS 主动证据采集",
            state="shadow",
            requested=bool(active_evidence.get("requested")),
            applicable=bool(active_evidence.get("applicable")),
            executed=int(_number(active_evidence.get("plan_count"), 0)) > 0,
            headline=(
                f"已生成 {int(_number(active_evidence.get('plan_count'), 0))} 个候选无关采集计划"
            ),
            metrics={
                "plan_count": int(_number(active_evidence.get("plan_count"), 0)),
                "planned_actions": int(
                    _number(active_evidence.get("planned_actions"), 0)
                ),
                "followups_linked": int(
                    _number(active_evidence.get("followups_linked"), 0)
                ),
            },
            ranking_effect="forbidden",
            interpretation="只建议下一次采什么，不把预测证据当成已观测证据，也不改变当前结果。",
        ),
    ]

    observable_effects = [
        effect["module_id"] for effect in effects
        if effect["executed"] and effect["module_id"] != "core_spectral_ranker"
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "task_id": str(task.get("id") or ""),
        "core_ranking_changed": bool(top1_changed or production_ranking_changed),
        "production_ranking_changed_records": production_ranking_changed,
        "observable_nonranking_effects": observable_effects,
        "why_core_looks_unchanged": [
            "影子模块按治理契约禁止改变检索、Top-1、置信度和释放。",
            "IIMN、连续 MS1 特征和 Expected Compounds 产生独立证据层，不并入普通 MS2 命中数。",
            "嵌合去卷积只在已选身份内受限增分；其 Top-1 不变是强制不变量。",
            "EviGate 与 EviAcquire 分别管理当前声称和下一次采集，不回写历史排序。",
        ],
        "effects": effects,
        "legacy_results_mutated": False,
    }
