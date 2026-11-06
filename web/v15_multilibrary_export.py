"""Stable long-table projection of V15 multi-library candidate evidence."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "mssystem-v15-multilibrary-export/1.0"

MULTILIBRARY_COLUMNS = [
    "Entity_ID",
    "Decision_Status",
    "Decision_Scope",
    "Candidate_Rank",
    "Is_Winner",
    "Candidate_Name",
    "Formula",
    "InChIKey",
    "Identity_Key",
    "Candidate_Entity_Score",
    "Legacy_V14_Entity_Score",
    "Collision_Energy_Bin",
    "Library_Source_Family",
    "Library_Search_Status",
    "Source_Participation_Explicit",
    "Reference_Record_Count",
    "Condition_Compatible_Record_Count",
    "Unique_Content_Profile_Count",
    "Independent_Single_Source_Profile_Count",
    "Best_Library_Score",
    "Best_Condition_Compatible_Score",
    "Representative_Spectrum_ID",
    "Representative_Source_Accession",
    "Reference_Instrument_Family",
    "Reference_Activation",
    "Reference_Collision_Energy",
    "Reference_Collision_Energy_Unit",
    "Collision_Energy_Delta",
    "Condition_Tier",
    "Rank_Influence_Eligible",
    "Rank_Influence_Allowed",
    "Rank_Influence_Basis",
    "Legacy_Best_Reference_Score",
    "Conditioned_Consensus_Score",
    "Conditioned_Multilibrary_Score",
    "Display_Only_All_Condition_Median",
    "Unconditioned_Consensus_Affects_Ranking",
    "Score_Is_Probability",
]

_LABELS_ZH = {
    "Entity_ID": "候选独立实体 ID",
    "Decision_Status": "结论状态",
    "Decision_Scope": "结论层级",
    "Candidate_Rank": "实体内候选排名",
    "Is_Winner": "是否当前第一候选",
    "Candidate_Name": "候选名称",
    "Formula": "分子式",
    "InChIKey": "InChIKey",
    "Identity_Key": "二维身份键",
    "Candidate_Entity_Score": "V15 实体排序分",
    "Legacy_V14_Entity_Score": "V14 审计实体分",
    "Collision_Energy_Bin": "样本 CE 分组",
    "Library_Source_Family": "谱库",
    "Library_Search_Status": "该谱库检索状态",
    "Source_Participation_Explicit": "已明确记录该谱库参与检索",
    "Reference_Record_Count": "该谱库参考谱记录数",
    "Condition_Compatible_Record_Count": "采集条件可比较参考谱数",
    "Unique_Content_Profile_Count": "去镜像后内容谱型数",
    "Independent_Single_Source_Profile_Count": "该谱库独有内容谱型数",
    "Best_Library_Score": "该谱库最佳谱图分",
    "Best_Condition_Compatible_Score": "该谱库最佳条件相容谱图分",
    "Representative_Spectrum_ID": "代表馆藏谱图编号",
    "Representative_Source_Accession": "代表来源 accession",
    "Reference_Instrument_Family": "参考谱仪器家族",
    "Reference_Activation": "参考谱碎裂方式",
    "Reference_Collision_Energy": "参考谱 CE",
    "Reference_Collision_Energy_Unit": "参考谱 CE 单位",
    "Collision_Energy_Delta": "样本与参考谱 CE 差",
    "Condition_Tier": "采集条件相容层级",
    "Rank_Influence_Eligible": "该谱库证据可进入条件共识",
    "Rank_Influence_Allowed": "本 CE 条件共识已达到改分门槛",
    "Rank_Influence_Basis": "本 CE 排序分依据",
    "Legacy_Best_Reference_Score": "本 CE 传统最佳参考谱分",
    "Conditioned_Consensus_Score": "本 CE 条件相容共识分",
    "Conditioned_Multilibrary_Score": "本 CE 条件化多库分",
    "Display_Only_All_Condition_Median": "全采集条件中位数（仅展示）",
    "Unconditioned_Consensus_Affects_Ranking": "全采集条件中位数参与排序",
    "Score_Is_Probability": "分数是否为概率",
}


def multilibrary_column_labels(language: str = "zh") -> list[str]:
    if language == "en":
        return list(MULTILIBRARY_COLUMNS)
    return [_LABELS_ZH[column] for column in MULTILIBRARY_COLUMNS]


def _entities_from_path(path: Path) -> Iterable[Mapping[str, Any]]:
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        yield from (
            entity
            for entity in payload.get("entities") or []
            if isinstance(entity, Mapping)
        )
        return
    with path.open("rb") as handle:
        yield from ijson.items(handle, "entities.item", use_float=True)


def multilibrary_evidence_rows(
    entities: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for entity in entities:
        decision = entity.get("decision") or {}
        winner_key = str(((entity.get("winner") or {}).get("identity_key") or ""))
        for candidate in entity.get("candidates") or []:
            if not isinstance(candidate, Mapping):
                continue
            display = candidate.get("candidate") or {}
            candidate_key = str(candidate.get("identity_key") or "")
            for ce_record in candidate.get("collision_energy_evidence") or []:
                if not isinstance(ce_record, Mapping):
                    continue
                evidence = ce_record.get("conditioned_multilibrary_evidence") or {}
                if not isinstance(evidence, Mapping):
                    continue
                for library in evidence.get("library_evidence_matrix") or []:
                    if not isinstance(library, Mapping):
                        continue
                    condition = library.get("representative_condition") or {}
                    rows.append({
                        "Entity_ID": entity.get("entity_id"),
                        "Decision_Status": decision.get("status"),
                        "Decision_Scope": decision.get("scope"),
                        "Candidate_Rank": candidate.get("entity_rank"),
                        "Is_Winner": bool(candidate_key and candidate_key == winner_key),
                        "Candidate_Name": display.get("compound_name"),
                        "Formula": display.get("formula"),
                        "InChIKey": display.get("inchikey"),
                        "Identity_Key": candidate_key,
                        "Candidate_Entity_Score": candidate.get("entity_score"),
                        "Legacy_V14_Entity_Score": candidate.get(
                            "legacy_v14_entity_score"
                        ),
                        "Collision_Energy_Bin": ce_record.get(
                            "collision_energy_bin"
                        ),
                        "Library_Source_Family": library.get("source_family"),
                        "Library_Search_Status": library.get("search_status"),
                        "Source_Participation_Explicit": library.get(
                            "source_participation_explicit"
                        ),
                        "Reference_Record_Count": library.get(
                            "reference_record_count"
                        ),
                        "Condition_Compatible_Record_Count": library.get(
                            "condition_compatible_record_count"
                        ),
                        "Unique_Content_Profile_Count": library.get(
                            "unique_content_profile_count"
                        ),
                        "Independent_Single_Source_Profile_Count": library.get(
                            "independent_single_source_profile_count"
                        ),
                        "Best_Library_Score": library.get("best_score"),
                        "Best_Condition_Compatible_Score": library.get(
                            "best_condition_compatible_score"
                        ),
                        "Representative_Spectrum_ID": library.get(
                            "representative_spectrum_id"
                        ),
                        "Representative_Source_Accession": library.get(
                            "representative_source_accession"
                        ),
                        "Reference_Instrument_Family": condition.get(
                            "reference_instrument_family"
                        ),
                        "Reference_Activation": condition.get(
                            "reference_activation"
                        ),
                        "Reference_Collision_Energy": condition.get(
                            "reference_collision_energy"
                        ),
                        "Reference_Collision_Energy_Unit": condition.get(
                            "collision_energy_unit"
                        ),
                        "Collision_Energy_Delta": condition.get(
                            "collision_energy_delta"
                        ),
                        "Condition_Tier": condition.get("tier"),
                        "Rank_Influence_Eligible": library.get(
                            "rank_influence_eligible"
                        ),
                        "Rank_Influence_Allowed": evidence.get(
                            "rank_influence_allowed"
                        ),
                        "Rank_Influence_Basis": evidence.get(
                            "rank_influence_basis"
                        ),
                        "Legacy_Best_Reference_Score": evidence.get(
                            "legacy_best_reference_score"
                        ),
                        "Conditioned_Consensus_Score": evidence.get(
                            "provenance_consensus_score"
                        ),
                        "Conditioned_Multilibrary_Score": evidence.get(
                            "conditioned_multilibrary_score"
                        ),
                        "Display_Only_All_Condition_Median": evidence.get(
                            "display_only_unconditioned_consensus_score"
                        ),
                        "Unconditioned_Consensus_Affects_Ranking": evidence.get(
                            "unconditioned_consensus_affects_ranking", False
                        ),
                        "Score_Is_Probability": evidence.get(
                            "score_is_probability", False
                        ),
                    })
    return rows


def read_multilibrary_evidence(path: Path) -> dict[str, Any]:
    path = Path(path)
    return {
        "schema_version": SCHEMA_VERSION,
        "source_artifact": str(path.resolve()),
        "rows": multilibrary_evidence_rows(_entities_from_path(path)),
        "score_is_probability": False,
        "unconditioned_consensus_affects_ranking": False,
    }


__all__ = [
    "MULTILIBRARY_COLUMNS",
    "SCHEMA_VERSION",
    "multilibrary_column_labels",
    "multilibrary_evidence_rows",
    "read_multilibrary_evidence",
]
