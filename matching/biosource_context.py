"""Conservative biosynthetic context for user-facing interpretation.

Identity evidence and biosynthetic context answer different questions.  This
module never changes retrieval, ranking or confidence.  It only turns a
reliably assigned chemical class into an explicitly class-level pathway
context, while refusing to fabricate a candidate-specific pathway when the
class or knowledge coverage is inadequate.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
import re
from typing import Any


_CLASS_PATHWAY_CONTEXT = {
    "flavonoid": "莽草酸途径 → 苯丙素 → 查耳酮 → 黄酮各亚型",
    "alkaloid": "氨基酸衍生途径（具体前体氨基酸与骨架仍需亚类证据）",
    "terpenoid": "MVA / MEP 异戊二烯途径（具体支路取决于萜类亚型）",
    "iridoid": "MEP 途径 → 单萜 → 环烯醚萜支路",
    "saponin": "MVA 途径 → 三萜／甾体骨架 + 后续糖基化",
    "phenolic_acid": "莽草酸／苯丙素途径",
    "coumarin": "莽草酸 → 苯丙素 → 香豆素支路",
    "lignan": "莽草酸 → 苯丙素二聚",
    "quinone": "醋酸-丙二酸聚酮途径或莽草酸支路（需醌亚类区分）",
    "steroid": "MVA 途径 → 角鲨烯 → 羊毛甾醇／甾醇支路",
    "xanthone": "莽草酸与聚酮混合来源（需山酮亚类区分）",
    "stilbene": "莽草酸／苯丙素 + 聚酮延伸（芪合成酶支路）",
    "diarylheptanoid": "莽草酸／苯丙素来源并伴随聚酮延伸（类别级背景）",
    "sugar": "中心碳代谢与糖代谢（类别级背景，不能区分具体单糖异构体）",
    "fatty_acyl": "脂肪酸生物合成、链延长与修饰（类别级背景）",
    "lipid": "脂肪酸／甘油脂／磷脂代谢（需脂质亚类与 sn 位点证据）",
    "amino_acid": "氨基酸代谢（具体前体与支路待亚类证据）",
    "nucleoside": "核苷酸／核苷代谢（具体碱基与糖基支路未解析）",
}

_PARTIAL_CLASS_CONTEXT = {
    "glycoside": "已支持糖基化修饰；苷元／母核的生物合成途径尚未解析",
}


def _text(value: Any) -> str:
    return " ".join(str(value or "").split())


def _unique(values: Iterable[Any], *, limit: int = 6) -> list[str]:
    output: list[str] = []
    for value in values:
        text = re.sub(
            r"^(?:生源|biosynthetic source)\s*[:：]\s*",
            "",
            _text(value),
            flags=re.IGNORECASE,
        )
        if text and text not in output:
            output.append(text)
        if len(output) >= limit:
            break
    return output


def biosource_context(
    match: Mapping[str, Any],
    supplied_pathways: Iterable[Any] = (),
) -> dict[str, Any]:
    """Return honest pathway context with explicit scope and coverage.

    ``supplied_pathways`` normally comes from the existing class rule engine.
    Even when a compound name is known, these labels remain chemical-class
    context unless a separately curated candidate-specific record exists.
    """

    supplied = _unique(supplied_pathways)
    classification = _text(
        match.get("display_class") or match.get("classification")
    ).lower()
    reliable = bool(match.get("classification_reliable"))
    compound_role = _text(match.get("compound_role")).lower()
    external_background_roles = {
        "external_contaminant",
        "curated_contaminant",
        "pesticide_screening_target",
        "environmental_contaminant",
        "industrial_contaminant",
        "transformation_product",
    }
    if compound_role in external_background_roles:
        return {
            "status": "not_applicable_external_background",
            "pathways": [],
            "reason": (
                "natural_product_biosynthesis_not_applicable_to_"
                "external_background"
            ),
            "reason_zh": (
                "该候选已归为潜在外源背景物；天然产物生物合成途径"
                "推断不适用"
            ),
            "classification": classification or None,
            "classification_reliable": reliable,
            "evidence_scope": "external_background_role_context",
            "candidate_specific": False,
            "knowledge_coverage": "not_applicable",
            "identity_rank_influence": False,
        }
    if (
        reliable
        and classification in _PARTIAL_CLASS_CONTEXT
        and supplied
        and all(
            value == _PARTIAL_CLASS_CONTEXT[classification]
            or "尚未解析" in value
            for value in supplied
        )
    ):
        return {
            "status": "partial_class_context",
            "pathways": supplied,
            "reason": "generic_class_does_not_resolve_core_biosynthesis",
            "reason_zh": (
                "当前可靠类别只说明修饰类型，不能据此编造母核生物合成途径"
            ),
            "classification": classification,
            "classification_reliable": True,
            "evidence_scope": "chemical_class_context_partial",
            "candidate_specific": False,
            "knowledge_coverage": "generic_class_only",
            "identity_rank_influence": False,
        }
    if supplied:
        return {
            "status": "class_context",
            "pathways": supplied,
            "reason": None,
            "reason_zh": None,
            "classification": classification or None,
            "classification_reliable": reliable,
            "evidence_scope": "chemical_class_context",
            "candidate_specific": False,
            "knowledge_coverage": "covered_class_rule",
            "identity_rank_influence": False,
        }
    if reliable and classification in _CLASS_PATHWAY_CONTEXT:
        return {
            "status": "class_context",
            "pathways": [_CLASS_PATHWAY_CONTEXT[classification]],
            "reason": "reliable_class_context_without_candidate_specific_rule",
            "reason_zh": (
                "基于可靠化学类别给出生物合成背景；不是候选特异生源证明"
            ),
            "classification": classification,
            "classification_reliable": True,
            "evidence_scope": "chemical_class_context",
            "candidate_specific": False,
            "knowledge_coverage": "covered_class_context",
            "identity_rank_influence": False,
        }
    if reliable and classification in _PARTIAL_CLASS_CONTEXT:
        return {
            "status": "partial_class_context",
            "pathways": [_PARTIAL_CLASS_CONTEXT[classification]],
            "reason": "generic_class_does_not_resolve_core_biosynthesis",
            "reason_zh": (
                "当前可靠类别只说明修饰类型，不能据此编造母核生物合成途径"
            ),
            "classification": classification,
            "classification_reliable": True,
            "evidence_scope": "chemical_class_context_partial",
            "candidate_specific": False,
            "knowledge_coverage": "generic_class_only",
            "identity_rank_influence": False,
        }
    reason = (
        "classification_not_reliable_for_biosource_context"
        if classification and not reliable
        else "biosource_knowledge_not_covered_for_current_class"
    )
    reason_zh = (
        "当前类别归属尚不可靠，不能据此生成生源途径"
        if classification and not reliable
        else "当前受控知识库尚未覆盖该候选的生源背景"
    )
    return {
        "status": "unavailable",
        "pathways": [],
        "reason": reason,
        "reason_zh": reason_zh,
        "classification": classification or None,
        "classification_reliable": reliable,
        "evidence_scope": "unavailable",
        "candidate_specific": False,
        "knowledge_coverage": "not_covered",
        "identity_rank_influence": False,
    }


__all__ = ["biosource_context"]
