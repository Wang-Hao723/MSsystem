"""One source of truth for user-facing confidence and score interpretation.

The guide is deliberately presentation-only.  It explains saved results but
never rescans spectra, changes candidate ranking, or promotes an identification
claim.  HTML, PDF, Excel, the result page and the integrated report all consume
the same versioned policy payload.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from html import escape
import json
from pathlib import Path
from typing import Any


CONFIDENCE_GUIDE_SCHEMA_VERSION = "1.0.0"
CONFIDENCE_GUIDE_POLICY_ID = "mssystem-confidence-interpretation-1.0.0"
CONFIDENCE_GUIDE_HTML_ID = "mssystem-confidence-guide"


_ZH_PAYLOAD: dict[str, Any] = {
    "language": "zh",
    "title": "怎样理解置信度与各项评分",
    "subtitle": "先看证据等级，再看可报告身份范围，最后看是否达到正式确证条件。",
    "read_first": (
        "最重要的边界：High 不等于已经确证该化合物；置信度和各项分数也不是“该化合物存在的概率”。"
        "只有正式结论状态显示“已确证／已发布”时，才可使用确定性表述。"
    ),
    "threshold_note": (
        "High／Medium／Low 是系统按本任务保存的质量误差、谱图相似度、碎片覆盖、候选竞争、"
        "谱库来源、异构体歧义和证据门槛等共同作出的规则判定，不是只看某一个分数。"
        "具体阈值与权重以本任务参数和系统版本为准；参数或版本不同的分数不宜直接横向比较。"
    ),
    "levels": [
        {
            "key": "confirmed",
            "label": "已确证／正式结论已发布",
            "meaning": "已通过当前项目规定的正式结论条件；通常需要同方法标准品及项目要求的正交证据。",
            "use": "可在限定的样品、方法和批次范围内表述为“确认检出该化合物”。",
            "avoid": "不要把该结论外推到未测样品，也不要超出已发布结论所声明的结构与方法范围。",
            "next_step": "保留标准品、原始谱图、保留时间及发布记录，供审计和复现。",
            "tone": "confirmed",
        },
        {
            "key": "high",
            "label": "High（高可信证据）",
            "meaning": "在页面注明的身份范围内具有较强 MS/MS／谱库支持，但仍属于谱学鉴定，不自动等同标准品确证。",
            "use": "可作为主要结果、优先报告对象和标准品确认清单；表述时必须带身份范围。",
            "avoid": "不能单凭 High 写成“已确定／已确证该化合物”，也不能自动声明立体化学、取代位置或唯一异构体。",
            "next_step": "重要结论优先补同方法标准品保留时间与 MS/MS，必要时增加离子淌度、NMR 或其他正交证据。",
            "tone": "high",
        },
        {
            "key": "medium",
            "label": "Medium（中等证据候选）",
            "meaning": "候选具有一定一致性，但至少一个关键证据轴仍不充分或存在竞争与歧义。",
            "use": "适合候选短名单、靶向复采、通路／类别探索和后续验证设计。",
            "avoid": "不能作为最终化合物鉴定，也不宜据此给出确定的生物学或定量结论。",
            "next_step": "优先复核原始 MS/MS、增加重复针和诊断碎片，并补标准品或其他独立证据。",
            "tone": "medium",
        },
        {
            "key": "low",
            "label": "Low／证据不足",
            "meaning": "只有弱、局部或不完整的谱库线索；系统保留它是为了不丢失潜在线索，而不是把它列为已鉴定结果。",
            "use": "仅用于探索、人工复核、补谱库和下一轮采集规划。",
            "avoid": "不要计入已鉴定化合物数量，不要写成“检出该命名化合物”。",
            "next_step": "重新采集更高质量 MS/MS，核对加合物、质量误差、空白和候选谱库来源。",
            "tone": "low",
        },
        {
            "key": "unknown",
            "label": "Unknown／未匹配",
            "meaning": "当前数据和规则下没有可辩护的命名候选通过门槛。",
            "use": "可作为未知特征、分子网络节点或后续非靶向发现对象。",
            "avoid": "未匹配不等于该化合物不存在，也不等于该信号没有研究价值。",
            "next_step": "改善谱图质量、扩充合适极性／仪器类型的谱库，或采用分子式、类别和网络层面的注释。",
            "tone": "unknown",
        },
        {
            "key": "background",
            "label": "空白／背景相关",
            "meaning": "存在谱库候选，但信号也与任务空白、携带污染或系统背景相关。",
            "use": "用于污染、背景和方法学审计；核验批次及色谱方法可比性后再决定是否排除。",
            "avoid": "默认不计作样品支持的主要化合物；背景标记也不等于系统已执行跨针谱图相减。",
            "next_step": "检查空白角色、批次、保留时间、样品／空白强度比和重复针一致性。",
            "tone": "background",
        },
    ],
    "identity_scopes": [
        {
            "label": "完整身份／确切结构",
            "meaning": "包括系统声明范围内的全部结构身份信息；若涉及立体化学，必须有相应证据。",
            "report": "只有该层级本身达到要求时，才可报告完整化合物身份。",
        },
        {
            "label": "单一二维连接结构",
            "meaning": "原子连接关系可支持到一个二维结构，但立体化学或其他完整身份信息可能尚未解决。",
            "report": "建议写“与该二维连接结构一致的 High 级候选／检出”。",
        },
        {
            "label": "异构体组",
            "meaning": "证据支持同一分子式的一组竞争候选，当前 MS/MS 不能可靠点名其中一个。",
            "report": "只能报告“该异构体组中的某一成员”，不能把首位候选当成唯一结构。",
        },
        {
            "label": "化合物类别",
            "meaning": "只能支持类别或骨架层面，不能支持具体命名化合物。",
            "report": "用于类别组成和探索性解释，不能替代结构鉴定。",
        },
        {
            "label": "推定候选／谱图线索／文献目标线索",
            "meaning": "分别表示候选级、仅相似谱图级或文献目标级证据，结构确定性依次受限。",
            "report": "必须保留“候选／线索”措辞，不得升级为确定检出。",
        },
    ],
    "scores": [
        {
            "key": "cosine",
            "label": "Cosine／余弦相似度（0–1）",
            "meaning": "比较样品谱与库谱的碎片整体形状；越高表示越相似。",
            "boundary": "不是身份概率；异构体或共享骨架也可能得到高分。",
        },
        {
            "key": "entropy",
            "label": "Entropy／谱图熵相似度（0–1）",
            "meaning": "从碎片强度分布和信息量角度补充衡量两张谱的相似性。",
            "boundary": "与 Cosine 相关，属于同一谱图证据家族，不能当作第二份完全独立的确认。",
        },
        {
            "key": "reverse_match",
            "label": "Reverse match／反向覆盖",
            "meaning": "查看库谱中的特征峰有多少被样品谱覆盖；低值提示库谱关键峰缺失。",
            "boundary": "覆盖高仍不能单独区分近缘异构体。",
        },
        {
            "key": "mass_error",
            "label": "质量误差（ppm／Da）与质量得分",
            "meaning": "衡量样品前体与候选理论／库前体质量是否一致；绝对误差通常越接近 0 越好。",
            "boundary": "精确质量相同的异构体很多，质量吻合本身不能确认结构。",
        },
        {
            "key": "fragment_support",
            "label": "匹配碎片数、碎片得分与强度得分",
            "meaning": "分别反映直接匹配离子数量／覆盖和匹配峰强度关系。",
            "boundary": "数量多通常更有支持力，但常见碎片并不具有唯一性；少碎片会限制置信度。",
        },
        {
            "key": "final_score",
            "label": "Final score／综合匹配分",
            "meaning": "按当前任务保存的权重综合谱图、质量、碎片和强度等证据，用于评价候选匹配质量。",
            "boundary": "不是后验概率；不同参数、版本或数据质量下不宜直接比较。",
        },
        {
            "key": "rank_score",
            "label": "Rank score／候选排序分",
            "meaning": "在当前候选集合中安排先后顺序的分数；可能包含已声明的重排或一致性信息。",
            "boundary": "第一名只是当前候选中的最优者，不代表第一名已经被确证。",
        },
        {
            "key": "candidate_margin",
            "label": "Candidate margin／候选差距",
            "meaning": "首位与次位候选的分数间隔；差距越小，候选竞争越难解决。",
            "boundary": "大差距只说明在当前库和规则下领先，不排除谱库外候选。",
        },
        {
            "key": "confidence_label",
            "label": "High／Medium／Low 置信标签",
            "meaning": "结合多项门槛、证据来源和身份范围作出的规则等级。",
            "boundary": "不是由一个固定分数直接换算，也不是 90%／70% 之类的概率。",
        },
        {
            "key": "formal_claim",
            "label": "项目 q 值／正式结论状态",
            "meaning": "项目完成校准并满足发布规则后，用于控制和发布正式结论。",
            "boundary": "缺少 q 值表示“未完成项目级校准”，不是 q=0；只有已发布状态才代表正式结论。",
        },
        {
            "key": "auxiliary_scores",
            "label": "分类概率、专家分、MSClassify 与靶向证据等级",
            "meaning": "分别回答类别预测、规则解释或精确质量／同位素靶向支持等旁证问题。",
            "boundary": "这些是独立解释轴，不能替代化合物身份置信度；专家分为未校准辅助摘要，靶向精确质量也不是结构确证。",
        },
    ],
    "safeguards": {
        "confidence_is_probability": False,
        "confidence_is_project_qvalue": False,
        "high_means_confirmed": False,
        "unmatched_means_absent": False,
        "auxiliary_scores_change_identity_rank_by_definition": False,
    },
}


_EN_PAYLOAD: dict[str, Any] = {
    "language": "en",
    "title": "How to interpret confidence and scores",
    "subtitle": "Read the evidence tier, then the reportable identity scope, and finally the formal confirmation status.",
    "read_first": (
        "Most important boundary: High does not mean that the compound has been confirmed, and neither confidence nor any displayed score is the probability that the compound is present. "
        "Use definitive wording only when the formal conclusion is explicitly released/confirmed."
    ),
    "threshold_note": (
        "High/Medium/Low is a rule-based decision using the saved task's mass error, spectral similarity, fragment coverage, candidate competition, library provenance, isomer ambiguity and evidence gates. "
        "It is not determined by one score. Thresholds and weights follow the saved task settings and system version; scores from different settings or versions should not be compared directly."
    ),
    "levels": [
        {
            "key": "confirmed", "label": "Confirmed / formal conclusion released",
            "meaning": "The current project's formal release criteria were met, normally including an authentic standard under the same method and any required orthogonal evidence.",
            "use": "Within the stated sample, method and batch scope, it may be reported as 'confirmed detected'.",
            "avoid": "Do not extrapolate to unmeasured samples or beyond the identity and method scope of the released claim.",
            "next_step": "Retain standard, raw spectra, retention-time and release records for audit and reproducibility.", "tone": "confirmed",
        },
        {
            "key": "high", "label": "High (strong evidence)",
            "meaning": "Strong MS/MS/library support at the stated identity scope, but not automatically authentic-standard confirmation.",
            "use": "Suitable as a main result, reporting priority and standard-confirmation shortlist; always retain the identity scope.",
            "avoid": "Do not write 'confirmed compound' from High alone or claim stereochemistry, substitution position or one unique isomer without supporting scope.",
            "next_step": "For important claims, add same-method standard RT and MS/MS and, where needed, ion mobility, NMR or other orthogonal evidence.", "tone": "high",
        },
        {
            "key": "medium", "label": "Medium (plausible candidate)",
            "meaning": "The candidate is plausible, but at least one key evidence axis is incomplete, competitive or ambiguous.",
            "use": "Useful for shortlisting, targeted reacquisition, pathway/class exploration and validation planning.",
            "avoid": "Do not treat it as a final identification or base definitive biological or quantitative claims on it.",
            "next_step": "Review raw MS/MS, add replicate and diagnostic-fragment evidence, and obtain a standard or another independent line of evidence.", "tone": "medium",
        },
        {
            "key": "low", "label": "Low / insufficient evidence",
            "meaning": "Only weak, partial or incomplete library evidence is present; it is retained as a lead, not as an identified result.",
            "use": "Use only for exploration, manual review, library expansion and acquisition planning.",
            "avoid": "Do not count it as identified or report the named compound as detected.",
            "next_step": "Acquire better MS/MS and check adduct, mass error, blanks and library provenance.", "tone": "low",
        },
        {
            "key": "unknown", "label": "Unknown / unmatched",
            "meaning": "No defensible named candidate passed the current data and rules.",
            "use": "Retain as an unknown feature, molecular-network node or discovery target.",
            "avoid": "Unmatched does not prove absence and does not mean the signal lacks research value.",
            "next_step": "Improve spectra, expand polarity/instrument-appropriate libraries or annotate at formula, class or network level.", "tone": "unknown",
        },
        {
            "key": "background", "label": "Blank/background related",
            "meaning": "A library candidate exists, but the signal is also associated with the task blank, carryover or system background.",
            "use": "Use for contamination/background and method audit after checking batch and method comparability.",
            "avoid": "Do not count it as a sample-supported main compound by default; the flag does not mean cross-injection spectral subtraction was performed.",
            "next_step": "Check blank role, batch, RT, sample/blank ratio and replicate consistency.", "tone": "background",
        },
    ],
    "identity_scopes": [
        {"label": "Full identity / exact structure", "meaning": "All identity information within the declared scope, including stereochemistry when claimed, is supported.", "report": "Report a full compound identity only when this scope itself meets the required evidence."},
        {"label": "Single 2D connectivity", "meaning": "One atom-connectivity graph is supportable, while stereochemistry or other full-identity details may remain unresolved.", "report": "Prefer 'High candidate/detection consistent with this 2D connectivity'."},
        {"label": "Isomer group", "meaning": "Evidence supports a competing same-formula group, but current MS/MS cannot reliably select one member.", "report": "Report only 'a member of this isomer group'; do not promote the lead candidate to a unique structure."},
        {"label": "Compound class", "meaning": "Only class or scaffold-level evidence is supportable.", "report": "Use for class composition and exploratory interpretation, not structural identification."},
        {"label": "Putative candidate / spectral lead / literature target lead", "meaning": "Candidate-level, similar-spectrum-only or literature-target evidence with progressively constrained structural certainty.", "report": "Retain candidate/lead wording and do not upgrade it to confirmed detection."},
    ],
    "scores": [
        {"key": "cosine", "label": "Cosine similarity (0–1)", "meaning": "Overall fragment-pattern similarity between sample and library spectra; higher is more similar.", "boundary": "Not an identity probability; isomers and shared scaffolds may also score highly."},
        {"key": "entropy", "label": "Spectral entropy similarity (0–1)", "meaning": "Complementary similarity based on fragment intensity distribution and information content.", "boundary": "Correlated with cosine and not a second fully independent confirmation."},
        {"key": "reverse_match", "label": "Reverse match / reference coverage", "meaning": "How much of the reference spectrum is covered by the sample spectrum; a low value flags missing reference features.", "boundary": "High coverage alone cannot resolve close isomers."},
        {"key": "mass_error", "label": "Mass error (ppm/Da) and mass score", "meaning": "Agreement between sample precursor and candidate theoretical/library mass; absolute error closer to zero is usually better.", "boundary": "Many isomers share exact mass, so mass agreement alone cannot confirm structure."},
        {"key": "fragment_support", "label": "Matched ions, fragment score and intensity score", "meaning": "Direct matched-ion count/coverage and matched-peak intensity agreement.", "boundary": "More support is generally better, but common fragments are not unique and sparse spectra constrain confidence."},
        {"key": "final_score", "label": "Final score / composite match score", "meaning": "A weighted summary of spectral, mass, fragment and intensity evidence under the saved task settings.", "boundary": "Not a posterior probability and not directly comparable across different settings, versions or data quality."},
        {"key": "rank_score", "label": "Rank score", "meaning": "Orders candidates within the current candidate set and may include declared reranking or consistency information.", "boundary": "Rank 1 is only the best current candidate, not a confirmed identity."},
        {"key": "candidate_margin", "label": "Candidate margin", "meaning": "Score gap between the leading and runner-up candidates; a small gap indicates unresolved competition.", "boundary": "A large gap only describes the current library and rules; it does not exclude out-of-library candidates."},
        {"key": "confidence_label", "label": "High/Medium/Low confidence label", "meaning": "A rule tier combining several gates, provenance and reportable identity scope.", "boundary": "Not a single-score conversion and not a 90%/70% probability."},
        {"key": "formal_claim", "label": "Project q-value / formal conclusion status", "meaning": "Used to control and release formal conclusions after project calibration and release rules are met.", "boundary": "A missing q-value means not calibrated, not q=0; only a released status is a formal conclusion."},
        {"key": "auxiliary_scores", "label": "Classification probability, expert score, MSClassify and targeted evidence tiers", "meaning": "Separate axes for class prediction, rule explanation or exact-mass/isotope targeted support.", "boundary": "They do not replace identity confidence. Expert score is an uncalibrated explanatory summary, and exact-mass targeting is not structure confirmation."},
    ],
    "safeguards": {
        "confidence_is_probability": False,
        "confidence_is_project_qvalue": False,
        "high_means_confirmed": False,
        "unmatched_means_absent": False,
        "auxiliary_scores_change_identity_rank_by_definition": False,
    },
}


def _language(language: str) -> str:
    return "en" if str(language or "").lower().startswith("en") else "zh"


@lru_cache(maxsize=1)
def _static_guide_payloads() -> dict[str, Any]:
    """Load the browser-readable policy file used as the operational source."""

    path = (
        Path(__file__).resolve().parent
        / "static" / "policies" / "confidence-guide.json"
    )
    try:
        payloads = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return {}
    if not isinstance(payloads, dict):
        return {}
    return payloads


def confidence_guide_payload(language: str = "zh") -> dict[str, Any]:
    """Return a detached, versioned payload safe for API consumers."""

    normalized = _language(language)
    static_payload = _static_guide_payloads().get(normalized)
    payload = deepcopy(
        static_payload
        if isinstance(static_payload, dict)
        else (_EN_PAYLOAD if normalized == "en" else _ZH_PAYLOAD)
    )
    payload.update({
        "schema_version": CONFIDENCE_GUIDE_SCHEMA_VERSION,
        "policy_id": CONFIDENCE_GUIDE_POLICY_ID,
    })
    return payload


def confidence_guide_table_rows(language: str = "zh") -> list[dict[str, str]]:
    """Flatten the guide for Excel/CSV-style report sheets."""

    payload = confidence_guide_payload(language)
    rows: list[dict[str, str]] = []
    if payload["language"] == "en":
        rows.extend([
            {"Section": "READ FIRST", "Name": "Core boundary", "Meaning": payload["read_first"], "Recommended use": "Read before interpreting any named compound.", "Do not infer": "High or a numeric score equals confirmation/probability.", "Next step": payload["threshold_note"]},
        ])
        for item in payload["levels"]:
            rows.append({"Section": "Evidence tier", "Name": item["label"], "Meaning": item["meaning"], "Recommended use": item["use"], "Do not infer": item["avoid"], "Next step": item["next_step"]})
        for item in payload["identity_scopes"]:
            rows.append({"Section": "Identity scope", "Name": item["label"], "Meaning": item["meaning"], "Recommended use": item["report"], "Do not infer": "Do not report beyond this scope.", "Next step": "Add orthogonal evidence if a narrower identity is required."})
        for item in payload["scores"]:
            rows.append({"Section": "Score/field", "Name": item["label"], "Meaning": item["meaning"], "Recommended use": "Use with the other saved evidence fields and the confidence tier.", "Do not infer": item["boundary"], "Next step": "Review task settings and raw evidence when the result is decision-critical."})
        return rows

    rows.extend([
        {"说明类型": "必读", "名称": "核心边界", "它代表什么": payload["read_first"], "可以怎样使用": "在解读任何命名化合物之前先阅读。", "不能怎样理解": "不能把 High 或数字分数直接理解成确证或概率。", "建议下一步": payload["threshold_note"]},
    ])
    for item in payload["levels"]:
        rows.append({"说明类型": "证据等级", "名称": item["label"], "它代表什么": item["meaning"], "可以怎样使用": item["use"], "不能怎样理解": item["avoid"], "建议下一步": item["next_step"]})
    for item in payload["identity_scopes"]:
        rows.append({"说明类型": "身份范围", "名称": item["label"], "它代表什么": item["meaning"], "可以怎样使用": item["report"], "不能怎样理解": "不能超出该身份范围报告。", "建议下一步": "如需更具体身份，补充相应正交证据。"})
    for item in payload["scores"]:
        rows.append({"说明类型": "评分／字段", "名称": item["label"], "它代表什么": item["meaning"], "可以怎样使用": "与其他保存的证据字段和置信等级一起判断。", "不能怎样理解": item["boundary"], "建议下一步": "结果用于关键决策时，回查任务参数和原始证据。"})
    return rows


def confidence_guide_html(language: str = "zh") -> str:
    """Render a report-first guide with compact overview and full disclosures."""

    payload = confidence_guide_payload(language)
    zh = payload["language"] == "zh"
    routes = (
        [
            ("01", "先看证据等级", "High／Medium／Low 回答证据强弱"),
            ("02", "再看身份范围", "二维连接、异构体组与类别不可混用"),
            ("03", "最后看正式状态", "只有已确证／已发布才可确定性表述"),
        ]
        if zh else
        [
            ("01", "Read the evidence tier", "High / Medium / Low describes evidence strength"),
            ("02", "Check identity scope", "2-D identity, isomer group and class differ"),
            ("03", "Check formal status", "Use definitive wording only after release"),
        ]
    )
    route_html = "".join(
        f"<article><span>{number}</span><b>{escape(title)}</b><small>{escape(note)}</small></article>"
        for number, title, note in routes
    )
    tier_map = "".join(
        "<article class='mcg-tier mcg-{tone}'><span></span><div><b>{label}</b>"
        "<p>{meaning}</p></div></article>".format(
            tone=escape(item["tone"]),
            label=escape(item["label"]),
            meaning=escape(item["meaning"]),
        )
        for item in payload["levels"]
    )
    level_details = "".join(
        "<article class='mcg-level mcg-{tone}'><h3>{label}</h3>"
        "<p><b>{means}</b>{meaning}</p><div class='mcg-level-actions'>"
        "<p><b>{use_label}</b>{use}</p><p class='mcg-avoid'><b>{avoid_label}</b>{avoid}</p>"
        "<p><b>{next_label}</b>{next_step}</p></div></article>".format(
            tone=escape(item["tone"]), label=escape(item["label"]),
            means="代表：" if zh else "Meaning: ", meaning=escape(item["meaning"]),
            use_label="可以用来：" if zh else "Use for: ", use=escape(item["use"]),
            avoid_label="不能理解成：" if zh else "Do not infer: ", avoid=escape(item["avoid"]),
            next_label="建议下一步：" if zh else "Next step: ", next_step=escape(item["next_step"]),
        )
        for item in payload["levels"]
    )
    scopes = "".join(
        "<article><span>{number:02d}</span><div><b>{label}</b><p>{meaning}</p>"
        "<em>{prefix}{report}</em></div></article>".format(
            number=index,
            label=escape(item["label"]),
            meaning=escape(item["meaning"]),
            prefix="建议表述：" if zh else "Recommended wording: ",
            report=escape(item["report"]),
        )
        for index, item in enumerate(payload["identity_scopes"], start=1)
    )
    scores = "".join(
        "<article><span>{number:02d}</span><div><b>{label}</b><p>{meaning}</p>"
        "<em>{prefix}{boundary}</em></div></article>".format(
            number=index,
            label=escape(item["label"]),
            meaning=escape(item["meaning"]),
            prefix="误读边界：" if zh else "Boundary: ",
            boundary=escape(item["boundary"]),
        )
        for index, item in enumerate(payload["scores"], start=1)
    )
    return f"""
<style>
#{CONFIDENCE_GUIDE_HTML_ID}{{--mcg-navy:#102a43;--mcg-teal:#16746f;margin:1.35rem auto;max-width:1180px;color:#203747;font:14px/1.62 -apple-system,BlinkMacSystemFont,'Segoe UI','Microsoft YaHei',Arial,sans-serif;box-sizing:border-box}}
#{CONFIDENCE_GUIDE_HTML_ID},#{CONFIDENCE_GUIDE_HTML_ID} *{{box-sizing:border-box}} #{CONFIDENCE_GUIDE_HTML_ID} p,#{CONFIDENCE_GUIDE_HTML_ID} h2,#{CONFIDENCE_GUIDE_HTML_ID} h3{{margin:0}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-shell{{overflow:hidden;border:1px solid #d7e2e9;border-radius:16px;background:#fff;box-shadow:0 12px 34px rgba(22,42,60,.08)}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-head{{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:1.25rem;align-items:center;padding:1.15rem 1.3rem;border-left:6px solid var(--mcg-teal);background:linear-gradient(110deg,#f5faf9 0%,#fff 70%)}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-kicker{{display:block;color:#188179;font-size:.68rem;font-weight:850;letter-spacing:.15em}} #{CONFIDENCE_GUIDE_HTML_ID} h2{{margin:.12rem 0 .18rem;color:var(--mcg-navy);font-size:1.35rem;line-height:1.3}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-sub{{color:#667b89;font-size:.82rem}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-boundaries{{display:flex;flex-wrap:wrap;justify-content:flex-end;gap:.42rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-boundaries span{{padding:.34rem .58rem;border:1px solid #cddde4;border-radius:999px;background:#fff;color:#49616f;font-size:.72rem;white-space:nowrap}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-boundaries b{{color:#163f57}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-body{{padding:1rem 1.2rem 1.05rem;background:#fbfcfd}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-alert{{display:grid;grid-template-columns:28px minmax(0,1fr);gap:.7rem;padding:.78rem .9rem;border:1px solid #ead7a7;border-radius:11px;background:#fffaf0;color:#684d18}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-alert>span{{display:grid;width:27px;height:27px;place-items:center;border-radius:8px;background:#f1c65e;color:#5c420d;font-weight:850}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-alert b{{display:block;margin-bottom:.12rem;font-size:.76rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-alert p{{font-size:.76rem}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-threshold{{margin:.55rem 0 .75rem;padding:.62rem .75rem;border-left:3px solid #7ba9a3;background:#edf6f4;color:#536b73;font-size:.72rem}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-route{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.5rem;margin-bottom:.85rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-route article{{position:relative;padding:.66rem .7rem .62rem 2.55rem;border:1px solid #dce6eb;border-radius:10px;background:#fff}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-route span{{position:absolute;left:.7rem;top:.68rem;color:#258078;font-size:.66rem;font-weight:850;letter-spacing:.08em}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-route b{{display:block;color:#294557;font-size:.74rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-route small{{display:block;margin-top:.08rem;color:#748592;font-size:.66rem;line-height:1.45}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-section-title{{display:flex;align-items:center;justify-content:space-between;gap:1rem;margin:.15rem 0 .45rem;color:#1d3a4d}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-section-title b{{font-size:.83rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-section-title span{{color:#81909b;font-size:.65rem}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier-map{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.48rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier{{display:grid;grid-template-columns:8px minmax(0,1fr);gap:.55rem;min-width:0;padding:.65rem .7rem;border:1px solid #dce5ea;border-radius:10px;background:#fff}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier>span{{width:8px;height:8px;margin-top:.26rem;border-radius:50%;background:var(--tier,#748191)}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier b{{color:#294758;font-size:.72rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier p{{margin-top:.1rem;color:#6a7b87;font-size:.65rem;line-height:1.48}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-confirmed{{--tier:#147d64}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-high{{--tier:#2474b8}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-medium{{--tier:#c28a20}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-low{{--tier:#bc5a31}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-background{{--tier:#7456a6}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-unknown{{--tier:#748191}}
#{CONFIDENCE_GUIDE_HTML_ID} details{{margin-top:.55rem;border:1px solid #d9e3e8;border-radius:10px;background:#fff}} #{CONFIDENCE_GUIDE_HTML_ID} summary{{display:flex;align-items:center;justify-content:space-between;gap:1rem;padding:.68rem .8rem;color:#294557;font-size:.73rem;font-weight:780;cursor:pointer;list-style:none}} #{CONFIDENCE_GUIDE_HTML_ID} summary::-webkit-details-marker{{display:none}} #{CONFIDENCE_GUIDE_HTML_ID} summary:after{{color:#738592;font-size:1rem;content:'+'}} #{CONFIDENCE_GUIDE_HTML_ID} details[open]>summary{{border-bottom:1px solid #e3eaee;background:#f7fafb}} #{CONFIDENCE_GUIDE_HTML_ID} details[open]>summary:after{{content:'−'}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-level-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.55rem;padding:.65rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level{{padding:.72rem;border:1px solid #dde6ea;border-top:3px solid var(--tier,#748191);border-radius:9px;background:#fff}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level h3{{color:#294557;font-size:.75rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level p{{margin-top:.24rem;color:#627581;font-size:.67rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level-actions{{display:grid;gap:.22rem;margin-top:.38rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level-actions p{{padding:.36rem .42rem;border-radius:6px;background:#f4f8f7}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level-actions .mcg-avoid{{background:#fff4ef;color:#91492e}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:0 .8rem;padding:.35rem .8rem .7rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list article{{display:grid;grid-template-columns:28px minmax(0,1fr);gap:.55rem;padding:.62rem 0;border-bottom:1px solid #e7edf0}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list article>span{{display:grid;width:26px;height:26px;place-items:center;border-radius:7px;background:#edf4f5;color:#47717a;font-size:.6rem;font-weight:800}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list b{{display:block;color:#2d495a;font-size:.7rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list p{{margin-top:.1rem;color:#687a86;font-size:.65rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list em{{display:block;margin-top:.2rem;color:#28746d;font-size:.63rem;font-style:normal}}
#{CONFIDENCE_GUIDE_HTML_ID} .mcg-score-list em{{padding:.28rem .35rem;border-radius:5px;background:#fff4ef;color:#93492d}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-policy{{padding:.6rem .15rem 0;color:#98a4ad;font:10px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace;text-align:right}}
@media(max-width:780px){{#{CONFIDENCE_GUIDE_HTML_ID} .mcg-head{{grid-template-columns:1fr}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-boundaries{{justify-content:flex-start}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-route,#{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier-map{{grid-template-columns:1fr}}}}
@media(max-width:560px){{#{CONFIDENCE_GUIDE_HTML_ID} .mcg-body{{padding:.8rem}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level-grid,#{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list{{grid-template-columns:1fr}}}}
@media print{{#{CONFIDENCE_GUIDE_HTML_ID}{{break-before:page}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-shell{{box-shadow:none}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-tier-map{{display:none}} #{CONFIDENCE_GUIDE_HTML_ID} details:not([open])>*:not(summary){{display:block!important}} #{CONFIDENCE_GUIDE_HTML_ID} details>summary{{break-after:avoid}} #{CONFIDENCE_GUIDE_HTML_ID} .mcg-level,#{CONFIDENCE_GUIDE_HTML_ID} .mcg-card-list article{{break-inside:avoid}}}}
</style>
<section id="{CONFIDENCE_GUIDE_HTML_ID}" data-policy-id="{CONFIDENCE_GUIDE_POLICY_ID}" data-schema-version="{CONFIDENCE_GUIDE_SCHEMA_VERSION}">
  <div class="mcg-shell">
    <header class="mcg-head"><div><span class="mcg-kicker">REPORT READING GUIDE</span><h2>{escape(payload['title'])}</h2><p class="mcg-sub">{escape(payload['subtitle'])}</p></div><div class="mcg-boundaries"><span><b>High</b> ≠ {'已确证' if zh else 'confirmed'}</span><span><b>{'分数' if zh else 'Score'}</b> ≠ {'存在概率' if zh else 'presence probability'}</span></div></header>
    <div class="mcg-body">
      <div class="mcg-alert"><span>!</span><div><b>{'阅读报告前先记住' if zh else 'Remember before reading'}</b><p>{escape(payload['read_first'])}</p></div></div>
      <p class="mcg-threshold">{escape(payload['threshold_note'])}</p>
      <div class="mcg-route">{route_html}</div>
      <div class="mcg-section-title"><b>{'证据等级速查' if zh else 'Evidence-tier quick map'}</b><span>{'先定位等级，再按需展开完整定义' if zh else 'Locate the tier first; open definitions when needed'}</span></div>
      <div class="mcg-tier-map">{tier_map}</div>
      <details class="mcg-detail"><summary><span>{'01 · 完整证据等级说明' if zh else '01 · Full evidence-tier definitions'}</span><small>{len(payload['levels'])} {'级' if zh else 'tiers'}</small></summary><div class="mcg-level-grid">{level_details}</div></details>
      <details class="mcg-detail"><summary><span>{'02 · 身份范围决定最多能把结构说到哪一步' if zh else '02 · Identity scope limits how specifically a structure may be named'}</span><small>{len(payload['identity_scopes'])} {'层' if zh else 'levels'}</small></summary><div class="mcg-card-list">{scopes}</div></details>
      <details class="mcg-detail"><summary><span>{'03 · 各项评分和辅助字段怎样读' if zh else '03 · How to read scores and auxiliary fields'}</span><small>{len(payload['scores'])} {'项' if zh else 'fields'}</small></summary><div class="mcg-card-list mcg-score-list">{scores}</div></details>
      <div class="mcg-policy">{escape(payload['policy_id'])} · {escape(payload['schema_version'])}</div>
    </div>
  </div>
</section>""".strip()


def inject_confidence_guide_html(html: str, language: str = "zh") -> str:
    """Insert the canonical guide once near the beginning of an HTML report."""

    source = str(html or "")
    if f'id="{CONFIDENCE_GUIDE_HTML_ID}"' in source or f"id='{CONFIDENCE_GUIDE_HTML_ID}'" in source:
        return source
    guide = confidence_guide_html(language)
    lower = source.lower()
    for marker in ("</header>", "</h1>"):
        index = lower.find(marker)
        if index >= 0:
            end = index + len(marker)
            return f"{source[:end]}\n{guide}\n{source[end:]}"
    body_index = lower.find("<body")
    if body_index >= 0:
        body_end = lower.find(">", body_index)
        if body_end >= 0:
            return f"{source[:body_end + 1]}\n{guide}\n{source[body_end + 1:]}"
    return f"{guide}\n{source}"


def confidence_guide_pdf_story(
    language: str,
    *,
    regular_font: str,
    bold_font: str,
) -> list[Any]:
    """Build a two-page, report-style reading guide for PDF exports."""

    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.platypus import (
        KeepTogether, PageBreak, Paragraph, Spacer, Table, TableStyle,
    )

    payload = confidence_guide_payload(language)
    zh = payload["language"] == "zh"
    styles = getSampleStyleSheet()
    normal = ParagraphStyle(
        "ConfidenceGuideNormalZh" if zh else "ConfidenceGuideNormalEn",
        parent=styles["Normal"], fontName=regular_font, fontSize=7.4,
        leading=10.1, textColor=colors.HexColor("#425866"),
        wordWrap="CJK", splitLongWords=True,
    )
    small = ParagraphStyle(
        "ConfidenceGuideSmallZh" if zh else "ConfidenceGuideSmallEn",
        parent=normal, fontSize=6.8, leading=9.1,
    )
    card = ParagraphStyle(
        "ConfidenceGuideCardZh" if zh else "ConfidenceGuideCardEn",
        parent=normal, fontSize=7.15, leading=9.55,
    )
    white_title = ParagraphStyle(
        "ConfidenceGuideWhiteTitleZh" if zh else "ConfidenceGuideWhiteTitleEn",
        parent=normal, fontName=bold_font, fontSize=14.5, leading=17,
        textColor=colors.white,
    )
    white_kicker = ParagraphStyle(
        "ConfidenceGuideWhiteKickerZh" if zh else "ConfidenceGuideWhiteKickerEn",
        parent=normal, fontName=bold_font, fontSize=6.5, leading=8,
        textColor=colors.HexColor("#9FE4DE"),
    )
    section = ParagraphStyle(
        "ConfidenceGuideSectionZh" if zh else "ConfidenceGuideSectionEn",
        parent=normal, fontName=bold_font, fontSize=9.2, leading=11,
        textColor=colors.HexColor("#173B50"), spaceBefore=1, spaceAfter=3,
    )
    policy = ParagraphStyle(
        "ConfidenceGuidePolicyZh" if zh else "ConfidenceGuidePolicyEn",
        parent=small, fontSize=5.8, leading=7, textColor=colors.HexColor("#8A98A3"),
        alignment=2,
    )
    tone_colors = {
        "confirmed": "#147D64", "high": "#2474B8", "medium": "#C28A20",
        "low": "#BC5A31", "background": "#7456A6", "unknown": "#748191",
    }

    hero = Table([
        [Paragraph("REPORT READING GUIDE", white_kicker)],
        [Paragraph(escape(payload["title"]), white_title)],
        [Paragraph(escape(payload["subtitle"]), ParagraphStyle(
            "ConfidenceGuideHeroSubZh" if zh else "ConfidenceGuideHeroSubEn",
            parent=normal, fontSize=7.2, leading=9.2,
            textColor=colors.HexColor("#D6E8EA"),
        ))],
    ], colWidths=[17.5 * cm])
    hero.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#123B50")),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor("#123B50")),
        ("LEFTPADDING", (0, 0), (-1, -1), 13),
        ("RIGHTPADDING", (0, 0), (-1, -1), 13),
        ("TOPPADDING", (0, 0), (-1, 0), 10),
        ("TOPPADDING", (0, 1), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, 1), 2),
        ("BOTTOMPADDING", (0, 2), (-1, 2), 10),
    ]))

    boundary_text = (
        "<b>High ≠ 已确证　　分数 ≠ 存在概率</b><br/>" if zh
        else "<b>High ≠ confirmed　　Score ≠ presence probability</b><br/>"
    ) + escape(payload["read_first"])
    boundary = Table([[Paragraph(boundary_text, card)]], colWidths=[17.5 * cm])
    boundary.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FFF8E8")),
        ("BOX", (0, 0), (-1, -1), 0.55, colors.HexColor("#E8C77A")),
        ("LINEBEFORE", (0, 0), (0, -1), 4, colors.HexColor("#D5A72D")),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
    ]))

    routes = (
        [
            ("01", "先看证据等级", "High／Medium／Low 回答证据强弱"),
            ("02", "再看身份范围", "二维连接、异构体组与类别不可混用"),
            ("03", "最后看正式状态", "只有已确证／已发布才可确定性表述"),
        ]
        if zh else
        [
            ("01", "Read the evidence tier", "High / Medium / Low describes evidence strength"),
            ("02", "Check identity scope", "2-D identity, isomer group and class differ"),
            ("03", "Check formal status", "Use definitive wording only after release"),
        ]
    )
    route_cells = [
        Paragraph(
            f"<font color='#258078'><b>{number}</b></font>　<b>{escape(title)}</b>"
            f"<br/><font color='#71828D'>{escape(note)}</font>",
            card,
        )
        for number, title, note in routes
    ]
    route_table = Table([route_cells], colWidths=[5.72 * cm] * 3, hAlign="LEFT")
    route_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F6F9FA")),
        ("BOX", (0, 0), (-1, -1), 0.45, colors.HexColor("#D7E3E8")),
        ("INNERGRID", (0, 0), (-1, -1), 0.45, colors.HexColor("#D7E3E8")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
    ]))

    tier_cells: list[Any] = []
    for item in payload["levels"]:
        tone = tone_colors.get(str(item.get("tone") or "unknown"), "#748191")
        text = (
            f"<font color='{tone}'><b>{escape(item['label'])}</b></font><br/>"
            f"{escape(item['meaning'])}<br/>"
            f"<font color='#24705F'><b>{'可用：' if zh else 'Use: '}</b>{escape(item['use'])}</font><br/>"
            f"<font color='#9A4C2D'><b>{'边界：' if zh else 'Boundary: '}</b>{escape(item['avoid'])}</font><br/>"
            f"<font color='#536F82'><b>{'下一步：' if zh else 'Next: '}</b>{escape(item['next_step'])}</font>"
        )
        tier_cells.append(Paragraph(text, card))
    tier_rows = [tier_cells[index:index + 2] for index in range(0, len(tier_cells), 2)]
    tier_table = Table(tier_rows, colWidths=[8.62 * cm, 8.62 * cm], hAlign="LEFT")
    tier_style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BOX", (0, 0), (-1, -1), 0.45, colors.HexColor("#D9E3E8")),
        ("INNERGRID", (0, 0), (-1, -1), 0.45, colors.HexColor("#D9E3E8")),
        ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.white, colors.HexColor("#F8FAFB")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
    ]
    for index, item in enumerate(payload["levels"]):
        row_index, column_index = divmod(index, 2)
        tier_style.append((
            "LINEABOVE", (column_index, row_index), (column_index, row_index), 2.2,
            colors.HexColor(tone_colors.get(str(item.get("tone") or "unknown"), "#748191")),
        ))
    tier_table.setStyle(TableStyle(tier_style))

    first_page = [
        hero, Spacer(1, 0.18 * cm), boundary, Spacer(1, 0.15 * cm),
        Paragraph(escape(payload["threshold_note"]), small),
        Spacer(1, 0.18 * cm), route_table, Spacer(1, 0.25 * cm),
        Paragraph("01 · " + ("证据等级速查" if zh else "EVIDENCE-TIER QUICK REFERENCE"), section),
        tier_table,
    ]

    scope_cells = []
    for index, item in enumerate(payload["identity_scopes"], start=1):
        scope_cells.append(Paragraph(
            f"<font color='#258078'><b>{index:02d}</b></font>　<b>{escape(item['label'])}</b><br/>"
            f"{escape(item['meaning'])}<br/><font color='#24705F'>"
            f"<b>{'建议表述：' if zh else 'Recommended wording: '}</b>{escape(item['report'])}</font>",
            card,
        ))
    scope_rows = [scope_cells[index:index + 2] for index in range(0, len(scope_cells), 2)]
    scope_table = Table(scope_rows, colWidths=[8.62 * cm, 8.62 * cm], hAlign="LEFT")
    scope_table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BOX", (0, 0), (-1, -1), 0.45, colors.HexColor("#D6E3E7")),
        ("INNERGRID", (0, 0), (-1, -1), 0.45, colors.HexColor("#D6E3E7")),
        ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.HexColor("#F7FAFA"), colors.white]),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
    ]))

    score_cells = []
    for index, item in enumerate(payload["scores"], start=1):
        score_cells.append(Paragraph(
            f"<font color='#476F82'><b>{index:02d}</b></font>　<b>{escape(item['label'])}</b><br/>"
            f"{escape(item['meaning'])}<br/><font color='#9A4C2D'>"
            f"<b>{'误读边界：' if zh else 'Boundary: '}</b>{escape(item['boundary'])}</font>",
            small,
        ))
    if len(score_cells) % 2:
        score_cells.append(Paragraph("", small))
    score_rows = [score_cells[index:index + 2] for index in range(0, len(score_cells), 2)]
    score_table = Table(score_rows, colWidths=[8.62 * cm, 8.62 * cm], hAlign="LEFT")
    score_table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BOX", (0, 0), (-1, -1), 0.4, colors.HexColor("#DCE5E9")),
        ("INNERGRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#DCE5E9")),
        ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.white, colors.HexColor("#F8FAFB")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    second_page = [
        Paragraph("02 · " + ("身份范围" if zh else "IDENTITY SCOPE"), section),
        scope_table, Spacer(1, 0.25 * cm),
        Paragraph("03 · " + ("各项评分与辅助字段" if zh else "SCORES AND AUXILIARY FIELDS"), section),
        score_table, Spacer(1, 0.15 * cm),
        Paragraph(
            f"{escape(payload['policy_id'])} · {escape(payload['schema_version'])}",
            policy,
        ),
        Spacer(1, 0.25 * cm),
    ]
    return [KeepTogether(first_page), PageBreak(), *second_page]


def format_confidence_guide_worksheet(worksheet: Any, language: str = "zh") -> None:
    """Replace the flat first Excel sheet with a readable report-style guide."""

    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.worksheet.page import PageMargins

    payload = confidence_guide_payload(language)
    zh = payload["language"] == "zh"
    worksheet.delete_rows(1, worksheet.max_row)
    for merged_range in list(worksheet.merged_cells.ranges):
        worksheet.unmerge_cells(str(merged_range))
    for row_dimension in worksheet.row_dimensions.values():
        row_dimension.height = None
    worksheet.sheet_format.defaultRowHeight = 18

    navy = "123B50"
    teal = "16746F"
    ink = "294557"
    muted = "6D7F8B"
    line = "D8E3E8"
    soft = "F6F9FA"
    amber = "FFF4D9"
    amber_line = "D6A62D"
    rust = "994B2E"
    white = "FFFFFF"
    tone_colors = {
        "confirmed": "147D64", "high": "2474B8", "medium": "C28A20",
        "low": "BC5A31", "background": "7456A6", "unknown": "748191",
    }
    thin = Side(style="thin", color=line)
    card_border = Border(left=thin, right=thin, top=thin, bottom=thin)

    def merged_value(
        cell_range: str,
        value: str,
        *,
        fill: str = white,
        font: Font | None = None,
        alignment: Alignment | None = None,
        border: Border | None = card_border,
    ) -> None:
        worksheet.merge_cells(cell_range)
        cell = worksheet[cell_range.split(":", 1)[0]]
        cell.value = value
        cell.fill = PatternFill("solid", fgColor=fill)
        cell.font = font or Font(name="Aptos", size=10, color=ink)
        cell.alignment = alignment or Alignment(
            horizontal="left", vertical="center", wrap_text=True,
        )
        if border:
            cell.border = border

    for column, width in zip("ABCDEF", (15, 23, 15, 23, 15, 23)):
        worksheet.column_dimensions[column].width = width

    merged_value(
        "A1:F1", payload["title"], fill=navy,
        font=Font(name="Aptos Display", size=18, bold=True, color=white),
        alignment=Alignment(horizontal="left", vertical="bottom"), border=None,
    )
    merged_value(
        "A2:F2", "REPORT READING GUIDE  ·  " + payload["subtitle"], fill=navy,
        font=Font(name="Aptos", size=9, color="D9EAEB"),
        alignment=Alignment(horizontal="left", vertical="top", wrap_text=True), border=None,
    )
    worksheet.row_dimensions[1].height = 34
    worksheet.row_dimensions[2].height = 29
    worksheet.row_dimensions[3].height = 8

    boundary_title = (
        "先记住：High ≠ 已确证　｜　分数 ≠ 存在概率"
        if zh else "Remember: High ≠ confirmed　|　Score ≠ presence probability"
    )
    merged_value(
        "A4:F4", f"{boundary_title}\n{payload['read_first']}", fill=amber,
        font=Font(name="Aptos", size=10, bold=True, color="6B4B13"),
        alignment=Alignment(horizontal="left", vertical="center", wrap_text=True),
        border=Border(left=Side(style="thick", color=amber_line), bottom=thin),
    )
    worksheet.row_dimensions[4].height = 52
    merged_value(
        "A5:F5", payload["threshold_note"], fill="EDF6F4",
        font=Font(name="Aptos", size=9, color="4D6970"),
        alignment=Alignment(horizontal="left", vertical="center", wrap_text=True),
    )
    worksheet.row_dimensions[5].height = 46
    worksheet.row_dimensions[6].height = 8

    def section_band(row: int, number: str, title: str, note: str = "") -> None:
        merged_value(
            f"A{row}:F{row}", f"{number}  {title}" + (f"   ·   {note}" if note else ""),
            fill=navy, font=Font(name="Aptos", size=10, bold=True, color=white),
            alignment=Alignment(horizontal="left", vertical="center"), border=None,
        )
        worksheet.row_dimensions[row].height = 24

    section_band(7, "01", "阅读路径" if zh else "READING PATH")
    routes = (
        [
            ("先看证据等级", "High／Medium／Low 回答证据强弱"),
            ("再看身份范围", "二维连接、异构体组与类别不可混用"),
            ("最后看正式状态", "只有已确证／已发布才可确定性表述"),
        ]
        if zh else
        [
            ("Read the evidence tier", "High / Medium / Low describes evidence strength"),
            ("Check identity scope", "2-D identity, isomer group and class differ"),
            ("Check formal status", "Use definitive wording only after release"),
        ]
    )
    for index, ((title, note), columns) in enumerate(zip(routes, ("AB", "CD", "EF")), start=1):
        merged_value(
            f"{columns[0]}8:{columns[1]}8", f"{index:02d}  {title}\n{note}", fill=soft,
            font=Font(name="Aptos", size=9, bold=True, color=ink),
            alignment=Alignment(horizontal="left", vertical="center", wrap_text=True),
        )
    worksheet.row_dimensions[8].height = 48
    worksheet.row_dimensions[9].height = 8

    section_band(
        10, "02", "证据等级速查" if zh else "EVIDENCE-TIER QUICK REFERENCE",
        "先定位等级，再看使用边界" if zh else "Locate the tier, then read its boundary",
    )
    pairs = ("AB", "CD", "EF")
    for index, item in enumerate(payload["levels"]):
        card_group = index // 3
        columns = pairs[index % 3]
        label_row = 11 + card_group * 2
        detail_row = label_row + 1
        tone = tone_colors.get(str(item.get("tone") or "unknown"), "748191")
        merged_value(
            f"{columns[0]}{label_row}:{columns[1]}{label_row}", item["label"], fill=tone,
            font=Font(name="Aptos", size=10, bold=True, color=white),
            alignment=Alignment(horizontal="left", vertical="center"), border=None,
        )
        detail = (
            f"{'代表' if zh else 'Meaning'}｜{item['meaning']}\n"
            f"{'可用' if zh else 'Use'}｜{item['use']}\n"
            f"{'边界' if zh else 'Boundary'}｜{item['avoid']}\n"
            f"{'下一步' if zh else 'Next'}｜{item['next_step']}"
        )
        merged_value(
            f"{columns[0]}{detail_row}:{columns[1]}{detail_row}", detail, fill="FBFCFD",
            font=Font(name="Aptos", size=8, color=ink),
            alignment=Alignment(horizontal="left", vertical="top", wrap_text=True),
        )
        worksheet.row_dimensions[label_row].height = 24
        worksheet.row_dimensions[detail_row].height = 104
    worksheet.row_dimensions[15].height = 8

    section_band(16, "03", "身份范围" if zh else "IDENTITY SCOPE")
    for offset, item in enumerate(payload["identity_scopes"], start=17):
        merged_value(
            f"A{offset}:B{offset}", item["label"], fill="E7F2F1",
            font=Font(name="Aptos", size=9, bold=True, color=teal),
            alignment=Alignment(horizontal="left", vertical="center", wrap_text=True),
        )
        merged_value(
            f"C{offset}:F{offset}",
            f"{item['meaning']}\n{'建议表述' if zh else 'Recommended wording'}｜{item['report']}",
            fill=white, font=Font(name="Aptos", size=8.5, color=ink),
            alignment=Alignment(horizontal="left", vertical="center", wrap_text=True),
        )
        worksheet.row_dimensions[offset].height = 47

    section_band(22, "04", "各项评分与辅助字段" if zh else "SCORES AND AUXILIARY FIELDS")
    score_rows = 6
    for index in range(score_rows):
        row = 23 + index
        for item_index, columns in ((index, "ABC"), (index + score_rows, "DEF")):
            if item_index >= len(payload["scores"]):
                continue
            item = payload["scores"][item_index]
            value = (
                f"{item_index + 1:02d}  {item['label']}\n{item['meaning']}\n"
                f"{'误读边界' if zh else 'Boundary'}｜{item['boundary']}"
            )
            merged_value(
                f"{columns[0]}{row}:{columns[-1]}{row}", value,
                fill=soft if index % 2 == 0 else white,
                font=Font(name="Aptos", size=8, color=ink),
                alignment=Alignment(horizontal="left", vertical="center", wrap_text=True),
            )
            worksheet[f"{columns[0]}{row}"].border = Border(
                left=Side(style="medium", color=teal), right=thin, top=thin, bottom=thin,
            )
        worksheet.row_dimensions[row].height = 70
    worksheet.row_dimensions[29].height = 8

    merged_value(
        "A30:F30", f"{payload['policy_id']}  ·  {payload['schema_version']}",
        fill=white, font=Font(name="Aptos", size=7, color="8B99A3"),
        alignment=Alignment(horizontal="right", vertical="center"), border=None,
    )
    worksheet.row_dimensions[30].height = 20
    worksheet.sheet_view.showGridLines = False
    worksheet.sheet_view.zoomScale = 85
    worksheet.freeze_panes = "A7"
    worksheet.sheet_properties.tabColor = teal
    worksheet.sheet_properties.pageSetUpPr.fitToPage = True
    worksheet.page_setup.orientation = "landscape"
    worksheet.page_setup.paperSize = worksheet.PAPERSIZE_A4
    worksheet.page_setup.fitToWidth = 1
    worksheet.page_setup.fitToHeight = 2
    worksheet.print_area = "A1:F30"
    worksheet.print_title_rows = None
    worksheet.page_margins = PageMargins(
        left=0.28, right=0.28, top=0.35, bottom=0.35, header=0.15, footer=0.15,
    )
    worksheet.oddFooter.center.text = "MSsystem · " + (
        "置信度解释规则" if zh else "Confidence interpretation policy"
    )
    worksheet.oddFooter.center.size = 8
    worksheet.oddFooter.center.color = "7D8A93"
