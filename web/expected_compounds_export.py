"""Export helpers for the independent Expected Compounds evidence layer.

Expected Compounds is a continuous-MS1 targeted XIC workflow.  Its rows must
remain visibly separate from ordinary MS2 spectral-library matches in every
report format.  This module provides one column contract and one set of
scientific boundary notes for CSV, Excel, HTML, and PDF exports.
"""

from __future__ import annotations

import hashlib
import html
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from matching.compound_names import resolve_compound_name


EVIDENCE_LABELS = {
    "zh": {
        "MS1_EXACT_MASS_CANDIDATE": "MS1 精确质量候选",
        "MS1_ISOTOPE_SUPPORTED": "MS1 + 同位素支持",
        "MS1_MS2_DIAGNOSTIC_PROBABLE": "MS1 + 诊断离子筛查支持",
        "EXPERIMENTAL_LIBRARY_CONFIRMED": "MS/MS 实验谱库候选支持",
    },
    "en": {
        "MS1_EXACT_MASS_CANDIDATE": "MS1 exact-mass candidate",
        "MS1_ISOTOPE_SUPPORTED": "MS1 + isotope supported",
        "MS1_MS2_DIAGNOSTIC_PROBABLE": "MS1 + diagnostic-ion screening support",
        "EXPERIMENTAL_LIBRARY_CONFIRMED": "MS/MS experimental-library candidate support",
    },
}

STATUS_LABELS = {
    "zh": {
        "ok": "已运行",
        "disabled": "本次分析已关闭",
        "not_run": "旧任务未运行",
        "failed": "运行失败",
        "unavailable_requires_mzml_with_ms1": "不可用：缺少连续 MS1",
        "unavailable_unknown_or_mixed_polarity": "不可用：极性未解析",
        "unavailable_no_selected_target_dataset": (
            "不可用：未绑定样本专属靶向库"
        ),
        "unavailable_no_presence_only_targets": "不可用：没有文献靶标",
        "unavailable_no_valid_expected_ions": "不可用：没有可搜索离子",
        "unavailable_no_ms1_scans": "不可用：mzML 中没有 survey MS1",
        "unavailable_no_usable_ms1_scans": (
            "不可用：mzML 中没有可用 survey MS1"
        ),
    },
    "en": {
        "ok": "Completed",
        "disabled": "Disabled for this analysis",
        "not_run": "Not run for this legacy task",
        "failed": "Failed",
        "unavailable_requires_mzml_with_ms1": (
            "Unavailable: continuous MS1 is absent"
        ),
        "unavailable_unknown_or_mixed_polarity": (
            "Unavailable: polarity is unresolved"
        ),
        "unavailable_no_selected_target_dataset": (
            "Unavailable: no sample-specific target catalogue was bound"
        ),
        "unavailable_no_presence_only_targets": (
            "Unavailable: no literature targets"
        ),
        "unavailable_no_valid_expected_ions": (
            "Unavailable: no searchable expected ions"
        ),
        "unavailable_no_ms1_scans": (
            "Unavailable: no survey MS1 scans in mzML"
        ),
        "unavailable_no_usable_ms1_scans": (
            "Unavailable: no usable survey MS1 scans in mzML"
        ),
    },
}

# Standard export contract: the main targeted table answers four questions only:
# what candidate, which target database, where the feature is, and how strong
# the evidence is.  Internal IDs, boolean flags and isomer-engine diagnostics
# remain available in the JSON/HTML audit layer rather than the routine table.
EXPECTED_FEATURE_COLUMNS = [
    "Evidence",
    "Candidate_Name",
    "Candidate_Count",
    "Target_Database",
    "Formula",
    "Adduct",
    "Related_Adducts",
    "Polarity",
    "Target_mz",
    "Observed_mz",
    "Mass_Error_ppm",
    "RT_min",
    "Peak_Area",
    "Peak_Quality",
    "OpenMS_Verification",
    "Background_Status",
    "Raw_Peak_Count",
    "Signal_to_Noise",
    "MS2_Diagnostic_Support",
    "Experimental_Library_Support",
]

EXPECTED_FEATURE_LABELS = {
    "zh": {
        "Evidence": "证据等级",
        "Candidate_Name": "主要候选",
        "Candidate_Count": "候选数",
        "Target_Database": "靶向数据库",
        "Formula": "分子式",
        "Adduct": "代表加合物",
        "Related_Adducts": "共洗脱加合物",
        "Polarity": "离子模式",
        "Target_mz": "理论 m/z",
        "Observed_mz": "实测 m/z",
        "Mass_Error_ppm": "误差 (ppm)",
        "RT_min": "主峰 RT (min)",
        "Peak_Area": "代表峰面积",
        "Peak_Quality": "峰质量 (0–10)",
        "OpenMS_Verification": "OpenMS 峰验证",
        "Background_Status": "背景状态",
        "Raw_Peak_Count": "原始峰段数",
        "Signal_to_Noise": "信噪比",
        "MS2_Diagnostic_Support": "诊断 MS2 谱图数",
        "Experimental_Library_Support": "实验谱库谱图数",
    },
    "en": {
        column: column for column in EXPECTED_FEATURE_COLUMNS
    },
}

# Full traceability export.  Unlike EXPECTED_FEATURE_COLUMNS this contract is
# deliberately one row per accepted XIC peak segment and exposes peak
# boundaries/quality checks.  It is a separate download so routine users are
# never asked to interpret internal peak-engine rows as compounds.
EXPECTED_RAW_FEATURE_COLUMNS = [
    "Feature_ID",
    "Evidence",
    "Candidate_Name",
    "Candidate_Count",
    "Target_Database",
    "Formula",
    "Adduct",
    "Polarity",
    "Target_mz",
    "Observed_mz",
    "Mass_Error_ppm",
    "Apex_RT_min",
    "Boundary_Start_RT_min",
    "Boundary_End_RT_min",
    "FWHM_s",
    "Peak_Height",
    "Peak_Area",
    "Signal_to_Noise",
    "Scan_Points",
    "Gap_Ratio",
    "Peak_Quality",
    "OpenMS_Verification",
    "Background_Status",
    "MS2_Diagnostic_Support",
    "Experimental_Library_Support",
]

EXPECTED_RAW_FEATURE_LABELS = {
    "zh": {
        "Feature_ID": "靶向峰 ID",
        "Evidence": "证据等级",
        "Candidate_Name": "候选名称",
        "Candidate_Count": "候选数",
        "Target_Database": "靶向数据库",
        "Formula": "分子式",
        "Adduct": "加合物",
        "Polarity": "离子模式",
        "Target_mz": "理论 m/z",
        "Observed_mz": "实测 m/z",
        "Mass_Error_ppm": "误差 (ppm)",
        "Apex_RT_min": "顶点 RT (min)",
        "Boundary_Start_RT_min": "起始 RT (min)",
        "Boundary_End_RT_min": "结束 RT (min)",
        "FWHM_s": "FWHM (s)",
        "Peak_Height": "峰高",
        "Peak_Area": "峰面积",
        "Signal_to_Noise": "信噪比",
        "Scan_Points": "有效扫描点",
        "Gap_Ratio": "缺口比例",
        "Peak_Quality": "峰质量 (0–10)",
        "OpenMS_Verification": "OpenMS 峰验证",
        "Background_Status": "背景状态",
        "MS2_Diagnostic_Support": "诊断 MS2 谱图数",
        "Experimental_Library_Support": "实验谱库谱图数",
    },
    "en": {
        column: column for column in EXPECTED_RAW_FEATURE_COLUMNS
    },
}

EXPECTED_SUMMARY_COLUMNS = [
    "Status_Label",
    "Selected_Datasets",
    "Target_Catalog_Rows",
    "Active_Target_Rows",
    "Quarantined_Target_Rows",
    "Expected_Ion_Hypotheses",
    "Detected_Features",
    "Supported_Features",
    "Mass_Only_Features",
    "MS1_MS2_Diagnostic_Probable",
    "Experimental_Library_Confirmed",
    "Warning",
]

ANALYSIS_WARNING_COLUMNS = [
    "Severity",
    "Modules",
    "Code",
    "Message",
    "Action",
]

ANALYSIS_WARNING_LABELS = {
    "zh": ["级别", "适用模块", "代码", "状态说明", "建议操作"],
    "en": ANALYSIS_WARNING_COLUMNS,
}


def _number(value: Any) -> float | int | None:
    if isinstance(value, bool):
        return int(value)
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed):
        return None
    return int(parsed) if parsed.is_integer() else parsed


def _join(values: Any, separator: str = " | ") -> str:
    if values is None:
        return ""
    if isinstance(values, str):
        return values
    if isinstance(values, Mapping):
        values = values.values()
    if not isinstance(values, Sequence):
        return str(values)
    return separator.join(
        str(value).strip()
        for value in values
        if str(value).strip()
    )


def _unique(values: Sequence[Any]) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        output.append(text)
    return output


def _display_sources(values: Any) -> str:
    """Make internal dataset tokens readable in a user-facing export."""
    displayed = []
    source_values = [values] if isinstance(values, str) else list(values or [])
    for value in _unique(source_values):
        text = value
        for prefix in ("user_contrib:", "literature:", "database:"):
            if text.startswith(prefix):
                text = text[len(prefix):]
                break
        displayed.append(text)
    return _join(displayed)


def _result_note(
    level: str,
    candidate_count: int,
    language: str,
) -> str:
    notes = {
        "zh": {
            "MS1_EXACT_MASS_CANDIDATE": "精确质量候选，需进一步证据确认",
            "MS1_ISOTOPE_SUPPORTED": "精确质量并有同位素证据支持",
            "MS1_MS2_DIAGNOSTIC_PROBABLE": "诊断离子筛查支持",
            "EXPERIMENTAL_LIBRARY_CONFIRMED": "MS/MS 实验谱库候选支持",
            "multiple": "同一离子特征对应多个数据库候选",
        },
        "en": {
            "MS1_EXACT_MASS_CANDIDATE": (
                "Exact-mass candidate; additional evidence is required"
            ),
            "MS1_ISOTOPE_SUPPORTED": (
                "Exact mass with supporting isotope evidence"
            ),
            "MS1_MS2_DIAGNOSTIC_PROBABLE": (
                "Screening support from reported MS2 diagnostic ions"
            ),
            "EXPERIMENTAL_LIBRARY_CONFIRMED": (
                "Supported as an MS/MS experimental-library candidate; not standard-confirmed"
            ),
            "multiple": (
                "One ion feature maps to multiple database candidates"
            ),
        },
    }["en" if language == "en" else "zh"]
    parts = [notes.get(level, level)]
    if candidate_count > 1:
        parts.append(notes["multiple"])
    return "；".join(parts) if language != "en" else "; ".join(parts)


def _identity_key(member: Mapping[str, Any]) -> str:
    """Return a full-structure identity key without collapsing stereoisomers."""

    stored = str(member.get("identity_key") or "").strip()
    if stored:
        return stored
    inchikey = str(member.get("inchikey") or "").strip().upper()
    if inchikey:
        return f"ik::{inchikey}"
    target_id = str(
        member.get("target_id") or member.get("catalog_row_id") or ""
    ).strip()
    if target_id:
        return f"target::{target_id}"
    return "name::" + "|".join(
        (
            str(member.get("compound_name") or "").strip().casefold(),
            str(member.get("formula") or "").strip(),
        )
    )


def _candidate_display(
    names: Sequence[Any],
    identity_count: int,
    language: str,
) -> str:
    """Return one normal compound name for routine reader-facing outputs.

    Candidate multiplicity remains in ``Candidate_Count`` and all names remain
    available in the detailed audit fields.
    """

    clean_names = _unique(list(names))
    if not clean_names:
        return "未命名靶标" if language != "en" else "Unnamed target"
    return clean_names[0]


def _clean_member_candidate_names(
    members: Sequence[Mapping[str, Any]],
    fallback_names: Sequence[Any] | None = None,
) -> list[str]:
    """Return reader-safe chemical names, preferring exact structure aliases.

    Presence-only catalogues can retain opaque source identifiers such as
    ``NP-000004`` for provenance.  Those identifiers remain in the target
    member audit, but they must not be rendered as candidate compound names.
    Exact InChIKey aliases are preferred because they can recover a real name
    (for example Schaftoside) from an opaque catalogue label.
    """

    ranked: list[tuple[int, int, str]] = []
    for index, member in enumerate(members):
        raw_name = str(member.get("compound_name") or "").strip()
        display, status = resolve_compound_name(
            raw_name, member.get("inchikey")
        )
        # Synthetic fixtures and user-authored catalogues sometimes use
        # ``Candidate A``/``Candidate B`` as an intentional reader label.  Keep
        # that exact label rather than rewriting it to a generic computational
        # placeholder; opaque registry/accession codes remain excluded below.
        if status == "computational_candidate" and re.fullmatch(
            r"Candidate\s+[A-Za-z0-9]+", raw_name, re.IGNORECASE
        ):
            display = raw_name
            status = "source_name"
        if (
            not display
            or status == "missing"
            or status.startswith("unsafe_")
            or display.casefold().startswith("unresolved compound")
        ):
            continue
        priority = 0 if status == "exact_inchikey_alias" else 1
        ranked.append((priority, index, display))

    if not ranked:
        for index, raw_name in enumerate(fallback_names or []):
            display, status = resolve_compound_name(raw_name, "")
            if (
                not display
                or status != "source_name"
                or display.casefold().startswith("unresolved compound")
            ):
                continue
            ranked.append((2, index, display))

    ranked.sort(key=lambda item: (item[0], item[1], item[2].casefold()))
    return _unique([item[2] for item in ranked])


def _normalized_analysis_warnings(
    warnings: Sequence[Any] | None,
) -> list[dict[str, str]]:
    """Return unique report-safe module warnings."""

    output: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for value in warnings or []:
        if isinstance(value, Mapping):
            modules = _join(value.get("modules") or [], " / ")
            record = {
                "severity": str(value.get("severity") or "warning"),
                "modules": modules,
                "code": str(value.get("code") or ""),
                "message": str(value.get("message") or ""),
                "action": str(value.get("action") or ""),
            }
        else:
            record = {
                "severity": "warning",
                "modules": "",
                "code": "",
                "message": str(value or ""),
                "action": "",
            }
        key = (
            record["code"],
            record["modules"],
            record["message"],
            record["action"],
        )
        if not record["message"] or key in seen:
            continue
        seen.add(key)
        output.append(record)
    return output


def analysis_warning_rows(
    warnings: Sequence[Any] | None,
) -> list[dict[str, str]]:
    """Return one row per unique module warning for audit exports."""

    return [
        {
            "Severity": warning["severity"],
            "Modules": warning["modules"],
            "Code": warning["code"],
            "Message": warning["message"],
            "Action": warning["action"],
        }
        for warning in _normalized_analysis_warnings(warnings)
    ]


def reportlab_font_names(language: str = "zh") -> tuple[str, str]:
    """Return ReportLab fonts that render Chinese and Latin text visibly."""

    if language == "en":
        # Do not rely on ReportLab's unembedded Type-1 Helvetica.  On the
        # production macOS/Poppler path it remains extractable but can render
        # as completely invisible text, leaving only table borders.  Embed a
        # TrueType face so PDF appearance is self-contained and reproducible.
        from reportlab import __file__ as reportlab_module_path
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont

        regular_name = "MSsystemLatin"
        bold_name = "MSsystemLatinBold"
        if (
            regular_name not in pdfmetrics.getRegisteredFontNames()
            or bold_name not in pdfmetrics.getRegisteredFontNames()
        ):
            bundled_fonts = Path(reportlab_module_path).resolve().parent / "fonts"
            candidates = [
                (
                    Path("/System/Library/Fonts/Supplemental/Arial.ttf"),
                    Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf"),
                ),
                (
                    Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
                    Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
                ),
                (
                    Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
                    Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
                ),
                (bundled_fonts / "Vera.ttf", bundled_fonts / "VeraBd.ttf"),
            ]
            font_pair = next(
                (
                    (regular_path, bold_path)
                    for regular_path, bold_path in candidates
                    if regular_path.is_file() and bold_path.is_file()
                ),
                None,
            )
            if font_pair is None:
                # This is an export-integrity failure, not a cosmetic fallback:
                # generating a potentially blank PDF would be misleading.
                raise RuntimeError(
                    "No embeddable TrueType font is available for PDF export"
                )
            regular_path, bold_path = font_pair
            if regular_name not in pdfmetrics.getRegisteredFontNames():
                pdfmetrics.registerFont(
                    TTFont(regular_name, str(regular_path))
                )
            if bold_name not in pdfmetrics.getRegisteredFontNames():
                pdfmetrics.registerFont(TTFont(bold_name, str(bold_path)))
            pdfmetrics.registerFontFamily(
                regular_name,
                normal=regular_name,
                bold=bold_name,
                italic=regular_name,
                boldItalic=bold_name,
            )
        return regular_name, bold_name

    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont

    font_name = "MSsystemUnicode"
    if font_name not in pdfmetrics.getRegisteredFontNames():
        candidates = [
            Path("/System/Library/Fonts/Supplemental/Arial Unicode.ttf"),
            Path("/Library/Fonts/Arial Unicode.ttf"),
            Path("/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttf"),
            Path("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"),
        ]
        font_path = next(
            (candidate for candidate in candidates if candidate.exists()),
            None,
        )
        if font_path is not None:
            pdfmetrics.registerFont(
                TTFont(font_name, str(font_path))
            )
            pdfmetrics.registerFontFamily(
                font_name,
                normal=font_name,
                bold=font_name,
                italic=font_name,
                boldItalic=font_name,
            )
        else:
            from reportlab.pdfbase.cidfonts import UnicodeCIDFont

            font_name = "STSong-Light"
            if font_name not in pdfmetrics.getRegisteredFontNames():
                pdfmetrics.registerFont(UnicodeCIDFont(font_name))
    return font_name, font_name


def _deduplicate_members(
    members: Sequence[Mapping[str, Any]],
) -> list[Mapping[str, Any]]:
    unique: dict[str, Mapping[str, Any]] = {}
    for member in members:
        key = _identity_key(member)
        if key and key not in unique:
            unique[key] = member
    return list(unique.values())


def _selected_target_members(
    feature: Mapping[str, Any],
) -> list[Mapping[str, Any]]:
    """Narrow an ambiguity group only when measured MS2 resolves a member."""

    members = _deduplicate_members([
        member
        for member in feature.get("target_members") or []
        if isinstance(member, Mapping)
    ])
    experimental = [
        item
        for item in feature.get("experimental_library_support") or []
        if isinstance(item, Mapping)
    ]
    resolving_support = experimental or [
        item
        for item in feature.get("ms2_diagnostic_support") or []
        if isinstance(item, Mapping)
    ]
    supported_target_ids = {
        str(item.get("target_id") or "").strip()
        for item in resolving_support
        if str(item.get("target_id") or "").strip()
    }
    supported_inchikeys = {
        str(item.get("target_inchikey") or "").strip().upper()
        for item in experimental
        if str(item.get("target_inchikey") or "").strip()
    }
    if not supported_target_ids and not supported_inchikeys:
        return members
    selected = [
        member
        for member in members
        if (
            str(member.get("target_id") or "").strip()
            in supported_target_ids
            or str(member.get("inchikey") or "").strip().upper()
            in supported_inchikeys
        )
    ]
    return selected or members


def _openms_label(feature: Mapping[str, Any], language: str) -> str:
    status = str(feature.get("openms_feature_link_status") or "").strip()
    labels = {
        "zh": {
            "linked_openms_feature": "已由独立 OpenMS 峰验证",
            "ambiguous_openms_features": "已关联多个 OpenMS 峰",
            "no_openms_feature_in_mz_rt_window": "未关联到 OpenMS 峰",
            "openms_features_unavailable": "OpenMS 峰表不可用",
            "target_feature_missing_mz_or_rt": "缺少 m/z 或 RT",
            "": "旧任务未计算",
        },
        "en": {
            "linked_openms_feature": "Linked to independent OpenMS feature",
            "ambiguous_openms_features": "Multiple OpenMS features linked",
            "no_openms_feature_in_mz_rt_window": "No OpenMS feature linked",
            "openms_features_unavailable": "OpenMS feature table unavailable",
            "target_feature_missing_mz_or_rt": "Target m/z or RT missing",
            "": "Not calculated for legacy task",
        },
    }[language]
    return labels.get(status, status)


def _background_label(feature: Mapping[str, Any], language: str) -> str:
    status = str(feature.get("blank_status") or "").strip()
    category = str(feature.get("blank_evidence_category") or "").strip()
    labels = {
        "zh": {
            "not_requested": "",
            "no_blank_feature": "空白中未检出",
            "sample_enriched": "样本富集",
            "blank_associated": "空白相关",
            "blank_dominant": "背景主导",
            "low_blank_prevalence": "空白中低频检出",
            "unevaluable": "空白证据不可判定",
            "": "",
        },
        "en": {
            "not_requested": "",
            "no_blank_feature": "Not detected in blank",
            "sample_enriched": "Sample enriched",
            "blank_associated": "Blank associated",
            "blank_dominant": "Background dominant",
            "low_blank_prevalence": "Low blank prevalence",
            "unevaluable": "Blank evidence unevaluable",
            "": "",
        },
    }[language]
    if category:
        return labels.get(category, category)
    return labels.get(status, status or labels[""])


def _unique_support_count(
    support: Sequence[Mapping[str, Any]],
    *,
    library: bool,
) -> int:
    keys = set()
    for item in support:
        spectrum_id = str(item.get("spectrum_id") or "").strip()
        library_id = (
            str(item.get("library_spectrum_id") or "").strip()
            if library else ""
        )
        keys.add((spectrum_id, library_id))
    return len(keys)


def status_label(payload: Mapping[str, Any], language: str = "zh") -> str:
    language = "en" if language == "en" else "zh"
    summary = payload.get("summary") or {}
    status = str(
        payload.get("status") or summary.get("status") or "not_run"
    )
    return STATUS_LABELS[language].get(status, status)


def expected_summary_row(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> dict[str, Any]:
    summary = payload.get("summary") or {}
    catalog = payload.get("catalog") or {}
    evidence = summary.get("evidence_counts") or {}
    status = str(
        payload.get("status") or summary.get("status") or "not_run"
    )
    selected = (
        summary.get("selected_dataset_sources")
        or catalog.get("selected_dataset_sources")
        or catalog.get("dataset_sources")
        or []
    )
    compound_rows = _expected_detail_rows(payload, language)
    compound_evidence_counts = {
        level: sum(
            str(row.get("Evidence_Level") or "") == level
            for row in compound_rows
        )
        for level in EVIDENCE_LABELS["en"]
    }
    quality_audit = (
        summary.get("target_catalog_quality_audit")
        or catalog.get("quality_audit")
        or {}
    )
    target_catalog_rows = int(
        summary.get("target_catalog_rows")
        or catalog.get("target_rows")
        or 0
    )
    quarantined_target_rows = int(
        summary.get("quarantined_target_rows")
        or quality_audit.get("quarantined_target_rows")
        or 0
    )
    active_target_rows = int(
        summary.get("active_target_rows")
        or quality_audit.get("active_target_rows")
        or max(0, target_catalog_rows - quarantined_target_rows)
    )
    return {
        "Status": status,
        "Status_Label": status_label(payload, language),
        "Requested": bool(payload.get("requested", summary.get("requested"))),
        "Method": str(
            payload.get("method")
            or "Expected Compounds / continuous-MS1 targeted XIC"
        ),
        "Target_Catalog_Rows": target_catalog_rows,
        "Active_Target_Rows": active_target_rows,
        "Quarantined_Target_Rows": quarantined_target_rows,
        "Target_Catalog_Quality_Audit": dict(quality_audit),
        "Expected_Ion_Hypotheses": int(
            summary.get("expected_ion_hypotheses")
            or catalog.get("expected_ion_hypotheses")
            or 0
        ),
        "Detected_Features": int(
            summary.get("detected_feature_count")
            or len(payload.get("features") or [])
        ),
        "Compound_Summary_Rows": len(compound_rows),
        "Supported_Compounds": sum(
            int(row.get("_Supported_Feature_Count") or 0) > 0
            for row in compound_rows
        ),
        "Mass_Only_Compounds": int(
            sum(
                int(row.get("_Supported_Feature_Count") or 0) == 0
                for row in compound_rows
            )
        ),
        "Compound_Evidence_Counts": compound_evidence_counts,
        "Supported_Features": int(summary.get("supported_features") or 0),
        "Mass_Only_Features": int(summary.get("mass_only_features") or 0),
        "Ambiguous_Features": int(
            summary.get("ambiguous_feature_count") or 0
        ),
        "Candidate_Target_Rows": int(
            summary.get("candidate_target_rows") or 0
        ),
        "Candidate_Target_Names": int(
            summary.get("candidate_target_names") or 0
        ),
        "MS1_Exact_Mass_Candidates": int(
            evidence.get("MS1_EXACT_MASS_CANDIDATE") or 0
        ),
        "MS1_Isotope_Supported": int(
            evidence.get("MS1_ISOTOPE_SUPPORTED") or 0
        ),
        "MS1_MS2_Diagnostic_Probable": int(
            evidence.get("MS1_MS2_DIAGNOSTIC_PROBABLE") or 0
        ),
        "Experimental_Library_Confirmed": int(
            evidence.get("EXPERIMENTAL_LIBRARY_CONFIRMED") or 0
        ),
        "MS1_Scan_Count": int(summary.get("ms1_scan_count") or 0),
        "Usable_MS1_Scan_Count": int(
            summary.get("usable_ms1_scan_count") or 0
        ),
        "Selected_Datasets": _join(selected),
        "Warning": str(summary.get("warning") or payload.get("warning") or ""),
        "Ordinary_MS2_Match_Increment": 0,
    }


def _expected_raw_detail_rows(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> list[dict[str, Any]]:
    language = "en" if language == "en" else "zh"
    labels = EVIDENCE_LABELS[language]
    rows: list[dict[str, Any]] = []
    for feature in payload.get("features") or []:
        if not isinstance(feature, Mapping):
            continue
        all_members = _deduplicate_members([
            member
            for member in feature.get("target_members") or []
            if isinstance(member, Mapping)
        ])
        members = _selected_target_members(feature)
        names = _clean_member_candidate_names(
            members,
            [
                *(feature.get("candidate_names") or []),
                feature.get("candidate_name"),
            ],
        )
        diagnostic = [
            item for item in feature.get("ms2_diagnostic_support") or []
            if isinstance(item, Mapping)
        ]
        library = [
            item for item in feature.get("experimental_library_support") or []
            if isinstance(item, Mapping)
        ]
        associated = [
            item for item in feature.get("associated_ms2") or []
            if isinstance(item, Mapping)
        ]
        matched_diagnostic = []
        for item in diagnostic:
            spectrum = str(item.get("spectrum_id") or "")
            ions = ",".join(
                f"{float(value):.6f}"
                for value in item.get("matched_diagnostic_ions") or []
                if _number(value) is not None
            )
            matched_diagnostic.append(
                f"{spectrum}:{ions}" if spectrum else ions
            )
        isotope = feature.get("formula_isotope_evidence") or {}
        level = str(
            feature.get("evidence_level") or "MS1_EXACT_MASS_CANDIDATE"
        )
        candidate_name = str(names[0] if names else "")
        target_databases = (
            feature.get("dataset_sources")
            or [member.get("dataset_source") for member in members]
        )
        source_plants = _unique([
            member.get("plant_source") for member in members
        ])
        target_inchikeys = _unique([
            member.get("inchikey") for member in members
        ])
        identity_keys = _unique([
            _identity_key(member) for member in members
        ])
        group_identity_keys = _unique([
            _identity_key(member).casefold() for member in all_members
        ])
        raw_group_candidate_names = _unique(
            list(feature.get("candidate_names") or [])
        )
        group_candidate_names = _clean_member_candidate_names(
            all_members, raw_group_candidate_names
        )
        if not group_candidate_names:
            group_candidate_names = list(names)
        candidate_count = (
            len(identity_keys)
            if diagnostic or library
            else int(
                feature.get("identity_count")
                or len(identity_keys)
                or len(names)
                or 0
            )
        )
        primary_candidate = (
            names[0] if names else candidate_name
        )
        candidate_display = (
            primary_candidate
            if candidate_count <= 1
            else (
                f"{primary_candidate} 等 {candidate_count} 个候选"
                if language == "zh"
                else f"{primary_candidate} + {candidate_count - 1} candidates"
            )
        )
        rows.append({
            "Evidence": labels.get(level, level),
            "Candidate_Name": candidate_name,
            "Candidate_Display": candidate_display,
            "Candidate_Count": candidate_count,
            "Target_Database": _display_sources(target_databases),
            "Source_Plant": _join(source_plants),
            "Formula": str(feature.get("formula") or ""),
            "Adduct": str(feature.get("adduct") or ""),
            "Polarity": str(feature.get("polarity") or ""),
            "Target_mz": _number(feature.get("target_mz")),
            "Observed_mz": _number(feature.get("observed_apex_mz")),
            "Mass_Error_ppm": _number(feature.get("mass_error_ppm")),
            "RT_min": _number(feature.get("apex_rt_minutes")),
            "Peak_Area": _number(
                feature.get("baseline_corrected_area_intensity_seconds")
                or feature.get("xic_area_intensity_seconds")
            ),
            "Peak_Height": _number(feature.get("apex_intensity")),
            "Signal_to_Noise": _number(feature.get("signal_to_noise")),
            "Scan_Points": int(feature.get("scan_points") or 0),
            "MS2_Diagnostic_Support": _unique_support_count(
                diagnostic, library=False
            ),
            "Experimental_Library_Support": _unique_support_count(
                library, library=True
            ),
            "Peak_Quality": _number(feature.get("peak_quality_score")),
            "Gap_Ratio": _number(feature.get("gap_ratio")),
            "OpenMS_Verification": _openms_label(feature, language),
            "Background_Status": _background_label(feature, language),
            "InChIKey": _join(target_inchikeys),
            "Result_Note": _result_note(
                level, candidate_count, language
            ),
            "Feature_ID": str(feature.get("feature_id") or ""),
            "All_Candidate_Names": _join(names),
            "Observed_Apex_mz": _number(feature.get("observed_apex_mz")),
            "Apex_RT_min": _number(feature.get("apex_rt_minutes")),
            "Boundary_Start_RT_min": _number(
                feature.get("boundary_start_rt_minutes")
            ),
            "Boundary_End_RT_min": _number(
                feature.get("boundary_end_rt_minutes")
            ),
            "FWHM_s": _number(feature.get("fwhm_seconds")),
            "Apex_Intensity": _number(feature.get("apex_intensity")),
            "XIC_Area_Intensity_s": _number(
                feature.get("baseline_corrected_area_intensity_seconds")
                or feature.get("xic_area_intensity_seconds")
            ),
            "Signal_to_Noise": _number(feature.get("signal_to_noise")),
            "Scan_Points": int(feature.get("scan_points") or 0),
            "XIC_Nonzero_Points": int(
                feature.get("xic_nonzero_points") or 0
            ),
            "Evidence_Level": level,
            "Evidence_Label": labels.get(level, level),
            "Orthogonal_Evidence_Supported": bool(
                feature.get("cd_pattern_accepted")
            ),
            "Isomer_Ambiguous": bool(feature.get("isomer_ambiguous")),
            "Identity_Count": candidate_count,
            "Identification_Scope": str(
                feature.get("identification_scope") or ""
            ),
            "Isotope_Evidence_Status": str(isotope.get("status") or ""),
            "Associated_MS2_Count": len(associated),
            "Diagnostic_MS2_Support_Count": len(diagnostic),
            "Matched_Diagnostic_Ions": _join(matched_diagnostic),
            "Experimental_Library_Support_Count": len(library),
            "Experimental_Library_Spectrum_IDs": _join(_unique([
                item.get("library_spectrum_id") for item in library
            ])),
            "Target_IDs": _join(_unique([
                member.get("target_id") for member in members
            ])),
            "Target_InChIKeys": _join(_unique([
                member.get("inchikey") for member in members
            ])),
            "Dataset_Sources": _join(
                feature.get("dataset_sources") or []
            ),
            "Generation_Sources": _join(
                feature.get("generation_sources") or []
            ),
            "Evidence_Interpretation": str(
                feature.get("evidence_interpretation") or ""
            ),
            "Counted_In_Ordinary_MS2_Matches": False,
            "_Identity_Keys": identity_keys,
            "_Candidate_Names": names,
            "_Group_Identity_Keys": group_identity_keys,
            "_Group_Candidate_Names": group_candidate_names,
            "_Raw_Group_Candidate_Names": raw_group_candidate_names,
            "_Target_Databases": _unique(
                [target_databases]
                if isinstance(target_databases, str)
                else list(target_databases or [])
            ),
            "_Diagnostic_Support": diagnostic,
            "_Library_Support": library,
            "_OpenMS_Status": str(
                feature.get("openms_feature_link_status") or ""
            ),
            "_Blank_Status": str(feature.get("blank_status") or ""),
            "_Blank_Category": str(
                feature.get("blank_evidence_category") or ""
            ),
            "_Blank_Flagged": bool(feature.get("blank_flagged")),
            "_CD_Pattern_Accepted": bool(
                feature.get("cd_pattern_accepted")
            ),
            "_OpenMS_Supported": bool(
                feature.get("openms_feature_supported")
            ),
            "_Expected_Ion_ID": str(
                feature.get("expected_ion_id") or ""
            ),
        })
    priority = {
        "EXPERIMENTAL_LIBRARY_CONFIRMED": 4,
        "MS1_MS2_DIAGNOSTIC_PROBABLE": 3,
        "MS1_ISOTOPE_SUPPORTED": 2,
        "MS1_EXACT_MASS_CANDIDATE": 1,
    }
    rows.sort(
        key=lambda row: (
            priority.get(str(row.get("Evidence_Level") or ""), 0),
            float(row.get("Peak_Height") or 0),
        ),
        reverse=True,
    )
    return rows


def _expected_detail_rows(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> list[dict[str, Any]]:
    """Return one readable row per compound/ambiguity group.

    Grouping intentionally matches the browser: formula plus the complete
    target-identity set.  A diagnostic result may resolve one member for its
    evidence label, but it must not split the parent ambiguity group into an
    extra displayed compound.  Different adducts and RT peaks remain in the
    raw audit layer.
    """

    language = "en" if language == "en" else "zh"
    raw_rows = _expected_raw_detail_rows(payload, language)
    evidence_priority = {
        "EXPERIMENTAL_LIBRARY_CONFIRMED": 4,
        "MS1_MS2_DIAGNOSTIC_PROBABLE": 3,
        "MS1_ISOTOPE_SUPPORTED": 2,
        "MS1_EXACT_MASS_CANDIDATE": 1,
    }

    def representative_rank(row: Mapping[str, Any]) -> tuple[Any, ...]:
        blank_rank = {
            "sample_enriched": 3,
            "not_detected": 2,
            "blank_associated": 1,
            "blank_dominant": 0,
        }
        return (
            evidence_priority.get(
                str(row.get("Evidence_Level") or ""), 0
            ),
            int(bool(row.get("_CD_Pattern_Accepted"))),
            len(row.get("_Library_Support") or []),
            len(row.get("_Diagnostic_Support") or []),
            int(str(row.get("Isotope_Evidence_Status") or "") == "supporting"),
            int(bool(row.get("_OpenMS_Supported"))),
            blank_rank.get(str(row.get("_Blank_Category") or ""), 0),
            float(row.get("Peak_Quality") or 0),
            float(row.get("Peak_Area") or 0),
            float(row.get("Peak_Height") or 0),
            -abs(float(row.get("Mass_Error_ppm") or 0)),
            -float(row.get("Apex_RT_min") or 0),
        )

    grouped: dict[tuple[str, tuple[str, ...]], list[dict[str, Any]]] = {}
    for row in raw_rows:
        identities = tuple(sorted(
            str(value).strip().casefold()
            for value in row.get("_Group_Identity_Keys") or []
            if str(value).strip()
        ))
        fallback_names = tuple(sorted(
            str(name).strip().casefold()
            for name in (
                row.get("_Group_Candidate_Names")
                or row.get("_Candidate_Names")
                or []
            )
            if str(name).strip()
        ))
        key = (
            str(row.get("Formula") or "").strip(),
            identities or fallback_names,
        )
        grouped.setdefault(key, []).append(row)

    output: list[dict[str, Any]] = []
    for (_formula, _group_tokens), rows in grouped.items():
        ranked_rows = sorted(
            rows, key=representative_rank, reverse=True
        )
        representative = ranked_rows[0]
        candidate_names = _unique([
            name
            for row in ranked_rows
            for name in (
                row.get("_Group_Candidate_Names")
                or row.get("_Candidate_Names")
                or []
            )
        ])
        if not candidate_names:
            candidate_names = _unique([
                name
                for row in ranked_rows
                for name in row.get("_Candidate_Names") or []
            ])
        identity_keys = _unique([
            identity
            for row in ranked_rows
            for identity in row.get("_Group_Identity_Keys") or []
        ])
        candidate_count = (
            len(identity_keys) or len(candidate_names) or 1
        )
        candidate_display = _candidate_display(
            candidate_names, candidate_count, language
        )
        databases = _unique([
            source
            for row in ranked_rows
            for source in row.get("_Target_Databases") or []
        ])
        diagnostic_support = [
            item
            for row in ranked_rows
            for item in row.get("_Diagnostic_Support") or []
        ]
        library_support = [
            item
            for row in ranked_rows
            for item in row.get("_Library_Support") or []
        ]
        adducts = sorted(_unique([
            row.get("Adduct") for row in ranked_rows
        ]))
        rt_values = [
            float(value)
            for row in ranked_rows
            for value in [_number(row.get("Apex_RT_min"))]
            if value is not None
        ]
        supported_feature_count = sum(
            bool(row.get("_CD_Pattern_Accepted"))
            for row in ranked_rows
        )
        ion_hypothesis_count = len(_unique([
            row.get("_Expected_Ion_ID") for row in ranked_rows
        ]))
        combined = dict(representative)
        combined.update({
            "Candidate_Name": candidate_display,
            "Candidate_Display": candidate_display,
            "Candidate_Count": candidate_count,
            "Target_Database": _display_sources(databases),
            "Related_Adducts": _join(adducts),
            "Raw_Peak_Count": len(ranked_rows),
            "Component_Peak_Count": len(ranked_rows),
            "MS2_Diagnostic_Support": _unique_support_count(
                diagnostic_support, library=False
            ),
            "Experimental_Library_Support": _unique_support_count(
                library_support, library=True
            ),
            "All_Candidate_Names": _join(candidate_names),
            "Dataset_Sources": _join(databases),
            "Identity_Count": candidate_count,
            "Isomer_Ambiguous": candidate_count > 1,
            "Multiple_Names_Same_Identity": (
                candidate_count == 1 and len(candidate_names) > 1
            ),
            "Result_Note": _result_note(
                str(representative.get("Evidence_Level") or ""),
                candidate_count,
                language,
            ),
            "_Identity_Keys": identity_keys,
            "_Candidate_Names": candidate_names,
            "_Raw_Rows": ranked_rows,
            "_Supported_Feature_Count": supported_feature_count,
            "_Ion_Hypothesis_Count": ion_hypothesis_count,
            "_RT_Min": min(rt_values) if rt_values else None,
            "_RT_Max": max(rt_values) if rt_values else None,
            "_Representative_Rank": representative_rank(representative),
        })
        output.append(combined)

    output.sort(
        key=lambda row: (
            row.get("_Representative_Rank") or (),
            int(row.get("_Supported_Feature_Count") or 0),
            int(row.get("Raw_Peak_Count") or 0),
        ),
        reverse=True,
    )
    return output


def expected_feature_rows(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> list[dict[str, Any]]:
    """Return only the focused columns used by standard Excel/CSV exports."""
    rows = []
    for row in _expected_detail_rows(payload, language):
        focused = {
            column: row.get(column, "")
            for column in EXPECTED_FEATURE_COLUMNS
        }
        focused["Candidate_Name"] = row.get(
            "Candidate_Display", row.get("Candidate_Name", "")
        )
        rows.append(focused)
    return rows


def expected_group_rows(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> list[dict[str, Any]]:
    """Return a compact JSON projection of the deduplicated candidate groups.

    The full Expected Compounds artifact intentionally retains every XIC peak
    and every target-member record for audit/export.  Sending that artifact to
    the browser made routine result viewing unnecessarily expensive.  This
    projection keeps the scientific evidence and group-level ambiguity while
    omitting ``_Raw_Rows`` and other internal trace payloads.  It therefore
    changes presentation only; no feature or evidence level is recomputed.
    """

    rows: list[dict[str, Any]] = []
    for detail in _expected_detail_rows(payload, language):
        identity_tokens = sorted(
            str(value).strip().casefold()
            for value in detail.get("_Identity_Keys") or []
            if str(value).strip()
        )
        candidate_names = [
            str(value).strip()
            for value in detail.get("_Candidate_Names") or []
            if str(value).strip()
        ]
        stable_basis = "|".join([
            str(detail.get("Formula") or "").strip(),
            *(identity_tokens or [name.casefold() for name in candidate_names]),
        ])
        group_id = "ecg_" + hashlib.sha256(
            stable_basis.encode("utf-8")
        ).hexdigest()[:16]
        rows.append({
            "group_id": group_id,
            "evidence_level": str(
                detail.get("Evidence_Level")
                or "MS1_EXACT_MASS_CANDIDATE"
            ),
            "evidence_label": str(detail.get("Evidence") or ""),
            "candidate_display": str(
                detail.get("Candidate_Display")
                or detail.get("Candidate_Name")
                or ""
            ),
            "candidate_names": candidate_names,
            "candidate_count": int(detail.get("Candidate_Count") or 0),
            "target_database": str(detail.get("Target_Database") or ""),
            "formula": str(detail.get("Formula") or ""),
            "adduct": str(detail.get("Adduct") or ""),
            "related_adducts": str(detail.get("Related_Adducts") or ""),
            "polarity": str(detail.get("Polarity") or ""),
            "target_mz": _number(detail.get("Target_mz")),
            "observed_mz": _number(detail.get("Observed_mz")),
            "mass_error_ppm": _number(detail.get("Mass_Error_ppm")),
            "rt_min": _number(detail.get("RT_min")),
            "rt_range_min": _number(detail.get("_RT_Min")),
            "rt_range_max": _number(detail.get("_RT_Max")),
            "peak_area": _number(detail.get("Peak_Area")),
            "peak_height": _number(detail.get("Peak_Height")),
            "signal_to_noise": _number(detail.get("Signal_to_Noise")),
            "peak_quality": _number(detail.get("Peak_Quality")),
            "openms_verification": str(
                detail.get("OpenMS_Verification") or ""
            ),
            "background_status": str(
                detail.get("Background_Status") or ""
            ),
            "raw_peak_count": int(detail.get("Raw_Peak_Count") or 0),
            "supported_feature_count": int(
                detail.get("_Supported_Feature_Count") or 0
            ),
            "ion_hypothesis_count": int(
                detail.get("_Ion_Hypothesis_Count") or 0
            ),
            "ms2_diagnostic_support": int(
                detail.get("MS2_Diagnostic_Support") or 0
            ),
            "experimental_library_support": int(
                detail.get("Experimental_Library_Support") or 0
            ),
            "isomer_ambiguous": bool(detail.get("Isomer_Ambiguous")),
            "ambiguity_scope": (
                "same_formula_or_mass_database_candidates"
                if bool(detail.get("Isomer_Ambiguous"))
                else "single_database_identity"
            ),
            "multiple_names_same_identity": bool(
                detail.get("Multiple_Names_Same_Identity")
            ),
            "result_note": str(detail.get("Result_Note") or ""),
        })
    return rows


def expected_raw_feature_rows(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> list[dict[str, Any]]:
    """Return the one-row-per-XIC-segment audit table."""

    rows = []
    for detail in _expected_raw_detail_rows(payload, language):
        rows.append({
            column: detail.get(column, "")
            for column in EXPECTED_RAW_FEATURE_COLUMNS
        })
    return rows


def expected_feature_column_labels(language: str = "zh") -> list[str]:
    language = "en" if language == "en" else "zh"
    labels = EXPECTED_FEATURE_LABELS[language]
    return [labels[column] for column in EXPECTED_FEATURE_COLUMNS]


def expected_raw_feature_column_labels(
    language: str = "zh",
) -> list[str]:
    language = "en" if language == "en" else "zh"
    labels = EXPECTED_RAW_FEATURE_LABELS[language]
    return [
        labels[column] for column in EXPECTED_RAW_FEATURE_COLUMNS
    ]


def expected_feature_groups(
    payload: Mapping[str, Any],
    language: str = "zh",
) -> dict[str, list[dict[str, Any]]]:
    """Split targeted results by evidence so weak candidates cannot bury hits."""

    groups = {
        "high": [],
        "isotope": [],
        "mass_only": [],
    }
    for detail in _expected_detail_rows(payload, language):
        focused = {
            column: detail.get(column, "")
            for column in EXPECTED_FEATURE_COLUMNS
        }
        focused["Candidate_Name"] = detail.get(
            "Candidate_Display", detail.get("Candidate_Name", "")
        )
        level = str(detail.get("Evidence_Level") or "")
        supported_count = int(
            detail.get("_Supported_Feature_Count") or 0
        )
        if supported_count and level in {
            "EXPERIMENTAL_LIBRARY_CONFIRMED",
            "MS1_MS2_DIAGNOSTIC_PROBABLE",
        }:
            groups["high"].append(focused)
        elif supported_count:
            groups["isotope"].append(focused)
        else:
            groups["mass_only"].append(focused)
    return groups


def append_expected_csv_section(
    writer: Any,
    payload: Mapping[str, Any],
    language: str = "zh",
) -> None:
    """Append a visibly separate Expected Compounds table to an existing CSV."""

    language = "en" if language == "en" else "zh"
    summary = expected_summary_row(payload, language)
    rows = expected_feature_rows(payload, language)
    title = (
        "=== Expected Compounds（独立连续 MS1 靶向层；不计入普通 MS2 匹配）==="
        if language == "zh"
        else "=== Expected Compounds (independent continuous-MS1 targeted layer; "
        "excluded from ordinary MS2 match counts) ==="
    )
    writer.writerow([])
    writer.writerow([title])
    writer.writerow(["Item", "Value", "Note"])
    summary_items = [
        ("Status", summary["Status_Label"], ""),
        ("Target_Database", _display_sources(
            str(summary.get("Selected_Datasets") or "").split(" | ")
        ), ""),
        ("Target_Catalog_Rows", summary["Target_Catalog_Rows"], ""),
        (
            "Active_Target_Rows",
            summary["Active_Target_Rows"],
            "通过运行时身份与化学元数据门控" if language == "zh"
            else "Passed runtime identity and chemistry metadata gates",
        ),
        (
            "Quarantined_Target_Rows",
            summary["Quarantined_Target_Rows"],
            "保留源库记录，但不生成生产离子假设" if language == "zh"
            else "Source rows retained but excluded from production ion hypotheses",
        ),
        ("Expected_Ion_Hypotheses", summary["Expected_Ion_Hypotheses"], ""),
        (
            "Compound_Summary_Rows",
            summary["Compound_Summary_Rows"],
            (
                "标准靶向表：每个去重靶向化合物仅一行；候选总数另列"
                if language == "zh"
                else "Routine table: one row per deduplicated targeted compound; candidate count is separate"
            ),
        ),
        (
            "Raw_XIC_Peak_Segments",
            summary["Detected_Features"],
            (
                "完整逐峰记录请下载靶向峰审计 CSV"
                if language == "zh"
                else "Download targeted peak-audit CSV for every segment"
            ),
        ),
        ("Supported_Compounds", summary["Supported_Compounds"], ""),
        ("Mass_Only_Compounds", summary["Mass_Only_Compounds"], ""),
        (
            "Counting_Note",
            "Independent",
            (
                "不计入普通 MS2 谱库匹配数"
                if language == "zh"
                else "Excluded from ordinary MS2 library-match counts"
            ),
        ),
    ]
    writer.writerows(summary_items)
    writer.writerow([])
    writer.writerow(expected_feature_column_labels(language))
    for row in rows:
        writer.writerow([row.get(column, "") for column in EXPECTED_FEATURE_COLUMNS])


def style_analysis_workbook(workbook: Any) -> None:
    """Keep each record on one line and make the intended reading order clear."""

    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.styles.numbers import BUILTIN_FORMATS_REVERSE
    from openpyxl.utils import get_column_letter
    import copy

    header_fill = PatternFill("solid", fgColor="17365D")
    header_font = Font(color="FFFFFF", bold=True, size=10)
    header_border = Border(
        bottom=Side(style="medium", color="9EADBC")
    )
    band_fill = PatternFill("solid", fgColor="F4F7FB")
    text_alignment = Alignment(
        horizontal="left", vertical="center", wrap_text=False
    )
    number_alignment = Alignment(
        horizontal="right", vertical="center", wrap_text=False
    )
    summary_alignment = Alignment(
        horizontal="left", vertical="center", wrap_text=True
    )
    width_overrides = {
        "证据等级": 23,
        "主要候选": 38,
        "候选数": 9,
        "靶向数据库": 24,
        "分子式": 14,
        "代表加合物": 14,
        "共洗脱加合物": 25,
        "离子模式": 10,
        "理论 m/z": 14,
        "实测 m/z": 14,
        "误差 (ppm)": 12,
        "主峰 RT (min)": 13,
        "代表峰面积": 15,
        "峰质量 (0–10)": 13,
        "OpenMS 峰验证": 24,
        "背景状态": 22,
        "原始峰段数": 12,
        "信噪比": 11,
        "诊断 MS2 谱图数": 16,
        "实验谱库谱图数": 16,
        "靶向峰 ID": 24,
        "候选名称": 38,
        "顶点 RT (min)": 13,
        "起始 RT (min)": 13,
        "结束 RT (min)": 13,
        "FWHM (s)": 11,
        "峰高": 15,
        "峰面积": 15,
        "有效扫描点": 12,
        "缺口比例": 11,
        "适用模块": 25,
        "代码": 34,
        "状态说明": 52,
        "建议操作": 58,
        "化合物": 38,
        "鉴定组": 10,
        "候选化合物": 38,
        "组内排名": 10,
        "候选角色": 29,
        "组内候选数": 11,
        "鉴定结论": 31,
        "鉴定组置信度": 14,
        "候选置信度": 14,
        "MS2 隔离目标 m/z": 17,
        "检索用实测/回退 m/z": 19,
        "检索质量轴来源": 29,
        "样本 CE": 11,
        "样本 CE 单位": 12,
        "样本仪器类型": 17,
        "样品离子形式": 16,
        "样品离子形式依据": 36,
        "候选谱库加合物": 18,
        "代表馆藏 CE": 13,
        "代表馆藏 CE 单位": 15,
        "代表馆藏仪器类型": 19,
        "候选 MS/MS 相似度": 19,
        "候选谱熵相似度": 18,
        "候选综合匹配分": 18,
        "候选质量误差 (ppm)": 20,
        "代表谱匹配碎片数": 19,
        "代表谱匹配碎片 m/z": 32,
        "代表谱匹配碎片 (m/z:强度)": 40,
        "样本 MS2 碎片离子（实际检索谱，m/z:强度）": 42,
        "实际检索 MS2 峰数": 18,
        "原始 MS2 碎片离子（已保存输入峰，m/z:强度）": 42,
        "原始 MS2 输入峰数": 18,
        "候选化学类别": 22,
        "鉴定组化学类别": 22,
        "分类置信度": 14,
        "分类依据": 24,
        "分类证据状态": 15,
        "分类可靠性说明": 36,
        "馆藏谱图编号": 24,
        "InChIKey": 29,
        "SMILES": 48,
        "样本谱图编号": 21,
        "备注": 40,
        "前体 m/z": 14,
        "置信度": 14,
        "Final Score": 12,
        "Cosine": 11,
        "化学类别": 19,
        "谱库来源": 24,
        "谱图编号": 17,
        "Spectrum": 17,
        "Compound_Name": 38,
        "Formula": 14,
        "Adduct": 14,
        "Precursor_mz": 14,
        "Retention_Time": 11,
        "Confidence": 14,
        "Final_Score": 12,
        "Cosine_Score": 11,
        "Mass_PPM": 12,
        "Classification": 19,
        "Dataset_Source": 24,
        "Peaks_Count": 11,
        "Status": 24,
        "Identification group": 18,
        "Candidate compound": 38,
        "Rank within group": 16,
        "Candidate role": 32,
        "Candidates in group": 18,
        "Identification conclusion": 34,
        "Group confidence": 18,
        "Candidate confidence": 20,
    }
    numeric_formats = {
        "理论 m/z": "0.000000",
        "实测 m/z": "0.000000",
        "前体 m/z": "0.0000",
        "MS2 隔离目标 m/z": "0.0000",
        "检索用实测/回退 m/z": "0.0000",
        "误差 (ppm)": "0.00",
        "主峰 RT (min)": "0.000",
        "代表峰面积": "0.00E+00",
        "峰质量 (0–10)": "0.00",
        "原始峰段数": "0",
        "信噪比": "0.00",
        "诊断 MS2 谱图数": "0",
        "实验谱库谱图数": "0",
        "顶点 RT (min)": "0.000",
        "起始 RT (min)": "0.000",
        "结束 RT (min)": "0.000",
        "FWHM (s)": "0.00",
        "峰高": "0.00E+00",
        "峰面积": "0.00E+00",
        "有效扫描点": "0",
        "缺口比例": "0.000",
        "候选数": "0",
        "Final Score": "0.000",
        "Cosine": "0.000",
        "候选 MS/MS 相似度": "0.0000",
        "候选谱熵相似度": "0.0000",
        "候选综合匹配分": "0.0000",
        "候选质量误差 (ppm)": "0.00",
        "Target_mz": "0.000000",
        "Observed_mz": "0.000000",
        "Precursor_mz": "0.0000",
        "Mass_Error_ppm": "0.00",
        "Mass_PPM": "0.00",
        "RT_min": "0.000",
        "Retention_Time": "0.000",
        "Peak_Area": "0.00E+00",
        "Signal_to_Noise": "0.00",
        "Cosine_Score": "0.000",
        "Final_Score": "0.000",
        "MS2_Diagnostic_Support": "0",
        "Experimental_Library_Support": "0",
        "Candidate_Count": "0",
        "Peaks_Count": "0",
    }
    evidence_fills = {
        "MS/MS 实验谱库候选支持": PatternFill("solid", fgColor="C6E0B4"),
        "MS1 + 诊断离子筛查支持": PatternFill(
            "solid", fgColor="BDD7EE"
        ),
        "MS/MS experimental-library candidate support": PatternFill(
            "solid", fgColor="C6E0B4"
        ),
        "MS1 + diagnostic-ion screening support": PatternFill(
            "solid", fgColor="BDD7EE"
        ),
        # Retain fills for archived V12/V14 rows carrying legacy labels.
        "严格实验谱库支持": PatternFill("solid", fgColor="C6E0B4"),
        "MS1 + 实测 MS2 诊断离子": PatternFill(
            "solid", fgColor="BDD7EE"
        ),
        "MS1 + 同位素支持": PatternFill("solid", fgColor="FFF2CC"),
        "MS1 精确质量候选": PatternFill("solid", fgColor="E7E6E6"),
        "Experimental-library confirmed": PatternFill(
            "solid", fgColor="C6E0B4"
        ),
        "MS1 + measured MS2 diagnostic ions": PatternFill(
            "solid", fgColor="BDD7EE"
        ),
        "MS1 + isotope supported": PatternFill(
            "solid", fgColor="FFF2CC"
        ),
        "MS1 exact-mass candidate": PatternFill(
            "solid", fgColor="E7E6E6"
        ),
    }

    for worksheet in workbook.worksheets:
        worksheet.sheet_view.showGridLines = False
        worksheet.sheet_view.zoomScale = 90
        worksheet.freeze_panes = (
            "A2" if worksheet.title.startswith("01_") else "B2"
        )
        if worksheet.title.startswith("00_"):
            tab_color = "147D64"
        elif worksheet.title.startswith("01_"):
            tab_color = "1F4E78"
        elif "高证据" in worksheet.title:
            tab_color = "70AD47"
        elif "同位素" in worksheet.title:
            tab_color = "FFC000"
        elif "仅质量" in worksheet.title:
            tab_color = "A5A5A5"
        elif "模块提示" in worksheet.title or "Warnings" in worksheet.title:
            tab_color = "FFC000"
        elif "峰审计" in worksheet.title or "Peak_Audit" in worksheet.title:
            tab_color = "C55A11"
        elif "MS2" in worksheet.title:
            tab_color = "5B9BD5"
        else:
            tab_color = "C55A11"
        worksheet.sheet_properties.tabColor = tab_color
        worksheet.sheet_properties.pageSetUpPr.fitToPage = True
        worksheet.page_setup.fitToWidth = 1
        worksheet.page_setup.fitToHeight = 0
        worksheet.print_title_rows = "1:1"
        worksheet.row_dimensions[1].height = 30

        headers = {
            cell.column: str(cell.value or "")
            for cell in worksheet[1]
        }
        for cell in worksheet[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.border = header_border
            cell.alignment = Alignment(
                horizontal="center",
                vertical="center",
                wrap_text=True,
            )

        if (
            "模块提示" in worksheet.title
            or "Module_Warnings" in worksheet.title
        ):
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for column_index, header in headers.items():
                worksheet.column_dimensions[
                    get_column_letter(column_index)
                ].width = width_overrides.get(header, {
                    "Severity": 12,
                    "Modules": 25,
                    "Code": 34,
                    "Message": 52,
                    "Action": 58,
                }.get(header, 22))
            for row_index in range(2, worksheet.max_row + 1):
                worksheet.row_dimensions[row_index].height = 42
                for cell in worksheet[row_index]:
                    cell.alignment = Alignment(
                        horizontal="left",
                        vertical="top",
                        wrap_text=True,
                    )
            continue

        if worksheet.title.startswith("00_"):
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            guide_widths = {
                "说明类型": 14,
                "名称": 32,
                "它代表什么": 58,
                "可以怎样使用": 58,
                "不能怎样理解": 58,
                "建议下一步": 58,
                "Section": 16,
                "Name": 32,
                "Meaning": 58,
                "Recommended use": 58,
                "Do not infer": 58,
                "Next step": 58,
            }
            for column_index, header in headers.items():
                worksheet.column_dimensions[
                    get_column_letter(column_index)
                ].width = guide_widths.get(header, 28)
            for row_index in range(2, worksheet.max_row + 1):
                worksheet.row_dimensions[row_index].height = (
                    120 if row_index == 2 else 84
                )
                for cell in worksheet[row_index]:
                    cell.alignment = Alignment(
                        horizontal="left",
                        vertical="top",
                        wrap_text=True,
                    )
                if row_index % 2 == 0:
                    for cell in worksheet[row_index]:
                        cell.fill = band_fill
            continue

        if worksheet.title.startswith("01_"):
            worksheet.column_dimensions["A"].width = 29
            worksheet.column_dimensions["B"].width = 32
            worksheet.column_dimensions["C"].width = 68
            for row_index in range(2, worksheet.max_row + 1):
                worksheet.row_dimensions[row_index].height = 27
                for cell in worksheet[row_index]:
                    cell.alignment = summary_alignment
                if row_index in {2, 3, 4, 5}:
                    worksheet.cell(row_index, 1).font = Font(
                        bold=True, color="1F4E78"
                    )
            continue

        worksheet.auto_filter.ref = worksheet.dimensions
        sample_limit = min(worksheet.max_row, 200)
        evidence_column = next(
            (
                index for index, header in headers.items()
                if header in {"证据等级", "Evidence"}
            ),
            None,
        )
        observed_widths = {
            index: len(header) for index, header in headers.items()
        }
        style_cache = {}

        # Register style components in the same order as the former
        # row-pass/column-pass/evidence-pass loops. This preserves even the
        # internal StyleArray identifiers while allowing repeated final styles
        # to be assigned in one operation.
        if worksheet.max_row >= 2:
            workbook._fills.add(band_fill)
            registered_alignments = set()
            for column_index in range(1, worksheet.max_column + 1):
                header = headers.get(column_index, "")
                number_format = numeric_formats.get(header)
                saw_number = False
                for row_index in range(2, worksheet.max_row + 1):
                    is_number = isinstance(
                        worksheet.cell(row_index, column_index).value,
                        (int, float),
                    )
                    alignment_name = "number" if is_number else "text"
                    if alignment_name not in registered_alignments:
                        workbook._alignments.add(
                            number_alignment
                            if is_number
                            else text_alignment
                        )
                        registered_alignments.add(alignment_name)
                    saw_number = saw_number or is_number
                    if (
                        len(registered_alignments) == 2
                        and (not number_format or saw_number)
                    ):
                        break
                if (
                    number_format
                    and saw_number
                    and number_format not in BUILTIN_FORMATS_REVERSE
                ):
                    workbook._number_formats.add(number_format)
            if evidence_column:
                registered_fills = set()
                for row_index in range(2, worksheet.max_row + 1):
                    value = str(
                        worksheet.cell(row_index, evidence_column).value or ""
                    )
                    if value in evidence_fills and value not in registered_fills:
                        workbook._fills.add(evidence_fills[value])
                        registered_fills.add(value)

        for row_index, row in enumerate(
            worksheet.iter_rows(
                min_row=2,
                max_row=worksheet.max_row,
                max_col=worksheet.max_column,
            ),
            start=2,
        ):
            worksheet.row_dimensions[row_index].height = 21
            for column_index, cell in enumerate(row, start=1):
                header = headers.get(column_index, "")
                if row_index <= sample_limit and cell.value is not None:
                    observed_widths[column_index] = max(
                        observed_widths[column_index], len(str(cell.value))
                    )
                is_number = isinstance(cell.value, (int, float))
                alignment_name = "number" if is_number else "text"
                number_format = (
                    numeric_formats.get(header) if is_number else None
                )
                evidence_value = (
                    str(cell.value or "")
                    if column_index == evidence_column
                    else ""
                )
                fill_name = (
                    evidence_value
                    if evidence_value in evidence_fills
                    else ("band" if row_index % 2 == 0 else "")
                )
                base_style = tuple(cell._style) if cell.has_style else None
                cache_key = (
                    base_style,
                    alignment_name,
                    number_format,
                    fill_name,
                )
                cached_style = style_cache.get(cache_key)
                if cached_style is not None:
                    cell._style = copy.copy(cached_style)
                    continue
                if fill_name == "band":
                    cell.fill = band_fill
                elif fill_name:
                    cell.fill = evidence_fills[fill_name]
                cell.alignment = (
                    number_alignment if is_number else text_alignment
                )
                if number_format:
                    cell.number_format = number_format
                style_cache[cache_key] = copy.copy(cell._style)

        for column_index in range(1, worksheet.max_column + 1):
            header = headers.get(column_index, "")
            width = width_overrides.get(
                header,
                min(max(observed_widths[column_index] + 2, 10), 22),
            )
            worksheet.column_dimensions[
                get_column_letter(column_index)
            ].width = width


def _fmt(value: Any, digits: int = 3, missing: str = "—") -> str:
    number = _number(value)
    if number is None:
        return missing
    return f"{float(number):.{digits}f}"


def render_expected_html_section(
    payload: Mapping[str, Any],
    language: str = "zh",
    *,
    include_navigation: bool = True,
    analysis_warnings: Sequence[Any] | None = None,
) -> str:
    """Render a focused targeted-library result section."""

    language = "en" if language == "en" else "zh"
    summary = expected_summary_row(payload, language)
    rows = _expected_detail_rows(payload, language)
    supported_rows = [
        row for row in rows
        if int(row.get("_Supported_Feature_Count") or 0) > 0
    ]
    mass_only_rows = [
        row for row in rows
        if int(row.get("_Supported_Feature_Count") or 0) == 0
    ]
    ok = summary["Status"] == "ok"
    labels = {
        "zh": {
            "title": "靶向库匹配结果",
            "subtitle": (
                "基于所选 Expected Compounds 靶向数据库提取连续 MS1 "
                "XIC；与普通 MS2 谱库鉴定分开统计。"
            ),
            "catalog": "靶向库记录",
            "active_catalog": "进入计算的有效记录",
            "quarantined_catalog": "运行时隔离记录",
            "hypotheses": "离子假设",
            "features": "靶向化合物（已去重）",
            "supported": "有补充证据的化合物",
            "raw_features": "原始 XIC 峰（审计）",
            "warning": "MS1 模块适用性提示",
            "click_expand": "点击展开",
            "empty": "本次没有形成可导出的靶向库匹配结果。",
            "boundary": (
                "结果边界：目标理论质量与实测质量相符时，只能形成候选证据，不能自动视为"
                "唯一结构确证；普通 MS2 谱库鉴定数量不会在本节重复计算。"
            ),
            "complete": "靶向数据库检出结果",
            "complete_note": (
                "主表每个去重目标显示一个规范化化合物名称；候选数单独列出。"
                "每行可展开全部原始 XIC 峰，其他候选名称保留在详细审计 CSV。"
            ),
            "columns": [
                "#", "化合物名称", "候选数", "靶向数据库",
                "分子式 / 检测加合物", "离子 / 原始峰审计", "代表 m/z",
                "代表 RT / 全部 RT", "代表强度 / 面积", "最高证据",
            ],
            "raw_columns": [
                "#", "靶向峰 ID", "加合物", "理论→实测 m/z",
                "ppm", "RT (min)", "峰高", "峰面积", "证据",
            ],
            "supported_heading": "有补充证据支持的靶向化合物",
            "mass_heading": "仅精确质量的靶向候选",
            "mass_note": (
                "这部分只有目标 m/z 与色谱峰形，默认收起；不能据此确认"
                "具体化合物。"
            ),
            "raw_expand": "展开原始 XIC 峰",
            "ions": "离子",
            "peaks": "峰",
            "supported_peaks": "峰有正交支持",
            "recommendation": "建议",
            "nav_targeted": "靶向库结果",
            "nav_ms2": "MS2 谱库鉴定",
            "quality_gate": "靶标库质量门控",
            "quality_gate_note": (
                "源库记录保留不变；只有通过身份、分子式和可执行离子信息检查的记录进入计算。"
            ),
            "full_name": "查看完整名称",
        },
        "en": {
            "title": "Targeted Library Matches",
            "subtitle": (
                "Continuous-MS1 XIC screening against the selected Expected "
                "Compounds target databases, reported separately from MS2 IDs."
            ),
            "catalog": "Target-library records",
            "active_catalog": "Active records used",
            "quarantined_catalog": "Runtime-quarantined records",
            "hypotheses": "Ion hypotheses",
            "features": "Deduplicated targeted compounds",
            "supported": "Compounds with supporting evidence",
            "raw_features": "Raw XIC peaks (audit)",
            "warning": "MS1 module applicability notes",
            "click_expand": "click to expand",
            "empty": "No exportable targeted-library matches were produced.",
            "boundary": (
                "Result boundary: exact mass alone is candidate evidence and does "
                "not establish a uniquely confirmed structure. Ordinary MS2 "
                "library identifications are not counted again here."
            ),
            "complete": "Detected targeted compounds",
            "complete_note": (
                "The main table shows one normalized compound name per deduplicated "
                "target and reports candidate count separately. Expand a row to audit "
                "every raw XIC peak; other candidate names remain in the audit CSV."
            ),
            "columns": [
                "#", "Compound name", "Candidate count", "Target database",
                "Formula / detected adducts", "Ions / raw-peak audit",
                "Representative m/z", "Representative RT / all RTs",
                "Representative intensity / area", "Highest evidence",
            ],
            "raw_columns": [
                "#", "Target peak ID", "Adduct", "Theoretical→observed m/z",
                "ppm", "RT (min)", "Peak height", "Peak area", "Evidence",
            ],
            "supported_heading": "Targeted compounds with supporting evidence",
            "mass_heading": "Exact-mass-only targeted candidates",
            "mass_note": (
                "These rows have target m/z and chromatographic peak shape only. "
                "They are collapsed by default and cannot establish a specific compound."
            ),
            "raw_expand": "Expand raw XIC peaks",
            "ions": "ions",
            "peaks": "peaks",
            "supported_peaks": "peaks with orthogonal support",
            "recommendation": "Recommended action",
            "nav_targeted": "Targeted library matches",
            "nav_ms2": "MS2 library identifications",
            "quality_gate": "Target-catalog quality gate",
            "quality_gate_note": (
                "Source rows are preserved; only records that pass identity, formula, "
                "and executable-ion checks enter computation."
            ),
            "full_name": "View full name",
        },
    }[language]
    warnings = _normalized_analysis_warnings(analysis_warnings)
    expected_warning = str(summary.get("Warning") or "").strip()
    if expected_warning:
        warnings = _normalized_analysis_warnings([
            *warnings,
            {
                "modules": ["Expected Compounds"],
                "code": "expected_compounds_status",
                "message": expected_warning,
            },
        ])
    evidence = summary.get("Compound_Evidence_Counts") or {}
    evidence_text = " · ".join(
        f"{html.escape(EVIDENCE_LABELS[language][key])} "
        f"{int(evidence.get(key) or 0)}"
        for key in EVIDENCE_LABELS[language]
    )
    quality_audit = summary.get("Target_Catalog_Quality_Audit") or {}
    reason_labels = {
        "zh": {
            "identity_status:rejected_identity_formula_mismatch": "已拒绝的名称/分子式身份冲突",
            "invalid_nonempty_molecular_formula": "无效分子式",
            "no_usable_formula_exact_mass_or_precursor": "缺少可用分子式、精确质量和前体",
        },
        "en": {
            "identity_status:rejected_identity_formula_mismatch": "rejected name/formula identity conflicts",
            "invalid_nonempty_molecular_formula": "invalid molecular formulae",
            "no_usable_formula_exact_mass_or_precursor": "no usable formula, exact mass, or precursor",
        },
    }[language]
    quality_reason_text = " · ".join(
        f"{reason_labels.get(str(reason), str(reason))} {int(count or 0)}"
        for reason, count in sorted(
            (quality_audit.get("quarantine_reason_counts") or {}).items(),
            key=lambda item: (-int(item[1] or 0), str(item[0])),
        )
        if int(count or 0) > 0
    )

    def raw_audit_table(row: Mapping[str, Any]) -> str:
        raw_rows = list(row.get("_Raw_Rows") or [])
        raw_body = []
        for raw_index, raw in enumerate(raw_rows, 1):
            raw_body.append(
                "<tr>"
                f'<td style="text-align:center">{raw_index}</td>'
                f"<td>{html.escape(str(raw.get('Feature_ID') or '—'))}</td>"
                f"<td>{html.escape(str(raw.get('Adduct') or '—'))}</td>"
                f'<td style="font-family:monospace">{_fmt(raw.get("Target_mz"), 6)}'
                f' → {_fmt(raw.get("Observed_Apex_mz"), 6)}</td>'
                f'<td style="font-family:monospace">{_fmt(raw.get("Mass_Error_ppm"), 2)}</td>'
                f'<td style="font-family:monospace">{_fmt(raw.get("Apex_RT_min"), 3)}</td>'
                f'<td style="font-family:monospace">{_fmt(raw.get("Peak_Height"), 1)}</td>'
                f'<td style="font-family:monospace">{_fmt(raw.get("Peak_Area"), 1)}</td>'
                f"<td>{html.escape(str(raw.get('Evidence_Label') or '—'))}</td>"
                "</tr>"
            )
        raw_headers = "".join(
            f"<th>{html.escape(column)}</th>"
            for column in labels["raw_columns"]
        )
        return (
            '<div style="overflow:auto;max-height:300px;margin-top:6px;'
            'border:1px solid #dbeafe;border-radius:6px;background:white">'
            '<table style="width:100%;border-collapse:collapse;font-size:10px;'
            'min-width:980px"><thead style="background:#eaf2ff;color:#1e3a8a">'
            f"<tr>{raw_headers}</tr></thead><tbody>"
            + "".join(raw_body)
            + "</tbody></table></div>"
        )

    def group_body(group_rows: Sequence[Mapping[str, Any]]) -> str:
        body_rows: list[str] = []
        for index, row in enumerate(group_rows, 1):
            names = [
                str(value).strip()
                for value in (row.get("_Candidate_Names") or [])
                if str(value).strip()
            ]
            name_text = names[0] if names else ""
            if not name_text:
                name_text = str(row.get("Candidate_Name") or "—")
            escaped_name = html.escape(name_text)
            if len(name_text) > 56:
                compact_name = html.escape(name_text[:53].rstrip()) + "…"
                name_cell = (
                    f'<b title="{html.escape(name_text, quote=True)}">'
                    f"{compact_name}</b>"
                    '<details class="expected-full-compound-name" '
                    'style="margin-top:4px;font-size:10px">'
                    '<summary style="cursor:pointer;color:#1d4ed8;'
                    'font-weight:600">'
                    f"{html.escape(labels['full_name'])}</summary>"
                    '<div style="margin-top:4px;padding:5px 6px;'
                    'border:1px solid #dbeafe;border-radius:4px;'
                    'background:#f8fbff;overflow-wrap:anywhere">'
                    f"{escaped_name}</div></details>"
                )
            else:
                name_cell = f"<b>{escaped_name}</b>"
            identity_count = int(row.get("Candidate_Count") or 1)
            evidence_level = str(row.get("Evidence_Level") or "")
            evidence_color = {
                "EXPERIMENTAL_LIBRARY_CONFIRMED": "#166534",
                "MS1_MS2_DIAGNOSTIC_PROBABLE": "#1d4ed8",
                "MS1_ISOTOPE_SUPPORTED": "#92400e",
                "MS1_EXACT_MASS_CANDIDATE": "#475569",
            }.get(evidence_level, "#475569")
            rt_range = (
                f"{_fmt(row.get('_RT_Min'), 3)}–"
                f"{_fmt(row.get('_RT_Max'), 3)}"
            )
            raw_count = int(row.get("Raw_Peak_Count") or 0)
            ion_count = int(row.get("_Ion_Hypothesis_Count") or 0)
            supported_count = int(
                row.get("_Supported_Feature_Count") or 0
            )
            raw_details = (
                '<details class="expected-raw-feature-audit" '
                f'data-raw-peak-count="{raw_count}" style="margin-top:4px">'
                '<summary style="cursor:pointer;color:#1d4ed8;font-size:10px;'
                'font-weight:700">'
                f"{html.escape(labels['raw_expand'])} ({raw_count})"
                "</summary>"
                f"{raw_audit_table(row)}</details>"
            )
            body_rows.append(
                "<tr>"
                f'<td style="text-align:center">{index}</td>'
                f'<td style="line-height:1.35">{name_cell}</td>'
                f'<td style="text-align:center">{identity_count}</td>'
                f"<td>{html.escape(str(row.get('Target_Database') or '—'))}</td>"
                f"<td>{html.escape(str(row.get('Formula') or '—'))}<br>"
                f"<small>{html.escape(str(row.get('Related_Adducts') or row.get('Adduct') or '—'))}</small></td>"
                f"<td><b>{ion_count}</b> {html.escape(labels['ions'])} / "
                f"<b>{raw_count}</b> {html.escape(labels['peaks'])}<br>"
                f"<small>{supported_count} {html.escape(labels['supported_peaks'])}</small>"
                f"{raw_details}</td>"
                f'<td style="font-family:monospace">{_fmt(row.get("Target_mz"), 6)}<br>'
                f'<small>→ {_fmt(row.get("Observed_Apex_mz"), 6)} · '
                f'{_fmt(row.get("Mass_Error_ppm"), 2)} ppm</small></td>'
                f'<td style="font-family:monospace">{_fmt(row.get("Apex_RT_min"), 3)}'
                f"<br><small>{rt_range}</small></td>"
                f'<td style="font-family:monospace">{_fmt(row.get("Peak_Height"), 1)}'
                f'<br><small>A={_fmt(row.get("Peak_Area"), 1)}</small></td>'
                '<td><span style="display:inline-block;padding:3px 7px;'
                'border-radius:999px;background:#fff;'
                f'color:{evidence_color};border:1px solid {evidence_color};'
                'font-weight:700;white-space:nowrap">'
                f"{html.escape(str(row.get('Evidence_Label') or '—'))}</span></td>"
                "</tr>"
            )
        return "".join(body_rows)

    def group_table(group_rows: Sequence[Mapping[str, Any]]) -> str:
        headers = "".join(
            f'<th style="padding:7px 6px;text-align:left">'
            f"{html.escape(column)}</th>"
            for column in labels["columns"]
        )
        return (
            '<div style="overflow:auto;max-height:900px;margin-top:10px">'
            '<table style="width:100%;border-collapse:collapse;font-size:12px;'
            'min-width:1380px;table-layout:fixed;overflow-wrap:anywhere">'
            '<colgroup><col style="width:52px"><col style="width:210px">'
            '<col style="width:72px"><col style="width:150px">'
            '<col style="width:126px"><col style="width:210px">'
            '<col style="width:158px"><col style="width:102px">'
            '<col style="width:140px"><col style="width:160px"></colgroup>'
            f'<thead style="background:#1d4ed8;color:white"><tr>{headers}</tr></thead>'
            '<tbody style="vertical-align:top">'
            + group_body(group_rows)
            + "</tbody></table></div>"
        )

    if rows:
        table_html = ""
        if supported_rows:
            table_html += (
                '<details open class="expected-supported-groups" '
                'style="margin-top:12px">'
                '<summary style="cursor:pointer;font-weight:700;color:#1e3a8a">'
                f"{html.escape(labels['supported_heading'])} "
                f"({len(supported_rows)})</summary>"
                f"{group_table(supported_rows)}</details>"
            )
        if mass_only_rows:
            table_html += (
                '<details class="expected-mass-only-groups" '
                'style="margin-top:12px;padding:9px 10px;border:1px solid '
                '#cbd5e1;border-radius:7px;background:white">'
                '<summary style="cursor:pointer;font-weight:700;color:#475569">'
                f"{html.escape(labels['mass_heading'])} "
                f"({len(mass_only_rows)}) — "
                f"{html.escape(labels['click_expand'])}</summary>"
                '<p style="margin:6px 0 0;color:#92400e;font-size:11px">'
                f"{html.escape(labels['mass_note'])}</p>"
                f"{group_table(mass_only_rows)}</details>"
            )
    else:
        table_html = (
            '<p style="color:#64748b;margin-top:12px">'
            f"{html.escape(labels['empty'])}</p>"
        )

    warning_html = ""
    if warnings:
        warning_rows = []
        for warning in warnings:
            module = warning.get("modules") or warning.get("code") or ""
            warning_rows.append(
                '<div style="padding:7px 0;border-top:1px solid #fde68a">'
                f"<b>{html.escape(module)}</b>"
                f"<div>{html.escape(warning['message'])}</div>"
                + (
                    '<div style="margin-top:2px;color:#a16207">'
                    f"{html.escape(labels['recommendation'])}："
                    f"{html.escape(warning['action'])}</div>"
                    if warning.get("action") else ""
                )
                + "</div>"
            )
        warning_html = (
            '<details id="ms1-module-warnings" style="margin-top:10px;'
            'padding:8px 10px;border-left:4px solid #d97706;'
            'background:#fffbeb;color:#92400e;font-size:12px;line-height:1.55">'
            '<summary style="cursor:pointer;font-weight:700">'
            f"⚠ {html.escape(labels['warning'])}（{len(warnings)}）"
            '<span style="margin-left:5px;font-weight:400;color:#a16207">'
            f"{html.escape(labels['click_expand'])}</span></summary>"
            '<div style="margin-top:4px">'
            + "".join(warning_rows)
            + "</div></details>"
        )

    status_color = "#166534" if ok else "#92400e"
    navigation = (
        f"""
  <nav aria-label="report sections"
       style="display:flex;gap:8px;flex-wrap:wrap;margin:0 0 14px">
    <a href="#expected-compounds-report"
       style="text-decoration:none;padding:6px 10px;border-radius:7px;background:#1d4ed8;color:white;font-weight:700">
      {html.escape(labels['nav_targeted'])}
    </a>
    <a href="#ms2-identification-report"
       style="text-decoration:none;padding:6px 10px;border-radius:7px;background:white;color:#1d4ed8;border:1px solid #93c5fd;font-weight:700">
      {html.escape(labels['nav_ms2'])}
    </a>
  </nav>"""
        if include_navigation else ""
    )
    return f"""
<section id="expected-compounds-report"
         style="margin:24px 0;padding:18px;border:1px solid #93c5fd;border-left:5px solid #2563eb;border-radius:10px;background:#f8fbff;color:#172033;">
  {navigation}
  <div style="display:flex;justify-content:space-between;gap:12px;align-items:flex-start;flex-wrap:wrap">
    <div>
      <h2 style="margin:0;color:#1e3a8a;font-size:20px">{html.escape(labels['title'])}</h2>
      <p style="margin:4px 0 0;color:#475569;font-size:13px">{html.escape(labels['subtitle'])}</p>
    </div>
    <span style="padding:4px 10px;border-radius:999px;background:#fff7ed;color:{status_color};font-weight:700">
      {html.escape(str(summary['Status_Label']))}
    </span>
  </div>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:9px;margin-top:14px">
    <div style="background:white;border:1px solid #dbeafe;border-radius:8px;padding:10px"><b>{summary['Target_Catalog_Rows']}</b><br><small>{html.escape(labels['catalog'])}</small></div>
    <div style="background:#f0fdf4;border:1px solid #86efac;border-radius:8px;padding:10px"><b>{summary['Active_Target_Rows']}</b><br><small>{html.escape(labels['active_catalog'])}</small></div>
    <div style="background:{'#fffbeb' if summary['Quarantined_Target_Rows'] else '#f0fdf4'};border:1px solid {'#fbbf24' if summary['Quarantined_Target_Rows'] else '#86efac'};border-radius:8px;padding:10px"><b>{summary['Quarantined_Target_Rows']}</b><br><small>{html.escape(labels['quarantined_catalog'])}</small></div>
    <div style="background:white;border:1px solid #dbeafe;border-radius:8px;padding:10px"><b>{summary['Expected_Ion_Hypotheses']}</b><br><small>{html.escape(labels['hypotheses'])}</small></div>
    <div style="background:#eff6ff;border:2px solid #2563eb;border-radius:8px;padding:10px"><b>{summary['Compound_Summary_Rows']}</b><br><small>{html.escape(labels['features'])}</small></div>
    <div style="background:white;border:1px solid #dbeafe;border-radius:8px;padding:10px"><b>{summary['Supported_Compounds']}</b><br><small>{html.escape(labels['supported'])}</small></div>
    <div style="background:white;border:1px solid #dbeafe;border-radius:8px;padding:10px"><b>{summary['Detected_Features']}</b><br><small>{html.escape(labels['raw_features'])}</small></div>
  </div>
  <p style="margin:12px 0 0;padding:9px 11px;background:{'#fffbeb' if summary['Quarantined_Target_Rows'] else '#f0fdf4'};border-left:4px solid {'#d97706' if summary['Quarantined_Target_Rows'] else '#16a34a'};border-radius:7px;font-size:12px;line-height:1.55">
    <b>{html.escape(labels['quality_gate'])}：</b>{html.escape(labels['quality_gate_note'])}
    {html.escape((' ' + quality_reason_text + '.') if quality_reason_text else '')}
  </p>
  <p style="margin:12px 0 0;padding:9px 11px;background:#eff6ff;border-radius:7px;font-size:12px;line-height:1.55">
    {evidence_text}. {html.escape(labels['boundary'])}
  </p>
  {warning_html}
  <h3 style="margin:15px 0 0;color:#1e3a8a;font-size:15px">{html.escape(labels['complete'])}</h3>
  <p style="margin:3px 0;color:#64748b;font-size:12px">{html.escape(labels['complete_note'])}</p>
  {table_html}
</section>
"""


def inject_expected_html(
    html_content: str,
    payload: Mapping[str, Any],
    language: str = "zh",
    *,
    analysis_warnings: Sequence[Any] | None = None,
) -> str:
    """Insert targeted results before the ordinary MS2 table, not at the end."""

    section = render_expected_html_section(
        payload,
        language,
        include_navigation=True,
        analysis_warnings=analysis_warnings,
    )
    ms2_anchor = '<div id="ms2-identification-report"></div>'
    lower = html_content.lower()
    for marker in (
        '<div class="toolbar"',
        '<section class="table-box"',
        '<div class="table-container"',
    ):
        position = lower.find(marker)
        if position >= 0:
            return (
                html_content[:position]
                + section
                + ms2_anchor
                + html_content[position:]
            )
    for marker in ("</main>", "</body>"):
        position = lower.rfind(marker)
        if position >= 0:
            return (
                html_content[:position]
                + section
                + ms2_anchor
                + html_content[position:]
            )
    return html_content + section + ms2_anchor


def render_expected_html_report(
    payload: Mapping[str, Any],
    source_name: str,
    language: str = "zh",
    *,
    analysis_warnings: Sequence[Any] | None = None,
) -> str:
    """Return a self-contained HTML report containing only targeted matches."""

    language = "en" if language == "en" else "zh"
    if language == "zh":
        page_title = "靶向库匹配报告"
        source_label = "分析文件"
        print_label = "打印 / 另存为 PDF"
        footer = (
            "本报告仅包含 Expected Compounds 靶向库检索结果；"
            "普通 MS2 谱库鉴定请下载相应 MS2 报告。"
        )
    else:
        page_title = "Targeted Library Match Report"
        source_label = "Analysis file"
        print_label = "Print / Save as PDF"
        footer = (
            "This report contains only Expected Compounds targeted-library "
            "results. Download the MS2 report for ordinary spectral-library IDs."
        )
    section = render_expected_html_section(
        payload,
        language,
        include_navigation=False,
        analysis_warnings=analysis_warnings,
    )
    return f"""<!doctype html>
<html lang="{'en' if language == 'en' else 'zh-CN'}">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>{html.escape(source_name)} — {html.escape(page_title)}</title>
  <style>
    *{{box-sizing:border-box}}
    body{{margin:0;background:#f3f6fa;color:#172033;font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}}
    main{{max-width:1700px;margin:auto;padding:24px}}
    header{{display:flex;justify-content:space-between;gap:18px;align-items:center;flex-wrap:wrap;padding:22px 26px;border-radius:12px;background:linear-gradient(135deg,#1e3a8a,#2563eb);color:white}}
    header h1{{margin:0 0 5px;font-size:25px}}
    header p{{margin:0;opacity:.9}}
    button{{border:1px solid #ffffff88;border-radius:7px;background:#ffffff1f;color:white;padding:8px 12px;font-weight:700;cursor:pointer}}
    td,th{{border-bottom:1px solid #dbeafe;padding:7px 6px;vertical-align:top}}
    footer{{padding:18px;color:#64748b;text-align:center}}
    @media print{{body{{background:white}}main{{max-width:none;padding:0}}button{{display:none}}details>*{{display:block!important}}summary{{list-style:none}}}}
  </style>
</head>
<body>
<main>
  <header>
    <div>
      <h1>{html.escape(page_title)}</h1>
      <p>{html.escape(source_label)}：{html.escape(source_name)}</p>
    </div>
    <button onclick="window.print()">{html.escape(print_label)}</button>
  </header>
  {section}
  <footer>{html.escape(footer)}</footer>
</main>
</body>
</html>"""


def expected_pdf_story(
    payload: Mapping[str, Any],
    language: str = "zh",
    *,
    max_rows: int = 200,
    analysis_warnings: Sequence[Any] | None = None,
) -> list[Any]:
    """Build a compact PDF story; full detail remains in Excel/CSV."""

    from reportlab.lib import colors
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import Paragraph, Spacer, Table, TableStyle

    language = "en" if language == "en" else "zh"
    summary = expected_summary_row(payload, language)
    rows = _expected_detail_rows(payload, language)
    styles = getSampleStyleSheet()
    regular_font, bold_font = reportlab_font_names(language)
    for style_name in ("Normal", "Title", "Heading1", "Heading2"):
        styles[style_name].fontName = regular_font
    note = ParagraphStyle(
        f"ExpectedNote_{language}",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=7,
        leading=9,
        textColor=colors.HexColor("#334155"),
    )
    cell_text = ParagraphStyle(
        f"ExpectedCell_{language}",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=6.2,
        leading=7.2,
        textColor=colors.HexColor("#172033"),
    )
    if language == "zh":
        title = "靶向库匹配结果"
        status_text = (
            f"状态：{summary['Status_Label']}；靶向库记录 "
            f"{summary['Target_Catalog_Rows']}（有效 "
            f"{summary['Active_Target_Rows']}，运行时隔离 "
            f"{summary['Quarantined_Target_Rows']}）；离子假设 "
            f"{summary['Expected_Ion_Hypotheses']}；原始 XIC 峰段 "
            f"{summary['Detected_Features']}；检出化合物 "
            f"{summary['Compound_Summary_Rows']}；具有补充证据 "
            f"{summary['Supported_Compounds']}。"
        )
        boundary = (
            "本节与普通 MS2 谱库鉴定分开统计。目标理论质量与实测质量相符时，只能形成候选"
            "证据，不能自动视为唯一结构确证。PDF 最多展示前 "
            f"{max_rows} 条；完整明细见 Excel 的分层靶向工作表或靶向 CSV。"
        )
        header = [
            "#", "化合物", "靶向数据库", "分子式/加合物", "目标 m/z",
            "ppm", "RT", "峰面积", "证据",
        ]
    else:
        title = "Targeted Library Matches"
        status_text = (
            f"Status: {summary['Status_Label']}; target-library records "
            f"{summary['Target_Catalog_Rows']} (active "
            f"{summary['Active_Target_Rows']}, runtime-quarantined "
            f"{summary['Quarantined_Target_Rows']}); ion hypotheses "
            f"{summary['Expected_Ion_Hypotheses']}; raw XIC segments "
            f"{summary['Detected_Features']}; detected compounds "
            f"{summary['Compound_Summary_Rows']}; with supporting evidence "
            f"{summary['Supported_Compounds']}."
        )
        boundary = (
            "This section is counted separately from ordinary MS2 library "
            "identifications. Exact mass alone does not establish a unique structure. "
            f"The PDF shows at most {max_rows} rows; complete detail is in the "
            "layered targeted Excel sheets or targeted CSV."
        )
        header = [
            "#", "Compound", "Target database", "Formula/adduct", "Target m/z",
            "ppm", "RT", "Peak area", "Evidence",
        ]

    story: list[Any] = [
        Spacer(1, 0.4 * cm),
        Paragraph(title, styles["Heading2"]),
        Paragraph(html.escape(status_text), note),
        Paragraph(html.escape(boundary), note),
    ]
    warnings = _normalized_analysis_warnings(analysis_warnings)
    expected_warning = str(summary.get("Warning") or "").strip()
    if expected_warning:
        warnings = _normalized_analysis_warnings([
            *warnings,
            {
                "modules": ["Expected Compounds"],
                "code": "expected_compounds_status",
                "message": expected_warning,
            },
        ])
    if warnings:
        warning_title = (
            "MS1 模块适用性提示"
            if language != "en"
            else "MS1 module applicability notes"
        )
        action_label = "建议" if language != "en" else "Action"
        warning_text = "；".join(
            " ".join(filter(None, [
                (
                    f"{warning['modules']}:"
                    if warning.get("modules") else ""
                ),
                warning["message"],
                (
                    f"{action_label}: {warning['action']}"
                    if warning.get("action") else ""
                ),
            ]))
            for warning in warnings
        )
        story.append(Paragraph(
            f"<b>{html.escape(warning_title)} "
            f"({len(warnings)})：</b> {html.escape(warning_text)}",
            note,
        ))
    if not rows:
        return story

    table_data = [header]
    for index, row in enumerate(rows[:max_rows], 1):
        primary_name = str(
            row.get("Candidate_Display")
            or row.get("Candidate_Name")
            or ""
        )
        table_data.append([
            str(index),
            Paragraph(html.escape(primary_name[:60]), cell_text),
            Paragraph(html.escape(str(
                row.get("Target_Database")
                or row.get("Dataset_Sources")
                or ""
            )[:36]), cell_text),
            Paragraph(
                f"{html.escape(str(row['Formula'] or '')[:18])}<br/>"
                f"{html.escape(str(row['Adduct'] or '')[:16])}",
                cell_text,
            ),
            _fmt(row["Target_mz"], 5),
            _fmt(row["Mass_Error_ppm"], 2),
            _fmt(row["Apex_RT_min"], 3),
            _fmt(row["Peak_Area"], 1),
            Paragraph(
                html.escape(str(row["Evidence_Label"])[:36]),
                cell_text,
            ),
        ])
    table = Table(
        table_data,
        colWidths=[
            0.55 * cm,
            2.7 * cm,
            2.2 * cm,
            2.2 * cm,
            1.7 * cm,
            1.0 * cm,
            1.2 * cm,
            1.5 * cm,
            3.3 * cm,
        ],
        repeatRows=1,
    )
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1D4ED8")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), bold_font),
        ("FONTNAME", (0, 1), (-1, -1), regular_font),
        ("FONTSIZE", (0, 0), (-1, -1), 6.5),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#CBD5E1")),
        (
            "ROWBACKGROUNDS",
            (0, 1),
            (-1, -1),
            [colors.white, colors.HexColor("#EFF6FF")],
        ),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.extend([Spacer(1, 0.18 * cm), table])
    return story
