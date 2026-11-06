#!/usr/bin/env python3
"""Generate a deterministic V12 publication-report QA artifact.

This is a rendering regression fixture, not a replacement for re-analysis of
the original raw file.  It reuses the immutable saved-peak view for the user's
MS2 586678 case, then renders it with the current V12 report code.  Assertions
fail closed when classification, ion-form provenance, product-structure data,
targeted-table geometry, mirror-spectrum colours, or terminology regress.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.spectrum import clean_peaks
from web.app import _gen_rich
from web.expected_compounds_export import inject_expected_html


DEFAULT_SOURCE_REPORT = (
    PROJECT_ROOT.parent
    / "mssystem_V13"
    / "reports"
    / "V13_MS2_586678_detailed_report.html"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "reports"
    / "MSSYSTEM_V12_PUBLICATION_QA_REPORT_20260820.html"
)
DEFAULT_RECEIPT = (
    PROJECT_ROOT
    / "docs"
    / "audits"
    / "MSSYSTEM_V12_PUBLICATION_QA_REPORT_20260820.json"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_text_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def _write_json_atomic(path: Path, payload: dict) -> None:
    _write_text_atomic(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
    )


def _extract_saved_query_peaks(source_report: Path) -> list[list[float]]:
    content = source_report.read_text(encoding="utf-8")
    start_marker = "样品 MS/MS 完整保存峰表"
    end_marker = "Library 参考库谱峰"
    start = content.find(start_marker)
    end = content.find(end_marker, start + len(start_marker))
    if start < 0 or end < 0:
        raise ValueError("The source report does not contain the saved query-peak table")
    section = content[start:end]
    rows = re.findall(
        r'<td class="mono">([0-9.]+)</td><td>([0-9.]+)</td>'
        r'<td>[0-9.]+%</td>',
        section,
    )
    peaks = [[float(mz), float(intensity)] for mz, intensity in rows]
    if len(peaks) != 182:
        raise ValueError(
            f"Expected 182 immutable saved peaks for MS2 586678, found {len(peaks)}"
        )
    return peaks


def _match() -> dict:
    product_smiles = "CCCCCCCCCCCC[N+](C)C"
    dimensions = {
        "neutral_loss": {
            "status": "no_evidence",
            "score": 0.0,
            "evidence_count": 0,
            "reason": "no_current_policy_neutral_loss_assignment",
        },
        "rule_coverage": {
            "status": "scored",
            "score": 0.50,
            "evidence_count": 1,
            "reason": "tertiary_amine_oxide_structure_rule_applied",
        },
        "fg_signature": {
            "status": "scored",
            "score": 0.60,
            "evidence_count": 1,
            "reason": "candidate_contains_tertiary_amine_oxide_headgroup",
        },
        "feature_library": {
            "status": "not_applicable",
            "score": 0.0,
            "evidence_count": 0,
            "reason": "no_curated_rule_for_ontology_class",
        },
        "structure_fragment": {
            "status": "scored",
            "score": 0.75,
            "evidence_count": 1,
            "reason": "candidate_bond_fragment_matches_observed_product",
        },
        "diagnostic_ion": {
            "status": "not_applicable",
            "score": 0.0,
            "evidence_count": 0,
            "reason": "no_curated_rule_for_ontology_class",
            "reason_label_cn": (
                "候选已由结构规则归类为三级胺氧化物；当前发布的规则库尚无"
                "该类别的特异诊断离子规则。参考库常见碎片在下方独立复核。"
            ),
        },
        "bond_energy": {
            "status": "unavailable",
            "score": 0.0,
            "evidence_count": 0,
            "reason": "no_released_bond_energy_model_for_this_candidate",
        },
    }
    candidate = {
        "compound_name": "N,N-Dimethyldodecylamine N-oxide",
        "source_compound_name": "N,N-Dimethyldodecylamine N-oxide",
        "formula": "C14H31NO",
        "smiles": "CCCCCCCCCCCC[N+](C)(C)[O-]",
        "inchikey": "SYELZBGXAIXKHU-UHFFFAOYSA-N",
        "classification": "Unknown",
        "structure_quality": "validated",
        "dataset_source": "MSMS-Public_experimentspectra-VS19-Pos",
        "source_tier": "experimental_public_library",
        "lib_spectrum_id": "MSMS-Public_experimentspectra-pos-VS19_3017342",
        "lib_source_accession": "MassBank:WANA005215",
        "adduct": "[M+H]+",
        "cosine": 0.9564,
        "final_score": 0.7005,
        "rank_score": 0.7005,
        "mass_ppm": 4.04,
        "mass_score": 0.92,
        "fragment_score": 0.81,
        "intensity_score": 0.73,
        "expert_score": 0.225,
        "expert_score_source": "base_8d_legacy_audit_only",
        "confidence": "High",
        "confidence_class": "High",
        "isomer_group_confidence": "Unknown",
        "lib_peaks": [
            [57.0700, 29.0],
            [58.0650, 506.0],
            [62.0600, 64.0],
            [71.0860, 19.0],
            [212.2370, 38.0],
        ],
        "pathway_json": {
            "n_spectra": 10,
            "diagnostic_ions": [
                {"mz": 58.0640, "freq": 0.50, "recurring": True},
                {"mz": 57.0697, "freq": 0.40, "recurring": True},
                {"mz": 62.0603, "freq": 0.30, "recurring": True},
                {"mz": 71.0856, "freq": 0.20, "recurring": True},
                {"mz": 212.2361, "freq": 0.80, "recurring": True},
                {"mz": 213.2400, "freq": 0.50, "recurring": True},
            ],
            "neutral_losses": [],
        },
        "expert_details": {
            "classification_provenance": {
                "reliable": True,
                "effective_classification": "amine_oxide",
                "analysis_class": "amine_oxide",
                "reliability_reason": "smiles_structure_rule:tertiary_amine_oxide",
                "input_confidence": 1.0,
            },
            "dimension_assessments": dimensions,
            "neutral_loss_analysis": {
                "best_pathway": [],
                "neutral_losses": [],
                "unassigned_mass_differences": [{
                    "from_mz": 230.248772,
                    "to_mz": 213.242535,
                    "mass_difference": 17.006237,
                    "reason": "not_assigned_as_NH3_under_current_accurate_mass_policy",
                }],
            },
            "diagnostic_ion_analysis": {"matched_ions": []},
            "structure_fragment_analysis": {
                "plausible_fragments": [{
                    "observed_mz": 213.242535,
                    "fragment_smiles": product_smiles,
                    "fragment_mass": 213.245101,
                    "expected_mz": 213.245101,
                    "delta": 0.002566,
                    "adduct": "[M]+",
                    "assignment_type": "candidate_structure_hypothesis",
                }],
            },
            "score_fusion": {
                "selected_source": "base_8d",
                "components": {
                    "base_8d": {
                        "available": True,
                        "score": 0.225,
                        "label": "旧版规则汇总索引（仅供审计）",
                    }
                },
            },
            "chemical_reasoning": (
                "候选结构含长链三级胺氧化物头基；m/z 213.24254 可由候选依赖的"
                "N–O 断裂产物结构解释。该结构是假设，不构成气相机理确证。"
            ),
        },
    }
    candidate["top_candidates"] = [{
        key: candidate[key]
        for key in (
            "compound_name", "source_compound_name", "formula", "smiles",
            "inchikey", "classification", "structure_quality", "cosine",
        )
    }]
    return candidate


def _targeted_payload() -> dict:
    compound_name = "N,N-Dimethyldodecylamine N-oxide"
    return {
        "schema_version": 1,
        "task_id": "publication-qa-ms2-586678",
        "method": "Expected Compounds / continuous-MS1 targeted XIC",
        "status": "ok",
        "requested": True,
        "catalog": {
            "target_rows": 1,
            "expected_ion_hypotheses": 2,
            "dataset_sources": ["qa:publication_layout_exact_case"],
        },
        "summary": {
            "status": "ok",
            "target_catalog_rows": 1,
            "expected_ion_hypotheses": 2,
            "detected_feature_count": 1,
            "supported_features": 1,
            "mass_only_features": 0,
            "ambiguous_feature_count": 0,
            "candidate_target_rows": 1,
            "candidate_target_names": 1,
            "ms1_scan_count": 1560,
            "usable_ms1_scan_count": 1558,
            "evidence_counts": {
                "MS1_EXACT_MASS_CANDIDATE": 0,
                "MS1_ISOTOPE_SUPPORTED": 0,
                "MS1_MS2_DIAGNOSTIC_PROBABLE": 0,
                "EXPERIMENTAL_LIBRARY_CONFIRMED": 1,
            },
            "counted_in_standard_spectral_library_matches": False,
            "standard_spectral_library_match_increment": 0,
        },
        "features": [{
            "feature_id": "LF0000838",
            "candidate_name": compound_name,
            "candidate_names": [compound_name],
            "formula": "C14H31NO",
            "adduct": "[M+H]+",
            "polarity": "POS",
            "target_mz": 230.248227,
            "observed_apex_mz": 230.248278,
            "mass_error_ppm": 0.218,
            "apex_rt_minutes": 9.777833,
            "boundary_start_rt_minutes": 9.70,
            "boundary_end_rt_minutes": 9.86,
            "fwhm_seconds": 7.8,
            "apex_intensity": 3963.546,
            "xic_area_intensity_seconds": 18650.0,
            "signal_to_noise": 18.4,
            "scan_points": 12,
            "xic_nonzero_points": 12,
            "evidence_level": "EXPERIMENTAL_LIBRARY_CONFIRMED",
            "cd_pattern_accepted": True,
            "isomer_ambiguous": False,
            "identity_count": 1,
            "identification_scope": "experimental_library_supported_candidate",
            "formula_isotope_evidence": {
                "status": "supporting",
                "tolerance_ppm": 5.0,
                "isotope_trace_count": 2,
            },
            "associated_ms2": [{"spectrum_id": "586678"}],
            "ms2_diagnostic_support": [],
            "experimental_library_support": [{
                "library_spectrum_id": (
                    "MSMS-Public_experimentspectra-pos-VS19_3017342"
                ),
                "high_accuracy_fragment_tolerance_da": 0.0025,
                "matched_fragment_pairs": 4,
                "minimum_matched_fragment_pairs": 3,
                "library_intensity_coverage": 1.0,
                "minimum_library_intensity_coverage": 0.5,
            }],
            "target_members": [{
                "target_id": "qa-lauramine-oxide",
                "inchikey": "SYELZBGXAIXKHU-UHFFFAOYSA-N",
            }],
            "dataset_sources": ["qa:publication_layout_exact_case"],
            "generation_sources": ["formula_exact_mass"],
            "evidence_interpretation": (
                "Continuous-MS1 feature with isotope support and a high-accuracy "
                "experimental-library MS2 comparison; this is not same-method "
                "authentic-standard confirmation."
            ),
            "ordinary_library_match_counted": False,
        }],
        "scientific_boundaries": {
            "ordinary_library_matched_count_unchanged": True,
            "same_method_authentic_standard_confirmation": False,
        },
    }


def _render(source_report: Path) -> tuple[str, dict]:
    raw_peaks = _extract_saved_query_peaks(source_report)
    query_peaks = clean_peaks(
        raw_peaks,
        230.24877166748,
        precursor_charge=1,
    )
    row = {
        "precursor_mz": 230.24877166748,
        "spectrum_id": "586678",
        "precursor_scan_ref": "scanId=586078",
        "peaks_count": len(raw_peaks),
        "retention_time": 9.777833333333,
        "retention_time_unit": "minute",
        "collision_energy": 30.0,
        "collision_energy_unit": "eV",
        "ion_mode": "POS",
        "precursor_charge": 1,
        "tic": 15354.3,
        "base_peak": 5085.603,
        "raw_peaks": raw_peaks,
        "raw_peaks_complete": True,
        "query_peaks": query_peaks,
        "lcms_adduct": "[M+H]+",
        "lcms_neutral_mass": 229.241445,
        "lcms_feature_link_status": "linked_openms_feature",
        "matches": [_match()],
    }
    html = _gen_rich(
        [row],
        source_name="historical sample · spectrum 586678 · V12 QA",
        total_spectra=1,
        matched_spectra=1,
        high_supported_compounds=1,
    )
    html = inject_expected_html(html, _targeted_payload(), "zh")
    banner = """
<section id="publication-qa-scope" style="margin:0 0 16px;padding:12px 14px;border:1px solid #0f766e;border-left:5px solid #0f766e;border-radius:8px;background:#f0fdfa;color:#134e4a;font-size:12px;line-height:1.65">
  <b>V12 出版报告回归验收工件</b>：使用 MS2 586678 的不可变保存峰表和当前 V12 代码重新渲染，用于检查类别、样品离子形式、规则证据维度、候选产物结构、镜像谱配色与靶向表列宽。它不是对已删除原任务的重新分析，也不替代同方法标准品验证。
</section>
"""
    marker = '<div class="container">'
    html = html.replace(marker, marker + banner, 1)
    expected_query_label = (
        f"质控后用于匹配 {len(query_peaks)} 个峰；原始工件保存 "
        f"{len(raw_peaks)} 个峰"
    )
    assertions = {
        "exact_case_present": "N,N-Dimethyldodecylamine N-oxide" in html,
        "reliable_amine_oxide_class_present": "amine_oxide" in html,
        "sample_ion_form_present": "样品离子形式：[M+H]+" in html,
        "sample_ion_source_present": "lcms_ion_identity" in html,
        "product_structure_section_present": "观测峰支持的候选产物结构" in html,
        "product_structure_smiles_present": "CCCCCCCCCCCC[N+](C)C" in html,
        "product_structure_is_hypothesis": "candidate_structure_hypothesis" in html,
        "query_library_colors_differ": (
            'fill="#007A78"' in html and 'fill="#6F4E9B"' in html
        ),
        "shared_axis_mirror_present": (
            'data-spectrum-layout="mirror-shared-axis"' in html
        ),
        "processed_and_saved_peak_counts_separated": expected_query_label in html,
        "targeted_name_column_is_bounded": (
            '<col style="width:52px"><col style="width:210px">' in html
        ),
        "targeted_section_precedes_ms2": (
            html.find('id="expected-compounds-report"')
            < html.find('id="ms2-identification-report"')
        ),
        "legacy_expert_scalar_demoted_to_audit": (
            "旧式聚合值仅保存在机器可读审计字段中，不是概率，也不参与身份排序"
            in html
            and "旧式未校准汇总值 0.225" not in html
        ),
        "structured_dimension_scores_rendered": all(
            value in html for value in ("0.50 · 1项", "0.60 · 1项", "0.75 · 1项")
        ),
        "dimension_status_counts_are_explicit": (
            "有支持证据 3/7；已评价但未观察 1；不适用/未建立 3" in html
        ),
        "legacy_expert_label_absent": "Expert: 0.216" not in html,
        "legacy_classification_wording_absent": "当前类别体系未覆盖" not in html,
        "legacy_decision_wording_absent": "暂缓确认" not in html,
        "legacy_score_wording_absent": "技术分数" not in html,
        "false_ammonia_assignment_absent": "氨（NH3）" not in html,
    }
    return html, {
        "raw_saved_peak_count": len(raw_peaks),
        "processed_query_peak_count": len(query_peaks),
        "assertions": assertions,
        "all_assertions_passed": all(assertions.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE_REPORT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()

    source_report = args.source_report.expanduser().resolve()
    output = args.output.expanduser().resolve()
    receipt_path = args.receipt.expanduser().resolve()
    html, audit = _render(source_report)
    if not audit["all_assertions_passed"]:
        failures = [
            key for key, passed in audit["assertions"].items() if not passed
        ]
        raise RuntimeError(f"Publication QA assertions failed: {failures}")
    _write_text_atomic(output, html)
    receipt = {
        "schema_version": "1.1.0",
        "audit_id": "mssystem-v12-publication-report-qa-20260820",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "deterministic report-rendering regression; not raw-file re-analysis "
            "or same-method authentic-standard confirmation"
        ),
        "source_report": {
            "path": str(source_report),
            "sha256": _sha256(source_report),
        },
        "output_report": {
            "path": str(output),
            "sha256": _sha256(output),
            "bytes": output.stat().st_size,
        },
        **audit,
    }
    _write_json_atomic(receipt_path, receipt)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
