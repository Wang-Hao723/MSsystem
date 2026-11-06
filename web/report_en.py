"""Standalone English HTML report renderer.

The Chinese standard/detailed reports grew independently and contain a mixture
of Chinese labels and free-text expert reasoning.  This module deliberately
renders from the structured, serialized evidence instead of translating those
HTML strings.  Numeric scores and candidate identities therefore remain
identical to the task result, while all system-authored prose is English and
scientific limitations are stated explicitly.
"""
from __future__ import annotations

import html
import json
import math
from datetime import datetime
from typing import Any, Iterable, Optional

from web.report_common import (
    _standard_compound_name,
    classification_evidence_status,
    exact_structure_confidence,
    final_claim_badge_html,
    honest_annotation,
    isomer_candidate_projection,
    reportable_identity_confidence,
    reportable_identity_label,
    render_isomer_candidate_competition_html,
    svg_for_smiles,
    two_dimensional_connectivity_confidence,
)
from web.raw_spectrum_locator import (
    raw_spectrum_locator,
    raw_spectrum_locator_text,
)


_CLASS_EN = {
    "黄酮": "Flavonoid",
    "黄酮类": "Flavonoid",
    "生物碱": "Alkaloid",
    "生物碱类": "Alkaloid",
    "萜类": "Terpenoid",
    "甾体": "Steroid",
    "皂苷": "Saponin",
    "酚酸": "Phenolic acid",
    "香豆素": "Coumarin",
    "木脂素": "Lignan",
    "醌": "Quinone",
    "糖": "Sugar",
    "糖苷": "Glycoside",
    "脂质": "Lipid",
    "脂肪酸": "Fatty acyl",
    "氨基酸": "Amino acid",
    "核苷": "Nucleoside",
}


def _esc(value: Any) -> str:
    return html.escape(str(value if value is not None else ""), quote=True)


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _fmt(value: Any, digits: int = 3, missing: str = "—") -> str:
    number = _finite_float(value)
    return f"{number:.{digits}f}" if number is not None else missing


def _match_of(row: dict) -> Optional[dict]:
    if row.get("match"):
        return row["match"]
    matches = row.get("matches") or []
    return matches[0] if matches else None


def _expert_dimension_audit(match: dict) -> dict[str, int]:
    details = match.get("expert_details") or {}
    assessments = details.get("dimension_assessments") or {}
    counts = {
        "supporting": 0,
        "no_evidence": 0,
        "not_applicable": 0,
        "unavailable": 0,
    }
    if not isinstance(assessments, dict):
        assessments = {}
    for value in assessments.values():
        if not isinstance(value, dict):
            continue
        status = str(value.get("status") or "unavailable")
        if status == "scored" and (
            int(value.get("evidence_count") or 0) > 0
            or (_finite_float(value.get("score")) or 0.0) > 0
        ):
            counts["supporting"] += 1
        elif status == "no_evidence" or status == "scored":
            counts["no_evidence"] += 1
        elif status == "not_applicable":
            counts["not_applicable"] += 1
        else:
            counts["unavailable"] += 1
    counts["total"] = sum(counts.values())
    return counts


def _identity_key(row: dict) -> str:
    match = _match_of(row) or {}
    return str(
        match.get("inchikey")
        or match.get("smiles")
        or match.get("compound_name")
        or f"mz:{_fmt(row.get('precursor_mz'), 5)}"
    )


def _deduplicate(rows: Iterable[dict]) -> tuple[list[dict], list[dict], int]:
    matched = sorted(
        (row for row in rows if _match_of(row)),
        key=lambda row: (
            _finite_float(
                (_match_of(row) or {}).get(
                    "rank_score",
                    (_match_of(row) or {}).get("final_score"),
                )
            ) or 0.0,
            _finite_float((_match_of(row) or {}).get("final_score")) or 0.0,
            _finite_float((_match_of(row) or {}).get("cosine")) or 0.0,
        ),
        reverse=True,
    )
    unique: list[dict] = []
    seen: set[str] = set()
    for row in matched:
        key = _identity_key(row)
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)
    unmatched = [row for row in rows if not _match_of(row)]
    return unique, unmatched, len(matched) - len(unique)


def class_name_en(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() in {"other", "unknown", "n/a", "none"}:
        return "Unclassified"
    if text in _CLASS_EN:
        return _CLASS_EN[text]
    return text.replace("_", " ").strip().title()


# Internal shorthand used throughout the renderer.
_class_en = class_name_en


def _classification_en(match: dict) -> str:
    evidence = classification_evidence_status(match)
    label = _class_en(evidence["label"])
    if evidence["status"] == "tentative":
        return f"Tentative class clue: {label}"
    if evidence["status"] == "unresolved":
        return "Classification unresolved"
    return label


def _confidence_badge(match: dict) -> str:
    compound = reportable_identity_confidence(match)
    label = reportable_identity_label(match, language="en")
    full_identity = exact_structure_confidence(match)
    css = {
        "High": "high",
        "RT-confirmed High": "high",
        "Medium": "medium",
        "Strong Candidate": "medium",
        "Low": "low",
        "Putative": "low",
        "Unknown": "unknown",
    }.get(compound, "unknown")
    return (
        f'<span class="badge {css}" title="Reportable-scope confidence; '
        f'full-identity confidence: {_esc(full_identity)}">'
        f"{_esc(label)}</span>"
        + (
            f'<small class="scope-note">Full identity: {_esc(full_identity)}</small>'
            if compound != full_identity else ""
        )
    )


def _parse_pathway(value: Any) -> Optional[dict]:
    if isinstance(value, dict):
        return value
    if not value:
        return None
    try:
        parsed = json.loads(value)
    except (TypeError, ValueError):
        return None
    return parsed if isinstance(parsed, dict) else None


def _observed(mz: Any, peaks: list, tolerance: float = 0.05) -> bool:
    target = _finite_float(mz)
    if target is None:
        return False
    for peak in peaks or []:
        try:
            if abs(float(peak[0]) - target) <= tolerance:
                return True
        except (TypeError, ValueError, IndexError):
            continue
    return False


def _pathway_annotation(item: dict, diagnostic: bool) -> str:
    if diagnostic:
        # Stored diagnostic-ion labels in old databases can be Chinese.  The
        # curated stable key is language-neutral and safer for an English report.
        key = str(item.get("diag_key") or "").strip()
        if key:
            return key.replace("_", " ").title()
    annotation = str(item.get("annotation") or "").strip()
    if annotation and not any("\u4e00" <= ch <= "\u9fff" for ch in annotation):
        return annotation
    key = str(item.get("nl_key") or item.get("diag_key") or "").strip()
    return key.replace("_", " ").title() if key else "Unassigned recurring feature"


def _render_pathway(match: dict, query_peaks: list) -> str:
    pathway = _parse_pathway(match.get("pathway_json"))
    if not pathway:
        return (
            '<section class="evidence"><h4>Fragmentation evidence</h4>'
            '<p class="muted">No compound-specific library fragmentation signature is '
            "stored for this candidate.</p></section>"
        )
    diagnostics = pathway.get("diagnostic_ions") or []
    losses = pathway.get("neutral_losses") or []
    chains = pathway.get("loss_chains") or []
    n_spectra = pathway.get("n_spectra", "unknown")
    parts = [
        '<section class="evidence"><h4>Library-derived fragmentation signature</h4>',
        "<p class=\"method-note\">This is an empirical signature aggregated from "
        f"<b>{_esc(n_spectra)}</b> reference MS/MS spectrum/s for the candidate. "
        "Observed product ions and neutral-loss sequences support a match, but they "
        "are not atom-mapped bond-cleavage mechanisms and do not by themselves "
        "establish a uniquely confirmed molecular structure.</p>",
        f"<p><b>Polarity applicability:</b> "
        f"{_esc(str(match.get('pathway_polarity_status') or 'not recorded').replace('_', ' '))}; "
        f"<b>adduct applicability:</b> "
        f"{_esc(str(match.get('pathway_adduct_status') or 'not recorded').replace('_', ' '))}.</p>",
    ]
    if diagnostics:
        rows = []
        for item in diagnostics[:20]:
            rows.append(
                "<tr>"
                f"<td>{_fmt(item.get('mz'), 4)}</td>"
                f"<td>{_esc(_pathway_annotation(item, True))}</td>"
                f"<td>{_fmt((_finite_float(item.get('freq')) or 0) * 100, 0)}%</td>"
                f"<td>{_fmt(item.get('rel_intensity'), 1)}</td>"
                f"<td>{'Yes' if item.get('recurring') else 'Curated only'}</td>"
                f"<td>{'Observed' if _observed(item.get('mz'), query_peaks) else 'Not observed'}</td>"
                "</tr>"
            )
        parts.extend(
            [
                "<h5>Diagnostic product ions</h5>",
                '<table class="mini"><thead><tr><th>m/z</th><th>Assignment</th>'
                "<th>Reference frequency</th><th>Maximum relative intensity</th>"
                "<th>Evidence source</th><th>Query</th></tr></thead><tbody>",
                "".join(rows),
                "</tbody></table>",
            ]
        )
    if losses:
        rows = []
        for item in losses[:20]:
            rows.append(
                "<tr>"
                f"<td>{_fmt(item.get('loss_da'), 4)}</td>"
                f"<td>{_fmt(item.get('fragment_mz'), 4)}</td>"
                f"<td>{_esc(_pathway_annotation(item, False))}</td>"
                f"<td>{_fmt((_finite_float(item.get('freq')) or 0) * 100, 0)}%</td>"
                f"<td>{'Observed' if _observed(item.get('fragment_mz'), query_peaks) else 'Not observed'}</td>"
                "</tr>"
            )
        parts.extend(
            [
                "<h5>Neutral losses</h5>",
                '<table class="mini"><thead><tr><th>Loss (Da)</th><th>Product m/z</th>'
                "<th>Assignment</th><th>Reference frequency</th><th>Query</th>"
                "</tr></thead><tbody>",
                "".join(rows),
                "</tbody></table>",
            ]
        )
    if chains:
        chain_lines = []
        for chain in chains[:8]:
            values = " &rarr; ".join(_fmt(value, 4) for value in (chain or []))
            if values:
                chain_lines.append(f"<li>{values}</li>")
        if chain_lines:
            parts.extend(
                [
                    "<h5>Observed reference product-ion sequences</h5>",
                    "<ul class=\"compact\">",
                    "".join(chain_lines),
                    "</ul>",
                ]
            )
    if not (diagnostics or losses or chains):
        parts.append('<p class="muted">The stored pathway record contains no trusted features.</p>')
    parts.append("</section>")
    return "".join(parts)


def _render_isotope(row: dict, match: dict) -> str:
    evidence = match.get("formula_isotope_evidence") or {}
    if not evidence:
        status = str(row.get("isotope_evidence_status") or "not available").replace("_", " ")
        return (
            '<section class="evidence"><h4>MS1 isotope evidence</h4>'
            f'<p class="muted">Formula-envelope comparison: {_esc(status)}. '
            "Isotope evidence requires a genuine adjacent MS1 scan; product-ion "
            "peaks are never treated as a precursor isotope envelope.</p></section>"
        )
    observed = evidence.get("observed") or {}
    theoretical = evidence.get("theoretical") or {}
    adjustment = match.get("isotope_rank_adjustment")
    return (
        '<section class="evidence"><h4>MS1 isotope evidence</h4>'
        f"<p><b>Status:</b> {_esc(evidence.get('status', 'unknown'))} &nbsp; "
        f"<b>Consistency score:</b> {_fmt(evidence.get('score'), 3)} &nbsp; "
        f"<b>Bounded rank adjustment:</b> {_fmt(adjustment, 5)}</p>"
        f"<p>M+1/M observed {_fmt(observed.get('m1_over_m'), 5)} vs theoretical "
        f"{_fmt(theoretical.get('m1_over_m'), 5)}; M+2/M observed "
        f"{_fmt(observed.get('m2_over_m'), 5)} vs theoretical "
        f"{_fmt(theoretical.get('m2_over_m'), 5)}.</p>"
        '<p class="method-note">This is formula-level supporting evidence with a '
        "maximum ranking influence of 0.01. Co-elution and isolation interference "
        "can distort isotope ratios; the final score and confidence tier are not "
        "inflated by this term.</p></section>"
    )


def _render_adduct(match: dict) -> str:
    query_adduct = str(match.get("query_adduct") or "").strip()
    library_adduct = str(match.get("adduct") or "").strip()
    candidate_hypothesis = str(
        match.get("candidate_adduct_hypothesis") or ""
    ).strip()
    reliable = bool(match.get("query_adduct_reliable"))
    basis = str(match.get("mass_match_basis") or "precursor m/z").replace("_", " ")
    cross = bool(match.get("cross_adduct_match"))
    statement = (
        "The query adduct was parsed from explicit input metadata and can support "
        "neutral-mass matching."
        if reliable
        else "No chemically specific query adduct was available. A single precursor "
        "m/z is insufficient to distinguish protonated, sodiated, ammoniated, or "
        "other adduct forms, so the system does not claim an automatic adduct identity."
    )
    return (
        '<section class="evidence"><h4>Adduct and neutral-mass evidence</h4>'
        f"<p><b>Query adduct:</b> {_esc(query_adduct or 'Not assigned')} &nbsp; "
        f"<b>Library adduct:</b> {_esc(library_adduct or 'Not recorded')} &nbsp; "
        f"<b>Candidate ion-form hypothesis:</b> "
        f"{_esc(candidate_hypothesis or 'Not inferred')} &nbsp; "
        f"<b>Mass-match basis:</b> {_esc(basis)} &nbsp; "
        f"<b>Cross-adduct match:</b> {'Yes' if cross else 'No'}</p>"
        f'<p class="method-note">{_esc(statement)} '
        "A candidate ion-form hypothesis, when shown, is inferred from that "
        "candidate's structure/formula and library precursor metadata; it is not "
        "independent sample adduct evidence.</p></section>"
    )


def _render_ce_profile(match: dict) -> str:
    evidence = match.get("ce_profile_evidence") or {}
    status = str(evidence.get("status") or "not available").replace("_", " ")
    best = evidence.get("best_profile") or {}
    if not best:
        return (
            '<section class="evidence"><h4>Multi-source CE evidence</h4>'
            f"<p><b>Status:</b> {_esc(status)}</p>"
            '<p class="method-note">No directly comparable, structure- and '
            "adduct-matched multi-collision-energy profile was available. Bare "
            "collision-energy numbers are not treated as eV.</p></section>"
        )
    trajectory_rows = []
    for item in best.get("matched_fragment_trajectories") or []:
        trajectory_rows.append(
            "<tr>"
            f"<td>{_fmt(item.get('fragment_mz'), 5)}</td>"
            f"<td>{_esc(item.get('trend') or 'unknown')}</td>"
            f"<td>{_esc(item.get('n_levels') or 0)}</td>"
            f"<td>{_fmt(item.get('apex_energy'), 2)}</td>"
            "</tr>"
        )
    trajectory_table = ""
    if trajectory_rows:
        trajectory_table = (
            '<table class="mini"><thead><tr><th>Fragment m/z</th>'
            "<th>CE trend</th><th>Energy levels</th><th>Apex CE</th>"
            f"</tr></thead><tbody>{''.join(trajectory_rows)}</tbody></table>"
        )
    return (
        '<section class="evidence"><h4>Multi-source CE evidence</h4>'
        f"<p><b>Status:</b> {_esc(status)} &nbsp; "
        f"<b>Query CE:</b> {_fmt(evidence.get('query_collision_energy'), 2)} "
        f"{_esc(evidence.get('query_collision_energy_unit') or '')} &nbsp; "
        f"<b>Nearest library CE:</b> {_fmt(best.get('closest_energy'), 2)} "
        f"{_esc(best.get('collision_energy_unit') or '')} &nbsp; "
        f"<b>CE difference:</b> {_fmt(best.get('energy_delta'), 2)}</p>"
        f"<p><b>Source dataset:</b> "
        f"{_esc(best.get('source_dataset') or 'not recorded')} &nbsp; "
        f"<b>Identity provenance:</b> "
        f"{_esc(best.get('identity_verification') or 'not recorded')}</p>"
        f"<p><b>Instrument domain:</b> {_esc(best.get('instrument_type') or 'unknown')} / "
        f"{_esc(best.get('instrument_model') or 'unknown')}<br>"
        f"<b>Profile:</b> {_esc(best.get('n_energy_levels') or 0)} energy levels, "
        f"{_esc(best.get('n_profile_spectra') or 0)} spectra, quality "
        f"{_esc(best.get('quality_grade') or 'unknown')}<br>"
        f"<b>Nearest-energy consensus cosine:</b> "
        f"{_fmt(best.get('consensus_cosine'), 4)} &nbsp; "
        f"<b>entropy:</b> {_fmt(best.get('consensus_entropy'), 4)}</p>"
        f"{trajectory_table}"
        '<p class="method-note">This evidence uses exact full InChIKey, adduct, '
        "polarity, source dataset, instrument and CE-unit isolated summaries. It is "
        "shadow-only: score contribution = 0 and it did not change final score, "
        "rank or confidence.</p></section>"
    )


def _render_isomer(match: dict) -> str:
    ann = honest_annotation(match)
    projection = isomer_candidate_projection(match)
    ambiguous = projection.get("status") != "not_applicable"
    candidates = projection.get("candidate_members") or []
    isomer_set = candidates[1:]
    detail = match.get("isomer_fragment_detail") or {}
    rank_applied = bool(
        match.get("isomer_fragment_rank_applied")
        or detail.get("rank_override_enabled")
    )
    rows = []

    def shadow_score(candidate: dict, key: str) -> str:
        metric = ((candidate.get("multimetric_shadow") or {}).get(key) or {})
        return _fmt(metric.get("score"), 4)

    for index, candidate in enumerate(candidates, 1):
        candidate_detail = candidate.get("isomer_fragment_detail") or {}
        diagnostic = (
            f"{candidate_detail.get('unique_supported', 0)}/"
            f"{candidate_detail.get('unique_total', 0)} product ions; "
            f"{candidate_detail.get('neutral_loss_supported', 0)}/"
            f"{candidate_detail.get('neutral_loss_total', 0)} neutral losses"
            if candidate_detail else "—"
        )
        rows.append(
            "<tr>"
            f"<td>{index}<br>{_esc('ranking lead (not uniquely resolved)' if index == 1 else 'competitive isomer')}</td>"
            f"<td>{_esc(candidate.get('compound_name') or 'Unnamed candidate')}</td>"
            f"<td>{_esc(candidate.get('formula') or '—')}</td>"
            f"<td class=\"mono\">{_esc(candidate.get('inchikey') or '—')}</td>"
            f"<td>{_fmt(candidate.get('cosine'), 4)}</td>"
            f"<td>{shadow_score(candidate, 'direct_cosine_hungarian')}</td>"
            f"<td>{shadow_score(candidate, 'modified_cosine_hungarian')}</td>"
            f"<td>{shadow_score(candidate, 'neutral_loss_cosine')}</td>"
            f"<td>{_fmt(candidate.get('final_score'), 4)}</td>"
            f"<td>{_fmt(candidate.get('rank_score'), 4)}</td>"
            f"<td>{_esc(candidate.get('dataset_source') or '—')} / "
            f"{_esc(candidate.get('lib_spectrum_id') or '—')}</td>"
            f"<td>{_esc(diagnostic)}</td>"
            "</tr>"
        )
    names = (
        ", ".join(
            str(item.get("compound_name") or "")
            if isinstance(item, dict) else str(item)
            for item in isomer_set
        )
        if isomer_set else "None recorded"
    )
    rank_text = (
        "Experimental isomer-fragment rank override was active for this result."
        if rank_applied
        else "Candidate-specific fragment evidence is annotation-only and did not "
        "override the production ranking."
    )
    observed_product_ions = ", ".join(
        _fmt(item.get("query_mz"), 5)
        for item in detail.get("candidate_specific_product_ions") or []
        if item.get("observed") and item.get("query_mz") is not None
    )
    observed_neutral_losses = ", ".join(
        f"−{_fmt(item.get('query_loss_da'), 5)}"
        for item in detail.get("candidate_specific_neutral_losses") or []
        if item.get("observed") and item.get("query_loss_da") is not None
    )
    literature = (
        ann.get("isomer_literature_evidence")
        or (ann.get("isomer_group_evidence") or {}).get(
            "literature_evidence"
        )
        or {}
    )
    literature_sources = literature.get("sources") or []
    literature_items = "".join(
        "<li>"
        f'<a href="{_esc(source.get("url") or "")}">'
        f'{_esc(source.get("title") or source.get("source_id") or "Source")}</a>'
        f" — {_esc(source.get('evidence_scope') or '')}"
        f" (PMID {_esc(source.get('pmid') or '—')}; "
        f"DOI {_esc(source.get('doi') or '—')})"
        "</li>"
        for source in literature_sources[:5]
        if source.get("url")
    )
    parts = [
        '<section class="evidence"><h4>Isomer assessment</h4>',
        f"<p><b>Competitive same-formula isomers:</b> {'Yes' if ambiguous else 'No'}"
        f" &nbsp; <b>Candidate count:</b> {_esc(projection.get('candidate_count', 0))}"
        f" &nbsp; <b>Competing candidates after the lead:</b> "
        f"{_esc(projection.get('competitive_candidate_count', 0))}<br>"
        f"<b>Competitive set after the ranking lead:</b> {_esc(names)}<br>"
        f"<b>Exact-structure confidence:</b> "
        f"{_esc(ann.get('exact_structure_confidence') or ann.get('compound_confidence') or 'Unknown')}"
        f" &nbsp; <b>2-D-connectivity confidence:</b> "
        f"{_esc(two_dimensional_connectivity_confidence(match))}"
        f" &nbsp; <b>Isomer-group confidence:</b> "
        f"{_esc(ann.get('isomer_group_confidence') or 'Unknown')}</p>",
        f"<p><b>Candidate-specific fragment score:</b> "
        f"{_fmt(match.get('isomer_fragment_score'), 4)} &nbsp; "
        f"<b>Supported unique ions:</b> {_esc(detail.get('unique_supported', 0))}/"
        f"{_esc(detail.get('unique_total', 0))} &nbsp; "
        f"<b>Supported candidate-specific neutral losses:</b> "
        f"{_esc(detail.get('neutral_loss_supported', 0))}/"
        f"{_esc(detail.get('neutral_loss_total', 0))}<br>"
        f"<b>Metadata comparability:</b> "
        f"{_esc(detail.get('metadata_comparability') or 'not available')}<br>"
        f"<b>Observed candidate-specific product m/z:</b> "
        f"{_esc(observed_product_ions or 'none')} &nbsp; "
        f"<b>Observed candidate-specific neutral losses:</b> "
        f"{_esc(observed_neutral_losses or 'none')}</p>",
        f'<p class="method-note">{_esc(rank_text)} LC retention time under the same '
        "method, ion mobility/CCS, authentic standards, or additional orthogonal "
        "evidence may be required for positional and stereochemical isomers.</p>",
        (
            "<p><b>Curated literature/method leads:</b></p>"
            f"<ul>{literature_items}</ul>"
            '<p class="method-note">Literature metadata is not a measured '
            "reference spectrum or a same-method authenticated standard and "
            "does not affect Final, Top-1 or confidence.</p>"
            if literature_items else ""
        ),
    ]
    isopair = match.get("isopair_evidence") or {}
    if isopair:
        gate_result = ann.get(
            "isopair_high_precision_gate_result", "not_applicable"
        )
        parts.append(
            "<p><b>IsoPair-MS formula-conditional decision:</b> "
            f"{_esc(isopair.get('status') or 'not evaluated')} &nbsp; "
            f"<b>Ranking changed:</b> "
            f"{'yes' if isopair.get('ranking_changed') else 'no'} &nbsp; "
            f"<b>High precision gate:</b> {_esc(gate_result)}<br>"
            f"<b>Baseline connectivity margin:</b> "
            f"{_fmt(isopair.get('baseline_margin'), 4)} &nbsp; "
            f"<b>Fused winner score/margin:</b> "
            f"{_fmt(isopair.get('fused_winner_score'), 4)} / "
            f"{_fmt(isopair.get('fused_winner_margin'), 4)}</p>"
            '<p class="method-note">IsoPair-MS operates only within the same '
            "molecular formula and adduct and predicts the first InChIKey block "
            "(2-D connectivity). A near-tie rerank is never promoted to High by "
            "the model alone, and no general stereochemical claim is made.</p>"
        )
    if rows:
        parts.extend(
            [
                '<table class="mini"><thead><tr><th>Rank / role</th><th>Candidate</th>'
                "<th>Formula</th><th>InChIKey</th><th>Production cosine</th><th>Direct</th>"
                "<th>Modified-Hungarian</th><th>Neutral loss</th><th>Final</th>"
                "<th>Rank score</th><th>Library source / spectrum</th>"
                "<th>Candidate-specific evidence</th></tr></thead><tbody>",
                "".join(rows),
                "</tbody></table>",
            ]
        )
    parts.append("</section>")
    return "".join(parts)


def _relative_peaks(peaks: list) -> list[tuple[float, float]]:
    parsed = []
    for peak in peaks or []:
        try:
            mz, intensity = float(peak[0]), abs(float(peak[1]))
        except (TypeError, ValueError, IndexError):
            continue
        if mz > 0 and intensity > 0:
            parsed.append((mz, intensity))
    maximum = max((intensity for _, intensity in parsed), default=0.0)
    if maximum <= 0:
        return []
    return [(mz, 100.0 * intensity / maximum) for mz, intensity in parsed]


def _mirror_spectrum_svg(query: list, library: list, precursor_mz: Any) -> str:
    q = _relative_peaks(query)
    lib = _relative_peaks(library)
    if not q and not lib:
        return '<p class="muted">No peak arrays are available for visualization.</p>'
    maximum_mz = max(
        [mz for mz, _ in q + lib]
        + [value for value in [_finite_float(precursor_mz)] if value is not None]
        + [1.0]
    )
    width, height, left, right, mid = 920, 260, 38, 12, 128
    plot_width = width - left - right

    def x(mz: float) -> float:
        return left + plot_width * mz / maximum_mz

    bars = []
    for mz, intensity in q:
        bars.append(
            f'<line x1="{x(mz):.2f}" y1="{mid}" x2="{x(mz):.2f}" '
            f'y2="{mid - intensity:.2f}" stroke="#007A78" stroke-width="1.2"/>'
        )
    for mz, intensity in lib:
        bars.append(
            f'<line x1="{x(mz):.2f}" y1="{mid}" x2="{x(mz):.2f}" '
            f'y2="{mid + intensity:.2f}" stroke="#6F4E9B" stroke-width="1.2"/>'
        )
    return (
        f'<svg class="spectrum" viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Mirror plot of query and library MS/MS spectra">'
        f'<line x1="{left}" y1="{mid}" x2="{width-right}" y2="{mid}" stroke="#64748b"/>'
        f'{"".join(bars)}'
        f'<text x="42" y="18" fill="#007A78" font-size="12">Query (upward) · {len(q)} peaks</text>'
        f'<text x="42" y="248" fill="#6F4E9B" font-size="12">Library (downward) · {len(lib)} peaks</text>'
        f'<text x="{width-right-80}" y="{mid-5}" fill="#475569" font-size="10">'
        f"{maximum_mz:.1f} m/z</text></svg>"
    )


def _render_evograph(graph: Any) -> str:
    """Render candidate-conditioned formula-composition hypotheses in English."""

    if not isinstance(graph, dict):
        return ""
    context = graph.get("context") or {}
    audit = graph.get("graph_audit") or {}
    nodes = [node for node in (graph.get("nodes") or []) if isinstance(node, dict)]
    node_by_id = {str(node.get("node_id") or ""): node for node in nodes}
    fragment_rows = []
    for node in nodes:
        if node.get("node_type") != "fragment":
            continue
        candidates = node.get("formula_candidates") or []
        best = candidates[0] if candidates else {}
        status = {
            "heuristically_supported": "heuristically supported",
            "ambiguous": "ambiguous",
            "unassigned": "unassigned",
        }.get(
            str(node.get("incoming_edge_status") or "unassigned"),
            str(node.get("incoming_edge_status") or "unassigned"),
        )
        ambiguity = (
            f"{len(candidates)} compositions" if len(candidates) > 1
            else "one retained composition" if best
            else "no composition within constraints"
        )
        fragment_rows.append(
            "<tr>"
            f"<td>{_fmt(node.get('mz'), 5)}</td>"
            f"<td>{_esc(best.get('formula') or '—')}</td>"
            f"<td>{_fmt(best.get('mass_error_ppm'), 2)}</td>"
            f"<td>{_esc(status)}</td><td>{_esc(ambiguity)}</td></tr>"
        )
    edge_rows = []
    for edge in graph.get("edges") or []:
        if not isinstance(edge, dict):
            continue
        if not edge.get("selected") or edge.get("kind") != "formula_conserving":
            continue
        parent = node_by_id.get(str(edge.get("parent_id") or ""), {})
        child = node_by_id.get(str(edge.get("child_id") or ""), {})
        prior = edge.get("common_loss_prior") or {}
        edge_rows.append(
            "<tr>"
            f"<td>{_fmt(parent.get('mz'), 5)} → {_fmt(child.get('mz'), 5)}</td>"
            f"<td>{_esc(edge.get('parent_formula') or '—')} → "
            f"{_esc(edge.get('child_formula') or '—')}</td>"
            f"<td>−{_fmt(edge.get('observed_loss_da'), 5)}</td>"
            f"<td>{_esc(edge.get('loss_formula') or '—')}</td>"
            f"<td>{_fmt(edge.get('mass_error_ppm'), 2)}</td>"
            f"<td>{_esc(prior.get('name') or edge.get('selection_status') or '—')}</td></tr>"
        )
    available = graph.get("formula_graph_available") is True
    reasons = ", ".join(str(value) for value in context.get("reasons") or [])
    status = (
        f"Checked {int(audit.get('fragment_node_count') or 0)} observed product ions; "
        f"retained {int(audit.get('selected_edge_count') or 0)} graph edges."
        if available else
        "Formula graph unavailable: "
        + (reasons or str(graph.get("reason") or "formula/adduct/charge context incomplete"))
    )
    tables = ""
    if fragment_rows:
        tables += (
            "<h5>Observed product ions and candidate neutral compositions</h5>"
            '<table class="mini"><thead><tr><th>Observed m/z</th><th>Candidate composition</th>'
            '<th>Error (ppm)</th><th>Incoming-edge state</th><th>Composition ambiguity</th>'
            f"</tr></thead><tbody>{''.join(fragment_rows)}</tbody></table>"
        )
    if edge_rows:
        tables += (
            "<h5>Formula-conserving fragmentation edges</h5>"
            '<table class="mini"><thead><tr><th>Observed path</th><th>Parent → child formula</th>'
            '<th>Observed loss</th><th>Loss composition</th><th>Error (ppm)</th><th>Annotation</th>'
            f"</tr></thead><tbody>{''.join(edge_rows)}</tbody></table>"
        )
    return (
        '<section class="evidence"><h4>Observed-fragment composition graph (EvoGraph-F)</h4>'
        f"<p><b>{_esc(status)}</b><br>Selected-ion/isolation target m/z "
        f"{_fmt(graph.get('selected_ion_mz'), 6)}; analytical search mass m/z "
        f"{_fmt(graph.get('search_precursor_mz', graph.get('precursor_mz')), 6)}; "
        f"parent formula mass error {_fmt(context.get('mass_consistency_error_ppm'), 2)} ppm.</p>"
        f"{tables}"
        '<p class="method-note">These are candidate-conditioned elemental-composition hypotheses under an explicit adduct and charge-retention model. They are not confirmed product structures or atom-mapped bond-cleavage mechanisms; heuristic weights are not probabilities. EvoGraph-F contributes 0 to retrieval, ranking, confidence and High promotion. CFM/RDKit product structures, when available, are a separate optional output.</p>'
        "</section>"
    )


def _matched_fragment_rows(query: list, library: list, tolerance: float = 0.05) -> str:
    q = _relative_peaks(query)
    lib = _relative_peaks(library)
    matched = []
    used_library: set[int] = set()
    for qmz, qint in sorted(q, key=lambda item: item[1], reverse=True):
        options = [
            (abs(qmz - lmz), index, lmz, lint)
            for index, (lmz, lint) in enumerate(lib)
            if index not in used_library and abs(qmz - lmz) <= tolerance
        ]
        if not options:
            continue
        delta, index, lmz, lint = min(options)
        used_library.add(index)
        matched.append((qmz, lmz, delta, qint, lint))
        if len(matched) >= 30:
            break
    if not matched:
        return '<p class="muted">No peak pairs were available within ±0.05 Da.</p>'
    rows = "".join(
        "<tr>"
        f"<td>{qmz:.5f}</td><td>{lmz:.5f}</td><td>{delta:.5f}</td>"
        f"<td>{qint:.1f}</td><td>{lint:.1f}</td></tr>"
        for qmz, lmz, delta, qint, lint in matched
    )
    return (
        '<table class="mini"><thead><tr><th>Query m/z</th><th>Library m/z</th>'
        "<th>|Δm/z|</th><th>Query relative intensity</th>"
        f"<th>Library relative intensity</th></tr></thead><tbody>{rows}</tbody></table>"
    )


def _detail_panel(index: int, row: dict, source_name: str = "") -> str:
    match = _match_of(row) or {}
    display_name = (
        _standard_compound_name(match.get("compound_name", ""), match)
        or "Unnamed candidate"
    )
    query = row.get("query_peaks") or row.get("peaks") or []
    library = row.get("lib_peaks") or match.get("lib_peaks") or []
    structure = svg_for_smiles(str(match.get("smiles") or ""), size=(250, 180))
    structure_html = structure if structure else '<p class="muted">Structure unavailable.</p>'
    expert = match.get("expert_details") or {}
    pathway_align = (expert.get("pathway_align") or {}).get("score")
    expert_fusion = expert.get("score_fusion") or {}
    expert_source = (
        expert_fusion.get("selected_source")
        or match.get("expert_score_source")
        or "none"
    )
    expert_source_label = {
        "base_8d": "Base 8 dimensions",
        "v3_kb": "V3 knowledge base",
        "compound_pathway": "Candidate pathway",
        "none": "No expert evidence",
    }.get(str(expert_source), str(expert_source).replace("_", " ").title())
    expert_counts = _expert_dimension_audit(match)
    shadow = match.get("contrastive_expert_shadow") or {}
    isolation = row.get("precursor_isolation_evidence") or {}
    locator = raw_spectrum_locator(row, source_name)
    locator_text = raw_spectrum_locator_text(
        row, source_name, language="en"
    )
    rt_text = (
        f"{locator['retention_time_min']:.3f} min "
        f"({locator['retention_time_clock']})"
        if locator["retention_time_min"] is not None
        else "not recorded"
    )
    selected_mz = _finite_float(row.get("precursor_mz"))
    search_mz = _finite_float(row.get("search_precursor_mz"))
    feature_mass_used = bool(row.get("lcms_feature_used_for_search"))
    mass_summary = f"selected-ion m/z {_fmt(selected_mz, 5)}"
    if feature_mass_used and search_mz is not None:
        mass_summary += f" / measured-MS1 search m/z {_fmt(search_mz, 5)}"
    return (
        f'<details class="compound-detail" id="compound-{index}">'
        f"<summary>#{index} {_esc(display_name)} "
        f"— RT {_esc(rt_text)} — MS2 {_esc(locator['spectrum_id'] or 'not recorded')} "
        f"— {_esc(mass_summary)} — final {_fmt(match.get('final_score'), 4)}</summary>"
        '<div class="detail-grid"><section><h4>Identity and provenance</h4>'
        f"<p><b>Name:</b> {_esc(display_name)}<br>"
        f"<b>Formula:</b> {_esc(match.get('formula') or '—')}<br>"
        f"<b>InChIKey:</b> <span class=\"mono\">{_esc(match.get('inchikey') or '—')}</span><br>"
        f"<b>SMILES:</b> <span class=\"mono wrap\">{_esc(match.get('smiles') or '—')}</span><br>"
        f"<b>Classification:</b> {_esc(_classification_en(match))}<br>"
        f"<b>Classification evidence:</b> "
        f"{_esc(classification_evidence_status(match)['status'])}; source "
        f"{_esc(str(classification_evidence_status(match)['source']).replace('_', ' '))}; "
        f"reason {_esc(classification_evidence_status(match)['reason'])}; "
        f"descriptive score {_fmt(match.get('display_classification_confidence', match.get('classification_confidence')), 3)}<br>"
        f"<b>Library source:</b> {_esc(match.get('dataset_source') or 'Not recorded')} "
        f"({_esc(match.get('source_tier') or 'unknown')})</p></section>"
        f'<section class="structure"><h4>Candidate structure</h4>{structure_html}</section></div>'
        '<section class="evidence"><h4>Identification evidence</h4>'
        f"<p>{_confidence_badge(match)}</p>"
        '<div class="metric-grid">'
        f"<div><b>Cosine</b><span>{_fmt(match.get('cosine'), 4)}</span></div>"
        f"<div><b>Entropy</b><span>{_fmt(match.get('entropy'), 4)}</span></div>"
        f"<div><b>Reverse coverage</b><span>{_fmt(match.get('reverse_match'), 4)}</span></div>"
        f"<div><b>Mass error</b><span>{_fmt(match.get('mass_ppm'), 2)} ppm</span></div>"
        f"<div><b>Final score</b><span>{_fmt(match.get('final_score'), 4)}</span></div>"
        f"<div><b>Rule dimensions with support</b><span>{expert_counts['supporting']}/{expert_counts['total'] or 8}</span></div>"
        f"<div><b>Evaluated with no observation</b><span>{expert_counts['no_evidence']}</span></div>"
        f"<div><b>Not applicable / unavailable</b><span>{expert_counts['not_applicable']} / {expert_counts['unavailable']}</span></div>"
        f"<div><b>Rule-evidence record source</b><span>{_esc(expert_source_label)}</span></div>"
        f"<div><b>Pathway alignment</b><span>{_fmt(pathway_align, 4)}</span></div>"
        f"<div><b>Isolation purity</b><span>{_fmt(isolation.get('purity'), 3)}</span></div>"
        "</div>"
        f"<p class=\"method-note\">Rule evidence is an applicability-aware interpretation, not an identity probability, and contributes 0 to identity ranking. Contrastive expert shadow: {_esc(shadow.get('status') or 'not available')}; "
        f"{len(shadow.get('matched_ions') or [])} ion/s and "
        f"{len(shadow.get('matched_losses') or [])} loss/es; score contribution = 0.</p>"
        "</section>"
        '<section class="evidence"><h4>Query–library mirror spectrum</h4>'
        f'<div class="raw-locator"><b>Raw MS2 locator</b><br>{_esc(locator_text)}<br>'
        '<span>Search by spectrum/native ID first. If the viewer cannot search that identifier, combine retention time and precursor m/z, then verify the query fragment peaks.</span></div>'
        f"{_mirror_spectrum_svg(query, library, row.get('search_precursor_mz', row.get('precursor_mz')))}"
        f"<h5>Matched product ions</h5>{_matched_fragment_rows(query, library)}</section>"
        f"{_render_evograph(row.get('evograph_fragmentation'))}"
        f"{_render_pathway(match, query)}"
        f"{_render_isotope(row, match)}"
        f"{_render_adduct(match)}"
        f"{_render_ce_profile(match)}"
        f"{_render_isomer(match) if isomer_candidate_projection(match).get('status') != 'not_applicable' else ''}"
        "</details>"
    )


_CSS = """
:root{--blue:#123c69;--cyan:#1f7a8c;--ink:#172033;--muted:#667085;--line:#d9e2ec;--bg:#f4f7fa}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.container{max-width:1600px;margin:auto;padding:24px}.hero{background:linear-gradient(135deg,#123c69,#1f7a8c);color:#fff;padding:26px 30px;border-radius:14px;display:flex;justify-content:space-between;gap:20px;align-items:center}
.hero h1{margin:0 0 5px;font-size:26px}.hero p{margin:2px 0;opacity:.92}.hero button{border:1px solid #ffffff88;background:#ffffff18;color:white;padding:9px 14px;border-radius:8px;cursor:pointer}
.cards{display:grid;grid-template-columns:repeat(5,minmax(130px,1fr));gap:10px;margin:14px 0}.card,.notice,.table-box,.compound-detail{background:#fff;border:1px solid var(--line);border-radius:10px;box-shadow:0 2px 5px #102a4310}
.card{padding:13px;text-align:center}.card strong{display:block;font-size:22px;color:var(--blue)}.card span{color:var(--muted);font-size:12px;text-transform:uppercase}
.notice{padding:14px 18px;margin:12px 0}.notice h2{font-size:17px;margin:0 0 7px}.notice p{margin:5px 0}.method-note{color:#475467;background:#f8fafc;border-left:3px solid #94a3b8;padding:8px 10px}
.toolbar{display:flex;gap:8px;flex-wrap:wrap;margin:12px 0}.toolbar input,.toolbar select{padding:8px;border:1px solid var(--line);border-radius:7px;background:#fff}
.table-box{overflow:auto}table{width:100%;border-collapse:collapse}th{background:var(--blue);color:#fff;text-align:left;font-size:12px;padding:8px;white-space:nowrap}td{padding:8px;border-bottom:1px solid #edf2f7;vertical-align:top}tbody tr:hover{background:#edf7fa}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.wrap{word-break:break-all}.badge{display:inline-block;border-radius:999px;padding:2px 8px;font-size:11px;font-weight:700}.badge.high{background:#d1fadf;color:#05603a}.badge.medium{background:#fef0c7;color:#93370d}.badge.low{background:#fee4e2;color:#912018}.badge.unknown{background:#eaecf0;color:#344054}.badge.outline{background:#fff;border:1px solid #98a2b3;color:#475467}.scope-note{display:block;margin-top:3px;color:#667085;font-size:10px}
.compound-detail{margin:12px 0;overflow:hidden}.compound-detail>summary{padding:12px 15px;font-weight:700;cursor:pointer;background:#f8fafc}.compound-detail>section,.compound-detail>.detail-grid{margin:14px}
.detail-grid{display:grid;grid-template-columns:2fr 1fr;gap:14px}.detail-grid section,.evidence{border:1px solid var(--line);border-radius:9px;padding:12px}.detail-grid h4,.evidence h4{margin:0 0 8px;color:var(--blue)}.evidence h5{margin:14px 0 6px;color:#344054}
.structure svg{max-height:180px}.metric-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:8px}.metric-grid>div{background:#f8fafc;padding:8px;border-radius:7px}.metric-grid b,.metric-grid span{display:block}.metric-grid span{font-size:17px;color:var(--cyan)}
.mini{font-size:12px}.mini th{background:#345d7e}.compact{margin:5px 0}.muted{color:var(--muted)}.spectrum{width:100%;max-height:300px;background:#fff}
.raw-locator{margin:0 0 10px;padding:9px 11px;border:1px solid #9cc9c6;border-left:4px solid #176b68;border-radius:6px;background:#f0f8f7;color:#214c50;font-size:12px;line-height:1.55}.raw-locator span{color:#52676d}
footer{text-align:center;color:var(--muted);padding:25px}.hidden-row{display:none}
@media(max-width:900px){.cards{grid-template-columns:repeat(2,1fr)}.detail-grid{grid-template-columns:1fr}.metric-grid{grid-template-columns:repeat(2,1fr)}.container{padding:10px}.hero{align-items:flex-start;flex-direction:column}}
@media print{body{background:#fff}.container{max-width:none;padding:0}.hero{print-color-adjust:exact}.toolbar,.hero button{display:none}.compound-detail{break-inside:avoid}.compound-detail>summary{list-style:none}.compound-detail:not([open])>*:not(summary){display:block}}
"""


_JS = """
function filterRows(){
  const q=(document.getElementById('q').value||'').toLowerCase();
  const c=document.getElementById('classFilter').value;
  let shown=0;
  document.querySelectorAll('#resultsBody tr').forEach(row=>{
    const ok=(!q||row.dataset.search.includes(q))&&(!c||row.dataset.class===c);
    row.classList.toggle('hidden-row',!ok); if(ok) shown++;
  });
  document.getElementById('shown').textContent=shown;
}
function sortRows(){
  const body=document.getElementById('resultsBody'), mode=document.getElementById('sort').value;
  const rows=Array.from(body.querySelectorAll('tr'));
  rows.sort((a,b)=>mode==='name'?a.dataset.name.localeCompare(b.dataset.name):
    parseFloat(b.dataset[mode]||0)-parseFloat(a.dataset[mode]||0));
  rows.forEach(row=>body.appendChild(row));
}
"""


def generate_english_html(
    results: list[dict],
    source_name: str = "Analysis",
    *,
    detailed: bool = False,
    total_spectra: int | None = None,
    matched_spectra: int | None = None,
    high_supported_compounds: int | None = None,
    high_confidence_scope_groups: int | None = None,
    high_confidence_exact_structures: int | None = None,
    high_confidence_isomer_groups: int | None = None,
) -> str:
    """Return a self-contained, English-only system-authored HTML report."""
    rows = list(results or [])
    matched, unmatched, duplicate_count = _deduplicate(rows)
    total_scan_count = (
        max(0, int(total_spectra))
        if total_spectra is not None
        else len(rows)
    )
    matched_scan_count = (
        max(0, min(total_scan_count, int(matched_spectra)))
        if matched_spectra is not None
        else sum(1 for row in rows if _match_of(row))
    )
    unmatched_scan_count = max(0, total_scan_count - matched_scan_count)
    unique_matched_count = len(matched)
    scan_match_rate = (
        100.0 * matched_scan_count / total_scan_count
        if total_scan_count else 0.0
    )
    high_supported = (
        max(0, int(high_supported_compounds))
        if high_supported_compounds is not None
        else len({
            _identity_key(row)
            for row in rows
            if _match_of(row)
            and reportable_identity_confidence(_match_of(row) or {}) == "High"
        })
    )
    high_exact = max(0, int(
        high_confidence_exact_structures
        if high_confidence_exact_structures is not None else high_supported
    ))
    high_isomer = max(0, int(
        high_confidence_isomer_groups
        if high_confidence_isomer_groups is not None else 0
    ))
    high_scope = max(0, int(
        high_confidence_scope_groups
        if high_confidence_scope_groups is not None
        else high_exact + high_isomer
    ))
    classes = sorted({_classification_en(_match_of(row) or {}) for row in matched})
    table_rows = []
    for index, row in enumerate(matched, 1):
        match = _match_of(row) or {}
        locator = raw_spectrum_locator(row, source_name)
        rt_text = (
            f"{locator['retention_time_min']:.3f} min"
            if locator["retention_time_min"] is not None
            else "not recorded"
        )
        cls = _classification_en(match)
        isomer_projection = isomer_candidate_projection(match)
        isomer_names = [
            str(candidate.get("compound_name") or "")
            for candidate in isomer_projection.get("competing_candidates") or []
            if str(candidate.get("compound_name") or "").strip()
        ]
        name = (
            _standard_compound_name(match.get("compound_name", ""), match)
            or "Unnamed candidate"
        )
        search_text = " ".join(
            [
                name,
                str(match.get("formula") or ""),
                str(match.get("inchikey") or ""),
                cls,
                " ".join(isomer_names),
            ]
        ).lower()
        name_cell = (
            f'<a href="#compound-{index}">{_esc(name)}</a>'
            if detailed
            else _esc(name)
        )
        competition_cell = (
            f'<b>{int(isomer_projection.get("candidate_count") or 0)} candidates</b><br>'
            f'Ranking lead: {_esc((isomer_projection.get("lead_candidate") or {}).get("compound_name") or name)}<br>'
            f'Competitors: {_esc("; ".join(isomer_names) or "details unavailable")}'
            if isomer_projection.get("status") != "not_applicable"
            else "—"
        )
        table_rows.append(
            f'<tr data-search="{_esc(search_text)}" data-class="{_esc(cls)}" '
            f'data-name="{_esc(name.lower())}" data-final="{_fmt(match.get("final_score"), 8, "0")}" '
            f'data-cosine="{_fmt(match.get("cosine"), 8, "0")}">'
            f"<td>{index}</td><td class=\"mono\"><b>selected m/z {_fmt(row.get('precursor_mz'), 5)}</b><br>"
            f"{'search m/z ' + _fmt(row.get('search_precursor_mz'), 5) + '<br>' if row.get('lcms_feature_used_for_search') else ''}"
            f"RT {_esc(rt_text)}<br>MS2 {_esc(locator['spectrum_id'] or 'not recorded')}</td>"
            f"<td>{name_cell}</td>"
            f"<td>{competition_cell}</td>"
            f"<td>{_esc(match.get('formula') or '—')}</td><td>{_esc(cls)}</td>"
            f"<td>{_confidence_badge(match)}</td>"
            f"<td>{final_claim_badge_html(match, language='en')}</td>"
            f"<td>{_fmt(match.get('cosine'), 4)}</td>"
            f"<td>{_fmt(match.get('entropy'), 4)}</td><td>{_fmt(match.get('reverse_match'), 4)}</td>"
            f"<td>{_fmt(match.get('mass_ppm'), 2)}</td><td>{_fmt(match.get('final_score'), 4)}</td>"
            f"<td>{_esc(match.get('adduct') or '—')}</td></tr>"
        )
    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    class_options = "".join(
        f'<option value="{_esc(value)}">{_esc(value)}</option>' for value in classes
    )
    details = (
        '<h2 id="evidence-details">Compound evidence details</h2>'
        + "".join(
            _detail_panel(index, row, source_name)
            for index, row in enumerate(matched, 1)
        )
        if detailed
        else ""
    )
    html_document = f"""<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_esc(source_name)} — LC–MS/MS Identification Report</title><style>{_CSS}</style></head>
<body><main class="container">
<header class="hero"><div><h1>LC–MS/MS Identification Report</h1>
<p><b>Source:</b> {_esc(source_name)}</p><p>Generated: {_esc(generated)}</p></div>
<button onclick="window.print()">Print / Save as PDF</button></header>
<section class="cards">
<div class="card"><strong>{total_scan_count}</strong><span>Total MS/MS spectra</span></div>
<div class="card"><strong>{matched_scan_count}</strong><span>MS/MS spectra with a library candidate ({scan_match_rate:.1f}%)</span></div>
<div class="card"><strong>{unique_matched_count}</strong><span>Unique matched compounds</span></div>
<div class="card"><strong>{high_scope}</strong><span>High evidence groups<br><small>single structures {high_exact} · isomer groups {high_isomer}</small></span></div>
<div class="card"><strong>{unmatched_scan_count}</strong><span>Unmatched MS/MS spectra</span></div>
</section>
<section class="notice"><h2>Interpretation and scope</h2>
<p>Scores and candidate order are copied from the saved analysis result; report generation does not rescore spectra.</p>
<p><b>Compound name</b> is the best library candidate retained by the saved analysis. The adjacent confidence label communicates the strength of that assignment without replacing the name with internal audit terminology.</p>
<p><b>Identification level:</b> a library-spectrum match is normally a putative Level 2 annotation. Level 1 confirmation requires an authentic standard measured under the same analytical method, including compatible retention time.</p>
<p><b>Fragmentation pathways:</b> stored pathways are empirical recurring ions and neutral-loss sequences mined from reference spectra. They are not automatically atom-mapped reaction mechanisms.</p>
<p><b>Adducts and isotopes:</b> a specific adduct is used only when supported by input metadata. Isotope envelopes come from adjacent MS1 scans, never from MS2 fragments.</p>
<p><b>Raw-data review:</b> each result keeps the observed query MS2 retention time, precursor m/z, and spectrum/native ID together. Search the source mzML by spectrum/native ID first; use RT plus precursor m/z as the fallback locator.</p>
<p><b>Classification:</b> “Unclassified” means the database or calibrated spectrum classifier did not provide a sufficiently specific class; it is not a chemical class.</p>
{f'<p>{duplicate_count} lower-scoring duplicate library assignments were collapsed by compound identity.</p>' if duplicate_count else ''}
</section>
<div class="toolbar"><input id="q" oninput="filterRows()" placeholder="Filter by name, formula, InChIKey, or class">
<select id="classFilter" onchange="filterRows()"><option value="">All classes</option>{class_options}</select>
<select id="sort" onchange="sortRows()"><option value="final">Final score, descending</option>
<option value="cosine">Cosine, descending</option><option value="name">Name, A–Z</option></select>
<span class="muted">Showing <b id="shown">{unique_matched_count}</b> compounds</span></div>
<section class="table-box"><table><thead><tr><th>#</th><th>Raw MS2 locator</th><th>Compound</th><th>Isomer competition</th>
<th>Formula</th><th>Classification</th><th>Evidence grade</th><th>Identification conclusion</th><th>Cosine</th><th>Entropy</th>
<th>Reverse coverage</th><th>Mass error (ppm)</th><th>Final score</th><th>Library adduct</th>
</tr></thead><tbody id="resultsBody">{''.join(table_rows)}</tbody></table></section>
{render_isomer_candidate_competition_html(matched, language="en")}
{details}
<footer>Generated by MSsystem V12. Automated annotations must be interpreted with the underlying evidence and analytical context.</footer>
</main><script>{_JS}</script></body></html>"""
    return html_document


__all__ = ["class_name_en", "generate_english_html"]
