"""Rescore every frozen candidate against every saved query spectrum.

The candidate pool is input, not reselected per scan. No cell is discarded by
Top-K truncation. Missing references/metadata remain distinct from measured
zero similarity. This computational building block is not a production ranker.
"""
from __future__ import annotations

import math
import re
from collections import Counter

from matching.conditioned_multilibrary import reference_condition, reference_content_fingerprint
from matching.evirank_retrieval import _polarity
from matching.library_record_store import semantic_sha256
from matching.scores import batch_modified_cosine_scores, entropy_score, reverse_match_ratio
from utils.spectrum import process_peaks, resolve_peak_processing_config

ALGORITHM_ID = "mssystem-v16-complete-candidate-query-panel-1.0.0"
QUERY_FIELDS = ("spectrum_id", "collision_energy", "collision_energy_unit", "collision_energy_mode",
    "instrument_type", "instrument", "activation_method", "instrument_domain_id", "collision_cell_region")


def project_query(row):
    result = {key: row.get(key) for key in QUERY_FIELDS}
    result.update({"precursor_mz": row.get("search_precursor_mz") or row.get("precursor_mz"),
        "polarity": _polarity(row.get("polarity") or row.get("ion_mode")),
        "charge": row.get("charge") or row.get("precursor_charge"),
        "adduct": row.get("query_adduct") or row.get("adduct"),
        "peaks": row.get("raw_peaks") or row.get("peaks") or row.get("query_peaks") or [],
        "raw_peaks_complete": row.get("raw_peaks_complete")})
    try:
        precursor = float(result["precursor_mz"])
    except (TypeError, ValueError):
        raise ValueError("Positive measured query precursor is required") from None
    if not math.isfinite(precursor) or precursor <= 0 or not result["polarity"]:
        raise ValueError("Finite positive precursor and explicit query polarity are required")
    result["precursor_mz"] = precursor
    return result


def _metrics(query_peaks, reference, processed, score, query, tolerance):
    return {"modified_cosine": round(score, 8),
        "spectral_entropy": round(float(entropy_score(query_peaks, processed, mz_tol=tolerance)), 8),
        "reverse_match": round(float(reverse_match_ratio(query_peaks, processed, mz_tol=tolerance,
            query_pmz=query["precursor_mz"], lib_pmz=reference["precursor_mz"])), 8),
        "reference_spectrum_id": str(reference.get("spectrum_id") or ""),
        "reference_content_fingerprint": reference_content_fingerprint(reference.get("peaks") or [])}


def rescore_complete_panel(source_rows, candidate_connectivities, retriever, *,
                           include_comparisons=False, baseline_connectivity=None):
    ids = sorted(set(str(value).split(":")[-1].split("-")[0].upper() for value in candidate_connectivities))
    if not ids or len(ids) > 40 or any(re.fullmatch(r"[A-Z]{14}", value) is None for value in ids):
        raise ValueError("A frozen panel must contain 1 to 40 distinct 2D connectivity keys")
    queries = [project_query(row) for row in source_rows]
    if not queries or len(queries) > 200:
        raise ValueError("A bounded panel requires 1 to 200 measured query spectra")
    if len({query["spectrum_id"] for query in queries}) != len(queries) or any(not q["spectrum_id"] for q in queries):
        raise ValueError("Unique explicit query spectrum IDs are required")
    query_fingerprints = {reference_content_fingerprint(query["peaks"]) for query in queries}
    views = retriever.policy["spectral_views"]
    if any(view["similarity"] != "matchms_modified_cosine_greedy" for view in views):
        raise ValueError("This panel version implements only explicit modified-cosine spectral views")
    minimum_peaks = int(retriever.policy["candidate_pool"]["minimum_library_fragment_peaks"])
    cells, window_audits, packets = [], [], []
    formulae = {identity: set() for identity in ids}
    for query in queries:
        records, window_audit = retriever._fetch_window(query["precursor_mz"], query["polarity"])
        window_audits.append({"query_id": query["spectrum_id"], **window_audit})
        records = [record for record in records if record["connectivity_identity"] in ids
                   and reference_content_fingerprint(record.get("peaks") or []) not in query_fingerprints]
        for record in records:
            formula = str(record.get("molecular_formula") or record.get("formula") or "").strip()
            if formula:
                formulae[record["connectivity_identity"]].add(formula)
        if include_comparisons:
            records.sort(key=lambda record: (record["connectivity_identity"],
                str(reference_content_fingerprint(record.get("peaks") or [])),
                str(record.get("spectrum_id") or "")))
        for view in views:
            config = resolve_peak_processing_config(view["peak_processing_preset"])
            tolerance = float(view["fragment_tolerance_da"])
            query_peaks, _ = process_peaks(query["peaks"], query["precursor_mz"], config=config,
                                           precursor_charge=query.get("charge"))
            ready, peak_lists, seen = [], [], set()
            for record in records:
                condition = reference_condition(query, record, near_delta_by_unit={"eV": 5.0, "NCE": 5.0})
                fingerprint = reference_content_fingerprint(record.get("peaks") or [])
                key = (record["connectivity_identity"], fingerprint,
                       bool(condition["eligible_for_conditioned_consensus"]))
                if include_comparisons:
                    # Keep acquisition strata needed by the paired comparison.
                    # Mirrors are collapsed before reference medians or votes.
                    key += (semantic_sha256({field: record.get(field) for field in (
                        *QUERY_FIELDS[1:], "dataset_source", "adduct", "charge", "ion_mode")}),)
                if key in seen:
                    continue
                seen.add(key)
                peaks, _ = process_peaks(record.get("peaks") or [], record["precursor_mz"],
                                         config=config, precursor_charge=record.get("charge"))
                if len(peaks) >= minimum_peaks:
                    ready.append((record, condition))
                    peak_lists.append(peaks)
            valid_query = len(query_peaks) >= minimum_peaks
            scores = batch_modified_cosine_scores(query_peaks, peak_lists, query["precursor_mz"],
                [float(record["precursor_mz"]) for record, _ in ready], mz_tol=tolerance) if ready and valid_query else []
            by_candidate = {identity: [] for identity in ids}
            for (record, condition), peaks, score in zip(ready, peak_lists, scores):
                if not math.isfinite(float(score)):
                    raise ValueError("Non-finite observed similarity")
                by_candidate[record["connectivity_identity"]].append((float(score), record, peaks, condition))
            if include_comparisons:
                packets.append({"query": query, "view": view["name"], "tolerance_da": tolerance,
                    "query_peaks": query_peaks,
                    "references": [{"record": record, "peaks": peaks, "score": score,
                        "condition": condition} for values in by_candidate.values()
                        for score, record, peaks, condition in values]})
            for identity in ids:
                eligible = by_candidate[identity]
                order = lambda value: (-value[0], str(value[1].get("spectrum_id") or ""))
                best = min(eligible, key=order) if eligible else None
                same_condition = [value for value in eligible if value[3]["eligible_for_conditioned_consensus"]]
                conditioned = min(same_condition, key=order) if same_condition else None
                project = lambda value: _metrics(query_peaks, value[1], value[2], value[0], query, tolerance) if value else None
                cells.append({"candidate_connectivity": identity, "query_id": query["spectrum_id"],
                    "spectral_view": view["name"], "query_condition": {key: query.get(key) for key in QUERY_FIELDS},
                    "status": "measured" if best else ("query_low_information" if not valid_query
                        else "no_reference_in_registered_window_after_quality_and_mirror_filters"),
                    "all_condition_reference_metrics": project(best),
                    "condition_matched_reference_metrics": project(conditioned),
                    "reference_profile_count": len(eligible), "condition_matched_profile_count": len(same_condition),
                    "same_condition_evidence_available": conditioned is not None,
                    "truncated_topk_controls_observation": False,
                    "raw_query_complete_declared": query.get("raw_peaks_complete") is True})
    expected = len(ids) * len(queries) * len(views)
    if len(cells) != expected:
        raise AssertionError("Complete candidate-by-query-by-view panel has missing cells")
    result = {"algorithm_id": ALGORITHM_ID, "candidate_connectivities": ids,
        "query_projection_sha256": semantic_sha256(queries), "library_snapshot_id": retriever.library_snapshot_id,
        "query_count": len(queries), "candidate_count": len(ids), "view_count": len(views),
        "expected_cells": expected, "cells": cells, "window_audits": window_audits,
        "cell_status_counts": dict(Counter(cell["status"] for cell in cells)),
        "condition_matched_cells": sum(cell["same_condition_evidence_available"] for cell in cells),
        "measured_zero_modified_cosine_cells": sum(
            cell["all_condition_reference_metrics"] is not None
            and cell["all_condition_reference_metrics"]["modified_cosine"] == 0 for cell in cells),
        "runtime_effect": False, "new_winner_selected": False, "score_is_probability": False,
        "query_content_mirrors_excluded": True, "reference_duplicates_create_extra_votes": False,
        "cross_instrument_or_unit_conversion_applied": False}
    result["near_collision_energy_delta_by_unit"] = {"eV": 5.0, "NCE": 5.0}
    if include_comparisons:
        from pipeline.v16_candidate_comparison import build_candidate_comparison
        result["comparison_evidence"] = build_candidate_comparison(
            ids, packets, candidate_formulae=formulae, baseline_connectivity=baseline_connectivity)
        result["algorithm_id"] = "mssystem-v16-complete-candidate-query-panel-1.1.0"
    result["semantic_sha256"] = semantic_sha256(result)
    return result
