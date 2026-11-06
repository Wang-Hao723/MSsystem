"""Read-only replay of the peak view used for a saved candidate comparison."""
from __future__ import annotations

import hashlib
import json
import math

import numpy as np
from matchms.similarity.spectrum_similarity_functions import (
    collect_peak_pairs,
    score_best_matches,
)


def _peak_array(peaks):
    array = np.asarray(peaks, dtype=np.float64)
    if array.size == 0:
        return np.empty((0, 2), dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("scoring peaks must have shape (N, 2)")
    if not np.isfinite(array).all() or (array[:, 0] <= 0).any() or (array[:, 1] < 0).any():
        raise ValueError("invalid scoring peak values")
    return array[np.argsort(array[:, 0], kind="mergesort")].copy()


def _fingerprint(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()).hexdigest()


def capture_spectral_evidence(candidate, *, query_precursor_mz, mz_tol=0.1):
    """Capture an existing comparison view without changing any candidate score.

    Pair indices address the sorted peak arrays saved in this object, not the
    raw acquisition array. Pairs are replayed with the same matchms primitives
    as ``cosine_score``; they do not claim a fragment-to-atom assignment.
    """
    q = _peak_array(candidate.get("_score_query_peaks") or [])
    l = _peak_array(candidate.get("_score_library_peaks") or [])
    if not len(q) or not len(l):
        return {"schema": "mssystem.spectral_evidence.v1", "status": "unavailable",
                "reason": "scoring_peak_view_missing", "ranking_influence": False}
    q_pmz = float(query_precursor_mz)
    l_pmz = float(candidate.get("lib_search_precursor_mz") or candidate.get("lib_precursor_mz"))
    tolerance = float(mz_tol)
    if not all(math.isfinite(x) and x > 0 for x in (q_pmz, l_pmz, tolerance)):
        raise ValueError("invalid precursor or tolerance")
    q_saved, l_saved = q.tolist(), l.tolist()
    q_max, l_max = float(q[:, 1].max()), float(l[:, 1].max())
    if q_max > 0:
        q[:, 1] /= q_max
    if l_max > 0:
        l[:, 1] /= l_max
    shift = q_pmz - l_pmz
    pair_sets, mode_sets = [], []
    for mode, offset in (("direct", 0.0), ("precursor_shifted", shift)):
        if mode == "precursor_shifted" and abs(shift) <= tolerance:
            continue
        pairs = collect_peak_pairs(q, l, tolerance, shift=offset)
        if pairs is not None and len(pairs):
            pair_sets.append(pairs)
            mode_sets.extend([mode] * len(pairs))
    recorded, replayed = [], 0.0
    if pair_sets and q_max > 0 and l_max > 0:
        pairs = np.concatenate(pair_sets, axis=0)
        order = np.argsort(pairs[:, 2], kind="mergesort")[::-1]
        ordered = pairs[order]
        replayed = float(score_best_matches(ordered, q, l, 0.0, 1.0)[0])
        used_q, used_l = set(), set()
        for position in order:
            qi, li = int(pairs[position, 0]), int(pairs[position, 1])
            if qi in used_q or li in used_l:
                continue
            used_q.add(qi)
            used_l.add(li)
            mode = mode_sets[position]
            qmz, lmz = q_saved[qi][0], l_saved[li][0]
            recorded.append({
                "query_peak_index": qi, "library_peak_index": li,
                "query_mz": qmz, "lib_mz": lmz,
                "match_mode": mode,
                "residual_da": qmz - lmz - (shift if mode == "precursor_shifted" else 0.0),
            })
    if not math.isfinite(replayed):
        raise ValueError("nonfinite cosine replay")
    comparison = candidate.get("raw_cosine", candidate.get("cosine"))
    comparison = float(comparison) if comparison is not None else None
    if comparison is not None and not math.isfinite(comparison):
        comparison = None
    output = {
        "schema": "mssystem.spectral_evidence.v1", "status": "captured",
        "score_scope": "raw_spectrum_candidate_comparison",
        "pair_origin": "deterministic_replay_of_saved_scoring_views",
        "index_basis": "zero_based_saved_scoring_peak_arrays",
        "algorithm": "ModifiedCosineGreedy", "mz_power": 0.0, "intensity_power": 1.0,
        "tolerance_da": tolerance, "query_precursor_mz": q_pmz,
        "library_precursor_mz": l_pmz, "precursor_shift_da": shift,
        "query_view": str(candidate.get("search_pass") or "primary_top15"),
        "query_peaks": q_saved, "library_peaks": l_saved,
        "query_peak_sha256": _fingerprint(q_saved), "library_peak_sha256": _fingerprint(l_saved),
        "reference": {"spectrum_id": candidate.get("lib_spectrum_id"),
                      "source_accession": candidate.get("lib_source_accession"),
                      "dataset_source": candidate.get("dataset_source")},
        "matched_pairs": recorded, "replayed_modified_cosine": replayed,
        "comparison_cosine": comparison,
        "replay_matches_comparison": comparison is not None and math.isclose(
            comparison, replayed, rel_tol=1e-12, abs_tol=1e-12),
        "ranking_influence": False, "fragment_atom_assignment": False,
    }
    output["view_sha256"] = _fingerprint(output)
    return output


def attach_spectral_evidence(best, sources, *, query_precursor_mz):
    """Attach only after selection; public candidate order and scores stay intact."""
    from matching.library_quality import compound_identity_key

    def key(candidate):
        return (compound_identity_key(candidate), str(candidate.get("lib_spectrum_id")),
                str(candidate.get("dataset_source")), str(candidate.get("lib_source_accession")))

    by_reference = {key(source): source for source in sources}
    by_reference[key(best)] = best
    captured = {}
    for candidate in [best, *(best.get("top_candidates") or [])]:
        reference_key = key(candidate)
        if reference_key in captured:
            candidate["spectral_evidence"] = captured[reference_key]
            continue
        source = by_reference.get(reference_key)
        try:
            evidence = capture_spectral_evidence(source or {}, query_precursor_mz=query_precursor_mz)
        except (TypeError, ValueError, IndexError) as exc:
            evidence = {"schema": "mssystem.spectral_evidence.v1", "status": "unavailable",
                        "reason": str(exc), "ranking_influence": False}
        candidate["spectral_evidence"] = evidence
        captured[reference_key] = evidence
