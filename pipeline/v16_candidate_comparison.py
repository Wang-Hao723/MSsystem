"""Paired, positive-ion evidence for a frozen complete candidate panel.

Shared product ions contribute zero to the fragment contrast. An unobserved
query ion never contributes a penalty. CE trajectories use only ions actually
observed at both energies, with references from the same acquisition/provenance
stratum. These are review scores, not calibrated probabilities or identities.
"""
from __future__ import annotations

import itertools
import math
import re
from collections import Counter, defaultdict
from statistics import median

from matching.conditioned_multilibrary import (
    normalize_source_family, reference_condition, reference_content_fingerprint,
)
from matching.evirank_retrieval import _polarity
from matching.library_quality import normalize_formula_text
from matching.library_record_store import semantic_sha256

ALGORITHM_ID = "mssystem-v16-paired-fragment-ce-comparison-1.0.0"
POLICY = {
    "reference_peak_prevalence_minimum": 0.5,
    "minimum_trajectory_energies": 2,
    "minimum_trajectory_ions": 2,
    "minimum_distinct_supporting_product_ions": 2,
    "minimum_reference_content_profiles_across_energies": 2,
    "trajectory_reference_ce_delta": 0.0,
    "reference_reducer": "median_of_unique_content_profiles",
    "query_replicate_reducer": "median_of_unique_query_content_within_condition",
    "production_rank_override": False,
}


def _compact(value):
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def _ion(value):
    charge = value.get("charge") or value.get("precursor_charge")
    try:
        numeric_charge = float(charge)
        charge = int(numeric_charge) if math.isfinite(numeric_charge) and numeric_charge.is_integer() else None
    except (TypeError, ValueError):
        charge = None
    return (_polarity(value.get("polarity") or value.get("ion_mode")), charge,
            str(value.get("adduct") or value.get("query_adduct") or "").replace(" ", "").lower())


def _condition_key(value, *, include_source=False):
    check = reference_condition(value, value, near_delta_by_unit={"eV": 0., "NCE": 0.})
    if (not check["eligible_for_conditioned_consensus"]
            or check["platform_anchor_scope"] == "analyzer_family_metadata_fallback"
            or not all(_ion(value))):
        return None
    platform = (str(value.get("instrument_domain_id") or "").strip()
                or _compact(value.get("instrument")))
    context = (platform, check["query_instrument_family"], check["query_activation"],
               _compact(value.get("collision_cell_region")), check["collision_energy_unit"], *_ion(value))
    if include_source:
        source = normalize_source_family(value.get("dataset_source"))
        if not source:
            return None
        context += (source,)
    return (context, float(check["query_collision_energy"]))


def _query_bin(query):
    key = _condition_key(query)
    if key is not None:
        return key
    # Incomplete metadata can group technical repeats for descriptive scores,
    # but it is never accepted by the trajectory/conditioned evidence gate.
    fields = ("collision_energy", "collision_energy_unit", "collision_energy_mode",
              "instrument", "instrument_type", "activation_method", "instrument_domain_id",
              "collision_cell_region", "polarity", "adduct", "charge")
    return ("incomplete:" + semantic_sha256({f: query.get(f) for f in fields}), None)


def _unit_area(peaks):
    total = sum(float(p[1]) for p in peaks)
    return [(float(p[0]), float(p[1]) / total) for p in peaks] if total > 0 else []


def _align(query_peaks, reference_peaks, tolerance):
    """Deterministic one-to-one direct product-ion alignment, never mass-shifted."""
    edges = sorted((abs(q[0] - r[0]), -q[1] * r[1], i, j)
                   for i, q in enumerate(query_peaks) for j, r in enumerate(reference_peaks)
                   if abs(q[0] - r[0]) <= tolerance)
    used_q, used_r, result = set(), set(), {}
    for error, _, i, j in edges:
        if i not in used_q and j not in used_r:
            used_q.add(i)
            used_r.add(j)
            result[i] = (reference_peaks[j][1], error, j)
    return result


def _unique_references(references):
    unique = {}
    for item in sorted(references, key=lambda x: str(x["record"].get("spectrum_id") or "")):
        fingerprint = reference_content_fingerprint(item["record"].get("peaks") or [])
        if fingerprint:
            unique.setdefault(fingerprint, item)
    return list(unique.values())


def _profile(references, query_peaks, tolerance):
    references = _unique_references(references)
    aligned = [_align(query_peaks, _unit_area(r["peaks"]), tolerance) for r in references]
    ions = []
    for index in range(len(query_peaks)):
        hits = [i for i, values in enumerate(aligned) if index in values]
        ions.append({"prevalence": len(hits) / len(references), "any_observed": bool(hits),
            "intensity": median(values.get(index, (0., None))[0] for values in aligned),
            "references": [{"spectrum_id": str(references[i]["record"].get("spectrum_id") or ""),
                "content_fingerprint": reference_content_fingerprint(references[i]["record"].get("peaks") or []),
                "mass_error_da": round(aligned[i][index][1], 7)} for i in hits]})
    return {"ions": ions, "count": len(references),
            "cosine": median(float(item["score"]) for item in references),
            "fingerprints": sorted(reference_content_fingerprint(item["record"]["peaks"]) for item in references)}


def _contrast(packet, a, b, *, stratum=None):
    query = packet["query"]
    peaks = _unit_area(packet["query_peaks"])
    supports = {"a": [], "b": []}
    shared, tracks = 0, []
    prevalence = POLICY["reference_peak_prevalence_minimum"]
    for (mz, intensity), pa, pb in zip(peaks, a["ions"], b["ions"]):
        if pa["any_observed"] and pb["any_observed"]:
            shared += 1
        for side, own, other in (("a", pa, pb), ("b", pb, pa)):
            if own["prevalence"] >= prevalence and not other["any_observed"]:
                supports[side].append({"query_mz": round(mz, 6),
                    "neutral_loss_da": round(float(query["precursor_mz"]) - mz, 6),
                    "query_fraction_of_processed_ion_current": round(intensity, 8),
                    "reference_prevalence": round(own["prevalence"], 8),
                    "positive_support": round(intensity * own["prevalence"], 8),
                    "reference_observations": own["references"],
                    "other_candidate_absence_proven": False})
        if pa["prevalence"] >= prevalence and pb["prevalence"] >= prevalence:
            tracks.append({"mz": mz, "query": intensity, "a": pa["intensity"], "b": pb["intensity"]})
    return {"query_id": str(query["spectrum_id"]), "query_fingerprint": reference_content_fingerprint(query["peaks"]),
        "query_bin": _query_bin(query), "view": packet["view"], "stratum": stratum,
        "tolerance_da": packet["tolerance_da"], "precursor_mz": query["precursor_mz"],
        "fragment_delta_a_minus_b": round(sum(p["positive_support"] for p in supports["a"])
                                          - sum(p["positive_support"] for p in supports["b"]), 8),
        "cosine_delta_a_minus_b": round(a["cosine"] - b["cosine"], 8),
        "a_unique_reference_profiles": a["count"], "b_unique_reference_profiles": b["count"],
        "a_fingerprints": a["fingerprints"], "b_fingerprints": b["fingerprints"],
        "shared_observed_product_ions_excluded": shared,
        "a_observed_contrast_ions": supports["a"], "b_observed_contrast_ions": supports["b"],
        "tracks": tracks}


def _common_tracks(items, tolerance):
    if not items:
        return []
    anchors = items[0]["tracks"]
    maps = [_align([(t["mz"], 1.) for t in anchors], [(t["mz"], 1.) for t in item["tracks"]], tolerance)
            for item in items[1:]]
    # Keep the actual one-to-one assignment. No unmatched ion becomes zero.
    output = []
    for index, anchor in enumerate(anchors):
        if not all(index in mapping for mapping in maps):
            continue
        matched = [anchor]
        for item, mapping in zip(items[1:], maps):
            matched.append(item["tracks"][mapping[index][2]])
        output.append({"mz": median(t["mz"] for t in matched),
                       **{key: median(t[key] for t in matched) for key in ("query", "a", "b")}})
    return output


def _collapse(observations):
    grouped = defaultdict(dict)
    for item in sorted(observations, key=lambda x: x["query_id"]):
        key = (item["query_bin"], item["view"], item["stratum"])
        grouped[key].setdefault(item["query_fingerprint"], item)
    bins = []
    for (query_bin, view, stratum), unique in sorted(grouped.items(), key=lambda x: repr(x[0])):
        items = list(unique.values())
        bins.append({"query_bin": query_bin, "view": view, "stratum": stratum,
            "unique_query_content_count": len(items), "query_ids": sorted(x["query_id"] for x in items),
            "fragment_delta_a_minus_b": round(median(x["fragment_delta_a_minus_b"] for x in items), 8),
            "cosine_delta_a_minus_b": round(median(x["cosine_delta_a_minus_b"] for x in items), 8),
            "a_fingerprints": sorted({f for x in items for f in x["a_fingerprints"]}),
            "b_fingerprints": sorted({f for x in items for f in x["b_fingerprints"]}),
            "tracks": _common_tracks(items, items[0]["tolerance_da"]),
            "precursor_min": min(x["precursor_mz"] for x in items),
            "precursor_max": max(x["precursor_mz"] for x in items),
            "tolerance_da": items[0]["tolerance_da"]})
    return bins


def _trajectories(bins):
    grouped = defaultdict(list)
    for item in bins:
        query_context, energy = item["query_bin"]
        if energy is not None and item["stratum"] is not None:
            reference_context, ref_energy = item["stratum"]
            if energy == ref_energy:
                grouped[(query_context, item["view"], reference_context)].append(item)
    result = []
    for (query_context, view, reference_context), items in sorted(grouped.items(), key=lambda x: repr(x[0])):
        items.sort(key=lambda x: x["query_bin"][1])
        if len(items) < POLICY["minimum_trajectory_energies"]:
            continue
        if max(x["precursor_max"] for x in items) - min(x["precursor_min"] for x in items) > 0.01:
            continue
        edges = []
        for low, high in zip(items, items[1:]):
            matched = []
            alignment = _align([(t["mz"], 1.) for t in low["tracks"]],
                               [(t["mz"], 1.) for t in high["tracks"]], low["tolerance_da"])
            for index, ion in enumerate(low["tracks"]):
                if index in alignment:
                    other = high["tracks"][alignment[index][2]]
                    changes = {key: other[key] - ion[key] for key in ("query", "a", "b")}
                    matched.append({"product_mz": round(ion["mz"], 6),
                        **{f"{key}_relative_intensity_change": round(value, 8) for key, value in changes.items()},
                        "a_absolute_change_error": round(abs(changes["query"] - changes["a"]), 8),
                        "b_absolute_change_error": round(abs(changes["query"] - changes["b"]), 8)})
            if len(matched) >= POLICY["minimum_trajectory_ions"]:
                edges.append({"low_ce": low["query_bin"][1], "high_ce": high["query_bin"][1],
                    "observed_ion_count": len(matched), "ions": matched,
                    "trajectory_error_delta_b_minus_a": round(sum(t["b_absolute_change_error"]
                        - t["a_absolute_change_error"] for t in matched) / len(matched), 8)})
        if edges:
            result.append({"query_context": query_context, "reference_context": reference_context,
                "view": view, "energies": [x["query_bin"][1] for x in items], "edges": edges,
                "trajectory_error_delta_b_minus_a": round(median(x["trajectory_error_delta_b_minus_a"] for x in edges), 8),
                "a_unique_reference_contents": len({f for x in items for f in x["a_fingerprints"]}),
                "b_unique_reference_contents": len({f for x in items for f in x["b_fingerprints"]})})
    return result


def build_candidate_comparison(ids, packets, *, candidate_formulae, baseline_connectivity=None):
    ids = sorted(ids)
    if baseline_connectivity is not None:
        baseline_connectivity = str(baseline_connectivity).split(":")[-1].split("-")[0].upper()
        if baseline_connectivity not in ids:
            raise ValueError("Baseline connectivity is outside the frozen pool")
    formulae = {key: sorted({normalized for f in candidate_formulae.get(key, [])
                            if (normalized := normalize_formula_text(f))}) for key in ids}
    indexed, strict_rejections = [], Counter()
    for packet in sorted(packets, key=lambda p: (str(p["query"]["spectrum_id"]), p["view"])):
        all_refs, strict = defaultdict(list), defaultdict(lambda: defaultdict(list))
        for item in packet["references"]:
            record = item["record"]
            identity = record["connectivity_identity"]
            all_refs[identity].append(item)
            query_key, reference_key = _condition_key(packet["query"]), _condition_key(record, include_source=True)
            condition = item["condition"]
            if (query_key is not None and reference_key is not None
                    and condition["eligible_for_conditioned_consensus"]
                    and condition["platform_anchor_scope"] != "analyzer_family_metadata_fallback"
                    and condition["collision_energy_delta"] == 0.0
                    and _ion(record) == _ion(packet["query"])):
                strict[identity][reference_key].append(item)
            else:
                strict_rejections[condition["tier"]] += 1
        peaks = _unit_area(packet["query_peaks"])
        indexed.append((packet,
            {key: _profile(refs, peaks, packet["tolerance_da"]) for key, refs in all_refs.items()},
            {key: {stratum: _profile(refs, peaks, packet["tolerance_da"]) for stratum, refs in groups.items()}
             for key, groups in strict.items()}))
    pairs, wins, fragments, cosines = [], Counter(), defaultdict(list), defaultdict(list)
    for a, b in itertools.combinations(ids, 2):
        pair = {"a": a, "b": b, "same_known_formula": len(formulae[a]) == len(formulae[b]) == 1 and formulae[a] == formulae[b]}
        all_observations, matched = [], []
        if pair["same_known_formula"]:
            for packet, all_refs, strict in indexed:
                if a not in all_refs or b not in all_refs:
                    continue
                all_observations.append(_contrast(packet, all_refs[a], all_refs[b]))
                for key in sorted(set(strict.get(a, {})) & set(strict.get(b, {})), key=repr):
                    matched.append(_contrast(packet, strict[a][key], strict[b][key], stratum=key))
        all_bins, matched_bins = _collapse(all_observations), _collapse(matched)
        trajectories = _trajectories(matched_bins)
        delta = median(x["fragment_delta_a_minus_b"] for x in all_bins) if all_bins else None
        cosine_delta = median(x["cosine_delta_a_minus_b"] for x in all_bins) if all_bins else None
        preferred = None
        required_views = {packet["view"] for packet in packets}
        for side, sign, identity in (("a", 1, a), ("b", -1, b)):
            observed = []
            tolerance = max((x["tolerance_da"] for x in matched), default=0.02)
            for mz in sorted(ion["query_mz"] for x in matched for ion in x[f"{side}_observed_contrast_ions"]):
                if not observed or mz - observed[-1] > tolerance:
                    observed.append(mz)
            if (matched_bins and trajectories and len(observed) >= POLICY["minimum_distinct_supporting_product_ions"]
                    and {x["view"] for x in trajectories} == required_views
                    and all(sign * x["fragment_delta_a_minus_b"] > 0 and sign * x["cosine_delta_a_minus_b"] > 0 for x in matched_bins)
                    and all(sign * x["trajectory_error_delta_b_minus_a"] > 0
                            and min(x["a_unique_reference_contents"], x["b_unique_reference_contents"]) >= 2 for x in trajectories)):
                preferred = identity
        if preferred:
            wins[preferred] += 1
        if delta is not None:
            fragments[a].append(delta)
            fragments[b].append(-delta)
            cosines[a].append(cosine_delta)
            cosines[b].append(-cosine_delta)
        detailed = (baseline_connectivity or ids[0]) in {a, b}
        pair.update({"status": "compared" if all_bins else "unavailable",
            "reason": None if all_bins else ("no_shared_measured_queries" if pair["same_known_formula"] else "formula_missing_conflicting_or_different"),
            "exploratory_fragment_delta_a_minus_b": round(delta, 8) if delta is not None else None,
            "exploratory_cosine_delta_a_minus_b": round(cosine_delta, 8) if cosine_delta is not None else None,
            "exact_condition_paired_bin_count": len(matched_bins),
            "exploratory_paired_bins": [{k: v for k, v in x.items() if k != "tracks"} for x in all_bins],
            "exact_condition_paired_bins": [{k: v for k, v in x.items() if k != "tracks"} for x in matched_bins],
            "ce_trajectories": trajectories, "strict_evidence_preferred_candidate": preferred,
            "detailed_observations_included": detailed,
            "observations": [{k: v for k, v in x.items() if k != "tracks"} for x in all_observations] if detailed else [],
            "condition_matched_observations": [{k: v for k, v in x.items() if k != "tracks"} for x in matched] if detailed else [],
            "primary_ranking_changed": False})
        pairs.append(pair)
    candidates = [{"connectivity": key, "formulae": formulae[key], "strict_pair_wins": wins[key],
        "compared_same_formula_candidates": len(fragments[key]),
        "median_positive_fragment_contrast": round(median(fragments[key]), 8) if fragments[key] else None,
        "median_paired_cosine_delta": round(median(cosines[key]), 8) if cosines[key] else None} for key in ids]
    candidates.sort(key=lambda x: (-x["strict_pair_wins"],
        -(x["median_positive_fragment_contrast"] if x["median_positive_fragment_contrast"] is not None else -2),
        -(x["median_paired_cosine_delta"] if x["median_paired_cosine_delta"] is not None else -2), x["connectivity"]))
    for rank, item in enumerate(candidates, 1):
        item["review_rank"] = rank
    result = {"algorithm_id": ALGORITHM_ID, "policy": dict(POLICY), "candidates": candidates, "pairs": pairs,
        "baseline_connectivity": baseline_connectivity, "primary_selected_connectivity": baseline_connectivity,
        "candidate_count": len(ids), "pair_count": len(pairs), "compared_pair_count": sum(x["status"] == "compared" for x in pairs),
        "exact_condition_pair_count": sum(x["exact_condition_paired_bin_count"] > 0 for x in pairs),
        "trajectory_pair_count": sum(bool(x["ce_trajectories"]) for x in pairs),
        "strict_preference_pair_count": sum(bool(x["strict_evidence_preferred_candidate"]) for x in pairs),
        "strict_reference_rejection_tiers": dict(strict_rejections),
        "primary_ranking_changed": False, "score_is_probability": False,
        "review_order_is_identity_assignment": False, "query_absence_penalties_applied": False,
        "product_ion_and_derived_loss_counted_twice": False, "spectral_views_are_independent_replicates": False,
        "conditionless_reference_scores_are_exploratory": True, "candidate_names_or_truth_used": False,
        "technical_repeats_collapsed": True, "reference_mirrors_collapsed": True,
        "cross_unit_conversion_applied": False}
    result["semantic_sha256"] = semantic_sha256(result)
    return result
