"""Evidence-aware molecular networking for MSsystem V7.

The V6 network treated every MS2 scan as a molecular family member, used one
fixed fragment binning as a proxy for modified cosine, and counted singleton
components as molecular families.  This module keeps the fast sparse candidate
search, but makes the scientific decisions on exact peak alignments:

* OpenMS chromatographic MS1 features become feature-consensus nodes whenever
  available; only unlinked spectra fall back to an explicitly labelled MS2
  precursor-event node;
* direct, precursor-shifted (modified-cosine), and neutral-loss evidence are
  scored separately;
* validated Ion Identity Molecular Network (IIMN) memberships are an explicit
  evidence layer, never silently inferred from precursor mass alone;
* mutual/either Top-K policy is explicit and reported;
* singleton components remain visible as singletons but are not called
  molecular families.

The public result remains backward compatible: ``edges`` is still a list of
``(source_result_index, target_result_index, score, matched_peaks)`` tuples.
The richer ``edge_details`` collection powers V7 UI and exports.
"""
from __future__ import annotations

import bisect
import math
from collections import Counter, OrderedDict, defaultdict
from typing import Any, Iterable

import numpy as np
import scipy.sparse as sp

from utils.spectrum import clean_peaks


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _collision_energy(value: Any) -> float | None:
    if isinstance(value, dict):
        value = value.get("value")
    return _number(value)


def _normalised_peaks(peaks: Iterable[Iterable[float]]) -> list[tuple[float, float]]:
    parsed: list[tuple[float, float]] = []
    for peak in peaks or []:
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if mz > 0 and intensity > 0 and math.isfinite(mz) and math.isfinite(intensity):
            parsed.append((mz, intensity))
    if not parsed:
        return []
    # Square-root intensity weighting is conventional for molecular networking:
    # it prevents one base peak from overwhelming reproducible lower peaks.
    weighted = [(mz, math.sqrt(intensity)) for mz, intensity in parsed]
    norm = math.sqrt(sum(intensity * intensity for _, intensity in weighted))
    if norm <= 0:
        return []
    return sorted((mz, intensity / norm) for mz, intensity in weighted)


def _consensus_peaks(
    members: list[dict],
    *,
    precursor_mz: float,
    top_n: int,
    rel_threshold: float,
    consensus_tolerance_da: float,
) -> tuple[list[tuple[float, float]], dict]:
    """Create a deterministic consensus from repeated precursor-event scans."""
    cleaned_members: list[list[tuple[float, float]]] = []
    raw_members = 0
    for row in members:
        raw = row.get("raw_peaks") or []
        source = raw or row.get("query_peaks") or []
        if raw:
            raw_members += 1
        cleaned = clean_peaks(
            source,
            precursor_mz,
            top_n=top_n,
            rel_threshold=rel_threshold,
            precursor_charge=row.get("precursor_charge"),
        )
        if cleaned:
            peak_max = max((float(intensity) for _, intensity in cleaned), default=0.0)
            if peak_max > 0:
                cleaned_members.append([
                    (float(mz), float(intensity) / peak_max)
                    for mz, intensity in cleaned
                ])

    if not cleaned_members:
        return [], {
            "raw_member_count": raw_members,
            "usable_member_count": 0,
            "peak_source": "missing",
        }
    if len(cleaned_members) == 1:
        return cleaned_members[0], {
            "raw_member_count": raw_members,
            "usable_member_count": 1,
            "peak_source": "raw" if raw_members else "search_fallback",
        }

    # Use nearest bins only to form a consensus. Exact inter-node comparison is
    # performed later and is not affected by this binning.
    per_bin: dict[int, list[tuple[float, float, int]]] = defaultdict(list)
    width = max(0.005, float(consensus_tolerance_da))
    for member_index, peaks in enumerate(cleaned_members):
        best_for_member: dict[int, tuple[float, float]] = {}
        for mz, intensity in peaks:
            bin_index = int(round(mz / width))
            previous = best_for_member.get(bin_index)
            if previous is None or intensity > previous[1]:
                best_for_member[bin_index] = (mz, intensity)
        for bin_index, (mz, intensity) in best_for_member.items():
            per_bin[bin_index].append((mz, intensity, member_index))

    minimum_occurrence = max(1, int(math.ceil(len(cleaned_members) * 0.30)))
    consensus: list[tuple[float, float]] = []
    for observations in per_bin.values():
        if len({item[2] for item in observations}) < minimum_occurrence:
            continue
        total_intensity = sum(item[1] for item in observations)
        mz = (
            sum(item[0] * item[1] for item in observations) / total_intensity
            if total_intensity > 0 else float(np.median([item[0] for item in observations]))
        )
        # Missing observations contribute zero. This rewards reproducible peaks
        # without deleting diagnostic peaks from two-scan stepped-energy data.
        intensity = total_intensity / len(cleaned_members)
        consensus.append((mz, intensity))
    consensus.sort(key=lambda peak: peak[1], reverse=True)
    consensus = sorted(consensus[:top_n])
    return consensus, {
        "raw_member_count": raw_members,
        "usable_member_count": len(cleaned_members),
        "peak_source": (
            "raw" if raw_members == len(members)
            else "mixed" if raw_members
            else "search_fallback"
        ),
    }


def prepare_feature_nodes(
    results: list[dict],
    *,
    top_n: int = 70,
    rel_threshold: float = 0.01,
    consensus_tolerance_da: float = 0.05,
) -> tuple[list[dict], dict]:
    """Build nodes from true LC-MS features with an explicit MS2 fallback."""
    groups: OrderedDict[str, list[tuple[int, dict]]] = OrderedDict()
    excluded_blank = 0
    missing_precursor = 0
    missing_peak_list = 0
    lcms_linked_spectra = 0
    fallback_spectra = 0
    for result_index, row in enumerate(results or []):
        if not isinstance(row, dict):
            continue
        if row.get("blank_excluded"):
            excluded_blank += 1
            continue
        precursor_mz = _number(row.get("precursor_mz"))
        if precursor_mz is None or precursor_mz <= 0:
            missing_precursor += 1
            continue
        if not (row.get("raw_peaks") or row.get("query_peaks")):
            missing_peak_list += 1
            continue
        lcms_feature_id = str(row.get("lcms_feature_id") or "")
        if lcms_feature_id:
            feature_id = f"lcms:{lcms_feature_id}"
            lcms_linked_spectra += 1
        else:
            feature_id = str(
                row.get("precursor_feature_id") or f"scan:{result_index}"
            )
            fallback_spectra += 1
        groups.setdefault(feature_id, []).append((result_index, row))

    nodes: list[dict] = []
    spectra_in_nodes = 0
    for feature_id, indexed_members in groups.items():
        member_indices = [item[0] for item in indexed_members]
        members = [item[1] for item in indexed_members]
        precursor_values = [
            value for value in (_number(row.get("precursor_mz")) for row in members)
            if value is not None
        ]
        if not precursor_values:
            continue
        precursor_mz = float(np.median(precursor_values))
        # RS artifacts are condition-specific measured representatives, not
        # feature-wide consensus spectra. Preserve the producer's explicitly
        # ranked artifact (RS:<feature>:1), regardless of saved row order.
        precomputed_rows = [
            (index, row) for index, row in indexed_members
            if row.get("lcms_consensus_peaks")
        ]
        representative_rows = [
            item for item in precomputed_rows
            if str(item[1].get("lcms_consensus_spectrum_id") or "").startswith("RS:")
        ]
        peak_representative = None
        if representative_rows:
            def artifact_order(item):
                artifact_id = str(item[1].get("lcms_consensus_spectrum_id") or "")
                suffix = artifact_id.rsplit(":", 1)[-1]
                return (int(suffix) if suffix.isdigit() else float("inf"),
                        artifact_id, str(item[1].get("spectrum_id") or ""))
            peak_representative = min(representative_rows, key=artifact_order)
            precomputed = peak_representative[1]["lcms_consensus_peaks"]
        else:
            precomputed = precomputed_rows[0][1]["lcms_consensus_peaks"] if precomputed_rows else None
        if precomputed:
            peaks = [
                (float(peak[0]), float(peak[1]))
                for peak in precomputed
                if len(peak) >= 2
            ]
            peak_meta = {
                "raw_member_count": sum(
                    bool(row.get("raw_peaks")) for row in members
                ),
                "usable_member_count": 1 if peak_representative else len(members),
                "peak_source": ("lcms_feature_representative" if peak_representative
                                else "lcms_feature_consensus"),
            }
        else:
            peaks, peak_meta = _consensus_peaks(
                members,
                precursor_mz=precursor_mz,
                top_n=top_n,
                rel_threshold=rel_threshold,
                consensus_tolerance_da=consensus_tolerance_da,
            )
        similarity_eligible = len(peaks) >= 3

        matches = [
            (row.get("match"), index)
            for index, row in indexed_members
            if isinstance(row.get("match"), dict) and row["match"].get("compound_name")
        ]
        best_match: dict | None = None
        best_match_index: int | None = None
        if matches:
            best_match, best_match_index = max(
                matches,
                key=lambda item: (
                    _number(item[0].get("final_score")) or 0.0,
                    _number(item[0].get("cosine")) or 0.0,
                    -item[1],
                ),
            )
        representative_index = (
            best_match_index if best_match_index is not None else member_indices[0]
        )
        representative = results[representative_index]
        retention_times = [
            value for value in (_number(row.get("retention_time")) for row in members)
            if value is not None
        ]
        energies = sorted({
            round(value, 4)
            for value in (_collision_energy(row.get("collision_energy")) for row in members)
            if value is not None
        })
        iimn_ids = [
            str(row.get("ion_identity_network_id"))
            for row in members if row.get("ion_identity_network_id")
        ]
        iimn_id = Counter(iimn_ids).most_common(1)[0][0] if iimn_ids else None
        linked_annotations = [
            row.get("ion_identity_linked_annotation")
            for row in members
            if isinstance(row.get("ion_identity_linked_annotation"), dict)
        ]
        linked_annotation = linked_annotations[0] if linked_annotations else None

        match = best_match or {}
        nodes.append({
            "idx": representative_index,
            "feature_id": feature_id,
            "feature_source": (
                "openms_lcms_feature"
                if feature_id.startswith("lcms:")
                else "ms2_precursor_event_fallback"
            ),
            "member_result_indices": member_indices,
            "spectrum_count": len(members),
            "spectrum_id": representative.get("spectrum_id"),
            "member_spectrum_ids": [row.get("spectrum_id") for row in members],
            "precursor_mz": precursor_mz,
            "retention_time": (
                float(np.median(retention_times)) if retention_times else None
            ),
            "retention_time_min": min(retention_times) if retention_times else None,
            "retention_time_max": max(retention_times) if retention_times else None,
            "peaks": peaks,
            "peaks_count": len(peaks),
            "peak_source": peak_meta["peak_source"],
            "peak_evidence": ({
                "artifact_type": "representative_spectrum",
                "spectrum_id": peak_representative[1].get("spectrum_id"),
                "result_index": peak_representative[0],
                "artifact_id": peak_representative[1].get("lcms_consensus_spectrum_id"),
                "collision_energy": peak_representative[1].get("collision_energy"),
                "collision_energy_unit": peak_representative[1].get("collision_energy_unit"),
                "collision_energy_mode": peak_representative[1].get("collision_energy_mode"),
                "activation_method": peak_representative[1].get("activation_method"),
                "precursor_mz": peak_representative[1].get("precursor_mz"),
                "retention_time": peak_representative[1].get("retention_time"),
                "ion_mode": peak_representative[1].get("ion_mode"),
                "precursor_charge": peak_representative[1].get("precursor_charge"),
                "instrument_configuration_id": peak_representative[1].get(
                    "instrument_configuration_id"
                ),
                "selection": "saved_candidate_independent_primary_artifact",
                "annotation_result_index": best_match_index,
            } if peak_representative else None),
            "network_eligibility": (
                "similarity_eligible"
                if similarity_eligible
                else "singleton_insufficient_fragments"
            ),
            "raw_member_count": peak_meta["raw_member_count"],
            "usable_member_count": peak_meta["usable_member_count"],
            "match": best_match,
            "compound_name": match.get("compound_name"),
            "classification": match.get("classification"),
            "confidence": match.get("confidence"),
            "evidence_tier": match.get("evidence_tier"),
            "final_score": match.get("final_score"),
            "cosine": match.get("cosine"),
            "entropy": match.get("entropy"),
            "mass_ppm": match.get("mass_ppm"),
            "smiles": match.get("smiles"),
            "inchikey": match.get("inchikey"),
            "formula": match.get("formula"),
            "adduct": match.get("adduct"),
            "dataset_source": match.get("dataset_source"),
            "source_tier": match.get("source_tier"),
            "ion_mode": representative.get("ion_mode"),
            "precursor_charge": representative.get("precursor_charge"),
            "collision_energy": energies[0] if len(energies) == 1 else None,
            "collision_energies": energies,
            "ion_identity_status": representative.get("ion_identity_status"),
            "ion_identity_network_id": iimn_id,
            "ion_identity": representative.get("ion_identity"),
            "ion_identity_linked_annotation": linked_annotation,
            "lcms_feature_id": representative.get("lcms_feature_id"),
            "lcms_feature_intensity": representative.get(
                "lcms_feature_intensity"
            ),
            "lcms_feature_charge": representative.get("lcms_feature_charge"),
            "lcms_isotope_trace_count": representative.get(
                "lcms_isotope_trace_count"
            ),
            "lcms_adduct": representative.get("lcms_adduct"),
            "lcms_neutral_mass": representative.get("lcms_neutral_mass"),
            "lcms_adduct_group_id": representative.get(
                "lcms_adduct_group_id"
            ),
            "lcms_consensus_spectrum_id": representative.get(
                "lcms_consensus_spectrum_id"
            ),
        })
        spectra_in_nodes += len(members)

    return nodes, {
        "input_spectra": len(results or []),
        "spectra_in_nodes": spectra_in_nodes,
        "excluded_blank_spectra": excluded_blank,
        "missing_precursor_spectra": missing_precursor,
        "missing_peak_list_spectra": missing_peak_list,
        "omitted_input_spectra": max(
            0,
            len(results or [])
            - excluded_blank
            - spectra_in_nodes,
        ),
        "input_node_coverage_pct": round(
            100.0
            * spectra_in_nodes
            / max(1, len(results or []) - excluded_blank),
            3,
        ),
        "feature_groups_considered": len(groups),
        "lcms_linked_spectra": lcms_linked_spectra,
        "fallback_spectra": fallback_spectra,
        "lcms_feature_nodes": sum(
            str(node.get("feature_source")) == "openms_lcms_feature"
            for node in nodes
        ),
        "fallback_precursor_event_nodes": sum(
            str(node.get("feature_source")) == "ms2_precursor_event_fallback"
            for node in nodes
        ),
        "consensus_nodes": sum(node["spectrum_count"] > 1 for node in nodes),
        "similarity_eligible_nodes": sum(
            node.get("network_eligibility") == "similarity_eligible"
            for node in nodes
        ),
        "insufficient_fragment_nodes": sum(
            node.get("network_eligibility")
            == "singleton_insufficient_fragments"
            for node in nodes
        ),
        "collapsed_spectra": max(0, spectra_in_nodes - len(nodes)),
    }


def _view_peaks(node: dict, view: str) -> list[tuple[float, float]]:
    if view == "fragment":
        return node["peaks"]
    precursor = float(node["precursor_mz"])
    return [
        (precursor - float(mz), float(intensity))
        for mz, intensity in node["peaks"]
        if precursor - float(mz) > 0.5
    ]


def _sparse_view(
    nodes: list[dict],
    *,
    view: str,
    bin_width: float,
    offset_fraction: float,
) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    offset = float(offset_fraction) * bin_width
    for row_index, node in enumerate(nodes):
        bin_max: dict[int, float] = {}
        for mz, intensity in _normalised_peaks(_view_peaks(node, view)):
            bin_index = int(math.floor((mz + offset) / bin_width))
            if intensity > bin_max.get(bin_index, 0.0):
                bin_max[bin_index] = intensity
        for bin_index, intensity in bin_max.items():
            rows.append(row_index)
            cols.append(bin_index)
            values.append(intensity)
    if not rows:
        empty = sp.csr_matrix((len(nodes), 1), dtype=np.float32)
        return empty, empty
    matrix = sp.coo_matrix(
        (
            np.asarray(values, dtype=np.float32),
            (np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)),
        ),
        shape=(len(nodes), max(cols) + 1),
    ).tocsr()
    norms = np.sqrt(np.asarray(matrix.multiply(matrix).sum(axis=1)).ravel())
    norms[norms == 0] = 1.0
    normalised = (sp.diags(1.0 / norms) @ matrix).tocsr()
    occupied = (matrix != 0).astype(np.float32).tocsr()
    return normalised, occupied


def _sparse_candidate_pairs(
    nodes: list[dict],
    *,
    mz_tolerance_da: float,
    precursor_window: float,
    score_floor: float,
    overlap_floor: int,
    block_size: int = 384,
) -> set[tuple[int, int]]:
    """Generate a recall-oriented candidate union from fragment/loss views."""
    n = len(nodes)
    if n < 2:
        return set()
    precursor = np.asarray([node["precursor_mz"] for node in nodes], dtype=float)
    candidates: set[tuple[int, int]] = set()

    # Exact enumeration is cheap for small tasks and guarantees that shifted
    # analogues are not missed by the approximate sparse prefilter.
    if n <= 250:
        for i in range(n):
            for j in range(i + 1, n):
                if abs(precursor[i] - precursor[j]) <= precursor_window:
                    candidates.add((i, j))
        return candidates

    for view in ("fragment", "neutral_loss"):
        for offset in (0.0, 1.0 / 3.0, 2.0 / 3.0):
            spectra, occupied = _sparse_view(
                nodes,
                view=view,
                bin_width=mz_tolerance_da,
                offset_fraction=offset,
            )
            spectra_t = spectra.T.tocsr()
            occupied_t = occupied.T.tocsr()
            for start in range(0, n, block_size):
                end = min(n, start + block_size)
                scores = (spectra[start:end] @ spectra_t).toarray()
                overlaps = (
                    occupied[start:end] @ occupied_t
                ).toarray().astype(np.int16)
                for local_index in range(end - start):
                    i = start + local_index
                    mask = (
                        (scores[local_index] >= score_floor)
                        & (overlaps[local_index] >= overlap_floor)
                    )
                    for j in np.flatnonzero(mask):
                        j = int(j)
                        if j <= i:
                            continue
                        if abs(precursor[i] - precursor[j]) <= precursor_window:
                            candidates.add((i, j))
    return candidates


def _alignment_candidates(
    peaks_a: list[tuple[float, float]],
    peaks_b: list[tuple[float, float]],
    *,
    shift_da: float,
    tolerance_da: float,
) -> list[tuple[float, int, int]]:
    mz_b = [peak[0] for peak in peaks_b]
    candidates: list[tuple[float, int, int]] = []
    for index_a, (mz_a, intensity_a) in enumerate(peaks_a):
        target = mz_a + shift_da
        left = bisect.bisect_left(mz_b, target - tolerance_da)
        right = bisect.bisect_right(mz_b, target + tolerance_da)
        for index_b in range(left, right):
            candidates.append((
                intensity_a * peaks_b[index_b][1],
                index_a,
                index_b,
            ))
    return candidates


def _greedy_alignment(
    peaks_a: list[tuple[float, float]],
    peaks_b: list[tuple[float, float]],
    *,
    shifts_da: tuple[float, ...],
    tolerance_da: float,
) -> tuple[float, int]:
    norm_a = _normalised_peaks(peaks_a)
    norm_b = _normalised_peaks(peaks_b)
    if not norm_a or not norm_b:
        return 0.0, 0
    candidates: dict[tuple[int, int], float] = {}
    for shift in shifts_da:
        for product, index_a, index_b in _alignment_candidates(
            norm_a,
            norm_b,
            shift_da=shift,
            tolerance_da=tolerance_da,
        ):
            key = (index_a, index_b)
            if product > candidates.get(key, 0.0):
                candidates[key] = product
    used_a: set[int] = set()
    used_b: set[int] = set()
    score = 0.0
    matched = 0
    for (index_a, index_b), product in sorted(
        candidates.items(), key=lambda item: item[1], reverse=True
    ):
        if index_a in used_a or index_b in used_b:
            continue
        used_a.add(index_a)
        used_b.add(index_b)
        score += product
        matched += 1
    return min(1.0, float(score)), matched


def _exact_edge(
    node_a: dict,
    node_b: dict,
    *,
    cosine_threshold: float,
    min_matched_peaks: int,
    mz_tolerance_da: float,
) -> dict | None:
    precursor_delta = float(node_b["precursor_mz"] - node_a["precursor_mz"])
    direct, direct_n = _greedy_alignment(
        node_a["peaks"],
        node_b["peaks"],
        shifts_da=(0.0,),
        tolerance_da=mz_tolerance_da,
    )
    modified, modified_n = _greedy_alignment(
        node_a["peaks"],
        node_b["peaks"],
        shifts_da=(0.0, precursor_delta),
        tolerance_da=mz_tolerance_da,
    )
    loss_a = _view_peaks(node_a, "neutral_loss")
    loss_b = _view_peaks(node_b, "neutral_loss")
    neutral_loss, neutral_loss_n = _greedy_alignment(
        loss_a,
        loss_b,
        shifts_da=(0.0,),
        tolerance_da=mz_tolerance_da,
    )

    edge_type: str | None = None
    score = 0.0
    matched = 0
    if direct >= cosine_threshold and direct_n >= min_matched_peaks:
        edge_type = "spectral"
        score, matched = direct, direct_n
    elif (
        abs(precursor_delta) > mz_tolerance_da
        and modified >= cosine_threshold
        and modified_n >= min_matched_peaks
        and max(direct, neutral_loss) >= max(0.30, cosine_threshold * 0.50)
    ):
        edge_type = "analog"
        score, matched = modified, modified_n
    elif (
        neutral_loss >= cosine_threshold
        and neutral_loss_n >= min_matched_peaks
        and modified >= max(0.35, cosine_threshold * 0.70)
    ):
        edge_type = "neutral_loss"
        score, matched = neutral_loss, neutral_loss_n
    if edge_type is None:
        return None
    return {
        "score": round(score, 8),
        "cosine": round(score, 8),
        "matched_peaks": int(matched),
        "edge_type": edge_type,
        "evidence_layers": [edge_type],
        "direct_cosine": round(direct, 8),
        "direct_matched_peaks": int(direct_n),
        "modified_cosine": round(modified, 8),
        "modified_matched_peaks": int(modified_n),
        "neutral_loss_cosine": round(neutral_loss, 8),
        "neutral_loss_matched_peaks": int(neutral_loss_n),
        "precursor_delta_da": round(abs(precursor_delta), 8),
    }


def _topk_edge_filter(
    edges: list[dict],
    *,
    topk: int,
    policy: str,
) -> list[dict]:
    if topk <= 0 or policy == "off" or not edges:
        return edges
    by_node: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for edge_index, edge in enumerate(edges):
        by_node[edge["source_pos"]].append((edge["score"], edge_index))
        by_node[edge["target_pos"]].append((edge["score"], edge_index))
    selected: dict[int, set[int]] = {}
    for node, values in by_node.items():
        values.sort(
            key=lambda item: (
                item[0],
                edges[item[1]]["matched_peaks"],
                -item[1],
            ),
            reverse=True,
        )
        selected[node] = {edge_index for _, edge_index in values[:topk]}
    kept: list[dict] = []
    for edge_index, edge in enumerate(edges):
        left = edge_index in selected.get(edge["source_pos"], set())
        right = edge_index in selected.get(edge["target_pos"], set())
        if (policy == "mutual" and left and right) or (
            policy == "either" and (left or right)
        ):
            kept.append(edge)
    return kept


class _DisjointSet:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left == root_right:
            return
        if self.size[root_left] < self.size[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        self.size[root_left] += self.size[root_right]


def _component_cap(
    edges: list[dict],
    *,
    node_count: int,
    max_component_size: int,
) -> list[dict]:
    if max_component_size <= 0:
        return edges
    union = _DisjointSet(node_count)
    kept: list[dict] = []
    for edge in sorted(
        edges,
        key=lambda item: (item["score"], item["matched_peaks"]),
        reverse=True,
    ):
        left, right = edge["source_pos"], edge["target_pos"]
        root_left, root_right = union.find(left), union.find(right)
        if root_left == root_right:
            kept.append(edge)
            continue
        if union.size[root_left] + union.size[root_right] <= max_component_size:
            union.union(left, right)
            kept.append(edge)
    return sorted(kept, key=lambda item: (item["source_pos"], item["target_pos"]))


def _iimn_edges(nodes: list[dict]) -> list[dict]:
    """Connect only IIMN groups already validated from continuous MS1 traces."""
    groups: dict[str, list[int]] = defaultdict(list)
    for position, node in enumerate(nodes):
        network_id = str(node.get("ion_identity_network_id") or "")
        if network_id and str(node.get("ion_identity_status") or "") == "networked":
            groups[network_id].append(position)
    edges: list[dict] = []
    for network_id, positions in groups.items():
        positions.sort()
        # A star keeps the graph interpretable while preserving connectivity.
        representative = max(
            positions,
            key=lambda position: (
                bool(nodes[position].get("match")),
                _number(nodes[position].get("final_score")) or 0.0,
                -position,
            ),
        )
        for position in positions:
            if position == representative:
                continue
            left, right = sorted((representative, position))
            edges.append({
                "source_pos": left,
                "target_pos": right,
                "score": 1.0,
                "cosine": 0.0,
                "matched_peaks": 0,
                "edge_type": "ion_identity",
                "evidence_layers": ["ion_identity"],
                "direct_cosine": 0.0,
                "direct_matched_peaks": 0,
                "modified_cosine": 0.0,
                "modified_matched_peaks": 0,
                "neutral_loss_cosine": 0.0,
                "neutral_loss_matched_peaks": 0,
                "precursor_delta_da": round(
                    abs(nodes[left]["precursor_mz"] - nodes[right]["precursor_mz"]),
                    8,
                ),
                "ion_identity_network_id": network_id,
            })
    return edges


def _merge_edge_layers(edges: list[dict]) -> list[dict]:
    merged: dict[tuple[int, int], dict] = {}
    for edge in edges:
        key = tuple(sorted((edge["source_pos"], edge["target_pos"])))
        if key not in merged:
            merged[key] = dict(edge)
            continue
        current = merged[key]
        current["evidence_layers"] = sorted(set(
            current.get("evidence_layers", [])
            + edge.get("evidence_layers", [])
        ))
        if edge["score"] > current["score"] and edge["edge_type"] != "ion_identity":
            for field in (
                "score", "cosine", "matched_peaks", "edge_type",
                "direct_cosine", "direct_matched_peaks",
                "modified_cosine", "modified_matched_peaks",
                "neutral_loss_cosine", "neutral_loss_matched_peaks",
                "precursor_delta_da",
            ):
                current[field] = edge[field]
        if edge.get("ion_identity_network_id"):
            current["ion_identity_network_id"] = edge["ion_identity_network_id"]
    return sorted(
        merged.values(), key=lambda item: (item["source_pos"], item["target_pos"])
    )


def _components(node_count: int, edges: list[dict]) -> list[list[int]]:
    union = _DisjointSet(node_count)
    for edge in edges:
        union.union(edge["source_pos"], edge["target_pos"])
    groups: dict[int, list[int]] = defaultdict(list)
    for position in range(node_count):
        groups[union.find(position)].append(position)
    return sorted(groups.values(), key=lambda values: (-len(values), values[0]))


def _identity_key(node: dict) -> str:
    inchikey = str(node.get("inchikey") or "").strip().upper()
    if inchikey:
        return inchikey
    return str(node.get("compound_name") or "").strip().lower()


def _annotate_nodes(
    nodes: list[dict],
    edges: list[dict],
    components: list[list[int]],
) -> list[dict]:
    degree = Counter()
    neighbours: dict[int, list[tuple[int, dict]]] = defaultdict(list)
    for edge in edges:
        left, right = edge["source_pos"], edge["target_pos"]
        degree[left] += 1
        degree[right] += 1
        neighbours[left].append((right, edge))
        neighbours[right].append((left, edge))

    connected_family_number = 0
    for component_number, members in enumerate(components):
        family_id: int | None = None
        if len(members) >= 2:
            family_id = connected_family_number
            connected_family_number += 1
        for position in members:
            nodes[position]["component_id"] = component_number
            nodes[position]["family_id"] = family_id
            nodes[position]["degree"] = int(degree[position])

    for position, node in enumerate(nodes):
        node["propagated_class"] = None
        node["propagated_compound"] = None
        node["propagation_evidence"] = None
        if node.get("match"):
            node["annotation_source"] = "library"
            continue
        linked = node.get("ion_identity_linked_annotation")
        if isinstance(linked, dict) and linked.get("compound_name"):
            node["annotation_source"] = "ion_identity"
            node["propagated_compound"] = linked.get("compound_name")
            node["propagation_evidence"] = (
                "validated MS1 co-elution/adduct link; not independent MS/MS identification"
            )
            continue

        annotated_neighbours = []
        for neighbour_position, edge in neighbours.get(position, []):
            neighbour = nodes[neighbour_position]
            if not neighbour.get("match") or edge["edge_type"] == "ion_identity":
                continue
            annotated_neighbours.append((neighbour, edge))
        if not annotated_neighbours:
            node["annotation_source"] = "unknown"
            continue
        classes = [
            item[0].get("classification")
            for item in annotated_neighbours
            if item[0].get("classification")
            and str(item[0]["classification"]).lower() != "unknown"
        ]
        if classes:
            top_class, top_count = Counter(classes).most_common(1)[0]
            if top_count / len(classes) >= 0.60:
                node["propagated_class"] = top_class
        identities = [
            _identity_key(item[0]) for item in annotated_neighbours
            if _identity_key(item[0])
        ]
        # A compound name is propagated only when two independent annotated
        # neighbours agree. One neighbour supports class-level interpretation.
        if identities:
            identity, count = Counter(identities).most_common(1)[0]
            if count >= 2:
                representative = next(
                    item[0] for item in annotated_neighbours
                    if _identity_key(item[0]) == identity
                )
                node["propagated_compound"] = representative.get("compound_name")
        if node["propagated_class"] or node["propagated_compound"]:
            node["annotation_source"] = "propagated"
            node["propagation_evidence"] = "direct evidence-bearing network neighbours"
        else:
            node["annotation_source"] = "unknown"
    return nodes


def build_evidence_network(
    results: list[dict],
    *,
    cos_thresh: float = 0.70,
    min_peaks: int = 6,
    topk: int = 10,
    topk_policy: str = "mutual",
    mz_tol: float = 0.10,
    precursor_window: float = 200.0,
    max_component_size: int = 0,
    top_n: int = 70,
    rel_threshold: float = 0.01,
) -> dict:
    """Build the V7 feature-consensus, multi-evidence molecular network."""
    policy = str(topk_policy or "either").lower()
    if policy not in {"mutual", "either", "off"}:
        raise ValueError("topk_policy must be mutual, either, or off")
    params = {
        "engine": "v7_evidence_network",
        "node_mode": "lcms_feature_consensus_with_precursor_event_fallback",
        "cos_thresh": float(cos_thresh),
        "min_peaks": int(min_peaks),
        "topk": int(topk),
        "topk_policy": policy,
        "mz_tol": float(mz_tol),
        "precursor_window": float(precursor_window),
        "max_component_size": int(max_component_size),
        "top_n": int(top_n),
        "rel_threshold": float(rel_threshold),
    }
    nodes, preparation = prepare_feature_nodes(
        results,
        top_n=top_n,
        rel_threshold=rel_threshold,
        consensus_tolerance_da=max(0.01, min(float(mz_tol), 0.10)),
    )
    node_count = len(nodes)
    if node_count < 2:
        if node_count == 1:
            nodes[0].update({
                "component_id": 0,
                "family_id": None,
                "degree": 0,
                "annotation_source": "library" if nodes[0].get("match") else "unknown",
                "propagated_class": None,
                "propagated_compound": None,
                "propagation_evidence": None,
            })
        stats = {
            **preparation,
            "n_spectra": len(results or []),
            "n_nodes": node_count,
            "n_edges": 0,
            "n_families": 0,
            "n_molecular_families": 0,
            "n_components_total": node_count,
            "connected_node_count": 0,
            "orphan_count": node_count,
            "singleton_count": node_count,
            "largest_family": 0,
            "n_matched": sum(bool(node.get("match")) for node in nodes),
            "n_propagated": 0,
            "n_unknown_pure": sum(not node.get("match") for node in nodes),
            "raw_peak_nodes": sum(
                node.get("peak_source") == "raw" for node in nodes
            ),
            "full_peak_nodes": sum(
                node.get("peak_source") in {"raw", "lcms_feature_consensus", "lcms_feature_representative"}
                for node in nodes
            ),
            "raw_peak_coverage_pct": round(
                100.0
                * sum(
                    node.get("peak_source")
                    in {"raw", "lcms_feature_consensus", "lcms_feature_representative"}
                    for node in nodes
                )
                / max(1, node_count), 2
            ),
            "full_peak_coverage_pct": round(
                100.0
                * sum(
                    node.get("peak_source")
                    in {"raw", "lcms_feature_consensus", "lcms_feature_representative"}
                    for node in nodes
                )
                / max(1, node_count), 2
            ),
            "warnings": [],
            "params": params,
        }
        return {
            "nodes": nodes,
            "edges": [],
            "edge_details": [],
            "families": [],
            "stats": stats,
            "params": params,
        }

    candidate_floor = max(0.18, min(0.48, float(cos_thresh) * 0.60))
    overlap_floor = max(2, int(min_peaks) - 3)
    candidate_pairs = _sparse_candidate_pairs(
        nodes,
        mz_tolerance_da=float(mz_tol),
        precursor_window=float(precursor_window),
        score_floor=candidate_floor,
        overlap_floor=overlap_floor,
    )
    exact_edges: list[dict] = []
    for source_pos, target_pos in sorted(candidate_pairs):
        detail = _exact_edge(
            nodes[source_pos],
            nodes[target_pos],
            cosine_threshold=float(cos_thresh),
            min_matched_peaks=int(min_peaks),
            mz_tolerance_da=float(mz_tol),
        )
        if detail:
            detail.update({
                "source_pos": source_pos,
                "target_pos": target_pos,
            })
            exact_edges.append(detail)

    topk_edges = _topk_edge_filter(
        exact_edges,
        topk=int(topk),
        policy=policy,
    )
    capped_edges = _component_cap(
        topk_edges,
        node_count=node_count,
        max_component_size=int(max_component_size),
    )
    iimn_edges = _iimn_edges(nodes)
    edges = _merge_edge_layers([*capped_edges, *iimn_edges])
    components = _components(node_count, edges)
    _annotate_nodes(nodes, edges, components)

    family_summaries: list[dict] = []
    for members in components:
        if len(members) < 2:
            continue
        family_id = nodes[members[0]]["family_id"]
        matched = [nodes[position] for position in members if nodes[position].get("match")]
        classes = [
            node.get("classification") for node in matched
            if node.get("classification")
            and str(node["classification"]).lower() != "unknown"
        ]
        best = max(
            matched,
            key=lambda node: _number(node.get("final_score")) or 0.0,
        ) if matched else None
        precursor_values = [nodes[position]["precursor_mz"] for position in members]
        family_summaries.append({
            "id": family_id,
            "size": len(members),
            "spectra_count": sum(nodes[position]["spectrum_count"] for position in members),
            "matched_count": len(matched),
            "consensus_class": (
                Counter(classes).most_common(1)[0][0] if classes else None
            ),
            "top_compound": best.get("compound_name") if best else None,
            "top_score": best.get("final_score") if best else None,
            "precursor_min": round(min(precursor_values), 4),
            "precursor_max": round(max(precursor_values), 4),
            "edge_types": sorted({
                edge["edge_type"] for edge in edges
                if edge["source_pos"] in members and edge["target_pos"] in members
            }),
        })
    family_summaries.sort(key=lambda family: (-family["size"], family["id"]))

    public_details: list[dict] = []
    public_edges: list[tuple[int, int, float, int]] = []
    for edge in edges:
        source = nodes[edge["source_pos"]]["idx"]
        target = nodes[edge["target_pos"]]["idx"]
        public = {
            key: value for key, value in edge.items()
            if key not in {"source_pos", "target_pos"}
        }
        public.update({"source": source, "target": target})
        public_details.append(public)
        public_edges.append((
            source,
            target,
            float(edge["score"]),
            int(edge["matched_peaks"]),
        ))

    connected_nodes = sum(len(component) for component in components if len(component) >= 2)
    singletons = sum(len(component) == 1 for component in components)
    largest_connected = max(
        (len(component) for component in components if len(component) >= 2),
        default=0,
    )
    raw_nodes = sum(node.get("peak_source") == "raw" for node in nodes)
    full_peak_nodes = sum(
        node.get("peak_source") in {"raw", "lcms_feature_consensus", "lcms_feature_representative"}
        for node in nodes
    )
    fallback_nodes = sum(
        node.get("peak_source") in {"search_fallback", "mixed"} for node in nodes
    )
    warnings: list[dict] = []
    if fallback_nodes:
        warnings.append({
            "code": "raw_peak_fallback",
            "message": (
                f"{fallback_nodes} node(s) use the reduced search-peak view; "
                "rerun legacy tasks in V7 for complete network evidence."
            ),
        })
    if not iimn_edges:
        warnings.append({
            "code": "iimn_not_available",
            "message": (
                "No validated MS1 ion-identity edges were available; spectral "
                "families are shown without adduct collapsing."
            ),
        })
    if not preparation.get("lcms_feature_nodes"):
        warnings.append({
            "code": "lcms_feature_fallback",
            "message": (
                "No MS2 spectrum was linked to an OpenMS chromatographic "
                "feature; nodes are explicit MS2 precursor-event fallbacks."
            ),
        })
    if preparation.get("insufficient_fragment_nodes"):
        warnings.append({
            "code": "insufficient_fragment_singletons",
            "message": (
                f"{preparation['insufficient_fragment_nodes']} nodes have fewer "
                "than three usable fragment peaks. They remain visible as "
                "singletons but cannot form spectral-similarity edges."
            ),
        })
    singleton_fraction = singletons / max(1, node_count)
    if node_count >= 10 and singleton_fraction >= 0.75:
        warnings.append({
            "code": "high_singleton_fraction",
            "message": (
                f"{singleton_fraction:.1%} of nodes are singletons at the "
                "current evidence thresholds. Review centroiding, feature "
                "linkage, fragment tolerance and acquisition depth before "
                "loosening the cosine gate."
            ),
        })
    if (
        int(max_component_size) == 0
        and largest_connected > max(100, int(math.ceil(node_count * 0.20)))
    ):
        warnings.append({
            "code": "possible_network_hairball",
            "message": (
                f"The largest family contains {largest_connected}/{node_count} "
                "nodes. Prefer mutual Top-K or stricter spectral evidence "
                "before applying an artificial component-size cap."
            ),
        })

    edge_type_counts = Counter(edge["edge_type"] for edge in edges)
    n_matched = sum(node.get("annotation_source") == "library" for node in nodes)
    n_propagated = sum(node.get("annotation_source") == "propagated" for node in nodes)
    n_ion_identity = sum(
        node.get("annotation_source") == "ion_identity" for node in nodes
    )
    stats = {
        **preparation,
        "n_spectra": len(results or []),
        "n_nodes": node_count,
        "n_edges": len(edges),
        "n_families": len(family_summaries),
        "n_molecular_families": len(family_summaries),
        "n_components_total": len(components),
        "connected_node_count": connected_nodes,
        "orphan_count": singletons,
        "singleton_count": singletons,
        "largest_family": largest_connected,
        "n_matched": n_matched,
        "n_propagated": n_propagated,
        "n_ion_identity_linked": n_ion_identity,
        "n_unknown_pure": node_count - n_matched - n_propagated - n_ion_identity,
        "raw_peak_nodes": raw_nodes,
        "full_peak_nodes": full_peak_nodes,
        "raw_peak_fallback_nodes": fallback_nodes,
        # Backward-compatible key now reflects any full measured/consensus
        # view, not only a single-scan raw list.
        "raw_peak_coverage_pct": round(
            100.0 * full_peak_nodes / max(1, node_count), 2
        ),
        "full_peak_coverage_pct": round(
            100.0 * full_peak_nodes / max(1, node_count), 2
        ),
        "connected_node_fraction_pct": round(
            100.0 * connected_nodes / max(1, node_count), 2
        ),
        "singleton_fraction_pct": round(
            100.0 * singletons / max(1, node_count), 2
        ),
        "candidate_pairs": len(candidate_pairs),
        "n_edges_threshold": len(exact_edges),
        "n_edges_after_topk": len(topk_edges),
        "n_edges_removed_topk": len(exact_edges) - len(topk_edges),
        "n_edges_removed_component_cap": len(topk_edges) - len(capped_edges),
        "component_cap_applied": bool(max_component_size > 0),
        "spectral_edges": edge_type_counts["spectral"],
        "analog_edges": edge_type_counts["analog"],
        "neutral_loss_edges": edge_type_counts["neutral_loss"],
        "ion_identity_edges": sum(
            "ion_identity" in edge.get("evidence_layers", []) for edge in edges
        ),
        "topk_policy": policy,
        "warnings": warnings,
        "params": params,
    }
    return {
        "nodes": nodes,
        "edges": public_edges,
        "edge_details": public_details,
        "families": family_summaries,
        "stats": stats,
        "params": params,
    }
