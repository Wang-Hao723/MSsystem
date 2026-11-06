"""GNPS 风格分子网络 (Feature-Based Molecular Networking)

把任务结果里的 N 个 query 谱图按 MS/MS 谱图相似度聚类成"分子家族 (molecular family)"，
不依赖标准品库匹配 —— 这正是解决"匹配率 ~39%、仍有 ~61% 未匹配"问题的核心：
未匹配的未知物若与已匹配节点谱图相似，会被划入同一家族，并把已知注释**传播**给家族邻居。

核心算法 (GNPS / FBMN 标准):
  1. 两两 Modified-Cosine 相似度 (这里用分箱向量化快速近似, 见下)
  2. 建边: cosine >= cos_thresh 且 匹配峰数 >= min_peaks
  3. 拓扑过滤: 仅保留双方互为 top-k 的边 (GNPS 默认 10)
  4. 大家族过滤: 优先保留高分边，把连通分量限制在指定节点数以内
  5. 连通分量 = 分子家族 (union-find)
  6. 注释传播: 家族内已知化合物的分类/名称传播到未知节点

性能方案 (N×N 余弦瓶颈):
  - 用 scipy.sparse 把每条谱图编码成 m/z 分箱强度向量 (bin_width=mz_tol),
    行 L2 归一化后, 余弦 = S @ S.T。
  - 分块 (block) 计算, 内存占用 O(block×N), 可扩展到上万谱图。
  - precursor m/z ± window 预过滤 (结构类似物亚结构差 < window Da)。

注意:
  - 分箱余弦是 GNPS Modified Cosine 的快速近似 (不显式做 precursor-shift 对齐),
    用于网络相似度排序足够; 高质量候选边可在调用端用 scores.cosine_score 精算。
  - 每个 spectrum_id = 一个节点 (这里不做结果去重, 去重针对的是重复库匹配, 与网络无关)。
"""
from __future__ import annotations

import math
import csv
import io
import json
from collections import Counter, defaultdict
from typing import Any, Optional

import numpy as np
import scipy.sparse as sp

try:  # GraphML 导出用; 不可用时降级为纯文本 edge list
    import networkx as nx
    _HAS_NX = True
except Exception:  # pragma: no cover
    _HAS_NX = False

from utils.spectrum import clean_peaks


# ──────────────────────────────────────────────────────────────────────────
# 数据准备: 从 result 列表抽取并清洗谱图
# ──────────────────────────────────────────────────────────────────────────

def prepare_nodes(results: list[dict],
                  top_n: int = 70,
                  rel_threshold: float = 0.01) -> list[dict]:
    """从任务结果列表构建网络节点。

    Args:
        results: _get_task_results() 返回的 result 列表, 每条含
                 query_peaks / precursor_mz / spectrum_id / retention_time / match
        top_n:    谱图清洗保留的最大峰数 (网络用比匹配的 top-15 更多峰)
        rel_threshold: 相对强度阈值, 过滤噪声峰

    Returns:
        节点字典列表, 每个含:
            idx, spectrum_id, precursor_mz, retention_time, peaks (清洗后),
            peaks_count, match (dict|None), compound_name, classification, confidence
        峰不足的谱图被跳过。
    """
    nodes: list[dict] = []
    for idx, r in enumerate(results):
        # Molecular networking needs the richer denoised MS/MS view.  The
        # library search intentionally stores only its Top-15 scoring view in
        # ``query_peaks``; using that here silently discarded low-intensity
        # diagnostic ions even though ``raw_peaks`` is saved for audit/reporting
        # for every new task.  Legacy tasks fall back to query_peaks.
        peaks = r.get("raw_peaks") or r.get("query_peaks") or []
        pmz = r.get("precursor_mz")
        if not peaks or pmz is None or pmz <= 0:
            continue
        cleaned = clean_peaks(peaks, pmz, top_n=top_n, rel_threshold=rel_threshold)
        if len(cleaned) < 3:  # 峰太少, 无法可靠比较
            continue

        m = r.get("match")
        is_match = isinstance(m, dict) and bool(m) and m.get("compound_name")
        nodes.append({
            "idx": idx,  # 原结果列表索引, 便于回溯
            "spectrum_id": r.get("spectrum_id"),
            "precursor_mz": float(pmz),
            "retention_time": r.get("retention_time"),
            "peaks": cleaned,
            "peaks_count": len(cleaned),
            "match": m if is_match else None,
            "compound_name": (m.get("compound_name") if is_match else None),
            "classification": (m.get("classification") if is_match else None),
            "confidence": (m.get("confidence") if is_match else None),
            "evidence_tier": (m.get("evidence_tier") if is_match else None),
            "final_score": (m.get("final_score") if is_match else None),
            "cosine": (m.get("cosine") if is_match else None),
            "entropy": (m.get("entropy") if is_match else None),
            "mass_ppm": (m.get("mass_ppm") if is_match else None),
            "smiles": (m.get("smiles") if is_match else None),
            "inchikey": (m.get("inchikey") if is_match else None),
            "formula": (m.get("formula") if is_match else None),
            "adduct": (m.get("adduct") if is_match else None),
            "dataset_source": (m.get("dataset_source") if is_match else None),
            "source_tier": (m.get("source_tier") if is_match else None),
            "ion_mode": r.get("ion_mode"),
            "precursor_charge": r.get("precursor_charge"),
            "collision_energy": r.get("collision_energy"),
        })
    return nodes


# ──────────────────────────────────────────────────────────────────────────
# 向量化: 谱图 → 稀疏分箱强度矩阵
# ──────────────────────────────────────────────────────────────────────────

def _bin_index(mz: float, bin_width: float) -> int:
    """m/z → 分箱索引 (同一 bin 内的峰视为可匹配)。"""
    return int(math.floor(mz / bin_width))


def build_sparse_vectors(nodes: list[dict],
                         mz_tol: float = 0.1) -> tuple[sp.csr_matrix, sp.csr_matrix, np.ndarray]:
    """把节点谱图编码成稀疏分箱矩阵。

    - 同一 bin 内多个峰取最大强度 (GNPS-style)。
    - S 的每行 L2 归一化 → 行点积 = 余弦相似度。
    - Sbool 为 0/1 占用矩阵 → 行点积 = 共享 bin 数 (匹配峰数近似)。

    Returns:
        S:      (N, B) csr, 行归一化强度
        Sbool:  (N, B) csr, 0/1 峰占用
        pmz:    (N,)   precursor m/z 数组
    """
    n = len(nodes)
    rows, cols, vals, bcols = [], [], [], []
    for i, node in enumerate(nodes):
        bin_max: dict[int, float] = {}
        for mz, inten in node["peaks"]:
            b = _bin_index(mz, mz_tol)
            if inten > bin_max.get(b, 0.0):
                bin_max[b] = float(inten)
        if not bin_max:
            continue
        for b, v in bin_max.items():
            rows.append(i)
            cols.append(b)
            vals.append(v)
            bcols.append(b)

    if not rows:
        empty = sp.csr_matrix((n, 1), dtype=np.float32)
        return empty, empty.astype(bool).tocsr(), np.zeros(n, dtype=np.float64)

    data = np.asarray(vals, dtype=np.float32)
    row_idx = np.asarray(rows, dtype=np.int64)
    col_idx = np.asarray(cols, dtype=np.int64)
    M = sp.coo_matrix((data, (row_idx, col_idx)), shape=(n, int(col_idx.max()) + 1)).tocsr()

    # 行 L2 归一化 (余弦)
    sq = np.asarray(M.multiply(M).sum(axis=1)).ravel()
    norm = np.sqrt(sq)
    norm[norm == 0] = 1.0
    S = sp.diags(1.0 / norm) @ M
    S = S.tocsr()

    Sbool = (M != 0).astype(np.float32).tocsr()

    pmz = np.array([node["precursor_mz"] for node in nodes], dtype=np.float64)
    return S, Sbool, pmz


# ──────────────────────────────────────────────────────────────────────────
# 两两相似度 (分块, 内存可控)
# ──────────────────────────────────────────────────────────────────────────

def pairwise_scores(S: sp.csr_matrix, Sbool: sp.csr_matrix, pmz: np.ndarray,
                    cos_thresh: float = 0.7, min_peaks: int = 6,
                    precursor_window: float = 200.0,
                    block: int = 512) -> list[tuple[int, int, float, int]]:
    """计算满足阈值的候选边 (i<j)。

    分块计算 S[block] @ S.T, 只保留:
        cosine >= cos_thresh AND matched_peaks >= min_peaks
        AND |precursor_i - precursor_j| <= precursor_window AND i<j

    Returns:
        [(i, j, cosine, matched_peaks), ...]
    """
    n = S.shape[0]
    edges: list[tuple[int, int, float, int]] = []
    if n < 2:
        return edges

    S_T = S.T.tocsr()          # (B, N)
    B_T = Sbool.T.tocsr()      # (B, N)
    order = np.argsort(pmz)    # 按 precursor 排序, 用于窗口预过滤
    pmz_sorted = pmz[order]

    for start in range(0, n, block):
        end = min(n, start + block)
        sub = S[start:end]            # (b, N) via S[block] @ S.T
        cos_block = (sub @ S_T).toarray()          # (b, N) float32
        ov_block = (Sbool[start:end] @ B_T).toarray().astype(np.int32)  # (b, N)

        for ki in range(end - start):
            i = start + ki
            ci = cos_block[ki]
            oi = ov_block[ki]
            # precursor 窗口: 二分定位 |pmz_i - pmz_j| <= window 的 j 范围
            # 二分查找必须作用于已排序数组；旧实现误在原始 pmz 上查找，
            # 当输入扫描并非按前体 m/z 排序时会漏边或纳入窗口外候选。
            lo = np.searchsorted(pmz_sorted, pmz[i] - precursor_window, side="left")
            hi = np.searchsorted(pmz_sorted, pmz[i] + precursor_window, side="right")
            cand = order[lo:hi]
            cand = cand[cand > i]                 # i<j
            if cand.size == 0:
                continue
            cvals = ci[cand]
            ovals = oi[cand]
            mask = (cvals >= cos_thresh) & (ovals >= min_peaks)
            for j, cos, matched in zip(cand[mask], cvals[mask], ovals[mask]):
                edges.append((int(i), int(j), float(cos), int(matched)))
    return edges


# ──────────────────────────────────────────────────────────────────────────
# 拓扑过滤: 每节点 top-k 邻居 (GNPS: 任一端点把它排进 top-k 则保留)
# ──────────────────────────────────────────────────────────────────────────

def topk_filter(edges: list[tuple[int, int, float, int]],
                n_nodes: int, topk: int = 10) -> list[tuple[int, int, float, int]]:
    """仅保留双方互为 top-k 的边（GNPS mutual Top-K 规则）。"""
    if topk <= 0 or not edges:
        return edges

    by_node: dict[int, list[tuple[float, int, tuple]]] = defaultdict(list)
    for e in edges:
        i, j, cos, matched = e
        by_node[i].append((cos, j, e))
        by_node[j].append((cos, i, e))

    top_neighbors: dict[int, set[int]] = {}
    for node, lst in by_node.items():
        lst.sort(key=lambda t: (t[0], t[2][3]), reverse=True)
        top_neighbors[node] = {neighbor for _, neighbor, _ in lst[:topk]}
    return [
        edge for edge in edges
        if edge[1] in top_neighbors.get(edge[0], set())
        and edge[0] in top_neighbors.get(edge[1], set())
    ]


# ──────────────────────────────────────────────────────────────────────────
# 连通分量 (union-find) = 分子家族
# ──────────────────────────────────────────────────────────────────────────

class _DSU:
    def __init__(self, n: int):
        self.p = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if self.size[ra] < self.size[rb]:
                ra, rb = rb, ra
            self.p[rb] = ra
            self.size[ra] += self.size[rb]


def component_size_filter(
    edges: list[tuple[int, int, float, int]],
    n_nodes: int,
    max_component_size: int = 100,
) -> list[tuple[int, int, float, int]]:
    """Break hairball networks while retaining the strongest alignments first.

    Edges are processed from high to low similarity. An edge that would merge
    two components beyond ``max_component_size`` is omitted; edges already
    inside an accepted component are retained. This is a deterministic,
    efficient equivalent of iteratively removing weak bridges from oversized
    components.
    """
    if max_component_size <= 0 or not edges:
        return edges
    dsu = _DSU(n_nodes)
    kept: list[tuple[int, int, float, int]] = []
    for edge in sorted(edges, key=lambda e: (e[2], e[3]), reverse=True):
        i, j, *_ = edge
        root_i, root_j = dsu.find(i), dsu.find(j)
        if root_i == root_j:
            kept.append(edge)
            continue
        if dsu.size[root_i] + dsu.size[root_j] <= max_component_size:
            dsu.union(root_i, root_j)
            kept.append(edge)
    # Stable public output independent of score-sort used during filtering.
    kept.sort(key=lambda e: (e[0], e[1]))
    return kept


def connected_components(n_nodes: int,
                         edges: list[tuple[int, int, float, int]]
                         ) -> tuple[np.ndarray, dict[int, list[int]]]:
    """返回 (family_id_per_node, {family_id: [node_idx, ...]})。"""
    dsu = _DSU(n_nodes)
    for i, j, *_ in edges:
        dsu.union(i, j)
    fam = np.array([dsu.find(i) for i in range(n_nodes)], dtype=np.int64)
    # 重映射为 0..K-1
    uniq, inv = np.unique(fam, return_inverse=True)
    families: dict[int, list[int]] = defaultdict(list)
    for node_idx, fid in enumerate(inv.tolist()):
        families[fid].append(node_idx)
    return inv, families


# ──────────────────────────────────────────────────────────────────────────
# 注释传播
# ──────────────────────────────────────────────────────────────────────────

def propagate_annotations(nodes: list[dict],
                          family_id: np.ndarray,
                          families: dict[int, list[int]]) -> None:
    """家族内传播: 已知节点 (库匹配) 的分类/名称传给同家族未知节点。

    每个 node 新增字段:
        family_id           : 所属家族 id
        annotation_source   : library | propagated | unknown
        propagated_class    : 家族共识分类 (仅 propagated 节点)
        propagated_compound : 家族最高分化合物名 (仅 propagated 节点)
    """
    # 先算每个家族的共识注释
    family_meta: dict[int, dict] = {}
    for fid, members in families.items():
        matched = [nodes[m] for m in members if nodes[m]["match"]]
        classes = [n["classification"] for n in matched
                   if n["classification"] and str(n["classification"]).lower() != "unknown"]
        consensus_class = (Counter(classes).most_common(1)[0][0]
                           if classes else None)
        # 最高 final_score 的化合物作为家族代表
        best = max(matched, key=lambda n: (n.get("final_score") or 0.0)) if matched else None
        family_meta[fid] = {
            "consensus_class": consensus_class,
            "rep_compound": best["compound_name"] if best else None,
            "rep_score": best.get("final_score") if best else None,
        }

    for i, node in enumerate(nodes):
        fid = int(family_id[i])
        node["family_id"] = fid
        meta = family_meta[fid]
        if node["match"]:
            node["annotation_source"] = "library"
            node["propagated_class"] = None
            node["propagated_compound"] = None
        elif meta["consensus_class"] or meta["rep_compound"]:
            node["annotation_source"] = "propagated"
            node["propagated_class"] = meta["consensus_class"]
            node["propagated_compound"] = meta["rep_compound"]
        else:
            node["annotation_source"] = "unknown"
            node["propagated_class"] = None
            node["propagated_compound"] = None


# ──────────────────────────────────────────────────────────────────────────
# 主入口
# ──────────────────────────────────────────────────────────────────────────

def build_network(results: list[dict],
                  cos_thresh: float = 0.7,
                  min_peaks: int = 6,
                  topk: int = 10,
                  mz_tol: float = 0.1,
                  precursor_window: float = 200.0,
                  max_component_size: int = 0,
                  top_n: int = 70,
                  rel_threshold: float = 0.01,
                  topk_policy: str = "mutual") -> dict:
    """构建 V7 多证据、特征共识分子网络。

    Args:
        results: 任务结果列表 (_get_task_results 返回)
        cos_thresh:      精确建边余弦阈值
        min_peaks:       建边最少精确匹配峰数
        topk:            每节点保留的最相似邻居数
        mz_tol:          精确峰匹配容差 (Da)
        precursor_window: 仅比较 precursor m/z 差 ≤ 此值的谱图 (Da)
        max_component_size: 单个分子家族最大节点数；0 表示不限制
        top_n:           谱图清洗保留峰数
        rel_threshold:   相对强度阈值
        topk_policy:     either / mutual / off；默认 mutual，遵循保守的
                         双向近邻证据，不用单向近邻制造巨型组件

    Returns:
        {nodes, edges, edge_details, families, stats, params}
    """
    from matching.evidence_network import build_evidence_network

    return build_evidence_network(
        results,
        cos_thresh=cos_thresh,
        min_peaks=min_peaks,
        topk=topk,
        topk_policy=topk_policy,
        mz_tol=mz_tol,
        precursor_window=precursor_window,
        max_component_size=max_component_size,
        top_n=top_n,
        rel_threshold=rel_threshold,
    )


# ──────────────────────────────────────────────────────────────────────────
# 序列化 / 导出
# ──────────────────────────────────────────────────────────────────────────

def _safe(v: Any) -> Any:
    """转成 JSON/GraphML 安全的标量。"""
    if v is None:
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return None if math.isnan(float(v)) else round(float(v), 4)
    return v


def node_public(node: dict) -> dict:
    """节点去掉 peaks 后的公开字段 (peaks 体积大, 单独按需给)。"""
    peak_evidence = node.get("peak_evidence")
    annotation_result_index = (
        peak_evidence.get("annotation_result_index")
        if isinstance(peak_evidence, dict)
        else None
    )
    peak_spectrum_id = (
        peak_evidence.get("spectrum_id")
        if isinstance(peak_evidence, dict)
        else node.get("spectrum_id")
    )
    peak_artifact_id = (
        peak_evidence.get("artifact_id")
        if isinstance(peak_evidence, dict)
        else node.get("lcms_consensus_spectrum_id")
    )
    return {
        "id": node["idx"],
        "feature_id": node.get("feature_id"),
        "feature_source": node.get("feature_source"),
        "lcms_feature_id": node.get("lcms_feature_id"),
        "lcms_feature_intensity": _safe(node.get("lcms_feature_intensity")),
        "lcms_feature_charge": _safe(node.get("lcms_feature_charge")),
        "lcms_isotope_trace_count": _safe(
            node.get("lcms_isotope_trace_count")
        ),
        "lcms_adduct": node.get("lcms_adduct"),
        "lcms_neutral_mass": _safe(node.get("lcms_neutral_mass")),
        "lcms_adduct_group_id": node.get("lcms_adduct_group_id"),
        "lcms_consensus_spectrum_id": peak_artifact_id,
        "annotation_lcms_consensus_spectrum_id": node.get(
            "lcms_consensus_spectrum_id"
        ),
        "member_result_indices": node.get("member_result_indices") or [node["idx"]],
        "spectrum_count": int(node.get("spectrum_count") or 1),
        "spectrum_id": _safe(peak_spectrum_id),
        "peak_spectrum_id": _safe(peak_spectrum_id),
        "peak_result_index": (
            peak_evidence.get("result_index")
            if isinstance(peak_evidence, dict)
            else node.get("idx")
        ),
        "annotation_spectrum_id": _safe(node.get("spectrum_id")),
        "precursor_mz": round(float(
            peak_evidence.get("precursor_mz")
            if isinstance(peak_evidence, dict)
            and peak_evidence.get("precursor_mz") is not None
            else node["precursor_mz"]
        ), 4),
        "retention_time": _safe(
            peak_evidence.get("retention_time")
            if isinstance(peak_evidence, dict)
            else node.get("retention_time")
        ),
        "feature_precursor_mz": round(float(node["precursor_mz"]), 4),
        "feature_retention_time": _safe(node.get("retention_time")),
        "peaks_count": node["peaks_count"],
        "peak_source": node.get("peak_source", "unknown"),
        "peak_evidence": peak_evidence,
        "annotation_result_index": annotation_result_index,
        "network_eligibility": node.get(
            "network_eligibility", "similarity_eligible"
        ),
        "raw_member_count": int(node.get("raw_member_count") or 0),
        "degree": node.get("degree", 0),
        "matched": node["match"] is not None,
        "compound_name": node.get("compound_name"),
        "classification": node.get("classification"),
        "confidence": node.get("confidence"),
        "evidence_tier": node.get("evidence_tier"),
        "final_score": _safe(node.get("final_score")),
        "cosine": _safe(node.get("cosine")),
        "entropy": _safe(node.get("entropy")),
        "mass_ppm": _safe(node.get("mass_ppm")),
        "smiles": node.get("smiles"),
        "inchikey": node.get("inchikey"),
        "formula": node.get("formula"),
        "adduct": node.get("adduct"),
        "dataset_source": node.get("dataset_source"),
        "source_tier": node.get("source_tier"),
        "ion_mode": (
            peak_evidence.get("ion_mode")
            if isinstance(peak_evidence, dict)
            else node.get("ion_mode")
        ),
        "precursor_charge": _safe(
            peak_evidence.get("precursor_charge")
            if isinstance(peak_evidence, dict)
            else node.get("precursor_charge")
        ),
        "collision_energy": _safe(
            peak_evidence.get("collision_energy")
            if isinstance(peak_evidence, dict)
            else node.get("collision_energy")
        ),
        "collision_energy_unit": (
            peak_evidence.get("collision_energy_unit")
            if isinstance(peak_evidence, dict)
            else None
        ),
        "collision_energy_mode": (
            peak_evidence.get("collision_energy_mode")
            if isinstance(peak_evidence, dict)
            else None
        ),
        "activation_method": (
            peak_evidence.get("activation_method")
            if isinstance(peak_evidence, dict)
            else None
        ),
        "instrument_configuration_id": (
            peak_evidence.get("instrument_configuration_id")
            if isinstance(peak_evidence, dict)
            else None
        ),
        "collision_energies": node.get("collision_energies") or [],
        "ion_identity_status": node.get("ion_identity_status"),
        "ion_identity_network_id": node.get("ion_identity_network_id"),
        "ion_identity": node.get("ion_identity"),
        "annotation_source": node.get("annotation_source", "unknown"),
        "propagated_class": node.get("propagated_class"),
        "propagated_compound": node.get("propagated_compound"),
        "propagation_evidence": node.get("propagation_evidence"),
        "family_id": node.get("family_id"),
        "component_id": node.get("component_id"),
    }


def edge_public(edge: tuple | dict) -> dict:
    if isinstance(edge, dict):
        return {
            "source": int(edge["source"]),
            "target": int(edge["target"]),
            "cosine": round(float(edge.get("cosine", edge.get("score", 0.0))), 4),
            "score": round(float(edge.get("score", edge.get("cosine", 0.0))), 4),
            "matched_peaks": int(edge.get("matched_peaks") or 0),
            "edge_type": edge.get("edge_type", "spectral"),
            "evidence_layers": edge.get("evidence_layers") or [
                edge.get("edge_type", "spectral")
            ],
            "direct_cosine": round(float(edge.get("direct_cosine") or 0.0), 4),
            "modified_cosine": round(float(edge.get("modified_cosine") or 0.0), 4),
            "neutral_loss_cosine": round(
                float(edge.get("neutral_loss_cosine") or 0.0), 4
            ),
            "precursor_delta_da": round(
                float(edge.get("precursor_delta_da") or 0.0), 4
            ),
            "ion_identity_network_id": edge.get("ion_identity_network_id"),
        }
    i, j, cos, matched = edge
    return {"source": i, "target": j, "cosine": round(cos, 4),
            "score": round(cos, 4), "matched_peaks": int(matched),
            "edge_type": "spectral", "evidence_layers": ["spectral"]}


def to_nodes_csv(network: dict) -> str:
    """Export a Cytoscape/Gephi-friendly node table with stable public IDs."""
    fields = [
        "id", "spectrum_id", "precursor_mz",
        "peak_spectrum_id", "peak_result_index",
        "annotation_spectrum_id", "annotation_result_index",
        "retention_time", "feature_precursor_mz",
        "feature_retention_time", "feature_id", "feature_source",
        "lcms_feature_id", "lcms_feature_intensity", "lcms_feature_charge",
        "lcms_isotope_trace_count", "lcms_adduct", "lcms_neutral_mass",
        "lcms_adduct_group_id", "lcms_consensus_spectrum_id",
        "annotation_lcms_consensus_spectrum_id", "spectrum_count", "ion_mode",
        "precursor_charge", "collision_energy", "collision_energy_unit",
        "collision_energy_mode", "activation_method",
        "instrument_configuration_id", "peaks_count", "family_id",
        "component_id", "peak_source", "network_eligibility",
        "peak_evidence",
        "raw_member_count", "degree",
        "annotation_source", "matched", "compound_name", "formula",
        "adduct", "inchikey", "smiles", "classification", "confidence",
        "evidence_tier", "cosine", "entropy", "final_score", "mass_ppm",
        "dataset_source", "source_tier", "propagated_class", "propagated_compound",
        "propagation_evidence", "ion_identity_status",
        "ion_identity_network_id", "ion_identity",
    ]
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    for node in network.get("nodes") or []:
        row = node_public(node)
        if isinstance(row.get("peak_evidence"), dict):
            row["peak_evidence"] = json.dumps(
                row["peak_evidence"], ensure_ascii=False, sort_keys=True,
                separators=(",", ":"),
            )
        writer.writerow(row)
    return buf.getvalue()


def to_edges_csv(network: dict) -> str:
    """Export an edge table enriched with node/family identifiers."""
    fields = [
        "source", "target", "cosine", "matched_peaks", "family_id", "score",
        "edge_type",
        "evidence_layers", "direct_cosine", "modified_cosine",
        "neutral_loss_cosine", "ion_identity_network_id",
        "source_spectrum_id", "target_spectrum_id", "source_precursor_mz",
        "target_precursor_mz", "precursor_delta_da",
    ]
    nodes = {node["idx"]: node for node in network.get("nodes") or []}
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=fields)
    writer.writeheader()
    source_edges = network.get("edge_details") or network.get("edges") or []
    for edge in source_edges:
        row = edge_public(edge)
        source = nodes.get(row["source"], {})
        target = nodes.get(row["target"], {})
        source_mz = _safe(source.get("precursor_mz"))
        target_mz = _safe(target.get("precursor_mz"))
        writer.writerow({
            **row,
            "family_id": source.get("family_id"),
            "source_spectrum_id": _safe(source.get("spectrum_id")),
            "target_spectrum_id": _safe(target.get("spectrum_id")),
            "source_precursor_mz": source_mz,
            "target_precursor_mz": target_mz,
            "precursor_delta_da": (
                row.get("precursor_delta_da")
                if row.get("precursor_delta_da") is not None else (
                    round(abs(float(source_mz) - float(target_mz)), 4)
                    if source_mz is not None and target_mz is not None else None
                )
            ),
            "evidence_layers": "|".join(row.get("evidence_layers") or []),
        })
    return buf.getvalue()


def to_cytoscape_elements(network: dict,
                          include_peaks: bool = False) -> dict:
    """生成 cytoscape.js elements {nodes, edges}。"""
    nodes = []
    for node in network["nodes"]:
        pub = node_public(node)
        if include_peaks:
            pub["peaks"] = node["peaks"]
        nodes.append({"data": pub})
    source_edges = network.get("edge_details") or network.get("edges") or []
    edges = [{"data": dict(edge_public(e), **{"id": f"e{ei}"})}
             for ei, e in enumerate(source_edges)]
    return {"nodes": nodes, "edges": edges}


def to_graphml(network: dict) -> Optional[str]:
    """导出 Cytoscape 兼容的 GraphML 字符串。networkx 不可用时返回 None。"""
    if not _HAS_NX:
        return None
    G = nx.Graph()
    for node in network["nodes"]:
        pub = node_public(node)
        # GraphML 属性需简单类型
        attrs = {
            k: (
                "" if v is None
                else "|".join(map(str, v)) if isinstance(v, (list, tuple, set))
                else json.dumps(
                    v, ensure_ascii=False, sort_keys=True, separators=(",", ":")
                ) if isinstance(v, dict)
                else v
            )
            for k, v in pub.items()
            if v is not None or k in
            ("compound_name", "classification", "annotation_source")
        }
        G.add_node(pub["id"], **attrs)
    source_edges = network.get("edge_details") or network.get("edges") or []
    for ei, edge in enumerate(source_edges):
        public = edge_public(edge)
        attrs = {
            key: ("|".join(value) if isinstance(value, list) else value)
            for key, value in public.items()
            if key not in {"source", "target"} and value is not None
        }
        G.add_edge(
            public["source"],
            public["target"],
            id=f"e{ei}",
            **attrs,
        )
    import io
    buf = io.BytesIO()
    nx.write_graphml(G, buf)
    return buf.getvalue().decode("utf-8")
