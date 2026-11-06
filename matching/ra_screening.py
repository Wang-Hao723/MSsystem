"""类风湿关节炎(RA)靶向筛查模块.

对"已完成分析任务的命中结果 + 分子网络"做 RA 领域解读, 不改核心检索流程.
三档 RA 相关性判定:
  direct     命中化合物(按 InChIKey14→化合物名)在 RA 参考库内 → 已证实抗 RA 活性
  propagated 节点所在分子家族含 RA 直接命中节点 → 候选结构类似物
  class      命中的化合物分类属 RA 高相关类 → RA 相关化学类别

参考库 data/ra_reference_compounds.csv 是领域知识载体(药材/机制/毒性/文献).
化学身份以化合物名为主匹配键(库命中名已标准化), InChIKey14 为增强键(由 SMILES
经 RDKit 推导; 复杂萜类 SMILES 留空则仅靠名称匹配).
"""
from __future__ import annotations

import csv
import logging
import re
from collections import defaultdict
from typing import Any, Optional

from config.settings import (
    RA_NAME_MATCH_MIN_LEN,
    RA_REFERENCE_CSV,
    RA_RELEVANT_CLASSES,
)

logger = logging.getLogger(__name__)

# ── 模块级缓存 (首次加载后只读, 线程安全) ──
_PANEL_CACHE: Optional[list[dict]] = None
_PANEL_BY_IKEY14: dict[str, dict] = {}
_PANEL_BY_NAME: dict[str, dict] = {}

# 名称归一化时剥离的常见后缀/盐/立体标记
_NAME_STRIP_RE = re.compile(r"\s+", flags=re.UNICODE)
_NAME_NOISE_TOKENS = (
    "(+)", "(-)", "(±)", "(rac)", "rac-", "dl-", "l-", "d-", "trans-", "cis-",
    "hydrate", "hydrochloride", "hcl", "mesylate", "sulfate",
)


def _normalize_name(name: str) -> str:
    """化合物名归一化: 小写 + 去空白/立体标记/盐后缀."""
    if not name:
        return ""
    s = name.strip().lower()
    for tok in _NAME_NOISE_TOKENS:
        s = s.replace(tok, "")
    s = _NAME_STRIP_RE.sub(" ", s).strip()
    return s


def _inchikey14(ik: Optional[str]) -> str:
    """取 InChIKey 前 14 位(连接层, 忽略立体/质子)."""
    if not ik:
        return ""
    ik = ik.strip()
    if len(ik) < 14 or ik[14:15] not in ("-", ""):
        return ""
    return ik[:14]


def _ikey14_from_smiles(smiles: Optional[str]) -> str:
    """由 SMILES 经 RDKit 推导 InChIKey14; 失败/为空返回 ''."""
    if not smiles:
        return ""
    try:
        from rdkit import Chem  # 延迟导入, 避免无 rdkit 环境硬依赖
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        ik = Chem.MolToInchiKey(mol)
        return _inchikey14(ik)
    except Exception as exc:  # noqa: BLE001
        logger.debug("RDKit 推导 InChIKey14 失败 (%s): %s", smiles, exc)
        return ""


def load_ra_panel(path: Optional[str] = None, *, force: bool = False) -> list[dict]:
    """加载 RA 参考化合物库. 返回条目列表(每条含 name/ikey14/mass + 元数据).

    结果缓存; force=True 强制重读. 单元测试可传入自定义 path.
    """
    global _PANEL_CACHE, _PANEL_BY_IKEY14, _PANEL_BY_NAME
    if _PANEL_CACHE is not None and not force and path is None:
        return _PANEL_CACHE

    csv_path = path or str(RA_REFERENCE_CSV)
    panel: list[dict] = []
    by_ikey14: dict[str, dict] = {}
    by_name: dict[str, dict] = {}

    try:
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name_en = (row.get("name_en") or "").strip()
                name_cn = (row.get("name_cn") or "").strip()
                if not name_en and not name_cn:
                    continue
                smiles = (row.get("smiles") or "").strip()
                # InChIKey14: CSV 提供优先, 否则由 SMILES 推导
                ikey14 = _inchikey14(row.get("inchikey")) or _ikey14_from_smiles(smiles)
                try:
                    exact_mass = float(row.get("exact_mass") or 0.0)
                except (TypeError, ValueError):
                    exact_mass = 0.0
                entry = {
                    "name_cn": name_cn,
                    "name_en": name_en,
                    "formula": (row.get("formula") or "").strip(),
                    "exact_mass": exact_mass,
                    "smiles": smiles,
                    "inchikey14": ikey14,
                    "source_herb": (row.get("source_herb") or "").strip(),
                    "source_herb_en": (row.get("source_herb_en") or "").strip(),
                    "mechanism": (row.get("mechanism") or "").strip(),
                    "activity_note": (row.get("activity_note") or "").strip(),
                    "toxicity_flag": (row.get("toxicity_flag") or "").strip(),
                    "reference": (row.get("reference") or "").strip(),
                    "category": (row.get("category") or "").strip() or "tcm",
                }
                panel.append(entry)
                if ikey14:
                    by_ikey14.setdefault(ikey14, entry)
                for nm in (name_en, name_cn):
                    nrm = _normalize_name(nm)
                    if len(nrm) >= RA_NAME_MATCH_MIN_LEN:
                        by_name.setdefault(nrm, entry)
    except FileNotFoundError:
        logger.warning("RA 参考库不存在: %s", csv_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("加载 RA 参考库失败: %s", exc)

    if path is None:  # 仅缓存默认路径的加载结果
        _PANEL_CACHE = panel
        _PANEL_BY_IKEY14 = by_ikey14
        _PANEL_BY_NAME = by_name
    return panel


def match_to_panel(compound_name: Optional[str],
                   inchikey: Optional[str],
                   *,
                   by_ikey14: Optional[dict] = None,
                   by_name: Optional[dict] = None) -> Optional[dict]:
    """把一个库命中化合物匹配到 RA 参考库条目.

    优先 InChIKey14(连接层精确匹配), 回退化合物名归一化匹配. 返回条目或 None.
    """
    if by_ikey14 is None:
        by_ikey14 = _PANEL_BY_IKEY14
    if by_name is None:
        by_name = _PANEL_BY_NAME

    ik = _inchikey14(inchikey)
    if ik and ik in by_ikey14:
        return by_ikey14[ik]

    nrm = _normalize_name(compound_name)
    if nrm and len(nrm) >= RA_NAME_MATCH_MIN_LEN:
        # 1) 精确归一化匹配
        hit = by_name.get(nrm)
        if hit:
            return hit
        # 2) 包含匹配(库命中名可能是 "Sinomenine (M+H)" 等)
        for key, entry in by_name.items():
            if nrm == key or key in nrm or nrm in key:
                return entry
    return None


def _classification_is_ra_relevant(classification: Optional[str]) -> bool:
    """命中的化合物分类是否属 RA 高相关类 (大小写不敏感子串匹配)."""
    if not classification:
        return False
    cls = classification.lower()
    return any(tok in cls for tok in RA_RELEVANT_CLASSES)


def _public_entry(entry: dict, *, with_meta: bool = True) -> dict:
    """裁剪 RA 条目用于 JSON 输出(去掉 smiles 等冗长字段)."""
    out = {
        "name_cn": entry["name_cn"],
        "name_en": entry["name_en"],
        "formula": entry["formula"],
        "exact_mass": round(entry["exact_mass"], 4) if entry["exact_mass"] else None,
        "source_herb": entry["source_herb"],
        "mechanism": entry["mechanism"],
        "toxicity_flag": entry["toxicity_flag"],
        "category": entry["category"],
    }
    if with_meta:
        out["activity_note"] = entry["activity_note"]
        out["reference"] = entry["reference"]
    return out


def screen(results: list[dict], net: dict, panel: list[dict],
           *, max_propagated: int = 300) -> dict:
    """对已完成分析的 results + 分子网络 net 做 RA 相关性解读.

    返回:
      direct_hits    直接命中(化合物在 RA 参考库)列表, 按 final_score 降序
      propagated     候选类似物(家族含 RA 直接命中)列表
      class_relevant 命中分类属 RA 相关类列表
      summary        {totals, by_mechanism, by_herb, by_toxicity}
      network_overlay {node_idx -> {ra_flag, mechanism}}
    """
    by_ikey14: dict[str, dict] = {}
    by_name: dict[str, dict] = {}
    for entry in panel:
        if entry["inchikey14"]:
            by_ikey14.setdefault(entry["inchikey14"], entry)
        for nm in (entry["name_en"], entry["name_cn"]):
            nrm = _normalize_name(nm)
            if len(nrm) >= RA_NAME_MATCH_MIN_LEN:
                by_name.setdefault(nrm, entry)

    direct_hits: list[dict] = []
    class_relevant: list[dict] = []
    direct_spec_ids: set[str] = set()

    # Pass 1 — 直接命中 + Pass 2 — 类别相关
    for r in results or []:
        m = r.get("match")
        if not isinstance(m, dict) or not m:
            continue
        cname = m.get("compound_name")
        entry = match_to_panel(cname, m.get("inchikey"),
                               by_ikey14=by_ikey14, by_name=by_name)
        spec_id = r.get("spectrum_id")
        base = {
            "spectrum_id": spec_id,
            "precursor_mz": r.get("precursor_mz"),
            "retention_time": r.get("retention_time"),
            "compound_name": cname,
            "final_score": m.get("final_score"),
            "cosine": m.get("cosine"),
            "classification": m.get("classification"),
            "inchikey": m.get("inchikey"),
        }
        if entry:
            row = {**base, "ra_flag": "direct",
                   "ra_entry": _public_entry(entry, with_meta=True)}
            direct_hits.append(row)
            if spec_id:
                direct_spec_ids.add(spec_id)
        elif _classification_is_ra_relevant(m.get("classification")):
            class_relevant.append({**base, "ra_flag": "class",
                                   "ra_entry": None})

    direct_hits.sort(key=lambda x: (x.get("final_score") or 0.0), reverse=True)
    class_relevant.sort(key=lambda x: (x.get("final_score") or 0.0), reverse=True)

    # Pass 3 — 分子家族传播: 家族含 RA 直接命中 → 其余节点为候选类似物
    propagated: list[dict] = []
    propagated_total = 0
    network_overlay: dict[str, dict] = {}

    nodes = net.get("nodes", []) if net else []
    fam_nodes: dict[Any, list[dict]] = defaultdict(list)
    node_by_spec: dict[str, dict] = {}
    for n in nodes:
        fid = n.get("family_id")
        if fid is not None:
            fam_nodes[fid].append(n)
        sid = n.get("spectrum_id")
        if sid:
            node_by_spec[sid] = n

    if direct_spec_ids and nodes:
        for fid, members in fam_nodes.items():
            has_direct = any(m.get("spectrum_id") in direct_spec_ids for m in members)
            if not has_direct:
                continue
            for n in members:
                sid = n.get("spectrum_id")
                is_direct = sid in direct_spec_ids
                mechanism = ""
                if is_direct:
                    mechanism = (match_to_panel(
                        (n.get("match") or {}).get("compound_name"),
                        (n.get("match") or {}).get("inchikey"),
                        by_ikey14=by_ikey14, by_name=by_name) or {}).get("mechanism", "")
                    network_overlay[str(n.get("idx"))] = {
                        "ra_flag": "direct", "mechanism": mechanism}
                    continue
                # 非直接但同家族 → 候选类似物
                propagated_total += 1
                network_overlay[str(n.get("idx"))] = {
                    "ra_flag": "propagated", "mechanism": mechanism}
                if len(propagated) < max_propagated:
                    m = n.get("match") or {}
                    propagated.append({
                        "node_idx": n.get("idx"),
                        "spectrum_id": sid,
                        "precursor_mz": n.get("precursor_mz"),
                        "retention_time": n.get("retention_time"),
                        "peaks_count": n.get("peaks_count"),
                        "annotation_source": n.get("annotation_source"),
                        "propagated_compound": n.get("propagated_compound"),
                        "propagated_class": n.get("propagated_class"),
                        "ra_flag": "propagated",
                    })

    # 给 class_relevant 节点也打 overlay 标记(若有对应网络节点)
    for r in class_relevant:
        sid = r.get("spectrum_id")
        n = node_by_spec.get(sid)
        if n is not None and str(n.get("idx")) not in network_overlay:
            network_overlay[str(n.get("idx"))] = {
                "ra_flag": "class", "mechanism": ""}

    # ── 汇总统计 ──
    by_mechanism: dict[str, int] = defaultdict(int)
    by_herb: dict[str, int] = defaultdict(int)
    by_toxicity: dict[str, int] = defaultdict(int)
    for d in direct_hits:
        e = d.get("ra_entry") or {}
        if e.get("mechanism"):
            by_mechanism[e["mechanism"]] += 1
        if e.get("source_herb"):
            by_herb[e["source_herb"]] += 1
        if e.get("toxicity_flag"):
            by_toxicity[e["toxicity_flag"]] += 1

    return {
        "direct_hits": direct_hits,
        "propagated": propagated,
        "propagated_total": propagated_total,
        "class_relevant": class_relevant,
        "summary": {
            "n_direct": len(direct_hits),
            "n_propagated": propagated_total,
            "n_propagated_shown": len(propagated),
            "n_class_relevant": len(class_relevant),
            "by_mechanism": dict(by_mechanism),
            "by_herb": dict(by_herb),
            "by_toxicity": dict(by_toxicity),
            "panel_size": len(panel),
        },
        "network_overlay": network_overlay,
    }
