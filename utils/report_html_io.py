"""raw_blank_analysis_v2 HTML 分析报告解析

外部前处理工具 (raw_blank_analysis_v2) 把 LC-MS 原始数据 (.d/.mzML) 加工成一份
HTML 峰级报告: 每个 UV 色谱峰一张卡片, 卡片里含该峰的 MS2 谱图 (前体 m/z + scan +
碎片 m/z/% + 离子模式)。原始 .d 数据常留在采集仪器 (Windows) 上拿不到, 所以这类报告
本身就是化合物鉴定的可用输入 —— MS/MS 谱图全部嵌在 HTML 里。

本模块把报告解析成与 utils.mgf_io.parse_mgf_file / utils.mzml_io.parse_mzml_file
相同的谱图记录 schema, 使其能无缝进入 web/server.py + web/app.py 的现有分析管线
(谱图匹配/打分/导出), 无需改动匹配层。

输出每条谱图:
    {
        "spectrum_id": int,          # 取 MS2 的 scan= id, 缺失则用峰序号
        "precursor_mz": float,        # 前体 m/z (前体 m/z=...)
        "peaks": [[mz, int], ...],    # MS2 碎片, int 用报告中相对强度 % (下游按 max 归一化)
        "title": str,                 # "UV峰#N RT=... prec=... scan=... POL"
        "retention_time": float|None, # 分钟
        "scan_number": int|None,
        "polarity": "POS"|"NEG"|None, # 来自卡片 离子模式 <strong>
        "polarity_source": "report_meta"|None,
    }
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

# 卡片 h3: "UV峰#1: RT=7.515min (7.4901-7.6714)"
_RT_RE = re.compile(r"RT\s*=\s*([\d.]+)")
# ms2-header 文本: "前体 m/z=729.2604 (scan=452672)"
_PREC_RE = re.compile(r"前体\s*m/z\s*=\s*([\d.]+)")
_SCAN_RE = re.compile(r"scan\s*=\s*(\d+)")
# 极性: ms1-meta 里 "离子模式: <strong>NEG</strong>"
_POL_RE = re.compile(r"离子模式[:：]\s*</?strong>?\s*(NEG|POS)", re.IGNORECASE)
# 相对强度: "100.0%" / "84.4%"
_PCT_RE = re.compile(r"(-?[\d.]+)\s*%")


def is_raw_report_html(file_path) -> bool:
    """探测文件是否为可解析的 raw_blank_analysis 风格 HTML 报告。

    判据: 是 .html/.htm 且文本里同时出现 "前体 m/z" 与 "MS2 碎片" —— 这是峰值卡片的
    硬签名, 既能命中 raw_blank_analysis_v2 这类报告, 又能挡住 mssystem 自己导出的
    结果 HTML (无 MS2 碎片表) 和任意网页。
    """
    try:
        p = Path(file_path)
        if p.suffix.lower() not in (".html", ".htm"):
            return False
    except Exception:
        return False
    try:
        # 只读文本前 ~512KB 足够看签名, 避免对超大文件全量读
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            head = f.read(512 * 1024)
        return ("前体" in head and "m/z" in head and "MS2" in head
                and ("碎片" in head or "前体 m/z" in head))
    except Exception:
        return False


def _parse_intensity_pct(cell_text: str) -> Optional[float]:
    """从单元格文本 (如 '100.0%') 抽相对强度数值, 失败返回 None。"""
    if not cell_text:
        return None
    m = _PCT_RE.search(cell_text)
    if not m:
        return None
    try:
        v = float(m.group(1))
    except ValueError:
        return None
    return v if v > 0 else None


def _row_mz(tr) -> Optional[float]:
    """从碎片表行抽 m/z: 优先 button[data-mz] (最干净), 兜底 mz-link span 文本。"""
    btn = tr.find("button", attrs={"data-mz": True})
    if btn:
        try:
            return float(btn["data-mz"])
        except (TypeError, ValueError):
            pass
    span = tr.find(class_="mz-link")
    if span:
        m = re.search(r"[\d.]+", span.get_text())
        if m:
            try:
                return float(m.group(0))
            except ValueError:
                pass
    return None


def _row_intensity_pct(tr) -> Optional[float]:
    """从碎片表行抽相对强度 %: 找文本里带 % 的单元格。"""
    for td in tr.find_all("td"):
        v = _parse_intensity_pct(td.get_text())
        if v is not None:
            return v
    return None


def parse_raw_report_html(file_path) -> list[dict]:
    """解析 raw_blank_analysis HTML 报告, 返回 MS2 谱图记录列表。

    每个 UV 峰卡片 → 一条 MS2 谱图 (卡片自带该峰的代表 MS2)。MS1 候选 / 加合离子检测
    表不提取 (匹配用不到, 且与 MS2 碎片表结构相同易混)。碎片强度用报告中相对强度 %
    作为绝对强度输入 —— 下游 clean_peaks 会按 max 归一化 + 按前体过滤, 行为与 MGF 一致。
    """
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        soup = BeautifulSoup(f.read(), "lxml")

    records: list[dict] = []
    cards = soup.find_all("section", class_="peak-card")
    if not cards:
        # 兜底: 某些报告可能用 <div class="peak-card"> 而非 section
        cards = soup.find_all(class_="peak-card")

    for idx, card in enumerate(cards, start=1):
        title_text = card.get_text(" ", strip=True)

        # ── MS2 谱图: 以 "MS2 碎片" 标题锚定, 避开 MS1 top5 / 加合检测表 ──
        ms2_h4 = None
        for h4 in card.find_all(["h4", "h3"]):
            if "MS2" in h4.get_text() and "碎片" in h4.get_text():
                ms2_h4 = h4
                break
        if ms2_h4 is None:
            logger.debug("卡片 #%d 无 MS2 碎片段, 跳过", idx)
            continue

        # 前体 m/z + scan: 卡片里 class=ms2-header 的 div (在 MS2 标题与碎片表之间)
        ms2_header = card.find(class_="ms2-header")
        header_txt = ms2_header.get_text(" ", strip=True) if ms2_header else ""
        if not header_txt:
            # 兜底: 取 MS2 标题之后第一个 div 的文本
            nxt = ms2_h4.find_next("div")
            header_txt = nxt.get_text(" ", strip=True) if nxt else ""

        pmz = None
        m = _PREC_RE.search(header_txt)
        if m:
            try:
                pmz = float(m.group(1))
            except ValueError:
                pmz = None

        # 离子模式: 卡片级 (ms1-meta 里 <strong>NEG/POS</strong>); 缺失再扫全文兜底
        polarity = None
        meta = card.find(class_="ms1-meta")
        strong = meta.find("strong") if meta else None
        strong_txt = strong.get_text(strip=True).upper() if strong else ""
        if strong_txt in {"NEG", "POS"}:
            polarity = strong_txt
        else:
            pm = _POL_RE.search(title_text)
            if pm:
                polarity = pm.group(1).upper()

        # scan id (spectrum_id 用), 缺失回退峰序号
        scan_id = None
        sm = _SCAN_RE.search(header_txt) or _SCAN_RE.search(title_text)
        if sm:
            try:
                scan_id = int(sm.group(1))
            except ValueError:
                scan_id = None
        if scan_id is None:
            # ms1-meta 的 scan_id: 451741~452075 是 MS1 范围, 非 MS2 scan; 仅在无 MS2 scan 时兜底取首个
            pm2 = re.search(r"scan_id[:：]\s*(\d+)", title_text)
            if pm2:
                try:
                    scan_id = int(pm2.group(1))
                except ValueError:
                    scan_id = None
        spec_id = scan_id if scan_id is not None else idx

        # RT (分钟)
        rt = None
        h3 = card.find(["h3"])
        h3_txt = h3.get_text(" ", strip=True) if h3 else title_text
        rm = _RT_RE.search(h3_txt)
        if rm:
            try:
                rt = float(rm.group(1))
            except ValueError:
                rt = None

        # ── MS2 碎片表: MS2 标题后的第一张 table ──
        frag_table = ms2_h4.find_next("table")
        peaks: list[list[float]] = []
        if frag_table is not None:
            for tr in frag_table.find_all("tr"):
                tds = tr.find_all("td")
                if len(tds) < 2:
                    continue  # 表头 / 空行
                mz = _row_mz(tr)
                pct = _row_intensity_pct(tr)
                if mz is None or pct is None:
                    continue
                peaks.append([float(mz), float(pct)])

        # 无前体或无碎片 → 无可匹配谱图, 跳过
        if pmz is None or not peaks:
            logger.info("卡片 #%d 无有效前体/碎片 (pmz=%s, n_frag=%d), 跳过",
                        idx, pmz, len(peaks))
            continue

        title = (f"UV峰#{idx} RT={rt}min prec={pmz:.4f}"
                 + (f" scan={scan_id}" if scan_id is not None else "")
                 + (f" {polarity}" if polarity else ""))

        records.append({
            "spectrum_id": spec_id,
            "precursor_mz": float(pmz),
            "peaks": peaks,
            "title": title,
            "retention_time": rt,
            "scan_number": scan_id,
            "polarity": polarity,
            "polarity_source": "report_meta" if polarity else None,
            # raw_blank_analysis 的报告只导出代表 MS2 的少量峰；当前实测
            # sample_report.html 每谱固定 5 峰。把限制传给结果层，避免把
            # HTML 摘要当成完整 mzML/MGF 谱图解释。
            "input_format": "raw_blank_analysis_html",
            "limited_fragment_list": len(peaks) <= 5,
            "input_evidence_note": (
                f"HTML 报告仅提供 {len(peaks)} 个 MS2 碎片；"
                "化合物名级鉴定可能因谱图被截断而无法成立"
            ) if len(peaks) <= 5 else "",
        })

    logger.info("HTML 报告解析完成: %s → %d 条 MS2 谱图 (共 %d 张峰值卡片)",
                Path(file_path).name, len(records), len(cards))
    return records


def records_to_mgf(records: list[dict], out_path) -> int:
    """把解析出的谱图记录写成 MGF 文件 (供离线复用/备份/其它工具)。

    MGF 用 CHARGE 锁极性 (-1 / 1), 与 utils.mgf_io.parse_mgf_file 的极性推导对齐,
    避免再导入时极性丢失导致负离子样本匹配出 [M+H]+ 假阳性。返回写出的谱图数。
    """
    n = 0
    with open(out_path, "w", encoding="utf-8") as f:
        for r in records:
            pmz = r["precursor_mz"]
            pol = r.get("polarity")
            charge = "1-" if pol == "NEG" else ("1+" if pol == "POS" else "1+"
                                                if pmz > 0 else "0")
            title = f"{r.get('title', '')} scanId={r.get('scan_number') or r.get('spectrum_id')}"
            rt = r.get("retention_time")
            f.write("BEGIN IONS\n")
            f.write(f"TITLE={title}\n")
            f.write(f"PEPMASS={pmz}\n")
            f.write(f"CHARGE={charge}\n")
            if pol:
                f.write(f"IONMODE={pol}\n")
            if rt is not None:
                # RTINSECONDS 单位为秒
                f.write(f"RTINSECONDS={rt * 60.0}\n")
            for mz, inten in r["peaks"]:
                f.write(f"{mz} {inten}\n")
            f.write("END IONS\n\n")
            n += 1
    return n
