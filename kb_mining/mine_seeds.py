"""
mssystem v3 — iridoid/saponin 种子挖掘 (补 NP 关键类经验规则)
==============================================================
DB 的 classification 不含 iridoid/saponin 标签 (无法用 mine_rules 挖), 故改用
**已知标志化合物名** 作种子, 从库中找它们的真实谱图, 挖诊断离子/中性丢失。
合并进 kb_v3/empirical_rules.json 的 iridoid / saponin 键。

cd /Users/zizi/Desktop/Mssystem/mssystem_v3
python3 kb_mining/mine_seeds.py
"""
from __future__ import annotations

import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

BIN = 0.025
MZ_MIN, MZ_MAX = 30.0, 1000.0
LOSS_MIN, LOSS_MAX = 1.0, 400.0
FREQ_MIN_FRAC = 0.20
TOPK = 20
MIN_N = 8  # 至少 8 条种子谱图才挖
OUT = ROOT / "kb_v3" / "empirical_rules.json"

# 种子标志化合物 (英文/中文关键词, 库 compound_name 子串匹配)
SEEDS = {
    "iridoid": ["catalpol", "geniposide", "genipin", "aucubin", "asperuloside",
                "loganin", "secologanin", "oleuropein", "gentiopicroside", "gentiopicros",
                "swer tiamarin", "sweroside", "morroniside", "paederoside", "picroside",
                "harpagide", "gardoside", "verbenalin", "monotropein", "geniposidic"],
    "saponin": ["ginsenoside", "glycyrrhiz", "astragaloside", "saikosaponin", "dioscin",
                "diosgenin", "platycodin", "gypenoside", "notoginsenoside", "protopanax",
                "oleanolic", "ursolic", "maslinic", "hederagenin", "escin", "diosgen",
                "polygalacic", "echinocystic", "quillaic", "gypsogenin"],
}


def polarity_of(adduct, charge):
    a = (adduct or "").strip()
    if "+" in a and "-" not in a: return "pos"
    if "-" in a and "+" not in a: return "neg"
    try:
        c = int(charge) if charge not in (None, "") else 0
        if c > 0: return "pos"
        if c < 0: return "neg"
    except Exception:
        pass
    return "unknown"


def mine_class(con, markers):
    """从 con 找 markers 对应化合物的谱图, 挖诊断离子/中性丢失。"""
    cur = con.cursor()
    diag = defaultdict(int); loss = defaultdict(int); n_spec = defaultdict(int)
    # 用 LIKE OR 找化合物
    like_sql = " OR ".join(["LOWER(c.compound_name) LIKE ?"] * len(markers))
    params = [f"%{m.lower()}%" for m in markers]
    rows = cur.execute(f"""
        SELECT s.precursor_mz, s.peaks_json, s.adduct, s.charge
        FROM Spectra s JOIN Mapping m ON m.spectrum_id=s.id
        JOIN Compounds c ON m.compound_id=c.id
        WHERE ({like_sql})
    """, params).fetchall()
    for r in rows:
        try:
            pmz = float(r["precursor_mz"]); peaks = json.loads(r["peaks_json"] or "[]")
        except Exception:
            continue
        if not peaks or pmz <= 0:
            continue
        pol = polarity_of(r["adduct"], r["charge"])
        if pol == "unknown":
            continue
        n_spec[pol] += 1
        sd = set(); sl = set()
        for p in peaks:
            try:
                mz, it = float(p[0]), float(p[1])
            except Exception:
                continue
            if it <= 0: continue
            if MZ_MIN <= mz <= MZ_MAX: sd.add(int(round(mz / BIN)))
            d = pmz - mz
            if LOSS_MIN <= d <= LOSS_MAX: sl.add(int(round(d / BIN)))
        for b in sd: diag[(pol, b)] += 1
        for b in sl: loss[(pol, b)] += 1
    out = {}
    for pol, ntot in n_spec.items():
        if ntot < MIN_N:
            continue
        thr = max(FREQ_MIN_FRAC * ntot, 2)
        di = sorted(((b * BIN + BIN / 2, c / ntot, c) for (pp, b), c in diag.items()
                     if pp == pol and c >= thr), key=lambda x: -x[1])[:TOPK]
        nl = sorted(((b * BIN + BIN / 2, c / ntot, c) for (pp, b), c in loss.items()
                     if pp == pol and c >= thr), key=lambda x: -x[1])[:TOPK]
        out[pol] = {"n_spectra": ntot,
                    "diagnostic_ions": [{"mz": round(m, 4), "freq": round(f, 3), "n_support": int(c)} for m, f, c in di],
                    "neutral_losses": [{"loss_da": round(m, 4), "freq": round(f, 3), "n_support": int(c)} for m, f, c in nl]}
    return out


def main():
    print("=== 种子挖掘 iridoid/saponin ===")
    emp = json.load(open(OUT, encoding="utf-8")) if OUT.exists() else {}
    for db_name, db_path in [("small", "data/local_library_small.db"), ("full", "data/local_library.db")]:
        con = sqlite3.connect(str(ROOT / db_path)); con.row_factory = sqlite3.Row
        for cls, markers in SEEDS.items():
            res = mine_class(con, markers)
            if not res:
                continue
            # 合并: 取 n_spectra 更大的
            cur = emp.get(cls, {})
            for pol, d in res.items():
                if pol not in cur or d["n_spectra"] > cur[pol].get("n_spectra", 0):
                    cur[pol] = d
            emp[cls] = cur
            for pol, d in res.items():
                top_nl = ", ".join(f"{x['loss_da']}" for x in d["neutral_losses"][:3])
                print(f"  [{db_name}] {cls} {pol}: {d['n_spectra']}谱, {len(d['diagnostic_ions'])}诊断离子, {len(d['neutral_losses'])}丢失 [{top_nl}]")
        con.close()
    # 中性丢失注释
    NL = {18.0106:"H2O",43.9898:"CO2",46.0055:"HCOOH",132.0423:"pentose",146.0579:"deoxyhexose",
          162.0528:"hexose",176.0321:"glucuronic_acid",180.0634:"hexose+H2O",246.0618:"acetylhexose",308.0905:"rutinose"}
    for cls in ("iridoid","saponin"):
        for pol,d in emp.get(cls,{}).items():
            for nl in d["neutral_losses"]:
                for k,v in NL.items():
                    if abs(nl["loss_da"]-k)<0.03: nl["annotation"]=v
    json.dump(emp, open(OUT,"w",encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\n✅ 合并到 {OUT}: iridoid={list(emp.get('iridoid',{}))} saponin={list(emp.get('saponin',{}))}")


if __name__ == "__main__":
    main()
