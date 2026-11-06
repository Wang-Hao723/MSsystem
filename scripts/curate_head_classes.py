"""
mssystem v3 — P2·2a 头部化合物分类策展 (curated)
====================================================
对头部高频化合物 (RA 参考物 + 17 标准品 + 常见 TCM/植化明星物) 用**策展名字→类**
字典 (文献/药植学常识) 替换 ML 猜测的 classification_v3, 标 source='curated', conf=1.0。

**只策展类无歧义的明星化合物** (策展质量 > 数量; 名字唯一指向一个 v3 类)。
含糊的 (苷+萜 双重归属等) 不强标, 留给人工/NP Classifier。
仅写 P2 副本; S5 红线回归通过后才换入 live。

用法:
  cd /Users/zizi/Desktop/Mssystem/mssystem_v3
  python3 scripts/curate_head_classes.py --db data/local_library_p2.db
  python3 scripts/curate_head_classes.py --db data/local_library_p2_full.db
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# ── 策展名字 → v3 类 (高置信, 明星化合物; 含同义词) ──
# 仅收录类无歧义者。多归属 (苷+萜/苷+黄酮) 取最具体的 NP 骨架类。
CURATED = {
    # 黄酮 (aglycone + 苷)
    "quercetin": "flavonoid", "kaempferol": "flavonoid", "luteolin": "flavonoid",
    "apigenin": "flavonoid", "naringenin": "flavonoid", "chrysin": "flavonoid",
    "catechin": "flavonoid", "epicatechin": "flavonoid", "myricetin": "flavonoid",
    "isorhamnetin": "flavonoid", "galangin": "flavonoid", "diosmetin": "flavonoid",
    "genistein": "flavonoid", "daidzein": "flavonoid", "biochanin": "flavonoid",
    "formononetin": "flavonoid", "hesperetin": "flavonoid", "eriodictyol": "flavonoid",
    "taxifolin": "flavonoid", "rutin": "flavonoid", "hesperidin": "flavonoid",
    "naringin": "flavonoid", "hyperoside": "flavonoid", "isoquercitrin": "flavonoid",
    "astragalin": "flavonoid", "quercitrin": "flavonoid", "vitexin": "flavonoid",
    "orientin": "flavonoid", "baicalin": "flavonoid", "baicalein": "flavonoid",
    "wogonin": "flavonoid", "puerarin": "flavonoid", "liquiritin": "flavonoid",
    "liquiritigenin": "flavonoid", "isoliquiritin": "flavonoid",
    # 生物碱
    "berberine": "alkaloid", "palmatine": "alkaloid", "jatrorrhizine": "alkaloid",
    "coptisine": "alkaloid", "magnoflorine": "alkaloid", "sinomenine": "alkaloid",
    "ephedrine": "alkaloid", "pseudoephedrine": "alkaloid", "atropine": "alkaloid",
    "hyoscyamine": "alkaloid", "scopolamine": "alkaloid", "morphine": "alkaloid",
    "codeine": "alkaloid", "papaverine": "alkaloid", "colchicine": "alkaloid",
    "nicotine": "alkaloid", "pilocarpine": "alkaloid", "reserpine": "alkaloid",
    "yohimbine": "alkaloid", "ajmaline": "alkaloid", "vincristine": "alkaloid",
    "vinblastine": "alkaloid", "matrine": "alkaloid", "oxymatrine": "alkaloid",
    "sophocarpine": "alkaloid", "cytisine": "alkaloid", "tetrandrine": "alkaloid",
    "fangchinoline": "alkaloid", "aconitine": "alkaloid", "mesaconitine": "alkaloid",
    "hypaconitine": "alkaloid", "lycorine": "alkaloid", "galantamine": "alkaloid",
    # 黄嘌呤生物碱
    "caffeine": "xanthine", "theobromine": "xanthine", "theophylline": "xanthine",
    # 酚酸
    "chlorogenic acid": "phenolic_acid", "caffeic acid": "phenolic_acid",
    "ferulic acid": "phenolic_acid", "p-coumaric acid": "phenolic_acid",
    "coumaric acid": "phenolic_acid", "sinapic acid": "phenolic_acid",
    "gallic acid": "phenolic_acid", "protocatechuic acid": "phenolic_acid",
    "vanillic acid": "phenolic_acid", "syringic acid": "phenolic_acid",
    "salicylic acid": "phenolic_acid", "rosmarinic acid": "phenolic_acid",
    "neochlorogenic acid": "phenolic_acid", "isochlorogenic acid": "phenolic_acid",
    "cryptochlorogenic acid": "phenolic_acid", "3-caffeoylquinic": "phenolic_acid",
    "5-caffeoylquinic": "phenolic_acid", "4-caffeoylquinic": "phenolic_acid",
    # 香豆素
    "coumarin": "coumarin", "umbelliferone": "coumarin", "aesculetin": "coumarin",
    "esculetin": "coumarin", "scopoletin": "coumarin", "psoralen": "coumarin",
    "bergapten": "coumarin", "xanthyletin": "coumarin", "osthol": "coumarin",
    "imperatorin": "coumarin", "phellopterin": "coumarin", "isopimpinellin": "coumarin",
    # 木脂素
    "pinoresinol": "lignan", "pinoresinol diglucoside": "lignan",
    "sesamin": "lignan", "magnolol": "lignan", "honokiol": "lignan",
    "schisandrin": "lignan", "schizandrin": "lignan", "podophyllotoxin": "lignan",
    "arctigenin": "lignan", "arctiin": "lignan", "enterolactone": "lignan",
    "secoisolariciresinol": "lignan", "matairesinol": "lignan",
    # 醌
    "emodin": "quinone", "rhein": "quinone", "chrysophanol": "quinone",
    "physcion": "quinone", "aloe-emodin": "quinone", "alizarin": "quinone",
    "purpurin": "quinone", "shikonin": "quinone", "alkannin": "quinone",
    "juglone": "quinone", "plumbagin": "quinone", "lapachol": "quinone",
    "tanshinone": "quinone", "cryptotanshinone": "quinone",
    # 二苯乙烯
    "resveratrol": "stilbene", "pterostilbene": "stilbene",
    "oxyresveratrol": "stilbene", "polydatin": "stilbene", "piceid": "stilbene",
    # 山酮
    "mangiferin": "xanthone", "gentisin": "xanthone", "bellidifolin": "xanthone",
    "swertianin": "xanthone", "decussatin": "xanthone",
    # 甾体 (含强心苷)
    "cholesterol": "steroid", "stigmasterol": "steroid", "beta-sitosterol": "steroid",
    "β-sitosterol": "steroid", "sitosterol": "steroid", "campesterol": "steroid",
    "ergosterol": "steroid", "diosgenin": "steroid", "sarsasapogenin": "steroid",
    "tigogenin": "steroid", "ecdysone": "steroid", "ecdysterone": "steroid",
    "20-hydroxyecdysone": "steroid", "cyasterone": "steroid", "estradiol": "steroid",
    "estrone": "steroid", "progesterone": "steroid", "testosterone": "steroid",
    "cortisol": "steroid", "cortisone": "steroid", "digoxin": "steroid",
    "digitoxin": "steroid", "ouabain": "steroid",
    # 萜
    "artemisinin": "terpenoid", "artemether": "terpenoid", "dihydroartemisinin": "terpenoid",
    "atractylon": "terpenoid", "atractylenolide": "terpenoid", "bisabolol": "terpenoid",
    "caryophyllene": "terpenoid", "eugenol": "terpenoid", "thymol": "terpenoid",
    "carvacrol": "terpenoid", "menthol": "terpenoid", "camphor": "terpenoid",
    "geraniol": "terpenoid", "loganing": "terpenoid", "loganin": "terpenoid",
    "triptolide": "terpenoid", "tripdiolide": "terpenoid", "triptonide": "terpenoid",
    "celastrol": "terpenoid", "andrographolide": "terpenoid", "ginkgolide": "terpenoid",
    "bilobalide": "terpenoid", "parthenolide": "terpenoid", "costunolide": "terpenoid",
    # 皂苷
    "ginsenoside": "saponin", "notoginsenoside": "saponin", "astragaloside": "saponin",
    "gypenoside": "saponin", "saikosaponin": "saponin", "dioscin": "saponin",
    "protodioscin": "saponin", "glycyrrhizic acid": "saponin", "glycyrrhizin": "saponin",
    "glycyrrhetinic acid": "saponin", "glycyrrhetinate": "saponin",
    "asiaticoside": "saponin", "madecassoside": "saponin", "hederagenin": "saponin",
    "oleanolic acid": "saponin", "ursolic acid": "saponin", "boswellic acid": "saponin",
    # 核苷
    "adenosine": "nucleoside", "guanosine": "nucleoside", "uridine": "nucleoside",
    "cytidine": "nucleoside", "thymidine": "nucleoside", "inosine": "nucleoside",
}


def _norm(name: str) -> str:
    if not name:
        return ""
    s = str(name).strip().lower()
    s = re.sub(r"[*_\-]+", " ", s)        # 分隔符归一
    s = re.sub(r"\s+", " ", s)
    return s


# 名字含这些标记 = 异构体对/干扰物标注残渣 (如 "Ornithine_Estrone (known isomers...)")
# 不是单化合物, 策展会误标 (estrone 出现在氨基酸名里 ≠ 甾体化合物) → 一律跳过。
_ARTIFACT = re.compile(r"known isomer|isobaric|isomers:|co-elut|peak in run", re.I)


def is_artifact_name(name: str) -> bool:
    return bool(_ARTIFACT.search(name or ""))


def predict_curated(name: str) -> str | None:
    """返回策展类或 None (无线索)。最长键优先 (e.g. 'glycyrrhizic acid' > 'acid')。"""
    n = _norm(name)
    if not n:
        return None
    for key in sorted(CURATED, key=len, reverse=True):   # 长键先, 防短键误吞
        k = _norm(key)
        if not k:
            continue
        # 整词匹配: n 等于 k, 或以 k 起头(后接空格), 或含 ' (k' / 'k-' 边界
        if n == k or n.startswith(k + " ") or (" " + k) in n or re.search(
                r"\b" + re.escape(k) + r"\b", n):
            return CURATED[key]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(ROOT / "data" / "local_library_p2.db"))
    ap.add_argument("--dry", action="store_true", help="只统计不写")
    args = ap.parse_args()
    db = str(Path(args.db).resolve())
    print(f"=== P2·2a 头部策展 ({'DRY' if args.dry else 'WRITE'}) ===\n  DB: {db}")

    con = sqlite3.connect(db)
    cur = con.cursor()
    rows = cur.execute(
        "SELECT id, compound_name, classification_v3, classification_source "
        "FROM Compounds").fetchall()
    print(f"  化合物: {len(rows):,} | 策展字典: {len(CURATED)} 键")

    updates, changed, upgraded, same = [], 0, 0, 0
    n_artifact = 0
    change_log = Counter()   # (old_cls -> new_cls)
    for cid, name, cls_v3, src in rows:
        if is_artifact_name(name):
            n_artifact += 1
            continue
        cu = predict_curated(name)
        if cu is None:
            continue
        if cu != (cls_v3 or ""):
            changed += 1
            change_log[f"{cls_v3} -> {cu}"] += 1
        elif src != "curated":
            upgraded += 1
        else:
            same += 1
            continue
        updates.append((cu, "curated", 1.0, cid))

    print(f"\n  匹配策展: {len(updates)+same} 化合物 (跳过 {n_artifact} 个异构体对/残渣名)")
    print(f"    类被修正 (ML/旧标 → 策展): {changed}")
    print(f"    source 升级为 curated (类不变): {upgraded}")
    print(f"    已是 curated 跳过: {same}")
    print(f"\n  Top 类修正 (old -> new):")
    for k, v in change_log.most_common(20):
        print(f"    {v:5d}  {k}")

    if updates and not args.dry:
        cur.executemany(
            "UPDATE Compounds SET classification_v3=?, classification_source=?, "
            "classification_confidence=? WHERE id=?", updates)
        con.commit()
        print(f"\n  ✅ 写入 {len(updates)} 行 (source=curated, conf=1.0)")
    elif args.dry:
        print(f"\n  (dry-run, 未写入)")
    con.close()


if __name__ == "__main__":
    main()
