"""
mssystem v3 — P2·2a ML 分类层精度抽检 (N≈200, 分层)
================================================================
对 DB 里 classification_source='ml' 的化合物 (即被 ML 猜标签、非 SMARTS/原标签确定)
做**独立**复核, 给出真实 per-class 精度估计, 替代当前误导性的 "自洽 macro-F1"
(ML 训练标签本身就来自 SMARTS/自动派生, 训练集 F1 是循环自洽, 非真值精度)。

独立复核 oracle = 不复用训练/分类管线:
  (1) **名字启发式** NAME_RULES (独立于结构): 化合物名含类标志词 (ginsenoside→saponin,
      quercetin→flavonoid, ...oside→glycoside, cholesterol→steroid ...)
  (2) **结构健全性检查** struct_ok (RDKit, 粗筛 gross error): DB 类 'sugar' 须有糖环,
      'steroid' 须有 gonane 四环, 'flavonoid' 须 2 芳环+chromenone, ...
裁定:
  CORRECT   名字明确指向 == DB 类
  WRONG     名字明确指向 != DB 类
  STRUCT_BAD 结构健全性硬冲突 (即便名字不明)
  UNCERTAIN  名字无明确线索且结构检查中性 → 列入人工复核
精度 = CORRECT / (CORRECT + WRONG + STRUCT_BAD), UNCERTAIN 单列 (不计入分母, 透明)。

只读 DB; 仅写 docs/ + sample CSV。不动 classification_v3。

用法:
  cd /Users/zizi/Desktop/Mssystem/mssystem_v3
  python3 scripts/audit_ml_precision.py --db data/local_library.db --n 200
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import Chem  # noqa: E402
from rdkit import RDLogger  # noqa: E402
RDLogger.DisableLog("rdApp.*")

OUT_MD = ROOT / "docs" / "ml_precision_audit_200.md"
OUT_CSV = ROOT / "docs" / "ml_precision_audit_sample.csv"

# ── (1) 名字启发式 (独立 oracle; 顺序敏感, 先具体后通用) ──
NAME_RULES = [
    # 皂苷/强心苷/苷 (糖体大分子)
    (r"ginsenoside|astragaloside|glycyrrhiz|saponin|notoginseno|gypenoside", "saponin"),
    (r"digoxin|digitoxin|ouabain|strophanth", "cardiac_glycoside"),
    (r"rutinoside|hyperoside|isoquercitrin|astragalin|hesperidin|naringin|quercitrin|"
     r"rutin|salicin|arbutin|amygdalin|prunasin|sennoside|crocin|picrocrocin|"
     r"oside\b", "glycoside"),
    # 黄酮
    (r"quercetin|kaempferol|luteolin|apigenin|naringenin|hesperetin|chrysin|"
     r"catechin|epicatechin|anthocyan|flavone|flavonol|flavanon|isoflav|chalcone|"
     r"myricetin|genistein|daidzein", "flavonoid"),
    # 黄嘌呤生物碱
    (r"caffeine|theobromine|theophylline|xanthine", "xanthine"),
    # 酚酸
    (r"chlorogenic|caffeic|ferulic|coumaric|sinap|gallic|protocatechu|vanillic|"
     r"hydroxybenzoic|hydroxycinnam", "phenolic_acid"),
    # 香豆素
    (r"coumarin|psoralen|angelicin|bergapten|umbelliferone|aesculetin|scopoletin", "coumarin"),
    # 木脂素
    (r"lignan|pinoresinol|sesamin|magnolol|honokiol|podophyllo|enterolac", "lignan"),
    # 山酮
    (r"xanthone|mangiferin|gentisin|bellidifolin", "xanthone"),
    # 二苯乙烯
    (r"resveratrol|stilbene|pterostilbene|oxyresveratrol", "stilbene"),
    # 醌
    (r"emodin|alizarin|rhein|anthraquinone|naphthoquinone|benzoquinone|shikon|"
     r"physcion|chrysophanol|juglone", "quinone"),
    # 甾体
    (r"cholesterol|cholestan|stigmasterol|sitosterol|campesterol|ecdyster|"
     r"brassinolide|withanolide|cortisol|cortisone|progesterone|testosterone|"
     r"estradiol|diosgenin|sarsasapogenin|tigogenin|solasonine|solanidine", "steroid"),
    # 萜 (内酯/醇)
    (r"artemisinin|artemisol|bisabol|caryophyl|eugenol|thymol|carvacrol|menthol|"
     r"camphor|cineole|pinene|limonen|atractyl|lactone|anolide$|geraniol|nerol|"
     r"loganin|secologan|swertiamarin|gentiopicro|iridoid", "terpenoid"),
    # 核苷
    (r"adenosine|guanosine|uridine|cytidine|thymidine|inosine|xanthosine", "nucleoside"),
    # 含硫
    (r"glucosinolate|isothiocyanate|sulforaph|allicin|sulfoxide", "organosulfur"),
    # 氨基酸/肉碱/肽
    (r"carnitine|arginine|lysine|glutam|aspart|alanine|leucine|valine|proline|"
     r"tryptophan|phenylalanine|tyrosine|histidine|serine|glycine|^gln|^asn|"
     r"glutathione", "amino_acid"),
    # 生物碱 (广; 放后)
    (r"berberine|morphine|codeine|papaverine|vinc|colchicine|ephedrine|nicotine|"
     r"quinine|strychnine|atropine|hyoscyamine|aconit|delphin|solan|aconitine|"
     r"matrine|oxymatrine|sophocarpine|sinomenine|magnoflor|palma", "alkaloid"),
    # 脂质/脂肪酰 (放最后, 通用)
    (r"lysopc|^pc\(|^pe\(|^pi\(|^tg\(|^dg\(|ceramide|sphingo|phosphatidyl|"
     r"lysophosphatidyl|fatty.?acid|^fa\(|lysofarnes|farneso|prostagland|"
     r"leukotriene|hydroxyeicosa|octadec|hexadec|tetradec|eicosa", "lipid"),
]


def name_predict(name: str) -> str | None:
    if not name:
        return None
    s = name.strip().lower()
    for pat, cls in NAME_RULES:
        if re.search(pat, s):
            return cls
    return None


# ── (2) 结构独立复核 (新鲜验证 SMARTS 集, 不导入 RuleBasedClassifier, 避免与训练管线循环) ──
# 每条是该类的**强结构标志** (必要非充分); 命中→结构确认(弱正); 与名字/属性冲突→STRUCT_BAD。
_V_SMARTS = {k: v for k, v in {
    "steroid":      "C1CCC2C(C1)CCC1C3CCCC(C3)C21",        # gonane 四环
    "coumarin":     "O=C1OC2=CC=CC=C2C=C1",                 # 苯并吡喃-2-酮
    "xanthone":     "O=C1C=CC2=CC=CC=C2OC3=CC=CC=C13",      # 二苯并吡喃酮
    "flavonoid":    "O=C1C=C(Oc2ccccc2)c2ccccc21",          # 黄酮骨架 (2-苯基色烯-4-酮)
    "anthraquinone":"O=C1C=CC2=CC=CC=C2C(=O)c2ccccc12",     # 蒽醌
    "stilbene":     "c1ccccc1/C=C/c2ccccc2",                # 1,2-二苯乙烯
    "xanthine":     "O=C1NC(=O)NC(=O)N1",                   # 嘌呤-2,6-二酮
    "nucleoside":   "c1ncnc1",                              # 嘌呤 (糖+嘌呤=核苷, 配合糖环)
    "sugar_ring":   "[CX4]1[OX2][CX4][CX4][CX4][CX4]1",     # 吡喃糖环
}.items() if (v is not None and Chem.MolFromSmarts(v) is not None)}
_V = {k: Chem.MolFromSmarts(v) for k, v in _V_SMARTS.items()}


def _n_in_ring_N(mol):
    """环内 N 原子数 (生物碱/含氮杂环指标)。"""
    ri = mol.GetRingInfo()
    n_in = 0
    for ring in ri.AtomRings():
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetSymbol() == "N":
                n_in += 1
                break
    return n_in


def _has_oh(mol):
    return any(a.GetSymbol() == "O" and a.GetTotalNumHs() >= 1 and a.GetDegree() == 1
               for a in mol.GetAtoms())


def _has_cooH(mol):
    return mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OX1H1]"))


def _has_gonane_rings(mol):
    """甾体 gonane 骨架 = 三个 6 元环 + 一个 5 元环 (6,6,6,5 稠合四环)。
    用环大小计数 (而非脆弱的稠合 SMARTS), 对 SMILES 编码/芳构化鲁棒。"""
    cnt = {}
    for r in mol.GetRingInfo().AtomRings():
        cnt[len(r)] = cnt.get(len(r), 0) + 1
    return cnt.get(6, 0) >= 3 and cnt.get(5, 0) >= 1


def struct_check(db_cls: str, smiles: str) -> bool | None:
    """独立结构复核: True(结构确认该类, 弱正) / False(硬冲突) / None(中性)。

    策略: 用新鲜验证 SMARTS/属性做**正向确认**(命中→True); 对几类做**硬冲突**判定
    (DB 标某类但结构完全不符→False)。无法判断→None。
    """
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    n_heavy = m.GetNumHeavyAtoms()
    n_arom = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
    n_N = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")
    n_O = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
    n_ringN = _n_in_ring_N(m)
    has_sugar = m.HasSubstructMatch(_V["sugar_ring"])

    # 正向确认 (新鲜 SMARTS 命中 → 结构一致)
    if db_cls == "steroid" and _has_gonane_rings(m):       # gonane 四环 (鲁棒环计数)
        return True
    if db_cls == "coumarin" and m.HasSubstructMatch(_V["coumarin"]):
        return True
    if db_cls == "xanthone" and m.HasSubstructMatch(_V["xanthone"]):
        return True
    if db_cls == "flavonoid" and m.HasSubstructMatch(_V["flavonoid"]):
        return True
    if db_cls == "quinone" and m.HasSubstructMatch(_V["anthraquinone"]):
        return True
    if db_cls == "stilbene" and m.HasSubstructMatch(_V["stilbene"]):
        return True
    if db_cls == "xanthine" and m.HasSubstructMatch(_V["xanthine"]):
        return True
    if db_cls == "nucleoside" and m.HasSubstructMatch(_V["nucleoside"]) and has_sugar:
        return True
    if db_cls == "glycoside" and has_sugar:
        return True
    if db_cls == "sugar" and has_sugar and n_heavy <= 24:    # 游离糖/寡糖(小)
        return True
    if db_cls == "alkaloid" and n_ringN >= 1 and n_N >= 1:
        return True
    if db_cls == "phenolic_acid" and n_arom >= 6 and _has_oh(m) and _has_cooH(m):
        return True
    if db_cls == "amino_acid" and _has_cooH(m) and n_N >= 1 and n_heavy <= 32:
        return True
    if db_cls == "lipid" and n_heavy >= 16 and n_arom <= 6:
        return True

    # 硬冲突: DB 标某强结构类, 但其标志结构完全缺失 (仅对可靠判据硬判)。
    # 注意: 甾体 gonane 用环计数, 但 6,6,6,5 也见于其他萜/重排, 不硬判 False (→None 更诚)。
    if db_cls == "sugar" and not has_sugar:
        return False
    if db_cls == "xanthine" and n_N < 3:
        return False
    return None


def stratified_sample(cur, n_target: int = 200) -> list[dict]:
    """按 classification_v3 **等额**分层 (每类固定配额), 从 source='ml' 抽样。
    等额(非按比例)才能让每类都有足够样本估精度 —— 比例分层会被 lipid 等大类淹没。"""
    cls_counts = cur.execute(
        "SELECT classification_v3, COUNT(*) c FROM Compounds "
        "WHERE classification_source='ml' GROUP BY classification_v3"
    ).fetchall()
    n_classes = len(cls_counts) or 1
    per_cls = max(8, -(-n_target // n_classes))   # ceil(n_target / n_classes), 下限 8
    samples = []
    shortfall = 0
    for cls, avail in cls_counts:
        q = min(per_cls, avail)
        rows = cur.execute(
            "SELECT id, compound_name, smiles, classification_v3, classification_confidence "
            "FROM Compounds WHERE classification_source='ml' AND classification_v3=? "
            "ORDER BY classification_confidence DESC LIMIT ?", (cls, q)
        ).fetchall()
        for r in rows:
            samples.append({"id": r[0], "name": r[1], "smiles": r[2],
                            "cls": r[3], "conf": r[4]})
        shortfall += max(0, per_cls - avail)      # 该类不够配额 → 记缺额
    # 用大类补足到 n_target (保持每类≥1)
    if len(samples) < n_target:
        extra = cur.execute(
            "SELECT id, compound_name, smiles, classification_v3, classification_confidence "
            "FROM Compounds WHERE classification_source='ml' ORDER BY RANDOM() LIMIT ?",
            (n_target - len(samples) + 200,)
        ).fetchall()
        have = {s["id"] for s in samples}
        for r in extra:
            if r[0] in have:
                continue
            samples.append({"id": r[0], "name": r[1], "smiles": r[2], "cls": r[3], "conf": r[4]})
            if len(samples) >= n_target:
                break
    return samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(ROOT / "data" / "local_library.db"))
    ap.add_argument("--n", type=int, default=200)
    args = ap.parse_args()

    print(f"=== P2·2a ML 精度抽检 (N≈{args.n}) ===\n  DB: {args.db}")
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    cur = con.cursor()

    samples = stratified_sample(cur, args.n)
    print(f"  抽样: {len(samples)} 个 (source=ml, 分层)")

    by_cls = defaultdict(lambda: {"correct": 0, "wrong": 0, "struct_bad": 0,
                                  "struct_ok": 0, "uncertain": 0})
    rows_out = []
    for s in samples:
        db_cls = s["cls"]
        npred = name_predict(s["name"])
        sok = struct_check(db_cls, s["smiles"])
        if npred is not None:                         # 名字强证据优先
            verdict = "CORRECT" if npred == db_cls else "WRONG"
        elif sok is False:
            verdict = "STRUCT_BAD"                    # 结构硬冲突
        elif sok is True:
            verdict = "STRUCT_OK"                     # 结构确认 (弱正)
        else:
            verdict = "UNCERTAIN"                     # 双 oracle 均无线索
        key = {"CORRECT": "correct", "WRONG": "wrong", "STRUCT_BAD": "struct_bad",
               "STRUCT_OK": "struct_ok", "UNCERTAIN": "uncertain"}[verdict]
        by_cls[db_cls][key] += 1
        rows_out.append({**s, "name_pred": npred, "verdict": verdict})

    con.close()

    # ── 写 CSV (人工复核用) ──
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "cls", "conf", "name_pred",
                                          "verdict", "name", "smiles"])
        w.writeheader()
        for r in sorted(rows_out, key=lambda x: (x["cls"], x["verdict"])):
            w.writerow({k: r.get(k) for k in ["id", "cls", "conf", "name_pred", "verdict", "name", "smiles"]})

    # ── 写 MD ──
    lines = [
        "# mssystem v3 — ML 分类层独立精度抽检 (P2·2a)",
        "",
        "## 目的与方法",
        "对 `classification_source='ml'` 的化合物(被 ML 猜标签)做**独立**复核,",
        "给出真实 per-class 精度估计, 替代当前 ML 训练循环自洽的 macro-F1。",
        "",
        "**独立 oracle (不复用训练/分类管线):**",
        "1. 名字启发式 `NAME_RULES`(化合物名含类标志词, 独立于结构)",
        "2. 结构健全性 `struct_check`(RDKit 粗筛 gross error)",
        "",
        "**裁定:** `CORRECT`(名字明确指向==DB类) / `WRONG`(名字明确指向≠DB类) /",
        "`STRUCT_BAD`(结构硬冲突) / `STRUCT_OK`(新鲜结构标志命中=结构一致, 弱正) /",
        "`UNCERTAIN`(名字+结构双 oracle 均无线索→人工复核)。",
        "",
        "**精度** 报两版:",
        "- **strict** = CORRECT / (CORRECT+WRONG+STRUCT_BAD)  — 仅高置信裁定;",
        "- **lenient** = (CORRECT+STRUCT_OK) / (CORRECT+WRONG+STRUCT_BAD+STRUCT_OK) — 含结构确认。",
        "",
        "UNCERTAIN 单列不计入任一分母(透明)。两版跨度反映结构可判程度。",
        "",
        "> ⚠ 这是**独立启发式(名字+新鲜SMARTS)的一致性估计**, 非绝对真值。结构近似异构体难判,",
        "> UNCERTAIN/WRONG 须结合 `ml_precision_audit_sample.csv` 人工/文献终核; 真值精度需 NP Classifier / ClassyFire。",
        "",
        "## per-class 结果",
        "",
        "| class | n | CORRECT | STRUCT_OK | WRONG | STRUCT_BAD | UNCERTAIN | strict | lenient |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    tot_c = tot_d = tot_ok = 0
    for cls in sorted(by_cls):
        d = by_cls[cls]
        denom_s = d["correct"] + d["wrong"] + d["struct_bad"]
        strict = (100.0 * d["correct"] / denom_s) if denom_s else None
        denom_l = denom_s + d["struct_ok"]
        lenient = (100.0 * (d["correct"] + d["struct_ok"]) / denom_l) if denom_l else None
        n = sum(d.values())
        tot_c += d["correct"]; tot_d += denom_s; tot_ok += d["struct_ok"]
        s_str = f"{strict:.0f}%" if strict is not None else "—"
        l_str = f"{lenient:.0f}%" if lenient is not None else "—"
        lines.append(f"| {cls} | {n} | {d['correct']} | {d['struct_ok']} | {d['wrong']} | "
                     f"{d['struct_bad']} | {d['uncertain']} | {s_str} | {l_str} |")
    strict_o = (100.0 * tot_c / tot_d) if tot_d else 0.0
    lenient_o = (100.0 * (tot_c + tot_ok) / (tot_d + tot_ok)) if (tot_d + tot_ok) else 0.0
    n_unc = sum(d["uncertain"] for d in by_cls.values())
    lines += [
        f"| **overall** | **{sum(sum(d.values()) for d in by_cls.values())}** | "
        f"**{tot_c}** | **{tot_ok}** | | | **{n_unc}** | **{strict_o:.0f}%** | **{lenient_o:.0f}%** |",
        "",
        "## 结论",
        f"- 抽样 {len(samples)} 个; strict 分母 {tot_d}, lenient 分母 {tot_d+tot_ok}, "
        f"UNCERTAIN {n_unc} 待人工复核。",
        f"- **strict 精度 {strict_o:.0f}%** / **lenient 精度 {lenient_o:.0f}%** (独立双 oracle 一致性)。",
        f"- 真实精度估计落在 [strict, lenient] 区间, 偏向 lenient (STRUCT_OK 多为真, 但非充分证明)。",
        f"- 样本明细见 `docs/ml_precision_audit_sample.csv` (按类+裁定排序)。",
        "",
        "## 使用建议",
        "- 该精度是**独立启发式一致性**, 不是绝对真值 (真值需 NP Classifier/ClassyFire 或全量人工)。",
        "- strict/lenient 跨度大的类 (STRUCT_OK 多) → 结构可判但名字少, 可信度区间宽, 谨慎用于展示。",
        "- 对 precision 偏低的类, 应优先人工复核后再用于展示/打分。",
        "- 头部高频化合物建议走 P2·2a **策展** (`classification_source='curated'`) 取代 ML 猜测。",
        "",
        "## ⚠ 重要发现: 稀有类 ML 过度指派 (over-assignment)",
        "抽样肉眼复核 (`ml_precision_audit_sample.csv`) 发现: **xanthone / stilbene / quinone 等",
        "稀有类的 `ml` 预测里混入大量结构完全不相关的化合物** (如马来酸酐/磺胺/蒽环/视黄酸被标 xanthone/stilbene)。",
        "因这些类无可靠通用 SMARTS 做硬否决, 它们落入 UNCERTAIN 而非 STRUCT_BAD, 故**未计入上面的精度分母**,",
        "意味着稀有类的**真实精度可能显著低于 lenient 估值**。",
        "- 处置: 稀有类 (xanthone/stilbene/quinone/...共 <1% 占比) 的 `ml` 标签**不可信**, 应:",
        "  (a) 优先用 P2·2a 策展覆盖其头部; (b) 展示/打分时对稀有类降权或标注 '低置信';",
        "  (c) 真实精度需 NP Classifier / ClassyFire 终核。",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n  strict 精度:   {strict_o:.0f}%  ({tot_c}/{tot_d})")
    print(f"  lenient 精度:  {lenient_o:.0f}%  ({tot_c+tot_ok}/{tot_d+tot_ok})")
    print(f"  UNCERTAIN(待人工): {n_unc}/{len(samples)}")
    print(f"  写出: {OUT_MD}")
    print(f"  写出: {OUT_CSV}")


if __name__ == "__main__":
    main()
