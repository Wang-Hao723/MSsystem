"""
mssystem v3 Phase 3.2-3.3 — 构建统一 schema 的天然产物裂解规则库 v3
====================================================================
合并:
  (A) 经验挖掘 (kb_v3/empirical_rules.json) — 数值诊断离子/中性丢失 + 频率
  (B) 文献策展知识 (下方 CURATED) — 骨架/官能团/生源途径/植物科属 4 维 + 机制

统一 schema (全数值, 可机器匹配), 输出 kb_v3/fragmentation_rules_v3.json。
根治 v2 KB 的两类缺陷: schema 不一致 + 符号 m/z。

cd /Users/zizi/Desktop/Mssystem/mssystem_v3
python3 kb_v3/build_kb_v3.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EMPIRICAL = ROOT / "kb_v3" / "empirical_rules.json"
OUT = ROOT / "kb_v3" / "fragmentation_rules_v3.json"

# ─────────────────────────────────────────────────────────────────────
# (B) 文献策展知识 — 基于专业书籍/文献 + 8 篇本地 PDF + 化学规律
# 来源: 《中药质谱分析》(刘淑莹 2012); Cuyckens & Claeys 2004 (flavonoid MS);
#        ClassyFire 本体; GNPS/MassBank; Spectraverse (Anal. Chem. 2026);
#        质谱裂解特征/ PDF (β-carboline, Sceletium, NP structure determination)
# ─────────────────────────────────────────────────────────────────────
CURATED = {
    "flavonoid": {
        "class_en": "Flavonoids", "core_skeleton": "C6-C3-C6 (2-苯基色原酮)",
        "biosynthesis": "莽草酸途径 → 苯丙素 → 查耳酮 → 各亚型",
        "subtypes": {
            "flavone": "芹菜素/木犀草素/黄芩素 (C2=C3, C4=O)",
            "flavonol": "槲皮素/山柰酚/杨梅素 (C3-OH)",
            "flavanone": "柚皮素 (C2-C3 饱和)",
            "flavan-3-ol": "儿茶素 (无C4酮)",
            "isoflavone": "葛根素/大豆苷元 (B环在C3)",
            "chalcone": "C6-C3-C6 开环",
            "anthocyanidin": "花色素 (flavylium 阳离子)",
        },
        "mechanisms": [
            {"name": "RDA 逆Diels-Alder (C环1,3断裂)", "products": "1,3A⁻ (m/z 151 for 5,7-diOH) + 1,3B⁻",
             "note": "黄酮/黄酮醇诊断; A环碎片常为基峰", "ref": "Cuyckens & Claeys 2004; Fabre 2001"},
            {"name": "糖苷键断裂脱糖", "products": "[M-H-162]⁻ (hexose), [M-H-308]⁻ (rutinose)",
             "note": "O-苷先脱糖得苷元; C-苷难断裂", "ref": "March 2006"},
            {"name": "CO/CO2 丢失", "products": "[A-CO]⁻, [A-CO2]⁻", "ref": "—"},
        ],
        "functional_group_effects": {
            "3-OH": "黄酮醇 vs 黄酮区分 (3-OH 引导不同 RDA)",
            "5,7-diOH": "A环 m/z 151/153 诊断",
            "3',4'-diOH(B环)": "m/z 135/137 诊断",
            "O-glycoside": "脱 162(hexose)/146(deoxyhexose)/308(rutinose)",
            "C-glycoside": "抗脱糖, 丢失 120/90 (糖环内碎裂)",
            "prenyl": "异戊烯基黄酮, 丢失 69 (C5H9) / 56",
            "methoxy": "丢失 15 (CH3) / 32 (CH3OH)",
        },
        "plant_families": ["Asteraceae(菊科)", "Fabaceae(豆科)", "Lamiaceae(唇形科)", "Rutaceae(芸香科)"],
    },
    "alkaloid": {
        "class_en": "Alkaloids", "core_skeleton": "含氮碱性天然产物 (生源: 氨基酸衍生)",
        "biosynthesis": "氨基酸途径 (Trp/Tyr/Phe/Lys → 各类生物碱)",
        "subtypes": {
            "isoquinoline": "小檗碱/吗啡/青藤碱 (苄基异喹啉)",
            "indole": "β-卡波林/长春新碱/士的宁",
            "quinoline": "奎宁/喜树碱",
            "quinolizidine": "苦参碱/金雀花碱 (豆科)",
            "pyrrolizidine": "千里光碱 (菊科, 肝毒)",
            "tropane": "阿托品/东莨菪碱 (茄科)",
            "diterpene_alkaloid": "乌头碱/次乌头碱 (毛茛科, 极毒)",
            "purine(xanthine)": "咖啡因/可可碱 (嘌呤生物碱)",
            "steroidal_alkaloid": "茄碱/番茄碱 (茄科)",
        },
        "mechanisms": [
            {"name": "α-断裂 (N旁C-C)", "products": "亚胺离子/iminium (稳定)",
             "note": "胺类生物碱α断裂产生共振稳定碳正离子", "ref": "《中药质谱分析》2012"},
            {"name": "RDA (含N环)", "products": "依赖骨架", "ref": "—"},
            {"name": "脱水/脱甲基", "products": "[M+H-18]+, [M+H-15]+", "ref": "—"},
        ],
        "functional_group_effects": {
            "ester": "乌头碱型: 丢失 60 (AcOH) / 76 (benzoyl 相关)",
            "methoxy": "丢失 31/32",
            "N-oxide": "丢失 16 (O)",
        },
        "plant_families": ["Ranunculaceae(毛茛科=乌头)", "Papaveraceae(罂粟科)", "Solanaceae(茄科)", "Fabaceae(豆科=苦参)", "Rubiaceae(茜草科)"],
    },
    "terpenoid": {
        "class_en": "Terpenoids", "core_skeleton": "异戊二烯(C5)单元聚合",
        "biosynthesis": "MVA + MEP 途径 (C5 → C10/C15/C20/C30/C40)",
        "subtypes": {
            "monoterpene(C10)": "环烯醚萜前体/薄荷烷/龙脑 (丢失 69 C5H9)",
            "sesquiterpene(C15)": "倍半萜内酯(白术内酯/土木香内酯)",
            "diterpene(C20)": "丹参酮/穿心莲内酯/雷公藤内酯",
            "triterpene(C30)": "齐墩果烷/乌苏烷/达玛烷型 (皂苷元)",
            "tetraterpene(C40)": "类胡萝卜素",
        },
        "mechanisms": [
            {"name": "逆Diels-Alder", "products": "异戊二烯片段 m/z 91/105/119/131", "ref": "—"},
            {"name": "脱水 (含OH)", "products": "[M+H-18]+ / [M-H2O-H]-", "ref": "—"},
            {"name": "脱甲基", "products": "-15", "ref": "—"},
        ],
        "functional_group_effects": {
            "lactone": "倍半萜内酯: 丢失 44/56 (内酯开环)",
            "ester": "乌头酯型生物碱同",
            "OH": "脱水 18",
            "glycoside": "皂苷先脱糖",
        },
        "plant_families": ["Asteraceae(菊科=倍半萜内酯/白术)", "Lamiaceae(唇形科=丹参酮)", "Lauraceae(樟科=单萜)", "Apocynaceae(夹竹桃科)"],
    },
    "iridoid": {
        "class_en": "Iridoids & Secoiridoids", "core_skeleton": "环戊烷[c]吡喃 (monoterpene)",
        "biosynthesis": "MEP途径 → 环烯醚萜苷",
        "subtypes": {"iridoid_glycoside": "京尼平苷/梓醇/桃叶珊瑚苷/车叶草苷",
                     "secoiridoid": "橄榄苦苷/龙胆苦苷/番木鳖苷"},
        "mechanisms": [
            {"name": "脱糖得苷元", "products": "[M-H-162]⁻ (hexose)", "ref": "—"},
            {"name": "苷元脱 CO/CO2/H2O", "products": "m/z 145/127/109 系列", "ref": "《中药质谱分析》"},
        ],
        "functional_group_effects": {"glycoside": "脱 162", "ester": "脱 60 (AcOH) 京尼平苷", "lactone": "环戊吡喃开环"},
        "plant_families": ["Scrophulariaceae(玄参科=地黄/梓醇)", "Caprifoliaceae(忍冬科=金银花)",
                           "Gentianaceae(龙胆科=龙胆苦苷)", "Oleaceae(木犀科=橄榄苦苷)", "Rubiaceae(茜草科=车叶草苷)"],
    },
    "saponin": {
        "class_en": "Saponins (triterpene/steroid glycosides)", "core_skeleton": "三萜/甾体苷元 + 糖链",
        "biosynthesis": "MVA途径 (三萜) + 糖基化",
        "subtypes": {"triterpene_saponin": "人参皂苷/甘草酸/柴胡皂苷 (达玛烷/齐墩果烷型)",
                     "steroid_saponin": "薯蓣皂苷/重楼皂苷"},
        "mechanisms": [
            {"name": "逐糖级联丢失", "products": "[M-H-162]⁻→[-162]⁻→...→苷元 (m/z 455/470 齐墩果酸等)",
             "note": "糖链序列可由级联丢失质量读出", "ref": "—"},
        ],
        "functional_group_effects": {"hexose": 162, "deoxyhexose": 146, "glucuronic_acid": 176,
                                     "pentose": 132, "acyl": "乙酰糖 +42/+60"},
        "plant_families": ["Araliaceae(五加科=人参皂苷)", "Fabaceae(豆科=甘草酸)", "Liliaceae(百合科=重楼)",
                           "Caryophyllaceae(石竹科=牛膝皂苷)"],
    },
    "phenolic_acid": {
        "class_en": "Phenolic acids", "core_skeleton": "酚环 + 羧基 (苯甲酸/肉桂酸型)",
        "biosynthesis": "莽草酸途径",
        "subtypes": {"hydroxybenzoic": "没食子酸/原儿茶酸/水杨酸/香草酸",
                     "hydroxycinnamic": "咖啡酸/阿魏酸/芥子酸/p-香豆酸/绿原酸(奎宁酸酯)"},
        "mechanisms": [
            {"name": "脱羧 (neg)", "products": "[M-H-44]⁻ (CO2)", "ref": "—"},
            {"name": "脱水 (pos)", "products": "[M+H-18]+", "ref": "—"},
            {"name": "绿原酸酯键断裂", "products": "咖啡酸 m/z 179/180 + 奎宁酸 m/z 191", "ref": "Clifford 2005"},
        ],
        "functional_group_effects": {"methoxy": "阿魏酸/芥子酸 脱 CH3(15)/CH3OH(32)",
                                     "quinic_acid_ester": "绿原酸系列 (3-/4-/5-咖啡酰奎宁酸异构)"},
        "plant_families": ["Lamiaceae(唇形科=丹参酚酸)", "Asteraceae(菊科=绿原酸)", "Caprifoliaceae(忍冬科=绿原酸)"],
    },
    "coumarin": {
        "class_en": "Coumarins", "core_skeleton": "苯并吡喃-2-酮 (C6-C3 内酯)",
        "biosynthesis": "莽草酸途径 (phenylpropanoid)",
        "subtypes": {"simple": "伞形酮/七叶苷", "furanocoumarin": "补骨脂素/白芷素 (芸香科)",
                     "pyranocoumarin": "邪蒿素", "isoprenylated": "紫花前胡苷"},
        "mechanisms": [
            {"name": "脱 CO2 (neg)", "products": "[M-H-44]⁻", "ref": "—"},
            {"name": "内酯环碎裂", "products": "[M-28] CO 丢失", "ref": "—"},
            {"name": "脱糖 (苷)", "products": "-162", "ref": "—"},
        ],
        "functional_group_effects": {"prenyl": "异戊烯基香豆素 丢失 56/69", "methoxy": "脱 15"},
        "plant_families": ["Rutaceae(芸香科=呋喃香豆素/白芷)", "Apiaceae(伞形科=前胡/当归)", "Fabaceae"],
    },
    "lignan": {
        "class_en": "Lignans", "core_skeleton": "C6-C3 二聚 (8-8' 连接的苯丙素)",
        "biosynthesis": "莽草酸途径 → 苯丙素二聚",
        "subtypes": {"dibenzylbutyrolactone": "牛蒡苷元/罗汉松脂素",
                     "aryltetralin": "五味子素/鬼臼毒素", "furofuran": "松脂醇/松脂醇二葡萄糖苷"},
        "mechanisms": [
            {"name": "脱糖 (苷)", "products": "-162/-324", "ref": "—"},
            {"name": "脱甲基/亚甲二氧基", "products": "-15/-44", "ref": "—"},
            {"name": "苄基断裂", "products": "酚丙烷碎片", "ref": "—"},
        ],
        "functional_group_effects": {"methylenedioxy": "丢失 44 (CH2O+CO)", "methoxy": "脱 15"},
        "plant_families": ["Schisandraceae(五味子科)", "Asteraceae(菊科=牛蒡)", "Pinaceae"],
    },
    "quinone": {
        "class_en": "Quinones", "core_skeleton": "蒽醌/萘醌/苯醌",
        "biosynthesis": "醋酸-丙二酸途径 (聚酮) / 莽草酸",
        "subtypes": {"anthraquinone": "大黄素/芦荟大黄素 (蓼科)", "naphthoquinone": "紫草素/胡桃醌",
                     "benzoquinone": "辅酶Q类"},
        "mechanisms": [
            {"name": "脱 CO/CO2", "products": "[M-28]/[M-44]", "ref": "—"},
            {"name": "脱 H2O (含OH蒽醌)", "products": "-18", "ref": "—"},
            {"name": "脱甲基/羧基", "products": "大黄素型 -15/-CH3", "ref": "—"},
        ],
        "functional_group_effects": {"OH": "脱水", "COOH": "脱 CO2/CO 大黄酸型", "prenyl": "脱 56/69"},
        "plant_families": ["Polygonaceae(蓼科=大黄蒽醌)", "Boraginaceae(紫草科=紫草素萘醌)", "Leguminosae"],
    },
    "steroid": {
        "class_en": "Steroids", "core_skeleton": "环戊烷并全氢菲 (gonane, 6-6-6-5 四环)",
        "biosynthesis": "MVA途径 (角鲨烯 → 羊毛甾醇 → 甾醇)",
        "subtypes": {"phytosterol": "β-谷甾醇/豆甾醇", "ecdysteroid": "β-蜕皮甾酮",
                     "withanolide": "醉茄内酯 (茄科)", "cardiac_glycoside": "强心苷 (夹竹桃科, 甾体+丁烯内酯)"},
        "mechanisms": [
            {"name": "脱水", "products": "[M+H-18]+ (-H2O), 连续脱 -18/-18", "ref": "—"},
            {"name": "侧链断裂", "products": "甾体核碎片 m/z 255/271/273", "ref": "—"},
            {"name": "强心苷丁烯内酯丢失", "products": "特征 m/z 111 (butenolide)", "ref": "—"},
        ],
        "functional_group_effects": {"OH": "脱水 18 (多个OH多次脱)", "glycoside": "脱糖",
                                     "side_chain": "侧链断裂质量取决于C数"},
        "plant_families": ["Apocynaceae(夹竹桃科=强心苷)", "Solanaceae(茄科=甾体皂苷/醉茄内酯)", "Liliaceae"],
    },
    "xanthone": {
        "class_en": "Xanthones", "core_skeleton": "二苯并-γ-吡喃酮 (9H-xanthen-9-one)",
        "biosynthesis": "聚酮途径 (龙胆科) / 莽草酸",
        "subtypes": {"simple": "1,3,6,7-四羟基 (芒果苷元)", "prenylated": "异戊烯基山酮", "glycoside": "芒果苷"},
        "mechanisms": [{"name": "脱 CO", "products": "[M-28]", "ref": "—"}, {"name": "脱糖(苷)", "products": "-162", "ref": "—"},
                       {"name": "RDA (山酮环)", "products": "依赖取代", "ref": "—"}],
        "functional_group_effects": {"OH": "脱水", "prenyl": "脱 69"},
        "plant_families": ["Gentianaceae(龙胆科)", "Clusiaceae/Hypericaceae(金丝桃科)", "Mangifera(芒果)"],
    },
    "stilbene": {
        "class_en": "Stilbenes", "core_skeleton": "1,2-二苯乙烯 (C6-C2-C6)",
        "biosynthesis": "莽草酸+聚酮 (芪合成酶 STS)",
        "subtypes": {"resveratrol": "白藜芦醇/虎杖苷", "pterostilbene": "紫檀芪", "glycoside": "虎杖苷(polydatin)"},
        "mechanisms": [{"name": "脱糖(苷)", "products": "-162 → 苷元 m/z 227 (白藜芦醇)", "ref": "—"},
                       {"name": "脱 CO/H2O", "products": "[M-28]/[M-18]", "ref": "—"},
                       {"name": "逆Diels-Alder", "products": "m/z 143/107 (取代苯环)", "ref": "—"}],
        "functional_group_effects": {"OH": "脱水/重排", "glycoside": "脱 162"},
        "plant_families": ["Polygonaceae(蓼科=虎杖)", "Vitaceae(葡萄科=葡萄)", "Leguminosae(豆科)"],
    },
}

# ─────────────────────────────────────────────────────────────────────
# 植物科属化学特征 (3.3 dim 4) — 用户明确要求
# family → {characteristic_classes, marker_compounds, fragmentation_notes}
# ─────────────────────────────────────────────────────────────────────
PLANT_FAMILY = {
    "Ranunculaceae(毛茛科)": {"classes": ["alkaloid(二萜生物碱)"], "markers": ["乌头碱(aconitine)", "次乌头碱", "小檗碱(黄连)"],
                              "notes": "乌头碱型二萜生物碱极毒; 丢失 60(AcOH)/90/120(苯甲酰) 酯链; 炮制减毒"},
    "Apocynaceae(夹竹桃科)": {"classes": ["steroid(强心苷)", "alkaloid(吲哚)"], "markers": ["地高辛/洋地黄毒苷", "长春新碱"],
                              "notes": "强心苷=甾体+丁烯内酯, 特征 m/z 111; 长春花吲哚生物碱"},
    "Solanaceae(茄科)": {"classes": ["alkaloid(甾体/托品)", "steroid"], "markers": ["茄碱(solanine)", "东莨菪碱/阿托品"],
                          "notes": "甾体生物碱 (茄碱) + 托品生物碱 (莨菪烷型)"},
    "Rutaceae(芸香科)": {"classes": ["coumarin(呋喃香豆素)", "alkaloid(喹啉/呋喃)", "flavonoid"], "markers": ["补骨脂素", "白芷素", "黄柏碱(小檗碱)", "橙皮苷"],
                          "notes": "呋喃香豆素特征; 喹啉生物碱; 黄酮苷"},
    "Lamiaceae(唇形科)": {"classes": ["terpenoid(二萜)", "phenolic_acid(丹酚酸)"], "markers": ["丹参酮IIA", "丹酚酸B", "迷迭香酸"],
                           "notes": "丹参=二萜醌(丹参酮)+缩酚酸(丹酚酸); 酚酸丰富"},
    "Asteraceae(菊科)": {"classes": ["terpenoid(倍半萜内酯)", "phenolic_acid(绿原酸)"], "markers": ["白术内酯", "青蒿素", "绿原酸", "洋蓟素"],
                          "notes": "倍半萜内酯丰富; 绿原酸系列"},
    "Fabaceae(豆科)": {"classes": ["flavonoid(异黄酮)", "alkaloid(喹嗪啶)", "saponin(三萜皂苷)"], "markers": ["葛根素/大豆苷元", "苦参碱", "甘草酸"],
                        "notes": "异黄酮特征(豆科); 喹嗪啶生物碱(苦参); 三萜皂苷"},
    "Polygonaceae(蓼科)": {"classes": ["quinone(蒽醌)", "stilbene"], "markers": ["大黄素/大黄酸", "白藜芦醇/虎杖苷"],
                            "notes": "蒽醌(大黄); 二苯乙烯(虎杖)"},
    "Araliaceae(五加科)": {"classes": ["saponin(三萜皂苷)"], "markers": ["人参皂苷Rg1/Re/Rb1", "三七皂苷"],
                            "notes": "达玛烷型三萜皂苷, 逐糖级联丢失"},
    "Zingiberaceae(姜科)": {"classes": ["terpenoid(二芳基庚烷/二萜)"], "markers": ["姜黄素(curcuminoid)", "姜辣素"],
                             "notes": "二芳基庚烷类(姜黄素); 酚性二萜"},
    "Gentianaceae(龙胆科)": {"classes": ["iridoid(裂环环烯醚萜)", "xanthone"], "markers": ["龙胆苦苷", "芒果苷"],
                              "notes": "裂环环烯醚萜苷(龙胆苦苷); 山酮"},
    "Caprifoliaceae(忍冬科)": {"classes": ["iridoid", "phenolic_acid"], "markers": ["马钱苷/番木鳖苷", "绿原酸"],
                                "notes": "环烯醚萜苷 + 绿原酸 (金银花)"},
    "Scrophulariaceae(玄参科)": {"classes": ["iridoid"], "markers": ["梓醇(catalpol)", "胡黄连苷"],
                                   "notes": "环烯醚萜苷 (地黄/梓醇) — 用户杜仲苷/车叶草苷同类"},
    "Boraginaceae(紫草科)": {"classes": ["quinone(萘醌)"], "markers": ["紫草素(shikonin)"],
                              "notes": "萘醌类 (紫草), 异戊烯基侧链"},
    "Oleaceae(木犀科)": {"classes": ["iridoid(secoiridoid)"], "markers": ["橄榄苦苷(oleuropein)", "女贞苷"],
                          "notes": "裂环环烯醚萜苷 (女贞子/连翘)"},
    "Schisandraceae(五味子科)": {"classes": ["lignan"], "markers": ["五味子甲素/乙素/醇甲"],
                                  "notes": "木脂素 (五味子), 联苯环辛二烯型"},
}


def build():
    empirical = json.load(open(EMPIRICAL, encoding="utf-8")) if EMPIRICAL.exists() else {}
    classes = {}
    for cls, cur in CURATED.items():
        emp = empirical.get(cls, {})
        classes[cls] = {
            "class": cls,
            "class_en": cur["class_en"],
            "core_skeleton": cur["core_skeleton"],
            "biosynthesis_pathway": cur["biosynthesis"],
            "skeleton_subtypes": cur["subtypes"],                    # 维度1: 骨架亚型
            "functional_group_effects": cur["functional_group_effects"],  # 维度2: 官能团
            "fragmentation_mechanisms": cur["mechanisms"],            # 文献机制
            "empirical": {                                            # 维度数据: 数值规则
                pol: {
                    "n_spectra": d.get("n_spectra", 0),
                    "diagnostic_ions": d.get("diagnostic_ions", []),
                    "neutral_losses": d.get("neutral_losses", []),
                } for pol, d in emp.items()
            },
            "plant_families": cur["plant_families"],
        }
    kb = {
        "meta": {
            "title": "mssystem 天然产物裂解规则库 v3",
            "version": "3.0.0",
            "schema": "unified: class→{skeleton_subtypes, functional_group_effects, mechanisms, empirical{pos/neg{diagnostic_ions[mz/freq/n_support], neutral_losses[loss_da/freq/n_support/annotation]}}, plant_families}",
            "principle": "全数值 m/z/loss 可机器匹配; 经验挖掘(频率) ∪ 文献策展(机制/4维)",
            "sources": ["913k本地库经验挖掘", "《中药质谱分析》刘淑莹2012", "Cuyckens&Claeys2004",
                        "ClassyFire本体", "GNPS/MassBank",
                        "Spectraverse (Gupta et al., Anal. Chem. 2026; spectral-data curation)",
                        "质谱裂解特征/8篇PDF"],
            "n_classes": len(classes),
        },
        "classes": classes,
        "plant_family_chemistry": PLANT_FAMILY,
        "common_neutral_losses": {
            "18.0106": "H2O", "27.9949": "CO", "28.0061": "CH2=CH2", "31.9898": "S", "33.9877": "H2S",
            "43.9898": "CO2", "46.0055": "HCOOH", "57.0215": "CH3NCO", "59.9847": "COS",
            "79.9568": "SO3", "120.0423": "C7H4O2", "131.0375": "C8H3O?",
            "132.0423": "pentose", "146.0579": "deoxyhexose", "162.0528": "hexose",
            "176.0321": "glucuronic_acid", "180.0634": "hexose+H2O", "246.0618": "acetylhexose",
            "308.0905": "rutinose(hexose+deoxyhexose)",
        },
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(kb, f, ensure_ascii=False, indent=2)
    print(f"✅ v3 KB 已构建: {OUT}")
    print(f"   {len(classes)} 类 (4维: 骨架/官能团/生源/科属) + {len(PLANT_FAMILY)} 植物科属 + {len(kb['common_neutral_losses'])} 通用中性丢失")
    for cls, d in classes.items():
        emp = d["empirical"]
        pols = [f"{p}:{len(emp[p]['diagnostic_ions'])}di/{len(emp[p]['neutral_losses'])}nl" for p in emp if emp[p]["n_spectra"] > 0]
        print(f"   {cls:14s} {d['class_en']:30s} 经验[{', '.join(pols) or '无'}]")


if __name__ == "__main__":
    build()
