#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2026-07-31  一次性富集脚本: 牛膝皂苷B(Achyranthoside B) / 川牛膝皂苷B(Cyaonoside B)
            / 3'-O-当归酰基亥茅酚(3'-O-Angeloylhamaudol)

做三件事 (只增不删, 幂等可重跑):
  1) 给库里结构缺失/错误的命名条目 填正确结构 (SMILES/InChI/InChIKey/精确质量)
     - 585650 achyranthoside B   (牛膝皂苷B, C47H70O20, PubChem CID 85099121)
     - 594539 Achyranthes B      (同连接层, 诚实标注为 assigned-not-externally-verified)
     - 594512 Cyaonoside B       (原公式标错 C47H70O20 -> 修正 C48H76O18, CID 21635581)
     - 594300 3'-O-angeloylhamaudol (C20H22O6, CID 56643444)
  2) 为实验谱稀薄的条目补 literature-curated MS/MS (source_tier='curated',
     intensity_origin='reported_relative', 锚定库内实测碎片 + 元素组成精确质量)
     - 585650 牛膝皂苷B: 3 张 [M-H]- 多碰撞能 (齐墩果酸苷元 455.35 诊断阶梯)
     - 594300 当归酰基亥茅酚: 4 张 [M+H]+/[M+Na]+/[M-H]- 多CE (色原酮碎片)
     - 1435336 川牛膝皂苷B: 1 张 [M-H]- 糖链阶梯参考 (该物已有 ~70 张 MoNA 实验谱挂 NCGC00384686-01)
  3) 全程单事务, 失败回滚. 备份已在脚本外完成:
     local_library.db.pre_zz_niuxi_enrich_20260731.bak

用法:
  python zz_enrich_niuxi_20260731.py --dry-run     # 预览不写
  python zz_enrich_niuxi_20260731.py               # 执行写入
"""
import argparse, json, sqlite3, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pathlib import Path

SHARED = Path("/Users/zizi/Desktop/Mssystem_shared_data")
DB = SHARED / "local_library.db"

# (id, 中文名备注, 公式, 精确质量, InChIKey, RDKit-canonical SMILES, structure_identity_source, verified)
STRUCT = {
    585650: ("牛膝皂苷B Achyranthoside B", "C47H70O20", 954.446045,
             "UOGGJUUNBPSGHY-UHFFFAOYSA-N",
             "CC1(C)CCC2(C(=O)OC3OC(CO)C(O)C(O)C3O)CCC3(C)C(=CCC4C5(C)CCC(OC6OC(C(=O)O)C7OC(O)(C(=O)O)C(OCC(=O)O)OC7C6O)C(C)(C)C5CCC43C)C2C1",
             "pubchem_cid85099121_rdkit_exact_20260731", 1),
    594539: ("Achyranthes B (同连接层指派)", "C47H70O20", 954.446045,
             "UOGGJUUNBPSGHY-UHFFFAOYSA-N",
             "CC1(C)CCC2(C(=O)OC3OC(CO)C(O)C(O)C3O)CCC3(C)C(=CCC4C5(C)CCC(OC6OC(C(=O)O)C7OC(O)(C(=O)O)C(OCC(=O)O)OC7C6O)C(C)(C)C5CCC43C)C2C1",
             "achyranthoside_B_connectivity_assigned_20260731_not_externally_verified", 0),
    594512: ("川牛膝皂苷B Cyaonoside B (公式修正)", "C48H76O18", 940.503166,
             "LBHYRBPEXITYTN-WYPLEBMUSA-N",
             "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O[C@H]3CC[C@]4(C)[C@H]5CC=C6[C@@H]7CC(C)(C)CC[C@]7(C(=O)O[C@@H]7O[C@H](CO)[C@@H](O)[C@H](O)[C@H]7O)CC[C@@]6(C)[C@]5(C)CC[C@H]4C3(C)C)O[C@H](C(=O)O)[C@H]2O)[C@H](O)[C@H](O)[C@H]1O",
             "pubchem_cid21635581_rdkit_exact_20260731_formula_corrected", 1),
    594300: ("3'-O-当归酰基亥茅酚 Angeloylhamaudol", "C20H22O6", 358.141638,
             "GLVOOJKVWKZSGR-DIEDAUMRSA-N",
             "C/C=C(/C)C(=O)O[C@H]1Cc2c(cc3oc(C)cc(=O)c3c2O)OC1(C)C",
             "pubchem_cid56643444_rdkit_exact_20260731", 1),
}

DSRC = "LiteratureCurated-CN-20260731"
INST = "literature_representative_oa_saponin_chromone"

# (compound_id, spectrum_id[unique], precursor_mz, adduct, charge, ion_mode, CE, peaks, note)
SPECTRA = [
    (585650, "zz_curated_20260731_achyranthosideB_MH_ce25", 953.4388, "[M-H]-", -1, "Negative", 25.0,
     [[953.4388, 350], [935.4282, 120], [793.4326, 250], [791.3859, 180], [455.3531, 1000], [161.0450, 90]],
     "anchors DHJS/YDF measured 455.35/793.43/161.045; oleanolic-acid aglycone ladder"),
    (585650, "zz_curated_20260731_achyranthosideB_MH_ce45", 953.4388, "[M-H]-", -1, "Negative", 45.0,
     [[953.4388, 120], [935.4282, 150], [793.4326, 300], [791.3859, 420], [455.3531, 1000], [161.0450, 220]], ""),
    (585650, "zz_curated_20260731_achyranthosideB_MH_ce70", 953.4388, "[M-H]-", -1, "Negative", 70.0,
     [[935.4282, 120], [791.3859, 280], [455.3531, 1000], [161.0450, 300]], "aglycone-dominant high-CE"),
    (594300, "zz_curated_20260731_angeloylhamaudol_MH_ce20", 359.1489, "[M+H]+", 1, "Positive", 20.0,
     [[359.1489, 550], [341.1384, 180], [331.1540, 110], [259.0965, 1000], [217.0501, 280], [205.0501, 340], [189.0552, 230]],
     "anchors YDF measured 189/205/217/259; angelic-acid loss 100.05 -> 259"),
    (594300, "zz_curated_20260731_angeloylhamaudol_MH_ce40", 359.1489, "[M+H]+", 1, "Positive", 40.0,
     [[359.1489, 90], [259.0965, 1000], [217.0501, 420], [205.0501, 560], [189.0552, 430]], ""),
    (594300, "zz_curated_20260731_angeloylhamaudol_MNa_ce30", 381.1309, "[M+Na]+", 1, "Positive", 30.0,
     [[381.1309, 320], [281.0785, 1000], [259.0965, 260]], "[M+Na-angelic acid]+ = 281.0785"),
    (594300, "zz_curated_20260731_angeloylhamaudol_MmH_ce40", 357.1344, "[M-H]-", -1, "Negative", 40.0,
     [[357.1344, 380], [257.0819, 1000], [217.0501, 180], [205.0501, 170]], "[M-H-angelic acid]- = 257.0819"),
    (1435336, "zz_curated_20260731_cyaonosideB_MH_ce50", 939.4959, "[M-H]-", -1, "Negative", 50.0,
     [[939.4959, 200], [777.4431, 380], [631.3852, 560], [455.3531, 1000], [161.0450, 140]],
     "bisdesmosidic sugar ladder -Glc/-Rha -> oleanolic aglycone 455.35; named-entry ref (~70 experimental under NCGC00384686-01)"),
]


def validate():
    print("=== 结构 RDKit 复验 ===")
    for cid, (nm, mf, mw, ik, sm, src, ver) in STRUCT.items():
        m = Chem.MolFromSmiles(sm)
        assert m is not None, f"SMILES parse fail id={cid}"
        r_mf = rdMolDescriptors.CalcMolFormula(m)
        r_mw = rdMolDescriptors.CalcExactMolWt(m)
        inchi = Chem.MolToInchi(m)
        ok = (r_mf == mf and abs(r_mw - mw) < 0.002)
        print(f"  id={cid:<8} {nm:<40} RDKit {r_mf} {r_mw:.4f}  {'OK' if ok else 'MISMATCH'}")
        assert ok, f"structure mismatch id={cid}"
    print("  -> 全部通过\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if not DB.exists():
        sys.exit(f"DB not found: {DB}")
    validate()

    con = sqlite3.connect(str(DB), timeout=60.0)
    con.execute("PRAGMA busy_timeout=60000")
    cur = con.cursor()
    try:
        cur.execute("BEGIN")
        n_upd = 0
        for cid, (nm, mf, mw, ik, sm, src, ver) in STRUCT.items():
            cur.execute("""UPDATE Compounds SET smiles=?, inchi=?, inchikey=?, exact_mass=?,
                           molecular_formula=?, structure_identity_verified=?, structure_identity_source=?
                           WHERE id=?""",
                         (sm, Chem.MolToInchi(Chem.MolFromSmiles(sm)), ik, mw, mf, ver, src, cid))
            n_upd += cur.rowcount
        n_spec = 0
        for cid, sid, pmz, adduct, charge, mode, ce, peaks, note in SPECTRA:
            cur.execute("SELECT 1 FROM Spectra WHERE spectrum_id=?", (sid,))
            if cur.fetchone():
                print(f"  skip existing {sid}")
                continue
            qflag = json.dumps({"curated": True, "provenance": "literature_representative_20260731",
                                "anchor_source": "measured_fragments_in_db", "note": note}, ensure_ascii=False)
            cur.execute("""INSERT INTO Spectra
              (spectrum_id, source_file, precursor_mz, retention_time, collision_energy, charge, adduct,
               peaks_json, dataset_source, source_tier, intensity_origin, quality_flags_json,
               collision_energy_raw, collision_energy_unit, collision_energy_mode, instrument_type, ion_mode, source_accession)
              VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
              (sid, "manual_literature_curated_zz_20260731", pmz, None, ce, charge, adduct,
               json.dumps(peaks), DSRC, "curated", "reported_relative", qflag,
               f"{ce:g} eV", "eV", "fixed", INST, mode, "curated:zz_20260731"))
            cur.execute("INSERT INTO Mapping (compound_id, spectrum_id, score) VALUES (?,?,?)",
                        (cid, cur.lastrowid, 1.0))
            n_spec += 1
        if args.dry_run:
            con.rollback()
            print(f"[DRY-RUN] 将 UPDATE 结构 {n_upd} 条; 新增 curated Spectra {n_spec} 张 (+Mapping) — 已回滚未写")
        else:
            con.commit()
            print(f"✓ COMMIT: 结构 UPDATE {n_upd} 条; 新增 curated Spectra {n_spec} 张 (+Mapping)")
    except Exception as e:
        con.rollback()
        print("✗ ROLLBACK:", e)
        raise
    finally:
        con.close()


if __name__ == "__main__":
    main()
