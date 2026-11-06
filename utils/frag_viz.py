"""碎裂可视化工具: RDKit 质量匹配法断键推断(无需 Docker, 离线方案)。

当 CFM-ID 产物结构缓存不可用时, 用 RDKit 在前体 2D 结构上标注"可能断裂的键":
枚举非环键 → 模拟断键 → 片段精确质量匹配中性丢失 Δm → 命中则高亮该键。
与 web/app.py::_fragmentation_svg 同源逻辑, 抽出来给 server.py 网页接口复用。
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw


def rdkit_bond_highlight_svg(smiles, losses, size=(300, 220)):
    """RDKit 质量匹配法碎裂可视化: 基于中性丢失 Δm 推断断键位置, 高亮前体上这些键。

    Args:
        smiles: 前体 SMILES
        losses: 中性丢失 Δm 列表 (float, Da)
        size: SVG 尺寸
    Returns:
        (svg_str, n_highlights) — svg 为空串表示无法生成
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "", 0
        clean = []
        for l in (losses or []):
            try:
                v = float(l)
            except Exception:
                continue
            if v > 0:
                clean.append(round(v, 2))
        clean = list(dict.fromkeys(clean))[:15]
        if not clean:
            return "", 0
        total = Descriptors.ExactMolWt(mol)
        hl = []
        for bond in mol.GetBonds():
            if bond.IsInRing():
                continue
            em = Chem.RWMol(mol)
            try:
                em.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                frags = Chem.GetMolFrags(em, asMols=True)
            except Exception:
                continue
            hit = False
            for f in frags:
                fm = Descriptors.ExactMolWt(f)
                for L in clean:
                    if abs(fm - L) <= 0.8 or abs(total - fm - L) <= 0.8:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                hl.append(bond.GetIdx())
        hl = list(dict.fromkeys(hl))
        AllChem.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DSVG(size[0], size[1])
        drawer.DrawMolecule(mol, highlightAtoms=[], highlightBonds=hl)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().split("?>", 1)[-1]
        svg = svg.replace("<svg ", '<svg preserveAspectRatio="xMidYMid meet" style="width:100%;height:auto" ', 1)
        return svg, len(hl)
    except Exception:
        return "", 0
