"""RDKit 结构式绘制：SMILES → 2D 分子图"""
import base64
import io
import logging

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, AllChem
    RDLogger.DisableLog("rdApp.*")
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available - structure images disabled")

# CoordGen 优先: 对复杂稠环/桥环体系(如吗啡烷类 morphinan、甾体、喹啉类生物碱)画出的 2D 朝向
# 更接近教科书标准; 默认 rdDepictor 常把这类多环画得扭曲别扭(用户反馈青藤碱结构"看着不对")。
try:
    from rdkit.Chem import rdCoordGen
    COORDGEN_AVAILABLE = True
except ImportError:
    COORDGEN_AVAILABLE = False


def _compute_depiction_coords(mol):
    """计算 2D 坐标: 优先 CoordGen (复杂稠环画法标准), 失败回退 rdDepictor。"""
    if COORDGEN_AVAILABLE:
        try:
            rdCoordGen.AddCoords(mol)
            return
        except Exception:
            pass
    AllChem.Compute2DCoords(mol)


def _sanitize_smiles(smiles):
    """Sanitize and canonicalize SMILES for reliable RDKit parsing.

    Handles common issues:
    - Invalid ring closure notation
    - Bad stereochemistry
    - Non-canonical forms that confuse some RDKit versions
    """
    if not smiles or not isinstance(smiles, str):
        return None
    s = smiles.strip()
    if not s or s.upper() in ("N/A", "NONE", "NULL", "-", "NA"):
        return None

    if not RDKIT_AVAILABLE:
        return s  # Return as-is; caller handles missing RDKit

    # First try: direct parse
    mol = Chem.MolFromSmiles(s)
    if mol is not None:
        try:
            return Chem.MolToSmiles(mol)  # Canonicalize
        except Exception:
            return s

    # Second try: parse without sanitization, then sanitize
    mol = Chem.MolFromSmiles(s, sanitize=False)
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol)
        except Exception:
            pass

    # Third try: remove stereochemistry and retry
    try:
        mol = Chem.MolFromSmiles(Chem.rdmolfiles.MolToSmiles(
            Chem.MolFromSmiles(s, sanitize=False), isomericSmiles=False
        ))
        if mol is not None:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass

    logger.debug("Could not sanitize SMILES: %s", s[:80])
    return None


def smiles_to_image_base64(smiles, size=(400, 300)):
    """SMILES → base64 PNG"""
    if not RDKIT_AVAILABLE or not smiles:
        return None

    sanitized = _sanitize_smiles(smiles)
    if sanitized is None:
        return None

    try:
        mol = Chem.MolFromSmiles(sanitized)
        if mol is None:
            return None
        _compute_depiction_coords(mol)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode()
    except Exception as e:
        logger.debug("Native PNG depiction unavailable: %s", e)
        try:
            # Some installed RDKit builds lack Cairo. Rasterize the same
            # RDKit SVG depiction with the available CairoSVG renderer.
            import cairosvg
            svg = smiles_to_svg(smiles, size=size)
            if not svg:
                return None
            data = cairosvg.svg2png(bytestring=svg.encode("utf-8"),
                                   output_width=size[0], output_height=size[1])
            return base64.b64encode(data).decode("ascii")
        except Exception as fallback_error:
            logger.debug("PNG rasterization failed: %s", fallback_error)
            return None


def smiles_to_data_uri(smiles, size=(400, 300)):
    """SMILES → 完整 data URI (PNG 优先, SVG 兜底).

    返回 'data:image/png;base64,...' 或 'data:image/svg+xml;base64,...',
    或 None (SMILES 无法解析时). RDKit 某些构建(如 base 环境 2025.03.6)
    未编译 Cairo → Draw.MolToImage(PNG) 抛 RuntimeError; 此时自动回退到
    MolDraw2DSVG(无需系统库), 浏览器/HTML 原生渲染 SVG. 因此本函数与运行
    环境是否有 Cairo 无关, 总能产出可用结构图.

    注意: smiles_to_image_base64 仅返回裸 PNG base64 (Cairo 缺失时为 None),
    供需要真实字节的调用方(如 st.image)使用; 需要"一定能渲染一张图"时用本函数.
    """
    if not RDKIT_AVAILABLE or not smiles:
        return None
    sanitized = _sanitize_smiles(smiles)
    if sanitized is None:
        return None
    try:
        mol = Chem.MolFromSmiles(sanitized)
        if mol is None:
            return None
        _compute_depiction_coords(mol)
        # PNG 优先 (需 Cairo)
        try:
            img = Draw.MolToImage(mol, size=size)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")
        except Exception:
            pass
        # SVG 兜底 (无系统依赖)
        drawer = Draw.MolDraw2DSVG(size[0], size[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return "data:image/svg+xml;base64," + base64.b64encode(svg.encode("utf-8")).decode("ascii")
    except Exception as e:
        logger.debug("Failed to generate image for SMILES %s: %s", smiles[:50], e)
        return None


def smiles_to_svg(smiles, size=(400, 300)):
    """SMILES → SVG 字符串"""
    if not RDKIT_AVAILABLE or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        _compute_depiction_coords(mol)
        drawer = Draw.MolDraw2DSVG(size[0], size[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def smiles_is_valid(smiles):
    """检查 SMILES 是否有效"""
    if not RDKIT_AVAILABLE or not smiles:
        return False
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False
