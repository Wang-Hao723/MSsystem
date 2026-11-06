"""CFM-ID 碎裂产物结构 — 通过 docker 调 fraggraph-gen, 解析碎片结构SMILES+精确质量, 带缓存。

用户已装 Docker + 拉了 kykeonanalytics/cfmid-arm64 (arm64原生, depth1约3s/分子)。
fraggraph-gen 输出: 每行 "id 质量 碎片SMILES结构" — 这就是质谱碎裂的产物结构(CFM-ID用
竞争碎裂模型+ML预测, 含rDA/重排等真机理, 不是RDKit单键断裂的H2O假象)。

用法:
  from utils.cfm_fragment import get_fragments
  frags = get_fragments("naringenin SMILES")  # [(mass, smiles), ...], 带缓存
"""
import json
import logging
import os
import subprocess
import uuid
from pathlib import Path

logger = logging.getLogger(__name__)

# arm64 原生镜像(快, ~3s/分子 depth1); 无则降级 wishartlab/cfmid(amd64慢~60s)
DOCKER_IMG = "kykeonanalytics/cfmid-arm64"
CACHE_DIR = Path(__file__).resolve().parent.parent / "output" / "cfm_fragments"


def run_fraggraph(smiles, depth=1, ionization="+", timeout=90):
    """调 docker fraggraph-gen, 返回 [(mass, frag_smiles), ...]。depth=1 快(3s), depth=2 慢(60s)。

    ionization: + (POS [M+H]) / - (NEG [M-H]) / * (EI)。
    返回按质量降序(前体在前)。
    """
    container_name = (
        f"mssystem-cfmid-{os.getpid()}-{uuid.uuid4().hex[:10]}"
    )
    try:
        out = subprocess.run(
            [
                "docker", "run", "--rm",
                "--name", container_name,
                "--label", "com.mssystem.component=cfmid-fraggraph",
                DOCKER_IMG, "fraggraph-gen",
             smiles, str(depth), ionization, "fragonly"],
            capture_output=True, text=True, timeout=timeout
        )
        if out.returncode != 0:
            logger.warning("fraggraph-gen failed (rc=%s): %s", out.returncode, out.stderr[:200])
            return []
        frags = []
        for line in out.stdout.splitlines():
            parts = line.split()
            # 格式: "id mass smiles [Intermediate Fragment]"
            if len(parts) >= 3 and parts[0].isdigit():
                try:
                    mass = float(parts[1])
                    smi = parts[2]
                    if smi and smi[0].isalpha():  # 合法SMILES开头
                        frags.append({"mass": mass, "smiles": smi})
                except (ValueError, IndexError):
                    continue
        frags.sort(key=lambda x: -x["mass"])  # 质量降序, 前体(最大)在前
        return frags
    except subprocess.TimeoutExpired:
        logger.warning("fraggraph-gen timeout for %s", smiles[:30])
        # Killing the ``docker run`` client does not stop the container.  The
        # old implementation consequently leaked hundreds of CPU-consuming
        # fraggraph containers.  Remove this exact, per-call named container.
        try:
            subprocess.run(
                ["docker", "rm", "-f", container_name],
                capture_output=True,
                text=True,
                timeout=15,
                check=False,
            )
        except Exception:
            logger.warning(
                "failed to remove timed-out CFM container %s",
                container_name,
            )
        return []
    except FileNotFoundError:
        logger.warning("docker not found — CFM碎裂不可用")
        return []
    except Exception as e:
        logger.warning("fraggraph-gen error: %s", e)
        return []


def get_fragments(smiles, depth=1, ionization="+", cache_only=False):
    """带磁盘缓存: 同 SMILES+depth+ionization 只算一次(3s), 后续读缓存(瞬时)。

    cache_only=True: **仅读缓存, 未缓存返回[]不跑docker** — 下载/批量导出路径必用,
    避免 detailed 报告逐化合物同步调 docker(~4s/个, 复杂的90s超时)阻塞单worker服务器。
    CFM 产物由单独的 precompute_cfm_fragments.py 预灌缓存(一次性, 可并行/后台)。
    返回 [{"mass":.., "smiles":..}, ...] 按质量降序。空/非法 SMILES 或失败返回 []。
    """
    if not smiles or not str(smiles).strip():
        return []
    smiles = str(smiles).strip()
    # 非法 SMILES 占位('N/A' 等)直接跳过 — 避免 fraggraph-gen segfault(rc=139)浪费时间
    if smiles.upper() in ("N/A", "NAN", "NONE", "NULL", "-"):
        return []
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    # 缓存键: 用 SMILES 的 md5(确定性, 跨进程稳定; Python 内置 hash() 因随机种子跨进程不一致会导致缓存永久 miss)
    import hashlib
    # v2_ 前缀: 极性派发修复后让旧(可能错极性)缓存自动失效, 无需手删 output/cfm_fragments/。
    cache_key = f"v2_d{depth}_{ionization}_{hashlib.md5(smiles.encode('utf-8')).hexdigest()[:16]}"
    cache_file = CACHE_DIR / f"{cache_key}.json"
    if cache_file.exists():
        try:
            cached = json.loads(cache_file.read_text())
            # 旧实现把 Docker 未启动/预测失败的 [] 也永久缓存，之后即使 Docker
            # 恢复仍永远返回空。空缓存只对 cache_only 读路径有效；交互请求应重试。
            if cached or cache_only:
                return cached
        except Exception:
            pass
    if cache_only:
        return []  # 下载/批量路径: 不跑 docker, 未缓存→空(RDKit兜底), 保持瞬时
    frags = run_fraggraph(smiles, depth=depth, ionization=ionization)
    if frags:
        try:
            cache_file.write_text(json.dumps(frags, ensure_ascii=False))
        except Exception:
            pass
    else:
        # 清掉历史负缓存，让 Docker 恢复后的下一次请求有机会重新预测。
        try:
            cache_file.unlink(missing_ok=True)
        except Exception:
            pass
    return frags


if __name__ == "__main__":
    # 自测: naringenin 碎裂产物结构
    naringenin = "C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O"
    import time
    t0 = time.time()
    frags = get_fragments(naringenin, depth=1, ionization="+")
    print(f"naringenin: {len(frags)} 个碎片产物 (耗时 {time.time()-t0:.1f}s)")
    print("前体 + top8 碎片结构:")
    for f in frags[:9]:
        print(f"  m/z {f['mass']:.3f}: {f['smiles'][:60]}")
