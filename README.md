# MSsystem

**English | [中文](#mssystem-1)**

Open-source LC–MS/MS annotation software for natural products. MSsystem matches every compatible MS/MS scan against experimental reference spectra, aggregates spectral support at the candidate-structure level across repeated scans and collision energies, and keeps the scan–reference–fragment links that support evidence review of candidate rankings and of the structural details reported.

This repository hosts the reviewed 2026-11-02 source snapshot released with the manuscript *MSsystem: natural product annotation using multiple tandem mass spectra*. The snapshot is byte-identical to Additional File 1 of the submission except for the public-release metadata changes listed in `RELEASE_NOTES.md`.

## Contents

| Path | Description |
| --- | --- |
| `matching/`, `pipeline/`, `acquisition/`, `calibration/`, `classifier/`, `governance/`, `kb_mining/`, `mssystem/`, `web/` | MSsystem source modules |
| `analyze_mzml_batch.py`, `build_from_mgf_msp.py` | Library construction and batch analysis entry points |
| `example/` | Executable example (Additional File 3): builds a local library from two MassBank references, searches two experimental queries, and verifies ranked candidates and matched peaks |
| `README_EN.md`, `README_ZH.md`, `PIPELINE.md`, `WEB_UI_GUIDE.md` | Snapshot documentation |
| `SOURCE_MANIFEST.json`, `SOURCE_PRIVACY_SCAN.json`, `CURRENT_LICENSE_SCOPE.json` | Snapshot integrity, privacy-scan and licence-scope records |
| `requirements.txt`, `requirements.lock` | Direct dependencies (lock records exact versions) |

## Installation

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.lock
```

The recorded runtime is CPython 3.13.5 on macOS arm64; the executable example was additionally validated in a fresh Python 3.12.3 environment on Linux x86-64 (`example/verification/`). Some packages are platform-sensitive; the lock lists direct dependencies rather than a complete transitive export. The full workbench also needs compatible shared data and a writable runtime directory (`MSSYSTEM_SHARED_DATA_ROOT`, `MSSYSTEM_RUNTIME_ROOT`, `MSSYSTEM_LIBRARY_WRITE_ENABLED`; see `README_EN.md`).

## Executable example

From the repository root:

```bash
uv venv --python python3.12 .venv
uv --no-cache pip install --python .venv/bin/python --only-binary :all: -r requirements.lock
bash example/run_demo.sh "$PWD/.venv/bin/python" "$PWD"
```

`example/run_demo.sh` takes the interpreter and the MSsystem source directory as arguments; inside this repository the source directory is the repository root. The standalone reviewer bundle described in `example/README.md` additionally contains the source archive itself. The four example spectra come from MassBank (CC BY-SA); attributions are listed in `example/README.md`.

## Licence

Project-owned source is released under the MIT License (root `LICENSE`). Vendored third-party bundles (Plotly.js, Cytoscape.js, SmilesDrawer, MassCube-derived fixtures) retain their original licences and notices — see `THIRD_PARTY_NOTICES.md` and `CURRENT_LICENSE_SCOPE.json`. Example spectra from MassBank retain CC BY-SA terms.

## Citation

See `CITATION.cff`. Please cite the released source snapshot together with the manuscript.

---

# MSsystem（中文）

面向天然产物的开源 LC–MS/MS 注释软件。MSsystem 将每张相容 MS/MS 扫描逐一与实验参考谱匹配，按候选结构在重复扫描与多个碰撞能量层面汇总谱学支持，并在候选记录中保留扫描—参考—碎片对应关系，供候选名次与可报告结构细节的复核。

本仓库随论文 *MSsystem: natural product annotation using multiple tandem mass spectra* 发布经审查的 2026-11-02 源码快照。除 `RELEASE_NOTES.md` 所列公开发布元数据调整外，快照与投稿附加文件 1 字节一致。可执行示例（附加文件 3）位于 `example/`，评价与统计复算包（附加文件 2）见 [Releases](../../releases)。

安装、运行示例与许可说明见上文英文部分；快照级中文文档见 `README_ZH.md` 与 `PIPELINE.md`。
