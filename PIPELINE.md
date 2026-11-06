# LC-MS/MS 质谱库匹配分析流程

## 概述

本项目实现了一个完整的 LC-MS/MS 数据分析和质谱库匹配流程，适用于中药复方等复杂天然产物样品的非靶向代谢组学分析。

**核心功能**：
- 从公开质谱数据库（GNPS、MassBank、MoNA 等）构建本地质谱库
- 空白溶剂峰扣除，去除背景污染
- 基于 Modified Cosine 相似度的 MS2 谱图匹配
- 自动化批量分析和结果报告生成

## 流程架构

```
┌─────────────────────────────────────────────────────────────┐
│ 1. 数据库构建 (build_from_mgf_msp.py)                       │
│    MGF/MSP 源文件 → 解析 → 去重 → SQLite 质谱库              │
├─────────────────────────────────────────────────────────────┤
│ 2. 空白扣除 (blank_subtraction.py)                           │
│    空白样品 MGF → BlankProfile → 样品峰扣除                   │
├─────────────────────────────────────────────────────────────┤
│ 3. 批量匹配 (analyze_mzml_batch.py)                          │
│    mzML/MGF 样品 + 质谱库 → Modified Cosine → 匹配结果       │
├─────────────────────────────────────────────────────────────┤
│ 4. 结果输出                                                  │
│    结果.csv (含 SELFIES/Confidence) + 报告.csv + 总体摘要.csv │
└─────────────────────────────────────────────────────────────┘
```

## 核心脚本

### 1. `build_from_mgf_msp.py` — 数据库构建

**功能**：从 MGF/MSP 源文件构建本地 SQLite 质谱库。

**输入**：`mgf_msp_database/` 目录下的 MGF 和 MSP 文件（支持 Pluskal Lab、GNPS、MoNA、MassBank 等多种格式）

**输出**：`output/local_library.db`（SQLite 数据库）

**关键特性**：
- 自动解析 MGF COMMENT 字段中的化合物名、SMILES、InChIKey、分子式等元数据
- 按 InChIKey 去重，优先保留具体分类标签（如 Flavonoid、Alkaloid）而非通用标签（如 MoNA-export-All_Spectra）
- 生成 Compounds、Spectra、Mapping 三表关联结构

**使用方法**：
```bash
python3 build_from_mgf_msp.py --mgf-dir /path/to/mgf_msp_database --output output
```

### 2. `blank_subtraction.py` — 空白溶剂峰扣除

**功能**：参考 MS-Dial 和 GNPS 策略，从空白样品聚合所有峰，对样品峰进行扣除。

**核心类**：
- `BlankProfile`：空白峰轮廓（m/z 数组 + 最大/平均/标准差强度）
- `load_blank_mgf()`：从空白 MGF 文件加载并聚合所有峰
- `subtract_blank()`：从样品峰中扣除空白贡献

**扣除模式**：
- `absolute`：直接强度减法（样品 - 空白）
- `relative`：按比例扣除
- `filter_only`：仅过滤低于阈值的峰

**使用方法**：
```python
from blank_subtraction import load_blank_mgf, subtract_blank
blank = load_blank_mgf("Blank-A.mgf")
clean_peaks = subtract_blank(sample_peaks, blank, mode="absolute", intensity_threshold_factor=0.1)
```

### 3. `analyze_mzml_batch.py` — 批量质谱匹配

**功能**：对 mzML/MGF 样品文件进行批量 MS2 谱图库匹配。

**输入**：
- 样品文件（.mzML 或 .mgf）
- 质谱库（local_library.db）
- 可选空白文件（自动按极性 A/P 匹配）

**输出**：每个样品生成：
- `<样品>_结果.csv`：匹配结果（含 SELFIES、Confidence 分级）
- `<样品>_报告.csv`：详细报告（含 Mass_PPM、各分项得分）
- `总体摘要.csv`：批量汇总统计

**关键参数**：
- 余弦相似度 ≥ 0.7（Modified Cosine）
- 综合得分 ≥ 0.5（cosine×0.4 + mass×0.2 + fragment×0.25 + intensity×0.15）
- m/z 容差：0.5 Da 或 10 ppm

**使用方法**：
```bash
# 无空白扣除
python3 analyze_mzml_batch.py --samples /path/to/samples --output results

# 含空白扣除（自动按极性匹配）
python3 analyze_mzml_batch.py \
  --samples /path/to/samples \
  --blank-a Blank-A.mgf \
  --blank-p Blank-P.mgf \
  --blank-mode absolute \
  --output results

# 仅处理特定文件
python3 analyze_mzml_batch.py --samples . --include "WH-1-A" "WH-1-P"
```

## 数据库来源

| 来源 | 谱图数 | <500 Da | 说明 |
|---:|---:|---:|---:|
| MoNA (MassBank of North America) | 794K | 2% | 公共质谱库，大分子为主 |
| MassBank Europe (NIST) | 139K | 84% | 精选小分子库 |
| MassBank Europe (RIKEN) | 167K | 84% | 天然产物/植物代谢物 |
| GNPS-LIBRARY | 21K | 63% | 社区精选天然产物库 |
| GNPS ALL_GNPS | 2.9M | — | GNPS 全库（下载中） |
| VF-NPL | 10K | 65% | 天然产物库 |
| HCD Natural Product | 1.5K | 59% | 天然产物 HCD 谱库 |
| MGF_classify_data | ~200K | 92% | 按分类整理的库（Flavonoid, Alkaloid 等） |
| **总计** | **~1.2M** | — | 去重后约 877K 条 |

## 匹配算法

### Modified Cosine 相似度
- 使用 matchms 库的 ModifiedCosine 算法
- 容忍度：m/z ±0.1 Da
- 强度归一化至 0-1

### 综合评分
```
Final_Score = Cosine_Score × 0.40 + Mass_Score × 0.20 + Fragment_Score × 0.25 + Intensity_Score × 0.15
```

### 置信度分级
- **High**：cos ≥ 0.9 且 final ≥ 0.5
- **Medium**：cos ≥ 0.8 且 final ≥ 0.4
- **Low**：cos ≥ 0.7
- **Unknown**：cos < 0.7

## 输出格式

### 结果.csv 字段
| 字段 | 说明 |
|---|---|
| Spectrum | 谱图编号（scan ID） |
| Precursor_mz | 前体离子 m/z |
| Match_Name | 匹配化合物名称 |
| Match_SMILES | SMILES 结构式 |
| Match_SELFIES | SELFIES 结构式 |
| Match_InChIKey | InChIKey 标识符 |
| Match_Formula | 分子式 |
| Final_Score | 综合得分 |
| Cosine_Score | 余弦相似度 |
| Confidence | 置信度（High/Medium/Low） |
| Class | 化合物分类 |

### 报告.csv 字段
含 Peaks_Count、Mass_PPM、Fragment_Score、Intensity_Score 等分项得分。

## 依赖环境

```
python >= 3.9
matchms >= 0.18.0
numpy
rdkit
selfies
sqlite3
```

安装：
```bash
conda create -n msnovelist python=3.9
conda activate msnovelist
pip install matchms numpy rdkit selfies
```
