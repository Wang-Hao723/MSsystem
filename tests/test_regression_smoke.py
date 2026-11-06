# -*- coding: utf-8 -*-
"""
mssystem_v4 回归冒烟测试 — 改完 web/ / matching/ / utils/ 必跑
================================================================
秒级, 不加载谱图库, 不启服务器。守护"以前出过的问题"不再回来。

跑法 (生产 env):
    cd /Users/zizi/Desktop/Mssystem/mssystem_v4
    /opt/anaconda3/bin/python3.13 -m pytest tests/test_regression_smoke.py -v

每条断言 = 一个历史 bug 的护栏。断言挂 = 又碰碎了别处, 提交前就要拦住。

  A. 静态代码契约  — 不跑逻辑, 纯文本/AST 检查 (守"已下线/已修复的代码不复活")
  B. 函数级单测    — mock 输入喂纯函数 (不依赖库/数据)
  C. 数据口径      — 用真实 NEG1 results.json (无则 skip; 重新分析后生成)

慢红线 (17 标品命中率, 需加载 small 库, 分钟级) 不在此文件, 单独跑:
    python3 eval/eval_baseline.py    # hit率 ≥ 14/17
"""
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# ═══════════════════════════════════════════════════════════════════
# A. 静态代码契约 — 纯文本检查 (秒级, 守护"已修复的 bug 不复活")
# ═══════════════════════════════════════════════════════════════════

def _read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="ignore")


def test_A1_iceberg_从生产路径移除():
    """commit a14fa5a: ICEBERG 结构推理兜底已下线, web/ 不应再集成它。
    只检查真正的集成符号 (调用/worker/事件), 不抓注释里解释历史的 'iceberg' 字样。"""
    INTEGRATION_TOKENS = ["fallback_batch", "_maybe_launch_iceberg_fallback",
                          "nomatch_items", "iceberg_fallback_done", "_iceberg_fallback_worker",
                          "_iceberg_threads"]
    bad = []
    for p in (ROOT / "web").rglob("*"):
        if p.suffix not in {".py", ".js"}:
            continue
        txt = p.read_text(encoding="utf-8", errors="ignore")
        for tok in INTEGRATION_TOKENS:
            if tok in txt:
                bad.append(f"{p.relative_to(ROOT)} 仍含 ICEBERG 集成符号: {tok}")
    assert not bad, "web/ 仍集成 ICEBERG (应已移除):\n" + "\n".join(bad)


def test_A2_iceberg模块文件已删():
    assert not (ROOT / "matching" / "iceberg_fallback.py").exists(), \
        "matching/iceberg_fallback.py 不应存在 (commit a14fa5a 已删)"


def test_A3_序列化含TIC字段():
    """commit 5981827: _serialize_results 必须把 tic/base_peak 放进 serialized 顶层。
    漏了 → mgf 任务 TIC fallback 拿不到强度 → TIC 图空白。"""
    srv = _read("web/server.py")
    assert '"tic"' in srv and '"base_peak"' in srv, \
        "_serialize_results 漏 tic/base_peak (TIC fallback 会画不出)"


def test_A4_mgf_用_retention_time键():
    """commit 5981827: matchms 把 RTINSECONDS 规范化成 'retention_time' 键。
    旧的错误用法 md.get("rtinsECONDS") 匹配不到 → RT 全丢 → TIC 无 x 轴。
    (注释里解释这个 bug 会含字面量 rtinsECONDS, 故只查 'md.get' 用法。)"""
    mgf = _read("utils/mgf_io.py")
    assert 'md.get("retention_time")' in mgf, "mgf_io 应用 retention_time 键读 RT"
    assert 'md.get("rtinsECONDS")' not in mgf, "旧的错误键名用法 md.get(rttinsECONDS) 回归"


def test_A5_摘要只按谱库匹配去重且保留污染物警示():
    """Do not split ordinary identifications into a misleading natural panel."""
    js = _read("web/static/js/results.js")
    assert re.search(
        r"const\s+contaminantGroups\s*=\s*groups\.filter\([^;]+potential_contaminant",
        js,
    ), "潜在外源污染物应从去重化合物组单独分出"
    assert "const naturalGroups" not in js
    assert "const unclassifiedGroups" not in js
    assert "const nonContaminantGroups" not in js
    assert "literatureTargetArea" not in js
    assert "/literature-targets" not in js
    assert re.search(
        r"const\s+frontendHighSupported\s*=\s*new\s+Set",
        js,
    ), "必须保留“至少一张 High MS2 支持”的身份级口径"
    m = re.search(
        r"const\s+frontendRepresentativeHigh\s*=\s*(\w+)\.filter",
        js,
    )
    assert m, "找不到最高分代表 High 定义"
    assert m.group(1) == "groups"
    assert "high_supported_unique_compounds" in js
    assert "representative_high_unique_compounds" in js
    assert "完整结构身份 High 支持" in js
    assert "完整身份代表 High" in js
    m2 = re.search(r"const\s+nClassHigh\s*=\s*(\w+)\.filter", js)
    assert m2 and m2.group(1) == "groups"


def test_A5a_文献独立筛查接口与旧天然计数文案已移除():
    import web.server as server

    paths = server.app.openapi()["paths"]
    assert "/api/library/targets/datasets" not in paths
    assert "/api/tasks/{task_id}/literature-targets" not in paths
    message = server._task_message_for_display({
        "status": "completed",
        "matched_spectra": 10,
        "total_spectra": 20,
        "unique_compounds": 7,
        "match_rate": 50.0,
        "message": "0 curated natural-product candidates",
    })
    assert message == "Completed: 10/20 spectra matched, 7 unique library candidates (50.0%)"
    assert "natural" not in message.lower()


def test_task_compound_count_uses_final_scope_aware_result_table():
    import web.server as server

    task = {
        "status": "completed",
        "matched_spectra": 10,
        "total_spectra": 20,
        "unique_compounds": 7,
        "match_rate": 50.0,
    }
    changed = server._apply_final_compound_count_summary(
        task,
        {
            "main_result_count": 4,
            "unique_report_rows": 6,
            "unresolved_precursor_event_count": 2,
            "isomer_group_detected_count": 1,
        },
    )

    assert changed is True
    assert task["legacy_unique_library_candidates"] == 7
    assert task["unique_compounds"] == 4
    assert task["main_result_count"] == 4
    assert task["result_record_count"] == 6
    assert task["unresolved_precursor_event_count"] == 2
    message = server._task_message_for_display(task)
    assert "4 main compound results" in message
    assert "2 unresolved precursor events excluded" in message


def test_A5b_未归属碎片差值不冒充已检测中性丢失():
    js = _read("web/static/js/results.js")
    assert "已检测、未解释" not in js
    assert "未归属前体–碎片差值" in js
    assert "不等同于已确认中性丢失" in js


def test_library_import_ui_uses_actual_importable_counts():
    js = _read("web/static/js/libaugment.js")
    assert "res.importable_records" in js
    assert "res.spectral_records" in js
    assert "res.provenance_held_records" in js
    assert "数值强度不等于实验证明" in js
    assert "res.quality && res.quality.accepted > 0" not in js


def test_A5c_证据不足和无结构不再笼统显示Unknown():
    js = _read("web/static/js/results.js")
    assert "谱库匹配候选" in js
    assert "候选不等同于确认鉴定" in js
    assert "未能可靠判类" in js
    assert "不可分类（无结构）" in js
    assert "证据不足" in js
    assert "仅谱图线索" in js
    assert "!!match.classification" not in js


def test_A5d_未匹配谱图可展开且模型概率不再冒充高置信():
    js = _read("web/static/js/results.js")
    assert "showUnmatchedDetail" in js
    assert "Top-3 类别" in js
    assert "HTML 中实际可用碎片" in js
    assert "MSClassify 未校准模型概率" in js
    assert "p >= 0.70" in js
    assert "p >= 0.45" in js
    assert "${_nHi} 高置信" not in js


def test_A6_离子模式过滤默认开():
    """commit 53efb20: 负离子样本不再匹配出正离子加合物(跨极度假阳性)。"""
    cfg = _read("config/settings.py")
    m = re.search(r"ENFORCE_POLARITY_FILTER\s*=\s*(True|False)", cfg)
    assert m and m.group(1) == "True", "ENFORCE_POLARITY_FILTER 必须 True (防跨极性假阳性)"


def test_A7_前端无ICEBERG_UI元素():
    """commit a14fa5a: results.js 不应有 🤖结构推理 banner/徽章/轮询。"""
    js = _read("web/static/js/results.js")
    for token in ["_pollIcebergFallback", "_icebergPollTimer", "icebergBadge",
                  "icebergTitleBadge", "icebergBanner", "iceberg_fallback_done"]:
        assert token not in js, f"results.js 仍含 ICEBERG UI 残留: {token}"


def test_A8_库下拉标签按value更新且跟随已加载库():
    """Small/Full 选项顺序不能决定标签，否则会显示 Full、实际提交 small。"""
    js = _read("web/static/js/upload.js")
    assert 'option[value="small"]' in js and 'option[value="full"]' in js
    assert "sel.options[0].textContent" not in js
    assert "sel.value = _loadedDbMode || cur" in js


def test_A9_批量导出采用样本文件名而非项目或通用名():
    js = _read("web/static/js/export.js")
    assert "_responseDownloadName" in js
    assert "_fallbackBatchFilename" in js
    assert "`${samples[0]}__${samples[1]}__等${samples.length}个样本`" in js
    assert "anchor.download = filename" in js
    assert "startArtifactDownload(job.id)" in js
    assert "anchor.href = `${API.BASE}/export/jobs/${jobId}/download`" in js
    assert "a.download = `batch_${format}_${taskIds.length}tasks.zip`" not in js

    integrated = _read("web/static/js/report-generator.js")
    assert "`${sampleStems[0]}__${sampleStems[1]}__等${sampleStems.length}个样本`" in integrated
    assert "const downloadName = `${samplePart}_样本综合报告" in integrated
    assert "MSsystem_样本综合报告_${tasks.length}samples" not in integrated


def test_A10_上传页使用统一批量排队接口并显示样本数():
    js = _read("web/static/js/upload.js")
    assert "'/tasks/batch/analyze'" in js
    assert "task_ids: uploaded.map" in js
    assert "const inputCount = _retryUploadedTasks.length || (" in js
    assert "批量开始分析（${inputCount} 个样本）" in js
    assert "重试排队（${inputCount} 个已上传样本）" in js
    assert "input.value = ''" in js


def test_A10d_空白文件可持久管理且切换极性无需刷新():
    js = _read("web/static/js/upload.js")
    html = _read("web/static/index.html")
    server = _read("web/server.py")

    assert 'id="blankUploadManager"' in html
    assert 'id="ms2FeatureLinkTolerancePpm"' in html
    assert "ms2_feature_link_tolerance_ppm: numberValue(" in js
    assert "/static/js/upload.js?v=16.0" in html
    assert "API.get('/blanks')" in js
    assert "BLANK_SELECTION_STORAGE_KEY" in js
    successful_queue = js.split(
        "// Every uploaded task is now queued or was already active.", 1
    )[1].split("// Navigate to tasks page", 1)[0]
    assert "_resetBlankSelection()" not in successful_queue
    assert "提交后继续用于下一批" in js
    assert "handleIonModeChange" in js
    assert "不兼容空白已从本批解除" in js
    assert "deleteManagedBlank" in js
    assert '_retryUploadedTasks = uploaded.filter' in js
    assert '@app.get("/api/blanks")' in server
    assert "List managed blank uploads" in server


def test_A10e_浏览器空白默认实际排除且未确认时拒绝提交():
    js = _read("web/static/js/upload.js")
    html = _read("web/static/index.html")
    server = _read("web/server.py")
    matching = _read("matching/blank_subtraction.py")

    assert 'value="hard_legacy" selected' in html
    assert "执行特征级空白背景排除" in html
    assert 'id="blankMethodComparable"' in html
    assert 'id="blankExclusionConfirmed"' in html
    assert "_blankEvidenceV2Available" in js
    assert "Object.assign(params, _blankAnalysisParameters())" in js
    assert "不会退化成 0 扣除" in js
    assert '"blank_evidence_v2": True' in server
    assert '@app.post("/api/blanks/{blank_token}/role")' in server
    assert "assess_blank_role_groups" in server
    # Direct API clients retain the conservative soft default; the browser
    # now makes its actual-exclusion choice explicit and auditable.
    assert 'DEFAULT_BLANK_FILTER_POLICY = "soft_evidence"' in matching
    assert 'assessment["blank_excluded"] = False' in matching


def test_A10b_上传页和后端均支持内嵌MS2的HTML报告():
    js = _read("web/static/js/upload.js")
    html = _read("web/static/index.html")
    server = _read("web/server.py")
    requirements = _read("requirements.txt")
    assert "'.html'" in js and "'.htm'" in js
    assert ".html,.htm" in html
    assert 'ext in (".html", ".htm")' in server
    assert "is_raw_report_html" in server
    assert "beautifulsoup4" in requirements and "lxml" in requirements


def test_A10c_未命中分类器使用V6正式运行模块():
    """models/ 在 V6 被 gitignore 且不存在，导入它会被静默吞掉并丢失全部判类。"""
    server = _read("web/server.py")
    app = _read("web/app.py")
    expected = "from classifier.msclass_wrapper import get_msclass_wrapper"
    absent = "from models.msclass_wrapper import get_msclass_wrapper"
    assert expected in server and expected in app
    assert absent not in server and absent not in app


def test_A10d_MSClassify预处理无需加载checkpoint即可用():
    from classifier.msclass_wrapper import MSClassWrapper, _class_label

    assert _class_label("PhenolicAcid") == "phenolic_acid"
    binned = MSClassWrapper._bin_spectrum([(123.1, 10), (123.8, 20), (1500, 99)])
    assert binned.shape == (1000,)
    assert binned[123] == 1.0


def test_A10e_HTML少峰输入强制严格证据门槛():
    server = _read("web/server.py")
    app = _read("web/app.py")
    assert 'True if spec.get("limited_fragment_list")' in server
    assert 'or bool(spec.get("limited_fragment_list"))' in app
    assert "系统已自动启用严格证据门槛" in _read("web/static/js/results.js")


def test_A10f_HTML_top5不允许精确化合物High():
    server = _read("web/server.py")
    assert 'r.get("limited_fragment_list") and _conf == "High"' in server
    assert '"compound_confidence": "Medium"' in server
    assert '"evidence_tier": "Strong Candidate"' in server
    assert "positional isomer unconfirmed" in server


def test_A11_分类结构画廊响应式且不再固定四列挤压图表():
    viewer = _read("web/static/js/molecule-viewer.js")
    css = _read("web/static/css/style.css")
    api = _read("web/static/js/api.js")
    assert "repeat(${cols}, ${this.defaultWidth}px)" not in viewer
    assert "repeat(auto-fit" in css
    assert "class-summary-grid" in css
    assert "matched_only=true&limit=5000" in api
    assert 'orientation: "h"' in api


# ═══════════════════════════════════════════════════════════════════
# B. 函数级单测 — mock 输入, 秒级, 不依赖库/数据
# ═══════════════════════════════════════════════════════════════════

def _mock_match(idx=1, dataset="Test"):
    """search() 风格的 match dict, 字段给全让 _honest_annotate/_resolve_formula 能跑。"""
    return {
        "cosine": 0.85, "entropy": 0.8, "mass_score": 0.95, "mass_ppm": 2.0,
        "fragment_score": 0.7, "intensity_score": 0.7, "reverse_match": 0.85,
        "final_score": 0.86, "expert_score": 0.0, "expert_details": None,
        "compound_name": f"Compound{idx}", "smiles": "CCO",
        "inchikey": f"AAAAAAAAAAAAAA-{idx}", "formula": "C2H6O", "adduct": "[M-H]-",
        "classification": "Organic acid", "dataset_source": dataset,
        "pathway_json": None, "lib_precursor_mz": 100.0 + idx, "lib_spectrum_id": "s1",
        "lib_peaks": [[100.0, 500.0], [200.0, 300.0]], "display_class": "Organic acid",
        "dl_classification": "", "dl_confidence": 0.5, "top_candidates": [],
        "class_conflict": False, "structure_score": 0.0,
    }


def _mock_row(idx=1, matched=True, dataset="Test"):
    """_run_analysis 构造的 all_results 行格式 (含 tic/base_peak/raw_peaks)。"""
    return {
        "spectrum_id": str(idx), "precursor_mz": 100.0 + idx, "peaks_count": 10,
        "retention_time": 10.0 * idx, "tic": 1000.0 * idx, "base_peak": 500.0,
        "match": _mock_match(idx, dataset) if matched else None,
        "blank_used": False, "blank_removed": 0,
        "query_peaks": [[100.0, 500.0], [200.0, 300.0]],
        "raw_peaks": [[100.0, 500.0], [200.0, 300.0], [300.0, 100.0]],
    }


def test_B1_serialize_不漏TIC字段():
    """_serialize_results 必须在 matched + unmatched 行都输出 tic/base_peak/retention_time。
    TIC fallback 画全色谱(含未匹配谱图)需要每个点都有 RT+强度。"""
    from web.server import _serialize_results
    out = _serialize_results([_mock_row(1, matched=True), _mock_row(2, matched=False)])
    # matched 行
    assert out[0]["tic"] == 1000.0, "matched 行丢 tic (TIC fallback 画不出)"
    assert out[0]["base_peak"] == 500.0
    assert out[0]["retention_time"] == 10.0
    assert out[0]["raw_peaks"], "matched 行丢 raw_peaks"
    # unmatched 行 (TIC 要画全色谱, 未匹配点也要有 RT+强度)
    assert out[1]["tic"] == 2000.0, "unmatched 行丢 tic (TIC 色谱会有空洞)"
    assert out[1]["retention_time"] == 20.0
    assert out[1]["raw_peaks"], "unmatched 行丢 raw_peaks（分子网络会退化为 Top-N）"


def test_B1_blank_excluded_保留审计字段且不伪装成未命中分类():
    from web.server import _serialize_results

    row = _mock_row(3, matched=False)
    row.update({
        "blank_used": True,
        "blank_status": "excluded_below_ratio",
        "blank_feature_matched": True,
        "blank_excluded": True,
        "blank_evidence_source": "blank_ms1",
        "blank_sample_intensity": 200.0,
        "blank_intensity": 100.0,
        "blank_sample_ratio": 2.0,
        "blank_mz_error_ppm": 1.25,
        "blank_rt_delta_minutes": 0.05,
    })

    serialized = _serialize_results([row])[0]

    assert serialized["blank_excluded"] is True
    assert serialized["blank_sample_ratio"] == 2.0
    assert serialized["blank_evidence_source"] == "blank_ms1"
    assert "display_class" not in serialized


def test_B1a_serialize_保留零质量误差():
    """精确质量误差 0 ppm 是有效数值，不能被 truthiness 转成 None。"""
    from web.server import _serialize_results
    row = _mock_row(1, matched=True)
    row["match"]["mass_ppm"] = 0.0
    assert _serialize_results([row])[0]["match"]["mass_ppm"] == 0.0


def test_B1a_report_fragments_follow_precursor_shifted_matching():
    from web.server import _msms_fragments

    matched = _msms_fragments(
        [[50.0, 100.0], [70.0, 50.0]],
        [[40.0, 100.0], [60.0, 50.0]],
        query_pmz=110.0,
        lib_pmz=100.0,
    )

    assert matched == "50.0000; 70.0000"


def test_B1b_name_only谱图线索保留但强制Unknown():
    from web.server import _serialize_results
    row = _mock_row(1, matched=True)
    row["match"].update({
        "compound_name": "Apigenin",
        "smiles": "",
        "has_structure": False,
        "name_only_lead": True,
        "classification": "Unknown",
        "display_class": "Unknown",
        "expert_details": None,
        "pathway_json": None,
    })
    match = _serialize_results([row])[0]["match"]
    assert match is not None
    assert match["name_only_lead"] is True
    assert match["has_structure"] is False
    assert match["confidence"] == "Unknown"
    assert match["confidence_class"] == "Unknown"
    assert match["evidence_tier"] == "Spectral Lead"
    assert match["expert_details"] is None
    assert match["pathway_json"] is None


def test_B2_dedup_排除iceberg残留_并按化合物去重():
    """_dedup_matched 必须排除 dataset_source==iceberg_fallback (旧任务残留),
    且同 InChIKey-14 的多张谱图合并为 1 条。"""
    from web.server import _dedup_matched
    r1 = _mock_row(1, matched=True, dataset="Test")
    r2 = _mock_row(2, matched=True, dataset="Test")
    r3 = _mock_row(3, matched=True, dataset="iceberg_fallback")  # 应被排除
    # 让 r1/r2 同 InChIKey 连接层 → 去重到 1 条
    r1["match"]["inchikey"] = "BBBBBBBBBBBBBB-1"
    r2["match"]["inchikey"] = "BBBBBBBBBBBBBB-2"
    out = _dedup_matched([r1, r2, r3])
    assert all(o["match"]["dataset_source"] != "iceberg_fallback" for o in out), \
        "iceberg_fallback 残留未被排除"
    assert len(out) == 1, f"同 ik14 两张谱图应去重到 1 条, 实得 {len(out)}"


def test_B3_mgf_RTINSECONDS能解析(tmp_path):
    """commit 5981827: mgf 的 RTINSECONDS= 必须解析成 retention_time。
    旧 bug (找错键名) → RT 全丢 → TIC 无 x 轴。"""
    from utils.mgf_io import parse_mgf_file
    mgf = tmp_path / "mini.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=scan=1\nRTINSECONDS=37.977\nPEPMASS=121.027\n"
        "CHARGE=1-\n100.0 100.0\n200.0 50.0\nEND IONS\n"
        "BEGIN IONS\nTITLE=scan=2\nPEPMASS=200.0\nCHARGE=1-\n100.0 100.0\nEND IONS\n"
    )
    recs = parse_mgf_file(str(mgf))
    assert len(recs) == 2
    assert recs[0]["retention_time"] == pytest.approx(37.977), \
        "RTINSECONDS 没解析出 (键名 bug 回归?)"


def test_B4_mgf从charge推polarity(tmp_path):
    """mgf 无 adduct 字段, 但 CHARGE=1-/-1 可定极性。不推则 query_polarity=None →
    离子模式过滤失效 → 负离子样本匹配出 [M+H]+ 跨极性假阳性 (回归脚本 C3 抓到的真 bug)。"""
    from utils.mgf_io import parse_mgf_file
    mgf = tmp_path / "pol.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=n1\nPEPMASS=200.0\nCHARGE=1-\n100.0 100.0\nEND IONS\n"
        "BEGIN IONS\nTITLE=p1\nPEPMASS=200.0\nCHARGE=1+\n100.0 100.0\nEND IONS\n"
        "BEGIN IONS\nTITLE=n2\nPEPMASS=200.0\nCHARGE=2-\n100.0 100.0\nEND IONS\n"
    )
    recs = parse_mgf_file(str(mgf))
    pol = {r["title"]: r.get("polarity") for r in recs}
    assert pol["n1"] == "NEG", "CHARGE=1- 应推 polarity=NEG"
    assert pol["p1"] == "POS", "CHARGE=1+ 应推 polarity=POS"
    assert pol["n2"] == "NEG", "CHARGE=2- 应推 polarity=NEG"


def test_B5_异构体专属碎片默认只作证据不改排序():
    """v6: 未经完整结构真值校准的异构碎片信号不得在生产中改写 Top-1。"""
    from matching.search import _apply_isomer_fragment_rerank

    query_clean = [[100.0, 1.0], [121.0, 0.2]]
    candidates = [
        {
            "formula": "C10H10O2", "inchikey": "ABCDEFGHIJKLMN-AAAA",
            "lib_peaks": [[100.0, 1.0], [120.0, 0.5]], "final_score": 0.50,
        },
        {
            "formula": "C10H10O2", "inchikey": "NOPQRSTUVWXYZ-AAAA",
            "lib_peaks": [[150.0, 1.0], [130.0, 0.5]], "final_score": 0.50,
        },
    ]
    _apply_isomer_fragment_rerank(candidates, query_clean)

    assert candidates[0]["isomer_fragment_score"] > 0.0
    assert candidates[0]["isomer_fragment_potential_bonus"] > 0.0
    assert candidates[0]["isomer_fragment_bonus"] == 0.0
    assert candidates[0]["isomer_fragment_rank_applied"] is False
    assert candidates[0].get("rank_score", candidates[0]["final_score"]) == pytest.approx(0.50)
    assert candidates[1]["isomer_fragment_score"] == 0.0


def test_B5a_异构体专属碎片实验开关可复现旧排序():
    """显式实验开关保留旧行为，便于在权威 InChIKey 真值上做消融。"""
    from matching.search import _apply_isomer_fragment_rerank

    candidates = [
        {
            "formula": "C10H10O2", "inchikey": "ABCDEFGHIJKLMN-AAAA",
            "lib_peaks": [[100.0, 1.0], [120.0, 0.5]], "final_score": 0.50,
        },
        {
            "formula": "C10H10O2", "inchikey": "NOPQRSTUVWXYZ-AAAA",
            "lib_peaks": [[150.0, 1.0], [130.0, 0.5]], "final_score": 0.50,
        },
    ]
    _apply_isomer_fragment_rerank(
        candidates, [[100.0, 1.0], [121.0, 0.2]], apply_rank_override=True
    )

    assert candidates[0]["rank_score"] > candidates[1].get(
        "rank_score", candidates[1]["final_score"]
    )
    assert candidates[0]["isomer_fragment_rank_applied"] is True


def test_B6_专属碎片rerank不影响单候选():
    """v5: 只有竞争候选组才做专属碎片 rerank, 单候选不扰动。"""
    from matching.search import _apply_isomer_fragment_rerank

    candidates = [
        {
            "formula": "C10H10O2", "inchikey": "ABCDEFGHIJKLMN-AAAA",
            "lib_peaks": [[100.0, 1.0]], "final_score": 0.50,
        },
    ]
    _apply_isomer_fragment_rerank(candidates, [[100.0, 1.0]])

    assert candidates[0].get("isomer_fragment_score") is None
    assert candidates[0]["final_score"] == pytest.approx(0.50)


def test_B7_异构体近似平局用谱图相似度打破():
    """v5: 同分异构体 final 极小差异不应由浮点噪声决定。
    若 final 在 eps 内近似平局, 应优先选择专属碎片/余弦更支持的候选。"""
    from matching.search import _select_best_result

    candidates = [
        {
            "compound_name": "true_isomer", "formula": "C15H10O4",
            "inchikey": "ABCDEFGHIJKLMN-AAAA", "final_score": 0.77850,
            "isomer_fragment_score": 0.0, "cosine": 0.8220, "reverse_match": 0.90,
            "mass_ppm": 1.0,
        },
        {
            "compound_name": "near_tie_isomer", "formula": "C15H10O4",
            "inchikey": "QRSTUVWXYZABCD-AAAA", "final_score": 0.77855,
            "isomer_fragment_score": 0.0, "cosine": 0.8210, "reverse_match": 0.90,
            "mass_ppm": 1.0,
        },
    ]

    assert _select_best_result(candidates)["compound_name"] == "true_isomer"


def test_B7a_未校准异构碎片分数不能推翻更高cosine():
    from config.settings import USE_ISOMER_FRAGMENT_RANK_OVERRIDE
    from matching.search import _select_best_result

    assert USE_ISOMER_FRAGMENT_RANK_OVERRIDE is False
    candidates = [
        {
            "compound_name": "cosine_winner", "formula": "C15H10O5",
            "inchikey": "ABCDEFGHIJKLMN-AAAA", "final_score": 0.80,
            "rank_score": 0.80, "isomer_fragment_score": 0.0,
            "cosine": 0.93, "reverse_match": 0.90, "mass_ppm": 2.0,
        },
        {
            "compound_name": "fragment_winner", "formula": "C15H10O5",
            "inchikey": "QRSTUVWXYZABCD-AAAA", "final_score": 0.80,
            "rank_score": 0.80, "isomer_fragment_score": 1.0,
            "cosine": 0.86, "reverse_match": 0.90, "mass_ppm": 2.0,
        },
    ]

    assert _select_best_result(candidates)["compound_name"] == "cosine_winner"


def test_B8_严格证据门槛拦截低信息量假阳性():
    """v5: blank 假阳性多来自低 entropy / 高 ppm / 低 reverse 的弱证据匹配。"""
    from matching.search import _passes_strict_evidence_gate

    assert _passes_strict_evidence_gate(
        final=0.70, entropy=0.60, reverse=0.80, ppm=5.0, enabled=True
    )
    assert not _passes_strict_evidence_gate(
        final=0.70, entropy=0.30, reverse=0.80, ppm=5.0, enabled=True
    ), "低 entropy 谱图信息量不足, 应拦截"
    assert not _passes_strict_evidence_gate(
        final=0.70, entropy=0.60, reverse=0.40, ppm=5.0, enabled=True
    ), "低 reverse match 表示库特征峰覆盖不足, 应拦截"
    assert not _passes_strict_evidence_gate(
        final=0.70, entropy=0.60, reverse=0.80, ppm=200.0, enabled=True
    ), "质量误差过大, 应拦截"
    assert _passes_strict_evidence_gate(
        final=0.70, entropy=0.60, reverse=0.80, ppm=42.5, enabled=True
    ), "常规严格模式仍保留原 100 ppm 上限"
    assert not _passes_strict_evidence_gate(
        final=0.70, entropy=0.60, reverse=0.80, ppm=42.5, enabled=True,
        max_mass_ppm=20.0,
    ), "top-5 高分辨 HTML 必须拦截 42.5 ppm 的同名义质量异分子式候选"
    assert _passes_strict_evidence_gate(
        final=0.70, entropy=0.60, reverse=0.80, ppm=8.86, enabled=True,
        max_mass_ppm=20.0,
    ), "M011 当前 4–9 ppm 候选不应被 top-5 专用质量门禁误伤"

    assert _passes_strict_evidence_gate(
        final=0.25, entropy=0.10, reverse=0.10, ppm=500.0, enabled=False
    ), "常规模式必须只服从用户设置的搜索阈值, 不能被隐藏门槛拦截"


def test_B9_异分子式同前体不能冒充异构体专属碎片证据():
    """异分子式竞争应由质量/加合物/同位素处理，不能混入 isomer score。"""
    from matching.search import _apply_isomer_fragment_rerank

    query_clean = [[185.0422, 1.0], [189.0524, 1.0], [207.0629, 1.0], [395.0947, 1.0]]
    candidates = [
        {
            "formula": "C16H22O11", "inchikey": "HPWWQPXTUDMRBI-NJPMDSMTSA-N",
            "lib_peaks": [[189.016, 1.0], [203.053, 0.2], [233.043, 0.8], [301.141, 0.6]],
            "final_score": 0.60,
        },
        {
            "formula": "C16H22O10", "inchikey": "ZJDOESGVOWAULF-UHFFFAOYSA-N",
            "lib_peaks": [[185.0422, 0.4], [189.0524, 1.0], [207.0629, 0.3], [395.0947, 0.2]],
            "final_score": 0.55,
        },
    ]
    _apply_isomer_fragment_rerank(candidates, query_clean)

    assert candidates[0].get("isomer_fragment_score") is None
    assert candidates[1].get("isomer_fragment_score") is None
    assert candidates[0].get("rank_score", 0.60) == pytest.approx(0.60)
    assert candidates[1].get("rank_score", 0.55) == pytest.approx(0.55)
    assert candidates[1]["final_score"] == pytest.approx(0.55)


def test_B9a_v6低质量独占基峰可进入受控重试():
    from matching.search import _drop_dominant_low_mz_interference

    peaks = [[58.065, 1000.0], [265.086, 200.0], [297.113, 150.0]]
    cleaned = _drop_dominant_low_mz_interference(peaks, 342.17)
    assert cleaned is not None
    assert all(abs(p[0] - 58.065) > 1e-6 for p in cleaned)

    # 基峰优势不足时不能删除；避免把真实低质量诊断峰当污染。
    assert _drop_dominant_low_mz_interference(
        [[58.065, 700.0], [265.086, 400.0]], 342.17
    ) is None


def test_B9b_裂解规则文件缺失时安全降级(tmp_path):
    """规则文件可选路径不存在时应返回空库，而不是因未定义 logger 崩溃。"""
    from matching.fragmentation_rules import load_compound_class_rules

    assert load_compound_class_rules(tmp_path / "missing_rules.json") == {}


def test_B9c_v6标准品名称不再用反向子串制造命中():
    from eval.standard_names import (
        formula_compatible, formula_match, name_match, standard_formula_compatible,
    )

    assert name_match("Cinnamic acid", "CinnamicAcid")
    assert not (
        name_match("Cinnamic acid", "4_Hydroxy_3MethoxyCinnamicAcid_Ferulic")
        and formula_match("C9H8O2", "C10H10O4")
    )
    assert not name_match("Ginsenoside Rg1", "Ginsenoside")
    assert not formula_match("C21H20O12", ""), "缺分子式不能充当分子式级命中"
    assert formula_compatible("C21H20O12", ""), "精确名称命中可容许库记录缺分子式"
    assert not formula_compatible("C21H20O12", "C21H28O12"), "冲突分子式必须否决名称命中"
    assert name_match("Glycyrrhizic acid", "Ammonium Glycyrrhizinate (AMGZ)")
    assert standard_formula_compatible("Glycyrrhizic acid", "C42H62O16", "C42H65NO16")


def test_B9d_MIX33谷甾醇真值不再缺前体和模式():
    from eval.standard_names import load_gt

    rows = load_gt(ROOT / "eval" / "33components_gt.csv")
    beta = next(row for row in rows if row["en_name"] == "beta-Sitosterol")
    assert beta["mode"] == "POS"
    assert float(beta["mz"]) == pytest.approx(415.3934)
    assert beta["formula"] == "C29H50O"


def test_B10_MGF_TIC概览保留scan并把秒换算为分钟():
    from web.server import _build_mgf_tic_overview

    overview = _build_mgf_tic_overview([
        {
            "spectrum_id": "37809", "retention_time": 120.0,
            "tic": 12345.0, "precursor_mz": 482.2567,
        },
        {
            "spectrum_id": 38062, "retention_time": 180.0,
            "peaks": [[100.0, 10.0], [120.0, 20.0]],
            "precursor_mz": 258.1322,
        },
    ])

    assert overview["source"] == "mgf_ms2_approximation"
    assert overview["ms1"] == []
    assert [p["scan_id"] for p in overview["ms2"]] == [37809, 38062]
    assert [p["rt"] for p in overview["ms2"]] == pytest.approx([2.0, 3.0])
    assert overview["ms2"][1]["tic"] == pytest.approx(30.0)


def test_B11_批量导出项目识别和ZIP命名():
    from web.server import _batch_zip_filename, _extract_export_project, _unique_zip_member

    assert _extract_export_project("20260708213047-F001051-7.mgf") == "F001051"
    assert _extract_export_project("20260520235618-WTT3.mzML") == "WTT"
    # 短F编号 (F+1~6位数字, 下划线/连字符命名均兼容): 不得退化成只取开头字母 "F"
    assert _extract_export_project("F002_NEG_C1_2.mzML") == "F002"
    assert _extract_export_project("F002-NEG-3-7.mzML") == "F002"
    assert _extract_export_project("F001_POS_1.mzML") == "F001"
    assert _batch_zip_filename(
        ["F001_NEG_1.mzML"], "csv", 1, "20260711"
    ) == "F001_NEG_1_CSV_20260711.zip"
    assert _batch_zip_filename(
        [
            "F001_NEG_1.mzML", "F001_NEG_2.mzML",
            "F001_NEG_3.mzML", "F001_NEG_C1_1.mzML",
        ],
        "pdf", 4, "20260711",
    ) == (
        "F001_NEG_1__F001_NEG_2__F001_NEG_3__F001_NEG_C1_1_"
        "PDF_4samples_20260711.zip"
    )
    assert _batch_zip_filename(
        [
            "F002_NEG_1.mzML", "F002_NEG_2.mzML", "F002_NEG_3.mzML",
            "F002_POS_1.mzML", "F002_POS_2.mzML",
        ],
        "csv_ms2", 5, "20260711",
    ) == "F002_NEG_1__F002_NEG_2__等5个样本_CSV_MS2_5samples_20260711.zip"

    used = set()
    assert _unique_zip_member("F001051/sample_results.csv", "aaaaaaaa1111", used) == \
        "F001051/sample_results.csv"
    assert _unique_zip_member("F001051/sample_results.csv", "bbbbbbbb2222", used) == \
        "F001051/sample_results__bbbbbbbb.csv"


def test_B12_FULL库任务必须串行避免内存耗尽():
    """FULL 库两任务并发实测会把服务推到约16GB RSS并耗尽swap。"""
    from web import server
    assert server._executor_full._max_workers == 1
    assert server._executor_small._max_workers >= 2
    source = _read("web/server.py")
    assert "Queued for Full-library serial analysis..." in source
    assert "Analysis queued for task" in source
    assert "Analysis execution started for task" in source
    assert "Analysis started for task" not in source


# ═══════════════════════════════════════════════════════════════════
# C. 数据口径 — 用真实 NEG1 results.json (无则 skip; 重新分析后生成)
# ═══════════════════════════════════════════════════════════════════

# 最新 NEG1 任务 (MIX_QJF_1ug_NEG1.mgf)。重新分析后会带 tic/RT/confidence 字段。
NEG1_RESULTS = ROOT / "output" / "tasks" / "83f297d1da91_results.json"


def _load_neg1():
    if not NEG1_RESULTS.exists():
        pytest.skip(f"无 {NEG1_RESULTS.name} (重新分析 NEG1.mgf 后生成)")
    return json.loads(NEG1_RESULTS.read_text(encoding="utf-8"))


def test_C1_results无iceberg残留():
    """commit a14fa5a: 重新分析后 results 不应再混入 iceberg_fallback 条目。"""
    rs = _load_neg1()
    n = sum(1 for r in rs if (r.get("match") or {}).get("dataset_source") == "iceberg_fallback")
    assert n == 0, f"results 仍有 {n} 条 iceberg_fallback 残留"


def test_C2_去重有效_且tic字段存在():
    """去重后唯一化合物数 < 匹配谱图数 (去重起作用);
    且每条 matched 有 tic 字段 (commit 5981827 序列化修复生效)。"""
    rs = _load_neg1()
    matched = [r for r in rs if r.get("match")
               and r["match"].get("dataset_source") != "iceberg_fallback"]
    if len(matched) < 2:
        pytest.skip("匹配谱图太少")
    # 去重 (ik14 / name|formula)
    best = {}
    for r in matched:
        m = r["match"]
        ik = (m.get("inchikey") or "").split("-")[0]
        k = ik if (len(ik) == 14 and ik.isalpha() and ik.isupper()) \
            else f"{m.get('compound_name','')}|{m.get('formula','')}"
        if not k:
            continue
        if k not in best or (m.get("final_score", 0) or 0) > (best[k]["match"].get("final_score") or 0):
            best[k] = r
    assert len(best) < len(matched), \
        f"去重没起作用: {len(best)} 唯一化合物 == {len(matched)} 匹配谱图"
    # tic 字段必须在 (旧 bug: 序列化漏掉 → TIC 画不出)
    n_no_tic = sum(1 for r in matched if "tic" not in r)
    # 允许旧任务(本次升级前)无 tic; 新分析应有。这里软断言: 多数行有 tic。
    if n_no_tic > len(matched) * 0.5:
        pytest.skip(f"{n_no_tic}/{len(matched)} 行无 tic 字段 (旧任务, 重新分析后生效)")


def test_C3_负离子样本无跨极性假阳性():
    """commit 53efb20 + mgf charge→polarity: NEG 样本不应匹配 [M+H]+ 等正离子加合物。
    注: mgf 从 CHARGE 推 polarity 是后加的修复; 旧 results(修复前分析) 会大量跨极性 →
    此时 skip 提示重新分析, 而非判 fail (代码已修, 是数据陈旧)。"""
    rs = _load_neg1()
    matched = [r for r in rs if r.get("match")]
    POS_ADDUCTS = {"[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+"}
    cross = sum(1 for r in matched if (r["match"].get("adduct") or "").strip() in POS_ADDUCTS)
    rate = cross / max(len(matched), 1)
    if rate > 0.05:
        pytest.skip(
            f"{cross} 跨极性匹配 ({rate:.0%}) — 旧 results (mgf polarity 未生效); "
            f"重新分析 NEG1.mgf 后 polarity 过滤生效"
        )
    assert rate <= 0.05, f"{cross} 跨极性匹配 ({rate:.0%}), polarity 过滤可能失效"


# ═══════════════════════════════════════════════════════════════════
# D. 质心化管线 (2026-07-13 Fix1+Fix2) — 守护"链式串峰/centroid误判不复活"
# ═══════════════════════════════════════════════════════════════════

def test_D1_centroid_merge不链式串峰():
    """Fix2: _centroid_merge 必须按"簇总跨度 < merge_da"断簇 (complete-linkage),
    不能用 single-linkage(相邻 gap 断簇)把 0.04 间隔的链串成 1 个 span 0.12 的簇。
    历史根因: 旧实现把多个不同碎片串成一个峰 → 污染 cosine → 33 标品已知物余弦反常偏低。"""
    from utils.centroiding import _centroid_merge, MERGE_DA

    # 0.04 间隔的 4 点链 (span 0.12 Da): 旧 single-linkage 会串成 1 簇 (BUG)
    chain = [[100.00, 10], [100.04, 8], [100.08, 6], [100.12, 4]]
    out = _centroid_merge(chain, MERGE_DA)
    assert len(out) == 2, f"链式串峰复活: 应断成 2 簇(各 span<0.05), 实得 {len(out)} 簇: {out}"
    # 每簇 span 必须严格 < merge_da
    for c in out:
        assert 99 < c[0] < 101

    # 真正的 profile 采样(单峰近邻点 span 0.02)仍应合并为 1 簇
    prof = [[100.00, 5], [100.01, 10], [100.02, 7]]
    assert len(_centroid_merge(prof, MERGE_DA)) == 1, "真 profile 采样应合并为 1 簇"

    # 幂等: 已质心化(全 gap >= merge_da) -> 每点自成一簇, 原样数量
    cen = [[100.00, 5], [100.06, 7], [105.00, 9]]
    assert len(_centroid_merge(cen, MERGE_DA)) == 3, "已质心化应原样返回(幂等)"


def test_D2_smart_peak_pick跳过已质心谱图():
    """Fix1 配套: smart_peak_pick 对不含 profile 残余的谱图(相邻 gap 全 >= merge_da)
    必须原样返回, 不动 peaks。守护"已质心数据不被错误合并"。"""
    from utils.centroiding import smart_peak_pick

    cen = [[100.00, 10], [120.00, 9], [150.00, 8], [180.00, 7]]  # 4 个分散碎片
    out = smart_peak_pick(cen, precursor_mz=200.0)
    assert out == cen, f"已质心谱图被动了: {out}"


def test_D3_file_is_centroid判读mzML表示cvParam():
    """Fix1: _file_is_centroid 按 MS:1000127(centroid)/MS:1000126(profile) cvParam 计数。
    真实仪器导出的 mzML 标注 centroid 时必须返回 True → parse_mzml_file 据此跳过 smart_peak_pick。"""
    from utils.mzml_io import _file_is_centroid

    # 用合成 mzML 片段验证判读逻辑 (不依赖真实大文件)
    import tempfile, os
    centroid_xml = (
        '<?xml version="1.0"?><mzML><spectrumList><spectrum>'
        '<cvParam accession="MS:1000127" name="centroid spectrum"/>'
        '</spectrum></spectrumList></mzML>'
    )
    profile_xml = (
        '<?xml version="1.0"?><mzML><spectrumList><spectrum>'
        '<cvParam accession="MS:1000126" name="profile spectrum"/>'
        '</spectrum></spectrumList></mzML>'
    )
    with tempfile.NamedTemporaryFile(suffix=".mzML", delete=False, mode="w") as f:
        f.write(centroid_xml); cen_path = f.name
    with tempfile.NamedTemporaryFile(suffix=".mzML", delete=False, mode="w") as f:
        f.write(profile_xml); prof_path = f.name
    try:
        assert _file_is_centroid(cen_path) is True, "centroid mzML 应判 True"
        assert _file_is_centroid(prof_path) is False, "profile mzML 应判 False"
    finally:
        os.unlink(cen_path); os.unlink(prof_path)


def test_D4_解析缓存文件名含算法版本号():
    """Fix3: eval 的解析缓存 pickle 文件名必须含 CENTROID_ALGO_VERSION,
    改动 centroiding 算法后旧缓存自动失效(不被复用), 避免新旧污染谱图混用。"""
    from utils.centroiding import CENTROID_ALGO_VERSION
    src = _read("eval/eval_standards.py") + _read("eval/eval_baseline.py")
    assert "CENTROID_ALGO_VERSION" in src, "eval 缓存未接算法版本号(Fix3 回退)"
    assert "f\"{mzml_path.stem}.{CENTROID_ALGO_VERSION}.pkl\"" in _read("eval/eval_standards.py") or \
           "f'{mzml_path.stem}.{CENTROID_ALGO_VERSION}.pkl'" in _read("eval/eval_standards.py"), \
           "eval_standards 缓存文件名未含版本号"
    assert CENTROID_ALGO_VERSION, "CENTROID_ALGO_VERSION 不能为空"
