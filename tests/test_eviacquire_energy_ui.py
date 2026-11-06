from pathlib import Path


def test_eviacquire_ui_exposes_sequential_energy_boundaries():
    text = Path("web/static/js/eviacquire.js").read_text(encoding="utf-8")

    assert "首轮 10/20/40 eV" in text
    assert "40 eV 后仍欠碎裂才补 60 eV" in text
    assert "不声称所有化合物“完全碎裂”" in text
    assert "RT 标品证据与碰撞能优化独立" in text
    assert "条件分支" in text
    assert "本次按无标品资源模式生成" in text
    assert "已从动作空间、排序和导出中移除" in text
    assert "只有任务明确登记候选标品或离子淌度/CCS 参比" in text
    assert "原方法 Medium 约 4 amu" in text
    assert "仪器原生 Narrow 约 1.3 amu" in text
    assert "同方法 RT 证据缺失（未登记标品不生成动作）" in text
    index = Path("web/static/index.html").read_text(encoding="utf-8")
    assert "/static/js/eviacquire.js?v=10.4" in index
