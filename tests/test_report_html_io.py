from pathlib import Path

import pytest

from utils.mgf_io import parse_mgf_file
from utils.report_html_io import (
    is_raw_report_html,
    parse_raw_report_html,
    records_to_mgf,
)


SYNTHETIC_REPORT = """<!doctype html><html><body>
<p><strong>样品:</strong> D:\\data\\sample.d</p>
<section class="peak-card">
  <h3>UV峰#1: RT=7.515min (7.49-7.67)</h3>
  <div class="ms1-meta">离子模式: <strong>NEG</strong> | scan_id: 1~2</div>
  <h4>过滤后 MS1 top5</h4>
  <table><tr><td>1</td><td><button data-mz="999.9">EIC</button></td><td>100%</td></tr></table>
  <h4>MS2 碎片</h4>
  <div class="ms2-header">前体 m/z=<span>729.2604</span> (scan=452672)</div>
  <table>
    <tr><th>#</th><th>m/z</th><th>%</th></tr>
    <tr><td>1</td><td><button data-mz="123.0451">EIC</button></td><td>100.0%</td></tr>
    <tr><td>2</td><td><button data-mz="557.1844">EIC</button></td><td>84.4%</td></tr>
  </table>
</section></body></html>"""


def test_report_detection_and_peak_card_parsing(tmp_path):
    report = tmp_path / "sample_report.html"
    report.write_text(SYNTHETIC_REPORT, encoding="utf-8")

    assert is_raw_report_html(report)
    records = parse_raw_report_html(report)

    assert len(records) == 1
    spec = records[0]
    assert spec["spectrum_id"] == 452672
    assert spec["scan_number"] == 452672
    assert spec["precursor_mz"] == pytest.approx(729.2604)
    assert spec["retention_time"] == pytest.approx(7.515)
    assert spec["polarity"] == "NEG"
    assert spec["polarity_source"] == "report_meta"
    assert spec["input_format"] == "raw_blank_analysis_html"
    assert spec["limited_fragment_list"] is True
    assert "仅提供 2 个 MS2 碎片" in spec["input_evidence_note"]
    assert spec["peaks"] == [[123.0451, 100.0], [557.1844, 84.4]]
    assert all(mz != 999.9 for mz, _ in spec["peaks"]), "MS1 表不能误当 MS2"


def test_arbitrary_html_is_rejected(tmp_path):
    page = tmp_path / "ordinary.html"
    page.write_text("<html><body><h1>普通网页</h1></body></html>", encoding="utf-8")
    assert not is_raw_report_html(page)
    assert parse_raw_report_html(page) == []


def test_report_to_mgf_roundtrip_preserves_negative_mode_and_rt(tmp_path):
    report = tmp_path / "sample_report.html"
    report.write_text(SYNTHETIC_REPORT, encoding="utf-8")
    mgf = tmp_path / "sample_report.mgf"

    assert records_to_mgf(parse_raw_report_html(report), mgf) == 1
    text = mgf.read_text(encoding="utf-8")
    assert "CHARGE=1-" in text
    assert "IONMODE=NEG" in text
    assert "RTINSECONDS=450.9" in text

    parsed = parse_mgf_file(mgf)
    assert len(parsed) == 1
    assert parsed[0]["polarity"] == "NEG"
    assert parsed[0]["retention_time"] == pytest.approx(450.9)


def test_real_m011_report_extracts_ten_negative_ms2_if_available():
    report = Path("/Users/zizi/Desktop/M011_results/sample_report.html")
    if not report.exists():
        pytest.skip("local M011 report not present")
    records = parse_raw_report_html(report)
    assert len(records) == 10
    assert all(r["polarity"] == "NEG" for r in records)
    assert all(r["limited_fragment_list"] is True for r in records)
    assert all("仅提供 5 个 MS2 碎片" in r["input_evidence_note"] for r in records)
    assert {r["spectrum_id"] for r in records} == {
        452672, 407106, 446918, 419299, 428251,
        351458, 399817, 438740, 377704, 481108,
    }
