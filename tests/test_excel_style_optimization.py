from __future__ import annotations

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from web.expected_compounds_export import style_analysis_workbook


HEADERS = [
    "证据等级", "主要候选", "候选数", "理论 m/z", "实测 m/z",
    "误差 (ppm)", "主峰 RT (min)", "代表峰面积", "信噪比",
    "Final Score", "样本谱图编号", "备注",
]
WIDTHS = {
    "证据等级": 23, "主要候选": 38, "候选数": 9,
    "理论 m/z": 14, "实测 m/z": 14, "误差 (ppm)": 12,
    "主峰 RT (min)": 13, "代表峰面积": 15, "信噪比": 11,
    "Final Score": 12, "样本谱图编号": 21, "备注": 40,
}
FORMATS = {
    "候选数": "0", "理论 m/z": "0.000000", "实测 m/z": "0.000000",
    "误差 (ppm)": "0.00", "主峰 RT (min)": "0.000",
    "代表峰面积": "0.00E+00", "信噪比": "0.00",
    "Final Score": "0.000",
}
EVIDENCE_COLORS = {
    "MS/MS 实验谱库候选支持": "C6E0B4",
    "MS1 + 诊断离子筛查支持": "BDD7EE",
    "MS1 + 同位素支持": "FFF2CC",
    "MS1 精确质量候选": "E7E6E6",
}


def make_workbook(rows=600):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "03_MS2谱库鉴定"
    sheet.append(HEADERS)
    evidence = list(EVIDENCE_COLORS)
    for index in range(rows):
        sheet.append([
            evidence[index % len(evidence)], f"candidate-{index % 37}",
            True if index == 17 else index % 7 + 1,
            100.0 + index / 10000, 100.0 + index / 9999,
            (index % 31 - 15) / 10, index / 1000, float(index * 101),
            (index % 97) / 3, (index % 1000) / 1000,
            f"scan={index + 1}", "=1+1" if index == 21 else "fixture",
        ])
    sheet.cell(3, 2).font = Font(name="Courier New", italic=True, color="7030A0")
    sheet.cell(4, 3).border = Border(left=Side(style="thin", color="FF0000"))
    return workbook


def legacy_style_oracle(workbook):
    sheet = workbook.active
    header_fill = PatternFill("solid", fgColor="17365D")
    header_font = Font(color="FFFFFF", bold=True, size=10)
    header_border = Border(bottom=Side(style="medium", color="9EADBC"))
    band_fill = PatternFill("solid", fgColor="F4F7FB")
    text_alignment = Alignment(horizontal="left", vertical="center", wrap_text=False)
    number_alignment = Alignment(horizontal="right", vertical="center", wrap_text=False)
    evidence_fills = {
        name: PatternFill("solid", fgColor=color)
        for name, color in EVIDENCE_COLORS.items()
    }
    sheet.sheet_view.showGridLines = False
    sheet.sheet_view.zoomScale = 90
    sheet.freeze_panes = "B2"
    sheet.sheet_properties.tabColor = "5B9BD5"
    sheet.sheet_properties.pageSetUpPr.fitToPage = True
    sheet.page_setup.fitToWidth = 1
    sheet.page_setup.fitToHeight = 0
    sheet.print_title_rows = "1:1"
    sheet.row_dimensions[1].height = 30
    headers = {cell.column: str(cell.value or "") for cell in sheet[1]}
    for cell in sheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.border = header_border
        cell.alignment = Alignment(
            horizontal="center", vertical="center", wrap_text=True
        )
    sheet.auto_filter.ref = sheet.dimensions
    sample_limit = min(sheet.max_row, 200)
    for row_index in range(2, sheet.max_row + 1):
        sheet.row_dimensions[row_index].height = 21
        if row_index % 2 == 0:
            for cell in sheet[row_index]:
                cell.fill = band_fill
    for column_index in range(1, sheet.max_column + 1):
        header = headers.get(column_index, "")
        letter = get_column_letter(column_index)
        if header in WIDTHS:
            width = WIDTHS[header]
        else:
            observed = len(header)
            for row_index in range(2, sample_limit + 1):
                value = sheet.cell(row_index, column_index).value
                if value is not None:
                    observed = max(observed, len(str(value)))
            width = min(max(observed + 2, 10), 22)
        sheet.column_dimensions[letter].width = width
        number_format = FORMATS.get(header)
        for row_index in range(2, sheet.max_row + 1):
            cell = sheet.cell(row_index, column_index)
            is_number = isinstance(cell.value, (int, float))
            cell.alignment = number_alignment if is_number else text_alignment
            if number_format and is_number:
                cell.number_format = number_format
    for row_index in range(2, sheet.max_row + 1):
        cell = sheet.cell(row_index, 1)
        fill = evidence_fills.get(str(cell.value or ""))
        if fill:
            cell.fill = fill


def assert_workbooks_equal(left, right):
    a, b = left.active, right.active
    assert a.max_row == b.max_row and a.max_column == b.max_column
    assert a.freeze_panes == b.freeze_panes
    assert a.auto_filter.ref == b.auto_filter.ref
    assert a.sheet_view.showGridLines == b.sheet_view.showGridLines
    assert a.sheet_view.zoomScale == b.sheet_view.zoomScale
    assert a.print_title_rows == b.print_title_rows
    assert a.page_setup.fitToWidth == b.page_setup.fitToWidth
    assert a.page_setup.fitToHeight == b.page_setup.fitToHeight
    for row_index in range(1, a.max_row + 1):
        assert a.row_dimensions[row_index].height == b.row_dimensions[row_index].height
        for column_index in range(1, a.max_column + 1):
            left_cell = a.cell(row_index, column_index)
            right_cell = b.cell(row_index, column_index)
            assert left_cell.value == right_cell.value
            assert tuple(left_cell._style) == tuple(right_cell._style)
    for column_index in range(1, a.max_column + 1):
        letter = get_column_letter(column_index)
        assert a.column_dimensions[letter].width == b.column_dimensions[letter].width


def test_cached_single_pass_styles_are_exactly_legacy_equivalent(tmp_path):
    legacy = make_workbook()
    optimized = make_workbook()
    legacy_style_oracle(legacy)
    style_analysis_workbook(optimized)
    assert_workbooks_equal(legacy, optimized)

    legacy_path = tmp_path / "legacy.xlsx"
    optimized_path = tmp_path / "optimized.xlsx"
    legacy.save(legacy_path)
    optimized.save(optimized_path)
    assert_workbooks_equal(
        load_workbook(legacy_path, data_only=False),
        load_workbook(optimized_path, data_only=False),
    )


def test_cached_style_arrays_remain_independently_editable():
    workbook = make_workbook(rows=20)
    style_analysis_workbook(workbook)
    sheet = workbook.active
    changed = sheet.cell(3, 11)
    unchanged = sheet.cell(3, 12)
    assert tuple(changed._style) == tuple(unchanged._style)
    assert changed._style is not unchanged._style
    unchanged_before = tuple(unchanged._style)
    changed.font = Font(name="Aptos", bold=True, color="FF0000")
    changed.number_format = "0.00000"
    assert tuple(unchanged._style) == unchanged_before
    assert tuple(changed._style) != unchanged_before
