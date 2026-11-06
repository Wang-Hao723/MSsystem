"""Export saved peak arrays without spreadsheet cell truncation or rounding.

The saved result-array position is the key, even when spectrum titles repeat.
Missing or historically incomplete arrays stay explicitly marked; this exporter
does not recreate peaks absent from the saved task.
"""
from __future__ import annotations

import hashlib
import json
import zlib
from collections import defaultdict
from collections.abc import Iterable, Mapping
from typing import Any

CHUNK_CHARACTERS = 30000
EXCEL_MAX_ROWS = 1048576
SCHEMA = "mssystem-saved-raw-peaks-1.0"


def attachment_disposition(filename: str) -> str:
    from urllib.parse import quote
    import re
    fallback = re.sub(r'[^\x20-\x7e]|["\\\\]', "_", filename)
    return f'attachment; filename="{fallback}"; filename*=UTF-8\'\'{quote(filename, safe="")}'


def canonical_peak_json(raw_peaks: Any) -> str:
    if raw_peaks is not None and not isinstance(raw_peaks, list):
        raise ValueError("saved raw_peaks must be an array or null")
    # ASCII also avoids invalid Excel XML control characters in JSON content.
    return json.dumps(raw_peaks, ensure_ascii=True, separators=(",", ":"), allow_nan=False)


def raw_peak_record(row: Mapping[str, Any], index: int, task_id: str,
                    source_results_sha256: str | None = None) -> dict[str, Any]:
    raw = row.get("raw_peaks")
    text = canonical_peak_json(raw)
    return {
        "schema": SCHEMA, "task_id": task_id, "result_index": index,
        "spectrum_id": row.get("spectrum_id"), "ion_mode": row.get("ion_mode"),
        "precursor_mz": row.get("precursor_mz"),
        "retention_time": row.get("retention_time"),
        "raw_peaks_present": "raw_peaks" in row,
        "raw_peaks_input_count": row.get("raw_peaks_input_count"),
        "raw_peaks_saved_count": row.get("raw_peaks_saved_count"),
        "raw_peaks_complete": row.get("raw_peaks_complete"),
        "peak_count": len(raw) if isinstance(raw, list) else 0,
        "raw_peaks_sha256": hashlib.sha256(text.encode("ascii")).hexdigest(),
        "source_results_sha256": source_results_sha256,
        "raw_peaks": raw,
    }


def iter_raw_peak_jsonl_gzip(rows: Iterable[Mapping[str, Any]], task_id: str,
                            source_results_sha256: str | None = None):
    """Yield a single gzip member, one complete saved spectrum per JSON line."""
    compressor = zlib.compressobj(level=6, wbits=31)
    for index, row in enumerate(rows):
        record = raw_peak_record(row, index, task_id, source_results_sha256)
        line = json.dumps(record, ensure_ascii=True, separators=(",", ":"),
                          allow_nan=False).encode("ascii") + b"\n"
        block = compressor.compress(line)
        if block:
            yield block
    yield compressor.flush()


def peak_sheet_names(language: str = "zh") -> tuple[str, str]:
    return (("14_Raw_Peak_Index", "15_Raw_Peak_JSON") if language == "en"
            else ("14_原始峰索引", "15_原始峰JSON"))


def prepare_reader_peak_cells(rows, *, task_id: str = "", results=None,
                              language: str = "zh"):
    """Shorten only oversized display cells, naming the complete download."""
    indices = defaultdict(list)
    for index, result in enumerate(results or []):
        indices[str(result.get("spectrum_id", ""))].append(index)
    index_sheet, chunk_sheet = peak_sheet_names(language)
    output = []
    for original in rows:
        row = dict(original)
        raw_text = str(row.get("Raw_MS_MS_Spectrum") or "")
        if len(raw_text) > CHUNK_CHARACTERS:
            indexes = indices.get(str(row.get("Spectrum", "")), [])
            index_note = ",".join(map(str, indexes[:20]))
            if len(indexes) > 20:
                index_note += ",…"
            if language == "en":
                note = (f"Saved peak preview only; complete arrays: {index_sheet} / "
                        f"{chunk_sheet}, or Reports & Exports > Complete raw peaks "
                        f"(JSONL.gz). Match Spectrum; result indices: {index_note or 'see index'}. "
                        "Repeated spectrum IDs retain separate result indices.")
            else:
                note = (f"此处为已保存峰的节选；完整数组见「{index_sheet}」「{chunk_sheet}」，"
                        "或在报告与导出下载「完整原始峰 JSONL.gz」。"
                        f"按样本谱图编号关联；结果数组索引：{index_note or '见索引表'}。"
                        "重复谱图编号须区分结果数组索引。")
            if task_id:
                note += f" /api/tasks/{task_id}/download/raw-peaks"
            row["Raw_MS_MS_Spectrum"] = raw_text[:2000] + " …[" + note + "]"
            old_note = str(row.get("Export_Note") or "")
            obsolete = ("raw peak list exceeds one spreadsheet cell; the CSV cell is explicitly truncated"
                        if language == "en" else "原始峰列表超过单个表格单元格容量；该单元格将明确截断")
            old_note = old_note.replace(obsolete, "").strip("; ")
            row["Export_Note"] = "; ".join(filter(None, [old_note, note]))
        output.append(row)
    return output


def add_raw_peak_sheets(workbook, rows: Iterable[Mapping[str, Any]], task_id: str,
                        *, language: str = "zh", source_results_sha256: str | None = None):
    """Append indexed JSON chunks; never pass oversized strings to openpyxl."""
    from openpyxl.styles import Alignment, Font, PatternFill

    index_name, chunk_name = peak_sheet_names(language)
    common = ["Task_ID", "Result_Index_0Based", "Spectrum", "Ion_Mode", "Saved_Peak_Count",
              "Input_Peak_Count", "Saved_Array_Complete", "Raw_Field_Present", "Raw_Peaks_SHA256"]
    index_headers = common + ["Chunk_Count", "Chunk_Sheet_First", "Chunk_Row_First",
                              "Chunk_Sheet_Last", "Chunk_Row_Last", "Source_Results_SHA256"]
    chunk_headers = common + ["Chunk_Number_1Based", "Chunk_Count", "JSON_Chunk"]
    pages = {"index": [], "chunk": []}

    def new_sheet(kind):
        base = index_name if kind == "index" else chunk_name
        suffix = "" if not pages[kind] else f"_{len(pages[kind]) + 1}"
        ws = workbook.create_sheet(base + suffix)
        ws.append(index_headers if kind == "index" else chunk_headers)
        ws.freeze_panes = "C2"
        ws.sheet_view.showGridLines = False
        for cell in ws[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="34495E")
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[1].height = 32
        from openpyxl.utils import get_column_letter
        for col in range(1, ws.max_column + 1):
            ws.column_dimensions[get_column_letter(col)].width = 22
        ws.column_dimensions["C"].width = 35
        ws.column_dimensions["I"].width = 68
        if kind == "chunk":
            ws.column_dimensions["L"].width = 80
        pages[kind].append(ws)
        return ws

    def append(kind, values):
        ws = pages[kind][-1] if pages[kind] else new_sheet(kind)
        if ws.max_row >= EXCEL_MAX_ROWS:
            ws = new_sheet(kind)
        if any(isinstance(value, str) and len(value) > CHUNK_CHARACTERS for value in values):
            raise ValueError("raw peak export metadata exceeds the safe Excel cell length")
        ws.append(values)
        # A scan title beginning '=' is an identifier, not a formula.
        for cell, value in zip(ws[ws.max_row], values):
            if isinstance(value, str):
                cell.data_type = "s"
        return ws.title, ws.max_row

    new_sheet("index")
    new_sheet("chunk")
    spectra = peaks = chunks = 0
    for index, row in enumerate(rows):
        record = raw_peak_record(row, index, task_id, source_results_sha256)
        text = canonical_peak_json(record["raw_peaks"])
        parts = [text[pos:pos + CHUNK_CHARACTERS] for pos in range(0, len(text), CHUNK_CHARACTERS)]
        values = [task_id, index, record["spectrum_id"], record["ion_mode"], record["peak_count"],
                  record["raw_peaks_input_count"], record["raw_peaks_complete"],
                  record["raw_peaks_present"], record["raw_peaks_sha256"]]
        first = last = None
        for number, part in enumerate(parts, 1):
            last = append("chunk", values + [number, len(parts), part])
            first = first or last
        append("index", values + [len(parts), *first, *last, source_results_sha256])
        spectra += 1
        peaks += record["peak_count"]
        chunks += len(parts)
    for sheets in pages.values():
        for ws in sheets:
            ws.auto_filter.ref = ws.dimensions
    return {"schema": SCHEMA, "spectra": spectra, "saved_peaks": peaks,
            "chunks": chunks, "chunk_characters": CHUNK_CHARACTERS,
            "index_sheets": [ws.title for ws in pages["index"]],
            "chunk_sheets": [ws.title for ws in pages["chunk"]]}
