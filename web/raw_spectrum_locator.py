"""Human-readable locators for returning from reports to the source MS/MS.

The locator deliberately describes the observed query spectrum, not the
library spectrum.  It keeps retention time, precursor m/z, and the mzML
spectrum/native identifier together so a reviewer can find the same event in
vendor software or an mzML viewer without interpreting internal result IDs.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping, Optional


def _finite_number(value: Any) -> Optional[float]:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def retention_time_minutes(row: Mapping[str, Any]) -> Optional[float]:
    """Return the observed MS2 retention time in minutes.

    Current result artifacts store RT in minutes.  Explicit second/hour units
    are nevertheless normalized here so legacy or imported rows cannot be
    silently mislabeled in a report.
    """

    value = _finite_number(row.get("retention_time"))
    if value is None:
        return None
    unit = str(row.get("retention_time_unit") or "minutes").strip().lower()
    if "millisecond" in unit or unit in {"ms", "msec", "msecs"}:
        return value / 60_000.0
    if "second" in unit or unit in {"s", "sec", "secs"}:
        return value / 60.0
    if "hour" in unit or unit in {"h", "hr", "hrs"}:
        return value * 60.0
    return value


def format_rt_clock(minutes: Optional[float]) -> str:
    """Format minutes as a chromatogram-friendly ``MM:SS.s`` value."""

    if minutes is None or not math.isfinite(minutes) or minutes < 0:
        return ""
    total_seconds = minutes * 60.0
    whole_minutes = int(total_seconds // 60)
    seconds = total_seconds - whole_minutes * 60
    return f"{whole_minutes:02d}:{seconds:04.1f}"


def raw_spectrum_locator(
    row: Mapping[str, Any], source_name: Optional[str] = None
) -> dict[str, Any]:
    """Build a machine- and reader-friendly locator for one observed MS2."""

    rt_minutes = retention_time_minutes(row)
    source_file = Path(str(source_name)).name if source_name else ""
    spectrum_id = str(row.get("spectrum_id") or "").strip()
    precursor_scan_ref = str(row.get("precursor_scan_ref") or "").strip()
    selected_ion_mz = _finite_number(row.get("precursor_mz"))
    search_precursor_mz = _finite_number(row.get("search_precursor_mz"))
    search_source = str(row.get("search_precursor_mz_source") or "").strip()
    feature_used = bool(row.get("lcms_feature_used_for_search"))
    if search_precursor_mz is None:
        search_precursor_mz = selected_ion_mz
    if not feature_used:
        feature_used = search_source == "measured_continuous_ms1_feature_centroid"
    return {
        "source_file": source_file,
        "retention_time_min": rt_minutes,
        "retention_time_clock": format_rt_clock(rt_minutes),
        # ``precursor_mz`` is retained as a compatibility alias for the mzML
        # selected-ion/isolation target.  It must not be relabelled as a
        # measured HRMS centroid when a continuous-MS1 feature was used.
        "precursor_mz": selected_ion_mz,
        "selected_ion_mz": selected_ion_mz,
        "search_precursor_mz": search_precursor_mz,
        "search_precursor_mz_source": search_source,
        "lcms_feature_used_for_search": feature_used,
        "search_precursor_delta_ppm": _finite_number(
            row.get("search_precursor_delta_ppm")
        ),
        "spectrum_id": spectrum_id,
        "precursor_scan_ref": precursor_scan_ref,
        "collision_energy": _finite_number(row.get("collision_energy")),
        "collision_energy_unit": str(
            row.get("collision_energy_unit") or ""
        ).strip(),
    }


def raw_spectrum_locator_text(
    row: Mapping[str, Any],
    source_name: Optional[str] = None,
    language: str = "zh",
) -> str:
    """Render the locator as plain text suitable for HTML/PDF/CSV cells."""

    locator = raw_spectrum_locator(row, source_name)
    english = str(language).lower().startswith("en")
    parts: list[str] = []
    if locator["source_file"]:
        parts.append(
            ("Source file " if english else "源文件 ")
            + str(locator["source_file"])
        )
    rt_minutes = locator["retention_time_min"]
    if rt_minutes is None:
        parts.append("RT not recorded" if english else "RT 未记录")
    else:
        clock = locator["retention_time_clock"]
        parts.append(f"RT {rt_minutes:.3f} min ({clock})")
    selected_ion_mz = locator["selected_ion_mz"]
    search_precursor_mz = locator["search_precursor_mz"]
    if locator["lcms_feature_used_for_search"]:
        if selected_ion_mz is not None:
            parts.append(
                (
                    "MS2 isolation target m/z "
                    if english else "MS2 隔离目标 m/z "
                )
                + f"{selected_ion_mz:.5f}"
            )
        if search_precursor_mz is not None:
            parts.append(
                (
                    "measured MS1 search m/z "
                    if english else "检索用实测 MS1 m/z "
                )
                + f"{search_precursor_mz:.5f}"
            )
    elif selected_ion_mz is not None:
        parts.append(f"precursor m/z {selected_ion_mz:.5f}")
    spectrum_id = locator["spectrum_id"]
    parts.append(
        (
            f"MS2 spectrum/native ID {spectrum_id}"
            if english
            else f"MS2 谱图编号/native ID {spectrum_id}"
        )
        if spectrum_id
        else ("MS2 spectrum ID not recorded" if english else "MS2 谱图编号未记录")
    )
    if locator["precursor_scan_ref"]:
        parts.append(
            ("linked preceding scan " if english else "关联前级扫描 ")
            + str(locator["precursor_scan_ref"])
        )
    return " · ".join(parts)


__all__ = [
    "format_rt_clock",
    "raw_spectrum_locator",
    "raw_spectrum_locator_text",
    "retention_time_minutes",
]
