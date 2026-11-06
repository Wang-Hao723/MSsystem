"""Auditable, configurable product-ion processing.

The legacy ``clean_peaks`` function is kept as a compatibility wrapper, but the
implementation now exposes every destructive stage.  This is important for
untargeted workflows: a Top-N view used for library matching must not be
confused with centroiding, feature detection, or the richer view used to build
molecular networks.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any, Iterable

import numpy as np

from config.settings import REL_THRESHOLD, TOP_N_PEAKS


@dataclass(frozen=True)
class PeakProcessingConfig:
    """Configuration for one symmetric query/library MS2 processing view."""

    preset: str = "library_match"
    relative_intensity_threshold: float = REL_THRESHOLD
    precursor_exclusion_da: float = 0.5
    top_n: int = TOP_N_PEAKS
    upper_mass_bound_enabled: bool = True

    def validated(self) -> "PeakProcessingConfig":
        if not 0.0 <= float(self.relative_intensity_threshold) <= 1.0:
            raise ValueError("relative_intensity_threshold must be between 0 and 1")
        if not 0.0 <= float(self.precursor_exclusion_da) <= 10.0:
            raise ValueError("precursor_exclusion_da must be between 0 and 10 Da")
        if int(self.top_n) < 3 or int(self.top_n) > 500:
            raise ValueError("top_n must be between 3 and 500")
        return replace(
            self,
            relative_intensity_threshold=float(self.relative_intensity_threshold),
            precursor_exclusion_da=float(self.precursor_exclusion_da),
            top_n=int(self.top_n),
            upper_mass_bound_enabled=bool(self.upper_mass_bound_enabled),
        )

    def as_dict(self) -> dict[str, Any]:
        return asdict(self.validated())

    @property
    def cache_key(self) -> tuple[Any, ...]:
        value = self.validated()
        return (
            round(value.relative_intensity_threshold, 8),
            round(value.precursor_exclusion_da, 8),
            value.top_n,
            value.upper_mass_bound_enabled,
        )

    @property
    def is_legacy_default(self) -> bool:
        return self.cache_key == PeakProcessingConfig().cache_key


PEAK_PROCESSING_PRESETS: dict[str, PeakProcessingConfig] = {
    # V6-compatible search view.  It remains the default so existing score
    # calibration is not silently changed.
    "library_match": PeakProcessingConfig(
        preset="library_match",
        relative_intensity_threshold=0.05,
        precursor_exclusion_da=0.5,
        top_n=15,
        upper_mass_bound_enabled=True,
    ),
    # Richer view for exploratory matching of sparse natural-product spectra.
    "sensitive": PeakProcessingConfig(
        preset="sensitive",
        relative_intensity_threshold=0.01,
        precursor_exclusion_da=0.5,
        top_n=40,
        upper_mass_bound_enabled=True,
    ),
    # The established V7 molecular-network view.
    "network_rich": PeakProcessingConfig(
        preset="network_rich",
        relative_intensity_threshold=0.01,
        precursor_exclusion_da=0.5,
        top_n=70,
        upper_mass_bound_enabled=True,
    ),
}


def resolve_peak_processing_config(
    preset: str | None = None,
    *,
    relative_intensity_threshold: float | None = None,
    precursor_exclusion_da: float | None = None,
    top_n: int | None = None,
    upper_mass_bound_enabled: bool | None = None,
) -> PeakProcessingConfig:
    """Resolve a named preset with optional explicit overrides."""

    key = str(preset or "library_match").strip().lower()
    if key not in PEAK_PROCESSING_PRESETS:
        raise ValueError(
            f"Unknown peak-processing preset {key!r}; "
            f"choose one of {sorted(PEAK_PROCESSING_PRESETS)}"
        )
    config = PEAK_PROCESSING_PRESETS[key]
    return replace(
        config,
        relative_intensity_threshold=(
            config.relative_intensity_threshold
            if relative_intensity_threshold is None
            else float(relative_intensity_threshold)
        ),
        precursor_exclusion_da=(
            config.precursor_exclusion_da
            if precursor_exclusion_da is None
            else float(precursor_exclusion_da)
        ),
        top_n=config.top_n if top_n is None else int(top_n),
        upper_mass_bound_enabled=(
            config.upper_mass_bound_enabled
            if upper_mass_bound_enabled is None
            else bool(upper_mass_bound_enabled)
        ),
    ).validated()


def _absolute_charge(value):
    """Best-effort absolute precursor charge for cleaning bounds."""
    if value is None or isinstance(value, bool):
        return 1
    text = str(value).strip()
    if not text:
        return 1
    try:
        charge = abs(int(float(text.rstrip("+-"))))
    except (TypeError, ValueError):
        return 1
    return max(1, charge)


def fragment_upper_mz(precursor_mz, precursor_charge=None):
    """Upper product-ion bound that respects multiply charged precursors."""
    try:
        precursor = float(precursor_mz)
    except (TypeError, ValueError):
        return 0.0
    return precursor * _absolute_charge(precursor_charge) + 1.5


def _stage(name: str, input_count: int, output_count: int) -> dict[str, Any]:
    return {
        "stage": name,
        "input_peaks": int(input_count),
        "retained_peaks": int(output_count),
        "retained_from_previous_pct": (
            round(100.0 * output_count / input_count, 3) if input_count else None
        ),
    }


def process_peaks(
    peaks: Iterable[Iterable[Any]] | None,
    precursor_mz: Any,
    *,
    config: PeakProcessingConfig | None = None,
    precursor_charge: Any = None,
) -> tuple[list[list[float]], dict[str, Any]]:
    """Apply the configured MS2 stages and return ``(peaks, audit)``.

    Centroid/profile conversion is intentionally not performed here.  This
    function receives a peak list after representation-aware vendor/OpenMS
    processing and only creates a reproducible comparison view.
    """

    resolved = (config or PeakProcessingConfig()).validated()
    source = list(peaks or [])
    stages: list[dict[str, Any]] = []

    valid: list[tuple[float, float]] = []
    for peak in source:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if (
            np.isfinite(mz)
            and np.isfinite(intensity)
            and mz > 0
            and intensity > 0
        ):
            valid.append((mz, intensity))
    stages.append(_stage("valid_numeric_positive", len(source), len(valid)))

    if not valid:
        return [], {
            "schema_version": 2,
            "config": resolved.as_dict(),
            "input_peaks": len(source),
            "output_peaks": 0,
            "retained_pct": 0.0 if source else None,
            "precursor_dominance_rescue_triggered": False,
            "precursor_dominance_rescue_reason": "no_valid_peaks",
            "stages": stages,
        }

    mz_arr = np.asarray([item[0] for item in valid], dtype=float)
    int_arr = np.asarray([item[1] for item in valid], dtype=float)
    # Keep an immutable validated view for the narrow precursor-dominance
    # rescue below.  Some DDA spectra contain an unfragmented precursor that
    # is orders of magnitude stronger than every product ion.  Applying a
    # relative threshold against that precursor and removing it afterwards
    # can otherwise erase a measured product-ion spectrum completely.
    valid_mz_arr = mz_arr.copy()
    valid_int_arr = int_arr.copy()
    max_int = float(np.max(int_arr))
    relative_keep = int_arr >= max_int * resolved.relative_intensity_threshold
    mz_arr = mz_arr[relative_keep]
    int_arr = int_arr[relative_keep]
    stages.append(
        _stage("relative_intensity", len(valid), int(mz_arr.size))
    )

    try:
        precursor = float(precursor_mz)
    except (TypeError, ValueError):
        precursor = 0.0
    before = int(mz_arr.size)
    if precursor > 0 and resolved.precursor_exclusion_da > 0:
        keep = np.abs(mz_arr - precursor) > resolved.precursor_exclusion_da
        mz_arr = mz_arr[keep]
        int_arr = int_arr[keep]
    stages.append(_stage("precursor_exclusion", before, int(mz_arr.size)))

    before = int(mz_arr.size)
    if precursor > 0 and resolved.upper_mass_bound_enabled:
        upper_mz = fragment_upper_mz(precursor, precursor_charge)
        keep = mz_arr <= upper_mz
        mz_arr = mz_arr[keep]
        int_arr = int_arr[keep]
    stages.append(_stage("upper_mass_bound", before, int(mz_arr.size)))

    precursor_dominance_rescue_triggered = False
    precursor_dominance_rescue_reason = "not_needed"
    # Guarded fallback: only reconsider the validated product ions when the
    # ordinary view has fewer than three fragments, at least three measured
    # product ions exist, and an excluded precursor-region peak dominates the
    # strongest product ion.  The same configured relative threshold and Top-N
    # limit still apply; the threshold reference alone becomes the strongest
    # product ion.  This does not lower the evidence/release gates.
    if precursor > 0 and int(mz_arr.size) < 3:
        product_keep = (
            np.abs(valid_mz_arr - precursor) > resolved.precursor_exclusion_da
        )
        if resolved.upper_mass_bound_enabled:
            product_keep &= valid_mz_arr <= fragment_upper_mz(
                precursor, precursor_charge
            )
        product_mz = valid_mz_arr[product_keep]
        product_int = valid_int_arr[product_keep]
        precursor_region_keep = (
            np.abs(valid_mz_arr - precursor)
            <= resolved.precursor_exclusion_da
        )
        precursor_region_int = valid_int_arr[precursor_region_keep]
        if product_int.size >= 3 and precursor_region_int.size:
            product_max = float(np.max(product_int))
            precursor_region_max = float(np.max(precursor_region_int))
            if product_max > 0 and precursor_region_max > product_max:
                rescue_keep = (
                    product_int
                    >= product_max * resolved.relative_intensity_threshold
                )
                rescued_mz = product_mz[rescue_keep]
                rescued_int = product_int[rescue_keep]
                if rescued_mz.size >= 3 and rescued_mz.size > mz_arr.size:
                    before_rescue = int(mz_arr.size)
                    mz_arr = rescued_mz
                    int_arr = rescued_int
                    precursor_dominance_rescue_triggered = True
                    precursor_dominance_rescue_reason = (
                        "excluded_precursor_dominated_relative_threshold"
                    )
                    stages.append(
                        _stage(
                            "precursor_dominance_rescue",
                            before_rescue,
                            int(mz_arr.size),
                        )
                    )

    before = int(mz_arr.size)
    if before > resolved.top_n:
        top_idx = np.argsort(int_arr, kind="stable")[::-1][: resolved.top_n]
        mz_arr = mz_arr[top_idx]
        int_arr = int_arr[top_idx]
    stages.append(_stage("top_n", before, int(mz_arr.size)))

    order = np.argsort(mz_arr, kind="stable")
    output = [
        [float(mz_arr[index]), float(int_arr[index])]
        for index in order
    ]
    return output, {
        "schema_version": 2,
        "config": resolved.as_dict(),
        "input_peaks": len(source),
        "output_peaks": len(output),
        "retained_pct": (
            round(100.0 * len(output) / len(source), 3) if source else None
        ),
        "precursor_dominance_rescue_triggered": (
            precursor_dominance_rescue_triggered
        ),
        "precursor_dominance_rescue_reason": (
            precursor_dominance_rescue_reason
        ),
        "stages": stages,
    }


def aggregate_peak_processing_audits(
    audits: Iterable[dict[str, Any]],
) -> dict[str, Any]:
    """Aggregate per-spectrum stage counts without averaging percentages."""

    audit_list = [item for item in audits if isinstance(item, dict)]
    stage_totals: dict[str, dict[str, int]] = {}
    for audit in audit_list:
        for stage in audit.get("stages") or []:
            name = str(stage.get("stage") or "")
            if not name:
                continue
            totals = stage_totals.setdefault(
                name, {"input_peaks": 0, "retained_peaks": 0}
            )
            totals["input_peaks"] += int(stage.get("input_peaks") or 0)
            totals["retained_peaks"] += int(stage.get("retained_peaks") or 0)

    stages = []
    for name in (
        "valid_numeric_positive",
        "relative_intensity",
        "precursor_exclusion",
        "upper_mass_bound",
        "precursor_dominance_rescue",
        "top_n",
    ):
        totals = stage_totals.get(name, {"input_peaks": 0, "retained_peaks": 0})
        stages.append(_stage(name, totals["input_peaks"], totals["retained_peaks"]))

    input_total = sum(int(item.get("input_peaks") or 0) for item in audit_list)
    output_total = sum(int(item.get("output_peaks") or 0) for item in audit_list)
    return {
        "schema_version": 2,
        "spectra_audited": len(audit_list),
        "input_peaks": input_total,
        "output_peaks": output_total,
        "retained_pct": (
            round(100.0 * output_total / input_total, 3) if input_total else None
        ),
        "config": audit_list[0].get("config") if audit_list else None,
        "precursor_dominance_rescue_spectra": sum(
            bool(item.get("precursor_dominance_rescue_triggered"))
            for item in audit_list
        ),
        "stages": stages,
    }


def clean_peaks(peaks, precursor_mz, top_n=None, rel_threshold=None,
                precursor_charge=None, precursor_exclusion_da=0.5,
                upper_mass_bound_enabled=True):
    """
    清洗峰列表：过滤噪音、保留前 N 个强度峰、去除 precursor 附近峰
    
    Args:
        peaks: [[mz, intensity], ...]
        precursor_mz: 前体离子 m/z
        top_n: 保留峰数 (默认 TOP_N_PEAKS)
        rel_threshold: 相对强度阈值 (默认 REL_THRESHOLD)
        precursor_charge: 前体电荷；多电荷前体允许合法的单电荷产物高于前体 m/z
    
    Returns:
        按 m/z 排序的清洗后峰列表
    """
    config = PeakProcessingConfig(
        preset="custom" if top_n is not None or rel_threshold is not None else "library_match",
        relative_intensity_threshold=(
            REL_THRESHOLD if rel_threshold is None else rel_threshold
        ),
        precursor_exclusion_da=precursor_exclusion_da,
        top_n=TOP_N_PEAKS if top_n is None else top_n,
        upper_mass_bound_enabled=upper_mass_bound_enabled,
    )
    cleaned, _ = process_peaks(
        peaks,
        precursor_mz,
        config=config,
        precursor_charge=precursor_charge,
    )
    return cleaned


def normalize_peaks(peaks):
    """强度归一化到 0-1"""
    if not peaks:
        return peaks
    max_int = max(p[1] for p in peaks)
    if max_int <= 0:
        return peaks
    return [[p[0], p[1] / max_int] for p in peaks]
