"""Evidence-gated deconvolution of chimeric MS/MS spectra.

The implementation follows the central statistical idea used by DecoID: an
observed spectrum is represented as a non-negative sparse combination of
candidate component spectra.  A positive L1 penalty is deliberately mandatory;
plain NNLS is not used because it tends to assign small coefficients to many
irrelevant library spectra.

This module is intentionally independent from the library search engine.  It
only solves and audits a supplied component set.  Candidate discovery and any
confidence decision remain the responsibility of the caller so that the
evidence trail stays explicit.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


Peak = Sequence[float]


@dataclass(frozen=True)
class _PreparedSystem:
    axis_mz: np.ndarray
    query: np.ndarray
    matrix: np.ndarray
    component_ids: tuple[str, ...]


def _finite_positive_peaks(peaks: Iterable[Peak] | None) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for peak in peaks or []:
        if not isinstance(peak, (list, tuple, np.ndarray)) or len(peak) < 2:
            continue
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError):
            continue
        if np.isfinite(mz) and np.isfinite(intensity) and mz > 0 and intensity > 0:
            out.append((mz, intensity))
    return out


def _mass_tolerance(mz: float, tolerance_da: float, tolerance_ppm: float) -> float:
    return max(float(tolerance_da), abs(float(mz)) * float(tolerance_ppm) * 1e-6)


def _build_axis(
    spectra: Sequence[list[tuple[float, float]]],
    *,
    tolerance_da: float,
    tolerance_ppm: float,
) -> np.ndarray:
    """Create a shared, deterministic mass axis without fixed-bin edge effects."""
    masses = sorted(mz for spectrum in spectra for mz, _ in spectrum)
    if not masses:
        return np.empty(0, dtype=float)
    clusters: list[list[float]] = [[masses[0]]]
    for mz in masses[1:]:
        centre = float(np.mean(clusters[-1]))
        if abs(mz - centre) <= _mass_tolerance(centre, tolerance_da, tolerance_ppm):
            clusters[-1].append(mz)
        else:
            clusters.append([mz])
    return np.asarray([float(np.mean(cluster)) for cluster in clusters], dtype=float)


def _vectorize(
    peaks: list[tuple[float, float]],
    axis: np.ndarray,
    *,
    tolerance_da: float,
    tolerance_ppm: float,
) -> np.ndarray:
    vector = np.zeros(axis.size, dtype=float)
    if not peaks or axis.size == 0:
        return vector
    # Spectral matching is scale invariant.  Scale before bin accumulation so
    # malformed/vendor-normalised spectra near float64 limits cannot overflow
    # when several peaks fall into one tolerance cluster.
    intensity_scale = max(intensity for _, intensity in peaks)
    if not np.isfinite(intensity_scale) or intensity_scale <= 0:
        return vector
    for mz, intensity in peaks:
        pos = int(np.searchsorted(axis, mz))
        candidates = [idx for idx in (pos - 1, pos) if 0 <= idx < axis.size]
        if not candidates:
            continue
        idx = min(candidates, key=lambda item: abs(float(axis[item]) - mz))
        if abs(float(axis[idx]) - mz) <= _mass_tolerance(
            float(axis[idx]), tolerance_da, tolerance_ppm
        ):
            vector[idx] += intensity / intensity_scale
    return vector


def _safe_unit_vector(vector: np.ndarray) -> np.ndarray | None:
    """Return an L2-unit vector without overflow/underflow propagation."""
    values = np.asarray(vector, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        return None
    scale = float(np.max(np.abs(values), initial=0.0))
    if not np.isfinite(scale) or scale <= 0:
        return None
    scaled = values / scale
    norm = float(np.linalg.norm(scaled))
    if not np.isfinite(norm) or norm <= np.finfo(float).eps:
        return None
    unit = scaled / norm
    return unit if np.all(np.isfinite(unit)) else None


def _prepare_system(
    query_peaks: Iterable[Peak],
    components: Sequence[Mapping],
    *,
    tolerance_da: float,
    tolerance_ppm: float,
) -> _PreparedSystem | None:
    query = _finite_positive_peaks(query_peaks)
    prepared_components: list[tuple[str, list[tuple[float, float]]]] = []
    seen: set[str] = set()
    for index, component in enumerate(components):
        identifier = str(component.get("id") or f"component_{index}")
        if identifier in seen:
            continue
        peaks = _finite_positive_peaks(component.get("peaks"))
        if len(peaks) < 2:
            continue
        seen.add(identifier)
        prepared_components.append((identifier, peaks))
    if len(query) < 3 or len(prepared_components) < 2:
        return None

    axis = _build_axis(
        [query, *[item[1] for item in prepared_components]],
        tolerance_da=tolerance_da,
        tolerance_ppm=tolerance_ppm,
    )
    y = _vectorize(
        query, axis, tolerance_da=tolerance_da, tolerance_ppm=tolerance_ppm
    )
    y = _safe_unit_vector(y)
    if y is None:
        return None

    columns = []
    ids = []
    for identifier, peaks in prepared_components:
        column = _vectorize(
            peaks, axis, tolerance_da=tolerance_da, tolerance_ppm=tolerance_ppm
        )
        column = _safe_unit_vector(column)
        if column is None:
            continue
        columns.append(column)
        ids.append(identifier)
    if len(columns) < 2:
        return None
    return _PreparedSystem(
        axis_mz=axis,
        query=y,
        matrix=np.column_stack(columns),
        component_ids=tuple(ids),
    )


def nonnegative_lasso(
    matrix: np.ndarray,
    observation: np.ndarray,
    *,
    penalty: float,
    max_iterations: int = 1000,
    tolerance: float = 1e-8,
) -> np.ndarray:
    """Solve ``||Xb-y||² + 2*penalty*||b||₁`` with ``b >= 0``.

    Cyclic coordinate descent is deterministic and has no optional solver
    dependency.  ``penalty`` must be strictly positive by design.
    """
    x = np.asarray(matrix, dtype=float)
    y = np.asarray(observation, dtype=float)
    if x.ndim != 2 or y.ndim != 1 or x.shape[0] != y.shape[0]:
        raise ValueError("matrix/observation dimensions are incompatible")
    if not np.isfinite(penalty) or penalty <= 0:
        raise ValueError("a strictly positive L1 penalty is required")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("matrix/observation must contain only finite values")
    if np.any(x < 0) or np.any(y < 0):
        raise ValueError("non-negative matrix/observation values are required")
    beta = np.zeros(x.shape[1], dtype=float)
    residual = y.copy()
    norms = np.sum(x * x, axis=0)
    if not np.all(np.isfinite(norms)):
        raise ValueError("component norms are numerically unstable")
    for _ in range(max(1, int(max_iterations))):
        largest_change = 0.0
        for column_index in range(x.shape[1]):
            norm = float(norms[column_index])
            if norm <= 0:
                continue
            old = float(beta[column_index])
            residual += x[:, column_index] * old
            rho = float(np.dot(x[:, column_index], residual))
            new = max(0.0, (rho - float(penalty)) / norm)
            if not np.isfinite(new):
                raise ValueError("coordinate descent produced a non-finite coefficient")
            beta[column_index] = new
            residual -= x[:, column_index] * new
            largest_change = max(largest_change, abs(new - old))
        if largest_change <= tolerance:
            break
    return beta


def _component_fractions(matrix: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    contributions = np.sum(matrix * coefficients[np.newaxis, :], axis=0)
    total = float(np.sum(contributions))
    if total <= 0:
        return np.zeros_like(coefficients)
    return contributions / total


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator > 0 else 0.0


def deconvolve_chimeric_spectrum(
    query_peaks: Iterable[Peak],
    components: Sequence[Mapping],
    target_id: str,
    *,
    penalty: float = 0.02,
    tolerance_da: float = 0.01,
    tolerance_ppm: float = 20.0,
    min_reconstruction_cosine: float = 0.75,
    min_explained_energy: float = 0.40,
    min_target_fraction: float = 0.10,
    max_target_fraction_span: float = 0.20,
    max_component_similarity: float = 0.97,
) -> dict:
    """Deconvolve and return an auditable target spectrum.

    ``purified_target_peaks`` contains only intensity actually observed in the
    query.  At each mass bin it allocates observed intensity in proportion to
    the fitted target contribution; it never returns the library column itself.
    This prevents circularly manufacturing a perfect library match.
    """
    base = {
        "method": "nonnegative_l1_sparse_regression",
        "runtime_influence": False,
        "score_contribution": 0.0,
        "target_id": str(target_id),
        "penalty": float(penalty),
    }
    prepared = _prepare_system(
        query_peaks,
        components,
        tolerance_da=tolerance_da,
        tolerance_ppm=tolerance_ppm,
    )
    if prepared is None:
        return {**base, "status": "insufficient_input", "accepted": False}
    if str(target_id) not in prepared.component_ids:
        return {**base, "status": "target_component_missing", "accepted": False}

    target_index = prepared.component_ids.index(str(target_id))
    coefficients = nonnegative_lasso(
        prepared.matrix, prepared.query, penalty=penalty
    )
    if not np.all(np.isfinite(coefficients)):
        return {**base, "status": "numerical_instability", "accepted": False}
    # These systems have only a handful of columns.  Explicit non-negative
    # contribution summation is both transparent and avoids false floating-point
    # flags observed in Apple's BLAS/MPS process after neural-classifier calls.
    contribution_matrix = prepared.matrix * coefficients[np.newaxis, :]
    prediction = np.sum(contribution_matrix, axis=1)
    if not np.all(np.isfinite(prediction)):
        return {**base, "status": "numerical_instability", "accepted": False}
    fractions = _component_fractions(prepared.matrix, coefficients)
    active_floor = max(1e-6, float(np.max(coefficients)) * 0.02)
    active = np.flatnonzero(coefficients >= active_floor)
    reconstruction_cosine = _cosine(prepared.query, prediction)
    residual_energy = float(np.sum((prepared.query - prediction) ** 2))
    explained_energy = min(max(1.0 - residual_energy, 0.0), 1.0)

    # Components that are effectively identical cannot be uniquely separated.
    target_column = prepared.matrix[:, target_index]
    other_similarities = [
        _cosine(target_column, prepared.matrix[:, idx])
        for idx in range(prepared.matrix.shape[1])
        if idx != target_index
    ]
    maximum_similarity = max(other_similarities, default=0.0)

    # Stability under a modest penalty perturbation is an explicit guard
    # against solutions that exist only at one arbitrary regularisation value.
    stability_fractions = []
    target_active_all = True
    for factor in (0.75, 1.25):
        trial = nonnegative_lasso(
            prepared.matrix, prepared.query, penalty=penalty * factor
        )
        trial_fractions = _component_fractions(prepared.matrix, trial)
        stability_fractions.append(float(trial_fractions[target_index]))
        target_active_all = target_active_all and bool(trial[target_index] > 1e-6)
    target_fraction_span = max(
        [float(fractions[target_index]), *stability_fractions]
    ) - min([float(fractions[target_index]), *stability_fractions])

    total_prediction = np.sum(contribution_matrix, axis=1)
    target_prediction = contribution_matrix[:, target_index]
    allocation = np.divide(
        target_prediction,
        total_prediction,
        out=np.zeros_like(target_prediction),
        where=total_prediction > 1e-12,
    )
    purified = prepared.query * np.clip(allocation, 0.0, 1.0)
    if float(np.max(purified, initial=0.0)) > 0:
        purified /= float(np.max(purified))
    purified_peaks = [
        [round(float(prepared.axis_mz[idx]), 6), round(float(value), 8)]
        for idx, value in enumerate(purified)
        if value >= 0.005
    ]

    component_rows = []
    for idx, identifier in enumerate(prepared.component_ids):
        component_rows.append({
            "id": identifier,
            "coefficient": round(float(coefficients[idx]), 8),
            "fraction": round(float(fractions[idx]), 6),
            "active": bool(idx in active),
        })

    # Every pre-score gate is serialized independently.  The caller adds the
    # final spectral-gain gate after rescoring the purified *observed* peaks.
    # This makes "all gates passed" a machine-auditable invariant rather than
    # an inference from a single status string.
    quality_gates = {
        "multi_component_solution": bool(len(active) >= 2),
        "target_component_supported": bool(target_index in active),
        "reconstruction_quality": bool(
            reconstruction_cosine >= min_reconstruction_cosine
            and explained_energy >= min_explained_energy
        ),
        "regularization_stability": bool(
            target_active_all and target_fraction_span <= max_target_fraction_span
        ),
        "target_fraction": bool(
            float(fractions[target_index]) >= min_target_fraction
        ),
        "non_redundancy": bool(maximum_similarity < max_component_similarity),
        "sufficient_purified_peaks": bool(len(purified_peaks) >= 3),
    }
    status = "accepted"
    if not quality_gates["multi_component_solution"]:
        status = "single_component_solution"
    elif not quality_gates["target_component_supported"]:
        status = "target_not_supported"
    elif not quality_gates["non_redundancy"]:
        status = "non_unique_components"
    elif reconstruction_cosine < min_reconstruction_cosine:
        status = "poor_reconstruction"
    elif explained_energy < min_explained_energy:
        status = "low_explained_energy"
    elif not quality_gates["target_fraction"]:
        status = "target_fraction_too_small"
    elif not quality_gates["regularization_stability"]:
        status = "regularisation_unstable"
    elif not quality_gates["sufficient_purified_peaks"]:
        status = "insufficient_purified_peaks"

    accepted = status == "accepted"
    return {
        **base,
        "status": status,
        "accepted": accepted,
        "runtime_influence": accepted,
        "reconstruction_cosine": round(reconstruction_cosine, 6),
        "explained_energy": round(explained_energy, 6),
        "target_fraction": round(float(fractions[target_index]), 6),
        "target_fraction_stability": [round(value, 6) for value in stability_fractions],
        "target_fraction_span": round(target_fraction_span, 6),
        "active_component_count": int(len(active)),
        "component_count": int(prepared.matrix.shape[1]),
        "target_max_component_similarity": round(maximum_similarity, 6),
        "components": component_rows,
        "purified_target_peaks": purified_peaks,
        "quality_gates": quality_gates,
        "all_pre_score_gates_passed": bool(all(quality_gates.values())),
        "gate_thresholds": {
            "min_reconstruction_cosine": float(min_reconstruction_cosine),
            "min_explained_energy": float(min_explained_energy),
            "max_target_fraction_span": float(max_target_fraction_span),
            "min_target_fraction": float(min_target_fraction),
            "max_component_similarity": float(max_component_similarity),
            "min_purified_peaks": 3,
        },
        "residual_fraction": round(max(0.0, 1.0 - explained_energy), 6),
        "coefficient_interpretation": (
            "Sparse reconstruction coefficients are not quantitative abundance estimates."
        ),
    }
