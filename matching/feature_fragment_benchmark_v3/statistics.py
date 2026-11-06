"""Conditional likelihood ratios, CE response models, and evidence gates."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from scipy.stats import fisher_exact


Z_975 = 1.959963984540054


def _katz_haldane(
    target_events: int,
    target_total: int,
    background_events: int,
    background_total: int,
) -> dict[str, Any]:
    """Haldane-Anscombe corrected log risk-ratio interval."""
    a = target_events + 0.5
    nt = target_total + 1.0
    c = background_events + 0.5
    nb = background_total + 1.0
    risk_target = a / nt
    risk_background = c / nb
    point = risk_target / risk_background
    variance = max(0.0, 1.0 / a - 1.0 / nt + 1.0 / c - 1.0 / nb)
    half = Z_975 * math.sqrt(variance)
    log_point = math.log(point)
    return {
        "point": point,
        "log_point": log_point,
        "lower_95": math.exp(log_point - half),
        "upper_95": math.exp(log_point + half),
        "log_half_width_95": half,
        "method": "Katz log risk-ratio with Haldane-Anscombe 0.5 correction",
    }


def _counts(rows: Sequence[Mapping[str, Any]]) -> tuple[int, int]:
    return sum(bool(row.get("present")) for row in rows), len(rows)


def _hierarchical_resample(
    rows: Sequence[Mapping[str, Any]], rng: np.random.Generator
) -> list[Mapping[str, Any]]:
    by_source: dict[str, dict[str, list[Mapping[str, Any]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in rows:
        by_source[str(row["source_id"])][str(row["laboratory_id"])].append(row)
    sources = sorted(by_source)
    sampled: list[Mapping[str, Any]] = []
    for source_index in rng.integers(0, len(sources), size=len(sources)):
        source = sources[int(source_index)]
        labs = sorted(by_source[source])
        for lab_index in rng.integers(0, len(labs), size=len(labs)):
            lab = labs[int(lab_index)]
            identities = by_source[source][lab]
            for identity_index in rng.integers(
                0, len(identities), size=len(identities)
            ):
                sampled.append(identities[int(identity_index)])
    return sampled


def _bootstrap_lr(
    target: Sequence[Mapping[str, Any]],
    background: Sequence[Mapping[str, Any]],
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    source_count = len(
        {str(row["source_id"]) for row in list(target) + list(background)}
    )
    laboratory_count = len(
        {str(row["laboratory_id"]) for row in list(target) + list(background)}
    )
    if source_count < 3 or laboratory_count < 2 or not target or not background:
        return {
            "status": "insufficient_source_or_laboratory_clusters",
            "replicates_requested": replicates,
            "replicates_used": 0,
            "source_count": source_count,
            "laboratory_count": laboratory_count,
            "positive": None,
            "negative": None,
        }
    rng = np.random.default_rng(seed)
    positive_values: list[float] = []
    negative_values: list[float] = []
    for _ in range(replicates):
        bt = _hierarchical_resample(target, rng)
        bb = _hierarchical_resample(background, rng)
        if not bt or not bb:
            continue
        tp, nt = _counts(bt)
        bp, nb = _counts(bb)
        positive_values.append(_katz_haldane(tp, nt, bp, nb)["point"])
        negative_values.append(
            _katz_haldane(nt - tp, nt, nb - bp, nb)["point"]
        )
    if not positive_values:
        return {
            "status": "no_valid_bootstrap_replicates",
            "replicates_requested": replicates,
            "replicates_used": 0,
            "source_count": source_count,
            "laboratory_count": laboratory_count,
            "positive": None,
            "negative": None,
        }

    def interval(values: Sequence[float]) -> dict[str, float]:
        lower, upper = np.percentile(np.asarray(values), [2.5, 97.5])
        return {"lower_95": float(lower), "upper_95": float(upper)}

    return {
        "status": "estimated",
        "replicates_requested": replicates,
        "replicates_used": len(positive_values),
        "source_count": source_count,
        "laboratory_count": laboratory_count,
        "positive": interval(positive_values),
        "negative": interval(negative_values),
        "resampling_order": ["source", "laboratory", "identity"],
    }


def _wider_interval(
    analytic: Mapping[str, Any], bootstrap: Mapping[str, Any] | None
) -> dict[str, Any]:
    if not bootstrap:
        return {**analytic, "reported_interval_basis": "analytic_only_bootstrap_unavailable"}
    lower = min(float(analytic["lower_95"]), float(bootstrap["lower_95"]))
    upper = max(float(analytic["upper_95"]), float(bootstrap["upper_95"]))
    point = float(analytic["point"])
    half = max(math.log(point) - math.log(lower), math.log(upper) - math.log(point))
    return {
        **analytic,
        "analytic_lower_95": analytic["lower_95"],
        "analytic_upper_95": analytic["upper_95"],
        "bootstrap_lower_95": bootstrap["lower_95"],
        "bootstrap_upper_95": bootstrap["upper_95"],
        "lower_95": lower,
        "upper_95": upper,
        "log_half_width_95": half,
        "reported_interval_basis": "union_of_analytic_and_hierarchical_bootstrap_intervals",
    }


def estimate_conditional_lr(
    target: Sequence[Mapping[str, Any]],
    background: Sequence[Mapping[str, Any]],
    *,
    bootstrap_replicates: int = 2000,
    seed: int = 151503,
) -> dict[str, Any]:
    """Estimate LR+ and LR- from one exact acquisition domain.

    LR+ uses records marked ``positive_eligible`` (or all rows when that field
    is absent for a prevalidated caller).  LR- has an independent denominator
    and uses only rows marked ``negative_eligible``.  One input row is one
    identity vote.
    """
    positive_target = [
        row for row in target if row.get("positive_eligible", True) is True
    ]
    positive_background = [
        row for row in background if row.get("positive_eligible", True) is True
    ]
    tp, nt = _counts(positive_target)
    bp, nb = _counts(positive_background)
    if not nt or not nb:
        return {
            "status": "unavailable",
            "reason": "empty_positive_evidence_denominator",
            "lr_positive": None,
            "lr_negative": None,
        }
    analytic_positive = _katz_haldane(tp, nt, bp, nb)
    fisher_p = float(fisher_exact([[tp, nt - tp], [bp, nb - bp]])[1])
    positive_bootstrap = _bootstrap_lr(
        positive_target,
        positive_background,
        replicates=bootstrap_replicates,
        seed=seed,
    )
    lr_positive = _wider_interval(
        analytic_positive, positive_bootstrap.get("positive")
    )

    negative_target = [
        row for row in target if row.get("negative_eligible", False) is True
    ]
    negative_background = [
        row for row in background if row.get("negative_eligible", False) is True
    ]
    ntp, ntt = _counts(negative_target)
    nbp, nbt = _counts(negative_background)
    if ntt and nbt:
        analytic_negative = _katz_haldane(ntt - ntp, ntt, nbt - nbp, nbt)
        negative_bootstrap = _bootstrap_lr(
            negative_target,
            negative_background,
            replicates=bootstrap_replicates,
            seed=seed + 1,
        )
        lr_negative = _wider_interval(
            analytic_negative, negative_bootstrap.get("negative")
        )
    else:
        negative_bootstrap = {
            "status": "unavailable",
            "reason": "empty_verified_isolation_denominator",
        }
        lr_negative = None
    return {
        "status": "estimated",
        "target_present": tp,
        "target_total": nt,
        "background_present": bp,
        "background_total": nb,
        "target_frequency": tp / nt,
        "background_frequency": bp / nb,
        "target_absence_total_verified": ntt,
        "background_absence_total_verified": nbt,
        "fisher_exact_two_sided_p": fisher_p,
        "lr_positive": lr_positive,
        "lr_negative": lr_negative,
        "positive_bootstrap": positive_bootstrap,
        "negative_bootstrap": negative_bootstrap,
        "identity_vote_policy": "one_two_dimensional_connectivity_identity_one_vote",
    }


def benjamini_hochberg(p_values: Sequence[float | None]) -> list[float | None]:
    """Return monotone Benjamini-Hochberg adjusted p-values."""
    valid = [(index, float(value)) for index, value in enumerate(p_values) if value is not None]
    adjusted: list[float | None] = [None] * len(p_values)
    running = 1.0
    for rank_from_end, (index, value) in enumerate(
        sorted(valid, key=lambda item: item[1], reverse=True), start=1
    ):
        rank = len(valid) - rank_from_end + 1
        running = min(running, value * len(valid) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def holm_adjust(p_values: Sequence[float | None]) -> list[float | None]:
    """Return step-down Holm family-wise adjusted p-values."""
    valid = sorted(
        ((index, float(value)) for index, value in enumerate(p_values) if value is not None),
        key=lambda item: item[1],
    )
    adjusted: list[float | None] = [None] * len(p_values)
    running = 0.0
    total = len(valid)
    for rank, (index, value) in enumerate(valid):
        running = max(running, (total - rank) * value)
        adjusted[index] = min(1.0, running)
    return adjusted


def evaluate_primary_gate(
    statistics: Mapping[str, Any],
    *,
    target_formula_count: int,
    background_formula_count: int,
    target_sources: Sequence[str],
    background_sources: Sequence[str],
    target_laboratories: Sequence[str],
    background_laboratories: Sequence[str],
    holm_adjusted_p: float | None,
    require_negative_evidence: bool = False,
) -> dict[str, Any]:
    """Apply every preregistered V2 threshold plus V3 isolation gates."""
    failures: list[str] = []
    if statistics.get("status") != "estimated":
        return {
            "passed": False,
            "failures": ["likelihood_ratio_unavailable"],
            "evidence_grade": "G0",
        }
    if int(statistics["target_total"]) < 30:
        failures.append("test_target_identities_lt_30")
    if int(statistics["background_total"]) < 30:
        failures.append("test_background_identities_lt_30")
    if target_formula_count < 10:
        failures.append("test_target_formulas_lt_10")
    if background_formula_count < 10:
        failures.append("test_background_formulas_lt_10")
    if len(set(target_sources)) < 3 or len(set(background_sources)) < 3:
        failures.append("test_sources_lt_3")
    if len(set(target_laboratories)) < 2 or len(set(background_laboratories)) < 2:
        failures.append("test_laboratories_lt_2")
    lr_positive = statistics.get("lr_positive")
    if not lr_positive or float(lr_positive["lower_95"]) <= 1.0:
        failures.append("lr_positive_interval_not_above_1")
    if not lr_positive or float(lr_positive["log_half_width_95"]) > math.log(2.0):
        failures.append("lr_positive_interval_too_wide")
    if holm_adjusted_p is None or holm_adjusted_p > 0.05:
        failures.append("holm_family_wise_p_gt_0_05")
    lr_negative = statistics.get("lr_negative")
    negative_passed = bool(
        lr_negative
        and float(lr_negative["upper_95"]) < 1.0
        and float(lr_negative["log_half_width_95"]) <= math.log(2.0)
    )
    if require_negative_evidence and not negative_passed:
        failures.append("lr_negative_gate_failed_or_unavailable")
    positive_passed = not [
        reason for reason in failures if reason != "lr_negative_gate_failed_or_unavailable"
    ]
    grade = "G4" if positive_passed and negative_passed else "G3" if positive_passed else "G2"
    return {
        "passed": not failures,
        "positive_gate_passed": positive_passed,
        "negative_gate_passed": negative_passed,
        "failures": failures,
        "evidence_grade": grade,
    }


def partial_pooling_gate(domains: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Check whether already-valid exact domains permit secondary pooling."""
    failures = []
    if len(domains) < 2:
        failures.append("fewer_than_two_domains")
    if any(not row.get("primary_gate_passed") for row in domains):
        failures.append("one_or_more_exact_domains_failed_primary_gate")
    for field in ("activation", "ce_unit", "ce_mode", "polarity", "adduct", "charge"):
        if len({row.get(field) for row in domains}) > 1:
            failures.append(f"mixed_{field}")
    if min((int(row.get("shared_identity_count", 0)) for row in domains), default=0) < 20:
        failures.append("shared_identities_lt_20")
    if min((int(row.get("shared_ce_setting_count", 0)) for row in domains), default=0) < 3:
        failures.append("shared_ce_settings_lt_3")
    if any(not row.get("source_direction_consistent") for row in domains):
        failures.append("source_direction_inconsistent")
    if min((float(row.get("leave_one_domain_out_lower_95", 0.0)) for row in domains), default=0.0) <= 1.0:
        failures.append("leave_one_domain_out_interval_not_above_1")
    log_points = [math.log(float(row["lr_positive_point"])) for row in domains if float(row.get("lr_positive_point", 0.0)) > 0]
    tau = float(np.std(log_points, ddof=1)) if len(log_points) >= 2 else None
    if tau is None or tau >= math.log(1.5):
        failures.append("between_domain_log_lr_sd_too_large")
    return {
        "allowed": not failures,
        "may_rescue_failed_domain": False,
        "between_domain_log_lr_sd": tau,
        "failures": sorted(set(failures)),
    }


def _shared_identity_rows(
    observations: Sequence[Mapping[str, Any]], settings: Sequence[float]
) -> tuple[set[str], list[Mapping[str, Any]]]:
    by_identity: dict[str, set[float]] = defaultdict(set)
    for row in observations:
        by_identity[str(row["identity_id"])].add(float(row["ce_value"]))
    required = set(settings)
    shared = {identity for identity, seen in by_identity.items() if required <= seen}
    return shared, [row for row in observations if str(row["identity_id"]) in shared]


def fit_ce_response(
    observations: Sequence[Mapping[str, Any]],
    *,
    ce_unit: str,
    minimum_shared_identities: int = 20,
) -> dict[str, Any]:
    """Fit the preregistered fixed-CE model without cross-unit conversion."""
    if not observations:
        return {"status": "unavailable", "reason": "no_observations"}
    settings = sorted({float(row["ce_value"]) for row in observations})
    shared, rows = _shared_identity_rows(observations, settings)
    rates = []
    for setting in settings:
        at_setting = [row for row in rows if float(row["ce_value"]) == setting]
        rates.append(
            {
                "ce_value": setting,
                "present": sum(bool(row.get("present")) for row in at_setting),
                "total": len(at_setting),
                "frequency": (
                    sum(bool(row.get("present")) for row in at_setting) / len(at_setting)
                    if at_setting
                    else None
                ),
            }
        )
    base = {
        "ce_unit": ce_unit,
        "cross_unit_conversion_applied": False,
        "ce_setting_count": len(settings),
        "shared_identity_count": len(shared),
        "minimum_shared_identities": minimum_shared_identities,
        "rates": rates,
    }
    if len(shared) < minimum_shared_identities:
        bins: dict[float, list[Mapping[str, Any]]] = defaultdict(list)
        for row in observations:
            bins[5.0 * math.floor(float(row["ce_value"]) / 5.0)].append(row)
        return {
            **base,
            "status": "insufficient_shared_identities",
            "model_type": "five_native_unit_bins",
            "bins": [
                {
                    "minimum": lower,
                    "maximum_exclusive": lower + 5.0,
                    "identity_count": len({str(row["identity_id"]) for row in bucket}),
                    "requires_all_primary_gates": True,
                }
                for lower, bucket in sorted(bins.items())
            ],
        }
    if len(settings) == 3:
        return {
            **base,
            "status": "fitted",
            "model_type": "categorical_fixed_ce",
            "parameters": {
                str(row["ce_value"]): row["frequency"] for row in rates
            },
        }
    if len(settings) < 4:
        return {
            **base,
            "status": "insufficient_ce_settings",
            "model_type": "five_native_unit_bins",
        }

    try:
        from patsy import dmatrix
        from scipy.optimize import minimize
        from scipy.special import expit

        ce_values = np.asarray([float(row["ce_value"]) for row in rows])
        response = np.asarray([1.0 if row.get("present") else 0.0 for row in rows])
        design = dmatrix(
            "cr(ce, df=3)", {"ce": ce_values}, return_type="dataframe"
        )
        extra_names = []
        extra_columns = []
        if ce_unit == "NCE":
            precursor = np.asarray([float(row["precursor_mz"]) for row in rows])
            charge = np.asarray([float(row["charge"]) for row in rows])
            precursor_sd = float(precursor.std()) or 1.0
            charge_sd = float(charge.std()) or 1.0
            extra_columns = [
                (precursor - precursor.mean()) / precursor_sd,
                (charge - charge.mean()) / charge_sd,
            ]
            extra_names = ["precursor_mz_standardized", "charge_standardized"]
        matrix = np.column_stack([design.to_numpy(), *extra_columns])
        def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
            linear = matrix @ parameters
            probabilities = expit(linear)
            epsilon = 1e-12
            loss = -float(
                np.sum(
                    response * np.log(probabilities + epsilon)
                    + (1.0 - response) * np.log(1.0 - probabilities + epsilon)
                )
            )
            ridge = 1e-8 * float(parameters @ parameters)
            gradient = matrix.T @ (probabilities - response) + 2e-8 * parameters
            return loss + ridge, gradient

        fitted = minimize(
            lambda parameters: objective(parameters)[0],
            np.zeros(matrix.shape[1], dtype=float),
            jac=lambda parameters: objective(parameters)[1],
            method="BFGS",
        )
        names = list(design.design_info.column_names) + extra_names
        return {
            **base,
            "status": "fitted",
            "model_type": "binomial_natural_cubic_spline_df3",
            "parameter_names": names,
            "parameters": [float(value) for value in fitted.x],
            "converged": bool(fitted.success),
            "optimizer": "SciPy BFGS with 1e-8 L2 stabilization",
            "nce_covariates": ["precursor_mz", "charge"] if ce_unit == "NCE" else [],
        }
    except Exception as exc:
        return {
            **base,
            "status": "model_fit_failed",
            "model_type": "binomial_natural_cubic_spline_df3",
            "error_type": type(exc).__name__,
        }


def residual_phi(left: Sequence[bool], right: Sequence[bool]) -> float | None:
    """Compute phi for two frozen binary fragment pathways."""
    if len(left) != len(right) or not left:
        return None
    counts = Counter(zip(map(bool, left), map(bool, right)))
    a = counts[(True, True)]
    b = counts[(True, False)]
    c = counts[(False, True)]
    d = counts[(False, False)]
    denominator = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    return (a * d - b * c) / denominator if denominator else None
