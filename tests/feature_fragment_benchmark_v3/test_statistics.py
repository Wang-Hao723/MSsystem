from __future__ import annotations

from matching.feature_fragment_benchmark_v3.statistics import (
    estimate_conditional_lr,
    evaluate_primary_gate,
    fit_ce_response,
    holm_adjust,
    partial_pooling_gate,
)


def _rows(*, target: bool):
    rows = []
    for index in range(90):
        source = f"source-{index % 3}"
        laboratory = f"lab-{index % 2}"
        present = index % 10 < (8 if target else 2)
        rows.append(
            {
                "identity_id": ("target" if target else "background") + str(index),
                "formula_id": f"C{index + 10}H{index + 20}O2",
                "source_id": source,
                "laboratory_id": laboratory,
                "present": present,
                "positive_eligible": True,
                "negative_eligible": True,
            }
        )
    return rows


def test_conditional_lr_uses_wider_cluster_interval_and_passes_full_gate() -> None:
    target = _rows(target=True)
    background = _rows(target=False)
    stats = estimate_conditional_lr(
        target, background, bootstrap_replicates=300, seed=42
    )

    assert stats["status"] == "estimated"
    assert stats["lr_positive"]["point"] > 1.0
    assert stats["lr_positive"]["lower_95"] <= stats["lr_positive"]["analytic_lower_95"]
    assert stats["lr_positive"]["upper_95"] >= stats["lr_positive"]["analytic_upper_95"]
    assert stats["positive_bootstrap"]["resampling_order"] == [
        "source",
        "laboratory",
        "identity",
    ]
    gate = evaluate_primary_gate(
        stats,
        target_formula_count=90,
        background_formula_count=90,
        target_sources=["source-0", "source-1", "source-2"],
        background_sources=["source-0", "source-1", "source-2"],
        target_laboratories=["lab-0", "lab-1"],
        background_laboratories=["lab-0", "lab-1"],
        holm_adjusted_p=stats["fisher_exact_two_sided_p"],
    )
    assert gate["positive_gate_passed"] is True
    assert gate["negative_gate_passed"] is True
    assert gate["evidence_grade"] == "G4"


def test_holm_is_monotone_and_partial_pooling_cannot_rescue_failed_domain() -> None:
    adjusted = holm_adjust([0.001, 0.02, 0.04])
    assert adjusted == sorted(adjusted)
    domains = [
        {
            "primary_gate_passed": True,
            "activation": "HCD",
            "ce_unit": "NCE",
            "ce_mode": "fixed",
            "polarity": "positive",
            "adduct": "[M+H]+",
            "charge": 1,
            "shared_identity_count": 25,
            "shared_ce_setting_count": 4,
            "source_direction_consistent": True,
            "leave_one_domain_out_lower_95": 1.2,
            "lr_positive_point": 3.0,
        },
        {
            "primary_gate_passed": False,
            "activation": "HCD",
            "ce_unit": "NCE",
            "ce_mode": "fixed",
            "polarity": "positive",
            "adduct": "[M+H]+",
            "charge": 1,
            "shared_identity_count": 25,
            "shared_ce_setting_count": 4,
            "source_direction_consistent": True,
            "leave_one_domain_out_lower_95": 1.2,
            "lr_positive_point": 3.2,
        },
    ]
    result = partial_pooling_gate(domains)
    assert result["allowed"] is False
    assert result["may_rescue_failed_domain"] is False
    assert "one_or_more_exact_domains_failed_primary_gate" in result["failures"]


def test_ce_models_use_categorical_three_settings_and_spline_for_four() -> None:
    categorical = []
    spline = []
    for identity in range(25):
        for setting_index, ce in enumerate((10.0, 20.0, 30.0, 40.0)):
            row = {
                "identity_id": f"identity-{identity}",
                "ce_value": ce,
                "present": (identity + setting_index) % 5 < setting_index + 1,
                "precursor_mz": 250.0 + identity,
                "charge": 1 + identity % 2,
            }
            spline.append(row)
            if ce != 40.0:
                categorical.append(row)
    categorical_result = fit_ce_response(categorical, ce_unit="eV")
    spline_result = fit_ce_response(spline, ce_unit="NCE")

    assert categorical_result["model_type"] == "categorical_fixed_ce"
    assert categorical_result["status"] == "fitted"
    assert spline_result["model_type"] == "binomial_natural_cubic_spline_df3"
    assert spline_result["status"] == "fitted"
    assert spline_result["nce_covariates"] == ["precursor_mz", "charge"]
    assert spline_result["cross_unit_conversion_applied"] is False
