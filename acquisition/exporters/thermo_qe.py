"""Thermo QE-oriented advisory profile (CSV, not a method file)."""

from __future__ import annotations

from collections.abc import Mapping

from .vendor_neutral import export_plan_csv


def export_thermo_qe_advisory(plan: Mapping) -> str:
    return export_plan_csv(plan, profile="thermo-qe-advisory")
