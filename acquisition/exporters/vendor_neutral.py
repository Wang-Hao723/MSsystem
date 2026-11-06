"""Vendor-neutral advisory CSV export."""

from __future__ import annotations

import csv
import io
from collections.abc import Mapping


def _energy_stages(design: Mapping) -> str:
    values = []
    for stage in design.get("stages") or []:
        energies = ";".join(
            str(value) for value in stage.get("collision_energies") or []
        )
        values.append(f"{stage.get('stage_id')}={energies or 'method-owner-defined'}")
    return " | ".join(values)


def _profile_compatible(profile: str, design: Mapping) -> bool:
    design_profile = str(design.get("profile_id") or "")
    if profile == "thermo-qe-advisory" and design_profile.startswith("agilent_qtof"):
        return False
    return True


def export_plan_csv(plan: Mapping, *, profile: str = "vendor-neutral") -> str:
    output = io.StringIO()
    output.write(
        "# EviAcquire-MS advisory list; manual method review required; "
        "not an instrument method file\n"
    )
    output.write(f"# profile={profile}; module_state=shadow; plan_id={plan.get('plan_id')}\n")
    writer = csv.writer(output)
    writer.writerow([
        "priority", "action_kind", "precursor_mz", "retention_time_min",
        "rt_window_min", "polarity", "collision_energies",
        "collision_energy_unit", "isolation_width_da",
        "isolation_width_profile_id", "isolation_width_evidence_status",
        "max_injection_time_ms", "diagnostic_ions", "utility",
        "utility_low", "utility_high", "feasibility", "energy_profile_id",
        "energy_evidence_status", "energy_stages", "energy_response_metrics",
        "automatic_energy_ceiling", "universal_complete_fragmentation_assumed",
        "cross_instrument_conversion_applied", "method_owner_review_required",
        "vendor_profile_compatible", "advisory_only", "required_resources",
        "resource_status", "resource_identity_keys", "resource_declaration_id",
    ])
    for index, action in enumerate(plan.get("ranked_actions") or [], start=1):
        interval = action.get("utility_interval") or [None, None]
        design = action.get("collision_energy_design") or {}
        writer.writerow([
            index,
            action.get("action_kind"),
            action.get("precursor_mz"),
            action.get("retention_time_min"),
            action.get("rt_window_min"),
            action.get("polarity"),
            ";".join(str(value) for value in action.get("collision_energies") or []),
            action.get("collision_energy_unit"),
            action.get("isolation_width_da"),
            action.get("isolation_width_profile_id"),
            action.get("isolation_width_evidence_status"),
            action.get("max_injection_time_ms"),
            ";".join(str(value) for value in action.get("diagnostic_ions") or []),
            action.get("utility"),
            interval[0], interval[1],
            action.get("feasibility"),
            design.get("profile_id"),
            design.get("evidence_status"),
            _energy_stages(design),
            ";".join(design.get("response_metrics") or []),
            design.get("automatic_ceiling"),
            design.get("universal_complete_fragmentation_assumed", False),
            design.get("cross_instrument_conversion_applied", False),
            design.get("method_owner_review_required", bool(design)),
            _profile_compatible(profile, design),
            True,
            ";".join(action.get("required_resources") or []),
            action.get("resource_status", "not_required"),
            ";".join(action.get("resource_identity_keys") or []),
            action.get("resource_declaration_id"),
        ])
    return output.getvalue()
