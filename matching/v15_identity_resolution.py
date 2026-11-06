"""Conservative identity-scope resolution for the V15 ranked shortlist.

The V15 ranker is allowed to order candidates, but an ordering alone does not
establish that the first structure is uniquely identified.  This module keeps
the ranking byte-for-byte stable and determines the narrowest identity scope
that the competitive shortlist can support:

* one competitive connectivity may retain the upstream decision;
* multiple competitive connectivities with one molecular formula form an
  unresolved molecular-formula isomer group; and
* competing connectivities with conflicting or incomplete formula assignments
  form an unresolved precursor event and cannot be counted as a compound.

Candidate names, sample names, ground truth and biological annotations are not
read by this module.
"""

from __future__ import annotations

import copy
from collections import Counter
from collections.abc import Mapping
from typing import Any

from matching.library_quality import normalize_formula_text


MODULE_ID = "v15_identity_scope_resolution"
ALGORITHM_ID = "mssystem-v15-identity-scope-resolution-1.2.0"


def _identity(candidate: Mapping[str, Any] | None) -> str:
    return str((candidate or {}).get("identity_key") or "").strip()


def _formula(candidate: Mapping[str, Any]) -> str:
    display = candidate.get("candidate") or {}
    return normalize_formula_text(
        display.get("formula") if isinstance(display, Mapping) else None
    )


def _ranking_snapshot(payload: Mapping[str, Any]) -> list[tuple[Any, ...]]:
    return [
        (
            str(entity.get("entity_id") or ""),
            _identity(entity.get("winner")),
            _identity(entity.get("leading_candidate_when_abstained")),
            tuple(
                (
                    _identity(candidate),
                    candidate.get("entity_rank"),
                    candidate.get("entity_score"),
                )
                for candidate in (entity.get("candidates") or [])
                if isinstance(candidate, Mapping)
            ),
        )
        for entity in (payload.get("entities") or [])
        if isinstance(entity, Mapping)
    ]


def _competitive_candidates(
    entity: Mapping[str, Any], config: Mapping[str, Any]
) -> tuple[list[Mapping[str, Any]], dict[str, Any]]:
    """Return the complete bounded production shortlist.

    V15 has already bounded this list upstream (currently five candidates per
    precursor entity) and has already applied chemistry vetoes.  A second,
    fixed score threshold here used to silently discard lower-ranked live
    hypotheses.  That made the reported identity scope depend on a threshold
    whose error rate did not transport between the mixed-standard cohorts.

    The resolver is a demotion-only reporting guard, not a ranker.  It must
    therefore inspect every non-vetoed candidate that the production ranker
    retained.  Candidate order and scores remain untouched.
    """

    candidates = [
        candidate
        for candidate in (entity.get("candidates") or [])
        if isinstance(candidate, Mapping) and not bool(candidate.get("vetoed"))
    ]
    configured_maximum = int(config.get("maximum_ranked_shortlist_size") or 5)
    if configured_maximum < 1:
        raise ValueError("maximum_ranked_shortlist_size must be positive")
    if len(candidates) > configured_maximum:
        raise ValueError(
            "V15 entity contains more non-vetoed candidates than the locked "
            "identity-resolution shortlist"
        )
    return candidates, {
        "competition_set_rule": "all_non_vetoed_v15_ranked_shortlist_candidates",
        "maximum_ranked_shortlist_size": configured_maximum,
        "observed_non_vetoed_candidate_count": len(candidates),
        "fixed_score_threshold_used": False,
        "score_ratio_threshold_used": False,
        "collision_energy_threshold_used_for_candidate_exclusion": False,
        "ranking_lead_always_included_in_conflict_resolution": True,
    }


def apply_v15_identity_resolution(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve reportable scope without changing candidate order or scores."""

    output = copy.deepcopy(dict(payload))
    policy = output.get("policy") or {}
    config = policy.get("identity_resolution") or {}
    if config.get("enabled") is not True:
        raise ValueError("V15 identity-scope resolution must be enabled")
    if str(config.get("algorithm_id") or "") != ALGORITHM_ID:
        raise ValueError("V15 identity-scope resolution algorithm is not locked")
    before = _ranking_snapshot(output)
    resolution_statuses: Counter[str] = Counter()

    for entity in output.get("entities") or []:
        if not isinstance(entity, dict):
            continue
        decision = dict(entity.get("decision") or {})
        original_status = str(decision.get("status") or "")
        original_scope = str(decision.get("scope") or "")
        competitive, thresholds = _competitive_candidates(entity, config)
        formula_groups: dict[str, list[str]] = {}
        missing_formula_identities = []
        for candidate in competitive:
            identity = _identity(candidate)
            formula = _formula(candidate)
            if formula:
                formula_groups.setdefault(formula, []).append(identity)
            else:
                missing_formula_identities.append(identity)

        competitive_ids = [_identity(candidate) for candidate in competitive]
        formulae = sorted(formula_groups)
        selected = entity.get("winner") or entity.get(
            "leading_candidate_when_abstained"
        )
        selected_identity = _identity(selected)
        status = "single_competitive_connectivity"
        reportable_scope = original_scope
        reportable_formula = _formula(selected) if isinstance(selected, Mapping) else ""

        if len(competitive) > 1 and missing_formula_identities:
            status = "unresolved_incomplete_formula_event"
            reportable_scope = "unresolved_precursor_event"
            reportable_formula = ""
            decision.update(
                {
                    "status": "unresolved_formula_event",
                    "scope": reportable_scope,
                    "blockers": list(
                        dict.fromkeys(
                            [
                                *(decision.get("blockers") or []),
                                "competitive_connectivity_formula_assignment_incomplete",
                            ]
                        )
                    ),
                }
            )
        elif len(competitive) > 1 and len(formulae) > 1:
            status = "unresolved_cross_formula_event"
            reportable_scope = "unresolved_precursor_event"
            reportable_formula = ""
            decision.update(
                {
                    "status": "unresolved_formula_event",
                    "scope": reportable_scope,
                    "blockers": list(
                        dict.fromkeys(
                            [
                                *(decision.get("blockers") or []),
                                "competitive_connectivities_have_conflicting_formulae",
                            ]
                        )
                    ),
                }
            )
        elif len(competitive) > 1 and len(formulae) == 1:
            status = "unresolved_same_formula_isomer_group"
            reportable_scope = "molecular_formula_isomer_group"
            reportable_formula = formulae[0]
            decision.update(
                {
                    "status": "selected_isomer_group",
                    "scope": reportable_scope,
                    "blockers": list(
                        dict.fromkeys(
                            [
                                *(decision.get("blockers") or []),
                                "competitive_same_formula_connectivities",
                            ]
                        )
                    ),
                }
            )

        decision.update(
            {
                "reportable_identity_key": (
                    selected_identity
                    if reportable_scope == "2d_connectivity_candidate"
                    and decision.get("status") == "selected"
                    else None
                ),
                "reportable_formula": reportable_formula or None,
                "ranking_lead_identity_key": selected_identity or None,
                "ranking_lead_is_reportable_single_connectivity": bool(
                    reportable_scope == "2d_connectivity_candidate"
                    and decision.get("status") == "selected"
                ),
                "identity_resolution_module_id": MODULE_ID,
                "identity_resolution_algorithm_id": ALGORITHM_ID,
            }
        )
        entity["decision"] = decision
        entity["identity_resolution"] = {
            "module_id": MODULE_ID,
            "algorithm_id": ALGORITHM_ID,
            "status": status,
            "upstream_decision_status": original_status,
            "upstream_decision_scope": original_scope,
            "reportable_scope": reportable_scope,
            "reportable_formula": reportable_formula or None,
            "ranking_lead_identity_key": selected_identity or None,
            "competitive_candidate_identity_keys": competitive_ids,
            "competitive_candidate_count": len(competitive_ids),
            "candidate_formulae": formulae,
            "formula_groups": {
                formula: identities for formula, identities in sorted(formula_groups.items())
            },
            "missing_formula_candidate_identity_keys": missing_formula_identities,
            "thresholds": thresholds,
            "candidate_names_used": False,
            "ground_truth_used": False,
            "ranking_changed": False,
            "score_changed": False,
        }
        if isinstance(selected, dict):
            selected["reporting_semantics"] = (
                "ranking_lead_not_reportable_compound_identity"
                if reportable_scope == "unresolved_precursor_event"
                else "ranking_lead_within_unresolved_isomer_group"
                if reportable_scope == "molecular_formula_isomer_group"
                else "selected_candidate_at_upstream_scope"
            )
        resolution_statuses[status] += 1

    after = _ranking_snapshot(output)
    if after != before:
        raise RuntimeError("Identity-scope resolution changed V15 ranking or scores")

    decision_statuses = Counter(
        str((entity.get("decision") or {}).get("status") or "")
        for entity in output.get("entities") or []
        if isinstance(entity, Mapping)
    )
    decision_scopes = Counter(
        str((entity.get("decision") or {}).get("scope") or "")
        for entity in output.get("entities") or []
        if isinstance(entity, Mapping)
    )
    output.setdefault("metadata", {}).update(
        {
            "identity_resolution_status": "available",
            "identity_resolution_module_id": MODULE_ID,
            "identity_resolution_algorithm_id": ALGORITHM_ID,
            "identity_resolution_status_counts": dict(
                sorted(resolution_statuses.items())
            ),
            "decision_status_counts": dict(sorted(decision_statuses.items())),
            "decision_scope_counts": dict(sorted(decision_scopes.items())),
            "unresolved_cross_formula_entity_count": int(
                resolution_statuses.get("unresolved_cross_formula_event", 0)
            ),
            "unresolved_incomplete_formula_entity_count": int(
                resolution_statuses.get("unresolved_incomplete_formula_event", 0)
            ),
            "unresolved_same_formula_isomer_entity_count": int(
                resolution_statuses.get("unresolved_same_formula_isomer_group", 0)
            ),
        }
    )
    output.setdefault("invariants", {}).update(
        {
            "one_reportable_identity_scope_per_precursor_entity": True,
            "same_formula_competitive_connectivities_release_as_isomer_group": True,
            "cross_formula_competition_releases_compound_identity": False,
            "identity_resolution_changes_top1": False,
            "identity_resolution_changes_top5": False,
            "identity_resolution_changes_entity_score": False,
            "identity_resolution_uses_candidate_names": False,
            "identity_resolution_uses_ground_truth": False,
            "identity_resolution_uses_fixed_score_threshold": False,
            "all_non_vetoed_ranked_shortlist_candidates_enter_scope_resolution": True,
        }
    )
    return output


__all__ = [
    "ALGORITHM_ID",
    "MODULE_ID",
    "apply_v15_identity_resolution",
]
