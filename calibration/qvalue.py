"""Tie-aware empirical q-values for one-call-per-entity competitions."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Sequence

from .models import (
    CompetitionPair,
    CompetitionRole,
    QValuePoint,
    WinnerCall,
)


def compete_pairs(pairs: Iterable[CompetitionPair]) -> tuple[WinnerCall, ...]:
    """Resolve every target/decoy pair, awarding exact ties to the decoy."""
    winners = []
    seen_competitions: set[str] = set()
    seen_entities: set[str] = set()
    for pair in pairs:
        if pair.competition_id in seen_competitions:
            raise ValueError(f"Duplicate competition_id: {pair.competition_id}")
        if pair.entity_id in seen_entities:
            raise ValueError(f"Entity appears in multiple competitions: {pair.entity_id}")
        seen_competitions.add(pair.competition_id)
        seen_entities.add(pair.entity_id)
        role = (
            CompetitionRole.TARGET
            if pair.target_score > pair.decoy_score
            else CompetitionRole.DECOY
        )
        winners.append(
            WinnerCall(
                competition_id=pair.competition_id,
                entity_id=pair.entity_id,
                unit_level=pair.unit_level,
                claim_scope=pair.claim_scope,
                target_record_id=pair.target_record_id,
                role=role,
                score=max(pair.target_score, pair.decoy_score),
                polarity=pair.polarity,
                instrument_family=pair.instrument_family,
                source_tier=pair.source_tier,
            )
        )
    return tuple(winners)


def estimate_qvalue_curve(
    calls: Sequence[WinnerCall],
    *,
    decoy_scale: float,
    pseudocount: float = 1.0,
) -> tuple[QValuePoint, ...]:
    """Estimate a monotone q-value curve with score ties handled as one block."""
    if decoy_scale <= 0 or not math.isfinite(decoy_scale):
        raise ValueError("decoy_scale must be a finite positive number")
    if pseudocount < 0 or not math.isfinite(pseudocount):
        raise ValueError("pseudocount must be a finite non-negative number")
    grouped: dict[float, list[CompetitionRole]] = defaultdict(list)
    for call in calls:
        grouped[float(call.score)].append(call.role)
    cumulative_target = 0
    cumulative_decoy = 0
    rows: list[dict] = []
    for threshold in sorted(grouped, reverse=True):
        roles = grouped[threshold]
        cumulative_target += sum(role == CompetitionRole.TARGET for role in roles)
        cumulative_decoy += sum(role == CompetitionRole.DECOY for role in roles)
        estimated_false_targets = (cumulative_decoy + pseudocount) * decoy_scale
        raw_fdr = min(1.0, estimated_false_targets / max(cumulative_target, 1))
        rows.append(
            {
                "score_threshold": threshold,
                "target_winners": cumulative_target,
                "decoy_winners": cumulative_decoy,
                "raw_fdr": raw_fdr,
            }
        )
    running_min = 1.0
    qvalues = [1.0] * len(rows)
    for index in range(len(rows) - 1, -1, -1):
        running_min = min(running_min, rows[index]["raw_fdr"])
        qvalues[index] = running_min
    return tuple(
        QValuePoint(**row, q_value=qvalues[index])
        for index, row in enumerate(rows)
    )


def qvalue_for_score(curve: Sequence[QValuePoint], score: float) -> float:
    """Return the conservative curve value applicable to a target score."""
    if not curve:
        raise ValueError("Cannot assign a q-value from an empty curve")
    for point in curve:
        if score >= point.score_threshold or math.isclose(
            score, point.score_threshold, rel_tol=0.0, abs_tol=1e-12
        ):
            return point.q_value
    return curve[-1].q_value
