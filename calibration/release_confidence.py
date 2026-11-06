"""Fail-closed shadow qualification for V15 release-confidence changes.

This module never changes a production decision.  It provides strict contracts for
partition isolation, grouped cross-fitting, a fixed margin-only shadow rule,
paired inference helpers, and the final review gate.  Correctness labels are not
part of :class:`ShadowCandidate`, so the shadow rule cannot read them.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from enum import Enum
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


RELEASE_CONFIDENCE_SCHEMA_VERSION = "mssystem-v15-release-confidence-shadow/1.0.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    """Return a deterministic SHA-256 for a JSON-compatible object."""
    return hashlib.sha256(_canonical(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PartitionRole(str, Enum):
    DEVELOPMENT = "development"
    TUNING = "tuning"
    SEALED_VALIDATION = "sealed_validation"


class IsolationRecord(BaseModel):
    """Label-free identifiers used to prove partition separation."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    record_id: str
    partition: PartitionRole
    connectivity_id: str
    molecular_formula: str
    formula_family_id: str
    source_dataset_id: str
    source_record_ids: tuple[str, ...]
    source_content_fingerprints: tuple[str, ...]
    laboratory_id: str
    instrument_instance_id: str
    instrument_domain_id: str
    acquisition_batch_id: str
    polarity: str
    adduct: str

    @model_validator(mode="after")
    def validate_source_identifiers(self) -> "IsolationRecord":
        if len(set(self.source_record_ids)) != len(self.source_record_ids):
            raise ValueError("source_record_ids must be unique within a record")
        if len(set(self.source_content_fingerprints)) != len(
            self.source_content_fingerprints
        ):
            raise ValueError(
                "source_content_fingerprints must be unique within a record"
            )
        return self


class IsolationAudit(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    eligible: bool
    reason_codes: tuple[str, ...]
    partition_counts: dict[str, int]
    overlaps: dict[str, tuple[str, ...]]
    coarse_instrument_domain_overlap: tuple[str, ...]
    transfer_scope: str


_UNKNOWN_TOKENS = {"", "unknown", "none", "null", "na", "n/a", "unavailable"}


def _is_unknown(value: str) -> bool:
    return str(value).strip().lower() in _UNKNOWN_TOKENS


def _values_by_partition(
    records: Sequence[IsolationRecord],
    field: str,
    *,
    flatten: bool = False,
) -> dict[PartitionRole, set[str]]:
    output = {role: set() for role in PartitionRole}
    for record in records:
        value = getattr(record, field)
        if flatten:
            output[record.partition].update(str(item) for item in value)
        else:
            output[record.partition].add(str(value))
    return output


def _cross_partition_overlap(values: Mapping[PartitionRole, set[str]]) -> set[str]:
    seen: dict[str, set[PartitionRole]] = defaultdict(set)
    for role, role_values in values.items():
        for value in role_values:
            seen[value].add(role)
    return {value for value, roles in seen.items() if len(roles) > 1}


def audit_partition_isolation(
    records: Iterable[IsolationRecord],
) -> IsolationAudit:
    """Audit identity, formula, source, laboratory, instrument and batch isolation.

    Physical instrument instances must be disjoint.  A coarse acquisition domain
    may recur on a different physical instrument in a different laboratory; that
    supports only same-domain transfer, not a new-domain claim.
    """

    materialized = tuple(records)
    record_ids = [record.record_id for record in materialized]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("Isolation manifest record_id values must be unique")

    counts = Counter(record.partition for record in materialized)
    reasons: list[str] = []
    for role in PartitionRole:
        if counts[role] == 0:
            reasons.append(f"missing_{role.value}_partition")

    required_scalar_fields = (
        "connectivity_id",
        "molecular_formula",
        "formula_family_id",
        "source_dataset_id",
        "laboratory_id",
        "instrument_instance_id",
        "instrument_domain_id",
        "acquisition_batch_id",
        "polarity",
        "adduct",
    )
    for field in required_scalar_fields:
        if any(_is_unknown(getattr(record, field)) for record in materialized):
            reasons.append(f"missing_or_unknown_{field}")
    for field in ("source_record_ids", "source_content_fingerprints"):
        if any(
            not getattr(record, field)
            or any(_is_unknown(item) for item in getattr(record, field))
            for record in materialized
        ):
            reasons.append(f"missing_or_unknown_{field}")

    blocking_fields = {
        "connectivity_id": False,
        "formula_family_id": False,
        "source_dataset_id": False,
        "source_record_ids": True,
        "source_content_fingerprints": True,
        "laboratory_id": False,
        "instrument_instance_id": False,
        "acquisition_batch_id": False,
    }
    overlaps: dict[str, tuple[str, ...]] = {}
    for field, flatten in blocking_fields.items():
        overlap = sorted(
            _cross_partition_overlap(
                _values_by_partition(materialized, field, flatten=flatten)
            )
        )
        if overlap:
            overlaps[field] = tuple(overlap)
            reasons.append(f"cross_partition_{field}_overlap")

    domain_overlap = sorted(
        _cross_partition_overlap(
            _values_by_partition(materialized, "instrument_domain_id")
        )
    )
    if domain_overlap:
        transfer_scope = "same_domain_new_laboratory_and_instrument_instance"
    else:
        transfer_scope = "new_domain_new_laboratory_and_instrument_instance"

    reasons = list(dict.fromkeys(reasons))
    return IsolationAudit(
        eligible=not reasons,
        reason_codes=tuple(reasons),
        partition_counts={role.value: counts[role] for role in PartitionRole},
        overlaps=overlaps,
        coarse_instrument_domain_overlap=tuple(domain_overlap),
        transfer_scope=transfer_scope,
    )


class CrossfitAssignment(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    record_id: str
    fold: int = Field(ge=0)
    leakage_component_id: str


class CrossfitAudit(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    eligible: bool
    reason_codes: tuple[str, ...]
    fold_counts: dict[int, int]
    leakage_identifiers_crossing_folds: dict[str, tuple[str, ...]]


def _crossfit_identifier_values(record: IsolationRecord) -> dict[str, tuple[str, ...]]:
    return {
        "connectivity_id": (record.connectivity_id,),
        "formula_family_id": (record.formula_family_id,),
        "source_record_ids": record.source_record_ids,
        "source_content_fingerprints": record.source_content_fingerprints,
    }


def build_crossfit_assignments(
    records: Iterable[IsolationRecord],
    *,
    n_folds: int,
    seed: int,
) -> tuple[CrossfitAssignment, ...]:
    """Assign leakage-linked components to deterministic, blocked folds."""

    materialized = tuple(records)
    if n_folds < 2:
        raise ValueError("Cross-fitting requires at least two folds")
    if len({record.record_id for record in materialized}) != len(materialized):
        raise ValueError("Cross-fit record_id values must be unique")
    if not materialized:
        raise ValueError("Cross-fit records cannot be empty")

    parent = list(range(len(materialized)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    owners: dict[tuple[str, str], int] = {}
    for index, record in enumerate(materialized):
        for field, values in _crossfit_identifier_values(record).items():
            for value in values:
                key = (field, value)
                previous = owners.get(key)
                if previous is None:
                    owners[key] = index
                else:
                    union(index, previous)

    components: dict[int, list[int]] = defaultdict(list)
    for index in range(len(materialized)):
        components[find(index)].append(index)
    if len(components) < n_folds:
        raise ValueError(
            "Leakage-linked components are fewer than requested cross-fit folds"
        )

    component_rows: list[tuple[str, list[int], Counter[tuple[str, str, str]]]] = []
    for indices in components.values():
        record_ids = sorted(materialized[index].record_id for index in indices)
        component_id = f"component:{semantic_sha256(record_ids)[:24]}"
        strata = Counter(
            (
                materialized[index].instrument_domain_id,
                materialized[index].polarity,
                materialized[index].adduct,
            )
            for index in indices
        )
        component_rows.append((component_id, indices, strata))
    component_rows.sort(
        key=lambda item: (
            -len(item[1]),
            semantic_sha256({"seed": seed, "component_id": item[0]}),
        )
    )

    fold_sizes = [0] * n_folds
    fold_strata: list[Counter[tuple[str, str, str]]] = [
        Counter() for _ in range(n_folds)
    ]
    output: list[CrossfitAssignment] = []
    for component_id, indices, strata in component_rows:
        fold = min(
            range(n_folds),
            key=lambda candidate: (
                sum(
                    fold_strata[candidate][stratum] * count
                    for stratum, count in strata.items()
                ),
                fold_sizes[candidate],
                semantic_sha256(
                    {
                        "seed": seed,
                        "component_id": component_id,
                        "fold": candidate,
                    }
                ),
            ),
        )
        fold_sizes[fold] += len(indices)
        fold_strata[fold].update(strata)
        output.extend(
            CrossfitAssignment(
                record_id=materialized[index].record_id,
                fold=fold,
                leakage_component_id=component_id,
            )
            for index in indices
        )
    if any(size == 0 for size in fold_sizes):
        raise ValueError("Cross-fit assignment produced an empty fold")
    return tuple(sorted(output, key=lambda item: item.record_id))


def audit_crossfit_assignments(
    records: Iterable[IsolationRecord],
    assignments: Iterable[CrossfitAssignment],
) -> CrossfitAudit:
    materialized = tuple(records)
    assigned = tuple(assignments)
    record_lookup = {record.record_id: record for record in materialized}
    assignment_lookup = {item.record_id: item for item in assigned}
    reasons: list[str] = []
    if len(record_lookup) != len(materialized):
        reasons.append("duplicate_crossfit_record_id")
    if len(assignment_lookup) != len(assigned):
        reasons.append("duplicate_crossfit_assignment")
    if set(record_lookup) != set(assignment_lookup):
        reasons.append("crossfit_assignment_record_mismatch")

    crossing: dict[str, tuple[str, ...]] = {}
    for field in (
        "connectivity_id",
        "formula_family_id",
        "source_record_ids",
        "source_content_fingerprints",
    ):
        identifier_folds: dict[str, set[int]] = defaultdict(set)
        for record_id, record in record_lookup.items():
            assignment = assignment_lookup.get(record_id)
            if assignment is None:
                continue
            values = _crossfit_identifier_values(record)[field]
            for value in values:
                identifier_folds[value].add(assignment.fold)
        invalid = sorted(
            value for value, folds in identifier_folds.items() if len(folds) > 1
        )
        if invalid:
            crossing[field] = tuple(invalid)
            reasons.append(f"crossfit_{field}_leakage")

    fold_counts = Counter(item.fold for item in assigned)
    if len(fold_counts) < 2:
        reasons.append("crossfit_has_fewer_than_two_nonempty_folds")
    reasons = list(dict.fromkeys(reasons))
    return CrossfitAudit(
        eligible=not reasons,
        reason_codes=tuple(reasons),
        fold_counts=dict(sorted(fold_counts.items())),
        leakage_identifiers_crossing_folds=crossing,
    )


class ShadowRule(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    algorithm_id: str
    base_minimum_winner_margin: float = Field(ge=0.0)
    current_minimum_entity_score: float = Field(ge=0.0)
    current_minimum_winner_margin: float = Field(ge=0.0)
    current_minimum_ce_support: int = Field(ge=1)
    only_allowed_current_gate_failure: str
    minimum_independent_acquisition_support: int = Field(ge=2)
    maximum_same_formula_connectivities: int = Field(ge=1)
    minimum_unique_compatible_reference_profiles: int = Field(ge=1)
    minimum_compatible_reference_provenance_strata: int = Field(ge=1)
    minimum_development_domain_cases: int = Field(ge=1)
    historical_error_sentinel_record_ids: tuple[str, ...] = ()

    @model_validator(mode="after")
    def validate_margin_interval(self) -> "ShadowRule":
        if self.base_minimum_winner_margin >= self.current_minimum_winner_margin:
            raise ValueError(
                "The base release margin must be lower than the current gate margin"
            )
        if len(set(self.historical_error_sentinel_record_ids)) != len(
            self.historical_error_sentinel_record_ids
        ):
            raise ValueError("Historical error sentinel IDs must be unique")
        return self


class ShadowCandidate(BaseModel):
    """Label-free evidence presented to the fixed shadow rule."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    record_id: str
    base_formal_2d_release: bool
    current_formal_2d_release: bool
    current_gate_failed_checks: tuple[str, ...]
    winner_entity_score: float = Field(allow_inf_nan=False)
    winner_margin: float = Field(ge=0.0, allow_inf_nan=False)
    winner_ce_support: int = Field(ge=0)
    query_quality_status: str | None = None
    query_quality_failed_checks: tuple[str, ...] = ()
    same_formula_connectivity_count: int | None = Field(default=None, ge=0)
    unique_compatible_reference_profiles: int | None = Field(default=None, ge=0)
    compatible_reference_provenance_strata: int | None = Field(default=None, ge=0)
    independent_acquisition_support_count: int | None = Field(default=None, ge=0)
    independent_acquisition_top1_consistent: bool | None = None
    development_domain_case_count: int | None = Field(default=None, ge=0)
    formula_conflict: bool | None = None
    adduct_conflict: bool | None = None
    precursor_conflict: bool | None = None
    sample_source_compatibility_conflict: bool | None = None
    historical_error_sentinel: bool = False
    eligibility_metadata_complete: bool = False


class ShadowDecision(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = RELEASE_CONFIDENCE_SCHEMA_VERSION
    record_id: str
    algorithm_id: str
    eligible: bool
    status: str
    reason_codes: tuple[str, ...]
    would_restore_2d_connectivity_release: bool
    runtime_influence: bool = False
    production_release_changed: bool = False
    top1_changed: bool = False
    top5_changed: bool = False
    entity_score_changed: bool = False
    claim_scope: Literal["2d_connectivity_candidate"] = "2d_connectivity_candidate"
    authentic_standard_confirmation: bool = False
    q_value: None = None
    probability_correct: None = None
    score_is_probability: bool = False


def evaluate_shadow_candidate(
    candidate: ShadowCandidate,
    rule: ShadowRule,
) -> ShadowDecision:
    """Evaluate a candidate without changing the current release decision."""

    reasons: list[str] = []
    if candidate.historical_error_sentinel or (
        candidate.record_id in rule.historical_error_sentinel_record_ids
    ):
        reasons.append("historical_error_sentinel_locked")
    if not candidate.base_formal_2d_release:
        reasons.append("not_an_existing_base_2d_release")
    if candidate.current_formal_2d_release:
        reasons.append("already_released_by_current_gate")
    if candidate.current_gate_failed_checks != (
        rule.only_allowed_current_gate_failure,
    ):
        reasons.append("not_a_margin_only_current_gate_failure")
    if candidate.winner_entity_score < rule.current_minimum_entity_score:
        reasons.append("winner_score_below_current_gate")
    if candidate.winner_ce_support < rule.current_minimum_ce_support:
        reasons.append("winner_ce_support_below_current_gate")
    if not (
        rule.base_minimum_winner_margin
        <= candidate.winner_margin
        < rule.current_minimum_winner_margin
    ):
        reasons.append("winner_margin_outside_fixed_shadow_interval")
    if candidate.query_quality_status != "passed":
        reasons.append("query_quality_not_passed")
    if candidate.query_quality_failed_checks:
        reasons.append("query_quality_has_failed_checks")
    if not candidate.eligibility_metadata_complete:
        reasons.append("eligibility_metadata_incomplete")

    count_requirements = (
        (
            "same_formula_connectivity_count",
            candidate.same_formula_connectivity_count,
            lambda value: 1 <= value <= rule.maximum_same_formula_connectivities,
        ),
        (
            "unique_compatible_reference_profiles",
            candidate.unique_compatible_reference_profiles,
            lambda value: value >= rule.minimum_unique_compatible_reference_profiles,
        ),
        (
            "compatible_reference_provenance_strata",
            candidate.compatible_reference_provenance_strata,
            lambda value: value >= rule.minimum_compatible_reference_provenance_strata,
        ),
        (
            "independent_acquisition_support_count",
            candidate.independent_acquisition_support_count,
            lambda value: value >= rule.minimum_independent_acquisition_support,
        ),
        (
            "development_domain_case_count",
            candidate.development_domain_case_count,
            lambda value: value >= rule.minimum_development_domain_cases,
        ),
    )
    for field, value, predicate in count_requirements:
        if value is None:
            reasons.append(f"missing_{field}")
        elif not predicate(value):
            reasons.append(f"failed_{field}")
    if candidate.independent_acquisition_top1_consistent is not True:
        reasons.append("independent_acquisition_top1_not_consistent")

    for field in (
        "formula_conflict",
        "adduct_conflict",
        "precursor_conflict",
        "sample_source_compatibility_conflict",
    ):
        value = getattr(candidate, field)
        if value is None:
            reasons.append(f"missing_{field}_audit")
        elif value:
            reasons.append(field)

    reasons = list(dict.fromkeys(reasons))
    eligible = not reasons
    return ShadowDecision(
        record_id=candidate.record_id,
        algorithm_id=rule.algorithm_id,
        eligible=eligible,
        status="eligible_shadow_only" if eligible else "withheld_shadow",
        reason_codes=tuple(reasons),
        would_restore_2d_connectivity_release=eligible,
    )


def zero_event_one_sided_upper_bound(
    n: int,
    *,
    alpha: float = 0.05,
) -> float:
    """Exact one-sided Clopper-Pearson upper bound when zero events occur."""

    if n < 1:
        raise ValueError("n must be positive")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    return 1.0 - alpha ** (1.0 / n)


def minimum_zero_event_sample_size(
    maximum_event_rate: float,
    *,
    alpha: float = 0.05,
) -> int:
    """Minimum n for the zero-event upper bound to meet a fixed rate margin."""

    if not 0.0 < maximum_event_rate < 1.0:
        raise ValueError("maximum_event_rate must be in (0, 1)")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    return math.ceil(math.log(alpha) / math.log(1.0 - maximum_event_rate))


def exact_two_sided_mcnemar_pvalue(left_only: int, right_only: int) -> float:
    """Exact two-sided McNemar p-value from the two discordant counts."""

    if left_only < 0 or right_only < 0:
        raise ValueError("Discordant counts cannot be negative")
    discordant = left_only + right_only
    if discordant == 0:
        return 1.0
    tail = sum(
        math.comb(discordant, index) for index in range(min(left_only, right_only) + 1)
    ) / (2**discordant)
    return min(1.0, 2.0 * tail)


def exact_one_sided_mcnemar_pvalue(
    shadow_gain_only: int,
    current_gain_only: int,
) -> float:
    """Exact McNemar p-value for the prespecified shadow-gain alternative.

    The alternative is that the shadow rule has more favorable discordant
    outcomes than the current gate. Under the null, either member of a
    discordant pair is favorable with probability 0.5.
    """

    if shadow_gain_only < 0 or current_gain_only < 0:
        raise ValueError("Discordant counts cannot be negative")
    discordant = shadow_gain_only + current_gain_only
    if discordant == 0:
        return 1.0
    return sum(
        math.comb(discordant, index)
        for index in range(shadow_gain_only, discordant + 1)
    ) / (2**discordant)


class PairedReleaseOutcome(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    record_id: str
    connectivity_cluster_id: str
    formula_family_id: str
    current_correct_release: bool
    shadow_correct_release: bool
    current_wrong_release: bool
    shadow_wrong_release: bool


class BootstrapInterval(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    estimate: float
    lower: float
    upper: float
    confidence_level: float
    cluster_field: str
    resamples: int
    seed: int


def _quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot compute a quantile from no values")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def paired_cluster_bootstrap_difference(
    outcomes: Iterable[PairedReleaseOutcome],
    *,
    endpoint: Literal["correct_release", "wrong_release"],
    cluster_field: Literal["connectivity_cluster_id", "formula_family_id"],
    resamples: int = 10_000,
    seed: int = 20260827,
    confidence_level: float = 0.95,
) -> BootstrapInterval:
    """Cluster bootstrap the mean paired shadow-minus-current difference."""

    materialized = tuple(outcomes)
    if not materialized:
        raise ValueError("Paired outcomes cannot be empty")
    if resamples < 100:
        raise ValueError("At least 100 bootstrap resamples are required")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be in (0, 1)")

    current_field = f"current_{endpoint}"
    shadow_field = f"shadow_{endpoint}"
    clusters: dict[str, list[PairedReleaseOutcome]] = defaultdict(list)
    for row in materialized:
        clusters[getattr(row, cluster_field)].append(row)
    cluster_ids = sorted(clusters)

    def difference(rows: Sequence[PairedReleaseOutcome]) -> float:
        return sum(
            int(getattr(row, shadow_field)) - int(getattr(row, current_field))
            for row in rows
        ) / len(rows)

    estimate = difference(materialized)
    rng = random.Random(seed)
    samples: list[float] = []
    for _ in range(resamples):
        sampled_rows: list[PairedReleaseOutcome] = []
        for _ in cluster_ids:
            sampled_rows.extend(clusters[rng.choice(cluster_ids)])
        samples.append(difference(sampled_rows))
    tail = (1.0 - confidence_level) / 2.0
    return BootstrapInterval(
        estimate=estimate,
        lower=_quantile(samples, tail),
        upper=_quantile(samples, 1.0 - tail),
        confidence_level=confidence_level,
        cluster_field=cluster_field,
        resamples=resamples,
        seed=seed,
    )


class QualificationEvidence(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    split_isolation_passed: bool
    crossfit_audit_passed: bool
    rule_frozen_before_sealed_labels: bool
    external_validation_complete: bool
    multiplicity_plan_followed: bool
    historical_error_sentinels_rescued: int = Field(ge=0)
    sealed_evaluable_units: int = Field(ge=0)
    sealed_new_wrong_releases: int = Field(ge=0)
    sealed_correct_releases_recovered: int = Field(ge=0)
    sealed_correct_releases_lost: int = Field(default=0, ge=0)
    project_target_decoy_ready: bool = False


class QualificationDecision(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    eligible_for_runtime_review: bool
    reason_codes: tuple[str, ...]
    sealed_wrong_release_upper_bound: float | None
    sealed_recovery_one_sided_mcnemar_pvalue: float | None
    required_zero_event_sample_size: int
    minimum_correct_gains_when_no_losses: int
    qvalues_must_remain_null: bool
    production_activation_automatic: bool = False


def evaluate_runtime_review_qualification(
    evidence: QualificationEvidence,
    *,
    maximum_added_wrong_release_rate: float = 0.025,
    alpha: float = 0.05,
) -> QualificationDecision:
    """Apply the ordered safety-then-recovery review gate."""

    reasons: list[str] = []
    if not evidence.split_isolation_passed:
        reasons.append("split_isolation_failed")
    if not evidence.crossfit_audit_passed:
        reasons.append("crossfit_audit_failed")
    if not evidence.rule_frozen_before_sealed_labels:
        reasons.append("rule_not_frozen_before_sealed_labels")
    if not evidence.external_validation_complete:
        reasons.append("external_validation_incomplete")
    if not evidence.multiplicity_plan_followed:
        reasons.append("multiplicity_plan_not_followed")
    if evidence.historical_error_sentinels_rescued:
        reasons.append("historical_error_sentinel_regression")
    if evidence.sealed_new_wrong_releases:
        reasons.append("sealed_new_wrong_release_detected")
    if evidence.sealed_correct_releases_recovered < 1:
        reasons.append("no_sealed_correct_release_recovery")
    if evidence.sealed_correct_releases_lost:
        reasons.append("shadow_removed_current_correct_release")

    recovery_pvalue = None
    recovery_discordant = (
        evidence.sealed_correct_releases_recovered
        + evidence.sealed_correct_releases_lost
    )
    if recovery_discordant:
        recovery_pvalue = exact_one_sided_mcnemar_pvalue(
            evidence.sealed_correct_releases_recovered,
            evidence.sealed_correct_releases_lost,
        )
        if recovery_pvalue > alpha:
            reasons.append("sealed_recovery_mcnemar_not_significant")

    required_n = minimum_zero_event_sample_size(
        maximum_added_wrong_release_rate, alpha=alpha
    )
    upper_bound = None
    if evidence.sealed_evaluable_units:
        if evidence.sealed_new_wrong_releases == 0:
            upper_bound = zero_event_one_sided_upper_bound(
                evidence.sealed_evaluable_units, alpha=alpha
            )
            if upper_bound > maximum_added_wrong_release_rate:
                reasons.append("sealed_safety_interval_too_wide")
        else:
            reasons.append("zero_event_safety_bound_not_applicable")
    else:
        reasons.append("no_sealed_evaluable_units")

    reasons = list(dict.fromkeys(reasons))
    minimum_correct_gains = math.ceil(math.log(alpha) / math.log(0.5))
    return QualificationDecision(
        eligible_for_runtime_review=not reasons,
        reason_codes=tuple(reasons),
        sealed_wrong_release_upper_bound=upper_bound,
        sealed_recovery_one_sided_mcnemar_pvalue=recovery_pvalue,
        required_zero_event_sample_size=required_n,
        minimum_correct_gains_when_no_losses=minimum_correct_gains,
        # Release qualification cannot authorize q-value assignment. A separate
        # entity-level competition receipt and registry binding are required.
        qvalues_must_remain_null=True,
    )


def load_shadow_protocol(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != RELEASE_CONFIDENCE_SCHEMA_VERSION:
        raise ValueError("Unsupported release-confidence shadow protocol")
    return payload


def shadow_rule_from_protocol(protocol: Mapping[str, Any]) -> ShadowRule:
    payload = dict(protocol["fixed_shadow_rule"])
    payload["historical_error_sentinel_record_ids"] = tuple(
        protocol["historical_frozen_audit"]["historical_error_sentinel_record_ids"]
    )
    return ShadowRule.model_validate(payload)


def audit_qvalue_registry(path: Path) -> dict[str, Any]:
    """Return the registry state without assigning any per-record q-value."""

    registry = json.loads(path.read_text(encoding="utf-8"))
    production_id = registry.get("production_calibration_id")
    if production_id is None:
        return {
            "status": "unavailable",
            "production_calibration_id": None,
            "q_value": None,
            "qvalues_must_remain_null": True,
            "reason_codes": [
                "no_production_project_target_decoy_calibration_registered"
            ],
        }
    entry = next(
        (
            item
            for item in registry.get("entries", [])
            if item.get("calibration_id") == production_id
        ),
        None,
    )
    if entry is None or entry.get("production_eligible") is not True:
        return {
            "status": "unavailable",
            "production_calibration_id": production_id,
            "q_value": None,
            "qvalues_must_remain_null": True,
            "reason_codes": ["production_calibration_registry_entry_not_eligible"],
        }
    return {
        "status": "registered_assignment_still_requires_exact_project_binding",
        "production_calibration_id": production_id,
        "q_value": None,
        "qvalues_must_remain_null": True,
        "reason_codes": ["no_record_specific_competition_assignment_supplied"],
    }
