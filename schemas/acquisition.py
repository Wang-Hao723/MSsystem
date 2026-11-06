"""Versioned contracts for shadow-mode active evidence acquisition."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


ACQUISITION_SCHEMA_VERSION = "1.2.0"


class AcquisitionResourceDeclaration(BaseModel):
    """Explicit, fail-closed resources available to one planning request.

    Identity keys may be the planner's canonical ``ik::...`` keys or full
    InChIKeys.  Merely owning an instrument or having an unspecified standards
    cabinet is not sufficient: a positive action budget and candidate-specific
    reference identities are both required.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    declaration_id: str | None = None
    authentic_standard_identity_keys: tuple[str, ...] = ()
    authentic_standard_action_budget: int = Field(default=0, ge=0, le=200)
    ion_mobility_available: bool = False
    ccs_reference_identity_keys: tuple[str, ...] = ()
    ion_mobility_action_budget: int = Field(default=0, ge=0, le=200)

    @field_validator(
        "authentic_standard_identity_keys", "ccs_reference_identity_keys"
    )
    @classmethod
    def normalize_nonempty_keys(cls, values: tuple[str, ...]) -> tuple[str, ...]:
        cleaned = tuple(dict.fromkeys(
            value.strip() for value in values if str(value).strip()
        ))
        if len(cleaned) != len(values):
            raise ValueError("resource identity keys must be non-empty and unique")
        return cleaned

    @field_validator("declaration_id")
    @classmethod
    def normalize_declaration_id(cls, value: str | None) -> str | None:
        if value is None:
            return None
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("declaration_id cannot be blank")
        return cleaned

    @model_validator(mode="after")
    def enforce_fail_closed_resources(self) -> "AcquisitionResourceDeclaration":
        if (
            self.authentic_standard_action_budget > 0
            and not self.authentic_standard_identity_keys
        ):
            raise ValueError(
                "authentic-standard budget requires candidate-specific identity keys"
            )
        if self.ion_mobility_action_budget > 0 and not self.ion_mobility_available:
            raise ValueError(
                "ion-mobility budget requires explicit instrument availability"
            )
        if (
            self.ion_mobility_action_budget > 0
            and not self.ccs_reference_identity_keys
        ):
            raise ValueError(
                "ion-mobility budget requires candidate-specific CCS references"
            )
        return self


class AmbiguityCandidate(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    identity_key: str = Field(min_length=1)
    compound_name: str = Field(min_length=1)
    formula: str | None = None
    inchikey: str | None = None
    smiles: str | None = None
    adduct: str | None = None
    rank_score: float | None = None
    source_tier: str | None = None


class AmbiguityEdge(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    edge_id: str = Field(min_length=1)
    left: AmbiguityCandidate
    right: AmbiguityCandidate
    reasons: tuple[str, ...]
    unresolved: bool = True


class AcquisitionNeed(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    need_id: str = Field(min_length=1)
    task_id: str = Field(min_length=1)
    entity_id: str = Field(min_length=1)
    source_result_index: int = Field(ge=0)
    spectrum_id: str | None = None
    precursor_mz: float = Field(gt=0)
    retention_time_min: float | None = Field(default=None, ge=0)
    polarity: str | None = None
    candidate_independent_linkage_fields: tuple[str, ...]
    candidates: tuple[AmbiguityCandidate, ...] = Field(min_length=2)
    ambiguity_edges: tuple[AmbiguityEdge, ...] = Field(min_length=1)
    missing_evidence: tuple[str, ...]
    applicability: str = "applicable"
    baseline_snapshot: dict[str, Any] = Field(default_factory=dict)


class CollisionEnergyStage(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    stage_id: str = Field(min_length=1)
    collision_energies: tuple[float, ...] = ()
    acquire_when: str = Field(min_length=1)
    objective: str = Field(min_length=1)

    @model_validator(mode="after")
    def validate_energies(self) -> "CollisionEnergyStage":
        if any(value <= 0 for value in self.collision_energies):
            raise ValueError("collision energies must be positive")
        if tuple(sorted(set(self.collision_energies))) != self.collision_energies:
            raise ValueError("collision energies must be unique and sorted")
        return self


class CollisionEnergyDesign(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    profile_id: str = Field(min_length=1)
    strategy: str = Field(min_length=1)
    instrument_domain: str = Field(min_length=1)
    collision_energy_unit: str = Field(min_length=1)
    existing_collision_energy: float | None = Field(default=None, ge=0)
    observed_collision_energies: tuple[float, ...] = ()
    evidence_status: str = Field(min_length=1)
    stages: tuple[CollisionEnergyStage, ...] = Field(min_length=1)
    response_metrics: tuple[str, ...] = ()
    automatic_ceiling: float | None = Field(default=None, gt=0)
    numeric_schedule_available: bool
    universal_complete_fragmentation_assumed: bool = False
    cross_instrument_conversion_applied: bool = False
    method_owner_review_required: bool = True

    @model_validator(mode="after")
    def enforce_scientific_boundaries(self) -> "CollisionEnergyDesign":
        if self.universal_complete_fragmentation_assumed:
            raise ValueError("a universal complete-fragmentation energy is forbidden")
        if self.cross_instrument_conversion_applied:
            raise ValueError("cross-instrument collision-energy conversion is forbidden")
        if tuple(sorted(set(self.observed_collision_energies))) != self.observed_collision_energies:
            raise ValueError("observed collision energies must be unique and sorted")
        return self


class AcquisitionAction(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    action_id: str = Field(min_length=1)
    action_kind: str = Field(min_length=1)
    label: str = Field(min_length=1)
    precursor_mz: float = Field(gt=0)
    retention_time_min: float | None = Field(default=None, ge=0)
    rt_window_min: float = Field(default=0.15, gt=0)
    polarity: str | None = None
    collision_energies: tuple[float, ...] = ()
    collision_energy_unit: str | None = None
    collision_energy_design: CollisionEnergyDesign | None = None
    isolation_width_da: float | None = Field(default=None, gt=0)
    isolation_width_profile_id: str | None = None
    isolation_width_evidence_status: str | None = None
    max_injection_time_ms: float | None = Field(default=None, gt=0)
    diagnostic_ions: tuple[float, ...] = ()
    pairwise_separation: dict[str, float] = Field(default_factory=dict)
    minimum_pair_separation: float = Field(ge=0, le=1)
    measurement_success_factor: float = Field(ge=0, le=1)
    transferability_factor: float = Field(ge=0, le=1)
    purity_factor: float = Field(ge=0, le=1)
    duty_cycle_cost: float = Field(ge=0, le=1)
    risk_penalty: float = Field(ge=0, le=1)
    utility: float = Field(ge=0, le=1)
    utility_interval: tuple[float, float]
    rationale: tuple[str, ...]
    feasibility: str
    required_resources: tuple[str, ...] = ()
    resource_status: str = "not_required"
    resource_identity_keys: tuple[str, ...] = ()
    resource_declaration_id: str | None = None
    advisory_only: bool = True

    @model_validator(mode="after")
    def validate_interval(self) -> "AcquisitionAction":
        low, high = self.utility_interval
        if not 0 <= low <= self.utility <= high <= 1:
            raise ValueError("utility interval must contain utility inside [0, 1]")
        resource_gated = self.action_kind in {
            "authentic_standard_rt", "ion_mobility_ccs"
        }
        if resource_gated and self.resource_status != "explicitly_available":
            raise ValueError(
                "resource-gated actions require explicit resource availability"
            )
        if resource_gated and not self.required_resources:
            raise ValueError("resource-gated actions must declare required resources")
        if resource_gated and not self.resource_identity_keys:
            raise ValueError(
                "resource-gated actions require candidate-specific reference identities"
            )
        return self


class AcquisitionPlan(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = ACQUISITION_SCHEMA_VERSION
    plan_id: str = Field(min_length=1)
    policy_id: str = Field(min_length=1)
    task_id: str = Field(min_length=1)
    entity_id: str = Field(min_length=1)
    generated_at: str
    runtime_snapshot_id: str | None = None
    library_snapshot_id: str | None = None
    baseline_result_sha256: str
    module_state: str = "shadow"
    need: AcquisitionNeed
    ranked_actions: tuple[AcquisitionAction, ...] = Field(min_length=1)
    unsupported_pairs: tuple[str, ...] = ()
    provenance: dict[str, Any] = Field(default_factory=dict)
    production_ranking_changed: bool = False
    predicted_evidence_confirmatory: bool = False
    legacy_results_mutated: bool = False
    utility_semantics: str = "bounded_priority_not_probability"

    @model_validator(mode="after")
    def enforce_shadow_invariants(self) -> "AcquisitionPlan":
        if self.module_state != "shadow":
            raise ValueError("EviAcquire-MS plans must remain shadow")
        if self.production_ranking_changed:
            raise ValueError("An acquisition plan cannot change production ranking")
        if self.predicted_evidence_confirmatory:
            raise ValueError("Predicted evidence cannot be treated as confirmatory")
        return self


class EvidenceDelta(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    delta_id: str
    plan_id: str
    baseline_task_id: str
    followup_task_id: str
    linked_result_indices: tuple[int, ...]
    linkage_fields_used: tuple[str, ...]
    candidate_fields_used_for_linkage: tuple[str, ...] = ()
    baseline_candidate_keys: tuple[str, ...]
    followup_candidate_keys: tuple[str, ...]
    resolved_candidate_keys: tuple[str, ...]
    evidence_change: str
    observed_at: str
    production_ranking_changed: bool = False
    automatic_claim_release: bool = False

    @model_validator(mode="after")
    def enforce_candidate_independent_linkage(self) -> "EvidenceDelta":
        if self.candidate_fields_used_for_linkage:
            raise ValueError("Follow-up linkage must be candidate-independent")
        return self


class AcquisitionPlanRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    max_plans: int = Field(default=25, ge=1, le=200)
    max_actions_per_plan: int = Field(default=6, ge=1, le=12)
    refresh: bool = False


class FollowupLinkRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    followup_task_id: str = Field(min_length=1)
    mz_tolerance_ppm: float = Field(default=10.0, gt=0, le=100)
    rt_tolerance_min: float = Field(default=0.20, gt=0, le=5)
