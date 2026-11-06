"""Immutable contracts for the EviDiscrim-MS shadow sidecar.

The contract deliberately separates candidate-conditioned forward predictions
from measured evidence.  Nothing represented here is a posterior probability,
an independent confirmation, or permission to change the production result.
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


EVIDISCRIM_SCHEMA_VERSION = "2.0.0"
EVIDISCRIM_MODULE_ID = "evidiscrim_candidate_contrast"


class DiagnosticIon(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    mz: float = Field(gt=0)
    left_predicted_intensity: float = Field(ge=0, le=1)
    right_predicted_intensity: float = Field(ge=0, le=1)
    left_prediction_uncertainty: float = Field(ge=0, le=1)
    right_prediction_uncertainty: float = Field(ge=0, le=1)
    absolute_prediction_delta: float = Field(ge=0, le=1)
    reliability: float = Field(ge=0, le=1)
    diagnostic_weight: float = Field(ge=0, le=1)
    observed_relative_intensity: float = Field(ge=0, le=1)
    observed_present: bool
    absence_penalty_weight: float = Field(ge=0, le=1)
    predicted_direction: str


class EnergyPairSeparation(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    collision_energy: float = Field(gt=0)
    collision_energy_unit: str = Field(min_length=1)
    raw_spectral_distribution_divergence: float = Field(ge=0, le=1)
    uncertainty_reliability_factor: float = Field(ge=0, le=1)
    spectral_distribution_divergence: float = Field(ge=0, le=1)
    diagnostic_ion_count: int = Field(ge=0)
    diagnostic_ions: tuple[float, ...] = ()
    prediction_replicates_left: int = Field(ge=1)
    prediction_replicates_right: int = Field(ge=1)
    uncertainty_method: str = Field(min_length=1)
    architecture_families_left: tuple[str, ...] = ()
    architecture_families_right: tuple[str, ...] = ()
    common_architecture_families: tuple[str, ...] = ()
    common_architecture_family_count: int = Field(default=0, ge=0)
    prediction_providers_left: tuple[str, ...] = ()
    prediction_providers_right: tuple[str, ...] = ()
    ensemble_qualified: bool = False
    instrument_domain_qualified: bool = False
    between_architecture_uncertainty: float = Field(default=0, ge=0, le=1)
    within_architecture_uncertainty: float = Field(default=0, ge=0, le=1)


class ContrastivePair(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    pair_key: str = Field(min_length=1)
    left_identity_key: str = Field(min_length=1)
    right_identity_key: str = Field(min_length=1)
    left_compound_name: str
    right_compound_name: str
    eligible: bool
    status: str = Field(min_length=1)
    scope_reasons: tuple[str, ...] = ()
    abstention_reason: str | None = None
    observed_collision_energy: float | None = Field(default=None, gt=0)
    collision_energy_unit: str | None = None
    diagnostic_ions: tuple[DiagnosticIon, ...] = ()
    diagnostic_ion_count: int = Field(default=0, ge=0)
    total_diagnostic_weight: float = Field(default=0, ge=0)
    left_diagnostic_fit: float | None = Field(default=None, ge=0, le=1)
    right_diagnostic_fit: float | None = Field(default=None, ge=0, le=1)
    signed_fit_margin_left_minus_right: float | None = None
    counterfactual_winner_identity_key: str | None = None
    energy_separations: tuple[EnergyPairSeparation, ...] = ()
    ensemble_qualified: bool = False
    common_architecture_family_count: int = Field(default=0, ge=0)
    common_architecture_families: tuple[str, ...] = ()
    score_semantics: str = (
        "uncalibrated_candidate_conditioned_diagnostic_fit_not_probability"
    )
    runtime_influence: bool = False

    @model_validator(mode="after")
    def enforce_pair_boundaries(self) -> "ContrastivePair":
        if self.runtime_influence:
            raise ValueError("EviDiscrim-MS pairs cannot influence production")
        if not self.eligible and self.counterfactual_winner_identity_key:
            raise ValueError("An ineligible pair cannot declare a winner")
        if self.status != "resolved_counterfactual" and self.counterfactual_winner_identity_key:
            raise ValueError("Only a resolved counterfactual pair may declare a winner")
        if self.common_architecture_family_count != len(
            self.common_architecture_families
        ):
            raise ValueError(
                "common_architecture_family_count must match the family list"
            )
        if self.diagnostic_ion_count != len(self.diagnostic_ions):
            raise ValueError("diagnostic_ion_count must match diagnostic_ions")
        return self


class CandidateContrastSupport(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    identity_key: str = Field(min_length=1)
    compound_name: str
    production_top1: bool = False
    eligible_pairs: int = Field(ge=0)
    resolved_pairs: int = Field(ge=0)
    pair_wins: int = Field(ge=0)
    pair_losses: int = Field(ge=0)
    pair_abstentions: int = Field(ge=0)
    net_diagnostic_support: float
    counterfactual_rank: int | None = Field(default=None, ge=1)
    score_semantics: str = "signed_pairwise_fit_margin_not_probability"


class EnergySelection(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    status: str = Field(min_length=1)
    objective: str = "maximize_worst_pair_spectral_distribution_divergence"
    recommended_collision_energy: float | None = Field(default=None, gt=0)
    collision_energy_unit: str | None = None
    minimum_pair_separation: float | None = Field(default=None, ge=0, le=1)
    eligible_pair_count: int = Field(ge=0)
    covered_pair_count: int = Field(ge=0)
    energy_grid: tuple[float, ...] = ()
    per_energy: tuple[dict[str, Any], ...] = ()
    planning_only: bool = True
    universal_complete_fragmentation_assumed: bool = False
    cross_instrument_conversion_applied: bool = False

    @model_validator(mode="after")
    def enforce_energy_boundaries(self) -> "EnergySelection":
        if not self.planning_only:
            raise ValueError("Predicted CE selection must remain planning-only")
        if self.universal_complete_fragmentation_assumed:
            raise ValueError("A universal complete-fragmentation energy is forbidden")
        if self.cross_instrument_conversion_applied:
            raise ValueError("Cross-instrument energy conversion is forbidden")
        return self


class EviDiscrimSpectrum(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    row_index: int = Field(ge=0)
    source_row_sha256: str = Field(min_length=64, max_length=64)
    spectrum_id: str
    query_precursor_mz: float | None = Field(default=None, gt=0)
    observed_collision_energy: float | None = Field(default=None, gt=0)
    collision_energy_unit: str | None = None
    production_top1_identity_key: str | None = None
    candidate_count: int = Field(ge=0)
    eligible_pair_count: int = Field(ge=0)
    resolved_pair_count: int = Field(ge=0)
    abstained_pair_count: int = Field(ge=0)
    ineligible_pair_count: int = Field(ge=0)
    pairwise_comparisons: tuple[ContrastivePair, ...] = ()
    candidate_support: tuple[CandidateContrastSupport, ...] = ()
    counterfactual_top1_identity_key: str | None = None
    counterfactual_candidate_order_changed: bool = False
    spectrum_status: str = Field(min_length=1)
    spectrum_abstention_reason: str | None = None
    energy_selection: EnergySelection
    runtime_influence: bool = False

    @model_validator(mode="after")
    def enforce_spectrum_boundaries(self) -> "EviDiscrimSpectrum":
        if self.runtime_influence:
            raise ValueError("EviDiscrim-MS spectra cannot influence production")
        if self.counterfactual_candidate_order_changed and not self.counterfactual_top1_identity_key:
            raise ValueError("A changed counterfactual order requires a counterfactual Top-1")
        return self


class EviDiscrimArtifact(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = EVIDISCRIM_SCHEMA_VERSION
    task_id: str = Field(min_length=1)
    module_id: str = EVIDISCRIM_MODULE_ID
    policy_id: str = Field(min_length=1)
    status: str = Field(min_length=1)
    deployment_state: str = "shadow"
    runtime_influence: bool = False
    ranking_changed: bool = False
    confidence_changed: bool = False
    claim_scope_changed: bool = False
    high_promotions: int = Field(default=0, ge=0, le=0)
    candidate_conditioned: bool = True
    independent_measurement: bool = False
    source_results_sha256: str = Field(min_length=64, max_length=64)
    source_results_sha256_after: str = Field(min_length=64, max_length=64)
    source_results_mutated: bool = False
    model: dict[str, Any] = Field(default_factory=dict)
    runtime: dict[str, Any] = Field(default_factory=dict)
    spectra: tuple[EviDiscrimSpectrum, ...] = ()
    skipped: tuple[dict[str, Any], ...] = ()
    summary: dict[str, Any] = Field(default_factory=dict)
    interpretation: str

    @model_validator(mode="after")
    def enforce_shadow_invariants(self) -> "EviDiscrimArtifact":
        if self.deployment_state != "shadow":
            raise ValueError("EviDiscrim-MS must remain shadow")
        if any((self.runtime_influence, self.ranking_changed, self.confidence_changed,
                self.claim_scope_changed, self.source_results_mutated)):
            raise ValueError("EviDiscrim-MS cannot mutate production evidence")
        if self.high_promotions:
            raise ValueError("EviDiscrim-MS cannot promote High")
        if not self.candidate_conditioned or self.independent_measurement:
            raise ValueError("Forward predictions must remain candidate-conditioned")
        if self.source_results_sha256 != self.source_results_sha256_after:
            raise ValueError("Source result digest changed while building shadow artifact")
        if self.schema_version.startswith("2."):
            invalid = [
                pair.pair_key
                for spectrum in self.spectra
                for pair in spectrum.pairwise_comparisons
                if pair.status == "resolved_counterfactual"
                and not pair.ensemble_qualified
            ]
            if invalid:
                raise ValueError(
                    "Schema v2 resolved pairs require independent-architecture "
                    f"ensembles: {invalid[:3]}"
                )
            invalid_energy = []
            for spectrum in self.spectra:
                selection = spectrum.energy_selection
                if selection.status != "recommended_shadow":
                    continue
                matching_energy = [
                    separation
                    for pair in spectrum.pairwise_comparisons
                    if pair.eligible
                    for separation in pair.energy_separations
                    if separation.collision_energy
                    == selection.recommended_collision_energy
                    and separation.collision_energy_unit
                    == selection.collision_energy_unit
                ]
                if (
                    not matching_energy
                    or len(matching_energy) != selection.eligible_pair_count
                    or any(
                        not separation.ensemble_qualified
                        or not separation.instrument_domain_qualified
                        for separation in matching_energy
                    )
                ):
                    invalid_energy.append(spectrum.row_index)
            if invalid_energy:
                raise ValueError(
                    "Schema v2 CE recommendations require complete, domain-"
                    "qualified independent ensembles: "
                    f"{invalid_energy[:3]}"
                )
        return self


class EviDiscrimRunRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    collision_energies: tuple[float, ...] = ()
    max_candidates_per_spectrum: int = Field(default=5, ge=2, le=8)
    max_spectra: int = Field(default=25, ge=1, le=200)
    refresh: bool = False

    @field_validator("collision_energies")
    @classmethod
    def validate_energy_grid(cls, values: tuple[float, ...]) -> tuple[float, ...]:
        if any(value <= 0 for value in values):
            raise ValueError("collision energies must be positive")
        if tuple(sorted(set(values))) != values:
            raise ValueError("collision energies must be sorted and unique")
        return values


__all__ = [
    "CandidateContrastSupport",
    "ContrastivePair",
    "DiagnosticIon",
    "EnergyPairSeparation",
    "EnergySelection",
    "EVIDISCRIM_MODULE_ID",
    "EVIDISCRIM_SCHEMA_VERSION",
    "EviDiscrimArtifact",
    "EviDiscrimRunRequest",
    "EviDiscrimSpectrum",
]
