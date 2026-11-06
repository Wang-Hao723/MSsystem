"""Typed contracts for project-local target-decoy calibration."""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, model_validator

from schemas.decision import ClaimScope


CALIBRATION_SCHEMA_VERSION = "1.0.0"


class CalibrationUnitLevel(str, Enum):
    PRECURSOR_EVENT = "precursor_event"
    CHROMATOGRAPHIC_FEATURE = "chromatographic_feature"
    ION_IDENTITY_GROUP = "ion_identity_group"
    NEUTRAL_ANALYTE = "neutral_analyte"


class CompetitionRole(str, Enum):
    TARGET = "target"
    DECOY = "decoy"


class ProjectCalibrationStatus(str, Enum):
    CALIBRATED = "calibrated"
    PARTIAL = "partial"
    UNAVAILABLE = "unavailable"


class CalibrationMethod(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = CALIBRATION_SCHEMA_VERSION
    method_id: str
    method_version: str
    score_field: str
    higher_score_is_better: bool = True
    winner_tie_policy: str = "decoy"
    pseudocount: float = Field(default=1.0, ge=0.0)
    min_competitions: int = Field(default=100, ge=1)
    min_target_winners: int = Field(default=50, ge=1)
    min_decoy_winners: int = Field(default=10, ge=1)
    allowed_unit_levels: tuple[CalibrationUnitLevel, ...]
    allowed_claim_scopes: tuple[ClaimScope, ...]
    stratum_fallback: tuple[tuple[str, ...], ...]
    max_target_decoy_pool_ratio: float = Field(default=4.0, ge=1.0)
    qvalue_release_thresholds: tuple[float, ...] = (0.01, 0.05, 0.10)

    @model_validator(mode="after")
    def validate_predeclared_method(self) -> "CalibrationMethod":
        if not self.higher_score_is_better:
            raise ValueError("V9 calibration v1 supports higher-is-better scores only")
        if self.winner_tie_policy != "decoy":
            raise ValueError("Target-decoy ties must fail conservatively to decoy")
        allowed_fields = {
            "unit_level",
            "claim_scope",
            "polarity",
            "instrument_family",
            "source_tier",
        }
        if not self.stratum_fallback:
            raise ValueError("At least one predeclared stratum is required")
        for fields in self.stratum_fallback:
            if not fields or not set(fields) <= allowed_fields:
                raise ValueError(f"Unsupported calibration stratum: {fields}")
            if "unit_level" not in fields or "claim_scope" not in fields:
                raise ValueError(
                    "Every fallback level must preserve unit level and claim scope"
                )
        if len(set(self.stratum_fallback)) != len(self.stratum_fallback):
            raise ValueError("Calibration fallback strata must be unique")
        if any(not 0.0 < value <= 1.0 for value in self.qvalue_release_thresholds):
            raise ValueError("q-value thresholds must be in (0, 1]")
        return self


class CalibrationProject(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    project_id: str
    method_id: str
    runtime_snapshot_id: str
    library_snapshot_id: str
    candidate_space_id: str
    score_field: str
    target_pool_size: int = Field(ge=1)
    decoy_pool_size: int = Field(ge=1)
    decoy_design_id: str
    decoy_generation_seed: int
    method_frozen_before_data: bool
    independent_from_threshold_tuning: bool
    complete_target_decoy_competition: bool
    decoy_identity_overlap_count: int = Field(default=0, ge=0)
    instrument_domain: str
    notes: tuple[str, ...] = ()


class CompetitionPair(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    competition_id: str
    entity_id: str
    unit_level: CalibrationUnitLevel
    claim_scope: ClaimScope
    target_record_id: str
    target_score: float = Field(allow_inf_nan=False)
    decoy_score: float = Field(allow_inf_nan=False)
    polarity: str
    instrument_family: str
    source_tier: str


class WinnerCall(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    competition_id: str
    entity_id: str
    unit_level: CalibrationUnitLevel
    claim_scope: ClaimScope
    target_record_id: str
    role: CompetitionRole
    score: float = Field(allow_inf_nan=False)
    polarity: str
    instrument_family: str
    source_tier: str


class QValuePoint(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    score_threshold: float = Field(allow_inf_nan=False)
    target_winners: int = Field(ge=0)
    decoy_winners: int = Field(ge=0)
    raw_fdr: float = Field(ge=0.0, le=1.0)
    q_value: float = Field(ge=0.0, le=1.0)


class CalibrationReadiness(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    eligible: bool
    reason_codes: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    competition_count: int = Field(ge=0)
    target_winner_count: int = Field(ge=0)
    decoy_winner_count: int = Field(ge=0)
    decoy_scale: float | None = Field(default=None, gt=0.0)


class StratumCalibration(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    model_id: str
    fields: tuple[str, ...]
    values: tuple[str, ...]
    competition_count: int = Field(ge=1)
    target_winner_count: int = Field(ge=0)
    decoy_winner_count: int = Field(ge=0)
    decoy_scale: float = Field(gt=0.0)
    curve: tuple[QValuePoint, ...]


class CalibrationAssignment(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    competition_id: str
    entity_id: str
    target_record_id: str
    unit_level: CalibrationUnitLevel
    claim_scope: ClaimScope
    score: float = Field(allow_inf_nan=False)
    q_value: float = Field(ge=0.0, le=1.0)
    score_field: str
    method_id: str
    project_id: str
    stratum_model_id: str
    stratum_fields: tuple[str, ...]
    stratum_values: tuple[str, ...]
    transfer_scope: str = "project_local_exact_runtime_library_candidate_space"


class ProjectCalibrationArtifact(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = CALIBRATION_SCHEMA_VERSION
    status: ProjectCalibrationStatus
    method_id: str
    method_semantic_sha256: str
    project: CalibrationProject
    readiness: CalibrationReadiness
    strata: tuple[StratumCalibration, ...] = ()
    assignments: tuple[CalibrationAssignment, ...] = ()
    target_winners_unassigned: int = Field(default=0, ge=0)
    probability_available: bool = False
    local_fdr_available: bool = False
    qvalue_semantics: str = (
        "empirical project-local target-decoy q-value; not a posterior probability"
    )
    coverage_guarantee: bool = False
    production_ranking_changed: bool = False
    legacy_confidence_overwritten: bool = False

    @model_validator(mode="after")
    def enforce_fail_closed_artifact(self) -> "ProjectCalibrationArtifact":
        if self.status == ProjectCalibrationStatus.UNAVAILABLE and self.assignments:
            raise ValueError("Unavailable calibration artifacts cannot assign q-values")
        if self.probability_available or self.local_fdr_available:
            raise ValueError(
                "Target-decoy q-values alone cannot create probabilities or local FDR"
            )
        if self.production_ranking_changed or self.legacy_confidence_overwritten:
            raise ValueError("Calibration is a post-ranking, non-destructive layer")
        return self
