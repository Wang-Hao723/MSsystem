"""Calibrated claim decisions, kept separate from retrieval ranking."""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, model_validator


class ClaimScope(str, Enum):
    AUTHENTIC_STANDARD = "authentic_standard"
    CONNECTIVITY = "2d_connectivity"
    ISOMER_GROUP = "isomer_group"
    MOLECULAR_FORMULA = "molecular_formula"
    CHEMICAL_CLASS = "chemical_class"
    ANALOG = "analog"
    UNRESOLVED = "unresolved"
    ABSTAIN = "abstain"


class CalibrationStatus(str, Enum):
    CALIBRATED = "calibrated"
    UNAVAILABLE = "unavailable"
    OUT_OF_DOMAIN = "out_of_domain"
    FAILED = "failed"


class CalibrationVector(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    status: CalibrationStatus = CalibrationStatus.UNAVAILABLE
    method_id: str | None = None
    project_id: str | None = None
    stratum_model_id: str | None = None
    unit_level: str | None = None
    score_field: str | None = None
    transfer_scope: str | None = None
    q_connectivity: float | None = Field(default=None, ge=0.0, le=1.0)
    q_isomer_group: float | None = Field(default=None, ge=0.0, le=1.0)
    q_formula: float | None = Field(default=None, ge=0.0, le=1.0)
    q_class: float | None = Field(default=None, ge=0.0, le=1.0)
    probability: float | None = Field(default=None, ge=0.0, le=1.0)
    local_fdr: float | None = Field(default=None, ge=0.0, le=1.0)
    reason: str | None = None

    @model_validator(mode="after")
    def prohibit_unearned_q_values(self) -> "CalibrationVector":
        values = (
            self.q_connectivity,
            self.q_isomer_group,
            self.q_formula,
            self.q_class,
            self.probability,
            self.local_fdr,
        )
        if self.status != CalibrationStatus.CALIBRATED and any(
            value is not None for value in values
        ):
            raise ValueError("Uncalibrated records cannot expose q-values or probabilities")
        if self.status == CalibrationStatus.CALIBRATED:
            if not any(value is not None for value in values):
                raise ValueError("Calibrated records must expose a calibrated quantity")
            if not all((self.method_id, self.project_id, self.stratum_model_id)):
                raise ValueError(
                    "Calibrated quantities require method, project and stratum provenance"
                )
        return self


class ClaimDecision(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    policy_id: str
    status: str
    scope: ClaimScope
    released: bool = False
    eligible_before_replication: bool = False
    requires_independent_replication: bool = True
    reason_codes: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    independent_support_groups: tuple[str, ...] = ()
    independent_injection_count: int = Field(default=0, ge=0)
    ambiguity_identity_keys: tuple[str, ...] = ()
    coverage_guarantee: bool = False
    ranking_changed: bool = False
    confidence_field_overwritten: bool = False
    calibration: CalibrationVector = Field(default_factory=CalibrationVector)

    @model_validator(mode="after")
    def require_calibration_for_released_connectivity(self) -> "ClaimDecision":
        if (
            self.released
            and self.scope == ClaimScope.CONNECTIVITY
            and self.calibration.status != CalibrationStatus.CALIBRATED
        ):
            raise ValueError(
                "A released connectivity claim requires project-level calibration"
            )
        return self
