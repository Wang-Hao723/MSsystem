"""Evidence objects with explicit applicability, provenance and permitted use."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator


class ModuleState(str, Enum):
    OFF = "off"
    SHADOW = "shadow"
    QUALIFIED = "qualified"
    PRODUCTION = "production"


class EvidenceApplicability(str, Enum):
    APPLICABLE = "applicable"
    NOT_APPLICABLE = "not_applicable"
    UNAVAILABLE = "unavailable"
    FAILED = "failed"
    CONFLICTED = "conflicted"


class EvidenceStatus(str, Enum):
    SUPPORTING = "supporting"
    INSUFFICIENT = "insufficient"
    CONTRADICTORY = "contradictory"
    SCOPE_LIMITED = "scope_limited"
    INFORMATIONAL = "informational"
    NOT_RUN = "not_run"
    UNAVAILABLE = "unavailable"
    FAILED = "failed"
    CONFLICTED = "conflicted"
    UNKNOWN = "unknown"


class EvidenceProvenance(BaseModel):
    model_config = ConfigDict(extra="allow", frozen=True)

    source_type: str = "legacy_result"
    source_id: str | None = None
    task_id: str | None = None
    spectrum_id: str | None = None
    library_spectrum_id: str | None = None
    library_accession: str | None = None
    method_id: str | None = None
    instrument_family: str | None = None
    runtime_snapshot_id: str | None = None


class EvidenceNode(BaseModel):
    """One auditable evidence family, never an implicit duplicate vote."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    node_id: str = Field(min_length=1)
    label: str = Field(min_length=1)
    status: EvidenceStatus
    direction: str
    applicability: EvidenceApplicability
    reason: str
    value: Any = None
    provenance: EvidenceProvenance
    independence_group: str
    independent_measurement: bool = False
    used_for: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    runtime_influence: bool = False
    candidate_conditioned: bool = False

    @model_validator(mode="after")
    def protect_candidate_conditioned_evidence(self) -> "EvidenceNode":
        if not self.candidate_conditioned:
            return self
        prohibited = {"ranking", "confidence", "claim_scope", "release_gate"}
        overlap = prohibited.intersection(self.used_for)
        if self.runtime_influence or overlap:
            raise ValueError(
                "Candidate-conditioned evidence must remain shadow-only and cannot "
                f"affect {sorted(overlap)}"
            )
        return self
