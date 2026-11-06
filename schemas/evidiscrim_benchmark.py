"""Contracts for identity-disjoint EviDiscrim-MS blind validation."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


EVIDISCRIM_BENCHMARK_SCHEMA_VERSION = "1.0.0"
BENCHMARK_METHODS = (
    "production_baseline",
    "whole_spectrum_fiora",
    "evidiscrim_full",
    "evidiscrim_no_uncertainty",
    "evidiscrim_no_abstention",
)
BenchmarkMethod = Literal[
    "production_baseline",
    "whole_spectrum_fiora",
    "evidiscrim_full",
    "evidiscrim_no_uncertainty",
    "evidiscrim_no_abstention",
]


class BenchmarkCandidate(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    candidate_id: str = Field(min_length=1)
    smiles: str = Field(min_length=1)
    formula: str = Field(min_length=1)
    inchikey: str = ""
    source: str = Field(min_length=1)


class AuthenticStandardSpectrum(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    standard_id: str = Field(min_length=1)
    blind_code: str = Field(min_length=1)
    replicate_id: str = Field(min_length=1)
    spectrum_id: str = Field(min_length=1)
    source_dataset: str = Field(min_length=1)
    source_record_id: str = Field(min_length=1)
    truth_provenance: str = Field(min_length=1)
    authentic_standard_lot: str = Field(min_length=1)
    true_candidate_id: str = Field(min_length=1)
    true_smiles: str = Field(min_length=1)
    true_inchikey: str = ""
    formula: str = Field(min_length=1)
    adduct: str = Field(min_length=1)
    polarity: Literal["positive", "negative"]
    precursor_mz: float = Field(gt=0)
    collision_energy: float = Field(gt=0)
    collision_energy_unit: str = Field(min_length=1)
    instrument_vendor: str = Field(min_length=1)
    instrument_model: str = Field(min_length=1)
    analyzer_type: str = Field(min_length=1)
    activation_type: str = Field(min_length=1)
    acquisition_batch: str = Field(min_length=1)
    run_order: int = Field(ge=1)
    peaks: tuple[tuple[float, float], ...] = Field(min_length=2)
    candidates: tuple[BenchmarkCandidate, ...] = Field(min_length=2)

    @field_validator("peaks")
    @classmethod
    def validate_peaks(
        cls, peaks: tuple[tuple[float, float], ...]
    ) -> tuple[tuple[float, float], ...]:
        if any(mz <= 0 or intensity <= 0 for mz, intensity in peaks):
            raise ValueError("all authentic-standard peaks must be positive")
        return peaks

    @model_validator(mode="after")
    def validate_candidate_ids(self) -> "AuthenticStandardSpectrum":
        candidate_ids = [candidate.candidate_id for candidate in self.candidates]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError("candidate IDs must be unique within a blind case")
        if self.true_candidate_id not in candidate_ids:
            raise ValueError("the authentic identity must occur in the candidate set")
        return self


class LockedMethodPrediction(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    blind_code: str = Field(min_length=1)
    method: BenchmarkMethod
    ranked_candidate_ids: tuple[str, ...] = ()
    abstained: bool = False
    released: bool = False
    release_score: float | None = None
    architecture_families: tuple[str, ...] = ()
    architecture_model_artifact_sha256: dict[str, str] = Field(
        default_factory=dict
    )
    prediction_lock_sha256: str = Field(min_length=64, max_length=64)

    @model_validator(mode="after")
    def enforce_prediction_consistency(self) -> "LockedMethodPrediction":
        if self.abstained and self.released:
            raise ValueError("an abstained result cannot be released")
        if not self.abstained and not self.ranked_candidate_ids:
            raise ValueError("a non-abstained result needs a ranked candidate")
        families = set(self.architecture_families)
        receipts = set(self.architecture_model_artifact_sha256)
        if families != receipts:
            raise ValueError(
                "every declared architecture family needs exactly one model-artifact hash"
            )
        if any(
            len(value) != 64
            or any(character not in "0123456789abcdef" for character in value.lower())
            for value in self.architecture_model_artifact_sha256.values()
        ):
            raise ValueError("model-artifact hashes must be 64 hexadecimal characters")
        return self


class ProjectFdrObservation(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    evaluation_unit_id: str = Field(min_length=1)
    method: BenchmarkMethod
    released: bool
    false_discovery: bool
    source_artifact_sha256: str = Field(
        min_length=64,
        max_length=64,
        description=(
            "SHA-256 of the shared locked project evaluation/truth cohort; "
            "paired methods must reference the same hash per unit"
        ),
    )

    @model_validator(mode="after")
    def enforce_fdr_consistency(self) -> "ProjectFdrObservation":
        if self.false_discovery and not self.released:
            raise ValueError("a non-released result cannot be a false discovery")
        return self


__all__ = [
    "AuthenticStandardSpectrum",
    "BENCHMARK_METHODS",
    "BenchmarkCandidate",
    "BenchmarkMethod",
    "EVIDISCRIM_BENCHMARK_SCHEMA_VERSION",
    "LockedMethodPrediction",
    "ProjectFdrObservation",
]
