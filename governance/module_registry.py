"""Machine-enforced deployment states for analytical modules."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator

from schemas.evidence import ModuleState


REGISTRY_SCHEMA_VERSION = "1.0.0"
SHADOW_SAFE_EFFECTS = {
    "audit", "display", "counterfactual", "acquisition_recommendation", "export"
}


class RuntimeAssertion(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    setting: str
    expected: Any


class ValidationArtifact(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    path: str
    sha256: str | None = None
    required: bool = True
    interpretation: str


class ModuleDefinition(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    module_id: str
    display_name: str
    state: ModuleState
    scientific_role: str
    allowed_effects: tuple[str, ...]
    forbidden_effects: tuple[str, ...]
    applicability_requirements: tuple[str, ...] = ()
    runtime_assertions: tuple[RuntimeAssertion, ...] = ()
    validation_artifacts: tuple[ValidationArtifact, ...] = ()
    promotion_requirements: tuple[str, ...] = ()
    demotion_triggers: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @model_validator(mode="after")
    def enforce_state_contract(self) -> "ModuleDefinition":
        allowed = set(self.allowed_effects)
        forbidden = set(self.forbidden_effects)
        if allowed & forbidden:
            raise ValueError(
                f"Module {self.module_id} both allows and forbids {allowed & forbidden}"
            )
        if self.state == ModuleState.OFF and allowed:
            raise ValueError("Off modules cannot declare runtime effects")
        if self.state == ModuleState.SHADOW and not allowed <= SHADOW_SAFE_EFFECTS:
            raise ValueError(
                f"Shadow module {self.module_id} has unsafe effects: "
                f"{sorted(allowed - SHADOW_SAFE_EFFECTS)}"
            )
        if self.state in {ModuleState.QUALIFIED, ModuleState.PRODUCTION} and not (
            self.validation_artifacts or self.promotion_requirements
        ):
            raise ValueError(
                f"{self.state.value} module {self.module_id} lacks validation governance"
            )
        return self


class ModuleRegistry(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = REGISTRY_SCHEMA_VERSION
    registry_id: str
    modules: tuple[ModuleDefinition, ...]
    global_invariants: tuple[str, ...]

    @model_validator(mode="after")
    def unique_modules(self) -> "ModuleRegistry":
        module_ids = [module.module_id for module in self.modules]
        if len(module_ids) != len(set(module_ids)):
            raise ValueError("Module registry IDs must be unique")
        return self


def load_module_registry(path: Path) -> ModuleRegistry:
    return ModuleRegistry.model_validate_json(path.read_text(encoding="utf-8"))


def registry_semantic_sha256(registry: ModuleRegistry) -> str:
    payload = json.dumps(
        registry.model_dump(mode="json"),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def audit_module_registry(
    registry: ModuleRegistry,
    *,
    project_root: Path,
    runtime_settings: Any,
) -> dict:
    violations = []
    modules = []
    for module in registry.modules:
        assertions = []
        for assertion in module.runtime_assertions:
            actual = getattr(runtime_settings, assertion.setting, None)
            passed = actual == assertion.expected
            assertions.append(
                {
                    "setting": assertion.setting,
                    "expected": assertion.expected,
                    "actual": actual,
                    "passed": passed,
                }
            )
            if not passed:
                violations.append(
                    f"{module.module_id}:runtime_setting_mismatch:{assertion.setting}"
                )
        artifacts = []
        for artifact in module.validation_artifacts:
            path = project_root / artifact.path
            exists = path.is_file()
            actual_hash = hashlib.sha256(path.read_bytes()).hexdigest() if exists else None
            hash_ok = artifact.sha256 is None or actual_hash == artifact.sha256
            passed = (exists or not artifact.required) and hash_ok
            artifacts.append(
                {
                    **artifact.model_dump(mode="json"),
                    "exists": exists,
                    "actual_sha256": actual_hash,
                    "passed": passed,
                }
            )
            if not passed:
                violations.append(f"{module.module_id}:validation_artifact_failed:{artifact.path}")
        modules.append(
            {
                "module_id": module.module_id,
                "state": module.state.value,
                "allowed_effects": module.allowed_effects,
                "runtime_assertions": assertions,
                "validation_artifacts": artifacts,
                "passed": all(item["passed"] for item in assertions + artifacts),
            }
        )
    return {
        "schema_version": REGISTRY_SCHEMA_VERSION,
        "registry_id": registry.registry_id,
        "registry_semantic_sha256": registry_semantic_sha256(registry),
        "passed": not violations,
        "violation_count": len(violations),
        "violations": violations,
        "modules": modules,
        "global_invariants": registry.global_invariants,
    }
