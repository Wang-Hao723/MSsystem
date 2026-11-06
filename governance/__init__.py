"""Module lifecycle governance and final claim compilation."""

from .module_registry import (
    ModuleDefinition,
    ModuleRegistry,
    audit_module_registry,
    load_module_registry,
)

__all__ = [
    "ModuleDefinition",
    "ModuleRegistry",
    "audit_module_registry",
    "load_module_registry",
]
