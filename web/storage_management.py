"""FastAPI installation boundary for conservative storage maintenance."""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from typing import Any, Callable, Literal

from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, ConfigDict

from runtime.storage_maintenance import (
    StorageMaintenanceError,
    StorageMaintenanceManager,
)


logger = logging.getLogger(__name__)


class StoragePolicyRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)

    recycle_automatic_cleanup_enabled: bool
    recycle_retention_days: Literal[7, 30, 90]
    cache_automatic_cleanup_enabled: bool
    cache_retention_days: Literal[1, 7, 30]
    confirm_auto_delete: bool = False


class StoragePreviewRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)

    target: Literal["recycle", "cache"]


class StorageCleanRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)

    token: str
    confirmed: bool


def _raise_http(exc: StorageMaintenanceError) -> None:
    raise HTTPException(status_code=exc.status_code, detail=exc.detail) from exc


async def _automatic_worker(
    manager: StorageMaintenanceManager,
    stop_event: asyncio.Event,
    *,
    interval_seconds: float = 3600.0,
) -> None:
    """Wait one full hour before the first check; the manager enforces daily cadence."""

    while not stop_event.is_set():
        try:
            await asyncio.wait_for(stop_event.wait(), timeout=interval_seconds)
        except asyncio.TimeoutError:
            try:
                await asyncio.to_thread(manager.run_auto_cleanup)
            except Exception:
                # Automatic maintenance must never terminate the web service.
                logger.exception("automatic storage maintenance failed")
                continue


def install_storage_management(
    app: FastAPI,
    *,
    output_dir: Path,
    tasks_snapshot: Callable[[], Any],
    recycle_inventory_unlocked: Callable[[], dict[str, Any]],
    recycle_lock: Any,
    cache_lock: Any = None,
    append_audit: Callable[[dict[str, Any]], None] | None = None,
) -> StorageMaintenanceManager:
    """Install the storage API and delayed automatic-maintenance worker.

    ``tasks_snapshot`` must return a stable shallow copy of the task registry.
    ``recycle_inventory_unlocked`` is invoked only while ``recycle_lock`` is
    held. ``cache_lock`` must be the report/export lock and expose non-blocking
    ``acquire`` and ``release``. ``append_audit`` receives completed cleanup
    audit dictionaries; cleanup receipts remain durable even if it fails.
    """

    if getattr(app.state, "storage_maintenance_manager", None) is not None:
        raise RuntimeError("storage management is already installed")
    manager = StorageMaintenanceManager(
        output_dir=Path(output_dir),
        tasks_snapshot=tasks_snapshot,
        recycle_inventory_unlocked=recycle_inventory_unlocked,
        recycle_lock=recycle_lock,
        cache_lock=cache_lock,
        append_audit=append_audit,
    )
    app.state.storage_maintenance_manager = manager
    app.state.storage_maintenance_stop_event = None
    app.state.storage_maintenance_worker = None

    @app.get("/api/storage-management")
    async def storage_management_overview():
        try:
            return await asyncio.to_thread(manager.overview)
        except StorageMaintenanceError as exc:
            _raise_http(exc)

    @app.get("/api/storage-management/history")
    async def storage_management_history(limit: int = Query(default=10, ge=1, le=20)):
        try:
            return await asyncio.to_thread(manager.cleanup_history, limit)
        except StorageMaintenanceError as exc:
            _raise_http(exc)

    @app.put("/api/storage-management/policy")
    async def update_storage_management_policy(payload: StoragePolicyRequest):
        try:
            return await asyncio.to_thread(
                manager.save_policy,
                recycle_enabled=payload.recycle_automatic_cleanup_enabled,
                recycle_retention_days=payload.recycle_retention_days,
                cache_enabled=payload.cache_automatic_cleanup_enabled,
                cache_retention_days=payload.cache_retention_days,
                confirm_auto_delete=payload.confirm_auto_delete,
            )
        except StorageMaintenanceError as exc:
            _raise_http(exc)

    @app.post("/api/storage-management/preview")
    async def preview_storage_cleanup(payload: StoragePreviewRequest):
        try:
            return await asyncio.to_thread(manager.preview, payload.target)
        except StorageMaintenanceError as exc:
            _raise_http(exc)

    @app.post("/api/storage-management/clean")
    async def clean_storage(payload: StorageCleanRequest):
        try:
            return await asyncio.to_thread(
                manager.clean,
                payload.token,
                confirmed=payload.confirmed,
            )
        except StorageMaintenanceError as exc:
            _raise_http(exc)

    async def startup() -> None:
        stop_event = asyncio.Event()
        app.state.storage_maintenance_stop_event = stop_event
        app.state.storage_maintenance_worker = asyncio.create_task(
            _automatic_worker(manager, stop_event),
            name="mssystem-storage-maintenance",
        )

    async def shutdown() -> None:
        stop_event = app.state.storage_maintenance_stop_event
        worker = app.state.storage_maintenance_worker
        if stop_event is not None:
            stop_event.set()
        if worker is not None:
            try:
                await worker
            except asyncio.CancelledError:
                # A caller outside this shutdown hook may cancel the worker.
                # The normal shutdown path waits for an in-flight to_thread
                # cleanup so its pre-delete receipt can be finalized.
                raise
        app.state.storage_maintenance_worker = None

    app.add_event_handler("startup", startup)
    app.add_event_handler("shutdown", shutdown)
    return manager


__all__ = ["install_storage_management"]
