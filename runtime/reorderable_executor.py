"""A serial executor whose not-yet-started work can be reordered safely.

The standard :class:`concurrent.futures.ThreadPoolExecutor` deliberately keeps
its internal work queue private.  MSsystem's Full-library analyses must remain
strictly serial for memory safety, but operators also need to move an urgent
sample ahead of other *waiting* samples.  This module exposes only that narrow
operation; the running callable is never interrupted or moved.
"""

from __future__ import annotations

import threading
from concurrent.futures import Future
from dataclasses import dataclass
from typing import Any, Callable, Literal


QueueAction = Literal["front", "up", "down"]


@dataclass
class _WorkItem:
    key: str
    future: Future
    function: Callable[..., Any]
    args: tuple[Any, ...]
    kwargs: dict[str, Any]


class ReorderableSerialExecutor:
    """A Future-compatible, one-worker executor with an inspectable queue.

    ``queue_key`` defaults to the first positional argument, matching
    MSsystem's ``submit(_run_analysis, task_id)`` call.  Keys must be unique
    while running or waiting so UI actions cannot accidentally target two
    different jobs.
    """

    _max_workers = 1

    def __init__(self, *, thread_name_prefix: str = "serial") -> None:
        self._thread_name_prefix = str(thread_name_prefix or "serial")
        self._condition = threading.Condition()
        self._queue: list[_WorkItem] = []
        self._running_key: str | None = None
        self._worker: threading.Thread | None = None
        self._shutdown = False

    def submit(
        self,
        function: Callable[..., Any],
        /,
        *args: Any,
        queue_key: str | None = None,
        **kwargs: Any,
    ) -> Future:
        """Queue one callable and return a standard ``Future``."""
        if queue_key is None:
            if not args:
                raise ValueError("queue_key is required when no positional key exists")
            queue_key = str(args[0])
        key = str(queue_key).strip()
        if not key:
            raise ValueError("queue_key cannot be empty")

        future: Future = Future()
        item = _WorkItem(key, future, function, tuple(args), dict(kwargs))
        with self._condition:
            if self._shutdown:
                raise RuntimeError("cannot schedule new futures after shutdown")
            if key == self._running_key or any(row.key == key for row in self._queue):
                raise ValueError(f"queue_key is already active: {key}")
            self._queue.append(item)
            self._ensure_worker_locked()
            self._condition.notify_all()
        return future

    def snapshot(self) -> dict[str, Any]:
        """Return an atomic copy of the running key and waiting order."""
        with self._condition:
            return {
                "running": self._running_key,
                "queued": [item.key for item in self._queue],
            }

    def reorder(self, queue_key: str, action: QueueAction) -> dict[str, Any]:
        """Move a waiting item; the currently running item is immutable."""
        key = str(queue_key).strip()
        if action not in {"front", "up", "down"}:
            raise ValueError(f"unsupported queue action: {action}")
        with self._condition:
            index = next(
                (i for i, item in enumerate(self._queue) if item.key == key),
                None,
            )
            if index is None:
                raise KeyError(key)
            if action == "front":
                destination = 0
            elif action == "up":
                destination = max(0, index - 1)
            else:
                destination = min(len(self._queue) - 1, index + 1)
            if destination != index:
                item = self._queue.pop(index)
                self._queue.insert(destination, item)
            return {
                "running": self._running_key,
                "queued": [item.key for item in self._queue],
            }

    def cancel(self, queue_key: str) -> bool:
        """Remove a waiting item and cancel its Future, if still pending."""
        key = str(queue_key).strip()
        with self._condition:
            index = next(
                (i for i, item in enumerate(self._queue) if item.key == key),
                None,
            )
            if index is None:
                return False
            item = self._queue.pop(index)
            item.future.cancel()
            self._condition.notify_all()
            return True

    def shutdown(
        self,
        wait: bool = True,
        *,
        cancel_futures: bool = False,
    ) -> None:
        """Stop accepting work and optionally cancel jobs not yet started."""
        worker: threading.Thread | None
        with self._condition:
            self._shutdown = True
            if cancel_futures:
                pending, self._queue = self._queue, []
                for item in pending:
                    item.future.cancel()
            worker = self._worker
            self._condition.notify_all()
        if wait and worker is not None and worker is not threading.current_thread():
            worker.join()

    def _ensure_worker_locked(self) -> None:
        if self._worker is not None and self._worker.is_alive():
            return
        self._worker = threading.Thread(
            target=self._worker_loop,
            name=f"{self._thread_name_prefix}_0",
            daemon=True,
        )
        self._worker.start()

    def _worker_loop(self) -> None:
        while True:
            with self._condition:
                while not self._queue and not self._shutdown:
                    self._condition.wait()
                if not self._queue and self._shutdown:
                    return
                item = self._queue.pop(0)
                self._running_key = item.key

            if item.future.set_running_or_notify_cancel():
                try:
                    result = item.function(*item.args, **item.kwargs)
                except BaseException as exc:
                    item.future.set_exception(exc)
                else:
                    item.future.set_result(result)

            with self._condition:
                self._running_key = None
                self._condition.notify_all()

