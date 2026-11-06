"""Lazy MSClassify inference for spectra without a library hit.

V6 previously imported ``models.msclass_wrapper`` even though the ``models``
directory is intentionally excluded from the project checkout.  The import
failure was swallowed by the result serializer, so unmatched spectra silently
lost all class-level output.  This wrapper uses the installed MSClassify model
and checkpoint directly and keeps the optional dependency failure explicit.
"""

from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

_DEFAULT_ROOT = Path("/Users/zizi/Desktop/MSClassify")
_ROOT = Path(os.environ.get("MSCLASSIFY_ROOT", _DEFAULT_ROOT))
_CHECKPOINT = Path(
    os.environ.get(
        "MSCLASSIFY_CHECKPOINT",
        _ROOT / "checkpoints" / "best_model_small.pt",
    )
)
_WRAPPER: Optional["MSClassWrapper"] = None
_LOAD_ATTEMPTED = False


def _class_label(value: str) -> str:
    """Convert checkpoint labels such as PhenolicAcid to phenolic_acid."""
    return re.sub(r"([A-Z])", r"_\1", str(value)).lower().lstrip("_")


class MSClassWrapper:
    """Small inference adapter around the separately installed MSClassify."""

    def __init__(self, root: Path = _ROOT, checkpoint: Path = _CHECKPOINT):
        root = Path(root).resolve()
        checkpoint = Path(checkpoint).resolve()
        if not (root / "models" / "classifier.py").exists():
            raise FileNotFoundError(f"MSClassify model source not found: {root}")
        if not checkpoint.exists():
            raise FileNotFoundError(f"MSClassify checkpoint not found: {checkpoint}")

        import torch
        import torch.nn.functional as F

        # MSClassify is a separate local project.  Import its package without
        # making Mssystem_v6 depend on the historical V5 models directory.
        root_s = str(root)
        inserted = root_s not in sys.path
        if inserted:
            sys.path.insert(0, root_s)
        try:
            from models.classifier import MSClassifier
        finally:
            if inserted and sys.path and sys.path[0] == root_s:
                sys.path.pop(0)

        if torch.cuda.is_available():
            self.device = torch.device("cuda")
        elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
            self.device = torch.device("mps")
        else:
            self.device = torch.device("cpu")

        ckpt = torch.load(checkpoint, map_location="cpu", weights_only=False)
        state = ckpt.get("model_state_dict", ckpt)
        self.class_names = list(ckpt.get("class_names") or [])
        if not self.class_names:
            raise ValueError("MSClassify checkpoint has no class_names")
        sub_weight = state.get("subclass_head.net.4.weight")
        n_subclasses = int(sub_weight.shape[0]) if sub_weight is not None else len(self.class_names) * 10

        self.model = MSClassifier(
            input_dim=1000,
            num_classes=len(self.class_names),
            num_subclasses=n_subclasses,
            latent_dim=256,
            dropout=0.4,
            encoder_type="hybrid",
        )
        self.model.load_state_dict(state)
        self.model.to(self.device)
        self.model.eval()
        self._labels = [_class_label(x) for x in self.class_names]
        raw_class_acc = ckpt.get("per_class_acc") or {}
        self.class_acc = {
            _class_label(key): float(value)
            for key, value in raw_class_acc.items()
            if value is not None
        }
        self._torch = torch
        self._functional = F
        logger.info(
            "MSClassify loaded for unmatched spectra: %d classes, device=%s",
            len(self.class_names), self.device,
        )

    @staticmethod
    def _usable_peaks(peaks):
        result = []
        for peak in (() if peaks is None else peaks):
            try:
                mz, intensity = float(peak[0]), float(peak[1])
            except (TypeError, ValueError, IndexError):
                continue
            if np.isfinite(mz) and np.isfinite(intensity) and 0 <= mz < 1000 and intensity > 0:
                result.append((mz, intensity))
        return result

    @staticmethod
    def _bin_spectrum(peaks) -> np.ndarray:
        # A peak list is unordered. Select the 200 strongest usable observations
        # deterministically before binning; scale first to avoid overflow.
        usable = sorted(MSClassWrapper._usable_peaks(peaks), key=lambda p: (-p[1], p[0]))[:200]
        binned = np.zeros(1000, dtype=np.float64)
        if usable:
            scale = usable[0][1]
            for mz, intensity in usable:
                binned[int(mz)] += intensity / scale
            binned /= binned.max()
        return binned.astype(np.float32)

    def predict(self, peaks, precursor_mz=None):
        """Return a broad class prediction, or None for unusable spectra."""
        peaks = self._usable_peaks(peaks)
        if len(peaks) < 3:
            return None
        binned = self._bin_spectrum(peaks)
        if not np.any(binned):
            return None
        tensor = self._torch.frombuffer(
            bytearray(binned.tobytes()), dtype=self._torch.float32
        ).clone().unsqueeze(0).to(self.device)
        with self._torch.no_grad():
            logits = self.model(tensor)["class_logits"]
            probs = self._functional.softmax(logits, dim=-1).squeeze(0).cpu()
        if not bool(self._torch.isfinite(probs).all()) or len(probs) != len(self._labels):
            logger.warning("MSClassify rejected nonfinite or incompatible model output")
            return None
        pred = int(probs.argmax())
        top_idx = self._torch.topk(probs, min(3, len(self._labels))).indices.tolist()
        return {
            "class": self._labels[pred],
            "confidence": float(probs[pred]),
            "top3": [(self._labels[i], float(probs[i])) for i in top_idx],
            "probs": {self._labels[i]: float(probs[i]) for i in range(len(self._labels))},
        }


def get_msclass_wrapper() -> Optional[MSClassWrapper]:
    """Load once; return None when the optional external model is unavailable."""
    global _WRAPPER, _LOAD_ATTEMPTED
    if _WRAPPER is not None:
        return _WRAPPER
    if _LOAD_ATTEMPTED:
        return None
    _LOAD_ATTEMPTED = True
    try:
        _WRAPPER = MSClassWrapper()
    except Exception as exc:
        logger.warning("MSClassify unavailable; unmatched classification disabled: %s", exc)
        _WRAPPER = None
    return _WRAPPER
