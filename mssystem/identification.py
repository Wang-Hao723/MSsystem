from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import numpy as np
from matchms import Spectrum
try:
    from matchms.similarity import ModifiedCosineGreedy as ModifiedCosine
except ImportError:
    from matchms.similarity import ModifiedCosine

from .classification import FingerprintClassifier, RuleBasedClassifier

logger = logging.getLogger(__name__)


class IdentificationClassifier:
    """
    MS/MS spectrum identification against a spectral library.

    Optimised for large libraries (660k+ records):
      - Library records are pre-indexed by precursor_mz into Da-bucketed buckets.
        At query time only candidates within the m/z tolerance window are fetched,
        avoiding a full-table cosine scan.
      - Cosine scores are computed with ModifiedCosineGreedy (handles missing
        peaks and isotope patterns better than plain CosineGreedy).
    """

    # Number of Dalton per bucket for the m/z index.
    # Smaller bucket = faster lookup but more dict overhead. 0.1 Da is good for HRMS.
    _MZ_BUCKET_DALTON = 0.1

    def __init__(
        self,
        library_records: List[Dict[str, Any]],
        ml_classifier: FingerprintClassifier | None = None,
        rule_classifier: RuleBasedClassifier | None = None,
        cosine_threshold: float = 0.7,
        mz_tolerance_da: float = 0.5,
        mz_tolerance_ppm: float = 10.0,
    ) -> None:
        self.library_records = library_records
        self.ml_classifier = ml_classifier or FingerprintClassifier()
        self.rule_classifier = rule_classifier or RuleBasedClassifier()
        self.cosine_threshold = cosine_threshold

        # m/z tolerance: use whichever is larger (Da or ppm)
        self.mz_tolerance_da = mz_tolerance_da
        self.mz_tolerance_ppm = mz_tolerance_ppm

        # Pre-index library by precursor_mz for fast candidate retrieval
        self._mz_index: Dict[int, List[int]] = {}
        self._build_mz_index()

        # Shared ModifiedCosineGreedy scorer
        self._cosine = ModifiedCosine(tolerance=0.1)

    # ── m/z index ──────────────────────────────────────────────────────────────

    def _build_mz_index(self) -> None:
        """Build a Da-bucketed index of record indices keyed by precursor_mz."""
        bucket_size = self._MZ_BUCKET_DALTON
        self._mz_index.clear()
        min_mz = float("inf")
        max_mz = float("-inf")

        for idx, rec in enumerate(self.library_records):
            pmz = self._safe_float(rec.get("precursor_mz"))
            if pmz is None or pmz <= 0:
                continue
            bucket = int(pmz // bucket_size)
            self._mz_index.setdefault(bucket, []).append(idx)
            if pmz < min_mz:
                min_mz = pmz
            if pmz > max_mz:
                max_mz = pmz

        # Sort each bucket's index list by precursor_mz for binary search
        for bucket in self._mz_index:
            self._mz_index[bucket].sort(
                key=lambda i: self._safe_float(
                    self.library_records[i].get("precursor_mz")
                )
                or 0.0
            )

        logger.info(
            "m/z index built: %d buckets, mz range [%.2f, %.2f], %d records indexed",
            len(self._mz_index),
            min_mz,
            max_mz,
            sum(len(v) for v in self._mz_index.values()),
        )

    def _mz_candidates(self, query_mz: float) -> List[int]:
        """Return record indices within the m/z tolerance window."""
        da_lo, da_hi = self._mz_window_da(query_mz)

        # Collect buckets that overlap the window
        bucket_size = self._MZ_BUCKET_DALTON
        start_bucket = int(da_lo // bucket_size)
        end_bucket = int(da_hi // bucket_size)

        candidates: List[int] = []
        seen: set[int] = set()

        for bucket in range(start_bucket, end_bucket + 1):
            for idx in self._mz_index.get(bucket, []):
                if idx in seen:
                    continue
                seen.add(idx)
                pmz = self._safe_float(self.library_records[idx].get("precursor_mz"))
                if pmz is not None and da_lo <= pmz <= da_hi:
                    candidates.append(idx)

        return candidates

    def _mz_window_da(self, query_mz: float) -> tuple[float, float]:
        """Compute the m/z search window (lower, upper) in Da."""
        if query_mz <= 0:
            return (0.0, float("inf"))

        # Da-based window
        da_lo = query_mz - self.mz_tolerance_da
        da_hi = query_mz + self.mz_tolerance_da

        # ppm-based window (whichever is wider)
        ppm_delta = query_mz * self.mz_tolerance_ppm * 1e-6
        ppm_lo = query_mz - ppm_delta
        ppm_hi = query_mz + ppm_delta

        # Use the union (larger window)
        return (min(da_lo, ppm_lo), max(da_hi, ppm_hi))

    # ── Identification ─────────────────────────────────────────────────────────

    def identify(
        self, query_record: Dict[str, Any], top_k: int = 5
    ) -> List[Dict[str, Any]]:
        query_mz = self._safe_float(query_record.get("precursor_mz"))
        if query_mz is None or query_mz <= 0:
            logger.warning("Query record has no valid precursor_mz — returning empty")
            return []

        # Step 1: m/z pre-filtering (fast index scan, not full table)
        candidate_indices = self._mz_candidates(query_mz)
        if not candidate_indices:
            return []

        # Step 2: Build query matchms Spectrum (once)
        query_spec = self._to_matchms_spectrum(query_record)
        # Verify precursor_mz is in metadata (required by ModifiedCosine)
        # This can fail if precursor_mz is NaN, empty string, or other invalid values
        if query_spec.metadata.get("precursor_mz") is None:
            logger.warning(
                "Query spectrum missing precursor_mz in metadata after conversion. "
                "Original value: %s (type: %s)",
                query_record.get("precursor_mz"),
                type(query_record.get("precursor_mz"))
            )
            return []

        # Step 3: Cosine scoring (only on candidates)
        candidates: List[Dict[str, Any]] = []
        candidates_mass_only: List[Dict[str, Any]] = []  # For records without peaks

        for idx in candidate_indices:
            rec = self.library_records[idx]
            rec_peaks = rec.get("peaks", [])

            # Check if library record has fragment peaks
            has_lib_peaks = rec_peaks and len(rec_peaks) > 0

            # Build lib spectrum on-demand
            lib_spec = self._to_matchms_spectrum(rec)

            # ModifiedCosine requires precursor_mz on both sides
            lib_pmz = lib_spec.metadata.get("precursor_mz")
            if lib_pmz is None:
                logger.debug("Skipping library record without precursor_mz: %s", rec.get("name"))
                continue

            # ModifiedCosine score (handles neutral losses & fragment shifts)
            cosine_result = self._cosine.pair(query_spec, lib_spec)
            cosine_score = self._extract_score(cosine_result)

            # Precursor m/z match (ppm-based, high-resolution MS standard)
            mass_score = self._mass_match_score(query_record, rec)

            # If no library peaks but acceptable mass match, still consider as candidate
            # with lower confidence (mass-only match)
            # Allow wider mass tolerance for mass-only matches (up to 100 ppm)
            if not has_lib_peaks and mass_score > 0.01:  # Allows up to ~100 ppm
                logger.debug("Mass-only candidate: %s, mass_score=%.4f", rec.get("name"), mass_score)
                candidates_mass_only.append({
                    "compound_name": rec.get("name"),
                    "inchikey": rec.get("inchikey"),
                    "smiles": rec.get("smiles"),
                    "inchi": rec.get("inchi"),
                    "formula": rec.get("formula"),
                    "exact_mass": rec.get("exact_mass"),
                    "compound_class": rec.get("compound_class"),
                    "classification_rule": rec.get("compound_class") or rec.get("classification") or "Unknown",
                    "classification_ml": "Unknown",
                    "cosine_score": 0.0,  # No peaks to compare
                    "modified_cosine": 0.0,
                    "mass_score": mass_score,
                    "fragment_score": 0.0,
                    "intensity_score": 0.0,
                    "mass_ppm": round(abs(query_mz - rec.get("precursor_mz", 0)) / rec.get("precursor_mz", 1) * 1e6, 2) if rec.get("precursor_mz") else 999.0,
                    "isotope_score": 0.0,
                    "final_score": mass_score * 0.6,  # Lower score for mass-only match
                    "matched_precursor_mz": rec.get("precursor_mz"),
                    "matched_peaks": [],
                    "match_type": "mass_only",  # Tag as mass-only match
                })
                continue

            if cosine_score < self.cosine_threshold:
                continue

            # Fragment ion match rate: shared peaks / query peaks
            fragment_score = self._fragment_match_score(query_record, rec)

            # Spectral intensity correlation (Pearson on matched peaks)
            intensity_score = self._intensity_correlation_score(query_record, rec)

            # Composite score per analytical chemistry practice:
            #   - Spectral similarity (cosine) is the primary identifier
            #   - Precursor mass accuracy gates high-confidence matches
            #   - Fragment match rate confirms structural similarity
            #   - Intensity correlation validates peak patterns
            # Level 1 (mass < 5 ppm): cosine-dominant scoring
            # Level 2 (mass 5-20 ppm): penalty applied
            # Level 3 (mass > 20 ppm): significant penalty
            q_mz = self._safe_float(query_record.get("precursor_mz"))
            r_mz = self._safe_float(rec.get("precursor_mz"))
            ppm = abs(q_mz - r_mz) / r_mz * 1e6 if (q_mz and r_mz and r_mz > 0) else 999.0

            if ppm <= 5.0:
                # High mass accuracy: trust cosine heavily
                final_score = (
                    0.50 * cosine_score
                    + 0.20 * mass_score
                    + 0.20 * fragment_score
                    + 0.10 * intensity_score
                )
            elif ppm <= 20.0:
                # Medium mass accuracy: balance cosine with fragment evidence
                final_score = (
                    0.40 * cosine_score
                    + 0.20 * mass_score
                    + 0.25 * fragment_score
                    + 0.15 * intensity_score
                )
            else:
                # Low mass accuracy: rely more on fragment evidence
                final_score = (
                    0.30 * cosine_score
                    + 0.10 * mass_score
                    + 0.40 * fragment_score
                    + 0.20 * intensity_score
                )

            # Rule-based classification, fallback to DB compound_class
            rule_cls = self.rule_classifier.predict(
                smiles=rec.get("smiles"), inchi=rec.get("inchi")
            )
            if rule_cls == "Unknown":
                db_cls = rec.get("compound_class") or rec.get("classification") or ""
                if db_cls and db_cls not in ("Unknown", "All_Spectra", "N/A", "NA", "None", "CLASS1", ""):
                    rule_cls = db_cls

            candidates.append(
                {
                    "compound_name": rec.get("name"),
                    "inchikey": rec.get("inchikey"),
                    "smiles": rec.get("smiles"),
                    "inchi": rec.get("inchi"),
                    "formula": rec.get("formula"),
                    "exact_mass": rec.get("exact_mass"),
                    "compound_class": rec.get("compound_class"),
                    "classification_rule": rule_cls,
                    "classification_ml": self.ml_classifier.predict(
                        rec.get("smiles", "")
                    )
                    if rec.get("smiles")
                    else "Unknown",
                    "cosine_score": cosine_score,
                    "modified_cosine": cosine_score,
                    "mass_score": mass_score,
                    "fragment_score": fragment_score,
                    "intensity_score": intensity_score,
                    "mass_ppm": round(ppm, 2),
                    "isotope_score": intensity_score,
                    "final_score": final_score,
                    "matched_precursor_mz": rec.get("precursor_mz"),
                    "matched_peaks": rec.get("peaks", []),
                    "match_type": "full",  # Tag as full match (with peaks)
                }
            )

        # Merge mass-only candidates at the end (lower priority)
        candidates.extend(candidates_mass_only)

        candidates.sort(key=lambda x: x["final_score"], reverse=True)

        # Deduplicate by compound name - keep only the best match per compound
        seen_names = set()
        unique_candidates = []
        for c in candidates:
            name = c.get("compound_name", "")
            if name not in seen_names:
                seen_names.add(name)
                unique_candidates.append(c)
                if len(unique_candidates) >= top_k:
                    break

        return unique_candidates

    # ── Helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _safe_float(val: Any) -> Optional[float]:
        if val is None:
            return None
        try:
            f = float(val)
            return f if np.isfinite(f) else None
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _to_matchms_spectrum(record: Dict[str, Any]) -> Spectrum:
        peaks = record.get("peaks", [])
        if not peaks:
            return Spectrum(mz=np.array([]), intensities=np.array([]), metadata={})

        mz = np.array([float(p[0]) for p in peaks], dtype=float)
        intensities = np.array([float(p[1]) for p in peaks], dtype=float)

        # Normalize intensities to 0-1 range (required by matchms)
        max_int = intensities.max()
        if max_int > 0:
            intensities = intensities / max_int

        # matchms requires mz to be sorted in ascending order
        order = np.argsort(mz)
        mz = mz[order]
        intensities = intensities[order]

        metadata: Dict[str, Any] = {}
        pmz = record.get("precursor_mz")
        if pmz is not None:
            try:
                pmz_float = float(pmz)
                # Only set precursor_mz if it's a valid finite number
                if np.isfinite(pmz_float):
                    metadata["precursor_mz"] = pmz_float
            except (TypeError, ValueError):
                pass
        rt = record.get("retention_time")
        if rt is not None:
            try:
                rt_float = float(rt)
                if np.isfinite(rt_float):
                    metadata["retention_time"] = rt_float
            except (TypeError, ValueError):
                pass

        return Spectrum(mz=mz, intensities=intensities, metadata=metadata)

    @staticmethod
    def _extract_score(result: Any) -> float:
        """Extract a scalar score from matchms similarity results.

        Handles: None, numpy 0-d structured array {score, matches} (ModifiedCosine),
        numpy 1-d array, (score, count) tuple/list, and objects with a .score attribute.
        """
        if result is None:
            return 0.0
        score = result
        # matchms ModifiedCosine/CosineGreedy: 0-d structured array [('score', f8), ('matches', i8)]
        if isinstance(score, np.ndarray) and score.ndim == 0:
            score = score.tolist()[0]
        # matchms array results
        elif isinstance(score, np.ndarray) and score.size > 0:
            arr = score.tolist()
            score = arr[0] if isinstance(arr, (list, tuple)) else float(arr)
        # (score, count) tuple/list
        elif isinstance(score, (tuple, list)) and len(score) > 0:
            score = score[0]
        # matchms Score objects
        elif hasattr(score, "score"):
            score = score.score
        try:
            return float(score)
        except (TypeError, ValueError):
            return 0.0

    @staticmethod
    def _mass_match_score(query: Dict[str, Any], ref: Dict[str, Any]) -> float:
        """Precursor m/z match score based on ppm error.

        Analytical chemistry standard:
          - ≤ 2 ppm: excellent (TOF/Orbitrap high-res)
          - ≤ 5 ppm: good (Q-TOF)
          - ≤ 10 ppm: acceptable
          - ≤ 20 ppm: tolerable (low-res)
          - > 20 ppm: unlikely same compound
        """
        q = IdentificationClassifier._safe_float(query.get("precursor_mz"))
        r = IdentificationClassifier._safe_float(ref.get("precursor_mz"))
        if q is None or r is None or r <= 0:
            return 0.0
        ppm = abs(q - r) / r * 1e6
        # Gaussian-like decay centered at 0 ppm, σ = 10 ppm
        return float(np.exp(-0.5 * (ppm / 10.0) ** 2))

    @staticmethod
    def _fragment_match_score(query: Dict[str, Any], ref: Dict[str, Any],
                              tol: float = 0.02) -> float:
        """Fragment ion match score: shared peaks / total query peaks.

        Uses ±0.02 Da tolerance for fragment matching (high-resolution MS standard).
        Computes bidirectional matches: how many query fragment ions are
        explained by the library. Returns the larger count normalized by query
        size to avoid penalizing libraries with extra peaks.
        """
        q_peaks = query.get("peaks", [])
        r_peaks = ref.get("peaks", [])
        if not q_peaks or not r_peaks:
            return 0.0

        matched = 0
        used_ref = set()
        for qp in q_peaks:
            qmz = qp[0] if len(qp) >= 2 else 0.0
            for j, rp in enumerate(r_peaks):
                if j in used_ref:
                    continue
                rmz = rp[0] if len(rp) >= 2 else 0.0
                if abs(qmz - rmz) <= tol:
                    matched += 1
                    used_ref.add(j)
                    break

        # Bidirectional: also check library peaks matched by query
        matched_rev = 0
        used_q = set()
        for rp in r_peaks:
            rmz = rp[0] if len(rp) >= 2 else 0.0
            for j, qp in enumerate(q_peaks):
                if j in used_q:
                    continue
                qmz = qp[0] if len(qp) >= 2 else 0.0
                if abs(qmz - rmz) <= tol:
                    matched_rev += 1
                    used_q.add(j)
                    break

        # Use Dice coefficient: 2*|A∩B| / (|A| + |B|)
        total_q = len(q_peaks)
        total_r = len(r_peaks)
        if total_q + total_r == 0:
            return 0.0
        dice = 2.0 * max(matched, matched_rev) / (total_q + total_r)
        return min(1.0, dice)

    @staticmethod
    def _intensity_correlation_score(query: Dict[str, Any], ref: Dict[str, Any],
                                     tol: float = 0.02) -> float:
        """Pearson correlation of matched fragment ion intensities.

        For each pair of matched fragment ions, compare their relative
        intensities. High correlation indicates similar fragmentation
        patterns, supporting structural similarity.
        """
        q_peaks = query.get("peaks", [])
        r_peaks = ref.get("peaks", [])
        if not q_peaks or not r_peaks:
            return 0.0

        q_ints = []
        r_ints = []
        used_ref = set()
        for qp in q_peaks:
            qmz = qp[0] if len(qp) >= 2 else 0.0
            qi = qp[1] if len(qp) >= 2 else 0.0
            for j, rp in enumerate(r_peaks):
                if j in used_ref:
                    continue
                rmz = rp[0] if len(rp) >= 2 else 0.0
                ri = rp[1] if len(rp) >= 2 else 0.0
                if abs(qmz - rmz) <= tol:
                    q_ints.append(qi)
                    r_ints.append(ri)
                    used_ref.add(j)
                    break

        if len(q_ints) < 3:
            return 0.0

        q_arr = np.array(q_ints, dtype=float)
        r_arr = np.array(r_ints, dtype=float)

        # Pearson correlation
        q_mean = q_arr.mean()
        r_mean = r_arr.mean()
        q_centered = q_arr - q_mean
        r_centered = r_arr - r_mean

        denom = np.sqrt(float(np.sum(q_centered ** 2)) * float(np.sum(r_centered ** 2)))
        if denom < 1e-12:
            return 0.0

        pearson = float(np.sum(q_centered * r_centered)) / denom
        # Map [-1, 1] to [0, 1]
        return max(0.0, (pearson + 1.0) / 2.0)
