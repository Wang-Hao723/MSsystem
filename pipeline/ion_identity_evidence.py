"""Query-side IIMN input and chromatographic-evidence adapter.

The existing matching module is kept byte-identical so library preprocessing
identity need not change. Only the query-file reader and XIC shape gate are
replaced in a private function namespace. No globals are monkey-patched.
"""
from __future__ import annotations

import bisect
from collections import Counter
from pathlib import Path
from types import FunctionType

import numpy as np
import matching.ion_identity as base
from utils.spectrum_contract import is_auxiliary_detector_spectrum


def _scan_polarity(spectrum, *, mzxml=False):
    if mzxml:
        polarity = str(spectrum.get('polarity') or '').strip()
        return {'+': 'POS', '-': 'NEG'}.get(polarity, '')
    positive = 'positive scan' in spectrum
    negative = 'negative scan' in spectrum
    if positive == negative:
        return ''
    return 'POS' if positive else 'NEG'


def _rt_minutes(spectrum, *, mzxml=False):
    if not mzxml:
        return base._scan_rt_minutes(spectrum)
    raw = spectrum.get('retentionTime')
    value = base._number(raw)
    if value is None:
        return None
    # pyteomics mzXML converts ISO-8601 durations to unitfloat minutes.
    unit = str(getattr(raw, 'unit_info', '') or '').lower()
    if 'second' in unit:
        return value / 60.0
    if 'minute' in unit:
        return value
    # Do not guess the unit of an unsupported parser representation.
    return None


def extract_targeted_ms1_traces(file_path, targets, *, rt_half_window_min=.20,
                                mz_tolerance_da=.001, mz_tolerance_ppm=10.0):
    path = Path(file_path)
    suffix = path.suffix.lower()
    if suffix not in {'.mzml', '.mzxml', '.xml'}:
        return {}, dict(status='continuous_ms1_unavailable_for_format', ms1_scan_count=0)
    try:
        from pyteomics import mzml, mzxml
        if suffix == '.xml':
            from lxml import etree
            _, element = next(etree.iterparse(str(path), events=('start',)))
            local_name = etree.QName(element).localname.lower()
            if local_name not in {'mzml', 'indexedmzml', 'mzxml'}:
                return {}, dict(status='unsupported_xml_root', ms1_scan_count=0)
            is_mzxml = local_name == 'mzxml'
        else:
            is_mzxml = suffix == '.mzxml'
        reader_factory = mzxml.read if is_mzxml else mzml.read
    except Exception as exc:
        return {}, dict(status='ms1_reader_unavailable', error=str(exc), ms1_scan_count=0)

    valid = []
    for target in targets:
        identifier = str(target.get('id') or '')
        target_mz = base._number(target.get('mz'))
        target_rt = base._number(target.get('rt'))
        polarity = str(target.get('polarity') or '').strip().upper()
        if identifier and target_mz and target_mz > 0 and target_rt is not None:
            valid.append((target_rt, identifier, target_mz, polarity))
    valid.sort()
    centres = [item[0] for item in valid]
    traces = {identifier: dict(rt=[], intensity=[], target_mz=mz, polarity=polarity)
              for _, identifier, mz, polarity in valid}
    ms1_scan_count = 0
    counts = Counter()
    try:
        with reader_factory(str(path)) as reader:
            for spectrum in reader:
                if is_auxiliary_detector_spectrum(spectrum):
                    continue
                if spectrum.get('msLevel' if is_mzxml else 'ms level') != 1:
                    continue
                rt = _rt_minutes(spectrum, mzxml=is_mzxml)
                mzs = spectrum.get('m/z array')
                intensities = spectrum.get('intensity array')
                if rt is None or mzs is None or intensities is None:
                    continue
                mzs = np.asarray(mzs, dtype=float)
                intensities = np.asarray(intensities, dtype=float)
                if mzs.ndim != 1 or intensities.ndim != 1 or len(mzs) != len(intensities):
                    raise ValueError('MS1 arrays have incompatible dimensions')
                if not np.isfinite(mzs).all() or not np.isfinite(intensities).all():
                    raise ValueError('MS1 arrays contain nonfinite values')
                if len(mzs) > 1 and np.any(np.diff(mzs) < 0):
                    order = np.argsort(mzs, kind='mergesort')
                    mzs, intensities = mzs[order], intensities[order]
                polarity = _scan_polarity(spectrum, mzxml=is_mzxml)
                counts[polarity or 'UNKNOWN'] += 1
                ms1_scan_count += 1
                start = bisect.bisect_left(centres, rt - rt_half_window_min)
                end = bisect.bisect_right(centres, rt + rt_half_window_min)
                for _, identifier, target_mz, target_polarity in valid[start:end]:
                    # Only explicit same-polarity scans can support an ion-form
                    # relationship. Missing/conflicting polarity is unqualified.
                    if target_polarity not in {'POS', 'NEG'} or polarity != target_polarity:
                        continue
                    traces[identifier]['rt'].append(float(rt))
                    traces[identifier]['intensity'].append(base._xic_intensity(
                        mzs, intensities, target_mz, tolerance_da=mz_tolerance_da,
                        tolerance_ppm=mz_tolerance_ppm))
    except Exception as exc:
        return {}, dict(status='ms1_trace_extraction_failed', error=str(exc), ms1_scan_count=ms1_scan_count)

    usable = 0
    for trace in traces.values():
        intensities = np.asarray(trace['intensity'], dtype=float)
        trace['points'] = len(trace['rt'])
        trace['nonzero_points'] = int(np.count_nonzero(intensities > 0))
        trace['apex_intensity'] = float(np.max(intensities, initial=0.0))
        if trace['apex_intensity'] > 0 and trace['points'] >= 5:
            trace['apex_rt'] = float(trace['rt'][int(np.argmax(intensities))])
            trace['status'] = 'usable'
            usable += 1
        else:
            trace['apex_rt'] = None
            trace['status'] = 'insufficient_trace'
    return traces, dict(status='ok' if ms1_scan_count else 'no_ms1_scans',
                        ms1_scan_count=ms1_scan_count, target_count=len(valid),
                        usable_trace_count=usable, scan_polarity_counts=dict(counts),
                        backend='pyteomics_same_polarity_streaming_xic',
                        source_format='mzXML' if is_mzxml else 'mzML')


def chromatographic_shape_evidence(left, right, *, min_points=5, min_overlap=.50,
                                   min_correlation=.85, max_apex_delta_min=.10):
    result = base.chromatographic_shape_evidence(
        left, right, min_points=min_points, min_overlap=min_overlap,
        min_correlation=min_correlation, max_apex_delta_min=max_apex_delta_min)
    if not result['accepted']:
        return result
    mappings = [{round(float(rt), 8): float(value)
                 for rt, value in zip(trace.get('rt') or [], trace.get('intensity') or [])}
                for trace in (left, right)]
    common = sorted(set(mappings[0]).intersection(mappings[1]))
    support = []
    for mapping in mappings:
        values = np.asarray([mapping[rt] for rt in common], dtype=float)
        apex = int(np.argmax(values))
        support.append(dict(nonzero_points=int(np.count_nonzero(values > 0)),
                            nonzero_left=int(np.count_nonzero(values[:apex] > 0)),
                            nonzero_right=int(np.count_nonzero(values[apex + 1:] > 0))))
    sufficient = all(item['nonzero_points'] >= min_points and item['nonzero_left'] >= 2
                     and item['nonzero_right'] >= 2 for item in support)
    return {**result, 'accepted': sufficient,
            'status': result['status'] if sufficient else 'insufficient_nonzero_peak_shape_points',
            'nonzero_peak_shape_support': support}


_namespace = dict(base.annotate_ion_identity_networks.__globals__)
_namespace.update(extract_targeted_ms1_traces=extract_targeted_ms1_traces,
                  chromatographic_shape_evidence=chromatographic_shape_evidence)
annotate_ion_identity_networks = FunctionType(
    base.annotate_ion_identity_networks.__code__, _namespace,
    name=base.annotate_ion_identity_networks.__name__,
    argdefs=base.annotate_ion_identity_networks.__defaults__,
    closure=base.annotate_ion_identity_networks.__closure__)
annotate_ion_identity_networks.__kwdefaults__ = dict(
    base.annotate_ion_identity_networks.__kwdefaults__ or {})
