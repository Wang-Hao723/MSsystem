"""Query-only mzML reading; the accepted library-normalization module stays intact."""
from types import FunctionType
import re

from utils import mzml_io
from pipeline.query_charge_guard import guard_query_record


def load_from_mzml(filename, ms_level=2, metadata_harmonization=True, *, metadata_callback=None):
    """Keep valid empty spectra in the stream to preserve MS2 ordinal alignment."""
    import numpy as np
    from pyteomics.mzml import MzML
    from matchms import Spectrum
    from matchms.importing.parsing_utils import parse_mzml_mzxml_metadata, sort_by_mz

    class QueryMzML(MzML):
        # Retain the original XML text before the installed Pyteomics numeric
        # converter changes integer zero to None. This is a query-local reader;
        # the library reader and global Pyteomics classes are untouched.
        def _get_info_smart(self, element, **kwargs):
            terms = []
            if str(element.tag).rsplit("}", 1)[-1] == "selectedIon":
                terms = [child for child in element
                         if child.get("accession") == "MS:1000041"]
                if len(terms) > 1:
                    raise ValueError("mzML selected ion has multiple charge-state declarations")
                if terms:
                    raw_charge = terms[0].get("value", "")
                    # PSI-MS MS:1000041 has xsd:int value type. Validate the
                    # original text before Pyteomics can truncate 1.5 to 1.
                    if raw_charge.strip() and not re.fullmatch(
                            r"[+-]?[0-9]+", raw_charge.strip()):
                        raise ValueError(
                            "mzML charge state must be an integer literal; "
                            "empty or zero denotes unknown charge")
            value = super()._get_info_smart(element, **kwargs)
            if terms:
                value["_query_charge_state_raw"] = raw_charge
            return value

    def omit_missing_charge(value):
        # matchms calls int() for a present charge-state key. A missing value
        # means unknown, so omit only that key in its private metadata view.
        # NumPy peak arrays are passed through by reference, never duplicated.
        if isinstance(value, dict):
            return {key: omit_missing_charge(item) for key, item in value.items()
                    if not (key == "charge state" and item is None)}
        if isinstance(value, list):
            return [omit_missing_charge(item) for item in value]
        return value

    with QueryMzML(str(filename), dtype=dict, read_schema=False, iterative=True,
                   use_index=False, huge_tree=False, decode_binary=True) as reader:
        ordinal = 0
        for raw in reader:
            if raw.get("ms level") != ms_level:
                continue
            mz, intensities = raw.get("m/z array"), raw.get("intensity array")
            if mz is None or intensities is None:
                declared_empty = str(raw.get("defaultArrayLength")) == "0"
                contains_points = any(value is not None and len(value) for value in (mz, intensities))
                if not declared_empty or contains_points:
                    raise ValueError(f"mzML spectrum {raw.get('id', '')!r} is missing required peak arrays")
                mz, intensities = [], []
            mz, intensities = sort_by_mz(
                mz=np.asarray(mz, dtype=float), intensities=np.asarray(intensities, dtype=float)
            )
            metadata = parse_mzml_mzxml_metadata(omit_missing_charge(raw))
            if metadata_callback is not None:
                metadata.update(metadata_callback(raw, ordinal))
            ordinal += 1
            yield Spectrum(mz=mz, intensities=intensities,
                           metadata=metadata,
                           metadata_harmonization=metadata_harmonization)


def parse_mzml_file(file_path, isolation_window_width_da=None,
                    isolation_window_width_source=None, allow_unknown_representation=False):
    """Reuse the exact parser with a private reader binding, never a global patch.

    The library cache binds every utils module by content. Keeping this sample
    adapter in pipeline avoids changing or weakening that accepted identity.
    A fresh globals dictionary isolates concurrent calls and preserves the
    parser's existing metadata, isotope, polarity, and representation rules.
    """
    from pipeline.mzml_metadata import acquisition_metadata, read_instrument_header

    header = read_instrument_header(file_path)
    evidence = {}

    def extract(*args, **kwargs):
        evidence.update(mzml_io._extract_mzml_precursor_evidence(*args, **kwargs))
        return evidence

    def enrich(raw, ordinal):
        values = acquisition_metadata(raw, header)
        record = evidence.get(("ms2_ordinal", ordinal))
        if record is None:
            raise ValueError("mzML precursor metadata could not be linked to the query scan")
        record.update(values)
        precursors = (raw.get("precursorList") or {}).get("precursor") or []
        ions = ((precursors[0].get("selectedIonList") or {}).get("selectedIon")
                or []) if precursors else []
        if ions and "_query_charge_state_raw" in ions[0]:
            charge_raw = ions[0]["_query_charge_state_raw"]
            record["reported_precursor_charge"] = charge_raw
            if re.fullmatch(r"[+-]?\d+", charge_raw.strip()):
                magnitude = abs(int(charge_raw))
                record["precursor_charge"] = magnitude
                record["precursor_charge_status"] = (
                    "unknown_zero" if magnitude == 0 else "explicit")
            elif not charge_raw.strip():
                record["precursor_charge"] = None
                record["precursor_charge_status"] = "missing_value"
        metadata = {key: values[key] for key in ("instrument", "instrument_type") if values.get(key)}
        if values.get("scan_polarity"):
            metadata["ionmode"] = "positive" if values["scan_polarity"] == "POS" else "negative"
        return metadata

    def reader(filename):
        return load_from_mzml(filename, metadata_callback=enrich)

    original = mzml_io.parse_mzml_file
    namespace = dict(original.__globals__, load_from_mzml=reader,
                     _extract_mzml_precursor_evidence=extract,
                     # Each scan resolves its own configuration above. An
                     # unknown IC must not inherit the first IC's analyzer.
                     extract_mzml_instrument_metadata=lambda _: {})
    parser = FunctionType(original.__code__, namespace, original.__name__,
                          original.__defaults__, original.__closure__)
    parser.__kwdefaults__ = original.__kwdefaults__
    records = parser(file_path, isolation_window_width_da=isolation_window_width_da,
                     isolation_window_width_source=isolation_window_width_source,
                     allow_unknown_representation=allow_unknown_representation)
    for record in records:
        meta = evidence.get(record["spectrum_id"], {})
        if meta.get("instrument_configuration_id"):
            for key in ("instrument", "instrument_type"):
                record[key] = meta.get(key, "")
        for key in ("collision_energy_quantities", "instrument_configuration_id",
                    "acquisition_metadata_source", "ion_injection_time",
                    "mass_resolving_power", "preset_scan_configuration", "scan_filter_string",
                    "reported_precursor_charge", "precursor_charge_status"):
            if key in meta:
                record[key] = meta[key]
    return [guard_query_record(record) for record in records]
