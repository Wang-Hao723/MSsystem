"""Source-verified, read-only corrections above the immutable spectral cache.

The metadata registry binds every original normalized record and all source files.
Its metadata wrapper preserves peaks, structures, masses and reference membership.
The full-library loader first applies the separate source-bound compound mass
repair overlay, including corrected neutral-mass index positions. Peaks,
structures, precursor masses and the precursor-mass index remain unchanged.
"""
from collections.abc import Sequence
from pathlib import Path
import hashlib, json
import threading
from runtime.source_file_verification import verify_registered_sources
FIELDS = {'charge', 'source_ion_mode', 'adduct', 'source_accession',
          'instrument_type', 'instrument', 'source_ms_level'}
REGISTRY_PATH = Path(__file__).resolve().parents[1] / 'config/library_metadata_repairs.json'
_loader_lock = threading.RLock()

def load_library_with_verified_repairs(db_path=None):
    """Load the accepted full cache, then restore exactly source-bound fields.

    Small libraries remain unchanged. Full-library binding failures abort loading.
    The search module and its callers share the same records object so neutral
    mass and source-family lookup retain their existing behavior.
    """
    import importlib
    from config.settings import DB_PATH
    search_module = importlib.import_module('matching.search')
    requested = Path(db_path or DB_PATH).resolve()
    with _loader_lock:
        if requested != Path(DB_PATH).resolve():
            records, index = search_module.load_library(requested)
            return records, index
        records = getattr(search_module, '_lib_records', None)
        if records is not None and getattr(search_module, '_lib_db_path', None) == requested:
            index = search_module._mz_index
        elif getattr(search_module.load_library, '__module__', None) != 'matching.search':
            # Test/integration harnesses may provide a bound loader with their
            # own source database and cache fixture.
            records, index = search_module.load_library(requested)
        else:
            from runtime.library_cache_identity import load_runtime_compatible_record_store
            cached = load_runtime_compatible_record_store(requested, project_root=Path(__file__).resolve().parents[1])
            if cached is None:
                records, index = search_module.load_library(requested)
            else:
                records, index, neutral_index = cached
                search_module._lib_records = records
                search_module._mz_index = index
                search_module._nm_index = neutral_index
                search_module._lib_db_path = requested
                search_module._lib_source_families = records.source_families
        registry = validate_registry(json.loads(REGISTRY_PATH.read_text()))
        if isinstance(records, VerifiedMetadataRecords):
            if records.registry['semantic_sha256'] != registry['semantic_sha256']:
                raise ValueError('Loaded metadata registry differs; restart required')
            from runtime.source_compound_mass_repairs import get_source_compound_mass_state, source_compound_mass_status
            mass_state = get_source_compound_mass_state()
            if mass_state is not None and not source_compound_mass_status(records, True).get('normalized_overlay_applied'):
                raise ValueError('Loaded compound mass generation differs; restart required')
            return records, index
        try:
            from runtime.source_compound_mass_repairs import apply_source_compound_mass_records
            mass_records, neutral_index, mass_identity = apply_source_compound_mass_records(
                records, search_module._nm_index, bucket_size=search_module.MZ_BUCKET_DALTON)
            wrapped = VerifiedMetadataRecords(mass_records, registry)
            from runtime.source_compound_mass_repairs import mark_source_compound_mass_published
            mark_source_compound_mass_published(wrapped)
        except Exception:
            # Never leave a failed full-library load exposed as an accepted base.
            search_module._lib_records = None
            search_module._mz_index = None
            search_module._nm_index = None
            search_module._lib_db_path = None
            search_module._lib_source_families = ()
            raise
        # Publish both corrected record values and neutral positions together.
        search_module._lib_records = wrapped
        search_module._nm_index = neutral_index
        return wrapped, index

def library_metadata_repair_status(records, loaded):
    """Report actual application, rather than presence of a configuration file."""
    if not loaded:
        return {'status': 'not_loaded', 'applied_reference_count': 0}
    if isinstance(records, VerifiedMetadataRecords):
        return {'status': 'applied', **records.repair_status()}
    return {'status': 'not_applied', 'applied_reference_count': 0}

def metadata_repair_preflight(project_root):
    """Check registry, base manifest and source files before reporting readiness."""
    root = Path(project_root)
    path = root / 'config/library_metadata_repairs.json'
    if not path.is_file():
        return {'compatible': True, 'status': 'not_configured'}
    try:
        registry = validate_registry(json.loads(path.read_text()))
        import os
        config = json.loads((root / 'config/library_record_store.json').read_text())
        override = os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE', '').strip()
        if override.lower() in {'off', '0', 'false'}:
            raise ValueError('Verified metadata repairs require their accepted base cache')
        manifest = json.loads((Path(override or config['cache_dir']) / 'manifest.json').read_text())
        for base_key, registry_key in [('normalized_records_sha256', 'base_normalized_records_sha256'), ('record_count', 'base_record_count')]:
            if manifest[base_key] != registry[registry_key]:
                raise ValueError('Metadata repair base mismatch: ' + base_key)
        if manifest['source_db']['sha256'] != registry['source_database_sha256']:
            raise ValueError('Metadata repair database binding mismatch')
        verify_registered_sources(registry['source_files'])
        return {'compatible': True, 'status': 'source_bindings_verified',
                'registry_id': registry['registry_id'], 'registered_reference_count': len(registry['rows']),
                'original_record_integrity_verified_by': 'library_load'}
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        return {'compatible': False, 'status': 'invalid_binding', 'error': str(exc)}

def digest(x):
    return hashlib.sha256(json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=True).encode()).hexdigest()

def validate_registry(registry):
    if not registry['schema_version'] == 'mssystem.source-metadata-repairs/1.0.0':
        raise ValueError("Invalid source metadata repair binding: registry['schema_version'] == 'mssystem.source-metadata-repairs/1.0.0'")
    if not registry['semantic_sha256'] == digest({k: v for k, v in registry.items() if k != 'semantic_sha256'}):
        raise ValueError("Invalid source metadata repair binding: registry['semantic_sha256'] == digest({k: v for k, v in registry.items() if k != 'semantic_sha256'})")
    indices = set()
    for row in registry['rows']:
        if not (isinstance(row['record_index'], int) and row['record_index'] >= 0):
            raise ValueError("Invalid source metadata repair binding: isinstance(row['record_index'], int) and row['record_index'] >= 0")
        if not row['record_index'] not in indices:
            raise ValueError("Invalid source metadata repair binding: row['record_index'] not in indices")
        indices.add(row['record_index'])
        if not (set(row['updates']) <= FIELDS and set(row['updates'])):
            raise ValueError("Invalid source metadata repair binding: set(row['updates']) <= FIELDS and set(row['updates'])")
        if not row['source_binding']['structure_key_equal'] is True:
            raise ValueError("Invalid source metadata repair binding: row['source_binding']['structure_key_equal'] is True")
        if not (row['source_binding']['raw_peaks_equal'] or row['source_binding']['existing_importer_peak_transform_equal']):
            raise ValueError("Invalid source metadata repair binding: row['source_binding']['raw_peaks_equal'] or row['source_binding']['existing_importer_peak_transform_equal']")
    return registry

class VerifiedMetadataRecords(Sequence):

    def __init__(self, backing, registry, *, verify_sources=True):
        from matching.library_polarity import classify_spectrum_polarity
        self.backing = backing
        self.registry = validate_registry(registry)
        self.corrected = {}
        manifest = backing.manifest
        if not manifest['normalized_records_sha256'] == registry['base_normalized_records_sha256']:
            raise ValueError('Different base library')
        if not manifest['record_count'] == registry['base_record_count']:
            raise ValueError("Invalid source metadata repair binding: manifest['record_count'] == registry['base_record_count']")
        if not manifest['source_db']['sha256'] == registry['source_database_sha256']:
            raise ValueError("Invalid source metadata repair binding: manifest['source_db']['sha256'] == registry['source_database_sha256']")
        if verify_sources:
            verify_registered_sources(registry['source_files'])
        for repair in registry['rows']:
            i = repair['record_index']
            original = backing[i]
            if not digest(original) == repair['original_record_sha256']:
                raise ValueError('Original reference changed')
            if not str(original['spectrum_id']) == repair['spectrum_id']:
                raise ValueError("Invalid source metadata repair binding: str(original['spectrum_id']) == repair['spectrum_id']")
            fixed = {**original, **repair['updates']}
            fixed['polarity'] = classify_spectrum_polarity(fixed.get('adduct'), fixed.get('charge'), fixed.get('source_ion_mode'))
            if not fixed['polarity'] == repair['source_polarity']:
                raise ValueError("Invalid source metadata repair binding: fixed['polarity'] == repair['source_polarity']")
            fixed['source_metadata_repair'] = {'registry_id': registry['registry_id'], 'fields_restored': sorted(repair['updates']), 'source_block_sha256': repair['source_binding']['block_sha256'], 'spectral_peaks_and_structure_unchanged': True, 'original_mass_index_retained': True}
            self.corrected[i] = fixed

    def __len__(self):
        return len(self.backing)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[j] for j in range(*i.indices(len(self)))]
        if i < 0:
            i += len(self)
        if i < 0 or i >= len(self):
            raise IndexError(i)
        return dict(self.corrected[i]) if i in self.corrected else self.backing[i]

    def read_uncached(self, i):
        return dict(self.corrected[i]) if i in self.corrected else self.backing.read_uncached(i)

    def __getattr__(self, name):
        return getattr(self.backing, name)

    def repair_status(self):
        return {'registry_id': self.registry['registry_id'], 'applied_reference_count': len(self.corrected), 'registry_semantic_sha256': self.registry['semantic_sha256'], 'base_normalized_records_sha256': self.registry['base_normalized_records_sha256'], 'source_database_written': False, 'peaks_structures_and_mass_indices_changed': False}
