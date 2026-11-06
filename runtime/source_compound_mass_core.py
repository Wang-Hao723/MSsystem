"""Research-only, source-bound compound mass/formula overlay plus neutral index.

The immutable database/cache remain unchanged. This module deliberately does
not extend the unrelated production spectral-metadata repair field whitelist.
Apply records and the returned neutral index together in one loading step.
"""
from collections.abc import Sequence
from pathlib import Path
import hashlib,json,math
import numpy as np

SCHEMA='mssystem.source-compound-mass-repairs/1.0.0'
ENGINE_VERSION='2.0.0'
FIELDS={'exact_mass','formula','neutral_mass','neutral_mass_source'}
def digest(value):
 return hashlib.sha256(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=True).encode()).hexdigest()
def positive(value):
 return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and value>0

def validate_registry(registry):
 from rdkit import Chem
 from rdkit.Chem import rdMolDescriptors
 from matching.library_quality import normalize_formula_text,formula_exact_mass
 if registry.get('schema_version')!=SCHEMA or registry.get('semantic_sha256')!=digest({k:v for k,v in registry.items() if k!='semantic_sha256'}):
  raise ValueError('Invalid compound-mass registry identity')
 compounds={}
 for compound in registry['compounds']:
  key=compound['full_inchikey'];smi=compound['smiles'];formula=compound['formula'];mass=compound['source_exact_mass']
  if key in compounds or not positive(mass):raise ValueError('Duplicate identity or invalid source mass')
  mol=Chem.MolFromSmiles(smi)
  if mol is None or Chem.GetFormalCharge(mol) or len(Chem.GetMolFrags(mol))!=1 or any(a.GetIsotope() for a in mol.GetAtoms()):raise ValueError('Ineligible neutral structure')
  if Chem.MolToInchiKey(mol)!=key or normalize_formula_text(rdMolDescriptors.CalcMolFormula(mol))!=normalize_formula_text(formula):raise ValueError('Structure identity/formula mismatch')
  derived=rdMolDescriptors.CalcExactMolWt(mol);fm=formula_exact_mass(formula)
  if fm is None or abs(derived-fm)>1e-5 or abs(derived-mass)>1e-5:raise ValueError('Source mass differs from formula/structure')
  if not compound['original_source_evidence']:raise ValueError('Original source evidence missing')
  for ev in compound['original_source_evidence']:
   if ev['source_file'] not in registry['source_files'] or ev['formula']!=formula or ev['source_exact_mass']!=mass or ev['full_inchikey']!=key:raise ValueError('Source evidence differs')
  compounds[key]=compound
 seen=set()
 for row in registry['rows']:
  i=row['record_index'];u=row['updates']
  if isinstance(i,bool) or not isinstance(i,int) or not 0<=i<registry['base_record_count'] or i in seen:raise ValueError('Invalid/duplicate reference index')
  seen.add(i);compound=compounds[row['full_inchikey']]
  if set(u)!=FIELDS or u!={'exact_mass':compound['source_exact_mass'],'formula':compound['formula'],'neutral_mass':compound['source_exact_mass'],'neutral_mass_source':'compound_exact_mass'}:raise ValueError('Only explicitly verified compound fields may change')
  if row['original_values']['neutral_mass']!=row['old_neutral_mass'] or not positive(row['old_neutral_mass']):raise ValueError('Invalid old mass binding')
 if set(registry['excluded_quarantined_spectrum_ids']) & {r['spectrum_id'] for r in registry['rows']}:raise ValueError('Quarantined reference included')
 return compounds

def source_verification_timeout(expected):
 # Budget for a cold external drive; this changes waiting, never acceptance.
 # Clamp integers before division so malformed enormous sizes cannot overflow.
 size=expected.get('size_bytes')
 if isinstance(size,bool) or not isinstance(size,int) or size<0:
  raise ValueError('Invalid registered source size')
 rate=16*1024*1024
 return max(120.,60.+min(size,1740*rate)/rate)

def verify_source_evidence(registry,compounds):
 from runtime.source_file_verification import verify_source_file
 from rdkit import Chem
 # Full content, source blocks and literal fields are all still checked.
 # A fixed 120 seconds was insufficient for the registered 4.5/5.1GB files
 # during a real cold start. The per-file limit remains at most 1800 seconds.
 for name,expected in registry['source_files'].items():
  verify_source_file(name,expected,timeout_seconds=source_verification_timeout(expected),attempts=1)
 for compound in compounds.values():
  for ev in compound['original_source_evidence']:
   with Path(ev['source_file']).open('rb') as f:f.seek(ev['byte_offset']);raw=f.read(ev['byte_length'])
   if hashlib.sha256(raw).hexdigest()!=ev['block_sha256']:raise ValueError('Original evidence block changed')
   # This registry version explicitly accepts literal MGF formula+mass fields.
   fields={line.split('=',1)[0].strip().upper():line.split('=',1)[1].strip() for line in raw.decode('utf-8').splitlines() if '=' in line}
   if fields.get('FORMULA')!=ev['formula'] or float(fields.get('EXACTMASS','nan'))!=ev['source_exact_mass']:raise ValueError('Original source literal fields differ')
   mol=Chem.MolFromSmiles(fields.get('SMILES',''))
   if mol is None or Chem.MolToInchiKey(mol)!=ev['full_inchikey']:raise ValueError('Original source structure differs')

class _MassRecords(Sequence):
 def __init__(self,backing,registry,base_manifest,identity):
  self.backing=backing;self.registry=registry;self.corrected={};self.overlay_identity=identity
  self.manifest={**base_manifest,'effective_compound_mass_overlay':identity}
  for repair in registry['rows']:
   i=repair['record_index'];original=backing[i]
   if digest(original)!=repair['original_record_sha256'] or str(original['spectrum_id'])!=repair['spectrum_id'] or original.get('inchikey')!=repair['full_inchikey']:raise ValueError('Original record identity changed')
   if {k:original.get(k) for k in FIELDS}!=repair['original_values']:raise ValueError('Original mass/formula fields changed')
   fixed={**original,**repair['updates']}
   fixed['source_compound_mass_repair']={'registry_id':registry['registry_id'],'registry_semantic_sha256':registry['semantic_sha256'],
    'original_values':repair['original_values'],'source_compound_id':repair['sql_compound_id'],
    'mass_formula_evidence':'Source literal exact mass and formula agree with neutral single-component unlabeled full-key-identical structure',
    'charge_adduct_measurements_and_structure_unchanged':True}
   self.corrected[i]=fixed
 def __len__(self):return len(self.backing)
 def __getitem__(self,i):
  if isinstance(i,slice):return [self[j] for j in range(*i.indices(len(self)))]
  if i<0:i+=len(self)
  if i<0 or i>=len(self):raise IndexError(i)
  return dict(self.corrected[i]) if i in self.corrected else self.backing[i]
 def read_uncached(self,i):
  if i in self.corrected:return dict(self.corrected[i])
  if hasattr(self.backing,'read_uncached'):return self.backing.read_uncached(i)
  return self.backing[i]
 def __getattr__(self,name):return getattr(self.backing,name)

class SparseNeutralMassIndex:
 """Only seven neutral positions differ; preserve mass-then-index ordering."""
 def __init__(self,base,registry,identity):
  self.base=base;self.overlay_identity=identity
  self._changed={r['record_index']:r for r in registry['rows']}
  new_pairs=sorted((r['updates']['neutral_mass'],r['record_index']) for r in registry['rows'])
  self._new_masses=np.asarray([m for m,i in new_pairs],dtype=np.float64)
  self._new_indices=np.asarray([i for m,i in new_pairs],dtype=np.uint64)
  positions=np.flatnonzero(np.isin(base._record_indices,list(self._changed)))
  if len(positions)!=len(self._changed):raise ValueError('Neutral index missing or duplicating original references')
  for pos in positions:
   row=self._changed[int(base._record_indices[pos])]
   if float(base._masses[pos])!=row['old_neutral_mass']:raise ValueError('Neutral index old mass mismatch')
 def __len__(self):return len(self.base)
 def range_query(self,center,tolerance):
  center=float(center);tolerance=float(tolerance)
  if not math.isfinite(center) or not math.isfinite(tolerance) or tolerance<0:return []
  low,high=center-tolerance,center+tolerance
  left=int(np.searchsorted(self.base._masses,low,side='left'));right=int(np.searchsorted(self.base._masses,high,side='right'))
  values=[(float(m),int(i)) for m,i in zip(self.base._masses[left:right],self.base._record_indices[left:right]) if int(i) not in self._changed]
  new_left=int(np.searchsorted(self._new_masses,low,side='left'))
  new_right=int(np.searchsorted(self._new_masses,high,side='right'))
  values.extend((float(m),int(i)) for m,i in zip(self._new_masses[new_left:new_right],self._new_indices[new_left:new_right]))
  return [i for m,i in sorted(values)]

def apply_verified_mass_repairs(records,neutral_index,registry,*,base_manifest=None,bucket_size=None,verify_sources=True):
 """Return the record/neutral-index pair atomically; never alter the m/z index."""
 compounds=validate_registry(registry);base_manifest=base_manifest or getattr(records,'manifest',None)
 if not base_manifest:raise ValueError('Full base cache/database manifest required')
 if base_manifest['record_count']!=len(records) or len(records)!=registry['base_record_count'] or base_manifest['normalized_records_sha256']!=registry['base_normalized_records_sha256'] or base_manifest['source_db']['sha256']!=registry['source_database_sha256']:raise ValueError('Base library binding differs')
 if base_manifest.get('record_alignment_sha256')!=registry['base_record_alignment_sha256']:raise ValueError('Base record alignment differs')
 if verify_sources:verify_source_evidence(registry,compounds)
 identity={'schema_version':SCHEMA,'engine_version':ENGINE_VERSION,'registry_semantic_sha256':registry['semantic_sha256'],
  'base_normalized_records_sha256':registry['base_normalized_records_sha256'],
  'base_record_alignment_sha256':registry['base_record_alignment_sha256'],
  'effective_record_and_neutral_index_identity':digest({'base_record_alignment_sha256':registry['base_record_alignment_sha256'],'registry_semantic_sha256':registry['semantic_sha256']}),
  'corrected_reference_count':len(registry['rows']),'precursor_index_unchanged':True}
 wrapped=_MassRecords(records,registry,base_manifest,identity)
 if isinstance(neutral_index,dict):
  if not positive(bucket_size):raise ValueError('Explicit bucket size required for bucket index')
  patched=dict(neutral_index);affected=set()
  for row in registry['rows']:
   old=int(row['old_neutral_mass']/bucket_size);new=int(row['updates']['neutral_mass']/bucket_size);affected.update((old,new))
   if list(neutral_index.get(old,[])).count(row['record_index'])!=1:raise ValueError('Original bucket position missing/duplicated')
  changed={r['record_index'] for r in registry['rows']}
  for b in affected:patched[b]=[i for i in neutral_index.get(b,[]) if i not in changed]
  for row in registry['rows']:patched[int(row['updates']['neutral_mass']/bucket_size)].append(row['record_index'])
  for b in affected:
   patched[b].sort(key=lambda i:(wrapped[i]['neutral_mass'],i))
   if not patched[b]:del patched[b]
  return wrapped,patched,identity
 if not hasattr(neutral_index,'_masses') or not hasattr(neutral_index,'_record_indices'):raise ValueError('Unsupported neutral index representation')
 return wrapped,SparseNeutralMassIndex(neutral_index,registry,identity),identity
