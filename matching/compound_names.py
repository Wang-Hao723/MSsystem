"""Conservative compound-name resolution for public-library identifiers."""

from __future__ import annotations

import json
import re
from html import unescape
from functools import lru_cache
from pathlib import Path
from typing import Any


_ALIASES_PATH = Path(__file__).with_name("compound_name_aliases.json")
_PUBCHEM_ALIASES_PATH = Path(__file__).with_name(
    "compound_name_aliases_pubchem.json"
)
_CURATED_ALIASES_PATH = Path(__file__).with_name("compound_name_aliases_curated.json")
_LABELED_ID_RE = re.compile(
    r"^(?P<identifier>(?:NCGC|MLS)\d+(?:-\d+)?)(?P<separator>[!_])?(?P<payload>.*)$",
    re.IGNORECASE,
)
_FORMULA_PREFIX_RE = re.compile(r"^C\d+H\d+(?:[A-Z][a-z]?\d*)*_", re.IGNORECASE)
_CAS_RE = re.compile(r"^(?:CAS(?:RN)?\s*[:#]?\s*)?\d{2,7}-\d{2}-\d$", re.IGNORECASE)
_ACCESSION_RE = re.compile(
    r"^(?:DTXSID\d+|DTXCID\d+|HMDB\d+|LMFA\d+|LMSD\d+|CHEBI[:_]?\d+|"
    r"(?:S?CHEMBL)\d+|MLS\d+(?:-\d+)?|NPASS\d+|PUBCHEM\s*CID\s*:?\s*\d+|CID\s*:?\s*\d+|"
    r"AKOS\d+|Z\d{6,}|\d{4,})$",
    re.IGNORECASE,
)
_INCHIKEY_NAME_RE = re.compile(r"^[A-Z]{14}(?:-[A-Z]{10}-[A-Z])?$")
_PURE_NUMERIC_RE = re.compile(r"^[+-]?\d+(?:\.\d+)?$")
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
_FORMULA_CHARGE_RE = re.compile(r"\^?\d*[+-]$")
_ELEMENT_SYMBOLS = frozenset(
    """
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni
    Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I
    Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt
    Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr
    Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
    """.split()
)
# Compact library/development codes are provenance labels, not chemical names.
# The deliberately narrow 2–4 letter prefix avoids swallowing ordinary names
# such as ``ANILINE-2-SULFONATE`` or familiar digit-leading abbreviations such
# as ``6PPD``.  Exact full-InChIKey aliases are resolved before this guard, so
# a verified code such as AO4703 can still display its canonical chemical name.
_COMPACT_CODE_RE = re.compile(
    r"^(?P<prefix>[A-Z]{2,4})(?:[-_])?\d[A-Z0-9_-]{1,18}$"
)
_CHEMICAL_DESCRIPTOR_PREFIXES = {"ANTI", "BETA", "CIS", "SYN"}
_TRAILING_CAS_RE = re.compile(r"(?<=[A-Za-z)])\s*\d{2,7}-\d{2}-\d$")
_TRAILING_METADATA_RE = re.compile(
    r"\s+(?:CollisionEnergy\s*:\s*\S+|"
    r"(?:\[\d*M(?:[+-][A-Za-z0-9]+)*\]\d*[+-]?|"
    r"\d*M(?:[+-][A-Za-z0-9]+)+(?:\d*[+-])?|M[+-]?))\s*$",
    re.IGNORECASE,
)
_SPECTRAL_MATCH_TITLE_RE = re.compile(
    r"^Spectral\s+Match\s+to\s+(?P<name>.+?)\s+from\s+"
    r"(?P<library>NIST(?:\s*14)?|METLIN)$",
    re.IGNORECASE,
)
_QUOTED_COLLISION_TITLE_RE = re.compile(
    r'^[\"“](?P<name>.+?)\s+CollisionEnergy\s*:\s*[^\"”]+[\"”](?:\s+.*)?$',
    re.IGNORECASE,
)
_UNQUOTED_COLLISION_TITLE_RE = re.compile(
    r"^(?P<name>.+?)\s+CollisionEnergy\s*:\s*[0-9.,;/_-]+(?:\s+.*)?$",
    re.IGNORECASE,
)
_DEREPLICATOR_TITLE_RE = re.compile(
    r"^Dereplicator\s+Identification\s*-\s*(?P<name>.+)$",
    re.IGNORECASE,
)
_MISSING = {"", "nan", "none", "null", "na", "n/a", "n.a.", "-"}
_NOISE_RE = re.compile(r"^(?:\$:.*|spectrum[_\s-]*\d+|unknown|unnamed|anonymous)$", re.I)
_INTERNAL_ALIAS_RE = re.compile(
    r"^(?:MEG[XPMA-Z0-9]*|ACON\d*|UNK|TMP|LIB|SPEC)[_-]?\d{3,}$",
    re.IGNORECASE,
)
_FRAGMENT_ONLY_ALIAS_RE = re.compile(r"^\[[^\]]+\](?:[+-]|\d+[+-])?$")
_ALIAS_METADATA_RE = re.compile(
    r"(?:_major|_minor|putative|candidate|dereplicator\s+identification|"
    r"spectral\s+library\s+match|sirius\s+atomic\s+difference|"
    r"observed\s+from|collision\s*energy)",
    re.IGNORECASE,
)
_GENERIC_FRAGMENT_ALIASES = {
    "methyl", "ethyl", "propyl", "butyl", "hydroxyl", "phenyl",
    "aglycone", "glycoside", "sugar", "fragment", "unknown",
}
# --- Library-provenance extraction ---------------------------------------
# Source tools (GNPS dereplicator, spectral-library matches) prepend a
# provenance sentence to a real reference name, e.g.
# ``Dereplicator Identification - Cyclohomonervilasterol_Cyclohomonervilasterol``.
# Stripping these wrappers recovers the embedded chemical name instead of
# discarding the whole record as "Unresolved compound (...)".
_PROVENANCE_PREFIX_RE = re.compile(
    r"^(?:Dereplicator\s+Identification\s*-\s*|"
    r"Spectral\s+(?:Library\s+)?Match\s*(?:to\s+|[-:]\s*)?|"
    r"(?:MS2|MS/MS)\s+Match\s*(?:to\s+|[-:]\s*)?|"
    r"Library\s+(?:Hit|Match)\s*[-:]?\s*|"
    r"(?:ReSpect|GNPS|MassBank|NIST)\s*[:]\s*[A-Za-z0-9_]+\s+)",
    re.IGNORECASE,
)
# MassBank (EU + mirrors) imports concatenate the accession and several
# CH$NAME fields (preferred / class / systematic) with '|', e.g.
# ``MassbankEU:ETS00001 C10-LAS (TENTATIVE)|C10-linear alkylbenzyl sulfonate|
#   4-(decan-4-yl)benzenesulfonic acid``.  The whole concatenation is readable
# enough to pass ``alias_quality``, so without a dedicated handler it is shown
# verbatim.  We strip the accession prefix and keep only the most specific
# readable segment.  ``[A-Z]*`` (IGNORECASE) absorbs the mirror suffix (``EU``)
# without swallowing ordinary names; the ``:`` anchor keeps it specific.
_MASSBANK_PREFIX_RE = re.compile(
    r"^Massbank[A-Z]*\s*[:]\s*[A-Za-z0-9]{2,}(?:\s+|$)",
    re.IGNORECASE,
)
# Parenthetical caveats that wrap a real name, e.g.
# ``Aflatoxin B1 (putative_Observed from A. flavus)``.
_PAREN_PROVENANCE_RE = re.compile(
    r"\s*\([^)]*(?:putative|candidate|observed\s+from|collision|"
    r"major|minor|library\s+alias)[^)]*\)\s*",
    re.IGNORECASE,
)
# Bare trailing caveat words: ``Sagitol B putative``, ``Surfactin Unk``.
_TRAILING_CAVEAT_RE = re.compile(
    r"\s+(?:putative|candidate|unk(?:nown)?)\s*$", re.IGNORECASE
)
# Unconfirmed computational formula annotations (SIRIUS ``Candidate ... delta
# mass / sirius atomic difference``).  These MUST NOT be presented as confirmed
# compound names; they receive an honest placeholder label instead.  ``Candidate``
# is matched as a standalone word (it never appears inside a real chemical name)
# so bracketed SIRIUS class tags like ``[BA_core:...] Candidate ...`` are caught
# even when they carry no explicit delta-mass provenance.
_COMPUTATIONAL_CANDIDATE_RE = re.compile(
    r"(?:sirius\s+atomic\s+difference|delta\s+mass|"
    r"literature\s+atomic\s+difference|\bCandidate\s+\S)",
    re.IGNORECASE,
)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in _MISSING else text


def _identifier_core(value: Any) -> str:
    text = _clean(value).strip(" _!\"'“”‘’")
    previous = None
    while text and text != previous:
        previous = text
        text = _TRAILING_METADATA_RE.sub("", text).strip(" _!\"'“”‘’")
    return text


def is_molecular_formula_name(value: Any) -> bool:
    """Return True when a source title is only a simple molecular formula.

    Formula-only library labels are useful metadata, but they are not compound
    names. Validate every token against the periodic table instead of relying
    on a loose ``C...H...`` prefix so readable names containing digits are not
    accidentally hidden.
    """
    text = _identifier_core(value).replace(" ", "")
    if len(text) >= 2 and text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
    text = _FORMULA_CHARGE_RE.sub("", text)
    if not text or len(text) > 120:
        return False

    tokens = list(_FORMULA_TOKEN_RE.finditer(text))
    if len(tokens) < 2 or "".join(match.group(0) for match in tokens) != text:
        return False
    if any(match.group(1) not in _ELEMENT_SYMBOLS for match in tokens):
        return False
    if any(match.group(2) == "0" for match in tokens):
        return False
    return any(match.group(2) for match in tokens) or len(
        {match.group(1) for match in tokens}
    ) >= 2


def is_identifier_name(value: Any) -> bool:
    """True for a registry/accession label that is not a chemical name."""
    text = _identifier_core(value)
    if not text:
        return False
    labeled = _LABELED_ID_RE.match(text)
    if labeled:
        return not labeled.group("payload").strip(" _!")
    return bool(
        _PURE_NUMERIC_RE.fullmatch(text)
        or _CAS_RE.fullmatch(text)
        or _ACCESSION_RE.fullmatch(text)
        or _INCHIKEY_NAME_RE.fullmatch(text)
        or is_molecular_formula_name(text)
        or is_compact_code_name(text)
    )


def is_compact_code_name(value: Any) -> bool:
    """Return True for an opaque short code that must not pose as a name."""
    match = _COMPACT_CODE_RE.fullmatch(_identifier_core(value))
    return bool(
        match
        and match.group("prefix") not in _CHEMICAL_DESCRIPTOR_PREFIXES
    )


def is_registry_labeled_name(value: Any) -> bool:
    """True when the source label starts with a supported registry identifier."""
    text = _identifier_core(value)
    return bool(
        _LABELED_ID_RE.match(text)
        or _CAS_RE.fullmatch(text)
        or _ACCESSION_RE.fullmatch(text)
        or is_compact_code_name(text)
    )


def is_readable_compound_name(value: Any) -> bool:
    """Reject missing, noisy, or accession-only names for name-only leads."""
    text = _clean(value)
    return bool(text and len(text) <= 240 and not is_identifier_name(text) and not _NOISE_RE.match(text))


def alias_quality(value: Any) -> tuple[bool, str]:
    """Return whether a generated alias is safe to present as a full compound."""

    text = _clean(value)
    if not text:
        return False, "missing"
    # Long systematic/IUPAC names are still chemical names.  Reject only
    # pathological source payloads; the UI wraps long names and retains the
    # source label verbatim for audit.  Metadata/candidate phrases are checked
    # separately below and remain unsafe regardless of length.
    if len(text) > 2000:
        return False, "too_long"
    if is_identifier_name(text) or _INTERNAL_ALIAS_RE.fullmatch(text):
        return False, "identifier_or_internal_code"
    if text.casefold() in _GENERIC_FRAGMENT_ALIASES:
        return False, "generic_substituent_or_fragment"
    if _FRAGMENT_ONLY_ALIAS_RE.fullmatch(text):
        return False, "bracketed_fragment_only"
    if "?" in text or "\ufffd" in text:
        return False, "unresolved_stereochemistry_or_corrupt_text"
    if _ALIAS_METADATA_RE.search(text):
        return False, "putative_or_collision_metadata"
    if _NOISE_RE.match(text):
        return False, "noise_label"
    return True, "safe"


def _unresolved_name(inchikey: Any) -> str:
    key = _clean(inchikey).upper()
    token = key.split("-")[0] if re.fullmatch(r"[A-Z]{14}(?:-[A-Z]{10}-[A-Z])?", key) else ""
    return f"Unresolved compound ({token})" if token else "Unresolved library compound"


def _is_unresolved_display(value: Any) -> bool:
    text = _clean(value)
    return text == "Unresolved library compound" or (
        text.startswith("Unresolved compound (") and text.endswith(")")
    )


def _clean_dereplicator_title(value: Any) -> str:
    """Recover only the embedded Dereplicator label, without upgrading evidence.

    GNPS titles commonly duplicate a name as ``Name_Name`` and use underscores
    as spaces.  This parser is deliberately anchored to the exact title prefix:
    derivative records such as ``Suspect related to ...`` and ``Conjugate of
    ...`` remain unresolved rather than being promoted to the parent name.
    """

    text = _identifier_core(value)
    match = _DEREPLICATOR_TITLE_RE.fullmatch(text)
    if not match:
        return ""
    payload = match.group("name").strip(" _!\"'“”‘’")
    # Decode exporter HTML entities (some records are double encoded).
    for _ in range(2):
        decoded = unescape(payload)
        if decoded == payload:
            break
        payload = decoded
    # Collapse an exact duplicated payload before underscores become spaces.
    for separator_index, character in enumerate(payload):
        if character == "_" and payload[:separator_index] == payload[separator_index + 1 :]:
            payload = payload[:separator_index]
            break
    return re.sub(r"_+", " ", payload).strip(" _!\"'“”‘’")


@lru_cache(maxsize=1)
def load_compound_name_aliases() -> dict[str, str]:
    """Load generated aliases plus the small manually verified override file.

    The generated file can be rebuilt at any time.  Reading the curated file
    separately prevents a rebuild from silently dropping corrections for a
    chemically readable but wrong source name (for example ``AO1680``).
    Curated entries win, and all matches remain exact full-InChIKey matches.
    """
    merged: dict[str, str] = {}
    for path in (
        _ALIASES_PATH,
        _PUBCHEM_ALIASES_PATH,
        _CURATED_ALIASES_PATH,
    ):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            continue
        aliases = payload.get("aliases", payload) if isinstance(payload, dict) else {}
        if not isinstance(aliases, dict):
            continue
        for key, value in aliases.items():
            clean_key = _clean(key).upper()
            # Alias files are assembled from heterogeneous public-source titles.
            # Some otherwise exact InChIKey aliases still carry a spectrum ion
            # suffix (for example ``Schaftoside [M+H]+``).  An exact key must not
            # bypass the same presentation cleanup applied to source titles:
            # compound identity and observed/library ion form are separate data.
            clean_value = _identifier_core(value)
            safe, _ = alias_quality(clean_value)
            if clean_key and clean_value and safe:
                merged[clean_key] = clean_value
    return merged


def _collapse_duplicate_aliases(text: str) -> str:
    """Collapse GNPS dereplicator ``Name|Name`` / ``Name_Name`` artefacts.

    Only collapses when the separated halves match, so a genuine name that
    happens to contain ``_`` (rare in chemical names) is preserved.
    """
    for sep in ("|", ";"):
        if sep in text:
            first = text.split(sep, 1)[0].strip(" _!\"'“”‘’")
            if first:
                return first
    if "_" in text:
        # ``X_Y`` where ``Y == X`` (X itself may contain underscores), e.g.
        # ``Citrusin_VI_Citrusin_VI`` -> ``Citrusin_VI``.
        idx = 0
        while True:
            nxt = text.find("_", idx)
            if nxt == -1:
                break
            prefix, rest = text[:nxt], text[nxt + 1 :]
            if prefix and rest == prefix:
                return prefix
            idx = nxt + 1
    return text


def _computational_candidate_label(inchikey: Any) -> str:
    """Honest placeholder for an unconfirmed computational formula annotation."""
    key = _clean(inchikey).upper()
    token = (
        key.split("-")[0]
        if re.fullmatch(r"[A-Z]{14}(?:-[A-Z]{10}-[A-Z])?", key)
        else ""
    )
    return (
        f"Computational candidate ({token})" if token else "Computational candidate"
    )


def _extract_name_from_provenance(
    raw_name: Any, inchikey: Any
) -> tuple[str, str] | None:
    """Recover a real chemical name wrapped in library/dereplicator provenance.

    Returns ``(display, status)`` when a readable name can be safely extracted,
    a computational-candidate placeholder for unconfirmed SIRIUS-style
    annotations, or ``None`` when no safe name is recoverable.
    """
    text = _identifier_core(raw_name)
    if not text:
        return None
    if _COMPUTATIONAL_CANDIDATE_RE.search(text):
        return _computational_candidate_label(inchikey), "computational_candidate"
    text = _PROVENANCE_PREFIX_RE.sub("", text, count=1)
    text = _PAREN_PROVENANCE_RE.sub(" ", text)
    text = _TRAILING_CAVEAT_RE.sub("", text)
    text = _collapse_duplicate_aliases(text)
    text = _identifier_core(text).strip(" _!\"'“”‘’")
    if not text or _is_unresolved_display(text):
        return None
    if is_identifier_name(text) or is_molecular_formula_name(text):
        return None
    safe, _ = alias_quality(text)
    if not safe:
        return None
    if _COMPUTATIONAL_CANDIDATE_RE.search(text):
        return _computational_candidate_label(inchikey), "computational_candidate"
    lowered = _clean(raw_name).lower()
    if "putative" in lowered or "candidate" in lowered:
        return text, "putative_extracted"
    if "dereplicator" in lowered or "spectral" in lowered or "library" in lowered:
        return text, "library_match_extracted"
    return text, "metadata_extracted"


def _extract_massbank_multi_name(
    raw_name: Any, inchikey: Any
) -> tuple[str, str] | None:
    """Recover one chemical name from a MassBank '|'-concatenated title.

    MassBank records carry several CH$NAME fields (preferred / class /
    systematic); the library importer joins them with '|'.  The whole string is
    readable enough to pass ``alias_quality``, so it would otherwise display
    verbatim.  Strip the ``Massbank[X]:<accession>`` prefix and keep the most
    specific readable segment.  Returns ``(display, status)`` or ``None`` when
    the title is not this format (so non-MassBank names are untouched).
    """
    text = _identifier_core(raw_name)
    if "|" not in text or not _MASSBANK_PREFIX_RE.match(text):
        return None
    body = _MASSBANK_PREFIX_RE.sub("", text, count=1).strip(" _!\"'“”‘’")
    segments = [seg.strip(" _!\"'“”‘’") for seg in body.split("|")]
    safe_segs: list[str] = []
    for seg in segments:
        if not seg or is_identifier_name(seg) or is_molecular_formula_name(seg):
            continue
        ok, _ = alias_quality(seg)
        if ok:
            safe_segs.append(seg)
    if not safe_segs:
        return None
    # MassBank lists names from common -> systematic; the last readable segment
    # is the most specific (systematic/IUPAC) name and is always a real compound.
    best = safe_segs[-1]
    lowered = text.lower()
    if "tentative" in lowered or "putative" in lowered:
        return best, "putative_extracted"
    return best, "massbank_multi_name"


def resolve_compound_name(raw_name: Any, inchikey: Any = "") -> tuple[str, str]:
    """Return ``(display_name, status)`` while retaining deterministic provenance."""
    raw = _clean(raw_name)
    key = _clean(inchikey).upper()
    alias = load_compound_name_aliases().get(key)
    if alias:
        # Defensive normalization also protects long-running processes whose
        # alias cache was populated before the loader-level cleanup above.
        cleaned_alias = _identifier_core(alias)
        safe, _ = alias_quality(cleaned_alias)
        if cleaned_alias and safe:
            return cleaned_alias, "exact_inchikey_alias"

    match = _LABELED_ID_RE.match(_identifier_core(raw))
    if match:
        payload = match.group("payload").strip(" _!")
        if match.group("separator") == "_" and payload:
            payload = _FORMULA_PREFIX_RE.sub("", payload).strip(" _!")
        # MoNA/NCATS labels often append a CAS number directly to a readable name,
        # e.g. ``MLS...!Sinomenine115-53-7``.  The CAS suffix is provenance, not
        # part of the chemical name.  A payload that is only a registry number is
        # still unresolved and must not leak into the display.
        payload = _TRAILING_CAS_RE.sub("", payload).strip(" _!")
        if payload:
            if is_identifier_name(payload):
                return _unresolved_name(key), "identifier_unresolved"
            safe, reason = alias_quality(payload)
            if not safe:
                return _unresolved_name(key), f"unsafe_registry_payload_{reason}"
            return payload, "registry_payload"
        return _unresolved_name(key), "identifier_unresolved"

    # MoNA records imported from NIST-derived collections can store a search
    # provenance sentence as the compound title, for example
    # ``Spectral Match to Dibutyl phthalate from NIST14 M+H``.  The trailing
    # adduct is removed by ``_identifier_core`` above; retain the full raw title
    # in ``source_compound_name`` but show only the chemical name to readers.
    spectral_match = _SPECTRAL_MATCH_TITLE_RE.match(_identifier_core(raw))
    if spectral_match:
        display = spectral_match.group("name").strip(" _!\"'“”‘’")
        if display and not is_identifier_name(display):
            safe, reason = alias_quality(display)
            if safe:
                library = spectral_match.group("library").casefold()
                status = (
                    "nist_spectral_match_title"
                    if library.startswith("nist")
                    else "metlin_spectral_match_title"
                )
                return display, status
            return _unresolved_name(key), f"unsafe_spectral_title_{reason}"

    # Some GNPS/MoNA exports wrap a valid chemical title and the collision
    # energy in quotes, followed by a free-form adduct (including ``2M+Na``).
    # Parse the anchored exporter format instead of allowing the metadata
    # sentence to turn an otherwise readable name into an InChIKey placeholder.
    collision_title = (
        _QUOTED_COLLISION_TITLE_RE.fullmatch(raw)
        or _UNQUOTED_COLLISION_TITLE_RE.fullmatch(raw)
    )
    if collision_title:
        display = collision_title.group("name").strip(" _!\"'“”‘’")
        safe, reason = alias_quality(display)
        if safe:
            return display, "collision_energy_title"
        return _unresolved_name(key), f"unsafe_collision_title_{reason}"

    # DEREPLICATOR is an annotation method, not identity confirmation.  Show
    # its embedded human-readable label while retaining a distinct status and
    # the untouched source title for provenance.  Downstream confidence gates
    # therefore remain free to treat the annotation as tentative.
    dereplicator_title = _clean_dereplicator_title(raw)
    if dereplicator_title:
        safe, reason = alias_quality(dereplicator_title)
        if safe:
            return dereplicator_title, "dereplicator_annotation_title"
        return _unresolved_name(key), f"unsafe_dereplicator_title_{reason}"

    # MassBank imports join accession + several name fields with '|'.  Recover
    # the single most-specific readable name instead of showing the concatenation.
    massbank_multi = _extract_massbank_multi_name(raw, key)
    if massbank_multi is not None:
        return massbank_multi

    identifier_core = _identifier_core(raw)
    if (
        _PURE_NUMERIC_RE.fullmatch(identifier_core)
        or _CAS_RE.fullmatch(identifier_core)
        or _ACCESSION_RE.fullmatch(identifier_core)
        or _INCHIKEY_NAME_RE.fullmatch(identifier_core)
        or is_molecular_formula_name(identifier_core)
    ):
        return _unresolved_name(key), "identifier_unresolved"
    if is_compact_code_name(raw):
        return _unresolved_name(key), "code_identifier_unresolved"
    if not raw:
        return _unresolved_name(key), "missing"
    cleaned = _identifier_core(raw)
    safe, reason = alias_quality(cleaned)
    if not safe:
        # Source tools (GNPS dereplicator, spectral-library matches) wrap a real
        # reference name in provenance/caveat metadata.  Recover that name before
        # falling back to the unresolved placeholder; unconfirmed SIRIUS-style
        # computational candidates stay flagged, never posed as confirmed names.
        if reason == "putative_or_collision_metadata":
            extracted = _extract_name_from_provenance(raw, key)
            if extracted is not None:
                return extracted
        return _unresolved_name(key), f"unsafe_source_name_{reason}"
    return cleaned, "source_name"   # 剥尾部 adduct/collision-energy (如 "X M+H"→"X")


def display_compound_name(raw_name: Any, inchikey: Any = "") -> str:
    return resolve_compound_name(raw_name, inchikey)[0]


def normalize_match_compound_names(match: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(match, dict):
        return match

    replacement_candidates: dict[str, set[str]] = {}
    display_names: set[str] = set()
    visited: set[int] = set()
    visited_items: list[dict[str, Any]] = []

    def _normalize(item: dict[str, Any]) -> None:
        current = _clean(item.get("compound_name"))
        source = _clean(item.get("source_compound_name"))
        # A previous presentation placeholder is not source provenance.  Keep
        # genuine registry/library labels, but remove cached
        # ``Unresolved compound (...)`` strings so they cannot leak through
        # tooltips after an exact structure alias becomes available.
        if _is_unresolved_display(source):
            item.pop("source_compound_name", None)
            source = ""
        raw = source or current
        display, status = resolve_compound_name(raw, item.get("inchikey"))
        item["name_resolution_status"] = status
        if display:
            display_names.add(display)
        if display and display != current:
            if current:
                replacement_candidates.setdefault(current, set()).add(display)
            if source or (raw and not _is_unresolved_display(raw)):
                item["source_compound_name"] = raw
            item["compound_name"] = display
        if display and raw and display != raw:
            replacement_candidates.setdefault(raw, set()).add(display)

    def _walk(item: dict[str, Any]) -> None:
        marker = id(item)
        if marker in visited:
            return
        visited.add(marker)
        visited_items.append(item)
        _normalize(item)
        for key in (
            "top_candidates",
            "isomer_set",
            "same_formula_alternatives",
        ):
            for child in item.get(key) or []:
                if isinstance(child, dict):
                    _walk(child)

    _walk(match)

    # Legacy task JSON also cached an honesty sentence assembled from the old
    # raw names.  Normalising only ``compound_name`` left identifiers visible in
    # the subtitle/tooltips even though the main name was clean.  Apply all
    # unambiguous replacements in one regex pass: sequential ``str.replace``
    # used to rewrite an identifier inside the replacement it had just
    # inserted, producing ``Unresolved compound (Unresolved compound (...))``
    # whenever cached results were prepared more than once.
    replacements = {
        old: next(iter(values))
        for old, values in replacement_candidates.items()
        if old and len(values) == 1 and old not in values
    }
    pattern = (
        re.compile(
            "|".join(
                re.escape(old)
                for old in sorted(replacements, key=len, reverse=True)
            )
        )
        if replacements else None
    )
    for item in visited_items:
        headline = _clean(item.get("honest_headline"))
        if not headline:
            continue
        if pattern is not None:
            headline = pattern.sub(
                lambda found: replacements[found.group(0)],
                headline,
            )
        # Repair legacy strings that were already nested before this fix.
        for display in sorted(display_names, key=len, reverse=True):
            wrapper = f"Unresolved compound ({display})"
            while wrapper in headline:
                headline = headline.replace(wrapper, display)
        item["honest_headline"] = headline
    return match


def normalize_result_compound_names(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    for result in results or []:
        if isinstance(result, dict):
            normalize_match_compound_names(result.get("match"))
    return results
