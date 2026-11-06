/**
 * Read-only fallback for completed tasks whose current compound/flagship
 * artifacts cannot be materialized. It projects saved spectrum-level results
 * into the bundle shape consumed by SampleOverview and EvidenceReview without
 * upgrading any saved candidate to a current entity conclusion.
 */
const SavedResultFallback = (() => {
  const PAGE_LIMIT = 5000;

  function savedIdentity(match, row) {
    const supplied = String(match?.identity_key || '').trim();
    if (supplied) return supplied;
    const inchikey = String(match?.inchikey || '')
      .trim().toUpperCase().replace(/^INCHIKEY=/, '').replace(/\s+/g, '');
    if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(inchikey)) return `ik::${inchikey}`;
    const smiles = String(match?.smiles || '').trim();
    if (smiles && !['none', 'null', 'nan', 'n/a'].includes(smiles.toLowerCase())) return `sm::${smiles}`;
    const name = String(match?.compound_name || '')
      .trim().toLocaleLowerCase().replace(/[^\p{L}\p{N}]+/gu, '');
    const formula = String(match?.formula || '').trim();
    if (!name && !formula) {
      const index = resultIndex(row);
      return `result::${index == null ? String(row?.spectrum_id || 'unknown') : index}`;
    }
    return `nf::${name}|${formula}`;
  }

  function finiteScore(value, fallback = -Infinity) {
    if (value == null || String(value).trim() === '') return fallback;
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : fallback;
  }

  function resultIndex(row) {
    if (row?.result_index == null || String(row.result_index).trim() === '') return null;
    const parsed = Number(row?.result_index);
    return Number.isInteger(parsed) ? parsed : null;
  }

  function preferRepresentative(candidate, current) {
    if (!current) return true;
    const candidateMatch = candidate.match || {};
    const currentMatch = current.match || {};
    const finalDelta = finiteScore(candidateMatch.final_score)
      - finiteScore(currentMatch.final_score);
    if (finalDelta) return finalDelta > 0;
    const cosineDelta = finiteScore(candidateMatch.cosine)
      - finiteScore(currentMatch.cosine);
    if (cosineDelta) return cosineDelta > 0;
    const candidateIndex = resultIndex(candidate);
    const currentIndex = resultIndex(current);
    if (candidateIndex == null) return false;
    return currentIndex == null || candidateIndex < currentIndex;
  }

  function failureMetadata(error) {
    if (!error) return null;
    const status = Number(error.status ?? error.statusCode ?? error.response?.status);
    return {
      name: String(error.name || 'Error'),
      message: String(error.message || error),
      status: Number.isFinite(status) ? status : null,
      code: error.code == null ? null : String(error.code),
    };
  }

  async function fetchSavedRows(taskId) {
    const rows = [];
    let offset = 0;
    while (true) {
      const payload = await API.get(
        `/tasks/${encodeURIComponent(taskId)}/results?projection=summary&limit=${PAGE_LIMIT}&offset=${offset}`
      );
      const pageRows = Array.isArray(payload) ? payload
        : Array.isArray(payload?.results) ? payload.results : [];
      rows.push(...pageRows);
      const hasMore = Array.isArray(payload) ? pageRows.length === PAGE_LIMIT
        : payload?.has_more === true;
      if (!hasMore) break;
      if (!pageRows.length) {
        throw new Error('Saved result pagination stopped before all pages were returned.');
      }
      offset += pageRows.length;
    }
    return rows;
  }

  function groupSavedCandidates(rows) {
    const byIdentity = new Map();
    const matched = [];
    const unmatched = [];
    rows.forEach((sourceRow) => {
      const row = { ...sourceRow };
      if (sourceRow?.match && typeof sourceRow.match === 'object') {
        row.match = { ...sourceRow.match };
      }
      if (!row.match) {
        unmatched.push(row);
        return;
      }
      matched.push(row);
      const key = savedIdentity(row.match, row);
      let entry = byIdentity.get(key);
      if (!entry) {
        entry = { key, rows: [], rep: null };
        byIdentity.set(key, entry);
      }
      entry.rows.push(row);
      if (preferRepresentative(row, entry.rep)) entry.rep = row;
    });

    const groups = [...byIdentity.values()].map((entry, index) => {
      const rep = entry.rep;
      const match = rep.match || {};
      const candidateNames = [...new Set(entry.rows
        .map(row => String(row.match?.compound_name || '').trim())
        .filter(Boolean))];
      const status = {
        key: 'saved_candidate',
        conclusionKey: 'saved_candidate',
        label: '已保存谱级候选',
        note: '当前实体结论不可用；本组仅按保存候选身份归并谱级结果。',
        hasHigh: false,
        hasIsomerHigh: false,
        countedAsSampleCompound: false,
        countedAsPlantCandidate: false,
      };
      const compound = {
        compound_id: `saved::${entry.key}`,
        identity_key: entry.key,
        record_type: 'saved_spectrum_candidate',
        display_name: String(match.compound_name || 'Unnamed saved candidate'),
        candidate_names: candidateNames,
        spectra_count: entry.rows.length,
        representative: rep,
        conclusion: {
          key: 'saved_candidate',
          label: status.label,
          note: status.note,
          counted_as_sample_compound: false,
          counted_as_plant_candidate: false,
        },
      };
      return {
        key: compound.compound_id,
        compound,
        conclusion: compound.conclusion,
        rep,
        rows: entry.rows,
        spectra: entry.rows.length,
        peakClusters: null,
        leadCandidate: match,
        candidateMembers: [match],
        candidateMemberCount: 1,
        candidateNames,
        sampleEvidence: {},
        status,
        index,
        entity: null,
        flagship: null,
      };
    });

    groups.sort((left, right) => (
      finiteScore(right.rep?.match?.final_score) - finiteScore(left.rep?.match?.final_score)
      || (resultIndex(left.rep) ?? Number.MAX_SAFE_INTEGER)
        - (resultIndex(right.rep) ?? Number.MAX_SAFE_INTEGER)
    ));
    groups.forEach((group, index) => { group.index = index; });
    return { groups, matched, unmatched };
  }

  async function load(taskId, task = {}, error = null) {
    if (!taskId) throw new Error('A task id is required to load saved results.');
    const rows = await fetchSavedRows(taskId);
    const grouped = groupSavedCandidates(rows);
    const sourceFailure = failureMetadata(error);
    const sourceWarning = '当前实体结论不可用，先查看已保存谱级候选';
    const metadata = {
      saved_only: true,
      conclusion_source: 'saved_spectrum_results',
      source_warning: sourceWarning,
      source_failure: sourceFailure,
      saved_result_count: rows.length,
      saved_matched_spectrum_count: grouped.matched.length,
      saved_unmatched_spectrum_count: grouped.unmatched.length,
      saved_candidate_group_count: grouped.groups.length,
      main_result_count: 0,
      high_confidence_scope_group_count: 0,
      high_confidence_2d_connectivity_count: 0,
      high_confidence_isomer_group_count: 0,
      confirmed_count: 0,
      unresolved_precursor_event_count: 0,
    };
    const compounds = {
      metadata,
      compounds: grouped.groups.map(group => group.compound),
      governance_blocked_candidates: [],
      governance_blocked_count: 0,
      governance_blocked_has_more: false,
    };
    return {
      task: { ...((task && typeof task === 'object') ? task : {}) },
      compounds,
      governanceBlockedCandidates: [],
      governanceBlockedCount: 0,
      governanceBlockedHasMore: false,
      evitwin: null,
      confidenceGuide: null,
      replicateConsistency: { status: 'unavailable' },
      identificationRecords: [],
      matched: grouped.matched,
      unmatched: grouped.unmatched,
      unannotatedRecords: grouped.unmatched,
      flagship: { metadata: { saved_only: true, source_warning: sourceWarning }, entities: [] },
      legacyGroups: grouped.groups,
      featureGroups: grouped.groups,
      groups: grouped.groups,
      savedOnly: true,
      sourceWarning,
      sourceFailure,
    };
  }

  return { load };
})();

window.SavedResultFallback = SavedResultFallback;
