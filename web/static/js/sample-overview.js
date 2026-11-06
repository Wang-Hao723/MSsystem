/**
 * User-facing sample result overview.
 *
 * This layer answers the first questions a sample owner asks: which compounds
 * were detected, what can actually be claimed, and what remains unresolved.
 * Scan/entity artifacts stay available in the detail view, but are not used as
 * the customer-facing result unit.
 */
const SampleResultData = (() => {
  const cache = new Map();

  function interfaceLanguage() {
    return window.MSSLanguage?.getLanguage?.() === 'en' ? 'en' : 'zh';
  }

  function bundleCacheKey(taskId, language = interfaceLanguage()) {
    return `${taskId}::${language}`;
  }

  const number = value => {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : 0;
  };

  function retentionTimeMinutes(row) {
    const parsed = Number(row?.retention_time);
    if (!Number.isFinite(parsed)) return null;
    const unit = String(row?.retention_time_unit || 'minutes').trim().toLowerCase();
    if (unit.includes('millisecond') || ['ms', 'msec', 'msecs'].includes(unit)) return parsed / 60000;
    if (unit.includes('second') || ['s', 'sec', 'secs'].includes(unit)) return parsed / 60;
    if (unit.includes('hour') || ['h', 'hr', 'hrs'].includes(unit)) return parsed * 60;
    return parsed;
  }

  function identityKey(match) {
    const supplied = String(match?.identity_key || '').trim();
    if (supplied) return supplied;
    const key = String(match?.inchikey || '').trim().toUpperCase().replace(/^INCHIKEY=/, '').replace(/\s+/g, '');
    if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(key)) return `ik::${key}`;
    const smiles = String(match?.smiles || '').trim();
    if (smiles) return `sm::${smiles}`;
    const name = String(match?.compound_name || '').toLowerCase().replace(/[^a-z0-9]+/g, '');
    return `nf::${name}|${String(match?.formula || '').trim()}`;
  }

  function evidenceStatus(rows) {
    const matches = rows.map(row => row.match || {});
    const hasHigh = matches.some(match => (
      String(match.two_dimensional_connectivity_confidence || match.confidence || '') === 'High'
      && String(match.identification_scope || '') !== 'isomer_group'
    ));
    const hasIsomerHigh = matches.some(match => String(match.identification_scope || '') === 'isomer_group'
      && String(match.isomer_group_confidence || '') === 'High');
    if (hasIsomerHigh) {
      return { key: 'isomer', label: '高可信异构体组检出', note: 'High 级证据支持该分子式异构体组，不能点名单一异构体', hasHigh, hasIsomerHigh };
    }
    if (hasHigh) {
      return { key: 'high', label: '高可信单一二维连接结构检出', note: '单一二维连接结构具有 High 级谱库支持；完整立体身份仍按证据边界报告', hasHigh, hasIsomerHigh };
    }
    if (matches.some(match => match.confidence === 'Medium' || match.evidence_tier === 'Strong Candidate')) {
      return { key: 'medium', label: '中等证据候选', note: '具有中等谱库支持，身份范围维持为候选', hasHigh, hasIsomerHigh };
    }
    return { key: 'low', label: '证据不足', note: '仅作为待核查谱库线索，不进入主要结论', hasHigh, hasIsomerHigh };
  }

  function statusFromConclusion(conclusion) {
    const raw = String(conclusion?.key || 'insufficient');
    const savedEntityGrade = String(conclusion?.entity_confidence
      || conclusion?.entity_evidence_tier || conclusion?.entity_grade || '').trim().toLowerCase();
    const entityWasEvaluated = conclusion?.entity_evaluated === true
      || String(conclusion?.entity_evaluation_status || '').toLowerCase() === 'evaluated';
    const explicitMedium = raw === 'candidate'
      && savedEntityGrade === 'medium' && entityWasEvaluated;
    const known = new Set([
      'confirmed', 'high', 'isomer', 'background', 'unresolved_event',
      'unannotated', 'saved_candidate', 'low',
    ]);
    const key = explicitMedium ? 'medium'
      : (raw === 'candidate' ? 'candidate'
        : (raw === 'insufficient' ? 'low' : (known.has(raw) ? raw : 'ungraded')));
    const defaultLabel = key === 'candidate' ? '保存候选（实体未分级）'
      : (key === 'ungraded' ? '未分级' : '证据不足');
    const rawLabel = String(conclusion?.label || '').trim();
    return {
      key,
      conclusionKey: raw,
      label: key === 'candidate' && /中等|medium/i.test(rawLabel)
        ? defaultLabel : (rawLabel || defaultLabel),
      note: key === 'candidate'
        ? '已保存候选，未保存经过实体级评估的 Medium 等级'
        : String(conclusion?.note || ''),
      hasHigh: ['confirmed', 'high'].includes(raw),
      hasIsomerHigh: raw === 'isomer'
        && !!conclusion?.high_confidence_group_detected,
      countedAsSampleCompound: !!conclusion?.counted_as_sample_compound,
      countedAsPlantCandidate: !!(conclusion?.counted_as_plant_candidate
        ?? conclusion?.counted_as_sample_compound),
    };
  }

  function candidateIdentityTokens(candidate, identityScope = '') {
    const tokens = new Set();
    const supplied = String(candidate?.identity_key || '').trim();
    if (supplied) {
      tokens.add(`identity:${supplied}`);
      const connectivity = supplied.match(/(?:std2d:|ik::)?([A-Z]{14})(?:[-|:]|$)/i);
      if (connectivity) tokens.add(`connectivity:${connectivity[1].toUpperCase()}`);
    }
    const inchikey = String(candidate?.inchikey || '').trim().toUpperCase().replace(/^INCHIKEY=/, '');
    if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(inchikey)) {
      tokens.add(`inchikey:${inchikey}`);
      if (/2d|connectivity|isomer_group|std2d/i.test(identityScope)) {
        tokens.add(`connectivity:${inchikey.slice(0, 14)}`);
      }
    }
    const smiles = String(candidate?.smiles || '').trim();
    if (smiles) tokens.add(`smiles:${smiles}`);
    const name = String(candidate?.compound_name || '').trim().toLocaleLowerCase();
    const formula = String(candidate?.formula || '').trim();
    if (name && formula) tokens.add(`name-formula:${name}|${formula}`);
    return tokens;
  }

  function sameSavedRepresentativeCandidate(compound, rep, rawMatch, leadCandidate) {
    const repIndex = rep?.result_index;
    const storedRepIndex = compound?.representative_result_index;
    const leadIndex = leadCandidate?.result_index ?? leadCandidate?.representative_result_index;
    if (storedRepIndex != null && repIndex != null
        && Number(storedRepIndex) !== Number(repIndex)) return false;
    if (leadIndex != null && repIndex != null && Number(leadIndex) !== Number(repIndex)) return false;
    const leadSpectrum = String(leadCandidate?.spectrum_id || '').trim();
    const repSpectrum = String(rep?.spectrum_id || '').trim();
    if (leadSpectrum && repSpectrum && leadSpectrum !== repSpectrum) return false;
    const scope = String(compound?.identity_scope || compound?.identity_key || '');
    const representativeTokens = candidateIdentityTokens({
      ...rawMatch,
      identity_key: rawMatch?.identity_key || compound?.identity_key,
      formula: rawMatch?.formula || compound?.formula,
    }, scope);
    const leadTokens = candidateIdentityTokens(leadCandidate, scope);
    const structural = token => /^(identity|connectivity|inchikey|smiles):/.test(token);
    const representativeStructural = [...representativeTokens].filter(structural);
    const leadStructural = [...leadTokens].filter(structural);
    if (representativeStructural.length && leadStructural.length) {
      return leadStructural.some(token => representativeStructural.includes(token));
    }
    return [...leadTokens].some(token => representativeTokens.has(token));
  }

  function fullInchikey(value) {
    const key = String(value || '').trim().toUpperCase().replace(/^INCHIKEY=/, '');
    return /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(key) ? key : '';
  }

  function boundSpectralEvidence(rep, rawMatch, leadCandidate) {
    const binding = rawMatch?.spectral_evidence_binding || {};
    if (binding.schema !== 'same-scan-selected-reference/1'
        || binding.status !== 'bound_saved_reference') return null;
    if (String(binding.spectrum_id || '') !== String(rep?.spectrum_id || '')) return null;
    if (binding.result_index == null || rep?.result_index == null
        || Number(binding.result_index) !== Number(rep.result_index)) return null;
    const bindingKey = fullInchikey(binding.candidate_inchikey);
    const currentKeys = [rawMatch?.inchikey, leadCandidate?.inchikey]
      .map(fullInchikey).filter(Boolean);
    if (!bindingKey || !currentKeys.length || currentKeys.some(key => key !== bindingKey)) return null;
    const bindingReference = String(binding.lib_spectrum_id || '').trim();
    const currentReferences = [rawMatch?.lib_spectrum_id, leadCandidate?.lib_spectrum_id]
      .map(value => String(value || '').trim()).filter(Boolean);
    if (!bindingReference || !currentReferences.length
        || currentReferences.some(value => value !== bindingReference)) return null;
    return binding;
  }

  function projectBoundLeadScores(compound, rep, rawMatch, leadCandidate) {
    const projected = { ...rawMatch };
    const spectralBinding = boundSpectralEvidence(rep, rawMatch, leadCandidate);
    if (spectralBinding) {
      const copied = [];
      for (const field of ['cosine', 'mass_ppm']) {
        delete projected[field];
        const value = spectralBinding.metrics?.[field];
        if (value != null && String(value).trim() !== '' && Number.isFinite(Number(value))) {
          projected[field] = Number(value);
          copied.push(field);
        }
      }
      projected.native_spectrum_confidence = String(
        spectralBinding.saved_confidence || projected.native_spectrum_confidence || ''
      );
      if (copied.length) projected.saved_spectral_score_projection = {
        source: 'spectral_evidence_binding.metrics', fields: copied,
        schema: spectralBinding.schema,
        identity_reference_and_representative_binding_verified: true,
      };
      return projected;
    }
    if (rawMatch?.spectral_evidence_binding?.status === 'bound_saved_reference') {
      projected.unbound_saved_display_values = {
        cosine: rawMatch.cosine ?? null,
        mass_ppm: rawMatch.mass_ppm ?? null,
        reason: 'spectral_evidence_binding_failed_frontend_contract',
      };
      delete projected.cosine;
      delete projected.mass_ppm;
      return projected;
    }
    if (!sameSavedRepresentativeCandidate(compound, rep, rawMatch, leadCandidate)) return projected;
    const copied = [];
    for (const field of ['cosine', 'mass_ppm']) {
      const source = Number(leadCandidate?.[field]);
      const hasCurrent = projected[field] != null && String(projected[field]).trim() !== ''
        && Number.isFinite(Number(projected[field]));
      const hasSource = leadCandidate?.[field] != null
        && String(leadCandidate[field]).trim() !== '' && Number.isFinite(source);
      if (!hasCurrent && hasSource) {
        projected[field] = source;
        copied.push(field);
      }
    }
    if (copied.length) projected.saved_lead_score_projection = {
      source: 'lead_candidate', fields: copied,
      identity_and_representative_binding_verified: true,
    };
    return projected;
  }

  const EXTERNAL_COMPOUND_ROLES = new Set([
    'external_contaminant', 'curated_contaminant', 'pesticide_screening_target',
  ]);

  function explicitExternalReview(group) {
    const roles = [
      group?.leadCandidate?.compound_role,
      group?.rep?.match?.compound_role,
      group?.compound?.compound_role,
      group?.compound?.plant_report?.role,
    ].map(value => String(value || '').trim());
    return group?.compound?.plant_report?.excluded_from_plant_candidate_table === true
      || group?.rep?.match?.potential_contaminant === true
      || roles.some(role => EXTERNAL_COMPOUND_ROLES.has(role));
  }

  function localizedCandidateName(member) {
    const name = String(member?.compound_name || '').trim();
    const chinese = String(member?.compound_name_cn || '').trim();
    if (!name) return chinese;
    if (interfaceLanguage() === 'en' || !chinese
        || name.toLocaleLowerCase() === chinese.toLocaleLowerCase()) return name;
    return `${name}（${chinese}）`;
  }

  function groupCompounds(payload) {
    return (payload?.compounds || []).map((compound, index) => {
      const rep = { ...(compound.representative || {}) };
      const rawMatch = rep.match || {};
      const storedCandidateNames = Array.from(new Set(
        (Array.isArray(compound.candidate_names) ? compound.candidate_names : [])
          .map(value => String(value || '').trim()).filter(Boolean)
      ));
      const suppliedMembers = Array.isArray(compound.candidate_members)
        ? compound.candidate_members : [];
      const suppliedMemberNames = Array.from(new Set(
        suppliedMembers
          .map(member => String(member?.compound_name || '').trim())
          .filter(Boolean)
      ));
      // ``candidate_names`` is a broad provenance/search list in some saved
      // artifacts and may include candidates outside the retained competition
      // set.  The structured ``candidate_members`` rows and the server's
      // explicit member count are authoritative for what the report/export
      // calls an isomer-group competitor.
      const candidateNames = suppliedMemberNames.length
        ? suppliedMemberNames
        : storedCandidateNames;
      const leadCandidate = compound.lead_candidate || {};
      let candidateMemberCount = Math.max(
        Number(compound.candidate_member_count || 0),
        Number(rawMatch.isomer_group_size || 0),
        suppliedMembers.length,
      );
      if (!candidateMemberCount) candidateMemberCount = candidateNames.length;
      const isIsomerGroup = String(compound.conclusion?.key || '') === 'isomer';
      const isUnresolvedEvent = String(compound.conclusion?.key || '') === 'unresolved_event'
        || String(compound.record_type || '') === 'unresolved_precursor_event';
      const isUnannotated = String(compound.conclusion?.key || '') === 'unannotated'
        || String(compound.record_type || '') === 'unannotated_ion_feature';
      const leadName = String(leadCandidate.compound_name
        || rawMatch.compound_name || candidateNames[0] || 'Unnamed candidate');
      // Scientific identity and interface prose stay in separate fields.  Old
      // sidecars may contain a translated ``display_name``; never reuse it as
      // the compound identity shown in tables, charts, exports or evidence UI.
      const isomerNames = suppliedMembers.map(localizedCandidateName)
        .filter(Boolean);
      const eventMz = Number(rep.precursor_mz ?? rawMatch.precursor_mz);
      const eventRt = retentionTimeMinutes(rep);
      const eventLocator = [
        Number.isFinite(eventMz) ? `m/z ${eventMz.toFixed(5)}` : '',
        eventRt == null ? '' : `RT ${eventRt.toFixed(3)} min`,
      ].filter(Boolean).join(', ');
      const conciseDisplayName = isUnannotated
        ? `${interfaceLanguage() === 'en'
          ? 'Unannotated ion feature' : '未注释离子特征'}${eventLocator ? `（${eventLocator}）` : ''}`
        : isUnresolvedEvent
        ? `${interfaceLanguage() === 'en'
          ? 'Unresolved precursor event' : '未解析母离子事件'}${eventLocator ? `（${eventLocator}）` : ''}`
        : (isIsomerGroup && isomerNames.length > 1
          ? `${isomerNames.join(' / ')}${interfaceLanguage() === 'en'
            ? ' (isomer group)' : '（异构体组）'}`
          : localizedCandidateName(leadCandidate) || leadName);
      const projectedCompound = {
        ...compound,
        display_name: conciseDisplayName,
        candidate_member_count: candidateMemberCount,
      };
      const match = {
        ...projectBoundLeadScores(compound, rep, rawMatch, leadCandidate),
        identity_key: compound.identity_key,
        compound_display_name: conciseDisplayName,
        auxiliary_interpretation: compound.auxiliary_interpretation
          || rep.match?.auxiliary_interpretation || {},
      };
      rep.match = match;
      const occurrences = Array.isArray(compound.occurrences)
        ? compound.occurrences.map(item => ({ ...item, match })) : [];
      return {
        key: compound.compound_id || compound.identity_key,
        compound: projectedCompound,
        conclusion: compound.conclusion || {},
        rep,
        rows: occurrences.length ? occurrences : [rep],
        spectra: Number(compound.spectra_count || 0),
        peakClusters: Number(compound.peak_cluster_count || occurrences.length || 0),
        leadCandidate,
        candidateMembers: suppliedMembers,
        candidateMemberCount,
        candidateNames,
        sampleEvidence: compound.sample_evidence || {},
        status: statusFromConclusion(compound.conclusion),
        index: Number.isInteger(Number(compound.frozen_compound_index))
          ? Number(compound.frozen_compound_index)
          : (Number.isInteger(Number(compound.index)) ? Number(compound.index) : index),
      };
    });
  }

  function deduplicateCompoundGroups(featureGroups) {
    const statusRank = {
      confirmed: 9, high: 8, isomer: 7, medium: 6,
      candidate: 5, ungraded: 4, low: 3, background: 2,
      unresolved_event: 1, unannotated: 0,
    };
    const finite = (value, fallback = -1) => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed : fallback;
    };
    const memberKey = member => {
      const key = String(member?.inchikey || '').trim().toUpperCase();
      if (key) return `ik:${key}`;
      const smiles = String(member?.smiles || '').trim();
      if (smiles) return `smiles:${smiles}`;
      return `name:${String(member?.compound_name || '').trim().toLocaleLowerCase()}|formula:${String(member?.formula || '').trim()}`;
    };
    const memberIdentityKey = member => {
      const formula = String(member?.formula || '').trim();
      const key = String(member?.inchikey || '').trim().toUpperCase()
        .replace(/^INCHIKEY=/, '');
      if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(key)) {
        return `std2d:${key.slice(0, 14)}|formula:${formula}`;
      }
      const smiles = String(member?.smiles || '').trim();
      if (smiles) return `smiles:${smiles}|formula:${formula}`;
      const name = String(member?.compound_name || '').trim().toLocaleLowerCase();
      return name ? `name:${name}|formula:${formula}` : '';
    };
    const bucketKey = group => {
      if (['unresolved_event', 'unannotated'].includes(group?.status?.key)) {
        return `${group.status.key}:${String(group.compound?.result_id
          || group.compound?.feature_id || group.key || group.index)}`;
      }
      if (group?.status?.key === 'isomer') {
        const members = [...new Set((group.candidateMembers || [])
          .map(memberIdentityKey).filter(Boolean))].sort();
        const formula = String(group.compound?.formula
          || group.rep?.match?.formula || '').trim();
        if (members.length) return `isomer-group:${formula}|${members.join(';')}`;
      }
      return `identity:${String(group.compound?.identity_key
        || group.compound?.compound_id || group.key || group.index)}`;
    };
    const representativePriority = group => {
      const match = group?.rep?.match || {};
      return [
        group?.status?.countedAsPlantCandidate ? 1 : 0,
        statusRank[group?.status?.key] || 0,
        finite(match.rank_score ?? match.final_score),
        finite(match.cosine),
        -Math.abs(finite(match.mass_ppm, 999999)),
      ];
    };
    const betterRepresentative = (left, right) => {
      const a = representativePriority(left);
      const b = representativePriority(right);
      for (let index = 0; index < a.length; index += 1) {
        if (a[index] !== b[index]) return a[index] > b[index];
      }
      return Number(left?.index || 0) < Number(right?.index || 0);
    };
    const buckets = new Map();
    (featureGroups || []).forEach(group => {
      const key = bucketKey(group);
      if (!buckets.has(key)) buckets.set(key, []);
      buckets.get(key).push(group);
    });
    const merged = [...buckets.values()].map(featureRows => {
      let representative = featureRows[0];
      featureRows.slice(1).forEach(group => {
        if (betterRepresentative(group, representative)) representative = group;
      });
      const candidateNames = [...new Set(featureRows.flatMap(
        group => group.candidateNames || []
      ).map(value => String(value || '').trim()).filter(Boolean))];
      const memberMap = new Map();
      featureRows.flatMap(group => group.candidateMembers || []).forEach(member => {
        const key = memberKey(member);
        const current = memberMap.get(key);
        if (!current || finite(member.rank_score ?? member.final_score)
            > finite(current.rank_score ?? current.final_score)) memberMap.set(key, member);
      });
      const occurrenceMap = new Map();
      featureRows.flatMap(group => group.rows || []).forEach(row => {
        const rowKey = [row?.spectrum_id, row?.precursor_mz,
          retentionTimeMinutes(row)].map(value => String(value ?? '')).join('|');
        if (!occurrenceMap.has(rowKey)) occurrenceMap.set(rowKey, row);
      });
      const featureIds = [...new Set(featureRows.map(group => String(
        group.compound?.feature_id || group.rep?.lcms_feature_id
        || group.rep?.precursor_feature_id || ''
      )).filter(Boolean))];
      const linkedFeatureCount = featureIds.filter(value => value.startsWith('lcms:')).length;
      const curatedNatural = featureRows.some(group => (
        group.compound?.plant_report?.curated_natural_role === true
        || group.leadCandidate?.compound_role === 'curated_natural'
      ));
      const sampleEvidence = featureRows.reduce((summary, group) => {
        const evidence = group.sampleEvidence || {};
        Object.entries(evidence.category_counts || {}).forEach(([key, value]) => {
          summary.category_counts[key] = Number(summary.category_counts[key] || 0)
            + Number(value || 0);
        });
        summary.sample_supported_spectra += Number(evidence.sample_supported_spectra || 0);
        summary.background_related_spectra += Number(evidence.background_related_spectra || 0);
        summary.blank_evidence_available = summary.blank_evidence_available
          || evidence.blank_evidence_available === true;
        const ratio = Number(evidence.maximum_sample_to_blank_ratio);
        if (Number.isFinite(ratio)) summary.maximum_sample_to_blank_ratio = Math.max(
          Number(summary.maximum_sample_to_blank_ratio ?? -Infinity), ratio
        );
        return summary;
      }, {
        status: 'sample_supported', category_counts: {},
        sample_supported_spectra: 0, background_related_spectra: 0,
        maximum_sample_to_blank_ratio: null, blank_evidence_available: false,
      });
      if (!sampleEvidence.sample_supported_spectra) sampleEvidence.status = 'background_related';
      return {
        ...representative,
        key: bucketKey(representative),
        compound: {
          ...representative.compound,
          feature_ids: featureIds,
          feature_result_count: featureRows.length,
        },
        rows: [...occurrenceMap.values()],
        spectra: featureRows.reduce((total, group) => total + Number(group.spectra || 0), 0),
        peakClusters: featureRows.reduce(
          (total, group) => total + Number(group.peakClusters || 0), 0
        ),
        precursorEventCount: featureIds.length || featureRows.length,
        featureResultCount: featureRows.length,
        linkedFeatureCount,
        fallbackEventCount: Math.max(0, featureRows.length - linkedFeatureCount),
        curatedNatural,
        candidateNames,
        candidateMembers: [...memberMap.values()],
        candidateMemberCount: Math.max(
          memberMap.size,
          ...featureRows.map(group => Number(group.candidateMemberCount || 0)),
        ),
        sampleEvidence,
        featureGroups: featureRows,
      };
    });
    merged.sort((left, right) => {
      const leftNamed = !['unresolved_event', 'unannotated'].includes(left.status?.key);
      const rightNamed = !['unresolved_event', 'unannotated'].includes(right.status?.key);
      return Number(rightNamed) - Number(leftNamed)
        || Number(right.curatedNatural) - Number(left.curatedNatural)
        || Number(right.linkedFeatureCount || 0) - Number(left.linkedFeatureCount || 0)
        || Number(right.featureResultCount || 0) - Number(left.featureResultCount || 0)
        || (statusRank[right.status?.key] || 0) - (statusRank[left.status?.key] || 0)
        || finite(right.rep?.match?.rank_score ?? right.rep?.match?.final_score)
          - finite(left.rep?.match?.rank_score ?? left.rep?.match?.final_score)
        || String(left.compound?.display_name || '').localeCompare(
          String(right.compound?.display_name || '')
        );
    });
    merged.forEach((group, index) => { group.index = index; });
    return merged;
  }

  function groupHasIsomerCompetition(group) {
    const match = group?.rep?.match || {};
    const finalCompound = group?.compound;
    const hasFinalCompound = finalCompound
      && typeof finalCompound === 'object'
      && Object.keys(finalCompound).length > 0;
    const finalScope = String(finalCompound?.identity_scope || '').toLowerCase();
    if (group?.status?.key === 'isomer'
      || group?.conclusion?.isomer_group_only === true
      || [
        'isomer_group',
        'molecular_formula_isomer_group',
        'stereoisomer_group',
      ].includes(finalScope)) return true;
    if (hasFinalCompound) return false;
    return match.isomer_ambiguous === true
      && Number(group?.candidateMemberCount || 0) > 1;
  }

  function candidateIdentitySummaries(groups, identificationRecords = []) {
    const summaries = new Map();
    const tierRank = {
      confirmed: 8, high: 7, isomer: 6, medium: 5,
      candidate: 4, ungraded: 3, low: 2, background: 1,
    };
    const finite = value => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed : null;
    };
    const identityKeyForMember = member => {
      const formula = String(member?.formula || '').trim();
      const inchikey = String(member?.inchikey || '').trim().toUpperCase().replace(/^INCHIKEY=/, '');
      if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(inchikey)) {
        return `std2d:${inchikey.slice(0, 14)}|formula:${formula}`;
      }
      const name = String(member?.compound_name || '').trim().toLocaleLowerCase();
      return `name:${name}|formula:${formula}`;
    };
    (groups || []).forEach(group => {
      if (!group || ['unannotated', 'unresolved_event'].includes(group.status?.key)) return;
      const supplied = Array.isArray(group.candidateMembers)
        ? group.candidateMembers.filter(member => String(member?.compound_name || '').trim()) : [];
      const fallback = group.leadCandidate && String(group.leadCandidate.compound_name || '').trim()
        ? [group.leadCandidate] : [];
      const members = supplied.length ? supplied : fallback;
      const seenInFeature = new Set();
      members.forEach((member, memberIndex) => {
        const identityKey = identityKeyForMember(member);
        if (seenInFeature.has(identityKey)) return;
        seenInFeature.add(identityKey);
        const statusRank = tierRank[group.status?.key] || 0;
        const score = finite(member.rank_score ?? member.flagship_entity_score ?? member.final_score);
        const cosine = finite(member.cosine);
        const massPpm = finite(member.mass_ppm);
        const rt = retentionTimeMinutes(group.rep || {});
        const featureId = String(group.compound?.feature_id || '');
        const isLead = member.is_ranking_lead === true
          || String(member.role || '') === 'ranking_lead'
          || (!member.role && memberIndex === 0);
        if (!summaries.has(identityKey)) {
          summaries.set(identityKey, {
            identityKey,
            compoundName: String(member.compound_name || '').trim(),
            compoundNameCn: String(member.compound_name_cn || '').trim(),
            formula: String(member.formula || group.compound?.formula || '').trim(),
            inchikey: String(member.inchikey || '').trim(),
            featureResultCount: 0,
            rankingLeadCount: 0,
            ambiguousFeatureCount: 0,
            linkedFeatureCount: 0,
            fallbackEventCount: 0,
            supportingSpectra: 0,
            fragmentEvidenceScanCount: 0,
            bestScore: null,
            bestCosine: null,
            bestEntropy: null,
            bestReverseMatch: null,
            bestAbsMassPpm: null,
            rtMin: null,
            rtMax: null,
            bestTierRank: -1,
            bestStatus: null,
            curatedNatural: false,
            representativeGroup: group,
          });
        }
        const summary = summaries.get(identityKey);
        summary.featureResultCount += 1;
        summary.rankingLeadCount += isLead ? 1 : 0;
        summary.ambiguousFeatureCount += members.length > 1 ? 1 : 0;
        summary.linkedFeatureCount += featureId.startsWith('lcms:') ? 1 : 0;
        summary.fallbackEventCount += featureId.startsWith('lcms:') ? 0 : 1;
        summary.supportingSpectra += Number(group.spectra || 0);
        if (score !== null && (summary.bestScore === null || score > summary.bestScore)) {
          summary.bestScore = score;
          summary.representativeGroup = group;
        }
        if (cosine !== null && (summary.bestCosine === null || cosine > summary.bestCosine)) {
          summary.bestCosine = cosine;
        }
        if (massPpm !== null && (
          summary.bestAbsMassPpm === null || Math.abs(massPpm) < Math.abs(summary.bestAbsMassPpm)
        )) summary.bestAbsMassPpm = massPpm;
        if (rt !== null) {
          summary.rtMin = summary.rtMin === null ? rt : Math.min(summary.rtMin, rt);
          summary.rtMax = summary.rtMax === null ? rt : Math.max(summary.rtMax, rt);
        }
        if (statusRank > summary.bestTierRank) {
          summary.bestTierRank = statusRank;
          summary.bestStatus = group.status;
        }
        summary.curatedNatural = summary.curatedNatural
          || member.compound_role === 'curated_natural'
          || group.compound?.plant_report?.curated_natural_role === true;
      });
    });
    (identificationRecords || []).forEach(record => {
      const candidate = record?.candidate || {};
      const summary = summaries.get(identityKeyForMember(candidate));
      if (!summary) return;
      const scores = record?.scores || {};
      const cosine = finite(scores.cosine);
      const entropy = finite(scores.entropy);
      const reverseMatch = finite(scores.reverse_match);
      const massPpm = finite(scores.mass_ppm);
      summary.fragmentEvidenceScanCount += 1;
      if (cosine !== null && (summary.bestCosine === null || cosine > summary.bestCosine)) {
        summary.bestCosine = cosine;
      }
      if (entropy !== null && (summary.bestEntropy === null || entropy > summary.bestEntropy)) {
        summary.bestEntropy = entropy;
      }
      if (reverseMatch !== null && (
        summary.bestReverseMatch === null || reverseMatch > summary.bestReverseMatch
      )) summary.bestReverseMatch = reverseMatch;
      if (massPpm !== null && (
        summary.bestAbsMassPpm === null || Math.abs(massPpm) < Math.abs(summary.bestAbsMassPpm)
      )) summary.bestAbsMassPpm = massPpm;
    });
    return [...summaries.values()].sort((left, right) => (
      Number(right.curatedNatural) - Number(left.curatedNatural)
      || right.linkedFeatureCount - left.linkedFeatureCount
      || right.rankingLeadCount - left.rankingLeadCount
      || right.featureResultCount - left.featureResultCount
      || right.bestTierRank - left.bestTierRank
      || Number(right.bestScore ?? -1) - Number(left.bestScore ?? -1)
      || left.compoundName.localeCompare(right.compoundName)
    ));
  }

  async function load(taskId, refresh = false) {
    const language = interfaceLanguage();
    const key = bundleCacheKey(taskId, language);
    if (refresh) cache.delete(key);
    if (!cache.has(key)) {
      cache.set(key, (async () => {
        const [task, confidenceGuide] = await Promise.all([
          API.get(`/tasks/${encodeURIComponent(taskId)}/status`),
          window.ConfidenceGuide
            ? window.ConfidenceGuide.load(language)
            : Promise.resolve(null),
        ]);
        // Read saved decisions first. Browsing a result must not rebuild an
        // identity lane or require remounting the original reference resources.
        const flagship = await API.get(
          `/tasks/${encodeURIComponent(taskId)}/flagship?limit=5000&materialize=false`
        ).catch(error => {
          if (task.completed_batch) throw error;
          return { metadata: {}, entities: [] };
        });
        if (task.completed_batch) {
          while (flagship.has_more) {
            const page = await API.get(`/tasks/${encodeURIComponent(taskId)}/flagship?limit=5000&offset=${flagship.entities.length}&materialize=false`);
            if (!page.entities?.length) throw new Error('原生实体证据分页不完整。');
            flagship.entities.push(...page.entities);
            flagship.has_more = !!page.has_more;
          }
          if (Number.isInteger(flagship.entities_total) && flagship.entities.length !== flagship.entities_total) {
            throw new Error('原生实体证据数量与冻结文件不一致。');
          }
        }
        const compoundUrl = task.completed_batch
          ? `/tasks/${encodeURIComponent(taskId)}/compounds?limit=500&materialize=false&result_view=native`
          : `/tasks/${encodeURIComponent(taskId)}/compounds?limit=2000&materialize=false`;
        const [compounds, identificationRecords, evitwin] = await Promise.all([
          API.get(compoundUrl)
            .catch(async error => {
              if (task.completed_batch) throw error;
              if (typeof SavedResultFallback === 'undefined') throw error;
              return { _savedResultBundle: await SavedResultFallback.load(taskId, task, error) };
            }),
          API.get(`/tasks/${encodeURIComponent(taskId)}/identification-records?limit=5000&materialize=false`)
            .catch(() => ({ records: [] })),
          API.get(`/tasks/${encodeURIComponent(taskId)}/evitwin-ms?materialize=false`)
            .catch(() => null),
        ]);
        if (compounds._savedResultBundle) return compounds._savedResultBundle;
        while (compounds.has_more) {
          const offset = (compounds.compounds || []).length;
          const page = await API.get(`/tasks/${encodeURIComponent(taskId)}/compounds?limit=${task.completed_batch ? 500 : 2000}&offset=${offset}&materialize=false${task.completed_batch ? '&result_view=native' : ''}`);
          if (!page.compounds?.length) throw new Error('原生结论分页不完整，请刷新后重试。');
          compounds.compounds.push(...page.compounds);
          compounds.has_more = !!page.has_more;
        }
        if (task.completed_batch && compounds.compounds.length !== compounds.native_total) {
          throw new Error('原生结论读取数量与冻结文件不一致。');
        }
        if (task.completed_batch && compounds.compounds.some((row, index) => row.frozen_compound_index !== index)) {
          throw new Error('原生结论索引缺失或重复，请刷新后重试。');
        }
        while (compounds.governance_blocked_has_more) {
          const nextOffset = (compounds.governance_blocked_candidates || []).length;
          const page = await API.get(
            `/tasks/${encodeURIComponent(taskId)}/compounds?limit=1&governance_offset=${nextOffset}&materialize=false`
          );
          const rows = page.governance_blocked_candidates || [];
          if (!rows.length) throw new Error('治理未放行候选分页未完成，请刷新结果页。');
          compounds.governance_blocked_candidates.push(...rows);
          compounds.governance_blocked_has_more = !!page.governance_blocked_has_more;
        }
        const replicateConsistency = await API.get(
          `/tasks/${encodeURIComponent(taskId)}/replicate-consistency`
        ).catch(() => ({ status: 'unavailable' }));
        const featureGroups = groupCompounds(compounds);
        const replicateByIdentity = new Map(
          (replicateConsistency?.compounds || []).map(item => [item.identity_key, item])
        );
        featureGroups.forEach(group => {
          group.replicateConsistency = replicateByIdentity.get(
            group.compound?.identity_key
          ) || null;
        });
        const flagshipByResultIndex = new Map();
        (flagship?.entities || []).forEach(entity => {
          (entity?.member_result_indices || []).forEach(index => {
            flagshipByResultIndex.set(Number(index), entity);
          });
        });
        featureGroups.forEach(group => {
          const resultIndex = Number(group?.rep?.result_index);
          const entity = Number.isInteger(resultIndex)
            ? flagshipByResultIndex.get(resultIndex) : null;
          group.entity = entity || null;
          group.flagship = entity
            ? (entity.winner || entity.leading_candidate_when_abstained || null)
            : null;
        });
        const groups = deduplicateCompoundGroups(featureGroups);
        return {
        task,
        compounds,
        governanceBlockedCandidates: Array.isArray(compounds?.governance_blocked_candidates)
          ? compounds.governance_blocked_candidates : [],
        governanceBlockedCount: Number(compounds?.governance_blocked_count || 0),
        governanceBlockedHasMore: !!compounds?.governance_blocked_has_more,
        evitwin,
        confidenceGuide,
        replicateConsistency,
        identificationRecords: identificationRecords?.records || [],
        // The user-facing groups are identity-deduplicated, while the audit
        // export must retain every measured feature supporting an identity.
        matched: featureGroups.map(group => group.rep),
        flagship: flagship || { metadata: task.flagship_summary || {}, entities: [] },
        legacyGroups: featureGroups,
        featureGroups,
        groups,
        };
      })().catch(error => {
        cache.delete(key);
        throw error;
      }));
    }
    return cache.get(key);
  }

  function clear(taskId) {
    if (!taskId) {
      cache.clear();
      return;
    }
    const prefix = `${taskId}::`;
    Array.from(cache.keys()).forEach(key => {
      if (key.startsWith(prefix)) cache.delete(key);
    });
  }

  function confidenceAxes(bundle) {
    const validCount = value => value !== null && value !== undefined && value !== ''
      && Number.isInteger(Number(value)) && Number(value) >= 0 ? Number(value) : null;
    const raw = bundle.task?.compound_confidence_summary || {};
    const single = validCount(raw.representative_2d_connectivity_high_unique_compounds);
    const isomer = validCount(raw.representative_isomer_group_high_unique_compounds);
    const supplied = validCount(raw.representative_high_reportable_scope_unique_compounds
      ?? raw.high_reportable_scope_unique_compounds);
    const libraryHigh = supplied ?? (single !== null && isomer !== null ? single + isomer : null);
    const flagship = bundle.flagship || {};
    const entities = flagship.entities || [];
    const total = validCount(flagship.entities_total ?? flagship.metadata?.entity_count);
    const completeEntities = total !== null && total === entities.length && !flagship.has_more;
    const matched = entities.filter(entity => Number(entity.candidate_count) > 0);
    const blocked = matched.filter(entity => (
      entity.decision?.blockers || []
    ).includes('insufficient_distinct_collision_energies'));
    return {
      libraryHigh, librarySingleHigh: single, libraryIsomerHigh: isomer,
      breakdownConsistent: single !== null && isomer !== null && single + isomer === libraryHigh,
      entityHigh: validCount(bundle.compounds?.metadata?.high_confidence_scope_group_count
        ?? bundle.task?.compound_conclusion_summary?.high_confidence_scope_group_count),
      matchedEntities: completeEntities ? matched.length : null,
      singleEnergyEntities: completeEntities ? matched.filter(
        entity => Number(entity.distinct_collision_energy_bins) === 1
      ).length : null,
      energyBlockedEntities: completeEntities ? blocked.length : null,
    };
  }

  function confidenceAxesHtml(bundle, language = interfaceLanguage()) {
    const a = confidenceAxes(bundle);
    const en = language === 'en';
    const n = value => Number(value).toLocaleString(en ? 'en-US' : 'zh-CN');
    if (a.libraryHigh === null) return '';
    const breakdown = a.breakdownConsistent ? (en
      ? ` (${n(a.librarySingleHigh)} single-connectivity candidates; ${n(a.libraryIsomerHigh)} isomer-group candidates)`
      : `（${n(a.librarySingleHigh)} 个单一连接候选、${n(a.libraryIsomerHigh)} 个异构体组候选）`) : '';
    const energy = a.matchedEntities !== null && a.energyBlockedEntities > 0 ? (en
      ? `${n(a.singleEnergyEntities)} of ${n(a.matchedEntities)} entities with candidates have one energy bin; ${n(a.energyBlockedEntities)} carry an insufficient-energy blocker under the current entity rule.`
      : `${n(a.matchedEntities)} 个有候选前体实体中，${n(a.singleEnergyEntities)} 个只有一个碰撞能量分组；${n(a.energyBlockedEntities)} 个被当前实体规则标记为碰撞能量条件不足。`) : '';
    const explanation = en
      ? 'The single-spectrum High label describes library-match evidence. Entity decisions separately assess competing structures and additional measurement conditions. These counts use different deduplication rules and cannot be added or treated as confirmed compound counts.'
      : '单谱 High 表示谱库匹配证据；实体判定另外检查竞争结构和附加测量条件。两项使用不同去重规则，不能相加或当作确证化合物数。';
    return `<div class="sample-confidence-axes" data-confidence-axes="separate-evidence-and-identity" style="display:block;margin:12px 0;padding:12px 16px;background:#f2f8f7;border-left:3px solid #176b68;font-size:14px;line-height:1.6"><b>${en ? 'Single-spectrum library High' : '单谱谱库 High'}：${n(a.libraryHigh)} ${en ? 'candidate identities' : '个候选身份'}</b><p>${en ? 'Representative spectra meet the saved single-spectrum High rules' : '代表谱达到保存的单谱 High 条件'}${breakdown}${en ? '.' : '。'}</p>${energy ? `<p>${energy}</p>` : ''}<p>${explanation}</p></div>`;
  }

  function entityDecisionNote(entity, language = interfaceLanguage(), highEvidence = null) {
    if (highEvidence?.tier === 'High' && highEvidence?.supporting_spectra?.length) {
      const isomer = highEvidence.scope === 'isomer_group';
      return language === 'en'
        ? (isomer ? 'Same-candidate measured MS/MS supports High for the retained isomer group; individual isomers remain unresolved.'
          : 'Same-candidate measured MS/MS supports 2-D connectivity High; standard confirmation remains separate.')
        : (isomer ? '同一候选的实测 MS/MS 支持异构体组 High；组内单体仍未区分。'
          : '同一候选的实测 MS/MS 支持二维连接候选 High；标准品确证另列。');
    }
    const labels = {
      insufficient_distinct_collision_energies: ['未通过实体 High 的多能量条件', 'Entity High multi-energy condition not met'],
      competitive_same_formula_connectivities: ['保留了多个同分子式连接结构候选', 'Multiple same-formula connectivity candidates retained'],
      competitive_connectivities_have_conflicting_formulae: ['保留候选的分子式存在冲突', 'Retained candidate formulae conflict'],
      competitive_connectivity_formula_assignment_incomplete: ['保留候选的分子式信息不完整', 'Retained candidate formula information incomplete'],
      insufficient_entity_margin_or_energy_coverage: ['候选分数差或能量支持未通过实体条件', 'Candidate margin or energy support below entity conditions'],
      entity_score_below_minimum: ['实体分数未通过阈值', 'Entity score below threshold'],
      all_candidates_failed_candidate_independent_chemistry_gate: ['候选均未通过化学条件检查', 'All candidates failed chemistry checks'],
      candidate_independent_ion_identity_conflicts_with_mass_consistent_library_neutral_loss_adduct: ['实测离子身份与谱库中性丢失加合物存在歧义', 'Observed ion identity and library neutral-loss adduct are ambiguous'],
    };
    return [...new Set(entity?.decision?.blockers || [])].map(reason => (
      labels[reason]?.[language === 'en' ? 1 : 0] || reason
    )).join(language === 'en' ? '; ' : '；');
  }

  function singleSpectrumHighRepresentatives(rows) {
    const finite = (value, fallback = 0) => value !== null && value !== undefined
      && Number.isFinite(Number(value)) ? Number(value) : fallback;
    const sortKey = row => {
      const m = row.match || {};
      return [m.potential_contaminant === true ? 1 : 0,
        -finite(m.rank_score ?? m.final_score), -finite(m.final_score), -finite(m.cosine),
        Math.abs(finite(m.mass_ppm, 999999)), identityKey({ ...m, identity_key: '' }),
        String(m.dataset_source || ''), String(m.lib_spectrum_id || ''), String(row.spectrum_id || '')];
    };
    const compare = (left, right) => {
      const a = sortKey(left); const b = sortKey(right);
      for (let i = 0; i < a.length; i += 1) {
        if (a[i] !== b[i]) return a[i] < b[i] ? -1 : 1;
      }
      return 0;
    };
    const representatives = new Map();
    rows.forEach(row => {
      const m = row?.match;
      if (!m || m.dataset_source === 'iceberg_fallback'
        || ['', 'nan', 'none', 'null'].includes(String(m.smiles || '').trim().toLowerCase())) return;
      const key = identityKey({ ...m, identity_key: '' });
      if (!representatives.has(key) || compare(row, representatives.get(key)) < 0) representatives.set(key, row);
    });
    return [...representatives.values()].sort(compare).filter(row => {
      const m = row.match;
      const scope = String(m.identification_scope || '');
      const isomer = scope === 'isomer_group' || (!scope && m.isomer_ambiguous === true);
      return isomer ? m.isomer_group_confidence === 'High'
        : (m.two_dimensional_connectivity_confidence || m.confidence) === 'High';
    });
  }

  return {
    load, clear, identityKey, evidenceStatus, number, retentionTimeMinutes,
    groupCompounds, deduplicateCompoundGroups, statusFromConclusion,
    explicitExternalReview,
    groupHasIsomerCompetition,
    candidateIdentitySummaries,
    confidenceAxes, confidenceAxesHtml, singleSpectrumHighRepresentatives, entityDecisionNote,
  };
})();

const SampleOverview = (() => {
  let currentTaskId = null;
  let currentBundle = null;
  let currentChromatogram = null;
  let currentExpectedCompounds = null;
  let currentReferenceStandardQC = null;
  let currentReferenceStandardSets = [];
  let currentDualPolarityTasks = [];
  let currentDualPolarityShadow = null;
  let dualPolarityShadowBusy = false;
  let currentFilter = 'all';
  let currentAnalysisFilter = 'all';
  let currentSearch = '';
  let currentVisualView = 'cards';
  let currentVisualPage = 0;
  let currentVisualSort = 'original';
  const esc = value => Tasks.escapeHtml(String(value == null ? '' : value));
  const integer = value => Number(value || 0).toLocaleString();
  const pluralizedCandidate = value => (
    Number(value) === 1 ? 'candidate' : 'candidates'
  );

  function interfaceLanguage() {
    return window.MSSLanguage?.getLanguage?.() === 'en' ? 'en' : 'zh';
  }

  function isEnglishInterface() {
    return interfaceLanguage() === 'en';
  }

  function exactDualPolarity(task) {
    const mode = String(task?.ion_mode || '').trim().toUpperCase();
    return mode === 'POS' || mode === 'NEG' ? mode : '';
  }

  function isCompletedDualPolarityTask(task) {
    return String(task?.status || '').trim().toLowerCase() === 'completed'
      && !!exactDualPolarity(task);
  }

  function eligibleOppositePolarityTasks(currentTask) {
    if (!isCompletedDualPolarityTask(currentTask)) return [];
    const currentMode = exactDualPolarity(currentTask);
    const oppositeMode = currentMode === 'POS' ? 'NEG' : 'POS';
    return currentDualPolarityTasks.filter(task => (
      String(task?.id || '') !== String(currentTask?.id || '')
      && isCompletedDualPolarityTask(task)
      && exactDualPolarity(task) === oppositeMode
    ));
  }

  function unresolvedEventKind(group) {
    const reason = String(group?.compound?.unresolved_event_reason || '');
    const statuses = Array.isArray(group?.compound?.precursor_event_resolution_statuses)
      ? group.compound.precursor_event_resolution_statuses : [];
    return reason === 'incomplete_candidate_formula_assignment'
      || statuses.includes('formula_assignment_incomplete')
      ? 'incomplete_formula'
      : 'formula_conflict';
  }

  function englishInterfaceCopy(value) {
    const source = String(value || '').trim();
    if (!source) return '';
    let translated = window.MSS_I18N_CATALOG?.legacy?.['en-US']?.[source]
      || window.MSSLanguage?.translateText?.(source)
      || source;
    if (!/[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]/.test(translated)) {
      return translated;
    }
    const pathwayTerms = [
      ['莽草酸途径', 'Shikimate pathway'],
      ['莽草酸', 'shikimate'],
      ['醋酸-丙二酸途径', 'Acetate–malonate pathway'],
      ['乙酸-丙二酸途径', 'Acetate–malonate pathway'],
      ['氨基酸途径', 'Amino-acid pathway'],
      ['糖酵解相关途径', 'Glycolysis-related pathway'],
      ['其他生源背景', 'Other biosource context'],
    ];
    for (const [sourceTerm, targetTerm] of pathwayTerms) {
      translated = translated.split(sourceTerm).join(targetTerm);
    }
    return translated;
  }

  function score(value, digits = 4) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
  }

  function ppm(value) {
    const parsed = Number(value);
    if (!Number.isFinite(parsed)) return '—';
    return `${parsed > 0 ? '+' : ''}${parsed.toFixed(2)} ppm`;
  }

  function classificationPresentation(match) {
    const english = isEnglishInterface();
    const display = String(match?.display_class || '').trim();
    const raw = String(match?.classification || '').trim();
    const displayUsable = display && !['unknown', 'other', 'unclassified'].includes(display.toLowerCase());
    const rawUsable = raw && !['unknown', 'other', 'unclassified'].includes(raw.toLowerCase());
    const label = displayUsable ? display : (rawUsable ? raw : '');
    if (!label) {
      const identityHigh = [
        match?.confidence,
        match?.exact_structure_confidence,
        match?.two_dimensional_connectivity_confidence,
      ].some(value => String(value || '').trim().toLowerCase() === 'high');
      return {
        label: english
          ? (identityHigh
            ? 'No reliable class assignment (library-candidate evidence is High)'
            : 'No reliable class assignment')
          : (identityHigh
            ? '尚无可靠类别归属（谱库候选证据为 High）'
            : '尚无可靠类别归属'),
        raw: '',
        state: 'unresolved',
      };
    }
    const source = String(
      (displayUsable ? match?.display_classification_source : '')
      || match?.classification_source || ''
    ).trim().toLowerCase();
    const provenance = match?.expert_details?.classification_provenance || {};
    const provenanceMatches = provenance?.reliable === true
      && String(provenance?.analysis_class || '').trim().toLowerCase() === label.toLowerCase();
    const trusted = ['curated', 'manual', 'smarts', 'structure', 'rule', 'database_curated']
      .some(prefix => source === prefix || source.startsWith(`${prefix}_`))
      || source.startsWith('runtime_structure_rule');
    const untrusted = !source || ['unknown', 'other', 'existing', 'ml', 'msclass', 'msclass_spectrum', 'fragment', 'fragment_inference', 'spectrum'].includes(source)
      || ['ml_', 'model_', 'msclass_', 'fragment_'].some(prefix => source.startsWith(prefix));
    if (!untrusted && (trusted || match?.classification_reliable === true || provenanceMatches)) {
      return { label, raw: label, state: 'reliable' };
    }
    const modelSource = source === 'ml' || source === 'msclass_spectrum'
      || source.startsWith('ml_') || source.startsWith('model_') || source.startsWith('msclass_');
    const readableLabel = english ? label.replace(/_/g, ' ') : label;
    return {
      label: english
        ? (modelSource
          ? `Class requires review (model suggestion: ${readableLabel})`
          : `Class lead requiring review: ${readableLabel}`)
        : (modelSource
          ? `类别待核查（模型提示：${label}）`
          : `类别线索（待核查）：${label}`),
      raw: label,
      state: 'tentative',
    };
  }

  function classificationLabel(match) {
    return classificationPresentation(match).label;
  }

  function rawClassification(match) {
    return classificationPresentation(match).raw;
  }

  function biosourceReasonEnglish(analysis) {
    const code = String(analysis?.biosourceReason || '').trim();
    const labels = {
      natural_product_biosynthesis_not_applicable_to_external_background:
        'This candidate is classified as potential external background; natural-product biosynthetic inference is not applicable',
      generic_class_does_not_resolve_core_biosynthesis:
        'The current reliable class describes only the modification type and cannot establish the biosynthesis of the core scaffold',
      reliable_class_context_without_candidate_specific_rule:
        'Biosynthetic context is based on the reliable chemical class and is not candidate-specific source evidence',
      classification_not_reliable_for_biosource_context:
        'The current class assignment is not reliable enough to infer a biosynthetic pathway',
      biosource_knowledge_not_covered_for_current_class:
        'The curated knowledge base does not yet cover the biosource context of this candidate',
    };
    if (labels[code]) return labels[code];
    const translated = window.MSSLanguage?.translateText?.(
      String(analysis?.biosourceReasonZh || '').trim()
    );
    return translated || '';
  }

  function classificationState(match) {
    return classificationPresentation(match).state;
  }

  function biosyntheticFamily(value) {
    const text = String(value || '');
    if (text.includes('莽草酸')) return '莽草酸途径';
    if (text.includes('醋酸-丙二酸') || text.includes('乙酸-丙二酸')) return '醋酸-丙二酸途径';
    if (/MVA|MEP/i.test(text)) return 'MVA / MEP 途径';
    if (text.includes('氨基酸')) return '氨基酸途径';
    if (text.includes('糖酵解')) return '糖酵解相关途径';
    return text.replace(/^.*?生源\s*[:：]\s*/, '').slice(0, 42) || '其他生源背景';
  }

  function groupAnalysis(group) {
    const matches = (group?.rows || []).map(row => row.match || {});
    const representative = group?.rep?.match || matches[0] || {};
    const auxiliary = group?.compound?.auxiliary_interpretation
      || representative.auxiliary_interpretation || {};
    const expert = auxiliary.expert || {};
    const adduct = auxiliary.adduct_assignment || {};
    const biosource = auxiliary.biosource || {};
    const fragmentation = auxiliary.fragmentation || {};
    const flags = matches.map(match => match.display_evidence_flags || {});
    const biosyntheticPathways = [...new Set([
      ...(Array.isArray(biosource.pathways) ? biosource.pathways : []),
      ...flags.flatMap(flag => Array.isArray(flag.biosynthetic_pathways)
        ? flag.biosynthetic_pathways : []),
    ].map(value => String(value || '').trim()
      .replace(/^生源\s*[:：]\s*/, '')).filter(Boolean))];
    const maximum = key => Math.max(0, ...flags.map(flag => Number(flag[key]) || 0));
    const expertSources = [...new Set(matches.map(match => String(match.expert_score_source || '').trim())
      .filter(source => source && source !== 'none'))];
    const expertScores = matches.map(match => Number(match.expert_score)).filter(Number.isFinite);
    const pathwayStatuses = [...new Set(flags.map(flag => String(flag.pathway_status || '')).filter(Boolean))];
    const fragmentationStepDetails = Array.isArray(fragmentation.steps)
      ? fragmentation.steps : [];
    const fragmentationSteps = Number(fragmentation.step_count || 0)
      || maximum('fragmentation_step_count');
    const productStructures = Number(fragmentation.product_structure_count || 0)
      || maximum('product_structure_count');
    const hasLibraryPathwayReference = (
      fragmentation.candidate_specific_library_pathway?.status === 'available'
      || flags.some(flag => !!flag.has_library_pathway_reference)
    );
    const expertReasoning = Array.isArray(expert.reasoning) ? expert.reasoning : [];
    return {
      biosyntheticPathways,
      biosourceReason: String(biosource.reason || ''),
      biosourceReasonZh: String(biosource.reason_zh || ''),
      biosourceStatus: String(biosource.status || 'unavailable'),
      biosourceEvidenceScope: String(biosource.evidence_scope || 'unavailable'),
      hasReasoning: expertReasoning.length > 0 || flags.some(flag => !!flag.has_reasoning),
      expertReasoning,
      hasFragmentation: fragmentationSteps > 0 || productStructures > 0
        || flags.some(flag => !!flag.has_fragmentation),
      hasFragmentationEvidence: fragmentationSteps > 0 || productStructures > 0
        || Number(fragmentation.neutral_loss_count || 0) > 0
        || Number(fragmentation.diagnostic_ion_count || 0) > 0
        || flags.some(flag => !!flag.has_fragmentation_evidence),
      hasLibraryPathwayReference,
      hasNeutralLoss: Number(fragmentation.neutral_loss_count || 0) > 0
        || flags.some(flag => !!flag.has_neutral_loss),
      fragmentationSteps,
      fragmentationStepDetails,
      fragmentationReason: String(fragmentation.reason || ''),
      neutralLosses: Number(fragmentation.neutral_loss_count || 0)
        || maximum('neutral_loss_count'),
      diagnosticIons: Number(fragmentation.diagnostic_ion_count || 0)
        || maximum('diagnostic_ion_count'),
      productStructures,
      pathwayStatuses,
      pathwayMismatch: pathwayStatuses.includes('mismatch'),
      expertSources,
      hasExpert: ['evaluated', 'partially_evaluated'].includes(expert.status),
      expertPositive: expert.positive_support === true,
      expertCoverage: expert.coverage || {},
      sampleIonForm: String(adduct.sample_ion_form || ''),
      sampleIonFormSource: String(adduct.sample_ion_form_source || ''),
      libraryAdduct: String(adduct.library_adduct || representative.adduct || ''),
      candidateIonFormHypothesis: String(
        adduct.candidate_ion_form_hypothesis
        || representative.candidate_adduct_hypothesis || ''
      ),
    };
  }

  function analysisCoverage(bundle) {
    const classes = new Map();
    const pathways = new Map();
    const coverage = {
      total: bundle.groups.length,
      classReliable: 0, classTentative: 0, classUnresolved: 0,
      biosynthetic: 0,
      expertPositive: 0, expertNoEvidence: 0, expertUnavailable: 0,
      fragmentation: 0, neutralLoss: 0, libraryPathwayReference: 0,
      pathwayMismatch: 0,
    };
    bundle.groups.forEach(group => {
      const match = group.rep?.match || {};
      const className = rawClassification(match);
      const state = classificationState(match);
      if (state === 'reliable') coverage.classReliable += 1;
      else if (state === 'tentative') coverage.classTentative += 1;
      else coverage.classUnresolved += 1;
      if (className) classes.set(className, (classes.get(className) || 0) + 1);
      const analysis = groupAnalysis(group);
      if (analysis.biosyntheticPathways.length) coverage.biosynthetic += 1;
      [...new Set(analysis.biosyntheticPathways.map(biosyntheticFamily))].forEach(family => {
        pathways.set(family, (pathways.get(family) || 0) + 1);
      });
      if (analysis.expertPositive) coverage.expertPositive += 1;
      else if (analysis.hasExpert) coverage.expertNoEvidence += 1;
      else coverage.expertUnavailable += 1;
      if (analysis.hasFragmentation) coverage.fragmentation += 1;
      if (analysis.hasNeutralLoss) coverage.neutralLoss += 1;
      if (analysis.hasLibraryPathwayReference) coverage.libraryPathwayReference += 1;
      if (analysis.pathwayMismatch) coverage.pathwayMismatch += 1;
    });
    const top = (mapping, limit = 3) => [...mapping.entries()]
      .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0])).slice(0, limit);
    return { ...coverage, topClasses: top(classes), topPathways: top(pathways) };
  }

  function analysisChip(filter, label, count) {
    return `<button onclick="SampleOverview.setAnalysisFilter('${esc(filter)}')"><span>${esc(label)}</span><b>${integer(count)}</b></button>`;
  }

  function encodeFilterValue(value) {
    return encodeURIComponent(String(value || '')).replace(/'/g, '%27');
  }

  function analysisCoverageCard(kind, title, headline, description, chips, footnote) {
    return `<article class="sample-analysis-card analysis-${kind}">
      <div><span>${esc(kind.toUpperCase())}</span><strong>${esc(title)}</strong></div>
      <b>${esc(headline)}</b>
      <p>${esc(description)}</p>
      <div class="analysis-card-breakdown">${chips.join('')}</div>
      <em>${esc(footnote)}</em>
    </article>`;
  }

  function sampleAnalysisFunnel(task, coverage, conclusionMetadata = {}) {
    const funnel = task.identification_funnel || {};
    const workflow = task.untargeted_workflow || {};
    const moduleAvailability = task.module_availability || {};
    const featureModule = moduleAvailability.ms1_feature_detection || {};
    const ms1 = Number(task.ms1_scans || workflow.ms1_spectra || 0);
    const ms2 = Number(task.total_spectra || funnel.ms2_scans || 0);
    const auxiliary = Number(task.spectrum_representation_processing?.auxiliary_spectra_excluded_from_peak_picking || 0);
    const features = Number(workflow.feature_count || workflow.raw_features || 0);
    const linkedFeatures = Number(featureModule.lcms_features_with_ms2 || workflow.ms2_attachment?.feature_count || 0);
    const linkedMs2 = Number(featureModule.linked_ms2_spectra || workflow.ms2_attachment?.linked_ms2_spectra || 0);
    const matched = Number(task.matched_spectra || funnel.matched_ms2_scans || 0);
    const conclusionSummary = Object.keys(conclusionMetadata || {}).length
      ? conclusionMetadata : (task.compound_conclusion_summary || {});
    const candidates = Number(conclusionSummary.unique_report_rows || coverage.total || 0);
    const mainResults = Number(conclusionSummary.main_result_count || 0);
    const networks = Number(moduleAvailability.iimn?.network_count || task.ion_identity_networks || 0);
    const unmatched = Number(funnel.unmatched_ms2_scans || Math.max(0, ms2 - matched));
    return `<section class="sample-analysis-funnel" aria-label="鉴定结果概览">
      <div class="sample-analysis-funnel-head"><div><span class="panel-eyebrow">RESULT OVERVIEW</span><h3>鉴定结果概览</h3></div><p>化合物数量已经去重；MS2 扫描数、色谱特征数和化合物数是不同统计单位。</p></div>
      <div class="sample-funnel-grid">
        <div><span>输入数据</span><b>${integer(ms1)} MS1 · ${integer(ms2)} MS2</b><small>${auxiliary ? `${integer(auxiliary)} 张 DAD/UV 辅助谱已隔离` : '仅统计质谱扫描'}</small></div>
        <div><span>色谱特征</span><b>${integer(features)} 个 OpenMS 特征</b><small>${integer(linkedFeatures)} 个特征关联 ${integer(linkedMs2)} 张 MS2</small></div>
        <div><span>去重化合物结果</span><b>${integer(mainResults)} 个主要结果</b><small>${integer(candidates)} 个化合物／异构体组线索，含背景与低证据项</small></div>
        <div><span>化合物解释</span><b>${integer(coverage.classReliable)} 可靠分类 · ${integer(coverage.biosynthetic)} 类别级生源背景</b><small>${integer(coverage.fragmentation)} 个结果重建了样品碎裂；${integer(coverage.libraryPathwayReference)} 个有谱库参考</small></div>
        <div><span>离子网络</span><b>${integer(networks)} 个 IIMN</b><small>网络关系不冒充独立化合物鉴定</small></div>
        <div><span>未匹配保留</span><b>${integer(unmatched)} 张 MS2</b><small>${integer(funnel.unmatched_with_nonzero_query_peaks || 0)} 张仍有可用碎片</small></div>
      </div>
    </section>`;
  }

  function coreAnalysisSummary(bundle) {
    const coverage = analysisCoverage(bundle);
    const task = bundle.task || {};
    const funnel = task.identification_funnel || {};
    const unmatchedMs2 = Number(funnel.unmatched_ms2_scans || 0);
    const unmatchedWithPeaks = Number(funnel.unmatched_with_nonzero_query_peaks || 0);
    const classChips = [
      analysisChip('class-state:reliable', '可靠来源', coverage.classReliable),
      analysisChip('class-state:tentative', '待核查线索', coverage.classTentative),
      analysisChip('class-state:unresolved', '未可靠判类', coverage.classUnresolved),
      ...coverage.topClasses.map(([label, count]) => analysisChip(`class:${encodeFilterValue(label)}`, label, count)),
    ];
    const pathwayChips = [
      analysisChip('biosource:any', '有类别级生源背景', coverage.biosynthetic),
      analysisChip('biosource:none', '知识未覆盖／类别不足', coverage.total - coverage.biosynthetic),
      ...coverage.topPathways.map(([label, count]) => analysisChip(`biosource:${encodeFilterValue(label)}`, label, count)),
    ];
    const expertChips = [
      analysisChip('expert:positive', '至少一张谱有正向支持', coverage.expertPositive),
      analysisChip('expert:no_evidence', '已运行无正向支持', coverage.expertNoEvidence),
      analysisChip('expert:unavailable', '不可用／无来源', coverage.expertUnavailable),
    ];
    const fragmentChips = [
      analysisChip('fragment:path', '已重建样品碎裂', coverage.fragmentation),
      analysisChip('fragment:reference', '有谱库碎裂参考', coverage.libraryPathwayReference),
      analysisChip('fragment:neutral', '有中性丢失结果', coverage.neutralLoss),
      analysisChip('fragment:mismatch', '馆藏路径极性不符', coverage.pathwayMismatch),
    ];
    return `${sampleAnalysisFunnel(task, coverage, bundle.compounds?.metadata || {})}<section class="sample-core-analysis" aria-label="化合物解释信息">
      <div class="sample-core-analysis-head">
        <div><span class="panel-eyebrow">COMPOUND INTERPRETATION</span><h3>化合物解释信息</h3></div>
        <p>分类、生源、结构规则证据和推测碎裂途径均会展示；无法推理时明确说明原因。这些解释不改变化合物身份排序。</p>
      </div>
      <div class="sample-analysis-card-grid">
        ${analysisCoverageCard('class', '分类分析', `${integer(coverage.classReliable)} 可靠 · ${integer(coverage.classTentative)} 待核查`, '分类可靠性与谱库身份匹配等级独立评估；身份为 High 时，类别仍可能因当前体系未覆盖而缺失。', classChips, '点击类别或可靠性状态筛选候选')}
        ${analysisCoverageCard('biosource', '生源背景', `${integer(coverage.biosynthetic)} / ${integer(coverage.total)} 候选`, '只在化学类别可靠时展示类别级生物合成背景；不是候选特异证明，也不影响谱库身份等级。', pathwayChips, '未覆盖会显示知识边界，不再写成“候选没有触发规则”')}
        ${analysisCoverageCard('expert', '规则证据评估', `${integer(coverage.expertPositive)} 个结果有正向规则支持`, '逐维显示规则是否适用、是否观察到支持以及证据来源；它不直接决定身份。', expertChips, '每个维度明确标注有支持、未观察、不适用或未建立')}
        ${analysisCoverageCard('fragment', '样品碎裂重建', `${integer(coverage.fragmentation)} 个结果有样品步骤/产物`, '样品 MS/MS 重建、候选谱库参考和纯质量差证据分开统计；没有可支持路径时明确报告不可用。', fragmentChips, '候选谱库路径不等于本样品已重建，推测路径也不等同已确认气相机理')}
      </div>
      <div id="activeAnalysisFilter" class="active-analysis-filter" hidden></div>
      <div class="sample-unmatched-analysis-note"><div><b>未匹配谱图没有被删除：</b>${integer(unmatchedMs2)} 张当前没有谱库候选${unmatchedWithPeaks ? `，其中 ${integer(unmatchedWithPeaks)} 张仍有可用碎片` : ''}。未匹配只表示当前谱库与阈值下没有足够相似的参比谱。</div><button class="btn btn-secondary btn-sm" onclick="SampleOverview.openUnmatched()">进入未匹配谱图工作台 →</button></div>
    </section>`;
  }

  function contaminantPresentation(match) {
    if (!match?.potential_contaminant) return null;
    const role = String(match.compound_role || 'external_contaminant');
    const category = String(match.contaminant_category || '').trim();
    const label = role === 'pesticide_screening_target'
      ? '农药筛查对象'
      : '潜在外源背景';
    const details = [category, match.pesticide_use, match.screening_domain]
      .map(value => String(value || '').trim()).filter(Boolean);
    return {
      label: details.length ? `${label} · ${details[0]}` : label,
      note: String(match.contaminant_note || match.compound_role_note || '').trim(),
    };
  }


  function init() {
    document.getElementById('resultTaskSelect')?.addEventListener('change', loadForTask);
    document.addEventListener('mssystem:languagechange', () => {
      if (!currentTaskId) return;
      SampleResultData.clear(currentTaskId);
      if (document.getElementById('page-results')?.classList.contains('active')) {
        loadForTask();
      }
    });
  }

  function taskConfigurationLabel(task) {
    const english = isEnglishInterface();
    const blankCount = Number(task?.blank_summary?.assigned_blank_replicate_count
      ?? task?.blank_summary?.blank_replicate_count ?? 0);
    const selection = task?.expected_compounds_summary?.target_dataset_selection || {};
    const source = String(selection.selection_source
      || task?.params?.expected_compound_dataset_selection_source || 'none');
    const labels = [];
    const blankRoles = task?.blank_summary?.blank_roles || {};
    const unspecified = Number(blankRoles.unspecified || 0);
    labels.push(blankCount > 0
      ? (task?.blank_summary?.blank_operation === 'comparison_only'
        ? (english
          ? `Blank-background comparison: ${blankCount} (0 excluded${unspecified > 0 ? '; method pairing not verified' : ''})`
          : `空白背景比较 ${blankCount}（0 条扣除${unspecified > 0 ? '，方法配对未核验' : ''}）`)
        : (unspecified > 0
          ? (english
            ? `Untyped same-polarity blanks: ${blankCount} (method pairing not verified; 0 excluded)`
            : `同极性未定角色空白 ${blankCount}（方法配对未核验，0 条扣除）`)
          : (english
            ? `Confirmed background blanks: ${blankCount} (${Number(task?.blank_summary?.excluded_spectra || 0)} excluded)`
            : `已确认背景空白 ${blankCount}（扣除 ${Number(task?.blank_summary?.excluded_spectra || 0)} 条）`)))
      : (english ? 'No blank evidence' : '无空白证据'));
    if (source === 'filename_autoselect') labels.push(
      english ? 'Legacy filename-based target selection' : '旧版文件名自动靶向'
    );
    else if (source === 'user_explicit') labels.push(
      english ? 'Explicit targeted library' : '显式靶向库'
    );
    else labels.push(english ? 'No targeted library' : '未绑定靶向库');
    return labels.join(' · ');
  }

  function configurationBanner(task) {
    const english = isEnglishInterface();
    const label = taskConfigurationLabel(task);
    const selection = task?.expected_compounds_summary?.target_dataset_selection || {};
    const selectionSource = String(selection.selection_source
      || task?.params?.expected_compound_dataset_selection_source || 'none');
    const legacyAutoselect = selectionSource === 'filename_autoselect';
    const blank = task?.blank_summary || {};
    const hardExclusion = blank.blank_operation === 'confirmed_presearch_background_exclusion'
      || blank.filter_policy === 'hard_legacy';
    const blankOperation = english
      ? (hardExclusion
        ? `The selected blanks were used for feature-level background exclusion; ${integer(blank.excluded_spectra || 0)} MS2 scans were excluded. Exclusion records remain in the audit layer, and fragment intensities are not subtracted peak by peak across injections.`
        : 'Selected blanks are background-comparison inputs, not sample analyses. “Blank-related” is a comparison flag and does not mean that subtraction was performed.')
      : (hardExclusion
        ? `所选空白已用于特征级背景排除，共排除 ${integer(blank.excluded_spectra || 0)} 张 MS2；排除记录仍保存在审计层，不进行跨针碎片强度逐峰相减。`
        : '所选空白仅作为背景比较输入，不作为样品任务分析；“空白相关”是比较标记，不代表已经执行扣除。');
    const comparisonBoundary = english
      ? ' Runs with different blank or targeted-library configurations should not be compared directly; reprojection does not alter the original search results.'
      : '不同空白或靶向库配置的任务不可直接比较；重新投影不会篡改当时的原始搜索结果。';
    return `<div class="sample-result-scope" style="margin-bottom:.7rem;border-left-color:${legacyAutoselect ? '#d97706' : '#64748b'}"><b>${english ? 'Saved-run configuration' : '本次保存结果的运行配置'}</b><span>${esc(label)}${english ? '. ' : '。'}${esc(blankOperation)}${esc(comparisonBoundary)}</span></div>`;
  }

  async function refreshSelect(preferredTaskId = null) {
    const select = document.getElementById('resultTaskSelect');
    if (!select) return;
    try {
      const data = await API.get('/tasks?status=completed');
      const previous = preferredTaskId || select.value;
      const tasks = data.tasks || [];
      const latestSelection = select.value;
      select.innerHTML = '<option value="">-- 选择已完成的任务 --</option>' +
        tasks.map(task => {
          const stamp = task.upload_time
            ? String(task.upload_time).replace('T', ' ').replace(/\.\d+.*$/, '').slice(0, 16)
            : '时间未记录';
          return `<option value="${esc(task.id)}">${esc(task.filename)} · ${esc(stamp)} · ${esc(taskConfigurationLabel(task))} · ${integer(task.matched_spectra)}/${integer(task.total_spectra)} MS2 有候选</option>`;
        }).join('');
      const target = preferredTaskId || latestSelection || previous || currentTaskId;
      if (target && tasks.some(task => task.id === target)) select.value = target;
    } catch (error) {
      console.error('Failed to refresh sample result selector:', error);
    }
  }

  function summaryNumbers(bundle) {
    const task = bundle.task || {};
    const funnel = task.identification_funnel || {};
    const metadata = bundle.compounds?.metadata || {};
    const eventAggregation = metadata.precursor_event_aggregation || {};
    const derived = key => bundle.groups.filter(group => group.status.key === key).length;
    const namedGroups = bundle.groups.filter(group => ![
      'unannotated', 'unresolved_event',
    ].includes(group.status.key));
    const unresolved = Number(metadata.unresolved_precursor_event_count
      ?? derived('unresolved_event'));
    const unresolvedIncomplete = Number(
      eventAggregation.incomplete_formula_unresolved_event_count
      ?? bundle.groups.filter(group => (
        group.status.key === 'unresolved_event'
        && unresolvedEventKind(group) === 'incomplete_formula'
      )).length
    );
    const unresolvedFormulaConflict = Number(
      eventAggregation.cross_formula_unresolved_event_count
      ?? Math.max(0, unresolved - unresolvedIncomplete)
    );
    return {
      compounds: namedGroups.length,
      main: bundle.groups.filter(group => group.status.countedAsSampleCompound).length,
      plantMain: bundle.groups.filter(group => group.status.countedAsPlantCandidate
          && !explicitExternalReview(group)).length,
      externalReview: bundle.groups.filter(explicitExternalReview).length,
      confirmed: derived('confirmed'),
      high: derived('high'),
      isomer: derived('isomer'),
      isomerHigh: bundle.groups.filter(group => group.status.hasIsomerHigh).length,
      highScopeGroups: derived('high')
        + bundle.groups.filter(group => group.status.hasIsomerHigh).length,
      spectralHighScopeGroups: Number(metadata.spectral_high_confidence_scope_group_count
        ?? metadata.high_confidence_scope_group_count ?? 0),
      backgroundHigh: Number(metadata.background_high_confidence_scope_group_count ?? 0),
      unverifiedBlank: Number(metadata.unverified_blank_warning_count ?? 0),
      unverifiedBlankHigh: Number(metadata.unverified_blank_warning_high_count ?? 0),
      medium: derived('medium'),
      candidate: derived('candidate'),
      unresolved,
      unresolvedFormulaConflict,
      unresolvedIncomplete,
      insufficient: derived('low'),
      unannotated: Number(
        metadata.unannotated_feature_count ?? derived('unannotated')
      ),
      background: derived('background'),
      legacyHigh: Number(metadata.legacy_high_supported_identity_count ?? 0),
      unmatched: Number(funnel.unmatched_ms2_scans || 0),
      contaminants: Number(task.contaminant_unique_compounds
        ?? bundle.groups.filter(explicitExternalReview).length),
      matchedSpectra: Number(task.matched_spectra ?? bundle.matched.length),
      totalSpectra: Number(task.total_spectra ?? funnel.ms2_scans ?? 0),
    };
  }

  function resultHeadline(bundle, counts) {
    const task = bundle.task || {};
    const english = isEnglishInterface();
    const mode = String(
      task.ion_mode || task.params?.ion_mode || (english ? 'Unknown' : '未知')
    ).toUpperCase();
    const isomerNonHigh = counts.isomer - counts.isomerHigh;
    const plural = (value, singular, pluralForm = `${singular}s`) => (
      Number(value) === 1 ? singular : pluralForm
    );
    const summary = english
      ? `<p><b>${integer(counts.matchedSpectra)}</b> MS2 scans produced library candidates and were consolidated into <b>${integer(counts.compounds)}</b> deduplicated candidate ${plural(counts.compounds, 'identity', 'identities')}. The sample-supported primary results comprise <b>${integer(counts.main)}</b> ${plural(counts.main, 'identity', 'identities')}: <b>${integer(counts.highScopeGroups)}</b> pass the current entity High conditions (<b>${integer(counts.high)}</b> single 2-D ${plural(counts.high, 'connectivity', 'connectivities')} and <b>${integer(counts.isomerHigh)}</b> unresolved isomer ${plural(counts.isomerHigh, 'group')}), with <b>${integer(isomerNonHigh)}</b> additional moderate-evidence isomer ${plural(isomerNonHigh, 'group')}, <b>${integer(counts.medium)}</b> explicitly evaluated entity-Medium ${plural(counts.medium, 'candidate')}, and <b>${integer(counts.candidate)}</b> saved but entity-ungraded ${plural(counts.candidate, 'candidate')}; <b>${integer(counts.unresolved)}</b> precursor ${plural(counts.unresolved, 'event')} are unresolved and excluded from the independent compound count (<b>${integer(counts.unresolvedFormulaConflict)}</b> with conflicting known formulae; <b>${integer(counts.unresolvedIncomplete)}</b> with incomplete candidate-formula assignments); <b>${integer(counts.insufficient)}</b> candidate ${plural(counts.insufficient, 'identity', 'identities')} have insufficient evidence. Formal confirmations are reported separately: <b>${integer(counts.confirmed)}</b>.</p>`
      : `<p><b>${integer(counts.matchedSpectra)}</b> 张 MS2 获得谱库候选，按候选二维连接与分子式归并为 <b>${integer(counts.compounds)}</b> 个去重候选身份。样品支持的主要结果为 <b>${integer(counts.main)}</b> 个身份：<b>${integer(counts.highScopeGroups)}</b> 个通过当前实体 High 条件（<b>${integer(counts.high)}</b> 个单一二维连接、<b>${integer(counts.isomerHigh)}</b> 个未解异构体组），另有 <b>${integer(isomerNonHigh)}</b> 个中等证据异构体组、<b>${integer(counts.medium)}</b> 个明确保存并完成实体评估的 Medium 候选，以及 <b>${integer(counts.candidate)}</b> 个已保存但实体未分级的候选；<b>${integer(counts.unresolved)}</b> 个母离子事件未解析且不计作独立化合物（<b>${integer(counts.unresolvedFormulaConflict)}</b> 个已知分子式冲突，<b>${integer(counts.unresolvedIncomplete)}</b> 个候选分子式记录不完整）；<b>${integer(counts.insufficient)}</b> 个候选身份证据不足。正式确证结果单列为 <b>${integer(counts.confirmed)}</b>。</p>`;
    const blankBoundary = !counts.unverifiedBlank ? '' : (english
      ? `<div class="sample-result-scope"><b>Blank-evidence boundary</b><span><b>${integer(counts.unverifiedBlank)}</b> identification ${plural(counts.unverifiedBlank, 'group')} also contain signal in an untyped same-polarity blank; <b>${integer(counts.unverifiedBlankHigh)}</b> reach a High reportable scope. Because batch and chromatographic-method pairing have not been verified, these signals are retained as warnings and no sample candidate is removed on this basis.</span></div>`
      : `<div class="sample-result-scope"><b>空白证据边界</b><span><b>${integer(counts.unverifiedBlank)}</b> 个鉴定组在同极性未定角色空白中亦有信号，其中 <b>${integer(counts.unverifiedBlankHigh)}</b> 个达到 High 可报告范围；因批次／色谱方法配对未经核验，系统只保留警示，不据此剔除样本候选。</span></div>`);
    return `
      <section class="sample-result-hero">
        <div class="sample-result-hero-copy">
          <div class="sample-result-summary-heading">
            <span class="sample-result-status-dot" aria-hidden="true"></span>
            <div><span class="sample-result-kicker">ANALYSIS COMPLETE</span><h3>样本结果摘要</h3></div>
            <span class="sample-result-mode">${esc(mode)} 离子模式</span>
          </div>
          ${SampleResultData.confidenceAxesHtml(bundle)}
          ${summary}
          <div class="sample-result-scope"><b>${english ? 'Plant candidate table and role review' : '正文候选与外源复核'}</b><span>${english
            ? `The default table contains ${integer(counts.plantMain)} sample-supported candidate groups after excluding reviewed external candidates and mixed unresolved groups. ${integer(counts.externalReview)} groups are retained in role review with their original spectra, winners and evidence grades. Unclassified origin does not establish plant origin; the High counts above describe spectral evidence before role filtering.`
            : `默认正文候选表保留 ${integer(counts.plantMain)} 个样品支持候选组，已排除明确外源候选及含外源候选的未解混合组；${integer(counts.externalReview)} 个组保留在外源复核中，原谱、首选和证据等级不变。来源未核验不等于植物内源；上方 High 计数表示来源筛选前的谱库证据。`}</span></div>
          ${blankBoundary}
          <div class="sample-result-analysis-tags" aria-label="结果处理内容"><span>化合物去重</span><span>二维连接／异构体分层</span><span>空白证据分层</span><span>谱库证据分级</span><span>解释信息保留</span></div>
          <div class="sample-result-scope"><b>怎样理解 High 级鉴定组</b><span>系统按能够诚实报告的身份范围分开计数：单一二维连接结构与未能区分的同分子式异构体组绝不混写；立体化学／完整结构定义层另行标注。High 级谱库支持仍不等同于同方法标准品 Level 1 确证。</span></div>
        </div>
        <div class="sample-result-hero-actions">
          <button class="btn btn-primary btn-sm" onclick="SampleOverview.openReview()">鉴定详情</button>
          <button class="btn btn-secondary btn-sm" onclick="SampleOverview.openReport()">综合报告</button>
          <button class="btn btn-secondary btn-sm" onclick="SampleOverview.openFollowup()">确认实验</button>
        </div>
      </section>`;
  }

  function targetedResultsSection(payload) {
    if (typeof Results === 'undefined'
        || typeof Results.expectedCompoundsCard !== 'function') return '';
    return `<section id="sampleTargetedResults" class="sample-targeted-results" aria-label="靶向数据库结果">${Results.expectedCompoundsCard(payload)}</section>`;
  }

  function evitwinSection(payload) {
    const compounds = Array.isArray(payload?.compounds) ? payload.compounds : [];
    const actualCeShadow = payload?.actual_ce_shadow || {};
    const actualCeFeatures = Array.isArray(actualCeShadow.features)
      ? actualCeShadow.features : [];
    if (!compounds.length && !actualCeFeatures.length) return '';
    const summary = payload.summary || {};
    const rows = compounds.slice(0, 30).map(feature => {
      const winner = feature.evitwin_winner || {};
      const publicWinner = feature.public_library_winner || {};
      const decision = feature.decision || {};
      const rt = Number(feature.representative_retention_time_min);
      const similarity = Number(winner.median_similarity);
      const margin = Number(winner.winner_margin);
      const changed = feature.display_rank_changed_from_public_library === true;
      const competing = decision.reportable_scope === 'isomer_group'
        ? (Array.isArray(feature.candidate_scores) ? feature.candidate_scores : [])
          .filter(item => Number(item?.median_similarity) >= 0.5)
          .slice(0, 4)
          .map(item => `${esc(item.compound_name || item.compound_code || '候选')} ${Number(item.median_similarity).toFixed(3)}`)
        : [];
      const statusLabels = {
        qualified_authentic_rank: '标品谱优先候选',
        qualified_retention_order_rank: '标品谱＋保留顺序联合候选',
        authentic_isomer_ranking: '标品谱异构体组先导',
      };
      const scope = decision.reportable_scope === 'isomer_group'
        ? '异构体组；不释放单体确证'
        : '二维连接候选；尚未自动升级 Level 1';
      return `<tr>
        <td><b>${esc(winner.compound_name || winner.compound_code || '未命名标品')}</b><small style="display:block;color:var(--text-secondary);">${esc(winner.compound_code || '')}${winner.formula ? ` · ${esc(winner.formula)}` : ''}</small>${competing.length ? `<small style="display:block;margin-top:.25rem;color:#9a3412;">组内竞争：${competing.join('；')}</small>` : ''}</td>
        <td><b>${Number.isFinite(similarity) ? similarity.toFixed(3) : '—'}</b><small style="display:block;color:var(--text-secondary);">${integer(winner.collision_energy_count)} 档 CE · ${integer(winner.collision_energy_votes)} 票${Number.isFinite(margin) ? ` · margin ${margin.toFixed(3)}` : ''}</small></td>
        <td>${Number.isFinite(rt) ? `${rt.toFixed(3)} min` : 'RT 未记录'}</td>
        <td><span class="sample-evidence-badge evidence-${decision.reportable_scope === 'isomer_group' ? 'isomer' : 'high'}">${esc(statusLabels[decision.status] || '标品谱候选')}</span><small style="display:block;margin-top:.3rem;color:var(--text-secondary);">${esc(scope)}</small></td>
        <td>${changed ? `<b style="color:#0f766e;">已调整首位</b><small style="display:block;color:var(--text-secondary);">公开谱库原先：${esc(publicWinner.compound_name || '无候选')}</small>` : `<span>与公开谱库一致</span>`}</td>
      </tr>`;
    }).join('');
    const actualCeRows = actualCeFeatures.slice(0, 30).map(feature => {
      const winner = feature.shadow_winner || {};
      const publicWinner = feature.public_library_winner || {};
      const decision = feature.decision || {};
      const energies = Array.isArray(feature.collision_energy_conditions)
        ? feature.collision_energy_conditions.map(value => `${Number(value).toFixed(2)} eV`).join(' / ')
        : '—';
      const score = Number(winner.median_similarity);
      const margin = Number(winner.winner_margin);
      const differs = feature.shadow_lead_differs_from_public_library === true;
      return `<tr>
        <td><b>${esc(winner.compound_name || winner.compound_code || '未命名标品候选')}</b><small style="display:block;color:var(--text-secondary);">${esc(winner.compound_code || '')}${winner.formula ? ` · ${esc(winner.formula)}` : ''}</small></td>
        <td><b>${Number.isFinite(score) ? score.toFixed(3) : '—'}</b><small style="display:block;color:var(--text-secondary);">实际 CE ${esc(energies)}${Number.isFinite(margin) ? ` · margin ${margin.toFixed(3)}` : ''}</small></td>
        <td>${integer(winner.collision_energy_count)} 个查询 CE · ${integer(winner.collision_energy_votes)} 票</td>
        <td><span class="sample-evidence-badge evidence-isomer">候选组 shadow</span><small style="display:block;margin-top:.3rem;color:var(--text-secondary);">${esc((decision.reason_codes || []).join('；'))}</small></td>
        <td>${differs ? `<b style="color:#9a3412;">影子先导不同</b><small style="display:block;color:var(--text-secondary);">正式首位保持：${esc(publicWinner.compound_name || '无候选')}</small>` : '<span>影子先导与正式首位一致</span>'}</td>
      </tr>`;
    }).join('');
    return `<section class="card" style="border-left:4px solid #0d9488;">
      <div class="card-header"><div><div class="card-title">EviTwin-MS 标品谱联合排序</div><div class="card-subtitle">优先回答“分析到了哪些化合物”：真实单标 MS/MS 先于公开谱库候选显示；竞争异构体仍按组级边界报告</div></div></div>
      <section class="sample-result-metrics" aria-label="EviTwin-MS 排序摘要" style="margin-top:.7rem;">
        ${statusCard(summary.priority_compound_count, '优先标品谱候选', '按二维连接去重', 'high')}
        ${statusCard(summary.display_rank_change_count, '首位顺序已调整', '真实标品谱证据支持重新排序', 'confirmed')}
        ${statusCard(summary.retention_order_assignment_count, '保留顺序联合分配', '同质量异构体在同一混标内一对一排序', 'isomer')}
      </section>
      <div class="sample-result-scope" style="margin-top:.65rem;"><b>证据边界</b><span>这里改变的是标品参考通道的候选显示顺序，不删除公开谱库结果。跨方法绝对 RT 未锁定前不自动增加 Level 1 确证；同二维连接立体异构体不点名单体。</span></div>
      ${compounds.length ? `<div class="table-wrap" style="margin-top:.65rem;"><table><thead><tr><th>优先候选</th><th>标品谱证据</th><th>样品 RT</th><th>可报告范围</th><th>相对公开谱库</th></tr></thead><tbody>${rows}</tbody></table></div>` : ''}
      ${actualCeFeatures.length ? `<div class="sample-result-scope" style="margin-top:.8rem;border-color:#fed7aa;background:#fff7ed;"><b>已有实际 CE 谱图比较（shadow）</b><span>只使用文件中记录的查询 CE；参考 CE 之间只内插相似度数值，不内插谱图、不向参考范围外外推。单个查询 CE 仍只算一票，因此本表不改正式 Top1、不增加确证化合物数。</span></div><div class="table-wrap" style="margin-top:.55rem;"><table><thead><tr><th>真品谱候选组先导</th><th>实际 CE 相似度</th><th>查询证据数</th><th>当前状态</th><th>与正式首位比较</th></tr></thead><tbody>${actualCeRows}</tbody></table></div>` : ''}
    </section>`;
  }

  function replicateConsistencySection(payload) {
    if (payload?.status !== 'available') return '';
    const english = isEnglishInterface();
    const summary = payload.summary || {};
    const injections = Array.isArray(payload.injections) ? payload.injections : [];
    const injectionRows = injections.map(item => `<tr>
      <td><b data-scientific-value>${esc(item.filename || item.task_id)}</b><small style="display:block;color:var(--text-secondary);">${english ? 'Technical injection' : '重复针'} ${integer(item.replicate_number)}</small></td>
      <td>${integer(item.compound_group_count)}</td>
      <td>${integer(item.main_group_count)}</td>
      <td><b>${integer(item.high_group_count)}</b></td>
    </tr>`).join('');
    const interpretation = english
      ? String(payload.interpretation || 'Technical-repeat robustness is a separate cohort projection. It does not rescore spectra, change per-injection ranking, or constitute independent identity evidence.')
      : String(payload.interpretation_zh || '重复针稳健性是独立展示层：不重算谱图、不改动单针候选排序或置信等级，也不等同于独立的结构确证证据。');
    return `<section class="card" style="border-left:4px solid #4f46e5;">
      <div class="card-header"><div><div class="card-title">${english ? 'Technical-injection consistency' : '技术重复针一致性'}</div><div class="card-subtitle">${english ? 'Grouped by terminal filename replicate number and identical analysis parameters; spectra are not rescored and no per-injection Top1 or High classification is changed' : '按文件名末尾重复编号及完全相同的分析参数自动成组；不重算谱图、不改动任何单针 Top1 或 High'}</div></div></div>
      <section class="sample-result-metrics" aria-label="${english ? 'Technical-injection consistency metrics' : '重复针一致性指标'}" style="margin-top:.7rem;">
        ${statusCard(injections.map(item => integer(item.high_group_count)).join(' / '), english ? 'High groups per injection' : '各针 High', injections.map(item => item.filename).join(' / '), 'total')}
        ${statusCard(summary.robust_high_count, english ? 'Repeat-robust High' : '重复稳健 High', english ? `High in at least ${integer(payload.robustness_threshold)} / ${integer(payload.injection_count)} injections` : `至少 ${integer(payload.robustness_threshold)} / ${integer(payload.injection_count)} 针达到 High`, 'high')}
        ${statusCard(summary.all_injection_high_count, english ? 'High in every injection' : '全部重复针均 High', english ? 'High in all technical injections' : '所有技术重复均达到 High', 'confirmed')}
        ${statusCard(summary.single_injection_high_count, english ? 'High in one injection only' : '仅单针 High', english ? 'Indicates threshold or candidate-ranking instability; not repeat-robust' : '提示阈值或候选排序不稳定，不能称为重复稳健', 'isomer')}
      </section>
      <div class="sample-result-scope" style="margin-top:.65rem;"><b>${english ? 'Counting boundary' : '统计边界'}</b><span>${esc(interpretation)} ${english ? '“Repeat-robust High” is reported separately from the original per-injection High classifications.' : '“重复稳健 High”与每针原始 High 分开呈现。'}</span></div>
      <details style="margin-top:.65rem;"><summary style="cursor:pointer;font-weight:700;color:#3730a3;">${english ? 'View counts by injection' : '查看每针数量'}</summary><div class="table-wrap" style="margin-top:.5rem;"><table><thead><tr><th>${english ? 'File' : '文件'}</th><th>${english ? 'Identification groups' : '鉴定组'}</th><th>${english ? 'Primary results' : '主要结果'}</th><th>${english ? 'High groups' : 'High 组'}</th></tr></thead><tbody>${injectionRows}</tbody></table></div></details>
    </section>`;
  }

  function dualPolarityShadowPayload(response) {
    if (!response || typeof response !== 'object') return {};
    return response.sidecar || response.artifact || response.payload || response;
  }

  function finiteSummaryNumber(value) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : null;
  }

  function dualPolarityNotEvaluableFallback(summary) {
    return ['POS', 'NEG'].reduce((total, polarity) => {
      const extraction = summary?.anchor_extraction?.[polarity] || {};
      return total
        + Number(extraction.rows_without_reliable_anchor || 0)
        + Number(extraction.rejected_conflict_groups || 0)
        + Number(extraction.rejected_inconsistent_neutral_mass_groups || 0)
        + Number(extraction.rejected_missing_rt_groups || 0);
    }, 0);
  }

  function dualPolarityShadowResult(response) {
    const english = isEnglishInterface();
    if (!response) {
      return `<div class="v16-dual-result-empty">${english
        ? 'After an explicit run, this panel will show accepted links, ambiguous anchors, non-evaluable evidence and connectivity support.'
        : '显式运行后，这里会列出已接受连接、歧义锚点、不可评价证据和二维连接支持。'}</div>`;
    }
    const payload = dualPolarityShadowPayload(response);
    const summary = payload.summary || response.summary || {};
    const accepted = finiteSummaryNumber(summary.accepted_neutral_entity_count)
      ?? finiteSummaryNumber(summary.candidate_independent_accepted_edge_count)
      ?? finiteSummaryNumber(summary.neutral_entity_count) ?? 0;
    const compatible = finiteSummaryNumber(
      summary.candidate_independent_compatible_edge_count
    );
    const ambiguous = finiteSummaryNumber(summary.ambiguous_anchor_count)
      ?? (compatible == null ? 0 : Math.max(0, compatible - accepted));
    const notEvaluable = finiteSummaryNumber(summary.not_evaluable_anchor_count)
      ?? dualPolarityNotEvaluableFallback(summary);
    const support = finiteSummaryNumber(summary.dual_polarity_supported_connectivity_count)
      ?? finiteSummaryNumber(summary.dual_supported_connectivity_count) ?? 0;
    const invariants = payload.invariants || response.invariants || {};
    const invariantKeys = [
      'production_top1_changed',
      'production_top5_changed',
      'production_score_changed',
      'production_decision_changed',
      'production_confidence_changed',
      'production_compound_count_changed',
    ];
    const productionUnchanged = invariantKeys.every(key => invariants[key] === false);
    const pairId = String(response.pair_id || payload.metadata?.pair_id || '').trim();
    const state = String(summary.status || response.status || payload.status || '').trim();
    return `<div class="v16-dual-result-summary" aria-live="polite">
      <div class="v16-dual-result-heading"><b>${english ? 'Shadow evidence result' : '影子证据结果'}</b>${pairId ? `<span class="mono">${esc(pairId)}</span>` : ''}${state ? `<span>${esc(state)}</span>` : ''}</div>
      <div class="v16-dual-metrics">
        <div><strong>${integer(accepted)}</strong><span>${english ? 'Accepted neutral links' : '已接受中性实体连接'}</span></div>
        <div><strong>${integer(ambiguous)}</strong><span>${english ? 'Ambiguous anchors' : '歧义锚点'}</span></div>
        <div><strong>${integer(notEvaluable)}</strong><span>${english ? 'Not evaluable' : '不可评价'}</span></div>
        <div><strong>${integer(support)}</strong><span>${english ? 'Dual-polarity connectivity support' : '双极性二维连接支持'}</span></div>
      </div>
      <div class="v16-dual-invariant ${productionUnchanged ? 'is-verified' : 'is-incomplete'}">
        <b>${productionUnchanged
          ? (english ? 'Production invariants verified' : '生产不变量已核验')
          : (english ? 'Production invariant receipt is incomplete' : '生产不变量回执不完整')}</b>
        <span>${english
          ? 'This shadow result does not change V15 Top1/Top5, scores, decisions, confidence or compound counts. Cross-polarity MS2 similarity remains auxiliary only.'
          : '该影子结果不会改变 V15 Top1/Top5、分数、决策、置信度或化合物计数；跨极性 MS2 相似度仅作辅助注释。'}</span>
      </div>
    </div>`;
  }

  function dualPolarityShadowSection(task) {
    const english = isEnglishInterface();
    const currentMode = exactDualPolarity(task);
    const eligible = eligibleOppositePolarityTasks(task);
    const taskReady = isCompletedDualPolarityTask(task);
    const optionLabel = english
      ? 'Select the explicitly confirmed opposite-polarity task…'
      : '手动选择已确认的相反极性任务…';
    const options = eligible.map(candidate => {
      const stamp = candidate.upload_time
        ? String(candidate.upload_time).replace('T', ' ').replace(/\.\d+.*$/, '').slice(0, 16)
        : (english ? 'time not recorded' : '时间未记录');
      return `<option value="${esc(candidate.id)}">${esc(candidate.filename || candidate.id)} · ${esc(exactDualPolarity(candidate))} · ${esc(stamp)}</option>`;
    }).join('');
    const unavailable = !taskReady
      ? (english
        ? 'This tool requires a completed task whose saved ion_mode is exactly POS or NEG. MIXED, unknown and unfinished tasks are rejected.'
        : '该工具只接受状态为 completed 且保存 ion_mode 精确为 POS 或 NEG 的任务；MIXED、未知极性和未完成任务均拒绝。')
      : (!eligible.length
        ? (english
          ? `No completed ${currentMode === 'POS' ? 'NEG' : 'POS'} task is available for explicit selection.`
          : `当前没有可供手动选择的已完成 ${currentMode === 'POS' ? 'NEG' : 'POS'} 任务。`)
        : '');
    return `<section id="v16DualPolarityShadowCard" class="card v16-dual-card" data-current-task-id="${esc(task?.id || '')}" data-current-polarity="${esc(currentMode)}">
      <div class="v16-dual-card-header">
        <div><div class="card-title">${english ? 'V16 dual-polarity shadow evidence' : 'V16 双极性影子证据'}</div><div class="card-subtitle">${english
          ? 'Candidate-independent neutral-mass and retention-time links are frozen before candidate support is attached.'
          : '先按候选无关的中性质量与保留时间冻结连接，再附加候选支持。'}</div></div>
        <span class="v16-dual-shadow-badge">SHADOW · READ ONLY</span>
      </div>
      <div class="v16-dual-boundary"><b>${english ? 'No production effect' : '不影响生产结果'}</b><span>${english
        ? 'This analysis never auto-pairs filenames and cannot change V15 Top1, scores, confidence or compound counts.'
        : '系统不会按文件名自动配对；本分析不会改变 V15 Top1、分数、置信度和化合物计数。'}</span></div>
      ${unavailable ? `<div class="v16-dual-unavailable">${esc(unavailable)}</div>` : `<div class="v16-dual-form">
        <label class="v16-dual-field"><span>${english ? `Opposite-polarity task (${currentMode === 'POS' ? 'NEG' : 'POS'})` : `相反极性任务（${currentMode === 'POS' ? 'NEG' : 'POS'}）`}</span><select id="v16DualPolarityTask" class="form-select" onchange="SampleOverview.updateDualPolarityReadiness()"><option value="">${optionLabel}</option>${options}</select><small>${english ? 'Nothing is preselected; filename similarity is not pairing evidence.' : '默认不选中任何任务；文件名相似不构成配对证据。'}</small></label>
        <label class="v16-dual-field"><span>${english ? 'Pair label (optional)' : '配对标签（可选）'}</span><input id="v16DualPolarityPairLabel" class="form-input" maxlength="120" placeholder="${english ? 'e.g. sample-01' : '例如：样本-01'}"></label>
        <div class="v16-dual-confirmations">
          <label><input id="v16DualSameSampleConfirmed" type="checkbox" onchange="SampleOverview.updateDualPolarityReadiness()"><span>${english ? 'I confirm both tasks are injections of the same sample.' : '我确认两个任务来自同一样本。'}</span></label>
          <label><input id="v16DualSameLcMethodConfirmed" type="checkbox" onchange="SampleOverview.updateDualPolarityReadiness()"><span>${english ? 'I confirm both tasks used the same LC method.' : '我确认两个任务使用同一 LC 方法。'}</span></label>
        </div>
        <div class="v16-dual-actions"><button id="v16DualPolarityRun" class="btn btn-primary btn-sm" disabled onclick="SampleOverview.runDualPolarityShadow()">${english ? 'Build shadow evidence graph' : '构建影子证据图'}</button><span id="v16DualPolarityReadiness">${english ? 'Select one task and confirm both conditions.' : '请选择一个任务并确认两项条件。'}</span></div>
      </div>`}
      <div id="v16DualPolarityResult">${dualPolarityShadowResult(currentDualPolarityShadow)}</div>
    </section>`;
  }

  function updateDualPolarityReadiness() {
    const currentTask = currentBundle?.task || {};
    const selectedId = String(document.getElementById('v16DualPolarityTask')?.value || '');
    const selectedTask = eligibleOppositePolarityTasks(currentTask)
      .find(task => String(task.id) === selectedId);
    const sameSample = document.getElementById('v16DualSameSampleConfirmed')?.checked === true;
    const sameLcMethod = document.getElementById('v16DualSameLcMethodConfirmed')?.checked === true;
    const ready = !!selectedTask && sameSample && sameLcMethod && !dualPolarityShadowBusy;
    const button = document.getElementById('v16DualPolarityRun');
    if (button) {
      button.disabled = !ready;
      button.textContent = dualPolarityShadowBusy
        ? (isEnglishInterface() ? 'Building shadow graph…' : '正在构建影子图…')
        : (isEnglishInterface() ? 'Build shadow evidence graph' : '构建影子证据图');
    }
    const status = document.getElementById('v16DualPolarityReadiness');
    if (status) {
      status.textContent = ready
        ? (isEnglishInterface() ? 'Ready: pairing is explicit and both confirmations are recorded.' : '可以运行：任务已显式选择，两项确认均已记录。')
        : (isEnglishInterface() ? 'Select one task and confirm both conditions.' : '请选择一个任务并确认两项条件。');
    }
  }

  async function runDualPolarityShadow() {
    const currentTask = currentBundle?.task || {};
    const currentMode = exactDualPolarity(currentTask);
    const selectedId = String(document.getElementById('v16DualPolarityTask')?.value || '');
    const selectedTask = eligibleOppositePolarityTasks(currentTask)
      .find(task => String(task.id) === selectedId);
    const sameSample = document.getElementById('v16DualSameSampleConfirmed')?.checked === true;
    const sameLcMethod = document.getElementById('v16DualSameLcMethodConfirmed')?.checked === true;
    if (!isCompletedDualPolarityTask(currentTask) || !selectedTask) {
      Toast.show(isEnglishInterface()
        ? 'Select a completed task with the opposite saved polarity.'
        : '请手动选择一个已完成且保存极性相反的任务。', 'warning');
      return;
    }
    if (!sameSample || !sameLcMethod) {
      Toast.show(isEnglishInterface()
        ? 'Confirm both the sample identity and LC method before running.'
        : '运行前必须同时确认同一样本和同一 LC 方法。', 'warning');
      return;
    }
    const pairLabelInput = String(
      document.getElementById('v16DualPolarityPairLabel')?.value || ''
    ).trim();
    const payload = {
      positive_task_id: currentMode === 'POS' ? String(currentTask.id) : String(selectedTask.id),
      negative_task_id: currentMode === 'NEG' ? String(currentTask.id) : String(selectedTask.id),
      same_sample_confirmed: true,
      same_lc_method_confirmed: true,
      pair_label: pairLabelInput || `${String(currentTask.id)} + ${String(selectedTask.id)}`,
    };
    dualPolarityShadowBusy = true;
    updateDualPolarityReadiness();
    try {
      currentDualPolarityShadow = await API.post(
        '/workflows/dual-polarity-shadow',
        payload,
      );
      const result = document.getElementById('v16DualPolarityResult');
      if (result) result.innerHTML = dualPolarityShadowResult(currentDualPolarityShadow);
      Toast.show(isEnglishInterface()
        ? 'V16 shadow evidence graph completed; V15 production results are unchanged.'
        : 'V16 影子证据图已完成；V15 生产结果保持不变。', 'success');
    } catch (error) {
      Toast.show(`${isEnglishInterface() ? 'Dual-polarity shadow analysis failed' : '双极性影子分析失败'}：${error.message}`, 'error');
    } finally {
      dualPolarityShadowBusy = false;
      updateDualPolarityReadiness();
    }
  }

  function referenceStandardQCSection(payload, sets) {
    const available = payload?.reference_standard_set && payload?.summary;
    if (!available) {
      const options = (sets || []).map(item => (
        `<option value="${esc(item.id)}">${esc(item.label)} · ${integer(item.truth_identity_count)} 个</option>`
      )).join('');
      return `<section class="card" style="border-left:4px solid #0f766e;">
        <div class="card-title">混标／标准品回收审计</div>
        ${isEnglishInterface()
          ? '<p style="font-size:.78rem;line-height:1.55;margin:.45rem 0;color:var(--text-secondary);">No QC ground-truth set is linked to this task. The system does not infer one from the filename. Select a set explicitly to audit the stored blind-search results; <b>no reanalysis is required, and candidates, scores, and rankings remain unchanged</b>.</p>'
          : '<p style="font-size:.78rem;line-height:1.55;margin:.45rem 0;color:var(--text-secondary);">本任务尚未绑定质控真值集。系统不会按文件名猜测；可在这里显式选择并对已保存的盲检索结果做后验审计，<b>无需重新分析，也不会改变候选、分数或排序</b>。</p>'}
        <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;">
          <select id="resultReferenceStandardSet" class="form-select" style="max-width:520px;"><option value="">选择真实混标集…</option>${options}</select>
          <button class="btn btn-primary btn-sm" onclick="SampleOverview.applyReferenceStandardQC()">应用质控集（不重跑）</button>
        </div>
      </section>`;
    }
    const s = payload.summary || {};
    const set = payload.reference_standard_set || {};
    const notObservable = Math.max(0, Number(s.declared_truth_identities || 0)
      - Number(s.single_injection_compatible_precursor_identities || 0));
    const rows = (payload.per_standard || []).map(item => {
      const best = item.best_correct_top1 || {};
      const observable = !!item.observable_in_injection;
      return `<tr>
        <td><b>${esc(item.name_en || item.name_zh || '—')}</b><small style="display:block;color:var(--text-secondary);">${esc(item.name_zh || '')}</small></td>
        <td>${observable ? '可观测' : '本针无兼容 MS2 事件'}</td>
        <td>${!observable ? '不适用' : (item.saved_topk_connectivity ? `已找回 · 最佳第 ${integer(item.best_saved_rank_connectivity)}` : '未进入已保存候选')}</td>
        <td>${!observable ? '不适用' : (item.top1_connectivity_any_scan ? '是' : '否')}</td>
        <td>${!observable ? '不适用' : (item.reported_scope_supported_any_scan ? '支持' : '不支持')}</td>
        <td>${!observable ? '不适用' : (item.high_reportable_scope_supported_any_scan ? 'High' : '未达 High')}</td>
        <td>${best.candidate_name ? `${esc(best.candidate_name)} · cos ${score(best.cosine, 3)}` : '—'}</td>
      </tr>`;
    }).join('');
    const declared = Number(s.declared_truth_identities || 0);
    const observable = Number(s.single_injection_compatible_precursor_identities || 0);
    const ratio = value => `${integer(value)}/${integer(observable)}`;
    return `<section class="card" style="border-left:4px solid #0f766e;">
      <div class="card-header"><div><div class="card-title">${esc(set.label || set.id)} · 混标回收审计</div><div class="card-subtitle">后验质控真值不参与候选生成、评分或排序；${esc(payload.scope_zh || payload.scope || '')}</div></div></div>
      <section class="sample-result-metrics" aria-label="混标回收指标" style="margin-top:.7rem;">
        ${statusCard(declared, '声明标准品', '混标真值总数', 'total')}
        ${statusCard(`${integer(observable)}/${integer(declared)}`, '本针可观测', `${notObservable} 个在当前极性/选中 MS2 中不可观测`, 'total')}
        ${statusCard(ratio(s.observable_saved_topk_connectivity), '候选层找回', `${score(s.observable_saved_topk_connectivity_pct, 1)}% / 可观测标准品`, 'high')}
        ${statusCard(ratio(s.observable_top1_connectivity), 'Top1 找回', `${score(s.top1_connectivity_pct_observable, 1)}% / 可观测标准品`, 'high')}
        ${statusCard(ratio(s.observable_reported_scope_supported), '可报告范围找回', `${score(s.reported_scope_supported_pct_observable, 1)}% / 可观测标准品`, 'confirmed')}
        ${statusCard(ratio(s.observable_high_reportable_scope_supported), '达到 High 范围', `${score(s.high_reportable_scope_supported_pct_observable, 1)}% / 可观测标准品`, 'isomer')}
      </section>
      <div class="sample-result-scope" style="margin-top:.65rem;"><b>分母说明</b><span>${esc(payload.denominator_warning_zh || payload.denominator_warning || '')}</span></div>
      <details style="margin-top:.7rem;"><summary style="cursor:pointer;font-weight:700;color:#115e59;">查看 ${integer(s.declared_truth_identities)} 个标准品逐项结果</summary>
        <div class="table-wrap" style="margin-top:.55rem;"><table><thead><tr><th>标准品</th><th>本针可观测</th><th>候选层</th><th>Top1</th><th>可报告范围</th><th>High 范围</th><th>最佳正确候选</th></tr></thead><tbody>${rows}</tbody></table></div>
      </details>
    </section>`;
  }

  function statusCard(value, label, description, tone) {
    const displayValue = typeof value === 'string' ? esc(value) : integer(value);
    return `<div class="sample-result-metric metric-${tone}"><strong>${displayValue}</strong><span>${esc(label)}</span><small>${esc(description)}</small></div>`;
  }

  function candidateIdentitySummarySection(bundle) {
    const summaries = SampleResultData.candidateIdentitySummaries(
      bundle.featureGroups || bundle.groups || [], bundle.identificationRecords || []
    );
    if (!summaries.length) return '';
    const visible = summaries.slice(0, 50);
    const rows = visible.map((item, index) => {
      const localizedName = interfaceLanguage() === 'en' || !item.compoundNameCn
        ? item.compoundName
        : `${item.compoundName}（${item.compoundNameCn}）`;
      const status = item.bestStatus || { key: 'low', label: '待复核候选' };
      const rt = item.rtMin === null
        ? 'RT 未记录'
        : (item.rtMax !== null && Math.abs(item.rtMax - item.rtMin) > 0.001
          ? `RT ${item.rtMin.toFixed(3)}–${item.rtMax.toFixed(3)} min`
          : `RT ${item.rtMin.toFixed(3)} min`);
      const linkedText = item.linkedFeatureCount
        ? `${integer(item.linkedFeatureCount)} 个连续 MS1 特征已关联`
        : '无连续 MS1 特征关联';
      const fallbackText = item.fallbackEventCount
        ? `${integer(item.fallbackEventCount)} 个未关联 MS1 的 MS2 采集事件`
        : '无未关联 MS1 事件';
      return `<tr>
        <td><span class="candidate-summary-rank">${index + 1}</span></td>
        <td class="candidate-summary-identity"><b data-scientific-value>${esc(localizedName || 'Unnamed candidate')}</b><small>${esc(item.formula || '分子式未记录')}</small><small class="mono">${esc(item.inchikey || 'InChIKey 未记录')}</small>${item.curatedNatural ? '<span class="candidate-summary-curated">已核验天然产物身份记录</span>' : ''}</td>
        <td><span class="sample-evidence-badge evidence-${esc(status.key || 'low')}">${esc(status.label || '待复核候选')}</span><small>作为本针排序先导 ${integer(item.rankingLeadCount)} / ${integer(item.featureResultCount)} 条</small><small>${integer(item.ambiguousFeatureCount)} 条位于未解竞争候选组</small></td>
        <td><b>${esc(linkedText)}</b><small>${esc(fallbackText)}</small><small>这里汇总的是候选出现记录，不把 MS2 重复采集写成独立检出次数。</small></td>
        <td class="candidate-summary-score"><b>碎片余弦 ${item.bestCosine === null ? '—' : item.bestCosine.toFixed(3)}</b><small>谱熵相似度 ${item.bestEntropy === null ? '—' : item.bestEntropy.toFixed(3)} · 反向覆盖 ${item.bestReverseMatch === null ? '—' : item.bestReverseMatch.toFixed(3)}</small><small>${integer(item.fragmentEvidenceScanCount)} 条扫描级碎片证据 · 实体排序分 ${item.bestScore === null ? '—' : item.bestScore.toFixed(3)}</small><small>最小 |质量误差| ${item.bestAbsMassPpm === null ? '—' : `${Math.abs(item.bestAbsMassPpm).toFixed(2)} ppm`}</small></td>
        <td><b>${esc(rt)}</b><small>${integer(item.supportingSpectra)} 张组级关联 MS2；竞争组内不能全部归给单一候选。</small></td>
      </tr>`;
    }).join('');
    return `<section class="card candidate-identity-summary" aria-label="候选身份汇总">
      <div class="sample-compound-header"><div><div class="card-title">候选身份汇总（同一二维连接不重复列名）</div><div class="card-subtitle">按候选二维连接与分子式汇总，并保留每个身份关联的离子特征数、MS2 数和保留时间范围。CE 只显示采集条件，不参与这里的排序或扣分；碎片谱相似度单列。</div></div><span class="candidate-summary-count">${integer(visible.length)} / ${integer(summaries.length)} 个候选身份</span></div>
      <div class="candidate-summary-boundary"><b>结果边界</b><span>连续 MS1 特征已关联时可定位色谱峰；仅有未关联 MS1 的 MS2 采集事件时仍保留谱库候选，但不能把每次采集写成一次独立化合物检出。</span></div>
      <div class="table-wrap candidate-summary-table-wrap"><table class="candidate-summary-table"><thead><tr><th>#</th><th>候选身份</th><th>最高证据范围</th><th>色谱／采集支持</th><th>碎片匹配与排序</th><th>保留时间范围</th></tr></thead><tbody>${rows}</tbody></table></div>
      ${summaries.length > visible.length ? `<p class="sample-table-count">这里只列出前 ${integer(visible.length)} 个汇总身份；下方去重表保留全部候选身份，并把关联特征作为支持证据计数。</p>` : ''}
    </section>`;
  }

  function renderOverview(bundle) {
    const area = document.getElementById('sampleOverviewArea');
    if (!area) return;
    const task = bundle.task || {};
    const counts = summaryNumbers(bundle);
    const matchedRate = counts.totalSpectra ? (counts.matchedSpectra / counts.totalSpectra * 100).toFixed(1) : '0.0';
    const warningCount = Array.isArray(task.analysis_warnings) ? task.analysis_warnings.length : 0;
    const chromatogramHtml = window.SampleReport
      && typeof window.SampleReport.summaryChromatograms === 'function'
      ? window.SampleReport.summaryChromatograms(bundle.groups, currentChromatogram)
      : '';
    const detailedSummaryHtml = bundle.savedOnly
      ? `<div class="vr-source-note">${isEnglishInterface() ? 'This view reads saved scan-level candidates. Current entity conclusions and their High/confirmation counts have not been computed.' : '本视图读取已保存的谱级候选；当前实体结论及其 High／确证数量尚未计算。'}</div><p class="vr-plot-note">${esc(task.filename || '')} · ${integer(counts.totalSpectra)} MS2 · ${integer(counts.matchedSpectra)} ${isEnglishInterface() ? 'scans with candidates' : '张谱有候选'}</p>`
      : `
      ${configurationBanner(task)}
      ${resultHeadline(bundle, counts)}
      ${window.ConfidenceGuide
        ? window.ConfidenceGuide.render(bundle.confidenceGuide, {
          mode: 'web', id: 'sample-overview-confidence-guide',
          language: window.MSSLanguage?.getLanguage?.() || 'zh',
        })
        : ''}
      ${evitwinSection(bundle.evitwin)}
      ${replicateConsistencySection(bundle.replicateConsistency)}
      <section class="sample-result-metrics" aria-label="样本鉴定摘要">
        ${statusCard(SampleResultData.confidenceAxes(bundle).libraryHigh ?? (isEnglishInterface() ? 'Unavailable' : '未提供'), isEnglishInterface() ? 'Single-spectrum library High' : '单谱谱库 High 候选', isEnglishInterface() ? 'Deduplicated candidate identities; confirmation assessed separately' : '按单谱候选身份去重；结构确证另列', 'high')}
        ${statusCard(counts.confirmed, '已确证存在', '满足正式结论条件的化合物', 'confirmed')}
        ${statusCard(counts.highScopeGroups, '实体 High 级鉴定组', '通过当前实体规则的身份组数', 'high')}
        ${statusCard(counts.high, 'High 单一二维连接', '可报告到单一二维连接结构', 'high')}
        ${statusCard(counts.isomerHigh, 'High 异构体组', '只能报告到同分子式候选组', 'isomer')}
        ${statusCard(counts.isomer - counts.isomerHigh, '中等证据异构体组', '实体规则未通过 High；单谱证据另列', 'total')}
        ${statusCard(counts.medium, '实体 Medium', '仅统计已保存并完成实体评估的 Medium 等级', 'total')}
        ${statusCard(counts.candidate, '保存候选（实体未分级）', '不借用单谱等级推高实体等级', 'total')}
        ${statusCard(counts.background, '空白／背景相关', '谱库命中存在，但不计作样本化合物', 'contaminant')}
        ${statusCard(counts.unresolved, '未解析母离子事件', isEnglishInterface()
          ? `${integer(counts.unresolvedFormulaConflict)} known-formula conflicts; ${integer(counts.unresolvedIncomplete)} incomplete candidate-formula assignments; all are excluded from the independent compound count`
          : `${integer(counts.unresolvedFormulaConflict)} 个已知分子式冲突；${integer(counts.unresolvedIncomplete)} 个候选分子式记录不完整；均不计作独立化合物`, 'low')}
        ${statusCard(counts.insufficient, '证据不足', '保留线索，不进入主要结果', 'low')}
        ${statusCard(counts.unannotated, '未注释特征', '已检出离子特征，但普通谱库没有可用结构候选', 'low')}
      </section>
      <section class="sample-result-context">
        <div><span>样本文件</span><b>${esc(task.filename || '—')}</b></div>
        <div><span>MS2 谱图</span><b>${integer(counts.totalSpectra)}</b></div>
        <div><span>获得候选的 MS2</span><b>${integer(counts.matchedSpectra)}${isEnglishInterface() ? ` (${matchedRate}%)` : `（${matchedRate}%）`}</b></div>
        <div><span>分析提示</span><b>${warningCount
          ? (isEnglishInterface()
            ? `${integer(warningCount)} analysis notice${warningCount === 1 ? '' : 's'}`
            : `${integer(warningCount)} 项需注意`)
          : (isEnglishInterface() ? 'No additional warnings' : '无额外警告')}</b></div>
      </section>
      ${dualPolarityShadowSection(task)}
      ${coreAnalysisSummary(bundle)}
      ${referenceStandardQCSection(currentReferenceStandardQC, currentReferenceStandardSets)}
      ${targetedResultsSection(currentExpectedCompounds)}
      ${candidateIdentitySummarySection(bundle)}
      ${Results.renderGovernanceBlockedCard(
        bundle.governanceBlockedCandidates || [],
        String(task.id || currentTaskId || ''),
        bundle.governanceBlockedCount,
        bundle.governanceBlockedHasMore
      )}
    `;
    area.innerHTML = `
      ${typeof FrozenTaskUI !== 'undefined' ? FrozenTaskUI.panel(task) : ''}
      ${!task.completed_batch && typeof VisualResults !== 'undefined' ? VisualResults.summary(bundle, counts, isEnglishInterface()) : ''}
      ${chromatogramHtml}
      <section class="card sample-compound-card">
        <div class="sample-compound-header">
          <div><div class="card-title">${task.completed_batch ? '原生MS/MS候选、未解析事件与QC证据' : '去重候选身份'}</div><div class="card-subtitle">${task.completed_batch ? '已读取全部保存结论；按身份和前体事件展示，XBT联合筛选另列。来源未核实的候选保留为谱库假设。' : (isEnglishInterface() ? 'Candidate structures, spectral similarity and measured locations' : '候选结构、碎片相似度与实测位置')}</div></div>
          <div class="sample-compound-search"><input id="sampleCompoundSearch" class="form-input" type="search" placeholder="搜索数据库名称、分子式或完整 InChIKey" oninput="SampleOverview.setSearch(this.value)"></div>
        </div>
        <div class="vr-controls">
          <div class="vr-view-switch" role="group" aria-label="结果展示方式">
            <button data-visual-view="cards" aria-pressed="${currentVisualView === 'cards'}" onclick="SampleOverview.setVisualView('cards')">${isEnglishInterface() ? 'Structure cards' : '结构卡片'}</button>
            <button data-visual-view="table" aria-pressed="${currentVisualView === 'table'}" onclick="SampleOverview.setVisualView('table')">${isEnglishInterface() ? 'Compact table' : '紧凑表格'}</button>
            <button data-visual-view="audit" aria-pressed="${currentVisualView === 'audit'}" onclick="SampleOverview.setVisualView('audit')">${isEnglishInterface() ? 'All fields' : '详细字段'}</button>
          </div>
          <label class="vr-sort">${isEnglishInterface() ? 'Display order' : '展示顺序'} <select onchange="SampleOverview.setVisualSort(this.value)">
            <option value="original" ${currentVisualSort === 'original' ? 'selected' : ''}>${isEnglishInterface() ? 'Saved result order' : '原结果顺序'}</option>
            <option value="cosine" ${currentVisualSort === 'cosine' ? 'selected' : ''}>${isEnglishInterface() ? 'Cosine ↓' : '碎片余弦 ↓'}</option>
            <option value="rt" ${currentVisualSort === 'rt' ? 'selected' : ''}>RT ↑</option>
            <option value="mz" ${currentVisualSort === 'mz' ? 'selected' : ''}>m/z ↑</option>
          </select></label>
        </div>
        <div class="sample-filter-bar" role="group" aria-label="化合物结论筛选">
          ${bundle.savedOnly ? `
          <button class="sample-filter active" data-overview-filter="all" onclick="SampleOverview.setFilter('all')">${isEnglishInterface() ? 'All saved candidates' : '全部保存候选'} <span>${integer(bundle.groups.length)}</span></button>
          <button class="sample-filter" data-overview-filter="contaminant" onclick="SampleOverview.setFilter('contaminant')">${isEnglishInterface() ? 'External / role review' : '外源／身份复核'} <span>${integer(counts.externalReview)}</span></button>
          ` : `
          ${task.completed_batch ? '' : `<button class="sample-filter" data-overview-filter="main" onclick="SampleOverview.setFilter('main')">${isEnglishInterface() ? 'Report candidates' : '正文候选'} <span>${integer(counts.plantMain)}</span></button>`}
          <button class="sample-filter" data-overview-filter="confirmed" onclick="SampleOverview.setFilter('confirmed')">已确证 <span>${integer(counts.confirmed)}</span></button>
          <button class="sample-filter" data-overview-filter="high" onclick="SampleOverview.setFilter('high')">High 单一二维连接 <span>${integer(counts.high)}</span></button>
          <button class="sample-filter" data-overview-filter="isomer" onclick="SampleOverview.setFilter('isomer')">High 异构体组 <span>${integer(counts.isomerHigh)}</span></button>
          ${bundle.replicateConsistency?.status === 'available' ? `<button class="sample-filter" data-overview-filter="replicate-high" onclick="SampleOverview.setFilter('replicate-high')">${isEnglishInterface() ? 'Repeat-robust High' : '重复稳健 High'} <span>${integer(bundle.replicateConsistency?.summary?.robust_high_count)}</span></button>` : ''}
          <button class="sample-filter" data-overview-filter="medium" onclick="SampleOverview.setFilter('medium')">实体 Medium <span>${integer(counts.medium)}</span></button>
          <button class="sample-filter" data-overview-filter="candidate" onclick="SampleOverview.setFilter('candidate')">保存候选（实体未分级） <span>${integer(counts.candidate)}</span></button>
          <button class="sample-filter" data-overview-filter="background" onclick="SampleOverview.setFilter('background')">空白／背景 <span>${integer(counts.background)}</span></button>
          <button class="sample-filter" data-overview-filter="unresolved_event" onclick="SampleOverview.setFilter('unresolved_event')">未解析事件 <span>${integer(counts.unresolved)}</span></button>
          <button class="sample-filter" data-overview-filter="low" onclick="SampleOverview.setFilter('low')">证据不足 <span>${integer(counts.insufficient)}</span></button>
          <button class="sample-filter" data-overview-filter="unannotated" onclick="SampleOverview.setFilter('unannotated')">未注释特征 <span>${integer(counts.unannotated)}</span></button>
          <button class="sample-filter active" data-overview-filter="all" onclick="SampleOverview.setFilter('all')">全部去重身份／事件 <span>${integer(bundle.groups.length)}</span></button>
          <button class="sample-filter" data-overview-filter="contaminant" onclick="SampleOverview.setFilter('contaminant')">${isEnglishInterface() ? 'External / identity review' : '外源／身份复核'} <span>${integer(counts.externalReview)}</span></button>
          `}
        </div>
        <div id="sampleCompoundTable"></div>
      </section>
      <details class="vr-disclosure"><summary>${isEnglishInterface() ? 'All counts, analysis context & additional results' : '全部计数、分析条件与附加结果'}</summary><div class="vr-disclosure-body">${task.completed_batch ? targetedResultsSection(currentExpectedCompounds) : detailedSummaryHtml}</div></details>`;
    renderCompoundTable();
  }

  function matchesAnalysisFilter(group) {
    const filter = String(currentAnalysisFilter || 'all');
    if (filter === 'all') return true;
    const match = group.rep?.match || {};
    const analysis = groupAnalysis(group);
    if (filter.startsWith('class-state:')) {
      return classificationState(match) === filter.split(':')[1];
    }
    if (filter.startsWith('class:')) {
      return rawClassification(match) === decodeURIComponent(filter.slice(6));
    }
    if (filter === 'biosource:any') return analysis.biosyntheticPathways.length > 0;
    if (filter === 'biosource:none') return analysis.biosyntheticPathways.length === 0;
    if (filter.startsWith('biosource:')) {
      const family = decodeURIComponent(filter.slice(10));
      return analysis.biosyntheticPathways.some(value => biosyntheticFamily(value) === family);
    }
    if (filter === 'expert:positive') return analysis.expertPositive;
    if (filter === 'expert:no_evidence') return analysis.hasExpert && !analysis.expertPositive;
    if (filter === 'expert:unavailable') return !analysis.hasExpert;
    if (filter === 'fragment:path') return analysis.hasFragmentation;
    if (filter === 'fragment:reference') return analysis.hasLibraryPathwayReference;
    if (filter === 'fragment:neutral') return analysis.hasNeutralLoss;
    if (filter === 'fragment:mismatch') return analysis.pathwayMismatch;
    return true;
  }

  function analysisFilterLabel(filter) {
    const labels = {
      'class-state:reliable': '分类：可靠来源',
      'class-state:tentative': '分类：待核查线索',
      'class-state:unresolved': '分类：未可靠判类',
      'biosource:any': '生源：有命名规则',
      'biosource:none': '生源：暂无命名规则',
      'expert:positive': '规则证据：候选组至少一张谱有正向支持',
      'expert:no_evidence': '规则证据：已评价但无正向支持',
      'expert:unavailable': '规则证据：不可用／无来源',
      'fragment:path': '碎裂：已重建样品步骤/产物',
      'fragment:reference': '碎裂：有候选谱库参考',
      'fragment:neutral': '碎裂：有中性丢失结果',
      'fragment:mismatch': '碎裂：馆藏路径极性不符',
    };
    if (labels[filter]) return labels[filter];
    if (filter.startsWith('class:')) return `分类：${decodeURIComponent(filter.slice(6))}`;
    if (filter.startsWith('biosource:')) return `生源：${decodeURIComponent(filter.slice(10))}`;
    return '多层分析筛选';
  }

  function fragmentationStepLabel(step) {
    const english = isEnglishInterface();
    const from = Number(step?.from_mz);
    const to = Number(step?.to_mz);
    const transition = Number.isFinite(from) && Number.isFinite(to)
      ? `${from.toFixed(4)} → ${to.toFixed(4)}` : '';
    const loss = String(
      (english ? step?.neutral_loss_en : step?.neutral_loss)
      || (english ? step?.neutral_loss : step?.neutral_loss_en) || step?.formula
      || step?.mechanism || step?.type || ''
    ).trim();
    return [transition, loss
      ? (english ? `proposed loss/transformation: ${loss}` : `推测丢失／转化 ${loss}`)
      : ''].filter(Boolean).join(english ? ', ' : '，');
  }

  function renderCompoundTable() {
    const container = document.getElementById('sampleCompoundTable');
    if (!container || !currentBundle) return;
    const search = currentSearch.trim().toLowerCase();
    const filtered = currentBundle.groups.filter(group => {
      if (currentFilter === 'main' && (!group.status.countedAsPlantCandidate
          || SampleResultData.explicitExternalReview(group))) return false;
      if (currentFilter === 'confirmed' && group.status.key !== 'confirmed') return false;
      if (currentFilter === 'high' && group.status.key !== 'high') return false;
      if (currentFilter === 'isomer' && !group.status.hasIsomerHigh) return false;
      if (currentFilter === 'isomer_all' && group.status.key !== 'isomer') return false;
      if (currentFilter === 'replicate-high' && !group.replicateConsistency?.robust_high) return false;
      if (['medium', 'candidate', 'ungraded', 'background', 'unresolved_event', 'low', 'unannotated', 'saved_candidate'].includes(currentFilter)
          && group.status.key !== currentFilter) return false;
      if (currentFilter === 'contaminant'
          && !SampleResultData.explicitExternalReview(group)) return false;
      if (!matchesAnalysisFilter(group)) return false;
      if (!search) return true;
      const match = group.rep.match || {};
      const analysis = groupAnalysis(group);
      return CompoundSearch.matches({
        names: [match.compound_name, match.source_compound_name, ...(group.candidateNames || [])],
        terms: [
          ...(group.compound?.candidate_formulas || []), match.formula, match.inchikey,
          match.display_class, match.classification, match.contaminant_category,
          match.compound_role, match.screening_domain, match.pesticide_use,
          ...analysis.biosyntheticPathways,
        ],
        mzs: [group.rep?.precursor_mz],
      }, search);
    });
    const active = document.getElementById('activeAnalysisFilter');
    if (active) {
      active.hidden = currentAnalysisFilter === 'all';
      active.innerHTML = currentAnalysisFilter === 'all' ? ''
        : `<span>当前分析筛选：<b>${esc(analysisFilterLabel(currentAnalysisFilter))}</b></span><button onclick="SampleOverview.setAnalysisFilter('all')">清除筛选 ×</button>`;
    }
    const sorted = [...filtered];
    if (typeof VisualResults !== 'undefined') {
      const finite = VisualResults.finite;
      const numeric = group => currentVisualSort === 'rt'
        ? VisualResults.rt(group.rep)
        : finite(currentVisualSort === 'mz' ? group.rep?.precursor_mz : group.rep?.match?.cosine);
      if (currentVisualSort !== 'original') sorted.sort((a, b) => {
        const av = numeric(a), bv = numeric(b);
        if (av == null) return bv == null ? 0 : 1;
        if (bv == null) return -1;
        return currentVisualSort === 'cosine' ? bv - av : av - bv;
      });
    }
    if (currentVisualView !== 'audit' && typeof VisualResults !== 'undefined') {
      currentVisualPage = Math.max(0, Math.min(currentVisualPage, Math.ceil(sorted.length / 24) - 1));
      container.innerHTML = VisualResults.catalog(sorted, {
        taskId: currentTaskId, mode: currentVisualView, page: currentVisualPage,
        pageSize: 24, english: isEnglishInterface(),
      });
      VisualResults.mountStructures('sampleCompoundTable');
      return;
    }
    const visible = sorted.slice(0, 300);
    const taskBlank = currentBundle.task?.blank_summary || {};
    const taskHardBlankExclusion = (
      taskBlank.blank_operation === 'confirmed_presearch_background_exclusion'
      || taskBlank.filter_policy === 'hard_legacy'
    );
    const taskBlankExcluded = integer(taskBlank.excluded_spectra || 0);
    const english = isEnglishInterface();
    const rows = visible.map(group => {
      const row = group.rep;
      const match = row.match || {};
      const hasIsomerCompetition = SampleResultData.groupHasIsomerCompetition(group);
      const isUnresolvedEvent = group.status.key === 'unresolved_event';
      const unresolvedKind = isUnresolvedEvent ? unresolvedEventKind(group) : '';
      const incompleteFormulaEvent = unresolvedKind === 'incomplete_formula';
      const leadName = String(group.leadCandidate?.compound_name
        || match.compound_name || 'Unnamed candidate');
      const name = String(
        group.compound?.display_name
        || (isUnresolvedEvent ? '未解析母离子事件' : leadName)
      );
      const scopeLabels = {
        confirmed: '单一化合物，已满足正式确证条件',
        high: '单一二维连接结构，High 级检出',
        isomer: 'High 级同分子式异构体组；证据仅支持组级，不释放单体名称',
        medium: '候选结构，仍需补充证据',
        background: '空白背景比较提示（未执行跨针 MS/MS 扣除）',
        unresolved_event: english
          ? (incompleteFormulaEvent
            ? 'Candidate structures disagree and at least one candidate lacks a usable molecular formula; excluded from the independent compound count'
            : 'Known candidate molecular formulae conflict; excluded from the independent compound count')
          : (incompleteFormulaEvent
            ? '候选结构不一致且至少一个候选缺少可用分子式；不计作独立化合物'
            : '已知候选分子式互相冲突的母离子事件；不计作独立化合物'),
        unannotated: english
          ? 'Measured ion feature with no usable library candidate'
          : '已检出离子特征，普通谱库没有可用结构候选',
        candidate: english
          ? 'Saved library candidate; entity grade not recorded'
          : '已保存谱库候选；实体等级未记录',
        ungraded: english ? 'Entity grade not recorded' : '实体等级未记录',
        low: '低证据谱库线索',
      };
      const scope = group.status.key === 'isomer'
        ? (group.status.hasIsomerHigh
          ? 'High 级同分子式异构体组；证据仅支持组级，不释放单体名称'
          : '中等证据同分子式异构体组；证据仅支持组级，不释放单体名称')
        : (scopeLabels[group.status.key] || '谱库候选');
      const statusNote = isUnresolvedEvent && english
        ? (incompleteFormulaEvent
          ? 'Candidate structures disagree and at least one candidate lacks a usable molecular formula. The candidates are retained for review, but the event is excluded from the independent compound count.'
          : 'Known candidate molecular formulae conflict. The candidates are retained for review, but the event is excluded from the independent compound count.')
        : group.status.note;
      const rtMinutes = SampleResultData.retentionTimeMinutes(row);
      const rt = rtMinutes == null ? 'RT 未记录' : `RT ${rtMinutes.toFixed(3)} min`;
      const spectrumId = String(row.spectrum_id || '').trim();
      const contaminant = contaminantPresentation(match);
      const analysis = groupAnalysis(group);
      const rawBiosource = analysis.biosyntheticPathways[0]
        || (english ? biosourceReasonEnglish(analysis) : analysis.biosourceReasonZh);
      const biosource = (english ? englishInterfaceCopy(rawBiosource) : rawBiosource)
        || (english
          ? 'Current biosource knowledge is not covered (library-identity level unchanged)'
          : '当前生源知识未覆盖（不影响谱库身份等级）');
      const fragmentLabels = analysis.fragmentationStepDetails
        .map(fragmentationStepLabel).filter(Boolean).slice(0, 2);
      const fragmentation = analysis.productStructures
        ? `${integer(analysis.productStructures)} 个观测峰支持的候选产物结构${analysis.fragmentationSteps ? `；${integer(analysis.fragmentationSteps)} 个质量差步骤` : ''}`
        : fragmentLabels.length
        ? fragmentLabels.join(english ? '; ' : '；')
        : (analysis.hasFragmentation
          ? `已从本次样品重建 ${integer(analysis.fragmentationSteps)} 个碎裂步骤`
          : (analysis.hasLibraryPathwayReference
            ? '未重建样品碎裂步骤；仅有候选谱库碎裂特征参考'
            : (english
              ? 'Unavailable: current evidence is insufficient to propose a candidate-specific fragmentation pathway'
              : '未能推理：当前证据不足以形成候选特异碎裂途径')));
      const expertCoverage = analysis.expertCoverage || {};
      const expert = analysis.hasExpert
        ? `有支持证据 ${integer(expertCoverage.scored || 0)} 维；已评价未观察 ${integer(expertCoverage.no_evidence || 0)} 维；不适用 ${integer(expertCoverage.not_applicable || 0)} 维；未建立 ${integer(expertCoverage.unavailable || 0)} 维`
        : '规则证据评估不可用（没有逐维证据记录）';
      const sampleEvidence = group.sampleEvidence || {};
      const unverifiedBlank = !!group.conclusion?.unverified_blank_association;
      const blankText = unverifiedBlank
        ? (english ? 'Same-polarity blank warning (0 excluded)' : '同极性空白警示（0 扣除）')
        : taskHardBlankExclusion
        ? (english
          ? `Background exclusion applied; this result was retained (${taskBlankExcluded} task-wide exclusions)`
          : `已执行背景排除；本条保留（任务共扣除 ${taskBlankExcluded} 条）`)
        : group.status.key === 'background'
        ? (english ? 'Background-comparison flag (no cross-injection subtraction)' : '背景比较提示（未跨针相减）')
        : (sampleEvidence.blank_evidence_available
          ? (english ? 'Blank background compared (0 excluded)' : '空白背景已比较（0 扣除）')
          : (english ? 'Not evaluated (no blank evidence)' : '未评估（无空白证据）'));
      const blankDetail = unverifiedBlank
        ? (english
          ? 'Signal was also observed in an untyped same-polarity blank. Batch and chromatographic-method pairing have not been verified, so this is a warning only and does not remove the candidate.'
          : String(group.conclusion?.blank_warning || '同极性未定角色空白中亦观察到信号；批次／色谱方法配对未经核验，仅作警示，不据此剔除候选。'))
        : taskHardBlankExclusion
        ? (english
          ? `Before library searching, this task excluded ${taskBlankExcluded} blank-dominant MS2 scans using same-polarity m/z, RT and intensity evidence; the current result was retained. Raw exclusion records remain in the blank-audit layer, and fragment intensities were not subtracted peak by peak across injections.`
          : `本任务已在检索前按同极性 m/z/RT/强度证据排除 ${taskBlankExcluded} 张空白占优 MS2；当前这条结果未被排除。原始排除记录保留在空白审计层，未进行跨针碎片强度逐峰相减。`)
        : group.status.key === 'background'
        ? (english
          ? 'The current signal is blank-related in the background comparison and is therefore excluded from primary sample results. MS/MS was not subtracted across injections, and the saved candidate was not deleted.'
          : '当前信号在背景比较中与空白相关，因此不进入主要样品结果；系统没有跨针相减 MS/MS，也没有删除保存的候选。')
        : (sampleEvidence.blank_evidence_available
          ? (english
            ? 'Blank-background comparison is complete for this task; the category and sample-to-blank ratio are available in the evidence view. No cross-injection MS/MS subtraction was performed.'
            : '已完成当前任务空白背景比较；具体类别与样品／空白比见证据详情。本任务未执行跨针 MS/MS 扣除。')
          : (english
            ? 'No usable blank was provided; background has not been evaluated. This status does not mean that the blank check was passed.'
            : '未提供可用空白，背景尚未评估；该状态不等于已通过空白检查。'));
      const candidateCount = Math.max(
        Number(group.candidateMemberCount || 0),
        Number(match.isomer_group_size || 0),
      );
      const alternativeNames = isUnresolvedEvent
        ? (english
          ? `<small><b>Ranking lead retained for review:</b> <span data-scientific-value>${esc(leadName)}</span> (not used as a compound name)</small><small>${incompleteFormulaEvent ? `${integer(candidateCount || 1)} candidates disagree structurally and at least one lacks a usable formula assignment; known formulae: ${esc((group.compound?.candidate_formulas || []).join(' / ') || 'none recorded')}` : `${integer(candidateCount || 1)} ${pluralizedCandidate(candidateCount || 1)} ${Number(candidateCount || 1) === 1 ? 'has' : 'have'} conflicting known formulae: ${esc((group.compound?.candidate_formulas || []).join(' / ') || 'none recorded')}`}</small>`
          : `<small><b>仅供复核的排序先导：</b><span data-scientific-value>${esc(leadName)}</span>（不作为化合物名称）</small><small>${incompleteFormulaEvent ? `${integer(candidateCount || 1)} 个候选结构不一致，且至少一个候选缺少可用分子式；已知分子式：${esc((group.compound?.candidate_formulas || []).join(' / ') || '无')}` : `${integer(candidateCount || 1)} 个候选的已知分子式互相冲突：${esc((group.compound?.candidate_formulas || []).join(' / ') || '无')}`}</small>`)
        : (hasIsomerCompetition
          ? (english
            ? `<small><b>Top-ranked candidate in this injection:</b> <span data-scientific-value>${esc(leadName)}</span> (does not mean that one isomer has been resolved)</small><small>${integer(candidateCount || 1)} spectrally competing ${pluralizedCandidate(candidateCount || 1)} in this group; remaining names are available in the detailed view</small>`
            : `<small><b>本针 Top1 排序先导：</b><span data-scientific-value>${esc(leadName)}</span>（不代表单体已解析）</small><small>本组 ${integer(candidateCount || 1)} 个谱学竞争候选；其余名称不再挤在概览中</small>`)
          : '');
      const replicate = group.replicateConsistency || {};
      const replicateNote = replicate.injection_count
        ? `<small class="${replicate.robust_high ? 'sample-scope-badge' : ''}">${english ? 'Technical injections: detected' : '重复针：检出'} ${integer(replicate.detected_injection_count)} / ${integer(replicate.injection_count)} · High ${integer(replicate.high_injection_count)} / ${integer(replicate.injection_count)}${replicate.robust_high ? (english ? ' · Repeat-robust High' : ' · 重复稳健 High') : ''}</small>`
        : '';
      const cosine = Number(match.cosine);
      const unrated = english ? 'Unrated' : '未分级';
      const confidenceDisplay = isUnresolvedEvent
        ? (english
          ? (incompleteFormulaEvent
            ? 'Unresolved event (incomplete candidate-formula assignment; excluded from the compound count)'
            : 'Unresolved event (conflicting known candidate formulae; excluded from the compound count)')
          : (incompleteFormulaEvent
            ? '未解析事件（候选分子式记录不完整；不计作化合物）'
            : '未解析事件（已知候选分子式冲突；不计作化合物）'))
        : hasIsomerCompetition
        ? (english
          ? `${match.isomer_group_confidence || unrated} (isomer group) / ${match.exact_structure_confidence || match.confidence || unrated} (full-structure identity)`
          : `${match.isomer_group_confidence || unrated}（异构体组） / ${match.exact_structure_confidence || match.confidence || unrated}（完整结构身份）`)
        : (['candidate', 'ungraded'].includes(group.status.key)
          ? (english ? 'Unrated (entity)' : '未分级（实体）')
          : group.status.key === 'high'
          ? (english
            ? `${match.two_dimensional_connectivity_confidence || match.confidence || unrated} (2-D connectivity) / ${match.exact_structure_confidence || match.confidence || unrated} (full-structure identity)`
            : `${match.two_dimensional_connectivity_confidence || match.confidence || unrated}（二维连接） / ${match.exact_structure_confidence || match.confidence || unrated}（完整结构身份）`)
          : (match.two_dimensional_connectivity_confidence || match.exact_structure_confidence || match.confidence || unrated));
      const precursorGateNote = match.precursor_high_evidence_gate_required
        ? `High 已降为 Medium：${match.precursor_high_evidence_gate_reason_zh || match.precursor_high_evidence_gate_reason || '碎片配对稀疏且缺少独立前体支持'}`
        : '';
      const ceStatus = String(match.ce_comparability?.status || (
        group.status.key === 'unannotated' ? 'no_candidate' : 'not_recorded'
      ));
      const fragmentationStatus = String(
        row.fragmentation_status || row.fragmentation_annotation?.status
        || 'not_recorded'
      );
      const backgroundReliability = row.background_evidence || {};
      const featureId = String(group.compound?.feature_id || row.lcms_feature_id
        || row.precursor_feature_id || '未记录');
      const eventText = Number(group.precursorEventCount || 1) > 1
        ? `${integer(group.precursorEventCount)} 个前体事件 · ` : '';
      const peakText = group.peakClusters > 1
        ? `${integer(group.peakClusters)} 个色谱峰区` : `${integer(group.peakClusters || 1)} 个色谱峰区`;
      const formulaDisplay = isUnresolvedEvent
        ? ((group.compound?.candidate_formulas || []).join(' / ')
          || (english ? 'No usable candidate formula recorded' : '无可用候选分子式'))
        : (group.compound?.formula || match.formula || '分子式未记录');
      return `<tr onclick="SampleOverview.openEvidence(${group.index})" title="打开该候选身份的代表特征证据">
        <td class="sample-compound-name"><b${isUnresolvedEvent ? '' : ' data-scientific-value'}>${esc(name)}</b>${contaminant ? `<small class="compound-caution" title="${esc(contaminant.note || contaminant.label)}">${esc(contaminant.label)}</small>` : ''}<small>${esc(formulaDisplay)}</small>${alternativeNames}${replicateNote}<small>${esc(eventText)}<b>${integer(group.spectra)}</b> 张 MS2，已合并为 ${esc(peakText)}</small><button class="sample-candidate-detail-action" onclick="event.stopPropagation();SampleOverview.openEvidence(${group.index})">${english ? (isUnresolvedEvent ? `View all ${integer(candidateCount || 1)} unresolved candidates and evidence →` : (hasIsomerCompetition ? `View all ${integer(candidateCount || 1)} competing candidates and evidence →` : 'View evidence details →')) : (isUnresolvedEvent ? `查看全部 ${integer(candidateCount || 1)} 个未解析候选及证据 →` : (hasIsomerCompetition ? `查看全部 ${integer(candidateCount || 1)} 个竞争候选及证据 →` : '查看证据详情 →'))}</button></td>
        <td class="sample-identification-scope"><span class="sample-evidence-badge evidence-${group.status.key}">${esc(group.status.label)}</span><small class="sample-evidence-note">${esc(statusNote)}</small><small class="sample-evidence-note" data-entity-decision-reasons>${esc(SampleResultData.entityDecisionNote(group.entity, undefined, match.entity_high_evidence))}</small><span class="sample-scope-label">能够报告到</span><span class="sample-scope-badge">${esc(scope)}</span></td>
        <td class="sample-analysis-summary">
          <div><b>分类</b><span title="分类来源：${esc(match.classification_source || '未记录')}；分类可靠性与谱库身份匹配等级独立评估">${esc(classificationLabel(match))}</span></div>
          <div><b>生源背景</b><span title="${esc(english ? `${biosource}. Class-level context is independent of identity ranking.` : `${biosource}；类别级背景与身份排序独立`)}">${esc(biosource)}</span></div>
          <div><b>碎裂途径</b><span title="${esc(fragmentation)}">${esc(fragmentation)}</span></div>
          <div><b>CE 注释（不计分）</b><span title="只记录采集条件；不乘入碎片相似度、综合分、排序或 EviGate">${esc(ceStatus)}</span></div>
          <div><b>碎片信息</b><span>${esc(fragmentationStatus)} · ${integer(row.effective_fragment_count || 0)} 个有效峰</span></div>
          <div><b>背景三态</b><span title="新规则未排除或降权任何峰">B1 ${esc(backgroundReliability.B1_cross_precursor || 'unknown')} · B2 ${esc(backgroundReliability.B2_incoherent_chromatography || 'unknown')} · B3 ${esc(backgroundReliability.B3_blank_confirmed || 'unknown')} · ${esc(backgroundReliability.action || 'retain')}</span></div>
          <div><b>规则证据</b><span title="${esc(expert)}">${esc(expert)}</span></div>
        </td>
        <td class="sample-key-scores">
          <span><small>谱库证据级别（按可报告范围）</small><b>${esc(confidenceDisplay)}</b>${precursorGateNote ? `<small class="compound-caution" title="${esc(precursorGateNote)}">${esc(precursorGateNote)}</small>` : ''}</span>
          <span><small>谱图相似度</small><b>${Number.isFinite(cosine) ? cosine.toFixed(3) : '—'}</b></span>
          <span><small>质量误差</small><b>${esc(ppm(match.mass_ppm))}</b></span>
          <span class="sample-key-score-note" title="${esc(blankDetail)}"><small>样本／空白</small><b>${esc(blankText)}</b></span>
          <details onclick="event.stopPropagation()"><summary>评分构成（非概率）</summary><small>综合排序分 ${score(match.rank_score ?? match.final_score, 3)}；CE 注释、规则证据、分类、生源和碎裂解释均不参与身份排序。</small></details>
        </td>
        <td class="sample-raw-locator"><span class="mono">m/z ${Number(row.precursor_mz || 0).toFixed(5)}</span><small>${esc(rt)}</small><small class="mono">Feature ID ${esc(featureId)}</small><small>样品离子形式 ${esc(analysis.sampleIonForm || '未可靠判定')}${analysis.sampleIonFormSource ? ` · ${esc(analysis.sampleIonFormSource)}` : ''}</small>${analysis.candidateIonFormHypothesis ? `<small title="由当前候选结构/分子式与谱库前体反推；不是样品独立证据">候选离子假设 ${esc(analysis.candidateIonFormHypothesis)}（非独立证据）</small>` : ''}<small>${esc(peakText)} · ${integer(group.spectra)} 张 MS2</small><small>${spectrumId ? `代表谱 ${esc(spectrumId)}` : '代表谱编号未记录'}</small></td>
      </tr>`;
    }).join('');
    container.innerHTML = filtered.length ? `
      <div class="sample-table-count">显示 ${integer(visible.length)} / ${integer(filtered.length)} 个候选身份／未解析事件${filtered.length > visible.length ? '（请搜索或筛选以缩小范围）' : ''}</div>
      <div class="table-wrap sample-compound-table-wrap"><table class="sample-compound-table"><colgroup><col class="candidate-col"><col class="scope-col"><col class="analysis-col"><col class="score-col"><col class="locator-col"></colgroup><thead><tr><th>候选身份／未解析事件</th><th>鉴定结论</th><th>分类 / 生源 / 碎裂 / 规则证据</th><th>关键证据</th><th>检测位置（代表 MS2 原始定位）</th></tr></thead><tbody>${rows}</tbody></table></div>`
      : `<div class="sample-no-results">${isEnglishInterface() ? 'No compound results match the current filters.' : '当前筛选条件下没有对应的化合物结果。'}</div>`;
  }

  async function loadForTask() {
    const taskId = document.getElementById('resultTaskSelect')?.value;
    const area = document.getElementById('sampleOverviewArea');
    if (!taskId) {
      currentTaskId = null;
      currentBundle = null;
      currentChromatogram = null;
      currentExpectedCompounds = null;
      currentReferenceStandardQC = null;
      currentReferenceStandardSets = [];
      currentDualPolarityTasks = [];
      currentDualPolarityShadow = null;
      dualPolarityShadowBusy = false;
      if (area) area.innerHTML = '<div class="empty-state"><div class="empty-icon">🔬</div><div class="empty-title">选择样本以查看鉴定结果</div><div class="empty-desc">结果首页会优先列出化合物、关键评分与鉴定范围</div></div>';
      return;
    }
    currentTaskId = taskId;
    currentBundle = null;
    currentChromatogram = null;
    currentExpectedCompounds = null;
    currentReferenceStandardQC = null;
    currentReferenceStandardSets = [];
    currentDualPolarityTasks = [];
    currentDualPolarityShadow = null;
    dualPolarityShadowBusy = false;
    currentFilter = 'all';
    currentAnalysisFilter = 'all';
    currentSearch = '';
    currentVisualPage = 0;
    if (area) area.innerHTML = '<div class="sample-loading"><span class="spinner"></span><b>正在生成去重化合物结果…</b><small>合并重复 MS2，并保留分类、生源、规则证据与碎裂解释</small></div>';
    try {
      const [bundle, chromatogram, expectedCompounds, referenceStandardQC, referenceSets, completedTasks] = await Promise.all([
        SampleResultData.load(taskId),
        API.get(`/tasks/${encodeURIComponent(taskId)}/tic_overview`).catch(() => null),
        API.get(`/tasks/${encodeURIComponent(taskId)}/expected-compounds/groups?limit=1000`)
          .catch(error => ({ _request_error: error.message })),
        API.get(`/tasks/${encodeURIComponent(taskId)}/reference-standard-qc`)
          .catch(error => ({ _request_error: error.message })),
        API.get('/reference-standard-qc/sets')
          .then(payload => payload.sets || []).catch(() => []),
        API.get('/tasks?status=completed&projection=minimal')
          .catch(error => ({ tasks: [], _request_error: error.message })),
      ]);
      if (currentTaskId !== taskId) return;
      currentBundle = bundle;
      currentChromatogram = chromatogram;
      currentExpectedCompounds = expectedCompounds;
      currentReferenceStandardQC = referenceStandardQC;
      currentReferenceStandardSets = referenceSets;
      currentDualPolarityTasks = Array.isArray(completedTasks?.tasks)
        ? completedTasks.tasks.filter(task => (
          String(task?.status || '').trim().toLowerCase() === 'completed'
        ))
        : [];
      renderOverview(currentBundle);
    } catch (error) {
      if (currentTaskId !== taskId) return;
      if (area) area.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">结果加载失败</div><div class="empty-desc">${esc(error.message || error)}</div></div>`;
    }
  }

  async function applyReferenceStandardQC() {
    if (!currentTaskId) return;
    const setId = String(document.getElementById('resultReferenceStandardSet')?.value || '').trim();
    if (!setId) {
      Toast.show('请先明确选择这份文件真实对应的混标质控集', 'warning');
      return;
    }
    try {
      Toast.show('正在对已保存结果计算混标回收指标；不会重新检索…', 'info');
      currentReferenceStandardQC = await API.post(
        `/tasks/${encodeURIComponent(currentTaskId)}/reference-standard-qc`,
        { set_id: setId, mass_tolerance_ppm: 10.0 },
      );
      if (currentBundle?.task) {
        currentBundle.task.reference_standard_qc_summary = {
          status: 'available',
          set_id: setId,
          ...(currentReferenceStandardQC.summary || {}),
        };
      }
      renderOverview(currentBundle);
      Toast.show('混标回收审计已完成；盲检索结果未改变', 'success');
    } catch (error) {
      Toast.show(`混标回收审计失败：${error.message}`, 'error');
    }
  }

  function setFilter(filter) {
    currentFilter = filter || 'all';
    currentVisualPage = 0;
    document.querySelectorAll('[data-overview-filter]').forEach(button => {
      button.classList.toggle('active', button.dataset.overviewFilter === currentFilter);
    });
    renderCompoundTable();
  }

  function setSearch(value) {
    currentSearch = String(value || '');
    currentVisualPage = 0;
    renderCompoundTable();
  }

  function setAnalysisFilter(filter) {
    currentAnalysisFilter = String(filter || 'all');
    currentVisualPage = 0;
    renderCompoundTable();
    document.getElementById('sampleCompoundTable')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  function setVisualView(view) {
    currentVisualView = ['cards', 'table', 'audit'].includes(view) ? view : 'cards';
    currentVisualPage = 0;
    document.querySelectorAll('[data-visual-view]').forEach(button => {
      button.setAttribute('aria-pressed', String(button.dataset.visualView === currentVisualView));
    });
    renderCompoundTable();
  }

  function setVisualSort(sort) {
    currentVisualSort = ['original', 'cosine', 'rt', 'mz'].includes(sort) ? sort : 'original';
    currentVisualPage = 0;
    renderCompoundTable();
  }

  function changeVisualPage(delta) {
    currentVisualPage = Math.max(0, currentVisualPage + Number(delta || 0));
    renderCompoundTable();
    document.getElementById('sampleCompoundTable')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  async function openTask(taskId) {
    App.navigateTo('results');
    await refreshSelect(taskId);
    const select = document.getElementById('resultTaskSelect');
    if (!select) return;
    select.value = taskId;
    await loadForTask();
  }

  async function openEvidence(groupIndex) {
    const group = currentBundle?.groups.find(item => item.index === Number(groupIndex));
    const rawResultIndex = group?.rep?.result_index;
    const resultIndex = rawResultIndex == null || String(rawResultIndex).trim() === ''
      ? null : Number(rawResultIndex);
    if (!Number.isInteger(resultIndex) || resultIndex < 0) {
      Toast.show(isEnglishInterface() ? 'This saved record has no result index for locating its spectrum.' : '该保存记录没有可定位原始谱图的结果索引。', 'warning');
      return;
    }
    if (typeof EvidenceReview !== 'undefined' && currentTaskId) {
      await EvidenceReview.open(
        currentTaskId,
        'identification',
        '',
        Number.isInteger(resultIndex) ? resultIndex : null,
      );
    }
  }

  function openReview() {
    if (currentTaskId && typeof EvidenceReview !== 'undefined') EvidenceReview.open(currentTaskId, 'identification');
  }

  function openFollowup() {
    if (currentTaskId && typeof EvidenceReview !== 'undefined') EvidenceReview.open(currentTaskId, 'followup');
  }

  function openUnmatched() {
    if (currentTaskId && typeof EvidenceReview !== 'undefined') EvidenceReview.open(currentTaskId, 'unmatched');
  }

  function openReport() {
    if (currentTaskId && typeof SampleReport !== 'undefined') SampleReport.open(currentTaskId);
  }

  return {
    init, refreshSelect, loadForTask, openTask, setFilter, setSearch, setAnalysisFilter,
    setVisualView, setVisualSort, changeVisualPage,
    openEvidence, openReview, openFollowup, openUnmatched, openReport,
    applyReferenceStandardQC,
    updateDualPolarityReadiness, runDualPolarityShadow,
    get currentTaskId() { return currentTaskId; },
  };
})();
