/**
 * Integrated sample report.
 *
 * The report is a user-facing synthesis of the saved task. It never rescans,
 * rescores or changes candidate order. Full technical evidence remains
 * available through the server-authored standard and detailed reports.
 */
const SampleReport = (() => {
  let currentTaskId = null;
  let currentBundle = null;
  let currentPlans = { plans: [], metadata: {} };
  let currentSystem = {};
  let currentTic = null;
  let currentBlankExcludedRows = [];
  let currentExpected = null;
  let currentReferenceQC = null;
  const blankExcludedRowsCache = new Map();
  const esc = value => Tasks.escapeHtml(String(value == null ? '' : value));
  const isEnglishInterface = () => window.MSSLanguage?.getLanguage?.() === 'en';
  const integer = value => Number(value || 0).toLocaleString();
  const decimal = (value, digits = 2) => {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
  };
  const percentage = (value, total, digits = 1) => {
    const numerator = Number(value);
    const denominator = Number(total);
    return Number.isFinite(numerator) && Number.isFinite(denominator) && denominator > 0
      ? `${(numerator / denominator * 100).toFixed(digits)}%`
      : '—';
  };

  function englishBiosourceCopy(value) {
    const source = String(value || '').trim();
    if (!source) return '';
    let translated = window.MSSLanguage?.translateText?.(source) || source;
    const pathwayTerms = [
      ['莽草酸途径', 'Shikimate pathway'],
      ['莽草酸', 'shikimate'],
      ['醋酸-丙二酸途径', 'Acetate–malonate pathway'],
      ['乙酸-丙二酸途径', 'Acetate–malonate pathway'],
      ['氨基酸途径', 'Amino-acid pathway'],
      ['糖酵解相关途径', 'Glycolysis-related pathway'],
    ];
    for (const [sourceTerm, targetTerm] of pathwayTerms) {
      translated = translated.split(sourceTerm).join(targetTerm);
    }
    return /[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]/.test(translated)
      ? 'Biosource pathway context is not available in the current English catalog'
      : translated;
  }

  const zipEncoder = new TextEncoder();
  let zipCrcTable = null;

  function zipCrc32(bytes) {
    if (!zipCrcTable) {
      zipCrcTable = new Uint32Array(256);
      for (let index = 0; index < 256; index += 1) {
        let value = index;
        for (let bit = 0; bit < 8; bit += 1) {
          value = (value & 1) ? (0xedb88320 ^ (value >>> 1)) : (value >>> 1);
        }
        zipCrcTable[index] = value >>> 0;
      }
    }
    let crc = 0xffffffff;
    for (const byte of bytes) crc = zipCrcTable[(crc ^ byte) & 0xff] ^ (crc >>> 8);
    return (crc ^ 0xffffffff) >>> 0;
  }

  function zipDateTime(date = new Date()) {
    const year = Math.max(1980, date.getFullYear());
    return {
      date: ((year - 1980) << 9) | ((date.getMonth() + 1) << 5) | date.getDate(),
      time: (date.getHours() << 11) | (date.getMinutes() << 5) | Math.floor(date.getSeconds() / 2),
    };
  }

  function zipHeader(size, writer) {
    const bytes = new Uint8Array(size);
    writer(new DataView(bytes.buffer));
    return bytes;
  }

  async function zipPayload(bytes) {
    if (typeof CompressionStream !== 'function') return { bytes, method: 0 };
    try {
      const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream('deflate-raw'));
      const compressed = new Uint8Array(await new Response(stream).arrayBuffer());
      return compressed.byteLength < bytes.byteLength
        ? { bytes: compressed, method: 8 }
        : { bytes, method: 0 };
    } catch (_) {
      return { bytes, method: 0 };
    }
  }

  async function buildZip(entries, onProgress = null) {
    if (entries.length > 0xffff) throw new Error('批量报告数超过 ZIP 格式上限');
    const prepared = [];
    const timestamp = zipDateTime();
    for (let index = 0; index < entries.length; index += 1) {
      const entry = entries[index];
      const nameBytes = zipEncoder.encode(entry.name);
      const sourceBytes = typeof entry.content === 'string'
        ? zipEncoder.encode(entry.content)
        : entry.content;
      if (sourceBytes.byteLength > 0xffffffff) throw new Error(`文件过大，无法写入 ZIP：${entry.name}`);
      const compressed = await zipPayload(sourceBytes);
      prepared.push({
        nameBytes,
        sourceBytes,
        compressedBytes: compressed.bytes,
        method: compressed.method,
        crc: zipCrc32(sourceBytes),
      });
      onProgress?.({ phase: 'zip', completed: index + 1, total: entries.length, name: entry.name });
    }

    const localParts = [];
    const centralParts = [];
    let localOffset = 0;
    let centralSize = 0;
    prepared.forEach(entry => {
      const localHeader = zipHeader(30, view => {
        view.setUint32(0, 0x04034b50, true);
        view.setUint16(4, 20, true);
        view.setUint16(6, 0x0800, true);
        view.setUint16(8, entry.method, true);
        view.setUint16(10, timestamp.time, true);
        view.setUint16(12, timestamp.date, true);
        view.setUint32(14, entry.crc, true);
        view.setUint32(18, entry.compressedBytes.byteLength, true);
        view.setUint32(22, entry.sourceBytes.byteLength, true);
        view.setUint16(26, entry.nameBytes.byteLength, true);
        view.setUint16(28, 0, true);
      });
      localParts.push(localHeader, entry.nameBytes, entry.compressedBytes);

      const centralHeader = zipHeader(46, view => {
        view.setUint32(0, 0x02014b50, true);
        view.setUint16(4, 20, true);
        view.setUint16(6, 20, true);
        view.setUint16(8, 0x0800, true);
        view.setUint16(10, entry.method, true);
        view.setUint16(12, timestamp.time, true);
        view.setUint16(14, timestamp.date, true);
        view.setUint32(16, entry.crc, true);
        view.setUint32(20, entry.compressedBytes.byteLength, true);
        view.setUint32(24, entry.sourceBytes.byteLength, true);
        view.setUint16(28, entry.nameBytes.byteLength, true);
        view.setUint16(30, 0, true);
        view.setUint16(32, 0, true);
        view.setUint16(34, 0, true);
        view.setUint16(36, 0, true);
        view.setUint32(38, 0, true);
        view.setUint32(42, localOffset, true);
      });
      centralParts.push(centralHeader, entry.nameBytes);
      centralSize += centralHeader.byteLength + entry.nameBytes.byteLength;
      localOffset += localHeader.byteLength + entry.nameBytes.byteLength + entry.compressedBytes.byteLength;
    });

    const endRecord = zipHeader(22, view => {
      view.setUint32(0, 0x06054b50, true);
      view.setUint16(4, 0, true);
      view.setUint16(6, 0, true);
      view.setUint16(8, prepared.length, true);
      view.setUint16(10, prepared.length, true);
      view.setUint32(12, centralSize, true);
      view.setUint32(16, localOffset, true);
      view.setUint16(20, 0, true);
    });
    return new Blob([...localParts, ...centralParts, endRecord], { type: 'application/zip' });
  }

  function compactNumber(value) {
    const parsed = Number(value);
    if (!Number.isFinite(parsed)) return '—';
    if (Math.abs(parsed) >= 1e9) return `${(parsed / 1e9).toFixed(1)}G`;
    if (Math.abs(parsed) >= 1e6) return `${(parsed / 1e6).toFixed(1)}M`;
    if (Math.abs(parsed) >= 1e3) return `${(parsed / 1e3).toFixed(1)}k`;
    return parsed.toFixed(0);
  }

  function retentionTimeMinutes(row) {
    const parsed = Number(row?.retention_time);
    if (!Number.isFinite(parsed)) return null;
    const unit = String(row?.retention_time_unit || 'minutes').trim().toLowerCase();
    if (unit.includes('millisecond') || ['ms', 'msec', 'msecs'].includes(unit)) return parsed / 60000;
    if (unit.includes('second') || ['s', 'sec', 'secs'].includes(unit)) return parsed / 60;
    if (unit.includes('hour') || ['h', 'hr', 'hrs'].includes(unit)) return parsed * 60;
    return parsed;
  }

  function retentionTimeClock(minutes) {
    if (!Number.isFinite(minutes) || minutes < 0) return '—';
    const secondsTotal = minutes * 60;
    const minutesWhole = Math.floor(secondsTotal / 60);
    const seconds = secondsTotal - minutesWhole * 60;
    return `${String(minutesWhole).padStart(2, '0')}:${seconds.toFixed(1).padStart(4, '0')}`;
  }

  function rawSpectrumLocator(row) {
    const rt = retentionTimeMinutes(row);
    const selectedIonMz = Number(row?.precursor_mz);
    const searchPrecursorMz = Number(row?.search_precursor_mz);
    const searchSource = String(row?.search_precursor_mz_source || '');
    const featureUsed = row?.lcms_feature_used_for_search === true
      || searchSource === 'measured_continuous_ms1_feature_centroid';
    return {
      rt,
      clock: retentionTimeClock(rt),
      precursorMz: selectedIonMz,
      selectedIonMz,
      searchPrecursorMz: Number.isFinite(searchPrecursorMz)
        ? searchPrecursorMz : selectedIonMz,
      searchSource,
      featureUsed,
      searchDeltaPpm: Number(row?.search_precursor_delta_ppm),
      spectrumId: String(row?.spectrum_id || '').trim(),
      precursorScanRef: String(row?.precursor_scan_ref || '').trim(),
    };
  }

  function rawSpectrumLocatorCell(row) {
    const locator = rawSpectrumLocator(row);
    const rtText = locator.rt == null ? 'RT 未记录' : `RT ${locator.rt.toFixed(3)} min`;
    const mzText = Number.isFinite(locator.selectedIonMz)
      ? `${locator.featureUsed ? '隔离目标' : '前体'} m/z ${locator.selectedIonMz.toFixed(5)}`
      : 'm/z 未记录';
    const searchText = locator.featureUsed && Number.isFinite(locator.searchPrecursorMz)
      ? `<small>检索用实测 MS1 m/z ${esc(locator.searchPrecursorMz.toFixed(5))}${Number.isFinite(locator.searchDeltaPpm) ? ` · 与隔离目标 Δ ${esc(locator.searchDeltaPpm.toFixed(2))} ppm` : ''}</small>`
      : '';
    const spectrumText = locator.spectrumId ? `MS2 ${esc(locator.spectrumId)}` : 'MS2 编号未记录';
    return `<div class="report-raw-locator"><b>${esc(rtText)}</b><span>${esc(locator.clock)}</span><small>${esc(mzText)} · ${spectrumText}</small>${searchText}${locator.precursorScanRef ? `<small>关联前级扫描 ${esc(locator.precursorScanRef)}</small>` : ''}</div>`;
  }

  function nearestIntensity(points, retentionTime) {
    if (!points.length || !Number.isFinite(retentionTime)) return 0;
    let low = 0;
    let high = points.length - 1;
    while (low < high) {
      const middle = Math.floor((low + high) / 2);
      if (points[middle].rt < retentionTime) low = middle + 1;
      else high = middle;
    }
    const right = points[low];
    const left = points[Math.max(0, low - 1)];
    return Math.abs(Number(left.rt) - retentionTime) <= Math.abs(Number(right.rt) - retentionTime)
      ? Number(left.tic || 0)
      : Number(right.tic || 0);
  }

  function withinPpm(observedMz, referenceMz, tolerancePpm) {
    const observed = Number(observedMz);
    const reference = Number(referenceMz);
    const tolerance = Number(tolerancePpm);
    if (!Number.isFinite(observed) || !Number.isFinite(reference)
        || observed <= 0 || reference <= 0 || !Number.isFinite(tolerance)) {
      return false;
    }
    return Math.abs(observed - reference) / reference * 1e6 <= tolerance;
  }

  function downsampleChartPoints(points, maximum = 900) {
    if (points.length <= maximum) return points;
    const stride = Math.ceil(points.length / maximum);
    const reduced = [];
    for (let start = 0; start < points.length; start += stride) {
      const bucket = points.slice(start, Math.min(points.length, start + stride));
      const peak = bucket.reduce((best, point) => Number(point.tic || 0) > Number(best.tic || 0) ? point : best, bucket[0]);
      reduced.push(peak);
    }
    return reduced.sort((a, b) => Number(a.rt) - Number(b.rt));
  }

  function prominentChromatogramPeaks(points, maximumIntensity, limit = 12) {
    if (points.length < 7 || !Number.isFinite(maximumIntensity) || maximumIntensity <= 0) return [];
    const radius = Math.max(2, Math.min(5, Math.round(points.length / 250)));
    const intensityFloor = maximumIntensity * 0.04;
    const candidates = [];
    for (let index = radius; index < points.length - radius; index += 1) {
      const intensity = Number(points[index].tic || 0);
      if (intensity < intensityFloor) continue;
      let localMaximum = true;
      for (let offset = -radius; offset <= radius; offset += 1) {
        if (offset && Number(points[index + offset].tic || 0) > intensity) {
          localMaximum = false;
          break;
        }
      }
      if (localMaximum) candidates.push(points[index]);
    }
    const rtSpan = Number(points[points.length - 1].rt) - Number(points[0].rt);
    const minimumSeparation = Math.max(0.04, Math.min(0.32, rtSpan / 55));
    const selected = [];
    [...candidates].sort((a, b) => Number(b.tic || 0) - Number(a.tic || 0)).forEach(point => {
      if (selected.length >= limit) return;
      if (selected.every(existing => Math.abs(Number(existing.rt) - Number(point.rt)) >= minimumSeparation)) {
        selected.push(point);
      }
    });
    return selected.sort((a, b) => Number(a.rt) - Number(b.rt));
  }

  function matchedRowsFromGroups(groups) {
    return (groups || []).flatMap(group => (
      Array.isArray(group?.rows) && group.rows.length
        ? group.rows
        : (group?.rep ? [group.rep] : [])
    ));
  }

  function rowDisplayName(row) {
    const match = row?.match || {};
    return match.compound_name || 'Unnamed candidate';
  }

  function candidateEventMarkers(groups, source, minRt, maxRt) {
    const markerLabels = {
      confirmed: '已确证存在',
      high: 'High 单一二维连接检出',
      isomer: 'High 异构体组检出',
      medium: '中等证据候选',
    };
    const priority = { medium: 1, isomer: 2, high: 3, confirmed: 4 };
    const raw = (groups || []).flatMap(group => {
      const rows = Array.isArray(group?.rows) && group.rows.length
        ? group.rows
        : (group?.rep ? [group.rep] : []);
      const status = group?.status || SampleResultData.evidenceStatus(rows);
      if (!markerLabels[status?.key]) return [];
      return rows.map(row => {
        const locator = rawSpectrumLocator(row);
        if (!Number.isFinite(locator.rt) || locator.rt < minRt || locator.rt > maxRt) return null;
        const match = row.match || {};
        return {
          key: status.key,
          label: status.label || markerLabels[status.key],
          name: rowDisplayName(row),
          identityKey: SampleResultData.identityKey(match),
          rt: locator.rt,
          precursorMz: locator.precursorMz,
          spectrumId: locator.spectrumId,
          rankScore: Number(match.flagship_entity_score ?? match.rank_score ?? match.final_score ?? 0),
        };
      }).filter(Boolean);
    }).sort((a, b) => a.rt - b.rt || priority[b.key] - priority[a.key] || b.rankScore - a.rankScore);

    // Several DDA scans can occur within the same visual pixel. Merge only
    // genuinely local events; never merge the same identity across distant RTs.
    const clusterToleranceMinutes = 0.025;
    const clusters = [];
    raw.forEach(marker => {
      const cluster = clusters[clusters.length - 1];
      if (!cluster || marker.rt - cluster.lastRt > clusterToleranceMinutes) {
        clusters.push({
          best: marker,
          lastRt: marker.rt,
          spectra: 1,
          identities: new Set([marker.identityKey]),
          names: new Set([marker.name]),
        });
        return;
      }
      cluster.lastRt = marker.rt;
      cluster.spectra += 1;
      cluster.identities.add(marker.identityKey);
      cluster.names.add(marker.name);
      if (priority[marker.key] > priority[cluster.best.key]
          || (priority[marker.key] === priority[cluster.best.key] && marker.rankScore > cluster.best.rankScore)) {
        cluster.best = marker;
      }
    });
    return clusters.map(cluster => {
      const names = [...cluster.names];
      const primaryName = String(names[0] || '未命名候选');
      const compactPrimaryName = primaryName.length > 96
        ? `${primaryName.slice(0, 93)}…` : primaryName;
      return {
        ...cluster.best,
        name: names.length > 1
          ? (isEnglishInterface()
            ? `${compactPrimaryName}; ${names.length - 1} other identification ${names.length - 1 === 1 ? 'group' : 'groups'} in the same RT cluster`
            : `${compactPrimaryName}；同一 RT 簇另有 ${names.length - 1} 个鉴定组`)
          : compactPrimaryName,
        candidateCount: cluster.identities.size,
        spectra: cluster.spectra,
        intensity: nearestIntensity(source, cluster.best.rt),
      };
    }).sort((a, b) => priority[a.key] - priority[b.key] || a.rt - b.rt);
  }

  function chromatogramChart(groups, metric = 'tic', ticOverview = currentTic, compact = false) {
    const isBpc = metric === 'bpc';
    const valueKey = isBpc ? 'bpc' : 'tic';
    const chartName = isEnglishInterface()
      ? (isBpc ? 'Base-peak chromatogram (BPC)' : 'Total ion chromatogram (TIC)')
      : (isBpc ? '基峰色谱图（BPC）' : '总离子流色谱图（TIC）');
    const shortName = isBpc ? 'BPC' : 'TIC';
    const ms1 = Array.isArray(ticOverview?.ms1) ? ticOverview.ms1 : [];
    const ms2 = Array.isArray(ticOverview?.ms2) ? ticOverview.ms2 : [];
    const ticFromSourceChromatogram = ticOverview?.tic_source === 'mzml_tic_chromatogram_ms_level_array';
    const validMs1 = ms1
      .filter(point => point?.[valueKey] !== null && point?.[valueKey] !== undefined && point?.[valueKey] !== '')
      .map(point => ({
        rt: Number(point.rt),
        tic: Number(point[valueKey]),
        basePeakMz: Number(point.base_peak_mz),
        scanId: point.scan_id,
      }))
      .filter(point => Number.isFinite(point.rt) && Number.isFinite(point.tic) && point.tic >= 0)
      .sort((a, b) => a.rt - b.rt);
    const validMs2 = ms2.map(point => ({ ...point, rt: Number(point.rt), tic: Number(point.tic || 0) }))
      .filter(point => Number.isFinite(point.rt) && Number.isFinite(point.tic) && point.tic >= 0)
      .sort((a, b) => a.rt - b.rt);
    const source = validMs1.length >= 2 ? validMs1 : (isBpc ? [] : validMs2);
    if (source.length < 2) {
      const emptyMessage = isBpc
        ? (ticFromSourceChromatogram
          ? '源 mzML 只保存了 MS1 TIC 总强度，未保存 MS1 峰表或 BPC chromatogram；TIC 不能反推出真实 BPC。'
          : '当前任务没有可用的逐扫描 MS1 Base Peak Intensity；系统不会用 MS2 基峰冒充真实 BPC。')
        : '已保存结果没有可绘制的 RT/TIC 强度序列。';
      return `<div class="report-chart-empty">${emptyMessage}</div>`;
    }

    const points = downsampleChartPoints(source);
    const width = 1080;
    const height = 360;
    const left = 72;
    const right = 24;
    const top = 26;
    const bottom = 52;
    const plotWidth = width - left - right;
    const plotHeight = height - top - bottom;
    const minRt = Number(source[0].rt);
    const maxRt = Number(source[source.length - 1].rt);
    const maxTic = source.reduce((maximum, point) => Math.max(maximum, Number(point.tic || 0)), 1);
    const xFor = rt => left + ((Number(rt) - minRt) / Math.max(maxRt - minRt, 1e-9)) * plotWidth;
    const yFor = tic => top + (1 - Math.max(0, Number(tic || 0)) / maxTic) * plotHeight;
    const linePoints = points.map(point => `${xFor(point.rt).toFixed(2)},${yFor(point.tic).toFixed(2)}`).join(' ');
    const areaPath = `M ${xFor(points[0].rt).toFixed(2)} ${top + plotHeight} L ${linePoints.replaceAll(' ', ' L ')} L ${xFor(points[points.length - 1].rt).toFixed(2)} ${top + plotHeight} Z`;
    const xTicks = Array.from({ length: 7 }, (_, index) => {
      const fraction = index / 6;
      const value = minRt + (maxRt - minRt) * fraction;
      const x = left + plotWidth * fraction;
      return `<line x1="${x}" y1="${top}" x2="${x}" y2="${top + plotHeight}" class="tic-grid"/><text x="${x}" y="${height - 18}" text-anchor="middle" class="tic-axis-label">${value.toFixed(1)}</text>`;
    }).join('');
    const yTicks = Array.from({ length: 5 }, (_, index) => {
      const fraction = index / 4;
      const value = maxTic * (1 - fraction);
      const y = top + plotHeight * fraction;
      return `<line x1="${left}" y1="${y}" x2="${left + plotWidth}" y2="${y}" class="tic-grid"/><text x="${left - 12}" y="${y + 4}" text-anchor="end" class="tic-axis-label">${esc(compactNumber(value))}</text>`;
    }).join('');
    const scanStride = Math.max(1, Math.ceil(validMs2.length / 220));
    const scanRug = validMs2.filter((_, index) => index % scanStride === 0).map(point => {
      const x = xFor(point.rt).toFixed(2);
      return `<line x1="${x}" y1="${top + plotHeight + 3}" x2="${x}" y2="${top + plotHeight + 10}" class="tic-ms2-rug"><title>MS2 scan ${esc(point.scan_id ?? '—')} · RT ${point.rt.toFixed(3)} min</title></line>`;
    }).join('');
    const candidateMarkers = candidateEventMarkers(groups, source, minRt, maxRt);
    const markerCounts = candidateMarkers.reduce((counts, marker) => {
      counts[marker.key] += 1;
      return counts;
    }, { confirmed: 0, high: 0, isomer: 0, medium: 0 });
    const highConfidenceLocations = markerCounts.confirmed + markerCounts.high;
    const markers = candidateMarkers.map(marker => {
      const x = xFor(marker.rt);
      const y = yFor(marker.intensity);
      const title = `${marker.label}：${marker.name} · RT ${marker.rt.toFixed(3)} min · m/z ${Number.isFinite(marker.precursorMz) ? marker.precursorMz.toFixed(5) : '—'} · ${marker.spectra} 张局部 MS2`;
      if (marker.key === 'isomer') {
        const size = 7;
        return `<path d="M ${x.toFixed(2)} ${(y - size).toFixed(2)} L ${(x + size).toFixed(2)} ${y.toFixed(2)} L ${x.toFixed(2)} ${(y + size).toFixed(2)} L ${(x - size).toFixed(2)} ${y.toFixed(2)} Z" class="tic-isomer-marker"><title>${esc(title)}</title></path>`;
      }
      if (marker.key === 'medium') {
        return `<rect x="${(x - 4.5).toFixed(2)}" y="${(y - 4.5).toFixed(2)}" width="9" height="9" rx="1.5" class="tic-medium-marker"><title>${esc(title)}</title></rect>`;
      }
      return `<circle cx="${x.toFixed(2)}" cy="${y.toFixed(2)}" r="5.8" class="tic-high-marker"><title>${esc(title)}</title></circle>`;
    }).join('');
    const coverageTolerance = Math.max(0.08, Math.min(0.18, (maxRt - minRt) / 150));
    const bpcMzTolerancePpm = Math.max(1, Math.min(
      100,
      Number(currentBundle?.task?.params?.ms2_feature_link_tolerance_ppm || 20)
    ));
    const majorPeaks = prominentChromatogramPeaks(source, maxTic);
    const markerCoversPeak = (peak, marker) => {
      if (Math.abs(marker.rt - Number(peak.rt)) > coverageTolerance) return false;
      if (!isBpc || !Number.isFinite(Number(peak.basePeakMz))) return true;
      return withinPpm(peak.basePeakMz, marker.precursorMz, bpcMzTolerancePpm);
    };
    const uncoveredPeaks = majorPeaks.filter(
      peak => !candidateMarkers.some(marker => markerCoversPeak(peak, marker))
    );
    const allMatchedRows = matchedRowsFromGroups(groups);
    const uncoveredDiagnostics = uncoveredPeaks.map(peak => {
      const peakRt = Number(peak.rt);
      const peakMz = Number(peak.basePeakMz);
      const nearbyMs2Rows = validMs2.filter(
        point => Math.abs(point.rt - peakRt) <= coverageTolerance
      );
      const ms2Nearby = nearbyMs2Rows.length;
      const precursorMatchedMs2 = isBpc && Number.isFinite(peakMz)
        ? nearbyMs2Rows.filter(point => withinPpm(
          peakMz, point.precursor_mz, bpcMzTolerancePpm
        )).length
        : ms2Nearby;
      const blankExcludedTargeted = isBpc && Number.isFinite(peakMz)
        ? currentBlankExcludedRows.filter(row => {
          const rowRt = retentionTimeMinutes(row);
          return Number.isFinite(rowRt)
            && Math.abs(rowRt - peakRt) <= coverageTolerance
            && withinPpm(peakMz, Number(row.precursor_mz), bpcMzTolerancePpm)
            && row.blank_excluded === true;
        })
        : [];
      const matchedNearby = allMatchedRows.filter(row => {
        const locator = rawSpectrumLocator(row);
        if (!Number.isFinite(locator.rt)
            || Math.abs(locator.rt - peakRt) > coverageTolerance) return false;
        return !isBpc || !Number.isFinite(peakMz)
          || withinPpm(peakMz, locator.precursorMz, bpcMzTolerancePpm);
      }).length;
      if (!ms2Nearby) {
        return { reason: 'no_ms2', detail: '附近未采集 MS2；高 MS1 响应本身不能提供精确结构身份' };
      }
      if (isBpc && Number.isFinite(peakMz) && !precursorMatchedMs2) {
        return {
          reason: 'bpc_precursor_not_selected',
          detail: `同一 RT 窗口采集 ${ms2Nearby} 张 MS2，但没有一张隔离目标与该 BPC 基峰 m/z ${peakMz.toFixed(5)} 在 ${bpcMzTolerancePpm.toFixed(0)} ppm 内一致`,
        };
      }
      if (blankExcludedTargeted.length
          && blankExcludedTargeted.length === precursorMatchedMs2) {
        const ratios = blankExcludedTargeted
          .map(row => Number(row.blank_sample_ratio))
          .filter(Number.isFinite)
          .sort((a, b) => a - b);
        const ratioLabel = ratios.length
          ? (Math.abs(ratios[0] - ratios[ratios.length - 1]) < 1e-9
            ? ratios[0].toFixed(2)
            : `${ratios[0].toFixed(2)}–${ratios[ratios.length - 1].toFixed(2)}`)
          : '未记录';
        const ratioThreshold = Number(
          currentBundle?.task?.params?.blank_sample_ratio || 3
        );
        return {
          reason: 'blank_excluded',
          detail: `该峰离子的 ${precursorMatchedMs2} 张 MS2 均在谱库检索前被空白硬扣除（样本/空白比 ${ratioLabel}；排除阈值 ${ratioThreshold.toFixed(2)}）`,
        };
      }
      if (!matchedNearby) {
        const blankNote = blankExcludedTargeted.length
          ? `；另有 ${blankExcludedTargeted.length} 张在检索前被空白硬扣除`
          : '';
        return {
          reason: 'no_library_candidate',
          detail: isBpc && Number.isFinite(peakMz)
            ? `该峰离子附近采集 ${precursorMatchedMs2} 张 MS2，但没有候选通过当前谱库与分数门槛${blankNote}`
            : `附近采集 ${ms2Nearby} 张 MS2，但没有候选通过当前谱库与分数门槛`,
        };
      }
      return {
        reason: 'below_priority_threshold',
        detail: isBpc && Number.isFinite(peakMz)
          ? `该峰离子附近采集 ${precursorMatchedMs2} 张 MS2，其中 ${matchedNearby} 张获得候选，但最高仅为 Low/Unknown 线索`
          : `附近采集 ${ms2Nearby} 张 MS2，其中 ${matchedNearby} 张获得候选，但最高仅为 Low/Unknown 线索`,
      };
    });
    const uncoveredReasonCounts = uncoveredDiagnostics.reduce((counts, item) => {
      counts[item.reason] = (counts[item.reason] || 0) + 1;
      return counts;
    }, {});
    const uncoveredMarkers = uncoveredPeaks.map((peak, index) => {
      const x = xFor(peak.rt).toFixed(2);
      const y = yFor(peak.tic).toFixed(2);
      const mzLabel = isBpc && Number.isFinite(Number(peak.basePeakMz))
        ? ` · base-peak m/z ${Number(peak.basePeakMz).toFixed(5)}` : '';
      return `<g class="tic-unassigned-group"><circle cx="${x}" cy="${y}" r="7" class="tic-unassigned-marker"/><text x="${x}" y="${(Number(y) + 3.2).toFixed(2)}" text-anchor="middle" class="tic-unassigned-label">?</text><title>主要 ${shortName} 峰 ${index + 1} · RT ${Number(peak.rt).toFixed(3)} min${mzLabel} · ${uncoveredDiagnostics[index].detail}；不代表该峰没有化学信号</title></g>`;
    }).join('');
    const coveredMajorPeaks = majorPeaks.length - uncoveredPeaks.length;
    const english = isEnglishInterface();
    const uncoveredPeakLabel = english
      ? (isBpc
        ? 'major peaks without MS2 evidence for the base-peak ion'
        : 'major peaks without a nearby primary result')
      : (isBpc
        ? '未建立该基峰离子 MS2 结果的主要峰'
        : '暂无重点结果邻近的主要峰');
    const sourceLabel = english
      ? (isBpc
        ? (ticOverview?.bpc_source === 'mzml_basepeak_chromatogram'
          ? 'Source mzML BPC' : 'Measured MS1 BPC')
        : (validMs1.length >= 2
          ? (ticFromSourceChromatogram ? 'Source mzML MS1 TIC' : 'Measured MS1 TIC')
          : 'Approximate MS2 TIC'))
      : (isBpc
        ? (ticOverview?.bpc_source === 'mzml_basepeak_chromatogram' ? '源 mzML BPC' : '实测 MS1 BPC')
        : (validMs1.length >= 2
          ? (ticFromSourceChromatogram ? '源 mzML MS1 TIC' : '实测 MS1 TIC')
          : 'MS2 TIC 近似'));
    const gradientId = compact
      ? (isBpc ? 'sampleSummaryBpcFill' : 'sampleSummaryTicFill')
      : (isBpc ? 'sampleBpcFill' : 'sampleTicFill');
    const traceColor = isBpc ? '#6f5aa8' : '#1b7f79';
    const traceClass = isBpc ? 'bpc-main-line' : 'tic-main-line';
    const legendClass = isBpc ? 'legend-bpc' : 'legend-tic';
    const axisTitle = isBpc ? 'Base peak intensity' : 'TIC intensity';
    const chartAriaLabel = english
      ? `${chartName}, retention time ${minRt.toFixed(2)}–${maxRt.toFixed(2)} min; ${highConfidenceLocations} confirmed or high-confidence locations, ${markerCounts.isomer} isomer-group locations, ${markerCounts.medium} moderate-evidence candidate locations and ${uncoveredPeaks.length} ${uncoveredPeakLabel}`
      : `${chartName}，保留时间 ${minRt.toFixed(2)} 至 ${maxRt.toFixed(2)} 分钟；标出 ${highConfidenceLocations} 个已确证或高可信检测位置、${markerCounts.isomer} 个异构体组检测位置、${markerCounts.medium} 个中等证据候选检测位置和 ${uncoveredPeaks.length} 个${uncoveredPeakLabel}`;

    return `<div class="report-tic-frame">
      <svg class="report-tic-svg report-${metric}-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${chartAriaLabel}">
        <defs><linearGradient id="${gradientId}" x1="0" y1="0" x2="0" y2="1"><stop offset="0%" stop-color="${traceColor}" stop-opacity=".34"/><stop offset="100%" stop-color="${traceColor}" stop-opacity=".03"/></linearGradient></defs>
        ${xTicks}${yTicks}
        <path d="${areaPath}" fill="url(#${gradientId})"/>
        <polyline points="${linePoints}" class="${traceClass}"/>
        ${scanRug}${uncoveredMarkers}${markers}
        <line x1="${left}" y1="${top + plotHeight}" x2="${left + plotWidth}" y2="${top + plotHeight}" class="tic-axis"/>
        <text x="${left + plotWidth / 2}" y="${height - 1}" text-anchor="middle" class="tic-axis-title">Retention time (min)</text>
        <text x="17" y="${top + plotHeight / 2}" text-anchor="middle" transform="rotate(-90 17 ${top + plotHeight / 2})" class="tic-axis-title">${axisTitle}</text>
      </svg>
      ${compact
        ? `<div class="sample-chromatogram-meta"><span>${sourceLabel} · ${integer(source.length)} scans</span><span>${english ? `${integer(highConfidenceLocations + markerCounts.isomer + markerCounts.medium)} primary-result acquisition locations · ${integer(uncoveredPeaks.length)} ${uncoveredPeakLabel}` : `主要结果检测位置 ${integer(highConfidenceLocations + markerCounts.isomer + markerCounts.medium)} · ${uncoveredPeakLabel} ${integer(uncoveredPeaks.length)}`}</span></div>`
        : `<div class="report-chart-legend"><span><i class="${legendClass}"></i>${sourceLabel} · ${integer(source.length)} scans</span><span><i class="legend-ms2"></i>MS2 扫描位置 · ${integer(validMs2.length)}</span><span><i class="legend-high"></i>已确证／高可信检测位置 · ${integer(highConfidenceLocations)}</span><span><i class="legend-isomer"></i>异构体组检测位置 · ${integer(markerCounts.isomer)}</span><span><i class="legend-medium"></i>中等证据候选检测位置 · ${integer(markerCounts.medium)}</span><span><i class="legend-unassigned"></i>${uncoveredPeakLabel} · ${integer(uncoveredPeaks.length)}</span></div>
          <div class="report-tic-coverage"><b>${shortName} 主要峰说明</b><span>自动识别 ${integer(majorPeaks.length)} 个主要 ${shortName} 峰，其中 ${integer(coveredMajorPeaks)} 个${isBpc ? `同时满足 RT ±${coverageTolerance.toFixed(2)} min 与隔离目标 m/z ${bpcMzTolerancePpm.toFixed(0)} ppm ` : `在 ±${coverageTolerance.toFixed(2)} min 内 `}关联重点候选的 MS2 采集事件，${integer(uncoveredPeaks.length)} 个以灰色问号提示；灰色问号不应解释为“没有化合物”。</span><small>彩色符号按每个局部 High、异构体组 High 或 Medium MS2 事件投影，并仅合并 0.025 min 内重叠符号；不会再用化合物的全样本唯一代表 RT 隐藏其他事件。问号原因：附近未采 MS2 ${integer(uncoveredReasonCounts.no_ms2 || 0)} 个${isBpc ? `、同一 RT 有 MS2 但该基峰离子未被选作隔离目标 ${integer(uncoveredReasonCounts.bpc_precursor_not_selected || 0)} 个、因空白硬扣除未进入谱库检索 ${integer(uncoveredReasonCounts.blank_excluded || 0)} 个` : ''}、已采该峰离子但无候选过门槛 ${integer(uncoveredReasonCounts.no_library_candidate || 0)} 个、仅 Low/Unknown ${integer(uncoveredReasonCounts.below_priority_threshold || 0)} 个。</small></div>`}
    </div>`;
  }

  function summaryChromatograms(groups, ticOverview) {
    const ms1TraceCount = Number(ticOverview?.ms1_tic_point_count
      ?? (Array.isArray(ticOverview?.ms1) ? ticOverview.ms1.length : 0));
    const ms2Count = Array.isArray(ticOverview?.ms2) ? ticOverview.ms2.length : 0;
    const bpcAvailable = ticOverview?.bpc_available === true;
    const ticFromSourceChromatogram = ticOverview?.tic_source === 'mzml_tic_chromatogram_ms_level_array';
    const hasMs1Trace = ms1TraceCount >= 2;
    const ticLabel = ms1TraceCount >= 2
      ? (ticFromSourceChromatogram
        ? `${integer(ms1TraceCount)} 个源 MS1 TIC 采集点`
        : `${integer(ms1TraceCount)} 张连续 MS1`)
      : `${integer(ms2Count)} 张 MS2 近似`;
    const bpcLabel = bpcAvailable
      ? (ticOverview?.bpc_source === 'mzml_basepeak_chromatogram'
        ? `${integer(ms1TraceCount)} 个源 BPC 采集点`
        : `${integer(ms1TraceCount)} 张连续 MS1`)
      : (ticFromSourceChromatogram
        ? '源文件仅有 MS1 TIC，无法反推 BPC'
        : '当前数据未提供连续 MS1 基峰强度');
    const scopeText = hasMs1Trace
      ? '读取源 mzML 的 MS1 色谱记录并标出重点候选对应的 MS2 采集位置；DAD/UV 辅助光学谱不进入这两张图。'
      : '当前保存结果为 MS2-only 采集事件；TIC 仅按 MS2 谱图强度近似展示，不属于 MS1 色谱，BPC 不作推算。';
    return `<section class="sample-result-chromatograms" aria-label="样本 TIC 与 BPC 摘要">
      <header class="sample-chromatogram-header">
        <div><span>CHROMATOGRAM OVERVIEW</span><h3>TIC / BPC 色谱摘要</h3><p>${scopeText}</p></div>
        <button class="text-action" onclick="SampleOverview.openReport()">查看完整色谱与峰说明 →</button>
      </header>
      <div class="sample-chromatogram-grid">
        <article class="sample-chromatogram-card chromatogram-tic"><div class="sample-chromatogram-title"><div><b>总离子流色谱图</b><span>TIC · ${ticLabel}</span></div><small>Retention time</small></div>${chromatogramChart(groups || [], 'tic', ticOverview, true)}</article>
        <article class="sample-chromatogram-card chromatogram-bpc"><div class="sample-chromatogram-title"><div><b>基峰色谱图</b><span>BPC · ${bpcLabel}</span></div><small>Retention time</small></div>${chromatogramChart(groups || [], 'bpc', ticOverview, true)}</article>
      </div>
    </section>`;
  }

  function init() {
    document.getElementById('reportTaskSelect')?.addEventListener('change', event => load(event.target.value));
    document.addEventListener('mssystem:languagechange', () => {
      if (!currentTaskId) return;
      SampleResultData.clear(currentTaskId);
      if (document.getElementById('page-report')?.classList.contains('active')) {
        load(currentTaskId);
      }
    });
  }

  async function refreshSelect(preferredTaskId = null) {
    const select = document.getElementById('reportTaskSelect');
    if (!select) return;
    try {
      const data = await API.get('/tasks?status=completed');
      const previous = preferredTaskId || select.value;
      const tasks = data.tasks || [];
      const latestSelection = select.value;
      select.innerHTML = '<option value="">-- 选择已完成的任务 --</option>' +
        tasks.map(task =>
          `<option value="${esc(task.id)}">${esc(task.filename)} · ${integer(task.matched_spectra)}/${integer(task.total_spectra)} MS2 有候选</option>`
        ).join('');
      const target = preferredTaskId || latestSelection || previous || currentTaskId;
      if (target && tasks.some(task => task.id === target)) select.value = target;
    } catch (error) {
      console.error('Failed to refresh report selector:', error);
    }
  }

  function counts(bundle) {
    const task = bundle.task || {};
    const funnel = task.identification_funnel || {};
    const metadata = bundle.compounds?.metadata || {};
    const eventAggregation = metadata.precursor_event_aggregation || {};
    const derived = key => bundle.groups.filter(group => group.status.key === key).length;
    const namedGroups = bundle.groups.filter(group => ![
      'unannotated', 'unresolved_event',
    ].includes(group.status.key));
    const unresolved = derived('unresolved_event');
    const unresolvedIncomplete = Number(
      eventAggregation.incomplete_formula_unresolved_event_count
      ?? bundle.groups.filter(group => (
        group.status.key === 'unresolved_event'
        && unresolvedEventKind(group) === 'incomplete_formula'
      )).length
    );
    return {
      compounds: namedGroups.length,
      main: bundle.groups.filter(group => group.status.countedAsSampleCompound).length,
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
      unresolved,
      unresolvedFormulaConflict: Number(
        eventAggregation.cross_formula_unresolved_event_count
        ?? Math.max(0, unresolved - unresolvedIncomplete)
      ),
      unresolvedIncomplete,
      insufficient: derived('low'),
      background: derived('background'),
      legacyHigh: Number(metadata.legacy_high_supported_identity_count ?? 0),
      matchedSpectra: Number(task.matched_spectra ?? funnel.matched_ms2_scans ?? bundle.matched.length),
      totalSpectra: Number(task.total_spectra ?? funnel.ms2_scans ?? 0),
      unmatchedSpectra: Number(funnel.unmatched_ms2_scans ?? Math.max(0, Number(task.total_spectra || 0) - Number(task.matched_spectra || 0))),
      matchedEvents: Number(funnel.matched_ms2_precursor_events ?? task.matched_precursor_features ?? 0),
    };
  }

  function displayName(group) {
    const match = group.rep.match || {};
    return match.compound_display_name || group.compound?.display_name
      || match.compound_name || 'Unnamed candidate';
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

  function displayClassification(match) {
    const projected = String(match?.display_class || '').trim();
    const raw = String(match?.classification || '').trim();
    const projectedUsable = projected && !['unknown', 'other', 'unclassified'].includes(projected.toLowerCase());
    const rawUsable = raw && !['unknown', 'other', 'unclassified'].includes(raw.toLowerCase());
    const label = projectedUsable ? projected : (rawUsable ? raw : '');
    if (!label) {
      const identityHigh = [
        match?.confidence,
        match?.exact_structure_confidence,
        match?.two_dimensional_connectivity_confidence,
      ].some(value => String(value || '').trim().toLowerCase() === 'high');
      return identityHigh
        ? '尚无可靠类别归属（谱库候选证据为 High）'
        : '尚无可靠类别归属';
    }
    const source = String(
      (projectedUsable ? match?.display_classification_source : '')
      || match?.classification_source || ''
    ).trim().toLowerCase();
    const trusted = ['curated', 'manual', 'smarts', 'structure', 'rule', 'database_curated']
      .some(prefix => source === prefix || source.startsWith(`${prefix}_`))
      || source.startsWith('runtime_structure_rule');
    const untrusted = !source || ['unknown', 'other', 'existing', 'ml', 'msclass', 'msclass_spectrum', 'fragment', 'fragment_inference', 'spectrum'].includes(source)
      || ['ml_', 'model_', 'msclass_', 'fragment_'].some(prefix => source.startsWith(prefix));
    return !untrusted && (trusted || match?.classification_reliable === true)
      ? label : `待核查类别线索：${label}`;
  }

  function claimConclusion(group) {
    const conclusion = group?.conclusion || group?.compound?.conclusion || {};
    const key = String(conclusion.key || 'insufficient');
    if (key === 'confirmed') return { label: '已确证存在', css: 'eligible', note: String(conclusion.note || '已通过当前正式结论条件') };
    if (key === 'high') return { label: 'High 单一二维连接检出', css: 'eligible', note: String(conclusion.note || '单一二维连接结构具有 High 级谱库支持；完整结构身份按证据边界另行报告') };
    if (key === 'isomer') {
      const isHigh = !!conclusion.high_confidence_group_detected;
      return {
        label: isHigh ? 'High 异构体组检出' : '中等证据异构体组',
        css: 'candidate',
        note: String(conclusion.note || (isHigh
          ? 'High 级证据只支持异构体组，不能释放单一异构体名称'
          : '中等证据只支持异构体组，不能释放单一异构体名称')),
      };
    }
    if (key === 'candidate') return { label: '中等证据候选', css: 'candidate', note: String(conclusion.note || '身份范围维持为候选') };
    if (key === 'background') return { label: '空白／背景相关', css: 'unassessed', note: String(conclusion.note || '不计作样本化合物') };
    if (key === 'unresolved_event') {
      const incompleteFormula = unresolvedEventKind(group) === 'incomplete_formula';
      return {
        label: isEnglishInterface() ? 'Unresolved precursor event' : '未解析母离子事件',
        css: 'unassessed',
        note: isEnglishInterface()
          ? (incompleteFormula
            ? 'Candidate structures disagree and at least one candidate lacks a usable molecular formula. The candidates are retained for review, but the event is excluded from the independent compound count.'
            : 'Known candidate molecular formulae conflict. The candidates are retained for review, but the event is excluded from the independent compound count.')
          : String(conclusion.note || (incompleteFormula
            ? '候选结构不一致且至少一个候选缺少可用分子式；不计作独立化合物'
            : '已知候选分子式互相冲突；不计作独立化合物')),
      };
    }
    const match = group?.rep?.match || {};
    const decision = match?.final_claim_decision || {};
    if (!decision.policy_id) return { label: '证据不足', css: 'unassessed', note: String(conclusion.note || '仅保留为低证据线索') };
    if (decision.released) return { label: '可作结论', css: 'eligible', note: '已通过当前结论条件' };
    const reasonLabels = {
      project_level_qvalue_unavailable: '缺项目级 q 值',
      project_level_target_decoy_qvalue_not_yet_available: '缺项目级 q 值',
      insufficient_independent_injections: '独立进样不足',
      insufficient_independent_injection_groups: '独立进样不足',
      same_method_authentic_standard_rt_required: '需同方法标品 RT',
      project_qvalue_above_release_threshold: '项目 q 值未过门槛',
      blank_background_conflict: '空白背景冲突',
    };
    const reasons = [...new Set((decision.reason_codes || []).map(String))]
      .map(code => reasonLabels[code] || '').filter(Boolean);
    return { label: '未达到当前结论条件', css: 'candidate', note: reasons.slice(0, 2).join('；') || '尚缺结论条件' };
  }

  function contaminantLabel(match) {
    if (!match?.potential_contaminant) return '';
    const prefix = match.compound_role === 'pesticide_screening_target'
      ? '农药筛查对象' : '潜在外源背景';
    return match.contaminant_category
      ? `${prefix} · ${match.contaminant_category}` : prefix;
  }

  function massErrorPenalty(value) {
    const ppm = Math.abs(Number(value));
    if (!Number.isFinite(ppm) || ppm <= 50) return 1;
    if (ppm <= 200) return 1 - 0.25 * (ppm - 50) / 150;
    return Math.max(0.35, 0.75 - 0.40 * Math.min(ppm - 200, 300) / 300);
  }

  function expertWeight() {
    return 0;
  }

  function expertWeightedTerm(match) {
    const score = Number(match?.expert_score);
    if (!Number.isFinite(score) || score <= 0) return 0;
    return score * expertWeight() * massErrorPenalty(match?.mass_ppm);
  }

  function median(values) {
    const sorted = values.filter(Number.isFinite).sort((a, b) => a - b);
    if (!sorted.length) return null;
    const middle = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
  }

  function expertImpactSection(bundle) {
    return `<div class="report-expert-impact">
      <div class="report-expert-impact-head"><div><b>分类、生源、规则证据与碎裂解释</b><span>逐化合物保留展示，不直接参与身份排序</span></div><strong>解释信息</strong></div>
      <p>这些模块会在下表中给出实际推理结果；无法形成候选特异结论时明确标为不可用。推测碎裂途径不等同于已经实验确认的机理。</p>
      <small>身份排序与解释输出分开记录，避免缺失解释被误当作没有运行。</small>
    </div>`;
  }

  function structuredRuleEvidenceLabel(expert) {
    const english = isEnglishInterface();
    const dimensions = Array.isArray(expert?.dimensions)
      ? expert.dimensions : [];
    if (!dimensions.length) {
      return english
        ? 'No auditable dimension-level rule evidence was saved; no rule-support claim is made'
        : '未保存可审计的逐维规则证据；不作规则支持声明';
    }
    const label = dimension => String(
      (english ? dimension.label_en : dimension.label)
      || dimension.label_en || dimension.label || dimension.key
      || (english ? 'Unnamed dimension' : '未命名维度')
    ).trim();
    const supported = dimensions.filter(dimension => (
      dimension.status === 'scored'
      && Number(dimension.evidence_count || 0) > 0
    ));
    const noEvidence = dimensions.filter(dimension => (
      dimension.status === 'no_evidence'
    )).length;
    const notApplicable = dimensions.filter(dimension => (
      dimension.status === 'not_applicable'
    )).length;
    const unavailable = dimensions.filter(dimension => (
      dimension.status === 'unavailable'
    )).length;
    if (english) {
      const parts = supported.length
        ? [`Observed support: ${supported.slice(0, 3).map(dimension => `${label(dimension)} (${integer(dimension.evidence_count)})`).join(', ')}`]
        : ['All dimensions were evaluated; no positive rule support was observed'];
      if (noEvidence) parts.push(`${integer(noEvidence)} dimensions without observed support`);
      if (notApplicable) parts.push(`${integer(notApplicable)} not applicable`);
      if (unavailable) parts.push(`${integer(unavailable)} unavailable`);
      return parts.join('; ');
    }
    const parts = supported.length
      ? [`观测支持：${supported.slice(0, 3).map(dimension => `${label(dimension)}（${integer(dimension.evidence_count)}）`).join('、')}`]
      : ['已逐维评价，未观察到正向规则支持'];
    if (noEvidence) parts.push(`未观察 ${integer(noEvidence)} 维`);
    if (notApplicable) parts.push(`不适用 ${integer(notApplicable)} 维`);
    if (unavailable) parts.push(`未建立 ${integer(unavailable)} 维`);
    return parts.join('；');
  }

  function interpretationCell(group) {
    const english = isEnglishInterface();
    const match = group?.rep?.match || {};
    const auxiliary = group?.compound?.auxiliary_interpretation
      || match.auxiliary_interpretation || {};
    const biosource = auxiliary.biosource || {};
    const fragmentation = auxiliary.fragmentation || {};
    const expert = auxiliary.expert || {};
    const pathways = Array.isArray(biosource.pathways) ? biosource.pathways : [];
    const rawBiosourceLabel = String(pathways[0] || biosource.reason_zh
      || '当前生源知识未覆盖（不影响谱库身份等级）')
      .replace(/^生源\s*[:：]\s*/, '');
    const biosourceLabel = english
      ? englishBiosourceCopy(rawBiosourceLabel)
      : rawBiosourceLabel;
    const steps = Array.isArray(fragmentation.steps) ? fragmentation.steps : [];
    const firstStep = steps[0] || {};
    const from = Number(firstStep.from_mz);
    const to = Number(firstStep.to_mz);
    const loss = String((english ? firstStep.neutral_loss_en : firstStep.neutral_loss)
      || firstStep.neutral_loss_en || firstStep.neutral_loss
      || firstStep.formula || firstStep.mechanism || '').trim();
    const productStructureCount = Number(fragmentation.product_structure_count || 0);
    const hasLibraryPathwayReference = (
      fragmentation.candidate_specific_library_pathway?.status === 'available'
    );
    const fragmentText = Number.isFinite(from) && Number.isFinite(to)
      ? `${from.toFixed(4)} → ${to.toFixed(4)}${loss ? (english
        ? `; proposed loss/transformation: ${loss}`
        : `，推测丢失／转化 ${loss}`) : ''}`
      : productStructureCount > 0
        ? (english
          ? `${integer(productStructureCount)} observed peaks support candidate product structures in this sample`
          : `本次样品有 ${integer(productStructureCount)} 个观测峰支持的候选产物结构`)
        : Number(fragmentation.step_count || steps.length) > 0
          ? (english
            ? `${integer(fragmentation.step_count || steps.length)} fragmentation steps were reconstructed from this sample`
            : `已从本次样品重建 ${integer(fragmentation.step_count || steps.length)} 个步骤`)
          : hasLibraryPathwayReference
            ? (english
              ? 'No sample fragmentation steps were reconstructed; a fragmentation-feature reference is available only for the library candidate'
              : '未重建样品碎裂步骤；候选谱库有碎裂特征参考')
            : Number(fragmentation.neutral_loss_count || 0) > 0
              ? (english
                ? `${integer(fragmentation.neutral_loss_count)} mass-difference or neutral-loss observations are available, but no continuous pathway was reconstructed`
                : `仅有 ${integer(fragmentation.neutral_loss_count)} 项质量差/中性丢失证据，未重建连续途径`)
              : (english
                ? 'Current evidence is insufficient to form a candidate-specific fragmentation pathway'
                : '当前证据不足，未形成候选特异途径');
    return `<div><b>分类：</b>${esc(displayClassification(match))}</div>
      <small><b>生源背景：</b>${esc(biosourceLabel)}</small>
      <small><b>碎裂：</b>${esc(fragmentText)}</small>
      <small><b>规则证据：</b>${esc(structuredRuleEvidenceLabel(expert))}</small>`;
  }

  function candidateMembersForGroup(group) {
    const suppliedMembers = Array.isArray(group.candidateMembers)
      ? group.candidateMembers : [];
    const leadName = String(group.leadCandidate?.compound_name
      || group.rep?.match?.compound_name || '');
    if (suppliedMembers.length) return suppliedMembers;
    return (group.candidateNames || []).map((name, index) => ({
      compound_name: name,
      display_order: index + 1,
      is_ranking_lead: String(name) === leadName,
      name_only_fallback: true,
    }));
  }

  function groupHasIsomerCompetition(group) {
    const match = group?.rep?.match || {};
    const finalCompound = group?.compound;
    const hasFinalCompound = finalCompound
      && typeof finalCompound === 'object'
      && Object.keys(finalCompound).length > 0;
    const finalScope = String(finalCompound?.identity_scope || '').toLowerCase();
    if (group?.status?.key === 'isomer') return true;
    if (group?.conclusion?.isomer_group_only === true) return true;
    if ([
      'isomer_group',
      'molecular_formula_isomer_group',
      'stereoisomer_group',
    ].includes(finalScope)) return true;
    // The final compound conclusion is authoritative.  Do not resurrect a
    // wider, stale scan-level Top-K ambiguity after it resolved to one 2D
    // connectivity.
    if (hasFinalCompound) return false;
    return match.isomer_ambiguous === true
      && candidateMembersForGroup(group).length > 1;
  }

  function isomerCandidateTableForGroup(group) {
    const members = candidateMembersForGroup(group);
    const unresolvedEvent = group?.status?.key === 'unresolved_event';
    const incompleteFormulaEvent = unresolvedEvent
      && unresolvedEventKind(group) === 'incomplete_formula';
    const english = isEnglishInterface();
    const leadName = String(group.leadCandidate?.compound_name
      || group.rep?.match?.compound_name || '');
    const rows = members.map((member, memberIndex) => {
      const isLead = member.is_ranking_lead
        || String(member.compound_name || '') === leadName;
      const role = isLead
        ? (unresolvedEvent
          ? (english
            ? '<b class="report-isomer-lead">Event ranking lead</b><small>Not a compound identity</small>'
            : '<b class="report-isomer-lead">事件排序先导</b><small>不作为化合物身份</small>')
          : (english
            ? '<b class="report-isomer-lead">Injection Top-1 ranking lead</b><small>Does not resolve the individual isomer</small>'
            : '<b class="report-isomer-lead">本针 Top1 排序先导</b><small>不等同于单体已解析</small>'))
        : (member.name_only_fallback
          ? (english ? 'Candidate name retained from a legacy result' : '旧结果保留的候选名称')
          : (unresolvedEvent
            ? (english
              ? (incompleteFormulaEvent
                ? 'Candidate in an event with incomplete formula assignment'
                : 'Formula-conflict candidate')
              : (incompleteFormulaEvent
                ? '分子式记录不完整事件中的候选'
                : '分子式冲突候选'))
            : (english ? 'Spectrally competing candidate' : '谱学竞争候选')));
      return `<tr>
        <td><b>${integer(member.rank || member.display_order || memberIndex + 1)}</b></td>
        <td>${role}</td>
        <td class="report-isomer-identity"><b data-scientific-value>${esc(member.compound_name || 'Unnamed candidate')}</b><small>${esc(member.formula || group.compound?.formula || (english ? 'Formula not recorded' : '分子式未记录'))}</small><small class="mono">InChIKey ${esc(member.inchikey || '—')}</small><small class="mono">SMILES ${esc(member.smiles || '—')}</small></td>
        <td>${esc(member.adduct || '—')}</td>
        <td class="mono"><b>Cosine ${member.cosine == null ? '—' : decimal(member.cosine, 3)}</b><small>Entropy ${member.entropy == null ? '—' : decimal(member.entropy, 3)}</small></td>
        <td class="mono"><b>${english ? 'Ranking score' : '排序分'} ${member.rank_score == null ? '—' : decimal(member.rank_score, 3)}</b><small>Final ${member.final_score == null ? '—' : decimal(member.final_score, 3)}</small></td>
        <td>${esc(member.relationship || '—')}</td>
        <td>${esc(member.dataset_source || '—')}<small>${esc(member.lib_spectrum_id || (english ? 'Spectrum ID not recorded' : '谱图 ID 未记录'))}</small></td>
      </tr>`;
    }).join('');
    return {
      count: members.length,
      hasNameOnlyFallback: members.some(member => member.name_only_fallback),
      html: rows
        ? `<div class="table-wrap report-isomer-table-wrap"><table class="sample-report-table report-isomer-candidate-table"><colgroup><col style="width:5%"><col style="width:13%"><col style="width:29%"><col style="width:8%"><col style="width:11%"><col style="width:11%"><col style="width:10%"><col style="width:13%"></colgroup><thead><tr><th>#</th><th>${english ? 'Role within group' : '组内角色'}</th><th>${english ? 'Candidate identity and structural identifiers' : '候选身份与结构标识'}</th><th>${english ? 'Adduct' : '加合物'}</th><th>${english ? 'Spectral similarity' : '谱图相似性'}</th><th>${english ? 'Ranking evidence' : '排序证据'}</th><th>${english ? 'Candidate relationship' : '候选关系'}</th><th>${english ? 'Library source / spectrum' : '谱库来源／谱图'}</th></tr></thead><tbody>${rows}</tbody></table></div>`
        : `<p class="report-empty">${english ? 'This legacy task did not retain expandable candidate-member records.' : '该旧任务没有保存可展开的候选成员记录。'}</p>`,
    };
  }

  function compoundCards(groups, limit = 30) {
    return groups.slice(0, limit).map((group, index) => {
      const english = isEnglishInterface();
      const row = group.rep;
      const match = row.match || {};
      const hasIsomerCompetition = groupHasIsomerCompetition(group);
      const isUnresolvedEvent = group.status.key === 'unresolved_event';
      const incompleteFormulaEvent = isUnresolvedEvent
        && unresolvedEventKind(group) === 'incomplete_formula';
      const hasCandidateConflict = hasIsomerCompetition || isUnresolvedEvent;
      const scopeLabels = {
        confirmed: '单一化合物，已确证',
        high: '单一二维连接结构，High 级检出',
        isomer: 'High 级分子式异构体组',
        medium: '中等证据结构候选',
        background: '空白／背景相关',
        unresolved_event: english
          ? (incompleteFormulaEvent
            ? 'Candidate structures disagree and at least one candidate lacks a usable molecular formula; excluded from the independent compound count'
            : 'Known candidate molecular formulae conflict; excluded from the independent compound count')
          : (incompleteFormulaEvent
            ? '候选结构不一致且至少一个候选缺少可用分子式；不计作独立化合物'
            : '已知候选分子式互相冲突的母离子事件；不计作独立化合物'),
        low: '低证据谱库线索',
      };
      const scope = group.status.key === 'isomer'
        ? (group.status.hasIsomerHigh
          ? 'High 级分子式异构体组'
          : '中等证据分子式异构体组')
        : (scopeLabels[group.status.key] || '谱库候选');
      const claim = claimConclusion(group);
      const contaminant = contaminantLabel(match);
      const adductInterpretation = group.compound?.auxiliary_interpretation?.adduct_assignment
        || match.auxiliary_interpretation?.adduct_assignment || {};
      const sampleIon = String(adductInterpretation.sample_ion_form || '').trim();
      const candidateIon = String(
        adductInterpretation.candidate_ion_form_hypothesis
        || match.candidate_adduct_hypothesis || ''
      ).trim();
      const confidenceDisplay = isUnresolvedEvent
        ? (isEnglishInterface()
          ? (incompleteFormulaEvent
            ? 'Unresolved precursor event (incomplete candidate-formula assignment; excluded from the compound count)'
            : 'Unresolved precursor event (conflicting known candidate formulae; excluded from the compound count)')
          : (incompleteFormulaEvent
            ? '未解析母离子事件（候选分子式记录不完整；不计作独立化合物）'
            : '未解析母离子事件（已知候选分子式冲突；不计作独立化合物）'))
        : hasIsomerCompetition
        ? (isEnglishInterface()
          ? `${match.isomer_group_confidence || 'Unassigned'} (isomer group) / ${match.exact_structure_confidence || match.confidence || 'Unassigned'} (full-structure identity)`
          : `${match.isomer_group_confidence || '未分级'}（异构体组） / ${match.exact_structure_confidence || match.confidence || '未分级'}（完整结构身份）`)
        : (group.status.key === 'high'
          ? (isEnglishInterface()
            ? `${match.two_dimensional_connectivity_confidence || match.confidence || 'Unassigned'} (2-D connectivity) / ${match.exact_structure_confidence || match.confidence || 'Unassigned'} (full-structure identity)`
            : `${match.two_dimensional_connectivity_confidence || match.confidence || '未分级'}（二维连接） / ${match.exact_structure_confidence || match.confidence || '未分级'}（完整结构身份）`)
          : (match.two_dimensional_connectivity_confidence || match.exact_structure_confidence || match.confidence || '未分级'));
      const precursorGateNote = match.precursor_high_evidence_gate_required
        ? `High 已降为 Medium：${match.precursor_high_evidence_gate_reason_zh || match.precursor_high_evidence_gate_reason || '碎片配对稀疏且缺少独立前体支持'}`
        : '';
      const leadName = String(group.leadCandidate?.compound_name
        || match.compound_name || '未命名候选');
      const candidateCount = Math.max(
        Number(group.candidateMemberCount || 0),
        Number(match.isomer_group_size || 0),
      );
      const candidateSummary = isUnresolvedEvent
        ? (english
          ? `<small><b>Event ranking lead for review only:</b> <span data-scientific-value>${esc(leadName)}</span><span> (not used as a compound name)</span></small><small>${incompleteFormulaEvent ? `${integer(candidateCount || 1)} candidates disagree structurally and at least one lacks a usable molecular formula` : `${integer(candidateCount || 1)} candidates have conflicting known molecular formulae`}; the complete list is shown below</small>`
          : `<small><b>仅供复核的事件排序先导：</b><span data-scientific-value>${esc(leadName)}</span><span>（不作为化合物名称）</span></small><small>${incompleteFormulaEvent ? `${integer(candidateCount || 1)} 个候选结构不一致，且至少一个候选缺少可用分子式` : `${integer(candidateCount || 1)} 个候选的已知分子式互相冲突`}；完整清单见下表</small>`)
        : (hasIsomerCompetition
          ? (english
            ? `<small><b>Injection Top-1 ranking lead:</b> <span data-scientific-value>${esc(leadName)}</span><span> (not a resolved individual isomer)</span></small><small>${integer(candidateCount || 1)} spectrally competing candidates; see the complete isomer-candidate ranking details</small>`
            : `<small><b>本针 Top1 排序先导：</b><span data-scientific-value>${esc(leadName)}</span><span>（不是已解析单体）</span></small><small>${integer(candidateCount || 1)} 个谱学竞争候选；完整清单见“异构体候选排序明细”</small>`)
          : '');
      const replicate = group.replicateConsistency || {};
      const replicateSummary = replicate.injection_count
        ? `<small>重复针：检出 ${integer(replicate.detected_injection_count)}/${integer(replicate.injection_count)} · High ${integer(replicate.high_injection_count)}/${integer(replicate.injection_count)}${replicate.robust_high ? ' · 重复稳健 High' : ''}</small>` : '';
      const conflictCandidates = hasCandidateConflict
        ? isomerCandidateTableForGroup(group) : null;
      const formulaDisplay = isUnresolvedEvent
        ? ((group.compound?.candidate_formulas || []).join(' / ')
          || (english ? 'No usable candidate formula recorded' : '无可用候选分子式'))
        : (group.compound?.formula || match.formula || (english ? 'Formula not recorded' : '分子式未记录'));
      const candidateConflictSummary = !hasCandidateConflict ? '' : (english
        ? `<small>${integer(candidateCount || 1)} ${isUnresolvedEvent ? 'unresolved' : 'conflicting'} candidates; ${isUnresolvedEvent ? (incompleteFormulaEvent ? 'candidate-formula assignment is incomplete and the event is excluded from the independent compound count' : 'known formulae conflict and the event is excluded from the independent compound count') : 'the injection Top-1 ranking lead is shown'}</small>`
        : `<small>${integer(candidateCount || 1)} 个${isUnresolvedEvent ? '未解析' : '冲突'}候选；${isUnresolvedEvent ? (incompleteFormulaEvent ? '候选分子式记录不完整，不计作独立化合物' : '已知分子式互相冲突，不计作独立化合物') : '当前显示本针 Top1 排序先导'}</small>`);
      return `<details class="report-compound-card" data-report-disclosure data-report-kind="compound">
        <summary>
          <span class="report-compound-number">${String(index + 1).padStart(2, '0')}</span>
          <span class="report-compound-identity"><b${isUnresolvedEvent ? '' : ' data-scientific-value'}>${esc(displayName(group))}</b><small>${esc(formulaDisplay)} · ${english ? 'sample ion form' : '样品离子形式'} ${esc(sampleIon || (english ? 'not reliably assigned' : '未可靠判定'))}</small>${candidateConflictSummary}</span>
          <span class="report-compound-scope"><span class="sample-evidence-badge evidence-${group.status.key}">${esc(group.status.label)}</span><small>${esc(scope)}</small><small data-entity-decision-reasons>${esc(SampleResultData.entityDecisionNote(group.entity, undefined, match.entity_high_evidence))}</small></span>
          <span class="report-compound-metrics"><b>${integer(group.spectra)} 张 MS2</b><small>Cosine ${decimal(match.cosine, 3)} · ${decimal(match.mass_ppm, 2)} ppm</small></span>
          <span class="report-disclosure-action" aria-hidden="true"><span class="when-closed">展开详情</span><span class="when-open">收起详情</span></span>
        </summary>
        <div class="report-compound-body">
          <div class="report-compound-evidence-grid">
            <div><span>鉴定范围与结论边界</span><b>${esc(confidenceDisplay)}</b><p>${esc(claim.note)}</p>${precursorGateNote ? `<small class="report-inline-warning">${esc(precursorGateNote)}</small>` : ''}</div>
            <div><span>分类、生源、碎裂与规则证据</span>${interpretationCell(group)}</div>
            <div><span>支持与代表 MS2 原始定位</span><div class="report-support-count"><b>${integer(group.spectra)} 张 MS2</b><small>${integer(group.peakClusters || 0)} 个色谱峰区</small></div>${rawSpectrumLocatorCell(row)}</div>
            <div><span>离子身份与排序证据</span><b>综合排序分 ${decimal(match.rank_score ?? match.final_score, 3)}</b><small>谱图相似度 ${decimal(match.cosine, 3)}</small><small>质量误差 ${decimal(match.mass_ppm, 2)} ppm</small><small>样品离子形式：${esc(sampleIon || '未可靠判定')}</small>${candidateIon ? `<small title="由候选结构/分子式与谱库前体反推；不是样品独立证据">候选离子假设：${esc(candidateIon)}（非独立证据）</small>` : ''}${contaminant ? `<small class="report-inline-warning" title="${esc(match.contaminant_note || contaminant)}">${esc(contaminant)}</small>` : ''}${replicateSummary}</div>
          </div>
          ${candidateSummary ? `<div class="report-isomer-inline-head"><div><b>${english ? (isUnresolvedEvent ? (incompleteFormulaEvent ? 'All candidates in this incomplete-formula precursor event' : 'All formula-conflict candidates for this precursor event') : 'Complete competing candidates in this isomer group') : (isUnresolvedEvent ? (incompleteFormulaEvent ? '本分子式记录不完整母离子事件的全部候选' : '本母离子事件的全部分子式冲突候选') : '本异构体组的完整竞争候选')}</b>${candidateSummary}</div><span>${integer(conflictCandidates?.count || candidateCount || 0)} ${english ? 'candidates' : '个'}</span></div>${conflictCandidates?.html || ''}` : ''}
        </div>
      </details>`;
    }).join('');
  }

  function compoundTable(groups, title, note, limit = 30, collapsible = false) {
    const defaultOpen = groups.length > 0 && !collapsible;
    return `<section class="sample-report-section report-tier-section${collapsible ? ' sample-report-secondary' : ''}">
      <details class="sample-report-details report-tier-details" data-report-disclosure data-report-kind="tier" data-report-primary="${defaultOpen ? 'true' : 'false'}"${defaultOpen ? ' open' : ''}>
        <summary class="report-section-summary"><span><b>${esc(title)}</b><small>${esc(note)}</small></span><span class="report-section-summary-meta"><strong>${integer(groups.length)} 个</strong><em class="when-closed">展开</em><em class="when-open">收起</em></span></summary>
        <div class="sample-report-details-body">${groups.length ? `<div class="report-compound-list">${compoundCards(groups, limit)}</div>` : '<p class="report-empty">本样本没有符合该层级的候选。</p>'}${groups.length > limit ? `<p class="report-table-note">本综合报告列出前 ${limit} 个；全部 ${integer(groups.length)} 个去重结果及逐谱证据请见“鉴定详情”。</p>` : ''}</div>
      </details>
    </section>`;
  }

  function isomerCandidateAppendix(groups) {
    let hasNameOnlyFallback = false;
    let candidateTotal = 0;
    const groupDetails = groups.map((group, groupIndex) => {
      const leadName = String(group.leadCandidate?.compound_name
        || group.rep?.match?.compound_name || '');
      const candidateTable = isomerCandidateTableForGroup(group);
      candidateTotal += candidateTable.count;
      if (candidateTable.hasNameOnlyFallback) hasNameOnlyFallback = true;
      return `<details class="report-isomer-group" data-report-disclosure data-report-kind="isomer-group">
        <summary><span class="report-compound-number">${String(groupIndex + 1).padStart(2, '0')}</span><span><b>${esc(group.compound?.formula || group.rep?.match?.formula || '分子式未记录')}</b><small>本针 Top1 排序先导：<span data-scientific-value>${esc(leadName || 'Unnamed candidate')}</span></small></span><span class="report-section-summary-meta"><strong>${integer(candidateTable.count)} 个候选</strong><em class="when-closed">展开候选</em><em class="when-open">收起候选</em></span></summary>
        <div class="report-isomer-group-body">${candidateTable.html}</div>
      </details>`;
    });
    if (!groupDetails.length) return '';
    return `<section class="sample-report-section sample-report-secondary">
      <details class="sample-report-details report-tier-details report-isomer-appendix" data-report-disclosure data-report-kind="appendix" data-report-primary="false">
        <summary class="report-section-summary"><span><b>3.9 异构体候选排序明细</b><small>按异构体组展开；第 1 名只是本针排序先导，组内歧义未解析时不能作为单体确证名称。</small></span><span class="report-section-summary-meta"><strong>${integer(groups.length)} 组 / ${integer(candidateTotal)} 个候选</strong><em class="when-closed">展开</em><em class="when-open">收起</em></span></summary>
        <div class="sample-report-details-body"><p class="report-section-intro">候选逐行展示，不再把多个名称拼接到化合物名称列。${hasNameOnlyFallback ? '旧结果未保存全部候选级分值的行以“—”显示；可在鉴定详情 Top-K 回查，重新投影后会补齐候选成员证据。' : ''}</p><div class="report-isomer-group-list">${groupDetails.join('')}</div></div>
      </details>
    </section>`;
  }

  function evidenceOverview(bundle) {
    const groups = bundle.groups || [];
    const statusOrder = [
      ['confirmed', '已确证存在', '#0f6a4e', group => group.status.key === 'confirmed'],
      ['high', 'High 单一二维连接', '#1b7f5a', group => group.status.key === 'high'],
      ['isomer-high', 'High 异构体组', '#c9871c', group => group.status.key === 'isomer' && group.status.hasIsomerHigh],
      ['isomer-non-high', '中等证据异构体组', '#d4a44c', group => group.status.key === 'isomer' && !group.status.hasIsomerHigh],
      ['medium', '中等证据候选', '#3977b8', group => group.status.key === 'medium'],
      ['background', '空白背景比较提示', '#a76c35', group => group.status.key === 'background'],
      ['unresolved-event', '未解析母离子事件', '#b45309', group => group.status.key === 'unresolved_event'],
      ['low', '证据不足', '#9aa8ad', group => group.status.key === 'low'],
    ];
    const evidenceCounts = new Map(statusOrder.map(([key, , , predicate]) => [key, groups.filter(predicate).length]));
    const evidenceRows = statusOrder.map(([key, label, color]) => {
      const value = evidenceCounts.get(key) || 0;
      const width = groups.length ? Math.max(2, value / groups.length * 100) : 0;
      return `<div class="report-bar-row"><span>${esc(label)}</span><div><i style="width:${width.toFixed(2)}%;background:${color}"></i></div><b>${integer(value)}</b></div>`;
    }).join('');
    const classes = new Map();
    groups.forEach(group => {
      const match = group.rep?.match || {};
      const label = displayClassification(match);
      const key = label === '类别尚未可靠确定'
        || label.startsWith('尚无可靠类别归属')
        ? '类别未覆盖／尚未可靠判类' : label;
      classes.set(key, (classes.get(key) || 0) + 1);
    });
    const topClasses = [...classes.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0])).slice(0, 6);
    const largestClass = Math.max(...topClasses.map(([, value]) => value), 1);
    const classRows = topClasses.map(([label, value]) => `<div class="report-bar-row"><span>${esc(label)}</span><div><i style="width:${(value / largestClass * 100).toFixed(2)}%;background:#5c7f89"></i></div><b>${integer(value)}</b></div>`).join('');
    return `<div class="report-evidence-overview">
      <div class="report-overview-card"><div class="report-overview-card-head"><b>去重化合物结论构成</b><span>重复 MS2 不重复计算，异构体未强行点名</span></div>${evidenceRows}</div>
      <div class="report-overview-card"><div class="report-overview-card-head"><b>化合物类别分布</b><span>描述性解释，不直接决定身份</span></div>${classRows || '<p>当前没有可汇总的类别信息。</p>'}</div>
    </div>`;
  }

  function chromatographySection(bundle) {
    const task = bundle.task || {};
    const workflow = task.untargeted_workflow || {};
    const attachment = workflow.ms2_attachment || {};
    const consensus = workflow.consensus_spectra || {};
    const representation = task.spectrum_representation_processing || {};
    const taskMs1Count = Number(task.ms1_scans);
    const ms1SpectrumCount = task.ms1_scans !== null
      && task.ms1_scans !== undefined
      && Number.isFinite(taskMs1Count)
      ? taskMs1Count
      : Number(currentTic?.ms1_spectrum_count
        ?? (currentTic?.ms1_spectra_available === false ? 0 : currentTic?.ms1?.length)
        ?? 0);
    const ms1TraceCount = Number(currentTic?.ms1_tic_point_count
      ?? currentTic?.ms1?.length
      ?? 0);
    const ms2Count = Number(task.total_spectra || currentTic?.ms2?.length || 0);
    const bpcAvailable = currentTic?.bpc_available === true;
    const ticFromSourceChromatogram = currentTic?.tic_source === 'mzml_tic_chromatogram_ms_level_array';
    const linked = Number(attachment.linked_ms2_spectra || 0);
    const audit = representation.representation_audit
      || representation.representation_audit_before || {};
    const nestedAuxiliary = Number(audit?.counts?.auxiliary?.spectra || 0);
    const excludedAuxiliary = Math.max(
      Number(representation.auxiliary_spectra_excluded_from_output || 0),
      nestedAuxiliary,
    );
    const picked = Number(representation.profile_spectra_picked || 0);
    const rtValues = [...(currentTic?.ms1 || []), ...(currentTic?.ms2 || [])]
      .map(point => Number(point.rt)).filter(Number.isFinite);
    const rtBounds = rtValues.reduce((bounds, value) => ({
      minimum: Math.min(bounds.minimum, value),
      maximum: Math.max(bounds.maximum, value),
    }), { minimum: Infinity, maximum: -Infinity });
    const rtRange = rtValues.length ? `${rtBounds.minimum.toFixed(2)}–${rtBounds.maximum.toFixed(2)} min` : '—';
    return `<section class="sample-report-section sample-report-chromatography" id="sample-report-chromatography">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">02 · DATA QUALITY</span><h4>数据质量与色谱概览</h4><p>先判断这份数据是否可用，再解读化合物候选。TIC 展示源文件记录的 MS1 总离子响应，BPC 展示逐扫描最强离子响应；系统不会用 TIC 或 MS2 基峰伪造缺失的 MS1 BPC。</p></div><span>${esc(rtRange)}</span></div>
      <div class="report-tic-card"><div class="report-chart-title"><div><b>总离子流色谱图（TIC）</b><span>${ms1TraceCount >= 2 ? (ticFromSourceChromatogram ? `来自源 mzML TIC chromatogram 的 ${integer(ms1TraceCount)} 个 MS1 采集点；不含 MS1 峰表` : '基于连续 MS1 survey scans') : '输入缺少连续 MS1，使用 MS2 近似曲线'}</span></div><small>线性强度坐标</small></div>${chromatogramChart(bundle.groups || [], 'tic')}</div>
      <div class="report-tic-card report-bpc-card"><div class="report-chart-title"><div><b>基峰色谱图（BPC）</b><span>${bpcAvailable ? (currentTic?.bpc_source === 'mzml_basepeak_chromatogram' ? '来自源 mzML 的独立 BPC chromatogram' : '每张连续 MS1 survey scan 的 Base Peak Intensity') : (ticFromSourceChromatogram ? '源文件未提供 MS1 峰表或 BPC chromatogram，不能由 TIC 反推' : '仅在含连续 MS1 Base Peak Intensity 的 mzML 数据中提供')}</span></div><small>线性强度坐标</small></div>${chromatogramChart(bundle.groups || [], 'bpc')}</div>
      <div class="report-qc-grid">
        <div><span>采集谱图</span><b>${integer(ms1SpectrumCount)} MS1 <i>/</i> ${integer(ms2Count)} MS2</b><small>${ticFromSourceChromatogram ? `另恢复 ${integer(ms1TraceCount)} 个源 MS1 TIC 点；没有 MS1 峰表` : `实测 RT 范围 ${esc(rtRange)}`}</small></div>
        <div><span>LC–MS 色谱特征</span><b>${integer(workflow.feature_count || 0)}</b><small>OpenMS 质量轨迹与色谱峰</small></div>
        <div><span>MS2→MS1 特征归属</span><b>${integer(linked)} <i>/</i> ${integer(ms2Count)}</b><small>${percentage(linked, ms2Count)} 唯一归属；${integer(attachment.ambiguous_ms2_spectra || 0)} 张歧义未强制指派</small></div>
        <div><span>特征级共识 MS2</span><b>${integer(consensus.actual_consensus_spectra || 0)}</b><small>${integer(consensus.multi_spectrum_consensus || 0)} 个由多张 MS2 形成</small></div>
        <div><span>候选无关代表谱</span><b>${integer(consensus.representative_spectra || 0)}</b><small>${integer(consensus.multi_spectrum_representative_groups || 0)} 个特征含多张 MS2，当前仍选代表谱</small></div>
        <div><span>质心化与表征核验</span><b>${integer(picked || ms1SpectrumCount + ms2Count)}</b><small>仅对实际谱图执行 centroid/profile 合同检查</small></div>
        <div><span>辅助检测器谱隔离</span><b>${integer(excludedAuxiliary)}</b><small>${excludedAuxiliary ? 'DAD/UV 光学谱未混入质谱分析' : '未记录需隔离的辅助检测器谱'}</small></div>
      </div>
    </section>`;
  }

  function expectedSection(bundle) {
    const task = bundle.task || {};
    const summary = task.expected_compounds_summary || currentExpected?.summary || {};
    if (!summary.requested && !task.params?.enable_expected_compounds) {
      return `<section class="sample-report-section" id="sample-report-targeted"><div class="sample-report-section-head report-major-head"><div><span class="report-kicker">04 · TARGETED LAYER</span><h4>独立靶向化合物层</h4><p>本任务未启用 Expected Compounds 连续 MS1 靶向层。</p></div></div></section>`;
    }
    const groups = Array.isArray(currentExpected?.groups) ? currentExpected.groups : [];
    const rows = groups.slice(0, 8).map((group, index) => `<tr>
      <td>${index + 1}</td>
      <td class="report-candidate-cell"><b>${esc(group.candidate_display || '未命名候选')}</b><small>${esc(group.formula || '—')} · ${esc(group.adduct || '—')} · ${esc(group.target_database || '来源未记录')}</small></td>
      <td><span class="sample-evidence-badge evidence-targeted">${esc(group.evidence_label || 'MS1 线索')}</span><small>${esc(group.result_note || '')}</small></td>
      <td class="mono">RT ${decimal(group.rt_min, 3)} min<small>m/z ${decimal(group.observed_mz, 5)} · ${Number(group.mass_error_ppm) > 0 ? '+' : ''}${decimal(group.mass_error_ppm, 2)} ppm</small></td>
      <td>${esc(group.openms_verification || '未记录')}<small>${group.isomer_ambiguous ? '同质量数据库候选 ' + integer(group.candidate_count || 1) : '候选数 ' + integer(group.candidate_count || 1)}</small></td>
    </tr>`).join('');
    const supportedGroups = Number(currentExpected?.supported_group_total || groups.length || 0);
    const catalogAudit = summary.target_catalog_quality_audit
      || currentExpected?.catalog?.quality_audit || {};
    const targetRows = Number(
      summary.target_catalog_rows
      ?? currentExpected?.catalog?.target_rows
      ?? 0
    );
    const quarantinedRows = Number(
      summary.quarantined_target_rows
      ?? catalogAudit.quarantined_target_rows
      ?? 0
    );
    const activeRows = Number(
      summary.active_target_rows
      ?? catalogAudit.active_target_rows
      ?? Math.max(0, targetRows - quarantinedRows)
    );
    const quarantineReasons = Object.entries(
      catalogAudit.quarantine_reason_counts || {}
    ).filter(([, count]) => Number(count) > 0)
      .sort((left, right) => Number(right[1]) - Number(left[1]))
      .map(([reason, count]) => `${reason} ${integer(count)}`)
      .join('；');
    return `<section class="sample-report-section sample-report-targeted" id="sample-report-targeted">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">04 · TARGETED LAYER</span><h4>独立靶向化合物线索</h4><p>Expected Compounds 用精确质量和 XIC 峰形检查预设目标；与普通 MS2 谱库鉴定分层统计，不增加上文 ${integer((bundle.groups || []).length)} 个谱库候选。</p></div><span>独立结果层</span></div>
      <div class="report-targeted-kpis"><div><b>${integer(targetRows)}</b><span>源库目标行</span></div><div><b>${integer(activeRows)}</b><span>进入计算的有效行</span></div><div><b>${integer(quarantinedRows)}</b><span>运行时隔离行</span></div><div><b>${integer(summary.detected_feature_count || 0)}</b><span>检出 XIC 特征</span></div><div><b>${integer(summary.supported_features || 0)}</b><span>同位素支持特征</span></div><div><b>${integer(supportedGroups)}</b><span>同位素支持候选组</span></div><div><b>${integer(summary.openms_feature_validation?.linked_targeted_features || 0)}</b><span>OpenMS 交叉验证</span></div><div><b>${integer(summary.candidate_target_names || 0)}</b><span>涉及候选名称</span></div></div>
      <div class="report-boundary-note" style="border-color:${quarantinedRows ? '#d97706' : '#16a34a'};background:${quarantinedRows ? '#fffbeb' : '#f0fdf4'};"><b>靶标库质量门控</b><span>源库记录保留不变；仅化学身份与质量元数据通过检查的 ${integer(activeRows)}/${integer(targetRows)} 行进入计算。${quarantinedRows ? `隔离 ${integer(quarantinedRows)} 行${quarantineReasons ? `（${esc(quarantineReasons)}）` : ''}。` : '没有记录被隔离。'}</span></div>
      <div class="report-boundary-note"><b>证据边界</b><span>这些结果目前为 MS1 精确质量/同位素支持线索；本任务中没有实验谱库确认或实测诊断 MS2 升级，不能单凭此区分异构体。</span></div>
      ${rows ? `<div class="table-wrap"><table class="sample-report-table report-targeted-table"><colgroup><col style="width:6%"><col style="width:20%"><col style="width:22%"><col style="width:26%"><col style="width:26%"></colgroup><thead><tr><th>#</th><th>候选组</th><th>当前证据</th><th>实测坐标</th><th>交叉验证</th></tr></thead><tbody>${rows}</tbody></table></div>` : '<p class="report-empty">当前没有可列出的同位素支持候选组；仅质量候选保留在靶向审计导出中。</p>'}
    </section>`;
  }

  function methodSection(bundle) {
    const task = bundle.task || {};
    const params = task.params || {};
    const workflow = task.untargeted_workflow || {};
    const processing = task.peak_processing_summary || {};
    const expected = task.expected_compounds_summary || {};
    const libraryCount = currentSystem.database?.loaded_spectra;
    const libraryLabel = params.db_mode === 'full' ? '完整运行谱库' : `${params.db_mode || '当前'}运行谱库`;
    const ms1Available = task.module_availability?.continuous_ms1_available !== false;
    const expectedStatus = String(expected.status || '').toLowerCase();
    const vendorConversion = task.vendor_conversion || {};
    const acquisition = vendorConversion.source?.acquisition || {};
    const sourceDescription = task.source_raw_filename && task.source_raw_filename !== task.filename
      ? `源数据 ${esc(task.source_raw_filename)} 经 ${esc(vendorConversion.format || '标准化流程')} 转换为 ${esc(task.filename || 'mzML')}`
      : `输入文件为 ${esc(task.filename || '—')}`;
    const acquisitionFacts = [
      acquisition.instrument ? `仪器 ${esc(acquisition.instrument)}` : '',
      acquisition.acquisition_mode ? `采集模式 ${esc(acquisition.acquisition_mode)}` : '',
      acquisition.method_name ? `方法 ${esc(acquisition.method_name)}` : '',
      Number.isFinite(Number(acquisition.isolation_width_da ?? params.isolation_window_width_da))
        ? `隔离窗 ${decimal(acquisition.isolation_width_da ?? params.isolation_window_width_da, 2)} Da`
        : '',
      Array.isArray(acquisition.collision_energies) && acquisition.collision_energies.length
        ? `碰撞能 ${acquisition.collision_energies.map(value => decimal(value, 1)).join('/')} eV`
        : '',
    ].filter(Boolean).join('，');
    const featureLinkPpm = decimal(params.ms2_feature_link_tolerance_ppm ?? 20, 1);
    const searchMassMethod = ms1Available
      ? `MS2 隔离目标使用 RT、电荷及 ${featureLinkPpm} ppm（低质量端 0.003 Da 下限）的窗口关联连续 MS1 特征。关联结果唯一且通过一致性检查时，谱库检索使用实测 MS1 特征质心；mzML selected-ion m/z 仅保留为隔离设定值和原始谱定位。未关联、歧义、电荷冲突或质量差异常时，系统明确回退到 selected-ion m/z。`
      : '没有连续 MS1 时，谱库检索只能使用 mzML selected-ion m/z，并在结果中标记为回退质量轴。';
    const msmsSearchMethod = isEnglishInterface()
      ? `Query and library spectra used the same peak-processing strategy (${esc(params.peak_processing_preset || processing.preset || 'library_match')}), retaining at most ${integer(params.top_n_peaks || processing.top_n || 15)} major fragment peaks. The selected search-mass axis used an HRMS candidate window of ${decimal(params.mz_tolerance_da, 3)} Da / ${decimal(params.mz_tolerance_ppm, 1)} ppm; the minimum Cosine and composite-score thresholds were ${decimal(params.cosine_threshold, 2)} and ${decimal(params.final_score_threshold, 2)}, respectively. The ${params.db_mode === 'full' ? 'complete runtime library' : `${esc(params.db_mode || 'current')} runtime library`}${libraryCount ? `, containing ${integer(libraryCount)} searchable reference spectra,` : ''} was searched by ion polarity.`
      : `查询谱与馆藏谱采用一致的峰处理策略（${esc(params.peak_processing_preset || processing.preset || 'library_match')}），检索保留最多 ${integer(params.top_n_peaks || processing.top_n || 15)} 个主要碎片峰。对上述选定的检索质量轴应用 ${decimal(params.mz_tolerance_da, 3)} Da / ${decimal(params.mz_tolerance_ppm, 1)} ppm HRMS 候选门限；最低 Cosine 与综合分阈值分别为 ${decimal(params.cosine_threshold, 2)} 和 ${decimal(params.final_score_threshold, 2)}。使用${libraryLabel}${libraryCount ? `（当前载入 ${integer(libraryCount)} 张可检索参比谱）` : ''}，按离子极性检索。`;
    const compoundConclusionMethod = isEnglishInterface()
      ? 'Scan-level library searching first generates candidates and projects them to measured ion features. The report then consolidates named results by candidate 2-D connectivity plus formula, listing each identity once while retaining its feature count, retention-time range, and linked MS2 evidence. Different structures supported by one feature remain an unresolved isomer group or unresolved precursor event according to their formula evidence; features with no usable candidate remain separate unannotated features. Collision-energy comparability, fragmentation information, and B1/B2/B3 background states are audit annotations only and do not change the stored score, rank, or EviGate decision.'
      : '扫描级谱库检索先生成候选并归属到实测离子特征。报告再按候选二维连接与分子式归并具名结果：每个身份只列一次，同时保留其关联特征数、保留时间范围和 MS2 证据。同一特征支持多个结构时，根据分子式证据保留为未解析异构体组或未解析母离子事件；没有可用候选的特征仍作为未注释特征单列。碰撞能可比性、碎片信息量及 B1/B2/B3 背景三态只作审计标注，不改变已保存的评分、排序或 EviGate 结论。';
    return `<section class="sample-report-section" id="sample-report-method">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">07 · METHODS</span><h4>分析方法与参数</h4><p>用于复现和审计的技术说明；不干扰前文的结果阅读。</p></div></div>
      <div class="report-method-grid">
        <div><b>采集上下文与数据表征</b><p>${sourceDescription}。${acquisitionFacts ? `原始采集记录：${acquisitionFacts}。` : ''}系统识别为 ${esc(task.ion_mode || params.ion_mode || '未知')} 离子模式，共读取 ${integer(task.ms1_scans || 0)} 张 MS1 和 ${integer(task.total_spectra || 0)} 张 MS2。分析前核验质谱层级、极性及 centroid/profile 表征；对明确为 profile 的 MS1/MS2 逐谱进行高分辨峰提取，辅助光学检测器谱不进入质谱分析。</p></div>
        <div><b>LC–MS 特征与 MS2 归属</b><p>${ms1Available ? `连续 MS1 数据使用 ${esc(workflow.backend || 'OpenMS')} 完成质量轨迹、色谱峰和同位素/加合物假设处理；形成 ${integer(workflow.feature_count || 0)} 个 LC–MS 特征，并将 ${integer(workflow.ms2_attachment?.linked_ms2_spectra || 0)} 张 MS2 归属到色谱特征。` : '本文件没有可用的连续 MS1，因此结果以 MS2 前体事件为回退单位；依赖连续 MS1 的模块不作有效性声明。'} ${searchMassMethod}加合物和中性质量仅作为离子身份假设，不直接等同于化合物鉴定。</p></div>
        <div><b>MS/MS 预处理与谱库检索</b><p>${msmsSearchMethod}</p></div>
        <div><b>候选身份归并与结论分级</b><p>${compoundConclusionMethod}</p></div>
        <div><b>独立靶向层</b><p>${params.enable_expected_compounds ? `系统${expectedStatus === 'ok' || expectedStatus === 'available' ? '已' : '尝试'}运行 Expected Compounds 连续 MS1 靶向层，用精确质量和 XIC 峰形检查用户明确选择的文献目标。源库中身份冲突、无效分子式或无可执行离子信息的记录在运行时隔离，但保留原始审计记录。该层与普通 MS2 谱库鉴定分开统计，只有精确质量时不能据此区分异构体。` : '本任务未启用 Expected Compounds 连续 MS1 靶向层。'}</p></div>
        <div><b>鉴定结论条件检查</b><p>系统将证据等级与“可报告结论/尚未达到确认条件”分开。作出结论还需项目级校准、适用范围、空白背景，以及独立进样或同方法标准品保留时间等条件；报告生成不重新评分，也不改变保存的候选顺序。</p></div>
      </div>
    </section>`;
  }

  function qualitySection(bundle, summary) {
    const task = bundle.task || {};
    const funnel = task.identification_funnel || {};
    const blank = task.blank_summary || {};
    const claims = task.claim_contract_summary || {};
    const warnings = Array.isArray(task.analysis_warnings) ? task.analysis_warnings : [];
    const modules = task.module_availability || {};
    const limitationItems = [];
    limitationItems.push(`${integer(summary.unmatchedSpectra)} 张 MS2 没有谱库候选通过当前检索与分数门槛；这不等同于样本中不存在相应 LC–MS 特征。`);
    if (summary.isomer) limitationItems.push(`${integer(summary.isomer)} 个候选当前更稳妥地报告到异构体组，具体结构仍需同方法标品 RT、差异碎片或 CCS 等正交证据。`);
    if (!blank.requested) limitationItems.push('本任务未提供适用的空白样，无法把样品候选与溶剂、耗材或系统背景作完整分层。');
    else if (blank.applicable === false) limitationItems.push('已提供的空白与样本不适用或极性不匹配，因此空白评估未执行。');
    else if (blank.blank_operation === 'comparison_only') limitationItems.push(`所选空白仅用于背景比较；本任务扣除 ${integer(blank.excluded_spectra || 0)} 张 MS2，未执行跨针 MS/MS 碎片相减。${blank.method_comparability_verified ? '' : '批次／色谱方法可比性尚未确认，因此背景标记只能作为警示或结果分层。'}`);
    if (modules.continuous_ms1_available === false) limitationItems.push('输入缺少连续 MS1，IIMN、XIC 和依赖隔离纯度的部分分析模块不适用。');
    if (!Number(claims.released_count || 0)) limitationItems.push('当前候选均尚未达到正式确认条件；主要限制通常是单次进样、缺少项目内 q 值校准或尚无标品/独立重复，而不是“系统未分析出结果”。');
    warnings.forEach(warning => limitationItems.push(String(warning?.message || warning)));
    const blockerLabels = {
      competitive_same_formula_isomer: '同分子式竞争异构体',
      cosine_below_0_70: '谱图相似度不足',
      cosine_below_0_80: '谱图相似度未达到较强证据',
      cosine_below_high_threshold: '谱图相似度未达到 High 门槛',
      final_score_below_high_threshold: '综合证据分不足',
      reverse_coverage_below_0_60: '馆藏碎片反向覆盖不足',
      precursor_formula_adduct_conflict: '前体质量、分子式或加合物不一致',
      mass_error_above_50_ppm: '前体质量误差偏大',
      ion_identity_conflict: '离子身份网络存在冲突',
      class_consensus_below_0_50: '类别证据不足',
    };
    const blockers = Object.entries(funnel.non_high_primary_blockers || {})
      .sort((a, b) => Number(b[1]) - Number(a[1])).slice(0, 5);
    return `<section class="sample-report-section" id="sample-report-quality">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">05 · LIMITATIONS</span><h4>证据边界与结果局限</h4><p>明确区分“已观测证据”、“证据等级”与“鉴定结论”。</p></div></div>
      <div class="report-quality-grid">
        <div class="report-quality-main"><b>本报告需要怎样理解</b><ul>${limitationItems.map(item => `<li>${esc(item)}</li>`).join('')}</ul></div>
        <div class="report-quality-side"><b>未达到较强证据的主要原因</b>${blockers.length ? blockers.map(([key, value]) => `<div><span>${esc(blockerLabels[key] || key)}</span><strong>${integer(value)}</strong></div>`).join('') : '<p>当前任务未记录阻断原因分布。</p>'}</div>
      </div>
    </section>`;
  }

  function practicalAction(plan) {
    const actions = Array.isArray(plan.ranked_actions) ? plan.ranked_actions : [];
    return actions.find(action => [
      'targeted_ms2_multi_ce', 'narrow_isolation_ms2',
      'high_sensitivity_targeted_ms2', 'opposite_polarity_targeted_ms2',
    ].includes(action.action_kind)) || actions[0] || {};
  }

  function actionParameters(action) {
    const values = [];
    if (Array.isArray(action.collision_energies) && action.collision_energies.length) {
      values.push(`CE ${action.collision_energies.map(value => decimal(value, 0)).join('/')} ${action.collision_energy_unit || ''}`.trim());
    }
    if (action.isolation_width_da != null) values.push(`隔离窗 ${decimal(action.isolation_width_da, 2)} Da`);
    if (action.max_injection_time_ms != null) values.push(`最大注入 ${decimal(action.max_injection_time_ms, 0)} ms`);
    return values.join('；') || '参数需由仪器负责人结合平台能力确定';
  }

  function followupSection() {
    const plans = Array.isArray(currentPlans.plans) ? currentPlans.plans : [];
    if (!plans.length) {
      return `<section class="sample-report-section" id="sample-report-followup"><div class="sample-report-section-head report-major-head"><div><span class="report-kicker">06 · NEXT STEPS</span><h4>建议的确认实验</h4><p>将当前候选转化为可执行的确认路径。</p></div></div><p class="report-empty">当前尚未生成针对候选歧义的确认实验建议。可从“鉴定详情 → 确认实验建议”生成。</p></section>`;
    }
    const rows = plans.slice(0, 8).map((plan, index) => {
      const need = plan.need || {};
      const candidates = (need.candidates || []).map(candidate => candidate.compound_name).filter(Boolean);
      const action = practicalAction(plan);
      const strongest = (plan.ranked_actions || [])[0] || {};
      return `<tr><td>${index + 1}</td><td><b data-scientific-value>${esc(candidates.slice(0, 3).join(' / ') || 'Unnamed candidate group')}</b>${candidates.length > 3 ? `<small>另有 ${candidates.length - 3} 个近邻候选</small>` : ''}</td><td class="mono">${decimal(need.precursor_mz, 5)}<small>RT ${decimal(need.retention_time_min, 3)} min · ${esc(need.polarity || '—')}</small></td><td><b>${esc(action.label || '待人工设计')}</b><small>${esc(actionParameters(action))}</small></td><td>${strongest.label && strongest.action_id !== action.action_id ? `${esc(strongest.label)}<small>需要额外资源或标品</small>` : '同上'}</td><td>${esc((action.rationale || [])[0] || '需结合仪器能力人工复核')}</td></tr>`;
    }).join('');
    return `<section class="sample-report-section" id="sample-report-followup">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">06 · NEXT STEPS</span><h4>建议的确认实验</h4><p>这些建议用于解决同分异构体或近邻候选，尚未执行，不能当作本次已观测证据。</p></div><span>${integer(currentPlans.metadata?.plan_count || plans.length)} 组计划</span></div>
      <div class="table-wrap"><table class="sample-report-table report-followup-table"><thead><tr><th>#</th><th>待区分候选</th><th>目标坐标</th><th>优先可执行方案</th><th>证据最强方案</th><th>建议依据</th></tr></thead><tbody>${rows}</tbody></table></div>
      <p class="report-table-note">采集参数必须由仪器负责人审核。预测信息不作为确认性证据；未检出也不能自动作为候选反证。</p>
    </section>`;
  }

  function provenanceSection(bundle) {
    const task = bundle.task || {};
    const claims = task.claim_contract_summary || {};
    const context = claims.context || {};
    const injectionCount = new Set(context.independent_injection_groups || []).size;
    const minimumInjections = Number(context.min_independent_injections || 2);
    return `<section class="sample-report-section sample-report-provenance" id="sample-report-traceability">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">08 · TRACEABILITY</span><h4>报告范围与追溯</h4><p>保留任务、时间、运行实例和鉴定结论条件，便于归档。</p></div></div>
      <dl>
        <div><dt>任务 ID</dt><dd>${esc(task.id || currentTaskId || '—')}</dd></div>
        <div><dt>分析完成时间</dt><dd>${esc(task.end_time || '未记录')}</dd></div>
        <div><dt>报告生成时间</dt><dd>${esc(new Date().toLocaleString('zh-CN', { hour12: false }))}</dd></div>
        <div><dt>运行实例</dt><dd>${esc(currentSystem.service_generation || 'V11')} · ${esc(currentSystem.runtime_snapshot_id || '运行快照未记录')}</dd></div>
        <div><dt>鉴定结论条件</dt><dd>项目级校准 ${integer(claims.calibrated_count || 0)}/${integer(claims.record_count || 0)}；独立进样 ${integer(injectionCount)}/${integer(minimumInjections)}</dd></div>
        <div><dt>排序声明</dt><dd>报告读取保存的分析结果，不重新评分，不改变候选顺序。</dd></div>
      </dl>
    </section>`;
  }

  function referenceStandardQCReportSection() {
    const payload = currentReferenceQC;
    if (!payload?.reference_standard_set || !payload?.summary) return '';
    const s = payload.summary;
    const set = payload.reference_standard_set;
    const declared = Number(s.declared_truth_identities || 0);
    const observable = Number(s.single_injection_compatible_precursor_identities || 0);
    const ratio = value => `${integer(value)}/${integer(observable)}`;
    const rows = (payload.per_standard || []).map(item => {
      const isObservable = !!item.observable_in_injection;
      return `<tr><td><b data-scientific-value>${esc(item.name_en || item.name_zh || '—')}</b></td><td>${isObservable ? '可观测' : '无兼容 MS2 事件'}</td><td>${!isObservable ? '不适用' : (item.saved_topk_connectivity ? `已找回（最佳第 ${integer(item.best_saved_rank_connectivity)}）` : '未找回')}</td><td>${!isObservable ? '不适用' : (item.top1_connectivity_any_scan ? '是' : '否')}</td><td>${!isObservable ? '不适用' : (item.reported_scope_supported_any_scan ? '支持' : '不支持')}</td><td>${!isObservable ? '不适用' : (item.high_reportable_scope_supported_any_scan ? 'High' : '未达 High')}</td></tr>`;
    }).join('');
    return `<section class="sample-report-section" id="sample-report-reference-qc">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">QC · REFERENCE STANDARDS</span><h4>${esc(set.label || set.id)} · 混标回收审计</h4><p>真值仅用于盲检索完成后的质量控制，不参与候选生成、评分、排序或置信度校准。</p></div><span>${integer(s.declared_truth_identities)} 个声明标准品</span></div>
      <div class="sample-report-kpis"><div><strong>${integer(declared)}</strong><span>声明标准品</span></div><div><strong>${integer(observable)}/${integer(declared)}</strong><span>本针可观测</span></div><div><strong>${ratio(s.observable_saved_topk_connectivity)}</strong><span>候选层找回</span></div><div><strong>${ratio(s.observable_top1_connectivity)}</strong><span>Top1 找回</span></div><div><strong>${ratio(s.observable_reported_scope_supported)}</strong><span>可报告范围找回</span></div><div><strong>${ratio(s.observable_high_reportable_scope_supported)}</strong><span>达到 High 范围</span></div></div>
      <p class="report-section-intro"><b>可观测分母：</b>候选层 ${decimal(s.observable_saved_topk_connectivity_pct, 1)}%，Top1 ${decimal(s.top1_connectivity_pct_observable, 1)}%，可报告范围 ${decimal(s.reported_scope_supported_pct_observable, 1)}%，High 范围 ${decimal(s.high_reportable_scope_supported_pct_observable, 1)}%。${esc(payload.denominator_warning_zh || payload.denominator_warning || '')}</p>
      <div class="table-wrap"><table class="sample-report-table"><thead><tr><th>标准品</th><th>本针可观测</th><th>候选层</th><th>Top1</th><th>可报告范围</th><th>High 范围</th></tr></thead><tbody>${rows}</tbody></table></div>
    </section>`;
  }

  function replicateConsistencyReportSection(bundle) {
    const payload = bundle?.replicateConsistency || {};
    if (payload.status !== 'available') return '';
    const summary = payload.summary || {};
    const rows = (payload.injections || []).map(item => `<tr>
      <td><b>${esc(item.filename || item.task_id)}</b></td>
      <td>${integer(item.replicate_number)}</td>
      <td>${integer(item.compound_group_count)}</td>
      <td>${integer(item.main_group_count)}</td>
      <td><b>${integer(item.high_group_count)}</b></td>
    </tr>`).join('');
    return `<section class="sample-report-section" id="sample-report-replicates">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">REPEATABILITY · TECHNICAL INJECTIONS</span><h4>技术重复针一致性</h4><p>单针证据与重复稳健性分开报告；不通过重复统计篡改单针 Top1、分数或置信等级。</p></div><span>${integer(payload.injection_count)} 针</span></div>
      <div class="sample-report-kpis"><div><strong>${(payload.injections || []).map(item => integer(item.high_group_count)).join(' / ')}</strong><span>各针 High</span></div><div><strong>${integer(summary.robust_high_count)}</strong><span>重复稳健 High（≥${integer(payload.robustness_threshold)}/${integer(payload.injection_count)}）</span></div><div><strong>${integer(summary.all_injection_high_count)}</strong><span>全部重复针均 High</span></div><div><strong>${integer(summary.single_injection_high_count)}</strong><span>仅单针 High</span></div></div>
      <p class="report-section-intro">${esc(payload.interpretation_zh || '')} 成组依据为文件名末尾重复编号和完全相同的分析参数。</p>
      <div class="table-wrap"><table class="sample-report-table"><thead><tr><th>文件</th><th>重复编号</th><th>鉴定组</th><th>主要结果</th><th>High 组</th></tr></thead><tbody>${rows}</tbody></table></div>
    </section>`;
  }

  function authenticRankingReportSection(bundle) {
    const payload = bundle?.evitwin || {};
    const compounds = Array.isArray(payload.compounds) ? payload.compounds : [];
    const actualCeFeatures = Array.isArray(payload?.actual_ce_shadow?.features)
      ? payload.actual_ce_shadow.features : [];
    if (!compounds.length && !actualCeFeatures.length) return '';
    const summary = payload.summary || {};
    const rows = compounds.slice(0, 40).map(feature => {
      const winner = feature.evitwin_winner || {};
      const publicWinner = feature.public_library_winner || {};
      const decision = feature.decision || {};
      const similarity = Number(winner.median_similarity);
      const rt = Number(feature.representative_retention_time_min);
      const changed = feature.display_rank_changed_from_public_library === true;
      const competing = decision.reportable_scope === 'isomer_group'
        ? (Array.isArray(feature.candidate_scores) ? feature.candidate_scores : [])
          .filter(item => Number(item?.median_similarity) >= 0.5)
          .slice(0, 4)
          .map(item => `<span data-scientific-value>${esc(item.compound_name || item.compound_code || 'Unnamed candidate')}</span> ${Number(item.median_similarity).toFixed(3)}`)
        : [];
      const scope = decision.reportable_scope === 'isomer_group'
        ? '异构体组（不点名单一立体结构）'
        : '二维连接候选（未自动升级 Level 1）';
      return `<tr>
        <td><b data-scientific-value>${esc(winner.compound_name || winner.compound_code || 'Unnamed reference candidate')}</b><small>${esc(winner.compound_code || '')}${winner.formula ? ` · ${esc(winner.formula)}` : ''}</small>${competing.length ? `<small>组内竞争：${competing.join('；')}</small>` : ''}</td>
        <td><b>${Number.isFinite(similarity) ? similarity.toFixed(3) : '—'}</b><small>${integer(winner.collision_energy_count)} 档碰撞能 · ${integer(winner.collision_energy_votes)} 票</small></td>
        <td>${Number.isFinite(rt) ? `${rt.toFixed(3)} min` : 'RT 未记录'}</td>
        <td>${esc(scope)}</td>
        <td>${changed ? `已重新排序<small>公开谱库原首位：<span data-scientific-value>${esc(publicWinner.compound_name || 'No candidate')}</span></small>` : '与公开谱库首位一致'}</td>
      </tr>`;
    }).join('');
    const actualCeRows = actualCeFeatures.slice(0, 40).map(feature => {
      const winner = feature.shadow_winner || {};
      const publicWinner = feature.public_library_winner || {};
      const score = Number(winner.median_similarity);
      const margin = Number(winner.winner_margin);
      const energies = (feature.collision_energy_conditions || [])
        .map(value => `${Number(value).toFixed(2)} eV`).join(' / ');
      return `<tr>
        <td><b data-scientific-value>${esc(winner.compound_name || winner.compound_code || 'Unnamed reference candidate')}</b><small>${esc(winner.compound_code || '')}${winner.formula ? ` · ${esc(winner.formula)}` : ''}</small></td>
        <td><b>${Number.isFinite(score) ? score.toFixed(3) : '—'}</b><small>实际 CE ${esc(energies || '—')}${Number.isFinite(margin) ? ` · margin ${margin.toFixed(3)}` : ''}</small></td>
        <td>${integer(winner.collision_energy_count)} 个查询 CE · ${integer(winner.collision_energy_votes)} 票</td>
        <td>未解析候选组；shadow 不改排序</td>
        <td>${feature.shadow_lead_differs_from_public_library === true ? `影子先导不同<small>正式首位保持：<span data-scientific-value>${esc(publicWinner.compound_name || 'No candidate')}</span></small>` : '影子先导与正式首位一致'}</td>
      </tr>`;
    }).join('');
    return `<section class="sample-report-section" id="sample-report-authentic-ranking">
      <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">AUTHENTIC MS/MS · PRIORITY CANDIDATES</span><h4>真实标品谱优先候选</h4><p>EviTwin-MS 以同仪器、多碰撞能真实单标 MS/MS 对公开谱库已召回的候选进行独立联合排序，并在同质量异构体具备足够分离峰时加入一对一保留顺序约束。</p></div><span>${integer(compounds.length)} 个去重候选</span></div>
      <div class="sample-report-kpis"><div><strong>${integer(summary.priority_compound_count)}</strong><span>标品谱优先候选</span></div><div><strong>${integer(summary.display_rank_change_count)}</strong><span>首位显示已调整</span></div><div><strong>${integer(summary.retention_order_assignment_count)}</strong><span>保留顺序联合分配</span></div><div><strong>0</strong><span>自动新增 Level 1</span></div></div>
      <p class="report-section-intro"><b>如何理解：</b>这张表优先回答样本中由本实验室真实标品谱支持的候选。它会调整标品参考通道的显示顺序，但保留公开谱库结果供人工核对；跨方法绝对 RT 未锁定前不自动确证，同二维连接的立体异构体仅报告为异构体组。</p>
      ${compounds.length ? `<div class="table-wrap"><table class="sample-report-table"><thead><tr><th>优先候选</th><th>真实标品谱证据</th><th>样品 RT</th><th>当前可报告范围</th><th>相对公开谱库</th></tr></thead><tbody>${rows}</tbody></table></div>` : ''}
      ${actualCeFeatures.length ? `<h5 style="margin-top:1rem;">已有实际 CE 谱图比较（shadow）</h5><p class="report-section-intro">只使用查询文件中记录的 CE。参考 CE 之间只内插相似度数值，不内插谱图、不向参考范围外外推；单个查询 CE 仍只算一票。下表不改正式 Top1，也不增加确证化合物数。</p><div class="table-wrap"><table class="sample-report-table"><thead><tr><th>真品谱候选组先导</th><th>实际 CE 相似度</th><th>查询证据数</th><th>当前范围</th><th>与正式首位比较</th></tr></thead><tbody>${actualCeRows}</tbody></table></div>` : ''}
    </section>`;
  }

  function reportViewControls() {
    return `<div class="sample-report-view-tools" role="toolbar" aria-label="报告展开与收起控制">
      <div><b>交互阅读</b><span>鉴定层级和每个化合物均可展开；异构体组内可直接查看完整竞争候选。打印或另存 PDF 时会自动包含全部详情。</span></div>
      <div class="report-view-actions">
        <button type="button" data-report-action="primary">仅显示主要结果</button>
        <button type="button" data-report-action="expand">全部展开</button>
        <button type="button" data-report-action="collapse">全部收起</button>
        <span data-report-disclosure-status aria-live="polite"></span>
      </div>
    </div>`;
  }

  function updateDisclosureStatus(root) {
    if (!root?.querySelectorAll) return;
    const disclosures = [...root.querySelectorAll('details[data-report-disclosure]')];
    const opened = disclosures.filter(item => item.open).length;
    const status = root.querySelector('[data-report-disclosure-status]');
    if (status) status.textContent = `已展开 ${integer(opened)} / ${integer(disclosures.length)}`;
  }

  function setDisclosureMode(mode, root = document.getElementById('sampleReportArea')) {
    if (!root?.querySelectorAll) return;
    const disclosures = [...root.querySelectorAll('details[data-report-disclosure]')];
    disclosures.forEach(item => {
      if (mode === 'expand') item.open = true;
      else if (mode === 'collapse') item.open = false;
      else if (mode === 'primary') {
        item.open = item.dataset.reportKind === 'tier'
          && item.dataset.reportPrimary === 'true';
      }
    });
    updateDisclosureStatus(root);
  }

  function bindReportInteractions(root) {
    if (!root?.querySelectorAll || root.dataset.reportInteractionsBound === 'true') return;
    root.dataset.reportInteractionsBound = 'true';
    root.querySelectorAll('[data-report-action]').forEach(button => {
      button.addEventListener('click', () => setDisclosureMode(button.dataset.reportAction, root));
    });
    root.querySelectorAll('details[data-report-disclosure]').forEach(item => {
      item.addEventListener('toggle', () => updateDisclosureStatus(root));
    });
    updateDisclosureStatus(root);
  }

  function standaloneInteractionScript() {
    return `(function(){var root=document.getElementById('sampleReportDocument');if(!root)return;var printState=[];var all=function(){return Array.prototype.slice.call(root.querySelectorAll('details[data-report-disclosure]'));};var update=function(){var items=all();var opened=items.filter(function(item){return item.open;}).length;var status=root.querySelector('[data-report-disclosure-status]');if(status)status.textContent='已展开 '+opened+' / '+items.length;};root.querySelectorAll('[data-report-action]').forEach(function(button){button.addEventListener('click',function(){var mode=button.getAttribute('data-report-action');all().forEach(function(item){if(mode==='expand')item.open=true;else if(mode==='collapse')item.open=false;else if(mode==='primary')item.open=item.getAttribute('data-report-kind')==='tier'&&item.getAttribute('data-report-primary')==='true';});update();});});root.querySelectorAll('details[data-report-disclosure]').forEach(function(item){item.addEventListener('toggle',update);});window.addEventListener('beforeprint',function(){printState=all().map(function(item){return item.open;});all().forEach(function(item){item.open=true;});});window.addEventListener('afterprint',function(){all().forEach(function(item,index){item.open=Boolean(printState[index]);});update();});update();})();`;
  }

  function candidateIdentitySummarySection(bundle) {
    const summaries = SampleResultData.candidateIdentitySummaries(
      bundle.featureGroups || bundle.groups || [], bundle.identificationRecords || []
    );
    if (!summaries.length) return '';
    const visible = summaries.slice(0, 50);
    const rows = visible.map((item, index) => {
      const localizedName = isEnglishInterface() || !item.compoundNameCn
        ? item.compoundName
        : `${item.compoundName}（${item.compoundNameCn}）`;
      const status = item.bestStatus || { key: 'low', label: '待复核候选' };
      const rt = item.rtMin === null
        ? 'RT 未记录'
        : (item.rtMax !== null && Math.abs(item.rtMax - item.rtMin) > 0.001
          ? `RT ${item.rtMin.toFixed(3)}–${item.rtMax.toFixed(3)} min`
          : `RT ${item.rtMin.toFixed(3)} min`);
      return `<tr>
        <td>${index + 1}</td>
        <td class="report-candidate-cell"><b data-scientific-value>${esc(localizedName || 'Unnamed candidate')}</b><small>${esc(item.formula || '分子式未记录')}</small><small class="mono">${esc(item.inchikey || 'InChIKey 未记录')}</small>${item.curatedNatural ? '<small class="report-curated-natural">已核验天然产物身份记录</small>' : ''}</td>
        <td><span class="sample-evidence-badge evidence-${esc(status.key || 'low')}">${esc(status.label || '待复核候选')}</span><small>排序先导 ${integer(item.rankingLeadCount)} / ${integer(item.featureResultCount)} 条</small><small>未解竞争组 ${integer(item.ambiguousFeatureCount)} 条</small></td>
        <td><b>${integer(item.linkedFeatureCount)} 个连续 MS1 特征已关联</b><small>${integer(item.fallbackEventCount)} 个未关联 MS1 的 MS2 采集事件</small><small>重复采集不作为独立化合物检出次数。</small></td>
        <td class="mono"><b>碎片余弦 ${item.bestCosine === null ? '—' : item.bestCosine.toFixed(3)}</b><small>谱熵相似度 ${item.bestEntropy === null ? '—' : item.bestEntropy.toFixed(3)} · 反向覆盖 ${item.bestReverseMatch === null ? '—' : item.bestReverseMatch.toFixed(3)}</small><small>${integer(item.fragmentEvidenceScanCount)} 条扫描级碎片证据 · 实体排序分 ${item.bestScore === null ? '—' : item.bestScore.toFixed(3)}</small><small>最小 |质量误差| ${item.bestAbsMassPpm === null ? '—' : `${Math.abs(item.bestAbsMassPpm).toFixed(2)} ppm`}</small></td>
        <td><b>${esc(rt)}</b><small>${integer(item.supportingSpectra)} 张组级关联 MS2；竞争组内证据不能全部归给单一候选。</small></td>
      </tr>`;
    }).join('');
    return `<section class="sample-report-section report-candidate-identity-summary">
      <div class="sample-report-section-head"><div><span class="report-kicker">CANDIDATE IDENTITY INDEX</span><h4>候选身份索引（同一二维连接不重复列名）</h4><p>该索引把逐离子特征结果按候选二维连接与分子式汇总，便于找到低于 High 门槛但有碎片匹配的候选。原始特征作为支持证据保留，归并不会提高证据等级。</p></div><span>${integer(visible.length)} / ${integer(summaries.length)} 个身份</span></div>
      <div class="report-boundary-note"><b>评分规则</b><span>碎片谱相似度和质量证据单列；CE 只记录采集条件可比性，不乘入分数、不改变排序。仅有未关联 MS1 的 MS2 采集事件时保留候选，但不能写成独立色谱峰检出。</span></div>
      <div class="table-wrap"><table class="sample-report-table report-candidate-identity-table"><colgroup><col style="width:4%"><col style="width:24%"><col style="width:16%"><col style="width:23%"><col style="width:17%"><col style="width:16%"></colgroup><thead><tr><th>#</th><th>候选身份</th><th>最高证据范围</th><th>色谱／采集支持</th><th>碎片匹配与排序</th><th>保留时间范围</th></tr></thead><tbody>${rows}</tbody></table></div>
      ${summaries.length > visible.length ? `<p class="report-table-note">此索引列出前 ${integer(visible.length)} 个身份；后续分级章节列出各去重身份，并在身份内保留关联特征与扫描证据。</p>` : ''}
    </section>`;
  }

  function reportDocumentHtml(bundle) {
    const task = bundle.task || {};
    const summary = counts(bundle);
    const confirmed = bundle.groups.filter(group => group.status.key === 'confirmed');
    const high = bundle.groups.filter(group => group.status.key === 'high');
    const isomer = bundle.groups.filter(group => group.status.key === 'isomer');
    const allIsomerCompetition = bundle.groups.filter(
      groupHasIsomerCompetition
    );
    const isomerHigh = isomer.filter(group => group.status.hasIsomerHigh);
    const isomerNonHigh = isomer.filter(group => !group.status.hasIsomerHigh);
    const medium = bundle.groups.filter(group => group.status.key === 'medium');
    const background = bundle.groups.filter(group => group.status.key === 'background');
    const unresolvedEvents = bundle.groups.filter(
      group => group.status.key === 'unresolved_event'
    );
    const low = bundle.groups.filter(group => group.status.key === 'low');
    const unannotated = bundle.groups.filter(
      group => group.status.key === 'unannotated'
    );
    const topNames = [...confirmed, ...high, ...isomerHigh].slice(0, 8).map(displayName);
    const matchedRate = summary.totalSpectra ? (summary.matchedSpectra / summary.totalSpectra * 100).toFixed(1) : '0.0';
    const claims = task.claim_contract_summary || {};
    const replicate = bundle.replicateConsistency || {};
    const replicateSummary = replicate.summary || {};
    const blank = task.blank_summary || {};
    const reportLead = isEnglishInterface()
      ? `<p class="report-lead"><b>${integer(summary.matchedSpectra)}</b> MS2 scans produced library candidates (${matchedRate}% of all MS2 scans), consolidated into <b>${integer(summary.compounds)}</b> deduplicated candidate identities. The sample-supported primary set contains <b>${integer(summary.main)}</b> identities: <b>${integer(summary.highScopeGroups)}</b> pass the current entity High conditions, including <b>${integer(high.length)}</b> single 2-D connectivity results and <b>${integer(isomerHigh.length)}</b> unresolved isomer groups; an additional <b>${integer(isomerNonHigh.length)}</b> isomer groups and <b>${integer(medium.length)}</b> candidates have moderate evidence. <b>${integer(unresolvedEvents.length)}</b> precursor events remain unresolved, and <b>${integer(unannotated.length)}</b> ion features have no usable library candidate. <b>${integer(summary.insufficient)}</b> candidate identities have insufficient evidence. Authentic-standard confirmations are reported separately: <b>${integer(confirmed.length)}</b>.</p>`
      : `<p class="report-lead"><b>${integer(summary.matchedSpectra)}</b> 张 MS2 获得谱库候选（占全部 MS2 的 ${matchedRate}%），按候选二维连接与分子式归并为 <b>${integer(summary.compounds)}</b> 个去重候选身份。样品支持的主要结果为 <b>${integer(summary.main)}</b> 个身份：<b>${integer(summary.highScopeGroups)}</b> 个通过当前实体 High 条件，其中 <b>${integer(high.length)}</b> 个可报告到单一二维连接、<b>${integer(isomerHigh.length)}</b> 个只能报告到未解异构体组；另有 <b>${integer(isomerNonHigh.length)}</b> 个中等证据异构体组和 <b>${integer(medium.length)}</b> 个中等证据候选。<b>${integer(unresolvedEvents.length)}</b> 个母离子事件未解析且不计作独立化合物，<b>${integer(unannotated.length)}</b> 个离子特征没有可用谱库候选；<b>${integer(summary.insufficient)}</b> 个候选身份证据不足。正式标准品确证单列为 <b>${integer(confirmed.length)}</b> 个。</p>`;
    return `<article class="sample-report-document" id="sampleReportDocument">
      <header class="sample-report-cover">
        <div><span>MSsystem · I.D.E.A platform</span><h3>LC–MS/MS 样本综合分析报告</h3><p>${esc(task.filename || '未命名样本')}</p><em>从原始色谱与数据质量，到化合物候选、证据边界和确认建议</em></div>
        <div class="sample-report-cover-meta"><span>离子模式 <b>${esc(task.ion_mode || task.params?.ion_mode || '—')}</b></span><span>采集规模 <b>${integer(task.ms1_scans || 0)} MS1 / ${integer(task.total_spectra || 0)} MS2</b></span><span>谱库 <b>${task.params?.db_mode === 'full' ? 'Full' : esc(task.params?.db_mode || '当前库')}</b></span><span>任务 <b>${esc(String(task.id || currentTaskId).slice(0, 12))}</b></span></div>
      </header>
      <nav class="sample-report-toc" aria-label="报告目录"><b>阅读路径</b><a href="#sample-report-conclusion">01 结论</a><a href="#sample-report-confidence-guide">先读懂置信度</a>${((bundle.evitwin?.compounds || []).length || (bundle.evitwin?.actual_ce_shadow?.features || []).length) ? '<a href="#sample-report-authentic-ranking">标品谱候选比较</a>' : ''}${currentReferenceQC?.reference_standard_set ? '<a href="#sample-report-reference-qc">QC 混标回收</a>' : ''}${replicate.status === 'available' ? '<a href="#sample-report-replicates">重复针一致性</a>' : ''}<a href="#sample-report-chromatography">02 TIC 与数据质量</a><a href="#sample-report-results">03 鉴定结果</a><a href="#sample-report-targeted">04 靶向线索</a><a href="#sample-report-quality">05 局限</a><a href="#sample-report-followup">06 确认建议</a><a href="#sample-report-method">07 方法</a><a href="#sample-report-traceability">08 追溯</a></nav>
      ${reportViewControls()}
      <section class="sample-report-section sample-report-conclusion" id="sample-report-conclusion">
        <div class="sample-report-section-head report-major-head"><div><span class="report-kicker">01 · EXECUTIVE SUMMARY</span><h4>样本总体结论</h4><p>先给出可以带走的结论，再逐层展开原始数据、候选和局限。</p></div><span class="report-status-pill">分析完成</span></div>
        ${SampleResultData.confidenceAxesHtml(bundle)}
        ${reportLead}
        <div class="sample-report-kpis"><div><strong>${SampleResultData.confidenceAxes(bundle).libraryHigh === null ? "—" : integer(SampleResultData.confidenceAxes(bundle).libraryHigh)}</strong><span>${isEnglishInterface() ? "Single-spectrum library High" : "单谱谱库 High 候选身份"}</span></div><div><strong>${integer(summary.highScopeGroups)}</strong><span>实体 High 级鉴定组</span></div><div><strong>${integer(high.length)}</strong><span>High 单一二维连接</span></div><div><strong>${integer(isomerHigh.length)}</strong><span>High 异构体组</span></div><div><strong>${integer(isomerNonHigh.length)}</strong><span>中等证据异构体组</span></div><div><strong>${integer(medium.length)}</strong><span>中等证据候选</span></div><div><strong>${integer(unresolvedEvents.length)}</strong><span>未解析事件</span></div><div><strong>${integer(confirmed.length)}</strong><span>正式确证</span></div><div><strong>${integer(summary.unverifiedBlank)}</strong><span>未核验空白警示</span></div></div>
        ${summary.unverifiedBlank ? `<div class="report-inline-warning"><b>空白证据边界：</b>${integer(summary.unverifiedBlank)} 个鉴定组在同极性未定角色空白中亦有信号，其中 ${integer(summary.unverifiedBlankHigh)} 个达到 High 可报告范围；因批次／色谱方法配对未经核验，这些信号仅作背景警示，没有被自动移出样品结果。</div>` : ''}
        ${blank.requested ? `<div class="report-inline-warning"><b>${isEnglishInterface() ? 'Blank operation:' : '空白操作：'}</b>${blank.blank_operation === 'comparison_only' ? (isEnglishInterface() ? `Background comparison only; ${integer(blank.excluded_spectra || 0)} MS2 scans were excluded. No cross-injection subtraction of MS/MS fragment intensities was performed, and the selected blanks were not analyzed as sample tasks.` : `本任务仅作背景比较，扣除 ${integer(blank.excluded_spectra || 0)} 张 MS2；未执行跨针 MS/MS 碎片相减，所选空白也未作为样品任务分析。`) : (isEnglishInterface() ? `Background exclusion removed ${integer(blank.excluded_spectra || 0)} MS2 scans before searching; method comparability and exclusion authorization are documented in Traceability.` : `背景排除 ${integer(blank.excluded_spectra || 0)} 张 MS2；方法可比性和排除授权状态见追溯信息。`)}</div>` : ''}
        ${replicate.status === 'available' ? `<div class="report-inline-warning"><b>重复针稳健性：</b>各针 High 为 ${(replicate.injections || []).map(item => integer(item.high_group_count)).join(' / ')}；至少 ${integer(replicate.robustness_threshold)}/${integer(replicate.injection_count)} 针达到 High 的重复稳健组为 ${integer(replicateSummary.robust_high_count)} 个，全部重复针均 High 为 ${integer(replicateSummary.all_injection_high_count)} 个，仅单针 High 为 ${integer(replicateSummary.single_injection_high_count)} 个。</div>` : ''}
        ${(bundle.evitwin?.compounds || []).length ? `<div class="report-inline-warning"><b>真实标品谱联合排序：</b>本样本得到 ${integer(bundle.evitwin.summary?.priority_compound_count)} 个优先候选，其中 ${integer(bundle.evitwin.summary?.display_rank_change_count)} 个相对公开谱库首位重新排序；具体化合物、相似度、碰撞能覆盖和样品 RT 见下一节。</div>` : ''}
        ${topNames.length ? `<div class="report-top-candidates"><b>主要 High 级结果</b><div>${topNames.map(name => `<span data-scientific-value>${esc(name)}</span>`).join('')}</div><small>High 级鉴定组按实际身份范围报告，仍不等同于同方法标准品 Level 1 确证。</small></div>` : ''}
        <div class="report-conclusion-grid"><div><b>结果应怎样读</b><p>High 单一二维连接与 High 异构体组分开报告；二维连接 High 不自动声明立体化学或完整结构定义已解析。只有满足正式项目结论条件时才升级为确证。</p></div><div class="report-formal-status"><b>解释信息</b><p>分类、生源、规则证据和推测碎裂途径逐项展示；无法推理时明确说明。它们不直接改变化合物身份排序。</p></div></div>
      </section>
      ${window.ConfidenceGuide
        ? window.ConfidenceGuide.render(bundle.confidenceGuide, {
          mode: 'report', id: 'sample-report-confidence-guide',
          language: window.MSSLanguage?.getLanguage?.() || 'zh',
        })
        : ''}
      ${authenticRankingReportSection(bundle)}
      ${referenceStandardQCReportSection()}
      ${replicateConsistencyReportSection(bundle)}
      ${chromatographySection(bundle)}
      <section class="sample-report-section" id="sample-report-results"><div class="sample-report-section-head report-major-head"><div><span class="report-kicker">03 · IDENTIFICATION</span><h4>去重候选身份与检测证据</h4><p>同一候选身份只列一次；关联的离子特征、保留时间峰区和重复 MS2 作为支持记录合并展示，未解析母离子事件仍单列。</p></div><span>${integer(bundle.groups.length)} 个身份／事件</span></div>
        ${evidenceOverview(bundle)}
        ${expertImpactSection(bundle)}
        <div class="report-raw-review-guide"><b>如何回查原始 mzML</b><p>下表的<strong>最佳 MS2 原始定位</strong>把每个化合物对应的代表查询谱之保留时间、MS2 隔离目标和谱图编号/native ID 放在同一格；可靠关联连续 MS1 特征时，另列真正用于 HRMS 检索的实测 MS1 质心。源文件为 <span class="mono">${esc(task.filename || '未记录')}</span>。人工复核时优先按谱图编号/native ID 精确定位；若查看软件不支持编号检索，再联合使用 RT 与隔离目标查找，并以查询碎片峰核对。这里的 RT 是样本实测 MS2 的采集时间，不是库谱 RT。</p></div>
      </section>
      ${candidateIdentitySummarySection(bundle)}
      ${compoundTable(confirmed, '3.1 已确证存在', '只有满足当前项目正式结论条件的化合物列在这里。', 30)}
      ${compoundTable(high, '3.2 High 单一二维连接检出', '样本信号在单一二维连接结构范围具有 High 级 MS/MS 谱库支持；立体化学／完整结构定义层与标准品确证状态另行报告。', 30)}
      ${compoundTable(isomerHigh, '3.3 High 异构体组检出', 'High 级证据支持分子式异构体组，但不足以稳妥区分具体异构体。', 30)}
      ${compoundTable(isomerNonHigh, '3.4 中等证据异构体组', '实体规则将本组列为 Medium，候选异构体尚未排除；该等级不覆盖原单谱谱库证据。', 30)}
      ${compoundTable(medium, '3.5 中等证据候选', '实体规则将本记录列为 Medium；单谱谱库证据与结构确证分别报告。', 25)}
      ${compoundTable(background, '3.6 空白背景比较提示（默认收起）', '谱库命中存在，但背景比较提示该信号与任务空白相关；不进入主要样品结果。该标记不等于跨针相减或删除候选。', 25, true)}
      ${compoundTable(unresolvedEvents,
        isEnglishInterface() ? '3.7 Unresolved precursor events (collapsed by default)' : '3.7 未解析母离子事件（默认收起）',
        isEnglishInterface()
          ? 'Candidate structures from one precursor event could not be consolidated into a reportable identity: known molecular formulae may conflict, or at least one candidate may lack a usable formula. Complete candidates are retained for review, but the event is excluded from the independent compound count.'
          : '同一母离子事件的候选结构无法归并为一个可报告身份：已知分子式可能互相冲突，或至少一个候选缺少可用分子式。完整候选保留供复核，但该事件不计作独立化合物。',
        25, true)}
      ${compoundTable(low, '3.8 证据不足线索（默认收起）', '保留低证据谱库线索，但不进入样本主要鉴定结果。', 25, true)}
      ${compoundTable(unannotated, '3.9 未注释离子特征', '这些色谱离子特征已经检出，但普通谱库检索没有可用结构候选；原始定位与 Feature ID 保留供后续复核。', 25)}
      ${isomerCandidateAppendix(allIsomerCompetition)}
      ${expectedSection(bundle)}
      ${qualitySection(bundle, summary)}
      ${followupSection()}
      ${methodSection(bundle)}
      ${provenanceSection(bundle)}
    </article>`;
  }

  function render(bundle) {
    const area = document.getElementById('sampleReportArea');
    const toolbar = document.getElementById('sampleReportToolbar');
    if (!area) return;
    area.innerHTML = reportDocumentHtml(bundle);
    bindReportInteractions(area);
    toolbar?.classList.remove('hidden');
  }

  async function fetchReportPayload(taskId, systemSnapshot = null) {
    const [bundle, plans, system, tic, expected, referenceQC] = await Promise.all([
      SampleResultData.load(taskId),
      API.get(`/tasks/${encodeURIComponent(taskId)}/acquisition-plans?limit=25`).catch(() => ({ plans: [], metadata: {} })),
      systemSnapshot
        ? Promise.resolve(systemSnapshot)
        : API.get('/system/info').catch(() => ({})),
      API.get(`/tasks/${encodeURIComponent(taskId)}/tic_overview`).catch(() => null),
      API.get(`/tasks/${encodeURIComponent(taskId)}/expected-compounds/groups?evidence_view=supported&limit=8`).catch(() => null),
      API.get(`/tasks/${encodeURIComponent(taskId)}/reference-standard-qc`).catch(() => null),
    ]);
    let blankExcludedRows = [];
    if (Number(bundle?.task?.blank_summary?.excluded_spectra || 0) > 0) {
      if (!blankExcludedRowsCache.has(taskId)) {
        blankExcludedRowsCache.set(taskId, API.get(
          `/tasks/${encodeURIComponent(taskId)}/results?projection=summary&blank_view=flagged&limit=5000`
        ).then(payload => (
          Array.isArray(payload?.results)
            ? payload.results.filter(row => row?.blank_excluded === true)
            : []
        )).catch(() => []));
      }
      blankExcludedRows = await blankExcludedRowsCache.get(taskId);
    }
    return {
      taskId, bundle, plans, system, tic, expected, referenceQC,
      blankExcludedRows,
    };
  }

  function activateReportPayload(payload) {
    currentTaskId = payload.taskId;
    currentBundle = payload.bundle;
    currentPlans = payload.plans || { plans: [], metadata: {} };
    currentSystem = payload.system || {};
    currentTic = payload.tic || null;
    currentBlankExcludedRows = payload.blankExcludedRows || [];
    currentExpected = payload.expected || null;
    currentReferenceQC = payload.referenceQC || null;
  }

  async function load(taskId) {
    const area = document.getElementById('sampleReportArea');
    const toolbar = document.getElementById('sampleReportToolbar');
    currentTaskId = taskId || null;
    if (!taskId) {
      currentBundle = null;
      currentPlans = { plans: [], metadata: {} };
      currentTic = null;
      currentBlankExcludedRows = [];
      currentExpected = null;
      currentReferenceQC = null;
      toolbar?.classList.add('hidden');
      if (area) area.innerHTML = '<div class="empty-state"><div class="empty-icon">▤</div><div class="empty-title">选择样本以生成综合报告</div><div class="empty-desc">报告将先给出总体结论，再说明方法、局限和建议</div></div>';
      return;
    }
    currentBundle = null;
    currentPlans = { plans: [], metadata: {} };
    currentTic = null;
    currentBlankExcludedRows = [];
    currentExpected = null;
    currentReferenceQC = null;
    if (area) area.innerHTML = '<div class="sample-loading"><span class="spinner"></span><b>正在形成样本综合报告…</b><small>整理结论、分析方法、局限和确认实验建议</small></div>';
    toolbar?.classList.add('hidden');
    try {
      const task = await API.get(`/tasks/${encodeURIComponent(taskId)}/status`);
      if (currentTaskId !== taskId) return;
      if (task.completed_batch) {
        if (area) area.innerHTML = `${FrozenTaskUI.panel(task)}<iframe title="当前样本全量MS/MS与MS1证据报告" src="${API.BASE}/tasks/${encodeURIComponent(taskId)}/report" style="width:100%;height:80vh;min-height:640px;border:1px solid #dce4ee;border-radius:10px"></iframe>`;
        return;
      }
      const payload = await fetchReportPayload(taskId);
      if (currentTaskId !== taskId) return;
      activateReportPayload(payload);
      render(payload.bundle);
    } catch (error) {
      if (area) area.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">报告生成失败</div><div class="empty-desc">${esc(error.message || error)}</div></div>`;
    }
  }

  async function open(taskId) {
    App.navigateTo('report');
    await refreshSelect(taskId);
    const select = document.getElementById('reportTaskSelect');
    if (!select) return;
    select.value = taskId;
    await load(taskId);
  }

  function printableHTML(bundle = currentBundle) {
    const documentHtml = bundle ? reportDocumentHtml(bundle) : '';
    const title = bundle?.task?.filename || currentTaskId || 'sample';
    return `<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MSsystem 样本综合报告 · ${esc(title)}</title><style>${standaloneStylesV2()}</style></head><body>${documentHtml}<footer>MSsystem · 本报告读取保存的分析结果，不重新评分或改变候选排序。</footer><script>${standaloneInteractionScript()}</script></body></html>`;
  }

  function legacyStandaloneStyles() {
    return `*{box-sizing:border-box}body{margin:0;background:#edf2f4;color:#17252a;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;line-height:1.65}.sample-report-document{max-width:1120px;margin:24px auto;background:#fff;border:1px solid #d8e1e5}.sample-report-cover{padding:36px 42px;background:#123f43;color:#fff;display:flex;justify-content:space-between;gap:24px}.sample-report-cover span{font-size:12px;letter-spacing:.08em;opacity:.8}.sample-report-cover h3{font-size:28px;margin:6px 0}.sample-report-cover p{margin:0}.sample-report-cover-meta{display:grid;gap:7px;align-content:center}.sample-report-section{padding:28px 42px;border-bottom:1px solid #e1e8eb}.sample-report-section h4{font-size:20px;margin:0 0 12px}.report-lead{font-size:17px}.sample-report-kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:18px 0}.sample-report-kpis div{padding:16px;border:1px solid #cfe0df;background:#f2f8f7}.sample-report-kpis strong{display:block;font-size:25px;color:#176b68}.sample-report-kpis span{font-size:12px}.report-conclusion-grid,.report-method-grid,.report-quality-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}.report-method-grid{grid-template-columns:repeat(2,1fr)}.report-conclusion-grid>div,.report-method-grid>div,.report-quality-grid>div{padding:15px;border:1px solid #dde6e8;background:#fbfcfc}.sample-report-section-head{display:flex;justify-content:space-between;gap:20px;align-items:start}.sample-report-section-head>span{background:#e4f2ef;color:#176b68;padding:3px 10px;border-radius:999px;font-size:12px}.report-raw-review-guide{margin-top:14px;padding:12px 14px;border:1px solid #9cc9c6;border-left:4px solid #176b68;background:#f0f8f7}.report-raw-review-guide p{margin:3px 0 0}.report-raw-locator{min-width:170px;line-height:1.4}.report-raw-locator b,.report-raw-locator span,.report-raw-locator small{display:block}.report-raw-locator span{color:#718087;font-size:10px}.table-wrap{overflow-x:auto}table{width:100%;border-collapse:collapse;font-size:12px}th{background:#eaf2f2;text-align:left;padding:8px;border-bottom:2px solid #bed1d1}td{padding:8px;border-bottom:1px solid #e4eaec;vertical-align:top}small{display:block;color:#63757c}.sample-evidence-badge{display:inline-block;padding:2px 7px;border-radius:999px;font-size:11px;font-weight:700}.evidence-high{background:#dff4e8;color:#17613a}.evidence-isomer{background:#fff0cf;color:#8a5800}.evidence-medium{background:#e1edfb;color:#245e9d}.evidence-low{background:#edf0f2;color:#647077}.report-quality-side div{display:flex;justify-content:space-between;padding:6px 0;border-bottom:1px solid #e5eaec}.sample-report-provenance dl{display:grid;grid-template-columns:repeat(2,1fr);gap:8px}.sample-report-provenance dl div{display:flex;gap:12px}.sample-report-provenance dt{color:#64777e;min-width:110px}.sample-report-provenance dd{margin:0;word-break:break-word}.report-table-note,.report-empty,.report-section-intro{font-size:12px;color:#64777e}footer{text-align:center;padding:18px;font-size:11px;color:#718087}@media(max-width:760px){.sample-report-cover{display:block}.sample-report-cover-meta{margin-top:20px}.sample-report-kpis,.report-conclusion-grid,.report-method-grid,.report-quality-grid,.sample-report-provenance dl{grid-template-columns:1fr}.sample-report-section{padding:22px}}@media print{body{background:#fff}.sample-report-document{margin:0;max-width:none;border:0}.sample-report-section{break-inside:auto}.sample-report-cover,.sample-report-kpis div,.report-conclusion-grid>div,.report-method-grid>div{-webkit-print-color-adjust:exact;print-color-adjust:exact}footer{display:none}}`;
  }

  function standaloneStylesV2() {
    return `
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;background:#edf2f4;color:#17252a;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;font-size:16px;line-height:1.68}
.sample-report-document{max-width:1180px;margin:28px auto;background:#fff;border:1px solid #d5dfe2;box-shadow:0 16px 45px rgba(35,66,70,.12)}
.sample-report-cover{padding:25px 31px;background:#fff;color:#213a40;border-top:5px solid #1d756f;border-bottom:1px solid #dce5e7;display:flex;justify-content:space-between;align-items:center;gap:30px}
.sample-report-cover>div:first-child>span{font-size:11px;font-weight:800;letter-spacing:.12em;color:#376c6a}
.sample-report-cover h3{font-size:26px;line-height:1.25;margin:6px 0 3px;color:#18383e}
.sample-report-cover p{margin:0;font-size:15px;color:#394f56;overflow-wrap:anywhere}
.sample-report-cover em{display:block;margin-top:8px;color:#52666d;font-size:12px;font-style:normal}
.sample-report-cover-meta{display:grid;min-width:265px;gap:7px;align-content:center;padding-left:22px;border-left:1px solid #dce5e7}
.sample-report-cover-meta span{display:flex;justify-content:space-between;gap:24px;color:#566970;font-size:12px}
.sample-report-cover-meta b{color:#263f45;text-align:right}
.sample-report-toc{display:flex;flex-wrap:wrap;align-items:center;gap:8px;padding:15px 46px;border-bottom:1px solid #dfe8ea;background:#f7faf9}
.sample-report-toc b{margin-right:6px;color:#31565a;font-size:13px}
.sample-report-toc a{padding:6px 10px;border:1px solid #d6e3e2;border-radius:999px;background:#fff;color:#315d60;font-size:12px;text-decoration:none}
.sample-report-section{position:relative;padding:34px 46px;border-bottom:1px solid #e1e8eb}
.sample-report-section:last-child{border-bottom:0}
.sample-report-section h4{margin:2px 0 7px;color:#173f43;font-size:23px;line-height:1.35}
.report-kicker{display:block;color:#1a7771;font-size:11px;font-weight:800;letter-spacing:.11em;text-transform:uppercase}
.sample-report-section-head{display:flex;justify-content:space-between;gap:26px;align-items:flex-start}
.sample-report-section-head>div>p{margin:0;color:#65777c;font-size:14px}
.sample-report-section-head>span{flex:0 0 auto;padding:5px 10px;border-radius:999px;background:#e4f2ef;color:#176b68;font-size:12px;font-weight:800}
.report-status-pill{background:#dff4e8!important;color:#17613a!important}
.report-lead{max-width:1050px;margin:18px 0;color:#2d4d51;font-size:18px;line-height:1.72}
.sample-report-kpis{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin:20px 0}
.sample-report-kpis div,.report-targeted-kpis div{padding:17px 18px;border:1px solid #d4e4e2;border-radius:10px;background:#f3f9f8}
.sample-report-kpis strong,.report-targeted-kpis b{display:block;color:#176b68;font-size:27px;line-height:1.2}
.sample-report-kpis span,.report-targeted-kpis span{display:block;margin-top:5px;color:#5e7378;font-size:12px}
.report-top-candidates{margin:18px 0;padding:17px 19px;border:1px solid #dce7e8;border-radius:10px;background:#fbfcfc}
.report-top-candidates>b{color:#254b4f;font-size:14px}
.report-top-candidates>div{display:flex;flex-wrap:wrap;gap:7px;margin:10px 0}
.report-top-candidates>div span{padding:5px 9px;border-radius:6px;background:#e9f4f2;color:#185e5a;font-size:12px;font-weight:700}
.report-top-candidates small{color:#667a80;font-size:12px}
.report-conclusion-grid,.report-method-grid,.report-quality-grid,.report-evidence-overview{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}
.report-conclusion-grid>div,.report-method-grid>div,.report-quality-main,.report-quality-side,.report-overview-card{padding:18px 20px;border:1px solid #dfe7e9;border-radius:10px;background:#fbfcfc}
.report-conclusion-grid b,.report-method-grid b,.report-quality-grid b{color:#23494c;font-size:14px}
.report-conclusion-grid p,.report-method-grid p{margin:6px 0 0;color:#5a6f74;font-size:14px;line-height:1.68}
.report-formal-status{border-left:4px solid #c28a2d!important;background:#fffbf3!important}
.report-tic-card{margin-top:20px;padding:18px 20px 12px;border:1px solid #d9e5e7;border-radius:12px;background:#fff;box-shadow:0 5px 18px rgba(40,77,82,.07)}
.report-bpc-card{border-color:#ded9ed}.report-bpc-card .report-tic-coverage{background:#f7f5fb}
.report-chart-title{display:flex;justify-content:space-between;gap:20px;align-items:flex-start}
.report-chart-title b{display:block;color:#173f43;font-size:17px}
.report-chart-title span,.report-chart-title small{color:#6b7d82;font-size:12px}
.report-tic-frame{margin-top:7px}.report-tic-svg{display:block;width:100%;height:auto;min-height:290px}
.tic-grid{stroke:#e6edef;stroke-width:1}.tic-axis{stroke:#6d8186;stroke-width:1.2}.tic-main-line{fill:none;stroke:#176b68;stroke-width:2.2;stroke-linejoin:round;stroke-linecap:round}.bpc-main-line{fill:none;stroke:#6f5aa8;stroke-width:2.2;stroke-linejoin:round;stroke-linecap:round}.tic-ms2-rug{stroke:#93a5aa;stroke-width:1;opacity:.55}.tic-high-marker{fill:#df5a43;stroke:#fff;stroke-width:2}.tic-isomer-marker{fill:#d58a12;stroke:#fff;stroke-width:2}.tic-medium-marker{fill:#3977b8;stroke:#fff;stroke-width:1.7}.tic-unassigned-marker{fill:rgba(255,255,255,.92);stroke:#718287;stroke-width:2;stroke-dasharray:2 2}.tic-unassigned-label{fill:#617479;font-size:9px;font-weight:800;pointer-events:none}.tic-axis-label{fill:#718287;font-size:12px}.tic-axis-title{fill:#4f656a;font-size:13px;font-weight:700}
.report-chart-legend{display:flex;flex-wrap:wrap;gap:18px;padding:3px 8px 8px;color:#5f7277;font-size:12px}.report-chart-legend span{display:flex;align-items:center;gap:7px}.report-chart-legend i{display:inline-block;width:18px;height:3px;border-radius:2px}.legend-tic{background:#176b68}.legend-bpc{background:#6f5aa8}.legend-ms2{background:#93a5aa}.legend-high{width:9px!important;height:9px!important;border-radius:50%!important;background:#df5a43}.legend-isomer{width:9px!important;height:9px!important;border-radius:1px!important;background:#d58a12;transform:rotate(45deg)}.legend-medium{width:9px!important;height:9px!important;border-radius:1px!important;background:#3977b8}.legend-unassigned{width:10px!important;height:10px!important;border:2px dashed #718287;border-radius:50%!important;background:#fff}
.report-tic-coverage{display:grid;grid-template-columns:auto minmax(0,1fr);gap:3px 13px;margin:4px 8px 6px;padding:12px 14px;border-radius:8px;background:#f4f8f8;color:#536b70}.report-tic-coverage b{grid-row:1/span 2;color:#234d50;font-size:13px}.report-tic-coverage span{font-size:12px;line-height:1.5}.report-tic-coverage small{color:#718287;font-size:11px;line-height:1.45}
.report-chart-empty{padding:42px;text-align:center;color:#7a898d}
.evidence-unresolved_event{background:#ffedd5;color:#9a3412}
.report-qc-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin-top:16px}
.report-qc-grid>div{padding:16px 17px;border:1px solid #dee7e9;border-radius:9px;background:#f9fbfb}.report-qc-grid span{display:block;color:#63777c;font-size:12px}.report-qc-grid b{display:block;margin:4px 0;color:#183f43;font-size:20px}.report-qc-grid b i{color:#9aabad;font-style:normal}.report-qc-grid small{color:#687b80;font-size:12px}
.report-evidence-overview{margin:20px 0}.report-overview-card-head{display:flex;justify-content:space-between;gap:12px;margin-bottom:12px}.report-overview-card-head b{font-size:15px;color:#23494c}.report-overview-card-head span{color:#718287;font-size:11px;text-align:right}
.report-bar-row{display:grid;grid-template-columns:minmax(120px,1fr) minmax(140px,2fr) 42px;gap:10px;align-items:center;margin:10px 0;color:#536a70;font-size:13px}.report-bar-row>div{height:9px;overflow:hidden;border-radius:999px;background:#edf1f2}.report-bar-row i{display:block;height:100%;border-radius:999px}.report-bar-row>b{text-align:right;color:#284b4f}
.report-expert-impact{margin:18px 0;padding:18px 20px;border:1px solid #d7e1ee;border-left:4px solid #496f9b;border-radius:10px;background:#f7f9fc}.report-expert-impact-head{display:flex;justify-content:space-between;gap:18px;align-items:flex-start}.report-expert-impact-head b{display:block;color:#294d70;font-size:15px}.report-expert-impact-head span{display:block;margin-top:2px;color:#708090;font-size:12px}.report-expert-impact-head>strong{padding:4px 9px;border-radius:999px;background:#e2ebf6;color:#315f8c;font-size:12px}.report-expert-kpis{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:9px;margin:14px 0}.report-expert-kpis>div{padding:11px 12px;border:1px solid #dde5ef;border-radius:8px;background:#fff}.report-expert-kpis b{display:block;color:#2c587f;font-size:19px}.report-expert-kpis span{display:block;color:#6d7e8d;font-size:10px;line-height:1.35}.report-expert-impact p{margin:8px 0;color:#4f6577;font-size:13px}.report-expert-impact>small{color:#6d7d8b;font-size:11px}
.report-raw-review-guide,.report-boundary-note{margin-top:18px;padding:15px 17px;border:1px solid #9cc9c6;border-left:4px solid #176b68;border-radius:8px;background:#f0f8f7;color:#214c50}.report-raw-review-guide p{margin:5px 0 0;font-size:14px}.report-boundary-note{display:flex;gap:15px;margin:16px 0;font-size:13px}.report-boundary-note b{flex:0 0 auto}
.table-wrap{overflow-x:auto;margin-top:14px}table{width:100%;border-collapse:collapse}th{padding:11px 10px;border-bottom:2px solid #b9cecf;background:#eaf2f2;color:#3f575c;text-align:left;font-size:13px}td{padding:12px 10px;border-bottom:1px solid #e2eaec;color:#334e52;vertical-align:top;font-size:14px;line-height:1.5}.sample-report-table{min-width:900px}.report-identification-table,.report-isomer-candidate-table{width:100%;min-width:100%;table-layout:fixed}.report-identification-table th,.report-identification-table td,.report-isomer-candidate-table th,.report-isomer-candidate-table td{overflow-wrap:anywhere}.report-identification-table .report-candidate-cell,.report-identification-table .report-raw-locator,.report-identification-table .report-score-cell{min-width:0}.report-identification-conclusion .report-conclusion-note{margin-top:7px;color:#50676d}.report-support-location .report-support-count{margin-bottom:7px;padding-bottom:6px;border-bottom:1px solid #dfe8e9;text-align:left}.report-support-location .report-support-count>b{font-size:13px}.report-support-location .report-raw-locator{min-width:0}.sample-report-table tbody tr:nth-child(even){background:#fbfcfc}.sample-report-table tbody tr:hover{background:#f0f7f6}.sample-report-table small,.report-raw-locator span,.report-raw-locator small{display:block;margin-top:3px;color:#708187;font-size:12px}.report-candidate-cell{min-width:230px}.report-candidate-cell>b{display:block;color:#183f43;font-size:15px}.report-support-count{text-align:center}.report-support-count>b{font-size:18px;color:#176b68}.report-raw-locator{min-width:175px}.report-raw-locator>b{color:#173f43}.report-score-cell{min-width:145px}.report-score-cell>b{color:#214f53}.sample-evidence-badge{display:inline-block;padding:3px 8px;border-radius:999px;font-size:11px;font-weight:800}.evidence-confirmed{background:#ccebdd;color:#0b573e}.evidence-high{background:#dff4e8;color:#17613a}.evidence-isomer{background:#fff0cf;color:#8a5800}.evidence-medium{background:#e1edfb;color:#245e9d}.evidence-background{background:#f7e7d5;color:#8a4b18}.evidence-low{background:#edf0f2;color:#647077}.evidence-targeted{background:#ede8fb;color:#64439b}.report-inline-warning{color:#a14d18!important;font-weight:700}
.report-table-note,.report-empty,.report-section-intro{color:#667a80;font-size:13px}.sample-report-details>summary:not(.report-section-summary){display:flex;justify-content:space-between;gap:18px;padding:4px 0;cursor:pointer;color:#23494c;font-size:17px}.sample-report-details>summary b{color:#687b80;font-size:12px}.sample-report-details-body{padding-top:18px}.sample-report-details-body>.sample-report-section-head{margin-bottom:0}
.report-targeted-kpis{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:10px;margin:18px 0}.report-targeted-kpis b{font-size:24px}.report-targeted-table{width:100%;min-width:100%;table-layout:fixed}.report-targeted-table th,.report-targeted-table td{overflow-wrap:anywhere;line-height:1.45}.report-targeted-table .report-candidate-cell{min-width:0}.report-targeted-table .report-candidate-cell>b{font-size:14px;line-height:1.35}
.report-quality-main ul{margin:10px 0 0 20px;padding:0;color:#536a70;font-size:14px}.report-quality-main li{margin:8px 0}.report-quality-side>div{display:flex;justify-content:space-between;gap:15px;padding:9px 0;border-bottom:1px solid #e3eaec;color:#60757a;font-size:13px}.report-quality-side>div:last-child{border-bottom:0}.report-quality-side strong{color:#a14d18}
.sample-report-provenance dl{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:5px 22px;margin:15px 0 0}.sample-report-provenance dl div{display:flex;gap:14px;padding:9px 0;border-bottom:1px solid #e8edef;font-size:13px}.sample-report-provenance dt{flex:0 0 115px;color:#718287}.sample-report-provenance dd{margin:0;color:#334e52;overflow-wrap:anywhere}.mono{font-family:"SFMono-Regular",Consolas,monospace}
.sample-report-view-tools{display:flex;justify-content:space-between;align-items:center;gap:16px;padding:14px 46px;border-bottom:1px solid #dce7e7;background:#f2f8f7}.sample-report-view-tools>div:first-child{display:grid;min-width:0;gap:2px}.sample-report-view-tools>div:first-child b{color:#214b4e;font-size:13px}.sample-report-view-tools>div:first-child span{color:#65797e;font-size:12px;line-height:1.45}.report-view-actions{display:flex;flex:0 0 auto;flex-wrap:wrap;justify-content:flex-end;align-items:center;gap:7px}.report-view-actions button{padding:7px 11px;color:#28585a;border:1px solid #b9d2d0;border-radius:7px;background:#fff;font:inherit;font-size:12px;font-weight:750;cursor:pointer}.report-view-actions button:hover{border-color:#568f8b;background:#e8f5f3;color:#124f4c}.report-view-actions button:focus-visible,details[data-report-disclosure]>summary:focus-visible{outline:3px solid rgba(46,133,126,.24);outline-offset:2px}.report-view-actions [data-report-disclosure-status]{min-width:112px;color:#61767b;font-size:11px;text-align:right}
.report-tier-section{padding:0}.report-tier-details>summary,.report-compound-card>summary,.report-isomer-group>summary{list-style:none}.report-tier-details>summary::-webkit-details-marker,.report-compound-card>summary::-webkit-details-marker,.report-isomer-group>summary::-webkit-details-marker{display:none}.report-section-summary{display:grid;grid-template-columns:minmax(0,1fr) max-content;justify-content:space-between;align-items:center;gap:20px;width:100%;max-width:100%;min-height:82px;padding:17px 46px;box-sizing:border-box;cursor:pointer;background:#fff;transition:background .15s ease}.sample-report-secondary .report-section-summary{background:#fafcfc}.report-section-summary:hover{background:#f2f8f7}.report-section-summary>span:first-child{display:grid;min-width:0;gap:4px}.report-section-summary>span:first-child>b{color:#173f43;font-size:18px;line-height:1.35}.report-section-summary>span:first-child>small{max-width:850px;color:#667a80;font-size:12px;line-height:1.5;overflow-wrap:anywhere}.report-section-summary-meta{display:grid;min-width:104px;gap:2px;padding-left:7px;text-align:right}.report-section-summary-meta strong{color:#1b6965;font-size:13px;white-space:nowrap}.report-section-summary-meta em{color:#6b7e83;font-size:11px;font-style:normal;font-weight:700;white-space:nowrap}details[open]>summary .when-closed,details:not([open])>summary .when-open{display:none}.report-tier-details[open]>.report-section-summary{border-bottom:1px solid #dfe8e9;background:#f8fbfa}.report-tier-details>.sample-report-details-body{padding:16px 46px 24px;background:#f7faf9}
.report-compound-list,.report-isomer-group-list{display:grid;gap:10px}.report-compound-card,.report-isomer-group{overflow:hidden;border:1px solid #d9e4e5;border-radius:10px;background:#fff;transition:border-color .15s ease,box-shadow .15s ease}.report-compound-card[open],.report-isomer-group[open]{border-color:#a9cbc8;box-shadow:0 5px 16px rgba(37,78,81,.08)}.report-compound-card>summary{display:grid;grid-template-columns:40px minmax(190px,2fr) minmax(145px,1fr) minmax(135px,.8fr) 86px;gap:14px;align-items:start;min-height:72px;padding:12px 14px;cursor:pointer;transition:background .15s ease}.report-compound-card>summary>*{min-width:0}.report-compound-card>summary:hover{background:#f2f8f7}.report-compound-card[open]>summary{background:#edf7f5}.report-compound-number{display:inline-grid;width:32px;height:32px;place-items:center;color:#517176;border:1px solid #d1e0e0;border-radius:50%;background:#f7fafa;font-size:11px;font-weight:800}.report-compound-identity,.report-compound-scope,.report-compound-metrics{display:grid;min-width:0;gap:3px}.report-compound-identity>b{display:block;color:#173b40;font-size:14px;line-height:1.55;white-space:normal;overflow-wrap:anywhere;word-break:normal}.report-compound-card summary small,.report-isomer-group summary small{display:block;color:#6b7f84;font-size:11px;line-height:1.55;white-space:normal;overflow-wrap:anywhere;word-break:normal}.report-compound-scope{justify-items:start}.report-compound-metrics>b{color:#245b5d;font-size:12px;line-height:1.5}.report-disclosure-action{color:#276d69;font-size:11px;font-weight:800;line-height:1.5;text-align:right;white-space:normal}.report-disclosure-action:after{content:"⌄";display:inline-block;margin-left:5px;font-size:13px;transition:transform .15s ease}.report-compound-card[open] .report-disclosure-action:after{transform:rotate(180deg)}.report-compound-body{padding:16px;border-top:1px solid #dce7e7;background:#fbfcfc}.report-compound-evidence-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:11px}.report-compound-evidence-grid>div{min-width:0;padding:13px 14px;border:1px solid #e0e8ea;border-radius:8px;background:#fff}.report-compound-evidence-grid>div>span:first-child{display:block;margin-bottom:6px;color:#6c8085;font-size:10px;font-weight:800;letter-spacing:.04em}.report-compound-evidence-grid>div>b,.report-compound-evidence-grid>div>div>b{color:#244c50;font-size:12px}.report-compound-evidence-grid p{margin:5px 0 0;color:#50686d;font-size:12px;line-height:1.55}.report-compound-evidence-grid small{display:block;margin-top:3px;color:#6d8085;font-size:11px;line-height:1.45;overflow-wrap:anywhere}.report-compound-evidence-grid .report-support-count{margin:0 0 6px;text-align:left}.report-compound-evidence-grid .report-raw-locator{min-width:0}
.report-isomer-inline-head{display:flex;justify-content:space-between;align-items:flex-start;gap:16px;margin-top:14px;padding:12px 14px;border:1px solid #ecd7a4;border-radius:8px 8px 0 0;background:#fff9ea}.report-isomer-inline-head>div>b{display:block;color:#765119;font-size:13px}.report-isomer-inline-head small{display:block;margin-top:3px;color:#7c6b4c;font-size:11px}.report-isomer-inline-head>span{flex:0 0 auto;color:#8a641f;font-size:11px;font-weight:800}.report-isomer-table-wrap{margin-top:0;border:1px solid #e5dfce;border-top:0;border-radius:0 0 8px 8px;background:#fff}.report-isomer-candidate-table{min-width:920px;table-layout:fixed}.report-isomer-candidate-table th,.report-isomer-candidate-table td{overflow-wrap:anywhere}.report-isomer-identity>b{color:#183f43}.report-isomer-lead{color:#8a5800}.report-isomer-group>summary{display:grid;grid-template-columns:40px minmax(0,1fr) auto;gap:12px;align-items:center;padding:12px 14px;cursor:pointer}.report-isomer-group>summary:hover,.report-isomer-group[open]>summary{background:#fff9ec}.report-isomer-group-body{padding:0 14px 14px;border-top:1px solid #ece5d6;background:#fffdf8}.report-isomer-group-body .report-isomer-table-wrap{margin-top:12px}.report-isomer-appendix>.sample-report-details-body>.report-section-intro{margin:0 0 14px}.report-candidate-identity-table{width:100%;min-width:1040px;table-layout:fixed}.report-candidate-identity-table th,.report-candidate-identity-table td{line-height:1.55;vertical-align:top;white-space:normal;overflow-wrap:anywhere;word-break:normal}.report-candidate-identity-table b,.report-candidate-identity-table small{line-height:1.55}.report-curated-natural{color:#27643d!important;font-weight:800}
footer{text-align:center;padding:22px;color:#718087;font-size:12px}
@media(max-width:980px){.sample-report-view-tools{align-items:flex-start;flex-direction:column}.report-view-actions{justify-content:flex-start}.report-compound-card>summary{grid-template-columns:40px minmax(0,1fr) auto}.report-compound-scope{grid-column:2}.report-compound-metrics{grid-column:3;grid-row:1}.report-disclosure-action{grid-column:3;grid-row:2}}
@media(max-width:850px){body{font-size:15px}.sample-report-cover{display:block;padding:30px}.sample-report-cover-meta{margin-top:24px;min-width:0}.sample-report-toc,.sample-report-section{padding-left:24px;padding-right:24px}.sample-report-view-tools,.report-section-summary,.report-tier-details>.sample-report-details-body{padding-left:24px;padding-right:24px}.sample-report-kpis,.report-targeted-kpis,.report-expert-kpis{grid-template-columns:repeat(2,1fr)}.report-qc-grid,.report-conclusion-grid,.report-method-grid,.report-quality-grid,.report-evidence-overview,.sample-report-provenance dl,.report-compound-evidence-grid{grid-template-columns:1fr}.report-tic-svg{min-height:230px}.report-tic-coverage{grid-template-columns:1fr}.report-tic-coverage b{grid-row:auto}}
@media(max-width:650px){.report-compound-card>summary{grid-template-columns:36px minmax(0,1fr);gap:9px}.report-compound-scope,.report-compound-metrics,.report-disclosure-action{grid-column:2;grid-row:auto}.report-isomer-group>summary{grid-template-columns:36px minmax(0,1fr)}.report-isomer-group .report-section-summary-meta{grid-column:2;text-align:left}}
@media(max-width:520px){.report-section-summary{grid-template-columns:minmax(0,1fr);gap:10px}.report-section-summary-meta{display:flex;min-width:0;gap:10px;padding-left:0;text-align:left}}
@media print{@page{size:A4 landscape;margin:10mm}body{background:#fff;font-size:11pt;line-height:1.5}.sample-report-document{margin:0;max-width:none;border:0;box-shadow:none}.sample-report-cover{padding:24px 28px}.sample-report-cover h3{font-size:25px}.sample-report-toc,.sample-report-view-tools{display:none!important}.sample-report-section{padding:22px 28px}.report-tier-section{padding:0}.sample-report-section h4{font-size:18px}.report-lead{font-size:13px}.sample-report-kpis div,.report-qc-grid>div,.report-overview-card,.report-method-grid>div,.report-quality-main,.report-quality-side{break-inside:avoid}.report-tic-card,.sample-report-targeted,.report-compound-card{break-inside:avoid}.report-tic-svg{min-height:220px}.sample-report-details:not([open])>.sample-report-details-body,details[data-report-disclosure]:not([open])>.report-compound-body,details[data-report-disclosure]:not([open])>.report-isomer-group-body{display:block!important}.sample-report-details>summary,details[data-report-disclosure]>summary{list-style:none}.report-disclosure-action,.report-section-summary-meta em{display:none!important}.sample-report-table{font-size:9px}th{font-size:9px;padding:7px}td{font-size:9px;padding:7px}.sample-report-table small,.report-raw-locator span,.report-raw-locator small{font-size:8px}a{color:inherit;text-decoration:none}footer{display:none}}
`;
  }

  function csvCell(value) {
    const text = String(value == null ? '' : value);
    return /[",\n\r]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  }

  function downloadIsomerCandidateCsv() {
    if (!currentBundle) {
      Toast.show('请先选择样本', 'warning');
      return;
    }
    const groups = (currentBundle.groups || []).filter(
      groupHasIsomerCompetition
    );
    const header = [
      '代表谱图编号', '前体 m/z', 'RT (min)', '鉴定层级',
      '异构体组置信度', '异构体解析状态', '异构体歧义类型',
      '组内候选总数', '候选记录状态', '候选排名', '候选角色',
      '是否当前最佳候选', '候选化合物名称', '分子式', 'InChIKey',
      'SMILES', '候选谱库加合物', 'MS/MS 相似度', '谱熵相似度',
      '反向匹配覆盖', '质量误差 (ppm)', '综合匹配分', '候选排序分',
      '候选置信度', '化学类别', '与最佳候选关系', '谱库来源',
      '馆藏谱图编号',
    ];
    const rows = [];
    groups.forEach(group => {
      const representative = group.rep || {};
      const match = representative.match || {};
      const members = candidateMembersForGroup(group);
      const leadName = String(group.leadCandidate?.compound_name
        || match.compound_name || members[0]?.compound_name || '');
      const candidateCount = Math.max(
        Number(group.candidateMemberCount || 0),
        Number(match.isomer_group_size || 0),
        members.length
      );
      members.forEach((member, index) => {
        const isLead = member.is_ranking_lead === true
          || String(member.compound_name || '') === leadName
          || index === 0;
        rows.push([
          representative.spectrum_id || '',
          representative.precursor_mz ?? '',
          retentionTimeMinutes(representative) ?? '',
          match.identification_scope || 'isomer_group',
          match.isomer_group_confidence || '',
          match.isomer_resolution_status || '',
          match.isomer_ambiguity_type || '',
          candidateCount,
          member.name_only_fallback
            ? '旧任务仅保存候选名称'
            : 'available',
          member.rank || member.display_order || index + 1,
          isLead
            ? '当前最佳候选（排名首位，未代表唯一解析）'
            : '竞争异构体候选',
          isLead ? '是' : '否',
          member.compound_name || '',
          member.formula || group.compound?.formula || match.formula || '',
          member.inchikey || '',
          member.smiles || '',
          member.adduct || '',
          member.cosine ?? '',
          member.entropy ?? '',
          member.reverse_match ?? '',
          member.mass_ppm ?? '',
          member.final_score ?? '',
          member.rank_score ?? '',
          member.confidence || '',
          member.classification || '',
          member.relationship || '',
          member.dataset_source || '',
          member.lib_spectrum_id || '',
        ]);
      });
    });
    const csv = '\ufeff' + [header, ...rows]
      .map(values => values.map(csvCell).join(','))
      .join('\r\n') + '\r\n';
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    const stem = safeReportStem(
      currentBundle.task?.filename || 'sample', currentTaskId
    );
    anchor.href = url;
    anchor.download = `${stem}_isomer_competing_candidates.csv`;
    anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    Toast.show(
      `已导出 ${integer(groups.length)} 个异构体组、${integer(rows.length)} 个逐候选条目`,
      'success',
      7000
    );
  }

  function downloadReviewIndex() {
    if (!currentBundle) {
      Toast.show('请先选择样本', 'warning');
      return;
    }
    const groupByResultIndex = new Map(
      (currentBundle.featureGroups || currentBundle.groups).map(group => [
        String(group.rep?.result_index ?? ''), group,
      ])
    );
    const rows = [...currentBundle.matched].sort((a, b) => {
      const art = retentionTimeMinutes(a);
      const brt = retentionTimeMinutes(b);
      return (art == null ? Infinity : art) - (brt == null ? Infinity : brt)
        || Number(a.precursor_mz || 0) - Number(b.precursor_mz || 0)
        || String(a.spectrum_id || '').localeCompare(String(b.spectrum_id || ''));
    });
    const header = [
      '源文件', 'Feature_ID', '结果记录_ID', '化合物候选', '分子式', '证据状态', 'RT_min', 'RT_mm_ss',
      '异构体组候选总数', '当前最佳候选', '其余异构体竞争候选', '异构体候选完整明细_JSON',
      'MS2_selected_ion_mz', 'Search_precursor_mz', 'Search_precursor_source',
      'Selected_to_search_delta_ppm', 'MS2_spectrum_native_ID', '关联前级扫描', '结果索引',
      '碰撞能量', '碰撞能量单位', '可用碎片数', '展示类别',
      'Cosine', 'Final_rank_score', 'Legacy_rule_aggregate_uncalibrated_audit', 'Rule_aggregate_source', 'Legacy_rule_weighted_term_audit',
    ];
    const table = [header, ...rows.map(row => {
      const match = row.match || {};
      const group = groupByResultIndex.get(String(row.result_index ?? ''));
      const locator = rawSpectrumLocator(row);
      const candidateMembers = Array.isArray(group?.candidateMembers)
        ? group.candidateMembers : [];
      const leadCandidate = candidateMembers.find(member => member.is_ranking_lead)
        || group?.leadCandidate || candidateMembers[0] || {};
      const competingCandidates = candidateMembers.filter(member => (
        !member.is_ranking_lead
        && String(member.compound_name || '') !== String(leadCandidate.compound_name || '')
      ));
      const candidateDetails = candidateMembers.map(member => ({
        rank: member.rank ?? member.display_order ?? null,
        role: member.is_ranking_lead ? 'ranking_lead' : 'competitive_isomer',
        compound_name: member.compound_name || '',
        formula: member.formula || '',
        inchikey: member.inchikey || '',
        smiles: member.smiles || '',
        adduct: member.adduct || '',
        cosine: member.cosine ?? null,
        entropy: member.entropy ?? null,
        final_score: member.final_score ?? null,
        rank_score: member.rank_score ?? null,
        relationship: member.relationship || '',
        dataset_source: member.dataset_source || '',
        lib_spectrum_id: member.lib_spectrum_id || '',
      }));
      return [
        currentBundle.task?.filename || '',
        group?.compound?.feature_id || row.lcms_feature_id
          || row.precursor_feature_id || '',
        group?.compound?.result_id || '',
        displayName(group || { rep: row }),
        match.formula || '',
        group?.status?.label || '',
        locator.rt == null ? '' : locator.rt.toFixed(6),
        locator.clock === '—' ? '' : locator.clock,
        groupHasIsomerCompetition(group)
          ? Math.max(Number(group?.candidateMemberCount || 0), candidateMembers.length)
          : 0,
        groupHasIsomerCompetition(group) ? (leadCandidate.compound_name || match.compound_name || '') : '',
        groupHasIsomerCompetition(group)
          ? competingCandidates.map(member => member.compound_name).filter(Boolean).join('; ')
          : '',
        groupHasIsomerCompetition(group) ? JSON.stringify(candidateDetails) : '',
        Number.isFinite(locator.selectedIonMz) ? locator.selectedIonMz.toFixed(8) : '',
        Number.isFinite(locator.searchPrecursorMz) ? locator.searchPrecursorMz.toFixed(8) : '',
        locator.searchSource,
        Number.isFinite(locator.searchDeltaPpm) ? locator.searchDeltaPpm.toFixed(6) : '',
        locator.spectrumId,
        locator.precursorScanRef,
        row.result_index ?? '',
        row.collision_energy ?? '',
        row.collision_energy_unit || '',
        row.peaks_count ?? '',
        displayClassification(match),
        match.cosine ?? '',
        match.final_score ?? '',
        match.expert_score ?? '',
        match.expert_score_source || '',
        expertWeightedTerm(match).toFixed(6),
      ];
    })];
    const csv = '\ufeff' + table.map(values => values.map(csvCell).join(',')).join('\r\n');
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    const stem = String(currentBundle.task?.filename || 'sample').replace(/\.[^.]+$/, '').replace(/[^\w\u4e00-\u9fff-]+/g, '_');
    anchor.href = url;
    anchor.download = `${stem}_MS2人工复核索引.csv`;
    anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 60000);
    Toast.show(`已生成 ${integer(rows.length)} 张 MS2 的人工复核索引`, 'success');
  }

  async function generate(taskId = null) {
    const target = taskId || currentTaskId || document.getElementById('reportTaskSelect')?.value;
    if (!target) {
      Toast.show('请先选择样本', 'warning');
      return;
    }
    const task = await API.get(`/tasks/${encodeURIComponent(target)}/status`);
    if (task.completed_batch) {
      FrozenTaskUI.download(`${API.BASE}/tasks/${encodeURIComponent(target)}/download/html`);
      return;
    }
    Toast.show('正在整理样本综合报告，请稍候…', 'info', 6000);
    if (!currentBundle || currentTaskId !== target) await load(target);
    if (!currentBundle) return;
    const blob = new Blob([printableHTML()], { type: 'text/html;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    const stem = String(currentBundle.task?.filename || 'sample').replace(/\.[^.]+$/, '').replace(/[^\w\u4e00-\u9fff-]+/g, '_');
    anchor.href = url;
    anchor.download = `${stem}_样本综合报告.html`;
    anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 60000);
    Toast.show('样本综合报告已生成', 'success');
  }

  function currentReportContext() {
    return {
      taskId: currentTaskId,
      bundle: currentBundle,
      plans: currentPlans,
      system: currentSystem,
      tic: currentTic,
      blankExcludedRows: currentBlankExcludedRows,
      expected: currentExpected,
    };
  }

  function htmlForPayload(payload) {
    const previous = currentReportContext();
    activateReportPayload(payload);
    try {
      return printableHTML(payload.bundle);
    } finally {
      activateReportPayload(previous);
    }
  }

  function safeReportStem(filename, taskId) {
    const stem = String(filename || taskId || 'sample')
      .replace(/\.[^.]+$/, '')
      .replace(/[<>:"/\\|?*\u0000-\u001f]+/g, '_')
      .replace(/\s+/g, ' ')
      .replace(/^\.+|\.+$/g, '')
      .trim()
      .slice(0, 140);
    return stem || `sample_${String(taskId || '').slice(0, 8) || 'unknown'}`;
  }

  function uniqueReportArchiveName(filename, taskId, usedNames) {
    const stem = safeReportStem(filename, taskId);
    let candidate = `${stem}_样本综合报告.html`;
    let serial = 2;
    while (usedNames.has(candidate.toLocaleLowerCase())) {
      const suffix = String(taskId || '').slice(0, 8) || serial;
      candidate = `${stem}__${suffix}${serial > 2 ? `_${serial}` : ''}_样本综合报告.html`;
      serial += 1;
    }
    usedNames.add(candidate.toLocaleLowerCase());
    return candidate;
  }

  async function generateBatch(taskIds, onProgress = null) {
    const uniqueTaskIds = [...new Set((taskIds || []).map(value => String(value || '').trim()).filter(Boolean))];
    if (!uniqueTaskIds.length) throw new Error('请至少选择一个已完成样本');

    const sharedSystem = await API.get('/system/info').catch(() => ({}));
    const successful = [];
    const failed = [];
    const usedNames = new Set();
    for (let index = 0; index < uniqueTaskIds.length; index += 1) {
      const taskId = uniqueTaskIds[index];
      try {
        const payload = await fetchReportPayload(taskId, sharedSystem);
        const sourceFilename = payload.bundle?.task?.filename || taskId;
        const archiveName = uniqueReportArchiveName(sourceFilename, taskId, usedNames);
        successful.push({
          taskId,
          sourceFilename,
          archiveName,
          content: htmlForPayload(payload),
        });
      } catch (error) {
        failed.push({ taskId, error: String(error?.message || error || '未知错误') });
      }
      onProgress?.({
        phase: 'report',
        completed: index + 1,
        total: uniqueTaskIds.length,
        exported: successful.length,
        failed: failed.length,
      });
    }
    if (!successful.length) {
      throw new Error(`所选 ${uniqueTaskIds.length} 个样本均未能生成综合报告`);
    }

    const exportedAt = new Date().toISOString();
    const manifestRows = [
      ['Source_File', 'Task_ID', 'Output_File', 'Status', 'Error', 'Exported_At'],
      ...successful.map(item => [
        item.sourceFilename, item.taskId, item.archiveName, 'exported', '', exportedAt,
      ]),
      ...failed.map(item => ['', item.taskId, '', 'failed', item.error, exportedAt]),
    ];
    const manifest = '\ufeff' + manifestRows
      .map(row => row.map(csvCell).join(','))
      .join('\r\n');
    const entries = [
      { name: '00_批量报告清单.csv', content: manifest },
      ...successful.map(item => ({ name: item.archiveName, content: item.content })),
    ];
    const zipBlob = await buildZip(entries, onProgress);
    const now = new Date();
    const stamp = `${now.getFullYear()}${String(now.getMonth() + 1).padStart(2, '0')}${String(now.getDate()).padStart(2, '0')}`;
    const sampleStems = [...new Set(successful.map(item =>
      safeReportStem(item.sourceFilename, item.taskId)
    ))];
    let samplePart = 'MSsystem';
    if (sampleStems.length <= 4 && sampleStems.length) {
      samplePart = sampleStems.join('__');
    } else if (sampleStems.length > 4) {
      samplePart = `${sampleStems[0]}__${sampleStems[1]}__等${sampleStems.length}个样本`;
    }
    samplePart = samplePart.slice(0, 170);
    const countToken = successful.length > 1 ? `_${successful.length}samples` : '';
    const downloadName = `${samplePart}_样本综合报告${countToken}_${stamp}.zip`;
    const url = URL.createObjectURL(zipBlob);
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = downloadName;
    anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    return {
      downloadName,
      exported: successful.length,
      failed: failed.length,
      total: uniqueTaskIds.length,
    };
  }

  function openBatchExport() {
    App.navigateTo('export');
    const selector = document.getElementById('batchExportType');
    if (selector) {
      selector.value = window.MSSLanguage?.getLanguage() === 'en'
        ? 'html_standard|auto'
        : 'html_integrated|zh';
    }
    setTimeout(() => {
      document.getElementById('batchExportCard')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }, 0);
    Toast.show('请勾选多个已完成样本，然后点击“下载所选单层 / 报告”', 'info', 6500);
  }

  function printReport() {
    if (!currentBundle) {
      Toast.show('请先选择样本', 'warning');
      return;
    }
    const popup = window.open('', '_blank');
    if (!popup) {
      Toast.show('浏览器阻止了打印窗口，请允许弹出窗口后重试', 'warning');
      return;
    }
    popup.opener = null;
    popup.document.open();
    popup.document.write(printableHTML());
    popup.document.close();
    popup.addEventListener('load', () => {
      popup.focus();
      popup.print();
    }, { once: true });
  }

  function downloadServerReport(formatType = 'standard') {
    if (!currentTaskId) {
      Toast.show('请先选择样本', 'warning');
      return;
    }
    const language = window.MSSLanguage?.getLanguage() || 'zh';
    return ExportMod.queueSingleExport('html', {format_type:formatType, language}, currentTaskId);
  }

  return {
    init, refreshSelect, load, open, summaryChromatograms,
    generate, generateBatch, openBatchExport, printReport,
    downloadReviewIndex, downloadIsomerCandidateCsv, downloadServerReport,
    buildZipEntries: buildZip,
  };
})();

// Backward-compatible name for legacy callers; now produces the integrated report.
window.SampleReport = SampleReport;
window.ReportGen = SampleReport;
