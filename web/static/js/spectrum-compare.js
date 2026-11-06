/**
 * SpectrumComparator — mirror plot for a query and a saved library spectrum.
 * The component only displays recorded matches and scores; it does not infer them.
 */
class SpectrumComparator {
    constructor() {
        this.containerId = null;
    }

    static savedEvidence(match) {
        const evidence = match && match.spectral_evidence;
        return evidence && evidence.schema === 'mssystem.spectral_evidence.v1'
            && evidence.status === 'captured'
            && Array.isArray(evidence.query_peaks) && Array.isArray(evidence.library_peaks)
            && Array.isArray(evidence.matched_pairs) ? evidence : null;
    }

    /**
     * Render a mirror plot comparing two spectra.
     * @param {string} containerId
     * @param {Object|null} query - {peaks: [[mz, intensity], ...], name}
     * @param {Object|null} library - {peaks: [[mz, intensity], ...], name}
     * @param {Array} matchedPairs - [[query_mz, lib_mz], ...] or
     *   [{query_mz, lib_mz}, ...]
     * @param {Object} options - {height, title, score, scoreLabel, text}
     */
    render(containerId, query, library, matchedPairs = [], options = {}) {
        this.containerId = containerId;
        const container = document.getElementById(containerId);
        if (!container || typeof Plotly === 'undefined') return;

        const settings = options && typeof options === 'object' ? options : {};
        const textOptions = settings.text && typeof settings.text === 'object'
            ? settings.text : {};
        const labels = {
            query: textOptions.query || 'Query',
            library: textOptions.library || 'Library match',
            matchedQuery: textOptions.matchedQuery || 'Query — recorded match',
            otherQuery: textOptions.otherQuery || 'Query — other peak',
            matchedLibrary: textOptions.matchedLibrary || 'Reference — recorded match',
            otherLibrary: textOptions.otherLibrary || 'Reference — other peak',
            locatedPairs: textOptions.locatedPairs || 'Located recorded peak pairs',
            referenceMissing: textOptions.referenceMissing || 'Reference spectrum unavailable',
            empty: textOptions.empty || 'No valid spectrum peaks available',
            mz: textOptions.mz || 'm/z',
            relativeIntensity: textOptions.relativeIntensity || 'Relative intensity (%)',
        };

        const qPeaks = this._validPeaks(query && query.peaks);
        const lPeaks = this._validPeaks(library && library.peaks);
        const qNorm = this._normalize(qPeaks);
        const lNorm = this._normalize(lPeaks);
        const recordedPairs = Array.isArray(matchedPairs) ? matchedPairs : [];
        const locatedPairs = this._locatePairs(recordedPairs, qPeaks, lPeaks);
        const matchedQ = new Set(locatedPairs.map(pair => pair.queryIndex));
        const matchedL = new Set(locatedPairs.map(pair => pair.libraryIndex));

        const queryName = this._escapeHtml(this._displayName(query && query.name, labels.query));
        const libraryName = this._escapeHtml(this._displayName(library && library.name, labels.library));
        const traces = [];
        const addTrace = (peaks, indices, direction, name, color, spectrumLabel) => {
            if (indices.length === 0) return;
            const x = [];
            const y = [];
            const customdata = [];
            const markerSizes = [];
            indices.forEach(index => {
                const peak = peaks[index];
                const plottedIntensity = direction * peak.relativeIntensity;
                const hoverData = [peak.originalMz, peak.relativeIntensity];
                x.push(peak.mz, peak.mz, null);
                y.push(0, plottedIntensity, null);
                customdata.push(hoverData, hoverData, null);
                markerSizes.push(0, 5, 0);
            });
            traces.push({
                x,
                y,
                customdata,
                type: 'scatter',
                mode: 'lines+markers',
                name: this._escapeHtml(name),
                line: { color, width: 1.5 },
                marker: { color, size: markerSizes },
                hovertemplate: `${this._escapeHtml(spectrumLabel)} m/z: %{customdata[0]}`
                    + `<br>${this._escapeHtml(labels.relativeIntensity)}: %{customdata[1]:.1f}%<extra></extra>`,
                connectgaps: false,
            });
        };

        const qMatchedIndices = qNorm.map((_, index) => index).filter(index => matchedQ.has(index));
        const qOtherIndices = qNorm.map((_, index) => index).filter(index => !matchedQ.has(index));
        const lMatchedIndices = lNorm.map((_, index) => index).filter(index => matchedL.has(index));
        const lOtherIndices = lNorm.map((_, index) => index).filter(index => !matchedL.has(index));
        addTrace(qNorm, qOtherIndices, 1, labels.otherQuery, '#3b82f6', labels.query);
        addTrace(qNorm, qMatchedIndices, 1, labels.matchedQuery, '#22c55e', labels.query);
        addTrace(lNorm, lOtherIndices, -1, labels.otherLibrary, '#ef4444', labels.library);
        addTrace(lNorm, lMatchedIndices, -1, labels.matchedLibrary, '#16a34a', labels.library);

        const allMz = qNorm.concat(lNorm).map(peak => peak.mz);
        const maxMz = allMz.length > 0 ? Math.max(...allMz) * 1.05 : 1;
        const hasReference = lNorm.length > 0;
        const hasAnyPeaks = qNorm.length > 0 || hasReference;
        const annotations = [];
        if (!hasAnyPeaks) {
            annotations.push(this._statusAnnotation(labels.empty, 0.5));
        } else {
            if (qNorm.length > 0) {
                annotations.push({
                    x: 0.02, y: 0.95, xref: 'paper', yref: 'paper',
                    text: `<b>${queryName}</b>`, showarrow: false,
                    font: { size: 11, color: '#3b82f6' },
                });
            }
            if (hasReference) {
                annotations.push({
                    x: 0.02, y: 0.05, xref: 'paper', yref: 'paper',
                    text: `<b>${libraryName}</b>`, showarrow: false,
                    font: { size: 11, color: '#ef4444' },
                });
            } else {
                annotations.push(this._statusAnnotation(labels.referenceMissing, 0.12));
            }
        }
        annotations.push({
            x: 0.98, y: 0.5, xref: 'paper', yref: 'paper', xanchor: 'right',
            text: `${this._escapeHtml(labels.locatedPairs)}: ${locatedPairs.length}/${recordedPairs.length}`,
            showarrow: false,
            font: { size: 10, color: locatedPairs.length > 0 ? '#15803d' : '#6b7280' },
        });

        const baseTitle = this._displayName(settings.title, 'MS/MS spectrum comparison');
        const hasSavedScore = typeof settings.score === 'number' && Number.isFinite(settings.score);
        const scoreLabel = this._displayName(settings.scoreLabel, 'Saved cosine');
        const precision = Number.isInteger(settings.scorePrecision)
            && settings.scorePrecision >= 0 && settings.scorePrecision <= 10 ? settings.scorePrecision : null;
        const scoreText = hasSavedScore && precision !== null
            ? settings.score.toFixed(precision) : String(settings.score);
        const title = hasSavedScore
            ? `${this._escapeHtml(baseTitle)} — ${this._escapeHtml(scoreLabel)}: ${scoreText}`
            : this._escapeHtml(baseTitle);
        const height = typeof settings.height === 'number' && Number.isFinite(settings.height)
            && settings.height > 0 ? settings.height : 450;

        Plotly.newPlot(container, traces, {
            title: { text: title, font: { size: 14 } },
            xaxis: {
                title: this._escapeHtml(labels.mz),
                range: [0, maxMz],
                showgrid: true,
                gridcolor: '#f0f0f0',
            },
            yaxis: {
                title: this._escapeHtml(labels.relativeIntensity),
                range: hasReference ? [-115, 115] : [0, 115],
                showgrid: true,
                gridcolor: '#f0f0f0',
            },
            shapes: [{
                type: 'line', x0: 0, x1: maxMz, y0: 0, y1: 0,
                line: { color: '#d1d5db', width: 1 },
            }],
            annotations,
            margin: { t: 50, b: 50, l: 60, r: 20 },
            height,
            legend: { x: 0.5, y: 1.12, orientation: 'h' },
            hovermode: 'closest',
        }, {
            responsive: true,
            modeBarButtonsToAdd: ['toImage'],
        });
    }

    _finiteNumber(value) {
        if (typeof value === 'number') return Number.isFinite(value) ? value : null;
        if (typeof value !== 'string' || value.trim() === '') return null;
        const numeric = Number(value);
        return Number.isFinite(numeric) ? numeric : null;
    }

    _validPeaks(peaks) {
        if (!Array.isArray(peaks)) return [];
        const valid = [];
        peaks.forEach((peak, sourceIndex) => {
            if (!Array.isArray(peak) || peak.length < 2) return;
            const mz = this._finiteNumber(peak[0]);
            const intensity = this._finiteNumber(peak[1]);
            if (mz === null || intensity === null || mz <= 0 || intensity < 0) return;
            valid.push({ mz, intensity, originalMz: String(peak[0]), sourceIndex });
        });
        return valid;
    }

    _normalize(peaks) {
        const maxIntensity = peaks.length > 0
            ? Math.max(...peaks.map(peak => peak.intensity)) : 0;
        return peaks.map(peak => ({
            ...peak,
            relativeIntensity: maxIntensity > 0 ? (peak.intensity / maxIntensity) * 100 : 0,
        }));
    }

    _pairEndpoints(pair) {
        if (Array.isArray(pair) && pair.length >= 2) return [pair[0], pair[1]];
        if (pair && typeof pair === 'object') return [pair.query_mz, pair.lib_mz];
        return [null, null];
    }

    _uniquePeakIndex(peaks, endpoint) {
        const mz = this._finiteNumber(endpoint);
        if (mz === null || mz <= 0) return null;
        const indices = [];
        peaks.forEach((peak, index) => {
            if (Math.abs(peak.mz - mz) <= 1e-9) indices.push(index);
        });
        return indices.length === 1 ? indices[0] : null;
    }

    _locatePairs(pairs, queryPeaks, libraryPeaks) {
        const located = [];
        pairs.forEach(pair => {
            const [queryMz, libraryMz] = this._pairEndpoints(pair);
            const hasIndices = pair && !Array.isArray(pair) && typeof pair === 'object'
                && ('query_peak_index' in pair || 'library_peak_index' in pair);
            const indexedPeak = (peaks, sourceIndex, endpoint) => {
                const mz = this._finiteNumber(endpoint);
                if (!Number.isInteger(sourceIndex) || sourceIndex < 0 || mz === null) return null;
                const index = peaks.findIndex(peak => peak.sourceIndex === sourceIndex);
                return index >= 0 && Math.abs(peaks[index].mz - mz) <= 1e-9 ? index : null;
            };
            const queryIndex = hasIndices
                ? indexedPeak(queryPeaks, pair.query_peak_index, queryMz)
                : this._uniquePeakIndex(queryPeaks, queryMz);
            const libraryIndex = hasIndices
                ? indexedPeak(libraryPeaks, pair.library_peak_index, libraryMz)
                : this._uniquePeakIndex(libraryPeaks, libraryMz);
            if (queryIndex !== null && libraryIndex !== null) {
                located.push({ queryIndex, libraryIndex });
            }
        });
        return located;
    }

    _displayName(value, fallback) {
        if (value === null || value === undefined) return fallback;
        const text = String(value).trim();
        return text === '' ? fallback : text;
    }

    _escapeHtml(value) {
        return String(value)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;');
    }

    _statusAnnotation(text, y) {
        return {
            x: 0.5, y, xref: 'paper', yref: 'paper',
            text: this._escapeHtml(text), showarrow: false,
            font: { size: 11, color: '#6b7280' },
        };
    }

    clear() {
        if (this.containerId) {
            const el = document.getElementById(this.containerId);
            if (el && typeof Plotly !== 'undefined') Plotly.purge(el);
        }
    }
}

window.SpectrumComparator = SpectrumComparator;
