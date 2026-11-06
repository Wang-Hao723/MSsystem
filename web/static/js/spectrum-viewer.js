/**
 * SpectrumViewer — Interactive Plotly.js mass spectrum viewer
 * GNPS-inspired: bar/stem chart, zoom/pan, m/z hover, top-N peak annotations
 */
class SpectrumViewer {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.currentPeaks = null;
    }

    /**
     * Render an interactive mass spectrum
     * @param {Array} peaks - [[mz, intensity], ...]
     * @param {Object} options - {color, mzRange, annotateCount, title, height}
     */
    render(peaks, options = {}) {
        if (!this.container || !peaks || peaks.length === 0) return;
        this.currentPeaks = peaks;

        const mz = peaks.map(p => p[0]);
        const intensities = peaks.map(p => p[1]);
        const maxInt = Math.max(...intensities);
        if (maxInt <= 0) return;

        const relInt = intensities.map(i => (i / maxInt) * 100);
        const color = options.color || '#3b82f6';
        const annotateCount = options.annotateCount || 5;

        // Top-N peaks for annotation
        const sorted = peaks.map((p, i) => ({mz: p[0], int: p[1], idx: i}))
                            .sort((a, b) => b.int - a.int)
                            .slice(0, annotateCount);
        const annotations = sorted.map(p => ({
            x: p.mz,
            y: (p.int / maxInt) * 100 + 3,
            text: p.mz.toFixed(2),
            showarrow: false,
            font: { size: 10, color: '#374151' },
        }));

        // 真正的质谱棒状图 (stem chart): 每个峰画 [底部0, 峰顶intensity, null] 用 null 分隔多段 line
        // 而非 Plotly bar (柱状图, 有宽度, m/z范围大时挤成一片, 这是用户报告"图形错误"的根因)
        const stemX = [];
        const stemY = [];
        peaks.forEach(p => {
            stemX.push(p[0], p[0], null);
            stemY.push(0, (p[1] / maxInt) * 100, null);
        });

        const stemTrace = {
            x: stemX,
            y: stemY,
            type: 'scatter',
            mode: 'lines',
            line: { color: color, width: 1.5 },
            hoverinfo: 'skip',
            showlegend: false,
            name: options.title || 'Spectrum',
        };

        // hover锚点: 在峰顶放透明 marker, 否则纯 line 段难悬浮取信息
        const hoverTrace = {
            x: mz,
            y: relInt,
            type: 'scatter',
            mode: 'markers',
            marker: { color: color, size: 6, opacity: 0.001 },
            hovertemplate: 'm/z: %{x:.4f}<br>Intensity: %{y:.1f}%<extra></extra>',
            showlegend: false,
        };

        const layout = {
            title: options.title ? { text: options.title, font: { size: 14 } } : undefined,
            xaxis: {
                title: { text: 'm/z', font: { size: 12 } },
                range: options.mzRange || [0, Math.max(...mz) * 1.05],
                showgrid: true,
                gridcolor: '#f0f0f0',
            },
            yaxis: {
                title: { text: 'Relative Intensity (%)', font: { size: 12 } },
                range: [0, 115],
                showgrid: true,
                gridcolor: '#f0f0f0',
            },
            annotations: annotations,
            hovermode: 'x unified',
            margin: { t: options.title ? 40 : 20, b: 50, l: 60, r: 20 },
            dragmode: 'zoom',
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            height: options.height || 380,
        };

        Plotly.newPlot(this.container, [stemTrace, hoverTrace], layout, {
            responsive: true,
            modeBarButtonsToAdd: ['toImage'],
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        });
    }

    /**
     * Render TIC chromatogram or all-scan overview.
     *
     * 两种调用形态:
     * 1) 旧 fallback: renderTIC([{rt, tic, matched, name}, ...]) — 来自 results 兜底, 仅 MS2 stem.
     * 2) 全扫描概览: renderTIC({ms1: [{rt,tic}], ms2: [{rt,tic,scan_id,...}], matched_scan_ids: Set, scan_name_map: Map})
     *    MS1 渲染为连续色谱图(实线), MS2 匹配点红色 marker, MS2 未匹配点灰色小点.
     *
     * rt 单位: 分钟
     */
    renderTIC(payload, options = {}) {
        if (!this.container || !payload) return;

        // 形态 1: 旧数组形式(无 tic_overview.json 的旧任务 fallback)
        if (Array.isArray(payload)) {
            return this._renderTICStems(payload, options);
        }

        // 形态 2: 标准 TIC 总离子流色谱图
        // 主曲线 = MS1 TIC 连续色谱曲线(质谱分析标准画法: RT vs 总离子流, 连续实线).
        // 文件无 MS1 时退化为 MS2 TIC 按 RT 连线. 匹配化合物在曲线上标红点, 不画淹没曲线的密集散点.
        const ms1 = (payload.ms1 || []).slice().sort((a, b) => a.rt - b.rt);
        const ms2 = (payload.ms2 || []).slice().sort((a, b) => a.rt - b.rt);
        const matchedIds = payload.matched_scan_ids instanceof Set
            ? payload.matched_scan_ids
            : new Set(payload.matched_scan_ids || []);
        const nameMap = payload.scan_name_map || new Map();

        const cleanPoint = p => ({
            ...p,
            rt: Number(p.rt),
            tic: Number(p.tic || 0),
        });
        const validMs1 = ms1.map(cleanPoint).filter(p => Number.isFinite(p.rt) && Number.isFinite(p.tic));
        const validMs2 = ms2.map(cleanPoint).filter(p => Number.isFinite(p.rt) && Number.isFinite(p.tic));
        const usefulMs1 = validMs1.filter(p => p.tic > 0);
        const usefulMs2 = validMs2.filter(p => p.tic > 0);

        // 色谱曲线数据源: 有有效 MS1 survey scan 用 MS1(标准 TIC); 纯 MS2 文件或 MS1 TIC 全 0 时用 MS2 分箱聚合
        // 近似, 避免数千稀疏 DDA 点连成锯齿线(用户报告"看不到 TIC 曲线"的主因)。
        let main, srcLabel, isMs2Approx = false;
        if (usefulMs1.length >= 2) {
            main = validMs1; srcLabel = 'MS1';
        } else if (usefulMs2.length) {
            main = this._binTIC(usefulMs2); srcLabel = 'MS2'; isMs2Approx = true;
        } else {
            main = [];
        }

        const traces = [];
        const yAt = rt => {
            if (!main.length) return 0;
            let best = main[0], bd = Infinity;
            for (const m of main) {
                const d = Math.abs(m.rt - rt);
                if (d < bd) { bd = d; best = m; }
            }
            return best ? (best.tic || 0) : 0;
        };

        // ① TIC 连续色谱曲线(主体); MS2 近似时用 spline 平滑插值
        if (main.length) {
            traces.push({
                x: main.map(p => p.rt), y: main.map(p => p.tic),
                type: 'scatter', mode: 'lines',
                line: { color: '#1f77b4', width: 1.4, shape: isMs2Approx ? 'spline' : 'linear' },
                name: `${srcLabel} TIC${isMs2Approx ? ' (分箱近似)' : ''} (${main.length} ${isMs2Approx ? 'bins' : 'scans'})`,
                hovertemplate: 'RT: %{x:.2f} min<br>TIC: %{y:,.0f}<extra>' + srcLabel + ' TIC</extra>',
            });
        }

        // ② 全部 MS2 scan 位置: 灰色小点贴在 TIC 曲线上，避免只显示一条线看不到 scan 分布。
        if (usefulMs2.length && main.length) {
            const scanPoints = this._downsamplePoints(usefulMs2, 5000);
            traces.push({
                x: scanPoints.map(p => p.rt),
                y: scanPoints.map(p => yAt(p.rt)),
                type: 'scatter', mode: 'markers',
                marker: { color: '#94a3b8', size: 4, opacity: 0.38 },
                text: scanPoints.map(p => {
                    const scan = p.scan_id != null ? `Scan ${p.scan_id}` : 'MS2 scan';
                    const mz = p.precursor_mz ? `<br>m/z ${Number(p.precursor_mz).toFixed(4)}` : '';
                    return scan + mz;
                }),
                name: `MS2 scans (${usefulMs2.length})`,
                hovertemplate: 'RT: %{x:.2f} min<br>%{text}<extra>MS2 scan</extra>',
            });
        }

        // ③ 匹配化合物: 在曲线对应 RT 处标红点(贴着 TIC 曲线)
        const idHits = id => matchedIds.has(id) || matchedIds.has(String(id)) || matchedIds.has(Number(id));
        const matched = usefulMs2.filter(p => idHits(p.scan_id));
        if (matched.length && main.length) {
            traces.push({
                x: matched.map(p => p.rt), y: matched.map(p => yAt(p.rt)),
                type: 'scatter', mode: 'markers',
                marker: { color: '#dc2626', size: 8, line: { color: '#7f1d1d', width: 0.5 } },
                text: matched.map(p => {
                    const name = nameMap.get(p.scan_id) || nameMap.get(String(p.scan_id)) || nameMap.get(Number(p.scan_id)) || '';
                    return `Scan ${p.scan_id}${name ? `<br>${name}` : ''}`;
                }),
                name: `匹配化合物 (${matched.length})`,
                hovertemplate: 'RT: %{x:.2f} min<br>%{text}<extra>匹配</extra>',
            });
        }

        if (!traces.length) {
            this.container.innerHTML = '<div style="padding:2rem;text-align:center;color:#999;font-size:0.85rem;">无可绘制的 RT/TIC 数据</div>';
            return;
        }

        const titleText = usefulMs1.length >= 2
            ? `总离子流色谱图 (TIC) — ${validMs1.length} MS1 scans / ${usefulMs2.length} MS2 scans`
            : `总离子流色谱图 (TIC · MS2 近似) — 原始 ${usefulMs2.length} 个 MS2 scans 已按 RT 分箱聚合`;
        const layout = {
            title: { text: options.title || titleText, font: { size: 13 } },
            xaxis: { title: { text: 'Retention Time (min)', font: { size: 12 } }, showgrid: true, gridcolor: '#eee' },
            yaxis: { title: { text: 'Intensity', font: { size: 12 } }, showgrid: true, gridcolor: '#eee',
                     type: options.logY ? 'log' : 'linear' },
            margin: { t: 40, b: 50, l: 70, r: 20 },
            height: options.height || 300,
            paper_bgcolor: 'white', plot_bgcolor: 'white',
            legend: { orientation: 'h', y: 1.12 },
            hovermode: 'closest',
        };
        Plotly.newPlot(this.container, traces, layout, {
            responsive: true, modeBarButtonsToAdd: ['toImage'], displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        });
    }

    _downsamplePoints(points, maxN) {
        if (!points || points.length <= maxN) return points || [];
        const step = Math.ceil(points.length / maxN);
        const out = [];
        for (let i = 0; i < points.length; i += step) out.push(points[i]);
        return out;
    }

    /**
     * 纯 MS2 文件(无 MS1 survey scan)的 TIC 近似: 把稀疏的 DDA MS2 扫描按 RT 分箱聚合,
     * 每 bin 取该窗口内所有 MS2 扫描 TIC 之和 ≈ 该 RT 窗口的总离子流, 画出平滑色谱曲线。
     * bin 宽自适应(RT 跨度/400, 钳到 0.02~0.1 min), 避免几千点锯齿。
     */
    _binTIC(points) {
        if (!points || !points.length) return [];
        let rMin = Infinity, rMax = -Infinity;
        for (const p of points) { if (p.rt < rMin) rMin = p.rt; if (p.rt > rMax) rMax = p.rt; }
        const span = Math.max(rMax - rMin, 1e-6);
        let binW = span / 400;
        binW = Math.max(0.02, Math.min(0.1, binW));
        const bins = new Map();
        for (const p of points) {
            const b = Math.round((p.rt - rMin) / binW);
            const cur = bins.get(b);
            const t = p.tic || 0;
            if (cur) { cur.tic += t; cur.n += 1; }
            else bins.set(b, { tic: t, n: 1 });
        }
        const out = [];
        for (const [b, v] of bins) out.push({ rt: rMin + b * binW, tic: v.tic });
        return out.sort((a, b) => a.rt - b.rt);
    }

    /** 旧 fallback: 仅 MS2 stem chart, 用于无 tic_overview.json 的旧任务 */
    _renderTICStems(points, options = {}) {
        if (!points.length) return;
        const matched = points.filter(p => p.matched);

        const allStemX = []; const allStemY = [];
        points.forEach(p => { allStemX.push(p.rt, p.rt, null); allStemY.push(0, p.tic, null); });
        const mStemX = []; const mStemY = [];
        matched.forEach(p => { mStemX.push(p.rt, p.rt, null); mStemY.push(0, p.tic, null); });

        const traces = [
            { x: allStemX, y: allStemY, type: 'scatter', mode: 'lines',
              line: { color: '#cbd5e0', width: 1 }, hoverinfo: 'skip',
              name: `全部 scan (${points.length})`, showlegend: true },
            { x: mStemX, y: mStemY, type: 'scatter', mode: 'lines',
              line: { color: '#3b82f6', width: 1.5 }, hoverinfo: 'skip',
              name: `匹配 (${matched.length})`, showlegend: true },
            { x: matched.map(p => p.rt), y: matched.map(p => p.tic),
              type: 'scatter', mode: 'markers',
              marker: { color: '#3b82f6', size: 6, opacity: 0.6 },
              text: matched.map(p => p.name || ''),
              hovertemplate: 'RT: %{x:.2f} min<br>TIC: %{y:,.0f}<br>%{text}<extra></extra>',
              showlegend: false },
        ];
        Plotly.newPlot(this.container, traces, {
            title: { text: options.title || '总离子流图 (TIC) — MS2 扫描概览 (旧任务 fallback)', font: { size: 13 } },
            xaxis: { title: { text: 'Retention Time (min)', font: { size: 12 } }, showgrid: true, gridcolor: '#f0f0f0' },
            yaxis: { title: { text: 'TIC (sum intensity)', font: { size: 12 } }, showgrid: true, gridcolor: '#f0f0f0' },
            margin: { t: 40, b: 50, l: 70, r: 20 },
            height: options.height || 280,
            paper_bgcolor: 'white', plot_bgcolor: 'white',
            legend: { orientation: 'h', y: 1.12 },
        }, { responsive: true, modeBarButtonsToAdd: ['toImage'], displayModeBar: true,
             modeBarButtonsToRemove: ['lasso2d', 'select2d'] });
    }

    /**
     * Render multiple spectra overlaid
     * @param {Array} spectraList - [{peaks, name, color}, ...]
     */
    renderOverlay(spectraList, options = {}) {
        if (!this.container || !spectraList.length) return;

        const traces = spectraList.map(spec => {
            const mz = spec.peaks.map(p => p[0]);
            const maxInt = Math.max(...spec.peaks.map(p => p[1]));
            const relInt = spec.peaks.map(p => (p[1] / maxInt) * 100);
            return {
                x: mz, y: relInt,
                type: 'bar',
                name: spec.name || 'Spectrum',
                marker: { color: spec.color || '#3b82f6', opacity: 0.6 },
                hovertemplate: `${spec.name || ''} m/z: %{x:.4f}<br>Intensity: %{y:.1f}%<extra></extra>`,
            };
        });

        Plotly.newPlot(this.container, traces, {
            title: options.title ? { text: options.title, font: { size: 14 } } : undefined,
            barmode: 'overlay',
            xaxis: { title: 'm/z' },
            yaxis: { title: 'Relative Intensity (%)', range: [0, 115] },
            margin: { t: 40, b: 50, l: 60, r: 20 },
            height: options.height || 380,
        }, { responsive: true, modeBarButtonsToAdd: ['toImage'] });
    }

    /** Clear the viewer */
    clear() {
        if (this.container) Plotly.purge(this.container);
        this.currentPeaks = null;
    }
}

// Global instance
window.SpectrumViewer = SpectrumViewer;
