/**
 * MSsystem — Unified API Layer
 * Handles all HTTP requests + SSE real-time updates
 */
const API = (() => {
  const BASE = '/api';
  let _sseSource = null;
  let _sseListeners = [];

  // ── HTTP Helpers ──────────────────────────────────────────────────
  async function get(path) {
    const res = await fetch(BASE + path);
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(err.detail || '请求失败');
    }
    return res.json();
  }

  async function post(path, body) {
    const isFormData = body instanceof FormData;
    const res = await fetch(BASE + path, {
      method: 'POST',
      headers: isFormData ? {} : { 'Content-Type': 'application/json' },
      body: isFormData ? body : JSON.stringify(body),
    });
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(err.detail || '请求失败');
    }
    return res.json();
  }

  async function del(path, body = undefined) {
    const options = { method: 'DELETE' };
    if (body !== undefined) {
      options.headers = { 'Content-Type': 'application/json' };
      options.body = JSON.stringify(body);
    }
    const res = await fetch(BASE + path, options);
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(err.detail || '删除失败');
    }
    return res.json();
  }

  /**
   * Upload files with progress tracking
   * @param {File[]} files
   * @param {function} onProgress — callback(progressPercent)
   * @returns {Promise<object>} — { uploaded: [...], count: N }
   */
  function upload(files, onProgress) {
    return new Promise((resolve, reject) => {
      const formData = new FormData();
      for (const f of files) formData.append('files', f);

      const xhr = new XMLHttpRequest();
      xhr.open('POST', BASE + '/upload');

      xhr.upload.addEventListener('progress', e => {
        if (e.lengthComputable && onProgress) {
          onProgress(Math.round(e.loaded / e.total * 100));
        }
      });

      xhr.addEventListener('load', () => {
        if (xhr.status >= 200 && xhr.status < 300) {
          try { resolve(JSON.parse(xhr.responseText)); }
          catch (e) { reject(new Error('服务器返回格式错误')); }
        } else {
          try {
            const err = JSON.parse(xhr.responseText);
            reject(new Error(err.detail || `上传失败 (${xhr.status})`));
          } catch { reject(new Error(`上传失败 (${xhr.status})`)); }
        }
      });

      xhr.addEventListener('error', () => reject(new Error('无法连接到服务器，请检查服务是否启动')));
      xhr.addEventListener('abort', () => reject(new Error('上传已取消')));

      xhr.send(formData);
    });
  }

  /**
   * Upload one Agilent MassHunter .d directory and wait for validated mzML.
   * The browser exposes a directory as individual files; the relative-path
   * manifest lets the server reconstruct it inside a managed staging area.
   */
  function convertVendorFolder(files, options, onProgress) {
    return new Promise((resolve, reject) => {
      const selected = Array.from(files || []);
      const relativePaths = selected.map(f => f.webkitRelativePath || f.name);
      const formData = new FormData();
      formData.append('manifest', JSON.stringify(relativePaths));
      formData.append(
        'accept_vendor_licenses',
        options?.acceptVendorLicenses ? 'true' : 'false'
      );
      formData.append(
        'allow_container_pull',
        options?.allowContainerPull ? 'true' : 'false'
      );
      formData.append('backend', options?.backend || 'auto');
      formData.append(
        'peak_picking_strategy',
        options?.peakPickingStrategy || 'auto'
      );
      for (const file of selected) formData.append('files', file, file.name);

      const xhr = new XMLHttpRequest();
      xhr.open('POST', BASE + '/vendor/convert');
      xhr.upload.addEventListener('progress', event => {
        if (event.lengthComputable && onProgress) {
          onProgress(Math.round(event.loaded / event.total * 100));
        }
      });
      xhr.addEventListener('load', () => {
        if (xhr.status >= 200 && xhr.status < 300) {
          try { resolve(JSON.parse(xhr.responseText)); }
          catch (error) { reject(new Error('转换服务返回格式错误')); }
        } else {
          try {
            const error = JSON.parse(xhr.responseText);
            reject(new Error(error.detail || `原始数据转换失败 (${xhr.status})`));
          } catch {
            reject(new Error(`原始数据转换失败 (${xhr.status})`));
          }
        }
      });
      xhr.addEventListener(
        'error',
        () => reject(new Error('无法连接到原始数据转换服务'))
      );
      xhr.addEventListener('abort', () => reject(new Error('原始数据转换已取消')));
      xhr.send(formData);
    });
  }

  /**
   * Download a file (blob) from the API
   * @param {string} path
   * @param {string} filename
   */
  async function download(path, filename) {
    try {
      const res = await fetch(BASE + path);
      if (!res.ok) {
        const err = await res.json().catch(() => ({ detail: res.statusText }));
        throw new Error(err.detail || '下载失败');
      }
      const blob = await res.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      Toast.show(`下载失败: ${e.message}`, 'error');
    }
  }

  // ── SSE Real-time Updates ─────────────────────────────────────────
  function connectSSE() {
    if (_sseSource) return; // already connected

    try {
      _sseSource = new EventSource(BASE + '/tasks/events');

      _sseSource.addEventListener('task_update', (event) => {
        try {
          const data = JSON.parse(event.data);
          _sseListeners.forEach(fn => fn(data));
        } catch (e) {
          console.warn('SSE parse error:', e);
        }
      });

      _sseSource.addEventListener('error', () => {
        console.warn('SSE connection lost, reconnecting in 3s...');
        if (_sseSource) { _sseSource.close(); _sseSource = null; }
        setTimeout(connectSSE, 3000);
      });

      console.log('SSE connected for real-time task updates');
    } catch (e) {
      console.warn('SSE not available, falling back to polling');
    }
  }

  function onTaskUpdate(callback) {
    _sseListeners.push(callback);
  }

  function removeTaskUpdateListener(callback) {
    _sseListeners = _sseListeners.filter(fn => fn !== callback);
  }

  return {
    get, post, del, upload, convertVendorFolder, download,
    connectSSE, onTaskUpdate, removeTaskUpdateListener,
    BASE,
  };
})();

// ═══════════════════════════════════════════════════════════
// Evidence-viewer API helpers — spectrum, classification and report.
// ═══════════════════════════════════════════════════════════

async function fetchSpectrum(taskId, spectrumIdx) {
    const resp = await fetch("/api/tasks/" + taskId + "/results");
    if (!resp.ok) return null;
    const results = await resp.json();
    if (!results || !results[spectrumIdx]) return null;
    const r = results[spectrumIdx];
    return {
        peaks: r.query_peaks || r.peaks || [],
        metadata: { precursor_mz: r.precursor_mz || r.precursormz || 0, name: r.name || r.compound_name || "" }
    };
}

async function fetchSpectrumComparison(taskId, spectrumIdx) {
    const resp = await fetch("/api/tasks/" + taskId + "/results");
    if (!resp.ok) return null;
    const results = await resp.json();
    if (!results || !results[spectrumIdx]) return null;
    const r = results[spectrumIdx];
    return {
        query: { peaks: r.query_peaks || r.peaks || [], name: "Query" },
        library: { peaks: r.lib_peaks || [], name: r.name || r.compound_name || "Library", smiles: r.smiles || "" },
        matchedPairs: r.matched_peaks || [],
        score: r.final_score || r.score || 0,
        cosine: r.cosine || 0
    };
}

async function fetchClassificationStats(taskId) {
    const resp = await fetch("/api/tasks/" + taskId + "/results");
    if (!resp.ok) return null;
    const results = await resp.json();
    if (!results) return null;
    const stats = { classCounts: {}, total: results.length };
    results.forEach(function(r) {
        const cls = r.classification || r.Classification || "Unknown";
        stats.classCounts[cls] = (stats.classCounts[cls] || 0) + 1;
    });
    return stats;
}

// ═══════════════════════════════════════════════════════════
// Spectrum Page Controller
// ═══════════════════════════════════════════════════════════

var SpectrumPage = {
    viewer: null,
    results: [],
    currentIdx: 0,

    init: function() {
        this.viewer = new SpectrumViewer("spectrumContainer");
    },

    loadTask: function(taskId) {
        if (!taskId) return;
        var self = this;
        fetch("/api/tasks/" + taskId + "/results").then(function(r) { return r.json(); }).then(function(data) {
            self.results = (data && data.results) || data || [];
            self.currentIdx = 0;
            document.getElementById("specControlsCard").classList.remove("hidden");
            self.showCurrent();
        });
    },

    showCurrent: function() {
        if (!this.results.length) return;
        var r = this.results[this.currentIdx];
        var m = r.match || {};
        var peaks = r.query_peaks || [];
        var name = m.compound_name || r.spectrum_id || "Spectrum " + this.currentIdx;
        if (peaks.length) this.viewer.render(peaks, { title: name });
        document.getElementById("specCounter").textContent = (this.currentIdx + 1) + "/" + this.results.length;
        document.getElementById("specMetadata").innerHTML =
            "Precursor m/z: " + (r.precursor_mz || "N/A").toFixed(4) +
            " | Score: " + ((m.final_score || 0).toFixed(3)) +
            " | Class: " + (m.classification || "Unknown");
    },

    prev: function() { if (this.currentIdx > 0) { this.currentIdx--; this.showCurrent(); } },
    next: function() { if (this.currentIdx < this.results.length - 1) { this.currentIdx++; this.showCurrent(); } }
};

// ═══════════════════════════════════════════════════════════
// Classification Page Controller
// ═══════════════════════════════════════════════════════════

var ClassPage = {
    dashboard: null,
    molViewer: null,
    comparator: null,
    results: [],

    init: function() {
        this.dashboard = new ClassificationDashboard();
        this.molViewer = new MoleculeViewer();
        this.comparator = new SpectrumComparator();
    },

    loadTask: function(taskId) {
        if (!taskId) {
            document.getElementById("classDashboard").classList.add("hidden");
            return;
        }
        var self = this;
        var dashboard = document.getElementById("classDashboard");
        dashboard.classList.remove("hidden");
        document.getElementById("molecule-grid").innerHTML = '<div class="class-empty">正在加载分类与结构数据…</div>';
        fetch("/api/tasks/" + taskId + "/results?matched_only=true&limit=5000").then(function(r) {
            if (!r.ok) throw new Error("分类数据加载失败");
            return r.json();
        }).then(function(data) {
            self.results = (data && data.results) || data || [];
            self.renderDashboard();
        }).catch(function(err) {
            self.results = [];
            document.getElementById("molecule-grid").innerHTML = '<div class="class-empty">' + err.message + '</div>';
        });
    },

    renderDashboard: function() {
        // Class accuracy — use match.classification (nested in match object)
        var classCounts = {};
        this.results.forEach(function(r) {
            var m = r.match || {};
            var cls = m.classification || "Unknown";
            classCounts[cls] = (classCounts[cls] || 0) + 1;
        });

        // Summary and a horizontal distribution chart remain readable with long class names.
        var sorted = Object.entries(classCounts).sort(function(a,b) { return b[1] - a[1]; });
        var structures = this.results.filter(function(r) {
            var m = r.match || {};
            return !!(m.has_structure !== false && m.smiles && !["N/A", "nan", "none", "null"].includes(String(m.smiles).trim()));
        });
        document.getElementById("classMatchedCount").textContent = this.results.length.toLocaleString();
        document.getElementById("classCategoryCount").textContent = sorted.length.toLocaleString();
        document.getElementById("classStructureCount").textContent = structures.length.toLocaleString();
        document.getElementById("classTopCategory").textContent = sorted.length ? sorted[0][0] : "—";

        var chartRows = sorted.slice(0, 24).reverse();
        var names = chartRows.map(function(x) { return x[0]; });
        var vals = chartRows.map(function(x) { return x[1]; });
        var total = Math.max(1, this.results.length);
        var chartHeight = Math.max(360, Math.min(680, 100 + chartRows.length * 25));
        document.getElementById("class-accuracy-chart").style.height = chartHeight + "px";

        Plotly.newPlot("class-accuracy-chart", [{
            x: vals, y: names, type: "bar", orientation: "h",
            marker: { color: "#3b82f6" },
            text: vals.map(function(v) { return v + " (" + (v/total*100).toFixed(1) + "%)"; }),
            textposition: "auto",
            hovertemplate: "%{y}<br>%{x} 条谱图<extra></extra>",
        }], {
            xaxis: { title: "匹配谱图数", rangemode: "tozero", gridcolor: "#edf2f7" },
            yaxis: { automargin: true, tickfont: { size: 11 } },
            margin: { t: 12, b: 52, l: 170, r: 28 }, height: chartHeight,
            paper_bgcolor: "white", plot_bgcolor: "white",
        }, { responsive: true, displaylogo: false });

        // Molecule gallery — deduplicate structures before rendering.
        var seenStructures = new Set();
        var withSmiles = structures.filter(function(r) {
            var m = r.match || {};
            var fullInChIKey = String(m.inchikey || "").trim().toUpperCase();
            if (!/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(fullInChIKey)) fullInChIKey = "";
            // Preserve stereochemical identities: the first 14 characters only
            // encode connectivity and are not a complete compound identifier.
            var key = fullInChIKey || m.smiles || [m.compound_name, m.formula].join("|");
            if (seenStructures.has(key)) return false;
            seenStructures.add(key);
            return true;
        }).slice(0, 12);
        var moleculeGrid = document.getElementById("molecule-grid");
        document.getElementById("classMoleculeCount").textContent = withSmiles.length + " 个";
        if (withSmiles.length > 0) {
            this.molViewer.renderGrid("molecule-grid", withSmiles.map(function(r) {
                var m = r.match || {};
                return { smiles: m.smiles, name: m.compound_name || "Unknown" };
            }), { minCardWidth: 190 });
        } else {
            moleculeGrid.innerHTML = '<div class="class-empty">本任务没有可展示的结构候选</div>';
        }

        // Populate compare selector — only matched results
        var matched = this.results.filter(function(r) { return r.match; }).slice(0, 30);
        var sel = document.getElementById("compareSelect");
        sel.innerHTML = '<option value="">-- 选择匹配项 --</option>';
        matched.forEach(function(r, i) {
            var m = r.match || {};
            var opt = document.createElement("option");
            opt.value = i;
            opt.textContent = (i+1) + ". " + (m.compound_name || "Unknown") + " (score:" + ((m.final_score||0).toFixed(3)) + ")";
            sel.appendChild(opt);
        });
        // Store matched for comparison lookup
        this._matchedForCompare = matched;
        if (window.Plotly && Plotly.Plots) {
            requestAnimationFrame(function() { Plotly.Plots.resize(document.getElementById("class-accuracy-chart")); });
        }
    },

    loadComparison: function(idx) {
        if (idx === "") return;
        var matched = this._matchedForCompare || this.results;
        var r = matched[parseInt(idx)];
        if (!r) return;
        var m = r.match || {};
        var queryPeaks = r.query_peaks || [];
        var libPeaks = r.lib_peaks || [];
        if (queryPeaks.length > 0) {
            this.comparator.render("compareContainer",
                { peaks: queryPeaks, name: "Query" },
                { peaks: libPeaks.length > 0 ? libPeaks : queryPeaks, name: m.compound_name || "Library", smiles: m.smiles },
                m.matched_peaks || []);
        }
    }
};
