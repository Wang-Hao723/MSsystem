/**
 * MSsystem — Data Upload Module
 * Handles file selection, drag & drop, validation, upload progress, analysis start
 */
const Upload = (() => {
  let _stagedFiles = [];
  let _vendorFolderFiles = [];
  let _vendorRootName = '';
  let _vendorLocalSelections = [];
  let _vendorConverterReady = false;
  let _isUploading = false;
  let _retryUploadedTasks = []; // 上传已成功但排队失败时复用，避免再次上传/刷新
  let _blankPaths = [];   // 已上传空白的受控服务器 token，可包含同批重复空白
  let _blankUploads = []; // 当前选中的空白元数据；跨样本批次保留，直到用户主动解除
  let _availableBlankUploads = []; // 服务器已保存的空白；刷新后仍可选择或删除
  let _blankManagerDeleteSelection = new Set(); // 空白管理中的批量删除勾选
  let _loadedDbMode = null;
  let _blankEvidenceV2Available = false;
  let _expectedDatasets = [];
  let _expectedDatasetsLoaded = false;

  const ALLOWED_EXTS = ['.mgf', '.mzml', '.mzxml', '.xml', '.html', '.htm'];
  const MAX_FILE_SIZE_MB = 500;
  const MAX_BATCH_FILES = 100;
  const BLANK_SELECTION_STORAGE_KEY = 'mssystem.v7.blank-selection';
  const VENDOR_LICENSE_STORAGE_KEY = 'mssystem.v7.vendor-license-accepted.v1';
  const BLANK_ROLE_LABELS = {
    unspecified: '类型未注明（仅作证据）',
    background_blank: '背景空白（通用）',
    solvent_blank: '溶剂空白',
    method_blank: '方法空白',
    extraction_blank: '提取空白',
    carryover_blank: '进样后残留 / Carryover',
  };
  const BACKGROUND_EXCLUSION_ROLES = new Set([
    'background_blank', 'solvent_blank', 'method_blank', 'extraction_blank',
  ]);

  // ── Initialization ────────────────────────────────────────────────
  function init() {
    const zone = document.getElementById('uploadZone');
    if (!zone) return;

    // 动态拉取真实库规模填到下拉选项 (修复旧版 355K/894K 写死不更新)
    refreshDbCounts();
    refreshRuntimeIdentity();
    void refreshBlankUploads();
    void refreshVendorConverterStatus();

    // Drag & drop
    ['dragenter', 'dragover'].forEach(evt => {
      zone.addEventListener(evt, e => { e.preventDefault(); zone.classList.add('dragover'); });
    });
    ['dragleave', 'drop'].forEach(evt => {
      zone.addEventListener(evt, e => { e.preventDefault(); zone.classList.remove('dragover'); });
    });
    zone.addEventListener('drop', e => handleFiles(e.dataTransfer.files));

    // 普通多选 + 文件夹多选。每次选完即清空 input，移除后可再选同一文件。
    ['fileInput', 'folderInput'].forEach(id => {
      const input = document.getElementById(id);
      if (!input) return;
      input.addEventListener('change', e => {
        handleFiles(e.target.files);
        e.target.value = '';
      });
    });
    const vendorInput = document.getElementById('vendorRawInput');
    if (vendorInput) {
      vendorInput.addEventListener('change', event => {
        handleVendorFolder(event.target.files);
        event.target.value = '';
      });
    }
    restoreVendorLicenseAcceptance();
    const ionMode = document.getElementById('ionMode');
    if (ionMode) {
      ionMode.addEventListener('change', () => {
        void handleIonModeChange(ionMode.value);
      });
    }
  }

  function restoreVendorLicenseAcceptance() {
    const checkbox = document.getElementById('vendorLicenseAccepted');
    if (!checkbox) return;
    try {
      checkbox.checked = window.localStorage.getItem(
        VENDOR_LICENSE_STORAGE_KEY
      ) === 'accepted';
    } catch {
      checkbox.checked = false;
    }
    checkbox.addEventListener('change', () => {
      try {
        if (checkbox.checked) {
          window.localStorage.setItem(
            VENDOR_LICENSE_STORAGE_KEY,
            'accepted'
          );
        } else {
          window.localStorage.removeItem(VENDOR_LICENSE_STORAGE_KEY);
        }
      } catch {
        // Private browsing/storage denial: keep the current in-page choice.
      }
    });
  }

  // ── File Handling ─────────────────────────────────────────────────
  function fileIdentity(file) {
    const relativeName = file.webkitRelativePath || file.name;
    return `${relativeName}\u0000${file.size}\u0000${file.lastModified || 0}`;
  }

  function handleFiles(fileList) {
    _retryUploadedTasks = [];
    let unsupported = 0;
    let oversized = 0;
    let duplicates = 0;
    let limitReached = false;
    for (const f of Array.from(fileList || [])) {
      if (_stagedFiles.length >= MAX_BATCH_FILES) {
        limitReached = true;
        break;
      }
      const ext = '.' + f.name.split('.').pop().toLowerCase();
      if (!ALLOWED_EXTS.includes(ext)) {
        unsupported += 1;
        continue;
      }
      const sizeMB = f.size / (1024 * 1024);
      if (sizeMB > MAX_FILE_SIZE_MB) {
        oversized += 1;
        continue;
      }
      // 同一文件不重复加入；不同子目录的同名文件仍可分别分析。
      if (_stagedFiles.some(u => fileIdentity(u) === fileIdentity(f))) {
        duplicates += 1;
        continue;
      }
      _stagedFiles.push(f);
    }
    if (limitReached) {
      Toast.show(`单批最多选择 ${MAX_BATCH_FILES} 个样本，其余文件未加入`, 'warning');
    }
    if (unsupported) {
      Toast.show(`已忽略 ${unsupported} 个不支持的文件（仅接受 MGF / mzML / mzXML / 可解析 HTML）`, 'warning');
    }
    if (oversized) {
      Toast.show(`已忽略 ${oversized} 个超过 ${MAX_FILE_SIZE_MB} MB 上传上限的文件`, 'warning');
    }
    if (duplicates) {
      Toast.show(`已忽略 ${duplicates} 个重复选择的文件`, 'info');
    }
    renderFileList();
    updateUI();
  }

  function handleVendorFolder(fileList) {
    _retryUploadedTasks = [];
    const files = Array.from(fileList || []);
    const status = document.getElementById('vendorRawStatus');
    _vendorFolderFiles = [];
    _vendorRootName = '';
    _vendorLocalSelections = [];
    if (!files.length) {
      renderVendorFolder();
      updateUI();
      return;
    }
    const relativePaths = files.map(file => file.webkitRelativePath || file.name);
    const roots = [...new Set(relativePaths.map(path => String(path).split('/')[0]))];
    const lowerPaths = relativePaths.map(path => String(path).toLowerCase());
    const required = [
      '/acqdata/msscan.bin',
      '/acqdata/mspeak.bin',
      '/acqdata/acqmethod.xml',
    ];
    if (
      roots.length !== 1
      || !roots[0].toLowerCase().endsWith('.d')
      || !required.every(suffix => lowerPaths.some(path => path.endsWith(suffix)))
    ) {
      if (status) {
        status.innerHTML = '<span style="color:var(--danger);">请选择一个完整的 Agilent MassHunter .d 文件夹（必须含 AcqData/MSScan.bin、MSPeak.bin、AcqMethod.xml）。</span>';
      }
      Toast.show('所选文件夹不是完整的 Agilent MassHunter .d 原始数据', 'error');
      updateUI();
      return;
    }
    _vendorFolderFiles = files;
    _vendorRootName = roots[0];
    _setVendorPeakPicking(
      'auto',
      '等待方法检测 · centroid 将保留，profile 才质心化'
    );
    _setVendorAcquisitionNote('unknown');
    renderVendorFolder();
    void refreshVendorConverterStatus();
    updateUI();
  }

  async function _selectLocalVendorSources({
    endpoint,
    buttonId,
    openingText,
    sourceLabel,
    defaultButtonLabel,
  }) {
    const status = document.getElementById('vendorRawStatus');
    const button = document.getElementById(buttonId);
    const previousLabel = button?.textContent || '';
    if (button) {
      button.disabled = true;
      button.textContent = '正在打开本机目录选择器…';
    }
    if (status) {
      status.innerHTML = `<span style="color:var(--text-secondary);">${escapeHtml(openingText)}</span>`;
    }
    try {
      const selection = await API.post(endpoint, {});
      if (selection.cancelled) {
        renderVendorFolder();
        Toast.show(`已取消选择 ${sourceLabel}`, 'info');
        return;
      }
      _retryUploadedTasks = [];
      _vendorFolderFiles = [];
      const rawSelections = (
        Array.isArray(selection.selections) && selection.selections.length
          ? selection.selections
          : [selection]
      );
      _vendorLocalSelections = rawSelections.map(item => ({
        selection_token: String(item.selection_token || ''),
        raw_name: String(item.raw_name || ''),
        source_kind: String(item.source_kind || 'agilent_d'),
        source_format: String(item.source_format || ''),
        relative_source: String(item.relative_source || item.raw_name || ''),
        intake_role: String(item.intake_role || 'sample'),
        expected_polarities: Array.isArray(item.expected_polarities)
          ? item.expected_polarities.map(mode => String(mode).toUpperCase())
          : [],
        expires_in_seconds: Number(
          item.expires_in_seconds
          || selection.expires_in_seconds
          || 0
        ),
        acquisition: item.acquisition || {},
        spectrum_representation: String(
          item.spectrum_representation || 'unknown'
        ),
        recommended_peak_picking_strategy: String(
          item.recommended_peak_picking_strategy || 'auto'
        ),
      })).filter(item => item.selection_token && item.raw_name);
      if (!_vendorLocalSelections.length) {
        throw new Error(`本机选择器没有返回可读取的 ${sourceLabel}`);
      }
      _vendorRootName = String(
        selection.selected_root_name
        || _vendorLocalSelections[0].raw_name
      );
      const representation = _vendorLocalSelections.length === 1
        ? _vendorLocalSelections[0].spectrum_representation
        : 'mixed';
      _setVendorPeakPicking(
        _vendorLocalSelections.length === 1
          ? _vendorLocalSelections[0].recommended_peak_picking_strategy
          : 'auto',
        _vendorLocalSelections.length === 1
          ? representation === 'centroid'
            ? 'AcqMethod.xml 已确认 centroid · 原样保留'
            : representation === 'profile'
            ? 'AcqMethod.xml 检出 profile · 自动质心化'
            : representation === 'mixed'
            ? 'AcqMethod.xml 检出混合表征 · 逐谱自动处理'
            : '方法表征未知 · 转换后逐谱检查'
          : `${_vendorLocalSelections.length} 个厂商原始样本 · 逐个转换并核验`
      );
      _setVendorAcquisitionNote(
        _vendorLocalSelections.some(item => item.source_kind === 'thermo_raw')
          ? 'thermo'
          : _vendorLocalSelections.length === 1 ? representation : 'mixed'
      );
      renderVendorFolder();
      await refreshVendorConverterStatus();
      updateUI();
      const invalidCount = (
        Array.isArray(selection.invalid_directories)
          ? selection.invalid_directories.length
          : Array.isArray(selection.invalid_files)
          ? selection.invalid_files.length
          : 0
      );
      Toast.show(
        _vendorLocalSelections.length === 1
          ? representation === 'centroid'
            ? `已选择 ${_vendorRootName}：centroid 峰表将原样保留`
            : `已选择：${_vendorRootName}`
          : `已批量选择 ${_vendorLocalSelections.length} 个厂商原始样本`
            + `${invalidCount ? `；另有 ${invalidCount} 个无效输入已跳过` : ''}`,
        invalidCount ? 'warning' : 'success'
      );
    } catch (error) {
      renderVendorFolder();
      Toast.show(`本机 ${sourceLabel} 选择失败：${error.message}`, 'error');
    } finally {
      if (button) {
        button.disabled = false;
        button.textContent = previousLabel || defaultButtonLabel;
      }
    }
  }

  async function selectLocalVendorFolder() {
    return _selectLocalVendorSources({
      endpoint: '/vendor/select-local-folder',
      buttonId: 'btnSelectLocalVendorFolder',
      openingText: '本机目录选择器已打开，请选中整个 xxx.d 或其批次父目录后点“选择”。',
      sourceLabel: 'Agilent .d',
      defaultButtonLabel: '选择 Agilent .d / 父目录',
    });
  }

  async function selectLocalThermoFolder() {
    return _selectLocalVendorSources({
      endpoint: '/vendor/select-local-thermo-folder',
      buttonId: 'btnSelectLocalThermoFolder',
      openingText: '本机目录选择器已打开，请选含 Thermo .raw 的目录；支持 POS/NEG 一层子目录。',
      sourceLabel: 'Thermo RAW',
      defaultButtonLabel: '选择 Thermo QE .raw / 批次目录',
    });
  }

  async function refreshVendorConverterStatus() {
    const area = document.getElementById('vendorConverterStatus');
    const badge = document.getElementById('vendorConverterBadge');
    const setup = document.getElementById('vendorConverterSetup');
    if (!area) return;
    area.textContent = '正在检查 ProteoWizard 转换后端…';
    _setConverterBadge(badge, '正在检测…', 'is-checking');
    try {
      const status = await API.get('/vendor/converter/status');
      if (status.local_msconvert?.available) {
        _vendorConverterReady = true;
        area.innerHTML = '✅ 已找到本机 msconvert；转换完全在本机完成。';
        _setConverterBadge(badge, '本机就绪 · 断网可用', 'is-ready');
      } else if (status.docker?.daemon_available && status.docker?.image_available) {
        _vendorConverterReady = true;
        area.innerHTML = '✅ ProteoWizard Docker 镜像已缓存在本机；后续转换不需要网络。';
        _setConverterBadge(badge, '已缓存 · 断网可用', 'is-ready');
      } else if (status.docker?.daemon_available) {
        _vendorConverterReady = false;
        area.innerHTML = '⚠️ Docker 已启动，但本机没有转换镜像；请先单独点击“准备 / 下载”，分析过程不会联网下载。';
        _setConverterBadge(badge, '首次需联网准备', 'is-action');
      } else if (status.docker?.cli_available) {
        _vendorConverterReady = false;
        area.innerHTML = '⚠️ 已安装 Docker，但 Docker Desktop 尚未启动';
        _setConverterBadge(badge, '请启动 Docker', 'is-action');
      } else {
        _vendorConverterReady = false;
        area.innerHTML = '❌ 未找到 msconvert 或 Docker；请先安装本机可用的 ProteoWizard 或 Docker Desktop';
        _setConverterBadge(badge, '转换器不可用', 'is-error');
      }
      if (
        setup
        && (_vendorFolderFiles.length || _vendorLocalSelections.length)
        && !_vendorConverterReady
      ) {
        setup.open = true;
      }
    } catch (error) {
      _vendorConverterReady = false;
      area.textContent = `转换器状态检查失败：${error.message}`;
      _setConverterBadge(badge, '状态检查失败', 'is-error');
      if (setup && (_vendorFolderFiles.length || _vendorLocalSelections.length)) {
        setup.open = true;
      }
    }
  }

  function _setConverterBadge(element, text, stateClass) {
    if (!element) return;
    element.textContent = text;
    element.className = `converter-status-badge ${stateClass}`;
  }

  function _setVendorPeakPicking(strategy, summaryText) {
    const selector = document.getElementById('vendorPeakPicking');
    const summary = document.getElementById('vendorRepresentationSummary');
    if (
      selector
      && (
        !selector.options
        || Array.from(selector.options).some(option => option.value === strategy)
      )
    ) {
      selector.value = strategy;
    }
    if (summary) summary.textContent = summaryText;
  }

  function _setVendorAcquisitionNote(representation = 'unselected') {
    const title = document.getElementById('vendorAcquisitionPolicyTitle');
    const detail = document.getElementById('vendorAcquisitionPolicyDetail');
    if (!title || !detail) return;
    const copy = {
      centroid: [
        'AcqMethod.xml 已确认：centroid',
        '仅转换格式并原样保留峰表，不做二次质心化；MS1 与 MS2 均保留。',
      ],
      profile: [
        'AcqMethod.xml 检出：profile',
        '由 ProteoWizard 厂商算法完成质心化；转换后仍逐谱核验表征。',
      ],
      mixed: [
        'AcqMethod.xml 检出：centroid / profile 混合',
        'centroid 谱原样保留，仅对 profile 谱执行质心化。',
      ],
      unknown: [
        '转换时自动识别 centroid / profile',
        'centroid 原样保留；只有检测到 profile 时才执行质心化。',
      ],
      thermo: [
        'Thermo RAW：转换后逐谱核验',
        '读取精确仪器型号、极性、MS1/MS2、分辨率与隔离窗；profile 由厂商算法质心化，已有 centroid 不重复处理。',
      ],
      unselected: [
        'centroid 数据不二次质心化',
        '选择目录后读取 AcqMethod.xml；centroid 仅转换格式并保留峰表，MS1 与 MS2 均保留。',
      ],
    };
    [title.textContent, detail.textContent] = copy[representation] || copy.unknown;
  }

  async function prepareVendorConverter() {
    const accepted = document.getElementById('vendorLicenseAccepted')?.checked || false;
    const allowPull = document.getElementById('vendorAllowPull')?.checked || false;
    const area = document.getElementById('vendorConverterStatus');
    const button = document.getElementById('btnPrepareVendorConverter');
    if (!accepted) {
      Toast.show('请先由你本人阅读并勾选厂商读取器许可', 'warning');
      return;
    }
    const previousLabel = button?.textContent || '';
    if (button) {
      button.disabled = true;
      button.textContent = '正在准备（首次下载约 2.2 GB）…';
    }
    if (area) {
      area.innerHTML = allowPull
        ? '⏳ 正在检查；若镜像尚不存在，将开始首次下载。请保持 Docker Desktop 运行…'
        : '⏳ 正在检查本机已有转换器…';
    }
    try {
      const result = await API.post('/vendor/converter/prepare', {
        accept_vendor_licenses: true,
        allow_container_pull: allowPull,
      });
      if (result.downloaded) {
        Toast.show('ProteoWizard 镜像下载完成，转换器已就绪', 'success');
      } else {
        Toast.show('ProteoWizard 转换器已就绪，无需重复下载', 'success');
      }
    } catch (error) {
      Toast.show(`ProteoWizard 准备失败：${error.message}`, 'error');
    } finally {
      if (button) {
        button.disabled = false;
        button.textContent = previousLabel || '准备 / 下载 ProteoWizard 转换器';
      }
      await refreshVendorConverterStatus();
    }
  }

  function renderVendorFolder() {
    const status = document.getElementById('vendorRawStatus');
    const options = document.getElementById('vendorRawOptions');
    if (!_vendorFolderFiles.length && !_vendorLocalSelections.length) {
      if (status) status.innerHTML = '';
      if (options) {
        options.classList.add('hidden');
        options.open = false;
      }
      return;
    }
    if (status) {
      if (_vendorLocalSelections.length) {
        const selectionRows = _vendorLocalSelections.map((item, index) => {
          const acquisition = item.acquisition || {};
          const representation = item.spectrum_representation;
          const representationText = representation === 'centroid'
            ? 'centroid 保留'
            : representation === 'profile'
            ? 'profile 自动质心化'
            : representation === 'mixed'
            ? '混合表征逐谱处理'
            : '转换后逐谱核验';
          const polarities = Array.isArray(acquisition.polarities)
            && acquisition.polarities.length
            ? acquisition.polarities.join('/')
            : item.expected_polarities.length
            ? item.expected_polarities.join('/')
            : '极性待核验';
          const kind = item.source_kind === 'thermo_raw'
            ? 'Thermo RAW'
            : 'Agilent .d';
          const roleSelector = ` <select class="form-select form-select-inline"`
            + ` aria-label="${escapeHtml(item.raw_name)} 导入角色"`
            + ` onchange="Upload.updateVendorSourceRole(${index}, this.value)">`
            + `<option value="sample"${item.intake_role === 'sample' ? ' selected' : ''}>作为样本分析</option>`
            + `<option value="blank"${item.intake_role === 'blank' ? ' selected' : ''}>作为空白背景</option>`
            + '</select>';
          return `<li><strong>${escapeHtml(item.relative_source || item.raw_name)}</strong>`
            + ` · ${escapeHtml(kind)} · ${escapeHtml(polarities)} · ${escapeHtml(representationText)}`
            + `${roleSelector}</li>`;
        }).join('');
        status.innerHTML = `
          <div class="vendor-selection-row">
            <span>✅ <strong>已选择 ${_vendorLocalSelections.length} 个厂商原始样本</strong> · 本机只读、逐样本转换</span>
            <button type="button" class="btn btn-secondary btn-sm" onclick="Upload.clearVendorFolder()">移除全部</button>
          </div>
          <details${_vendorLocalSelections.length <= 6 ? ' open' : ''}>
            <summary>查看本批厂商原始数据清单</summary>
            <ul style="margin:.5rem 0 0 1.25rem;">${selectionRows}</ul>
          </details>`;
      } else {
        const totalMB = _vendorFolderFiles.reduce((sum, file) => sum + file.size, 0) / 1024 / 1024;
        status.innerHTML = `
          <div class="vendor-selection-row">
            <span>✅ <strong>${escapeHtml(_vendorRootName)}</strong> · ${_vendorFolderFiles.length} 个文件 · ${totalMB.toFixed(1)} MB · 转换时自动核验 centroid/profile</span>
            <button type="button" class="btn btn-secondary btn-sm" onclick="Upload.clearVendorFolder()">移除</button>
          </div>`;
      }
    }
    if (options) options.classList.remove('hidden');
  }

  function clearVendorFolder() {
    _retryUploadedTasks = [];
    _vendorFolderFiles = [];
    _vendorRootName = '';
    _vendorLocalSelections = [];
    _setVendorPeakPicking('auto', 'centroid 保留；仅在 profile 时质心化');
    _setVendorAcquisitionNote();
    const input = document.getElementById('vendorRawInput');
    if (input) input.value = '';
    renderVendorFolder();
    updateUI();
  }

  function updateVendorSourceRole(index, intakeRole) {
    const item = _vendorLocalSelections[index];
    if (!item) return;
    item.intake_role = intakeRole === 'blank' ? 'blank' : 'sample';
    renderVendorFolder();
    updateUI();
  }

  function removeFile(index) {
    _retryUploadedTasks = [];
    _stagedFiles.splice(index, 1);
    renderFileList();
    updateUI();
  }

  function clearFiles() {
    _retryUploadedTasks = [];
    _stagedFiles = [];
    _vendorFolderFiles = [];
    _vendorRootName = '';
    _vendorLocalSelections = [];
    _setVendorPeakPicking('auto', 'centroid 保留；仅在 profile 时质心化');
    _setVendorAcquisitionNote();
    ['fileInput', 'folderInput', 'vendorRawInput'].forEach(id => {
      const input = document.getElementById(id);
      if (input) input.value = '';
    });
    renderFileList();
    renderVendorFolder();
    updateUI();
  }

  function renderFileList() {
    const container = document.getElementById('fileList');
    if (!container) return;
    if (_stagedFiles.length === 0) {
      container.innerHTML = '';
      return;
    }
    const totalMB = _stagedFiles.reduce((sum, f) => sum + f.size, 0) / 1024 / 1024;
    container.innerHTML = `
      <div class="staged-files-toolbar">
        <div><strong>已选择 ${_stagedFiles.length} 个样本</strong><span> · 合计 ${totalMB.toFixed(1)} MB</span></div>
        <button type="button" class="btn btn-secondary btn-sm" onclick="Upload.clearFiles()">清空全部</button>
      </div>
      <div class="staged-files-list">${_stagedFiles.map((f, i) => `
      <div class="file-item">
        <span class="file-icon">📄</span>
        <span class="file-name">${escapeHtml(f.webkitRelativePath || f.name)}</span>
        <span class="file-size">${(f.size / 1024 / 1024).toFixed(2)} MB</span>
        <button class="file-remove" onclick="Upload.removeFile(${i})" title="移除">×</button>
      </div>`).join('')}</div>`;
  }

  // ── Library Management ────────────────────────────────────────────
  async function loadLibrary() {
    const btn = document.getElementById('btnLoadLib');
    const dbMode = document.getElementById('dbMode');
    if (!btn || !dbMode) return;

    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> 加载中...';
    try {
      const status = await API.post(`/library/load?db_mode=${dbMode.value}`);
      renderLibraryStatus(status);
      Toast.show(`数据库加载成功: ${status.spectra_count.toLocaleString()} 张谱图`, 'success');
    } catch (e) {
      Toast.show(`加载失败: ${e.message}`, 'error');
    } finally {
      btn.disabled = false;
      btn.textContent = '加载数据库';
    }
    updateUI();
  }

  async function checkLibraryStatus() {
    try {
      const status = await API.get('/library/status');
      renderLibraryStatus(status);
    } catch (e) {
      console.error('Library status check failed:', e);
    }
    updateUI();
  }

  function renderLibraryStatus(status) {
    const area = document.getElementById('libStatusArea');
    if (!area) return;
    if (status.loaded) {
      const dbMode = document.getElementById('dbMode');
      if (dbMode) {
        if (String(status.name || '').startsWith('Small')) _loadedDbMode = 'small';
        if (String(status.name || '').startsWith('Full')) _loadedDbMode = 'full';
        if (_loadedDbMode) dbMode.value = _loadedDbMode;
      }
      area.innerHTML = `
        <div class="lib-status loaded">
          <span>✅</span>
          <span><strong>${escapeHtml(status.name)}</strong> — ${status.spectra_count.toLocaleString()} 张谱图已加载 (${status.load_time_sec}s)</span>
        </div>`;
    } else {
      area.innerHTML = `
        <div class="lib-status not-loaded">
          <span>⚠️</span>
          <span>数据库未加载 — 请先选择数据库并点击"加载数据库"</span>
        </div>`;
    }
  }

  function numberValue(id, fallback) {
    const value = parseFloat(document.getElementById(id)?.value);
    return Number.isFinite(value) ? value : fallback;
  }

  function optionalNumberValue(id) {
    const raw = String(document.getElementById(id)?.value || '').trim();
    if (!raw) return null;
    const value = parseFloat(raw);
    return Number.isFinite(value) ? value : null;
  }

  function localizedCopy(zh, en) {
    return globalThis.I18n?.getLanguage?.() === 'en' ? en : zh;
  }

  function analysisParams(ionModeElementId = 'ionMode') {
    const params = {
      db_mode: document.getElementById('dbMode')?.value || 'full',
      cosine_threshold: numberValue('cosTh', 0.3),
      final_score_threshold: numberValue('finalTh', 0.25),
      mz_tolerance_da: numberValue('mzDa', 0.02),
      mz_tolerance_ppm: numberValue('mzPpm', 10),
      top_n_peaks: 15,
      peak_processing_preset: document.getElementById('peakProcessingPreset')?.value || 'library_match',
      peak_relative_intensity_threshold: optionalNumberValue('peakRelativeThreshold'),
      peak_precursor_exclusion_da: numberValue('peakPrecursorExclusion', 0.5),
      peak_top_n: optionalNumberValue('peakTopN'),
      peak_upper_mass_bound_enabled: document.getElementById('peakUpperMassBound')?.checked !== false,
      profile_peak_picking: document.getElementById('profilePeakPicking')?.value || 'auto',
      profile_peak_picker: document.getElementById('profilePeakPicker')?.value || 'hires',
      profile_peak_picker_signal_to_noise: numberValue('profilePeakPickerSnr', 1.0),
      enable_ms1_feature_detection: document.getElementById('enableMs1FeatureDetection')?.checked !== false,
      ms1_mass_error_ppm: numberValue('ms1MassErrorPpm', 10),
      ms2_feature_link_tolerance_ppm: numberValue('ms2FeatureLinkTolerancePpm', 20),
      ms1_noise_threshold_intensity: numberValue('ms1NoiseThreshold', 1000),
      ms1_chrom_peak_snr: numberValue('ms1ChromPeakSnr', 3),
      ms1_min_trace_length_seconds: numberValue('ms1MinTraceLength', 5),
      ms1_min_fwhm_seconds: numberValue('ms1MinFwhm', 1),
      ms1_max_fwhm_seconds: numberValue('ms1MaxFwhm', 60),
      ms1_charge_low: Math.round(numberValue('ms1ChargeLow', 1)),
      ms1_charge_high: Math.round(numberValue('ms1ChargeHigh', 1)),
      enable_adduct_grouping: document.getElementById('enableAdductGrouping')?.checked !== false,
      strict_evidence_gate: document.getElementById('strictEvidenceGate')?.checked || false,
      enable_chimeric_deconvolution: document.getElementById('enableChimericDeconvolution')?.checked !== false,
      enable_ion_identity_networking: document.getElementById('enableIonIdentityNetworking')?.checked !== false,
      ion_mode: document.getElementById(ionModeElementId)?.value || 'auto',
      enable_expected_compounds: document.getElementById('enableExpectedCompounds')?.checked === true,
      expected_compound_datasets: _selectedExpectedDatasets(),
      expected_mass_tolerance_ppm: numberValue('expectedMassTolerancePpm', 5),
      expected_min_peak_intensity: numberValue('expectedMinPeakIntensity', 1000),
      expected_min_snr: numberValue('expectedMinSnr', 1.5),
      expected_min_scan_points: Math.round(numberValue('expectedMinScanPoints', 5)),
      expected_max_gap_ratio: numberValue('expectedMaxGapRatio', 0.35),
      expected_max_fwhm_seconds: numberValue('expectedMaxFwhmSeconds', 60),
      expected_min_relative_valley_depth: numberValue('expectedMinValleyDepth', 0.10),
      expected_openms_rt_tolerance_minutes: numberValue('expectedOpenmsRtTolerance', 0.05),
      expected_background_sample_blank_ratio: numberValue('expectedBackgroundRatio', 5),
      expected_generate_common_adducts: document.getElementById('expectedGenerateCommonAdducts')?.checked !== false,
    };
    const rtMethodId = String(
      document.getElementById('rtMethodId')?.value || ''
    ).trim();
    if (rtMethodId) params.rt_method_id = rtMethodId;
    const isolationWidth = optionalNumberValue('isolationWindowWidth');
    if (isolationWidth != null && isolationWidth > 0) {
      params.isolation_window_width_da = isolationWidth;
    }
    if (params.ms1_min_fwhm_seconds > params.ms1_max_fwhm_seconds) {
      throw new Error('MS1 最小 FWHM 不能大于最大 FWHM');
    }
    if (params.ms1_charge_low > params.ms1_charge_high) {
      throw new Error('MS1 最小电荷不能大于最大电荷');
    }
    return params;
  }

  function analysisParamsForRun(ionModeElementId = 'ionMode') {
    const params = analysisParams(ionModeElementId);
    const blankEnabled = document.getElementById('enableBlank')?.checked === true;
    if (blankEnabled && !_blankPaths.length) {
      throw new Error(localizedCopy(
        '已启用空白背景评估，但本批没有选中的空白文件',
        'Blank-background assessment is enabled, but no blank is selected for this batch'
      ));
    }
    if (_blankPaths.length) {
      Object.assign(params, _blankAnalysisParameters());
    }
    return params;
  }

  // ── Analysis Start ────────────────────────────────────────────────
  async function startAnalysis() {
    if (
      _stagedFiles.length === 0
      && _vendorFolderFiles.length === 0
      && _vendorLocalSelections.length === 0
      && _retryUploadedTasks.length === 0
    ) {
      Toast.show('请先上传文件', 'warning');
      return;
    }
    if (_blankPaths.length) {
      try {
        _blankAnalysisParameters();
      } catch (error) {
        Toast.show(error.message, 'warning', 9000);
        return;
      }
    }
    if (_isUploading) return;
    _isUploading = true;

    const btn = document.getElementById('btnStartAnalysis');
    const progressWrap = document.getElementById('uploadProgressWrap');
    const progressBar = document.getElementById('uploadProgressBar');
    const progressText = document.getElementById('uploadProgressText');

    if (btn) { btn.disabled = true; btn.innerHTML = '<span class="spinner"></span> 上传中...'; }
    if (progressWrap) progressWrap.classList.add('active');

    try {
      let uploadResult = _retryUploadedTasks.length
        ? { uploaded: _retryUploadedTasks.slice(), failed: [] }
        : { uploaded: [], failed: [] };
      let importedBlankCount = 0;
      if (
        !_retryUploadedTasks.length
        && (_vendorFolderFiles.length || _vendorLocalSelections.length)
      ) {
        const accepted = document.getElementById('vendorLicenseAccepted')?.checked || false;
        if (!accepted) {
          throw new Error('使用厂商原始数据读取器前，请确认已阅读并接受 ProteoWizard 打包的厂商许可');
        }
        // Readiness is checked before any directory upload.  This call never
        // downloads because allow_container_pull is fixed to false.
        await API.post('/vendor/converter/prepare', {
          accept_vendor_licenses: true,
          allow_container_pull: false,
        });
        const peakPickingStrategy = document.getElementById('vendorPeakPicking')?.value || 'auto';
        if (_vendorLocalSelections.length) {
          const selections = _vendorLocalSelections.slice();
          const failedSelections = [];
          const convertedRows = [];
          const conversionResults = new Array(selections.length);
          const workerCount = Math.min(2, selections.length);
          let nextSelectionIndex = 0;
          let completedConversions = 0;
          const convertNextSelection = async () => {
            while (nextSelectionIndex < selections.length) {
              const index = nextSelectionIndex;
              nextSelectionIndex += 1;
              const selection = selections[index];
              if (progressText) {
                progressText.textContent = (
                  `正在本机转换 ${completedConversions}/${selections.length}；`
                  + `${workerCount} 路受控并行，每个文件独立校验 MS1/MS2…`
                );
              }
              try {
                conversionResults[index] = {
                  selection,
                  result: await API.post('/vendor/convert-local', {
                    selection_token: selection.selection_token,
                    accept_vendor_licenses: true,
                    allow_container_pull: false,
                    backend: 'auto',
                    intake_role: selection.intake_role === 'blank' ? 'blank' : 'sample',
                    blank_role: 'background_blank',
                    peak_picking_strategy: peakPickingStrategy,
                  }),
                };
              } catch (error) {
                conversionResults[index] = { selection, error };
              } finally {
                completedConversions += 1;
                if (progressBar) {
                  progressBar.style.width = `${Math.max(
                    5,
                    Math.round((completedConversions / selections.length) * 100)
                  )}%`;
                }
              }
            }
          };
          await Promise.all(
            Array.from({ length: workerCount }, () => convertNextSelection())
          );
          for (const item of conversionResults) {
            const selection = item.selection;
            if (item.error) {
              failedSelections.push(selection);
              uploadResult.failed.push({
                filename: selection.raw_name,
                error: item.error.message,
              });
              continue;
            }
            const vendorResult = item.result || {};
            const uploadedRows = vendorResult.uploaded || [];
            const blankRows = vendorResult.blank_uploads || [];
            importedBlankCount += blankRows.length;
            uploadResult.uploaded.push(...uploadedRows);
            uploadResult.failed.push(...(vendorResult.failed || []));
            convertedRows.push(
              ...uploadedRows,
              ...blankRows.map(row => ({
                filename: row.blank_filename,
                ms1_count: row.n_ms1_scans,
                ms2_count: row.n_ms2_contexts,
                intake_role: 'blank',
              }))
            );
            if (blankRows.length) {
              _mergeAvailableBlankUploads(blankRows);
              const selectedIdentities = new Set(
                _blankUploads.map(row => _blankContentIdentity(row))
              );
              for (const row of blankRows) {
                const identity = _blankContentIdentity(row);
                if (
                  _blankUploads.length < 10
                  && _blankSupportsMode(row)
                  && !selectedIdentities.has(identity)
                ) {
                  _blankUploads.push(row);
                  selectedIdentities.add(identity);
                }
              }
              _syncBlankPaths();
              _renderBlankStatuses();
            }
          }
          _vendorLocalSelections = failedSelections;
          _vendorRootName = failedSelections.length
            ? `${failedSelections.length} 个待重试厂商原始样本`
            : '';
          renderVendorFolder();
          if (progressBar) progressBar.style.width = '100%';
          if (convertedRows.length) {
            const totalMs1 = convertedRows.reduce(
              (sum, row) => sum + Number(row.ms1_count || 0),
              0
            );
            const totalMs2 = convertedRows.reduce(
              (sum, row) => sum + Number(row.ms2_count || 0),
              0
            );
            Toast.show(
              `厂商原始数据转换通过 ${convertedRows.length}/${selections.length}：`
                + `MS1 ${totalMs1.toLocaleString()} / MS2 ${totalMs2.toLocaleString()}`,
              failedSelections.length ? 'warning' : 'success'
            );
          }
        } else {
          if (progressText) progressText.textContent = '上传 Agilent .d 原始目录…';
          const vendorResult = await API.convertVendorFolder(
            _vendorFolderFiles,
            {
              acceptVendorLicenses: true,
              allowContainerPull: false,
              backend: 'auto',
              peakPickingStrategy,
            },
            pct => {
              if (progressBar) progressBar.style.width = pct + '%';
              if (progressText) {
                progressText.textContent = pct < 100
                  ? `Agilent .d 上传进度: ${pct}%`
                  : '原始目录已上传，ProteoWizard 正在转换并校验 MS1/MS2…';
              }
            }
          );
          uploadResult.uploaded.push(...(vendorResult.uploaded || []));
          uploadResult.failed.push(...(vendorResult.failed || []));
          const converted = (vendorResult.uploaded || [])[0];
          if (converted) {
            Toast.show(
              `.d 转换通过：MS1 ${Number(converted.ms1_count || 0).toLocaleString()} / MS2 ${Number(converted.ms2_count || 0).toLocaleString()}`
                + `${converted.isolation_width_da ? `；方法隔离窗 ${converted.isolation_width_da} Da` : ''}`,
              'success'
            );
          }
        }
      }
      // Convert the all-or-nothing proprietary directory first.  Ordinary
      // uploads follow only after it succeeds, avoiding stranded pending tasks.
      if (!_retryUploadedTasks.length && _stagedFiles.length) {
        const ordinaryResult = await API.upload(_stagedFiles, (pct) => {
          if (progressBar) progressBar.style.width = pct + '%';
          if (progressText) progressText.textContent = `开放格式上传进度: ${pct}%`;
        });
        uploadResult.uploaded.push(...(ordinaryResult.uploaded || []));
        uploadResult.failed.push(...(ordinaryResult.failed || []));
      }

      // Build one complete request for every entry point.  The queue page
      // uses the same helper, so a selected blank cannot disappear when a
      // pending task is started from "Run Queue" instead of this page.
      const params = analysisParamsForRun('ionMode');

      // One batch request gives the whole upload a single, explicit queue result.
      const uploaded = uploadResult.uploaded || [];
      const uploadFailed = uploadResult.failed || [];
      if (uploaded.length === 0) {
        if (importedBlankCount > 0 && uploadFailed.length === 0) {
          Toast.show(
            `已导入 ${importedBlankCount} 个厂商 RAW 空白；可继续选择样本并使用这些空白`,
            'success'
          );
          return;
        }
        throw new Error('服务器未接收任何可分析样本');
      }
      // Persist the server-side task IDs before queueing. If the batch request
      // itself fails (network/HTTP error), the next click can retry these IDs
      // without uploading the sample files again.
      _retryUploadedTasks = uploaded.slice();
      const batchResult = await API.post('/tasks/batch/analyze', {
        task_ids: uploaded.map(u => u.task_id),
        params,
      });
      const queued = batchResult.queued_count || 0;
      const skipped = (batchResult.skipped || []).length;
      const queueFailed = (batchResult.failed || []).length;
      const queueFailureItems = batchResult.failed || [];
      const failedNames = uploadFailed.slice(0, 3).map(item => item.filename).join('、');
      const queueFailureDetail = queueFailureItems
        .slice(0, 2)
        .map(item => item.error)
        .filter(Boolean)
        .join('；');
      Toast.show(
        `批量分析已排队：${queued}/${uploaded.length} 个样本`
          + `${_blankPaths.length && String(params.ion_mode || '').toUpperCase() === 'AUTO' ? '；正/负空白将按每个样本实际极性分配（同极性不等于同批次／同方法配对）' : ''}`
          + `${uploadFailed.length ? `，上传拒绝 ${uploadFailed.length} 个${failedNames ? `（${failedNames}）` : ''}` : ''}`
          + `${skipped ? `，跳过 ${skipped}` : ''}${queueFailed ? `，排队失败 ${queueFailed}` : ''}`
          + `${queueFailureDetail ? `：${queueFailureDetail}` : ''}`,
        (uploadFailed.length || queueFailed || queued === 0) ? 'warning' : 'success'
      );

      const acceptedTaskIds = new Set([
        ...(batchResult.queued || []),
        ...(batchResult.skipped || []).map(item => item.id),
      ]);
      const queueIncomplete = acceptedTaskIds.size < uploaded.length;
      _retryUploadedTasks = uploaded.filter(
        item => !acceptedTaskIds.has(item.task_id)
      );
      void refreshBlankUploads();
      if (queueFailed || queueIncomplete || uploadFailed.length) {
        if (uploadFailed.length && !_retryUploadedTasks.length) {
          const rejectedNames = new Set(uploadFailed.map(item => item.filename));
          _stagedFiles = _stagedFiles.filter(file => rejectedNames.has(file.name));
          renderFileList();
        }
        // Keep the form and current blank selection.  A mode/library correction
        // can now be retried immediately; already-uploaded task IDs are reused.
        updateUI();
        return;
      }

      // Every uploaded task is now queued or was already active.
      clearFiles();
      // Keep the explicit blank selection for the next sample batch.  Users
      // can still detach individual blanks or disable blank processing; a
      // successful queue must not silently change the next batch's evidence.
      // Navigate to tasks page
      App.navigateTo('tasks');

    } catch (e) {
      Toast.show(`上传或批量排队失败: ${e.message}`, 'error');
    } finally {
      _isUploading = false;
      if (btn) btn.disabled = false;
      if (progressWrap) { progressWrap.classList.remove('active'); }
      if (progressBar) { progressBar.style.width = '0%'; }
      updateUI();
    }
  }

  // ── UI State Management ───────────────────────────────────────────
  function updateUI() {
    const hasFiles = (
      _stagedFiles.length > 0
      || _vendorFolderFiles.length > 0
      || _vendorLocalSelections.length > 0
      || _retryUploadedTasks.length > 0
    );
    const settingsCard = document.getElementById('settingsCard');
    const btnStart = document.getElementById('btnStartAnalysis');
    const runReadiness = document.getElementById('analysisRunReadiness');
    const libLoaded = !!document.querySelector('.lib-status.loaded');

    if (settingsCard) {
      settingsCard.classList.remove('hidden');
      settingsCard.classList.toggle('is-waiting-for-sample', !hasFiles);
    }
    if (btnStart) btnStart.disabled = !(hasFiles && libLoaded);
    if (runReadiness) {
      runReadiness.textContent = !libLoaded
        ? '请先加载参考谱库；分析级别和全部阈值仍可预先设置。'
        : !hasFiles
        ? '分析参数已可设置；选择样本后即可开始分析。'
        : '样本与谱库已就绪；请确认分析级别后开始。';
      runReadiness.classList.toggle('is-ready', hasFiles && libLoaded);
    }
    if (btnStart && !_isUploading) {
      const inputCount = _retryUploadedTasks.length || (
        _stagedFiles.length
        + (_vendorFolderFiles.length ? 1 : _vendorLocalSelections.length)
      );
      btnStart.textContent = _retryUploadedTasks.length
        ? `↻ 重试排队（${inputCount} 个已上传样本）`
        : inputCount > 1
        ? `🚀 批量开始分析（${inputCount} 个样本）`
        : hasFiles ? '🚀 开始分析（1 个样本）' : '🚀 开始分析';
    }
    _updateWorkflowRail({ hasFiles, libLoaded });
  }

  function _updateWorkflowRail({ hasFiles, libLoaded }) {
    if (typeof document.querySelectorAll !== 'function') return;
    const steps = Array.from(
      document.querySelectorAll('.intake-step-rail-item')
    );
    if (steps.length < 4) return;
    steps.forEach(step => step.classList.remove('is-ready', 'is-current'));
    if (libLoaded) {
      steps[0].classList.add('is-ready');
    } else {
      steps[0].classList.add('is-current');
    }
    if (hasFiles) {
      steps[1].classList.add('is-ready');
    } else if (libLoaded) {
      steps[1].classList.add('is-current');
    }
    if (_blankUploads.length) {
      steps[2].classList.add('is-ready');
    }
    if (hasFiles && libLoaded) {
      steps[3].classList.add('is-current');
    }
  }

  // ── Helpers ───────────────────────────────────────────────────────
  function escapeHtml(str) {
    const div = document.createElement('div');
    div.textContent = str;
    return div.innerHTML;
  }

  // ── Expected Compounds target catalogues ──────────────────────────
  function toggleExpectedCompounds(enabled) {
    const controls = document.getElementById('expectedCompoundsControls');
    if (controls) controls.style.display = enabled ? 'block' : 'none';
    if (enabled && !_expectedDatasetsLoaded) {
      void loadExpectedCompoundDatasets();
    }
  }

  function _selectedExpectedDatasets() {
    if (typeof document.querySelectorAll !== 'function') return [];
    return Array.from(
      document.querySelectorAll('#expectedDatasetList input[name="expectedDataset"]:checked')
    ).map(input => input.value);
  }

  function updateExpectedDatasetStatus() {
    const status = document.getElementById('expectedDatasetStatus');
    if (!status) return;
    const selected = new Set(_selectedExpectedDatasets());
    const totalRows = _expectedDatasets.reduce(
      (sum, item) => sum + Number(item.target_rows || 0), 0
    );
    if (!selected.size) {
      status.textContent = _expectedDatasets.length
        ? `尚未选择目标库 · 靶向层不会运行（可用 ${_expectedDatasets.length} 个）`
        : '当前没有可用的 Expected Compounds 数据集';
      status.style.color = '#b45309';
      return;
    }
    const selectedRows = _expectedDatasets
      .filter(item => selected.has(String(item.dataset_source || '')))
      .reduce((sum, item) => sum + Number(item.target_rows || 0), 0);
    status.textContent = `已限定 ${selected.size}/${_expectedDatasets.length} 个数据集 · ${selectedRows.toLocaleString()} 条目标行`;
    status.style.color = '#1d4ed8';
  }

  function clearExpectedDatasetSelection() {
    if (typeof document.querySelectorAll !== 'function') return;
    document.querySelectorAll(
      '#expectedDatasetList input[name="expectedDataset"]'
    ).forEach(input => { input.checked = false; });
    updateExpectedDatasetStatus();
  }

  async function loadExpectedCompoundDatasets() {
    const list = document.getElementById('expectedDatasetList');
    const status = document.getElementById('expectedDatasetStatus');
    if (!list) return;
    if (status) {
      status.textContent = '正在读取可用数据集…';
      status.style.color = 'var(--text-secondary)';
    }
    try {
      const payload = await API.get('/expected-compounds/datasets');
      _expectedDatasets = Array.isArray(payload.datasets) ? payload.datasets : [];
      _expectedDatasetsLoaded = true;
      if (!_expectedDatasets.length) {
        list.innerHTML = '<span style="font-size:0.72rem;color:#b45309;">生产库中尚无 literature_presence_only 目标表。</span>';
        updateExpectedDatasetStatus();
        return;
      }
      list.innerHTML = _expectedDatasets.map((item, index) => {
        const source = String(item.dataset_source || `dataset_${index + 1}`);
        const rowCount = Number(item.target_rows || 0);
        const entityCount = Number(item.named_formula_entities || 0);
        return `
          <label style="display:flex;align-items:flex-start;gap:6px;cursor:pointer;min-width:0;">
            <input type="checkbox" name="expectedDataset" value="${escapeAttribute(source)}" style="width:auto;margin-top:2px;">
            <span style="min-width:0;font-size:0.72rem;line-height:1.35;">
              <strong data-i18n-skip style="display:block;overflow-wrap:anywhere;">${escapeHtml(source)}</strong>
              <span style="color:var(--text-secondary);">${rowCount.toLocaleString()} 行 · ${entityCount.toLocaleString()} 个名称/分子式实体</span>
            </span>
          </label>`;
      }).join('');
      list.querySelectorAll('input[name="expectedDataset"]').forEach(input => {
        input.addEventListener('change', updateExpectedDatasetStatus);
      });
      updateExpectedDatasetStatus();
    } catch (error) {
      _expectedDatasets = [];
      list.innerHTML = `<span style="font-size:0.72rem;color:#b91c1c;">数据集读取失败：${escapeHtml(error.message)}</span>`;
      if (status) {
        status.textContent = 'Expected Compounds 数据集不可用';
        status.style.color = '#b91c1c';
      }
    }
  }

  function escapeAttribute(str) {
    return escapeHtml(str)
      .replaceAll('"', '&quot;')
      .replaceAll("'", '&#39;');
  }

  async function refreshRuntimeIdentity() {
    const badge = document.getElementById('runtimeInstanceBadge');
    if (!badge) return;
    const currentPort = window.location.port || (window.location.protocol === 'https:' ? '443' : '80');
    try {
      const health = await API.get('/health');
      const generation = String(health.service_generation || '').toUpperCase();
      const expectedGeneration = 'V16';
      if (generation !== expectedGeneration) {
        badge.textContent = generation
          ? `⚠ 当前端口 ${currentPort} 运行的是 ${generation}，不是 MSsystem ${expectedGeneration}`
          : `⚠ 当前端口 ${currentPort} 不是可识别的 MSsystem 服务`;
        badge.style.background = '#fff5f5';
        badge.style.color = '#c53030';
        return;
      }
      _blankEvidenceV2Available = health.capabilities?.blank_evidence_v2 === true;
      const blankPolicy = document.getElementById('blankFilterPolicy');
      const blankPolicyStatus = document.getElementById('blankPolicyRuntimeStatus');
      if (blankPolicy) {
        blankPolicy.disabled = !_blankEvidenceV2Available;
      }
      if (blankPolicyStatus) {
        blankPolicyStatus.textContent = _blankEvidenceV2Available
          ? '✓ 当前服务支持实际特征级背景排除与仅比较两种模式'
          : '⚠ 当前进程仍是旧空白后端；带空白的新任务已暂时锁定，等待安全重启';
        blankPolicyStatus.style.color = _blankEvidenceV2Available ? '#276749' : '#b45309';
      }
      updateBlankPolicyUI();
      _renderBlankStatuses();
      const custom = String(currentPort) !== String(health.default_port || 8516);
      badge.textContent = `✓ 运行实例：MSsystem · I.D.E.A platform · 单一客户端/API 端口 ${currentPort}${custom ? '（自定义端口）' : ''}`;
      badge.style.background = '#f0fff4';
      badge.style.color = '#276749';
    } catch (e) {
      _blankEvidenceV2Available = false;
      badge.textContent = `⚠ 无法核验 MSsystem 服务（当前端口 ${currentPort}）`;
      badge.style.background = '#fffaf0';
      badge.style.color = '#975a16';
    }
  }

  // ── 动态库规模 (下拉选项真实条数) ──────────────────────────────────
  function _fmtCount(n) {
    if (n == null) return '?';
    if (n >= 1000000) return (n / 1000000).toFixed(2).replace(/0+$/, '').replace(/\.$/, '') + 'M';
    if (n >= 1000) return (n / 1000).toFixed(n >= 100000 ? 0 : 1).replace(/\.0$/, '') + 'K';
    return String(n);
  }

  function _dbCountLabel(prefix, entry) {
    if (!entry.available) return `${prefix} (未构建)`;
    const physical = entry.physical_spectra ?? entry.spectra;
    const loaded = entry.loaded_searchable_spectra;
    const eligible = entry.source_searchable_spectra;
    const projected = entry.runtime_searchable_spectra;
    const searchable = loaded != null ? loaded : (projected ?? eligible);
    const searchableNote = loaded != null ? '当前已载入' : '运行时可检索';
    return `${prefix} (${_fmtCount(searchable)} ${searchableNote} / ${_fmtCount(eligible)} 来源合格 / ${_fmtCount(physical)} 物理记录；${_fmtCount(entry.compounds)} 化合物)`;
  }

  async function refreshDbCounts() {
    const sel = document.getElementById('dbMode');
    if (!sel) return;
    try {
      const counts = await API.get('/library/counts');
      const small = counts.small || {};
      const full = counts.full || {};
      const sLabel = _dbCountLabel('Small <800 Da', small);
      const fLabel = _dbCountLabel('Full', full);
      // 保留当前选中值
      const cur = sel.value;
      const smallOption = sel.querySelector('option[value="small"]');
      const fullOption = sel.querySelector('option[value="full"]');
      if (smallOption) smallOption.textContent = sLabel;
      if (fullOption) fullOption.textContent = fLabel;
      sel.value = _loadedDbMode || cur;
    } catch (e) {
      console.warn('refreshDbCounts failed:', e);
    }
  }

  // ── Blank evidence workflow ────────────────────────────────────────
  function updateBlankPolicyUI() {
    const policy = String(
      document.getElementById('blankFilterPolicy')?.value || 'hard_legacy'
    );
    const hard = policy === 'hard_legacy';
    const confirmations = document.getElementById(
      'blankHardExclusionConfirmations'
    );
    const explanation = document.getElementById('blankPolicyExplanation');
    if (confirmations) confirmations.style.display = hard ? 'block' : 'none';
    if (explanation) {
      explanation.textContent = hard
        ? '实际背景排除下：样品/空白低于设定比值的同极性特征在检索前排除并计入“已扣除”；原始 MS/MS 与排除原因保留用于追溯。系统不跨针逐峰相减碎片强度。'
        : '仅比较模式明确不会扣除任何 MS2；系统只记录样品富集、空白相关、空白占优和 carryover 证据，结果中的“已扣除”将为 0。';
    }
    _renderBlankStatuses();
  }

  function _blankAnalysisParameters() {
    if (!_blankPaths.length) return {};
    if (!_blankEvidenceV2Available) {
      throw new Error(localizedCopy(
        '当前服务尚未启用新的空白背景处理后端，带空白任务不能排队。',
        'The active service does not expose the blank-evidence backend; a run with blanks cannot be queued.'
      ));
    }
    const policy = String(
      document.getElementById('blankFilterPolicy')?.value || 'hard_legacy'
    );
    const blankRoles = Object.fromEntries(
      _blankUploads.map(row => [
        String(row.blank_path),
        String(row.blank_role || 'background_blank'),
      ])
    );
    const hard = policy === 'hard_legacy';
    if (hard) {
      const invalid = _blankUploads.filter(
        row => !BACKGROUND_EXCLUSION_ROLES.has(
          String(row.blank_role || 'background_blank')
        )
      );
      if (invalid.length) {
        const names = invalid.slice(0, 4).map(
          row => row.blank_filename || row.blank_path
        ).join('、');
        throw new Error(localizedCopy(
          `实际空白排除要求把空白类型设为“背景空白/溶剂空白/方法空白/提取空白”；请先修改：${names}`,
          `Feature-level blank exclusion requires background, solvent, method, or extraction blank roles. Update: ${names}`
        ));
      }
      if (!document.getElementById('blankMethodComparable')?.checked) {
        throw new Error(localizedCopy(
          '已选择实际空白排除，但尚未确认空白与样品同批次、同色谱方法',
          'Feature-level blank exclusion is selected, but same-batch and same-LC-method comparability has not been confirmed'
        ));
      }
      if (!document.getElementById('blankExclusionConfirmed')?.checked) {
        throw new Error(localizedCopy(
          '已选择实际空白排除，但尚未确认从主检索中排除空白占优特征',
          'Feature-level blank exclusion is selected, but removal of blank-dominant features from the primary search has not been confirmed'
        ));
      }
    }
    return {
      blank_files: _blankPaths.slice(),
      blank_roles: blankRoles,
      blank_filter_policy: hard ? 'hard_legacy' : 'soft_evidence',
      blank_method_comparability_confirmed: hard,
      blank_background_exclusion_confirmed: hard,
      blank_sample_ratio: numberValue('blankSampleRatio', 3),
      blank_rt_tolerance_minutes: numberValue('blankRtTolerance', 0.30),
      blank_min_detection_fraction: numberValue(
        'blankMinDetectionFraction', 0.50
      ),
    };
  }

  function toggleBlank(enabled) {
    const ctrl = document.getElementById('blankControls');
    if (ctrl) ctrl.style.display = enabled ? 'block' : 'none';
    if (!enabled) {
      _resetBlankSelection();
    } else {
      updateBlankPolicyUI();
    }
  }

  function _currentIonMode() {
    return String(document.getElementById('ionMode')?.value || 'auto').toUpperCase();
  }

  function _blankPolarities(row) {
    return Array.isArray(row?.polarities)
      ? row.polarities.map(value => String(value).toUpperCase())
      : [];
  }

  function _blankContentIdentity(row) {
    const signatures = row?.vendor_conversion?.source?.required_file_signatures;
    if (signatures && typeof signatures === 'object') {
      const hashes = Object.keys(signatures).sort().map(key =>
        String(signatures[key]?.sha256 || '').trim().toLowerCase()
      ).filter(Boolean);
      if (hashes.length) return `vendor-source:${hashes.join(':')}`;
    }
    const directHash = String(
      row?.content_sha256
      || row?.sha256
      || row?.vendor_conversion?.output?.sha256
      || ''
    ).trim().toLowerCase();
    if (directHash) return `content:${directHash}`;
    return `token:${String(row?.blank_path || '')}`;
  }

  function _dedupeBlankUploads(rows) {
    const seen = new Set();
    return (rows || []).filter(row => {
      const identity = _blankContentIdentity(row);
      if (seen.has(identity)) return false;
      seen.add(identity);
      return true;
    });
  }

  function _blankSupportsMode(row, mode = _currentIonMode()) {
    const normalizedMode = String(mode || 'AUTO').toUpperCase();
    if (normalizedMode === 'AUTO') return true;
    return _blankPolarities(row).includes(normalizedMode);
  }

  function _blankRole(row) {
    const role = String(row?.blank_role || 'unspecified');
    return Object.prototype.hasOwnProperty.call(BLANK_ROLE_LABELS, role)
      ? role : 'unspecified';
  }

  function _blankRoleOptions(selectedRole) {
    return Object.entries(BLANK_ROLE_LABELS).map(([value, label]) =>
      `<option value="${value}"${value === selectedRole ? ' selected' : ''}>${escapeHtml(label)}</option>`
    ).join('');
  }

  function _persistBlankSelection() {
    try {
      sessionStorage.setItem(
        BLANK_SELECTION_STORAGE_KEY,
        JSON.stringify(_blankUploads.map(row => row.blank_path))
      );
    } catch (error) {
      console.warn('Blank selection could not be persisted:', error);
    }
  }

  function _savedBlankSelection() {
    try {
      const tokens = JSON.parse(sessionStorage.getItem(BLANK_SELECTION_STORAGE_KEY) || '[]');
      return Array.isArray(tokens) ? tokens.map(String) : [];
    } catch {
      return [];
    }
  }

  function _mergeAvailableBlankUploads(rows) {
    const byToken = new Map(
      _availableBlankUploads.map(row => [String(row.blank_path), row])
    );
    for (const row of rows || []) {
      if (row?.blank_path) byToken.set(String(row.blank_path), row);
    }
    _availableBlankUploads = Array.from(byToken.values());
  }

  async function refreshBlankUploads() {
    const manager = document.getElementById('blankUploadManager');
    try {
      const response = await API.get('/blanks');
      _availableBlankUploads = Array.isArray(response.blanks) ? response.blanks : [];
      const deletableTokens = new Set(
        _availableBlankUploads
          .filter(row => !row.in_use)
          .map(row => String(row.blank_path))
      );
      _blankManagerDeleteSelection = new Set(
        Array.from(_blankManagerDeleteSelection)
          .filter(token => deletableTokens.has(token))
      );
      const selectedTokens = _blankUploads.length
        ? _blankUploads.map(row => String(row.blank_path))
        : _savedBlankSelection();
      const availableByToken = new Map(
        _availableBlankUploads.map(row => [String(row.blank_path), row])
      );
      _blankUploads = selectedTokens
        .map(token => availableByToken.get(token))
        .filter(Boolean)
        .slice(0, 10);
      _syncBlankPaths();
      _renderBlankStatuses();
    } catch (error) {
      if (manager) {
        manager.innerHTML = `<span style="font-size:0.72rem;color:var(--danger);">空白文件清单读取失败：${escapeHtml(error.message)}</span>`;
      }
    }
  }

  async function _uploadBlankFiles(fileList) {
    const uploaded = [];
    try {
      for (const file of Array.from(fileList || [])) {
        const fd = new FormData();
        fd.append('file', file);
        uploaded.push(await API.post('/blank/upload', fd));
      }
    } catch (error) {
      // Keep one failed multi-file selection atomic from the user's
      // perspective; best-effort delete any files accepted earlier in it.
      await Promise.allSettled(
        uploaded.map(row => API.del(`/blanks/${encodeURIComponent(row.blank_path)}`))
      );
      throw error;
    }
    return uploaded;
  }

  function _syncBlankPaths() {
    _blankUploads = _dedupeBlankUploads(_blankUploads).slice(0, 10);
    _blankPaths = _blankUploads.map(row => row.blank_path);
    _persistBlankSelection();
  }

  function _resetBlankInputs() {
    ['blankFile', 'blankFileQuick'].forEach(id => {
      const input = document.getElementById(id);
      if (input) input.value = '';
    });
  }

  function _resetBlankSelection() {
    _blankUploads = [];
    _syncBlankPaths();
    _resetBlankInputs();
    const enabled = document.getElementById('enableBlank');
    if (enabled) enabled.checked = false;
    const controls = document.getElementById('blankControls');
    if (controls) controls.style.display = 'none';
    const quickStatus = document.getElementById('blankStatusQuick');
    const settingsStatus = document.getElementById('blankStatus');
    if (quickStatus) quickStatus.innerHTML = '';
    if (settingsStatus) settingsStatus.innerHTML = '';
    const clearButton = document.getElementById('btnClearBlankUploads');
    if (clearButton) clearButton.style.display = 'none';
    _renderBlankManager();
  }

  function _blankReadyHtml(uploaded) {
    const totalMs1 = uploaded.reduce((sum, row) => sum + Number(row.n_ms1_scans || 0), 0);
    const totalMs2 = uploaded.reduce((sum, row) => sum + Number(row.n_ms2_contexts || 0), 0);
    const mode = _currentIonMode();
    const hard = String(
      document.getElementById('blankFilterPolicy')?.value || 'hard_legacy'
    ) === 'hard_legacy';
    const pairingNote = mode === 'AUTO'
      ? '自动模式会逐个样本识别 POS/NEG，只把同极性空白分配给该样本。实际排除还要求确认同批次／同方法。'
      : `当前固定为 ${mode}，只会使用 ${mode} 空白。`;
    const rows = uploaded.map((row, index) => {
      const polarities = Array.isArray(row.polarities) && row.polarities.length
        ? row.polarities.join('/') : '极性未知';
      const role = _blankRole(row);
      return `
        <div style="display:flex;align-items:center;justify-content:space-between;gap:8px;margin-top:4px;padding:4px 7px;border:1px solid var(--border-color);border-radius:6px;background:#fff;">
          <span style="min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${escapeAttribute(row.blank_filename)}">
            ${escapeHtml(row.blank_filename)} <strong>[${escapeHtml(polarities)}]</strong>
            <span style="color:var(--text-secondary);"> · ${escapeHtml(BLANK_ROLE_LABELS[role])}</span>
          </span>
          <span style="display:flex;gap:5px;flex-shrink:0;">
            <button type="button" class="btn btn-secondary btn-sm" onclick="Upload.detachBlank(${index})">本批不用</button>
            <button type="button" class="btn btn-secondary btn-sm" onclick="Upload.removeBlank(${index})">删除文件</button>
          </span>
        </div>`;
    }).join('');
    return `
      <div>✓ 空白就绪：<strong>${uploaded.length} 个</strong>；MS1 ${totalMs1.toLocaleString()} 扫描 / MS2 ${totalMs2.toLocaleString()} 张</div>
      <div>${rows}</div>
      <div style="margin-top:4px;color:#276749;"><strong>极性配对：</strong>${escapeHtml(pairingNote)}</div>
      <div style="margin-top:4px;color:var(--text-secondary);">${hard
        ? '当前选择：实际特征级背景排除；确认项未勾选或空白类型不合格时系统将拒绝排队，不会退化成 0 扣除。'
        : '当前选择：仅背景比较，明确不扣除任何 MS2。'} carryover 空白只作残留证据。提交后继续用于下一批，直到你点击“本批不用”或关闭空白处理。</div>`;
  }

  function _renderBlankManager() {
    const manager = document.getElementById('blankUploadManager');
    if (!manager) return;
    if (!_availableBlankUploads.length) {
      _blankManagerDeleteSelection.clear();
      manager.innerHTML = '<span style="font-size:0.72rem;color:var(--text-secondary);">服务器当前没有可管理的空白文件。</span>';
      return;
    }
    const selectedTokens = new Set(_blankUploads.map(row => String(row.blank_path)));
    const mode = _currentIonMode();
    const deletableRows = _availableBlankUploads.filter(row => !row.in_use);
    const selectedDeleteCount = deletableRows.filter(
      row => _blankManagerDeleteSelection.has(String(row.blank_path))
    ).length;
    const allDeletableSelected = (
      deletableRows.length > 0
      && selectedDeleteCount === deletableRows.length
    );
    const rows = _availableBlankUploads.map((row, index) => {
      const token = String(row.blank_path || '');
      const selected = selectedTokens.has(token);
      const checkedForDelete = _blankManagerDeleteSelection.has(token);
      const compatible = _blankSupportsMode(row, mode);
      const polarities = _blankPolarities(row);
      const polarityLabel = polarities.length ? polarities.join('/') : '极性未知';
      const role = _blankRole(row);
      const inUse = !!row.in_use;
      const selectionButton = selected
        ? `<button type="button" class="btn btn-secondary btn-sm" onclick="Upload.detachManagedBlank(${index})">从本批移除</button>`
        : `<button type="button" class="btn btn-secondary btn-sm" onclick="Upload.selectManagedBlank(${index})"${compatible ? '' : ' disabled'}>用于本批</button>`;
      return `
        <div style="display:flex;align-items:center;justify-content:space-between;gap:8px;margin-top:5px;padding:6px 8px;border:1px solid var(--border-color);border-radius:7px;background:#fff;flex-wrap:wrap;">
          <span style="display:flex;align-items:flex-start;gap:7px;min-width:0;font-size:0.74rem;">
            <input type="checkbox" aria-label="选择删除 ${escapeAttribute(row.blank_filename || token)}"
                   onchange="Upload.toggleBlankDeleteSelection(${index}, this.checked)"
                   ${checkedForDelete ? 'checked' : ''}${inUse ? ' disabled' : ''}
                   title="${inUse ? '正在被任务使用，不能删除' : '加入批量删除'}"
                   style="width:auto;margin-top:2px;flex-shrink:0;">
            <span>
            <strong>${escapeHtml(row.blank_filename || token)}</strong>
            <span style="color:var(--text-secondary);"> · ${escapeHtml(polarityLabel)} · ${Number(row.size_mb || 0).toFixed(2)} MB</span>
            ${selected ? '<span style="color:#276749;"> · 本批已选</span>' : ''}
            ${inUse ? '<span style="color:#975a16;"> · 正被排队/运行任务使用</span>' : ''}
            ${!compatible && mode !== 'AUTO' ? `<span style="color:var(--danger);"> · 与 ${escapeHtml(mode)} 不兼容</span>` : ''}
            </span>
          </span>
          <span style="display:flex;gap:5px;flex-shrink:0;align-items:center;flex-wrap:wrap;">
            <select class="form-select" aria-label="空白类型"
                    onchange="Upload.updateBlankRole(${index}, this.value)"
                    ${_blankEvidenceV2Available ? '' : 'disabled'}
                    style="width:auto;max-width:190px;padding:2px 5px;font-size:.7rem;">
              ${_blankRoleOptions(role)}
            </select>
            ${selectionButton}
            <button type="button" class="btn btn-secondary btn-sm" onclick="Upload.deleteManagedBlank(${index})"${inUse ? ' disabled' : ''}>删除</button>
          </span>
        </div>`;
    }).join('');
    manager.innerHTML = `
      <div style="display:flex;align-items:center;justify-content:space-between;gap:8px;flex-wrap:wrap;">
        <div style="font-size:0.72rem;color:var(--text-secondary);">已上传空白文件（请注明类型；进样后 carryover 只作残留证据）</div>
        <div style="display:flex;align-items:center;gap:8px;flex-wrap:wrap;">
          <label style="font-size:0.72rem;display:flex;align-items:center;gap:4px;cursor:${deletableRows.length ? 'pointer' : 'not-allowed'};">
            <input type="checkbox" onchange="Upload.toggleSelectAllBlankDeletes(this.checked)"
                   ${allDeletableSelected ? 'checked' : ''}${deletableRows.length ? '' : ' disabled'}
                   style="width:auto;">
            ${allDeletableSelected ? '取消全选' : '全选可删除'}（${deletableRows.length}）
          </label>
          <span style="font-size:0.72rem;color:var(--text-secondary);">已选 ${selectedDeleteCount}</span>
          <button type="button" class="btn btn-danger btn-sm"
                  onclick="Upload.deleteSelectedBlanks()"${selectedDeleteCount ? '' : ' disabled'}>
            🗑 删除所选 (${selectedDeleteCount})
          </button>
        </div>
      </div>
      ${rows}`;
  }

  function toggleBlankDeleteSelection(index, checked) {
    const row = _availableBlankUploads[index];
    if (!row || row.in_use) return;
    const token = String(row.blank_path);
    if (checked) {
      _blankManagerDeleteSelection.add(token);
    } else {
      _blankManagerDeleteSelection.delete(token);
    }
    _renderBlankManager();
  }

  function toggleSelectAllBlankDeletes(checked) {
    if (checked) {
      _availableBlankUploads
        .filter(row => !row.in_use)
        .forEach(row => {
          _blankManagerDeleteSelection.add(String(row.blank_path));
        });
    } else {
      _blankManagerDeleteSelection.clear();
    }
    _renderBlankManager();
  }

  async function updateBlankRole(index, blankRole) {
    const row = _availableBlankUploads[index];
    if (!row) return;
    if (!_blankEvidenceV2Available) {
      Toast.show('新的空白证据后端等待安全重启，当前不能修改空白类型', 'warning');
      return;
    }
    try {
      const updated = await API.post(
        `/blanks/${encodeURIComponent(row.blank_path)}/role`,
        { blank_role: blankRole }
      );
      _availableBlankUploads[index] = updated;
      _blankUploads = _blankUploads.map(item =>
        item.blank_path === updated.blank_path ? updated : item
      );
      _syncBlankPaths();
      _renderBlankStatuses();
      Toast.show(`空白类型已设为：${BLANK_ROLE_LABELS[_blankRole(updated)]}`, 'success');
    } catch (error) {
      Toast.show(`空白类型更新失败：${error.message}`, 'error');
      await refreshBlankUploads();
    }
  }

  function _renderBlankStatuses() {
    if (!_blankUploads.length) {
      const quickStatus = document.getElementById('blankStatusQuick');
      const settingsStatus = document.getElementById('blankStatus');
      if (quickStatus) quickStatus.innerHTML = '';
      if (settingsStatus) settingsStatus.innerHTML = '';
      const clearButton = document.getElementById('btnClearBlankUploads');
      if (clearButton) clearButton.style.display = 'none';
      const enabled = document.getElementById('enableBlank');
      if (enabled) enabled.checked = false;
      const controls = document.getElementById('blankControls');
      if (controls) controls.style.display = 'none';
      _renderBlankManager();
      _updateOptionalQcState();
      updateUI();
      return;
    }
    const html = _blankReadyHtml(_blankUploads);
    const quickStatus = document.getElementById('blankStatusQuick');
    const settingsStatus = document.getElementById('blankStatus');
    if (quickStatus) quickStatus.innerHTML = html;
    if (settingsStatus) settingsStatus.innerHTML = html;
    const clearButton = document.getElementById('btnClearBlankUploads');
    if (clearButton) clearButton.style.display = 'inline-flex';
    const enabled = document.getElementById('enableBlank');
    if (enabled) enabled.checked = true;
    const controls = document.getElementById('blankControls');
    if (controls) controls.style.display = 'block';
    _renderBlankManager();
    _updateOptionalQcState();
    updateUI();
  }

  function _updateOptionalQcState() {
    const state = document.querySelector('.optional-qc-state');
    const panel = document.querySelector('.optional-qc-card');
    if (state) {
      state.textContent = _blankUploads.length
        ? `${_blankUploads.length} 个空白已用于本批`
        : '未启用';
    }
    if (panel) {
      panel.classList?.toggle?.('has-selection', _blankUploads.length > 0);
    }
  }

  async function _appendBlankFiles(fileList, status) {
    const selected = Array.from(fileList || []);
    if (!selected.length) return;
    if (_blankUploads.length + selected.length > 10) {
      _resetBlankInputs();
      throw new Error('每批最多使用 10 个空白文件');
    }
    if (status) status.textContent = `⏳ 上传并解析 ${selected.length} 个空白...`;
    const uploaded = await _uploadBlankFiles(selected);
    _mergeAvailableBlankUploads(uploaded);
    const mode = _currentIonMode();
    const compatible = uploaded.filter(row => _blankSupportsMode(row, mode));
    const incompatible = uploaded.filter(row => !_blankSupportsMode(row, mode));
    const selectedIdentities = new Set(
      _blankUploads.map(row => _blankContentIdentity(row))
    );
    for (const row of compatible) {
      const identity = _blankContentIdentity(row);
      if (!selectedIdentities.has(identity)) {
        _blankUploads.push(row);
        selectedIdentities.add(identity);
      }
    }
    _syncBlankPaths();
    _resetBlankInputs();
    _renderBlankStatuses();
    if (incompatible.length) {
      Toast.show(
        `已保存 ${incompatible.length} 个空白，但其极性与当前 ${mode} 不兼容，未加入本批；可在空白管理中切换后使用`,
        'warning'
      );
    }
  }

  function detachBlank(index) {
    const row = _blankUploads[index];
    if (!row) return;
    _blankUploads.splice(index, 1);
    _syncBlankPaths();
    _renderBlankStatuses();
    Toast.show(`本批已不再使用空白：${row.blank_filename}`, 'info');
  }

  function selectManagedBlank(index) {
    const row = _availableBlankUploads[index];
    if (!row) return;
    const mode = _currentIonMode();
    if (!_blankSupportsMode(row, mode)) {
      const polarities = _blankPolarities(row);
      Toast.show(
        `该空白为 ${polarities.join('/') || '极性未知'}，不能用于当前 ${mode} 分析`,
        'error'
      );
      return;
    }
    if (_blankUploads.length >= 10) {
      Toast.show('每批最多使用 10 个空白文件', 'warning');
      return;
    }
    const duplicateContent = _blankUploads.some(
      item => _blankContentIdentity(item) === _blankContentIdentity(row)
    );
    if (!duplicateContent) {
      _blankUploads.push(row);
      _syncBlankPaths();
      _renderBlankStatuses();
    } else if (!_blankUploads.some(item => item.blank_path === row.blank_path)) {
      Toast.show(localizedCopy(
        `该空白与本批已选文件内容相同，未重复计数：${row.blank_filename}`,
        `This blank is byte-identical to one already selected and was not counted twice: ${row.blank_filename}`
      ), 'info');
    }
  }

  function detachManagedBlank(index) {
    const row = _availableBlankUploads[index];
    if (!row) return;
    const selectedIndex = _blankUploads.findIndex(
      item => item.blank_path === row.blank_path
    );
    if (selectedIndex >= 0) detachBlank(selectedIndex);
  }

  async function _deleteBlankRow(row) {
    await API.del(`/blanks/${encodeURIComponent(row.blank_path)}`);
    _blankManagerDeleteSelection.delete(String(row.blank_path));
    _blankUploads = _blankUploads.filter(
      item => item.blank_path !== row.blank_path
    );
    _availableBlankUploads = _availableBlankUploads.filter(
      item => item.blank_path !== row.blank_path
    );
    _syncBlankPaths();
    _resetBlankInputs();
    _renderBlankStatuses();
  }

  async function removeBlank(index) {
    const row = _blankUploads[index];
    if (!row) return;
    try {
      await _deleteBlankRow(row);
      Toast.show(`已删除空白文件：${row.blank_filename}`, 'success');
    } catch (error) {
      Toast.show(`空白删除失败：${error.message}`, 'error');
    }
  }

  async function deleteManagedBlank(index) {
    const row = _availableBlankUploads[index];
    if (!row) return;
    if (row.in_use) {
      Toast.show('该空白正被排队或运行任务使用，暂不能删除', 'warning');
      return;
    }
    try {
      await _deleteBlankRow(row);
      Toast.show(`已删除空白文件：${row.blank_filename}`, 'success');
    } catch (error) {
      Toast.show(`空白删除失败：${error.message}`, 'error');
      await refreshBlankUploads();
    }
  }

  async function deleteSelectedBlanks() {
    const rows = _availableBlankUploads.filter(
      row => (
        !row.in_use
        && _blankManagerDeleteSelection.has(String(row.blank_path))
      )
    );
    if (!rows.length) {
      Toast.show('请先勾选要删除的空白文件', 'warning');
      return;
    }
    if (!confirm(`确认永久删除选中的 ${rows.length} 个空白文件？此操作不可恢复。`)) {
      return;
    }

    try {
      const result = await API.post('/blanks/batch/delete', {
        blank_tokens: rows.map(row => String(row.blank_path)),
      });
      const deletedTokens = new Set(
        (result.deleted || []).map(item =>
          String(item.blank_path || item)
        )
      );
      deletedTokens.forEach(token => {
        _blankManagerDeleteSelection.delete(token);
      });
      _blankUploads = _blankUploads.filter(
        row => !deletedTokens.has(String(row.blank_path))
      );
      _availableBlankUploads = _availableBlankUploads.filter(
        row => !deletedTokens.has(String(row.blank_path))
      );
      _syncBlankPaths();
      _resetBlankInputs();
      _renderBlankStatuses();
      await refreshBlankUploads();

      const deletedCount = Number(
        result.deleted_count ?? deletedTokens.size
      );
      const inUseCount = Number(
        result.skipped_in_use_count
        ?? (result.skipped_in_use || []).length
      );
      const failedCount = Number(
        result.failed_count ?? (result.failed || []).length
      );
      const message = `已批量删除 ${deletedCount} 个空白`
        + `${inUseCount ? `；${inUseCount} 个正被任务使用，已保留` : ''}`
        + `${failedCount ? `；${failedCount} 个删除失败` : ''}`
        + `${result.metadata_warning ? `；${result.metadata_warning}` : ''}`;
      Toast.show(
        message,
        (inUseCount || failedCount || result.metadata_warning)
          ? 'warning' : 'success'
      );
    } catch (error) {
      Toast.show(`批量删除空白失败：${error.message}`, 'error');
      await refreshBlankUploads();
    }
  }

  async function clearBlankUploads(deleteRemote = true) {
    const rows = _blankUploads.slice();
    if (!rows.length || !deleteRemote) {
      _resetBlankSelection();
      return;
    }
    try {
      const result = await API.post('/blanks/batch/delete', {
        blank_tokens: rows.map(row => String(row.blank_path)),
      });
      const deletedTokens = new Set(
        (result.deleted || []).map(
          item => String(item.blank_path || item)
        )
      );
      const failed = rows.filter(
        row => !deletedTokens.has(String(row.blank_path))
      );
      deletedTokens.forEach(token => {
        _blankManagerDeleteSelection.delete(token);
      });
      _availableBlankUploads = _availableBlankUploads.filter(
        row => !deletedTokens.has(String(row.blank_path))
      );
      _blankUploads = failed;
      _syncBlankPaths();
      _resetBlankInputs();
      if (failed.length) {
        _renderBlankStatuses();
        Toast.show(`${failed.length} 个空白仍被任务使用或删除失败`, 'warning');
      } else {
        _resetBlankSelection();
        Toast.show('本批所选空白文件已删除', 'success');
      }
    } catch (error) {
      Toast.show(`本批空白删除失败：${error.message}`, 'error');
      await refreshBlankUploads();
    }
  }

  async function handleIonModeChange(mode) {
    const normalizedMode = String(mode || 'AUTO').toUpperCase();
    if (normalizedMode === 'AUTO') {
      _renderBlankManager();
      return;
    }
    const incompatible = _blankUploads.filter(
      row => !_blankSupportsMode(row, normalizedMode)
    );
    if (!incompatible.length) {
      _renderBlankManager();
      return;
    }
    const incompatibleTokens = new Set(
      incompatible.map(row => String(row.blank_path))
    );
    _blankUploads = _blankUploads.filter(
      row => !incompatibleTokens.has(String(row.blank_path))
    );
    _syncBlankPaths();
    _renderBlankStatuses();
    Toast.show(
      `已切换为 ${normalizedMode}；${incompatible.length} 个不兼容空白已从本批解除，文件仍保留在空白管理中`,
      'warning'
    );
  }

  async function uploadBlank(fileList) {
    const status = document.getElementById('blankStatus');
    try {
      await _appendBlankFiles(fileList, status);
    } catch (e) {
      _resetBlankInputs();
      if (status) status.innerHTML = `✗ 空白上传失败：${escapeHtml(e.message)}`;
    }
  }

  // v3: 主上传区快捷空白入口 — 上传后自动启用扣除(勾选 enableBlank), 状态镜像到主区 + 分析参数区
  async function uploadBlankQuick(fileList) {
    const status = document.getElementById('blankStatusQuick');
    try {
      await _appendBlankFiles(fileList, status);
    } catch (e) {
      _resetBlankInputs();
      if (status) status.innerHTML = `✗ 空白上传失败：${escapeHtml(e.message)}`;
    }
  }

  return {
    init, handleFiles, removeFile, clearFiles,
    handleVendorFolder, clearVendorFolder,
    selectLocalVendorFolder, selectLocalThermoFolder, updateVendorSourceRole,
    prepareVendorConverter,
    loadLibrary, checkLibraryStatus, startAnalysis,
    analysisParams, analysisParamsForRun,
    refreshDbCounts, refreshRuntimeIdentity,
    toggleBlank, updateBlankPolicyUI, uploadBlank, uploadBlankQuick,
    refreshBlankUploads, handleIonModeChange,
    detachBlank, removeBlank, clearBlankUploads,
    selectManagedBlank, detachManagedBlank, deleteManagedBlank, updateBlankRole,
    toggleBlankDeleteSelection, toggleSelectAllBlankDeletes,
    deleteSelectedBlanks,
    loadExpectedCompoundDatasets, toggleExpectedCompounds,
    clearExpectedDatasetSelection, updateExpectedDatasetStatus,
    get stagedFiles() { return _stagedFiles; },
  };
})();
