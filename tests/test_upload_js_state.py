"""Behavioral regression tests for the browser upload-state controller."""

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def _run_node(script: str) -> None:
    completed = subprocess.run(
        ["node", "-e", script],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_blank_selection_survives_refresh_and_detaches_on_mode_switch():
    _run_node(
        r"""
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');

const elements = new Map();
function element(id) {
  if (!elements.has(id)) {
    elements.set(id, {
      id,
      value: id === 'ionMode' ? 'NEG' : '',
      checked: false,
      disabled: false,
      innerHTML: '',
      textContent: '',
      style: {},
      classList: { toggle() {}, add() {}, remove() {} },
      addEventListener() {},
      querySelector() { return null; },
    });
  }
  return elements.get(id);
}
const storage = new Map();
const deleted = [];
const toasts = [];
const blanks = [
  { blank_path: 'pos.mzML', blank_filename: 'pos.mzML', polarities: ['POS'], size_mb: 1, in_use: false },
  { blank_path: 'neg.mzML', blank_filename: 'neg.mzML', polarities: ['NEG'], size_mb: 1, in_use: false },
];
const context = {
  console,
  Map,
  Set,
  JSON,
  Promise,
  FormData: class {},
  sessionStorage: {
    getItem(key) { return storage.has(key) ? storage.get(key) : null; },
    setItem(key, value) { storage.set(key, String(value)); },
  },
  document: {
    getElementById: element,
    querySelector() { return {}; },
    createElement() { return { textContent: '', innerHTML: '' }; },
  },
  API: {
    async get(path) {
      assert.strictEqual(path, '/blanks');
      return { blanks, count: blanks.length };
    },
    async del(path) {
      deleted.push(path);
      return { deleted: path };
    },
  },
  Toast: { show(message, level) { toasts.push({ message, level }); } },
};
vm.createContext(context);
const source = fs.readFileSync('web/static/js/upload.js', 'utf8');
vm.runInContext(source + '\nthis.UploadUnderTest = Upload;', context);
const Upload = context.UploadUnderTest;

(async () => {
  await Upload.refreshBlankUploads();
  Upload.selectManagedBlank(1);
  assert.deepStrictEqual(
    JSON.parse(storage.get('mssystem.v7.blank-selection')),
    ['neg.mzML']
  );

  element('ionMode').value = 'POS';
  await Upload.handleIonModeChange('POS');
  assert.deepStrictEqual(
    JSON.parse(storage.get('mssystem.v7.blank-selection')),
    []
  );
  assert.ok(toasts.some(item => item.message.includes('不兼容空白已从本批解除')));

  Upload.selectManagedBlank(1);
  assert.deepStrictEqual(
    JSON.parse(storage.get('mssystem.v7.blank-selection')),
    []
  );

  element('ionMode').value = 'NEG';
  Upload.selectManagedBlank(1);
  await Upload.deleteManagedBlank(1);
  assert.deepStrictEqual(
    JSON.parse(storage.get('mssystem.v7.blank-selection')),
    []
  );
  assert.deepStrictEqual(deleted, ['/blanks/neg.mzML']);
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
"""
    )


def test_run_parameters_include_selected_blanks_for_upload_and_queue_pages():
    _run_node(
        r"""
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');

const elements = new Map();
function element(id) {
  if (!elements.has(id)) {
    const values = {
      ionMode: 'auto',
      batchIonMode: 'auto',
      dbMode: 'full',
      blankFilterPolicy: 'hard_legacy',
      blankSampleRatio: '3',
      blankRtTolerance: '0.30',
      blankMinDetectionFraction: '0.50',
    };
    elements.set(id, {
      id,
      value: values[id] || '',
      checked: false,
      disabled: false,
      innerHTML: '',
      textContent: '',
      style: {},
      classList: { toggle() {}, add() {}, remove() {} },
      addEventListener() {},
      querySelector() { return null; },
      querySelectorAll() { return []; },
    });
  }
  return elements.get(id);
}
const storage = new Map();
const blanks = [
  {
    blank_path: 'pos-token.mzML',
    blank_filename: 'AutoMS-Blank-pos.mzML',
    polarities: ['POS'],
    blank_role: 'background_blank',
    size_mb: 1,
    in_use: false,
  },
  {
    blank_path: 'neg-token.mzML',
    blank_filename: 'AutoMS-Blank-neg.mzML',
    polarities: ['NEG'],
    blank_role: 'background_blank',
    size_mb: 1,
    in_use: false,
  },
];
const context = {
  console,
  Map,
  Set,
  JSON,
  Promise,
  FormData: class {},
  sessionStorage: {
    getItem(key) { return storage.has(key) ? storage.get(key) : null; },
    setItem(key, value) { storage.set(key, String(value)); },
  },
  document: {
    getElementById: element,
    querySelector(selector) {
      return selector === '.lib-status.loaded' ? {} : null;
    },
    querySelectorAll() { return []; },
    createElement() { return { textContent: '', innerHTML: '' }; },
  },
  API: {
    async get(path) {
      if (path === '/health') {
        return {
          service_generation: 'V16',
          default_port: 8516,
          capabilities: { blank_evidence_v2: true },
        };
      }
      if (path === '/blanks') return { blanks, count: blanks.length };
      throw new Error(`unexpected GET ${path}`);
    },
  },
  Toast: { show() {} },
};
context.window = {
  location: { port: '8516', protocol: 'http:' },
};
vm.createContext(context);
const source = fs.readFileSync('web/static/js/upload.js', 'utf8');
vm.runInContext(source + '\nthis.UploadUnderTest = Upload;', context);
const Upload = context.UploadUnderTest;

(async () => {
  await Upload.refreshRuntimeIdentity();
  await Upload.refreshBlankUploads();
  Upload.selectManagedBlank(0);
  Upload.selectManagedBlank(1);
  element('blankMethodComparable').checked = true;
  element('blankExclusionConfirmed').checked = true;

  const uploadParams = Upload.analysisParamsForRun('ionMode');
  const queueParams = Upload.analysisParamsForRun('batchIonMode');
  for (const params of [uploadParams, queueParams]) {
    assert.deepStrictEqual(
      JSON.parse(JSON.stringify(params.blank_files)),
      ['pos-token.mzML', 'neg-token.mzML']
    );
    assert.deepStrictEqual(
      JSON.parse(JSON.stringify(params.blank_roles)),
      {
        'pos-token.mzML': 'background_blank',
        'neg-token.mzML': 'background_blank',
      }
    );
    assert.strictEqual(params.blank_filter_policy, 'hard_legacy');
    assert.strictEqual(params.blank_method_comparability_confirmed, true);
    assert.strictEqual(params.blank_background_exclusion_confirmed, true);
  }
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
"""
    )


def test_byte_identical_blank_uploads_are_not_counted_as_replicates():
    _run_node(
        r"""
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');

const elements = new Map();
function element(id) {
  if (!elements.has(id)) {
    elements.set(id, {
      id,
      value: id === 'ionMode' ? 'NEG' : '',
      checked: false,
      disabled: false,
      innerHTML: '',
      textContent: '',
      style: {},
      classList: { toggle() {}, add() {}, remove() {} },
      addEventListener() {},
      querySelector() { return null; },
      querySelectorAll() { return []; },
    });
  }
  return elements.get(id);
}
const storage = new Map();
const toasts = [];
const blanks = [
  {
    blank_path: 'first-neg.mzML',
    blank_filename: 'AutoMS-Blank-1-2-neg.mzML',
    polarities: ['NEG'],
    blank_role: 'background_blank',
    content_sha256: 'abc123',
    size_mb: 1,
    in_use: false,
  },
  {
    blank_path: 'second-neg.mzML',
    blank_filename: 'AutoMS-Blank-1-2-neg-copy.mzML',
    polarities: ['NEG'],
    blank_role: 'background_blank',
    content_sha256: 'abc123',
    size_mb: 1,
    in_use: false,
  },
];
const context = {
  console,
  Map,
  Set,
  JSON,
  Promise,
  FormData: class {},
  sessionStorage: {
    getItem(key) { return storage.has(key) ? storage.get(key) : null; },
    setItem(key, value) { storage.set(key, String(value)); },
  },
  document: {
    getElementById: element,
    querySelector() { return {}; },
    querySelectorAll() { return []; },
    createElement() { return { textContent: '', innerHTML: '' }; },
  },
  API: {
    async get(path) {
      assert.strictEqual(path, '/blanks');
      return { blanks, count: blanks.length };
    },
  },
  Toast: { show(message, level) { toasts.push({ message, level }); } },
};
vm.createContext(context);
const source = fs.readFileSync('web/static/js/upload.js', 'utf8');
vm.runInContext(source + '\nthis.UploadUnderTest = Upload;', context);
const Upload = context.UploadUnderTest;

(async () => {
  await Upload.refreshBlankUploads();
  Upload.selectManagedBlank(0);
  Upload.selectManagedBlank(1);
  assert.deepStrictEqual(
    JSON.parse(storage.get('mssystem.v7.blank-selection')),
    ['first-neg.mzML']
  );
  assert.ok(toasts.some(item => item.message.includes('未重复计数')));
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
"""
    )


def test_failed_queue_reuses_uploaded_task_without_page_refresh():
    _run_node(
        r"""
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');

const elements = new Map();
function element(id) {
  if (!elements.has(id)) {
    const values = {
      ionMode: 'NEG',
      dbMode: 'full',
      peakProcessingPreset: 'library_match',
    };
    elements.set(id, {
      id,
      value: values[id] || '',
      checked: false,
      disabled: false,
      innerHTML: '',
      textContent: '',
      style: {},
      classList: { toggle() {}, add() {}, remove() {} },
      addEventListener() {},
      querySelector() { return null; },
    });
  }
  return elements.get(id);
}
const storage = new Map();
const toasts = [];
let uploadCalls = 0;
let batchCalls = 0;
let navigations = 0;
const context = {
  console,
  Map,
  Set,
  JSON,
  Promise,
  FormData: class {},
  sessionStorage: {
    getItem(key) { return storage.has(key) ? storage.get(key) : null; },
    setItem(key, value) { storage.set(key, String(value)); },
  },
  document: {
    getElementById: element,
    querySelector(selector) {
      return selector === '.lib-status.loaded' ? {} : null;
    },
    createElement() { return { textContent: '', innerHTML: '' }; },
  },
  API: {
    async upload() {
      uploadCalls += 1;
      return {
        uploaded: [{ task_id: 'pending-1', filename: 'sample.mgf' }],
        failed: [],
      };
    },
    async post(path) {
      assert.strictEqual(path, '/tasks/batch/analyze');
      batchCalls += 1;
      if (batchCalls === 1) {
        throw new Error('temporary queue transport error');
      }
      if (batchCalls === 2) {
        return {
          queued: [],
          queued_count: 0,
          skipped: [],
          failed: [{ id: 'pending-1', error: 'blank polarity mismatch' }],
        };
      }
      return {
        queued: ['pending-1'],
        queued_count: 1,
        skipped: [],
        failed: [],
      };
    },
    async get(path) {
      assert.strictEqual(path, '/blanks');
      return { blanks: [], count: 0 };
    },
  },
  Toast: { show(message, level) { toasts.push({ message, level }); } },
  App: { navigateTo(page) { assert.strictEqual(page, 'tasks'); navigations += 1; } },
};
vm.createContext(context);
const source = fs.readFileSync('web/static/js/upload.js', 'utf8');
vm.runInContext(source + '\nthis.UploadUnderTest = Upload;', context);
const Upload = context.UploadUnderTest;
const sample = {
  name: 'sample.mgf',
  size: 100,
  lastModified: 1,
  webkitRelativePath: '',
};
Upload.handleFiles([sample]);

(async () => {
  await Upload.startAnalysis();
  assert.strictEqual(uploadCalls, 1);
  assert.strictEqual(batchCalls, 1);
  assert.strictEqual(navigations, 0);
  assert.ok(element('btnStartAnalysis').textContent.includes('重试排队'));

  element('ionMode').value = 'POS';
  await Upload.startAnalysis();
  assert.strictEqual(uploadCalls, 1, 'the already-uploaded task must be reused');
  assert.strictEqual(batchCalls, 2);
  assert.strictEqual(navigations, 0);
  assert.ok(element('btnStartAnalysis').textContent.includes('重试排队'));

  await Upload.startAnalysis();
  assert.strictEqual(uploadCalls, 1, 'queue validation failure must not re-upload');
  assert.strictEqual(batchCalls, 3);
  assert.strictEqual(navigations, 1);
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
"""
    )
