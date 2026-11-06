import json
from pathlib import Path
import subprocess

from fastapi.testclient import TestClient


ROOT = Path(__file__).resolve().parents[1]


def test_blank_batch_delete_removes_free_files_and_preserves_in_use(
    monkeypatch, tmp_path,
):
    import web.server as server

    blank_dir = tmp_path / "blanks"
    blank_dir.mkdir()
    blank_index = tmp_path / "blank_uploads_index.json"
    free_a = "11111111_free-a.mzML"
    in_use = "22222222_in-use.mzML"
    free_b = "33333333_free-b.mgf"
    missing = "44444444_missing.mzML"
    for token in (free_a, in_use, free_b):
        (blank_dir / token).write_bytes(b"blank")

    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(server, "BLANK_UPLOADS_FILE", blank_index)
    monkeypatch.setattr(server, "_tasks", {
        "running-task": {
            "id": "running-task",
            "status": "running",
            "params": {"blank_files": [in_use]},
        }
    })
    with server._blank_uploads_file_lock:
        server._write_blank_upload_index_unlocked({
            token: {
                "blank_path": token,
                "blank_filename": server._blank_original_filename(token),
                "blank_role": "unspecified",
            }
            for token in (free_a, in_use, free_b)
        })

    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/blanks/batch/delete",
            json={
                "blank_tokens": [
                    free_a, in_use, free_b, free_a, missing,
                ]
            },
        )
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["requested_count"] == 4
        assert body["deleted_count"] == 2
        assert {
            item["blank_path"] for item in body["deleted"]
        } == {free_a, free_b}
        assert body["skipped_in_use_count"] == 1
        assert body["skipped_in_use"] == [{
            "blank_path": in_use,
            "blank_filename": "in-use.mzML",
            "active_task_ids": ["running-task"],
        }]
        assert body["failed_count"] == 1
        assert body["failed"][0]["blank_path"] == missing

        assert not (blank_dir / free_a).exists()
        assert not (blank_dir / free_b).exists()
        assert (blank_dir / in_use).exists()
        saved = json.loads(blank_index.read_text(encoding="utf-8"))
        assert set(saved["uploads"]) == {in_use}
    finally:
        client.close()


def test_blank_manager_can_select_all_deletable_rows_and_delete_once():
    script = r"""
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');

const elements = new Map();
function element(id) {
  if (!elements.has(id)) {
    elements.set(id, {
      id,
      value: id === 'ionMode' ? 'auto' : '',
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

let blanks = [
  { blank_path: 'pos-a.mzML', blank_filename: 'pos-a.mzML', polarities: ['POS'], size_mb: 1, in_use: false },
  { blank_path: 'neg-a.mzML', blank_filename: 'neg-a.mzML', polarities: ['NEG'], size_mb: 1, in_use: false },
  { blank_path: 'active.mzML', blank_filename: 'active.mzML', polarities: ['POS'], size_mb: 1, in_use: true, active_task_ids: ['task-1'] },
];
const posted = [];
const toasts = [];
const context = {
  console,
  Map,
  Set,
  JSON,
  Promise,
  FormData: class {},
  confirm() { return true; },
  sessionStorage: { getItem() { return null; }, setItem() {} },
  document: {
    getElementById: element,
    querySelector() { return null; },
    querySelectorAll() { return []; },
    createElement() { return { textContent: '', innerHTML: '' }; },
  },
  API: {
    async get(path) {
      assert.strictEqual(path, '/blanks');
      return { blanks, count: blanks.length };
    },
    async post(path, payload) {
      posted.push({
        path,
        payload: JSON.parse(JSON.stringify(payload)),
      });
      const deleted = payload.blank_tokens.map(token => ({
        blank_path: token,
        blank_filename: token,
      }));
      const deletedTokens = new Set(payload.blank_tokens);
      blanks = blanks.filter(row => !deletedTokens.has(row.blank_path));
      return {
        requested_count: payload.blank_tokens.length,
        deleted,
        deleted_count: deleted.length,
        skipped_in_use: [],
        skipped_in_use_count: 0,
        failed: [],
        failed_count: 0,
      };
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
  Upload.toggleSelectAllBlankDeletes(true);
  assert.ok(
    element('blankUploadManager').innerHTML.includes('删除所选 (2)')
  );
  await Upload.deleteSelectedBlanks();
  assert.deepStrictEqual(posted, [{
    path: '/blanks/batch/delete',
    payload: { blank_tokens: ['pos-a.mzML', 'neg-a.mzML'] },
  }]);
  assert.ok(
    element('blankUploadManager').innerHTML.includes('active.mzML')
    || element('blankUploadManager').innerHTML.includes('正被排队/运行任务使用')
  );
  assert.ok(toasts.some(item => item.message.includes('已批量删除 2 个空白')));
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
"""
    completed = subprocess.run(
        ["node", "-e", script],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
