"""Browser-level contract for the presentation-only compound-table search."""

import ast
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SEARCH = ROOT / "web/static/js/compound-search.js"
PYTHON_PROJECTION = ROOT / "web/governance_blocked_projection.py"
INDEX = ROOT / "web/static/index.html"


def _python_aliases():
    tree = ast.parse(PYTHON_PROJECTION.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            getattr(target, "id", None) == "ZH_COMPOUND_ALIASES"
            for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise AssertionError("ZH_COMPOUND_ALIASES is missing from the Python projection")


def test_compound_search_reuses_the_governance_alias_dictionary_and_is_loaded_first():
    source = SEARCH.read_text(encoding="utf-8")
    literal = re.search(
        r"const ZH_COMPOUND_ALIASES = Object\.freeze\((\{.*?\})\);",
        source,
        flags=re.DOTALL,
    )
    assert literal, "compound-search.js must carry the presentation alias dictionary"
    assert json.loads(literal.group(1)) == _python_aliases()

    page = INDEX.read_text(encoding="utf-8")
    assert page.index("/static/js/compound-search.js") < page.index(
        "/static/js/sample-overview.js"
    )


def test_compound_search_projects_aliases_competing_names_and_mz_without_mutation():
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to execute the browser search helper")
    program = f"""
const fs = require('fs');
const vm = require('vm');
vm.runInThisContext(fs.readFileSync({json.dumps(str(SEARCH))}, 'utf8'));
const rows = [
  {{ id: 'sinomenine', names: ['sinomenine'], terms: ['C19H23NO4', 'SINOIK'], mzs: [330.169] }},
  {{ id: 'competition', names: ['coptisine', 'berberine'], terms: ['C20H18NO4'], mzs: [336.123] }},
];
const before = JSON.stringify(rows);
const select = query => rows.filter(row => CompoundSearch.matches(row, query)).map(row => row.id);
const result = {{
  empty: select(''),
  absent: select('no-such-compound'),
  english: select('sinomenine'),
  chinese: select('青藤碱'),
  competitor: select('小檗碱'),
  formula: select('C19H23NO4'),
  inchikey: select('SINOIK'),
  mzLowerBoundary: CompoundSearch.matches({{ mzs: [300.0] }}, '299.5'),
  mzUpperBoundary: CompoundSearch.matches({{ mzs: [300.0] }}, '300.5'),
  mzOutside: CompoundSearch.matches({{ mzs: [300.0] }}, '300.5001'),
  unchanged: before === JSON.stringify(rows),
}};
console.log(JSON.stringify(result));
"""
    run = subprocess.run([node, "-e", program], check=True, capture_output=True, text=True)
    assert json.loads(run.stdout) == {
        "empty": ["sinomenine", "competition"],
        "absent": [],
        "english": ["sinomenine"],
        "chinese": ["sinomenine"],
        "competitor": ["competition"],
        "formula": ["sinomenine"],
        "inchikey": ["sinomenine"],
        "mzLowerBoundary": True,
        "mzUpperBoundary": True,
        "mzOutside": False,
        "unchanged": True,
    }
