from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_task_page_exposes_explicit_cross_occurrence_action():
    html = (ROOT / "web/static/index.html").read_text(encoding="utf-8")
    javascript = (ROOT / "web/static/js/tasks.js").read_text(encoding="utf-8")

    assert 'id="btnCrossOccurrence"' in html
    assert "Tasks.crossOccurrenceRank()" in html
    assert "same_lc_method_confirmed: true" in javascript
    assert "activate: true" in javascript
    assert "至少勾选 3 个" in javascript
    assert "不会把结果标为标准品确证" in javascript


def test_server_registers_rank_activate_read_and_rollback_routes():
    from web.server import app

    routes = {
        (route.path, method)
        for route in app.routes
        for method in (getattr(route, "methods", None) or [])
    }
    assert ("/api/workflows/cross-occurrence-ranking", "POST") in routes
    assert ("/api/workflows/cross-occurrence-ranking/{cohort_id}", "GET") in routes
    assert ("/api/workflows/cross-occurrence-ranking/rollback", "POST") in routes
