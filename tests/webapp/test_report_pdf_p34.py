"""Plan 34 (post PDF/TeX removal) -- the webapp worked-steps report is the print-ready
HTML report + its JSON machine record. PDF/TeX report artifacts have been removed: the
worker never writes them, the job page never links them, and the download endpoint
refuses trace.pdf / trace.tex cleanly (404, whitelist).

Offline/fixture-based, in the idiom of tests/webapp/test_pages.py: the worker's
artifact production is exercised through the shared spec core (quiverlab.hpc.spec.run,
the same path the queued worker and instant tier use), and pages via TestClient.
"""
import json

from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.store import JobStore

# The runner fixture: a monomial GF(2) algebra with an HH computation + report artifact.
_SPEC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..2"],
    "artifacts": {"pdf": True, "tikz": False},
}


def _done_job(cfg, result_json='{"ok": true}'):
    store = JobStore(cfg.db_path)
    store.init_schema()
    jid = store.create_job({}, ip="h")
    store.claim_next()
    art = cfg.artifacts_dir / jid
    art.mkdir(parents=True, exist_ok=True)
    (art / "result.json").write_text(result_json)
    store.mark_done(jid, str(art))
    return jid, art


# --------------------------------------------------------------------------- #
# Worker artifact production (the spec core write_trace path): HTML report + JSON.
# --------------------------------------------------------------------------- #
def test_worker_produces_html_report_and_json(tmp_path):
    from quiverlab.hpc.spec import run as spec_run
    art = tmp_path / "art"
    art.mkdir()
    res = spec_run(_SPEC, art)
    assert (art / "trace_steps.html").exists(), "the worker must render the HTML report"
    assert (art / "trace.json").exists(), "the JSON machine record is persisted alongside"
    # PDF/TeX report artifacts are never produced:
    assert not (art / "trace.pdf").exists()
    assert not (art / "trace.tex").exists()
    assert res["meta"]["pdf"] == "worked steps in trace_steps.html"


def test_worker_html_report_is_print_ready(tmp_path):
    from quiverlab.hpc.spec import run as spec_run
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_SPEC, art)
    html = (art / "trace_steps.html").read_text()
    assert "@media print" in html and "<math" in html         # typeset + print CSS
    assert "Math is shown as TeX source" not in html          # no apology


# --------------------------------------------------------------------------- #
# Job page: honest HTML/JSON labels; no PDF/TeX links, no why-no-PDF note.
# --------------------------------------------------------------------------- #
def test_job_page_links_html_report_no_pdf_or_tex(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert f"/download/{jid}/trace_steps.html" in r.text
    assert "Report (HTML" in r.text                            # honest label
    # No PDF/TeX links, and no why-no-PDF apology (both removed):
    assert f"/download/{jid}/trace.pdf" not in r.text
    assert f"/download/{jid}/trace.tex" not in r.text
    assert "no LaTeX toolchain" not in r.text and "failed to compile" not in r.text
    assert "LaTeX source" not in r.text


def test_job_page_bilingual_html_label(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace_steps.html").write_text("<p>steps</p>")
    es = TestClient(create_app(cfg)).get(f"/es/job/{jid}")
    assert "Reporte (HTML" in es.text


def test_job_page_links_trace_json_when_present(tmp_path):
    # Plan 34: the JSON machine record is offered on the job page with an honest label.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace.json").write_text('{"quiverlab_trace_schema": 1, "events": []}')
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert f"/download/{jid}/trace.json" in r.text
    assert "Report data (JSON)" in r.text
    es = TestClient(create_app(cfg)).get(f"/es/job/{jid}")
    assert "Datos del reporte (JSON)" in es.text                 # bilingual label


def test_job_page_omits_trace_json_link_when_absent(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, _art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert f"/download/{jid}/trace.json" not in r.text           # not fabricated


def test_download_trace_json_served_as_application_json(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace.json").write_text('{"quiverlab_trace_schema": 1, "events": []}')
    r = TestClient(create_app(cfg)).get(f"/download/{jid}/trace.json")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("application/json")
    assert json.loads(r.content)["quiverlab_trace_schema"] == 1


def test_download_pdf_and_tex_are_refused(tmp_path):
    """The removed report artifacts are OFF the download whitelist: even if a stale
    file with that name sits in the artifact dir, the endpoint 404s (never a 500)."""
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace.pdf").write_bytes(b"%PDF-1.5\n%%EOF\n")
    (art / "trace.tex").write_text(r"\documentclass{article}")
    client = TestClient(create_app(cfg))
    for name in ("trace.pdf", "trace.tex"):
        r = client.get(f"/download/{jid}/{name}")
        assert r.status_code == 404, f"{name} must not be served"
