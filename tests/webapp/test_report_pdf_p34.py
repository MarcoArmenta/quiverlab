"""Plan 34 -- the webapp worked-steps report: a real trace.pdf when a LaTeX engine
is available, an HONEST HTML report otherwise, and a job page that labels each
artifact truthfully (no "PDF" tag on an HTML file) and says WHY a PDF is missing.

Offline/fixture-based, in the idiom of tests/webapp/test_pages.py: the worker's
artifact production is exercised through the shared spec core (quiverlab.hpc.spec.run,
the same path the queued worker and instant tier use), and pages via TestClient.
"""
import json
import pathlib

import pytest
from fastapi.testclient import TestClient

from quiverlab.trace.writer import have_latex
from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.store import JobStore

# The runner fixture: a monomial GF(2) algebra with an HH computation + PDF artifact.
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
# Worker artifact production (the spec core write_trace path).
# --------------------------------------------------------------------------- #
def test_worker_produces_trace_pdf_when_engine_present(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    from quiverlab.hpc.spec import run as spec_run
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")

    def fake_compile(tex, out_pdf, engine):
        pathlib.Path(out_pdf).write_bytes(b"%PDF-1.5 fake\n%%EOF\n")
        return 1

    monkeypatch.setattr(W, "_compile_pdf", fake_compile)
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_SPEC, art)
    assert (art / "trace.pdf").exists(), "the worker must render trace.pdf with an engine"
    assert (art / "trace.tex").exists(), "the LaTeX source is persisted alongside"
    assert not (art / "trace_steps.html").exists()


def test_worker_html_fallback_is_print_ready_when_no_engine(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    from quiverlab.hpc.spec import run as spec_run
    monkeypatch.setattr(W, "have_latex", lambda: None)          # no toolchain
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_SPEC, art)
    assert not (art / "trace.pdf").exists()
    steps = art / "trace_steps.html"
    assert steps.exists(), "the worker must fall back to a print-ready HTML report"
    html = steps.read_text()
    assert "@media print" in html and "<math" in html         # typeset + print CSS
    assert "Math is shown as TeX source" not in html          # no apology


# --------------------------------------------------------------------------- #
# Job page: honest labels + a why-no-PDF note.
# --------------------------------------------------------------------------- #
def test_job_page_links_pdf_as_report_pdf(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace.pdf").write_bytes(b"%PDF-1.5\n%%EOF\n")
    (art / "trace.tex").write_text(r"\documentclass{article}")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert f"/download/{jid}/trace.pdf" in r.text
    assert "Report (PDF)" in r.text
    # a real PDF present => no missing-PDF apology
    assert "no LaTeX toolchain" not in r.text and "failed to compile" not in r.text


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


def test_job_page_labels_html_report_honestly_no_fake_pdf(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: None)         # deterministic reason
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert f"/download/{jid}/trace_steps.html" in r.text
    assert "Report (HTML" in r.text                            # honest label
    assert f"/download/{jid}/trace.pdf" not in r.text          # NO fake PDF link
    assert "no LaTeX toolchain" in r.text                      # says WHY


def test_job_page_missing_pdf_reason_compile(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: "pdflatex")   # engine present, no pdf
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert "failed to compile" in r.text                       # honest, engine-aware
    assert "no LaTeX toolchain" not in r.text


def test_job_page_bilingual_pdf_reason(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: None)
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')
    (art / "trace_steps.html").write_text("<p>steps</p>")
    es = TestClient(create_app(cfg)).get(f"/es/job/{jid}")
    assert "Reporte (HTML" in es.text
    assert "cadena de herramientas LaTeX" in es.text           # Spanish reason


# --------------------------------------------------------------------------- #
# End to end with a REAL engine: artifact dir has trace.tex AND trace.pdf, and the
# page links the PDF.
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# MAJOR-3: the WORKER records why no PDF was produced, and the job page reports
# the RECORDED reason -- not a re-probe on a possibly-different web tier (the
# split-deployment bug: worker has LaTeX, web tier does not, or vice versa).
# --------------------------------------------------------------------------- #
def test_worker_records_pdf_fallback_reason_no_toolchain(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    from quiverlab.hpc.spec import run as spec_run
    monkeypatch.setattr(W, "have_latex", lambda: None)
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_SPEC, art)
    meta = json.loads((art / "result.json").read_text())["meta"]
    assert meta["pdf_fallback_reason"] == "no_toolchain"


def test_worker_records_pdf_fallback_reason_compile(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    from quiverlab.hpc.spec import run as spec_run
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")

    def boom(tex, out_pdf, engine):
        raise RuntimeError("compile aborted")

    monkeypatch.setattr(W, "_compile_pdf", boom)
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_SPEC, art)
    meta = json.loads((art / "result.json").read_text())["meta"]
    assert meta["pdf_fallback_reason"] == "compile"


def test_worker_pdf_present_records_no_fallback_reason(tmp_path, monkeypatch):
    """A produced PDF leaves the meta byte-stable -- NO pdf_fallback_reason key."""
    from quiverlab.trace import writer as W
    from quiverlab.hpc.spec import run as spec_run
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")

    def fake_compile(tex, out_pdf, engine):
        pathlib.Path(out_pdf).write_bytes(b"%PDF-1.5 fake\n%%EOF\n")
        return 1

    monkeypatch.setattr(W, "_compile_pdf", fake_compile)
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_SPEC, art)
    meta = json.loads((art / "result.json").read_text())["meta"]
    assert "pdf_fallback_reason" not in meta


def test_job_page_prefers_recorded_compile_over_probe(tmp_path, monkeypatch):
    """Worker RECORDED 'compile'; the web tier's probe would say 'no_toolchain'. The
    page must report the recorded reason -- the whole point of MAJOR-3."""
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: None)          # probe would say no_toolchain
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, json.dumps(
        {"quiverlab_version": "9.9.9", "results": {},
         "meta": {"pdf_fallback_reason": "compile"}}))
    (art / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert "failed to compile" in r.text                        # RECORDED reason wins
    assert "no LaTeX toolchain" not in r.text                   # NOT the probe's answer


def test_job_page_prefers_recorded_no_toolchain_over_probe(tmp_path, monkeypatch):
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: "pdflatex")    # probe would say compile
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, json.dumps(
        {"quiverlab_version": "9.9.9", "results": {},
         "meta": {"pdf_fallback_reason": "no_toolchain"}}))
    (art / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert "no LaTeX toolchain" in r.text                       # RECORDED reason wins
    assert "failed to compile" not in r.text


def test_job_page_legacy_row_without_recorded_reason_falls_back_to_probe(tmp_path, monkeypatch):
    """A pre-Plan-34 result.json (no meta.pdf_fallback_reason) still gets an honest
    note via the live probe -- the documented legacy fallback."""
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: None)
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid, art = _done_job(cfg, '{"quiverlab_version": "9.9.9", "results": {}}')  # NO meta
    (art / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert "no LaTeX toolchain" in r.text                       # probe fallback


@pytest.mark.skipif(have_latex() is None, reason="no LaTeX toolchain on PATH")
def test_end_to_end_real_pdf_and_page_link(tmp_path):
    from quiverlab.hpc.spec import run as spec_run
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    store = JobStore(cfg.db_path)
    store.init_schema()
    jid = store.create_job({}, ip="h")
    store.claim_next()
    art = cfg.artifacts_dir / jid
    art.mkdir(parents=True, exist_ok=True)
    spec_run(_SPEC, art)
    store.mark_done(jid, str(art))
    assert (art / "trace.tex").exists()
    assert (art / "trace.pdf").exists(), "a real engine must yield trace.pdf"
    assert (art / "trace.pdf").read_bytes()[:5] == b"%PDF-"
    r = TestClient(create_app(cfg)).get(f"/job/{jid}")
    assert f"/download/{jid}/trace.pdf" in r.text and "Report (PDF)" in r.text
