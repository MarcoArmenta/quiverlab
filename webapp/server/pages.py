"""Server-rendered bilingual pages (Task 11) + language-neutral artifact
downloads.

Every page is mounted twice -- under ``/…`` (English) and ``/es/…`` (Spanish) --
by the SAME handlers; the language is a closure constant, never a route
parameter, so FastAPI never exposes it. Every chrome string goes through
``i18n.t`` bound to the page's language; the header language toggle links to the
same page under the other prefix.

Downloads are strict: the job id is ULID-validated with Task 9's ``valid_ulid``,
the filename is a whitelist of the runner's REAL artifact names (note
``trace_steps.html`` -- the HTML worked-steps report -- NOT ``trace.html``),
the file is resolved and confirmed to live directly inside the job's artifact
dir (no traversal), and it is served ONLY for a ``done`` job with the correct
Content-Type."""
from __future__ import annotations

import json
from pathlib import Path

from fastapi import Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates

from webapp.server.catalog import build_catalog
from webapp.server.i18n import t as _t
from webapp.server.references import grouped_bibliography
from webapp.server.security import sanitize_error_string, valid_ulid

_TEMPLATES = Jinja2Templates(
    directory=str(Path(__file__).resolve().parent.parent / "templates"))

# The runner (webapp/server/runner.py) writes these exact names: result.json, the
# print-ready worked-steps HTML report trace_steps.html, the JSON machine record
# trace.json (Plan 34 -- the complete worked-steps event stream, always persisted),
# and tikz.tex (the quiver diagram source). PDF/TeX report artifacts have been
# removed. The whitelist is also the traversal defense -- a name outside it is a 404.
_CONTENT_TYPES = {
    "result.json": "application/json",
    "trace_steps.html": "text/html; charset=utf-8",
    "trace.json": "application/json; charset=utf-8",
    "tikz.tex": "text/plain; charset=utf-8",
}
_ALLOWED = frozenset(_CONTENT_TYPES)

# (url prefix, language) -- the two mount points every page gets.
_LANGS = (("", "en"), ("/es", "es"))


def _other_url(path: str) -> str:
    """The same page under the other language prefix (for the header toggle)."""
    if path == "/es" or path.startswith("/es/"):
        return path[3:] or "/"
    return "/es" + ("" if path == "/" else path)


def _ctx(request: Request, cfg, lang: str, prefix: str, **extra) -> dict:
    """The shared template context. ``t`` is pre-bound to this page's language so
    templates call ``t("key")`` with no lang argument.

    ``offline``/``caps`` come from ``app.state`` (stamped by
    ``webapp.server.offline.create_offline_app``); on the deployed server they are
    absent, so the offline footer block never renders there (Plan 28)."""
    state = request.app.state
    return {"lang": lang, "prefix": prefix,
            "t": (lambda k: _t(k, lang)),
            "other_url": _other_url(request.url.path),
            "docs_url": cfg.docs_url,                # "" ⇒ no Docs link (base.html)
            "big_max_cells": cfg.big_ops_threshold,  # big.reject {maxcells}
            "big_jobs_enabled": cfg.big_jobs_enabled,
            "offline": bool(getattr(state, "offline", False)),
            "caps": getattr(state, "caps", None),
            **extra}


def register(app, cfg, store) -> None:
    for prefix, lang in _LANGS:
        _mount_pages(app, cfg, store, prefix, lang)
    _mount_download(app, cfg, store)


def _mount_pages(app, cfg, store, prefix: str, lang: str) -> None:
    # Each call has its own prefix/lang locals; the closures below bind them
    # correctly (no loop-variable capture bug). lang/prefix are NOT route
    # parameters -- they are closure constants, so FastAPI does not expose them.
    @app.get(prefix or "/", response_class=HTMLResponse)
    def index(request: Request):
        return _TEMPLATES.TemplateResponse(
            request, "index.html",
            _ctx(request, cfg, lang, prefix, catalog=build_catalog()))

    @app.get(prefix + "/about", response_class=HTMLResponse)
    def about(request: Request):
        return _TEMPLATES.TemplateResponse(
            request, "about.html", _ctx(request, cfg, lang, prefix))

    # ----- the draw-a-quiver canvas (docs/gui ported; Plan 10 editor, server
    # compute). gui.js is vendored BYTE-IDENTICAL under static/gui/ (parity is
    # test-pinned) and loads its siblings by page-relative URLs ("gui/…"), so
    # every asset is mounted under BOTH language prefixes alongside the page.
    @app.get(prefix + "/draw", response_class=HTMLResponse)
    def draw(request: Request):
        return _TEMPLATES.TemplateResponse(
            request, "draw.html", _ctx(request, cfg, lang, prefix))

    @app.get(prefix + "/gui/manifest.json")
    def gui_manifest():
        # The Pyodide worker downloads a wheel named here; the server-backed
        # worker only reads the version. "server" is a truthful non-null wheel
        # marker (gui.js refuses a null wheel as "engine payload not built").
        import quiverlab
        return {"schema": 1, "wheel": "server",
                "quiverlab_version": getattr(quiverlab, "__version__", "unknown")}

    _GUI_STATIC = Path(__file__).resolve().parent.parent / "static" / "gui"
    _GUI_FILES = {"gui.js": "text/javascript; charset=utf-8",
                  "gui.css": "text/css; charset=utf-8",
                  "worker.js": "text/javascript; charset=utf-8",
                  "mathjax-katex-shim.js": "text/javascript; charset=utf-8",
                  "presets.json": "application/json"}

    @app.get(prefix + "/gui/{name}")
    def gui_asset(name: str):
        # Fixed whitelist (same posture as /download): unknown names never
        # touch the filesystem.
        if name not in _GUI_FILES:
            return JSONResponse(status_code=404, content={"message": "not found"})
        return FileResponse(_GUI_STATIC / name, media_type=_GUI_FILES[name])

    @app.get(prefix + "/literature", response_class=HTMLResponse)
    def literature(request: Request):
        return _TEMPLATES.TemplateResponse(
            request, "literature.html",
            _ctx(request, cfg, lang, prefix, groups=grouped_bibliography()))

    @app.get(prefix + "/job/{job_id}", response_class=HTMLResponse)
    def job_page(request: Request, job_id: str):
        # Validate the id shape BEFORE touching the store (Task 9 helper).
        if not valid_ulid(job_id):
            return HTMLResponse(_t("job.not_found", lang), status_code=404)
        job = store.get_job(job_id)
        if job is None:
            return HTMLResponse(_t("job.not_found", lang), status_code=404)
        art = cfg.artifacts_dir / job_id
        # The worked-steps artifact is the print-ready HTML report (Plan 34; PDF/TeX
        # output has been removed). The browser's Print -> Save as PDF exports a PDF.
        has_html_report = (art / "trace_steps.html").exists()
        # The JSON machine record of the worked steps (Plan 34 -- the complete event
        # stream) is offered alongside the HTML report.
        has_trace_json = (art / "trace.json").exists()
        has_tikz = (art / "tikz.tex").exists()
        reproduce = version = None
        references: list = []
        if job.status == "done" and (art / "result.json").exists():
            try:
                data = json.loads((art / "result.json").read_text(encoding="utf-8"))
                reproduce = data.get("reproduce")
                version = data.get("quiverlab_version")
                references = data.get("references", []) or []  # runner-resolved
            except (json.JSONDecodeError, OSError):
                pass
        # "Reproduce on a cluster" (Plan 28): the request spec as a runnable YAML
        # config, rendered next to the local-Python snippet. Only for a done job
        # whose stored spec is a real request (skips e.g. empty test fixtures); a
        # serialisation hiccup must never 500 the page.
        cluster_yaml = None
        if job.status == "done" and isinstance(job.spec, dict) and job.spec.get("algebra"):
            try:
                from webapp.server.clusterconfig import cluster_config_yaml
                cluster_yaml = cluster_config_yaml(job.spec)
            except Exception:
                cluster_yaml = None
        # Genericise the error at the READ boundary (same helper the JSON API
        # uses): the stored value stays raw for forensics, but no unexpected
        # internal type/message is rendered into the page.
        return _TEMPLATES.TemplateResponse(
            request, "job.html",
            _ctx(request, cfg, lang, prefix, job=job,
                 has_html_report=has_html_report,
                 has_trace_json=has_trace_json,
                 has_tikz=has_tikz, reproduce=reproduce,
                 version=version, references=references, cluster_yaml=cluster_yaml,
                 error=sanitize_error_string(job.error)))


def _mount_download(app, cfg, store) -> None:
    @app.get("/download/{job_id}/{name}")
    def download(job_id: str, name: str):
        # ULID + whitelist first: a malformed id or a non-artifact name never
        # reaches the store or the filesystem.
        if not valid_ulid(job_id) or name not in _ALLOWED:
            return JSONResponse(status_code=404, content={"message": "not found"})
        # Only completed jobs expose their artifacts.
        job = store.get_job(job_id)
        if job is None or job.status != "done":
            return JSONResponse(status_code=404, content={"message": "not found"})
        base = (cfg.artifacts_dir / job_id).resolve()
        path = (base / name).resolve()
        # Defense in depth: the resolved file must live directly inside the job's
        # artifact dir (the whitelist already forbids separators/traversal).
        if path.parent != base or not path.is_file():
            return JSONResponse(status_code=404, content={"message": "not found"})
        return FileResponse(path, media_type=_CONTENT_TYPES[name], filename=name)
