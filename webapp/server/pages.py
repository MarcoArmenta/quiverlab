"""Server-rendered bilingual pages (Task 11) + language-neutral artifact
downloads.

Every page is mounted twice -- under ``/…`` (English) and ``/es/…`` (Spanish) --
by the SAME handlers; the language is a closure constant, never a route
parameter, so FastAPI never exposes it. Every chrome string goes through
``i18n.t`` bound to the page's language; the header language toggle links to the
same page under the other prefix.

Downloads are strict: the job id is ULID-validated with Task 9's ``valid_ulid``,
the filename is a whitelist of the runner's REAL artifact names (note
``trace_steps.html`` -- the HTML worked-steps fallback -- NOT ``trace.html``),
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

# The runner (webapp/server/runner.py) writes these exact names: result.json,
# the worked-steps trace.pdf OR its HTML fallback trace_steps.html, and tikz.tex.
# The whitelist is also the traversal defense -- a name outside it is a 404.
_CONTENT_TYPES = {
    "result.json": "application/json",
    "trace.pdf": "application/pdf",
    "trace_steps.html": "text/html; charset=utf-8",
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
    templates call ``t("key")`` with no lang argument."""
    return {"lang": lang, "prefix": prefix,
            "t": (lambda k: _t(k, lang)),
            "other_url": _other_url(request.url.path),
            "docs_url": cfg.docs_url,                # "" ⇒ no Docs link (base.html)
            "big_max_cells": cfg.big_ops_threshold,  # big.reject {maxcells}
            "big_jobs_enabled": cfg.big_jobs_enabled,
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
        # The worked-steps artifact is either the PDF or its HTML fallback.
        trace_name = ("trace.pdf" if (art / "trace.pdf").exists()
                      else "trace_steps.html" if (art / "trace_steps.html").exists()
                      else None)
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
        # Genericise the error at the READ boundary (same helper the JSON API
        # uses): the stored value stays raw for forensics, but no unexpected
        # internal type/message is rendered into the page.
        return _TEMPLATES.TemplateResponse(
            request, "job.html",
            _ctx(request, cfg, lang, prefix, job=job, trace_name=trace_name,
                 has_tikz=has_tikz, reproduce=reproduce, version=version,
                 references=references, error=sanitize_error_string(job.error)))


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
