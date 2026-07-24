"""FastAPI app: the JSON API tier (server-rendered pages arrive in Task 11).
Small requests run synchronously under a hard wall-time net; larger ones become
jobs. Errors surface a client-safe error type (a QuiverlabError subclass name or
a runner tag); unexpected internal exceptions are genericised and their full
tracebacks stay in the server log."""
from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from webapp.server.catalog import build_catalog
from webapp.server.config import Config, get_config
from webapp.server.estimator import classify
from webapp.server.instant import run_with_timeout
from webapp.server.limits import check_can_queue
from webapp.server.runner import RunError, build_algebra
from webapp.server.schema import ComputeRequest
from webapp.server.security import (
    GENERIC_ERROR_MESSAGE, GENERIC_ERROR_TYPE, SECURITY_HEADERS, hash_ip,
    is_safe_error_type, sanitize_error, valid_ulid,
)
from webapp.server.store import JobStore

_log = logging.getLogger("quiverlab_web.app")
_STATIC = Path(__file__).resolve().parent.parent / "static"


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _stamp_security_headers(response):
    """Stamp the strict security headers on a response (idempotent via
    ``setdefault``). One source of truth (``SECURITY_HEADERS``) shared by the
    header middleware AND the catch-all exception handler -- the latter runs
    above the middleware, so it must stamp here too or an unhandled 500 ships
    bare."""
    for k, v in SECURITY_HEADERS.items():
        response.headers.setdefault(k, v)
    return response


def client_ip(request: Request) -> str:
    """The originating client address: the first X-Forwarded-For hop (behind the
    reverse proxy) or the direct peer. Hashed with a salt before it ever reaches
    the store or a log -- the raw value stays inside this function's callers only
    long enough to hash it."""
    fwd = request.headers.get("x-forwarded-for")
    if fwd:
        return fwd.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


def _build_or_error(spec):
    """Build the algebra, normalising every failure to a ``RunError``.

    ``build_algebra`` raises a mix of shapes: ``RunError`` (its FieldError /
    missing-builder tags), a raw webapp ``CatalogError`` (from
    ``validate_family``), a library ``QuiverlabError`` subclass, or -- for a bad
    param type -- an arbitrary exception from the builder. We tag each with its
    class name; ``sanitize_error`` then decides pass-through vs genericise."""
    try:
        return build_algebra(spec)
    except RunError:
        raise
    except Exception as exc:            # CatalogError / QuiverlabError / unexpected
        raise RunError(type(exc).__name__, str(exc))


def _error_response(error_type: str, message: str) -> JSONResponse:
    """A JSON error response with the error type sanitised. Safe (client-facing)
    errors are a 422; a genericised internal error is a 500, and its original
    type is logged server-side (never returned)."""
    safe = is_safe_error_type(error_type)
    etype, msg = sanitize_error(error_type, message)
    if not safe:
        _log.warning("genericised internal error (original type=%s)", error_type)
    return JSONResponse(status_code=(422 if safe else 500),
                        content={"error_type": etype, "message": msg})


def _reject_message(reason: str | None) -> str:
    if reason == "big_disabled":
        return ("computation exceeds the anonymous limit and the big-job tier is "
                "unavailable; run locally: pip install quiverlab")
    return "computation too large; narrow the degree range or run locally: pip install quiverlab"


def create_app(cfg: Config | None = None, mailer=None) -> FastAPI:
    cfg = cfg or get_config()
    store = JobStore(cfg.db_path)
    store.init_schema()
    app = FastAPI(title="quiverlab-web")
    if _STATIC.exists():
        app.mount("/static", StaticFiles(directory=_STATIC), name="static")

    @app.middleware("http")
    async def _security_headers(request: Request, call_next):
        return _stamp_security_headers(await call_next(request))

    @app.exception_handler(Exception)
    async def _unhandled_exception(request: Request, exc: Exception) -> JSONResponse:
        # A GENUINELY unhandled exception escapes the route and unwinds past the
        # header middleware to Starlette's ServerErrorMiddleware -- which would
        # otherwise ship a bare PlainTextResponse with NO security headers. This
        # catch-all runs there: it logs the type only (never the request body or
        # client IP), returns the same genericised 500 the honest path uses (no
        # class name or message leaks), and stamps the security headers itself.
        _log.error("unhandled server error (type=%s)", type(exc).__name__)
        return _stamp_security_headers(JSONResponse(
            status_code=500,
            content={"error_type": GENERIC_ERROR_TYPE, "message": GENERIC_ERROR_MESSAGE}))

    def _ip_hash(request: Request) -> str:
        # Hash immediately: the raw address must never reach the store or a log.
        return hash_ip(client_ip(request), cfg.ip_hash_salt)

    @app.get("/api/catalog")
    def api_catalog():
        return build_catalog()

    @app.post("/api/compute")
    def api_compute(req: ComputeRequest, request: Request):
        iph = _ip_hash(request)
        try:
            # Build the algebra to read its dimension. This runs BEFORE (outside)
            # the wall net, so its construction cost is bounded by catalog
            # validation (the degree/param caps in validate_family), not the
            # instant timeout.
            A = _build_or_error(req.algebra)   # validation + honest errors
            dim = A.dim                        # library attribute (never .dimension())
        except RunError as exc:
            return _error_response(exc.error_type, exc.message)

        info = classify(dim, req, cfg)
        tier = info["tier"]

        if tier == "instant":
            try:
                result = run_with_timeout(req, cfg)   # built algebra NOT retained
            except RunError as exc:
                return _error_response(exc.error_type, exc.message)
            if result is not None:
                return {"tier": "instant", "result": result}
            tier = "queued"                    # wall net overflowed -> queue anonymously

        if tier == "big":
            # NEVER a hard reject at the anonymous boundary: offer the big-job
            # tier with the honest numbers; the page reveals the email field.
            return JSONResponse(status_code=202,
                                content={"tier": "big", "estimate": info["estimate"]})

        if tier == "reject":
            return JSONResponse(status_code=422,
                                content={"error_type": "TooLarge", "reason": info["reason"],
                                         "estimate": info["estimate"],
                                         "message": _reject_message(info["reason"])})

        # tier == "queued"
        # Soft limit: this read-then-write (count check, then create_job) is NOT
        # atomic across concurrent requests, so the cap can be marginally exceeded
        # under a race -- acceptable for anonymous abuse control.
        refusal = check_can_queue(store, cfg, iph, _now_iso())
        if refusal:
            return JSONResponse(status_code=429,
                                content={"error_type": "RateLimited", "message": refusal})
        jid = store.create_job(req.model_dump(by_alias=True), ip=iph)
        return JSONResponse(status_code=202, content={"tier": "queued", "job_id": jid})

    @app.post("/api/jobs")
    def api_create_job(req: ComputeRequest, request: Request):
        iph = _ip_hash(request)
        # Soft limit: read-then-write is non-atomic across concurrent requests
        # (see api_compute); acceptable for anonymous abuse control.
        refusal = check_can_queue(store, cfg, iph, _now_iso())
        if refusal:
            return JSONResponse(status_code=429,
                                content={"error_type": "RateLimited", "message": refusal})
        jid = store.create_job(req.model_dump(by_alias=True), ip=iph)
        return JSONResponse(status_code=202, content={"job_id": jid})

    @app.get("/api/jobs/{job_id}")
    def api_get_job(job_id: str):
        # Validate the id shape BEFORE touching the store -- a malformed id can
        # never name a real job.
        if not valid_ulid(job_id):
            return JSONResponse(status_code=404, content={"message": "no such job"})
        job = store.get_job(job_id)
        if job is None:
            return JSONResponse(status_code=404, content={"message": "no such job"})
        return {"id": job.id, "status": job.status, "progress": job.progress,
                "error": job.error, "created_at": job.created_at,
                "finished_at": job.finished_at}

    _register_pages(app, cfg, store)   # Task 11
    # Task 12: feedback form page, JSON submit API, token-gated admin table.
    # Imported here (not at module top) so feedback.py can reuse this module's
    # client_ip/_now_iso helpers without an import cycle.
    from webapp.server.feedback import register_feedback
    register_feedback(app, cfg, store)
    # Task 13: big-job email magic-link tier (spec §17). Imported here (not at
    # module top) so bigjobs.py can reuse this module's _build_or_error /
    # _error_response helpers without an import cycle. The production worker
    # passes the real SMTP mailer; None here leaves the default transport, built
    # lazily only when a link is actually sent.
    from webapp.server.bigjobs import register_big_jobs
    register_big_jobs(app, cfg, store, mailer=mailer)
    return app


def _register_pages(app, cfg, store) -> None:
    # Task 11: server-rendered bilingual pages + artifact downloads.
    from webapp.server.pages import register
    register(app, cfg, store)
