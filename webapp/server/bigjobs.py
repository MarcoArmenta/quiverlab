"""Big-job magic-link flow (spec §17). A computation that overflows the anonymous
tier is offered as a *big job*: the requester leaves an email, receives a
single-use, spec-bound, expiring link, and clicking it queues the job. No
account, no session, no cookie.

The token is ``base64url(payload).base64url(HMAC-SHA256)``. It carries the
pending-big id, a hash of the exact spec, and an expiry -- all HMAC-signed with
``cfg.token_secret`` and verified in constant time. Single use is enforced by
atomically consuming the ``pending_big`` row (a second click finds nothing). The
plaintext email lives only in that row (and the queued job) and is deleted right
after the completion notice; a salted hash stays for rate-limiting. No email is
ever logged.
"""
from __future__ import annotations

import base64
import hashlib
import hmac
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

from fastapi import Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from pydantic import field_validator

from webapp.server.app import _build_or_error, _error_response
from webapp.server.estimator import classify
from webapp.server.i18n import t as _t
from webapp.server.mail import smtp_mailer
from webapp.server.runner import RunError
from webapp.server.schema import ComputeRequest
from webapp.server.security import hash_ip

_TEMPLATES = Jinja2Templates(
    directory=str(Path(__file__).resolve().parent.parent / "templates"))

# (url prefix, language) -- the verify page is mounted twice, same handler.
_LANGS = (("", "en"), ("/es", "es"))


class BigJobRequest(ComputeRequest):
    """A ``ComputeRequest`` plus the requester's email and preferred language.
    The email is shape-checked only (no delivery probe); the real validation is
    that the magic link is clicked from that inbox."""

    email: str
    lang: str = "en"

    @field_validator("email")
    @classmethod
    def _email_shape(cls, v: str) -> str:
        v = v.strip()
        if "@" not in v or not (3 <= len(v) <= 200):
            raise ValueError("a valid email is required for a big job")
        return v

    @field_validator("lang")
    @classmethod
    def _known_lang(cls, v: str) -> str:
        return v if v in ("en", "es") else "en"


# --------------------------------------------------------------------------- #
# Token: base64url(payload).base64url(HMAC-SHA256)
# --------------------------------------------------------------------------- #

def _b64(b: bytes) -> str:
    return base64.urlsafe_b64encode(b).decode().rstrip("=")


def _unb64(s: str) -> bytes:
    return base64.urlsafe_b64decode(s + "=" * (-len(s) % 4))


def _canonical(obj) -> str:
    """A byte-stable JSON encoding (sorted keys, no whitespace) so the same
    payload/spec always hashes identically across a JSON round-trip."""
    return json.dumps(obj, separators=(",", ":"), sort_keys=True)


def sign_token(secret: str, payload: dict) -> str:
    """Sign a payload dict into ``body.sig``. The body is a readable (base64url)
    encoding of the canonical payload; the signature is HMAC-SHA256 over it."""
    body = _b64(_canonical(payload).encode())
    sig = hmac.new(secret.encode(), body.encode(), hashlib.sha256).digest()
    return body + "." + _b64(sig)


def verify_token(secret: str, token: str, now: int) -> dict | None:
    """Return the payload iff the signature is valid (constant-time compare) and
    the token has not expired; otherwise ``None``. A malformed token, a bad
    signature, the wrong secret, or an ``exp`` at/after ``now`` all fail closed."""
    try:
        body, sig = token.split(".", 1)
    except ValueError:
        return None
    expected = _b64(hmac.new(secret.encode(), body.encode(), hashlib.sha256).digest())
    if not hmac.compare_digest(sig, expected):     # constant-time
        return None
    try:
        payload = json.loads(_unb64(body))
    except Exception:
        return None
    if int(payload.get("exp", 0)) < now:
        return None
    return payload


def _spec_hash(spec: dict) -> str:
    """A stable hash of the exact spec, binding a token to the computation it
    verifies (so a valid link cannot be repointed at a different job)."""
    return hashlib.sha256(_canonical(spec).encode()).hexdigest()


def _now() -> int:
    return int(datetime.now(timezone.utc).timestamp())


def _week_ago_iso() -> str:
    return (datetime.now(timezone.utc) - timedelta(days=7)).strftime(
        "%Y-%m-%dT%H:%M:%SZ")


def _email_hash(email: str, cfg) -> str:
    """Salted hash of the normalised address, used for rate-limiting and kept
    after the plaintext is deleted. The raw address never reaches a log or a
    long-lived column beyond the pending row / in-flight job."""
    return hash_ip(email.strip().lower(), cfg.token_secret)


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #

def register_big_jobs(app, cfg, store, mailer=None) -> None:
    """Mount ``POST /api/jobs/big`` (request a link) and ``GET /verify/{token}``
    plus its ``/es`` twin (consume a link). ``mailer`` overrides the default SMTP
    transport (tests inject a fake); it is ``None`` in production, where the
    worker/app wiring supplies the real ``smtp_mailer``."""

    @app.post("/api/jobs/big")
    def submit_big(req: BigJobRequest):
        # SMTP unconfigured => the whole tier is off; say "run locally" (spec §17).
        if not cfg.big_jobs_enabled:
            return JSONResponse(status_code=503, content={
                "error_type": "BigJobsDisabled",
                "message": "Large jobs are unavailable. Run locally: pip install quiverlab"})
        # Read the algebra dimension (library attribute -- never .dimension()) to
        # size the job; honest build errors surface as the same sanitised payload
        # the sync endpoint uses.
        try:
            A = _build_or_error(req.algebra)
            dim = A.dim
        except RunError as exc:
            return _error_response(exc.error_type, exc.message)
        info = classify(dim, req, cfg)
        if info["tier"] != "big":
            # Only genuinely over-anonymous-cap work belongs here.
            return JSONResponse(status_code=422, content={
                "error_type": "NotBig", "reason": info["reason"],
                "estimate": info["estimate"]})
        eh = _email_hash(req.email, cfg)
        # Per-email caps (on the HASH, not the address). Soft: the read-then-write
        # is non-atomic, acceptable for abuse control.
        if store.count_big_running_for_email_hash(eh) >= cfg.per_email_running_max:
            return JSONResponse(status_code=429, content={
                "error_type": "RateLimited",
                "message": "You already have a big job running."})
        if (store.count_big_since_for_email_hash(eh, _week_ago_iso())
                >= cfg.per_email_weekly_max):
            return JSONResponse(status_code=429, content={
                "error_type": "RateLimited", "message": "Weekly big-job limit reached."})
        spec = req.model_dump(by_alias=True, exclude={"email", "lang"})
        pid = store.create_pending_big(spec, req.email, eh, lang=req.lang)
        payload = {"pid": pid, "sh": _spec_hash(spec),
                   "exp": _now() + cfg.big_token_ttl_seconds}
        token = sign_token(cfg.token_secret, payload)
        prefix = "/es" if req.lang == "es" else ""
        url = cfg.public_base_url.rstrip("/") + prefix + "/verify/" + token
        send = mailer or smtp_mailer(cfg)
        send(req.email, _t("mail.verify_subject", req.lang),
             _t("mail.verify_body", req.lang).replace("{url}", url))
        return JSONResponse(status_code=202, content={"status": "sent"})

    for prefix, lang in _LANGS:
        _mount_verify(app, cfg, store, prefix, lang)


def _mount_verify(app, cfg, store, prefix: str, lang: str) -> None:
    # prefix/lang are closure constants (never route params) -- FastAPI does not
    # expose the language in the path, matching pages.py / feedback.py.
    @app.get(prefix + "/verify/{token}", response_class=HTMLResponse)
    def verify(request: Request, token: str):
        job_id = None
        payload = verify_token(cfg.token_secret, token, _now())
        # Check the big-queue cap BEFORE consuming, so a transiently full queue
        # does not burn a still-valid single-use link (the user can retry later).
        if payload is not None and store.count_big_pending() < cfg.big_queue_max:
            data = store.consume_pending_big(payload.get("pid", ""))   # single-use
            if (data is not None
                    and hmac.compare_digest(_spec_hash(data["spec"]),
                                            payload.get("sh", ""))):
                job_id = store.create_job(
                    data["spec"], ip=hash_ip("verified", cfg.ip_hash_salt),
                    tier="big", email=data["email"], email_hash=data["email_hash"],
                    wall_seconds=cfg.big_job_wall_seconds,
                    mem_bytes=cfg.big_job_mem_bytes, lang=data.get("lang", "en"))
        ctx = {"lang": lang, "prefix": prefix,
               "t": (lambda k: _t(k, lang)), "docs_url": cfg.docs_url,
               "other_url": (("/es" if lang == "en" else "") + "/verify/" + token),
               "job_id": job_id}
        return _TEMPLATES.TemplateResponse(request, "verify.html", ctx,
                                           status_code=(200 if job_id else 410))
