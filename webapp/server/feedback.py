"""Feedback: a bilingual form page, a JSON submit endpoint with honeypot +
rate limiting, and a token-gated admin table. No accounts, no emails (v1).

Abuse controls (spec global constraints): a filled honeypot is dropped SILENTLY
(a 201-shaped response, nothing stored); the message is 10-4000 chars; the daily
per-hashed-IP cap is enforced (429 over it); `job_ref` is optional but
ULID-validated when present. The client address is hashed (reusing the app's
`client_ip` + `security.hash_ip`) BEFORE any store call -- the raw address never
reaches the store or a log. The admin table exists only when a token is
configured and compares it in CONSTANT TIME (`hmac.compare_digest`); the only
contacts it shows are the feedback `contact` field -- never a job email.
"""
from __future__ import annotations

import hmac
import html
import json
from pathlib import Path
from typing import Literal

from fastapi import Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field, field_validator, model_validator

from webapp.server.app import _now_iso, client_ip
from webapp.server.i18n import t as _t
from webapp.server.limits import check_feedback_allowed
from webapp.server.security import hash_ip, valid_ulid

_TEMPLATES = Jinja2Templates(
    directory=str(Path(__file__).resolve().parent.parent / "templates"))

# (url prefix, language) -- the page is mounted twice, same handler.
_LANGS = (("", "en"), ("/es", "es"))
GITHUB_ISSUES_URL = "https://github.com/MarcoArmenta/quiverlab/issues"


class FeedbackRequest(BaseModel):
    category: Literal["problem", "feature", "literature"]
    message: str = Field(min_length=10, max_length=4000)
    contact: str | None = Field(default=None, max_length=200)
    job_ref: str | None = Field(default=None, max_length=64)
    # literature-category structured fields (required only for that category):
    reference: str | None = Field(default=None, max_length=2000)
    why_relevant: str | None = Field(default=None, max_length=2000)
    website: str = ""                       # honeypot: real users never fill this

    @field_validator("message")
    @classmethod
    def _nonblank(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("message must not be blank")
        return v

    @field_validator("contact")
    @classmethod
    def _blank_contact_is_none(cls, v: str | None) -> str | None:
        # A blank contact is stored as NULL, not an empty string.
        if v is None:
            return None
        return v.strip() or None

    @field_validator("job_ref")
    @classmethod
    def _job_ref_ulid(cls, v: str | None) -> str | None:
        # Optional, but a present job_ref must be a well-formed ULID -- a
        # malformed one can never name a real job (global constraint).
        if v is None:
            return None
        v = v.strip()
        if not v:
            return None
        if not valid_ulid(v):
            raise ValueError("job_ref must be a valid ULID when provided")
        return v

    @model_validator(mode="after")
    def _literature_fields(self):
        if self.category == "literature":
            if len((self.reference or "").strip()) < 4:
                raise ValueError("literature: reference must be 4-2000 chars")
            if len((self.why_relevant or "").strip()) < 10:
                raise ValueError("literature: why-relevant must be 10-2000 chars")
        return self


def register_feedback(app, cfg, store) -> None:
    for prefix, lang in _LANGS:
        _mount_page(app, cfg, prefix, lang)

    @app.post("/api/feedback")
    def api_feedback(req: FeedbackRequest, request: Request):
        if req.website.strip():             # honeypot: pretend success, store nothing
            return JSONResponse(status_code=201, content={"reference": "ignored"})
        # Hash the address immediately: the raw value never reaches the store.
        iph = hash_ip(client_ip(request), cfg.ip_hash_salt)
        refusal = check_feedback_allowed(store, cfg, iph, _now_iso())
        if refusal:
            return JSONResponse(status_code=429,
                                content={"error_type": "RateLimited", "message": refusal})
        extra = None
        if req.category == "literature":
            extra = json.dumps({"reference": req.reference,
                                "why_relevant": req.why_relevant})
        fid = store.create_feedback(req.category, req.message, req.contact,
                                    iph, req.job_ref, extra)
        return JSONResponse(status_code=201, content={"reference": fid})

    if cfg.admin_token:                     # route exists only when a token is configured
        @app.get("/admin/feedback", response_class=HTMLResponse)
        def admin_feedback(token: str = ""):
            # Constant-time compare so a wrong token leaks no timing signal.
            # 403 (not 401): the token is a query-param authorization, not an
            # HTTP auth scheme, so no WWW-Authenticate challenge is owed.
            if not hmac.compare_digest(token, cfg.admin_token):
                return HTMLResponse("forbidden", status_code=403)
            rows = store.list_feedback()
            out = ["<h1>Feedback</h1><table border=\"1\" cellpadding=\"4\">",
                   "<tr><th>when</th><th>category</th><th>message</th>"
                   "<th>contact</th><th>job</th><th>extra</th></tr>"]
            for r in rows:
                # Every cell escaped: feedback text is fully untrusted input.
                out.append(
                    "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                        html.escape(r["created_at"]), html.escape(r["category"]),
                        html.escape(r["message"]), html.escape(r.get("contact") or ""),
                        html.escape(r.get("job_ref") or ""), html.escape(r.get("extra") or "")))
            out.append("</table>")
            return HTMLResponse("".join(out))


def _mount_page(app, cfg, prefix: str, lang: str) -> None:
    # lang/prefix are closure constants, never route parameters -- FastAPI does
    # not expose the language in the path (same pattern as pages.py).
    @app.get(prefix + "/feedback", response_class=HTMLResponse)
    def feedback_page(request: Request, job: str = ""):
        # Only prefill a well-formed ULID (the API validates it anyway; this
        # keeps an arbitrary ?job=... out of the rendered field).
        job_ref = job if valid_ulid(job) else ""
        return _TEMPLATES.TemplateResponse(
            request, "feedback.html",
            {"lang": lang, "prefix": prefix,
             "t": (lambda k: _t(k, lang)),
             "other_url": ("/feedback" if lang == "es" else "/es/feedback"),
             "docs_url": cfg.docs_url,
             "job_ref": job_ref, "github": GITHUB_ISSUES_URL})
