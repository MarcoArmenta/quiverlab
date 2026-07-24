"""Task 15 -- full-stack acceptance, in-process via ``TestClient``.

Adapted from the task brief to the REAL delivered stack (the brief's inline code
was written against a stale surface):

  * Reference family is ``QuantumCI(q=1)`` -- the brief's ``truncated_polynomial``
    is NOT in ``quiverlab.families()`` and would never build (the same 2026-07-24
    amendment already applied to ``test_api``/``test_worker``/``test_bigjobs``).
    ``QuantumCI(q=1)`` has dim 4; over GF(2) ``hh_cohomology`` dims are
    [4, 8, 12, 16, ...].
  * Queueing is forced with ``QLWEB_INSTANT_MAX_DEGREE=2`` and a degree-4 request
    (the real estimator's degree cap), matching ``test_api``'s pattern.
  * Result artifact shape is ``result["results"]["hh_cohomology"]["dims"]`` (a
    list), the artifact download route is ``GET /download/{jid}/result.json``, and
    the job status route ``GET /api/jobs/{jid}`` returns ``{"status": ...}``.
  * The Spanish index lives at ``/es`` (canonical; ``/es/`` 307-redirects to it)
    and carries the invariant label "Polinomio de Coxeter".

This exercises the whole stack end-to-end WITHOUT Docker or a network: the app
factory, the SQLite job store, the resource-capped ``spawn`` worker child, the
instant wall-time net, the bilingual pages, the feedback round-trip, and the
big-job magic-link tier (with an injected mailer -- no SMTP is ever touched).

Explicitly marked ``fast``: this file's basename collides with the deep-bucket
``test_acceptance.py`` heuristic in ``tests/conftest.py`` (meant for the heavy
engine acceptance tests), so an explicit bucket marker is required to keep the
webapp acceptance test in the fast leg alongside its siblings. Its subprocess
computations are small QuantumCI(q=1) runs -- the same ones the other fast webapp
tests already spawn.

A real multi-process smoke (uvicorn + a worker loop as OS subprocesses, driven
over HTTP) lives at ``scripts/webapp_smoke.py``; it is not a pytest so that port
allocation and process teardown never make CI flaky.
"""
import json
import re

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.store import JobStore
from webapp.worker.worker import worker_tick

pytestmark = pytest.mark.fast

# QuantumCI(q=1): dim 4, nondegenerate over GF(2)/GF(3); the real reference
# family throughout the webapp suite.
_QCI = {"kind": "family", "family": "QuantumCI", "params": {"q": 1}}


def _gf(p=2):
    return dict(_QCI, field={"kind": "GF", "p": p, "n": 1})


class FakeMailer:
    """Captures (to, subject, body) instead of sending -- the network is never
    touched, and the captured tuples let the big-job test read back the link."""

    def __init__(self):
        self.sent = []

    def __call__(self, to, subject, body):
        self.sent.append((to, subject, body))


def test_end_to_end_queued_job(tmp_path):
    # INSTANT_MAX_DEGREE=2 forces a degree-4 request into the queued tier.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_MAX_DEGREE": "2"})
    client = TestClient(create_app(cfg))
    body = {"schema": 1, "algebra": _gf(2),
            "compute": ["hh_cohomology:0..4"],
            "artifacts": {"pdf": True, "tikz": False}}
    r = client.post("/api/compute", json=body)
    assert r.status_code == 202, r.text
    assert r.json()["tier"] == "queued"
    jid = r.json()["job_id"]

    # Run the pending job in a real resource-capped spawn child.
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg) is True

    status = client.get(f"/api/jobs/{jid}").json()
    assert status["status"] == "done", status.get("error")

    dl = client.get(f"/download/{jid}/result.json")
    assert dl.status_code == 200
    result = json.loads(dl.content)
    dims = result["results"]["hh_cohomology"]["dims"]
    assert isinstance(dims, list) and dims == [4, 8, 12, 16, 20]

    # The worked-steps artifact (pdf=True) is downloadable -- PDF when a LaTeX
    # toolchain is present, else the self-contained HTML fallback.
    trace = client.get(f"/download/{jid}/trace.pdf")
    if trace.status_code != 200:
        trace = client.get(f"/download/{jid}/trace_steps.html")
    assert trace.status_code == 200, "worked-steps artifact should be downloadable"


def test_instant_path_end_to_end(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    body = {"schema": 1, "algebra": _gf(3),
            "compute": ["coxeter_polynomial"],
            "artifacts": {"pdf": False, "tikz": False}}
    r = client.post("/api/compute", json=body)
    assert r.status_code == 200, r.text
    assert r.json()["tier"] == "instant"
    # The instant result carries the computed block inline (no job, no artifact).
    assert "coxeter_polynomial" in r.json()["result"]["results"]


def test_spanish_page_and_feedback_roundtrip(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    es = client.get("/es")                       # canonical (/es/ 307s here)
    assert es.status_code == 200 and "Polinomio de Coxeter" in es.text
    fb = client.post("/api/feedback", json={
        "category": "feature", "message": "Please add support for gentle algebras.",
        "contact": "", "job_ref": "", "website": ""})
    assert fb.status_code == 201 and fb.json()["reference"]


def test_literature_page_and_suggestion(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    lit = client.get("/literature")
    assert lit.status_code == 200                # renders the curated bibliography
    r = client.post("/api/feedback", json={
        "category": "literature", "message": "Consider adding this reference.",
        "reference": "arXiv:1610.05741",
        "why_relevant": "Volkov homotopy liftings, used by the Gerstenhaber bracket.",
        "contact": "", "job_ref": "", "website": ""})
    assert r.status_code == 201 and r.json()["reference"]


def test_big_job_end_to_end_mocked_smtp(tmp_path):
    # Force a TINY computation into the big tier by dropping every anonymous
    # threshold to 0, so the whole magic-link flow exercises end-to-end while the
    # actual job stays fast (no real 4-hour big computation in CI).
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org",
                           "QLWEB_PUBLIC_BASE_URL": "https://ql.example",
                           "QLWEB_INSTANT_OPS_THRESHOLD": "0",
                           "QLWEB_QUEUED_OPS_THRESHOLD": "0",
                           "QLWEB_INSTANT_MAX_DEGREE": "0",
                           "QLWEB_QUEUED_MAX_DEGREE": "0"})
    mail = FakeMailer()
    client = TestClient(create_app(cfg, mailer=mail))
    body = {"schema": 1, "algebra": _gf(2),
            "compute": ["hh_cohomology:0..3"], "artifacts": {"pdf": False, "tikz": False},
            "email": "user@example.org", "lang": "en"}
    assert client.post("/api/jobs/big", json=body).status_code == 202
    token = re.search(r"/verify/([^\s]+)", mail.sent[0][2]).group(1)
    assert client.get("/verify/" + token).status_code == 200

    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg, mailer=mail) is True      # runs the queued big job
    assert any(s[0] == "user@example.org" and "finished" in s[1].lower()
               for s in mail.sent[1:])                       # completion email sent


def test_big_tier_disabled_without_smtp(tmp_path):
    # The SMTP-off leg of the smoke: with no relay configured the whole big-job
    # tier is unavailable and the app says "run locally" (a clean 503), never a
    # hang or a crash.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    assert cfg.big_jobs_enabled is False
    client = TestClient(create_app(cfg))
    body = {"schema": 1, "algebra": {"kind": "family", "family": "QuantumCI",
                                     "params": {"q": 1}, "field": {"kind": "CC"}},
            "compute": ["hh_cohomology:0..30"], "artifacts": {"pdf": False, "tikz": False},
            "email": "user@example.org", "lang": "en"}
    r = client.post("/api/jobs/big", json=body)
    assert r.status_code == 503
    assert r.json()["error_type"] == "BigJobsDisabled"
