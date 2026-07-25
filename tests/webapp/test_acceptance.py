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
import sqlite3

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.cache import canonical_key, library_version
from webapp.server.config import Config
from webapp.server.store import JobStore
from webapp.worker.sweeper import sweep_cache_once, sweep_once
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


# --------------------------------------------------------------------------- #
# Plan 25 -- the result cache: never recompute a known example, across users and
# tiers. Email verification gates the COST of computing, not access to the maths.
# --------------------------------------------------------------------------- #

def _job_ids(cfg):
    conn = sqlite3.connect(cfg.db_path)
    try:
        return [r[0] for r in conn.execute("SELECT id FROM jobs").fetchall()]
    finally:
        conn.close()


def _pending_big_count(cfg):
    conn = sqlite3.connect(cfg.db_path)
    try:
        return conn.execute("SELECT COUNT(*) FROM pending_big").fetchone()[0]
    finally:
        conn.close()


def _queue_body(compute="hh_cohomology:0..4", pdf=True):
    return {"schema": 1, "algebra": _gf(2), "compute": [compute],
            "artifacts": {"pdf": pdf, "tikz": False}}


def test_queued_result_replayed_from_cache_across_users(tmp_path):
    # User A's queued computation, once finished, is replayed instantly for User B's
    # identical request -- zero recompute, no new job, no worker run needed.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_MAX_DEGREE": "2"})   # deg-4 -> queued
    client = TestClient(create_app(cfg))
    body = _queue_body()

    r1 = client.post("/api/compute", json=body)
    assert r1.status_code == 202 and r1.json()["tier"] == "queued"
    jid = r1.json()["job_id"]
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg) is True                    # runs + caches it
    assert len(_job_ids(cfg)) == 1

    # Identical request from "another user": a cache HIT pointing at the SAME job,
    # with no second job created and no worker tick required.
    r2 = client.post("/api/compute", json=body)
    assert r2.status_code == 200, r2.text
    assert r2.json()["tier"] == "cached"
    assert r2.json()["job_id"] == jid                         # the very same result
    assert r2.json()["cached_at"]                             # honest timestamp
    assert len(_job_ids(cfg)) == 1                            # NO duplicate job
    assert store.count_pending() == 0                         # nothing re-queued

    # The permalink page renders the bilingual "previously computed" note.
    page = client.get(f"/job/{jid}")
    assert page.status_code == 200
    assert "served instantly without recomputing" in page.text
    es = client.get(f"/es/job/{jid}")
    assert "sin recalcular" in es.text                        # ES parity


def test_instant_tier_checks_cache_before_computing(tmp_path):
    # A request that CLASSIFIES as instant is still served from the cache first,
    # skipping the instant compute entirely. Seed the cache via a queued run under
    # a tight instant cap, then re-ask under a generous cap (where it is instant):
    # the reply is the cached permalink, not a fresh instant result.
    body = _queue_body(compute="hh_cohomology:0..4", pdf=False)
    cfg_q = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                             "QLWEB_INSTANT_MAX_DEGREE": "2"})
    client_q = TestClient(create_app(cfg_q))
    jid = client_q.post("/api/compute", json=body).json()["job_id"]
    store = JobStore(cfg_q.db_path)
    assert worker_tick(store, cfg_q) is True

    cfg_i = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                             "QLWEB_INSTANT_MAX_DEGREE": "8"})   # now instant-classified
    client_i = TestClient(create_app(cfg_i))
    r = client_i.post("/api/compute", json=body)
    assert r.status_code == 200 and r.json()["tier"] == "cached"
    assert r.json()["job_id"] == jid
    assert len(_job_ids(cfg_i)) == 1                            # no fresh instant job either


def test_big_job_cache_hit_sends_no_email_and_mints_no_token(tmp_path):
    # THE big-job requirement: a big request whose canonical key is already cached
    # is served immediately, with NO email, NO token, NO pending row, NO email hash
    # -- across users. Only genuinely new big examples take the magic-link path.
    over = {"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_SMTP_HOST": "relay",
            "QLWEB_SMTP_FROM": "q@e.org", "QLWEB_PUBLIC_BASE_URL": "https://ql.example",
            "QLWEB_INSTANT_OPS_THRESHOLD": "0", "QLWEB_QUEUED_OPS_THRESHOLD": "0",
            "QLWEB_INSTANT_MAX_DEGREE": "0", "QLWEB_QUEUED_MAX_DEGREE": "0"}
    cfg = Config.from_env(over)
    body = {"schema": 1, "algebra": _gf(2), "compute": ["hh_cohomology:0..3"],
            "artifacts": {"pdf": False, "tikz": False},
            "email": "first@example.org", "lang": "en"}

    # First user: the full magic-link flow computes and caches the big example.
    mail1 = FakeMailer()
    client1 = TestClient(create_app(cfg, mailer=mail1))
    assert client1.post("/api/jobs/big", json=body).status_code == 202
    token = re.search(r"/verify/([^\s]+)", mail1.sent[0][2]).group(1)
    assert client1.get("/verify/" + token).status_code == 200
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg, mailer=mail1) is True
    first_job = _job_ids(cfg)[0]

    # Second user, FRESH mailer, same DB, identical big request: served from cache.
    mail2 = FakeMailer()
    client2 = TestClient(create_app(cfg, mailer=mail2))
    r = client2.post("/api/jobs/big", json=dict(body, email="second@example.org"))
    assert r.status_code == 200, r.text
    assert r.json()["status"] == "cached"
    assert r.json()["job_id"] == first_job                     # replays the SAME result
    # The crux, asserted via the fake mailer + the store:
    assert mail2.sent == []                                    # NO verification email
    assert _pending_big_count(cfg) == 0                        # NO token / pending row minted
    assert _job_ids(cfg) == [first_job]                        # NO duplicate big job

    # Privacy: neither the cached result nor the cache row leaks the first requester.
    stripped = {k: body[k] for k in ("schema", "algebra", "compute", "artifacts")}
    row = store.cache_row(canonical_key(stripped, library_version()))
    assert "first@example.org" not in json.dumps(row)
    result_json = (cfg.artifacts_dir / first_job / "result.json").read_text(encoding="utf-8")
    assert "first@example.org" not in result_json


def test_big_cache_hit_served_even_when_smtp_disabled(tmp_path):
    # Email gates COST, not access: with SMTP off (big tier "disabled"), an ALREADY
    # cached example is still served -- a cached big example needs no relay. Seed the
    # cache with an anonymous queued run of the same maths, then POST /api/jobs/big.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_MAX_DEGREE": "2"})     # no SMTP -> big off
    assert cfg.big_jobs_enabled is False
    client = TestClient(create_app(cfg))
    base = {"schema": 1, "algebra": _gf(2), "compute": ["hh_cohomology:0..4"],
            "artifacts": {"pdf": False, "tikz": False}}
    jid = client.post("/api/compute", json=base).json()["job_id"]
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg) is True                      # cached

    # A big request for the SAME maths would 503 (SMTP off) -- but it is cached, so
    # it is served instead. The cache check runs BEFORE the disabled-tier refusal.
    r = client.post("/api/jobs/big", json=dict(base, email="u@e.org", lang="en"))
    assert r.status_code == 200, r.text
    assert r.json()["status"] == "cached" and r.json()["job_id"] == jid
    assert _pending_big_count(cfg) == 0


def test_version_bump_invalidates_the_cache(tmp_path, monkeypatch):
    # A library-version bump changes every key, so a previously cached example is a
    # MISS and recomputes -- the cache never replays a stale-version result.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_MAX_DEGREE": "2"})
    client = TestClient(create_app(cfg))
    body = _queue_body(pdf=False)
    jid = client.post("/api/compute", json=body).json()["job_id"]
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg) is True                      # cached under the live version

    # Same request under the same version: a hit.
    assert client.post("/api/compute", json=body).json()["tier"] == "cached"

    # Bump the version: the key changes, so the same request now MISSES and queues
    # a brand-new job instead of replaying the (now stale-version) cached one.
    import quiverlab
    monkeypatch.setattr(quiverlab, "__version__", "99.0.0-bumped")
    r = client.post("/api/compute", json=body)
    assert r.status_code == 202 and r.json()["tier"] == "queued"
    assert r.json()["job_id"] != jid                            # a fresh computation


def test_cache_lru_sweep_evicts_and_retention_reclaims(tmp_path):
    # The LRU size-cap sweep runs alongside retention. Fill the cache past a cap of
    # 1, sweep: the least-recently-hit entry is evicted (unpinned), and the ordinary
    # retention sweep then reclaims that now-unpinned old job -- while the surviving
    # cached job's artifacts are kept.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_CACHE_MAX_ENTRIES": "1"})
    store = JobStore(cfg.db_path)
    store.init_schema()

    def _cached_done_job(tag):
        jid = store.create_job({"compute": [tag]}, ip="i")
        store.claim_next()
        store.mark_done(jid, str(cfg.artifacts_dir / jid))
        (cfg.artifacts_dir / jid).mkdir(parents=True, exist_ok=True)
        return jid

    cold = _cached_done_job("cartan")
    store.cache_put("cold", cold, library_version(), "2026-01-01T00:00:00Z")
    warm = _cached_done_job("dimension")
    store.cache_put("warm", warm, library_version(), "2026-06-01T00:00:00Z")

    assert sweep_cache_once(store, cfg) == 1                     # cap=1 -> evict the cold one
    assert store.cache_row("cold") is None and store.cache_row("warm") is not None

    # Retention now reclaims the unpinned (evicted, old) job; the still-cached job
    # AND its artifacts survive.
    removed = sweep_once(store, cfg, now_iso="2999-01-01T00:00:00Z")
    assert cold in removed and warm not in removed
    assert not (cfg.artifacts_dir / cold).exists()              # evicted job's artifacts gone
    assert (cfg.artifacts_dir / warm).exists()                  # cached job's artifacts kept


# --------------------------------------------------------------------------- #
# Plan 26 -- no-code module input: the module compute kinds served end-to-end on
# all three tiers (instant / queued+cached / big+cached), with the loud
# RelationError surfaced as a clean 4xx.
# --------------------------------------------------------------------------- #

# k[x]/(x^3) over GF(2): one vertex + a nilpotent loop -- a tiny module algebra.
_MOD_ALG = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
            "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}
_M2 = {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}


def _module_body(compute, module=None, **extra):
    body = {"schema": 2, "algebra": _MOD_ALG, "compute": compute,
            "artifacts": {"pdf": False, "tikz": False},
            "module": module if module is not None else _M2}
    body.update(extra)
    return body


def test_module_instant_compute(tmp_path):
    # A small module dispatches on the instant tier and carries its result + the
    # resolved references inline -- no job, no artifact.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute",
                    json=_module_body(["dimension_vector", "rad_top_soc", "tau"]))
    assert r.status_code == 200, r.text
    assert r.json()["tier"] == "instant"
    results = r.json()["result"]["results"]
    assert results["dimension_vector"]["dimvec"] == {"1": 2}
    assert results["rad_top_soc"]["socle"]["dimvec"] == {"1": 1}
    assert results["tau"]["dimvec"] == {"1": 2}
    # References carry through like every other kind.
    keys = {e["key"] for e in r.json()["result"]["references"]}
    assert "assem_book" in keys


def test_module_relation_violation_is_a_clean_4xx(tmp_path):
    # A module whose matrices break the relations (x acting invertibly => x^3 != 0)
    # gets the library's loud error as a clean 4xx -- never a 500, never a silent
    # wrong answer.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute",
                    json=_module_body(["dimension_vector"],
                                      module={"dims": {"1": 1}, "maps": {"x": [[1]]}}))
    assert r.status_code == 422, r.text
    assert "relation" in r.json()["message"].lower()


def test_module_ext_two_modules_instant(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute",
                    json=_module_body(["ext:0..3"],
                                      ext_target={"builtin": {"kind": "simple", "vertex": 1}}))
    assert r.status_code == 200, r.text
    assert r.json()["result"]["results"]["ext"]["dims"] == [1, 1, 1, 1]


def test_queued_module_replayed_from_cache_across_users(tmp_path):
    # A deep module resolution routes to the queued tier; once computed it is
    # replayed instantly for an identical request, no recompute, no second job.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_MAX_DEGREE": "2"})    # deg-4 -> queued
    client = TestClient(create_app(cfg))
    body = _module_body(["projective_resolution:0..4"])

    r1 = client.post("/api/compute", json=body)
    assert r1.status_code == 202 and r1.json()["tier"] == "queued", r1.text
    jid = r1.json()["job_id"]
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg) is True                     # runs + caches it

    result = json.loads((cfg.artifacts_dir / jid / "result.json").read_text(encoding="utf-8"))
    assert result["results"]["projective_resolution"]["top"] == 4

    r2 = client.post("/api/compute", json=body)
    assert r2.status_code == 200 and r2.json()["tier"] == "cached"
    assert r2.json()["job_id"] == jid                          # the very same result
    assert len(_job_ids(cfg)) == 1                             # no duplicate job


def test_big_module_cache_hit_sends_no_email(tmp_path):
    # The big-job requirement, for a MODULE request: a first user runs it through
    # the magic-link flow (which caches it); a second user's identical module
    # request is served from cache with NO email, NO token, NO pending row.
    over = {"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_SMTP_HOST": "relay",
            "QLWEB_SMTP_FROM": "q@e.org", "QLWEB_PUBLIC_BASE_URL": "https://ql.example",
            "QLWEB_INSTANT_OPS_THRESHOLD": "0", "QLWEB_QUEUED_OPS_THRESHOLD": "0",
            "QLWEB_INSTANT_MAX_DEGREE": "0", "QLWEB_QUEUED_MAX_DEGREE": "0"}
    cfg = Config.from_env(over)
    body = _module_body(["dimension_vector", "tau"], email="first@example.org", lang="en")

    mail1 = FakeMailer()
    client1 = TestClient(create_app(cfg, mailer=mail1))
    assert client1.post("/api/jobs/big", json=body).status_code == 202, "expected the magic-link flow"
    token = re.search(r"/verify/([^\s]+)", mail1.sent[0][2]).group(1)
    assert client1.get("/verify/" + token).status_code == 200
    store = JobStore(cfg.db_path)
    assert worker_tick(store, cfg, mailer=mail1) is True
    first_job = _job_ids(cfg)[0]

    # Second user, fresh mailer, identical module request: served from cache.
    mail2 = FakeMailer()
    client2 = TestClient(create_app(cfg, mailer=mail2))
    r = client2.post("/api/jobs/big", json=dict(body, email="second@example.org"))
    assert r.status_code == 200, r.text
    assert r.json()["status"] == "cached" and r.json()["job_id"] == first_job
    assert mail2.sent == []                                     # NO verification email
    assert _pending_big_count(cfg) == 0                        # NO token / pending row
    assert _job_ids(cfg) == [first_job]                        # NO duplicate big job

    # Privacy: the cached module result leaks nothing about the first requester.
    result_json = (cfg.artifacts_dir / first_job / "result.json").read_text(encoding="utf-8")
    assert "first@example.org" not in result_json
    assert json.loads(result_json)["results"]["tau"]["dimvec"] == {"1": 2}
