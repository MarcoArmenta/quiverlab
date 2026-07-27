"""Task 9 -- FastAPI app tests.

Adapted per the 2026-07-24 amendment and the REAL in-repo interfaces (the brief's
inline test code was written against a stale surface):
  * reference family is ``QuantumCI(q=1)`` -- ``truncated_polynomial`` is NOT in
    ``families()`` (amendment pt 3; ``QuantumCI`` q=1 is nondegenerate over GF(2)).
  * the estimator (Task 5) caps the anonymous ``queued`` band at degree 20 and
    turns an over-anonymous CC/degree-40 request into ``reject``/``big_disabled``
    (NOT ``queued``). So ``test_deep_request_becomes_job`` uses a genuine
    queued-band request (degree 12), exactly as Task 5's own
    ``test_deep_request_is_queued`` does.
  * the algebra dimension is the attribute ``A.dim`` (amendment pt 4).

Beyond the brief's 10 tests this file also pins the adjudicated decisions:
error-type genericisation, salted-IP storage, security headers on API AND error
responses, ULID rejection, and rate-limit 429s.
"""
import sqlite3

from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.security import (
    GENERIC_ERROR_MESSAGE, hash_ip, sanitize_error, valid_ulid,
)
from webapp.server.store import JobStore


FAMILY = "QuantumCI"


def _cfg(tmp_path, **extra):
    env = {"QLWEB_DATA_DIR": str(tmp_path)}
    env.update({k: str(v) for k, v in extra.items()})
    return Config.from_env(env)


def _client(tmp_path, **extra):
    return TestClient(create_app(_cfg(tmp_path, **extra)))


def _gf_body(compute, pdf=False):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": FAMILY,
                        "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": compute, "artifacts": {"pdf": pdf, "tikz": False}}


def _cc_body(compute, pdf=False):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": FAMILY,
                        "params": {"q": 1}, "field": {"kind": "CC"}},
            "compute": compute, "artifacts": {"pdf": pdf, "tikz": False}}


# --------------------------------------------------------------------------- #
# Brief's 10 tests (adapted to the real family + real estimator semantics)
# --------------------------------------------------------------------------- #

def test_catalog_endpoint(tmp_path):
    r = _client(tmp_path).get("/api/catalog")
    assert r.status_code == 200
    assert any(f["name"] == FAMILY for f in r.json()["families"])


def test_instant_compute(tmp_path):
    r = _client(tmp_path).post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r.status_code == 200, r.text
    assert r.json()["tier"] == "instant"
    assert "hh_cohomology" in r.json()["result"]["results"]


def test_deep_request_becomes_job(tmp_path):
    # Degree over the instant cap (8) but within the anonymous queued band (<=20)
    # -> a queued anonymous job. (The brief's original degree-40 CC request is
    # `reject`/`big_disabled` under the real estimator, not `queued`.)
    r = _client(tmp_path).post("/api/compute", json=_cc_body(["hh_cohomology:0..12"], pdf=True))
    assert r.status_code == 202
    assert r.json()["tier"] == "queued"
    assert r.json()["job_id"]


def test_failure_surfaces_error_type(tmp_path):
    body = {"schema": 1,
            "algebra": {"kind": "family", "family": "does_not_exist",
                        "params": {}, "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    r = _client(tmp_path).post("/api/compute", json=body)
    assert r.status_code == 422
    assert r.json()["error_type"] in ("CatalogError", "FieldError")


def test_wall_net_converts_instant_to_job(tmp_path):
    # instant_wall_seconds=0 makes every instant run "time out", exercising the
    # hard net: an under-estimated instant computation becomes a queued job.
    client = _client(tmp_path, QLWEB_INSTANT_WALL_SECONDS=0)
    r = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r.status_code == 202
    assert r.json()["tier"] == "queued"
    assert r.json()["job_id"]


def test_security_headers_present(tmp_path):
    r = _client(tmp_path).get("/api/catalog")
    csp = r.headers["content-security-policy"]
    assert "default-src 'self'" in csp
    assert "script-src 'self'" in csp          # no 'unsafe-inline'
    assert "'unsafe-inline'" not in csp.split("style-src")[0]  # scripts stay strict
    assert r.headers["x-content-type-options"] == "nosniff"
    assert r.headers["x-frame-options"] == "DENY"
    assert r.headers["referrer-policy"]


def test_stored_ip_is_hashed_not_raw(tmp_path):
    cfg = _cfg(tmp_path, QLWEB_INSTANT_MAX_DEGREE=0)
    client = TestClient(create_app(cfg))
    # TestClient's default client host is "testclient"; assert it is never stored raw.
    r = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..2"]))
    assert r.status_code == 202
    rows = sqlite3.connect(cfg.db_path).execute("SELECT ip FROM jobs").fetchall()
    assert rows and all(ip != "testclient" and len(ip) == 64 for (ip,) in rows)


def test_queue_cap_returns_429(tmp_path):
    client = _client(tmp_path, QLWEB_INSTANT_MAX_DEGREE=0,   # force queueing
                     QLWEB_GLOBAL_QUEUE_MAX=1, QLWEB_PER_IP_RUNNING_MAX=5)
    body = _gf_body(["hh_cohomology:0..2"])
    assert client.post("/api/compute", json=body).status_code == 202   # fills the queue
    r2 = client.post("/api/compute", json=body)
    assert r2.status_code == 429
    assert r2.json()["error_type"] == "RateLimited"


def test_over_anonymous_cap_offers_big_with_estimate(tmp_path):
    # SMTP set -> big-job tier enabled -> create_app's production-secret guard
    # (correction #4) needs a real (non-default) signing secret.
    client = _client(tmp_path, QLWEB_SMTP_HOST="relay", QLWEB_SMTP_FROM="q@e.org",
                     QLWEB_TOKEN_SECRET="test-secret-not-the-default")
    r = client.post("/api/compute", json=_cc_body(["hh_cohomology:0..30"]))
    assert r.status_code == 202          # NOT a hard reject
    body = r.json()
    assert body["tier"] == "big"
    assert body["estimate"]["cells"] > 0 and body["estimate"]["minutes"] >= 1


def test_beyond_big_cap_returns_structured_422(tmp_path):
    client = _client(tmp_path, QLWEB_SMTP_HOST="relay", QLWEB_SMTP_FROM="q@e.org",
                     QLWEB_TOKEN_SECRET="test-secret-not-the-default")
    r = client.post("/api/compute", json=_cc_body(["hh_cohomology:0..200"]))
    assert r.status_code == 422
    assert r.json()["reason"] == "beyond_big_cap"
    assert r.json()["estimate"]["cells"] > 0


# --------------------------------------------------------------------------- #
# Adjudicated-decision tests (required by the task instructions)
# --------------------------------------------------------------------------- #

def test_hash_ip_is_salted_sha256():
    import hashlib
    h = hash_ip("1.2.3.4", "pepper")
    assert h == hashlib.sha256(b"pepper:1.2.3.4").hexdigest()
    assert len(h) == 64
    # Salt changes the digest; the raw address is not recoverable from it.
    assert hash_ip("1.2.3.4", "other") != h
    assert "1.2.3.4" not in h


def test_sanitize_error_passes_safe_and_generices_unexpected():
    # QuiverlabError subclass names surface verbatim.
    assert sanitize_error("RelationError", "bad relation") == ("RelationError", "bad relation")
    assert sanitize_error("FieldError", "not a prime") == ("FieldError", "not a prime")
    # The runner's own honest domain tags surface verbatim.
    assert sanitize_error("CatalogError", "unknown family") == ("CatalogError", "unknown family")
    assert sanitize_error("ResultTooLarge", "too big")[0] == "ResultTooLarge"
    # An UNEXPECTED internal exception type name is genericised (nothing leaks).
    etype, msg = sanitize_error("KeyError", "'internal_dict_key'")
    assert etype == "InternalError"
    assert "KeyError" not in etype and "internal_dict_key" not in msg


def test_unexpected_error_is_genericised(tmp_path, monkeypatch):
    # The runner leaks type+str into RunError for server-log fidelity; the API
    # must NOT pass an unexpected type/message through. Force build_algebra to
    # raise a raw exception carrying a secret-looking message.
    import webapp.server.app as appmod

    def _boom(_spec):
        raise KeyError("SECRET_INTERNAL_DETAIL")

    monkeypatch.setattr(appmod, "build_algebra", _boom)
    r = _client(tmp_path).post("/api/compute", json=_gf_body(["hh_cohomology:0..2"]))
    assert r.status_code == 500
    assert r.json()["error_type"] == "InternalError"
    assert "KeyError" not in r.text
    assert "SECRET_INTERNAL_DETAIL" not in r.text


def test_security_headers_on_error_responses(tmp_path):
    client = _client(tmp_path)
    # 422 error response carries the headers.
    bad = {"schema": 1,
           "algebra": {"kind": "family", "family": "does_not_exist",
                       "params": {}, "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    r422 = client.post("/api/compute", json=bad)
    assert r422.status_code == 422
    assert "default-src 'self'" in r422.headers["content-security-policy"]
    # 404 error response carries the headers too.
    r404 = client.get("/api/jobs/" + "0" * 26)   # valid ULID shape, absent
    assert r404.status_code == 404
    assert r404.headers["x-content-type-options"] == "nosniff"


def test_get_job_rejects_non_ulid(tmp_path):
    client = _client(tmp_path)
    for bad in ("not-a-ulid", "O" * 26, "1" * 25, "I" * 26, "%2e%2e"):
        r = client.get(f"/api/jobs/{bad}")
        assert r.status_code == 404, (bad, r.status_code)
    # The validator itself is strict Crockford / length 26.
    assert not valid_ulid("not-a-ulid")
    assert not valid_ulid("O" * 26)      # 'O' is not in the Crockford alphabet
    assert valid_ulid("0" * 26)


def test_create_and_get_job_roundtrip(tmp_path):
    client = _client(tmp_path, QLWEB_INSTANT_MAX_DEGREE=0)
    r = client.post("/api/jobs", json=_gf_body(["hh_cohomology:0..2"]))
    assert r.status_code == 202
    jid = r.json()["job_id"]
    assert valid_ulid(jid)
    got = client.get(f"/api/jobs/{jid}")
    assert got.status_code == 200
    assert got.json()["id"] == jid
    assert got.json()["status"] == "pending"


def test_jobs_endpoint_rate_limited(tmp_path):
    client = _client(tmp_path, QLWEB_GLOBAL_QUEUE_MAX=1, QLWEB_PER_IP_RUNNING_MAX=5)
    body = _gf_body(["hh_cohomology:0..2"])
    assert client.post("/api/jobs", json=body).status_code == 202
    r2 = client.post("/api/jobs", json=body)
    assert r2.status_code == 429
    assert r2.json()["error_type"] == "RateLimited"


def test_get_missing_job_returns_404(tmp_path):
    r = _client(tmp_path).get("/api/jobs/" + "7" * 26)   # valid shape, never created
    assert r.status_code == 404


def test_unhandled_exception_gets_security_headers_and_no_leak(tmp_path, monkeypatch):
    # A GENUINELY unhandled exception -- one that escapes the route (not a
    # RunError caught into an honest error response) -- must still ship ALL the
    # security headers and leak neither the exception type nor its message.
    # Force the store's get_job to raise an arbitrary exception with a secret.
    from webapp.server.store import JobStore

    def _boom(self, job_id):
        raise RuntimeError("SECRET_INTERNAL_LEAK")

    monkeypatch.setattr(JobStore, "get_job", _boom)
    # raise_server_exceptions=False so the 500 response is observable.
    client = TestClient(create_app(_cfg(tmp_path)), raise_server_exceptions=False)
    r = client.get("/api/jobs/" + "0" * 26)   # valid ULID -> reaches get_job
    assert r.status_code == 500
    # All four security headers present on the genuinely-unhandled 500.
    assert "default-src 'self'" in r.headers["content-security-policy"]
    assert r.headers["x-content-type-options"] == "nosniff"
    assert r.headers["x-frame-options"] == "DENY"
    assert r.headers["referrer-policy"] == "no-referrer"
    # Genericised: neither the exception type nor its (secret) message leaks.
    assert r.json()["error_type"] == "InternalError"
    assert "RuntimeError" not in r.text
    assert "SECRET_INTERNAL_LEAK" not in r.text


def test_async_job_error_is_genericised_at_read_boundary(tmp_path):
    # The queued/async path STORES the raw "Type: message" for server-side
    # forensics; the READ boundary (GET /api/jobs AND the /job page) must
    # genericise an UNEXPECTED internal type exactly as the sync path does, while
    # an honest QuiverlabError/runner-tag failure still shows its type+message.
    cfg = _cfg(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    # An unexpected internal failure exactly as the worker would store it.
    bad = store.create_job(_gf_body(["hh_cohomology:0..2"]), ip="h")
    store.mark_failed(bad, "KeyError: 'SECRET_INTERNAL'")
    # An honest, client-safe refusal (a QuiverlabError subclass name).
    good = store.create_job(_gf_body(["hh_cohomology:0..2"]), ip="h")
    store.mark_failed(good, "RelationError: relation is not admissible")

    client = TestClient(create_app(cfg))

    # Unexpected internal error: nothing leaks -- on BOTH the API and the page.
    api = client.get(f"/api/jobs/{bad}")
    assert api.status_code == 200
    assert api.json()["error"] == "InternalError: " + GENERIC_ERROR_MESSAGE
    assert "KeyError" not in api.text and "SECRET_INTERNAL" not in api.text
    page = client.get(f"/job/{bad}")
    assert page.status_code == 200
    assert "KeyError" not in page.text and "SECRET_INTERNAL" not in page.text
    assert "InternalError" in page.text

    # Honest QuiverlabError: its type + message DO pass through verbatim.
    assert client.get(f"/api/jobs/{good}").json()["error"] == \
        "RelationError: relation is not admissible"
    assert "RelationError" in client.get(f"/job/{good}").text


def test_client_ip_trusts_last_xff_hop(tmp_path):
    # Behind the single trusted reverse proxy (Caddy APPENDS the real peer as the
    # LAST X-Forwarded-For hop), the stored/rate-limit key is that last hop --
    # never the leftmost, which a client behind an appending proxy can spoof.
    cfg = _cfg(tmp_path, QLWEB_INSTANT_MAX_DEGREE=0)      # force a queued row
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..2"]),
                    headers={"X-Forwarded-For": "9.9.9.9, 10.0.0.5"})
    assert r.status_code == 202
    (stored,) = sqlite3.connect(cfg.db_path).execute("SELECT ip FROM jobs").fetchone()
    assert stored == hash_ip("10.0.0.5", cfg.ip_hash_salt)   # Caddy-appended peer
    assert stored != hash_ip("9.9.9.9", cfg.ip_hash_salt)    # spoofable leftmost ignored
