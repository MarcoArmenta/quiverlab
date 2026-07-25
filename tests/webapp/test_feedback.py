"""Task 12 -- feedback: bilingual page, JSON submit API, token-gated admin view.

Adapted from the brief against the REAL in-repo interfaces:
  * the i18n `fb.*` keys and `nav.feedback` already ship in BOTH en.json/es.json
    (Task 11), so the page tests key off the real ES strings ("Comentarios",
    "Sugerir bibliografia");
  * IP hashing reuses the app's `client_ip` + `security.hash_ip` (never a raw
    address in the store);
  * `job_ref` is ULID-validated when present (the plan's global constraint the
    brief's inline code omitted) -- see the invalid/valid job_ref tests;
  * the admin token arrives in the ``X-Admin-Token`` HEADER (never a query
    string -- uvicorn logs query strings, so a ``?token=`` would leak the secret
    into access logs); a wrong token -> 401.

Beyond the brief this file pins the additional adjudicated coverage the task
requires: message length bounds at 9/10/4001 chars, the 6th same-day submission
over the DEFAULT cap of 5, a structural assertion that the admin compare is
constant-time (`hmac.compare_digest`), and that the admin table HTML-escapes
untrusted feedback text.
"""
import json

from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.security import valid_ulid
from webapp.server.store import JobStore

_ULID = "01AN4Z07BY79KA1307SR9X4MV3"        # a well-formed Crockford-base32 ULID


def _client(tmp_path, **env):
    e = {"QLWEB_DATA_DIR": str(tmp_path)}
    e.update(env)
    return TestClient(create_app(Config.from_env(e)))


def _body(**over):
    b = {"category": "problem", "message": "HH^3 dims look wrong for A5 over GF(2).",
         "contact": "", "job_ref": "", "website": ""}
    b.update(over)
    return b


# --------------------------------------------------------------------------- #
# Submit API: happy path, honeypot, validation, rate limit.
# --------------------------------------------------------------------------- #
def test_feedback_roundtrip(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/feedback", json=_body())
    assert r.status_code == 201, r.text
    ref = r.json()["reference"]
    rows = JobStore(cfg.db_path).list_feedback()
    assert rows and rows[0]["id"] == ref
    assert rows[0]["category"] == "problem"


def test_feedback_honeypot_dropped_silently(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/feedback", json=_body(website="http://spam"))
    assert r.status_code == 201                          # looks fine to the bot
    # Body shape is INDISTINGUISHABLE from a real success: a 26-char Crockford
    # ULID reference. A parsing bot cannot detect the drop.
    ref = r.json()["reference"]
    assert len(ref) == 26 and valid_ulid(ref)
    assert JobStore(cfg.db_path).list_feedback() == []   # but nothing stored


def test_feedback_honeypot_whitespace_dropped(tmp_path):
    # A whitespace-only honeypot value is still a bot -- any NON-EMPTY raw value
    # triggers the silent drop (no .strip() bypass).
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/feedback", json=_body(website="   "))
    assert r.status_code == 201                           # dropped silently
    assert JobStore(cfg.db_path).list_feedback() == []    # nothing stored


def test_feedback_too_short_is_422(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(message="hi"))
    assert r.status_code == 422


def test_feedback_message_9_chars_is_422(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(message="a" * 9))
    assert r.status_code == 422                           # min length is 10


def test_feedback_message_10_chars_ok(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(message="a" * 10))
    assert r.status_code == 201, r.text                  # exactly the lower bound


def test_feedback_message_4000_chars_ok(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(message="a" * 4000))
    assert r.status_code == 201, r.text                  # exactly the upper bound


def test_feedback_message_4001_chars_is_422(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(message="a" * 4001))
    assert r.status_code == 422                           # max length is 4000


def test_feedback_rate_limit_429(tmp_path):
    client = _client(tmp_path, QLWEB_FEEDBACK_DAILY_MAX="2")
    for _ in range(2):
        assert client.post("/api/feedback", json=_body()).status_code == 201
    assert client.post("/api/feedback", json=_body()).status_code == 429


def test_feedback_sixth_same_day_is_429(tmp_path):
    client = _client(tmp_path)                            # DEFAULT cap == 5
    for _ in range(5):
        assert client.post("/api/feedback", json=_body()).status_code == 201
    assert client.post("/api/feedback", json=_body()).status_code == 429


def test_feedback_invalid_job_ref_is_422(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(job_ref="not-a-ulid"))
    assert r.status_code == 422                           # ULID-validated when present


def test_feedback_valid_job_ref_is_stored(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/feedback", json=_body(job_ref=_ULID))
    assert r.status_code == 201, r.text
    assert JobStore(cfg.db_path).list_feedback()[0]["job_ref"] == _ULID


# --------------------------------------------------------------------------- #
# Pages (bilingual).
# --------------------------------------------------------------------------- #
def test_feedback_page_prefills_job_ref(tmp_path):
    r = _client(tmp_path).get("/feedback?job=" + _ULID)
    assert r.status_code == 200
    assert _ULID in r.text
    assert "Comentarios" not in r.text                   # English page


def test_feedback_page_es(tmp_path):
    r = _client(tmp_path).get("/es/feedback")
    assert r.status_code == 200
    assert "Comentarios" in r.text


def test_feedback_page_has_literature_option(tmp_path):
    r = _client(tmp_path).get("/es/feedback")
    assert "Sugerir bibliografía" in r.text          # literature category label (ES)


# --------------------------------------------------------------------------- #
# Admin view (only when a token is configured; constant-time compare).
# --------------------------------------------------------------------------- #
def test_admin_route_absent_without_token(tmp_path):
    r = _client(tmp_path).get("/admin/feedback", headers={"X-Admin-Token": "whatever"})
    assert r.status_code == 404


def test_admin_token_gate(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_ADMIN_TOKEN": "sekret"})
    client = TestClient(create_app(cfg))
    ref = client.post("/api/feedback", json=_body()).json()["reference"]
    wrong = client.get("/admin/feedback", headers={"X-Admin-Token": "wrong"})
    assert wrong.status_code == 401                       # 401 per plan (not 403)
    assert wrong.text == "unauthorized"
    # The token lives in the X-Admin-Token HEADER, never a query string (a
    # ?token=... would land the secret in uvicorn's access log): a query-param
    # token does NOT authenticate.
    assert client.get("/admin/feedback?token=sekret").status_code == 401
    ok = client.get("/admin/feedback", headers={"X-Admin-Token": "sekret"})
    assert ok.status_code == 200
    assert "HH^3 dims look wrong" in ok.text
    assert ref in ok.text                                 # the id column renders the reference


def test_admin_non_ascii_token_is_401_not_500(tmp_path):
    # A non-ASCII token header must NOT crash hmac.compare_digest (TypeError on
    # str) -- both sides are encoded to bytes, so it refuses cleanly. (Passed as
    # UTF-8 bytes because httpx rejects a non-ASCII str header value client-side;
    # the server still decodes it to a non-ASCII str, exercising the guard.)
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_ADMIN_TOKEN": "sekret"})
    client = TestClient(create_app(cfg))
    r = client.get("/admin/feedback", headers={"X-Admin-Token": "café".encode("utf-8")})
    assert r.status_code == 401                           # clean refusal, never a 500


def test_admin_escapes_untrusted_text(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_ADMIN_TOKEN": "t"})
    client = TestClient(create_app(cfg))
    client.post("/api/feedback", json=_body(message="<script>alert(1)</script> bug report"))
    r = client.get("/admin/feedback", headers={"X-Admin-Token": "t"})
    assert r.status_code == 200
    assert "<script>alert(1)</script>" not in r.text     # escaped, not raw
    assert "&lt;script&gt;" in r.text


def test_admin_compare_is_constant_time(tmp_path):
    # Structural, not a file-wide grep: inspect the ACTUAL admin handler the app
    # registered and assert its body uses hmac.compare_digest and never a plain
    # `== cfg.admin_token`.
    import inspect

    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_ADMIN_TOKEN": "t"})
    app = create_app(cfg)
    handler = next(r.endpoint for r in app.routes
                   if getattr(r, "path", None) == "/admin/feedback")
    src = inspect.getsource(handler)
    assert "hmac.compare_digest" in src                  # constant-time compare in the handler
    assert " == cfg.admin_token" not in src              # never a plain == token compare


# --------------------------------------------------------------------------- #
# Literature category: structured extra JSON.
# --------------------------------------------------------------------------- #
def test_literature_submission_roundtrip(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/feedback", json=_body(
        category="literature",
        message="This paper's resolution should be cited.",
        reference="arXiv:1406.2300",
        why_relevant="It is the Chouhy-Solotar resolution the engine implements."))
    assert r.status_code == 201, r.text
    row = JobStore(cfg.db_path).list_feedback()[0]
    assert row["category"] == "literature"
    extra = json.loads(row["extra"])
    assert extra["reference"] == "arXiv:1406.2300"
    assert "Chouhy-Solotar" in extra["why_relevant"]


def test_literature_requires_structured_fields(tmp_path):
    r = _client(tmp_path).post("/api/feedback", json=_body(
        category="literature", message="Please cite something relevant."))
    assert r.status_code == 422           # reference + why_relevant missing


def test_non_literature_stores_no_extra(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    client.post("/api/feedback", json=_body())          # category problem
    assert JobStore(cfg.db_path).list_feedback()[0]["extra"] is None
