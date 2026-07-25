"""Task 13 -- big-job email magic-link tier (spec §17).

Adapted from the task brief: the reference family is ``QuantumCI(q=1)`` (the
brief's ``truncated_polynomial`` is not in ``quiverlab.families()`` and would
never build -- the same 2026-07-24 amendment already applied to the worker
tests). ``QuantumCI(q=1)`` has dim 4; over CC with ``hh_cohomology:0..30`` the
degree (30, in (queued_max=20, big_max=40]) lands it in the big tier.

No test performs a real SMTP send: every path uses the injected ``FakeMailer``.
"""
import logging
import re
import smtplib
import sqlite3
import time

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app as create_app_
from webapp.server.bigjobs import _email_hash, _spec_hash, sign_token, verify_token
from webapp.server.config import Config
from webapp.server.security import hash_ip
from webapp.server.store import JobStore


class FakeMailer:
    """Captures (to, subject, body) instead of sending -- the network is never
    touched, and the captured tuples let a test read back the magic link."""

    def __init__(self):
        self.sent = []                       # list of (to, subject, body)

    def __call__(self, to, subject, body):
        self.sent.append((to, subject, body))


def _big_cfg(tmp_path, **over):
    e = {"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_SMTP_HOST": "relay",
         "QLWEB_SMTP_FROM": "quiverlab@example.org",
         "QLWEB_PUBLIC_BASE_URL": "https://ql.example"}
    e.update(over)
    return Config.from_env(e)


def _big_body(**over):
    b = {"schema": 1,
         "algebra": {"kind": "family", "family": "QuantumCI",
                     "params": {"q": 1}, "field": {"kind": "CC"}},
         "compute": ["hh_cohomology:0..30"], "artifacts": {"pdf": False, "tikz": False},
         "email": "user@example.org", "lang": "en"}
    b.update(over)
    return b


def _token_from(mailer):
    body = mailer.sent[-1][2]
    return re.search(r"/verify/([^\s]+)", body).group(1)


def _all_job_ids(store):
    conn = sqlite3.connect(store.db_path)
    try:
        return [r[0] for r in conn.execute("SELECT id FROM jobs").fetchall()]
    finally:
        conn.close()


def _pending_big_count(store):
    conn = sqlite3.connect(store.db_path)
    try:
        return conn.execute("SELECT COUNT(*) FROM pending_big").fetchone()[0]
    finally:
        conn.close()


# --------------------------------------------------------------------------- #
# Token: sign / verify / expiry / tamper
# --------------------------------------------------------------------------- #

def test_token_roundtrip_expiry_and_tamper():
    secret = "s3cret"
    tok = sign_token(secret, {"pid": "P", "sh": "H", "exp": 1000})
    assert verify_token(secret, tok, now=999)["pid"] == "P"
    assert verify_token(secret, tok, now=1001) is None            # expired
    assert verify_token(secret, tok + "x", now=999) is None       # tampered sig
    assert verify_token("other", tok, now=999) is None            # wrong secret
    assert verify_token(secret, "no-dot-here", now=999) is None   # malformed


# --------------------------------------------------------------------------- #
# SMTP-off => tier disabled end-to-end
# --------------------------------------------------------------------------- #

def test_big_disabled_returns_503(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})       # no SMTP
    assert cfg.big_jobs_enabled is False
    r = TestClient(create_app_(cfg)).post("/api/jobs/big", json=_big_body())
    assert r.status_code == 503


# --------------------------------------------------------------------------- #
# Request-big -> pending row + exactly one verification mail; link is single-use
# --------------------------------------------------------------------------- #

def test_submit_and_verify_end_to_end(tmp_path):
    cfg = _big_cfg(tmp_path)
    mail = FakeMailer()
    client = TestClient(create_app_(cfg, mailer=mail))
    r = client.post("/api/jobs/big", json=_big_body())
    assert r.status_code == 202 and r.json()["status"] == "sent"
    # Exactly one verification mail, to the requester, never logged (see the
    # caplog test); the store holds the pending row but no job yet.
    assert len(mail.sent) == 1 and mail.sent[0][0] == "user@example.org"
    store = JobStore(cfg.db_path)
    assert _all_job_ids(store) == []                              # not queued until verified

    token = _token_from(mail)
    v = client.get("/verify/" + token)
    assert v.status_code == 200
    jobs = [store.get_job(j) for j in _all_job_ids(store)]
    big = [j for j in jobs if j.tier == "big"]
    assert len(big) == 1
    assert big[0].wall_seconds == cfg.big_job_wall_seconds
    assert big[0].mem_bytes == cfg.big_job_mem_bytes
    assert big[0].email == "user@example.org"                    # kept until completion
    # Single use: the same link cannot queue a second job.
    v2 = client.get("/verify/" + token)
    assert "already used" in v2.text.lower() or v2.status_code == 410
    assert len([j for j in _all_job_ids(store)]) == 1            # still exactly one job


# --------------------------------------------------------------------------- #
# Spec binding: a tampered spec-hash never queues a job, and the link is spent
# --------------------------------------------------------------------------- #

def test_spec_binding_rejects_mismatched_hash(tmp_path):
    cfg = _big_cfg(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    pid = store.create_pending_big({"schema": 1, "compute": ["cartan"]},
                                   "u@e.org", "eh", lang="en")
    bad = sign_token(cfg.token_secret, {"pid": pid, "sh": "not-the-real-hash",
                                        "exp": int(time.time()) + 60})
    v = TestClient(create_app_(cfg, mailer=FakeMailer())).get("/verify/" + bad)
    assert "already used" in v.text.lower() or v.status_code == 410
    assert store.consume_pending_big(pid) is None                # consumed/rejected safely
    assert _all_job_ids(store) == []                             # no job created


def test_expired_token_is_a_clean_rejection(tmp_path):
    cfg = _big_cfg(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    pid = store.create_pending_big({"schema": 1, "compute": ["cartan"]},
                                   "u@e.org", "eh", lang="en")
    sh = _spec_hash({"schema": 1, "compute": ["cartan"]})
    expired = sign_token(cfg.token_secret,
                         {"pid": pid, "sh": sh, "exp": int(time.time()) - 1})
    v = TestClient(create_app_(cfg, mailer=FakeMailer())).get("/verify/" + expired)
    assert "already used" in v.text.lower() or v.status_code == 410
    # An expired link never consumes the pending row (nothing was built): the row
    # is still present and would only be reaped by the sweeper.
    assert store.consume_pending_big(pid) is not None
    assert _all_job_ids(store) == []


# --------------------------------------------------------------------------- #
# Per-email caps => 429
# --------------------------------------------------------------------------- #

def test_per_email_running_cap_429(tmp_path):
    cfg = _big_cfg(tmp_path, QLWEB_PER_EMAIL_RUNNING_MAX="1")
    store = JobStore(cfg.db_path)
    store.init_schema()
    eh = hash_ip("user@example.org", cfg.token_secret)
    store.create_job({}, ip="i", tier="big", email="user@example.org", email_hash=eh)
    r = TestClient(create_app_(cfg, mailer=FakeMailer())).post(
        "/api/jobs/big", json=_big_body())
    assert r.status_code == 429


def test_per_email_weekly_cap_429(tmp_path):
    cfg = _big_cfg(tmp_path, QLWEB_PER_EMAIL_WEEKLY_MAX="2",
                   QLWEB_PER_EMAIL_RUNNING_MAX="100")
    store = JobStore(cfg.db_path)
    store.init_schema()
    eh = hash_ip("user@example.org", cfg.token_secret)
    # Two completed big jobs this week already exhaust the weekly cap.
    for _ in range(2):
        jid = store.create_job({}, ip="i", tier="big",
                               email="user@example.org", email_hash=eh)
        store.mark_done(jid, str(cfg.artifacts_dir / jid))
    r = TestClient(create_app_(cfg, mailer=FakeMailer())).post(
        "/api/jobs/big", json=_big_body())
    assert r.status_code == 429


# --------------------------------------------------------------------------- #
# Queue cap: a transiently full big queue must NOT burn a still-valid link
# --------------------------------------------------------------------------- #

def test_big_queue_cap_does_not_burn_link(tmp_path):
    mail = FakeMailer()
    cfg0 = _big_cfg(tmp_path, QLWEB_BIG_QUEUE_MAX="0")
    client0 = TestClient(create_app_(cfg0, mailer=mail))
    assert client0.post("/api/jobs/big", json=_big_body()).status_code == 202
    token = _token_from(mail)
    # Queue is full (cap 0): the link is rejected but NOT consumed.
    v = client0.get("/verify/" + token)
    assert v.status_code == 410
    # Same DB, capacity restored: the untouched link now verifies and queues.
    cfg1 = _big_cfg(tmp_path, QLWEB_BIG_QUEUE_MAX="20")
    client1 = TestClient(create_app_(cfg1, mailer=mail))
    v2 = client1.get("/verify/" + token)
    assert v2.status_code == 200
    store = JobStore(cfg1.db_path)
    assert len([j for j in (store.get_job(i) for i in _all_job_ids(store))
                if j.tier == "big"]) == 1


# --------------------------------------------------------------------------- #
# Bilingual verify page
# --------------------------------------------------------------------------- #

def test_es_verify_page(tmp_path):
    cfg = _big_cfg(tmp_path)
    tok = sign_token(cfg.token_secret, {"pid": "nope", "sh": "x",
                                        "exp": int(time.time()) + 60})
    r = TestClient(create_app_(cfg, mailer=FakeMailer())).get("/es/verify/" + tok)
    assert "Trabajo grande" in r.text


# --------------------------------------------------------------------------- #
# Completion mail sent on done; plaintext email cleared, hash kept
# --------------------------------------------------------------------------- #

def test_email_deleted_after_completion(tmp_path):
    from webapp.worker.worker import worker_tick
    cfg = _big_cfg(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    eh = hash_ip("user@example.org", cfg.token_secret)
    jid = store.create_job(
        {"schema": 1, "algebra": {"kind": "family", "family": "QuantumCI",
         "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": ["hh_cohomology:0..3"], "artifacts": {"pdf": False, "tikz": False}},
        ip="i", tier="big", email="user@example.org", email_hash=eh,
        wall_seconds=cfg.big_job_wall_seconds, mem_bytes=cfg.big_job_mem_bytes, lang="en")
    mail = FakeMailer()
    worker_tick(store, cfg, mailer=mail)
    done = store.get_job(jid)
    assert done.status == "done", done.error
    assert done.email is None                    # plaintext deleted after notice
    assert done.email_hash == eh                 # hash kept for weekly rate-limiting
    assert mail.sent and "user@example.org" == mail.sent[-1][0]


# --------------------------------------------------------------------------- #
# Email is never in the admin feedback view, and never in any log line
# --------------------------------------------------------------------------- #

def test_email_never_in_admin_feedback(tmp_path):
    cfg = _big_cfg(tmp_path, QLWEB_ADMIN_TOKEN="tok")
    mail = FakeMailer()
    client = TestClient(create_app_(cfg, mailer=mail))
    client.post("/api/jobs/big", json=_big_body())               # pending big w/ email
    admin = client.get("/admin/feedback", headers={"X-Admin-Token": "tok"})
    assert admin.status_code == 200
    assert "user@example.org" not in admin.text                  # never shows job emails


def test_email_never_appears_in_logs(tmp_path, caplog):
    cfg = _big_cfg(tmp_path)
    mail = FakeMailer()
    client = TestClient(create_app_(cfg, mailer=mail))
    with caplog.at_level(logging.DEBUG):
        r = client.post("/api/jobs/big", json=_big_body())
        assert r.status_code == 202
        token = _token_from(mail)
        client.get("/verify/" + token)
        # Exercise the completion-mail path too -- the other place email flows.
        store = JobStore(cfg.db_path)
        job = next(j for j in (store.get_job(i) for i in _all_job_ids(store))
                   if j.tier == "big")
        from webapp.server.mail import notify_completion
        notify_completion(cfg, job, "done", mailer=mail)
    assert "user@example.org" not in caplog.text


# --------------------------------------------------------------------------- #
# Mail failures never reach the logs -- including the exception paths
# --------------------------------------------------------------------------- #

def test_submit_mail_failure_is_502_and_never_logs_address(tmp_path, caplog):
    cfg = _big_cfg(tmp_path)
    addr = "leak-submit@example.org"

    def boom(to, subject, body):
        # A recipient refusal whose str() carries the address -- a naive log would
        # leak it.
        raise smtplib.SMTPRecipientsRefused({addr: (550, b"x")})

    client = TestClient(create_app_(cfg, mailer=boom), raise_server_exceptions=False)
    with caplog.at_level(logging.DEBUG):
        r = client.post("/api/jobs/big", json=_big_body(email=addr))
    assert r.status_code == 502
    assert addr not in caplog.text                       # never logged, even on failure
    # No pending row leaked in a half-state, and no job was created.
    store = JobStore(cfg.db_path)
    assert _pending_big_count(store) == 0
    assert _all_job_ids(store) == []


def test_completion_mail_failure_keeps_job_done_and_never_logs_address(tmp_path, caplog):
    from webapp.worker.worker import worker_tick
    cfg = _big_cfg(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    addr = "leak-done@example.org"
    eh = hash_ip(addr, cfg.token_secret)
    jid = store.create_job(
        {"schema": 1, "algebra": {"kind": "family", "family": "QuantumCI",
         "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": ["hh_cohomology:0..3"], "artifacts": {"pdf": False, "tikz": False}},
        ip="i", tier="big", email=addr, email_hash=eh,
        wall_seconds=cfg.big_job_wall_seconds, mem_bytes=cfg.big_job_mem_bytes, lang="en")

    def boom(to, subject, body):
        raise RuntimeError("smtp down for " + to)        # str() carries the address

    with caplog.at_level(logging.DEBUG):
        worker_tick(store, cfg, mailer=boom)
    done = store.get_job(jid)
    assert done.status == "done", done.error             # the result is never lost
    assert done.email is None                            # email still cleared on failure
    assert done.email_hash == eh                         # hash kept for rate-limiting
    assert addr not in caplog.text                       # never logged, even on failure


# --------------------------------------------------------------------------- #
# Strict single-address shape check
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("bad", ["a@b.org\nBcc: x@y.org", "a b@c.org", "a@b"])
def test_email_shape_rejects_bad_addresses(tmp_path, bad):
    cfg = _big_cfg(tmp_path)
    client = TestClient(create_app_(cfg, mailer=FakeMailer()))
    r = client.post("/api/jobs/big", json=_big_body(email=bad))
    assert r.status_code == 422


# --------------------------------------------------------------------------- #
# Plus-alias cap bypass: +tag variants share one rate-limit bucket
# --------------------------------------------------------------------------- #

def test_plus_alias_shares_rate_bucket(tmp_path):
    cfg = _big_cfg(tmp_path, QLWEB_PER_EMAIL_RUNNING_MAX="1")
    store = JobStore(cfg.db_path)
    store.init_schema()
    # A running big job under user+1@x.org occupies the shared bucket ...
    eh = _email_hash("user+1@x.org", cfg)
    store.create_job({}, ip="i", tier="big", email="user+1@x.org", email_hash=eh)
    # ... so user+2@x.org hits the same running cap instead of slipping past.
    r = TestClient(create_app_(cfg, mailer=FakeMailer())).post(
        "/api/jobs/big", json=_big_body(email="user+2@x.org"))
    assert r.status_code == 429


# --------------------------------------------------------------------------- #
# Completion-mail permalink is localized to the job's language
# --------------------------------------------------------------------------- #

def test_completion_mail_localizes_url_for_es(tmp_path):
    from webapp.server.mail import notify_completion
    cfg = _big_cfg(tmp_path)

    class _Job:
        id = "JOB123"
        email = "user@example.org"
        lang = "es"

    mail = FakeMailer()
    notify_completion(cfg, _Job(), "done", mailer=mail)
    assert "/es/job/JOB123" in mail.sent[-1][2]
