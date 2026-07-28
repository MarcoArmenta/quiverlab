import time
from pathlib import Path

import yaml

from webapp.server.config import Config
from webapp.server.store import JobStore
from webapp.worker.run_loop import start_workers

DEPLOY = Path("webapp/deploy")


def test_compose_has_services_and_data_mount():
    compose = yaml.safe_load((DEPLOY / "docker-compose.yml").read_text(encoding="utf-8"))
    assert set(compose["services"]) >= {"app", "worker", "caddy"}
    app = compose["services"]["app"]
    assert app["environment"]["NUMBA_NUM_THREADS"] == 2 or app["environment"]["NUMBA_NUM_THREADS"] == "2"
    assert any("/data" in v for v in app["volumes"])


def test_compose_requires_prod_secrets_and_recoverable_shutdown():
    text = (DEPLOY / "docker-compose.yml").read_text(encoding="utf-8")
    compose = yaml.safe_load(text)
    # Both prod secrets are wired with the required-var syntax on BOTH tiers, so
    # `compose up` refuses to start when they are unset (no dev-default ships).
    for svc in ("app", "worker"):
        env = compose["services"][svc]["environment"]
        assert "${QLWEB_IP_HASH_SALT:?" in env["QLWEB_IP_HASH_SALT"]
        assert "${QLWEB_TOKEN_SECRET:?" in env["QLWEB_TOKEN_SECRET"]
    # A template exists and lists exactly the two required secrets (empty values).
    example = (DEPLOY / ".env.example").read_text(encoding="utf-8")
    assert "QLWEB_IP_HASH_SALT=" in example and "QLWEB_TOKEN_SECRET=" in example
    assert "openssl rand -hex 32" in example
    # Worker gets a stop grace so graceful drain has room before Docker SIGKILLs.
    assert "stop_grace_period" in compose["services"]["worker"]
    # App healthcheck + caddy waits for it healthy; Caddyfile mounted read-only.
    assert "healthcheck" in compose["services"]["app"]
    assert compose["services"]["caddy"]["depends_on"]["app"]["condition"] == "service_healthy"
    assert any(v.endswith(":ro") and "Caddyfile" in v
               for v in compose["services"]["caddy"]["volumes"])


def test_caddyfile_has_tls_site():
    text = (DEPLOY / "Caddyfile").read_text(encoding="utf-8")
    assert "reverse_proxy" in text
    assert "app:8000" in text


def test_provisioning_mentions_arbutus_persistent_instance():
    text = (DEPLOY / "PROVISIONING.md").read_text(encoding="utf-8").lower()
    assert "persistent" in text and "floating ip" in text and "volume" in text


def test_provisioning_encodes_security_baseline():
    text = (DEPLOY / "PROVISIONING.md").read_text(encoding="utf-8")
    low = text.lower()
    assert "never open ssh to `0.0.0.0/0`" in low     # SSH-CIDR warning (§10)
    assert "fail2ban" in low                          # host hardening step
    assert "unattended-upgrades" in low               # weekly security updates
    assert "ubuntu 24.04 lts" in low                  # current LTS image


def test_two_workers_drain_jobs(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_WORKER_PROCESSES": "2"})
    store = JobStore(cfg.db_path)
    store.init_schema()
    # A tiny, fast GF(2) job that two loops drain to `done`. Mirrors the real
    # family + compute of tests/webapp/test_worker.py (`QuantumCI(q=1)`); the
    # library exposes no `truncated_polynomial` family (see GET /api/catalog).
    spec = {"schema": 1,
            "algebra": {"kind": "family", "family": "QuantumCI",
                        "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": ["hh_cohomology:0..3"], "artifacts": {"pdf": False, "tikz": False}}
    jids = [store.create_job(spec, ip="h") for _ in range(2)]
    procs, _stop = start_workers(cfg, count=2)
    try:
        deadline = time.monotonic() + 180
        while time.monotonic() < deadline:
            if all(store.get_job(j).status in ("done", "failed") for j in jids):
                break
            time.sleep(1)
        for j in jids:
            assert store.get_job(j).status == "done", store.get_job(j).error
    finally:
        for p in procs:
            p.terminate()
        for p in procs:
            p.join()


def test_workers_stop_gracefully_on_flag(tmp_path):
    # With no pending jobs, each loop idles checking the shared stop event. Setting
    # it (the SIGTERM path in main(), without terminate) makes the loops exit on
    # their own with exit code 0 -- the true graceful drain.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path), "QLWEB_WORKER_PROCESSES": "2"})
    JobStore(cfg.db_path).init_schema()
    procs, stop = start_workers(cfg, count=2)
    try:
        stop.set()
        for p in procs:
            p.join(15)
        assert all(not p.is_alive() for p in procs)
        assert all(p.exitcode == 0 for p in procs)
    finally:
        for p in procs:
            if p.is_alive():
                p.terminate()
        for p in procs:
            p.join()
