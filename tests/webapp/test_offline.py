"""Plan 28 -- offline GUI mode tests (fast).

Exercised without a real uvicorn server: ``create_offline_app`` builds the same
FastAPI app ``serve_offline`` would serve, and a ``TestClient`` drives it. Covers
the offline config (big-job tier off), host-resource-derived defaults + display
(``detect_resources`` monkeypatched), the seed-cache first-run copy + a served
cache hit, the banner content, and the estimator memory field surfaced in the API
payload + both language pages.
"""
import importlib.util
import pathlib

import pytest
from fastapi.testclient import TestClient

from webapp.server import offline as offline_mod
from webapp.server.app import create_app
from webapp.server.cache import canonical_key, library_version
from webapp.server.config import Config
from webapp.server.estimator import classify, estimate_bytes, human_bytes
from webapp.server.offline import (
    _banner_lines, build_offline_config, create_offline_app, runtime_caps,
)
from webapp.server.schema import ComputeRequest
from webapp.server.store import JobStore

_SEED_CACHE_PY = (pathlib.Path(__file__).resolve().parents[2]
                  / "container" / "seed_cache.py")

FAKE_RES = {"cores": 4, "mem_bytes": 8 * 1024 ** 3, "gpus": 0, "gpus_used": False}
FAKE_RES_GPU = {"cores": 16, "mem_bytes": 64 * 1024 ** 3, "gpus": 2, "gpus_used": False}

# A tiny, fast request (Cartan only -- no HH blow-up) used for the seed bundle.
LOOP_CARTAN = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False},
}


def _load_seed_cache():
    spec = importlib.util.spec_from_file_location("seed_cache", _SEED_CACHE_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def patched_res(monkeypatch):
    """Deterministic host resources + no ambient seed/data env leaking in."""
    monkeypatch.setattr(offline_mod, "detect_resources", lambda: dict(FAKE_RES))
    monkeypatch.delenv("QUIVERLAB_SEED_CACHE", raising=False)
    monkeypatch.delenv("QUIVERLAB_DATA", raising=False)


def _app(tmp_path, env=None):
    return create_offline_app(data_dir=tmp_path, env=env if env is not None else {})


# --------------------------------------------------------------------------- #
# Offline config + host-resource-derived defaults
# --------------------------------------------------------------------------- #

def test_offline_disables_big_jobs_and_flags_state(patched_res, tmp_path):
    app, cfg, res = _app(tmp_path)
    assert cfg.big_jobs_enabled is False        # no SMTP configured
    assert app.state.offline is True
    assert app.state.resources == FAKE_RES


def test_worker_mem_cap_derives_from_detected_ram(patched_res, tmp_path):
    _app_, cfg, _res = _app(tmp_path)
    assert cfg.job_mem_bytes == FAKE_RES["mem_bytes"] * 4 // 5   # four fifths of RAM
    assert cfg.worker_processes == 1                             # laptop default


def test_explicit_env_overrides_detected_defaults(patched_res, tmp_path):
    _app_, cfg, _res = _app(tmp_path, env={"QLWEB_JOB_MEM_BYTES": "123456",
                                           "QLWEB_WORKER_PROCESSES": "3"})
    assert cfg.job_mem_bytes == 123456 and cfg.worker_processes == 3


def test_build_offline_config_is_isolated_from_environment(tmp_path):
    res = {"cores": 2, "mem_bytes": 2 * 1024 ** 3, "gpus": 0, "gpus_used": False}
    cfg = build_offline_config(tmp_path, res, env={})
    assert cfg.data_dir == tmp_path
    assert cfg.job_mem_bytes == res["mem_bytes"] * 4 // 5


# --------------------------------------------------------------------------- #
# Banner + footer content (detected resources alongside configured caps)
# --------------------------------------------------------------------------- #

def test_banner_shows_open_url_and_caps(patched_res, tmp_path):
    app, cfg, _res = _app(tmp_path)
    text = "\n".join(_banner_lines(8000, cfg, app.state.caps, True))
    assert "open  http://localhost:8000" in text
    assert "4 core(s)" in text and "RAM detected" in text
    assert "worker caps:" in text and "worker(s)" in text
    assert "big-job/email tier: disabled" in text
    assert "GPU" not in text                    # no GPU line when none detected


def test_banner_gpu_note_when_detected(monkeypatch, tmp_path):
    monkeypatch.setattr(offline_mod, "detect_resources", lambda: dict(FAKE_RES_GPU))
    monkeypatch.delenv("QUIVERLAB_SEED_CACHE", raising=False)
    app, cfg, _res = _app(tmp_path)
    text = "\n".join(_banner_lines(8000, cfg, app.state.caps, True))
    assert "2 GPU(s) detected -- not used" in text
    assert "exact CPU computation" in text


def test_footer_shows_resources_bilingual(monkeypatch, tmp_path):
    monkeypatch.setattr(offline_mod, "detect_resources", lambda: dict(FAKE_RES_GPU))
    monkeypatch.delenv("QUIVERLAB_SEED_CACHE", raising=False)
    app, _cfg, _res = _app(tmp_path)
    c = TestClient(app)
    en = c.get("/").text
    assert "Local limits:" in en and "CPU cores" in en and "RAM detected" in en
    assert "2 GPU(s) detected" in en and "not used" in en
    es = c.get("/es").text
    assert "Límites locales:" in es and "RAM detectada" in es
    assert "2 GPU detectada(s)" in es and "sin usar" in es
    # The ES page must not leak the English strings for any of the new keys.
    for eng in ("Local limits", "CPU cores", "RAM detected", "worker memory cap",
                "not used", "worker process(es)"):
        assert eng not in es, f"English leaked into /es footer: {eng!r}"


def test_footer_absent_on_deployed_server(tmp_path):
    # The deployed (non-offline) app never sets app.state.offline -> no footer block.
    c = TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))
    html = c.get("/").text
    assert "Local limits" not in html and "RAM detected" not in html


def test_runtime_caps_human_readable(patched_res, tmp_path):
    _app_, cfg, res = _app(tmp_path)
    caps = runtime_caps(cfg, FAKE_RES)
    assert caps["cores"] == 4
    assert caps["ram_human"] == "8.0 GiB"
    assert caps["worker_mem_human"].endswith("GiB")
    assert caps["worker_mem_bytes"] == cfg.job_mem_bytes


# --------------------------------------------------------------------------- #
# Seed cache: first-run copy honored + a cache hit served
# --------------------------------------------------------------------------- #

def test_seed_first_run_copy_and_cache_hit(monkeypatch, tmp_path):
    monkeypatch.setattr(offline_mod, "detect_resources", lambda: dict(FAKE_RES))
    # Build a seed bundle: DB + sibling artifacts/ (the shape offline copies).
    seed = _load_seed_cache()
    seed_db = tmp_path / "seedsrc" / "seed-cache.db"
    ok, _skipped, total = seed.seed([LOOP_CARTAN], seed_db)
    assert ok == 1 and total == 1
    monkeypatch.setenv("QUIVERLAB_SEED_CACHE", str(seed_db))

    data_dir = tmp_path / "data"           # fresh -> first run
    app, cfg, _res = create_offline_app(data_dir=data_dir, env={})
    assert app.state.seeded is True
    assert cfg.db_path.exists()            # seed DB copied in
    # artifacts copied under the data dir's artifacts/
    assert any(cfg.artifacts_dir.glob("*/result.json"))

    c = TestClient(app)
    r = c.post("/api/compute", json=LOOP_CARTAN)
    assert r.status_code == 200 and r.json()["tier"] == "cached", r.text
    jid = r.json()["job_id"]
    page = c.get(f"/job/{jid}")
    assert page.status_code == 200
    assert "Reproduce on a cluster" in page.text
    assert "quiverlab-hpc run this-file.yaml" in page.text


def test_offline_worker_computes_and_caches_a_queued_job(patched_res, tmp_path):
    # The embedded-worker path under the RAM-derived offline config: a queued job
    # computes in a resource-capped child and becomes replayable via the cache.
    from webapp.worker.worker import worker_tick
    _app_, cfg, _res = _app(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    dump = ComputeRequest.model_validate(LOOP_CARTAN).model_dump(by_alias=True)
    jid = store.create_job(dump, ip="")
    assert worker_tick(store, cfg) is True
    assert store.get_job(jid).status == "done"
    assert store.cache_row(canonical_key(dump, library_version())) is not None


def test_seed_not_recopied_when_db_exists(monkeypatch, tmp_path):
    monkeypatch.setattr(offline_mod, "detect_resources", lambda: dict(FAKE_RES))
    seed = _load_seed_cache()
    seed_db = tmp_path / "seedsrc" / "seed-cache.db"
    seed.seed([LOOP_CARTAN], seed_db)
    monkeypatch.setenv("QUIVERLAB_SEED_CACHE", str(seed_db))
    data_dir = tmp_path / "data"
    first, _c, _r = create_offline_app(data_dir=data_dir, env={})
    assert first.state.seeded is True
    second, _c2, _r2 = create_offline_app(data_dir=data_dir, env={})
    assert second.state.seeded is False        # DB already present -> no re-copy


# --------------------------------------------------------------------------- #
# Estimator memory field: exact ints, in the API payload, on both language pages
# --------------------------------------------------------------------------- #

def test_estimate_bytes_and_human_are_exact_ints():
    assert estimate_bytes(10, 4, "GF") == 10 * 10 * 5 * 8
    assert estimate_bytes(10, 4, "CC") == 10 * 10 * 5 * 64
    assert isinstance(estimate_bytes(10, 4, "GF"), int)
    assert human_bytes(0) == "0 B"
    assert human_bytes(1023) == "1023 B"
    assert human_bytes(1024) == "1.0 KiB"
    assert human_bytes(140000) == "136.7 KiB"


def test_classify_estimate_carries_memory(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = ComputeRequest.model_validate({
        "schema": 1, "algebra": {"kind": "family", "family": "QuantumCI",
                                 "params": {"q": 1}, "field": {"kind": "CC"}},
        "compute": ["hh_cohomology:0..30"], "artifacts": {"pdf": False, "tikz": False}})
    est = classify(300, req, cfg)["estimate"]
    assert isinstance(est["bytes"], int) and est["bytes"] > 0
    assert isinstance(est["mem_human"], str) and est["mem_human"]


def _cc_body(compute):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": "QuantumCI",
                        "params": {"q": 1}, "field": {"kind": "CC"}},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False}}


def test_memory_field_in_big_estimate_payload(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org",
                           # big-job tier enabled -> a real signing secret, else
                           # create_app's production-secret guard (correction #4) refuses.
                           "QLWEB_TOKEN_SECRET": "test-secret-not-the-default"})
    c = TestClient(create_app(cfg))
    r = c.post("/api/compute", json=_cc_body(["hh_cohomology:0..30"]))
    assert r.status_code == 202 and r.json()["tier"] == "big"
    est = r.json()["estimate"]
    assert isinstance(est["bytes"], int) and est["bytes"] > 0
    assert est["mem_human"]


def test_big_warn_template_shows_memory_both_languages(tmp_path):
    c = TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))
    en = c.get("/").text
    assert "~{mem} memory)" in en                 # EN big.warn carries the mem slot
    es = c.get("/es").text
    assert "~{mem} de memoria)" in es              # ES carries it, translated
    assert " memory)" not in es                    # and does not leak the English word


def test_estimate_field_in_reject_payload_when_big_disabled(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})     # no SMTP -> disabled
    c = TestClient(create_app(cfg))
    r = c.post("/api/compute", json=_cc_body(["hh_cohomology:0..30"]))
    assert r.status_code == 422
    assert isinstance(r.json()["estimate"]["bytes"], int)
