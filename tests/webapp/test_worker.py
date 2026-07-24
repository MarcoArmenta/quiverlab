"""Task 8 -- worker process tests.

Adapted from the task brief per the 2026-07-24 amendment: the reference family
is ``QuantumCI(q=1, field=GF(2))`` (the brief's ``truncated_polynomial`` is not
in ``families()``, so it would never build). The worker runs each job in a
``spawn`` child so the rlimit/thread caps set in the child never touch the parent
test process, and the memory-cap probe is Linux-only (macOS does not enforce a
hard RLIMIT_AS).
"""
import multiprocessing as mp
import sys

import pytest

from webapp.server.config import Config
from webapp.server.store import JobStore
from webapp.worker.worker import worker_tick
from webapp.worker.sweeper import sweep_once


_ALG = {"kind": "family", "family": "QuantumCI",
        "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}}


def _setup(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    store = JobStore(cfg.db_path)
    store.init_schema()
    return cfg, store


def test_worker_runs_a_small_job(tmp_path):
    cfg, store = _setup(tmp_path)
    spec = {"schema": 1, "algebra": _ALG,
            "compute": ["hh_cohomology:0..3"],
            "artifacts": {"pdf": False, "tikz": False}}
    jid = store.create_job(spec, ip="1.2.3.4")
    did = worker_tick(store, cfg)
    assert did is True
    job = store.get_job(jid)
    assert job.status == "done", job.error
    assert (tmp_path / "artifacts" / jid / "result.json").exists()


def test_worker_streams_progress(tmp_path):
    cfg, store = _setup(tmp_path)
    spec = {"schema": 1, "algebra": _ALG,
            "compute": ["dimension", "hh_cohomology:0..3", "cartan"],
            "artifacts": {"pdf": False, "tikz": False}}
    jid = store.create_job(spec, ip="1.2.3.4")
    worker_tick(store, cfg)
    job = store.get_job(jid)
    assert job.status == "done", job.error
    # The child streamed per-item checkpoints the parent drained into the store.
    assert "of" in job.progress and job.progress["of"] == 3


def test_worker_records_failure(tmp_path):
    cfg, store = _setup(tmp_path)
    spec = {"schema": 1,
            "algebra": {"kind": "family", "family": "does_not_exist",
                        "params": {}, "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    jid = store.create_job(spec, ip="1.2.3.4")
    worker_tick(store, cfg)
    job = store.get_job(jid)
    assert job.status == "failed"
    assert job.error


def test_sweeper_removes_old(tmp_path):
    cfg, store = _setup(tmp_path)
    jid = store.create_job({}, ip="1.1.1.1")
    store.claim_next()
    store.mark_done(jid, artifact_dir=str(tmp_path / "artifacts" / jid))
    (tmp_path / "artifacts" / jid).mkdir(parents=True, exist_ok=True)
    removed = sweep_once(store, cfg, now_iso="2099-01-01T00:00:00Z")
    assert jid in removed
    assert store.get_job(jid) is None


def _probe_caps(wall, mem, q):
    import resource

    from webapp.worker.worker import _apply_caps
    _apply_caps(wall, mem)
    q.put(resource.getrlimit(resource.RLIMIT_AS)[0])


@pytest.mark.skipif(sys.platform != "linux",
                    reason="RLIMIT_AS is only reliably enforced on Linux")
def test_memory_cap_applied_in_child(tmp_path):
    # Deterministic: assert the child sets RLIMIT_AS to the given per-job cap.
    # (Forcing a real OOM is environment-dependent and flaky; this proves the
    #  guard is wired without endangering the test process.)
    cap = 512 * 1024 ** 2
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=_probe_caps, args=(600, cap, q))
    p.start()
    p.join(10)
    assert q.get() == cap
