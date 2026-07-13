import os
import numpy as np

from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.scan3 import quantum_ci
from quiverlab.engine.deepen import deepen
from quiverlab.engine.resolutions_minimal import minimal_homology_dims


def test_deepen_matches_minimal_homology_dims(tmp_path):
    """deepen's HH_* equals the batch minimal_homology_dims up to the same degree."""
    A = truncated_polynomial(3)               # k[x]/x^3
    ck = str(tmp_path / "ck")
    out = deepen(A, ck, prime=32003, max_degree=5)
    ref = minimal_homology_dims(A, 5, primes=(32003,))[32003]
    assert out["HH"] == ref[:len(out["HH"])]
    assert out["HH"][0] == ref[0]


def test_deepen_resume(tmp_path):
    """A stop-early run then a continue run == one full run (resume correctness)."""
    A = truncated_polynomial(3)
    ck = str(tmp_path / "ck")
    deepen(A, ck, prime=32003, max_degree=2)                # stop early
    out = deepen(A, ck, prime=32003, max_degree=6)          # resume -> continue
    ref = minimal_homology_dims(A, 6, primes=(32003,))[32003]
    assert out["HH"] == ref[:len(out["HH"])]
    assert out["max_degree_reached"] >= 2                   # actually went past the cap


def test_deepen_memory_wall(tmp_path):
    """A tiny transient budget stops with stop_reason='memory' and exact partial HH."""
    A = quantum_ci(2)                                      # dim 4, non-trivial resolution
    ck = str(tmp_path / "ck")
    out = deepen(A, ck, prime=32003, max_transient_bytes=1)  # 1 byte -> trips at degree 1
    assert out["stop_reason"] == "memory"
    assert "wall_radK_bytes" in out and out["wall_radK_bytes"] > 1
    ref = minimal_homology_dims(A, 6, primes=(32003,), max_transient_bytes=1)[32003]
    assert out["HH"] == ref                                 # both exact to the same depth


def test_deepen_checkpoint_reload(tmp_path):
    """After a degree, _load_ckpt returns a payload at the expected degree."""
    from quiverlab.engine.deepen import _load_ckpt
    A = truncated_polynomial(3)
    ck = str(tmp_path / "ck")
    deepen(A, ck, prime=32003, max_degree=3)
    payload = _load_ckpt(ck)
    assert payload is not None and payload["n"] >= 2
    assert "cur" in payload and "rks" in payload


def test_save_ckpt_no_shared_tmp_race(tmp_path, monkeypatch):
    """Concurrent same-degree checkpointing must not crash on a shared temp file -- the
    os.replace FileNotFoundError that killed the doubly-launched job 64134356."""
    import os as _os
    import quiverlab.engine.deepen as D
    ck = str(tmp_path / "ck")
    payload = {"n": 7, "cur": None, "cur_r": 1, "rks": {0: 1},
               "last_gens": None, "HH": [1, 2], "per_degree": []}
    real_replace = _os.replace
    fired = {"v": False}

    def replace_hook(src, dst):
        # mid-save of "process A", let a second "process B" fully checkpoint the SAME
        # degree to the SAME dir, then let A finish its os.replace.
        if not fired["v"]:
            fired["v"] = True
            monkeypatch.setattr(_os, "getpid", lambda: 99999)   # become process B
            D._save_ckpt(ck, dict(payload))
            monkeypatch.setattr(_os, "getpid", lambda: 11111)   # back to process A
        return real_replace(src, dst)

    monkeypatch.setattr(_os, "getpid", lambda: 11111)
    monkeypatch.setattr(_os, "replace", replace_hook)
    D._save_ckpt(ck, dict(payload))                              # must NOT raise
    assert _os.path.exists(_os.path.join(ck, "ckpt_7.pkl"))


def test_save_ckpt_latest_is_monotonic(tmp_path):
    """A slower/stale concurrent writer must not regress latest.txt and lose deeper
    progress (so a resume always reloads the deepest checkpoint)."""
    import quiverlab.engine.deepen as D
    ck = str(tmp_path / "ck")
    base = {"cur": None, "cur_r": 1, "rks": {0: 1}, "last_gens": None,
            "HH": [], "per_degree": []}
    D._save_ckpt(ck, dict(base, n=30))
    D._save_ckpt(ck, dict(base, n=12))                          # a stale writer's old degree
    with open(os.path.join(ck, "latest.txt")) as f:
        assert int(f.read().strip()) == 30


def test_deepen_predictive_time_stop(tmp_path, monkeypatch):
    """deepen must stop cleanly (stop=time) BEFORE attempting a degree that won't finish in
    the walltime budget -- not overrun it (on the cluster an overrun = SIGKILL / TIMEOUT)."""
    import quiverlab.engine.deepen as D

    class _Clock:                       # a fake clock advanced only by _advance_resolution
        def __init__(self): self.t = 0.0
        def time(self): return self.t
    clock = _Clock()
    real_adv = D._advance_resolution

    def slow_adv(*a, **k):
        clock.t += 100.0                # each degree "costs" 100 s of walltime
        return real_adv(*a, **k)

    monkeypatch.setattr(D, "time", clock)
    monkeypatch.setattr(D, "_advance_resolution", slow_adv)

    A = truncated_polynomial(3)         # infinite (periodic) resolution -> never terminates early
    out = deepen(A, str(tmp_path / "ck"), prime=32003, time_limit_s=250)
    assert out["stop_reason"] == "time"
    assert clock.t <= 250               # stopped before overrunning the budget (old code -> 300)


def test_deepen_finalize_only_reads_checkpoint(tmp_path):
    """finalize_only rewrites the summary from the latest checkpoint WITHOUT advancing --
    the recovery path for a job that timed out before writing a fresh JSON."""
    A = truncated_polynomial(3)
    ck = str(tmp_path / "ck")
    full = deepen(A, ck, prime=32003, max_degree=4)
    fin = deepen(A, ck, prime=32003, finalize_only=True)
    assert fin["stop_reason"] == "checkpoint"
    assert fin["HH"] == full["HH"]
    assert fin["max_degree_reached"] == full["max_degree_reached"]


# not ported: test_deepen_job_finalize_only, test_deepen_job_runs_one_entry -- both shell
# out (subprocess) to hanlab/cluster/deepen_job.py and read hanlab/open_zoo_catalog_v2.json.
# cluster/ and the open-zoo catalog are EXCLUDED from Plan 02 (deferred to Plans 04/06).
