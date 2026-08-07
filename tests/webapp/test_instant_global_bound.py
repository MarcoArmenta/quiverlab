"""Fix M3 -- a process-wide ceiling on concurrent INSTANT children.

Each instant request spawns a resource-capped child. The per-IP rate gate bounds
one client's RATE but not the AGGREGATE number of children alive across many IPs,
so a wide flood could spawn unboundedly many at once. ``QLWEB_INSTANT_GLOBAL_MAX``
(config default 0 = unlimited = prior behaviour; cloud profile 8) bounds the live
instant children per process, enforced in ``webapp/server/instant.py`` via a
process-wide counter. Over the ceiling the instant tier returns a graceful 503
"busy, try again" -- never a 500, never a silent queue.

The gate is exercised directly (no real spawn children); the app wiring is checked
by monkeypatching ``run_with_timeout`` to raise ``InstantBusy``."""
import contextlib

import pytest
from fastapi.testclient import TestClient

from webapp.server import app as app_module
from webapp.server import instant as instant_mod
from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.instant import InstantBusy, _instant_slot

pytestmark = pytest.mark.fast

# QuantumCI(q=1) over GF(3) + coxeter_polynomial classifies INSTANT (see
# tests/webapp/test_acceptance.py::test_instant_path_end_to_end).
_INSTANT_BODY = {
    "schema": 1,
    "algebra": {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
                "field": {"kind": "GF", "p": 3, "n": 1}},
    "compute": ["coxeter_polynomial"],
    "artifacts": {"pdf": False, "tikz": False},
}


@pytest.fixture(autouse=True)
def _reset_inflight(monkeypatch):
    # The counter is a module global (genuinely process-wide); reset it so tests
    # never cross-contaminate.
    monkeypatch.setattr(instant_mod, "_instant_inflight", 0)


# --------------------------------------------------------------------------- #
# The config knob: additive, default 0 = unlimited.
# --------------------------------------------------------------------------- #

def test_instant_global_max_defaults_to_zero(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    assert cfg.instant_global_max == 0                 # unlimited = prior behaviour


def test_instant_global_max_is_overridable(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_GLOBAL_MAX": "8"})
    assert cfg.instant_global_max == 8


# --------------------------------------------------------------------------- #
# The gate: default (0) never bounds; a positive limit trips at N+1.
# --------------------------------------------------------------------------- #

def test_gate_disabled_at_zero_never_counts_or_trips():
    with contextlib.ExitStack() as stack:
        for _ in range(50):
            stack.enter_context(_instant_slot(0))      # 0 = disabled
        assert instant_mod._instant_inflight == 0      # disabled path counts nothing
    assert instant_mod._instant_inflight == 0


def test_gate_trips_at_n_plus_one_then_recovers():
    limit = 3
    with contextlib.ExitStack() as stack:
        for _ in range(limit):                         # fill exactly to the ceiling
            stack.enter_context(_instant_slot(limit))
        assert instant_mod._instant_inflight == limit
        with pytest.raises(InstantBusy):               # the (N+1)th is refused
            with _instant_slot(limit):
                pass
        # A refused acquire must NOT have incremented the counter.
        assert instant_mod._instant_inflight == limit
    # Every held slot is released when the stack unwinds -- no leak.
    assert instant_mod._instant_inflight == 0


def test_gate_releases_slot_even_when_body_raises():
    with pytest.raises(ValueError):
        with _instant_slot(2):
            assert instant_mod._instant_inflight == 1
            raise ValueError("boom")
    assert instant_mod._instant_inflight == 0          # released in finally


# --------------------------------------------------------------------------- #
# The app wiring: InstantBusy -> graceful 503 (not 500, not a queue).
# --------------------------------------------------------------------------- #

def test_instant_busy_becomes_a_clean_503(tmp_path, monkeypatch):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))

    def _boom(req, c):                                 # simulate a saturated ceiling
        raise InstantBusy("instant tier at capacity (8/8 children live)")

    monkeypatch.setattr(app_module, "run_with_timeout", _boom)
    r = client.post("/api/compute", json=_INSTANT_BODY)
    assert r.status_code == 503, r.text                # busy, not 500, not queued
    body = r.json()
    assert body["error_type"] == "Busy"
    assert "capacity" in body["message"].lower() or "busy" in body["message"].lower()
    # A capacity refusal must NOT leak an internal-exception 500 shape.
    assert body["error_type"] != "InternalError"


def test_instant_success_path_unaffected_when_gate_disabled(tmp_path):
    # Default cfg (global max 0) leaves the instant path byte-identical: a small
    # instant request still returns a 200 instant result (real spawn child).
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    assert cfg.instant_global_max == 0
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute", json=_INSTANT_BODY)
    assert r.status_code == 200, r.text
    assert r.json()["tier"] == "instant"
    assert "coxeter_polynomial" in r.json()["result"]["results"]
