"""Plan 78 / R13 -- both runners (the wheel's ``hpc.spec`` core and the Pyodide GUI twin)
emit the ``deformations`` block byte-for-byte identical, so the browser cannot drift from
the server report / cluster run. An algebra-level kind carrying a DIM BUDGET (the
``hh1_lie`` precedent), with the plan's OWN ``DEFORM_MAXDIM = 32`` -- NOT P70's 48, a
different cost law (HH-richness, not ``A.dim``)."""
import importlib.util
import json
import pathlib

import pytest

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

pytestmark = [pytest.mark.oracle_crossengine]

# The plan's OBSTRUCTED benchmark point: dim-4 `QuantumCI(-1)`, HH^2 = 5, measured at
# 43.5 s -- while dim-20 kZ_10/J^2 with HH^2 = 0 is 0.024 s. That inversion is the whole
# reason the estimator sizes on HH-RICHNESS and not on A.dim.
#
# MIND THE q CONVENTION (this bit me while writing the test). `families.quantum.QuantumCI`
# documents its relation as `xy + q*yx`, so `QuantumCI(-1)` is the COMMUTATOR
# `x*y - y*x` -- i.e. the COMMUTATIVE k[x,y]/(x^2, y^2), HH^* = [4, 4, 5, 6]. Under the
# more common literature convention `yx = q*xy`, "q = -1" would mean ANTI-commuting, which
# in this codebase is `QuantumCI(1)` (`x*y + y*x`) and has HH^* = [2, 4, 6, 8] -- a
# DIFFERENT algebra with HH^2 = 6. Spelling the relation out here (rather than naming the
# family) keeps the fixture honest about which algebra is actually pinned.
_QCI_QQ = {"schema": 2,
           "algebra": {"kind": "quiver", "vertices": [1],
                       "arrows": {"x": [1, 1], "y": [1, 1]},
                       "relations": ["x*x", "y*y", "x*y - y*x"],
                       "field": {"kind": "QQ"}},
           "compute": ["deformations"],
           "artifacts": {"pdf": False, "tikz": False}}

# kZ_3/J^2 -- rad^2 = 0, the UNOBSTRUCTED / dg-Lie-certified regime, and the cheap side of
# the cost law (HH^2 thin). Also exercises the explicit budget suffix.
_RADSQ_QQ = {"schema": 2,
             "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                         "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                         "relations": ["a*b", "b*c", "c*a"],
                         "field": {"kind": "QQ"}},
             "compute": ["deformations:32"],
             "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, spec):
    s = importlib.util.spec_from_file_location("gui_runner_deformations_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(spec))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["deformations"]


def _assert_twins(body, tmp_path):
    gui = _gui_block(body, body["compute"][0])
    ref = _server_block(body, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
    return ref


@pytest.mark.deep
def test_deformations_twin_obstructed_quantum_ci(tmp_path):
    # DEEP by explicit marker (conftest: an explicit bucket wins over the directory
    # auto-assignment, which would put tests/gui/ in `fast`). Measured 85 s -- a 25 % hit
    # on the whole fast suite, on EVERY OS x Python matrix cell, for one algebra. The
    # cheap rad^2 = 0 twin below (0.03 s) stays in `fast`, so the cross-runner contract is
    # still gated everywhere; only the expensive obstruction pin rides the deep leg. The
    # 85 s vs 0.03 s split between these two tests IS the HH-richness cost law the
    # estimator sizes on.
    b = _assert_twins(_QCI_QQ, tmp_path)
    assert b["kind"] == "deformations" and b["status"] == "complete"
    assert b["hh2_dim"] == 5 and b["hh3_dim"] == 6
    # OBSTRUCTED: some basis direction has [alpha, alpha] != 0, so it does not extend.
    assert b["unobstructed"] is False
    assert b["obstruction_witness"]


def test_deformations_twin_radsq_budget_suffix(tmp_path):
    b = _assert_twins(_RADSQ_QQ, tmp_path)
    assert b["kind"] == "deformations"
    # rad^2 = 0 => the B(A)[1] L-infinity companion is UNCONDITIONALLY dg-Lie (Plan 78's
    # guaranteed deliverable; no l_>=3 is ever shipped, the spike stays frozen).
    assert b["dg_lie"] is True
    assert b["l3_status"] is not None
