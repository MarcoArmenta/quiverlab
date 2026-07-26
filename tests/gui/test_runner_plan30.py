"""Plan 30 -- the GUI runner's new module surface (client/Pyodide copy): the
resolution `summands` field, `decompose`, `tor`, and the tau/tau^- input
certificate. Mirrors the webapp server-tier dispatch (tests/webapp/test_runner_plan30).

The decompose/tor engines are imported LAZILY, so the branch is green whether or
not they have landed: decompose tests skip if absent (present in this worktree);
the tor test asserts the honest "engine unavailable" error when tor is absent and
the real dims when present.
"""
import json

import pytest

try:
    import quiverlab.modules.decompose as _decompose_mod  # noqa: F401
    _HAVE_DECOMPOSE = True
except ImportError:
    _HAVE_DECOMPOSE = False

try:
    import quiverlab.modules.tor as _tor_mod  # noqa: F401
    _HAVE_TOR = True
except ImportError:
    _HAVE_TOR = False

_needs_decompose = pytest.mark.skipif(
    not _HAVE_DECOMPOSE, reason="Plan-30 Part-A decompose engine not present yet")

_A2 = {"schema": 1,
       "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                   "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}},
       "compute": [], "artifacts": {"pdf": False, "tikz": True}}
_LOOP = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": [], "artifacts": {"pdf": False, "tikz": True},
         "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}}
_S1S2 = {"dims": {"1": 1, "2": 1}, "maps": {"a": [[0]]}}


def _ready(runner, req):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out


def _one(runner, spec):
    return json.loads(runner.compute_one(spec))


# --------------------------------------------------------------------------- #
# Resolution summands (Marco #3)
# --------------------------------------------------------------------------- #

def test_resolution_blocks_carry_latex_summands(runner):
    _ready(runner, _LOOP)
    pr = _one(runner, "projective_resolution:0..3")["block"]
    assert pr["summands"] == ["P_{1}", "P_{1}", "P_{1}", "P_{1}"]
    assert "betti" in pr and "terms" in pr             # raw fields retained
    ir = _one(runner, "injective_resolution:0..2")["block"]
    assert ir["summands"] == ["I_{1}", "I_{1}", "I_{1}"]


def test_summands_latex_multiplicity_and_order(runner):
    assert runner._summands_latex([1, 1, 2], "P") == "P_{1}^{2} \\oplus P_{2}"
    assert runner._summands_latex([], "I") == "0"


# --------------------------------------------------------------------------- #
# decompose
# --------------------------------------------------------------------------- #

@_needs_decompose
def test_decompose_block(runner):
    _ready(runner, dict(_A2, module=_S1S2))
    b = _one(runner, "decompose")["block"]
    assert b["iso_classes"] == 2 and len(b["summands"]) == 2
    assert all(s["indecomposable"] for s in b["summands"])
    assert len(b["citations"][0]) == 2


@_needs_decompose
def test_decompose_uncertifiable_is_a_named_error(runner):
    _ready(runner, _LOOP)                              # char 2 <= dim 2
    out = _one(runner, "decompose")
    assert out["ok"] is False
    assert out["error"]["type"] != "InternalError"    # a NAMED library error


# --------------------------------------------------------------------------- #
# tau / tau^- input certificate
# --------------------------------------------------------------------------- #

@_needs_decompose
def test_tau_certifies_decomposable_input(runner):
    _ready(runner, dict(_A2, module=_S1S2))
    b = _one(runner, "tau")["block"]
    assert b["indecomposable"] is False
    assert b["note_key"] == "mod.tau_additive" and len(b["decomposition"]) == 2


@_needs_decompose
def test_tau_certifies_indecomposable_input(runner):
    _ready(runner, dict(_A2, module={"builtin": {"kind": "simple", "vertex": 2}}))
    b = _one(runner, "tau")["block"]
    assert b["indecomposable"] is True and "decomposition" not in b


def test_tau_omits_certificate_when_uncertifiable(runner):
    # char 2 <= dim 2: cert omitted (no false claim), block still computes.
    _ready(runner, _LOOP)
    b = _one(runner, "tau")["block"]
    assert "indecomposable" not in b and "note_key" not in b
    assert b["dimvec"] == {"1": 2}


# --------------------------------------------------------------------------- #
# tor (honest either way)
# --------------------------------------------------------------------------- #

def test_tor_dispatch_is_honest_about_the_engine(runner):
    req = dict(_A2, module={"dims": {"1": 1, "2": 1}, "maps": {"a": [[1]]}},
               tor_target={"dims": {"1": 1}, "side": "left"})
    _ready(runner, req)
    out = _one(runner, "tor:0..2")
    if _HAVE_TOR:
        assert out["ok"] is True
        assert out["block"]["kind"] == "tor" and "dims" in out["block"]
    else:
        assert out["ok"] is False
        assert out["error"]["type"] != "InternalError"
        assert "Tor" in out["error"]["message"] and "Plan 29" in out["error"]["message"]


# --------------------------------------------------------------------------- #
# Worked-steps bundle for module computes (Marco #5): trace_tex/.html populated
# --------------------------------------------------------------------------- #

def test_module_worked_steps_tex_when_report_requested(runner):
    # With the report requested (artifacts.pdf), a module compute emits the Part-C
    # step events, so trace_tex()/trace_html() cover the module -- the .tex the GUI
    # download button ships (Plan 30 C1).
    req = dict(_LOOP, compute=["projective_resolution:0..3"],
               artifacts={"pdf": True, "tikz": False})
    _ready(runner, req)
    assert _one(runner, "projective_resolution:0..3")["ok"]
    tex = runner.trace_tex()
    assert tex and "P_" in tex and "projective" in tex
    assert runner.trace_html()


def test_module_trace_absent_when_report_not_requested(runner):
    # No report requested -> no module events emitted -> trace stays empty (the
    # module compute pays no re-run cost). Byte-parity with the pre-Plan-30 GUI.
    req = dict(_LOOP, compute=["projective_resolution:0..3"],
               artifacts={"pdf": False, "tikz": True})
    _ready(runner, req)
    assert _one(runner, "projective_resolution:0..3")["ok"]
    assert runner.trace_tex() == "" and runner.trace_html() == ""
