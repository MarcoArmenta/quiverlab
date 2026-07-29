"""Module-block fixes from Marco's 2026-07-29 pass over the desktop app.

  * ``projective_dimension`` / ``injective_dimension`` carried NO ``latex``, so the
    draw page typeset a literal "undefined" under each (his example-a showed two);
  * an unresolved probe is reported as a certified LOWER BOUND, never a bare
    ``\\infty`` the engine did not prove;
  * a ``tau`` / ``tau_minus`` request now also carries the SECOND module N's
    translate, with its full per-arrow matrices.

Both runners are checked together: ``quiverlab.hpc.spec`` (server / container /
CLI) and its Pyodide twin ``docs/gui/runner.py`` must agree, or the draw page and
the desktop app disagree about the same computation.
"""
import json
import tempfile

import pytest

from quiverlab.hpc import spec

# NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
# Plan-32 audit requires the audited class counts to be environment-independent.
# These are contract tests for the two runners' block shapes.

_A2 ={"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 5}}
# k[x]/(x^3): a loop, so every non-projective module has INFINITE projective
# dimension -- the unresolved-probe branch.
_LOOP = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
         "relations": ["x*x*x"], "field": {"kind": "GF", "p": 5}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_m0729", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


# --------------------------------------------------------------------------- #
# pd / id: a display latex, and an honest unresolved branch
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("run", [_server, _pyodide])
def test_finite_homological_dimensions_state_the_equality(run):
    req = {"schema": 2, "algebra": _A2,
           "module": {"builtin": {"kind": "projective", "vertex": 1}},
           "compute": ["projective_dimension", "injective_dimension"]}
    res = run(req)
    assert res["projective_dimension"]["latex"] == r"\operatorname{pd} M = 0"
    assert res["projective_dimension"]["value"] == 0
    assert "note" not in res["projective_dimension"]      # nothing to qualify
    assert res["injective_dimension"]["latex"].startswith(r"\operatorname{id} M = ")


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_unresolved_probe_states_a_lower_bound_not_infinity(run):
    """The resolution simply did not terminate by the probed depth; claiming
    ``= \\infty`` would assert what was not computed."""
    req = {"schema": 2, "algebra": _LOOP,
           "module": {"builtin": {"kind": "simple", "vertex": 1}},
           "compute": ["projective_dimension"]}
    b = run(req)["projective_dimension"]
    assert b["value"] is None and b["finite"] is False
    assert b["latex"] == r"\operatorname{pd} M > 32" == \
        r"\operatorname{pd} M > %d" % b["bound"]
    assert r"\infty" not in b["latex"]
    assert "certified lower bound" in b["note"]


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_no_block_ever_lacks_its_display_latex(run):
    """The root cause of the two "undefined" lines: a renderer typesets
    ``block.latex``, so every latex-rendered kind must ship one."""
    req = {"schema": 2, "algebra": _A2,
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["dimension_vector", "tau", "tau_minus",
                       "projective_dimension", "injective_dimension"]}
    res = run(req)
    for kind in ("dimension_vector", "tau", "tau_minus",
                 "projective_dimension", "injective_dimension"):
        assert res[kind].get("latex"), kind


# --------------------------------------------------------------------------- #
# tau / tau^- of the SECOND module N
# --------------------------------------------------------------------------- #
_WITH_TARGET = {
    "schema": 2, "algebra": _A2,
    "module": {"builtin": {"kind": "injective", "vertex": 2}},
    "ext_target": {"builtin": {"kind": "simple", "vertex": 2}},
    "compute": ["tau", "tau_minus"],
}


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_tau_covers_the_second_module_with_full_matrices(run):
    res = run(_WITH_TARGET)
    for kind, sym in (("tau", r"\tau N"), ("tau_minus", r"\tau^{-} N")):
        targets = res[kind]["targets"]
        assert len(targets) == 1
        t = targets[0]
        assert t["name"] == "N" and t["role"] == "ext_target"
        assert sym in t["latex"]
        # a non-zero translate ships as a FULL representation, like M's
        if not t["is_zero"]:
            assert set(t["repr"]) >= {"dims", "maps"}


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_no_target_named_means_no_targets_key(run):
    """Byte-stability: a single-module request's tau block is unchanged."""
    req = {"schema": 2, "algebra": _A2,
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["tau"]}
    assert "targets" not in run(req)["tau"]


def test_a_target_that_refuses_is_reported_not_fatal(monkeypatch):
    """tau M is already computed and valid, so a loud refusal on N degrades to an
    honest entry rather than losing the whole block."""
    from quiverlab.errors import QuiverlabError

    def boom(mod, kind, name):
        if name == "N":
            raise QuiverlabError("no translate here")
        return {"name": name, "latex": "x", "is_zero": True, "dimvec": {}, "dim": 0}

    monkeypatch.setattr(spec, "_ar_translate", boom)
    got = spec._target_translates("tau", [("ext_target", object())])
    assert got == [{"name": "N", "error": "no translate here", "role": "ext_target"}]


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_ext_still_refuses_loudly_when_its_target_is_missing(run):
    """The opportunistic target build for tau must not turn a missing ext_target
    into a silent None fed to ext_dims."""
    req = {"schema": 2, "algebra": _A2,
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["ext:0..2"]}
    with pytest.raises(Exception):
        run(req)
