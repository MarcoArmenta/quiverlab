"""Plan 53 GUI/webapp cross-runner contract: the extended ``homological_profile``
block (phidim/psidim/phi_spectrum/lit additive keys) + the new ``fractional_cy``
compute kind. Both runners -- ``quiverlab.hpc.spec`` (server / container / CLI) and
its Pyodide twin ``docs/gui/runner.py`` -- must agree BYTE-FOR-BYTE, or the draw page
and the desktop app disagree about the same computation.

NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
Plan-32 audit requires the audited class counts to be environment-independent. These
are contract tests for the two runners' block shapes."""
import json
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec

# kA3/rad^2 (gl.dim 2, rep-finite): phidim = psidim = 2, spectrum [0,1,2], LIT gorenstein.
_KA3 = {"kind": "quiver", "vertices": [1, 2, 3],
        "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": ["a*b"],
        "field": {"kind": "GF", "p": 32003}}
# k[x]/(x^3): symmetric local Nakayama, stable CY dimension 1/1 (weakly 1-CY).
_KXX3 = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
         "relations": ["x*x*x"], "field": {"kind": "GF", "p": 32003}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p53", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def test_homological_profile_extended_shape():
    req = {"schema": 1, "algebra": _KA3, "compute": ["homological_profile"]}
    b = _server(req)["homological_profile"]
    # the four Plan-53 additive keys are present with values (kA3/rad^2 is rep-finite)
    assert b["phidim"]["value"] == 2 and b["phidim"]["exact"] is True
    assert b["psidim"]["value"] == 2
    assert b["phi_spectrum"]["values"] == [0, 1, 2]
    assert b["phi_spectrum"]["gaps"] == [] and b["phi_spectrum"]["complete"] is True
    assert b["lit"]["family"] == "gorenstein" and b["lit"]["findim_upper"] == 2


def test_fractional_cy_block_shape():
    req = {"schema": 1, "algebra": _KXX3, "compute": ["fractional_cy"]}
    b = _server(req)["fractional_cy"]
    assert b["cy_dimension"] == "1/1" and (b["m"], b["ell"]) == (1, 1)
    assert b["weakly_n_cy"] == 1 and b["status"] == "certified"
    assert b["tier"] == "weak-on-generators"        # the display never overstates it


@pytest.mark.parametrize("req", [
    {"schema": 1, "algebra": _KA3, "compute": ["homological_profile"]},
    {"schema": 1, "algebra": _KXX3, "compute": ["fractional_cy"]},
])
def test_runners_byte_identical(req):
    k = req["compute"][0]
    sv = _server(req)[k]
    py = _pyodide(req)[k.split(":")[0]]
    assert (json.dumps(sv, sort_keys=True, default=str)
            == json.dumps(py, sort_keys=True, default=str))


def test_fractional_cy_non_self_injective_is_error_entry():
    # kA3 (hereditary, NOT self-injective): a clean typed error entry, never a 500.
    kA3_lin = {"kind": "quiver", "vertices": [1, 2, 3],
               "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
               "field": {"kind": "GF", "p": 32003}}
    req = {"schema": 1, "algebra": kA3_lin, "compute": ["fractional_cy"]}
    b = _server(req)["fractional_cy"]
    assert "error" in b and "self-injective" in b["error"]
    assert "references" not in b            # no citations on the error-entry shape
