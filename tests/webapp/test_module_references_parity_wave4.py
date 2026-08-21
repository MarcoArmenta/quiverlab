"""Wave 4 (v1.0.1, 2026-08-21): the two runners must stamp the SAME ``references``
on every module-kind block.

The server core ``quiverlab.hpc.spec`` runs every module block through
``_with_refs`` -> ``block["references"] = list(_MOD_REFS[kind])`` (the server then
AGGREGATES those keys into the report's References section, spec.py ~L1367). Its
Pyodide twin ``docs/gui/runner.py`` computed the very same ``keys`` but only put
them in ``citations`` -- the ``references`` key was DROPPED for the default module
kinds (decompose, tau, tau_minus, rad_top_soc, dimension_vector, projective/
injective dimension, ext, tor, projective/injective resolution). Neither renderer
reads ``references``, so nothing was visibly wrong -- but the twins are meant to
emit byte-identical blocks, and only the ALGEBRA kinds + barcode/almost_split/
tilting/orbit carried it. The two ``_MOD_REFS`` tables are already byte-identical,
so parity is a one-line-per-block addition on the twin (no server change, so the
Plan-28 delegation golden is untouched).

NO oracle-class marker: tests/webapp/ collects only with the [web] extra (the
Plan-32 audited class counts must stay environment-independent); this is a
cross-runner block-shape contract test.
"""
import importlib.util
import json
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec

_A3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a": [1, 2], "b": [2, 3]},
       "relations": ["a*b"], "field": {"kind": "GF", "p": 5}}
_MOD = {"builtin": {"kind": "projective", "vertex": 1}}
_N = {"builtin": {"kind": "simple", "vertex": 2}}

# The default module kinds that were missing ``references`` on the twin.
_KINDS = ("decompose", "tau", "tau_minus", "rad_top_soc", "dimension_vector",
          "projective_dimension", "injective_dimension", "ext",
          "projective_resolution", "injective_resolution")


def _req(kind):
    r = {"schema": 2, "algebra": _A3, "module": _MOD}
    if kind == "ext":
        r["ext_target"] = _N
        r["compute"] = ["ext:0..3"]
    elif kind in ("projective_resolution", "injective_resolution"):
        r["compute"] = [kind + ":0..3"]
    else:
        r["compute"] = [kind]
    return r


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _load_twin():
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_refparity", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    return mod


def _pyodide(req):
    mod = _load_twin()
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


@pytest.mark.parametrize("kind", _KINDS)
def test_twin_stamps_the_same_references_as_the_server(kind):
    sb = _server(_req(kind))[kind]
    gb = _pyodide(_req(kind))[kind]
    # the server has always carried references; the twin must now match them.
    assert sb.get("references"), f"server {kind} block lost its references"
    assert gb.get("references") == sb.get("references"), (
        f"{kind}: twin references {gb.get('references')!r} != server "
        f"{sb.get('references')!r}")
    # citations were already identical -- keep that a regression guard.
    assert gb.get("citations") == sb.get("citations"), kind


def test_mod_refs_tables_are_byte_identical():
    twin = _load_twin()
    assert twin._MOD_REFS == spec._MOD_REFS
