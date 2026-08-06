"""The orbit_geometry module kind (Plan 49 / C8): schema-2 gated, served by
hpc.spec, mirrored by the Pyodide twin. Reports orbit dim + rigidity + honest
codim; the Kac canonical decomposition is present on hereditary Dynkin, refused
(None + note) off it. Cross-runner byte-identity via the shared builder."""
import importlib.util
import json
import pathlib
import tempfile

from quiverlab.hpc import spec

_kA3 = {"kind": "quiver", "vertices": [1, 2, 3],
        "arrows": {"a1": [1, 2], "a2": [2, 3]}, "relations": [],
        "field": {"kind": "CC"}}
_kxx = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
        "relations": ["x*x"], "field": {"kind": "CC"}}


def _server_block(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]["orbit_geometry"]


def _pyodide_block(req):
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_twin_p49", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    r = json.loads(mod.compute_one("orbit_geometry"))
    assert r["ok"], r
    return r["block"]


def test_orbit_geometry_block_shape():
    # kA3 over CC, module = builtin simple S_2, compute ["orbit_geometry"].
    req = {"schema": 2, "algebra": _kA3,
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["orbit_geometry"]}
    block = _server_block(req)
    assert isinstance(block["orbit_dim"], int)
    assert block["rigid"] is True                    # S_2 a real-root brick
    assert block["ext1_self"] == 0
    assert block["hereditary"] is True
    assert block["codim_semantics"] == "hereditary"
    cd = block["canonical_decomposition"]
    assert isinstance(cd, list) and cd               # non-empty: dim (0,1,0) = S_2 itself
    assert list(cd[0]["root"]) == [0, 1, 0] and cd[0]["name"] == "S_2"
    # the block carries the Voigt reference + resolved citations
    assert "voigt_rigidity" in block["references"]
    assert any(k == "Voigt1977" for k, _ in block["citations"])
    # the orbit-dim identity: dim O_M = group_dim - end_dim
    assert block["orbit_dim"] == block["group_dim"] - block["end_dim"]


def test_orbit_geometry_codim_gloss_general():
    # a NON-hereditary algebra (k[x]/(x^2), module S) => general codim semantics,
    # canonical decomposition refused loudly (None + a note mentioning "hereditary").
    req = {"schema": 2, "algebra": _kxx,
           "module": {"builtin": {"kind": "simple", "vertex": 1}},
           "compute": ["orbit_geometry"]}
    block = _server_block(req)
    assert block["hereditary"] is False
    assert block["codim_semantics"] == "general"
    assert block["canonical_decomposition"] is None
    assert "hereditary" in block["canonical_note"]
    assert block["ext1_self"] == 1 and block["rigid"] is False


def test_twin_parity():
    # both runners byte-identical via the shared orbit_geometry_block builder.
    for algebra, builtin in ((_kA3, {"kind": "simple", "vertex": 2}),
                             (_kxx, {"kind": "simple", "vertex": 1})):
        req = {"schema": 2, "algebra": algebra, "module": {"builtin": builtin},
               "compute": ["orbit_geometry"]}
        s = _server_block(req)
        p = _pyodide_block(req)
        assert json.dumps(s, sort_keys=True) == json.dumps(p, sort_keys=True)
