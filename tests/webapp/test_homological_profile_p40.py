"""The ``homological_profile`` scalar compute kind (Plan 40 / C6): the whole
homological-dimension family in ONE no-code block -- served by ``quiverlab.hpc.spec``
(server / container / CLI) and mirrored byte-for-byte by its Pyodide twin
``docs/gui/runner.py`` (the draw page / desktop app). If the two disagree, the same
computation reads differently in the two surfaces.

NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
Plan-32 audit requires the audited class counts to be environment-independent (this
is a cross-runner contract test, like test_module_blocks_m0729)."""
import json
import tempfile

from quiverlab.hpc import spec

# kA2 over GF(7): char 7 > dim of each simple, so the Igusa-Todorov phi/psi of the
# sum of the simples is a VALUE, not the char-caveat error entry.
_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 7}}
# k[x]/(x^3) over GF(2): SELF-INJECTIVE, so the devil's-advocate theorem guard
# serves phi = psi = 0 directly (see test_selfinjective_gf2_profile...); the
# per-entry-error contract now lives on the non-self-injective _RADSQ2_GF2.
_LOOP3 = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
          "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2}}


def _req(alg):
    return {"schema": 1, "algebra": alg, "compute": ["homological_profile"]}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]["homological_profile"]


def _pyodide(req):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p40", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    r = json.loads(mod.compute_one("homological_profile"))
    assert r["ok"], r
    return r["block"]


def test_profile_block_shape():
    block = _server(_req(_A2))
    assert block["kind"] == "homological_profile"
    assert block["global_dimension"]["value"] == 1
    assert block["global_dimension"]["exact"] is True
    assert block["finitistic"]["exact"] is True
    assert block["finitistic"]["lower"] == 1 and block["finitistic"]["upper"] == 1
    assert block["dominant"]["value"] == 1 and block["dominant"]["infinite"] is False
    assert block["gorenstein"]["is_gorenstein"] is True
    assert "igusa_todorov" in block
    it = block["igusa_todorov"]
    assert it["module"] == "(+)_v S_v" and it["phi"] == 1 and it["psi"] == 1


def test_char_caveat_reports_a_per_entry_error_not_a_crash(monkeypatch):
    # Devil's-advocate repoint (2026-08-05): k[x]/(x^3) over GF(2) is
    # SELF-INJECTIVE, so the theorem guard now serves the CORRECT phi=psi=0
    # (see the companion test). The per-entry-error CONTRACT -- an IT refusal
    # degrades to {"error": ...} and never kills the rest of the profile --
    # is exercised directly: any QuiverlabError from phi must be contained.
    # (A natural char-caveat fixture exists -- local radsq two loops over
    # GF(2) -- but its OTHER profile entries walk exponentially growing
    # syzygies to the depth guards, which is a minutes-long test; the
    # contract under test is containment, not the arithmetic trigger.)
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules import homdims as _hd

    def _refuse(M, budget=512, bound=64):
        raise QuiverlabError("decompose cannot certify over GF(2) (test)")

    monkeypatch.setattr(_hd, "igusa_todorov_phi", _refuse)
    block = _server(_req(_LOOP3))
    assert block["kind"] == "homological_profile"
    assert "error" in block["igusa_todorov"]
    assert block["dominant"]["infinite"] is True     # the rest still computed


def test_selfinjective_gf2_profile_reports_theorem_values():
    # the old fixture, now served CORRECTLY by the self-injective theorem
    # guard (phi = psi = 0), where it previously degraded to a per-entry
    # error; the dominant-dimension certificate is unchanged.
    block = _server(_req(_LOOP3))
    assert block["igusa_todorov"].get("phi") == 0
    assert block["igusa_todorov"].get("psi") == 0
    assert "error" not in block["igusa_todorov"]
    assert block["dominant"]["infinite"] is True and block["dominant"]["value"] is None


def test_twin_parity():
    # both runners byte-identical on the block (strip the server-only `invariant` tag).
    for alg in (_A2, _LOOP3):
        req = _req(alg)
        s = {k: v for k, v in _server(req).items() if k != "invariant"}
        p = {k: v for k, v in _pyodide(req).items() if k != "invariant"}
        assert json.dumps(s, sort_keys=True) == json.dumps(p, sort_keys=True)
