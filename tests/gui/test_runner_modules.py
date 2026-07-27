"""Plan 26 -- the GUI runner's module compute kinds (client/Pyodide copy).

Mirrors the webapp server-tier module dispatch: build a module from the request's
module block (explicit dims+maps, or a builtin S/P/I pick-list; right/left),
dispatch dimension_vector / rad_top_soc / tau / tau^- / Ext / projective &
injective resolutions / projective & injective dimension. Block shapes carry a
`latex` display and `citations`, like the algebra invariants; a relation-violating
module is a NAMED library error, never a bare InternalError.
"""
import json

# k[x]/(x^3) over GF(2), plus a length-2 nilpotent module M (x = [[0,0],[1,0]]).
_LOOP = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": [], "artifacts": {"pdf": False, "tikz": True},
         "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}}
# A_2: 1 --a--> 2, no relations -- a genuine multi-vertex algebra.
_A2 = {"schema": 1,
       "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                   "relations": [], "field": {"kind": "GF", "p": 2, "n": 1}},
       "compute": [], "artifacts": {"pdf": False, "tikz": True}}


def _ready(runner, req):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out
    return out


def _one(runner, spec):
    return json.loads(runner.compute_one(spec))


def test_module_goldens(runner):
    _ready(runner, _LOOP)
    dv = _one(runner, "dimension_vector")["block"]
    assert dv["dimvec"] == {"1": 2} and dv["dim"] == 2 and dv["side"] == "right"
    assert dv["latex"] and len(dv["citations"][0]) == 2       # [key, formatted]
    rts = _one(runner, "rad_top_soc")["block"]
    # Plan 34 (Marco): rad/top/soc are FULL representations -- dim VECTOR (`dims`)
    # + exact per-arrow matrices (`maps`); no redundant total-dim field.
    assert rts["radical"]["dims"] == {"1": 1}
    assert rts["socle"]["dims"] == {"1": 1}
    assert rts["radical"]["maps"] == {"x": [[0]]}
    assert "dim" not in rts["radical"] and "dimvec" not in rts["radical"]
    assert _one(runner, "tau")["block"]["dimvec"] == {"1": 2}
    assert _one(runner, "tau_minus")["block"]["dimvec"] == {"1": 2}
    pr = _one(runner, "projective_resolution:0..3")["block"]
    assert pr["top"] == 3 and len(pr["terms"]) >= 4 and "pd" in pr
    ir = _one(runner, "injective_resolution:0..3")["block"]
    assert "injective_dimension" in ir
    assert _one(runner, "projective_dimension")["block"]["value"] is None    # infinite
    assert _one(runner, "injective_dimension")["block"]["finite"] is False


def test_ext_with_two_modules(runner):
    req = dict(_LOOP, ext_target={"builtin": {"kind": "simple", "vertex": 1}})
    _ready(runner, req)
    ext = _one(runner, "ext:0..4")["block"]
    assert ext["dims"] == [1, 1, 1, 1, 1] and ext["target"]["dimvec"] == {"1": 1}


def test_builtin_pick_lists(runner):
    for kind, want in (("simple", 1), ("projective", 3), ("injective", 3)):
        req = dict(_LOOP, module={"builtin": {"kind": kind, "vertex": 1}})
        _ready(runner, req)
        assert _one(runner, "dimension_vector")["block"]["dim"] == want


def test_tau_of_projective_is_zero(runner):
    req = dict(_LOOP, module={"builtin": {"kind": "projective", "vertex": 1}})
    _ready(runner, req)
    b = _one(runner, "tau")["block"]
    assert b["is_zero"] is True and "0" in b["latex"]


def test_multi_vertex_block_matrix(runner):
    # arrow a:1->2 acts by the 1x1 block [[1]]; runner places it in the full action.
    req = dict(_A2, module={"dims": {"1": 1, "2": 1}, "maps": {"a": [[1]]}})
    _ready(runner, req)
    assert _one(runner, "dimension_vector")["block"]["dimvec"] == {"1": 1, "2": 1}


def test_left_module_side(runner):
    req = dict(_LOOP, module=dict(_LOOP["module"], side="left"))
    _ready(runner, req)
    b = _one(runner, "dimension_vector")["block"]
    assert b["side"] == "left" and b["dimvec"] == {"1": 2}


def test_relation_violation_is_named_not_internal(runner):
    req = dict(_LOOP, module={"dims": {"1": 1}, "maps": {"x": [[1]]}})
    _ready(runner, req)
    out = _one(runner, "dimension_vector")
    assert out["ok"] is False
    assert out["error"]["type"] != "InternalError"           # a NAMED library error
    assert "relation" in out["error"]["message"].lower()


def test_wrong_shape_is_named_error(runner):
    req = dict(_A2, module={"dims": {"1": 1, "2": 1}, "maps": {"a": [[1, 0]]}})
    _ready(runner, req)
    out = _one(runner, "dimension_vector")
    assert out["ok"] is False and out["error"]["type"] == "RequestError"


def test_module_snippet_reproduces(runner):
    req = dict(_LOOP, compute=["dimension_vector", "tau"])
    _ready(runner, req)
    _one(runner, "dimension_vector")
    src = runner.python_snippet()
    ns = {}
    exec(src, ns)                                             # runs end to end
    assert ns["M"].dim == 2


def test_gui_request_shape_with_string_entries(runner):
    # The EXACT request gui.js buildRequest() emits: string matrix entries (an empty
    # cell normalizes to "0"), string dims keys, a builtin ext_target -- all DATA.
    req = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                       "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["dimension_vector", "ext:0..2"],
           "artifacts": {"pdf": False, "tikz": True},
           "module": {"dims": {"1": 2}, "maps": {"x": [["0", "0"], ["1", "0"]]},
                      "side": "right"},
           "ext_target": {"builtin": {"kind": "simple", "vertex": 1}, "side": "right"}}
    _ready(runner, req)
    assert _one(runner, "dimension_vector")["block"]["dimvec"] == {"1": 2}
    assert _one(runner, "ext:0..2")["block"]["dims"] == [1, 1, 1]


def test_module_kinds_in_result_bundle_and_estimate(runner):
    _ready(runner, dict(_LOOP, compute=["dimension_vector"]))
    _one(runner, "dimension_vector")
    bundle = json.loads(runner.result_bundle())
    assert any(r.get("invariant") == "dimension_vector" for r in bundle["results"])
    est = json.loads(runner.estimate(1.0))
    assert est["ok"] and any(b["invariant"] == "dimension_vector" for b in est["breakdown"])
