"""Plan 26 -- module compute-kind runner tests (webapp server tier).

The runner builds a module from the request's `module` block (explicit dims+maps,
or a builtin S/P/I pick-list; right or left) and dispatches the module-level
invariants: dimension vector, rad/top/soc, Ext^n(M,N), tau/tau^-, projective &
injective resolutions, and projective/injective dimension. Matrices are exact
DATA (never evaluated, never floats); a module that violates the relations raises
the library's loud error, surfaced as a clean typed refusal (never a 500).
"""
import json

from webapp.server.runner import run_spec, RunError
from webapp.server.schema import ComputeRequest

# k[x]/(x^3) over GF(2): one vertex, a nilpotent loop -- a compact module zoo.
_LOOP = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
         "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}
# The A_2 quiver 1 --a--> 2, no relations: a genuine multi-vertex algebra.
_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 2, "n": 1}}
# a length-2 nilpotent module M over k[x]/(x^3): x = [[0,0],[1,0]].
_M2 = {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}


def _run(tmp_path, algebra, compute, module=None, ext_target=None):
    body = {"schema": 2, "algebra": algebra, "compute": compute,
            "artifacts": {"pdf": False, "tikz": False}}
    if module is not None:
        body["module"] = module
    if ext_target is not None:
        body["ext_target"] = ext_target
    return run_spec(ComputeRequest.model_validate(body), tmp_path)


def test_dimension_vector(tmp_path):
    r = _run(tmp_path, _LOOP, ["dimension_vector"], _M2)
    b = r["results"]["dimension_vector"]
    assert b["dimvec"] == {"1": 2} and b["dim"] == 2 and b["side"] == "right"
    # references + resolved citations attached, like every algebra invariant.
    assert b["references"] and len(b["citations"][0]) == 2


def test_rad_top_soc(tmp_path):
    b = _run(tmp_path, _LOOP, ["rad_top_soc"], _M2)["results"]["rad_top_soc"]
    # Plan 34 (Marco): each of rad/top/soc is a FULL representation -- the dim
    # VECTOR (`dims`) + the exact per-arrow action matrices (`maps`); the redundant
    # total-dim field is gone. The shape mirrors the module INPUT (feedable back in).
    assert b["radical"]["dims"] == {"1": 1}
    assert b["top"]["dims"] == {"1": 1}
    assert b["socle"]["dims"] == {"1": 1}
    assert b["radical"]["maps"] == {"x": [[0]]}
    assert b["top"]["maps"] == {"x": [[0]]}
    assert b["socle"]["maps"] == {"x": [[0]]}
    assert "dim" not in b["radical"] and "dimvec" not in b["radical"]


def test_tau_and_tau_minus(tmp_path):
    r = _run(tmp_path, _LOOP, ["tau", "tau_minus"], _M2)
    assert r["results"]["tau"]["dimvec"] == {"1": 2}
    assert r["results"]["tau_minus"]["dimvec"] == {"1": 2}
    assert r["results"]["tau"]["is_zero"] is False


def test_tau_of_projective_is_zero(tmp_path):
    # tau(P) = 0 for a projective (here P_1 of k[x]/x^3 has dim 3).
    b = _run(tmp_path, _LOOP, ["tau"],
             {"builtin": {"kind": "projective", "vertex": 1}})["results"]["tau"]
    assert b["is_zero"] is True and b["dim"] == 0


def test_projective_and_injective_resolution(tmp_path):
    r = _run(tmp_path, _LOOP, ["projective_resolution:0..3", "injective_resolution:0..3"], _M2)
    pr = r["results"]["projective_resolution"]
    assert pr["top"] == 3 and len(pr["terms"]) >= 4 and "pd" in pr and "betti" in pr
    ir = r["results"]["injective_resolution"]
    assert ir["top"] == 3 and "injective_dimension" in ir


def test_projective_and_injective_dimension(tmp_path):
    # S_1 of k[x]/x^3 has infinite projective AND injective dimension (periodic).
    r = _run(tmp_path, _LOOP, ["projective_dimension", "injective_dimension"],
             {"builtin": {"kind": "simple", "vertex": 1}})
    pd = r["results"]["projective_dimension"]
    idim = r["results"]["injective_dimension"]
    assert pd["value"] is None and pd["finite"] is False
    assert idim["value"] is None and idim["finite"] is False


def test_finite_projective_dimension_over_A2(tmp_path):
    # S_2 over the A_2 quiver is projective, so pd(S_2) = 0 (finite).
    b = _run(tmp_path, _A2, ["projective_dimension"],
             {"builtin": {"kind": "simple", "vertex": 2}})["results"]["projective_dimension"]
    assert b["value"] == 0 and b["finite"] is True


def test_ext_with_two_modules(tmp_path):
    r = _run(tmp_path, _LOOP, ["ext:0..4"], _M2,
             ext_target={"builtin": {"kind": "simple", "vertex": 1}})
    b = r["results"]["ext"]
    assert b["top"] == 4 and b["dims"] == [1, 1, 1, 1, 1]
    assert b["target"]["dimvec"] == {"1": 1}
    assert b["references"] == ["module_ext"]


def test_builtin_pick_lists(tmp_path):
    for kind, expect in (("simple", 1), ("projective", 3), ("injective", 3)):
        b = _run(tmp_path, _LOOP, ["dimension_vector"],
                 {"builtin": {"kind": kind, "vertex": 1}})["results"]["dimension_vector"]
        assert b["dim"] == expect


def test_multi_vertex_block_matrix(tmp_path):
    # arrow a:1->2 acts by the 1x1 block [[1]] (target x source dims); the runner
    # places it into the full vertex-ordered action matrix. M is then P_1 (dim 2).
    r = _run(tmp_path, _A2, ["dimension_vector"],
             {"dims": {"1": 1, "2": 1}, "maps": {"a": [[1]]}})
    assert r["results"]["dimension_vector"]["dimvec"] == {"1": 1, "2": 1}


def test_left_module_side(tmp_path):
    b = _run(tmp_path, _LOOP, ["dimension_vector"],
             dict(_M2, side="left"))["results"]["dimension_vector"]
    assert b["side"] == "left" and b["dimvec"] == {"1": 2}


def test_relation_violation_is_a_typed_refusal_not_a_500(tmp_path):
    # x acting invertibly ([[1]]) violates x^3 = 0 -> the library's loud module
    # error, surfaced as a QuiverlabError-typed RunError (a SAFE error type, so the
    # app returns a clean 4xx, never a 500, never a silent wrong answer).
    try:
        _run(tmp_path, _LOOP, ["dimension_vector"], {"dims": {"1": 1}, "maps": {"x": [[1]]}})
        assert False, "expected a loud refusal"
    except RunError as exc:
        from webapp.server.security import is_safe_error_type
        assert is_safe_error_type(exc.error_type), exc.error_type
        assert "relation" in exc.message.lower()


def test_wrong_shape_block_is_a_schema_refusal(tmp_path):
    try:
        _run(tmp_path, _A2, ["dimension_vector"],
             {"dims": {"1": 1, "2": 1}, "maps": {"a": [[1, 0]]}})
        assert False, "expected a shape refusal"
    except RunError as exc:
        assert exc.error_type == "SchemaError" and "1x1" in exc.message


def test_unknown_arrow_and_vertex_refused(tmp_path):
    for module in ({"dims": {"1": 1}, "maps": {"zzz": [[0]]}},
                   {"dims": {"9": 1}, "maps": {}}):
        try:
            _run(tmp_path, _LOOP, ["dimension_vector"], module)
            assert False, "expected a schema refusal"
        except RunError as exc:
            assert exc.error_type == "SchemaError"


def test_module_result_json_and_reproduce_execs(tmp_path):
    r = _run(tmp_path, _LOOP, ["dimension_vector", "tau"], _M2)
    on_disk = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert on_disk["results"]["dimension_vector"]["dimvec"] == {"1": 2}
    # The copy-paste reproduce snippet is runnable and rebuilds the module.
    ns = {}
    exec(r["reproduce"], ns)
    assert ns["M"].dim == 2


def test_top_level_references_aggregate_module_citations(tmp_path):
    r = _run(tmp_path, _LOOP, ["dimension_vector", "projective_resolution:0..2"], _M2)
    keys = {e["key"] for e in r["references"]}
    assert {"assem_book", "minimal_resolution"} <= keys
