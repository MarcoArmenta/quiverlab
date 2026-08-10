"""Plan 74 Task G1: the SkewGroupAlgebra construction family + skew_group_hh kind,
served by the webapp/hpc tier. Unmarked per the Plan-32 extras-gated ruling (a
cross-runner/server contract, not an oracle)."""
import json

from webapp.server.cache import canonical_key
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest


def _run(body, tmp_path):
    return run_spec(ComputeRequest.model_validate(body), tmp_path)


def _family(params, field, compute):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": "SkewGroupAlgebra",
                        "params": params, "field": field},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False}}


_FLAGSHIP_PARAMS = {"vertices": [1], "arrows": {"x": [1, 1]}, "relations": ["x*x"],
                    "generators": [{"vertex_perm": {"1": 1}, "arrow_perm": {"x": "x"},
                                    "arrow_scalars": {"x": -1}}]}


def test_flagship_block(tmp_path):
    body = _family(_FLAGSHIP_PARAMS, {"kind": "QQ"}, ["skew_group_hh:3"])
    block = _run(body, tmp_path)["results"]["skew_group_hh"]
    assert block["dim"] == 4
    assert block["group_order"] == 2
    assert block["dims"] == [1, 1, 1, 1]
    assert block["direct_dims"] == [1, 1, 1, 1]
    assert block["agrees"] is True
    assert block["status"] == "complete"
    assert block["char_ok"] is True
    assert "stefan_hopf_galois" in block["references"]
    assert "StefanHopfGalois1995" in [c[0] for c in block["citations"]]
    classes = {s["class"]: s["inv"] for s in block["decomposition"]}
    assert classes["e"] == [1, 1, 1, 1] and classes["g1"] == [0, 0, 0, 0]


def test_generator_order_canonical_key(tmp_path):
    # two orderings of an S3 action -> the SAME canonical key (order-normalized)
    g1 = {"vertex_perm": {"1": 2, "2": 3, "3": 1}, "arrow_perm": {}}
    g2 = {"vertex_perm": {"1": 2, "2": 1, "3": 3}, "arrow_perm": {}}
    base = {"vertices": [1, 2, 3], "arrows": {}, "relations": []}
    b_ab = _family({**base, "generators": [g1, g2]}, {"kind": "GF", "p": 7}, ["skew_group_hh:0"])
    b_ba = _family({**base, "generators": [g2, g1]}, {"kind": "GF", "p": 7}, ["skew_group_hh:0"])
    r_ab = ComputeRequest.model_validate(b_ab)
    r_ba = ComputeRequest.model_validate(b_ba)
    k_ab = canonical_key(r_ab.model_dump(by_alias=True), "0.1.0.dev0")
    k_ba = canonical_key(r_ba.model_dump(by_alias=True), "0.1.0.dev0")
    assert k_ab == k_ba
    # and the block is byte-identical
    blk_ab = _run(b_ab, tmp_path)["results"]["skew_group_hh"]
    blk_ba = _run(b_ba, tmp_path)["results"]["skew_group_hh"]
    assert json.dumps(blk_ab, sort_keys=True) == json.dumps(blk_ba, sort_keys=True)
    assert blk_ab["dim"] == 18 and blk_ab["dims"] == [2]


def test_modular_case_direct_only(tmp_path):
    # Z/2 vertex-swap on the 2-cycle Nakayama over GF(2): char 2 | |G| = 2.
    params = {"vertices": [1, 2], "arrows": {"a": [1, 2], "b": [2, 1]},
              "relations": ["a*b", "b*a"],
              "generators": [{"vertex_perm": {"1": 2, "2": 1},
                              "arrow_perm": {"a": "b", "b": "a"}}]}
    body = _family(params, {"kind": "GF", "p": 2}, ["skew_group_hh:2"])
    block = _run(body, tmp_path)["results"]["skew_group_hh"]
    assert block["char_ok"] is False
    assert block["status"] == "modular"
    assert block["decomposition"] is None
    assert block["dims"] == block["direct_dims"]     # direct HH of A|xG only


def test_reproduce_snippet_execs_and_twin_matches():
    # The SkewGroupAlgebra reproduce snippet is runnable Python that rebuilds A|xG,
    # and both runners emit byte-identical construction lines. (SkewGroupAlgebra is a
    # STRUCTURE-CONSTANT algebra with no quiver, so it is NOT in the quiver-presented
    # _FAMILIES parametrization of test_construction_families_input.py -- covered here.)
    import importlib.util
    import pathlib
    from quiverlab.hpc.spec import _synthetic_reproduce_lines
    params = dict(_FLAGSHIP_PARAMS)
    for ref in ("ql.GF(5)", "ql.CC", "QQ"):
        server = _synthetic_reproduce_lines("SkewGroupAlgebra", params, ref)
        rp = pathlib.Path(__file__).resolve().parents[1] / "gui" / "test_skew_group_runner_twin_p74.py"
        runner_path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
        spec = importlib.util.spec_from_file_location("gui_runner_repro", runner_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert server == mod._synthetic_reproduce_lines("SkewGroupAlgebra", params, ref)
    # the GF(5) form execs to the right dimension (2*2 = 4)
    ns: dict = {}
    exec("import quiverlab as ql\n" + "\n".join(
        _synthetic_reproduce_lines("SkewGroupAlgebra", params, "ql.GF(5)")), ns)
    assert ns["A"].dim == 4


def test_bad_action_is_clean_error(tmp_path):
    # x |-> y swap on an ASYMMETRIC base is not an automorphism -> clean 4xx-style
    # CatalogError at build, never a 500.
    params = {"vertices": [1], "arrows": {"x": [1, 1], "y": [1, 1]},
              "relations": ["x*x", "x*y", "y*x", "y*y*y"],
              "generators": [{"vertex_perm": {"1": 1}, "arrow_perm": {"x": "y", "y": "x"}}]}
    body = _family(params, {"kind": "QQ"}, ["skew_group_hh:2"])
    try:
        _run(body, tmp_path)
    except Exception as exc:                    # ComputeError / CatalogError family build
        assert "automorph" in str(exc).lower() or "relation" in str(exc).lower() \
            or "bijective" in str(exc).lower()
        return
    raise AssertionError("expected a loud refusal for the non-automorphism action")
