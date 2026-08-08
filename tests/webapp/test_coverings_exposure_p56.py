"""fundamental_group + simply_connected kinds (Plan 56): served by hpc.spec,
mirrored by the Pyodide twin, three-valued verdict rendered honestly, both runners
byte-identical on the block."""
import importlib.util
import json
import os
import tempfile

from quiverlab.hpc.spec import parse_request, run


def _square(rel, compute):
    return {
        "schema": 1,
        "algebra": {
            "kind": "quiver",
            "vertices": [1, 2, 3, 4],
            "arrows": {"a": [1, 2], "b": [2, 4], "c": [1, 3], "d": [3, 4]},
            "relations": rel,
            "field": {"kind": "GF", "p": 7, "n": 1},
        },
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False},
    }


def _zito(compute):
    return {
        "schema": 1,
        "algebra": {
            "kind": "quiver",
            "vertices": [1, 2, 3, 4, 5],
            "arrows": {"al": [1, 2], "be": [2, 3], "ga": [3, 5],
                       "de": [2, 4], "ep": [4, 5]},
            "relations": ["al*be*ga - al*de*ep"],
            "field": {"kind": "QQ"},
        },
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False},
    }


def _twin_block(req, kind):
    runner_path = os.path.join(os.getcwd(), "docs/gui/runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p56", runner_path)
    gr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gr)
    built = json.loads(gr.run_build(json.dumps(req)))
    assert built["ok"], built
    return json.loads(gr.compute_one(kind))["block"]


def test_fundamental_group_block_shape(tmp_path):
    out = run(parse_request(_square(["a*b - c*d"], ["fundamental_group"])), tmp_path)
    b = out["results"]["fundamental_group"]
    assert b["abelianization"]["free_rank"] == 0
    assert b["abelianization"]["invariant_factors"] == []
    assert "assem_delapena" in b["references"]
    assert "AssemDelaPena1996" in [k for k, _ in b["citations"]]   # resolved bibtex key
    assert "intrinsic_note" in b                       # honest-decidability documentation


def test_simply_connected_three_valued(tmp_path):
    # square WITHOUT relation: False, witness kind nontrivial_pi1ab
    out = run(parse_request(_square([], ["simply_connected"])), tempfile.mkdtemp())
    b = out["results"]["simply_connected"]
    assert b["verdict"] is False and b["witness"]["kind"] == "nontrivial_pi1ab"
    # the Zito example: simply connected (True) but NOT strongly (strongly.verdict False)
    outz = run(parse_request(_zito(["simply_connected"])), tempfile.mkdtemp())
    bz = outz["results"]["simply_connected"]
    assert bz["verdict"] is True
    assert bz["strongly"]["verdict"] is False and bz["strongly"]["witness"]["vertex"] == 2


def test_twin_parity(tmp_path):
    reqs = [
        (_square(["a*b - c*d"], ["fundamental_group"]), "fundamental_group"),
        (_zito(["simply_connected"]), "simply_connected"),
    ]
    for req, kind in reqs:
        sb = run(parse_request(req), tempfile.mkdtemp())["results"][kind]
        tb = _twin_block(req, kind)
        assert json.dumps(sb, sort_keys=True) == json.dumps(tb, sort_keys=True)
