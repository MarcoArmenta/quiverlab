"""The tame_wild kind through the spec runner (Plan 62 / R19).

``quiverlab.hpc.spec`` exposes ``parse_request`` + ``run``; this drives the real
``_dispatch`` branch. The webapp runner delegates to the same dispatch (pinned
byte-stable by ``tame_wild_a5_cc`` in ``tests/webapp/_runner_goldens.json``), and the
Pyodide twin (``docs/gui/runner.py``) is checked byte-identical in
``tests/webapp/test_tame_wild_exposure_p62.py::test_twin_parity``.
"""
from quiverlab.hpc.spec import parse_request
from quiverlab.hpc.spec import run as spec_run


def _linear(n, field=None):
    return {"kind": "quiver", "vertices": list(range(1, n + 1)),
            "arrows": {f"a{i}": [i, i + 1] for i in range(1, n)},
            "relations": [], "field": field or {"kind": "QQ"}}


def _star(arms, field=None):
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt)
            ars[f"e{ai}_{nxt}"] = [prev, nxt]
            prev = nxt
            nxt += 1
    return {"kind": "quiver", "vertices": verts, "arrows": ars,
            "relations": [], "field": field or {"kind": "QQ"}}


def _req(algebra, compute):
    return parse_request({"schema": 1, "algebra": algebra, "compute": compute,
                          "artifacts": {"pdf": False, "tikz": False}})


def test_tame_wild_dynkin_rep_finite(tmp_path):
    b = spec_run(_req(_linear(5), ["tame_wild"]), tmp_path)["results"]["tame_wild"]
    assert b["rep_type"] == "rep-finite"
    assert b["weakly_positive"] is True and b["strongly_simply_connected"] is True
    assert b["references"][0] == "bongartz_criterion" and b["citations"]
    assert b["is_unit_form"] is True and b["field_alg_closed"] is True


def test_tame_wild_euclidean_tame(tmp_path):
    b = spec_run(_req(_star([1, 2, 5]), ["tame_wild"]), tmp_path)["results"]["tame_wild"]
    assert b["rep_type"] == "tame"
    assert b["weakly_nonnegative"] is True and b["weakly_positive"] is False
    assert b["witness"] is not None                       # isotropic radical direction


def test_tame_wild_T237_wild_witness(tmp_path):
    b = spec_run(_req(_star([1, 2, 6]), ["tame_wild"]), tmp_path)["results"]["tame_wild"]
    assert b["rep_type"] == "wild"
    assert b["weakly_nonnegative"] is False and b["witness_value"] < 0
    assert len(b["gram"]) == 10                            # the 10-variable form


def test_tame_wild_non_triangular_is_error_block(tmp_path):
    # k[x]/(x^3): a loop -> no unit form -> {"error": ...}, never a crash / 500.
    loop = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
            "relations": ["x*x*x"], "field": {"kind": "QQ"}}
    b = spec_run(_req(loop, ["tame_wild"]), tmp_path)["results"]["tame_wild"]
    assert "error" in b and ("triangular" in b["error"] or "cycle" in b["error"])
    assert "bongartz_criterion" in b["references"]


def test_reproduce_snippet_names_the_method(tmp_path):
    out = spec_run(_req(_linear(5), ["tame_wild"]), tmp_path)
    assert "tame_wild_certificate" in out["reproduce"]
