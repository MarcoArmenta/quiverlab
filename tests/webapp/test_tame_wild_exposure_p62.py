"""tame_wild kind (Plan 62 / R19): served by hpc.spec, mirrored by the Pyodide twin,
the rep-finite / tame / wild verdict rendered honestly, both runners byte-identical
on the block."""
import importlib.util
import json
import os
import tempfile

from quiverlab.hpc.spec import parse_request, run


# The verdict layer is algebraically-closed-only (CC, or -- absent P61's field flag --
# any char-0 field by base change): CC has a verdict in BOTH regimes, so verdict-asserting
# requests use CC (over QQ post-P61 the flag refuses the verdict). GF(p) = the refusal.
def _linear(n, compute, field=None):
    field = field or {"kind": "CC"}
    return {
        "schema": 1,
        "algebra": {
            "kind": "quiver",
            "vertices": list(range(1, n + 1)),
            "arrows": {f"a{i}": [i, i + 1] for i in range(1, n)},
            "relations": [],
            "field": field,
        },
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False},
    }


def _star(arms, compute, field=None):
    field = field or {"kind": "CC"}
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt)
            ars[f"e{ai}_{nxt}"] = [prev, nxt]
            prev = nxt
            nxt += 1
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": verts, "arrows": ars,
                    "relations": [], "field": field},
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False},
    }


def _kron(m, compute, field=None):
    field = field or {"kind": "CC"}
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                    "arrows": {f"a{k}": [1, 2] for k in range(m)},
                    "relations": [], "field": field},
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False},
    }


def _twin_block(req, kind):
    runner_path = os.path.join(os.getcwd(), "docs/gui/runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p62", runner_path)
    gr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gr)
    built = json.loads(gr.run_build(json.dumps(req)))
    assert built["ok"], built
    return json.loads(gr.compute_one(kind))["block"]


def test_tame_wild_block_shape(tmp_path):
    out = run(parse_request(_linear(5, ["tame_wild"])), tmp_path)   # A5 tree
    b = out["results"]["tame_wild"]
    assert b["rep_type"] == "rep-finite"
    assert b["weakly_positive"] is True and b["is_unit_form"] is True
    assert b["strongly_simply_connected"] is True
    assert "bongartz_criterion" in b["references"]
    assert "Bongartz1984" in [k for k, _ in b["citations"]]         # resolved bibtex key
    assert "scope_note" in b and b["gram"]


def test_tame_wild_wild_witness(tmp_path):
    out = run(parse_request(_star([1, 2, 6], ["tame_wild"])), tempfile.mkdtemp())  # T237
    b = out["results"]["tame_wild"]
    assert b["rep_type"] == "wild"
    assert b["weakly_nonnegative"] is False
    assert b["witness"] is not None and b["witness_value"] < 0


def test_tame_wild_off_scope_keeps_form(tmp_path):
    # 2-Kronecker over CC: form computed (weakly nonnegative True), but NOT simply
    # connected, so the verdict is withheld (rep_type null) -- the honest layer split
    # (over CC the FIELD gate passes, so the withholding is purely the P56 gate).
    out = run(parse_request(_kron(2, ["tame_wild"])), tempfile.mkdtemp())
    b = out["results"]["tame_wild"]
    assert b["weakly_nonnegative"] is True and b["rep_type"] is None
    assert b["simply_connected"] is False and "simply connected" in b["reason"]


def test_tame_wild_non_CC_keeps_form(tmp_path):
    # GF(7): the form is field-free (computed), the verdict refused (not alg. closed).
    out = run(parse_request(_linear(5, ["tame_wild"], {"kind": "GF", "p": 7, "n": 1})),
              tempfile.mkdtemp())
    b = out["results"]["tame_wild"]
    assert b["weakly_positive"] is True and b["field_alg_closed"] is False
    assert b["rep_type"] is None


def test_twin_parity(tmp_path):
    reqs = [
        (_linear(5, ["tame_wild"]), "tame_wild"),        # rep-finite
        (_kron(2, ["tame_wild"]), "tame_wild"),          # off-scope (form only)
        (_star([1, 2, 6], ["tame_wild"]), "tame_wild"),  # wild
    ]
    for req, kind in reqs:
        sb = run(parse_request(req), tempfile.mkdtemp())["results"][kind]
        tb = _twin_block(req, kind)
        assert json.dumps(sb, sort_keys=True) == json.dumps(tb, sort_keys=True)
