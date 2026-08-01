"""docs/gui/runner.py serves the module Ext / Tor kinds SHAPE-IDENTICALLY to the server
tier, INCLUDING the Plan 35 wave-3a explicit representatives (basis_classes /
chain_basis / differentials). The Pyodide twin and the wheel's spec core cannot import
each other, so the contract is that the block one emits is key-for-key equal to the
other's -- and both carry the reps.
"""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
         "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}},
         "ext_target": {"builtin": {"kind": "simple", "vertex": 1}},
         # tor_target side omitted -> both tiers default it to LEFT (the Tor convention)
         "tor_target": {"builtin": {"kind": "simple", "vertex": 1}},
         "compute": ["ext:0..3", "tor:0..3"],
         "artifacts": {"pdf": False, "tikz": False}}
_KINDS = ("ext", "tor")


def _load_gui_runner():
    spec = importlib.util.spec_from_file_location("gui_runner", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _gui_blocks(body):
    gui = _load_gui_runner()
    assert json.loads(gui.run_build(json.dumps(body)))["ok"]
    blocks = {}
    for spec in body["compute"]:
        out = json.loads(gui.compute_one(spec))
        assert out["ok"], out
        blocks[out["block"]["kind"]] = out["block"]
    return blocks


def _server_blocks(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]


# The wave-3a fields (plus the dims/target they annotate) must be key-for-key equal
# across the two runners. The whole-block equality is NOT asserted: module blocks carry
# a pre-existing tier difference (the server's `_with_refs` adds a `references` key the
# Pyodide runner does not) that predates and is orthogonal to this wave.
_COMPARE = ("kind", "top", "dims", "target", "basis_classes", "chain_basis",
            "differentials")


def test_ext_tor_blocks_match_spec_runner(tmp_path):
    gui = _gui_blocks(_BODY)
    ref = _server_blocks(_BODY, tmp_path)
    for kind in _KINDS:
        for key in _COMPARE:
            assert gui[kind].get(key) == ref[kind].get(key), (kind, key)


def test_both_runners_carry_explicit_reps(tmp_path):
    """Both runners emit basis_classes / chain_basis / differentials for Ext and Tor
    (the key-for-key equality passes even if BOTH dropped them, so pin PRESENCE)."""
    gui = _gui_blocks(_BODY)
    ref = _server_blocks(_BODY, tmp_path)
    for kind in _KINDS:
        for field in ("basis_classes", "chain_basis", "differentials"):
            assert field in gui[kind], (kind, field, "gui")
            assert field in ref[kind], (kind, field, "server")
        # the degree keys are the string range 0..3
        assert sorted(gui[kind]["basis_classes"]) == ["0", "1", "2", "3"], kind
