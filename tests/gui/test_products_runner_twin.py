"""docs/gui/runner.py serves the HH product kinds (cup / cap / bracket /
connes_b, Plan 35) SHAPE-IDENTICALLY to the server tier.

The Pyodide runner twin (``docs/gui/runner.py``) and the wheel's spec core
(``quiverlab.hpc.spec``, exercised here through ``webapp.server.runner.run_spec``)
duplicate the algebra dispatch by design -- they cannot import each other. The
contract is that a product block one runner emits is key-for-key equal to the
other's. This is the deliverable: the twin's block == the server's block.

The runner twin is loaded exactly as ``tests/gui/conftest.py`` loads it (a fresh
module, the same file the browser gets). It is driven by ``run_build`` +
``compute_one`` (there is no ``run_request``); each compute returns ``{"block": ...}``
whose value we compare to the server's ``results[kind]``.
"""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")


def _load_gui_runner():
    spec = importlib.util.spec_from_file_location("gui_runner", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# k<x>/(x^3) over GF(2): a single-loop truncated polynomial algebra -- all four
# product kinds are defined on it (bracket is GF(p)-only + window-bounded).
_BODY = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1],
                     "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                     "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": ["cup:0..2", "cap:0..2", "bracket:0..2", "connes_b:0..2"],
         "artifacts": {"pdf": False, "tikz": False}}
_KINDS = ("cup", "cap", "bracket", "connes_b")


def _gui_blocks(body):
    gui = _load_gui_runner()
    built = json.loads(gui.run_build(json.dumps(body)))
    assert built["ok"], built
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


def test_product_blocks_match_spec_runner(tmp_path):
    gui = _gui_blocks(_BODY)
    ref = _server_blocks(_BODY, tmp_path)
    for kind in _KINDS:
        assert gui[kind] == ref[kind], kind


def test_both_runners_carry_explicit_reps(tmp_path):
    """Plan 35 explicit representatives: both runners emit basis_classes /
    chain_basis / differentials for every product kind (the key-for-key equality
    above passes even if BOTH dropped them, so pin their PRESENCE here)."""
    gui = _gui_blocks(_BODY)
    ref = _server_blocks(_BODY, tmp_path)
    for kind in _KINDS:
        for field in ("basis_classes", "chain_basis", "differentials"):
            assert field in gui[kind], (kind, field, "gui")
            assert field in ref[kind], (kind, field, "server")
        # every product block carries a homology or cohomology side (never empty)
        assert gui[kind]["basis_classes"], kind


def test_product_block_carries_references_and_citations(tmp_path):
    # The twin's product block keeps `references` (from `.blocks()`) AND the
    # resolved `citations` pairs, exactly like every other invariant block --
    # this is what the key-for-key equality above pins, spelled out.
    gui = _gui_blocks(_BODY)
    for kind in _KINDS:
        block = gui[kind]
        assert block["references"], kind
        assert block["citations"] and len(block["citations"][0]) == 2, kind


def test_cyclic_homology_block_matches_spec_runner(tmp_path):
    # Plan-35 follow-up: the cyclic-homology block the twin emits is key-for-key
    # equal to the server tier's (same kind/top/dims/engine/references/citations).
    body = dict(_BODY, compute=["cyclic_homology:0..4"])
    gui = _gui_blocks(body)
    ref = _server_blocks(body, tmp_path)
    # the GUI keys blocks by table.kind ("HC_"); the server keys by compute kind.
    assert gui["HC_"] == ref["cyclic_homology"]
    assert gui["HC_"]["references"] == ["cyclic"] and gui["HC_"]["citations"]


def test_python_snippet_covers_product_kinds():
    # A product compute must have a snippet entry, or python_snippet() KeyErrors
    # (the same failure mode the `dimension` correction fixed).
    gui = _load_gui_runner()
    assert json.loads(gui.run_build(json.dumps(_BODY)))["ok"]
    for spec in _BODY["compute"]:
        assert json.loads(gui.compute_one(spec))["ok"]
    snippet = gui.python_snippet()
    for call in ("A.cup_products(2)", "A.cap_products(2)",
                 "A.gerstenhaber_brackets(2)", "A.connes_differentials(2)"):
        assert call in snippet, call
