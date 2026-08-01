"""Marco 2026-07-31 (ADDENDUM 2): the ``json_guide`` documents, per computation, HOW to
recover every computed object from ``result.json``. Two things are pinned here:

  * THE GUIDE CANNOT LIE. Every recipe's ``path`` is WALKED against the actual result
    dict with a trivial dot / ``["key"]`` / ``[int]`` parser and must resolve to an
    existing object (the builder self-validates, this re-checks independently across the
    kinds an example actually computes).
  * The report appendix tolerates a pre-guide cache (no ``json_guide``) -- it renders
    the "Reading the JSON record" section only when a non-empty guide is passed.
"""
import re
import tempfile

from quiverlab.trace.json_guide import build_json_guide
from quiverlab.trace.render_html import render_html
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

# dot / ["str"] / [int] steps only (the guide's advertised, trivially parseable syntax).
_STEP = re.compile(r'\.([A-Za-z_][A-Za-z0-9_]*)|\["([^"]*)"\]|\[(\d+)\]')


def _resolve(result, path):
    """Walk ``path`` (rooted at ``results.<kind>...``) against the result dict; raises
    KeyError / IndexError / TypeError if any step is missing (i.e. the guide lied)."""
    assert path.startswith("results"), path
    obj = result["results"]
    for m in _STEP.finditer(path[len("results"):]):
        name, skey, idx = m.groups()
        key = name if name is not None else (skey if skey is not None else int(idx))
        obj = obj[key]
    return obj


def _run(compute, extra=None):
    body = {
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                    "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
                    "field": {"kind": "GF", "p": 7, "n": 1}},
        "module": {"builtin": {"kind": "simple", "vertex": 1, "side": "right"}},
        "compute": compute,
    }
    body.update(extra or {})
    req = ComputeRequest.model_validate(body)
    with tempfile.TemporaryDirectory() as td:
        return run_spec(req, td)


def _loop_run(compute, extra=None):
    """The x^3 loop over GF(2): the small example that carries products / cyclic."""
    body = {
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                    "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": compute,
    }
    body.update(extra or {})
    req = ComputeRequest.model_validate(body)
    with tempfile.TemporaryDirectory() as td:
        return run_spec(req, td)
def test_guide_paths_all_resolve_module_and_algebra_kinds():
    res = _run(["hh_cohomology:0..2", "cartan", "global_dimension",
                "ext:0..2", "projective_resolution:0..3", "rad_top_soc",
                "dimension_vector", "tau"],
               {"ext_target": {"builtin": {"kind": "simple", "vertex": 3,
                                           "side": "right"}}})
    guide = res["json_guide"]
    assert guide, "the guide must not be empty for a real computation"
    for g in guide:                       # every path resolves -> the guide cannot lie
        _resolve(res, g["path"])
        assert set(g) == {"object", "path", "note"}
    covered = {g["path"].split(".")[1].split("[")[0] for g in guide}
    for kind in ("hh_cohomology", "cartan", "ext", "projective_resolution",
                 "global_dimension", "rad_top_soc"):
        assert kind in covered, (kind, sorted(covered))
def test_guide_paths_all_resolve_products_cyclic_tor():
    res = _loop_run(["cup:0..2", "connes_b:0..2", "cyclic_homology:0..3", "tor:0..2"],
                    {"module": {"builtin": {"kind": "simple", "vertex": 1,
                                            "side": "right"}},
                     "tor_target": {"builtin": {"kind": "simple", "vertex": 1,
                                                "side": "left"}}})
    guide = res["json_guide"]
    assert guide
    for g in guide:
        _resolve(res, g["path"])
    covered = {g["path"].split(".")[1].split("[")[0] for g in guide}
    for kind in ("cup", "connes_b", "cyclic_homology", "tor"):
        assert kind in covered, (kind, sorted(covered))
    # a product recipe names the table lookup + the degree-keyed reps under a side
    assert any(g["path"] == "results.cup.tables" for g in guide)
    assert any(g["path"].startswith('results.cup.basis_classes["') for g in guide)
    assert any(g["path"] == 'results.connes_b.matrices["0"]' for g in guide)
def test_build_json_guide_accepts_runner_list_shape():
    """The Pyodide runner accumulates a LIST of ``invariant``-tagged blocks; the guide
    reads it via results_html.normalize, same as the {kind: block} server shape."""
    res = _run(["hh_cohomology:0..1", "cartan"])
    list_form = [dict(blk, invariant=kind) for kind, blk in res["results"].items()]
    guide = build_json_guide(list_form)
    assert guide
    assert {g["path"].split(".")[1].split("[")[0] for g in guide} >= {"hh_cohomology",
                                                                       "cartan"}
def test_guide_skips_error_blocks():
    guide = build_json_guide({"ext": {"error": {"type": "DepthLimitError",
                                                "message": "x"}}})
    assert guide == []


def test_appendix_renders_only_with_a_nonempty_guide():
    """Renderer tolerance: no guide (a pre-guide cache) -> no appendix; a real guide ->
    the 'Reading the JSON record' section with the object/path/note table."""
    events = []
    assert "Reading the JSON record" not in render_html(events, title="t")
    assert "Reading the JSON record" not in render_html(events, title="t", json_guide=[])
    guide = [{"object": "Ext dimensions", "path": "results.ext.dims",
              "note": "list; dim in degree n is dims[n]."}]
    html = render_html(events, title="t", json_guide=guide)
    assert "Reading the JSON record" in html
    assert "results.ext.dims" in html
