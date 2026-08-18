"""Plan 79 Task G1 -- the browser half of the `cluster_category` kind: the Pyodide twin's
OWN dispatch, byte-parity with the server tier, and the gui.js wiring in BOTH copies.

Both runners route through the single shared builder
`cluster.category.cluster_category_block`, so the blocks are byte-identical by
construction; this file proves it by driving the twin's own run_build + compute_one.

Unmarked per the Plan-32 extras-gated ruling.
"""
import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "docs" / "gui" / "runner.py"
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEB = ROOT / "webapp" / "static" / "gui" / "gui.js"

_A3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a": [1, 2], "b": [2, 3]},
       "relations": [], "field": {"kind": "QQ"}}


def _runner():
    spec = importlib.util.spec_from_file_location("gui_runner_p79", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _twin_block(compute="cluster_category", algebra=None):
    m = _runner()
    body = {"schema": 1, "algebra": algebra or _A3, "compute": [compute],
            "artifacts": {"pdf": False, "tikz": False}}
    built = json.loads(m.run_build(json.dumps(body)))
    assert built["ok"], built
    out = json.loads(m.compute_one(compute))
    assert out["ok"], out
    return out["block"]


def _server_block(compute="cluster_category", algebra=None):
    import tempfile

    from quiverlab.hpc.spec import parse_request
    from quiverlab.hpc.spec import run as spec_run
    body = {"schema": 1, "algebra": algebra or _A3, "compute": [compute],
            "artifacts": {"pdf": False, "tikz": False}}
    with tempfile.TemporaryDirectory() as td:
        return spec_run(parse_request(body), pathlib.Path(td))["results"]["cluster_category"]


def test_the_twin_serves_the_kind_through_its_own_dispatch():
    b = _twin_block()
    assert b["kind"] == "cluster_category"
    assert b["num_indec"] == 9
    assert b["num_cluster_tilting"]["count"] == 14


def test_twin_parity_is_byte_identical_to_the_server():
    twin, server = _twin_block(), _server_block()
    assert (json.dumps(twin, sort_keys=True, default=str)
            == json.dumps(server, sort_keys=True, default=str))


def test_twin_parity_with_an_explicit_budget():
    twin = _twin_block("cluster_category:64")
    server = _server_block("cluster_category:64")
    assert twin["budget_pairs"] == 64
    assert (json.dumps(twin, sort_keys=True, default=str)
            == json.dumps(server, sort_keys=True, default=str))


def test_the_twin_parses_the_budget_grammar_like_the_servers():
    m = _runner()
    assert m._parse_compute("cluster_category") == ("cluster_category", None)
    assert m._parse_compute("cluster_category:64") == ("cluster_category", 64)


def test_the_twin_eta_puts_the_kind_above_the_module_level_kinds():
    m = _runner()
    assert m.ETA_MODEL["scalars"]["cluster_category"] > m.ETA_MODEL["scalars"]["koszul"]


def test_both_gui_copies_are_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEB.read_bytes()


def test_gui_wires_every_cluster_category_touchpoint():
    js = GUI_DOCS.read_text(encoding="utf-8")
    assert 'id="qlgui-cluster_category"' in js
    assert 'id="qlgui-cluster_category-budget"' in js
    assert '"cluster_category", "cluster_category-budget",' in js
    assert 'compute.push("cluster_category:" + el["cluster_category-budget"].value)' in js
    assert 'name === "cluster_category"' in js
    assert 'cluster_category: { cb: "cluster_category", top: "cluster_category-budget"' in js
    assert '"wall_chamber", "cluster_category", "silting"' in js


def test_the_gui_renderer_never_prints_an_uncertified_count_as_a_number():
    js = GUI_DOCS.read_text(encoding="utf-8")
    start = js.index('name === "cluster_category"')
    branch = js[start:start + 3000]
    assert "ct.certified" in branch                     # the count is gated on it
    assert "NOT certified" in branch
    assert "cited, not computed" in branch              # the 2-CY scope is stated
    assert "b.error" in branch


def test_i18n_has_the_picker_label_in_all_four_locales():
    for lang in ("en", "es", "fr", "zh"):
        cat = json.loads(
            (ROOT / "webapp" / "server" / "i18n" / f"{lang}.json").read_text(
                encoding="utf-8"))
        assert cat["pick.kind.cluster_category"].strip()
    draw = (ROOT / "webapp" / "templates" / "draw.html").read_text(encoding="utf-8")
    assert "data-pick-kind-cluster_category" in draw
