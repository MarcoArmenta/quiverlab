"""Random representations (GitHub #3, Samuel Leblanc).

Covers the wheel core ``quiverlab.hpc.spec.random_module``, its byte-identical
Pyodide twin ``docs/gui/runner.py::random_module``, and the webapp route
``POST /api/gui/random-module``. Randomness is EXACT (GF(p): prime-subfield ints;
char 0: small ints), relation-aware (rejection sampling), and refuses char-0-with-
relations up front. Determinism in the seed is the contract that lets both runners
agree, so the parity test is the linchpin.
"""
import importlib.util
import json
import pathlib

import pytest

from quiverlab.hpc import spec
from quiverlab.hpc.spec import ModuleSpec, _build_module, build_algebra

_ROOT = pathlib.Path(__file__).resolve().parents[2]


def _load_gui_runner():
    """Import docs/gui/runner.py (the Pyodide twin) as a module."""
    path = _ROOT / "docs" / "gui" / "runner.py"
    sp = importlib.util.spec_from_file_location("gui_runner_twin", path)
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


_RUNNER = _load_gui_runner()

_KA3 = {"kind": "quiver", "vertices": [1, 2, 3],
        "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
        "field": {"kind": "GF", "p": 7, "n": 1}}
_LOOP2 = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
          "relations": ["x*x"], "field": {"kind": "GF", "p": 5, "n": 1}}
_SQUARE_GF2 = {"kind": "quiver", "vertices": [1, 2, 3, 4],
               "arrows": {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4]},
               "relations": ["a*c - b*d"], "field": {"kind": "GF", "p": 2, "n": 1}}


# --------------------------------------------------------------------------- #
# (a) cross-runner parity: spec core == Pyodide twin, byte-for-byte, per seed.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("alg,dims,side,seed", [
    (_KA3, {"1": 1, "2": 1, "3": 1}, "right", 42),
    (_KA3, {"1": 2, "2": 1, "3": 1}, "right", 1),
    (_KA3, {"1": 1, "2": 1, "3": 1}, "left", 7),
    (_LOOP2, {"1": 2}, "right", 11),
    (_SQUARE_GF2, {"1": 1, "2": 1, "3": 1, "4": 1}, "right", 7),
])
def test_spec_and_runner_agree_byte_for_byte(alg, dims, side, seed):
    core = spec.random_module(alg, dims, side=side, seed=seed, tries=200)
    twin = json.loads(_RUNNER.random_module(json.dumps(
        {"algebra": alg, "dims": dims, "side": side, "seed": seed, "tries": 200})))
    assert json.dumps(core, sort_keys=True) == json.dumps(twin, sort_keys=True)


# --------------------------------------------------------------------------- #
# (b) no relations, any field: one valid draw with correctly shaped entries.
# --------------------------------------------------------------------------- #
def test_no_relations_gf_first_draw_valid():
    out = spec.random_module(_KA3, {"1": 1, "2": 1, "3": 1}, seed=3)
    assert out["tries"] == 1 and set(out["maps"]) == {"a", "b"}
    for mat in out["maps"].values():
        for row in mat:
            for x in row:
                assert isinstance(x, int) and 0 <= x < 7


def test_no_relations_char0_integer_entries():
    ccf = {"kind": "quiver", "vertices": [1, 2, 3],
           "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
           "field": {"kind": "CC"}}
    out = spec.random_module(ccf, {"1": 1, "2": 1, "3": 1}, seed=3)
    assert out["tries"] == 1
    for mat in out["maps"].values():
        for row in mat:
            for x in row:
                assert isinstance(x, int) and -5 <= x <= 5


def test_shapes_follow_side_and_dims():
    # right: block is dim[target] x dim[source]; a:1->2 with dims 2,3 -> 3x2.
    out = spec.random_module(_KA3, {"1": 2, "2": 3, "3": 1}, side="right", seed=0)
    assert len(out["maps"]["a"]) == 3 and len(out["maps"]["a"][0]) == 2
    assert len(out["maps"]["b"]) == 1 and len(out["maps"]["b"][0]) == 3


# --------------------------------------------------------------------------- #
# (c) relations + GF(p): a seeded draw satisfies the relations and BUILDS.
# --------------------------------------------------------------------------- #
def test_relations_gf2_success_builds_a_real_module():
    dims = {"1": 1, "2": 1, "3": 1, "4": 1}
    out = spec.random_module(_SQUARE_GF2, dims, seed=7, tries=200)
    assert "maps" in out and out["tries"] >= 1
    A = build_algebra(_SQUARE_GF2)
    M = _build_module(A, ModuleSpec(dims=dims, maps=out["maps"], side="right"), "M")
    assert M.dim == sum(dims.values())         # built with no relation violation


def test_relations_success_is_deterministic_in_seed():
    dims = {"1": 1, "2": 1, "3": 1, "4": 1}
    a = spec.random_module(_SQUARE_GF2, dims, seed=7, tries=200)
    b = spec.random_module(_SQUARE_GF2, dims, seed=7, tries=200)
    assert a == b


# --------------------------------------------------------------------------- #
# (d) char-0 WITH relations: principled refusal.
# --------------------------------------------------------------------------- #
def test_char0_with_relations_refuses():
    cc = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
          "relations": ["x*x"], "field": {"kind": "CC"}}
    assert spec.random_module(cc, {"1": 2}, seed=1) == {"error": "char0"}


# --------------------------------------------------------------------------- #
# (e) budget exhaustion.
# --------------------------------------------------------------------------- #
def test_budget_exhaustion_reports_tries():
    # tries=1 with a seed whose first draw violates the relation.
    out = spec.random_module(_SQUARE_GF2, {"1": 1, "2": 1, "3": 1, "4": 1},
                             seed=7, tries=1)
    assert out == {"error": "budget", "tries": 1}


# --------------------------------------------------------------------------- #
# (g) zero-dim vertex -> empty blocks, no crash.
# --------------------------------------------------------------------------- #
def test_zero_dim_vertex_gives_empty_blocks():
    out = spec.random_module(_KA3, {"1": 1, "2": 0, "3": 1}, seed=5)
    assert out["maps"]["a"] == []          # rows = dim[2] = 0
    assert out["maps"]["b"] == [[]]        # rows = dim[3] = 1, cols = dim[2] = 0


def test_oversize_dims_raise_module_too_large():
    from quiverlab.hpc.spec import ComputeError
    with pytest.raises(ComputeError):
        spec.random_module(_KA3, {"1": 9999}, seed=1)


# --------------------------------------------------------------------------- #
# (f) the webapp route.
# --------------------------------------------------------------------------- #
def _client(tmp_path):
    from fastapi.testclient import TestClient

    from webapp.server.app import create_app
    from webapp.server.config import Config
    return TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))


def test_route_happy_path_and_seed_echo(tmp_path):
    c = _client(tmp_path)
    r = c.post("/api/gui/random-module",
               json={"algebra": _KA3, "dims": {"1": 1, "2": 1, "3": 1}, "seed": 42})
    assert r.status_code == 200
    body = r.json()
    assert body["seed"] == 42 and set(body["maps"]) == {"a", "b"}


def test_route_defaults_a_seed_and_returns_it(tmp_path):
    c = _client(tmp_path)
    r = c.post("/api/gui/random-module",
               json={"algebra": _KA3, "dims": {"1": 1, "2": 1, "3": 1}})
    assert r.status_code == 200 and isinstance(r.json()["seed"], int)


def test_route_refusals_and_errors_are_4xx_not_500(tmp_path):
    c = _client(tmp_path)
    cc = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
          "relations": ["x*x"], "field": {"kind": "CC"}}
    assert c.post("/api/gui/random-module",
                  json={"algebra": cc, "dims": {"1": 2}}).status_code == 422
    assert c.post("/api/gui/random-module",
                  json={"algebra": _SQUARE_GF2,
                        "dims": {"1": 1, "2": 1, "3": 1, "4": 1},
                        "seed": 7, "tries": 1}).json()["error"] == "budget"
    assert c.post("/api/gui/random-module",
                  json={"algebra": _KA3, "dims": {"1": 9999}}).status_code == 422
    # an infinite-dimensional algebra is a clean 422, never a 500
    inf = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
           "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}
    assert c.post("/api/gui/random-module",
                  json={"algebra": inf, "dims": {"1": 1}}).status_code == 422
    assert c.post("/api/gui/random-module", json={"dims": {"1": 1}}).status_code == 400


# --------------------------------------------------------------------------- #
# (h) draw.html carries the data-attrs under every language mount.
# --------------------------------------------------------------------------- #
def test_draw_page_carries_feature_attrs_in_every_language(tmp_path):
    from webapp.server.i18n import LANGS
    c = _client(tmp_path)
    attrs = ["data-rel-rad2", "data-rel-comm", "data-rel-comm-cycle",
             "data-rel-generated", "data-rel-none", "data-mod-random",
             "data-mod-random-dims", "data-mod-random-done",
             "data-mod-random-fail", "data-mod-random-char0"]
    for lang in LANGS:
        prefix = "" if lang == "en" else "/" + lang
        html = c.get(prefix + "/draw").text
        for a in attrs:
            assert a + "=" in html, (lang, a)
