"""The four product kinds through the spec runner (Plan 35).

``quiverlab.hpc.spec`` exposes the dict->request validator as ``parse_request``
and the runner as ``run`` (there is no ``parse_config``); this drives the real
``_dispatch`` branch for cup/cap/bracket/connes_b. The webapp runner
(``webapp.server.runner.run_spec``) delegates to this same dispatch -- the
byte-stable golden ``products_loop_gf2`` in ``tests/webapp/_runner_goldens.json``
pins that path.
"""
import pytest

from quiverlab.hpc.spec import run as spec_run, parse_request as parse_config


def _req(compute):
    return parse_config({
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1],
                    "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                    "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False}})


def test_cup_block(tmp_path):
    res = spec_run(_req(["cup:0..2"]), tmp_path)
    b = res["results"]["cup"]
    assert b["kind"] == "cup" and b["top"] == 2
    assert b["basis"].startswith(("bar/", "cs/"))
    assert b["citations"] and b["references"]
    degs = [t["degrees"] for t in b["tables"]]
    assert [0, 0] in degs and [1, 1] in degs


def test_all_four_kinds_serve(tmp_path):
    res = spec_run(_req(["cup:0..2", "cap:0..2", "bracket:0..2",
                         "connes_b:0..2"]), tmp_path)
    assert set(res["results"]) == {"cup", "cap", "bracket", "connes_b"}
    assert res["results"]["bracket"]["window"] == 2
    assert res["results"]["connes_b"]["ranks"].keys() == {"0", "1"}


def test_range_required(tmp_path):
    from quiverlab.hpc.spec import ComputeError
    with pytest.raises(Exception):
        spec_run(_req(["cup"]), tmp_path)     # ComputeError -> typed 4xx upstream
