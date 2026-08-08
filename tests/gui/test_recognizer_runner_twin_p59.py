"""The Pyodide GUI runner (docs/gui/runner.py) serves the Plan-59 `string_homological`
and `toupie` kinds SHAPE-IDENTICALLY to the server core (quiverlab.hpc.spec): the SAME
library block builders (modules.string_homological.string_homological_block /
families.toupie.toupie_block) + `references`->citations, so the recognizer verdicts,
HH row and sl_a line cannot drift between the browser GUI and the server report."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_KD4 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [0, 1, 2, 3],
                    "arrows": {"a": [1, 0], "b": [2, 0], "c": [3, 0]},
                    "relations": [], "field": {"kind": "QQ"}},
        "compute": ["string_homological"], "artifacts": {"pdf": False, "tikz": False}}
_KRON3 = {"schema": 1,
          "algebra": {"kind": "quiver", "vertices": [1, 2],
                      "arrows": {"a": [1, 2], "b": [1, 2], "c": [1, 2]},
                      "relations": [], "field": {"kind": "QQ"}},
          "compute": ["toupie"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, kind):
    spec = importlib.util.spec_from_file_location("gui_runner_recognizer_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(kind))
    assert out["ok"], out
    return out["block"]


def test_gui_string_homological_block_shape():
    b = _gui_block(_KD4, "string_homological")
    assert b["verdict"] == "not_string" and b["is_string"] is False
    assert b["witness"]["summand_count"] >= 3
    assert "suarez_alvarez" in b["references"]


def test_gui_toupie_block_shape():
    b = _gui_block(_KRON3, "toupie")
    assert b["is_toupie"] is True and b["direct_arrow_count"] == 3
    assert b["hh"][:2] == [1, 8] and b["sl_a"]["dim"] == 8
    assert "alsolotar_toupie" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    for body, kind in ((_KD4, "string_homological"), (_KRON3, "toupie")):
        server = run_spec(ComputeRequest.model_validate(body), tmp_path)["results"][kind]
        gui = _gui_block(body, kind)
        assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
