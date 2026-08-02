"""Plan 28 -- config-export tests (fast).

The webapp-side YAML helper (``webapp/server/clusterconfig.py``) is tested
thoroughly: its output round-trips through a YAML parser back into a schema-valid
compute request, carries the runnable header comment, and preserves module
matrices + exact fraction entries.

The ``docs/gui/gui.js`` ``configYaml()`` emitter MIRRORS the helper. It is a pure,
self-contained function, so we extract it textually and -- when node is available
-- run it on the same requests and assert the SAME round-trip. A pinning check
asserts both files carry the mirror comment so a future edit to one flags the
other.
"""
import json
import pathlib
import shutil
import subprocess

import pytest
import yaml

from webapp.server.clusterconfig import HEADER, cluster_config_yaml
from webapp.server.schema import ComputeRequest

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_GUI_JS = _ROOT / "docs" / "gui" / "gui.js"
_CLUSTERCONFIG = _ROOT / "webapp" / "server" / "clusterconfig.py"

# A representative set: a family request, a quiver request, and a no-code module
# request with a fraction entry + a left side (the hardest cases for the emitter).
REQUESTS = [
    {"schema": 1, "algebra": {"kind": "family", "family": "QuantumCI",
                              "params": {"q": 1}, "field": {"kind": "CC"}},
     "compute": ["hh_cohomology:0..4", "cartan"],
     "artifacts": {"pdf": True, "tikz": False}},
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1, 2],
                              "arrows": {"a": [1, 2]}, "relations": [],
                              "field": {"kind": "GF", "p": 3, "n": 1}},
     "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": True}},
    {"schema": 2, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                              "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
     "compute": ["dimension_vector", "projective_resolution:0..3"],
     "artifacts": {"pdf": False, "tikz": False},
     "module": {"dims": {"1": 2}, "maps": {"x": [["0", "1/2"], ["1", "0"]]},
                "side": "left"}},
]


def _dump(req):
    return ComputeRequest.model_validate(req).model_dump(by_alias=True)


# --------------------------------------------------------------------------- #
# Python helper: round-trip through YAML into a schema-valid request
# --------------------------------------------------------------------------- #

def test_helper_has_runnable_header():
    y = cluster_config_yaml(_dump(REQUESTS[0]))
    assert y.startswith(HEADER)
    assert "quiverlab-hpc run this-file.yaml -o result.json" in y


@pytest.mark.parametrize("req", REQUESTS)
def test_helper_round_trips_to_schema_valid_request(req):
    dump = _dump(req)
    loaded = yaml.safe_load(cluster_config_yaml(dump))     # YAML ignores the header comment
    revalidated = ComputeRequest.model_validate(loaded)    # must not raise
    assert revalidated.algebra.kind == req["algebra"]["kind"]
    assert revalidated.compute == req["compute"]


def test_helper_preserves_module_matrix_and_fraction():
    dump = _dump(REQUESTS[2])
    loaded = yaml.safe_load(cluster_config_yaml(dump))
    m = ComputeRequest.model_validate(loaded).module
    assert m.side == "left"
    assert m.dims == {"1": 2}
    # the fraction entry survives as an exact string, not a float
    assert m.maps["x"][0][1] == "1/2"


# --------------------------------------------------------------------------- #
# JS emitter mirror (node): same round-trip, so the two implementations agree
# --------------------------------------------------------------------------- #

def _extract_js_fn(name: str) -> str:
    """Textually extract a top-level ``function <name>(...) { ... }`` by matching
    braces. gui.js is an IIFE so the function is not exported at module scope --
    this is the same technique test_pages.py uses for ``isHttpUrl``."""
    src = _GUI_JS.read_text(encoding="utf-8")
    start = src.index(f"function {name}")
    depth, started, i = 0, False, start
    while i < len(src):
        c = src[i]
        if c == "{":
            depth += 1
            started = True
        elif c == "}":
            depth -= 1
            if started and depth == 0:
                return src[start:i + 1]
        i += 1
    raise AssertionError(f"function {name} not found / unbalanced in gui.js")


def test_gui_js_exposes_config_yaml():
    fn = _extract_js_fn("configYaml")
    assert "function keyScalar" in fn            # numeric dict keys stay strings
    src = _GUI_JS.read_text(encoding="utf-8")
    assert "window.QLGUI" in src and "configYaml: configYaml" in src


@pytest.mark.parametrize("req", REQUESTS)
def test_js_emitter_round_trips_like_the_helper(req):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")
    fn = _extract_js_fn("configYaml")
    harness = fn + "\nprocess.stdout.write(configYaml(" + json.dumps(req) + "));\n"
    proc = subprocess.run([node, "-e", harness], capture_output=True, text=True,
                          encoding="utf-8", timeout=30)
    assert proc.returncode == 0, proc.stderr
    loaded = yaml.safe_load(proc.stdout)         # JS-emitted YAML -> dict
    revalidated = ComputeRequest.model_validate(loaded)   # must be schema-valid
    assert revalidated.compute == req["compute"]
    if req.get("module"):
        assert revalidated.module.side == req["module"]["side"]
        assert revalidated.module.dims == req["module"]["dims"]


# --------------------------------------------------------------------------- #
# Pinning: the two emitters name each other so an edit to one flags the other
# --------------------------------------------------------------------------- #

def test_emitters_are_pinned_to_each_other():
    js = _GUI_JS.read_text(encoding="utf-8")
    py = _CLUSTERCONFIG.read_text(encoding="utf-8")
    assert "clusterconfig.py" in js, "gui.js must name the Python mirror"
    assert "gui.js" in py, "clusterconfig.py must name the JS mirror"
