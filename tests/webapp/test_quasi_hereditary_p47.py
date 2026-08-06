"""The quasi_hereditary algebra-scalar kind (Plan 47): schema-1, served by hpc.spec,
mirrored by the Pyodide twin, both byte-identical on the block."""
import json
import tempfile

from quiverlab.hpc.spec import parse_request, run


def _qh_request():
    return {
        "schema": 1,
        "algebra": {
            "kind": "quiver",
            "vertices": [1, 2, 3],
            "arrows": {"a1": [1, 2], "a2": [2, 3]},
            "relations": [],
            "field": {"kind": "GF", "p": 32003, "n": 1},
        },
        "compute": ["quasi_hereditary"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def test_quasi_hereditary_block_shape(tmp_path):
    out = run(parse_request(_qh_request()), tmp_path)
    b = out["results"]["quasi_hereditary"]
    assert b["is_quasi_hereditary"] is True
    assert b["gl_dim"] == {"value": 1, "exact": True}
    assert set(b["standard_dims"]) == {"1", "2", "3"}
    assert "DlabRingel1989" in [k for k, _ in b["citations"]]   # resolved bibtex key
    assert "order-dependent" in b["order_note"]


def test_twin_parity():
    import importlib.util
    import os
    server = run(parse_request(_qh_request()), tempfile.mkdtemp())
    sb = server["results"]["quasi_hereditary"]
    # run the same request through the Pyodide twin docs/gui/runner.py
    runner_path = os.path.join(os.getcwd(), "docs/gui/runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p47", runner_path)
    gr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gr)
    built = json.loads(gr.run_build(json.dumps(_qh_request())))
    assert built["ok"], built
    tb = json.loads(gr.compute_one("quasi_hereditary"))["block"]
    assert json.dumps(sb, sort_keys=True) == json.dumps(tb, sort_keys=True)
