"""The split_extension + arrow_removal kinds through the spec runner (Plan 72).

``quiverlab.hpc.spec`` exposes ``parse_request`` + ``run``; this drives the real
``_dispatch`` branches. The webapp runner delegates to the same dispatch, and the
Pyodide twin (``docs/gui/runner.py``) is checked byte-identical in
``tests/gui/test_split_arrow_runner_twin.py``.
"""
from quiverlab.hpc.spec import parse_request
from quiverlab.hpc.spec import run as spec_run


def _req(algebra, compute):
    return parse_request({"schema": 1, "algebra": algebra, "compute": compute,
                          "artifacts": {"pdf": False, "tikz": False}})


def _kA2(field=None):
    return {"kind": "quiver", "vertices": [1, 2], "arrows": {"a1": [1, 2]},
            "relations": [], "field": field or {"kind": "GF", "p": 7, "n": 1}}


def _P1(field=None):
    return {"kind": "quiver", "vertices": [1, 2],
            "arrows": {"x": [1, 1], "c": [1, 2]}, "relations": ["x^2"],
            "field": field or {"kind": "GF", "p": 7, "n": 1}}


def test_split_extension_block(tmp_path):
    b = spec_run(_req(_kA2(), ["split_extension:4"]), tmp_path)["results"]["split_extension"]
    assert b["status"] == "complete"
    assert b["assembled"] == [3, 1, 1, 1, 1] == b["direct"]
    assert b["agrees"] is True and b["exact"] is True
    assert b["flank_M"] == [2, 1, 0, 1, 2, 1]          # M-flank to top+1
    assert b["hh1_nonzero"] is True
    assert b["references"] == ["cmrs_split", "crs_trivial_ext_hh1"] and b["citations"]


def test_arrow_removal_block(tmp_path):
    b = spec_run(_req(_P1(), ["arrow_removal:4"]), tmp_path)["results"]["arrow_removal"]
    assert b["status"] == "complete"
    assert b["removed"] == ["c"]
    assert b["hom_A"] == b["hom_B"] == [3, 1, 1, 1, 1]
    assert b["hom_agrees"] is True and b["hom_low_agrees"] == [True, True]
    assert b["coh_low_delta"][0] == -2                 # disconnection at n=0
    assert b["references"] == ["clms_arrow_removal", "han_conjecture"] and b["citations"]


def test_unpresentable_extension_is_clean_error_not_500(tmp_path):
    # Over CC, TrivialExtension falls back to a structure-constant build with no
    # quiver -> the split-extension LES refuses CLEANLY (status='unsupported' + error),
    # never a 500 (the scope guard turns the AttributeError into a typed refusal).
    ccq = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a1": [1, 2]},
           "relations": [], "field": {"kind": "CC"}}
    b = spec_run(_req(ccq, ["split_extension:4"]), tmp_path)["results"]["split_extension"]
    assert b["status"] == "unsupported" and "error" in b
    assert b["references"] == ["cmrs_split", "crs_trivial_ext_hh1"]


def test_reproduce_snippets_name_the_methods(tmp_path):
    out = spec_run(_req(_kA2(), ["split_extension:4"]), tmp_path)
    assert "split_extension_cohomology" in out["reproduce"]
    out2 = spec_run(_req(_P1(), ["arrow_removal:4"]), tmp_path)
    assert "arrow_removal" in out2["reproduce"]
