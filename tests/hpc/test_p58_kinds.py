"""The coxeter_spectral compute kind through the spec runner (Plan 58 / R20).

``quiverlab.hpc.spec`` exposes the dict->request validator ``parse_request`` and the
runner ``run``; this drives the real ``_dispatch`` branch for ``coxeter_spectral``.
The webapp runner delegates to this same dispatch (byte-stable golden
``coxeter_spectral_3kronecker_qq`` in ``tests/webapp/_runner_goldens.json``), and the
Pyodide twin (``docs/gui/runner.py``) is checked byte-identical in
``tests/webapp/test_coxeter_spectral_p58.py``.
"""
from quiverlab.hpc.spec import parse_request
from quiverlab.hpc.spec import run as spec_run

_KA2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
        "relations": [], "field": {"kind": "QQ"}}
_KRON3 = {"kind": "quiver", "vertices": [1, 2],
          "arrows": {"a0": [1, 2], "a1": [1, 2], "a2": [1, 2]},
          "relations": [], "field": {"kind": "QQ"}}


def _req(algebra, compute):
    return parse_request({"schema": 1, "algebra": algebra, "compute": compute,
                          "artifacts": {"pdf": False, "tikz": False}})


def test_coxeter_spectral_ka2(tmp_path):
    out = spec_run(_req(_KA2, ["coxeter_spectral"]), tmp_path)
    b = out["results"]["coxeter_spectral"]
    assert b["kind"] == "coxeter_spectral"
    assert b["cyclotomic"] is True and b["quasi_unipotent"] is True
    assert [f["cyclotomic_index"] for f in b["factorization"]] == [3]
    assert b["coxeter_order"] == 3
    assert b["outside_unit_circle_count"] == 0
    assert b["spectral_radius"]["value"] == "1" and b["mahler_measure"]["value"] == "1"
    assert "dlPena2014mahler" in b["references"] and b["citations"]
    assert "A.coxeter_spectral()" in out["reproduce"]


def test_coxeter_spectral_three_kronecker(tmp_path):
    b = spec_run(_req(_KRON3, ["coxeter_spectral"]),
                 tmp_path)["results"]["coxeter_spectral"]
    assert b["cyclotomic"] is False and b["coxeter_order"] is None
    assert b["outside_unit_circle_count"] == 1
    assert b["spectral_radius"]["minpoly"] == [1, -7, 1] and b["spectral_radius"]["degree"] == 2
    assert b["spectral_radius"]["interval"] == ["6", "7"]
    assert b["lehmer_class_note"] and b["scope"]
    # references resolve to citation pairs (never a silent empty)
    assert len(b["citations"]) == len(set(tuple(c) for c in b["citations"]))
