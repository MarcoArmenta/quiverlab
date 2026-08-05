"""The ext_algebra + recognizers kinds through the spec runner (Plan 38).

``quiverlab.hpc.spec`` exposes the dict->request validator as ``parse_request``
and the runner as ``run``; this drives the real ``_dispatch`` branches for
ext_algebra / recognizers. The webapp runner (``webapp.server.runner.run_spec``)
delegates to this same dispatch -- the byte-stable goldens
``ext_algebra_exterior_gf7`` / ``recognizers_gentle_a3`` in
``tests/webapp/_runner_goldens.json`` pin that path, and the Pyodide twin
(``docs/gui/runner.py``) is checked byte-identical in
``tests/webapp/test_koszul_exposure_p38.py``.
"""
from quiverlab.hpc.spec import parse_request
from quiverlab.hpc.spec import run as spec_run

# the Koszul poster child: E(k<x,y>/(x^2, y^2, x*y+y*x)) is the symmetric algebra.
_EXTERIOR = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1], "y": [1, 1]},
             "relations": ["x*x", "y*y", "x*y+y*x"],
             "field": {"kind": "GF", "p": 7, "n": 1}}
_GENTLE = {"kind": "quiver", "vertices": [1, 2, 3],
           "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": ["a*b"],
           "field": {"kind": "GF", "p": 5, "n": 1}}
# k[x]/(x^2): infinite gl.dim, non-unimodular Cartan (det 2) -> no form type.
_LOOP2 = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
          "relations": ["x*x"], "field": {"kind": "GF", "p": 5, "n": 1}}


def _req(algebra, compute):
    return parse_request({"schema": 1, "algebra": algebra, "compute": compute,
                          "artifacts": {"pdf": False, "tikz": False}})


def test_ext_algebra_koszul(tmp_path):
    b = spec_run(_req(_EXTERIOR, ["ext_algebra:0..4"]),
                 tmp_path)["results"]["ext_algebra"]
    assert b["koszul"] is True
    assert b["graded_dims"][:3] == [1, 2, 3]            # E(A) = symmetric algebra
    assert b["obstruction"] is None
    assert b["certified_through_degree"] == 4
    assert b["generators_by_degree"] == {"1": 2}
    assert b["relations_by_degree"] == {"2": 1}
    assert b["references"] == ["priddy", "froberg_koszul", "polishchuk_positselski"]
    assert b["citations"]


def test_ext_algebra_default_top(tmp_path):
    # no range on the scalar kind => default top 6
    b = spec_run(_req(_EXTERIOR, ["ext_algebra"]),
                 tmp_path)["results"]["ext_algebra"]
    assert b["top"] == 6 and b["certified_through_degree"] == 6
    assert b["graded_dims"][:3] == [1, 2, 3]


def test_recognizers_gentle(tmp_path):
    b = spec_run(_req(_GENTLE, ["recognizers"]),
                 tmp_path)["results"]["recognizers"]
    f = b["flags"]
    assert f["is_gentle"] is True and f["is_string"] is True
    assert f["is_special_biserial"] is True and f["is_hereditary"] is False
    assert f["is_semisimple"] is False and f["is_nakayama"] is True
    assert b["dynkin_type"] == "A3" and b["form_type"] == "finite"
    assert b["references"] == ["assem_book"] and b["citations"]


def test_recognizers_form_type_null_when_cartan_not_unimodular(tmp_path):
    b = spec_run(_req(_LOOP2, ["recognizers"]),
                 tmp_path)["results"]["recognizers"]
    assert b["form_type"] is None                       # k[x]/(x^2): det C = 2
    assert b["dynkin_type"] is None                     # a loop is not a diagram
    # every flag is a plain bool here (no refusal on a presented algebra)
    assert all(isinstance(v, bool) for v in b["flags"].values())
    assert b["flags"]["is_selfinjective"] is True        # k[x]/(x^2) is self-injective
