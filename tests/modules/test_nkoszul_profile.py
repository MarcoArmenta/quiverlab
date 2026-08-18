# SPDX-License-Identifier: MIT
"""The assembled generalized-Koszulity profile (Plan 77 / R36): Plan 27's quadratic
verdict verbatim plus the four Plan-77 recognizers and the generation-degree table."""
import json

import pytest

from quiverlab import (PreprojectiveAlgebra, TruncatedPathAlgebra,
                       linear_path_algebra, truncated_polynomial)
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import koszul_profile, koszul_profile_block

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
def test_profile_kx3():
    b = koszul_profile_block(truncated_polynomial(3, field=QQ), top=8)
    assert b["kind"] == "koszul"
    assert b["quadratic_koszul"] is False              # Plan 27's verdict (not length 2)
    assert b["n_homogeneous"] == 3
    assert b["k2"]["generator_degrees"] == [1, 2]
    assert b["generation_degrees"]["1"][:4] == [0, 1, 3, 4]
    assert "berger_nonquadratic" in b["references"]
    assert b["almost_koszul"]["verdict"] is None       # q = 1: N-Koszul, not almost


@lit
def test_profile_preprojective_A3_almost():
    b = koszul_profile_block(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert b["quadratic_koszul"] is False
    assert b["almost_koszul"]["verdict"] is True
    assert b["almost_koszul"]["p"] == 2 and b["almost_koszul"]["q"] == 2
    assert b["k2"]["verdict"] is False                 # a degree-3 Yoneda generator
    assert b["n_homogeneous"] == 2                     # quadratic: defers to Plan 27


@xeng
def test_quadratic_fields_are_plan27_verbatim():
    for A in (linear_path_algebra(2, field=QQ), PreprojectiveAlgebra("A3", field=QQ),
              truncated_polynomial(4, field=QQ)):
        Y = A.ext_algebra(6)
        r = koszul_profile(A, top=6)
        assert r["quadratic_koszul"] is Y.koszul
        assert r["quadratic_reason"] == Y._koszul_reason
        expected = list(Y.koszul_obstruction) if Y.koszul_obstruction else None
        assert r["quadratic_obstruction"] == expected
        assert r["certified_through_degree"] == Y.certified_through_degree


@selfcert
def test_block_is_json_round_trippable():
    # The block crosses the wire (webapp) and the Pyodide bridge, so every key must be
    # JSON-safe: vertex-keyed maps stringified, no sets, no tuples.
    for A in (truncated_polynomial(3, field=QQ),
              TruncatedPathAlgebra("A4", 3, field=QQ),
              PreprojectiveAlgebra("A3", field=QQ)):
        b = koszul_profile_block(A, top=5)
        assert json.loads(json.dumps(b, sort_keys=True)) == b


@selfcert
def test_latex_line_never_claims_more_than_the_verdicts():
    # An unproved (window-bounded) claim must SAY "through degree W" in the summary.
    local = koszul_profile(truncated_polynomial(3, field=QQ), top=8)
    assert "through degree 8" in local["latex"]
    finite = koszul_profile(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert "3-Koszul" in finite["latex"] and "through degree" not in finite["latex"]
    almost = koszul_profile(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert "(2,2)" in almost["latex"].replace(" ", "")


@selfcert
def test_algebra_delegate_matches_the_module_function():
    A = truncated_polynomial(3, field=QQ)
    assert A.koszul_profile(6) == koszul_profile(A, 6)


@selfcert
def test_presentation_less_algebra_refuses_loudly():
    # A structure-constant algebra carries no path basis, so no internal degree and no
    # Yoneda presentation exist. The refusal is loud and typed (the runners turn it
    # into a clean 4xx, never a 500).
    from quiverlab import Algebra
    from quiverlab.errors import QuiverlabError
    one = QQ.one()
    sc = Algebra(QQ, [[[one]]], [one])                 # presentation-less: quiver None
    assert sc.quiver is None
    with pytest.raises(QuiverlabError):
        koszul_profile_block(sc, top=4)


@lit
def test_profile_kA4_J3_is_finite_gldim_3_koszul():
    b = koszul_profile_block(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert b["n_homogeneous"] == 3
    assert b["n_koszul"]["verdict"] is True and b["complete"] is True
    assert b["k2"]["verdict"] is True
    assert b["n_koszul"]["internal_degrees"]["1"] == [0, 1, 3]


@selfcert
@pytest.mark.parametrize("typ,expected", [
    ("A3", [0, 1, 2, 4, 5, 6, 8]),
    ("A4", [0, 1, 2, 5, 6, 7, 10]),
])
def test_a_QUADRATIC_algebra_still_gets_its_generation_degree_table(typ, expected):
    # REGRESSION. The generation-degree table is the HEADLINE primitive, so it must be
    # present whenever it is DEFINED (A length-graded), independently of which branch the
    # N-Koszul recognizer takes. It used to be read back off n_koszul["internal_degrees"],
    # which the N = 2 branch never fills -- it returns early to defer to Plan 27 -- so
    # EVERY quadratic algebra, including every Dynkin preprojective (exactly the showcase
    # almost-Koszul examples), rendered an EMPTY table in the GUI, the report and the
    # TikZ staircase while the data existed. Caught by the report-renderer test, whose
    # staircase drew nothing.
    b = koszul_profile_block(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert b["n_homogeneous"] == 2                     # the defer-to-Plan-27 branch
    assert b["generation_degrees"]                     # ... and the table is STILL there
    for v in b["generation_degrees"]:
        assert b["generation_degrees"][v] == expected
    # and it is the SAME data the almost-Koszul classifier read
    assert b["almost_koszul"]["break_hom_degree"] == 3


@selfcert
def test_every_simple_appears_in_the_table_not_just_the_first():
    b = koszul_profile_block(PreprojectiveAlgebra("A4", field=QQ), top=6)
    assert sorted(b["generation_degrees"]) == ["1", "2", "3", "4"]
