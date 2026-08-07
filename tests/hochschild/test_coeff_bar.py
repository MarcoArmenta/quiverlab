"""HH(A,M) via the bar reference. None short-circuit is byte-identical; regular(A)
through the general path == None (M=A oracle); exact path == ported GF(p) reference."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.errors import QuiverlabError


@pytest.mark.oracle_crossengine
def test_regular_coefficient_equals_no_coefficient():
    A = ql.truncated_polynomial(3, field=ql.CC)
    base_h = A.hochschild_homology(5, engine="bar").dims
    base_c = A.hochschild_cohomology(5, engine="bar").dims
    Mreg = Bimodule.regular(A)
    assert A.hochschild_homology(5, engine="bar", coefficients=Mreg).dims == base_h
    assert A.hochschild_cohomology(5, engine="bar", coefficients=Mreg).dims == base_c


@pytest.mark.oracle_crossengine
def test_exact_bar_matches_ported_gfp_reference():
    from quiverlab.engine import bimodule as eb
    from quiverlab.engine.adapter import to_engine
    A = ql.truncated_polynomial(3, field=ql.GF(32003))
    eng = to_engine(A.unit_adapted())
    Md = Bimodule.dual(A)
    ours = A.hochschild_cohomology(5, engine="bar", coefficients=Md).dims
    ref = eb.hochschild_cohomology_with_coefficients(eng, eb.dual_bimodule(eng), 5)[32003]
    assert ours == ref


@pytest.mark.oracle_crossengine
def test_exact_bar_homology_matches_ported_gfp_reference():
    from quiverlab.engine import bimodule as eb
    from quiverlab.engine.adapter import to_engine
    A = ql.truncated_polynomial(3, field=ql.GF(32003))
    eng = to_engine(A.unit_adapted())
    Md = Bimodule.dual(A)
    ours = A.hochschild_homology(5, engine="bar", coefficients=Md).dims
    ref = eb.hochschild_homology_with_coefficients(eng, eb.dual_bimodule(eng), 5)[32003]
    assert ours == ref


@pytest.mark.oracle_selfcert
def test_none_is_byte_identical_short_circuit():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    t0 = A.hochschild_cohomology(3)
    assert set(vars(t0)) == {"dims", "kind", "top", "algebra_repr", "engine", "references"}
    # coefficients provenance is only surfaced for a non-None coefficient
    t1 = A.hochschild_cohomology(3, coefficients=Bimodule.dual(A))
    assert getattr(t1, "coefficients", None) == "D(A)"


@pytest.mark.oracle_selfcert
def test_chain_accessor_dims_match_public_bar():
    # cross-check the P54 chain-level accessor against the public HH-with-coeff dims.
    from quiverlab.hochschild.coefficients import twisted_homology_classes
    A = ql.truncated_polynomial(3, field=ql.CC)
    M = Bimodule.twisted_by_nakayama(A)
    top = 4
    data = twisted_homology_classes(A, M, top)
    dims = A.hochschild_homology(top, engine="bar", coefficients=M).dims
    assert [len(d["classes"]) for d in data] == dims


def test_fast_refuses_coefficients():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    with pytest.raises(QuiverlabError):
        A.hochschild_homology(2, engine="fast", coefficients=Bimodule.dual(A))


def test_coefficient_algebra_mismatch_refused():
    A = ql.truncated_polynomial(3, field=ql.CC)
    B = ql.truncated_polynomial(3, field=ql.CC)   # a DIFFERENT object
    MB = Bimodule.dual(B)
    with pytest.raises(QuiverlabError):
        A.hochschild_cohomology(3, engine="bar", coefficients=MB)
