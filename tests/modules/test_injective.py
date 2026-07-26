"""Injective resolutions + injective dimension (Plan 23, Tier 1b item 2).

Theory oracles (no [qpa] extra):
  * inj.dim_A(M) = pd_{A^op}(DM); self-injective => 0 (projective) or inf (else);
  * hereditary kA_n: inj.dim of every module <= 1, max over simples = gl.dim = 1;
  * commutative square: max_v inj.dim S_v = gl.dim = 2;
  * injective envelope E^0 of a simple S_v is I_v.
"""
import pytest

from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial

pytestmark = [pytest.mark.oracle_literature]


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _cn(n, field=CC):
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = [f"a{i}*a{i % n + 1}" for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


@pytest.mark.parametrize("field", [CC, GF(2)])
def test_selfinjective_injdim_zero_iff_projective(field):
    A = truncated_polynomial(2, field=field)         # k[x]/(x^2), self-injective
    P = A.projective(1)                               # projective = injective here
    S = A.simple(1)                                   # non-projective
    assert P.injective_dimension() == 0
    assert S.injective_dimension(bound=8) is None     # infinite injective dimension
    # dual statement: projective dimension of the simple is likewise infinite
    assert S.projective_resolution(8).pd() is None


@pytest.mark.parametrize("field", [CC, GF(3)])
def test_hereditary_injdim_bounded_by_gldim(field):
    A = linear_path_algebra(3, field=field)           # gl.dim 1
    injdims = [A.simple(v).injective_dimension() for v in A.quiver.vertices]
    assert all(d in (0, 1) for d in injdims)
    assert max(injdims) == int(A.global_dimension())  # = 1
    assert A.simple(1).injective_dimension() == 0     # S_1 is injective in kA_3


def test_commutative_square_injdim_equals_gldim():
    A = _square()
    injdims = [A.simple(v).injective_dimension() for v in A.quiver.vertices]
    assert max(injdims) == 2                           # gl.dim(commutative square) = 2
    assert int(A.global_dimension()) == 2


@pytest.mark.parametrize("field", [CC, GF(5)])
def test_injective_envelope_of_simple_is_Iv(field):
    A = linear_path_algebra(3, field=field)
    for v in A.quiver.vertices:
        res = A.simple(v).injective_resolution(4)
        E0 = res.terms[0]
        assert E0.dimension_vector() == A.injective(v).dimension_vector()
        assert E0.is_isomorphic(A.injective(v))
        ok, why = E0.check_module()
        assert ok, why


def test_injective_resolution_matches_injdim():
    A = _square()
    for v in A.quiver.vertices:
        res = A.simple(v).injective_resolution(6)
        assert res.injective_dimension() == A.simple(v).injective_dimension()


def test_injective_resolution_terms_are_injective():
    # each E^n is a sum of I_u; its dimension vector is the sum of injective dimvecs.
    A = linear_path_algebra(3)
    res = A.simple(3).injective_resolution(4)
    for n, verts in enumerate(res.vertices):
        term = res.terms[n]
        if term is None:
            continue
        want = {w: 0 for w in A.quiver.vertices}
        for u in verts:
            for w, d in A.injective(u).dimension_vector().items():
                want[w] += d
        assert term.dimension_vector() == want
