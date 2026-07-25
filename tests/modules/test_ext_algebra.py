"""The Yoneda / Ext-algebra E(A) = Ext^*_A(A/J, A/J) as a graded quiver-with-relations
presentation (Plan 27). Worked anchors + the 7-oracle theory battery, the locked
left-to-right corner convention, byte-reproducible products, `as_algebra()` round-trips
and honest refusals, the Domain spread, and the monomial Anick gate.

Theory oracles (all run without [qpa]):
  1. k[x]/x^2            E = k[y],           Koszul
  2. k[x]/x^n (n>=3)     E = k[y,z]/(y^2),   not Koszul (non-quadratic), char-independent
  3. hereditary kQ       E = kQ/J^2,         E^{>=2}=0, finite, Koszul, self-hosting
  4. rad^2=0 kQ/J^2      E = kQ (free),      dim E^n = #paths_n
  5. quantum CI          E = quantum plane,  dim E^n = n+1, Koszul (G-quadratic)
  6. commutative square  gl.dim 2, dims (4,4,1), E ~= A, self-hosting
  7. semisimple          E = A in degree 0

Sources: Priddy, Trans. AMS 152 (1970) (G-quadratic => Koszul, key ``priddy``);
Froberg, LNPAM 205 (1999) (Hilbert-series criterion, key ``froberg_koszul``);
Green-Solberg-Zacharia, Trans. AMS 353 (2001) (minimal resolutions, key
``minimal_resolution``); Assem-Simson-Skowronski (2006), Ch. III (key ``assem_book``).
The Anick chain oracle reuses the Bardzell/CS S-sequence combinatorics
(arXiv:1406.2300; key ``chouhy_solotar``)."""
from collections import Counter

import pytest

from quiverlab import (CC, GF, Quiver, QuantumCI, RadicalSquareZero,
                       linear_path_algebra, truncated_polynomial)
from quiverlab.errors import QuiverlabError

# the koszul sibling (quadraticity / G-quadratic / Froberg) may or may not have landed;
# gate every verdict-dependent assertion on its availability.
try:
    import quiverlab.modules.koszul  # noqa: F401
    _HAS_KOSZUL = True
except ImportError:                                        # pragma: no cover
    _HAS_KOSZUL = False

_needs_koszul = pytest.mark.skipif(
    not _HAS_KOSZUL, reason="modules/koszul.py (Koszul verdict) not landed yet")


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def _square(field=CC):
    # commutative square with the single relation a*b - c*d (cf. tests/modules/test_ext.py)
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _gens(P):
    return {d: len(v) for d, v in P.generators_by_degree.items()}


def _rels(P):
    return {d: len(v) for d, v in P.relations_by_degree.items()}


def _anick_betti(A, n, max_degree):
    """Independent monomial oracle: for a monomial kQ/I, dim Ext^n(S_i,S_j) = #{Anick
    n-chains i->j}. Reuse the Bardzell/CS S-sequence the way the CS tests do (SSequence
    over reduction_system_of); count chains by (source, target) corner. Lives test-side:
    there is no public engine surface for it, and an independent oracle should not share
    the engine under test."""
    from quiverlab.resolutions_cs.ambiguities import SSequence
    from quiverlab.resolutions_cs.build import reduction_system_of
    ss = SSequence(reduction_system_of(A), max_degree=max_degree)
    return Counter((c.o, c.t) for c in ss.S(n))


def _engine_betti(P, n):
    eng = P._eng
    out = Counter()
    for i in eng.verts:
        for j in eng.verts:
            d = eng.ext_dim(i, n, j)
            if d:
                out[(i, j)] = d
    return out


# --------------------------------------------------------------------------------------
# oracle 1 -- k[x]/x^2  =>  E = k[y]
# --------------------------------------------------------------------------------------
def test_kx2_is_polynomial_ring():
    P = truncated_polynomial(2).ext_algebra(top=6)
    assert P.graded_dims_through(6) == [1] * 7          # dim E^n = 1 for all n
    assert _gens(P) == {1: 1}                            # a single degree-1 generator y
    assert P.relations_by_degree == {}                  # E = k[y] is free: no relations
    # infinite gl.dim (self-injective): not certified finite-dimensional, but Koszul =>
    # certified finitely generated (in degree 1).
    assert P.is_finite_dimensional is None
    assert P.is_finitely_generated_certified is True


@_needs_koszul
def test_kx2_is_koszul():
    P = truncated_polynomial(2).ext_algebra(top=4)
    assert P.koszul is True and P.koszul_obstruction is None


# --------------------------------------------------------------------------------------
# oracle 2 -- k[x]/x^n (n>=3)  =>  E = k[y,z]/(y^2), char-independent, NOT Koszul
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("n", [3, 4])
@pytest.mark.parametrize("field", [None, GF(2), GF(3), GF(32003)])
def test_kxn_is_kyz_mod_ysquared(n, field):
    A = truncated_polynomial(n, field=field)
    P = A.ext_algebra(top=4)
    assert P.graded_dims_through(4) == [1] * 5           # dims all 1 (one vertex)
    assert _gens(P) == {1: 1, 2: 1}                      # y (deg 1), z (deg 2)
    # the y^2 = 0 relation appears in degree 2 (char-independent: y_1 = .x^{n-2})
    assert len(P.relations_by_degree.get(2, [])) == 1
    rel = P.relations_by_degree[2][0]
    assert rel.degree == 2 and rel.word_length() == [2]  # a single length-2 word (y . y)


@_needs_koszul
@pytest.mark.parametrize("n", [3, 4])
@pytest.mark.parametrize("field", [None, GF(2), GF(3), GF(32003)])
def test_kxn_not_koszul(n, field):
    P = truncated_polynomial(n, field=field).ext_algebra(top=4)
    assert P.koszul is False
    assert P.koszul_obstruction is not None and P.koszul_obstruction[0] == 2
    assert P.is_finite_dimensional is None
    assert P.is_finitely_generated_certified is None      # deg-2 gen, higher not excluded


def test_ysquared_zero_iff_n_ge_3():
    # the discriminating product: y*y = 0 for n>=3, y*y != 0 for n=2 (E = k[y]).
    for n, expect_zero in [(2, False), (3, True), (5, True)]:
        P = truncated_polynomial(n).ext_algebra(top=3)
        eng, dom = P._eng, P.algebra.domain
        s = eng.corner_basis(1, 1, 1)[0]
        coords, _ = eng.product((1, 1, s), (1, 1, s))
        assert all(dom.is_zero(c) for c in coords) is expect_zero


# --------------------------------------------------------------------------------------
# oracle 3 -- hereditary kA_n  =>  E = kQ/J^2, E^{>=2}=0, self-hosting, Koszul
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("n", [2, 3])
def test_hereditary_ext_algebra(n):
    A = linear_path_algebra(n)
    P = A.ext_algebra(top=4)
    # gl.dim 1 (hereditary): complete, finite-dimensional.
    assert P.certified_through_degree == 1
    assert P.is_finite_dimensional is True
    assert P.is_finitely_generated_certified is True
    # E^{>=2} = 0
    hil = P.hilbert_matrix_through(3)
    assert all(x == 0 for deg in (2, 3) for row in hil[deg] for x in row)
    # generators = the arrows (n-1 of them); the ext-quiver equals Q (corner convention)
    assert _gens(P) == {1: n - 1}
    assert Counter(P.ext_quiver.arrows.values()) == Counter(A.quiver.arrows.values())
    # relations = J^2 (vanishing length-2 paths): for A_n exactly n-2 of them
    assert _rels(P).get(2, 0) == max(0, n - 2)


def test_hereditary_as_algebra_roundtrips():
    A = linear_path_algebra(3)
    E = A.ext_algebra(top=3).as_algebra()
    # E = kQ/J^2 on the same quiver: dim = #vertices + #arrows = 3 + 2 = 5, a genuine
    # associative unital Algebra (validated on build). Verify the defining J^2 structure
    # directly: the two arrow classes are orthogonal-idempotent-graded and compose to 0.
    assert E.dim == 5
    dom, idx = E.domain, {l: i for i, l in enumerate(E.basis_labels)}

    def bvec(label):
        v = [dom.zero()] * E.dim
        v[idx[label]] = dom.one()
        return v
    a, b = bvec("E1_1_2_0"), bvec("E1_2_3_0")             # the two degree-1 generators
    assert all(dom.is_zero(x) for x in E.multiply(a, b))  # J^2 = 0: arrows compose to 0
    e1 = bvec("E0_1_1_0")
    assert E.multiply(e1, a) == a                         # e_1 acts as the source idempotent


@_needs_koszul
def test_hereditary_is_koszul():
    assert linear_path_algebra(3).ext_algebra(top=3).koszul is True


# --------------------------------------------------------------------------------------
# oracle 4 -- rad^2 = 0  =>  E = kQ (free), dim E^n = #paths_n
# --------------------------------------------------------------------------------------
def test_radsq_one_loop_is_polynomial():
    # one loop, rad^2 = 0 is exactly k[x]/x^2: E = k[y], dims all 1.
    A = RadicalSquareZero(Quiver([1], {"x": (1, 1)}))
    P = A.ext_algebra(top=5)
    assert P.graded_dims_through(5) == [1] * 6
    assert _gens(P) == {1: 1} and P.relations_by_degree == {}


def test_radsq_two_loops_dim_is_2n():
    # two loops, rad^2 = 0: E = kQ free on 2 generators, dim E^n = 2^n.
    A = RadicalSquareZero(Quiver([1], {"x": (1, 1), "y": (1, 1)}))
    # dims 2^n straight off the Betti numbers of the periodic resolution (cheap, deep).
    res = A.simple(1).projective_resolution(6)
    assert [res.betti(n) for n in range(7)] == [2 ** n for n in range(7)]
    # the presentation is the free path algebra: 2 generators, no relations.
    P = A.ext_algebra(top=3)
    assert P.graded_dims_through(3) == [1, 2, 4, 8]
    assert _gens(P) == {1: 2} and P.relations_by_degree == {}


def test_radsq_line_worked_anchor():
    # 1 -> 2 -> 3, rad^2 = 0: E = kQ, dim E^n = #paths_n = [3, 2, 1]; the length-2 path
    # survives as y_a . y_b != 0 in corner (1, 3) (the locked left-to-right convention).
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}))
    P = A.ext_algebra(top=3)
    assert P.graded_dims_through(P.certified_through_degree) == [3, 2, 1]
    assert P.is_finite_dimensional is True                 # gl.dim 2
    assert _rels(P) == {}                                  # free: no relations
    eng, dom = P._eng, A.domain
    # Ext^1(S_1, S_2) and Ext^1(S_2, S_3) are 1-dim (the two arrows); their product lands
    # in Ext^2(S_1, S_3) and is nonzero.
    sa = eng.corner_basis(1, 1, 2)[0]
    sb = eng.corner_basis(2, 1, 3)[0]
    coords, k = eng.product((1, 1, sa), (2, 1, sb))
    assert k == 3 and any(not dom.is_zero(c) for c in coords)
    # ... while dim Ext^2(S_1, S_3) = 1 and every other corner in degree 2 vanishes.
    assert eng.ext_dim(1, 2, 3) == 1
    assert _engine_betti(P, 2) == Counter({(1, 3): 1})


# --------------------------------------------------------------------------------------
# oracle 5 -- quantum complete intersection  =>  E = quantum plane, dim E^n = n+1
# --------------------------------------------------------------------------------------
def test_quantum_ci_ext_algebra():
    A = QuantumCI("2")                                     # k<x,y>/(x^2, y^2, x*y + 2 y*x)
    P = A.ext_algebra(top=5)
    assert P.graded_dims_through(5) == [1, 2, 3, 4, 5, 6]  # dim E^n = n + 1 = |S_n|
    assert _gens(P) == {1: 2}                              # two degree-1 generators
    assert len(P.relations_by_degree.get(2, [])) == 1     # a single degree-2 relation
    assert P.relations_by_degree[2][0].word_length() == [2]  # quadratic (length-2 words)


@_needs_koszul
def test_quantum_ci_is_koszul():
    P = QuantumCI("2").ext_algebra(top=4)
    assert P.koszul is True                                # G-quadratic
    assert P.is_finitely_generated_certified is True


# --------------------------------------------------------------------------------------
# oracle 6 -- commutative square  =>  gl.dim 2, dims (4,4,1), E ~= A, self-hosting
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("field", [CC, GF(7)])
def test_commutative_square_ext_algebra(field):
    A = _square(field)
    P = A.ext_algebra(top=2)
    assert P.certified_through_degree == 2                 # gl.dim 2
    assert P.is_finite_dimensional is True
    assert P.graded_dims_through(2) == [4, 4, 1]           # 4 vertices, 4 arrows, 1 relation
    assert _gens(P) == {1: 4}                              # generated in degree 1 (Koszul)
    # the single degree-2 relation reproduces a*b - c*d (2 length-2 words in corner (1,4))
    rel = P.relations_by_degree[2][0]
    assert rel.source == 1 and rel.target == 4
    assert len(rel.terms) == 2 and rel.word_length() == [2]
    # E ~= A (Koszul self-dual): as_algebra() self-hosts with matching total dimension.
    E = P.as_algebra()
    assert E.dim == A.dim == 9


@_needs_koszul
def test_commutative_square_is_koszul():
    assert _square().ext_algebra(top=2).koszul is True


# --------------------------------------------------------------------------------------
# oracle 7 -- semisimple  =>  E = A in degree 0
# --------------------------------------------------------------------------------------
def test_semisimple_ext_algebra():
    A = Quiver([1, 2], {}).algebra(relations=[])           # k x k, no arrows
    P = A.ext_algebra(top=3)
    assert P.certified_through_degree == 0                 # gl.dim 0
    assert P.graded_dims_through(0) == [2]                 # E = R = k^{Q_0} in degree 0
    assert P.generators_by_degree == {} and P.relations_by_degree == {}
    assert P.is_finite_dimensional is True
    assert P.as_algebra().dim == 2


# --------------------------------------------------------------------------------------
# byte-reproducibility of the (canonicalised) Yoneda product
# --------------------------------------------------------------------------------------
def test_product_is_byte_reproducible():
    from quiverlab.modules.ext_algebra import ext_algebra
    A = truncated_polynomial(3)
    P1, P2 = ext_algebra(A, 4), ext_algebra(A, 4)
    s = P1._eng.corner_basis(1, 1, 1)[0]
    # a product past the bar window (y * z in degree 3) computed by two independent lifts
    z = P1._eng.corner_basis(1, 2, 1)[0]
    a = P1._eng.product((1, 1, s), (1, 2, z))
    b = P2._eng.product((1, 1, s), (1, 2, z))
    assert a == b                                          # identical bytes, both lifts


# --------------------------------------------------------------------------------------
# corner bookkeeping: dim Ext^1(S_i,S_j) = #arrows i->j; Ext^2 in the relation corner
# --------------------------------------------------------------------------------------
def test_corner_bookkeeping_matches_quiver():
    A = _square()
    P = A.ext_algebra(top=2)
    eng = P._eng
    # Ext^1 corners are exactly the arrows (one per arrow, correct source/target)
    assert _engine_betti(P, 1) == Counter(A.quiver.arrows.values())
    # Ext^2 sits only in corner (1, 4) -- the single minimal relation a*b - c*d
    assert _engine_betti(P, 2) == Counter({(1, 4): 1})
    # every degree-0 corner is the diagonal (E^0 = R = k^{Q_0})
    assert _engine_betti(P, 0) == Counter({(v, v): 1 for v in (1, 2, 3, 4)})


# --------------------------------------------------------------------------------------
# monomial Anick gate: the module engine equals the chain counts degreewise
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("name,A", [
    ("k[x]/x^3", truncated_polynomial(3)),
    ("rad2_line", RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}))),
    ("straddle", Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y*x"])),
])
def test_monomial_anick_gate(name, A):
    pytest.importorskip("quiverlab.groebner")
    P = A.ext_algebra(top=6)
    for n in range(7):
        assert _engine_betti(P, n) == _anick_betti(A, n, 6), f"{name} degree {n}"


# --------------------------------------------------------------------------------------
# honest refusals + truncation language
# --------------------------------------------------------------------------------------
def test_as_algebra_refuses_when_not_finite():
    # k[x]/x^2 is self-injective: E = k[y] is infinite-dimensional. as_algebra must refuse
    # loudly with truncation language, never fabricate a finite algebra.
    P = truncated_polynomial(2).ext_algebra(top=4)
    with pytest.raises(QuiverlabError) as exc:
        P.as_algebra()
    msg = str(exc.value).lower()
    assert "truncation" in msg or "finite" in msg
    assert str(P.certified_through_degree) in str(exc.value)


def test_hilbert_beyond_certified_is_loud_for_infinite():
    # for an infinite-gl.dim algebra, asking for dims beyond the certified degree must
    # raise rather than silently pad zeros.
    P = truncated_polynomial(3).ext_algebra(top=3)
    P.hilbert_matrix_through(3)                             # in range: fine
    with pytest.raises(QuiverlabError):
        P.hilbert_matrix_through(4)                         # out of certified range: loud


def test_structure_constant_algebra_refused():
    # no quiver presentation -> no path basis -> loud, uniform QuiverlabError (as with the
    # other module-homology surfaces).
    from quiverlab import Algebra
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        A.ext_algebra(3)


# --------------------------------------------------------------------------------------
# the public delegator
# --------------------------------------------------------------------------------------
def test_delegator_smoke():
    from quiverlab.modules.ext_algebra import YonedaPresentation
    P = truncated_polynomial(2).ext_algebra(top=3)
    assert isinstance(P, YonedaPresentation)
    assert repr(P)                                         # renders without error
