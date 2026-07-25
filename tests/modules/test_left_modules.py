"""Left modules alongside right, right as default (Plan 24, Tier 1b).

A left A-module IS a right A^op-module; the surface accepts ``side="right"|"left"``
with right the (byte-unchanged) default. These oracles pin the honest asymmetry
between the sides and the duality D that exchanges them.

Literature: Assem-Simson-Skowronski, *Elements of the Representation Theory of
Associative Algebras*, Vol. 1 (Cambridge, 2006), bib key ``ASS2006``. Left vs
right indecomposable projectives ``e_v A`` vs ``A e_v`` (Ch. II/III); the duality
``D = Hom_k(-,k): mod A <-> A-mod`` exchanging the sides (Ch. III); the
Auslander-Reiten translates and their duality (Ch. IV, Ch. VIII). Cited at the
chapter granularity we can verify -- no invented theorem numbers.
"""
import pytest

from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError


def _square(field=CC):
    # commutative square 1->2->4, 1->3->4 with a*b = c*d.
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _cn(n, field=CC):
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = [f"a{i}*a{i % n + 1}" for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _dv(M):
    return M.dimension_vector()


# --- side plumbing: default, attribute, base algebra, repr ------------------
def test_right_is_the_default_and_side_attribute():
    A = linear_path_algebra(2)
    P = A.projective(1)
    assert P.side == "right"
    assert P.base_algebra is A
    L = A.projective(1, side="left")
    assert L.side == "left"
    assert L.base_algebra is A                     # left A-module lives over A
    assert L.algebra is A.opposite()               # ... represented over A^op


def test_repr_mentions_side_only_when_left():
    A = linear_path_algebra(2)
    assert repr(A.projective(1)).startswith("P_1: right ")
    assert repr(A.projective(1, side="left")).startswith("P_1: left ")


def test_left_modules_are_genuine():
    A = _square(field=GF(5))
    for v in A.quiver.vertices:
        for X in (A.simple(v, side="left"), A.projective(v, side="left"),
                  A.injective(v, side="left")):
            ok, why = X.check_module()
            assert ok, why
            assert X.side == "left"


# --- the honest smallest asymmetry: kA_2 P(1) differs by side --------------
def test_ka2_projective_asymmetry():
    # ASS Ch. II: right P_v = e_v A (paths FROM v), left P_v = A e_v (paths TO v).
    A = linear_path_algebra(2)                     # 1 --a1--> 2
    assert _dv(A.projective(1)) == {1: 1, 2: 1}    # right e_1 A = {e_1, a1}
    assert _dv(A.projective(1, side="left")) == {1: 1, 2: 0}   # left A e_1 = {e_1}
    assert _dv(A.projective(2)) == {1: 0, 2: 1}    # right e_2 A = {e_2}
    assert _dv(A.projective(2, side="left")) == {1: 1, 2: 1}   # left A e_2 = {e_2, a1}


def test_ka2_full_spi_dimension_vectors_both_sides():
    A = linear_path_algebra(2)
    right = {"S": {v: _dv(A.simple(v)) for v in (1, 2)},
             "P": {v: _dv(A.projective(v)) for v in (1, 2)},
             "I": {v: _dv(A.injective(v)) for v in (1, 2)}}
    left = {"S": {v: _dv(A.simple(v, side="left")) for v in (1, 2)},
            "P": {v: _dv(A.projective(v, side="left")) for v in (1, 2)},
            "I": {v: _dv(A.injective(v, side="left")) for v in (1, 2)}}
    # simples agree (1-dim at v either side); P and I are Cartan-transposed mirrors.
    assert right["S"] == {1: {1: 1, 2: 0}, 2: {1: 0, 2: 1}}
    assert left["S"] == right["S"]
    assert right["P"] == {1: {1: 1, 2: 1}, 2: {1: 0, 2: 1}}
    assert left["P"] == {1: {1: 1, 2: 0}, 2: {1: 1, 2: 1}}
    assert right["I"] == {1: {1: 1, 2: 0}, 2: {1: 1, 2: 1}}
    assert left["I"] == {1: {1: 1, 2: 1}, 2: {1: 0, 2: 1}}


def test_square_spi_dimension_vectors_both_sides():
    # Commutative square: dim e_w P_v(right) = #paths v->w = C[v][w];
    # dim e_w P_v(left) = #paths w->v = C[w][v]. So left = Cartan-transposed.
    A = _square()
    C = A.cartan_matrix()
    verts = list(A.quiver.vertices)
    idx = {v: i for i, v in enumerate(verts)}
    for v in verts:
        rp = _dv(A.projective(v))
        lp = _dv(A.projective(v, side="left"))
        assert rp == {w: C[idx[v]][idx[w]] for w in verts}
        assert lp == {w: C[idx[w]][idx[v]] for w in verts}
    # I_v is D of the opposite projective; its dim vector is the transpose again.
    for v in verts:
        ri = _dv(A.injective(v))
        li = _dv(A.injective(v, side="left"))
        assert ri == {w: C[idx[w]][idx[v]] for w in verts}
        assert li == {w: C[idx[v]][idx[w]] for w in verts}


# --- D exchanges the two sides over the SAME algebra ------------------------
def test_D_exchanges_sides_over_same_algebra():
    # ASS Ch. III: D = Hom_k(-,k): mod A <-> A-mod is contravariant, D.D = id.
    A = linear_path_algebra(3)
    M = A.projective(1)                            # right A-module
    DM = M.dualize()
    assert DM.side == "left" and DM.base_algebra is A
    assert DM.dimension_vector() == M.dimension_vector()   # D preserves dim vectors
    DDM = DM.dualize()
    assert DDM.side == "right" and DDM.base_algebra is A
    assert M.is_isomorphic(DDM)                    # D.D = id, back to a right module


def test_injective_is_D_of_left_projective():
    # ASS Ch. III: I_v = D(A e_v). A e_v is the LEFT projective; D flips it to the
    # right injective. (This is the honest form of the Plan-23 identity.)
    for field in (CC, GF(5)):
        A = _square(field=field)
        for v in A.quiver.vertices:
            explicit = A.projective(v, side="left").dualize()   # right A-module
            implicit = A.injective(v)
            assert explicit.side == "right"
            assert explicit.dimension_vector() == implicit.dimension_vector()
            assert explicit.is_isomorphic(implicit), f"D disagreement at v={v}"


# --- tau on the left = tau over A^op, translated -----------------------------
def _action_eq(M, N):
    dom = M.domain
    if set(M.action) != set(N.action):
        return False
    for lab, a in M.action.items():
        b = N.action[lab]
        if len(a) != len(b) or any(len(ra) != len(rb) for ra, rb in zip(a, b)):
            return False
        if any(not dom.is_zero(dom.sub(x, y))
               for ra, rb in zip(a, b) for x, y in zip(ra, rb)):
            return False
    return True


@pytest.mark.parametrize("field", [CC, GF(3)])
def test_left_tau_is_right_tau_over_opposite(field):
    # A left A-module's tau is computed by the right-module tau over A^op; the
    # side translation is the identity on the representation, flipping only the tag.
    A = linear_path_algebra(3, field=field)
    Aop = A.opposite()
    for v in A.quiver.vertices:
        left_tau = A.simple(v, side="left").tau()          # left A-module
        right_tau = Aop.simple(v).tau()                    # right A^op-module
        assert left_tau.algebra is right_tau.algebra is Aop
        assert left_tau.side == "left" and right_tau.side == "right"
        assert left_tau.dimension_vector() == right_tau.dimension_vector()
        assert _action_eq(left_tau, right_tau)             # identical representation


def test_tau_preserves_side():
    A = linear_path_algebra(3)
    assert A.simple(1).tau().side == "right"
    assert A.simple(1, side="left").tau().side == "left"
    assert A.simple(1, side="left").tau_minus().side == "left"


# --- AR duality of translates: D exchanges tau and tau^- --------------------
@pytest.mark.parametrize("A", [linear_path_algebra(2), linear_path_algebra(3), _cn(3)])
def test_D_exchanges_tau_and_tauminus(A):
    # ASS Ch. IV/VIII (AR duality): tau(DM) ~= D(tau^- M) and D(tau M) ~= tau^-(DM).
    for v in A.quiver.vertices:
        M = A.simple(v)
        lhs, rhs = M.dualize().tau(), M.tau_minus().dualize()   # both left A-modules
        assert lhs.side == rhs.side == "left"
        assert lhs.dimension_vector() == rhs.dimension_vector()
        assert lhs.is_isomorphic(rhs)
        lhs2, rhs2 = M.tau().dualize(), M.dualize().tau_minus()
        assert lhs2.dimension_vector() == rhs2.dimension_vector()
        assert lhs2.is_isomorphic(rhs2)


# --- injective dimension of a left module = pd of its D-dual right module ----
@pytest.mark.parametrize("field", [CC, GF(3)])
def test_left_injdim_equals_pd_of_dual(field):
    # inj.dim of a left A-module M = pd of DM, where DM = D(M) is a RIGHT A-module.
    A = _square(field=field)
    for v in A.quiver.vertices:
        LM = A.simple(v, side="left")
        DM = LM.dualize()                          # right A-module over A
        assert DM.side == "right" and DM.base_algebra is A
        assert LM.injective_dimension(bound=12) == DM.projective_resolution(12).pd()


def test_left_injective_resolution_terms_are_left_injectives():
    A = linear_path_algebra(3)
    res = A.simple(3, side="left").injective_resolution(4)
    for n, term in enumerate(res.terms):
        if term is None:
            continue
        assert term.side == "left"                 # E^n of a left module is left-injective
        ok, why = term.check_module()
        assert ok, why


# --- local (single-vertex) algebras: sides agree ----------------------------
@pytest.mark.parametrize("n", [2, 3, 4])
@pytest.mark.parametrize("field", [CC, GF(2)])
def test_local_algebra_sides_agree(n, field):
    # k[x]/(x^n) is self-opposite; left and right invariants coincide.
    A = truncated_polynomial(n, field=field)
    for build in ("simple", "projective", "injective"):
        r = getattr(A, build)(1)
        el = getattr(A, build)(1, side="left")
        assert r.dimension_vector() == el.dimension_vector()
    # projective / injective dimension of the simple agree across sides
    assert (A.simple(1).projective_resolution(6).pd()
            == A.simple(1, side="left").projective_resolution(6).pd())
    assert (A.simple(1).injective_dimension(bound=6)
            == A.simple(1, side="left").injective_dimension(bound=6))


# --- comparing across sides is a category error, refused loudly -------------
def test_iso_refuses_across_sides():
    A = linear_path_algebra(2)
    right_s1 = A.simple(1)                          # dimvec {1:1, 2:0}
    left_s1 = A.simple(1, side="left")             # SAME dimvec {1:1, 2:0}
    assert right_s1.dimension_vector() == left_s1.dimension_vector()
    with pytest.raises(QuiverlabError):            # not False -- a category error
        right_s1.is_isomorphic(left_s1)
    with pytest.raises(QuiverlabError):
        A.hom(right_s1, left_s1)
    with pytest.raises(QuiverlabError):
        A.ext(right_s1, left_s1, 1)


def test_iso_refuses_across_different_base_algebras():
    A = linear_path_algebra(2)
    Aop = A.opposite()
    # both right modules, but over different (opposite) algebras: still a category error.
    with pytest.raises(QuiverlabError):
        A.simple(1).is_isomorphic(Aop.simple(1))


# --- A.module(...) constructor honours the side flag ------------------------
def test_algebra_module_constructor_both_sides():
    A = linear_path_algebra(2)
    # right module built explicitly = S_1 (+) S_2 impostor (arrow acts 0)
    R = A.module({1: 1, 2: 1}, {"a1": [[0, 0], [0, 0]]}, name="R")
    assert R.side == "right" and R.base_algebra is A
    ok, why = R.check_module()
    assert ok, why
    # left module: data is the A^op representation (arrow a1 reversed, acts 0)
    L = A.module({1: 1, 2: 1}, {"a1": [[0, 0], [0, 0]]}, side="left", name="L")
    assert L.side == "left" and L.base_algebra is A
    ok, why = L.check_module()
    assert ok, why
