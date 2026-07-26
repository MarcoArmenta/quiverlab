"""Ext and Tor where BOTH arguments are genuinely non-projective, non-injective,
AND non-simple modules (Plan 34 -- Marco's release feedback: "test Ext and Tor
with non-projective, non-injective, non-simple modules"; the standing batteries
lean on the S(v)/P(v)/I(v) builtins, which under-tests the machinery).

Every module in the menagerie is CERTIFIED non-trivial in-test (``_assert_interior``
/ ``_assert_interior_left``): non-simple, and -- being indecomposable -- projective
iff isomorphic to some ``P_v`` and injective iff isomorphic to some ``I_v``, both
refuted by the exact ``is_isomorphic`` certificate.  Every numeric pin is either a
literature/theory identity or a cross-engine/duality anchor -- never an uncited
number.

Menagerie:
  * kA_5 interior intervals [2,3], [3,4], [2,4] (built by EXPLICIT per-arrow
    matrices via ``A.module`` -- the Plan-26 no-code constructor path).  Over the
    equioriented A_n the indecomposables that are neither projective ([i,n]) nor
    injective ([1,j]) nor simple are exactly the intervals [i,j] with 1 < i < j < n.
    kA_5 is the SMALLEST equioriented A_n carrying >= 3 pairwise-distinct such
    interior indecomposables -- which the pairwise Euler/Ext battery below needs.
    (kA_3 has NONE; kA_4 has exactly ONE, namely [2,3] with dim vector (0,1,1,0) --
    itself interior, so "every non-simple indecomposable over kA_4 is projective or
    injective" is FALSE; kA_5 has three, [2,3],[3,4],[2,4].  The kA_4 counterexample
    is asserted in ``test_menagerie_modules_are_interior``.)
  * the commutative square kQ/(ab-cd): rad P_1 and P_1/soc P_1 (both interior).
  * a self-injective cyclic Nakayama kZ_3/J^3 (projective = injective, so a
    non-projective module is automatically non-injective): rad P_1, a length-2
    uniserial.

Pins:
  1. Hereditary Euler-form identity (ASS2006, Ch. III, ``assem_book``): over the
     hereditary kA_n, dim Hom(M,N) - dim Ext^1(M,N) = <dim M, dim N>, the Euler
     bilinear form <x,y> = x^T C^{-1} y from the repo Cartan matrix C
     (``A.cartan_matrix``).  The convention is pinned by computing the form TWO
     independent ways -- from C, and from the arrow combinatorial formula
     sum_v x_v y_v - sum_{a: i->j} x_i y_j -- and demanding they agree; and
     Ext^{>=2} = 0 (hereditary).
  2. Duality anchor (Plan 29, ``tensor_product`` / ``module_ext``): degreewise
     dim Tor_n^A(M,N) = dim Ext^n(M, DN) for a RIGHT M and a LEFT N (D side-aware:
     DN of a left module is a right module), cross-checked against the
     independently written Ext engine.
  3. Auslander-Reiten pin (ASS2006, Ch. IV/V, ``assem_book``): for a non-projective
     indecomposable M the AR sequence 0 -> tau M -> E -> M -> 0 exists, so
     Ext^1(M, tau M) != 0; indecomposability is certified via
     ``is_indecomposable`` (char 0 here, so the trace-form certificate is rigorous).
  4. Balance of Tor (Plan 29): resolving M vs resolving N over A^op agree.

Sources: Assem-Simson-Skowronski, "Elements of the Representation Theory of
Associative Algebras, Vol. 1", Cambridge Univ. Press (2006) (key ``assem_book``);
Cartan-Eilenberg, "Homological Algebra", Princeton Univ. Press (1956) (key
``tensor_product``); Green-Solberg-Zacharia, Trans. AMS 353 (2001) (key
``module_ext``).
"""
import sympy
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.tor import tor_dims


# --------------------------------------------------------------------------- #
# Menagerie constructors + non-triviality certificates
# --------------------------------------------------------------------------- #
def _square(field=QQ):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _cyclic_nakayama(n, loewy, field=QQ):
    """kZ_n / J^loewy: the n-cycle 1->2->...->n->1 with every length-``loewy`` path
    zero.  Self-injective (projective = injective) for every n, loewy (a cyclic
    Nakayama algebra); its non-simple length < loewy uniserials are non-projective."""
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = ["*".join(f"a{((i - 1 + t) % n) + 1}" for t in range(loewy))
            for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _p_mod_soc(P):
    """P / soc P as an explicit quotient module (soc = intersection of the arrow
    kernels), used to manufacture an interior module over a non-hereditary algebra."""
    from quiverlab.modules import linalg_mod as lm
    from quiverlab.modules.radtopsoc import _intersect, quotient
    dom = P.domain
    inter = None
    for a in P.algebra.quiver.arrows:
        ker = lm.kernel_columns(P.action[a], dom)
        inter = ker if inter is None else _intersect(inter, ker, dom)
    return quotient(P, inter or [], name=f"{P.name}/soc")


def _interval(A, i, j, n, name=None):
    """The uniserial interval module [i,j] over the equioriented kA_n, by EXPLICIT
    per-arrow matrices (Plan-26 no-code path): 1 at each vertex i..j, the connecting
    arrow a_k: k -> k+1 acting by the identity between consecutive slots."""
    dim = j - i + 1
    dimvec = {v: (1 if i <= v <= j else 0) for v in range(1, n + 1)}
    maps = {}
    for k in range(1, n):
        blk = [[0] * dim for _ in range(dim)]
        if i <= k <= j - 1:
            blk[(k + 1) - i][k - i] = 1        # a_k maps slot(k) -> slot(k+1)
        maps[f"a{k}"] = blk
    return A.module(dimvec, maps, name=name or f"[{i},{j}]")


def _assert_interior(A, M):
    """Certify M (an indecomposable RIGHT module) is non-simple AND non-projective AND
    non-injective: non-simple by dim; indecomposable so projective iff iso some P_v,
    injective iff iso some I_v -- both refuted by the exact is_isomorphic certificate."""
    assert M.dim > 1, f"{M.name} is simple (dim 1)"
    assert M.is_indecomposable(), f"{M.name} is decomposable"
    for v in A.quiver.vertices:
        P, I = A.projective(v), A.injective(v)
        assert not (M.dim == P.dim and M.is_isomorphic(P)), f"{M.name} == P_{v} (projective)"
        assert not (M.dim == I.dim and M.is_isomorphic(I)), f"{M.name} == I_{v} (injective)"


def _assert_interior_left(A, N):
    """The LEFT-module twin of ``_assert_interior``: N non-simple, and (indecomposable)
    not isomorphic to any left projective P_v nor any left injective I_v."""
    assert N.side == "left", f"{N.name} is not a left module"
    assert N.dim > 1, f"{N.name} is simple (dim 1)"
    assert N.is_indecomposable(), f"{N.name} is decomposable"
    for v in A.quiver.vertices:
        P, I = A.projective(v, side="left"), A.injective(v, side="left")
        assert not (N.dim == P.dim and N.is_isomorphic(P)), f"{N.name} == left P_{v}"
        assert not (N.dim == I.dim and N.is_isomorphic(I)), f"{N.name} == left I_{v}"


def _dimvec(A, M):
    dv = M.dimension_vector()
    return [dv[v] for v in A.quiver.vertices]


def _euler_from_cartan(A, M, N):
    """<dim M, dim N> = x^T C^{-1} y from the repo Cartan matrix C[i][j] = dim e_i A e_j.
    Derivation: for a finite-gl.dim algebra the class of M in the projective basis is
    g_M = C^{-T} dim M (since dim P_v = row v of C), so sum_i (-1)^i dim Ext^i(M,N) =
    g_M . dim N = (dim M)^T C^{-1} dim N."""
    C = sympy.Matrix(A.cartan_matrix())
    x = sympy.Matrix(_dimvec(A, M))
    y = sympy.Matrix(_dimvec(A, N))
    return (x.T * C.inv() * y)[0, 0]


def _euler_arrow(A, M, N):
    """The hereditary Euler form of a path algebra kQ (ASS2006, Ch. III):
    <x,y> = sum_v x_v y_v - sum_{a: i->j in Q_1} x_i y_j."""
    x, y = M.dimension_vector(), N.dimension_vector()
    s = sum(x[v] * y[v] for v in A.quiver.vertices)
    for _a, (src, tgt) in A.quiver.arrows.items():
        s -= x[src] * y[tgt]
    return s


def _tor_anchor(A, M, N, top):
    """Return tor_dims(A, M, N, top) after asserting the Plan-29 duality anchor
    dim Tor_n(M,N) = dim Ext^n(M, DN) degreewise against the Ext engine (D side-aware:
    N.dualize() of the LEFT module N is a RIGHT module)."""
    tors = tor_dims(A, M, N, top)
    exts = [A.ext(M, N.dualize(), n) for n in range(top + 1)]
    assert tors == exts, (M.name, N.name, tors, exts)
    return tors


# --------------------------------------------------------------------------- #
# Menagerie self-documentation: every module is non-P, non-I, non-simple
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_menagerie_modules_are_interior():
    A5 = linear_path_algebra(5, field=QQ)
    for ij in ((2, 3), (3, 4), (2, 4)):
        _assert_interior(A5, _interval(A5, ij[0], ij[1], 5))

    # Why kA_5 (not a smaller kA_n): it is the smallest equioriented A_n with >= 3
    # interior indecomposables.  The load-bearing counterexample to the tempting
    # false claim "over kA_4 every non-simple indecomposable is projective or
    # injective": kA_4's interval [2,3] (dim vector (0,1,1,0)) is itself INTERIOR.
    A4 = linear_path_algebra(4, field=QQ)
    _assert_interior(A4, _interval(A4, 2, 3, 4))

    S = _square()
    _assert_interior(S, S.projective(1).radical())
    _assert_interior(S, _p_mod_soc(S.projective(1)))

    B = _cyclic_nakayama(3, 3)
    assert B.is_selfinjective()                    # projective = injective here
    _assert_interior(B, B.projective(1).radical()) # non-projective => non-injective


# --------------------------------------------------------------------------- #
# Pin 1: hereditary Euler form  dim Hom - dim Ext^1 = <dim M, dim N>
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
@pytest.mark.oracle_crossengine
def test_hereditary_euler_form_kA5_interior_pairs():
    """kA_5 (hereditary): for every ordered pair of interior indecomposables
    [2,3],[3,4],[2,4], dim Hom(M,N) - dim Ext^1(M,N) equals the Euler form
    <dim M, dim N>, computed BOTH from the Cartan matrix and from the arrow formula
    (they must agree -- the convention pin), and Ext^{>=2} = 0 (ASS2006 III)."""
    A = linear_path_algebra(5, field=QQ)
    mods = [_interval(A, i, j, 5) for (i, j) in ((2, 3), (3, 4), (2, 4))]
    for M in mods:
        _assert_interior(A, M)
    for M in mods:
        for N in mods:
            euler_c = _euler_from_cartan(A, M, N)
            euler_a = _euler_arrow(A, M, N)
            assert euler_c == euler_a, (M.name, N.name, euler_c, euler_a)   # convention
            hom = A.hom(M, N)
            ext1 = A.ext(M, N, 1)
            assert hom - ext1 == euler_c, (M.name, N.name, hom, ext1, euler_c)
            assert A.ext(M, N, 2) == 0 and A.ext(M, N, 3) == 0             # hereditary


@pytest.mark.oracle_literature
def test_hereditary_euler_form_holds_on_a_decomposable_middle_term():
    """The Euler form is additive, so the identity survives on a DECOMPOSABLE module
    (here E = [2,3] (+) [3,4], a plain direct sum used purely to exercise Euler
    additivity -- NOT the middle term of an AR sequence).  E is
    non-projective and non-injective by support alone over kA_5: every projective
    summand P_i hits vertex 5 and every injective summand I_j hits vertex 1, but
    dim E = (0,1,2,1,0) vanishes at both ends."""
    A = linear_path_algebra(5, field=QQ)
    # E = [2,3] (+) [3,4]: dims v2=1, v3=2, v4=1; basis idx0=v2, idx1,idx2=v3, idx3=v4.
    dims = {2: 1, 3: 2, 4: 1}
    maps = {"a1": [[0] * 4 for _ in range(4)],
            "a2": [[0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],  # v2 -> v3 slot1
            "a3": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0]],  # v3 slot2 -> v4
            "a4": [[0] * 4 for _ in range(4)]}
    E = A.module(dims, maps, name="E=[2,3]+[3,4]")
    assert not E.is_indecomposable()                                   # genuinely a sum
    dvE = E.dimension_vector()
    assert dvE[1] == 0 and dvE[5] == 0                                 # non-proj & non-inj
    N = _interval(A, 2, 4, 5)
    _assert_interior(A, N)
    for M, X in ((E, N), (N, E)):
        euler_c = _euler_from_cartan(A, M, X)
        assert euler_c == _euler_arrow(A, M, X)
        assert A.hom(M, X) - A.ext(M, X, 1) == euler_c
        assert A.ext(M, X, 2) == 0


# --------------------------------------------------------------------------- #
# Pin 3: Auslander-Reiten -- Ext^1(M, tau M) != 0 for M non-projective indec
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_ar_ext_nonzero_kA5():
    """kA_5: M = [2,3] is non-projective indecomposable; tau M = [3,4] is ALSO interior
    (non-P, non-I, non-simple), and the AR sequence gives Ext^1(M, tau M) = 1 (ASS2006
    IV/V).  tau M is identified with the interval [3,4] by the exact iso certificate."""
    A = linear_path_algebra(5, field=QQ)
    M = _interval(A, 2, 3, 5)
    _assert_interior(A, M)
    tauM = M.tau()
    _assert_interior(A, tauM)                          # tau M is itself interior here
    assert tauM.is_isomorphic(_interval(A, 3, 4, 5))   # tau[2,3] = [3,4]
    assert A.ext(M, tauM, 1) >= 1                       # AR sequence exists => nonzero
    assert A.ext(M, tauM, 1) == 1


@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_ar_ext_nonzero_commutative_square():
    """The commutative square (non-hereditary, gl.dim 2): M = P_1/soc P_1 is a
    non-projective indecomposable and tau M = rad P_1 -- BOTH interior -- with
    Ext^1(M, tau M) = 1 (the AR sequence, ASS2006)."""
    S = _square()
    M = _p_mod_soc(S.projective(1))
    _assert_interior(S, M)
    tauM = M.tau()
    _assert_interior(S, tauM)
    assert tauM.is_isomorphic(S.projective(1).radical())   # tau(P_1/soc) = rad P_1
    assert S.ext(M, tauM, 1) >= 1                           # AR sequence exists => nonzero (literature)
    assert S.ext(M, tauM, 1) == 1                           # the engine value (selfcert refinement)


# --------------------------------------------------------------------------- #
# Pin 2 + 4: Tor duality anchor (both args interior) + balance
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_tor_duality_anchor_hereditary_kA5():
    """kA_5: a RIGHT interior M = [2,3] against a LEFT interior N = D([3,4]) (D of the
    interior right module [3,4] is an interior LEFT module).  dim Tor_n(M,N) =
    dim Ext^n(M, DN) degreewise; the only nonvanishing higher Tor is Tor_1 (hereditary),
    which fires (= 1) -- a non-vacuous pair."""
    A = linear_path_algebra(5, field=QQ)
    M = _interval(A, 2, 3, 5)
    N = _interval(A, 3, 4, 5).dualize()               # right [3,4] -> left D[3,4]
    _assert_interior(A, M)
    _assert_interior_left(A, N)
    tors = _tor_anchor(A, M, N, 3)
    assert tors == [0, 1, 0, 0]                        # nonvanishing, hereditary
    # The duality anchor (_tor_anchor) resolves M on BOTH sides (Tor and Ext^*(M,DN)),
    # so it rests on ONE engine path; resolving N instead (over A^op) is independent.
    assert tor_dims(A, M, N, 3, resolve="second") == tors


@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_tor_duality_anchor_commutative_square():
    """The commutative square: RIGHT interior M = rad P_1 against LEFT interior
    N = D(P_1/soc P_1).  dim Tor_0 = 2 (non-vacuous); the duality anchor holds
    degreewise."""
    S = _square()
    M = S.projective(1).radical()
    N = _p_mod_soc(S.projective(1)).dualize()
    _assert_interior(S, M)
    _assert_interior_left(S, N)
    tors = _tor_anchor(S, M, N, 3)
    assert tors == [2, 0, 0, 0]
    # Independent engine path: resolve N over A^op (the anchor resolves only M).
    assert tor_dims(S, M, N, 3, resolve="second") == tors


@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_tor_duality_anchor_selfinjective_higher_degrees():
    """A self-injective kZ_3/J^3 (projective = injective, gl.dim = infinity): RIGHT
    interior M = rad P_1 against LEFT interior N = D(rad P_1).  Higher Tor is periodic
    and FIRES -- dim Tor = [1,0,1,0,1] to degree 4 -- exercising the machinery past
    degree 1; the duality anchor holds at every degree."""
    B = _cyclic_nakayama(3, 3)
    assert B.is_selfinjective()
    M = B.projective(1).radical()
    N = B.projective(1).radical().dualize()
    _assert_interior(B, M)                             # non-projective => non-injective
    _assert_interior_left(B, N)
    tors = _tor_anchor(B, M, N, 4)
    assert tors == [1, 0, 1, 0, 1]


@pytest.mark.oracle_crossengine
def test_tor_balance_resolve_second_selfinjective():
    """Balance of Tor (Plan 29): resolving M vs resolving N over A^op computes the SAME
    groups.  Checked on the interior self-injective pair, where higher Tor is nonzero,
    so the two independent resolutions are genuinely exercised."""
    B = _cyclic_nakayama(3, 3)
    M = B.projective(1).radical()
    N = B.projective(1).radical().dualize()
    first = tor_dims(B, M, N, 4)
    second = tor_dims(B, M, N, 4, resolve="second")
    assert first == second == [1, 0, 1, 0, 1]
