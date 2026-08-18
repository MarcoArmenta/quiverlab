"""Complete resolution of A over A^e by splicing the minimal A^e-resolution with its
Nakayama-twisted A^e-dual (Plan 76 / R3; Bergh-Jorgensen JNCG 7 (2013)). GF(p);
self-injective scope; d.d=0 across the splice joint and in-window acyclicity are the
mandatory self-certs.

THE ARBITER IS ACYCLICITY, NOT d.d=0 (measured, Plan-75 lesson repeating).  The
negative differentials are carried as A^e-module maps in AMBIENT (A^e)^r coordinates,
where the composite d.d is computed WITHOUT reference to the corner tags -- so a wrong
Nakayama permutation pi leaves d.d=0 intact and is INVISIBLE to it (verified live below:
the pi=identity control on kZ3/J2 passes assert_dd_zero and still reports wrong dims).
What a wrong pi does break, loudly, is the CORNER TYPING of the dual half (the entries
stop living in eps_t . A^e . eps_s for their declared tags) and, downstream, exactness.
Both are asserted."""
import pytest

from quiverlab import NakayamaAlgebra, linear_path_algebra, truncated_polynomial
from quiverlab.engine.complete_resolution import complete_resolution
from quiverlab.errors import QuiverlabError
from quiverlab.fields import GF

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


@selfcert
@pytest.mark.parametrize("n,p", [(2, 32003), (2, 2), (3, 32003), (3, 3), (4, 32003)])
def test_dd_zero_and_acyclic_kxn(n, p):
    A = truncated_polynomial(n, field=GF(p))
    T = complete_resolution(A, 5, p)
    T.assert_dd_zero()                      # every joint, incl. d_0.d_1 and d_-1.d_0
    T.assert_acyclic(window=(-4, 4))


@selfcert
def test_dd_zero_multivertex_symmetric_nakayama():
    """The multi-vertex minimal engine must be fed UN-adapted (Plan-18 gotcha):
    unit_adapted() trips radical_basis / vertex detection on multi-vertex."""
    A = NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003))   # dim 6, symmetric
    T = complete_resolution(A, 4, 32003)
    T.assert_dd_zero()
    T.assert_acyclic(window=(-3, 3))
    assert T.term_rank(0) == T.term_rank(-1)     # T_-1 = dual of T_0 (swapped tags)


@selfcert
def test_kz3j2_nu_twist_REQUIRES_pi(monkeypatch):
    """M3 discriminating witness: kZ3/J2 is self-injective, NON-symmetric, vertex
    permutation pi = (1 2 3).  The pi-permuted dual is exact; forcing pi = identity
    (the plain swap) must FAIL loudly -- the value that discriminates the tag code."""
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))   # dim 6, non-symmetric
    T = complete_resolution(A, 4, 32003)                          # correct pi -> passes
    T.assert_dd_zero()
    T.assert_acyclic(window=(-3, 3))

    import quiverlab.engine.complete_resolution as cr
    monkeypatch.setattr(cr, "nakayama_permutation",
                        lambda A: {v: v for v in A.quiver.vertices})
    with pytest.raises(QuiverlabError):
        cr.complete_resolution(A, 4, 32003)


@selfcert
def test_dd_zero_alone_does_not_discriminate_pi(monkeypatch):
    """The measured negative result that justifies assert_acyclic (honest scope).

    With the corner-typing gate disabled, a WRONG pi still satisfies d.d=0 at every
    joint -- the composite is taken in ambient A^e coordinates, which carry no tag --
    and yet the reported cohomology is wrong.  Exactness is what sees it.  Pinning
    this stops a future refactor from 'simplifying' assert_acyclic away."""
    import quiverlab.engine.complete_resolution as cr
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))
    monkeypatch.setattr(cr, "nakayama_permutation",
                        lambda A: {v: v for v in A.quiver.vertices})
    T = cr.complete_resolution(A, 4, 32003, _skip_corner_gate=True)
    T.assert_dd_zero()                      # BLIND: passes with the wrong pi
    with pytest.raises(AssertionError):
        T.assert_acyclic(window=(-3, 3))    # sees it


@selfcert
def test_non_selfinjective_deferred():
    """M4: a Gorenstein-but-NOT-self-injective algebra (kA2, gldim 1) is DEFERRED --
    the D(P_n)-dual splice cannot build it -> loud refusal."""
    A = linear_path_algebra(2, field=GF(32003))
    with pytest.raises(QuiverlabError):
        complete_resolution(A, 4, 32003)


@selfcert
def test_not_selfinjective_radsq_refused():
    from quiverlab import Quiver
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}),
                          field=GF(32003))
    with pytest.raises(QuiverlabError):
        complete_resolution(A, 4, 32003)


@selfcert
def test_structure_constant_algebra_refused():
    """No quiver presentation -> no path-type basis -> loud refusal (never a guess)."""
    from quiverlab.core.algebra import Algebra
    A = truncated_polynomial(2, field=GF(32003))
    plain = Algebra(A.domain, A.T, A.unit)                 # presentation dropped
    with pytest.raises(QuiverlabError):
        complete_resolution(plain, 3, 32003)


@selfcert
def test_nakayama_permutation_witnesses():
    from quiverlab.engine.complete_resolution import nakayama_permutation
    sym = NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003))
    assert nakayama_permutation(sym) == {1: 1, 2: 2}              # symmetric -> identity
    non = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))
    assert nakayama_permutation(non) == {1: 2, 2: 3, 3: 1}        # live-verified 3-cycle
    assert not non.is_symmetric()


# ---------------------------------------------------------------------------
# Task I2 -- the Tate dimensions off the complete resolution
# ---------------------------------------------------------------------------

@lit
@pytest.mark.parametrize("n,p,expect", [
    (2, 32003, 1), (2, 2, 2), (3, 32003, 2), (3, 3, 3), (4, 32003, 3)])
def test_kxn_full_tate_ring_constant(n, p, expect):
    """k[x]/(x^n): dim HHhat^m = n-1 (char !| n) / n (char | n) for ALL m in Z.
    Derived from the explicit period-2 A^e-resolution (plan Live-verified facts)."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    A = truncated_polynomial(n, field=GF(p))
    dims = tate_cohomology_dims(A, 5, p)
    for m in range(-5, 6):
        assert dims[m] == expect, f"degree {m}"


@lit
def test_kxn_tate_hh0_differs_from_ordinary_hh0():
    """The headline HHhat^0 != HH^0 witness: the Tate degree-0 is the STABLE centre,
    a proper quotient of Z(A) = HH^0.  k[x]/(x^3) over GF(32003): 2 vs 3."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    A = truncated_polynomial(3, field=GF(32003))
    assert tate_cohomology_dims(A, 4, 32003)[0] == 2
    assert A.hochschild_cohomology(2).dims[0] == 3


@xeng
@pytest.mark.parametrize("n,p", [(3, 32003), (2, 2)])
def test_positive_agrees_with_ordinary_hh(n, p):
    """Bergh-Jorgensen threshold: HHhat^m == HH^m for m >= 1 (self-injective, d=0).
    The NEW splice engine against the SHIPPED bar/fast engine -- independent routes."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    A = truncated_polynomial(n, field=GF(p))
    tate = tate_cohomology_dims(A, 6, p)
    hh = A.hochschild_cohomology(6).dims
    for m in range(1, 7):
        assert tate[m] == hh[m], f"degree {m}"


@xeng
@pytest.mark.parametrize("spec", [(2, 3), (3, 2), (2, 2)])
def test_positive_agrees_multivertex_nakayama(spec):
    """Same threshold oracle on multi-vertex self-injective input (corner path),
    including the NON-symmetric kZ3/J2 whose pi is a 3-cycle."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    n, l = spec
    A = NakayamaAlgebra(n=n, l=l, cyclic=True, field=GF(32003))
    tate = tate_cohomology_dims(A, 4, 32003)
    hh = A.hochschild_cohomology(4).dims
    for m in range(1, 5):
        assert tate[m] == hh[m], f"degree {m}"


@selfcert
def test_symmetric_duality_self_consistency():
    """Bergh-Jorgensen, symmetric case (nu^2 = id): dim HHhat^{-j} == dim HHhat^{j-1}
    for all j -- a relation the engine is never told, computed on both halves."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    for A in (truncated_polynomial(3, field=GF(32003)),
              NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003))):
        assert A.is_symmetric()
        dims = tate_cohomology_dims(A, 4, 32003)
        for j in range(1, 5):
            assert dims[-j] == dims[j - 1], f"j={j} on {A}"


@selfcert
def test_nonsymmetric_negatives_are_not_the_mirror():
    """kZ3/J2 is Frobenius but NOT symmetric, so the symmetric mirror does NOT apply:
    the engine must not be silently mirroring the positive degrees."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))
    assert A.is_selfinjective() and not A.is_symmetric()
    dims = tate_cohomology_dims(A, 4, 32003)
    assert any(dims[-j] != dims[j - 1] for j in range(1, 5))


@selfcert
def test_tate_homology_frobenius_duality():
    """The HOMOLOGY duality holds for EVERY Frobenius algebra (Bergh-Jorgensen p.2):
    dim HHhat_n == dim HHhat_{-(n+1)} -- checked on the NON-symmetric witness too."""
    from quiverlab.engine.complete_resolution import tate_homology_dims
    for A in (truncated_polynomial(3, field=GF(32003)),
              NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))):
        dims = tate_homology_dims(A, 4, 32003)
        for n in range(0, 4):
            assert dims[n] == dims[-(n + 1)], f"n={n} on {A}"
