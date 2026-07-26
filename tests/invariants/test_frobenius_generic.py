"""Generic-Domain Frobenius / Nakayama / symmetry (Plan 19): socle-criterion
knowns over QQ / CC / GF(4), self-certification of the form and nu, the
inner-automorphism symmetry semantics, GF(p) parity with the engine.

Char-2 trap (global constraint): QuantumCI(q) degenerates when q == 0 in the
field — expectations here are per-field, never blanket."""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.families import QuantumCI
from quiverlab.families.zoo import build_from_record, load_catalog
from quiverlab.fields import QQ

pytestmark = [
    pytest.mark.oracle_literature,
    pytest.mark.oracle_selfcert,
    pytest.mark.oracle_crossengine,
]


def _rec(name):
    return next(r for r in load_catalog() if r.get("name") == name)


def _commxy(field):
    """k[x,y]/(x^2, y^2, xy): socle span(x, y) is 2-dim => NOT Frobenius."""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x^2", "y^2", "x*y", "y*x"], field=field)


def test_frobenius_knowns_any_domain():
    for f in (QQ, GF(4)):
        assert truncated_polynomial(4, field=f).is_frobenius() is True
        assert linear_path_algebra(2, field=f).is_frobenius() is False
        assert _commxy(f).is_frobenius() is False
    assert truncated_polynomial(2, field=CC).is_symmetric() is True


def test_exterior_frobenius_not_symmetric_qq():
    # Lambda(x, y): nu = diag(1, -1, -1, 1) is NOT inner (trivial vertex
    # permutation but no invertible u with nu(a)u = ua) — the case that
    # separates the inner test from the permutation shortcut.
    E2 = QuantumCI(1, field=QQ)
    assert E2.is_frobenius() is True
    assert E2.is_symmetric() is False


def test_commutative_ci_symmetric_qq():
    assert QuantumCI(-1, field=QQ).is_symmetric() is True


def test_quantum_ci_nakayama_qq():
    A = QuantumCI(2, field=QQ)
    dom = A.domain
    assert A.is_frobenius() is True
    assert A.is_symmetric() is False
    N = A.nakayama_automorphism()
    m = A.dim
    ix = A.basis_labels.index("x")
    iy = A.basis_labels.index("y")
    for i in range(m):
        for j in range(m):
            if i != j:
                assert dom.is_zero(N[i][j]), "nu of the quantum CI is diagonal"
    # quiverlab's QuantumCI(q) relation is x*y + q*y*x (so yx = -xy/q): the
    # Nakayama automorphism is nu = diag(1, -q, -1/q, 1) on (e, x, y, xy) —
    # NOT the engine docstring's yx - q*xy convention. Check: with the
    # xy-dual form, lambda(y nu(x)) = -q*lambda(yx) = -q*(-1/q) = 1 = lambda(xy).
    mq = dom.neg(dom.coerce(2))
    mqinv = dom.inv(mq)
    assert ((dom.eq(N[ix][ix], mq) and dom.eq(N[iy][iy], mqinv)) or
            (dom.eq(N[ix][ix], mqinv) and dom.eq(N[iy][iy], mq)))


def test_cn_3_2_frobenius_not_symmetric():
    # multi-vertex: kZ_3/rad^2 — Nakayama permutation is the 3-cycle, so
    # Frobenius but NOT symmetric (inner autos fix vertex classes)
    for f in (QQ, GF(4)):
        A = build_from_record(_rec("cn_3_2"), field=f)
        assert A.is_frobenius() is True
        assert A.is_symmetric() is False


def test_nakayama_self_certifies_across_domains():
    """The returned form and nu satisfy their defining equations EXACTLY:
    lambda(f_i f_j) = G[i][j] nondegenerate, lambda(ab) = lambda(b nu(a)),
    nu multiplicative, nu(1) = 1. No oracle needed — the axioms are the gate."""
    from quiverlab.fields import linalg
    from quiverlab.invariants.frobenius import (frobenius_form_generic,
                                                nakayama_automorphism_generic)
    zoo = [truncated_polynomial(3, field=QQ),
           truncated_polynomial(2, field=GF(4)),
           QuantumCI(2, field=QQ),
           QuantumCI(1, field=QQ),
           build_from_record(_rec("cn_3_2"), field=QQ)]
    for A in zoo:
        dom = A.domain
        m = A.dim
        lam, G = frobenius_form_generic(A)
        assert linalg.rank(G, dom) == m
        N = nakayama_automorphism_generic(A)

        def nu(vec):
            out = [dom.zero()] * m
            for j, c in enumerate(vec):
                if dom.is_zero(c):
                    continue
                for i in range(m):
                    out[i] = dom.add(out[i], dom.mul(c, N[i][j]))
            return out

        def lam_of(vec):
            acc = dom.zero()
            for t, c in enumerate(vec):
                acc = dom.add(acc, dom.mul(lam[t], c))
            return acc

        for a in range(m):
            fa = A._basis_vec(a)
            nua = nu(fa)
            for b in range(m):
                fb = A._basis_vec(b)
                # lambda(ab) = lambda(b nu(a))
                assert dom.eq(lam_of(A.multiply(fa, fb)),
                              lam_of(A.multiply(fb, nua)))
                # nu(ab) = nu(a) nu(b)
                got = nu(A.multiply(fa, fb))
                want = A.multiply(nua, nu(fb))
                assert all(dom.eq(x, y) for x, y in zip(got, want))
        one = nu(A.unit)
        assert all(dom.eq(x, y) for x, y in zip(one, A.unit))


def test_gfp_parity_with_engine():
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.coxeter import is_frobenius as eng_frob
    from quiverlab.invariants.frobenius import (is_frobenius_generic,
                                                is_symmetric_generic)
    zoo = [truncated_polynomial(3, field=GF(5)),
           QuantumCI(2, field=GF(5)),
           QuantumCI(1, field=GF(7)),
           _commxy(GF(5)),
           linear_path_algebra(2, field=GF(5))]
    for A in zoo:
        p = A.domain.p
        want = bool(eng_frob(to_engine(A.unit_adapted()), p))
        assert is_frobenius_generic(A) is want, repr(A)
    # symmetry parity where the engine's nu = id semantics and the inner
    # semantics provably coincide (nu diagonal / commutative cases)
    assert is_symmetric_generic(truncated_polynomial(3, field=GF(5))) is True
    assert is_symmetric_generic(QuantumCI(2, field=GF(5))) is False
    assert truncated_polynomial(3, field=GF(5)).is_symmetric() is True
    assert QuantumCI(2, field=GF(5)).is_symmetric() is False
