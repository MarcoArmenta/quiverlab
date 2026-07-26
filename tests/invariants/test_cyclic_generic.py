"""Generic-Domain cyclic homology (Plan 19): GF(p) parity with the engine,
mixed-complex identities over QQ, the char-0 Connes lambda-complex second
model (Loday Thm 2.1.5), and closed-form pins."""
import itertools

import pytest

from quiverlab import CC, GF, Quiver, truncated_polynomial
from quiverlab.families import QuantumCI
from quiverlab.fields import QQ, linalg

pytestmark = [
    pytest.mark.oracle_crossengine,
    pytest.mark.oracle_selfcert,
    pytest.mark.oracle_literature,
]


def _k(field):
    return Quiver([1], {}).algebra(relations=[], field=field)


def _kxk(field):
    return Quiver([1, 2], {}).algebra(relations=[], field=field)


def test_hc_generic_matches_engine_over_gfp():
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.cyclic import cyclic_homology_dims as engine_hc
    from quiverlab.hochschild.cyclic import cyclic_homology_dims as generic_hc
    cases = ((truncated_polynomial(2, field=GF(5)), 3),
             (truncated_polynomial(3, field=GF(3)), 3),
             (QuantumCI(2, field=GF(5)), 2))
    for A, top in cases:
        p = A.domain.p
        want = engine_hc(to_engine(A.unit_adapted()), top, primes=(p,))[p]
        got = generic_hc(A, top)
        assert got.kind == "HC_"
        assert got.dims == [int(d) for d in want]


def _matmul(X, Y, dom):
    """(a x b) @ (b x c) over dom, list-of-rows convention."""
    if not X or not Y:
        return []
    c = len(Y[0])
    out = [[dom.zero()] * c for _ in range(len(X))]
    for i, Xi in enumerate(X):
        for j, xv in enumerate(Xi):
            if dom.is_zero(xv):
                continue
            Yj = Y[j]
            for k in range(c):
                if not dom.is_zero(Yj[k]):
                    out[i][k] = dom.add(out[i][k], dom.mul(xv, Yj[k]))
    return out


def test_mixed_complex_identities_over_qq():
    from quiverlab.hochschild.bar import boundary_matrix
    from quiverlab.hochschild.cyclic import connes_B_matrix
    A = truncated_polynomial(2, field=QQ).unit_adapted()
    dom = A.domain
    b = {k: boundary_matrix(A, k, 10 ** 6)[0] for k in range(1, 5)}
    B = {k: connes_B_matrix(A, k, 10 ** 6)[0] for k in range(0, 4)}

    def zero(M):
        return all(dom.is_zero(x) for row in M for x in row)

    for n in range(0, 3):
        if n >= 1:
            assert zero(_matmul(b[n], b[n + 1], dom)), f"b^2 != 0 at {n}"
        assert zero(_matmul(B[n + 1], B[n], dom)), f"B^2 != 0 at {n}"
        bB = _matmul(b[n + 1], B[n], dom)
        if n >= 1:
            Bb = _matmul(B[n - 1], b[n], dom)
            S = [[dom.add(bB[i][j], Bb[i][j]) for j in range(len(bB[0]))]
                 for i in range(len(bB))]
        else:
            S = bB
        assert zero(S), f"bB + Bb != 0 at {n}"


def _lambda_hc(A, top):
    """Char-0 SECOND MODEL: HC_n = H_n(C^lambda) where C^lambda_n =
    A^{(x)(n+1)} / im(1 - t), t = (-1)^n cyclic rotation, b = the FULL cyclic
    Hochschild boundary (n+1 faces, wrap-around last).  Valid over fields
    containing Q (Loday, Cyclic Homology, Thm 2.1.5).  Independent of
    hochschild/bar.py AND hochschild/cyclic.py: unnormalized chains, no
    unit-adaptation, quotient model instead of bicomplex."""
    dom = A.domain
    m = A.dim

    def b_cols(n):
        """Columns of b: C_n -> C_{n-1} as vectors (rank is transpose-invariant)."""
        colws = list(itertools.product(range(m), repeat=n + 1))
        rowws = list(itertools.product(range(m), repeat=n))
        ridx = {w: i for i, w in enumerate(rowws)}
        cols = []
        for w in colws:
            v = [dom.zero()] * len(rowws)
            for i in range(n + 1):
                if i < n:
                    prod = A.T[w[i]][w[i + 1]]
                    keys = [w[:i] + (t,) + w[i + 2:] for t in range(m)]
                else:
                    prod = A.T[w[n]][w[0]]
                    keys = [(t,) + w[1:n] for t in range(m)]
                for t, cf in enumerate(prod):
                    if dom.is_zero(cf):
                        continue
                    val = cf if i % 2 == 0 else dom.neg(cf)
                    r = ridx[keys[t]]
                    v[r] = dom.add(v[r], val)
            cols.append(v)
        return cols

    def v_gens(n):
        """Spanning set of im(1 - t) on C_n."""
        ws = list(itertools.product(range(m), repeat=n + 1))
        idx = {w: i for i, w in enumerate(ws)}
        out = []
        for w in ws:
            rot = (w[-1],) + w[:-1]
            v = [dom.zero()] * len(ws)
            v[idx[w]] = dom.add(v[idx[w]], dom.one())
            s = dom.neg(dom.one()) if n % 2 == 0 else dom.one()   # -(-1)^n
            v[idx[rot]] = dom.add(v[idx[rot]], s)
            out.append(v)
        return out

    qdim, rkbar = {}, {0: 0}
    for n in range(top + 2):
        vg = v_gens(n)
        qdim[n] = m ** (n + 1) - (linalg.rank(vg, dom) if vg else 0)
        if n >= 1:
            vprev = v_gens(n - 1)
            rprev = linalg.rank(vprev, dom) if vprev else 0
            rkbar[n] = linalg.rank(b_cols(n) + vprev, dom) - rprev
    return [qdim[n] - rkbar[n] - rkbar[n + 1] for n in range(top + 1)]


def test_hc_lambda_complex_second_model_qq():
    for A, top in ((_k(QQ), 4), (_kxk(QQ), 3),
                   (truncated_polynomial(2, field=QQ), 3)):
        assert A.cyclic_homology(top).dims == _lambda_hc(A, top)


def test_hc_closed_forms_char0_and_gf4():
    for field in (QQ, GF(4)):
        assert _k(field).cyclic_homology(4).dims == [1, 0, 1, 0, 1]
        assert _kxk(field).cyclic_homology(4).dims == [2, 0, 2, 0, 2]


def test_hc_dual_numbers_cc_matches_qq():
    # dims are invariant under field extension QQ -> CC (flat base change)
    got = truncated_polynomial(2, field=CC).cyclic_homology(2)
    assert got.kind == "HC_"
    assert got.dims == truncated_polynomial(2, field=QQ).cyclic_homology(2).dims
