"""Plan 54 Task A -- the BV hypothesis gate (semisimple-nu / symmetric anchor).

Routing (per instance, in order): symmetric -> semisimple -> self-injective
Nakayama non-semisimple (BIKLZ-blocked, loud) -> loud refusal.
"""
import pytest

import quiverlab as ql
from quiverlab.hochschild.bv.hypothesis import classify_bv, nu_is_semisimple

pytestmark = pytest.mark.oracle_selfcert


def test_symmetric_kxy():
    # k[x,y]/(x^2,y^2) = QuantumCI(-1): symmetric, nu = id (inner).
    A = ql.QuantumCI(-1, field=ql.GF(5))
    h = classify_bv(A)
    assert h.applies is True
    assert h.route == "symmetric"
    assert h.nu_inner is True
    assert h.refusal is None
    assert "Tradler" in h.label


def test_semisimple_quantumci():
    # QuantumCI(q=2) over GF(5): Frobenius, NOT symmetric, diagonal semisimple nu.
    A = ql.QuantumCI(2, field=ql.GF(5))
    h = classify_bv(A)
    assert h.applies is True
    assert h.route == "semisimple"
    assert h.nu_semisimple is True
    assert h.nu_inner is False
    assert "semisimple" in h.label.lower()


def test_kx_xp_is_symmetric_not_biklz():
    # k[x]/(x^p) over GF(p): char | N, but it is SYMMETRIC (nu = id) -> symmetric
    # route, NOT the BIKLZ refusal (the e=1 special case is served by route 1).
    for p in (2, 3, 5):
        A = ql.truncated_polynomial(p, field=ql.GF(p))
        h = classify_bv(A)
        assert h.applies is True, (p, h.refusal)
        assert h.route == "symmetric", p


def test_biklz_blocked_nonsemisimple_selfinjective_nakayama():
    # kZ_2/J^4 over GF(2): Frobenius, self-injective Nakayama, NOT symmetric,
    # char | ord(nu) so nu is non-semisimple -> the general BIKLZ construction is
    # out of v1 scope: loud refusal.
    A = ql.NakayamaAlgebra(n=2, l=4, cyclic=True, field=ql.GF(2))
    h = classify_bv(A)
    assert h.applies is False
    assert h.route is None
    assert h.nu_semisimple is False
    assert "general" in h.refusal
    assert "is NOT implemented in v1" in h.refusal


def test_non_frobenius_refused():
    # kA_2 (linear path algebra): not Frobenius -> loud refusal.
    A = ql.linear_path_algebra(2, field=ql.GF(5))
    h = classify_bv(A)
    assert h.applies is False
    assert h.route is None
    assert "not Frobenius" in h.refusal


def test_nu_is_semisimple_direct():
    dom = ql.GF(5)
    # diagonal(1,3,2,1): distinct-ish, squarefree minpoly -> semisimple
    diag = [[1, 0, 0, 0], [0, 3, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]
    assert nu_is_semisimple(diag, dom) is True
    # a 2x2 Jordan block for eigenvalue 1 over GF(2): [[1,1],[0,1]] -> NOT semisimple
    dom2 = ql.GF(2)
    jordan = [[1, 1], [0, 1]]
    assert nu_is_semisimple(jordan, dom2) is False
