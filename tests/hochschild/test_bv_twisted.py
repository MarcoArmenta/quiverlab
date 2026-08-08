"""Plan 54 Task E -- the twisted-coefficient BV route (semisimple nu) over P52
twisted homology.

Delta on HH^*(A) for a NON-symmetric Frobenius algebra with semisimple Nakayama
automorphism (LZZ route), via P52's twisted Hochschild homology HH_*(A, {}_1A_nu),
the twisted Connes B_sigma (bv/twisted_connes.py) and the nu-twisted Frobenius
pairing (bv/transport.py::twisted_pairing_matrix). The flagship non-trivial-twist
instance is QuantumCI(q) (semisimple diagonal nu, NOT symmetric).

Certs exercised: (1) B_sigma DESCENDS to twisted homology (B of a cycle is a cycle
-- the paracyclic defect vanishes on homology; the loud descent self-cert inside
twisted_connes_class_matrix), (2) the perfect nu-twisted pairing (square +
invertible), (3) Delta^2 = 0, (4) the derived bracket == the INDEPENDENT
Gerstenhaber bracket in-window (the decisive cross-engine arbiter, which also
selects the twist direction).

DEVIATION recorded (vs spec MAJOR-2): for QuantumCI the two twist directions do NOT
give equal twisted-homology dims -- {}_1A_nu matches dim HH^* while {}_1A_{nu^{-1}}
does not -- so the perfect-pairing dim-match ALREADY rejects the wrong direction,
and the bracket arbiter then certifies the surviving one. Direction determination
is therefore even sharper than the spec assumed (dim-match + arbiter agree).
"""
import numpy as np
import pytest

import quiverlab as ql
from quiverlab.families import QuantumCI
from quiverlab.hochschild.bv.hypothesis import classify_bv
from quiverlab.hochschild.bv.twist import nu_inverse, twisted_homology_quotient
from quiverlab.hochschild.bv.twisted_connes import twisted_connes_class_matrix
from quiverlab.invariants.frobenius import nakayama_automorphism_generic


def _rank(M, p):
    from quiverlab.engine.coxeter import rref_mod_p
    A = np.array([[int(x) for x in r] for r in M], dtype=np.int64)
    return len(rref_mod_p(A % p, p)[1]) if A.size else 0


# ---------------------------------------------------------------------------
# routing: QuantumCI(q!=+-1) is a NON-symmetric semisimple-nu Frobenius algebra
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
def test_quantumci_routes_semisimple():
    A = QuantumCI(2, field=ql.GF(5))
    assert not A.is_symmetric() and A.is_frobenius()
    h = classify_bv(A)
    assert h.route == "semisimple" and h.nu_semisimple is True


# ---------------------------------------------------------------------------
# the twisted route computes and self-certifies (Delta^2 = 0 + descent + pairing)
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
def test_twisted_route_delta_squared_zero():
    A = QuantumCI(2, field=ql.GF(5))
    bv = A.bv_operator(3)
    assert "semisimple Nakayama automorphism" in bv.hypothesis
    assert "twist {}_1A_nu" in bv.hypothesis and bv.nakayama["inner"] is False
    assert bv.hh_dims == [2, 2, 1, 0]
    p = 5
    for n in range(2, bv.top + 1):
        Dn = np.array([[int(x) for x in r] for r in bv.matrices[n]], dtype=np.int64)
        Dp = np.array([[int(x) for x in r] for r in bv.matrices[n - 1]], dtype=np.int64)
        if Dn.size and Dp.size:
            assert not np.any((Dp @ Dn) % p), f"Delta^2 != 0 at {n}"


# ---------------------------------------------------------------------------
# the descent self-cert (MAJOR 4): B_sigma carries twisted cycles to cycles for
# the correct twist; twisted_connes_class_matrix raises loudly otherwise.
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
def test_twisted_connes_descends_on_correct_twist():
    A = QuantumCI(2, field=ql.GF(5))
    AU = A.unit_adapted()
    p = 5
    nu = nakayama_automorphism_generic(AU)
    nu_int = [[int(nu[i][j]) % p for j in range(AU.dim)] for i in range(AU.dim)]
    tw = twisted_homology_quotient(AU, nu_int, 3)
    twist = np.array(nu_int, dtype=np.int64)
    # descends (no raise) for every in-window degree
    for n in range(3):
        B = twisted_connes_class_matrix(AU, twist, tw, n, p)
        assert B.shape == (tw[n + 1]["dim"], tw[n]["dim"])


# ---------------------------------------------------------------------------
# twist-direction discrimination: only {}_1A_nu yields a perfect pairing (its
# twisted homology matches dim HH^*); {}_1A_{nu^{-1}} does not.
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
def test_twist_direction_discriminated():
    A = QuantumCI(2, field=ql.GF(5))
    AU = A.unit_adapted()
    p = 5
    nu = nakayama_automorphism_generic(AU)
    nui = nu_inverse(nu, A.domain)
    nu_int = [[int(nu[i][j]) % p for j in range(AU.dim)] for i in range(AU.dim)]
    nui_int = [[int(nui[i][j]) % p for j in range(AU.dim)] for i in range(AU.dim)]
    top = 2                                # top=2 already discriminates (nu=[2,2,1]
    d_nu = [twisted_homology_quotient(AU, nu_int, top)[n]["dim"] for n in range(top + 1)]
    d_nui = [twisted_homology_quotient(AU, nui_int, top)[n]["dim"] for n in range(top + 1)]
    hh = list(A.hochschild_cohomology(top).dims)
    assert d_nu == hh                      # {}_1A_nu is the perfect-pairing direction
    assert d_nui != hh                     # {}_1A_{nu^{-1}} is not (rejected)


# ---------------------------------------------------------------------------
# the ARBITER (oracle_crossengine): the QuantumCI derived bracket == the
# independent Gerstenhaber bracket in-window over ODD primes. THIS is the
# decisive correctness proof for the twisted route (no Frobenius input in the
# independent bracket).
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
def test_twisted_seven_term_relation():
    # (7T): Delta is a differential operator of order <= 2 wrt cup, on the twisted
    # route too. top=3 QuantumCI: the only in-window triple is (1,1,1).
    from quiverlab.hochschild.bv.bracket import seven_term_residual
    A = QuantumCI(2, field=ql.GF(5))
    top = 3
    bv = A.bv_operator(top)
    cup = A.cup_products(top)
    for pp in range(1, top + 1):
        for qq in range(1, top + 1):
            for rr in range(1, top + 1):
                if pp + qq + rr <= top:
                    assert seven_term_residual(A, bv, cup, pp, qq, rr) == 0, \
                        (pp, qq, rr)


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("q,prime", [(2, 5), (2, 32003), (3, 5)])
def test_quantumci_derived_bracket_equals_independent(q, prime):
    A = QuantumCI(q, field=ql.GF(prime))
    if classify_bv(A).route != "semisimple":
        pytest.skip("collapses to symmetric at this prime")
    top = 3
    bv = A.bv_operator(top)
    assert bv.bracket_check["agrees"] is True
    assert bv.bracket_check["pairs_checked"] >= 1
    derived = bv.derived_bracket
    independent = A.gerstenhaber_brackets(top)
    p = A.domain.p
    for key, dt in derived.tables.items():
        it = independent.tables[key]
        assert dt.dims == it.dims, key
        dl, dr, dout = dt.dims
        for k in range(dout):
            for i in range(dl):
                for j in range(dr):
                    assert int(dt.constants[k][i][j]) % p == \
                        int(it.constants[k][i][j]) % p, (key, k, i, j)
