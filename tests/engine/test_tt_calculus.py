"""The Tamarkin-Tsygan calculus substrate (PLAN item A step 3 / item D):
the cup product on HH^* and the cap-product module action HH^* (x) HH_* -> HH_*.

This is the substrate the *genuine* Auslander-Reiten action Theta is built on
(Armenta, arXiv:1908.02255: Theta = capping with the inverse-dualizing class).

Acceptance gates (all on the normalized bar complex, the reference oracle):

  (A) cochain/chain-level ("on the nose") identities -- these are the rigorous proof
      that the products DESCEND to (co)homology:
        * cup is associative;
        * cup Leibniz  δ(f ⌣ g) = δf ⌣ g + (-1)^p f ⌣ δg        (⌣ descends);
        * cup unit 1_A acts as the identity;
        * cap unit 1_A acts as the identity;
        * cap module law  f ∩ (g ∩ z) = (g ⌣ f) ∩ z;
        * cap Leibniz  b(f ∩ z) = (-1)^{p+1}(δf ∩ z) + (-1)^p (f ∩ b z)  (∩ descends).
  (B) induced (class-level) structure:
        * the representative machinery reproduces the engine's HH^* / HH_* dims;
        * HH^* is a graded-commutative unital ring; the cap action is a unital module;
        * concrete regression pins: HH^0(A) cup = the multiplication of Z(A); for
          k[x]/(x^2) the degree-1 class squares to zero while the degree-2 class does
          not (the known HH^*(k[x]/x^2) = k[u,v]/(u^2), deg u=1, deg v=2).
"""
import numpy as np
import pytest

# hanlab sys.path shim dropped: quiverlab uses absolute package imports.
from quiverlab.engine.hh_engine import (
    cn_basis,
    differential_matrix,
    truncated_polynomial,
    two_gen_local,
    hochschild_homology_dims,
)
from quiverlab.engine.scan3 import (
    cochain_basis,
    coboundary_matrix,
    quantum_ci,
    hochschild_cohomology_dims,
)
from quiverlab.engine import tt_calculus as tt

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims        # dim HH_n(A)
cohomology_dims = hochschild_cohomology_dims    # dim HH^n(A)
PRIME = 32003
P = PRIME


# --- panel of small algebras ------------------------------------------------
def _panel():
    return [
        truncated_polynomial(2),                                    # k[x]/(x^2)
        truncated_polynomial(3),                                    # k[x]/(x^3)
        two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "kxy"),
        quantum_ci(2),                                              # non-monomial
    ]


def _delta(alg, n, vec):
    idx = {g: i for i, g in enumerate(cochain_basis(alg, n + 1))}
    return coboundary_matrix(alg, n, cochain_basis(alg, n), idx) @ vec


def _bdy(alg, n, vec):
    idx = {g: i for i, g in enumerate(cn_basis(alg, n - 1))}
    return differential_matrix(alg, n, cn_basis(alg, n), idx) @ vec


def _basis_cochains(alg, n):
    L = len(cochain_basis(alg, n))
    for i in range(L):
        v = np.zeros(L, dtype=np.int64); v[i] = 1; yield v


def _basis_chains(alg, n):
    L = len(cn_basis(alg, n))
    for i in range(L):
        v = np.zeros(L, dtype=np.int64); v[i] = 1; yield v


# ===========================================================================
# (A) on-the-nose identities -- the descent oracles
# ===========================================================================
@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_cup_associative_on_cochains(alg):
    """(f ⌣ g) ⌣ h = f ⌣ (g ⌣ h) for all degree-1 basis cochains (exact)."""
    for f in _basis_cochains(alg, 1):
        for g in _basis_cochains(alg, 1):
            for h in _basis_cochains(alg, 1):
                lhs = tt.cup_cochain(alg, 2, 1, tt.cup_cochain(alg, 1, 1, f, g), h)
                rhs = tt.cup_cochain(alg, 1, 2, f, tt.cup_cochain(alg, 1, 1, g, h))
                assert np.array_equal(lhs, rhs)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
@pytest.mark.parametrize("pq", [(1, 1), (1, 2), (2, 1)])
def test_cup_leibniz(alg, pq):
    """δ(f ⌣ g) = δf ⌣ g + (-1)^p f ⌣ δg  (⌣ is a chain map => descends to HH^*)."""
    p, q = pq
    sign = (-1) ** p
    for f in _basis_cochains(alg, p):
        for g in _basis_cochains(alg, q):
            lhs = _delta(alg, p + q, tt.cup_cochain(alg, p, q, f, g)) % P
            t1 = tt.cup_cochain(alg, p + 1, q, _delta(alg, p, f), g)
            t2 = tt.cup_cochain(alg, p, q + 1, f, _delta(alg, q, g))
            rhs = (t1 + sign * t2) % P
            assert np.array_equal(lhs % P, rhs % P)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_cup_unit_on_cochains(alg):
    """1_A in C^0 is a two-sided cup unit:  1 ⌣ g = g  and  g ⌣ 1 = g."""
    u = tt.unit_cochain(alg)
    for q in (1, 2):
        for g in _basis_cochains(alg, q):
            assert np.array_equal(tt.cup_cochain(alg, 0, q, u, g), g)
            assert np.array_equal(tt.cup_cochain(alg, q, 0, g, u), g)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_cap_unit_on_chains(alg):
    """1_A in C^0 acts as the identity:  1 ∩ z = z  for every chain z."""
    u = tt.unit_cochain(alg)
    for n in (1, 2, 3):
        for z in _basis_chains(alg, n):
            assert np.array_equal(tt.cap_cochain(alg, 0, n, u, z), z)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
@pytest.mark.parametrize("pqn", [(1, 1, 3), (1, 2, 3), (2, 1, 3)])
def test_cap_module_law(alg, pqn):
    """f ∩ (g ∩ z) = (g ⌣ f) ∩ z  -- HH_* is a module over the cup ring (exact)."""
    p, q, n = pqn
    for f in _basis_cochains(alg, p):
        for g in _basis_cochains(alg, q):
            for z in _basis_chains(alg, n):
                lhs = tt.cap_cochain(alg, p, n - q, f, tt.cap_cochain(alg, q, n, g, z))
                rhs = tt.cap_cochain(alg, p + q, n, tt.cup_cochain(alg, q, p, g, f), z)
                assert np.array_equal(lhs, rhs)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
@pytest.mark.parametrize("pn", [(1, 3), (2, 3)])
def test_cap_leibniz(alg, pn):
    """b(f ∩ z) = (-1)^{p+1}(δf ∩ z) + (-1)^p (f ∩ b z)  (∩ descends to HH_*)."""
    p, n = pn
    s1 = (-1) ** (p + 1)
    s2 = (-1) ** p
    for f in _basis_cochains(alg, p):
        for z in _basis_chains(alg, n):
            lhs = _bdy(alg, n - p, tt.cap_cochain(alg, p, n, f, z)) % P
            t1 = tt.cap_cochain(alg, p + 1, n, _delta(alg, p, f), z)
            t2 = tt.cap_cochain(alg, p, n - 1, f, _bdy(alg, n, z))
            rhs = (s1 * t1 + s2 * t2) % P
            assert np.array_equal(lhs % P, rhs % P)


def test_cap_degree_guard():
    """The cap product needs p <= n; otherwise ValueError."""
    A = truncated_polynomial(2)
    f = np.zeros(len(cochain_basis(A, 2)), dtype=np.int64)
    z = np.zeros(len(cn_basis(A, 1)), dtype=np.int64)
    with pytest.raises(ValueError):
        tt.cap_cochain(A, 2, 1, f, z)


# ===========================================================================
# (B) induced (class-level) structure
# ===========================================================================
@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_class_machinery_matches_engine_dims(alg):
    """The representative quotients reproduce the engine's HH^* and HH_* dims."""
    N = 4
    coh = cohomology_dims(alg, N, primes=(P,))[P]
    hom = homology_dims(alg, N, primes=(P,))[P]
    for n in range(N + 1):
        assert tt.cohomology_classes(alg, n, P).dim == coh[n]
        assert tt.homology_classes(alg, n, P).dim == hom[n]


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_cup_unit_acts_as_identity_on_classes(alg):
    """The class of 1_A in HH^0 is a unit for the induced cup product."""
    for q in (1, 2):
        C, d0, dq, dout = tt.cup_product_matrix(alg, 0, q, P)
        if dq == 0:
            continue
        assert dout == dq
        # the unit class: coordinates of 1_A among the HH^0 reps
        H0 = tt.cohomology_classes(alg, 0, P)
        u = H0.coords(tt.unit_cochain(alg))
        # sum_i u_i C[:, i, j] must be the j-th basis vector of HH^q
        for j in range(dq):
            col = sum(int(u[i]) * C[:, i, j] for i in range(d0)) % P
            e = np.zeros(dout, dtype=np.int64); e[j] = 1
            assert np.array_equal(col, e)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
@pytest.mark.parametrize("pq", [(1, 1), (1, 2), (2, 1)])
def test_cup_graded_commutative_on_classes(alg, pq):
    """rep_i ⌣ rep_j = (-1)^{pq} rep_j ⌣ rep_i  in HH^{p+q}."""
    p, q = pq
    Cpq, dp, dq, _ = tt.cup_product_matrix(alg, p, q, P)
    Cqp, _, _, _ = tt.cup_product_matrix(alg, q, p, P)
    s = (-1) ** (p * q) % P
    for i in range(dp):
        for j in range(dq):
            assert np.array_equal(Cpq[:, i, j] % P, (s * Cqp[:, j, i]) % P)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_cup_associative_on_classes(alg):
    """(a ⌣ b) ⌣ c = a ⌣ (b ⌣ c) on HH^1 x HH^1 x HH^1 -> HH^3 classes."""
    C11, d1, _, d2 = tt.cup_product_matrix(alg, 1, 1, P)
    C21, _, _, d3a = tt.cup_product_matrix(alg, 2, 1, P)
    C11b, _, _, _ = tt.cup_product_matrix(alg, 1, 1, P)
    C12, _, _, d3b = tt.cup_product_matrix(alg, 1, 2, P)
    assert d3a == d3b
    for i in range(d1):
        for j in range(d1):
            for k in range(d1):
                # (i ⌣ j) ⌣ k
                left = np.zeros(d3a, dtype=np.int64)
                for a in range(d2):
                    left = (left + int(C11[a, i, j]) * C21[:, a, k]) % P
                # i ⌣ (j ⌣ k)
                right = np.zeros(d3b, dtype=np.int64)
                for b in range(d2):
                    right = (right + int(C11b[b, j, k]) * C12[:, i, b]) % P
                assert np.array_equal(left, right)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
def test_cap_unit_identity_on_homology_classes(alg):
    """The class of 1_A in HH^0 caps as the identity on HH_n."""
    H0 = tt.cohomology_classes(alg, 0, P)
    u = H0.coords(tt.unit_cochain(alg))
    for n in (1, 2, 3):
        C, d0, dn, dout = tt.cap_product_matrix(alg, 0, n, P)
        assert dout == dn
        for j in range(dn):
            col = sum(int(u[i]) * C[:, i, j] for i in range(d0)) % P
            e = np.zeros(dout, dtype=np.int64); e[j] = 1
            assert np.array_equal(col, e)


@pytest.mark.parametrize("alg", _panel(), ids=lambda a: a.name)
@pytest.mark.parametrize("pqn", [(1, 1, 3), (1, 2, 3)])
def test_cap_module_law_on_classes(alg, pqn):
    """f ∩ (g ∩ z) = (g ⌣ f) ∩ z, induced on HH^p, HH^q, HH_n classes."""
    p, q, n = pqn
    Ccup, dq_, dp_, dqp = tt.cup_product_matrix(alg, q, p, P)
    Ccap_g, _, dn, dnq = tt.cap_product_matrix(alg, q, n, P)
    Ccap_f, dpf, dnq2, dnpq = tt.cap_product_matrix(alg, p, n - q, P)
    Ccap_gf, dgf, dn2, dnpq2 = tt.cap_product_matrix(alg, p + q, n, P)
    assert dnpq == dnpq2 and dn == dn2
    Hp = tt.cohomology_classes(alg, p, P).dim
    Hq = tt.cohomology_classes(alg, q, P).dim
    for fi in range(Hp):
        for gj in range(Hq):
            for zk in range(dn):
                # f ∩ (g ∩ z)
                gz = Ccap_g[:, gj, zk]
                lhs = np.zeros(dnpq, dtype=np.int64)
                for a in range(dnq):
                    lhs = (lhs + int(gz[a]) * Ccap_f[:, fi, a]) % P
                # (g ⌣ f) ∩ z
                gf = Ccup[:, gj, fi]
                rhs = np.zeros(dnpq, dtype=np.int64)
                for b in range(dgf):
                    rhs = (rhs + int(gf[b]) * Ccap_gf[:, b, zk]) % P
                assert np.array_equal(lhs, rhs)


# ===========================================================================
# concrete regression pins
# ===========================================================================
def test_HH0_cup_is_center_multiplication():
    """HH^0(A) = Z(A) and the cup product there is the multiplication of Z(A).

    Since C^0 = A, a HH^0 representative IS an element of A; cup on HH^0 is the
    A-product.  For k[x]/(x^2) the center is all of A, so the table is A's own:
    1·1=1, 1·x=x, x·x=0."""
    A = truncated_polynomial(2)
    C, d0, d0b, dout = tt.cup_product_matrix(A, 0, 0, P)
    assert d0 == 2 and dout == 2
    H0 = tt.cohomology_classes(A, 0, P)
    # the two reps, as elements of A (C^0 = A), span {1, x}; verify the product
    # table reproduces e_0·e_0=e_0, e_0·e_1=e_1, e_1·e_1=0 once expressed in reps.
    # Convert: the reps are columns of H0.reps (length-2 A-vectors).
    reps = H0.reps % P
    # e_0 = 1 is rep-combination u0; e_1 = x is u1
    u0 = H0.coords(np.array([1, 0], dtype=np.int64))   # 1_A
    u1 = H0.coords(np.array([0, 1], dtype=np.int64))   # x
    def prod(a, b):
        out = np.zeros(dout, dtype=np.int64)
        for i in range(d0):
            for j in range(d0):
                out = (out + int(a[i]) * int(b[j]) * C[:, i, j]) % P
        return out
    assert np.array_equal(prod(u0, u0), u0)            # 1·1 = 1
    assert np.array_equal(prod(u0, u1), u1)            # 1·x = x
    assert np.array_equal(prod(u1, u1) % P, np.zeros(dout, dtype=np.int64))  # x·x = 0


def test_kx2_cup_ring_known_structure():
    """HH^*(k[x]/x^2) over char 0 is k[u,v]/(u^2), deg u = 1, deg v = 2:
    the degree-1 generator squares to zero, the degree-2 generator does not."""
    A = truncated_polynomial(2)
    C11, d1, _, d2 = tt.cup_product_matrix(A, 1, 1, P)
    assert d1 == 1 and d2 == 1
    assert np.array_equal(C11[:, 0, 0] % P, np.zeros(d2, dtype=np.int64))   # u^2 = 0
    C22, da, db, d4 = tt.cup_product_matrix(A, 2, 2, P)
    assert da == 1 and d4 == 1
    assert int(C22[0, 0, 0]) % P != 0                                       # v^2 != 0
