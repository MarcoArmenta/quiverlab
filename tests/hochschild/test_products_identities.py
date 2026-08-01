"""Plan 35 identity batteries. Table-level: the identities of a Gerstenhaber
algebra (graded-commutative associative cup, antisymmetric bracket, cup-Leibniz)
+ the cap module structure, on small algebras over several primes, plus
bar<->CS table equality in-window (the cross-engine gate). The induced B^2=0
identity is covered separately in tests/hochschild/test_connes_b.py."""
import itertools
import pytest

import quiverlab as ql

PRIMES = (32003, 2, 3, 5)


def _mat(table):
    """constants as int tensor [k][i][j]."""
    return [[[int(c) for c in row] for row in mat] for mat in table.constants]


def _compose(hp, p, q, r, prime):
    """(f_i ∪ g_j) ∪ h_l and f_i ∪ (g_j ∪ h_l) as coordinate tensors — equality
    is associativity at the table level."""
    Tpq, Tqr = hp.tables[(p, q)], hp.tables[(q, r)]
    Tpq_r, Tp_qr = hp.tables[(p + q, r)], hp.tables[(p, q + r)]
    dl, dm, dr_ = Tpq.dims[0], Tpq.dims[1], Tqr.dims[1]
    dout = Tpq_r.dims[2]
    left = [[[[0] * dout for _ in range(dr_)] for _ in range(dm)] for _ in range(dl)]
    right = [[[[0] * dout for _ in range(dr_)] for _ in range(dm)] for _ in range(dl)]
    A_, B_, C_, D_ = _mat(Tpq), _mat(Tpq_r), _mat(Tqr), _mat(Tp_qr)
    for i, j, l in itertools.product(range(dl), range(dm), range(dr_)):
        for k in range(dout):
            left[i][j][l][k] = sum(A_[m][i][j] * B_[k][m][l]
                                   for m in range(Tpq.dims[2])) % prime
            right[i][j][l][k] = sum(C_[m][j][l] * D_[k][i][m]
                                    for m in range(Tqr.dims[2])) % prime
    return left, right


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("prime", PRIMES)
def test_cup_graded_commutative_and_associative(prime):
    A = ql.truncated_polynomial(3, field=ql.GF(prime))
    hp = A.cup_products(3)
    for (p, q), t in hp.tables.items():
        if (q, p) not in hp.tables:
            continue
        s = hp.tables[(q, p)]
        sign = 1 if (p * q) % 2 == 0 else prime - 1
        M, N = _mat(t), _mat(s)
        for k in range(t.dims[2]):
            for i in range(t.dims[0]):
                for j in range(t.dims[1]):
                    assert M[k][i][j] % prime == (sign * N[k][j][i]) % prime
    for p, q, r in [(0, 0, 1), (0, 1, 1), (1, 1, 1)]:
        left, right = _compose(hp, p, q, r, prime)
        assert left == right, f"associativity fails at {(p, q, r)}"


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("prime", PRIMES)
def test_bracket_antisymmetry(prime):
    A = ql.truncated_polynomial(3, field=ql.GF(prime))
    hb = A.gerstenhaber_brackets(3)
    for (p, q), t in hb.tables.items():
        if (q, p) not in hb.tables:
            continue
        s, M, N = hb.tables[(q, p)], _mat(t), _mat(hb.tables[(q, p)])
        sign = 1 if ((p - 1) * (q - 1)) % 2 == 0 else prime - 1
        for k in range(t.dims[2]):
            for i in range(t.dims[0]):
                for j in range(t.dims[1]):
                    assert M[k][i][j] % prime == (-sign * N[k][j][i]) % prime, \
                        f"[f,g] != -(-1)^((p-1)(q-1))[g,f] at {(p, q)}"


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("prime", (7, 3))
def test_bar_vs_cs_cup_tables_in_window(prime):
    """The cross-engine gate: same DIMS and RANK-equivalent tables. Bases
    differ, so we compare basis-independent data: dims and, for each (p,q),
    the RANK of the flattened constants matrix (dl*dr x dout) mod p."""
    import numpy as np
    from quiverlab.engine.coxeter import rref_mod_p
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(prime))
    bar = A.cup_products(2, engine="bar")
    cs = A.cup_products(2, engine="cs")
    assert sorted(bar.tables) == sorted(cs.tables)
    for key in bar.tables:
        tb, tc = bar.tables[key], cs.tables[key]
        assert tb.dims == tc.dims
        def flat_rank(t):
            dl, dr, dout = t.dims
            if 0 in t.dims:
                return 0
            M = np.array([[int(t.constants[k][i][j]) for k in range(dout)]
                          for i in range(dl) for j in range(dr)], dtype=np.int64)
            _, piv = rref_mod_p(M % prime, prime)
            return len(piv)
        assert flat_rank(tb) == flat_rank(tc), f"table rank differs at {key}"


# ---------------------------------------------------------------------------
# Cap module identity: (z ∩ f) ∩ g = z ∩ (f ∪ g).
#
# Cap tables have degrees (p, n) meaning HH^p (x) HH_n -> HH_{n-p}, so
# constants[k][i][j] pairs a cohomology class (left index i, in HH^p) with a
# homology class (right index j, in HH_n) to give a homology class (out index
# k, in HH_{n-p}). Take f ∈ HH^p (index a), g ∈ HH^q (index c), z ∈ HH_n
# (index b); the shared output lives in HH_{n-p-q} (index k).
#
#   z ∩ f          = cap(f, z),        table Cap(p,   n)   [·][a][b]  -> HH_{n-p}
#   (z ∩ f) ∩ g    = cap(g, z∩f),      table Cap(q, n-p)   [·][c][·]  -> HH_{n-p-q}
#   f ∪ g          = cup(f, g),        table Cup(p,   q)   [·][a][c]  -> HH^{p+q}
#   z ∩ (f ∪ g)    = cap(f∪g, z),      table Cap(p+q, n)   [·][·][b]  -> HH_{n-p-q}
#
# Chaining (sum over the shared intermediate index m):
#   LHS[a][c][b][k] = Σ_m Cap(p,n)[m][a][b]   · Cap(q,n-p)[k][c][m]   (m ∈ HH_{n-p})
#   RHS[a][c][b][k] = Σ_m Cup(p,q)[m][a][c]   · Cap(p+q,n)[k][m][b]   (m ∈ HH^{p+q})
# The module axiom is exactly LHS == RHS (mod prime), no sign.
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
def test_cap_module_identity():
    prime = 7
    A = ql.truncated_polynomial(2, field=ql.GF(prime))
    cup = A.cup_products(3)
    cap = A.cap_products(3)
    for p, q, n in [(0, 0, 1), (0, 1, 2), (1, 1, 2), (1, 2, 3)]:
        Cpn = cap.tables[(p, n)]           # HH^p (x) HH_n     -> HH_{n-p}
        Cqm = cap.tables[(q, n - p)]       # HH^q (x) HH_{n-p} -> HH_{n-p-q}
        Upq = cup.tables[(p, q)]           # HH^p (x) HH^q     -> HH^{p+q}
        Cpqn = cap.tables[(p + q, n)]      # HH^{p+q} (x) HH_n -> HH_{n-p-q}
        M_Cpn, M_Cqm, M_Upq, M_Cpqn = _mat(Cpn), _mat(Cqm), _mat(Upq), _mat(Cpqn)
        da, dc, db = Cpn.dims[0], Upq.dims[1], Cpn.dims[1]
        dk = Cqm.dims[2]
        for a, c, b, k in itertools.product(range(da), range(dc), range(db), range(dk)):
            lhs = sum(M_Cpn[m][a][b] * M_Cqm[k][c][m]
                      for m in range(Cpn.dims[2])) % prime
            rhs = sum(M_Upq[m][a][c] * M_Cpqn[k][m][b]
                      for m in range(Upq.dims[2])) % prime
            assert lhs == rhs, \
                f"cap module identity fails at (p,q,n)={(p, q, n)} idx {(a, c, b, k)}"


# ---------------------------------------------------------------------------
# Bracket cup-Leibniz: [f, g ∪ h] = [f, g] ∪ h + (-1)^{(p-1)q} g ∪ [f, h].
#
# Bracket tables have degrees (p, q) meaning HH^p (x) HH^q -> HH^{p+q-1}, with
# left index = HH^p, right index = HH^q. Take f ∈ HH^p (index a), g ∈ HH^q
# (index c), h ∈ HH^r (index e); the shared output lives in HH^{p+q+r-1}
# (index k). Every piece chains through one shared intermediate index m:
#
#   g ∪ h        = cup(g, h),         Cup(q, r)       [m][c][e]  -> HH^{q+r}
#   [f, g ∪ h]   = bracket(f, g∪h),   Br(p, q+r)      [k][a][m]  -> HH^{p+q+r-1}
#     LHS[a][c][e][k] = Σ_m Cup(q,r)[m][c][e] · Br(p,q+r)[k][a][m]      (m ∈ HH^{q+r})
#
#   [f, g]       = bracket(f, g),     Br(p, q)        [m][a][c]  -> HH^{p+q-1}
#   [f, g] ∪ h   = cup([f,g], h),     Cup(p+q-1, r)   [k][m][e]  -> HH^{p+q+r-1}
#     T1[a][c][e][k]  = Σ_m Br(p,q)[m][a][c] · Cup(p+q-1,r)[k][m][e]    (m ∈ HH^{p+q-1})
#
#   [f, h]       = bracket(f, h),     Br(p, r)        [m][a][e]  -> HH^{p+r-1}
#   g ∪ [f, h]   = cup(g, [f,h]),     Cup(q, p+r-1)   [k][c][m]  -> HH^{p+q+r-1}
#     T2[a][c][e][k]  = Σ_m Br(p,r)[m][a][e] · Cup(q,p+r-1)[k][c][m]    (m ∈ HH^{p+r-1})
#
# The identity is LHS == T1 + (-1)^{(p-1)q} · T2 (mod prime).
#
# Fixture note (non-vacuity): on k[x]/(x^a) the cup(1,1) table is zero in odd
# characteristic, so on GF(7) EVERY term of the (1,1,1) identity would carry a
# vanishing cup(1,1) factor and the check would collapse to 0 == 0 — passing a
# wrong composite equally. GF(2) is the content-bearing case there
# (cup(1,1)/br(1,1) nonzero); its sign exponent (p-1)q = 0 (sign +1), but ±1
# also coincide mod 2, so GF(2) cannot distinguish the sign. The (2,1,1) case
# on k[x]/(x^3) over GF(3) carries a genuine sign: (p-1)q = 1 (sign -1 ≠ +1
# mod 3) and, although its LHS is zero there (cup(1,1)=0 and br(2,2)=0 in odd
# char), T1 and T2 are both nonzero, so the identity reduces to the non-trivial
# 0 == T1 - T2 — a wrong sign fails it. The two rows together cover both a
# content-bearing full identity (GF(2)) and a content-bearing sign (GF(3)).
# ---------------------------------------------------------------------------
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("a, prime, p, q, r", [(2, 2, 1, 1, 1), (3, 3, 2, 1, 1)])
def test_bracket_cup_leibniz(a, prime, p, q, r):
    A = ql.truncated_polynomial(a, field=ql.GF(prime))
    top = p + q + r - 1
    cup = A.cup_products(top)
    br = A.gerstenhaber_brackets(top)
    Uqr = cup.tables[(q, r)]              # g ∪ h
    Bpqr = br.tables[(p, q + r)]          # [f, g∪h]
    Bpq = br.tables[(p, q)]               # [f, g]
    Upqr = cup.tables[(p + q - 1, r)]     # [f,g] ∪ h
    Bpr = br.tables[(p, r)]               # [f, h]
    Uqpr = cup.tables[(q, p + r - 1)]     # g ∪ [f,h]
    M_Uqr, M_Bpqr = _mat(Uqr), _mat(Bpqr)
    M_Bpq, M_Upqr = _mat(Bpq), _mat(Upqr)
    M_Bpr, M_Uqpr = _mat(Bpr), _mat(Uqpr)
    sign = 1 if ((p - 1) * q) % 2 == 0 else prime - 1
    da, dc, de = Bpq.dims[0], Bpq.dims[1], Bpr.dims[1]
    dk = Bpqr.dims[2]
    for a_, c, e, k in itertools.product(range(da), range(dc), range(de), range(dk)):
        lhs = sum(M_Uqr[m][c][e] * M_Bpqr[k][a_][m]
                  for m in range(Uqr.dims[2])) % prime
        t1 = sum(M_Bpq[m][a_][c] * M_Upqr[k][m][e]
                 for m in range(Bpq.dims[2])) % prime
        t2 = sum(M_Bpr[m][a_][e] * M_Uqpr[k][c][m]
                 for m in range(Bpr.dims[2])) % prime
        rhs = (t1 + sign * t2) % prime
        assert lhs == rhs, \
            f"bracket cup-Leibniz fails at (p,q,r)={(p, q, r)} idx {(a_, c, e, k)}"
