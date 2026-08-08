"""Plan 54 Task D -- the bracket recovered from Delta via the BV relation, and the
in-window ARBITER against the independent Gerstenhaber bracket.

The BV relation expresses the Gerstenhaber bracket as the deviation of Delta from
being a cup-derivation:

    [alpha, beta] = eps(p,q) * ( Delta(alpha U beta) - Delta(alpha) U beta
                                 - (-1)^p alpha U Delta(beta) )                  (BV)

for alpha in HH^p, beta in HH^q. The overall sign ``eps`` is NOT assumed: it is
ARBITRATED ONCE (Plan-21 precedent) by demanding the derived bracket equal the
INDEPENDENT ``gerstenhaber_brackets`` (which quiverlab computes with NO Frobenius
input) in-window over odd primes. That arbitration fixes

    eps(p, q) = -(-1)^{(p-1) q}

(the record's named candidate; empirically ``-1`` when ``(p-1)q`` is even, ``+1``
when odd -- reverse-engineered on k[x,y]/(x^2,y^2), k[x]/(x^N) over GF(5)/GF(32003),
matching every pair including the (2,2) sign flip). The same fixed formula is
applied over EVERY field, GF(2) included (there ``-1 = +1`` and the relation becomes
``Delta(a U b) - Delta(a) U b - a U Delta(b) = [a,b]`` mod 2 -- still a genuine
constraint). A mismatch under this convention is a LOUD refusal (the correctness
gate for every non-symmetric route; the symmetric route it de-risks).
"""
import numpy as np

from quiverlab.errors import QuiverlabError


def _epsilon(p, q):
    """The arbitrated BV sign eps(p,q) = -(-1)^{(p-1)q}, as +1/-1."""
    return -1 if (((p - 1) * q) % 2 == 0) else 1


def _delta_tensor(bv, n):
    """Delta_n as a numpy int matrix (dim HH^{n-1} x dim HH^n), or None out of range."""
    if n < 1 or n > bv.top:
        return None
    rows = bv.matrices[n]
    if not rows:
        return np.zeros((bv.hh_dims[n - 1], bv.hh_dims[n]), dtype=np.int64)
    return np.array([[int(x) for x in row] for row in rows], dtype=np.int64)


def _cup_tensor(cup, p, q, hh):
    """The cup structure constants HH^p (x) HH^q -> HH^{p+q} as (dout, dl, dr), or
    a zero tensor of the right shape when the table is absent/empty."""
    t = cup.tables.get((p, q))
    dl, dr, dout = hh[p], hh[q], hh[p + q]
    C = np.zeros((dout, dl, dr), dtype=np.int64)
    if t is None:
        return C
    for k in range(dout):
        for i in range(dl):
            for j in range(dr):
                C[k, i, j] = int(t.constants[k][i][j])
    return C


def derived_bracket_tables(A, bv, cup, top):
    """Build ``[alpha, beta]`` from Delta and the cup product via (BV) for every
    pair ``p, q >= 1`` with ``p + q <= top`` (so Delta_{p+q}, Delta_p, Delta_q and
    the three cup tables are all in-window). Returns an ``HHProducts(kind="bracket")``
    serialized identically to ``gerstenhaber_brackets``."""
    from quiverlab.hochschild.products import HHProducts, ProductTable
    p_ = A.domain.p
    hh = bv.hh_dims
    tables = {}
    for p in range(1, top + 1):
        for q in range(1, top + 1):
            if p + q > top:
                continue
            out = p + q - 1
            dl, dr, dout = hh[p], hh[q], hh[out]
            Dpq = _delta_tensor(bv, p + q)                 # HH^{p+q} -> HH^{p+q-1}
            Dp = _delta_tensor(bv, p)                       # HH^p -> HH^{p-1}
            Dq = _delta_tensor(bv, q)                       # HH^q -> HH^{q-1}
            Cpq = _cup_tensor(cup, p, q, hh)               # HH^p x HH^q -> HH^{p+q}
            Cpm1q = _cup_tensor(cup, p - 1, q, hh)         # HH^{p-1} x HH^q -> HH^{out}
            Cpqm1 = _cup_tensor(cup, p, q - 1, hh)         # HH^p x HH^{q-1} -> HH^{out}
            T1 = np.einsum('lk,kij->lij', Dpq, Cpq) if (Dpq.size and Cpq.size) \
                else np.zeros((dout, dl, dr), np.int64)
            T2 = np.einsum('ai,laj->lij', Dp, Cpm1q) if (Dp.size and Cpm1q.size) \
                else np.zeros((dout, dl, dr), np.int64)
            T3 = np.einsum('bj,lib->lij', Dq, Cpqm1) if (Dq.size and Cpqm1.size) \
                else np.zeros((dout, dl, dr), np.int64)
            sgn_p = 1 if p % 2 == 0 else -1
            raw = (T1 - T2 - sgn_p * T3) % p_
            C = (_epsilon(p, q) * raw) % p_
            tables[(p, q)] = ProductTable(
                kind="bracket", degrees=(p, q), out_degree=out,
                dims=(dl, dr, dout),
                constants=tuple(tuple(tuple(str(int(C[k, i, j])) for j in range(dr))
                                      for i in range(dl)) for k in range(dout)))
    window = max((p + q - 1 for (p, q) in tables), default=0)
    return HHProducts(
        kind="bracket", top=top, tables=tables,
        engine="BV-derived (from Delta via the BV relation)",
        basis=bv.basis, window=window,
        references=["bracket", "gerstenhaber", "bv_lzz"])


def bracket_check(A, derived, independent, top):
    """Compare the derived bracket to the independent Gerstenhaber bracket
    table-for-table over the pairs both cover; returns
    ``{"engine", "window", "agrees"}``. Constants are compared as ints mod p (the
    class basis is shared -- the same tt cohomology classes both were built on)."""
    p_ = A.domain.p
    agrees = True
    checked = 0
    for key, dt in derived.tables.items():
        it = independent.tables.get(key)
        if it is None:
            continue
        if dt.dims != it.dims:
            agrees = False
            break
        checked += 1
        dl, dr, dout = dt.dims
        for k in range(dout):
            for i in range(dl):
                for j in range(dr):
                    if int(dt.constants[k][i][j]) % p_ != int(it.constants[k][i][j]) % p_:
                        agrees = False
                        break
                if not agrees:
                    break
            if not agrees:
                break
        if not agrees:
            break
    window = derived.window if agrees else 0
    return {"engine": independent.engine, "window": window,
            "agrees": bool(agrees and checked > 0), "pairs_checked": checked}


def attach_bracket_arbiter(A, bv, top, max_cells):
    """Attach ``bv.derived_bracket`` (built from Delta via (BV)) and
    ``bv.bracket_check`` (the in-window derived == independent equality). A
    mismatch is a LOUD refusal -- never a silent wrong Delta. For ``top < 2`` there
    is no bracket pair in-window; the arbiter is vacuous and records so (the Delta
    self-certs Delta^2=0 + perfect pairing still stand)."""
    if top < 2:
        bv.bracket_check = {"engine": "gerstenhaber", "window": 0,
                            "agrees": False, "pairs_checked": 0}
        return bv
    cup = A.cup_products(top, max_cells=max_cells)
    independent = A.gerstenhaber_brackets(top, max_cells=max_cells)
    derived = derived_bracket_tables(A, bv, cup, top)
    chk = bracket_check(A, derived, independent, top)
    if not chk["agrees"]:
        raise QuiverlabError(
            "BV transport does not reproduce the independent Gerstenhaber bracket "
            "in-window -- the hypothesis/convention does not certify for this "
            "instance (no silent wrong Delta)",
            hint="the algebra may be outside the certified BV scope")
    bv.bracket_check = chk
    bv.derived_bracket = derived
    return bv


def seven_term_residual(A, bv, cup, p, q, r):
    """The seven-term (BV / Koszul) relation residual for degrees (p, q, r) with
    p, q, r >= 1 and p+q+r <= bv.top, over F_p; a zero residual certifies Delta is a
    differential operator of order <= 2 (the BV-algebra characterization, checked
    directly, not merely invoked). Returns the max abs residual mod p (0 == holds).

        Delta(aUbUc) = Delta(aUb)Uc + (-1)^p aUDelta(bUc) + (-1)^{(p-1)q} bUDelta(aUc)
                       - Delta(a)UbUc - (-1)^p aUDelta(b)Uc - (-1)^{p+q} aUbUDelta(c).
    """
    p_ = A.domain.p
    hh = bv.hh_dims

    def C(x, y):
        return _cup_tensor(cup, x, y, hh)

    def D(n):
        return _delta_tensor(bv, n)

    dp, dq, dr = hh[p], hh[q], hh[r]
    # cup triple products, associative: (aUb)Uc == aU(bUc); use left assoc.
    # Uabc[l,i,j,k] = coeff of out_l in a_i U b_j U c_k, out in HH^{p+q+r}
    Cab = C(p, q)                       # (HH^{p+q}, p, q)
    Cab_c = C(p + q, r)                 # (HH^{p+q+r}, p+q, r)
    Uabc = np.einsum('mij,lmk->lijk', Cab, Cab_c)      # HH^{p+q+r}
    Cbc = C(q, r)
    Ca_bc = C(p, q + r)
    Uabc2 = np.einsum('mjk,lim->lijk', Cbc, Ca_bc)
    assert not np.any((Uabc - Uabc2) % p_), "cup not associative (bug)"

    sp = 1 if p % 2 == 0 else -1
    spq = 1 if (p + q) % 2 == 0 else -1
    spm1q = 1 if ((p - 1) * q) % 2 == 0 else -1

    # LHS: Delta(aUbUc): Delta_{p+q+r} applied to Uabc
    Dtot = D(p + q + r)
    LHS = np.einsum('ol,lijk->oijk', Dtot, Uabc) % p_   # HH^{p+q+r-1}
    dout = hh[p + q + r - 1]

    # RHS term A: Delta(aUb) U c
    DCab = np.einsum('ml,lij->mij', D(p + q), Cab) % p_   # Delta(aUb): HH^{p+q-1}
    Cpqm1_r = C(p + q - 1, r)
    tA = np.einsum('mij,omk->oijk', DCab, Cpqm1_r) % p_
    # RHS term B: (-1)^p a U Delta(bUc)
    DCbc = np.einsum('ml,ljk->mjk', D(q + r), Cbc) % p_   # HH^{q+r-1}
    Cp_qrm1 = C(p, q + r - 1)
    tB = (sp * np.einsum('mjk,oim->oijk', DCbc, Cp_qrm1)) % p_
    # RHS term C: (-1)^{(p-1)q} b U Delta(aUc)
    Cac = C(p, r)
    DCac = np.einsum('ml,lik->mik', D(p + r), Cac) % p_   # HH^{p+r-1}
    Cq_prm1 = C(q, p + r - 1)
    tC = (spm1q * np.einsum('mik,ojm->oijk', DCac, Cq_prm1)) % p_
    # RHS term D: - Delta(a) U b U c
    Da = D(p)
    aD_b = np.einsum('ai,laj->lij', Da, C(p - 1, q)) % p_  # Delta(a)Ub: HH^{p+q-1}
    tD = np.einsum('lij,olk->oijk', aD_b, C(p + q - 1, r)) % p_
    # RHS term E: - (-1)^p a U Delta(b) U c
    Db = D(q)
    aDb = np.einsum('bj,lib->lij', Db, C(p, q - 1)) % p_   # aUDelta(b): HH^{p+q-1}
    tE = (sp * np.einsum('lij,olk->oijk', aDb, C(p + q - 1, r))) % p_
    # RHS term F: - (-1)^{p+q} a U b U Delta(c)
    Dc = D(r)
    bDc = np.einsum('ck,ljc->ljk', Dc, C(q, r - 1)) % p_   # bUDelta(c): HH^{q+r-1}
    tF = (spq * np.einsum('mjk,oim->oijk', bDc, C(p, q + r - 1))) % p_

    RHS = (tA + tB + tC - tD - tE - tF) % p_
    return int(np.max(np.abs(((LHS - RHS) % p_).astype(np.int64)))) if LHS.size else 0
