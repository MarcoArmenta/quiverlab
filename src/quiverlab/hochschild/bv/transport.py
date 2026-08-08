"""Plan 54 Task B/E -- the Frobenius perfect pairing (dagger) and the BV transport
(double-dagger) giving Delta: HH^n -> HH^{n-1}.

The pairing (spec s2.1)

    <f, z>_n = lambda( a_0 . f(a_1 (x) ... (x) a_n) ),   f in HH^n, z in HH_n(A, M),

is evaluated by REUSING the Plan-35 cochain-on-chain machinery: with p = n the cap
product ``engine.tt_calculus.cap_cochain(E, n, n, f, z)`` lands in
``C_0 = A`` (the element ``a_0 . f(a_1..a_n)``), and lambda is the Frobenius
covector applied to it. Delta is the adjoint under (dagger):

    <Delta alpha, z>_{n-1} = <alpha, B_nu z>_n    for all z in HH_{n-1},           (double-dagger)

whose transpose bookkeeping over the class bases is (derived, see the Methodology
note in the plan)

    P_{n-1}^T . Delta_n = B_{n-1}^T . P_n^T   =>   Delta_n = (P_{n-1}^T)^{-1} B_{n-1}^T P_n^T,

with ``P_n[i][j] = <coh_i, hom_j>`` (dim HH^n x dim HH_n) and ``B_{n-1}`` the
Connes matrix HH_{n-1} -> HH_n (``connes_b_tables`` / twisted Connes). We solve the
system (no explicit inverse); the PERFECT-PAIRING self-cert requires every ``P_n``
square and invertible in-window, else a loud refusal. Any residual Koszul/pairing
sign is derivation-corrected and arbiter-pinned (``bv/bracket.py``), never assumed.

v1 is GF(p), in-window (the bracket arbiter that certifies correctness is
GF(p)-window-bounded); off GF(p) / past window refuse loudly at the public method.
"""
import numpy as np

from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.fields.primefield import PrimeField

BV_REFERENCES = ["bv_tradler", "bv_lzz", "bv_volkov", "bv_biklz",
                 "cyclic", "bracket", "cup"]


# ---------------------------------------------------------------------------
# the symmetric (trace-form) Frobenius covector, in the unit-adapted basis
# ---------------------------------------------------------------------------
def symmetric_frobenius_form(AU):
    """A VERIFIED-nondegenerate SYMMETRIC trace form ``lambda`` (lambda(ab) =
    lambda(ba)) of ``AU`` as a covector over AU's basis (Domain elements). Reuses
    the presentation-free trace-form space + witness search of
    ``invariants.frobenius`` (the same certificate ``is_symmetric`` uses). Loud if
    no nondegenerate symmetric form exists (i.e. AU is not symmetric)."""
    from quiverlab.fields.linalg import rank as _rank
    from quiverlab.invariants.frobenius import (
        _combine, _gram, _symmetric_witness_candidates, _trace_form_space)
    dom = AU.domain
    m = AU.dim
    S = _trace_form_space(AU)
    if not S:
        raise QuiverlabError("BV symmetric route: no trace form on this algebra")
    for coeffs in _symmetric_witness_candidates(S, dom):
        lam = _combine(S, coeffs, dom, m)
        if _rank(_gram(AU, lam), dom) == m:
            return lam
    raise QuiverlabError(
        "BV symmetric route: no nondegenerate symmetric trace form found "
        "(the algebra is not symmetric)")


# ---------------------------------------------------------------------------
# the pairing matrix P_n over GF(p) (symmetric route: ordinary homology classes)
# ---------------------------------------------------------------------------
def pairing_matrix(E, n, lam, coh_n, hom_n, p):
    """The (dagger) matrix ``P_n[i][j] = <coh_i, hom_j> = lambda(cap_cochain(...))``,
    shape ``dim HH^n x dim HH_n`` over F_p. ``coh_n``/``hom_n`` are engine
    ``_Quotient`` objects; ``lam`` is the Frobenius covector over AU's basis
    (length E.m)."""
    from quiverlab.engine.tt_calculus import cap_cochain
    dl, dr = coh_n.dim, hom_n.dim
    lam_i = [int(lam[t]) % p for t in range(E.m)]
    P = np.zeros((dl, dr), dtype=np.int64)
    for i in range(dl):
        f = coh_n.reps[:, i]
        for j in range(dr):
            z = hom_n.reps[:, j]
            cap = cap_cochain(E, n, n, f, z)          # over cn_basis(E, 0), length m
            acc = 0
            for t in range(E.m):
                ct = int(cap[t])
                if ct:
                    acc += lam_i[t] * ct
            P[i, j] = acc % p
    return P


def _assert_invertible(P, n, p):
    from quiverlab.engine.coxeter import rref_mod_p
    if P.shape[0] != P.shape[1]:
        raise QuiverlabError(
            f"BV perfect-pairing self-cert failed at degree {n}: the pairing "
            f"HH^{n} (x) HH_{n} -> k is {P.shape[0]}x{P.shape[1]}, not square "
            "(dim HH^n != dim HH_n) -- the coefficient bimodule is wrong",
            hint="this algebra's Frobenius duality does not hold in this degree")
    if P.shape[0] == 0:
        return
    _, piv = rref_mod_p(P % p, p)
    if len(piv) != P.shape[0]:
        raise QuiverlabError(
            f"BV perfect-pairing self-cert failed at degree {n}: the pairing "
            f"matrix is singular (rank {len(piv)} < {P.shape[0]}) -- not a "
            "perfect pairing for this coefficient bimodule",
            hint="wrong twist / not a Frobenius duality in this degree")


# ---------------------------------------------------------------------------
# the transport solve Delta_n = (P_{n-1}^T)^{-1} B_{n-1}^T P_n^T over F_p
# ---------------------------------------------------------------------------
def _delta_from_pairing_and_B(Pmats, Bmats, hh_dims, top, p):
    """Solve ``P_{n-1}^T Delta_n = B_{n-1}^T P_n^T`` for n = 1..top over F_p.
    ``Pmats[n]`` is P_n (already invertibility-certified); ``Bmats[n]`` is the
    Connes matrix HH_n -> HH_{n+1} (rows HH_{n+1}, cols HH_n), so ``Bmats[n-1]``
    is B_{n-1}: HH_{n-1} -> HH_n. Returns ``{n: numpy Delta_n}`` (hh^{n-1} x hh^n)."""
    from quiverlab.engine.coxeter import solve_mod_p
    deltas = {}
    for n in range(1, top + 1):
        kprev = hh_dims[n - 1]
        kn = hh_dims[n]
        if kprev == 0 or kn == 0:
            deltas[n] = np.zeros((kprev, kn), dtype=np.int64)
            continue
        Pprev = Pmats[n - 1] % p                       # hh^{n-1} x hh_{n-1} (square)
        Pn = Pmats[n] % p                              # hh^n x hh_n (square)
        Bprev = Bmats[n - 1] % p                       # hh_n x hh_{n-1}
        rhs = (Bprev.T @ Pn.T) % p                     # hh_{n-1} x hh_n
        PprevT = Pprev.T % p                            # hh_{n-1} x hh_{n-1}
        cols = []
        for c in range(rhs.shape[1]):
            x = solve_mod_p(PprevT, rhs[:, c], p)
            if x is None:
                raise QuiverlabError(
                    f"BV transport: the adjoint solve is inconsistent at degree "
                    f"{n} -- the pairing is not perfect (bug/hypothesis failure)")
            cols.append(x % p)
        deltas[n] = (np.array(cols, dtype=np.int64).T % p) if cols else \
            np.zeros((kprev, 0), dtype=np.int64)
    return deltas


def assert_delta_squared_zero(deltas, hh_dims, top, p):
    """Delta_{n-1} . Delta_n = 0 on HH for 2 <= n <= top (the square-zero
    self-cert). Loud on failure."""
    for n in range(2, top + 1):
        Dn = deltas[n] % p                              # hh^{n-1} x hh^n
        Dprev = deltas[n - 1] % p                       # hh^{n-2} x hh^{n-1}
        if Dn.size == 0 or Dprev.size == 0:
            continue
        prod = (Dprev @ Dn) % p                          # hh^{n-2} x hh^n
        if np.any(prod):
            raise QuiverlabError(
                f"BV self-cert Delta^2 = 0 FAILED at degree {n} "
                f"(Delta_{n-1} . Delta_{n} != 0) -- report this algebra")


def _rank_mod_p(M, p):
    from quiverlab.engine.coxeter import rref_mod_p
    if M.size == 0 or M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    _, piv = rref_mod_p(M % p, p)
    return len(piv)


# ---------------------------------------------------------------------------
# the symmetric route (Tradler; P52-free) -- the anchor
# ---------------------------------------------------------------------------
def bv_matrices_symmetric(A, top, max_cells=4_000_000):
    """Delta on HH^*(A) for a SYMMETRIC algebra (spec s2.2): the twist is trivial
    (nu inner), the pairing is the symmetric trace form and B_nu the ORDINARY
    Connes B. GF(p) only in v1. Returns a :class:`BVOperator` (the derived-bracket
    arbiter is attached by ``bv/bracket.py`` when wired)."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine import tt_calculus as TT
    from quiverlab.engine.scan3 import cochain_basis
    from quiverlab.hochschild.products import BVOperator, connes_b_tables
    dom = A.domain
    if not isinstance(dom, PrimeField):
        raise QuiverlabError(
            "BV operator v1 is GF(p) only (the bracket arbiter is "
            "GF(p)-window-bounded)", hint="compute over GF(p)")
    p = dom.p
    AU = A.unit_adapted()
    E = to_engine(AU)
    lam = symmetric_frobenius_form(AU)

    # guard the bar blow-up (cochain/chain pair cells), mirroring Plan 35.
    for n in range(top + 1):
        cells = len(cochain_basis(E, n))
        if cells * cells > max_cells:
            raise DepthLimitError(
                f"BV: the degree-{n} cochain basis pairs {cells * cells} cells "
                f"(> max_cells = {max_cells})", hint="raise max_cells or lower top")

    coh = {n: TT.cohomology_classes(E, n, p) for n in range(top + 1)}
    hom = {n: TT.homology_classes(E, n, p) for n in range(top + 1)}
    hh_dims = [coh[n].dim for n in range(top + 1)]

    Pmats = {}
    for n in range(top + 1):
        Pn = pairing_matrix(E, n, lam, coh[n], hom[n], p)
        _assert_invertible(Pn, n, p)
        Pmats[n] = Pn

    cb = connes_b_tables(A, top, max_cells=max_cells)     # ordinary Connes B
    Bmats = {n: np.array([[int(x) for x in row] for row in cb.matrices[n]],
                         dtype=np.int64).reshape(hh_dims[n + 1], hh_dims[n])
             for n in range(top)}

    deltas = _delta_from_pairing_and_B(Pmats, Bmats, hh_dims, top, p)
    assert_delta_squared_zero(deltas, hh_dims, top, p)

    matrices = {n: [[str(int(deltas[n][i, j])) for j in range(deltas[n].shape[1])]
                    for i in range(deltas[n].shape[0])]
                for n in range(1, top + 1)}
    ranks = {n: _rank_mod_p(deltas[n], p) for n in range(1, top + 1)}
    nak = {"matrix": [[str(1 if i == j else 0) for j in range(A.dim)]
                      for i in range(A.dim)],
           "semisimple": True, "order": 1, "inner": True}
    bv = BVOperator(
        top=top, hh_dims=hh_dims, matrices=matrices, ranks=ranks,
        hypothesis="symmetric (Tradler AIF 2008)", nakayama=nak,
        basis=f"bar/GF({p})", window=top, references=BV_REFERENCES)
    # attach the derived-bracket arbiter (Task D); a mismatch refuses loudly.
    from quiverlab.hochschild.bv.bracket import attach_bracket_arbiter
    attach_bracket_arbiter(A, bv, top, max_cells)
    return bv
