"""Plan 54 Task E -- the twisted Connes operator ``B_sigma`` on the sigma-twisted bar
complex, and its descent to twisted Hochschild homology.

The sigma-twisted chain complex ``C_n = M (x) Abar^{(x)n}`` with ``M = {}_1A_sigma``
(P52's convention, spec s3: right action twisted, twist on the FIRST face) is a
PARACYCLIC (not cyclic) module: the twisted cyclic operator

    t_sigma(a_0 (x) ... (x) a_n) = (-1)^n ( a_n (x) sigma(a_0) (x) a_1 (x) ... (x) a_{n-1} )

satisfies ``d_0 t = d_n`` against P52's faces, but ``t^{n+1}`` is a diagonal
``sigma``-action rather than the identity. Hence the naive Connes ``B = s o N``
(``N = sum_{i=0}^{n} t^i``, ``s`` inserting the unit in the M-slot) does NOT satisfy
``bB + Bb = 0`` at CHAIN level -- the defect is ``1 - T`` with ``T`` the paracyclic
``sigma``-action.

BUT for a SEMISIMPLE ``nu`` (the LZZ hypothesis) that action is trivial ON HOMOLOGY,
so ``B_sigma`` DESCENDS: it carries twisted CYCLES to twisted cycles (``b B z = 0``
whenever ``b z = 0``), well-defined modulo boundaries. This module builds the
chain-level ``B_sigma`` and PROJECTS it to the class basis, asserting the descent
per instance (``B(rep)`` is a cycle -- a loud refusal otherwise). The ultimate
correctness gate remains the bracket arbiter (``bv/bracket.py``).

``twist`` is the RIGHT twist of the coefficient itself (``twist = nu`` for
``{}_1A_nu``): the wrapping element carries exactly this matrix. Whether it is
``nu`` or its inverse is NOT assumed -- it is ARBITRATED per instance by the descent
self-cert (``b B z = 0`` on cycles) plus the bracket arbiter; the wrong power fails
the descent loudly (it never returns a silent wrong Delta).
"""
from itertools import product

import numpy as np

from quiverlab.errors import QuiverlabError


def _apply_twist(twist, vec, m, p):
    """The twist matrix applied to a slot-vector ``vec`` (a dict idx->coeff);
    ``twist[:, idx]`` is the image of basis element ``idx`` (columns = images)."""
    out = {}
    for idx, c in vec.items():
        col = twist[:, idx]
        for r in range(m):
            w = int(col[r]) % p
            if w:
                out[r] = (out.get(r, 0) + c * w) % p
    return out


def twisted_cyclic_tau(AU, twist, n, p):
    """The twisted cyclic operator ``t_sigma`` on ``C_n`` (bar basis, rows = cols =
    ``C_n``): ``t(a_0,...,a_n) = (-1)^n (a_n, twist(a_0), a_1,...,a_{n-1})``.
    ``t_0 = id``. ``twist`` is the coefficient's Nakayama matrix (columns = images)."""
    from quiverlab.hochschild import bar
    m = AU.dim
    cols = bar._cochain_basis(m, n, m)
    if n == 0:
        return np.eye(len(cols), dtype=np.int64)
    ri = {b: i for i, b in enumerate(cols)}
    T = np.zeros((len(cols), len(cols)), dtype=np.int64)
    sign = 1 if n % 2 == 0 else p - 1
    for ci, (s, J) in enumerate(cols):
        slots = [{si: 1} for si in ([s] + list(J))]
        v0, vlast = slots[0], slots[-1]
        new = [vlast, _apply_twist(twist, v0, m, p)] + slots[1:-1]
        keyed = [list(d.items()) for d in new]
        for combo in product(*keyed):
            cc = sign
            sidx = None
            Jout = []
            ok = True
            for pos, (idx, c) in enumerate(combo):
                cc = (cc * c) % p
                if pos == 0:
                    sidx = idx
                elif idx == 0:                        # unit in a bar slot -> degenerate
                    ok = False
                    break
                else:
                    Jout.append(idx)
            if not ok or cc == 0:
                continue
            r = ri[(sidx, tuple(Jout))]
            T[r, ci] = (T[r, ci] + cc) % p
    return T


def _s_extra_degeneracy(AU, n, p):
    """``s : C_n -> C_{n+1}``, ``s(a_0,...,a_n) = (1, a_0,...,a_n)`` (insert the unit
    in the M-slot; all n+1 old factors become bar factors, degenerate ones die)."""
    from quiverlab.hochschild import bar
    m = AU.dim
    cols = bar._cochain_basis(m, n, m)
    rows = bar._cochain_basis(m, n + 1, m)
    ri = {b: i for i, b in enumerate(rows)}
    S = np.zeros((len(rows), len(cols)), dtype=np.int64)
    for ci, (s, J) in enumerate(cols):
        newbar = (s,) + tuple(J)
        if 0 in newbar:                               # old M-slot was the unit
            continue
        S[ri[(0, newbar)], ci] = 1
    return S


def twisted_connes_matrix(AU, twist, n, p):
    """The chain-level twisted Connes ``B_sigma : C_n -> C_{n+1}`` (bar basis) =
    ``s o N``, ``N = sum_{i=0}^{n} t_sigma^i``. Descends to homology for semisimple
    nu (asserted downstream)."""
    T = twisted_cyclic_tau(AU, twist, n, p)
    N = np.eye(T.shape[0], dtype=np.int64)
    cur = np.eye(T.shape[0], dtype=np.int64)
    for _ in range(1, n + 1):
        cur = (cur @ T) % p
        N = (N + cur) % p
    S = _s_extra_degeneracy(AU, n, p)
    return (S @ N) % p


def twisted_connes_class_matrix(AU, twist, tw, n, p):
    """``B_{n} : HH_n(A, {}_1A_sigma) -> HH_{n+1}`` in the CLASS basis (rows =
    HH_{n+1} classes, cols = HH_n classes). Builds the chain ``B_sigma``, checks
    each ``B(rep_j)`` is a twisted CYCLE (the descent self-cert -- loud refusal on
    failure), and expresses it in the target class reps modulo boundaries.

    ``tw`` is the :func:`bv.twist.twisted_homology_quotient` output (``reps``,
    ``boundary`` per degree); ``twist`` is the coefficient's Nakayama matrix."""
    from quiverlab.engine.coxeter import solve_mod_p
    rows = tw[n + 1]["dim"]
    cols = tw[n]["dim"]
    out = np.zeros((rows, cols), dtype=np.int64)
    if rows == 0 or cols == 0:
        return out
    Bchain = twisted_connes_matrix(AU, twist, n, p)      # C_n -> C_{n+1}
    bnp1 = tw[n + 1]["boundary"]                             # b_{n+1}: C_{n+1} -> C_n
    bnp2 = tw[n + 2]["boundary"] if (n + 2) < len(tw) else None
    reps_next = tw[n + 1]["reps"]
    # column space of the class reps + the boundary image im b_{n+2} inside C_{n+1}
    basis_cols = [reps_next[j] for j in range(rows)]
    if bnp2 is not None and bnp2.size:
        basis_cols += [bnp2[:, c] for c in range(bnp2.shape[1])]
    Bmat = np.array(basis_cols, dtype=np.int64).T % p
    for j in range(cols):
        v = (Bchain @ tw[n]["reps"][j]) % p                  # B(rep_j) in C_{n+1}
        if bnp1 is not None and bnp1.size and np.any((bnp1 @ v) % p):
            raise QuiverlabError(
                f"BV twisted route: Connes B_sigma does not descend to homology at "
                f"degree {n} (B of a cycle is not a cycle) -- the Nakayama "
                f"automorphism does not act trivially on twisted homology; the "
                "semisimple-nu hypothesis is not certified for this instance",
                hint="only symmetric or semisimple-nu Frobenius algebras are served")
        x = solve_mod_p(Bmat, v % p, p)
        if x is None:
            raise QuiverlabError(
                f"BV twisted route: B_sigma(class) at degree {n} is not in the "
                "span of the target classes and boundaries (internal invariant) "
                "-- please report this algebra")
        out[:, j] = x[:rows] % p
    return out
