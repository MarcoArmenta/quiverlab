"""The Stefan conjugacy-class Hochschild decomposition of a skew group algebra
``A rtimes G`` (Plan 74 / R8).

Over ``char k does not divide |G|`` (Reynolds averaging),

    ``HH^n(A rtimes G) = ( (+)_{g in G} HH^n(A, {}_gA) )^G
                       = (+)_{[g] conj class} HH^n(A, {}_gA)^{Z(g)}``,

and likewise for homology. Each summand ``HH^n(A, {}_gA)`` is served by P52's
``Bimodule.twisted(A, phi=g_matrix)`` + the chain-level accessors
``twisted_cohomology_classes`` / ``twisted_homology_classes``. THIS module builds
the genuinely-new bridge P52 does not ship:

- :func:`_autom_action_on_classes` -- the chain map induced on the twisted bar
  (co)chain complex by an algebra automorphism ``h in Z(g)`` (the diagonal
  action ``m (x) a_1 (x) ... (x) a_n |-> h(m) (x) h(a_1) (x) ... (x) h(a_n)``,
  a chain map exactly because ``h`` commutes with ``g``), read on the P52 class
  representatives; and
- the Reynolds invariants ``(-)^{Z(g)}`` (idempotent ``e = (1/|Z(g)|) sum_h h``,
  ``dim invariants = rank e``) + the class-block assembly.

The DIRECT ``HH^*(A rtimes G)`` computed by the shipped engine on the constructed
algebra is the trustworthy ground truth; the decomposition is cross-checked
degreewise against it (``verify_direct=True``). Anchors: Stefan (JPAA 103, 1995,
the spectral-sequence origin), Shepler-Witherspoon (J. Algebra 351, 2012, the
additive decomposition + Z(g)-invariance).
"""
from dataclasses import dataclass

from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.families.skew_group import (
    GroupAction, _eye, _mat_apply, _mat_mul)
from quiverlab.fields.linalg import rank, solve
from quiverlab.hochschild.coefficients import (
    Bimodule, twisted_cohomology_classes, twisted_homology_classes)

_DECOMP_CITATIONS = ("stefan_hopf_galois", "shepler_witherspoon_group_actions",
                     "cibils_marcos_smash", "marcos_mv_invariants",
                     "witherspoon_gsm204")


@dataclass
class SkewGroupHH:
    """The assembled Stefan decomposition of ``HH^*(A rtimes G)`` (a data report;
    ``__bool__`` is intentionally undefined, mirroring ``HHTable``)."""
    algebra: object
    action: object
    side: str
    top: int
    dims: list
    summands: list          # per class rep: {"g": label, "hh": [...], "inv": [...]}
    direct_dims: object     # list[int] | None (when verify_direct)
    agrees: object          # bool | None
    status: str             # "complete" | "budget" | "modular" | "unsupported"
    note: str

    def describe(self):
        return f"HH_{self.side} of A|xG ({self.status})"


# --------------------------------------------------------------------------- #
# exact linear-algebra helpers on Domain matrices
# --------------------------------------------------------------------------- #
def _matinv(dom, M):
    n = len(M)
    cols = []
    for j in range(n):
        e = [dom.one() if i == j else dom.zero() for i in range(n)]
        x = solve(M, e, dom)
        if x is None:
            raise QuiverlabError("automorphism matrix is singular (not invertible)",
                                 hint="internal invariant -- the action is not by automorphisms")
        cols.append(x)
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def _transpose(M):
    return [[M[i][j] for i in range(len(M))] for j in range(len(M[0]))] if M else M


def _kron(dom, X, Y):
    nx, mx = len(X), len(X[0]) if X else 0
    ny, my = len(Y), len(Y[0]) if Y else 0
    out = [[dom.zero()] * (mx * my) for _ in range(nx * ny)]
    for i in range(nx):
        for j in range(mx):
            xij = X[i][j]
            if dom.is_zero(xij):
                continue
            for k in range(ny):
                for l in range(my):
                    ylk = Y[k][l]
                    if not dom.is_zero(ylk):
                        out[i * ny + k][j * my + l] = dom.mul(xij, ylk)
    return out


def _kpow(dom, X, n):
    R = X
    for _ in range(n - 1):
        R = _kron(dom, R, X)
    return R


def _to_unit_basis(A, Mh):
    """Express the automorphism matrix ``Mh`` (in A's own basis) in the unit-adapted
    basis ``B = A.unit_adapted()`` -- the basis the bar complex and the P52 classes
    live in. ``Mh^B = P^{-1} Mh P`` with ``P`` = the unit-adapting change of basis."""
    if A.is_unit_adapted:
        return [list(r) for r in Mh]
    dom = A.domain
    P = A._unit_adapting_change()
    Pinv = _matinv(dom, P)
    return _mat_mul(dom, Pinv, _mat_mul(dom, Mh, P))


def _hbar(dom, MhB):
    """The induced map on Abar = A / k.1: drop row 0 and column 0 of ``MhB`` (index
    0 is the unit in the unit-adapted basis)."""
    m = len(MhB)
    return [[MhB[i][j] for j in range(1, m)] for i in range(1, m)]


def _cols(M):
    return [[row[c] for row in M] for c in range(len(M[0]))] if (M and M[0]) else []


# --------------------------------------------------------------------------- #
# the bridge: h-action on the twisted HH classes
# --------------------------------------------------------------------------- #
def _autom_action_on_classes(dom, MhA, hbar, classes, n, side, image_cols):
    """The matrix (in the class-rep basis) of the automorphism ``h`` on
    ``HH_n(A, {}_gA)`` / ``HH^n(...)``, given ``h``'s matrix ``MhA`` on the
    COEFFICIENT slot, its induced map ``hbar`` on Abar, the P52 ``classes``
    (cycle reps), and the boundary/coboundary ``image_cols`` to quotient by.

    Homology transport: ``H_n = MhA (x) hbar^{(x)n}`` (the diagonal action).
    Cohomology transport: ``H^n = MhA (x) ((hbar^{-1})^T)^{(x)n}`` (precomposition
    by ``h^{-1}`` on the argument slots, post by ``h`` on the value slot).

    TWO DIFFERENT BASES, and mixing them is silently wrong on multi-vertex input.
    The bar slots live in the UNIT-ADAPTED basis (``hbar`` is cut from ``h``'s
    unit-adapted matrix), but the coefficient slot does NOT: P52's
    ``bar._coeff_in_unit_basis`` transports a coefficient bimodule by
    ``Bimodule.change_of_basis``, which by contract re-indexes only the A-element
    slot and leaves ``M``'s OWN basis alone -- so ``M`` stays in ``A``'s original
    basis. Hence ``MhA = action.matrix_of(h, A)`` (untransported) is what belongs
    here. The two coincide whenever ``A`` is already unit-adapted (every local
    algebra) or ``h`` fixes the unit-adapting vertex, which is why the error is
    invisible until a VERTEX-PERMUTING automorphism acts on a multi-vertex algebra;
    there it is loud (the transported class leaves the cycle span) rather than
    quietly wrong. Both conventions above are arbitrated non-vacuously against the
    chain-map property in ``tests/hochschild/test_skew_group_transport_basis.py``,
    on an ``hbar`` that is neither symmetric nor self-inverse (so ``hbar``,
    ``hbar^T``, ``hbar^{-1}``, ``(hbar^{-1})^T`` are four distinct matrices and
    exactly one works on each side)."""
    K = len(classes)
    if K == 0:
        return []
    if n == 0:
        Hn = MhA
    else:
        if side == "hom":
            Jblock = hbar
        else:
            Jblock = _transpose(_matinv(dom, hbar)) if (hbar and hbar[0]) else hbar
        Jpow = _kpow(dom, Jblock, n)
        Hn = _kron(dom, MhA, Jpow)
    # reduction basis: [class reps | image columns]
    reduce_cols = [list(c) for c in classes] + [list(c) for c in image_cols]
    dimC = len(reduce_cols[0])
    Mred = [[reduce_cols[j][i] for j in range(len(reduce_cols))] for i in range(dimC)]
    R = []
    for c in classes:
        cprime = _mat_apply(dom, Hn, c)
        x = solve(Mred, cprime, dom)
        if x is None:
            raise QuiverlabError(
                "the automorphism image of a homology class left the cycle span "
                "(the transport is not a chain map -- h does not commute with g?)",
                hint="internal invariant -- please report this action")
        R.append(x[:K])
    return [[R[j][i] for j in range(K)] for i in range(K)]


def _reynolds_invariant_dim(dom, mats, group_order):
    """``dim`` of the ``(-)^H`` invariants = ``rank`` of the Reynolds idempotent
    ``e = (1/|H|) sum_h R_h``. Certifies ``e^2 = e`` (loud otherwise)."""
    if not mats:
        return 0
    K = len(mats[0])
    if K == 0:
        return 0
    inv_order = dom.inv(dom.coerce(group_order))
    e = [[dom.zero()] * K for _ in range(K)]
    for M in mats:
        for i in range(K):
            ei, Mi = e[i], M[i]
            for j in range(K):
                ei[j] = dom.add(ei[j], Mi[j])
    for i in range(K):
        for j in range(K):
            e[i][j] = dom.mul(e[i][j], inv_order)
    e2 = _mat_mul(dom, e, e)
    if not all(dom.eq(e2[i][j], e[i][j]) for i in range(K) for j in range(K)):
        raise QuiverlabError(
            "the Reynolds average e = (1/|H|) sum h is not idempotent (e^2 != e) -- "
            "the group action on HH is not linear/consistent",
            hint="internal invariant -- please report this action")
    return rank(e, dom)


# --------------------------------------------------------------------------- #
# the assembly
# --------------------------------------------------------------------------- #
def stefan_decomposition(A, action, top, *, side="coh", max_cells=4_000_000,
                         verify_direct=False):
    """Assemble ``HH^*(A rtimes G)`` (``side="coh"``) or ``HH_*(...)``
    (``side="hom"``) along the conjugacy classes of ``G``, cross-checked against
    the DIRECT engine when ``verify_direct=True``. Char ``|`` ``|G|`` (modular)
    refuses LOUDLY (the clean decomposition fails; use the direct route)."""
    if side not in ("coh", "hom"):
        raise QuiverlabError(f"side must be 'coh' or 'hom', got {side!r}",
                             hint="'coh' for HH^*, 'hom' for HH_*")
    if not isinstance(action, GroupAction):
        raise QuiverlabError("action must be a GroupAction",
                             hint="build via GroupAction.cyclic / from_generators")
    action.check(A)
    dom = A.domain
    order = action.order
    char = dom.characteristic
    if char != 0 and order % char == 0:
        raise QuiverlabError(
            f"skew-group HH decomposition is not available in the modular case "
            f"(char {char} divides |G| = {order}): Maschke fails, the Reynolds "
            f"average is unavailable, and HH^*(A|xG) != ((+) HH^*(A,{{}}_gA))^G in "
            f"general (Shepler-Witherspoon 1905.09613). The DIRECT HH^*(A|xG) is "
            f"still computable on the constructed algebra.",
            hint="compute A.skew_group(action).hochschild_cohomology(top) directly")

    dims = [0] * (top + 1)
    summands = []
    status = "complete"
    note = ""
    # precompute each element's matrices (shared across summands). TWO BASES, on
    # purpose -- see _autom_action_on_classes: the COEFFICIENT slot stays in A's own
    # basis (Bimodule.change_of_basis never moves M's basis), while the BAR slots are
    # unit-adapted. They differ exactly when a vertex-permuting automorphism acts on a
    # multi-vertex algebra.
    MhA = [action.matrix_of(i, A) for i in range(order)]
    Hbar = [_hbar(dom, _to_unit_basis(A, MhA[i])) for i in range(order)]

    for cls in action.conjugacy_classes:
        gi = cls[0]
        Zg = action.centralizer(gi)
        g_mat = action.matrix_of(gi, A)
        Mg = Bimodule.twisted(A, phi=g_mat)
        try:
            if side == "hom":
                cd = twisted_homology_classes(A, Mg, top + 1, max_cells=max_cells)
            else:
                cd = twisted_cohomology_classes(A, Mg, top, max_cells=max_cells)
        except DepthLimitError as exc:
            status = "budget"
            note = f"summand for class {gi} exceeded max_cells: {exc}"
            summands.append({"g": _label(gi), "hh": None, "inv": None})
            continue
        hh = [len(cd[n]["classes"]) for n in range(top + 1)]
        invs = []
        for n in range(top + 1):
            if side == "hom":
                nb = cd[n + 1]["boundary"] if n + 1 < len(cd) else None
                image_cols = _cols(nb) if nb else []
            else:
                prev = cd[n - 1]["coboundary"] if n >= 1 else None
                image_cols = _cols(prev) if prev else []
            mats = []
            for hk in Zg:
                R = _autom_action_on_classes(dom, MhA[hk], Hbar[hk],
                                             cd[n]["classes"], n, side, image_cols)
                mats.append(R)
            invs.append(_reynolds_invariant_dim(dom, mats, len(Zg)))
        summands.append({"g": _label(gi), "hh": hh, "inv": invs})
        for n in range(top + 1):
            dims[n] += invs[n]

    direct_dims = None
    agrees = None
    if verify_direct:
        direct_dims, dstatus, dnote = _direct_dims(A, action, top, side, max_cells)
        if direct_dims is not None:
            if status != "complete":
                # decomposition capped but direct available -> can't compare full range
                agrees = None
            else:
                agrees = (list(dims) == list(direct_dims))
        else:
            if status == "complete":
                status = dstatus
                note = dnote

    return SkewGroupHH(algebra=A, action=action, side=side, top=top, dims=dims,
                       summands=summands, direct_dims=direct_dims, agrees=agrees,
                       status=status, note=note)


def _direct_dims(A, action, top, side, max_cells):
    from quiverlab.families.skew_group import skew_group_algebra
    S = skew_group_algebra(A, action)
    try:
        if side == "hom":
            return list(S.hochschild_homology(top, max_cells=max_cells).dims), "complete", ""
        return list(S.hochschild_cohomology(top, max_cells=max_cells).dims), "complete", ""
    except DepthLimitError as exc:
        return None, "budget", f"direct HH exceeded max_cells: {exc}"


def _label(i):
    return "e" if i == 0 else f"g{i}"


# --------------------------------------------------------------------------- #
# the three-tier compute-kind block (shared by both runners -> byte-identical)
# --------------------------------------------------------------------------- #
def skew_group_hh_block(A, budget, *, side="coh"):
    """The ``skew_group_hh`` compute-kind block: the Stefan conjugacy-class
    decomposition of ``HH^*(A rtimes G)`` with the DIRECT cross-check row. ``A`` must
    be a ``SkewGroupAlgebra`` (a smash carrying its base + action, stamped by
    :func:`skew_group_algebra`). ``budget`` is the top HH degree (default 3).

    Modular (char ``|`` ``|G|``): ``decomposition=None``, ``status="modular"``,
    ``dims`` = the DIRECT HH only (the constructor + direct HH stay available). A
    relation-violating action raises a clean error entry, never a 500."""
    refs = list(_DECOMP_CITATIONS)
    base = getattr(A, "_skew_base", None)
    action = getattr(A, "_skew_action", None)
    if base is None or action is None:
        return {"kind": "skew_group_hh",
                "error": "skew_group_hh needs a SkewGroupAlgebra input (a base algebra "
                         "with an explicit finite group action)",
                "references": refs}
    top = int(budget) if budget is not None else 3
    dom = base.domain
    order = action.order
    base_dim = base.dim
    char = dom.characteristic
    char_ok = not (char != 0 and order % char == 0)
    block = {"kind": "skew_group_hh", "base_dim": base_dim, "group_order": order,
             "dim": order * base_dim, "char_ok": char_ok, "side": side,
             "references": refs}
    if not char_ok:
        note = (f"modular case (char {char} divides |G| = {order}): the conjugacy-class "
                "decomposition is unavailable (Maschke fails); the DIRECT HH of A|xG is "
                "shown instead.")
        try:
            from quiverlab.families.skew_group import skew_group_algebra
            S = skew_group_algebra(base, action)
            method = (S.hochschild_cohomology if side == "coh"
                      else S.hochschild_homology)
            direct = list(method(top, verbose=False).dims)
            block.update({"decomposition": None, "dims": direct, "direct_dims": direct,
                          "agrees": None, "status": "modular", "note": note})
        except QuiverlabError as exc:
            block.update({"decomposition": None, "dims": None, "direct_dims": None,
                          "agrees": None, "status": "modular", "note": note + f" ({exc})"})
        return block
    try:
        rep = stefan_decomposition(base, action, top, side=side, verify_direct=True)
    except QuiverlabError as exc:
        block["error"] = str(exc)
        return block
    block.update({
        "decomposition": [{"class": s["g"], "hh": s["hh"], "inv": s["inv"]}
                          for s in rep.summands],
        "dims": list(rep.dims),
        "direct_dims": rep.direct_dims,
        "agrees": rep.agrees,
        "status": rep.status,
        "note": rep.note or None,
    })
    return block
