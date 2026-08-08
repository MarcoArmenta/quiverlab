"""Split-extension / trivial-extension Hochschild long exact sequence (Plan 72, R5).

For the split square-zero extension ``L = B ⋉ M`` (``M`` a two-sided ideal with
``M·M = 0``, ``L = B ⊕ M`` as ``k``-spaces), the short exact sequence of
``L``-bimodules ``0 → M → L → B → 0`` induces the **long exact sequence in
``HH^•(L, −)``** (Cibils–Marcos–Redondo–Solotar, ``math/0102194``)::

    ··· → HH^n(L,M) →^ι HH^n(L,L) →^π HH^n(L,B) →^{δ^n} HH^{n+1}(L,M) → ···

whose middle term ``HH^n(L,L) = HH^n(L)`` is assembled from the two flanks and the
connecting map ``δ`` (CMRS **Thm 4.1**: ``δ`` is the cup with the extension class),
and cross-checked against the direct ``HH^•(L)``. The flagship is the **trivial
extension** ``L = T(B) = B ⋉ D(B)`` (``M = D(B)``), the ``TrivialExtension(B)``
quiverlab already builds (Plan 31). Also: ``HH^1(B ⋉ M) ≠ 0`` for ``M ≠ 0`` (CMRS
**Thm 5.5**), the grading-derivation witness ``E(b+m) = m``; for directed ``B`` the
sharper ``HH^1(T(B)) = k ⊕ HH^1(B)``.

**How the LES is COMPUTED (design decision, Plan 72 fix round).** The card's snake
is realized on the **Chouhy–Solotar** ``L^e``-projective resolution rather than the
bar resolution. Reason: at the plan's own pinned depths the bar coboundary is
intractable (e.g. ``T(kD₄)``'s ``δ^4`` lands in ``HH^5(L,M)`` and the bar chain
group there has ``> 2·10^6`` rows), while CS is polynomial and reproduces every
live-verified pin exactly (``T(kD₄)`` flanks ``[4,1,0,0,1,4]``/``[1,0,0,1,4]``,
``δ = [0,0,0,1,2]``). The SES ``0 → M → L → B → 0`` of ``L``-bimodules makes the
coefficient space ``L = M ⊕ B`` split with ``M`` an **ideal (sub-bimodule)**, so the
regular CS Hom-complex is **block upper-triangular** in that split::

    δ_reg = [[δ_M,  X ],
             [ 0 , δ_B]]

(the ``B``-rows / ``M``-cols block is zero because ``M`` is a sub-bimodule). ``δ_M``
is ``C^•(L,M)``, ``δ_B`` is ``C^•(L,B)``, and the coupling ``X`` IS the snake
connecting map: for a ``B``-cocycle ``φ`` the lift ``(0, φ)`` has coboundary
``(Xφ, 0)`` whose ``M``-class is ``δ[φ]``. This is a genuine, computed ``δ`` (not a
dimension-implied rank), and ``assembled == direct`` is a real cross-engine oracle.
The ``HH^1`` grading witness is materialized on the (small, degrees 0–1–2) bar
complex.
"""
from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import nullspace, rank, solve


@dataclass
class SplitExtReport:
    base: object                       # B
    extension: object                  # L = B ⋉ M
    coeff_name: str                    # "D(B)" | the Bimodule.describe()
    side: str                          # "cohomology" | "homology"
    n_range: int                       # top
    flank_M: list                      # dim HH^n(L, M)  (n = 0..top+1 -- BOUNDARY RULE)
    flank_B: list                      # dim HH^n(L, B)  (n = 0..top)   [homology: 0..top+1]
    delta_ranks: list                  # rank δ^n        (computed snake ranks, n = 0..top)
    assembled: list                    # dim HH^n(L) from the LES  (n = 0..top)
    direct: list = None                # dim HH^n(L) directly (the cross-check; None if skipped)
    leading_B: list = field(default_factory=list)   # HH^n(B)  (informational; the p=0 term)
    leading_M: list = field(default_factory=list)   # H^n(B,M) (informational)
    hh1_nonzero: bool = None           # HH^1(L) != 0 (the grading-derivation witness)
    hh1_directed_formula: object = None  # HH^1(L) == 1 + HH^1(B) (directed B; None if not directed)
    exact: bool = None                 # the LES exactness self-cert holds
    agrees: object = None              # assembled == direct
    status: str = "complete"           # "complete" | "unsupported" | "budget" | "error"
    note: str = ""


# --------------------------------------------------------------------------- #
# constructor
# --------------------------------------------------------------------------- #
def split_extension(B, M=None, *, name=None):
    """The split square-zero extension ``L = B ⋉ M``. With ``M = None`` (default)
    returns ``TrivialExtension(B)`` (``M = D(B)``, presented, Plan 31), certified
    ``dim L = 2·dim B``. A general bimodule ``M`` is refused loudly (only the
    trivial-extension flagship, whose ``L`` is presentable via ``TrivialExtension``,
    is in scope — DD-A1)."""
    if M is not None:
        raise QuiverlabError(
            "split_extension currently supports only M = None (the trivial "
            "extension L = T(B) = B ⋉ D(B), whose L is presentable via "
            "TrivialExtension); a general bimodule M has no presented-L build here",
            hint="pass M=None for the trivial-extension LES")
    if B.quiver is None:
        raise QuiverlabError(
            "the split-extension LES needs a quiver presentation of B to inflate "
            "coefficients", hint="present B via Quiver(...).algebra(...)")
    from quiverlab.families.trivial_extension import TrivialExtension
    L = TrivialExtension(B)
    if L.dim != 2 * B.dim:                       # self-cert (mirrors Plan 31 D2)
        raise QuiverlabError(
            f"split_extension: dim L = {L.dim} must be 2·dim B = {2 * B.dim}",
            hint="the trivial-extension build failed its dimension certificate")
    if name is not None and L.basis_labels is not None:
        pass                                      # name is advisory; L already labelled
    return L


# --------------------------------------------------------------------------- #
# the M-block (ideal) of the coefficient space L
# --------------------------------------------------------------------------- #
def _mblock_indices(L, B):
    """The coordinate indices of ``L``'s basis spanning the ideal ``M = ⟨M-part
    arrows⟩`` (the kernel of ``π: L ↠ B``). ``M·M = 0`` and the presented ``T(B)``
    basis is corner-homogeneous, so ``M`` is spanned by a SUBSET of ``L``'s standard
    basis — asserted coordinate-aligned (loud otherwise)."""
    dom = L.domain
    m = L.dim
    vanishing = [a for a in L.quiver.arrows if a not in B.quiver.arrows]
    label_index = {lab: i for i, lab in enumerate(L.basis_labels)}
    mset = set(label_index[a] for a in vanishing)
    changed = True
    while changed:                                # close the two-sided ideal
        changed = False
        for i in list(mset):
            ei = L._basis_vec(i)
            for j in range(m):
                bj = L._basis_vec(j)
                for w in (L.multiply(bj, ei), L.multiply(ei, bj)):
                    for t, x in enumerate(w):
                        if not dom.is_zero(x) and t not in mset:
                            mset.add(t)
                            changed = True
    return mset


# --------------------------------------------------------------------------- #
# the CS block-partition snake
# --------------------------------------------------------------------------- #
def _cols(D):
    return [[row[c] for row in D] for c in range(len(D[0]))] if D and D[0] else []


def _rank(D, dom):
    return rank(D, dom) if D and D[0] else 0


def _reps_mod_image(cycles, image, dom):
    reps, base = [], list(image)
    base_rank = rank(base, dom) if base else 0
    for v in cycles:
        rr = rank(base + [v], dom)
        if rr > base_rank:
            reps.append(list(v))
            base = base + [list(v)]
            base_rank = rr
    return reps


def _apply(mat, vec, dom):
    if not mat or not mat[0]:
        return []
    return [sum((mat[r][c] * vec[c] for c in range(len(vec))), dom.zero())
            for r in range(len(mat))]


class _Blocks:
    """The CS regular Hom-complex block-partitioned by the coefficient split
    ``L = M ⊕ B`` (M = ideal). Holds, per degree, the ``δ_M`` / ``δ_B`` / ``X``
    (coh) or ``b_M`` / ``b_B`` / ``Y`` (hom) blocks and the block dims."""

    def __init__(self, L, B, top, side, max_cells):
        from quiverlab.resolutions_cs.build import reduction_system_of
        from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
        from quiverlab.resolutions_cs.homology import _require_admissible
        self.dom = L.domain
        self.side = side
        self.top = top
        mset = _mblock_indices(L, B)
        rs = reduction_system_of(L)
        _require_admissible(rs)
        # coh needs matrices to top+1 (M-flank -> top+1); hom needs to top+2
        # (B-flank -> top+1 reads b_{top+2}); build one degree past that.
        self.res = ChouhySolotarResolution(L, rs, max_degree=top + 2,
                                           max_cells=max_cells, coefficients=None)
        self.res.assert_dd_zero(upto=top + 1, side=side)
        self.res.assert_order_condition(upto=top + 1)
        self._mset = mset
        self._mask_cache = {}

    def _masks(self, n):
        if n not in self._mask_cache:
            basis = self.res._basis(n, self.side)
            mi = [i for i, (ch, j) in enumerate(basis) if j in self._mset]
            bi = [i for i, (ch, j) in enumerate(basis) if j not in self._mset]
            self._mask_cache[n] = (mi, bi)
        return self._mask_cache[n]

    def dim_M(self, n):
        return len(self._masks(n)[0])

    def dim_B(self, n):
        return len(self._masks(n)[1])

    def dim_L(self, n):
        return self.res.dim_C(n, self.side)

    @staticmethod
    def _sub(D, rows, cols):
        return [[D[r][c] for c in cols] for r in rows]

    def coh_blocks(self, n):
        """δ^n = [[δ_M, X],[0, δ_B]] with rows = C^{n+1}, cols = C^n. Returns
        (δ_M, δ_B, X) and asserts the (B-rows, M-cols) block is zero."""
        D = self.res.matrix(n, "coh")
        mr, br = self._masks(n + 1)
        mc, bc = self._masks(n)
        off = self._sub(D, br, mc)
        if any(not self.dom.is_zero(x) for row in off for x in row):
            raise QuiverlabError(
                f"split-extension LES: M is not a sub-bimodule at cochain degree "
                f"{n} (off-diagonal block nonzero) — the ideal partition is wrong",
                hint="internal invariant; please report this presentation")
        return (self._sub(D, mr, mc), self._sub(D, br, bc), self._sub(D, mr, bc))

    def hom_blocks(self, n):
        """b_n = [[b_M, Y],[0, b_B]] with rows = C_{n-1}, cols = C_n. Returns
        (b_M, b_B, Y) and asserts the (B-rows, M-cols) block is zero."""
        D = self.res.matrix(n, "hom")
        mr, br = self._masks(n - 1)
        mc, bc = self._masks(n)
        off = self._sub(D, br, mc)
        if any(not self.dom.is_zero(x) for row in off for x in row):
            raise QuiverlabError(
                f"split-extension homology: M is not a sub-bimodule at chain degree "
                f"{n} (off-diagonal block nonzero)",
                hint="internal invariant; please report this presentation")
        return (self._sub(D, mr, mc), self._sub(D, br, bc), self._sub(D, mr, bc))


def _cohomology_snake(L, B, top, max_cells):
    """Assemble HH^•(L) via the LES on the CS complex. Returns
    (flank_M[0..top+1], flank_B[0..top], delta[0..top], mid[0..top], assembled)."""
    bl = _Blocks(L, B, top, "coh", max_cells)
    dom = bl.dom
    DM, DB, X = {}, {}, {}
    for n in range(top + 2):
        dm, db, x = bl.coh_blocks(n)
        DM[n], DB[n], X[n] = dm, db, x
    # flank_M to top+1, flank_B to top
    flank_M, prev = [], 0
    for n in range(top + 2):
        r = _rank(DM[n], dom)
        flank_M.append(bl.dim_M(n) - r - prev)
        prev = r
    flank_B, prev = [], 0
    for n in range(top + 1):
        r = _rank(DB[n], dom)
        flank_B.append(bl.dim_B(n) - r - prev)
        prev = r
    mid, prev = [], 0
    for n in range(top + 1):
        r = _rank(bl.res.matrix(n, "coh"), dom)
        mid.append(bl.dim_L(n) - r - prev)
        prev = r
    # snake δ^n : HH^n(B) -> HH^{n+1}(M) via the coupling X[n]
    delta = []
    for n in range(top + 1):
        ker = (nullspace(DB[n], dom) if DB[n] and DB[n][0]
               else [[dom.one() if i == j else dom.zero() for j in range(bl.dim_B(n))]
                     for i in range(bl.dim_B(n))])
        im_prev = _cols(DB[n - 1]) if n > 0 else []
        b_classes = _reps_mod_image(ker, im_prev, dom)
        images = [_apply(X[n], phi, dom) for phi in b_classes]
        delta.append(_map_rank(images, _cols(DM[n]), dom))
    assembled = [(flank_M[n] - (delta[n - 1] if n > 0 else 0)) + (flank_B[n] - delta[n])
                 for n in range(top + 1)]
    exact = _ses_dim_identity(bl, "coh", top) and _delta_bounds(delta, flank_B, flank_M, top)
    return flank_M, flank_B, delta, mid, assembled, exact


def _homology_snake(L, B, top, max_cells):
    """Assemble HH_•(L) via the dual LES on the CS complex. Returns
    (flank_M[0..top], flank_B[0..top+1], delta[0..top+1], mid[0..top], assembled)."""
    bl = _Blocks(L, B, top, "hom", max_cells)
    dom = bl.dom
    bM, bB, Y = {}, {}, {}
    for n in range(1, top + 3):
        m_, b_, y_ = bl.hom_blocks(n)
        bM[n], bB[n], Y[n] = m_, b_, y_
    def rkM(n):
        return _rank(bM[n], dom) if n >= 1 else 0
    def rkB(n):
        return _rank(bB[n], dom) if n >= 1 else 0
    flank_M = [bl.dim_M(n) - rkM(n) - rkM(n + 1) for n in range(top + 1)]
    flank_B = [bl.dim_B(n) - rkB(n) - rkB(n + 1) for n in range(top + 2)]
    mid = []
    for n in range(top + 1):
        rb_n = _rank(bl.res.matrix(n, "hom"), dom) if n >= 1 else 0
        rb_n1 = _rank(bl.res.matrix(n + 1, "hom"), dom)
        mid.append(bl.dim_L(n) - rb_n - rb_n1)
    # snake δ_n : HH_n(B) -> HH_{n-1}(M) via the coupling Y[n]
    delta = [0] * (top + 2)
    for n in range(1, top + 2):
        ker = (nullspace(bB[n], dom) if bB[n] and bB[n][0]
               else [[dom.one() if i == j else dom.zero() for j in range(bl.dim_B(n))]
                     for i in range(bl.dim_B(n))])
        im_next = _cols(bB[n + 1]) if n + 1 <= top + 2 else []
        b_classes = _reps_mod_image(ker, im_next, dom)
        images = [_apply(Y[n], z, dom) for z in b_classes]
        delta[n] = _map_rank(images, _cols(bM[n]) if n >= 1 else [], dom)
    assembled = [(flank_M[n] - delta[n + 1]) + (flank_B[n] - delta[n])
                 for n in range(top + 1)]
    exact = _ses_dim_identity(bl, "hom", top) and _delta_bounds_hom(delta, flank_B, flank_M, top)
    return flank_M, flank_B, delta[: top + 1], mid, assembled, exact


def snake_lift_well_defined(B, top, *, M=None, max_cells=4_000_000):
    """Adversarial-lift self-cert (Task A2 Step 3): the snake ``δ^n[φ]`` is
    independent of the lift chosen — shifting the canonical lift ``(0, φ)`` by ANY
    ``M``-cochain ``m ∈ C^n(L,M) = ker π`` leaves ``δ^n[φ]``'s class unchanged
    (because ``δ_reg(m, φ) = (δ_M m + Xφ, 0)`` and ``δ_M m`` is an ``M``-coboundary).
    Returns True iff well-definedness holds at every degree with every basis shift."""
    L = split_extension(B, M)
    bl = _Blocks(L, B, top, "coh", max_cells)
    dom = bl.dom
    for n in range(top + 1):
        dm, db, x = bl.coh_blocks(n)
        ker = (nullspace(db, dom) if db and db[0]
               else [[dom.one() if i == j else dom.zero() for j in range(bl.dim_B(n))]
                     for i in range(bl.dim_B(n))])
        im_dm = _cols(dm)
        for phi in ker:
            base_img = _apply(x, phi, dom)
            r0 = rank(im_dm + [base_img], dom) if base_img else (rank(im_dm, dom) if im_dm else 0)
            # shift by every M-cochain basis vector e_k in C^n(L,M)
            for k in range(bl.dim_M(n)):
                shift = [dom.one() if t == k else dom.zero() for t in range(bl.dim_M(n))]
                dm_shift = _apply(dm, shift, dom)
                shifted = [base_img[t] + dm_shift[t] for t in range(len(base_img))] \
                    if base_img else dm_shift
                r1 = rank(im_dm + [shifted], dom) if shifted else r0
                if r1 != r0:
                    return False
    return True


def _map_rank(images, image_space, dom):
    """rank of the induced map to the quotient (target-space) / span(image_space)."""
    base = list(image_space)
    br = rank(base, dom) if base else 0
    got = 0
    for v in images:
        if not v:
            continue
        r2 = rank(base + [v], dom)
        if r2 > br:
            got += 1
            base.append(v)
            br = r2
    return got


def _ses_dim_identity(bl, side, top):
    """dim C^n(L,M) + dim C^n(L,B) == dim C^n(L,L) at every degree (the degreewise
    SES-of-complexes self-cert)."""
    for n in range(top + 1):
        if bl.dim_M(n) + bl.dim_B(n) != bl.dim_L(n):
            return False
    return True


def _delta_bounds(delta, flank_B, flank_M, top):
    for n in range(top + 1):
        if not (0 <= delta[n] <= min(flank_B[n], flank_M[n + 1])):
            return False
    return True


def _delta_bounds_hom(delta, flank_B, flank_M, top):
    for n in range(1, top + 2):
        hi = min(flank_B[n], flank_M[n - 1])
        if not (0 <= delta[n] <= hi):
            return False
    return True


# --------------------------------------------------------------------------- #
# leading pieces (p = 0) over the smaller algebra B
# --------------------------------------------------------------------------- #
def _leading_pieces(B, top, side):
    from quiverlab.hochschild.coefficients import Bimodule
    dualB = Bimodule.dual(B)
    if side == "coh":
        hhB = B.hochschild_cohomology(top, verbose=False).dims
        hBM = B.hochschild_cohomology(top, coefficients=dualB, verbose=False).dims
    else:
        hhB = B.hochschild_homology(top, verbose=False).dims
        hBM = B.hochschild_homology(top, coefficients=dualB, verbose=False).dims
    return hhB, hBM


# --------------------------------------------------------------------------- #
# the HH^1(L) != 0 grading-derivation witness (CMRS Thm 5.5)
# --------------------------------------------------------------------------- #
def hh1_grading_witness(B, *, M=None):
    """The ``HH^1(L) ≠ 0`` witness (CMRS Thm 5.5): the grading derivation
    ``E: L → L``, ``E(b+m) = m`` (``E|_B = 0``, ``E|_M = id``), materialized as an
    explicit degree-1 bar cocycle of ``L`` and certified **outer** (``[E] ≠ 0`` in
    ``HH^1(L)``). For **directed** ``B`` (acyclic quiver) the sharper
    ``HH^1(T(B)) = k ⊕ HH^1(B)`` is reported. Returns a dict with ``is_cocycle``,
    ``is_outer``, ``nonzero``, ``hh1_L``, ``hh1_B``, ``directed``, ``formula``."""
    from quiverlab.hochschild import bar
    L = split_extension(B, M)
    dom = L.domain
    m = L.dim
    mset = _mblock_indices(L, B)
    Lu = L.unit_adapted()
    P = L._unit_adapting_change()                 # columns = Lu basis in L coords
    # E in Lu's basis: projection onto the (transported) M-ideal.
    E = [[dom.zero()] * m for _ in range(m)]
    for j in range(m):
        pj = [P[r][j] for r in range(m)]
        masked = [pj[r] if r in mset else dom.zero() for r in range(m)]
        x = solve(P, masked, dom)
        if x is None:
            raise QuiverlabError(
                "hh1_grading_witness: unit-adapting change is singular",
                hint="internal invariant; please report this presentation")
        for r in range(m):
            E[r][j] = x[r]
    basis1 = bar._cochain_basis(m, 1, m)          # (s, (a,)), s in 0..m-1, a in 1..m-1
    e_vec = [E[s][a] for (s, (a,)) in basis1]
    d0, _, _ = bar.coboundary_matrix(Lu, 0, 4_000_000, None)
    d1, _, _ = bar.coboundary_matrix(Lu, 1, 4_000_000, None)
    prod = _apply(d1, e_vec, dom)
    is_cocycle = all(dom.is_zero(x) for x in prod)
    inner_cols = _cols(d0)                         # im δ^0 = inner derivations
    is_outer = solve([[inner_cols[k][r] for k in range(len(inner_cols))]
                      for r in range(len(e_vec))], e_vec, dom) is None if inner_cols \
        else any(not dom.is_zero(x) for x in e_vec)
    hh1_L = len(e_vec) - _rank(d1, dom) - _rank(d0, dom)
    hh1_B = B.hochschild_cohomology(1, verbose=False).dims[1]
    directed = bool(B.quiver.is_acyclic())
    formula = (hh1_L == 1 + hh1_B) if directed else None
    return {"is_cocycle": is_cocycle, "is_outer": is_outer,
            "nonzero": bool(is_cocycle and is_outer), "hh1_L": hh1_L,
            "hh1_B": hh1_B, "directed": directed, "formula": formula}


# --------------------------------------------------------------------------- #
# reports
# --------------------------------------------------------------------------- #
def _split_extension_report(B, top, side, M, max_cells):
    L = split_extension(B, M)
    coeff_name = "D(B)"
    if side == "coh":
        flank_M, flank_B, delta, mid, assembled, exact = _cohomology_snake(L, B, top, max_cells)
        leading_B, leading_M = _leading_pieces(B, top, "coh")
        side_name = "cohomology"
    else:
        flank_M, flank_B, delta, mid, assembled, exact = _homology_snake(L, B, top, max_cells)
        leading_B, leading_M = _leading_pieces(B, top, "hom")
        side_name = "homology"
    wit = hh1_grading_witness(B, M=M)
    agrees = assembled == mid
    return SplitExtReport(
        base=B, extension=L, coeff_name=coeff_name, side=side_name, n_range=top,
        flank_M=flank_M, flank_B=flank_B, delta_ranks=delta, assembled=assembled,
        direct=mid, leading_B=leading_B, leading_M=leading_M,
        hh1_nonzero=wit["nonzero"], hh1_directed_formula=wit["formula"],
        exact=bool(exact), agrees=bool(agrees), status="complete", note="")


def split_extension_cohomology(B, top, *, M=None, max_cells=4_000_000):
    """The split-extension cohomology LES assembly + the direct cross-check."""
    return _split_extension_report(B, top, "coh", M, max_cells)


def split_extension_homology(B, top, *, M=None, max_cells=4_000_000):
    """The split-extension homology twin (the dual SES ``0 → M → L → B → 0`` in
    ``HH_•(L, −)``)."""
    return _split_extension_report(B, top, "hom", M, max_cells)


# --------------------------------------------------------------------------- #
# GUI / webapp block builder (algebra-only, budget-carrying kind)
# --------------------------------------------------------------------------- #
_SPLIT_REFERENCES = ("cmrs_split", "crs_trivial_ext_hh1")


def split_extension_block(A, top=6, max_cells=4_000_000):
    """The ``split_extension`` compute-kind block (all three tiers). Interprets the
    drawn algebra ``A`` as ``B``, runs the cohomology LES, and returns the structured
    decomposition (flanks, leading pieces, assembled vs. direct, HH^1 witness). Loud
    refusals are turned into a clean ``status``/``error`` block (never a 500)."""
    block = {"kind": "split_extension", "references": list(_SPLIT_REFERENCES)}
    try:
        rep = split_extension_cohomology(A, top, max_cells=max_cells)
    except QuiverlabError as exc:
        block.update(status="unsupported", error=str(exc),
                     note="split-extension LES unavailable for this input")
        return block
    block.update(
        status=rep.status,
        base=repr(A).splitlines()[0],
        extension=repr(rep.extension).splitlines()[0],
        side=rep.side,
        top=rep.n_range,
        coeff=rep.coeff_name,
        flank_M=rep.flank_M, flank_B=rep.flank_B,
        delta_ranks=rep.delta_ranks,
        assembled=rep.assembled, direct=rep.direct,
        leading_B=rep.leading_B, leading_M=rep.leading_M,
        hh1_nonzero=rep.hh1_nonzero, hh1_directed_formula=rep.hh1_directed_formula,
        exact=rep.exact, agrees=rep.agrees, note=rep.note,
    )
    return block
