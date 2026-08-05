"""The four Plan-42 presets over the spectral-sequence engine.

* :func:`hochschild_bB_ss` -- the Hochschild ``(b, B)`` bicomplex (this file).
* :func:`radical_filtration_ss` -- the associated-graded / radical filtration
  (Task 6).
* :func:`grothendieck_double_complex` / :func:`cartan_eilenberg_ss` -- the
  Grothendieck / Cartan-Eilenberg change-of-rings SS (Task 7).

Every construction returns a pre-certified :class:`SpectralSequence` (the standing
``E_inf`` totals == total-homology self-certificate runs at ``__init__``). No
floats (``src/`` AST gate)."""
from quiverlab.errors import DepthLimitError
from quiverlab.modules import linalg_mod as lm
from quiverlab.specseq import _subspace as sub
from quiverlab.specseq.double import DoubleComplex
from quiverlab.specseq.filtered import FilteredComplex
from quiverlab.specseq.pages import SpectralSequence

_GUARD_HINT = ("the (b, B) bar bicomplex is exponential (dim C_n = m*(m-1)^n); "
               "raise max_cells or lower top only if you know what you are doing")


def _guard_bB_cells(m, top, max_cells):
    """Refuse LOUDLY, by length arithmetic (never enumeration), before any basis or
    matrix is materialized -- mirrors ``engine/cyclic._guard_cyclic_cells``."""
    mr = m - 1

    def cdim(c):
        return m * mr ** c

    maxc = top + 1
    for c in range(1, maxc + 1):                    # bar boundary b_c
        cells = cdim(c - 1) * cdim(c)
        if cells > max_cells:
            raise DepthLimitError(
                f"Hochschild (b,B) SS: bar boundary b_{c} would have {cells} cells "
                f"(> max_cells = {max_cells})", hint=_GUARD_HINT)
    for c in range(0, maxc):                          # Connes B_c
        cells = cdim(c + 1) * cdim(c)
        if cells > max_cells:
            raise DepthLimitError(
                f"Hochschild (b,B) SS: Connes B_{c} would have {cells} cells "
                f"(> max_cells = {max_cells})", hint=_GUARD_HINT)


def hochschild_bB_ss(A, top, max_cells=4_000_000):
    """The Hochschild ``(b, B)`` spectral sequence of ``A``, abutting to cyclic
    homology ``HC_*(A)``.

    The first-quadrant ``(b, B)`` bicomplex: ``D_{p,q} = C_{q-p}`` (``q >= p >= 0``,
    the normalized bar chains on the unit-adapted basis), vertical ``d_v = b``
    (``C_c -> C_{c-1}``, lowers ``q``), horizontal ``d_h = B`` (Connes, ``C_c ->
    C_{c+1}``, lowers ``p``). The two differentials ALREADY anticommute -- the
    mixed-complex identity ``bB + Bb = 0`` -- so no sign adjustment is needed (the
    DoubleComplex gate confirms it). Its total complex is ``engine/cyclic``'s
    ``Tot_n = (+)_p C_{n-2p}``, so ``E_inf`` total == ``HC_n`` (the cross-engine
    arbiter of the exact ``(p, q)`` layout).

    The exponential bar basis (``dim C_n = m*(m-1)^n``) is guarded up front by
    length arithmetic: an over-``max_cells`` request raises ``DepthLimitError``
    before any basis is materialized. Works over any exact Domain (the generic
    ``hochschild.cyclic`` / ``hochschild.bar`` matrices)."""
    from quiverlab.hochschild.bar import boundary_matrix
    from quiverlab.hochschild.cyclic import connes_B_matrix
    B0 = A.unit_adapted()
    dom = B0.domain
    m = B0.dim
    _guard_bB_cells(m, top, max_cells)
    maxc = top + 1
    cdim = {c: m * (m - 1) ** c for c in range(maxc + 1)}
    bmats = {c: boundary_matrix(B0, c, max_cells)[0] for c in range(1, maxc + 1)}
    Bmats = {c: connes_B_matrix(B0, c, max_cells)[0] for c in range(0, maxc)}
    terms = {}
    for p in range(0, maxc + 1):
        for q in range(p, maxc + 1):
            if p + q > maxc:
                continue
            d = cdim[q - p]
            if d:
                terms[(p, q)] = d
    d_h, d_v = {}, {}
    for (p, q) in terms:
        c = q - p
        if c >= 1 and (p, q - 1) in terms:            # b_c : D_{p,q} -> D_{p,q-1}
            d_v[(p, q)] = bmats[c]
        if p >= 1 and (p - 1, q) in terms:            # B_c : D_{p,q} -> D_{p-1,q}
            d_h[(p, q)] = Bmats[c]
    dc = DoubleComplex(terms, d_h, d_v, dom, check=True)
    return SpectralSequence(dc.column_filtration())


# --------------------------------------------------------------------------- #
# Preset 3 -- associated-graded / radical filtration (Task 6).
# --------------------------------------------------------------------------- #
def _rad_powers(M):
    """The descending radical chain ``[rad^0 M, rad^1 M, ..., rad^{L-1} M]`` as
    column-span bases INSIDE ``M``'s coordinates (``rad^0 = M`` whole, stopping just
    before the first zero power). ``rad^{i} M = sum_a action[a](rad^{i-1} M)`` -- the
    single-step ``radtopsoc.radical`` definition, iterated in place so every power
    stays in ``M``'s basis (``radical()`` alone returns a fresh submodule)."""
    dom = M.domain
    d = M.dim
    if d == 0:
        return [[]]
    arrows = list(M.algebra.quiver.arrows)
    ident = lm.identity(d, dom)
    full = [lm.col(ident, j) for j in range(d)]
    powers, cur = [full], full
    while arrows:
        gens = [lm.matvec(M.action[a], c, dom) for a in arrows for c in cur]
        nxt = sub.reduce_to_independent(gens, dom)
        if not nxt or sub.span_dim(nxt, dom) >= sub.span_dim(cur, dom):
            break                      # rad^i = 0 (nilpotent) -- or a defensive no-shrink
        powers.append(nxt)
        cur = nxt
    return powers


def radical_filtration_ss(X):
    """The associated-graded (radical-filtration) spectral sequence of a P39
    :class:`~quiverlab.modules.complexes.ChainComplex` ``X`` of ``A``-modules.

    ``F_p X_n = X_n . rad^{max(0,-p)}`` -- an INCREASING (in ``p``), exhaustive
    (``F_0 = X_n``), Hausdorff (``rad`` nilpotent so ``F_{-large} = 0``) subcomplex
    filtration (the differentials are module maps, hence preserve the radical
    powers). Converges to ``H_*(X)`` by the standing self-certificate. A complex of
    semisimple modules (``rad = 0``) gives the trivial one-step filtration
    (collapse at ``E_1``); for a Koszul algebra the minimal resolution is linear and
    the associated-graded sequence degenerates early -- the exact page is arbitrated
    per instance (never forced), see the verification page."""
    degs = X.degrees()
    powers = {n: _rad_powers(X.term(n)) for n in degs}
    maxL = max((len(powers[n]) for n in degs), default=1)
    lo = -(maxL - 1)
    filt = {}
    for n in degs:
        pw, Ln = powers[n], len(powers[n])
        levels = []
        for j in range(maxL):               # p = lo + j; radical exponent i = -p = maxL-1-j
            i = (maxL - 1) - j
            levels.append([list(c) for c in pw[i]] if i < Ln else [])
        filt[n] = levels
    return SpectralSequence(FilteredComplex.from_chain_complex(X, filt, lo=lo))


# --------------------------------------------------------------------------- #
# Presets 1 & 2 -- Grothendieck / Cartan-Eilenberg change-of-rings (Task 7).
# --------------------------------------------------------------------------- #
def _flat(mat):
    return [x for row in mat for x in row]


def _coords_in_hombasis(H, psi, dom):
    """Coordinates of the module map ``psi`` in the ``hom_space`` basis ``H`` (each a
    matrix of the same shape). Loud if ``psi`` is not in the span (not a module
    map / block-indexing bug)."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.fields import linalg as flinalg
    if not H:                                   # target Hom is 0-dimensional
        return []
    B = lm.cols_to_matrix([_flat(h) for h in H])
    x = flinalg.solve(B, _flat(psi), dom)
    if x is None:
        raise QuiverlabError(
            "Grothendieck double complex: a differential image is not expressible "
            "in the target Hom basis (block indexing bug)")
    return x


def _assert_change_of_rings(A, B):
    """``B = A/I'`` an admissible quotient: same quiver, ``rel(A) subset rel(B)``.
    The restriction of a ``B``-module to ``A`` (identical matrices) is valid exactly
    then (Design-decision 2; general change-of-rings is out of scope this release)."""
    from quiverlab.errors import QuiverlabError
    if A.quiver is None or B.quiver is None:
        raise QuiverlabError(
            "Cartan-Eilenberg SS needs quiver-presented A and B")
    if (list(A.quiver.vertices) != list(B.quiver.vertices)
            or A.quiver.arrows != B.quiver.arrows):
        raise QuiverlabError(
            "Cartan-Eilenberg SS: A and B must share the same quiver "
            "(B = A/I' an admissible quotient); full-generality change-of-rings for "
            "an arbitrary algebra map is out of scope this release")
    relA = {str(r) for r in (A.relations or [])}
    relB = {str(r) for r in (B.relations or [])}
    if not relA <= relB:
        raise QuiverlabError(
            "Cartan-Eilenberg SS: rel(A) must be a subset of rel(B) (B = kQ/I' with "
            "I subset of I'); the restriction of a B-module to A is valid only then")


def _restrict_to_A(M, A):
    """View a right ``B``-module ``M`` as a right ``A``-module (``A``, ``B`` share the
    quiver, ``rel(A) subset rel(B)``): rebuild from the COMMON arrow actions via
    ``from_arrow_action`` so ``A``'s longer paths are filled in by composition
    (``Module(A, M.dim, M.action)`` alone drops them and breaks ``A``-side resolutions
    -- e.g. ``a*b`` is nonzero in ``A`` but absent from a ``B = kQ/(a*b)``-module's
    action)."""
    from quiverlab.errors import QuiverlabError
    if M.side != "right":
        raise QuiverlabError(
            "Cartan-Eilenberg SS: a left-module M is out of scope this release "
            "(the change-of-rings preset takes a right B-module M)")
    arrow_action = {a: M.action[a] for a in A.quiver.arrows}
    return A.module(M.dimension_vector(), arrow_action, side="right", name=M.name)


def grothendieck_double_complex(M, U, N, p_len, q_len, max_term_dim=200000):
    """The Cartan-Eilenberg / Grothendieck change-of-rings double complex for the
    ``U = B`` case (the adjunction collapse
    ``Hom_B(Q_p, Hom_A(B, J^q)) = Hom_A(res_A Q_p, J^q)``).

    ``M`` is a ``B``-module, ``N`` an ``A``-module (``A = N.algebra``,
    ``B = M.algebra``); ``U`` must be ``B`` (the algebra itself -- the change-of-rings
    bimodule). ``D^{p,q} = Hom_A(res_A Q_p, J^q(N))`` where ``Q_.`` is the minimal
    ``B``-projective resolution of ``M`` (length ``p_len``) and ``J^.`` the minimal
    ``A``-injective coresolution of ``N`` (length ``q_len``). COHOMOLOGICAL, stored
    with NEGATED total degree (position ``(-p, -q)``, the P39 ``C^n := C_{-n}``
    discipline) so it lands in the homological :class:`DoubleComplex`; horizontal =
    pre-compose with ``d^Q`` (raises ``p``), vertical = ``(-1)^p .`` post-compose with
    ``d^J`` (raises ``q``, Koszul sign so the square anticommutes -- the gate confirms
    it). The abutment is certified per instance by :func:`cartan_eilenberg_ss`."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.hom import hom_space
    from quiverlab.modules.injective import injective_resolution
    from quiverlab.modules.module import Module
    from quiverlab.modules.resolution import minimal_resolution
    A, B = N.algebra, M.algebra
    if U is not B:
        raise QuiverlabError(
            "grothendieck_double_complex: only the change-of-rings U = B case is "
            "implemented this release (general (B,A)-bimodule U is deferred); pass "
            "U = M.algebra")
    _assert_change_of_rings(A, B)
    dom = A.domain
    terms_Q, dmats_Q = minimal_resolution(M, p_len, max_term_dim=max_term_dim)
    resA = {}                                    # p -> res_A Q_p (an A-module) or None
    for p, td in enumerate(terms_Q):
        Qp = td.module
        if Qp is not None and Qp.dim > 0:
            resA[p] = Module(A, Qp.dim, Qp.action, name=f"resA_Q{p}", side=Qp.side)
    inj = injective_resolution(N, q_len, max_term_dim=max_term_dim)
    Jt = {}                                      # q -> J^q (an A-module) or None
    for q, E in enumerate(inj.terms):
        if E is not None and E.dim > 0:
            Jt[q] = E
    # Hom basis per (p, q)
    Hb = {}
    for p in resA:
        for q in Jt:
            Hb[(p, q)] = hom_space(resA[p], Jt[q])
    terms, d_h, d_v = {}, {}, {}
    for (p, q), H in Hb.items():
        if H:
            terms[(-p, -q)] = len(H)             # negated storage
    for (p, q), H in Hb.items():
        if not H:
            continue
        # horizontal: pre-compose with d^Q_{p+1} : Q_{p+1} -> Q_p  (raises p)
        Ht = Hb.get((p + 1, q))
        if Ht and (p + 1) in resA:
            dQ = dmats_Q[p + 1] if (p + 1) < len(dmats_Q) else None
            if dQ and dQ[0]:
                cols = [_coords_in_hombasis(Ht, lm.matmul(phi, dQ, dom), dom)
                        for phi in H]
                mat = lm.cols_to_matrix(cols)
                if mat and mat[0]:
                    d_h[(-p, -q)] = mat
        # vertical: (-1)^p * post-compose with d^J_{q+1} : J^q -> J^{q+1} (raises q)
        Hd = Hb.get((p, q + 1))
        if Hd and (q + 1) in Jt:
            dJ = inj.differential(q + 1) if (q + 1) < len(inj.dmats) else None
            if dJ and dJ[0]:
                sign = dom.one() if p % 2 == 0 else dom.neg(dom.one())
                cols = []
                for phi in H:
                    c = _coords_in_hombasis(Hd, lm.matmul(dJ, phi, dom), dom)
                    cols.append([dom.mul(sign, x) for x in c])
                mat = lm.cols_to_matrix(cols)
                if mat and mat[0]:
                    d_v[(-p, -q)] = mat
    return DoubleComplex(terms, d_h, d_v, dom, check=True)


def cartan_eilenberg_ss(A, B, M, N, p_len=6, q_len=6):
    """The Cartan-Eilenberg change-of-rings spectral sequence for an ADMISSIBLE
    QUOTIENT ``B = A/I'`` (``A``, ``B`` share the quiver, ``rel(A) subset rel(B)``):
    ``E_2^{p,q} = Ext_B^p(M, Ext_A^q(B, N)) => Ext_A^{p+q}(M|_A, N)``.

    ``M`` is a ``B``-module, ``N`` an ``A``-module. Returns a pre-certified
    :class:`SpectralSequence` (the standing ``E_inf`` totals == total-homology
    self-cert runs at construction). The ABUTMENT is additionally CERTIFIED per
    instance against the module ``Ext`` cross-engine oracle over the window
    ``[0, min(p_len, q_len) - 1]``: if ``E_inf`` total (at the negated homological
    degree) does not equal ``Ext_A^n(M|_A, N)``, the per-instance change-of-rings
    hypothesis fails (or the truncation is too shallow) and it refuses loudly --
    NEVER a wrong abutment. Degenerate ``B = A`` collapses at ``E_2`` (``Ext_A^q(A,
    N) = N`` for ``q = 0``, else ``0``)."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.ext import ext_dims
    _assert_change_of_rings(A, B)
    dc = grothendieck_double_complex(M, B, N, p_len, q_len)
    ss = SpectralSequence(dc.column_filtration())
    # per-instance abutment certification (the CE hypothesis / truncation depth):
    window = min(p_len, q_len)
    M_over_A = _restrict_to_A(M, A)
    ext = ext_dims(A, M_over_A, N, max(window - 1, 0))
    totals = ss.total_homology_dims                       # keyed by NEGATED degree
    for n in range(window):
        got = totals.get(-n, 0)
        want = ext[n] if n < len(ext) else 0
        if got != want:
            raise QuiverlabError(
                "Cartan-Eilenberg SS: the change-of-rings acyclicity hypothesis is "
                f"not certified for this instance at degree {n} -- the abutment "
                f"E_inf total {got} != Ext_A^{n}(M|_A, N) = {want} (either "
                "Ext_B^{>0}(M, Hom_A(B, J^q)) != 0, or the truncation p_len/q_len is "
                "too shallow to certify this degree)",
                hint="raise p_len/q_len, or the instance falls outside the "
                     "admissible-quotient change-of-rings scope")
    return ss
