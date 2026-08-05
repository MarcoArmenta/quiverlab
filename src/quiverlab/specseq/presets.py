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
from quiverlab.specseq.double import DoubleComplex
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
