"""Finite repetitive-algebra slices (Plan 44 / C7).

The repetitive algebra ``hat(A)`` (Hughes-Waschbuesch) is the doubly-infinite matrix
algebra with ``A`` on the diagonal and the connecting bimodule ``D(A)`` on the
super-diagonal -- ``INFINITE`` dimensional. This ships certified FINITE slices only:
``copies`` copies of ``A`` (``copies >= 1``) joined by ``copies - 1`` copies of the
``D(A)`` connecting bimodule, ``dim == (2*copies - 1)*dim A``.

``copies == 1`` is ``A`` itself (byte-Cartan base case). For ``copies >= 2`` the target is
the block structure-constant algebra ``hat(A)_slice`` (``A`` on the diagonal blocks,
``D(A)`` on the connecting blocks between CONSECUTIVE copies, NOT wrapping the last copy
back to the first -- that wrap would be the trivial extension). It is CERTIFIED per instance
by ``check=True`` (associativity/unit) plus the dimension identity; for a string-representable
base over QQ / GF(p) it is additionally PRESENTED as a genuine ``kQ_slice/I`` via the
``TrivialExtension`` connecting-arrow idiom (so path-basis invariants serve it), falling back
to the structure-constant build otherwise. Float-free.

HONEST SCOPE: only finite slices are shipped; the full ``hat(A)`` and the Hughes-Waschbuesch
orbit identity ``T(A) ~ hat(A)/(nu)`` (a slice with periodic identification) are named
successors, not shipped here. Refuses ``copies < 1`` loudly.
"""
from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError

_CITATIONS = ("happel_trivial_extension", "hughes_waschbusche", "assem_book")


def _repetitive_structure_constants(A, copies):
    """``hat(A)_slice`` as raw structure constants (dim ``(2*copies - 1)*dim A``). Block
    layout: ``A_0, DA_0, A_1, DA_1, ..., A_{copies-1}`` (block ``2c`` = ``A_c``, block
    ``2c+1`` = ``DA_c``). ``DA_c`` is the ``(A_c, A_{c+1})``-bimodule ``D(A)`` (left action
    by copy ``c``, right by copy ``c+1``); the single-copy table comes from
    ``TrivialExtension``'s ``A |x D(A)`` build, spread across consecutive copies."""
    from quiverlab.families.trivial_extension import _trivial_extension_structure_constants
    dom = A.domain
    n = A.dim
    zero = dom.zero()
    Tsc = _trivial_extension_structure_constants(A)          # dim 2n: A-part 0..n-1, DA-part n..2n-1
    m = (2 * copies - 1) * n

    def offA(c):
        return (2 * c) * n

    def offDA(c):
        return (2 * c + 1) * n

    # DA_c connects copy c+1 -> copy c (the reversed connecting arrow), so it is the
    # (A_{c+1}, A_c)-bimodule: LEFT action by copy c+1, RIGHT action by copy c.
    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    for c in range(copies):                                  # A_c * A_c = A-mult
        for i in range(n):
            for k in range(n):
                pv = Tsc.T[i][k]
                for t in range(n):
                    T[offA(c) + i][offA(c) + k][offA(c) + t] = pv[t]
    for c in range(copies - 1):                              # A_{c+1} * DA_c = left action
        for i in range(n):
            for k in range(n):
                pv = Tsc.T[i][n + k]
                for t in range(n):
                    T[offA(c + 1) + i][offDA(c) + k][offDA(c) + t] = pv[n + t]
    for c in range(copies - 1):                              # DA_c * A_c = right action
        for i in range(n):
            for k in range(n):
                pv = Tsc.T[n + i][k]
                for t in range(n):
                    T[offDA(c) + i][offA(c) + k][offDA(c) + t] = pv[n + t]
    unit = [zero] * m
    for c in range(copies):
        for t in range(n):
            unit[offA(c) + t] = A.unit[t]
    R = Algebra.from_structure_constants(T, unit, field=dom, check=True)
    R._family_citations = _CITATIONS
    return R


def _presented_repetitive_slice(A, copies):
    """``kQ_slice/I`` presentation via the ``TrivialExtension`` connecting-arrow idiom:
    ``copies`` relabeled copies of ``Q_A`` + one connecting arrow per bimodule-socle vector
    of ``A`` running from copy ``c+1`` to copy ``c`` (``w in e_i A e_j -> (j, c+1) ->
    (i, c)``, reversed), between CONSECUTIVE copies. Certified by ``dim == (2*copies-1)*n``;
    raises ``QuiverlabError`` on a failed certificate (the caller falls back)."""
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.families._present import present_from_pi
    from quiverlab.families.trivial_extension import (_bimodule_socle, _corner_maps,
                                                      _socle_arrows)
    from quiverlab.invariants.pathbasis import path_type_basis
    dom = A.domain
    n = A.dim
    idem, rad, src, tgt = path_type_basis(A, "repetitive_slice")
    socle = _bimodule_socle(A, rad, dom)
    vsrc, vtgt = _corner_maps(A, idem, rad, src, tgt)
    arrows_info = _socle_arrows(A, socle, vsrc, vtgt, dom)   # [(i_vertex, j_vertex, phi)]

    R = _repetitive_structure_constants(A, copies)           # pi target (also the fallback)
    label_index = {lab: t for t, lab in enumerate(A.basis_labels)}

    def offA(c):
        return (2 * c) * n

    def offDA(c):
        return (2 * c + 1) * n

    def vlabel(v, c):
        return f"{v}__{c}"

    verts = [vlabel(v, c) for c in range(copies) for v in A.quiver.vertices]
    arrows = {}
    img = {}
    for c in range(copies):                                  # relabeled copies of Q_A
        for a, (s, t) in A.quiver.arrows.items():
            name = f"{a}__{c}"
            arrows[name] = (vlabel(s, c), vlabel(t, c))
            vec = [dom.zero()] * R.dim
            vec[offA(c) + label_index[a]] = dom.one()        # pi = length-1 path a in copy c
            img[name] = vec
    for c in range(copies - 1):                              # connecting arrows (c+1) -> c
        for s, (vi, vj, phi) in enumerate(arrows_info):
            name = f"te{s}__{c}"
            arrows[name] = (vlabel(vj, c + 1), vlabel(vi, c))
            vec = [dom.zero()] * R.dim
            for t in range(n):
                vec[offDA(c) + t] = phi[t]                    # pi = the D(A) covector in DA_c
            img[name] = vec

    Q = Quiver(verts, arrows)
    base_bound = A.loewy_length() + 2
    return present_from_pi(Q, img, R, dom, (2 * copies - 1) * n, base_bound,
                           citations=_CITATIONS)


def repetitive_slice(A, copies):
    """A finite slice of the repetitive algebra ``hat(A)``: ``copies`` copies of ``A``
    joined by ``copies - 1`` copies of ``D(A)``. Certified per instance by
    ``dim == (2*copies - 1)*dim A``. Refuses ``copies < 1`` loudly (the full repetitive
    algebra is infinite-dimensional)."""
    if not isinstance(copies, int) or isinstance(copies, bool) or copies < 1:
        raise QuiverlabError(
            f"repetitive_slice: copies must be an integer >= 1, got {copies!r}",
            hint="copies >= 1; the full repetitive algebra hat(A) is infinite-dimensional")
    if copies == 1:
        return A                                             # base case: A itself
    from quiverlab.errors import FieldError
    from quiverlab.families.trivial_extension import _is_string_representable_domain
    if _is_string_representable_domain(A.domain) and A.quiver is not None:
        try:
            return _presented_repetitive_slice(A, copies)
        except (FieldError, QuiverlabError):
            pass                                             # fall back to structure constants
    return _repetitive_structure_constants(A, copies)
