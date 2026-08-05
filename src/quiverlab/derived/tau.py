"""The derived AR translate ``tau_{D^b} = nu[-1]`` and its inverse ``tau^-_{D^b}``
on perfect complexes (Plan 43 / Happel).

The Serre functor of ``D^b(mod A)`` exists iff ``gl.dim A < infinity``; refuse loudly
otherwise (``k[x]/(x^2)`` is the pinned negative case). ``nu = D Hom_A(-, A)`` is applied
TERMWISE on a perfect (projective) complex: each projective term ``(+)P_v`` maps to
``nu`` of it (``= D Hom_A((+)P_v, A) = D((+)Ae_v) = (+)I_v``), each differential (a map
between projective terms) to its ``nu``-image between those injectives via the
corner-transpose, then the whole complex is shifted by ``-1``.

Basis discipline (the P39/P41 mismatch rule): the injective terms are built as
``dualize(Hom_A(term, A))`` -- exactly the objects ``nakayama_functor`` produces -- so
their k-basis is byte-consistent with the corner-transpose differentials; we never mix
a ``builders.injective`` basis with a corner-transpose basis. Correctness never rests on
the corner-transpose bookkeeping: the output self-certifies (``d.d = 0`` via
``ChainComplex(check=True)``) and is pinned against the trusted module ``tau`` (the
concentration + dimension-vector oracle) and the K0-action Coxeter identity."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import projective
from quiverlab.modules.complexes import ChainComplex
from quiverlab.modules.duality import dualize, _zero_module
from quiverlab.modules.module import _other_side
from quiverlab.modules.opposite import opposite_algebra
from quiverlab.modules.resolution import _direct_sum
from quiverlab.derived._corner import corner_transpose


def _require_finite_gldim(A):
    """Happel gate: no Serre functor / AR triangles / ``tau_{D^b}`` at infinite global
    dimension. ``GlobalDimension.exact`` is True exactly when every simple resolved
    within the bound (finite gl.dim); otherwise a certified lower bound -> refuse."""
    from quiverlab.modules.ext import global_dimension
    g = global_dimension(A)
    if not g.exact:                              # exact == resolved == finite here
        raise QuiverlabError(
            f"tau_Db: D^b(mod A) has no Serre functor -- global dimension is not "
            f"finite ({g!r}); the derived AR translate is undefined (Happel). "
            f"k[x]/(x^2) is the pinned negative case; needs finite gl.dim")
    return g.value


def _term_vertices(X, n, kind):
    """``[summand vertex list]`` of the perfect term ``X_n``: the provenance fast path
    (``X._proj_vertices``, set by ``from_projective_resolution`` / this module), else
    ``decompose`` + ``identify_standard`` (loud on a wrong summand; char scope
    QQ/GF(32003) for ``decompose``). ``kind`` is ``"projective"`` (nu input) or
    ``"injective"`` (nu^- input)."""
    # kind-tagged provenance first (devil's-advocate fix, 2026-08-05): the
    # fast path must never let one kind masquerade as the other.
    tagged = getattr(X, "_term_provenance", None)
    if tagged is not None:
        tkind, tmap = tagged
        if tkind != kind:
            raise QuiverlabError(
                f"this complex's terms are certified {tkind}, not {kind} -- "
                "tau_Db_minus consumes the INJECTIVE output shape of tau_Db; "
                "feed a projective complex to tau_Db instead.")
        if n in tmap:
            return list(tmap[n])
    prov = getattr(X, "_proj_vertices", None)
    if prov is not None and n in prov:
        if kind != "projective":
            raise QuiverlabError(
                "this complex carries PROJECTIVE term provenance but "
                f"{kind} terms are required -- tau_Db_minus refuses projective "
                "input loudly (the inverse translate needs injective terms).")
        return list(prov[n])
    from quiverlab.modules.decompose import decompose
    from quiverlab.modules.hom import identify_standard
    verts = []
    for s, _mult in decompose(X.term(n)):
        std = identify_standard(s)
        if not std or std[0] != kind:
            raise QuiverlabError(
                f"tau_Db: degree-{n} term is not {kind} -- X is not a perfect "
                f"{'projective' if kind == 'projective' else 'injective'} complex "
                f"(got summand {std})")
        verts.append(std[1])
    return verts


def _T(mat, rows, cols, dom):
    """Exact ``cols x rows`` transpose of a ``rows x cols`` matrix, empty-shape safe
    (``lm.transpose`` collapses a 0-row matrix to ``[]``, which breaks a ModuleHom of
    nonzero target dim)."""
    out = lm.zeros(cols, rows, dom)
    for i in range(rows):
        row = mat[i]
        for j in range(cols):
            out[j][i] = row[j]
    return out


def _corner_cover(A, verts, side):
    """``Hom_A((+)_v P_v, A) = (+)_v projective(A^op, v)``, tagged the contravariant
    side -- byte-identical to the object ``corner_transpose`` builds internally, so
    ``dualize`` of it shares the corner-transpose differentials' k-basis."""
    Aop = opposite_algebra(A)
    out_side = _other_side(side)
    mods = [projective(Aop, v) for v in verts]
    if not mods:
        return _zero_module(Aop, side=out_side)
    Q, _ = _direct_sum(mods, name="nu_cover", side=out_side)
    return Q


# --------------------------------------------------------------------------- #
# nu termwise (proj complex -> inj complex) and its inverse (inj -> proj).
# --------------------------------------------------------------------------- #
def _nu_map(d, src_verts, tgt_verts, A, side):
    """``nu(d)`` for a differential ``d: (+)P_{src} -> (+)P_{tgt}`` (rows=``(+)P_tgt``,
    cols=``(+)P_src``), as the matrix ``nu(X_src) -> nu(X_tgt)`` in the ``dualize``
    corner-cover bases. ``d^* = corner_transpose_A(d)`` (contravariant Hom(-, A)), then
    ``nu(d) = D(d^*)`` = its k-dual = the transpose in the injective (dualize) basis."""
    dom = A.domain
    # d^*: (+)Ae_tgt -> (+)Ae_src, rows=(+)Ae_src (=nu(X_src) dim), cols=(+)Ae_tgt.
    dstar = corner_transpose(d, from_verts=tgt_verts, to_verts=src_verts, A=A,
                             side=side)[3]
    rows = len(dstar)
    cols = len(dstar[0]) if (dstar and dstar[0]) else 0
    return _T(dstar, rows, cols, dom)            # D: nu(X_src) -> nu(X_tgt)


def _nu_inverse_map(g, hi_verts, lo_verts, A, side):
    """``nu^{-1}(g)`` for a differential ``g: (+)I_{hi} -> (+)I_{lo}`` between injective
    terms (rows=``(+)I_lo``, cols=``(+)I_hi``), as the matrix ``(+)P_hi -> (+)P_lo``.
    Inverts ``_nu_map``: ``nu(d) = D(corner_transpose_A(d))`` so
    ``nu^{-1}(g) = corner_transpose_{A^op}(D(g))`` -- the double-opposite corner
    transpose recovers the map between the A-projectives."""
    dom = A.domain
    Aop = opposite_algebra(A)
    rows = len(g)
    cols = len(g[0]) if (g and g[0]) else 0
    gt = _T(g, rows, cols, dom)                  # D(g): (+)proj(Aop,lo) -> (+)proj(Aop,hi)
    # corner_transpose over A^op: from=hi (rows/target of gt), to=lo (cols/source of gt)
    return corner_transpose(gt, from_verts=hi_verts, to_verts=lo_verts, A=Aop,
                            side=_other_side(side))[3]


def _apply_nu_perfect(X):
    A, dom, side = X.algebra, X.domain, X.side
    terms, dmats, vlists = {}, {}, {}
    for n in X.degrees():
        if X.term(n).dim == 0:
            continue
        vlists[n] = _term_vertices(X, n, "projective")
        terms[n] = dualize(_corner_cover(A, vlists[n], side))     # nu(X_n) = D Hom(X_n, A)
    for n in X.degrees():
        d = X._dmats.get(n)
        if d and d[0] and n in terms and (n - 1) in terms:
            dmats[n] = _nu_map(d, vlists[n], vlists[n - 1], A, side)
    out = ChainComplex(terms, dmats, check=True)          # d.d=0 self-cert
    # nu's OUTPUT terms are INJECTIVES -- tag correctly (devil's-advocate fix;
    # the old _proj_vertices label here was the masquerade that let a wrong
    # kind ride the fast path)
    out._term_provenance = ("injective", {n: list(vs) for n, vs in vlists.items()})
    return out


def _apply_nu_inverse_perfect(Y):
    A, dom, side = Y.algebra, Y.domain, Y.side
    terms, dmats, vlists = {}, {}, {}
    for n in Y.degrees():
        if Y.term(n).dim == 0:
            continue
        vlists[n] = _term_vertices(Y, n, "injective")
        blocks = [projective(A, v) for v in vlists[n]]            # nu^{-1}(I_v) = P_v
        Q, _ = _direct_sum(blocks, name="nuinv_term", side=side)
        terms[n] = Q
    for n in Y.degrees():
        g = Y._dmats.get(n)
        if g and g[0] and n in terms and (n - 1) in terms:
            dmats[n] = _nu_inverse_map(g, vlists[n], vlists[n - 1], A, side)
    out = ChainComplex(terms, dmats, check=True)
    out._proj_vertices = {n: list(vs) for n, vs in vlists.items()}
    out._term_provenance = ("projective", {n: list(vs) for n, vs in vlists.items()})
    return out


def tau_Db(X):
    """``tau_{D^b}(X) = (nu X)[-1]`` for a certified-perfect complex ``X`` over an
    algebra of finite global dimension (Happel). ``ChainComplex(check=True)`` inside
    ``_apply_nu_perfect`` is the ``d.d = 0`` self-cert."""
    if not X.is_perfect():
        raise QuiverlabError(
            "tau_Db: X must be a perfect complex; resolve it with "
            "projective_model(X, ...) first (P39).")
    _require_finite_gldim(X.algebra)
    return _apply_nu_perfect(X).shift(-1)


def tau_Db_minus(X):
    """``tau^-_{D^b}(X) = (nu^{-1} X)[+1]`` -- the inverse derived AR translate on a
    perfect complex of INJECTIVES (the output shape of :func:`tau_Db`). The net shift
    ``+1`` mirrors ``tau_Db``'s ``-1``; certified by the round-trip quasi-iso
    ``tau_Db_minus(tau_Db(X)) ~ X`` (homology dimension vectors agree degreewise)."""
    _require_finite_gldim(X.algebra)
    return _apply_nu_inverse_perfect(X).shift(1)
