"""tau-rigidity + g-vectors + the g-matrix of a support tau-tilting pair
(Plan 45 / C4, Adachi-Iyama-Reiten Compos. Math. 150 (2014)).

A right A-module M is tau-rigid iff Hom_A(M, tau M) = 0 (equivalently Ext^1(M, M) = 0
over a hereditary algebra). The g-vector g^M in K_0(proj A) = Z^{Q_0} is [P_0] - [P_1]
of the minimal projective presentation P_1 -> P_0 -> M -> 0 (read off the resolution's
summand-vertex multisets -- the exact idiom duality.transpose_module uses). g is additive
and g^{P_v} = e_v. For a support tau-tilting pair (M, P) the g-matrix has the module
summands' g-vectors together with -e_v for each killed projective P_v (the 2-term silting
shift P_v[1]); it is unimodular and determines the pair (AIR Thm 5.1 / g-vector
injectivity). Float-free; exact over the Domain."""
from __future__ import annotations

from collections import Counter

from quiverlab.modules.duality import tau
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.resolution import minimal_resolution


def is_tau_rigid(M):
    """True iff ``Hom_A(M, tau M) = 0`` (AIR Def 0.1). tau(projective) = 0 => every
    projective is tau-rigid; the check is one Hom-dimension call today."""
    return hom_dim(M, tau(M)) == 0


def g_vector(M):
    """The g-vector ``g^M`` in ``K_0(proj A) = Z^{Q_0}`` as a vertex-keyed dict:
    ``g^M_v = (mult of P_v in P_0) - (mult of P_v in P_1)`` for the minimal projective
    presentation ``P_1 -> P_0 -> M -> 0``. ``g^{P_v} = e_v``; additive on direct sums."""
    verts = list(M.algebra.quiver.vertices)          # K_0(proj A) basis
    terms, _dmats = minimal_resolution(M, 1)
    c0 = Counter(terms[0].vertices)                  # P_0 summand vertices
    c1 = Counter(terms[1].vertices) if len(terms) > 1 else Counter()   # P_1
    return {v: c0[v] - c1[v] for v in verts}


def g_columns(pair):
    """The ordered list of g-vectors (columns of the g-matrix), one per silting summand:
    the module summands' g-vectors (in summand order) followed by ``-e_v`` for each
    support vertex ``v`` (in vertex order). Each column is a list indexed by the vertex
    order of ``pair.algebra.quiver.vertices``."""
    verts = list(pair.algebra.quiver.vertices)
    idx = {v: i for i, v in enumerate(verts)}
    cols = []
    for Mi in pair.summands:                         # indecomposable module summands
        gv = g_vector(Mi)
        cols.append([gv[v] for v in verts])
    for v in sorted(pair.support, key=lambda w: idx[w]):   # killed projectives -> -e_v
        e = [0] * len(verts)
        e[idx[v]] = -1
        cols.append(e)
    return cols


def g_matrix(pair):
    """The ``n x n`` integer g-matrix of a support tau-tilting pair: rows = vertices,
    columns = :func:`g_columns` (module summands first, then killed projectives ``-e_v``).
    Unimodular (det +-1) by AIR; the SET of columns determines the pair (the canonical
    dedup key is :meth:`SupportTauTiltingPair.g_key`)."""
    cols = g_columns(pair)
    n = len(list(pair.algebra.quiver.vertices))
    return [[cols[c][r] for c in range(len(cols))] for r in range(n)]
