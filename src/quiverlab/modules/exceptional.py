"""Classical exceptional sequences of a hereditary algebra (Plan 65 / R28).

For a **hereditary** ``A = kQ`` (``Q`` acyclic, ``I = 0``, so ``gl.dim <= 1``):

- an **exceptional module** ``M`` has ``End_A(M)`` a division ring and
  ``Ext^{>=1}_A(M, M) = 0``. Over a hereditary algebra only ``Ext^1`` can be nonzero.
  We implement the **brick** criterion ``end_dim(M) = 1 and Ext^1(M, M) = 0`` (M3): a
  1-dimensional endomorphism ring is ``k`` -- a field, hence a division ring -- so this is
  never a FALSE positive; over a non-algebraically-closed field a module can be exceptional
  with a *larger* division-ring ``End`` (``end_dim > 1``), which the brick test returns
  ``False`` for (a scope-limited false negative that does NOT occur on the Dynkin/QQ
  battery -- every Dynkin indecomposable is a rigid brick). Over algebraically closed ``k``,
  brick <=> exceptional exactly.
- an **exceptional sequence** ``(E_1, ..., E_r)`` is a tuple of exceptional modules with, for
  every ``i < j``, ``Hom_A(E_j, E_i) = 0`` AND ``Ext^1_A(E_j, E_i) = 0`` (no maps or
  extensions BACKWARD, from a later term to an earlier one -- the Crawley-Boevey / Ringel
  convention; the opposite convention merely reverses each sequence, same count). Complete
  iff ``r = n = #simples``.
- the **braid group** ``B_n`` acts on complete exceptional sequences by mutation ``sigma_i``
  (Task B2), transitively (Crawley-Boevey Ottawa 1992; Ringel Contemp. Math. 171 (1994)).
- for **Dynkin** ``A`` the set of complete exceptional sequences is finite, of size
  ``#CES = n! * h^n / |W|`` (Obaid et al.); ``A_n = (n+1)^{n-1}``, ``D_4 = 162``.

Scope (contractual, loud typed refusals): **hereditary only** (non-hereditary refuses);
enumeration is **representation-finite (Dynkin) only** (rep-infinite hereditary has an
infinite braid orbit -> honest ``status`` cap). Over ``GF(p)`` the shipped
``decompose``/``is_isomorphic`` char caveat (char 0 or char > dim) applies. Exact only --
c-matrices are integer matrices, counts ``int``, verdicts ``bool``/``str``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.modules.ext import ext
from quiverlab.modules.hom import end_dim, hom_dim


def _require_hereditary(A):
    """Refuse loudly unless ``A`` is a hereditary path algebra ``kQ`` (``I = 0``). Delegates
    to :func:`invariants.recognizers.is_hereditary`, which also refuses a structure-constant-
    only algebra (no quiver)."""
    from quiverlab.invariants.recognizers import is_hereditary
    if not is_hereditary(A):
        raise QuiverlabError(
            "classical exceptional sequences require a hereditary algebra kQ (I = 0); "
            "this algebra has relations (gl.dim >= 2) or no quiver presentation",
            hint="the tau-exceptional surface (tau_exceptional_sequences) covers "
                 "non-hereditary tau-tilting-finite algebras")


# --------------------------------------------------------------------------- #
# recognizers
# --------------------------------------------------------------------------- #
def is_exceptional_module(A, M) -> bool:
    """True iff ``M`` is a brick with no self-extension (the BRICK criterion, M3):
    ``end_dim(M) == 1 and Ext^1_A(M, M) == 0``. Over a hereditary algebra ``gl.dim <= 1``, so
    higher self-Ext vanish automatically and ``Ext^1`` is the only obstruction. Coincides
    with "exceptional" over algebraically closed ``k`` and on the Dynkin/QQ battery."""
    _require_hereditary(A)
    return end_dim(M) == 1 and ext(A, M, M, 1) == 0


def is_exceptional_sequence(A, seq) -> bool:
    """True iff ``seq = (E_1, ..., E_r)`` is an exceptional sequence: each ``E_i`` is an
    exceptional module and, for every ``i < j``, ``Hom_A(E_j, E_i) = 0`` and
    ``Ext^1_A(E_j, E_i) = 0`` (no backward maps/extensions -- the fixed convention)."""
    _require_hereditary(A)
    seq = list(seq)
    for E in seq:
        if not (end_dim(E) == 1 and ext(A, E, E, 1) == 0):
            return False
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if hom_dim(seq[j], seq[i]) != 0:
                return False
            if ext(A, seq[j], seq[i], 1) != 0:
                return False
    return True


def c_matrix(A, seq):
    """The c-matrix of an exceptional sequence: rows = the dimension vectors of the terms in
    vertex order (integers). For a MODULE (unshifted) sequence every row is the non-negative
    dim-vector of ``E_i``; this is the hereditary exceptional-sequence reading of the P45
    c-vectors (wall normals / brick dim-vectors). Hereditary-only."""
    _require_hereditary(A)
    verts = list(A.quiver.vertices)
    rows = []
    for E in seq:
        dv = E.dimension_vector()
        rows.append([int(dv[v]) for v in verts])
    return rows
