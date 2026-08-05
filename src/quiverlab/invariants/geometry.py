"""Geometry of representations (Plan 49 / C8): orbit dimensions in the
representation variety, Voigt rigidity, and the Kac canonical decomposition.
Right modules; the orbit/rigidity layer is exact over EVERY Domain
(dim End + Ext are exact); the canonical decomposition is a hereditary-Dynkin
notion with a loud refusal off scope. Float-free (int/dict)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules.hom import end_dim


def group_dim(dimvec):
    """dim GL(d) = sum_v d_v^2 -- the base-change group acting on Rep(Q, d)."""
    return sum(int(n) * int(n) for n in dimvec.values())


def representation_variety_dim(algebra, dimvec):
    """dim Rep(Q, d) = sum over arrows a: i -> j of d_i * d_j (the AMBIENT quiver
    representation variety). For kQ/I the module variety mod_A(d) is the closed
    subvariety cut by the relations; this is its ambient space -- stated honestly
    by every codim consumer."""
    total = 0
    for (s, t) in algebra.quiver.arrows.values():
        total += int(dimvec[s]) * int(dimvec[t])
    return total


def orbit_dimension(M):
    """dim of the GL(d)-orbit of M in Rep(Q, d):
        dim O_M = dim GL(d) - dim_k End_A(M) = sum_v d_v^2 - dim End_A(M).
    Aut(M) is Zariski-open in the affine space End_A(M), so
    dim Stab(M) = dim End_A(M). Holds verbatim for kQ/I."""
    return group_dim(M.dimension_vector()) - end_dim(M)


def is_rigid(M):
    """Voigt: M is rigid iff Ext^1_A(M, M) = 0, and then O_M is OPEN in the
    module variety (rigid => open orbit, in general). One Ext call."""
    return M.algebra.ext(M, M, 1) == 0


def rigidity_codim(M):
    """dim Ext^1_A(M, M). HONEST semantics (stated by every caller):
      * HEREDITARY A: EQUALS codim of the orbit closure in Rep(Q, d) -- Voigt,
        since Rep(Q, d) is smooth: codim O_M = dim End(M) - <d,d> = dim Ext^1(M,M).
      * general kQ/I: only an UPPER BOUND -- codim_{mod_A(d)} O_M <= dim Ext^1(M,M)
        (equality iff mod_A(d) is smooth at M). NEVER claim equality off hereditary.
    """
    return M.algebra.ext(M, M, 1)
