"""Typed M2 comparisons. Reuses the QPA CrosscheckReport container (pure
dataclass; its ``qpa`` field is read as "the external system's value")."""
from __future__ import annotations

from itertools import combinations

from quiverlab.errors import QuiverlabError
from quiverlab.m2 import scripts, session
from quiverlab.qpa.crosscheck import CrosscheckReport


def _our_graded_dims(A, top):
    from quiverlab.modules.koszul import _algebra_graded_matrices
    mats = _algebra_graded_matrices(A, top)
    return [sum(sum(row) for row in mn) for mn in mats]


def crosscheck_graded_dims(A, top=6):
    ours = _our_graded_dims(A, top)
    theirs = scripts.parse_sentinels(
        session.run_script(scripts.graded_dims_script(A, top)))
    return CrosscheckReport(what=f"graded_dims:0..{top}", ours=ours,
                            qpa=theirs, agree=ours == theirs)


def _assert_commutative_presentation(A):
    # A.relations are parsed Relation objects; str() renders each back to source
    # form ("x*y - y*x"), which we normalize by stripping whitespace before the
    # textual commutator check.
    rels = {str(r).replace(" ", "") for r in A.relations}
    gens = sorted(A.quiver.arrows)
    for a, b in combinations(gens, 2):
        if f"{a}*{b}-{b}*{a}" not in rels and f"{b}*{a}-{a}*{b}" not in rels:
            raise QuiverlabError(
                f"not presented as commutative: missing commutator for "
                f"({a},{b}) -- the M2 side would be a different algebra.")


def crosscheck_commutative_ext(A, variables, relations, top=6):
    _assert_commutative_presentation(A)
    ours = A.ext_algebra(top=top).graded_dims_through(top)
    theirs = scripts.parse_sentinels(session.run_script(
        scripts.commutative_ext_script(A.domain.characteristic,
                                       variables, relations, top)))
    return CrosscheckReport(what=f"commutative_ext:0..{top}", ours=ours,
                            qpa=theirs, agree=ours == theirs)


_SUBJECTS = {
    "graded_dims": crosscheck_graded_dims,
    "commutative_ext": crosscheck_commutative_ext,
}


def crosscheck(algebra, what, *args, **kwargs):
    """Dispatch an M2 crosscheck. Valid subjects: graded_dims,
    commutative_ext. Anything else (e.g. hochschild) is refused loudly --
    M2 has no Hochschild theory (honest scope)."""
    if what not in _SUBJECTS:
        raise QuiverlabError(
            f"no M2 oracle for {what!r}; available: "
            f"{sorted(_SUBJECTS)} (M2 has no quiver or Hochschild support).")
    return _SUBJECTS[what](algebra, *args, **kwargs)
