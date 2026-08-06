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


def _koszul_double_complex(p):
    """The commutative Koszul double complex ``K(x) (x) K(y)`` over
    ``A = k[x,y]/(x^2, xy, y^2)`` (basis ``{1, x, y}``): both differentials are
    multiplication by ``x`` / ``y`` (a genuinely non-regular sequence on the
    maximal ideal, so the Koszul homology is nonzero in every degree). Its total
    complex is a small commutative complex M2 verifies independently."""
    from quiverlab import GF, Quiver
    from quiverlab.specseq.double import DoubleComplex
    dom = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(p)).domain
    mx = [[0, 0, 0], [1, 0, 0], [0, 0, 0]]        # mult by x on {1,x,y}: 1->x
    my = [[0, 0, 0], [0, 0, 0], [1, 0, 0]]        # mult by y: 1->y
    neg_my = [[-v for v in row] for row in my]     # Koszul sign (-1)^p on the p=1 col
    terms = {(0, 0): 3, (1, 0): 3, (0, 1): 3, (1, 1): 3}
    d_h = {(1, 0): mx, (1, 1): mx}
    d_v = {(0, 1): my, (1, 1): neg_my}
    return DoubleComplex(terms, d_h, d_v, dom)


def crosscheck_commutative_ss(p=7, top=4):
    """Plan-42 spectral-sequence crosscheck: our ``E_inf`` totals of the commutative
    Koszul double complex over ``ZZ/p`` vs M2's homology of the SAME total complex
    (the ``Complexes`` package). The compared quantity is the convention-robust
    ``E_inf`` totals (== total-complex homology, by strong convergence) -- M2 1.26's
    ``SpectralSequences`` package is unscriptable (it rides the removed
    ``ChainComplex`` type), so the ``E_2`` page grid is NOT compared; only the
    convergence target is (both systems must agree)."""
    from quiverlab.specseq.pages import SpectralSequence
    dc = _koszul_double_complex(p)
    ss = SpectralSequence(dc.column_filtration())
    Einf = ss.page(ss.convergence.e_infinity_page)
    tot = {}
    for (pp, qq) in Einf.spots:
        tot[pp + qq] = tot.get(pp + qq, 0) + Einf.dim(pp, qq)
    ours = [tot.get(n, 0) for n in range(top + 1)]
    _, dm, _ = dc.total()
    dmats = {n: [[int(x) % p for x in row] for row in dm[n]] for n in dm}
    theirs = scripts.parse_sentinels(
        session.run_script(scripts.commutative_ss_script(p, dmats, top)))
    return CrosscheckReport(what=f"commutative_ss_einf:0..{top}", ours=ours,
                            qpa=theirs, agree=ours == theirs)


_SUBJECTS = {
    "graded_dims": crosscheck_graded_dims,
    "commutative_ext": crosscheck_commutative_ext,
    # commutative_ss takes (p, top) not (algebra, ...) -- called directly by the
    # Plan-42 battery, not through the algebra-first crosscheck() dispatch.
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
