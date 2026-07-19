"""Validation battery (a): CS resolution == Bardzell resolution, exactly, on the
monomial algebra zoo (spec-§6 "CS specializes to Bardzell").

Two INDEPENDENT engines are compared over GF(32003):

  * CS  side: CSResolution(A) -> ChouhySolotarResolution built from
    reduction_system_of(A) (the Plan-03 Groebner reduction system of the quiverlab
    Algebra A); its differential is the CS leading map delta_n plus the
    order-condition correction (resolutions_cs/resolution.py).
  * Bardzell side: BardzellResolution(pres) built from the hand-written
    MonomialPresentation; its differential is Bardzell's alternating big/small
    formula (engine/resolutions_bardzell.py).

The differentials are constructed by disjoint code, so agreement of ranks and HH
dimensions is a genuine cross-check, not a tautology. (The associated-path /
chain enumeration IS shared -- CS's SSequence reuses MonomialPresentation's
associated_paths on the tip algebra's leading words -- so term-count agreement is
partly structural; the rank and HH agreement is the load-bearing part.)
"""
import pytest, numpy as np
from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation
from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.engine.adapter import to_engine
from quiverlab.resolutions_cs.engine_facade import CSResolution
pytest.importorskip("quiverlab.groebner")

# ---------------------------------------------------------------------------
# TRANSCRIPTION DEVIATION (flagged per brief): the brief's third case is
# ("local_radsq", 3, 10). local_radsq(3) is radical-square-zero, so dim C_n grows
# as 4*3^n: at N=10 the differential C_9 -> C_10 is 78732 x 236196 (~1.5e11 int64
# entries, ~148 GB dense), np.linalg.matrix_rank (float SVD) on it is infeasible,
# and the CS differential BUILD alone costs 146 s at n=6 and explodes past it
# (the HH test at depth N builds degree N+1). N=10 therefore fails its own
# transcribed test on memory AND time. Per the brief's minimal-fix clause, N is
# reduced to the largest feasible depth, N=5: the HH test then builds through
# n=6 (~2.5 min, ~200 MB), and CS==Bardzell is verified exactly there
# (term counts, ranks, and HH dims [4, 6, 14, 32, 72, 170] all agree). The other
# two cases and every other value are verbatim. See p04-task-8-report.md.
CASES = [("kx", 3, 12), ("cyclic_nakayama", (5, 3), 10), ("local_radsq", 3, 5)]


def _build(name, param):
    """(core.Algebra over GF(32003), matching MonomialPresentation) per case."""
    from quiverlab import Quiver, GF
    from quiverlab.families import truncated_polynomial
    F = GF(32003)
    if name == "kx":
        # k[x]/(x^param): one vertex, one loop, relation x^param.
        A = truncated_polynomial(param, field=F)
        pres = MonomialPresentation.truncated_polynomial(param)
        return A, pres
    if name == "cyclic_nakayama":
        # kZ_n / rad^ell: n-cycle 0->1->...->(n-1)->0, relations = all length-ell paths.
        n, ell = param
        Q = Quiver(vertices=list(range(n)),
                   arrows={f"a{i}": (i, (i + 1) % n) for i in range(n)})
        rels = ["*".join(f"a{(i + k) % n}" for k in range(ell)) for i in range(n)]
        A = Q.algebra(relations=rels, field=F)
        pres = MonomialPresentation.cyclic_nakayama(n, ell)
        return A, pres
    if name == "local_radsq":
        # g-loop local algebra, rad^2 = 0: one vertex, g loops, every length-2 path
        # is a relation.
        g = param
        Q = Quiver(vertices=[1], arrows={f"x{i}": (1, 1) for i in range(g)})
        rels = [f"x{i}*x{j}" for i in range(g) for j in range(g)]
        A = Q.algebra(relations=rels, field=F)
        pres = MonomialPresentation.local_radsq(g)
        return A, pres
    raise ValueError(f"unknown case {name!r}")


@pytest.mark.parametrize("name,param,N", CASES)
def test_cs_equals_bardzell_terms_and_ranks(name, param, N):
    A, pres = _build(name, param)
    E = to_engine(A.unit_adapted())
    cs, bd = CSResolution(A), BardzellResolution(pres)
    for n in range(N + 1):
        assert len(cs.term_basis(E, n)) == len(bd.term_basis(E, n))
    for n in range(1, N + 1):
        Mc = cs.differential_matrix(E, n, cs.term_basis(E, n),
                                    {g: i for i, g in enumerate(cs.term_basis(E, n - 1))})
        Mb = bd.differential_matrix(E, n, bd.term_basis(E, n),
                                    {g: i for i, g in enumerate(bd.term_basis(E, n - 1))})
        assert np.linalg.matrix_rank(Mc % 32003) == np.linalg.matrix_rank(Mb % 32003)


@pytest.mark.parametrize("name,param,N", CASES)
def test_cs_equals_bardzell_hh(name, param, N):
    A, pres = _build(name, param)
    E = to_engine(A.unit_adapted())
    hb = hochschild_homology_dims(E, N, primes=(32003,), resolution=BardzellResolution(pres))[32003]
    hc = hochschild_homology_dims(E, N, primes=(32003,), resolution=CSResolution(A))[32003]
    assert list(hc) == list(hb)
