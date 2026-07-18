"""Plan 02 exit criteria, executable.

Three things the engine port must demonstrate end-to-end:

  1. The deep monomial resolution reaches depths the bar oracle cannot, and the
     minimal (Bardzell) bimodule resolution certifies structural facts -- a finite
     global dimension shows up as vanishing generators past that dimension -- at a
     degree (30) where the normalized bar complex is combinatorially hopeless.
  2. On single-cycle monomial algebras the engine computes genuinely NONZERO deep
     Hochschild homology to depth 30, where the bar complex's term C_30 = a*(a-1)^30
     is astronomically out of reach (~8e14 cells already for k[x]/(x^4)).
  3. Where both the fast GF(p) engine and the pure bar oracle can run, they agree
     exactly through the public API.

The whole-suite pure-path check (QUIVERLAB_NO_NUMBA=1 sweep) is run by the task
checklist, not from inside a test.
"""
from quiverlab import GF, Quiver
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import (
    hochschild_homology_dims,
    truncated_polynomial,
)
from quiverlab.engine.resolutions_bardzell import (
    BardzellResolution,
    MonomialPresentation,
)

# hanlab __init__ alias, reproduced locally:
homology_dims = hochschild_homology_dims
CHAR0_PROXY = 32003  # large prime standing in for characteristic 0


def test_deep_monomial_homology_beyond_bar_reach():
    """kA_3 with rad^2 = 0 (arrows a: 1->2, b: 2->3, relation a*b): a monomial
    algebra of finite global dimension. The minimal (Bardzell) bimodule resolution
    walked to degree 30 certifies gl.dim <= 2 as VANISHING generators past degree 2,
    and the induced deep Hochschild homology (HH_0 = 3, HH_n = 0 for n >= 1 -- the
    quiver is acyclic, so no oriented cycle closes any associated path) is computed
    to depth 30 through the real engine algebra. The bar complex could never reach
    degree 30 on a three-vertex algebra; the engine does it instantly, with a
    certificate.
    """
    # Public route to the engine algebra: build kA_3/(a*b) over GF(5), lower it.
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=["a*b"], field=GF(5))
    E = to_engine(A.unit_adapted())

    # The matching monomial presentation for the Bardzell backend. Bank convention:
    # arrows = (arrow_id, source, target); relations = path-tuples of arrow_ids.
    pres = MonomialPresentation(
        vertices=[1, 2, 3],
        arrows=[("a", 1, 2), ("b", 2, 3)],
        relations=[("a", "b")],
    )
    res = BardzellResolution(pres)

    # (a) The finite-gl.dim certificate: the bimodule-resolution generators are the
    # associated paths AP^n. |AP^1| = 2 (arrows), |AP^2| = 1 (the single relation),
    # and AP^n = {} for every 3 <= n <= 30 -- the minimal resolution stops, walked to
    # depth 30 by the Anick n-chain recursion (something the bar complex cannot do).
    bound = res._path_bound(30)
    ap = {n: len(pres.associated_paths(n, bound)) for n in range(1, 31)}
    assert ap[1] == 2 and ap[2] == 1
    assert all(ap[n] == 0 for n in range(3, 31))

    # (b) The tensored-down complex terms C_n = term_basis(n) walked to degree 30:
    # C_0 > 0, and the deep term vanishes. (term_basis is algebra-independent for the
    # Bardzell backend -- it reads the presentation -- but we thread the real engine
    # algebra to keep the call end-to-end.)
    term_dims = [len(res.term_basis(E, n)) for n in range(31)]
    assert term_dims[0] > 0
    assert term_dims[-1] == 0

    # (c) The actual deep Hochschild homology through the engine, to depth 30.
    dims = homology_dims(E, 30, resolution=res)[5]
    assert len(dims) == 31
    assert dims == [3] + [0] * 30


def test_deep_truncated_homology_nonzero_beyond_bar_reach():
    """Single-cycle monomial algebras k[x]/(x^a): the Bardzell resolution has bounded
    term dim (a per degree) and reaches depth 30 instantly with genuinely NONZERO
    homology, while the bar complex's C_30 = a*(a-1)^30 (~8.2e14 for a = 4) is
    hopeless. Known closed forms, cross-checked against the bar oracle at shallow
    depth by test_bardzell_resolution, are asserted here at depth 30.
    """
    # k[x]/(x^2): HH_n = 2 for ALL n in characteristic 2 (a torsion pathology), and
    # [2, 1, 1, ...] in characteristic 0. Deep and nonzero throughout.
    A2 = truncated_polynomial(2)
    res2 = BardzellResolution(MonomialPresentation.truncated_polynomial(2))
    d2 = homology_dims(A2, 30, resolution=res2)
    assert d2[2] == [2] * 31                       # char 2: every HH_n = 2, nonzero
    assert all(v > 0 for v in d2[2])
    assert d2[CHAR0_PROXY] == [2] + [1] * 30       # char 0 proxy

    # k[x]/(x^4): [4] + [3]*30 in characteristic 0 (the depth-40 unlock, at depth 30).
    A4 = truncated_polynomial(4)
    res4 = BardzellResolution(MonomialPresentation.truncated_polynomial(4))
    d4 = homology_dims(A4, 30, resolution=res4)
    assert d4[CHAR0_PROXY] == [4] + [3] * 30
    assert all(v > 0 for v in d4[CHAR0_PROXY])


def test_engine_agrees_with_bar_on_truncated_x4():
    """k[x]/(x^4) over GF(5): the fast GF(p) engine and the pure bar oracle agree
    exactly through the public API, over every degree the (exponential) bar path can
    still reach (depth 5)."""
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^4"], field=GF(5))
    fast = A.hochschild_homology(5, engine="fast").dims
    bar = A.hochschild_homology(5, engine="bar").dims
    assert fast == bar == [4, 3, 3, 3, 3, 3]
