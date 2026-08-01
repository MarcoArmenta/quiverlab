import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import (
    cs_cohomology_dims, cs_homology_dims, cs_hh_basis)
pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_literature]


@pytest.mark.oracle_crossengine
@pytest.mark.oracle_selfcert
def test_cs_hh_basis_matches_dims_multivertex_zero_codomain():
    """Regression (Plan 35 wave 3d): ``cs_hh_basis`` must return EXACTLY
    ``cs_(co)homology_dims`` many representatives at every degree -- INCLUDING a degree
    whose outgoing differential lands in a 0-dimensional codomain (the whole C_n is then
    the kernel). Before the fix ``nullspace`` was handed a 0-row matrix and could not
    recover its width, silently returning [] -- so ``cs_hh_basis`` under-reported
    representatives whenever the top (co)chain differential vanished into a 0-dim space,
    e.g. every finite-global-dimension multi-vertex algebra at its top HH degree.

    ``kZ_3 / J^2`` (the 3-cycle with radical square zero) is the minimal witness: HH^1
    has dimension 1 while C^2 = 0, so the old code returned 0 cocycles there."""
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=["a*b", "b*c", "c*a"], field=GF(5))
    top = 5
    cdims = cs_cohomology_dims(A, top).dims
    hdims = cs_homology_dims(A, top).dims
    assert cdims[1] == 1                          # the witness: HH^1 = 1 with C^2 = 0
    for side, dims in (("coh", cdims), ("hom", hdims)):
        for n in range(top + 1):
            reps = cs_hh_basis(A, n, side)
            assert len(reps) == dims[n], (side, n, len(reps), dims[n])
            # self-cert: each rep annihilates its outgoing differential (b_0 excepted)
            if n == 0 and side == "hom":
                continue
            from quiverlab.resolutions_cs.build import reduction_system_of
            from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
            res = ChouhySolotarResolution(A, reduction_system_of(A), max_degree=n + 1)
            D = res.matrix(n, side)                # delta^n / b_n
            dom = A.domain
            for v in reps:
                for row in D:
                    s = 0
                    for j in range(len(v)):
                        s = dom.add(s, dom.mul(dom.coerce(row[j]), dom.coerce(v[j])))
                    assert dom.is_zero(s), (side, n)


def _A(field=CC, rels=("x*x",), arrows=None, verts=(1,)):
    return Quiver(list(verts), arrows or {"x": (1, 1)}).algebra(relations=list(rels), field=field)


def test_kx2_dims_char0_and_char2():
    assert cs_cohomology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    assert cs_homology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    assert cs_cohomology_dims(_A(field=GF(2)), 5).dims == [2, 2, 2, 2, 2, 2]


def test_square_dims():
    A = _A(rels=["a*b - c*d"], verts=(1, 2, 3, 4),
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    assert cs_cohomology_dims(A, 4).dims == [1, 0, 0, 0, 0]
    assert cs_homology_dims(A, 4).dims == [4, 0, 0, 0, 0]


def test_qci_homology_matches_bank_vector():
    A = _A(rels=["x*x", "y*y", "y*x - 2*x*y"], arrows={"x": (1, 1), "y": (1, 1)})
    assert cs_homology_dims(A, 12).dims == [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


@pytest.mark.oracle_selfcert
def test_hh_dims_unchanged_and_gate_green_after_go_loud(kx2_rs, qci_rs):
    """Regression for GO-LOUD (open-item #1): making SSequence.S(n>max_degree) raise
    must NOT change any HH dimension, and assert_dd_zero(upto=top+1, side="coh") must
    still pass through the public build (its internal S(max_degree+1) read is now the
    explicit empty cochain space, not the silent out-of-range []). Byte-identical dims."""
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    # k[x]/(x^2): known HH^* dims through top=6.
    assert cs_cohomology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    res = ChouhySolotarResolution(_A(), kx2_rs, max_degree=7)   # max_degree = top + 1
    res.assert_dd_zero(upto=7, side="coh")                       # reads S(8)=S(max_degree+1): explicit empty
    # quantum complete intersection over CC: known HH^* dims through top=6.
    Aq = _A(rels=["x*x", "y*y", "y*x - 2*x*y"], arrows={"x": (1, 1), "y": (1, 1)})
    assert cs_cohomology_dims(Aq, 6).dims == [2, 2, 1, 0, 0, 0, 0]
    resq = ChouhySolotarResolution(Aq, qci_rs(xi="2"), max_degree=7)
    resq.assert_dd_zero(upto=7, side="coh")


def test_engine_facade_is_resolution_protocol():
    from quiverlab.resolutions_cs.engine_facade import CSResolution
    from quiverlab.engine.resolutions import Resolution
    from quiverlab.engine.adapter import to_engine
    Ap = _A(field=GF(32003), rels=["a*b - c*d"], verts=(1, 2, 3, 4),
            arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    R = CSResolution(Ap)
    assert isinstance(R, Resolution)
    E = to_engine(Ap.unit_adapted())
    b2 = R.term_basis(E, 2)
    assert R.differential_matrix(E, 2, b2, {g: i for i, g in enumerate(R.term_basis(E, 1))}).shape[1] == len(b2)
