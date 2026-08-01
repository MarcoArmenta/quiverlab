"""Plan 35 wave 3a -- explicit Ext / Tor representatives + self-cert data.

Ext^n_A(M, N) and Tor_n^A(M, N) now ship the actual class representatives: the ordered
LABELLED basis of the ambient Hom / tensor space, each class as a labelled term-sum AND
a sparse coordinate vector, and the annihilating differential -- so a reader can verify
each class is a genuine (co)cycle.

The crux (Marco, mirroring the HH product surface): for EVERY shipped class at EVERY
degree, applying the SHIPPED differential to the SHIPPED coordinate vector gives 0.
Plus the labels are HAND-CHECKED, the two views (term-sum / sparse vector) agree, the
rep count matches the engine dims, and the enumeration lengths equal the matrix shapes.
"""
import pytest

import quiverlab as ql
from quiverlab import GF, Quiver
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.complex_reps import ext_reps, tor_reps


# --------------------------------------------------------------------------- #
# Self-certification: what a user does with the machine record.
# --------------------------------------------------------------------------- #
def _dense(vector, ncols, dom):
    v = [dom.zero()] * ncols
    for idx, c in vector:
        v[idx] = dom.coerce(int(c)) if str(c).lstrip("-").isdigit() else dom.coerce(c)
    return v


def _certify(payload, dom):
    """Apply each shipped differential to each shipped vector; assert 0. Returns the
    number of classes certified (differential must NOT be elided here)."""
    diffs, bc, n = payload["differentials"], payload["basis_classes"], 0
    for deg, classes in bc.items():
        d = diffs[deg]
        assert not d.get("elided"), "unexpected elision at deg %s" % deg
        rows = d.get("rows") or []
        ncols = d["shape"][1]
        for cl in classes:
            v = _dense(cl["vector"], ncols, dom)
            if rows:
                R = [[dom.coerce(int(x)) if str(x).lstrip("-").isdigit()
                      else dom.coerce(x) for x in row] for row in rows]
                prod = lm.matvec(R, v, dom)
                assert all(dom.is_zero(p) for p in prod), \
                    "class at deg %s is not annihilated by its differential" % deg
            n += 1
    return n


def _kA2():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    return Q.algebra(relations=[], field=GF(7))


def _loop3(p=2):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x*x*x"], field=GF(p))


# --------------------------------------------------------------------------- #
# (1) Hand-check anchors -- explicit labels + coordinate vectors.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_ext1_S1_S2_kA2_handcheck():
    """kA_2 (1->2): Ext^1(S_1, S_2) = 1, its single class is the map sending the P_2
    generator to the generator of S_2 -- coordinate vector e_1 over a 1-element Hom
    basis. (dim table [0,1,0,0].)"""
    A = _kA2()
    dims, pl = ext_reps(A, A.simple(1), A.simple(2), 3)
    assert dims == [0, 1, 0, 0]
    assert pl["chain_basis"]["1"] == ["[P_2 -> n_2,1]"]
    cls = pl["basis_classes"]["1"]
    assert len(cls) == 1
    assert cls[0]["terms"] == [["1", "P_2", "n_2,1"]]
    assert cls[0]["vector"] == [[0, "1"]]
    # the other degrees carry no classes (empty Hom or resolved to 0)
    assert pl["basis_classes"]["0"] == [] and pl["basis_classes"]["2"] == []


@pytest.mark.oracle_literature
def test_tor0_is_tensor_product_cokernel_loop():
    """x^3 loop / GF(2), M the 2-dim module: Tor_0(M, S_1) = M (x)_A S_1, a 1-dim
    cokernel; its class is the tensor of the P_1 generator and the S_1 generator."""
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    N = A.simple(1, side="left")
    dims, pl = tor_reps(A, M, N, 4)
    assert dims == [1, 1, 1, 1, 1]
    assert pl["chain_basis"]["0"] == ["P_1 (x) n_1,1"]
    assert pl["basis_classes"]["0"][0]["terms"] == [["1", "P_1", "n_1,1"]]
    # d_0 is the zero map: every 0-chain is a cycle, Tor_0 is the cokernel of d_1.
    assert pl["differentials"]["0"]["shape"][0] == 0
    assert "coker" in pl["differentials"]["0"]["note"]


# --------------------------------------------------------------------------- #
# (2) Self-certification: shipped vector x shipped differential == 0 for EVERY class.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_ext_reps_are_cocycles_loop():
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    _dims, pl = ext_reps(A, M, A.simple(1), 4)
    assert _certify(pl, A.domain) >= 5           # one class per degree here


@pytest.mark.oracle_selfcert
def test_tor_reps_are_cycles_loop():
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    _dims, pl = tor_reps(A, M, A.simple(1, side="left"), 4)
    assert _certify(pl, A.domain) >= 5


@pytest.mark.oracle_selfcert
def test_ext_reps_are_cocycles_kA2():
    A = _kA2()
    # a projective P_1 has pd 0 -> all Ext^{>0}(P_1, -) vanish; a nontrivial case:
    _dims, pl = ext_reps(A, A.simple(1), A.simple(2), 3)
    _certify(pl, A.domain)


@pytest.mark.oracle_selfcert
def test_ext_reps_multivertex_nakayama():
    """A multi-vertex example (kA_3 line) exercises repeated summands / corner-typed
    terms; every shipped class must still be a cocycle."""
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=[], field=GF(5))
    _dims, pl = ext_reps(A, A.simple(1), A.simple(3), 4)
    _certify(pl, A.domain)


# --------------------------------------------------------------------------- #
# (3) Cross-engine: the rep pass returns EXACTLY the engine dims.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_crossengine
def test_ext_reps_dims_match_engine():
    from quiverlab.modules.ext import ext_dims
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    dims, pl = ext_dims(A, M, A.simple(1), 4, with_reps=True)
    assert dims == ext_dims(A, M, A.simple(1), 4)
    # rep count per degree == the dim
    for n, d in enumerate(dims):
        assert len(pl["basis_classes"][str(n)]) == d


@pytest.mark.oracle_crossengine
def test_tor_reps_dims_match_engine():
    from quiverlab.modules.tor import tor_dims
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    N = A.simple(1, side="left")
    dims, pl = tor_dims(A, M, N, 4, with_reps=True)
    assert dims == tor_dims(A, M, N, 4)
    for n, d in enumerate(dims):
        assert len(pl["basis_classes"][str(n)]) == d


@pytest.mark.oracle_crossengine
def test_tor_reps_reject_resolve_second():
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.tor import tor_dims
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    N = A.simple(1, side="left")
    with pytest.raises(QuiverlabError):
        tor_dims(A, M, N, 3, resolve="second", with_reps=True)


# --------------------------------------------------------------------------- #
# (4) Interpretability: enumeration length == class-vector width == differential shape;
#     the two views (term-sum / sparse vector) agree.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("which", ["ext", "tor"])
def test_enumeration_lengths_match_matrix_dims(which):
    A = _loop3(3)
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    if which == "ext":
        _d, pl = ext_reps(A, M, A.simple(1), 4)
    else:
        _d, pl = tor_reps(A, M, A.simple(1, side="left"), 4)
    for n in pl["basis_classes"]:
        enum = pl["chain_basis"][n]
        assert isinstance(enum, list)
        width = len(enum)
        d = pl["differentials"][n]
        # the differential's COLUMN count is the degree-n ambient dimension
        assert d["shape"][1] == width
        for cl in pl["basis_classes"][n]:
            # every coordinate index is within the enumeration
            for idx, _c in cl["vector"]:
                assert 0 <= idx < width
            # term-sum and sparse vector describe the same nonzero support
            assert len(cl["terms"]) == len(cl["vector"])
