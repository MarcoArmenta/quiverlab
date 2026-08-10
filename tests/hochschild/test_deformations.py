"""Plan 78 Task 1: infinitesimal deformations (HH^2) + the primary obstruction
alpha |-> [alpha,alpha] in HH^3 (the targeted (2,2) CS Gerstenhaber self-bracket).

Pins are BASIS-INDEPENDENT only (Minor-2 / Plan-35 basis rule): we assert
`unobstructed`, the nonzero obstruction CLASS in HH^3, dims -- never the witness
string, the raw constants, or the diagonal-vanishing pattern (that pattern is
basis-DEPENDENT; live on this engine QuantumCI(0)'s diagonal happens to vanish, but
on another basis [a_1,a_1] != 0 -- so it is never pinned).
"""
import pytest

import quiverlab as ql
from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.hochschild.deformations import (
    infinitesimal_deformations, obstruction_map, _obstruction_data,
    self_bracket_class, _is_zero_vec)

pytestmark = pytest.mark.deep


def _kx(a, field=QQ):
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x^%d" % a], field=field)


def _kZ(n, ell, field=QQ):
    names = "abcdefghijklmnopqrstuvwxyz"
    verts = list(range(1, n + 1))
    arrows = {names[i]: (verts[i], verts[(i + 1) % n]) for i in range(n)}
    rels = []
    for start in range(n):
        word, cur = [], start
        for _ in range(ell):
            word.append(names[cur])
            cur = (cur + 1) % n
        rels.append("*".join(word))
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


# One expensive QuantumCI(0)/QQ obstruction (~18s) reused across assertions.
@pytest.fixture(scope="module")
def qci0_qq():
    return _obstruction_data(ql.QuantumCI(0, field=QQ))


# ---------------------------------------------------------------------------
@pytest.mark.oracle_literature
@pytest.mark.oracle_crossengine
def test_infinitesimal_kx2_tangent_is_one():
    """Pin 1: HH^2(k[x]/x^2)/QQ = 1 is the 1-dim deformation tangent."""
    A = _kx(2)
    info = infinitesimal_deformations(A)
    assert info["hh2_dim"] == 1
    assert info["basis"] == "cs/QQ"
    assert info["characteristic"] == 0
    assert info["char0_note"] is None
    # crossengine: == A.hochschild_cohomology(2)
    assert info["hh2_dim"] == A.hochschild_cohomology(2, engine="cs")[2]


@pytest.mark.oracle_literature
@pytest.mark.oracle_crossengine
def test_obstruction_quantumci0_is_real(qci0_qq):
    """Pin 2: QuantumCI(0)/QQ has a GENUINE second-order obstruction. Basis-independent
    pins ONLY: hh2/hh3 dims, unobstructed is False, and the witness produces a NONZERO
    class in HH^3 (degree 3, of the recorded HH^3 dim)."""
    A = ql.QuantumCI(0, field=QQ)
    data = qci0_qq
    assert data["hh2_dim"] == 3
    assert data["hh3_dim"] == 5
    # crossengine: dims == hochschild_cohomology
    table = A.hochschild_cohomology(3, engine="cs")
    assert (data["hh2_dim"], data["hh3_dim"]) == (table[2], table[3])

    ob = obstruction_map(A)
    assert ob.unobstructed is False
    assert ob.hh2_dim == 3 and ob.hh3_dim == 5
    assert ob.basis == "cs/QQ"
    assert ob.witness is not None
    assert "rrb_linfty_bardzell" in ob.references


@pytest.mark.oracle_selfcert
def test_obstruction_witness_class_is_nonzero_in_hh3(qci0_qq):
    """The witness direction's self-bracket [alpha,alpha] is a NONZERO class in HH^3
    (the basis-independent reason it is obstructed). Some pairwise bracket class is
    nonzero, and unobstructedness would require ALL of them to vanish."""
    from quiverlab.hochschild.deformations import _verdict
    data = qci0_qq
    dom = data["dom"]
    unobstructed, witness, coords = _verdict(data)
    assert unobstructed is False
    assert coords is not None
    cls = self_bracket_class(data, coords)
    assert not _is_zero_vec(cls, dom)          # nonzero class in HH^3
    assert len(cls) == data["hh3_dim"] == 5    # lives in degree 3, recorded dim
    # some pairwise bracket class is nonzero (whole quadratic form scanned)
    B, n = data["B"], data["hh2_dim"]
    assert any(not _is_zero_vec(B[i][j], dom) for i in range(n) for j in range(n))


@pytest.mark.oracle_selfcert
def test_obstruction_bracket_is_a_cocycle(qci0_qq):
    """[alpha_i, alpha_j] is a genuine HH^3 cocycle: delta^3 applied to the raw native
    bracket cochain vanishes (the NW theorem, verified on the load-bearing nonzero
    entry). Certifies the returned obstruction chain descends to HH^3."""
    from quiverlab.resolutions_cs.bracket import native_bracket
    data = qci0_qq
    res, dom, coh2, B, n = (data["res"], data["dom"], data["coh2"],
                            data["B"], data["hh2_dim"])
    # find a nonzero pairwise bracket (the load-bearing one)
    pair = next((i, j) for i in range(n) for j in range(n)
                if not _is_zero_vec(B[i][j], dom))
    i, j = pair
    vec = native_bracket(res, coh2[i], 2, coh2[j], 2)
    M = res.matrix(3, "coh")                   # delta^3 : C^3 -> C^4
    dvec = [sum((dom.mul(row[c], vec[c]) for c in range(len(vec))), dom.zero())
            for row in M]
    assert _is_zero_vec(dvec, dom)


@pytest.mark.oracle_literature
def test_gentle_nilpotent_examples_are_unobstructed():
    """Pin 3: small gentle / nilpotent examples are UNOBSTRUCTED (MRRS Thm 5.4 in
    action). kZ4/J^2 (rad^2=0, HH^3=0) and kZ3/J^3 (HH=[1,1,1,1], all self-brackets 0)."""
    for A in (_kZ(4, 2), _kZ(3, 3)):
        ob = obstruction_map(A)
        assert ob.unobstructed is True
        assert ob.witness is None
        assert ob.characteristic == 0


@pytest.mark.oracle_selfcert
def test_obstruction_charp_returns_field_general_block():
    """Char p: obstruction_map(QuantumCI(0)/GF(5)) returns the field-general HH^2/[a,a]
    block with a char0_note (the deformation INTERPRETATION is char-0 gated), not a crash.
    The bracket math is field-general so the obstruction is still detected."""
    A = ql.QuantumCI(0, field=GF(5))
    ob = obstruction_map(A)
    assert ob.characteristic == 5
    assert ob.char0_note is not None and "characteristic-0" in ob.char0_note
    assert ob.hh2_dim == 3 and ob.hh3_dim == 5
    assert ob.basis == "cs/GF(5)"
    assert ob.unobstructed is False


@pytest.mark.oracle_selfcert
def test_presentation_less_refuses_loudly():
    """A structure-constant algebra with no quiver presentation is refused loudly
    (the CS engine needs a presentation)."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.core.algebra import Algebra
    A = ql.QuantumCI(0, field=QQ)
    # a presentation-less structure-constant algebra: same T/unit, no quiver
    B = Algebra(A.domain, A.T, A.unit)
    assert B.quiver is None
    with pytest.raises(QuiverlabError):
        obstruction_map(B)
