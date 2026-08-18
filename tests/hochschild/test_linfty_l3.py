"""Plan 78 Task 4: the B(A)[1] L-infinity companion.

SPIKE OUTCOME (frozen / ledgered): v1 ships (a) the rad^2=0 dg-Lie certificate
(UNCONDITIONAL) + (b) the induced-l_2 == CS-bracket model-independence THEOREM statement
(cited, not computed -- no B(A) l_2 adapter); the general char-0 Bardzell l_3 is DEFERRED
(l3_bracket status='deferred'), l_4 is NEVER claimed zero. This is a complete honest v1.
"""
import pytest

import quiverlab as ql
from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild.deformations import (
    dg_lie_certificate, l3_bracket, _rad_square_zero, _is_monomial)

pytestmark = pytest.mark.deep


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


@pytest.mark.oracle_selfcert
def test_rad_square_zero_structural_check():
    """rad^2 = 0 iff dim A = |Q_0| + |Q_1| (the cheap structural certificate)."""
    assert _rad_square_zero(_kZ(4, 2)) is True        # dim 8 = 4 + 4
    assert _rad_square_zero(_kZ(3, 2)) is True         # dim 6 = 3 + 3
    assert _rad_square_zero(ql.QuantumCI(0, field=QQ)) is False   # dim 4, has yx


@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_dg_lie_certificate():
    """rad^2=0 monomial => dg-Lie (l_{>=3}=0, True); rad^2!=0 monomial => honest L-infinity
    (False); non-monomial => None (RRB's B(A) scope). UNCONDITIONAL (RRB)."""
    assert dg_lie_certificate(_kZ(4, 2)) is True       # rad^2=0 -> dg-Lie
    assert dg_lie_certificate(_kZ(5, 2)) is True
    assert dg_lie_certificate(ql.QuantumCI(0, field=QQ)) is False   # rad^2=(yx)!=0
    # non-monomial -> None
    comm = ql.QuantumCI(-1, field=QQ)                  # xy - yx is non-monomial
    assert not _is_monomial(comm)
    assert dg_lie_certificate(comm) is None


@pytest.mark.oracle_selfcert
def test_l3_bracket_is_deferred_honestly():
    """The general char-0 Bardzell l_3 is DEFERRED: l3_bracket returns status='deferred'
    with an honest note (never a fabricated value); l_4 is never claimed zero."""
    A = _kZ(4, 2)                                      # monomial, rad^2=0
    rep = l3_bracket(A)
    assert rep.status == "deferred"
    assert rep.monomial is True
    assert rep.dg_lie is True                          # rad^2=0 -> l_{>=3}=0 here
    assert "deferred" in rep.note.lower()
    assert "l_4" in rep.note.lower() or "l4" in rep.note.lower()
    assert "rrb_linfty_bardzell" in rep.references
    # a non-rad^2=0 monomial too
    rep2 = l3_bracket(ql.QuantumCI(0, field=QQ))
    assert rep2.status == "deferred"
    assert rep2.dg_lie is False


@pytest.mark.oracle_selfcert
def test_l3_bracket_refuses_non_monomial_and_char_p():
    """Scope boundary: l3_bracket refuses non-monomial (RRB's L-infinity is monomial) and
    char != 0 loudly."""
    with pytest.raises(QuiverlabError):
        l3_bracket(ql.QuantumCI(-1, field=QQ))         # non-monomial
    with pytest.raises(QuiverlabError):
        l3_bracket(_kZ(4, 2, field=GF(5)))             # char p


def _kA3_rad2(field=QQ):                               # gentle, acyclic, rad^2=0
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


@pytest.mark.oracle_selfcert
def test_deformation_structure_gentle_nilpotent():
    """deformation_structure on a gentle+acyclic+rad^2=0 algebra: nilpotent regime, dg-Lie
    certificate, l3_status='dg_lie', unobstructed. Cheap (HH^2=0 -> no bracket)."""
    from quiverlab.hochschild.deformations import deformation_structure
    D = deformation_structure(_kA3_rad2())
    assert D.characteristic == 0
    assert D.nilpotent_regime is True
    assert D.dg_lie is True
    assert D.l3_status == "dg_lie"
    assert D.unobstructed is True
    assert D.status == "complete"
    assert D.mc_description.startswith("MC = Z^2")


@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
def test_deformation_structure_quantumci0_full():
    """deformation_structure on QuantumCI(0)/QQ: obstructed, dg_lie False, l3 deferred, a
    canonical radical A_alpha populated with its Ext-algebra summary. The Algebra delegate
    agrees. (~18s -- the reference obstruction cost.)"""
    A = ql.QuantumCI(0, field=QQ)
    D = A.deformation_structure()                      # via the delegate
    assert D.hh2_dim == 3 and D.hh3_dim == 5
    assert D.unobstructed is False
    assert D.obstruction_witness is not None
    assert D.nilpotent_regime is None                  # non-gentle
    assert D.dg_lie is False                            # rad^2 != 0
    assert D.l3_status == "deferred"
    assert D.a_alpha is not None
    assert D.a_alpha["flat"] is True
    assert D.a_alpha["relations"]
    assert D.ext_algebra_summary is not None
    assert D.ext_algebra_summary.startswith("E(A_alpha):")
    assert D.base_change_note is not None
    assert "model-independence" in D.note.lower()


@pytest.mark.oracle_selfcert
def test_deformation_structure_charp_field_general():
    """Over char p deformation_structure returns the field-general HH block + a char0_note,
    no expensive bracket, unobstructed=None (char-0 notion)."""
    from quiverlab.hochschild.deformations import deformation_structure
    D = deformation_structure(ql.QuantumCI(0, field=GF(5)))
    assert D.characteristic == 5
    assert D.hh2_dim == 3 and D.hh3_dim == 5
    assert D.unobstructed is None
    assert D.dg_lie is None
    assert D.l3_status == "n/a (char p)"
    assert D.char0_note is not None
    assert D.status == "complete"
