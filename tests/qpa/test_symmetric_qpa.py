"""Plan 29 Part 0 — QPA crosscheck of ``is_symmetric`` / ``is_weakly_symmetric``.

quiverlab's symmetry verdicts vs GAP/QPA's ``IsSymmetricAlgebra`` /
``IsWeaklySymmetricAlgebra`` (``crosscheck_symmetric``). Reproduces the fixed
live bug: the multi-vertex symmetric Nakayama (Brauer star) algebras kZ_n/J^L
with n | (L-1) were wrongly reported non-symmetric by the former GF(p) engine
shortcut; QPA's ``IsSymmetricAlgebra`` returns ``true``.

qpa-marked: skips locally when GAP is absent, mandatory under
``QUIVERLAB_REQUIRE_QPA=1`` in CI.

Plan 31 update: ``TrivialExtension(A)`` is now a genuine ``kQ/rels``-presented
Algebra, so it CAN be fed to QPA's ``IsSymmetricAlgebra`` route; both that
crosscheck and QPA's native ``TrivialExtensionOfQuiverAlgebra`` construction
oracle live in ``tests/qpa/test_trivial_extension_qpa.py``. This file covers the
Nakayama / Brauer-star side (the original fixed-bug reproduction).
"""
import pytest

from quiverlab import GF
from quiverlab.families import NakayamaAlgebra, QuantumCI
from quiverlab.families.zoo import build_from_record, load_catalog
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _rec(name):
    return next(r for r in load_catalog() if r.get("name") == name)


# -- the headline: multi-vertex symmetric Brauer stars, QPA-confirmed ----------
@pytest.mark.parametrize("n,ell", [(2, 3), (3, 4), (4, 5)])
@pytest.mark.parametrize("field", [GF(3), QQ])
def test_brauer_stars_symmetric_agree_with_qpa(n, ell, field):
    # n | (ell - 1) => symmetric AND weakly symmetric; QPA IsSymmetricAlgebra=true.
    A = NakayamaAlgebra(n=n, l=ell, cyclic=True, field=field)
    rep = A.crosscheck("symmetric").assert_agree()
    assert rep.ours == [True, True]


# -- negatives: n does NOT divide (ell - 1) ------------------------------------
@pytest.mark.parametrize("n,ell", [(3, 3), (4, 7)])
@pytest.mark.parametrize("field", [GF(3), QQ])
def test_nakayama_non_dividing_agree_with_qpa(n, ell, field):
    A = NakayamaAlgebra(n=n, l=ell, cyclic=True, field=field)
    rep = A.crosscheck("symmetric").assert_agree()
    assert rep.ours == [False, False]


# -- weakly symmetric is STRICTLY weaker than symmetric (QPA agrees) -----------
@pytest.mark.parametrize("q", [1, 2])
def test_quantum_ci_weakly_but_not_symmetric_agree_with_qpa(q):
    # single vertex => identity Nakayama permutation => weakly symmetric, but nu
    # is not inner => NOT symmetric. QPA: IsWeaklySymmetricAlgebra=true,
    # IsSymmetricAlgebra=false.
    A = QuantumCI(q, field=QQ)
    rep = A.crosscheck("symmetric").assert_agree()
    assert rep.ours == [False, True]


def test_commutative_ci_symmetric_agree_with_qpa():
    QuantumCI(-1, field=QQ).crosscheck("symmetric").assert_agree()


# -- zoo members (multi-vertex, incl. the Plan-18 line_abc_cde) ----------------
@pytest.mark.parametrize("name", ["cn_3_2", "line_abc_cde"])
@pytest.mark.parametrize("field", [GF(3), QQ])
def test_zoo_members_symmetry_agree_with_qpa(name, field):
    A = build_from_record(_rec(name), field=field)
    A.crosscheck("symmetric").assert_agree()
