"""Plan 31 — QPA crosscheck of the presented TrivialExtension (design decision D6).

Two variants, both against the LIVE QPA backend:

* Variant A — feed OUR presented ``T(A)`` (a genuine ``kQ_T/I_T``) through the
  existing ``crosscheck_symmetric`` route (``A.crosscheck("symmetric")``): QPA
  rebuilds the presentation and confirms ``IsSymmetricAlgebra`` /
  ``IsWeaklySymmetricAlgebra`` — both ``true`` — matching our ``[True, True]``.
  (Impossible before Plan 31: ``T(A)`` was presentation-free and could not be fed
  to QPA's ``kQ/rels`` route.)

* Variant B — the NEW ``trivial_extension`` crosscheck (``A.crosscheck(
  "trivial_extension")``): QPA builds ``TrivialExtensionOfQuiverAlgebra(A)``
  NATIVELY and compares dim, arrow count, and the three self-injectivity
  predicates (``IsSymmetricAlgebra`` / ``IsWeaklySymmetricAlgebra`` /
  ``IsSelfinjectiveAlgebra``) against our presented build. QPA's dual-arrow labels
  differ (``te_a1_i_j``), so counts — not names — are compared.

Live-verified pins (2026-07-26): kA_2 → (dim 6, 2 arrows), kA_3 → (12, 3),
2-Kronecker → (8, 4); all three predicates ``true`` in every case.

qpa-marked: skips locally when GAP is absent, mandatory under
``QUIVERLAB_REQUIRE_QPA=1`` in CI. QPA fields are limited to QQ / prime GF(p).
"""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.families import TrivialExtension
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _kron(field):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=field)


def _comm_square(field):
    return Quiver([1, 2, 3, 4],
                  {"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4)}
                  ).algebra(relations=["a*c - b*d"], field=field)


# -- Variant A: our presented T(A) through QPA's IsSymmetricAlgebra route -------
@pytest.mark.parametrize("build", [
    lambda: linear_path_algebra(2, field=QQ),        # kA_2
    lambda: _comm_square(QQ),                        # commutative square
])
def test_presented_TA_symmetric_agrees_with_qpa(build):
    T = TrivialExtension(build())
    assert T.quiver is not None                       # presented -> feedable to QPA
    rep = T.crosscheck("symmetric").assert_agree()
    assert rep.ours == [True, True]                   # symmetric AND weakly symmetric


# -- Variant B: QPA's native TrivialExtensionOfQuiverAlgebra construction oracle -
@pytest.mark.parametrize("build,dim,arrows", [
    (lambda: linear_path_algebra(2, field=QQ), 6, 2),    # kA_2 -> T dim 6, 2 arrows
    (lambda: linear_path_algebra(3, field=QQ), 12, 3),   # kA_3 -> T dim 12, 3 arrows
    (lambda: _kron(QQ), 8, 4),                           # 2-Kronecker -> T dim 8, 4 arrows
])
def test_native_trivial_extension_construction_agrees_with_qpa(build, dim, arrows):
    base = build()
    # our presented build pinned to the QPA-live-verified dim + total arrow count ...
    T = TrivialExtension(base)
    assert T.dim == dim
    assert len(T.quiver.arrows) == arrows
    # ... and QPA's native TrivialExtensionOfQuiverAlgebra must agree on dim, arrow
    # count, and IsSymmetric / IsWeaklySymmetric / IsSelfinjective (all true).
    base.crosscheck("trivial_extension").assert_agree()
