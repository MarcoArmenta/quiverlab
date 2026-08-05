"""QPA (GAP) as the external oracle for the Plan-37 C1 categorical glue.

PROBE (documented by ``test_qpa_exposes_hom_glue_surface``): QPA exposes
``HomOverAlgebra(M, N)`` (a k-BASIS of Hom; ``Length`` = dim Hom) AND the
module-homomorphism kernel/image/cokernel surface of QPA Ch. 7 --
``KernelInclusion(f)`` / ``ImageInclusion(f)`` / ``CoKernelProjection(f)`` (there is
no ``KernelOfModuleHomomorphism``; the generic-mapping ``KernelInclusion`` is the
name QPA actually binds). So this battery covers BOTH dim Hom AND, for a canonical
hom, the kernel/image/cokernel dimension vectors -- not merely Hom dims.

For every ordered pair (M, N) of projectives/simples of kA_2, kA_3/(ab) and the
Plan-18 line algebra ``line_abc_cde``, ``crosscheck_hom_glue``:
  * compares ``len(A.hom_basis(M, N))`` with ``Length(HomOverAlgebra(M, N))``; and
  * when ``dim Hom = 1`` (the nonzero hom is unique up to a scalar, so its
    kernel/image/cokernel are well-defined invariants both bases agree on) compares
    the kernel/image/cokernel dimension vectors of that canonical hom against QPA.

Over GF(7) (Hom is field-generic; char 7 keeps every module small). qpa-marked:
skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.qpa import session
from quiverlab.qpa.crosscheck import crosscheck_hom_glue

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

_F = GF(7)


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=_F)


def _kA3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=_F)


def _line_abc_cde():
    # Plan-18 multi-vertex record: 1->2->...->6 with relations a*b*c and c*d*e.
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5), "e": (5, 6)})
    return Q.algebra(relations=["a*b*c", "c*d*e"], field=_F)


def _mods(A):
    verts = list(A.quiver.vertices)
    return [A.projective(v) for v in verts] + [A.simple(v) for v in verts]


def test_qpa_exposes_hom_glue_surface():
    """Document the probe the battery relies on: QPA binds HomOverAlgebra and the
    module-hom kernel/image/cokernel names (and NOT KernelOfModuleHomomorphism)."""
    lg = session.libgap_handle()
    for name in ("HomOverAlgebra", "KernelInclusion", "ImageInclusion",
                 "CoKernelProjection"):
        assert bool(lg.IsBoundGlobal(name)), name
    assert not bool(lg.IsBoundGlobal("KernelOfModuleHomomorphism"))


@pytest.mark.parametrize("build,min_pairs", [
    (_kA2, 16), (_kA3, 36), (_line_abc_cde, 144),
])
def test_hom_glue_battery(build, min_pairs):
    A = build()
    mods = _mods(A)
    seen_hom = 0
    seen_kernel = 0
    for M in mods:
        for N in mods:
            rep = crosscheck_hom_glue(A, M, N)
            rep.assert_agree()                 # dim Hom (+ ker/im/coker when dim=1)
            seen_hom += 1
            if "kernel" in rep.ours:           # a canonical dim-1 hom was compared
                seen_kernel += 1
    assert seen_hom == min_pairs
    assert seen_kernel >= 1                     # the kernel/image/cokernel path fired
