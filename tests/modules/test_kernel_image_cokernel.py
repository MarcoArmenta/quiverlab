"""Kernel/image/cokernel of a ModuleHom, self-certified by rank-nullity,
mono/epi flags, and the epi-mono factorization identity."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.morphism import hom_basis

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _all_homs(A):
    mods = [A.projective(v) for v in (1, 2, 3)] + [A.simple(v) for v in (1, 2, 3)]
    for M in mods:
        for N in mods:
            yield from hom_basis(M, N)


def test_rank_nullity_and_factorization_battery():
    A = _a3()
    seen = 0
    for f in _all_homs(A):
        K, iota = f.kernel()
        I, epi, mono = f.image()
        C, proj = f.cokernel()
        assert K.dim + f.rank() == f.src.dim              # rank-nullity
        assert I.dim == f.rank()
        assert C.dim == f.tgt.dim - f.rank()
        assert iota.is_mono() and epi.is_epi() and mono.is_mono() and proj.is_epi()
        assert epi.then(mono).matrix == f.matrix          # f = mono . epi
        assert iota.then(f).is_zero()                     # f . iota = 0
        assert f.then(proj).is_zero()                     # proj . f = 0
        seen += 1
    # Plan-37 threshold adjusted to reality: the 6-module x 6-module grid over this
    # Nakayama algebra kA3/(ab) has exactly 14 hom-basis elements total (verified
    # degreewise via A.hom; P3 is simple so P3 == S3). The plan's ">= 20" was an
    # overestimate; 12 keeps the "substantial nonempty battery" intent.
    assert seen >= 12                                     # battery is nonempty


def test_kernel_of_projective_cover_is_radical_syzygy():
    A = _a3()
    S1 = A.simple(1)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)
    from quiverlab.modules.morphism import ModuleHom
    f = ModuleHom(Q0, S1, d0)
    K, _ = f.kernel()
    # ker(P(S1) ->> S1) = rad P(S1)
    assert K.dimension_vector() == Q0.radical().dimension_vector()
