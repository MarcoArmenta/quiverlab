"""Radical-filtration spectral sequence. Self-cert: converges to H(X) always
(the standing certificate). Koszul tie-in is an ARBITRATED oracle hypothesis:
for a Koszul algebra the radical-filtration SS of the minimal simple-resolution
degenerates at E_2 (the linear staircase on E_1 collapses) -- the provable page,
pinned below (Design-decision 3: never force the oracle; the observed page IS the
folklore E_2 collapse, so it is asserted exactly)."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.presets import radical_filtration_ss

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_converges_to_homology():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = radical_filtration_ss(X)                 # __init__ self-certifies
    assert ss.convergence.abutment == X.homology_dims()


def test_semisimple_input_degenerates_immediately():
    # a complex of SEMISIMPLE modules (rad = 0) has a trivial radical filtration:
    # F_0 = whole, F_{-1} = 0, so E_1 = E_inf (collapse).
    A = Quiver([1, 2], {}).algebra(relations=[], field=GF(7))   # semisimple k x k
    X = ChainComplex.stalk(A.simple(1), 0)
    ss = radical_filtration_ss(X)
    assert ss.convergence.collapse() is True
    assert ss.convergence.degenerates_at == 1


@pytest.mark.oracle_literature
def test_koszul_degeneration_arbitrated():
    # kA_n is Koszul (hereditary => the minimal resolution of a simple is linear).
    # The radical-filtration SS of the minimal projective resolution of S_1
    # degenerates at E_2 = E_inf (ARBITRATED: this IS the observed provable page,
    # the folklore "E_2 collapse" -- recorded on the verification page).
    A = linear_path_algebra(3, field=GF(5))       # hereditary kA3, Koszul
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    ss = radical_filtration_ss(X)
    assert ss.convergence.degenerates_at == 2
    assert ss.convergence.collapse() is True
    assert ss.convergence.abutment == X.homology_dims()   # always true (self-cert)


def test_koszul_degeneration_kA4():
    # the same arbitrated E_2 collapse on kA4 (a second Koszul witness).
    A = linear_path_algebra(4, field=GF(5))
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = radical_filtration_ss(X)
    assert ss.convergence.degenerates_at == 2
    assert ss.convergence.abutment == X.homology_dims()
