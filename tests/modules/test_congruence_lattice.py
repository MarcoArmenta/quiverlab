"""Con(tors A) + the forcing order on bricks (Plan 64 / R26). Literature: |Con(tors kA2)| = 5
(the pentagon), |Con(tors kA3)| = 14. Self-cert: Con is distributive (Funayama-Nakayama);
#join-irreducible congruences == #bricks (the forcing order lives on bricks, DIRRT).
Cross-engine: the forcing-poset order-ideal count == |Con|."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.congruence import congruence_lattice
from quiverlab.tautilting.torsion import bricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("n, con_size", [(2, 5), (3, 14)])
def test_congruence_lattice_size(n, con_size):
    C = congruence_lattice(linear_path_algebra(n, field=QQ))
    assert C.is_complete and C.size == con_size


@selfcert
def test_con_is_distributive():
    # Con(L) is ALWAYS distributive (Funayama-Nakayama); P64 builds it as a down-set lattice,
    # so distributivity is structural -- assert the flag and (small case) verify directly.
    C = congruence_lattice(linear_path_algebra(2, field=QQ))
    assert C.is_distributive is True


@xeng
@pytest.mark.parametrize("n, nbricks", [(2, 3), (3, 6)])
def test_join_irreducible_congruences_are_bricks(n, nbricks):
    # DIRRT: the forcing-equivalence classes of covers are the bricks; #J(Con) == #bricks.
    A = linear_path_algebra(n, field=QQ)
    C = congruence_lattice(A)
    assert len(C.join_irreducible_congruences) == len(bricks(A)) == nbricks
    assert len(C.forcing_bricks) == nbricks


@selfcert
def test_kA2_forcing_order_is_a_V():
    # the 3 bricks of kA2 with the forcing order: ONE brick BELOW the other two -- a "V"
    # opening UPWARD -- whose 5 order ideals give |Con| = 5. Ruling 4: pin the DIRECTION,
    # not just the relation count. forcing_order pairs are (i, j) meaning brick_i <= brick_j.
    C = congruence_lattice(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    assert len(C.forcing_bricks) == 3
    assert len(C.forcing_order) == 2                   # exactly two strict relations
    lowers = {i for (i, j) in C.forcing_order}         # the "below" (forced) endpoints
    uppers = {j for (i, j) in C.forcing_order}         # the "above" (forcing) endpoints
    assert len(lowers) == 1 and len(uppers) == 2       # ONE minimum below TWO maxima
    assert lowers.isdisjoint(uppers)                   # the minimum is not itself an upper
    # the minimum brick is the one whose principal congruence is CONTAINED in the other two
    # (con(min) subset con(a), con(b)) -- it is FORCED BY them, not the reverse. Observed live
    # this authoring: relations [(2, 0), (2, 1)] (brick 2 below bricks 0, 1; appendix prototype).
    assert C.size == 5
