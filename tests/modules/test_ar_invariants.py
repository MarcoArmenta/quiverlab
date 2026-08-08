"""AR-component invariants (Plan 57 / R21). Literature: kA_n / D_4 are
representation-directed (Gamma_A acyclic, every indecomposable directing, no
regular modules); the NON-uniform cyclic Nakayama kupisch=[3,2,2] is knittable,
NOT self-injective, and NOT representation-directed (the reachable negative branch,
verified live) yet IS generalized standard (rep-finite => rad^inf = 0) -- the two
flags are independent. Self-cert: partition total; projective postprojective,
injective preinjective; rep-directed <=> acyclic; generalized_standard with the
nilpotency witness. Honest: self-injective (uniform cyclic) / rep-infinite refuse
loudly."""
import pytest

from quiverlab import Quiver, linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ
from quiverlab.families.dynkin import dynkin_quiver     # NOT a top-level export
from quiverlab.modules.ar_invariants import ar_invariants, ar_invariants_block

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_linear_a_n_is_representation_directed(n):
    inv = ar_invariants(linear_path_algebra(n, field=QQ))
    assert inv.is_complete and inv.is_representation_directed
    assert all(p in ("postprojective", "preinjective") for p in inv.partition.values())
    assert "regular" not in inv.partition.values()          # Dynkin: no regular
    assert len(inv.directing) == n * (n + 1) // 2           # every indecomposable directing


@lit
def test_d4_is_representation_directed():
    inv = ar_invariants(dynkin_quiver("D4").algebra(relations=[], field=QQ))
    assert inv.is_complete and inv.is_representation_directed and len(inv.directing) == 12


@selfcert
def test_partition_places_projectives_and_injectives():
    A = linear_path_algebra(3, field=QQ)
    inv = ar_invariants(A)
    # the partition is TOTAL (every knit vertex classified)
    assert set(inv.partition) == set(range(len(inv.names)))
    # a directed component: every orbit reaches a projective, so all postprojective
    assert all(v == "postprojective" for v in inv.partition.values())


@selfcert
def test_generalized_standard_true_with_nilpotency_witness():
    inv = ar_invariants(linear_path_algebra(3, field=QQ))
    assert inv.generalized_standard is True and inv.nilpotency_index == 3


@lit
def test_cyclic_nakayama_is_knittable_but_NOT_representation_directed():
    # LIVE-VERIFIED WITNESS (the recognizer's reachable negative branch): the
    # non-uniform cyclic Nakayama kupisch=[3,2,2] on Z_3 is NOT self-injective, the
    # knit COMPLETES (7 indecomposables), and Gamma_A has an oriented cycle.
    A = NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)
    assert A.is_selfinjective() is False                     # verified live
    inv = ar_invariants(A)
    assert inv.is_complete and inv.status == "complete"      # knit accepts it
    assert len(inv.directing) < 7                            # some module lies on a cycle
    assert inv.is_representation_directed is False           # the negative branch
    # the two flags are INDEPENDENT: rep-finite => rad^inf = 0 => generalized standard
    assert inv.generalized_standard is True


@selfcert
def test_non_directed_partition_carries_the_honesty_caveat():
    # On a NON-directed component the τ-partition is not the clean trichotomy, so the
    # block/note must carry the plan-mandated caveat (never claim the clean partition).
    A = NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)
    inv = ar_invariants(A)
    assert inv.is_representation_directed is False
    assert "buckets may overlap" in inv.note                 # library-level note
    b = ar_invariants_block(A)
    assert b["representation_directed"] is False
    assert b.get("partition_note") == (
        "τ-orbit classification; non-directed component — buckets may overlap")
    # a DIRECTED component carries NO caveat (the clean trichotomy is honest there)
    bd = ar_invariants_block(linear_path_algebra(3, field=QQ))
    assert bd["representation_directed"] is True and "partition_note" not in bd


@selfcert
def test_self_injective_and_wild_refuse_loudly():
    si = ar_invariants(NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ))  # uniform => self-inj
    assert si.is_complete is False and si.status == "unsupported"
    # a small budget_dim trips the wild knit cap fast (preprojectives explode in dim)
    wild = ar_invariants(Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)})
                         .algebra(relations=[], field=QQ),
                         budget_modules=40, budget_dim=20)
    assert wild.is_complete is False and wild.status in ("budget", "error")
    assert wild.is_representation_directed is False          # never a false "directed"
    assert wild.generalized_standard is False                # no rad^inf verdict off-scope
