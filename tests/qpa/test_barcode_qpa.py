"""QPA scope for R33 (Plan 69). QPA has NO persistence/barcode/commutative-ladder verb;
the crosschecks are (1) a standing guard that FAILS if a persistence surface ever appears,
(2) DecomposeModuleWithMultiplicities parity on A_n -- barcode's interval bars carry the
SAME summand dim-vectors as QPA's decomposition (the barcode IS the interval decomposition),
and (3) the CL(3) diagram DECOMPOSE closer -- the diagram's summands (our decompose, the
engine that feeds the AR-index matching) agree with QPA on a concrete CL(3) module. The 29
whole-AR-quiver vertex COUNT is NOT a QPA verb (QPA has no AR-quiver enumeration), so it
stays oracle_selfcert (reconciled against the Escolar-Hiraoka figure, Task 6 Step 1a)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.families.commutative_ladder import CommutativeLadder
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.barcode import barcode
from quiverlab.modules.module import Module
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

_F = GF(7)                                       # char 7 > dim of every A_5 module below
FWD = {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5)}


def _direct_sum(mods, name="DS"):
    A = mods[0].algebra
    dom = A.domain
    labels = set().union(*(m.action.keys() for m in mods))
    dims = [m.dim for m in mods]
    N = sum(dims)
    offs, o = [], 0
    for d in dims:
        offs.append(o)
        o += d
    action = {}
    for lab in labels:
        Mb = lm.zeros(N, N, dom)
        for m, off in zip(mods, offs):
            blk = m.action.get(lab)
            if blk is None:
                continue
            for i in range(m.dim):
                for j in range(m.dim):
                    Mb[off + i][off + j] = blk[i][j]
        action[lab] = Mb
    return Module(A, N, action, name=name, side=mods[0].side)


def _multiset(A, pairs):
    """Sorted multiset of per-summand dim-vector tuples (quiver-vertex order), expanded
    by multiplicity -- the same reduction crosscheck_decompose uses."""
    verts = list(A.quiver.vertices)
    return sorted(tuple(dv.get(v, 0) for v in verts)
                  for dv, mult in pairs for _ in range(mult))


def test_qpa_has_no_persistence_surface():
    lg = session.libgap_handle()
    for name in ("PersistenceDiagram", "CommutativeLadder", "Barcode"):
        assert not bool(lg.eval('IsBoundGlobal("%s")' % name)), \
            "QPA now ships %s -- add a real crosscheck (honest scope changed)" % name


def test_barcode_parity_with_qpa_decompose_on_a5():
    # A forward A_5 interval-sum module; its barcode bars == our decompose summands ==
    # QPA's DecomposeModuleWithMultiplicities summands.
    A = Quiver([1, 2, 3, 4, 5], FWD).algebra(relations=[], field=_F)
    M = _direct_sum([A.projective(1), A.simple(3), A.injective(5), A.projective(2)],
                    "P1+S3+I5+P2")
    assert M.check_module()[0]
    A.crosscheck("decompose", M).assert_agree()      # our decompose == QPA
    bars = _multiset(A, [(bar.dimvec, bar.multiplicity) for bar in barcode(M).bars])
    dec = _multiset(A, [(s.dimension_vector(), m) for s, m in M.decompose()])
    assert bars == dec                               # barcode bars == decompose (== QPA)


def test_cl3_diagram_decompose_closer_qpa():
    # The CL(3) generalized diagram's summands come from OUR decompose (then matched to
    # AR-quiver vertices). Close THAT engine across engines on a concrete CL(3) module:
    # our decompose == QPA DecomposeModuleWithMultiplicities. Char > dim (GF(23) > 18) so
    # both engines decide. (The whole-AR-quiver 29-vertex COUNT is not a QPA verb.)
    A = CommutativeLadder(3, field=GF(23))
    M = A.projective(A.quiver.vertices[0])
    A.crosscheck("decompose", M).assert_agree()
