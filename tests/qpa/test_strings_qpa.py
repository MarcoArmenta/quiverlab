"""QPA as the oracle for the string family (Plan 46).

QPA has recognizers (IsGentleAlgebra / IsSpecialBiserialAlgebra) but NO string/band
enumerator and NO AG invariant -- so the crosschecks are: recognizer PARITY on the
string zoo, and MODULE-level decompose of a sum of string modules (validating our
materialisation at module level). The NamesGVars()-style probe is the standing guard
that FAILS if QPA ever ships a string surface (the Plan-35 precedent). qpa-marked."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.invariants import recognizers as rec
from quiverlab.modules.morphism import direct_sum
from quiverlab.qpa import scripts, session
from quiverlab.strings.modules import string_module
from quiverlab.strings.walks import enumerate_strings

pytestmark = [pytest.mark.qpa,
              pytest.mark.skipif(session.should_skip_qpa(),
                                 reason="[qpa] backend not installed")]


def _gentle_a3(field=GF(32003)):
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


def test_qpa_has_no_string_enumeration_surface():
    for name in ("StringModules", "BandModules", "AGInvariant"):
        bound = bool(session.run('IsBoundGlobal("%s");' % name))
        assert not bound, ("QPA now ships %s -- add a real crosscheck (honest scope "
                           "changed)" % name)


@pytest.mark.parametrize("name,build", [
    ("gentle_a3", lambda: _gentle_a3(GF(5))),
    ("two_cycle", lambda: Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=GF(5))),
    ("kronecker", lambda: Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(
        relations=[], field=GF(5))),
])
def test_recognizer_parity_on_string_zoo(name, build):
    A = build()
    decl = scripts.quiver_and_algebra_script(A)          # binds GAP variable `A`
    ours = (rec.is_special_biserial(A), rec.is_gentle(A))
    qpa = (bool(session.run(decl + "\nIsSpecialBiserialAlgebra(A);")),
           bool(session.run(decl + "\nIsGentleAlgebra(A);")))
    assert ours == qpa


def test_string_module_sum_decompose_vs_qpa():
    # a direct sum of a few string modules must decompose (QPA) back to the same
    # dimension-vector multiset. Over GF(32003) so char > dim (the decompose caveat).
    A = _gentle_a3(GF(32003))
    walks = list(enumerate_strings(A, max_length=6).walks)[:4]
    mods = [string_module(A, w) for w in walks]
    D = mods[0]
    for M in mods[1:]:
        D, _, _ = direct_sum(D, M)
    A.crosscheck("decompose", D).assert_agree()           # DecomposeModuleWithMultiplicities
