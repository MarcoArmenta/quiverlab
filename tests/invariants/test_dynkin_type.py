"""Type detection, round-tripped against the families/dynkin generators
and cross-checked against the Tits-form signature (Gabriel).

The generator (families/dynkin.py) tabulates finite ADE and Euclidean ~A
only (~A_n = a cycle on n+1 vertices); ~D / ~E it refuses, and A1 / ~A1
degenerate (empty quiver / dropped double edge). So the round-trips use
the generator where it is faithful and build ~D / ~E (and A1, Kronecker
= ~A1) explicitly -- these are still genuine detection oracles."""
import pytest

from quiverlab import Quiver
from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.invariants.dynkin_type import dynkin_type, is_connected
from quiverlab.invariants.forms import form_type

pytestmark = pytest.mark.oracle_crossengine


# --- explicit affine ~D / ~E builders (generator refuses them) --------------
def _affine_D(n):
    """The Euclidean D~_n diagram (n+1 vertices): two 2-leaf forks joined by a
    spine (n>=5); the 4-leaf degree-4 star for n=4."""
    if n == 4:
        return Quiver([0, 1, 2, 3, 4],
                      {"a": (0, 1), "b": (0, 2), "c": (0, 3), "d": (0, 4)})
    spine = list(range(3, 3 + (n - 5)))            # internal degree-2 spine
    forkB = 3 + (n - 5)
    arrows = {"la1": (0, 1), "la2": (0, 2)}
    chain = [0] + spine + [forkB]
    for k, (u, v) in enumerate(zip(chain, chain[1:])):
        arrows[f"s{k}"] = (u, v)
    arrows["lb1"] = (forkB, forkB + 1)
    arrows["lb2"] = (forkB, forkB + 2)
    return Quiver(list(range(n + 1)), arrows)      # D~_n has n+1 vertices (0..n)


def _affine_E(n):
    """The Euclidean E~_n diagram: a single degree-3 center with arms
    (2,2,2) / (1,3,3) / (1,2,5) for n = 6 / 7 / 8."""
    arms = {6: (2, 2, 2), 7: (1, 3, 3), 8: (1, 2, 5)}[n]
    verts, arrows, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for step in range(length):
            verts.append(nxt)
            arrows[f"a{ai}_{step}"] = (prev, nxt)
            prev, nxt = nxt, nxt + 1
    return Quiver(verts, arrows)


FINITE = [("A", n) for n in range(2, 9)] + [("D", n) for n in range(4, 9)] \
         + [("E", 6), ("E", 7), ("E", 8)]


@pytest.mark.parametrize("typ,n", FINITE, ids=[f"{t}{n}" for t, n in FINITE])
def test_roundtrip_finite(typ, n):
    Q = dynkin_quiver(f"{typ}{n}")
    assert dynkin_type(Q) == (typ, n)


def test_a1_single_vertex():
    # the generator makes A1 an empty quiver; A1 is one isolated node
    assert dynkin_type(Quiver([1], {})) == ("A", 1)


AFFINE_A = [("~A", n) for n in range(2, 7)]


@pytest.mark.parametrize("typ,n", AFFINE_A, ids=[f"~A{n}" for _, n in AFFINE_A])
def test_roundtrip_affine_A(typ, n):
    Q = dynkin_quiver(f"~A{n}")                     # generator: cycle on n+1 vertices
    assert dynkin_type(Q) == (typ, n)


AFFINE_D = [("~D", n) for n in range(4, 8)]


@pytest.mark.parametrize("typ,n", AFFINE_D, ids=[f"~D{n}" for _, n in AFFINE_D])
def test_detect_affine_D(typ, n):
    assert dynkin_type(_affine_D(n)) == (typ, n)


AFFINE_E = [("~E", 6), ("~E", 7), ("~E", 8)]


@pytest.mark.parametrize("typ,n", AFFINE_E, ids=[f"~E{n}" for _, n in AFFINE_E])
def test_detect_affine_E(typ, n):
    assert dynkin_type(_affine_E(n)) == (typ, n)


def test_kronecker_is_affine_a1():
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    assert dynkin_type(Q) == ("~A", 1)


def test_agreement_with_tits_signature():
    # Dynkin <=> finite form; affine <=> tame form (hereditary, Gabriel)
    for typ, n in [("A", 4), ("D", 5), ("E", 6)]:
        A = dynkin_quiver(f"{typ}{n}").algebra()
        assert dynkin_type(A.quiver) is not None
        assert dynkin_type(A.quiver)[0] in ("A", "D", "E")
        assert form_type(A) == "finite"


def test_wild_star_is_none():
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (2, 1), "b": (3, 1), "c": (4, 1), "d": (5, 1), "e": (6, 1)})
    assert dynkin_type(Q) is None          # 5-star: not ADE, not affine


def test_loop_is_none():
    Q = Quiver([1], {"x": (1, 1)})
    assert dynkin_type(Q) is None          # Jordan quiver is ~A0, not tabulated


def test_disconnected_is_none_and_flagged():
    Q = Quiver([1, 2, 3], {"a": (1, 2)})
    assert is_connected(Q) is False
    assert dynkin_type(Q) is None
    assert Q.is_connected() is False       # the thin delegating method
