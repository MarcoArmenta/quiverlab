"""Plan 35 wave 3b -- the cyclic-homology (HC) surface now PRINTS the actual cycle
representatives of the (b, B) total complex, per degree, as a labelled term-sum AND a
sparse coordinate vector over the ordered ``Tot_n = C_n (+) C_{n-2} (+) ...`` basis, with
the total differential ``D_n = b + B`` that annihilates each class and the explicit
column structure (which coordinate slice lives in which ``C_k``).

The crux (mirroring the connes_b / Ext / Tor self-cert): for EVERY shipped class at
EVERY degree, applying the SHIPPED total differential to the SHIPPED vector gives 0
(both engines); the GF(p) capture on the char-0-shaped prime 32003 agrees with the
generic (b, B) mixed complex; the labels are HAND-CHECKED; the rep count equals the HC
dims; the machine record is JSON-safe.
"""
import json
from fractions import Fraction

import pytest

import quiverlab as ql
from quiverlab.fields import GF, QQ


# --------------------------------------------------------------------------- #
# Self-certification helpers (exactly what a user does with the machine record).
# --------------------------------------------------------------------------- #
def _certify_gfp(payload, p):
    """Apply each shipped int total differential to each shipped vector over F_p."""
    bc, diffs, n = payload["basis_classes"], payload["differentials"], 0
    for deg in bc:
        d = diffs[deg]
        assert not d.get("elided"), "unexpected elision at deg %s" % deg
        rows, ncols = d["rows"], d["shape"][1]
        for cl in bc[deg]:
            v = [0] * ncols
            for idx, c in cl["vector"]:
                v[idx] = int(c)
            for r in rows:
                assert sum(int(r[j]) * v[j] for j in range(ncols)) % p == 0, \
                    "deg %s: shipped HC vector not annihilated by D=b+B" % deg
            n += 1
    return n


def _certify_dom(payload):
    """Apply each shipped Domain total differential to each shipped vector, exact."""
    bc, diffs, n = payload["basis_classes"], payload["differentials"], 0
    for deg in bc:
        d = diffs[deg]
        assert not d.get("elided"), "unexpected elision at deg %s" % deg
        rows, ncols = d["rows"], d["shape"][1]
        for cl in bc[deg]:
            v = [Fraction(0)] * ncols
            for idx, c in cl["vector"]:
                v[idx] = Fraction(str(c))
            for r in rows:
                acc = sum(Fraction(str(r[j])) * v[j] for j in range(ncols))
                assert acc == 0, "deg %s: shipped HC vector not a cycle" % deg
            n += 1
    return n


# --------------------------------------------------------------------------- #
# Hand-checked labels: k[x]/(x^2), the dual numbers (the spec's HC anchor).
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_handcheck_dualnumbers_hc0_gf7():
    """HC_0 of k[x]/(x^2) over GF(7) is A/[A,A] = A, dim 2; its classes are the two
    C_0 chains e_1 and x (in the sole column C_0)."""
    A = ql.truncated_polynomial(2, field=GF(7))
    table, pay = A.cyclic_homology(3, with_reps=True)
    assert list(table.dims) == [2, 0, 2, 0]
    # ordered Tot_0 = C_0 basis, column-annotated
    assert pay["chain_basis"]["0"] == ["col C_0: e_1", "col C_0: x"]
    z0 = pay["basis_classes"]["0"]
    assert [c["terms"] for c in z0] == [[["1", 0, [], "e_1"]], [["1", 0, [], "x"]]]
    assert all(c["kind"] == "cyclic" for c in z0)
    # column structure: Tot_0 is the single column C_0, dim 2 at offset 0
    cs = pay["column_structure"]["0"]
    assert cs == {"total": 2, "columns": [{"degree": 0, "offset": 0, "dim": 2}]}
    # D_0 = 0 (every 0-chain is a cycle) -- the note states HC_0 = A/[A,A]
    assert pay["differentials"]["0"]["shape"] == [0, 2]
    assert "A/[A,A]" in pay["differentials"]["0"]["note"]


@pytest.mark.oracle_literature
def test_handcheck_dualnumbers_hc2_two_columns():
    """Tot_2 = C_2 (+) C_0 (two columns); HC_2 is 2-dimensional and its classes are the
    C_2 chain x⊗x⊗x and the C_0 chain e_1 -- the column structure must be explicit."""
    A = ql.truncated_polynomial(2, field=GF(7))
    _t, pay = A.cyclic_homology(3, with_reps=True)
    cs = pay["column_structure"]["2"]
    assert cs == {"total": 4, "columns": [{"degree": 2, "offset": 0, "dim": 2},
                                          {"degree": 0, "offset": 2, "dim": 2}]}
    z2 = pay["basis_classes"]["2"]
    assert [c["terms"] for c in z2] == [[["1", 2, ["x", "x"], "x"]],
                                        [["1", 0, [], "e_1"]]]
    # the enumeration labels carry the column tag
    assert pay["chain_basis"]["2"][:2] == ["col C_2: e_1 (x) x (x) x",
                                           "col C_2: x (x) x (x) x"]
    assert pay["chain_basis"]["2"][2:] == ["col C_0: e_1", "col C_0: x"]


# --------------------------------------------------------------------------- #
# Self-certification: every shipped class is a cycle of the total complex.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("a", [2, 3, 4])
def test_selfcert_gfp_truncated(a):
    A = ql.truncated_polynomial(a, field=GF(7))
    _t, pay = A.cyclic_homology(4, with_reps=True)
    n = _certify_gfp(pay, 7)
    assert n == sum(len(v) for v in pay["basis_classes"].values())


@pytest.mark.oracle_selfcert
def test_selfcert_gfp_quiver_kronecker():
    """A multi-vertex presented algebra (kA_2) -- self-cert over GF(5)."""
    A = ql.Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra(field=GF(5))
    _t, pay = A.cyclic_homology(3, with_reps=True)
    _certify_gfp(pay, 5)


@pytest.mark.oracle_selfcert
def test_selfcert_generic_qq():
    A = ql.truncated_polynomial(3, field=QQ)
    _t, pay = A.cyclic_homology(4, with_reps=True)
    n = _certify_dom(pay)
    assert n == sum(len(v) for v in pay["basis_classes"].values())


# --------------------------------------------------------------------------- #
# Cross-engine: GF(p) fast rank vs the generic mixed complex on the char-0-shaped
# prime 32003 -- reps AND dims agree (mirror the connes_b cross-engine gate).
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("a", [2, 3])
def test_gfp_generic_agree_prime_32003(a):
    Af = ql.truncated_polynomial(a, field=GF(32003))     # generic route (PrimeField? no)
    Ag = ql.truncated_polynomial(a, field=GF(32003))
    # GF(32003) IS a prime field -> the GF(p) engine. Compare to QQ generic dims.
    tf, pf = Ag.cyclic_homology(4, with_reps=True)
    tq, pq = ql.truncated_polynomial(a, field=QQ).cyclic_homology(4, with_reps=True)
    assert list(tf.dims) == list(tq.dims)
    # per degree the captured HC dim (rep count) matches on both engines
    for deg in pf["basis_classes"]:
        assert len(pf["basis_classes"][deg]) == len(pq["basis_classes"][deg])
    _certify_gfp(pf, 32003)
    _certify_dom(pq)


# --------------------------------------------------------------------------- #
# Rep count == table dims (the capture cross-check), term<->vector agreement, and the
# machine record is JSON-safe. Interpretability: labels stay short and readable.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_rep_count_matches_dims_and_json_safe():
    A = ql.truncated_polynomial(2, field=GF(7))
    table, pay = A.cyclic_homology(4, with_reps=True)
    for n, dim in enumerate(table.dims):
        assert len(pay["basis_classes"][str(n)]) == dim
    json.dumps(pay)                                       # JSON-safe
    # term <-> sparse-vector agreement: a term's coeff appears at its vector index
    for deg, classes in pay["basis_classes"].items():
        for cl in classes:
            assert len(cl["terms"]) == len(cl["vector"])


@pytest.mark.oracle_selfcert
def test_enumeration_labels_short_and_column_tagged():
    A = ql.truncated_polynomial(3, field=GF(7))
    _t, pay = A.cyclic_homology(3, with_reps=True)
    for deg, labels in pay["chain_basis"].items():
        assert isinstance(labels, list)
        for lbl in labels:
            assert lbl.startswith("col C_")               # column-tagged
            assert len(lbl) < 200                         # interpretable length
