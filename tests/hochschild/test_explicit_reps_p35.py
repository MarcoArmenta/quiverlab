"""Plan 35 explicit representatives -- the product surface now PRINTS the actual
(co)cycle that produced each structure constant, as a labeled term-sum AND a sparse
coordinate vector, alongside the ordered (co)chain enumeration and the annihilating
differential, so a reader can verify each class is a genuine (co)cycle.

The crux (Marco): for EVERY shipped class at EVERY degree, applying the SHIPPED (or
rebuilt-when-elided) differential to the SHIPPED vector gives 0 -- exactly what a user
would do. Plus the labels are HAND-CHECKED, the two views (term-sum / sparse vector)
agree, the rep count matches the table dims, and the shipped reps reproduce the
shipped constants through the raw engine product.
"""
import numpy as np
import pytest

import quiverlab as ql
from quiverlab.fields import QQ
from quiverlab.hochschild import basis_reps as BR


# --------------------------------------------------------------------------- #
# Shared self-certification helpers (what a user does with the machine record).
# --------------------------------------------------------------------------- #
def _dense_int(vector, ncols, p):
    v = [0] * ncols
    for idx, c in vector:
        v[idx] = int(c) % p
    return v


def _certify_gfp(block, p):
    """Apply each shipped int differential to each shipped vector over F_p; assert 0.
    Returns the number of classes certified (differential must NOT be elided here)."""
    diffs, bc, n = block["differentials"], block["basis_classes"], 0
    for side in bc:
        for deg in bc[side]:
            d = diffs[side][deg]
            assert not d.get("elided"), "unexpected elision at %s deg %s" % (side, deg)
            rows, ncols = d["rows"], d["shape"][1]
            for cl in bc[side][deg]:
                v = _dense_int(cl["vector"], ncols, p)
                for r in rows:
                    assert sum(int(r[j]) * v[j] for j in range(ncols)) % p == 0, \
                        "%s deg %s: shipped vector is not annihilated by delta/b" % (side, deg)
                n += 1
    return n


def _certify_dom(block, dom):
    diffs, bc, n = block["differentials"], block["basis_classes"], 0
    for side in bc:
        for deg in bc[side]:
            d = diffs[side][deg]
            assert not d.get("elided"), "unexpected elision at %s deg %s" % (side, deg)
            rows, ncols = d["rows"], d["shape"][1]
            for cl in bc[side][deg]:
                v = [dom.zero()] * ncols
                for idx, c in cl["vector"]:
                    v[idx] = dom.coerce(c)
                for r in rows:
                    acc = dom.zero()
                    for j in range(ncols):
                        acc = dom.add(acc, dom.mul(dom.coerce(r[j]), v[j]))
                    assert dom.is_zero(acc), "%s deg %s not a (co)cycle" % (side, deg)
                n += 1
    return n


# --------------------------------------------------------------------------- #
# Hand-checked labels: k[x]/(x^2). The spec's anchor.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_handcheck_kx2_gf7_cohomology_class():
    """HH^1 of k[x]/(x^2) over GF(7) is 1-dimensional; its class is the explicit
    cochain [x -> x] with coefficient 1 (spec hand-check anchor)."""
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    b = A.cup_products(2, engine="bar").blocks()
    cls = b["basis_classes"]["coh"]["1"]
    assert len(cls) == 1
    assert cls[0]["terms"] == [["1", ["x"], "x"]]         # alpha^1_1 = [x -> x]
    assert cls[0]["kind"] == "cochain" and cls[0]["degree"] == 1
    assert b["chain_basis"]["coh"]["1"] == ["[x -> e_1]", "[x -> x]"]


@pytest.mark.oracle_literature
def test_handcheck_kx2_gf2_unit_valued_term():
    """Over GF(2) the same HH^1 is 2-dimensional and one class is the unit-valued
    [x -> e_1] ('[x -> 1]'-type), the other [x -> x]."""
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    b = A.cup_products(1, engine="bar").blocks()
    cls = b["basis_classes"]["coh"]["1"]
    assert [c["terms"] for c in cls] == [[["1", ["x"], "e_1"]], [["1", ["x"], "x"]]]


@pytest.mark.oracle_literature
def test_handcheck_kx2_gf7_homology_cycle():
    """HH_1 of k[x]/(x^2) over GF(7): the cycle z^1_1 = (e_1 (x) x); HH_0 = {e_1, x}."""
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    b = A.cap_products(2, engine="bar").blocks()
    z1 = b["basis_classes"]["hom"]["1"]
    assert len(z1) == 1 and z1[0]["terms"] == [["1", ["x"], "e_1"]]   # e_1 (x) x
    assert z1[0]["kind"] == "chain"
    z0 = b["basis_classes"]["hom"]["0"]
    assert [c["terms"] for c in z0] == [[["1", [], "e_1"]], [["1", [], "x"]]]


@pytest.mark.oracle_literature
def test_handcheck_cs_kx3_labels():
    """CS route, k[x]/(x^3) over QQ: HH^1 classes [x -> x], [x -> x*x]; the value
    labels come from the algebra's basis ['e_1','x','x*x']."""
    A = ql.truncated_polynomial(3, field=QQ)
    b = A.cup_products(2, engine="cs").blocks()
    cls = b["basis_classes"]["coh"]["1"]
    assert [c["terms"] for c in cls] == [[["1", ["x"], "x"]], [["1", ["x"], "x*x"]]]
    assert b["basis"] == "cs/QQ"


# --------------------------------------------------------------------------- #
# Self-certification: delta @ vec = 0 / b @ vec = 0 from the SHIPPED record.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("kind", ["cup", "cap", "bracket"])
def test_selfcert_gfp_products(kind):
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    method = {"cup": A.cup_products, "cap": A.cap_products,
              "bracket": A.gerstenhaber_brackets}[kind]
    b = method(3, engine="bar").blocks() if kind != "bracket" else method(3).blocks()
    assert _certify_gfp(b, 7) > 0


@pytest.mark.oracle_selfcert
def test_selfcert_gfp_connes():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    assert _certify_gfp(A.connes_differentials(3).blocks(), 7) > 0


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("kind", ["cup", "cap"])
def test_selfcert_cs_products(kind):
    A = ql.truncated_polynomial(3, field=QQ)
    b = (A.cup_products if kind == "cup" else A.cap_products)(2, engine="cs").blocks()
    assert _certify_dom(b, QQ) > 0


@pytest.mark.oracle_selfcert
def test_selfcert_generic_connes():
    A = ql.truncated_polynomial(2, field=QQ)
    b = A.connes_differentials(2).blocks()
    assert "generic" in b["engine"]
    assert _certify_dom(b, QQ) > 0


# --------------------------------------------------------------------------- #
# The two views agree: sparse vector <-> term-sum <-> ordered enumeration.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("call", [
    lambda A: A.cup_products(2, engine="bar"),
    lambda A: A.cap_products(2, engine="bar"),
    lambda A: A.connes_differentials(2),
])
def test_sparse_terms_enumeration_consistent(call):
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    b = call(A).blocks()
    for side in b["basis_classes"]:
        for deg in b["basis_classes"][side]:
            enum = b["chain_basis"][side][deg]
            for cl in b["basis_classes"][side][deg]:
                assert len(cl["vector"]) == len(cl["terms"])
                for (idx, c), (tc, w, v) in zip(cl["vector"], cl["terms"]):
                    assert c == tc
                    assert BR.element_label(tuple(w), v, cl["kind"]) == enum[idx]


# --------------------------------------------------------------------------- #
# Rep count == table dims, and the shipped reps reproduce the shipped constants.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_crossengine
def test_rep_count_matches_table_dims():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    obj = A.cup_products(3, engine="bar")
    b = obj.blocks()
    for t in b["tables"]:
        p, q = t["degrees"]
        dl, dr, dout = t["dims"]
        assert len(b["basis_classes"]["coh"][str(p)]) == dl
        assert len(b["basis_classes"]["coh"][str(q)]) == dr
        assert len(b["basis_classes"]["coh"][str(p + q)]) == dout


@pytest.mark.oracle_crossengine
def test_shipped_reps_reproduce_shipped_constants():
    """Feed the SHIPPED cup reps through the raw engine cup product and read off the
    class coordinates; they must equal the SHIPPED structure constants (reps <->
    constants coherence -- the reps really are the ones that made the table)."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine import tt_calculus as TT
    from quiverlab.engine.scan3 import cochain_basis
    A = ql.truncated_polynomial(3, field=ql.GF(5))
    p = 5
    E = to_engine(A.unit_adapted())
    obj = A.cup_products(2, engine="bar")
    b = obj.blocks()
    bc = b["basis_classes"]["coh"]

    def dense(deg, cl):
        v = np.zeros(len(cochain_basis(E, deg)), dtype=np.int64)
        for idx, c in cl["vector"]:
            v[idx] = int(c) % p
        return v

    for t in b["tables"]:
        dp, dq = t["degrees"]
        out = dp + dq
        Hout = TT.cohomology_classes(E, out, p)
        for i, ci in enumerate(bc[str(dp)]):
            for j, cj in enumerate(bc[str(dq)]):
                h = TT.cup_cochain(E, dp, dq, dense(dp, ci), dense(dq, cj))
                coords = [int(x) for x in Hout.coords(h)]
                got = [int(t["constants"][k][i][j]) for k in range(t["dims"][2])]
                assert coords == got, (dp, dq, i, j, coords, got)


# --------------------------------------------------------------------------- #
# Elision: at qci tops the bar cochain spaces blow past the 250k cap; the
# cap+note path MUST fire, and rebuilding from the note certifies the class.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_elided_differential_note_rebuild_certifies():
    Q = ql.QuantumCI(2, 2, 2, field=ql.GF(5))       # dim 4; delta^4 is 972 x 324
    b = Q.cup_products(4, engine="bar").blocks()
    d4 = b["differentials"]["coh"]["4"]
    assert d4.get("elided") is True and d4["shape"] == [972, 324]
    assert "scan3.coboundary_matrix" in d4["note"]
    assert BR.MATRIX_CELL_CAP == 250_000 and 972 * 324 > BR.MATRIX_CELL_CAP
    # Rebuild delta^4 from the note's documented API and certify every deg-4 class.
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.scan3 import cochain_basis, coboundary_matrix
    E = to_engine(Q.unit_adapted())
    idx = {g: i for i, g in enumerate(cochain_basis(E, 5))}
    M = coboundary_matrix(E, 4, cochain_basis(E, 4), idx)
    cls = b["basis_classes"]["coh"]["4"]
    assert cls, "expected degree-4 classes to certify"
    for cl in cls:
        v = np.array(_dense_int(cl["vector"], 324, 5), dtype=np.int64)
        assert ((M @ v) % 5 == 0).all()


# --------------------------------------------------------------------------- #
# JSON-safe + tolerance for legacy blocks without the explicit-reps fields.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_blocks_json_safe_and_have_fields():
    import json
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    b = A.cap_products(2, engine="bar").blocks()
    json.dumps(b)                                    # no non-JSON objects
    for key in ("basis_classes", "chain_basis", "differentials"):
        assert key in b and "coh" in b[key] and "hom" in b[key]


@pytest.mark.oracle_selfcert
def test_legacy_block_without_reps_is_tolerated():
    """A block constructed without reps (a cached/old result) omits the fields
    cleanly -- the constructor and blocks() stay backward compatible."""
    from quiverlab.hochschild.products import HHProducts, ProductTable
    t = ProductTable(kind="cup", degrees=(0, 0), out_degree=0, dims=(1, 1, 1),
                     constants=(((("1",),),),))
    hp = HHProducts(kind="cup", top=0, tables={(0, 0): t},
                    engine="x", basis="bar/GF(7)", window=None, references=[])
    b = hp.blocks()
    assert "basis_classes" not in b and "chain_basis" not in b
    assert "differentials" not in b
