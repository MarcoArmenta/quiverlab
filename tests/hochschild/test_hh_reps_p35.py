"""Plan 35 wave 3d -- the PLAIN hh_cohomology / hh_homology blocks now carry the
explicit per-degree representatives, so the classical dictionary read-offs (HH^0's
central elements, HH^1's derivations D(arrow)=value, HH^2's deformation cochain, HH_0's
commutator residues) are backed by actual (co)cycles, not just framing prose.

The crux (mirroring the product / Ext-Tor / HC batteries): for EVERY shipped class at
EVERY captured degree, applying the SHIPPED differential to the SHIPPED coordinate
vector gives 0 -- over BOTH routes (the GF(p) bar route and the Chouhy-Solotar route).
Plus the k[x]/(x^2) HH^1 hand-check ties the derivation display to the captured
[x -> x], the rep count matches the table dims, and the term-sum / sparse-vector views
agree.
"""
from fractions import Fraction

import pytest

import quiverlab as ql
from quiverlab.fields import QQ
from quiverlab.fields.primefield import PrimeField
from quiverlab.hochschild.hh_reps import hh_reps_blocks
from quiverlab.trace import interpretations as I


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _capture(A, kind, top):
    tbl = (A.hochschild_cohomology if kind == "hh_cohomology"
           else A.hochschild_homology)(top, verbose=False)
    blk = hh_reps_blocks(A, kind, top, list(tbl.dims), tbl.engine)
    return list(tbl.dims), tbl.engine, blk


def _dense(vector, ncols):
    v = [Fraction(0)] * ncols
    for idx, c in vector:
        v[idx] = Fraction(str(c))
    return v


def _certify(block, p=None):
    """Apply each shipped differential to each shipped vector; assert 0 (mod p when a
    prime field). Returns the number of classes certified."""
    diffs, bc, n = block["differentials"], block["basis_classes"], 0
    for deg in bc:
        d = diffs[deg]
        if d.get("elided"):
            continue
        rows = d.get("rows") or []
        ncols = (d.get("shape") or [0, 0])[1]
        for cl in bc[deg]:
            v = _dense(cl["vector"], ncols)
            for r in rows:
                s = sum((Fraction(str(r[j])) * v[j] for j in range(ncols)), Fraction(0))
                if p is not None:
                    assert int(s.numerator) % p == 0, (deg, s)
                else:
                    assert s == 0, (deg, s)
            n += 1
    return n


# --------------------------------------------------------------------------- #
# Self-certification -- both routes.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("a,p", [(2, 7), (3, 2), (4, 5)])
@pytest.mark.parametrize("kind", ["hh_cohomology", "hh_homology"])
def test_gfp_route_self_certifies(a, p, kind):
    """k[x]/(x^a) over GF(p): the F_p fast dims path -> the GF(p) bar reps route; every
    shipped vector annihilates the shipped differential."""
    A = ql.truncated_polynomial(a, field=PrimeField(p))
    dims, engine, blk = _capture(A, kind, 4)
    assert "hanlab" in engine.lower(), engine          # the GF(p) route is exercised
    assert blk is not None
    assert _certify(blk, p) >= 1
    for deg, classes in blk["basis_classes"].items():
        assert len(classes) == dims[int(deg)]          # rep count == table dims


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("a", [3, 4])
@pytest.mark.parametrize("kind", ["hh_cohomology", "hh_homology"])
def test_cs_route_self_certifies(a, kind):
    """k[x]/(x^a) over QQ: no F_p route -> the Chouhy-Solotar reps route; every shipped
    vector annihilates the shipped differential over the Domain."""
    A = ql.truncated_polynomial(a, field=QQ)
    dims, engine, blk = _capture(A, kind, 4)
    assert "hanlab" not in engine.lower()
    assert blk is not None
    assert _certify(blk, None) >= 1
    for deg, classes in blk["basis_classes"].items():
        assert len(classes) == dims[int(deg)]


# --------------------------------------------------------------------------- #
# Hand-check: k[x]/(x^2) HH^1 = [x -> x], read off as D(x) = x (the derivation display
# tied to the captured representative). The reference anchor of the whole wave.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_dual_numbers_HH1_is_derivation_x_maps_to_x():
    A = ql.truncated_polynomial(2, field=PrimeField(7))
    dims, engine, blk = _capture(A, "hh_cohomology", 3)
    assert dims[1] == 1
    cls = blk["basis_classes"]["1"]
    assert len(cls) == 1
    # the single HH^1 class is the term-sum [x -> x]
    assert cls[0]["terms"] == [["1", ["x"], "x"]]
    # the classical dictionary reads it as the derivation D(x) = x
    assert I.element_readoff("hh_cohomology", 1, cls[0]["terms"]) == ["D(x) = x"]
    # dual numbers are commutative -> no inner derivations: dim Inn = rank delta^0 = 0
    assert blk["inner_dims"]["1"] == 0


@pytest.mark.oracle_literature
def test_dual_numbers_HH0_lists_the_centre():
    """HH^0 = Z(A) = A (commutative): the two central elements are e_1 and x."""
    A = ql.truncated_polynomial(2, field=PrimeField(7))
    _dims, _engine, blk = _capture(A, "hh_cohomology", 2)
    reads = [I.central_element(c["terms"]) for c in blk["basis_classes"]["0"]]
    assert reads == ["e_1", "x"]


@pytest.mark.oracle_literature
def test_dual_numbers_HH0_homology_is_commutator_quotient():
    """HH_0 = A/[A,A]: for a commutative algebra every element survives; the residues
    are e_1 and x."""
    A = ql.truncated_polynomial(2, field=PrimeField(7))
    _dims, _engine, blk = _capture(A, "hh_homology", 2)
    reads = [I.commutator_residue(c["terms"]) for c in blk["basis_classes"]["0"]]
    assert reads == ["e_1", "x"]


@pytest.mark.oracle_literature
def test_cs_route_dual_numbers_x3_HH1_read_off():
    """k[x]/(x^3) over QQ (CS route): HH^1 has the two derivations D(x) = x and
    D(x) = x*x -- read off the captured CS representatives."""
    A = ql.truncated_polynomial(3, field=QQ)
    _dims, _engine, blk = _capture(A, "hh_cohomology", 2)
    reads = []
    for c in blk["basis_classes"]["1"]:
        reads.extend(I.element_readoff("hh_cohomology", 1, c["terms"]))
    assert reads == ["D(x) = x", "D(x) = x*x"]


# --------------------------------------------------------------------------- #
# Term-sum / coordinate-vector coherence + honest no-route path.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_term_sum_and_vector_agree():
    A = ql.truncated_polynomial(3, field=PrimeField(5))
    _dims, _engine, blk = _capture(A, "hh_cohomology", 3)
    for deg, classes in blk["basis_classes"].items():
        enum = blk["chain_basis"][deg]
        for cl in classes:
            # one nonzero term per nonzero coordinate, indices agree, in the same order
            assert len(cl["terms"]) == len(cl["vector"])
            for (coeff_t, word, value), (idx, coeff_v) in zip(cl["terms"], cl["vector"]):
                assert coeff_t == coeff_v
                if not isinstance(enum, dict):          # not display-elided
                    assert value in enum[idx] and (not word or word[0] in enum[idx])


@pytest.mark.oracle_selfcert
def test_no_route_keeps_dims_only():
    """A structure-constants algebra off a prime field (no path basis for the GF(p) bar
    route, no quiver presentation for the CS route) has NEITHER reps route ->
    hh_reps_blocks returns None (the block stays dims-only, honestly)."""
    from quiverlab.core.algebra import Algebra
    from quiverlab.fields import GF
    # k[x]/(x^2) as a 2-dim structure-constants algebra over GF(4) = GF(2^2), built
    # WITHOUT a quiver presentation (A.quiver is None) and off a prime field.
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, [1, 0], field=GF(4),
                                         basis_labels=["e_1", "x"])
    tbl = A.hochschild_cohomology(2, verbose=False)
    blk = hh_reps_blocks(A, "hh_cohomology", 2, list(tbl.dims), tbl.engine)
    assert blk is None
