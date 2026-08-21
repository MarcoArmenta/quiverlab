"""Wave 5 (v1.0.1): ``A.enveloping()`` must not crash on GF(p^n) or CC.

Finding 5b: over GF(9) a raw ``TypeError`` (a GF(p^n) tuple coeff hit
``str(int(coeff))`` in the relation-string emitter); over CC a
``RelationError: malformed term ''``. The presentation path is gated with the
same ``_is_string_representable_domain`` check ``TrivialExtension`` uses, and off
that gate ``enveloping`` returns the structure-constant ``A (x) A^op`` build
(the certified iso target), so ``A^e`` computes over every field class.
"""
import pytest

import quiverlab
from quiverlab import fields
from quiverlab.families.extension import enveloping_algebra

selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_enveloping_gf_pn_no_raw_typeerror():
    B = quiverlab.linear_path_algebra(3, field=quiverlab.GF(9))
    Be = B.enveloping()
    assert Be.dim == B.dim ** 2 == 36     # A^e = A (x) A^op, dim (dim A)^2


@selfcert
def test_enveloping_cc_default_field():
    B = quiverlab.linear_path_algebra(3)  # CC, the DEFAULT field
    Be = B.enveloping()
    assert Be.dim == B.dim ** 2 == 36


@selfcert
def test_enveloping_gfp_and_qq_unchanged_bound_quiver():
    # DO NOT BREAK: over string-representable fields enveloping stays a bound-quiver
    # algebra with the same dim certificate.
    for dom in (quiverlab.GF(5), fields.QQ):
        B = quiverlab.linear_path_algebra(3, field=dom)
        Be = enveloping_algebra(B)
        assert Be.dim == B.dim ** 2 == 36
        assert Be.quiver is not None   # presented (bound-quiver) path preserved
