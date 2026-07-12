import pytest
import sympy
from quiverlab.errors import FieldError
from quiverlab.fields import GF
from quiverlab.fields.conway import CONWAY
from quiverlab.fields.finitefield import FiniteField, poly_is_irreducible


def test_dispatcher():
    assert GF(7).characteristic == 7          # PrimeField
    F4 = GF(4)
    assert isinstance(F4, FiniteField)
    assert F4.characteristic == 2 and F4.q == 4
    with pytest.raises(FieldError):
        GF(6)
    with pytest.raises(FieldError):
        GF(1)


def test_gf4_arithmetic():
    F4 = GF(4)
    x = F4.gen()
    x2 = F4.mul(x, x)
    # modulus x^2 + x + 1: x^2 = x + 1
    assert x2 == F4.add(x, F4.one())
    assert F4.mul(x, F4.inv(x)) == F4.one()
    # Frobenius sanity: a^4 = a for all a in GF(4)
    for a in F4.elements():
        assert F4.mul(F4.mul(a, a), F4.mul(a, a)) == a


def test_every_bundled_conway_entry_is_irreducible_and_primitive():
    for (p, n), coeffs in CONWAY.items():
        assert poly_is_irreducible(coeffs, p), f"CONWAY[{(p, n)}] not irreducible: replace from Lübeck's table"
        F = GF(p**n)
        x = F.gen()
        q = p**n
        # primitivity: order of x is exactly q - 1
        for ell in sympy.factorint(q - 1):
            assert F.pow(x, (q - 1) // ell) != F.one(), \
                f"CONWAY[{(p, n)}] generator not primitive: replace from Lübeck's table"
        assert F.pow(x, q - 1) == F.one()


def test_user_modulus_and_bad_modulus():
    # x^2 + 1 is irreducible over GF(3)
    F9 = GF(9, modulus=[1, 0, 1])
    assert F9.mul(F9.gen(), F9.gen()) == F9.neg(F9.one())
    with pytest.raises(FieldError):
        GF(9, modulus=[2, 0, 1])  # x^2 + 2 = (x+1)(x+2) mod 3: reducible


def test_beyond_table_needs_modulus():
    with pytest.raises(FieldError):
        GF(101**3)
