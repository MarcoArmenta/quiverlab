"""QQi: the Gaussian rationals Q(i), exactly (a first-class public field-spec).

Unlike CC (whose working domain is whatever exact subfield of C its entries
generate), QQi is *pinned* to Q(i) = sympy's QQ_I regardless of the entries.
This gives an explicit, named "exact Q(i)" field while reusing CC's proven exact
machinery: parsing goes through CC.parse_entry (accepts 'i', '1+2*i', rationals;
rejects floats/complex loudly), and arithmetic runs in SympyExactDomain(QQ_I) --
the same exact Gaussian-rationals domain CC hands back when 'i' is among the
entries. Everything is sympy.I / sympy.Rational: never a Python complex, 1j, or
float literal (the loud-exactness contract, spec D3/§4.2)."""
import sympy
from sympy.polys.constructor import construct_domain

from quiverlab.errors import FieldError
from quiverlab.fields.complexfield import CC, SympyExactDomain

# The Gaussian-rationals domain QQ_I, obtained exactly the way CC obtains it when
# 'i' appears among entries. Pinned once: Q(i) is a fixed field, not entry-derived.
_QI_DOMAIN, _ = construct_domain([sympy.I], field=True, extension=True)


class GaussianRationalField:
    """Field spec pinned to Q(i). Like ComplexField, it is not itself a Domain:
    make_domain hands back the exact working Domain (SympyExactDomain(QQ_I))."""

    name = "QQ(i)"
    characteristic = 0

    def __repr__(self):
        return self.name

    def parse_entry(self, x):
        # Reuse CC's exact-complex parser: 'i', '1+2*i', ints, Fractions, exact
        # sympy Exprs pass; floats/complex raise ExactnessError. Out-of-field
        # values (e.g. sqrt(2)) parse here but are rejected at coerce time,
        # keeping the working field exactly Q(i) rather than a broader closure.
        return CC.parse_entry(x)

    def make_domain(self, entries):
        for e in entries:
            self.parse_entry(e)          # loud exactness gate; values are pinned-out
        if not (_QI_DOMAIN.is_Field and getattr(_QI_DOMAIN, "is_Numerical", False)
                and _QI_DOMAIN.is_Exact):
            raise FieldError(f"pinned Q(i) domain {_QI_DOMAIN} is not an exact numerical field")
        dom = SympyExactDomain(_QI_DOMAIN)
        dom.name = f"QQ(i) (computing exactly in {_QI_DOMAIN})"
        return dom


QQi = GaussianRationalField()
