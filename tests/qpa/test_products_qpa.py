"""QPA (GAP) probe for a Hochschild PRODUCT surface (Plan 35, Task 7).

Plan 35 exposes the cup product / cap action / Gerstenhaber bracket / Connes
differential on HH^*(A). This file asks the honest question: does QPA 1.37 ship
any surface we could crosscheck those products against?

It does NOT. QPA 1.37 has NO Hochschild-cohomology-ring surface: no CupProduct,
no HochschildCohomologyRing(*), no HochschildCohomology(*), no HochschildProduct,
and NamesGVars() contains zero names matching "Hochschild" or "Cup". What it DOES
expose -- ExtAlgebraGenerators and YonedaProduct -- is the MODULE Ext algebra
Ext^*_A(M, M) with the Yoneda product (the Plan-27 surface), over the algebra A
itself. That is a different object from HH^*(A) = Ext^*_{A^e}(A, A): repurposing
it for Hochschild would require constructing the enveloping algebra A^e as a QPA
algebra and A as an A^e-module, a construction QPA does not provide. So there is
no Hochschild PRODUCT surface to compare against.

Consequently this probe SKIPS with an honest message (the expected outcome named
in the Task-7 brief). The covering oracle for the product surface is the Task-6
identity battery (Gerstenhaber/cup-associativity/cap-Leibniz/Connes-B^2=0), plus
the literature pins in tests/hochschild/test_products_literature.py. The
verification page (Task 13) records this honest-scope entry.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a Hochschild PRODUCT surface.
_HH_PRODUCT_NAMES = (
    "CupProduct",
    "HochschildCohomologyRing",
    "HochschildCohomology",
    "HochschildCohomologyRingInDegrees",
    "HochschildCohomologyRingInfo",
    "HochschildProduct",
    "MultiplicativeStructureOfExtAlgebra",
)


def test_qpa_exposes_no_hochschild_product_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST -- for no name mentioning "Hochschild"
    #     or "Cup". This must precede any IsBoundGlobal query below: in GAP,
    #     IsBoundGlobal("Foo") REGISTERS "Foo" into NamesGVars() (as an unbound
    #     known name), so scanning after the queries would echo them back and
    #     falsely "find" a surface. Verified live 2026-08-01.
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    hochschild_like = sorted(n for n in gvar_names
                             if "hochschild" in n.lower() or "cup" in n.lower())

    # (2) None of the named Hochschild-product entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
             for name in _HH_PRODUCT_NAMES}
    present = [name for name, ok in bound.items() if ok]

    # (3) The Ext-algebra/Yoneda surface that DOES exist is module Ext
    #     (Ext^*_A(M,M)), not the Hochschild HH^*(A) cup product -- record that it
    #     is present so the skip is honest about what QPA offers instead.
    module_ext_surface = {
        name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
        for name in ("ExtAlgebraGenerators", "YonedaProduct")
    }

    if not present and not hochschild_like:
        pytest.skip(
            "QPA 1.37 exposes no Hochschild product surface "
            "(no CupProduct / HochschildCohomologyRing*; NamesGVars has no "
            "Hochschild/Cup name). Its ExtAlgebraGenerators/YonedaProduct is the "
            f"MODULE Ext algebra, not HH^*(A). module_ext_surface={module_ext_surface}. "
            "Covering oracle: the Task-6 product identity batteries + the "
            "literature pins in tests/hochschild/test_products_literature.py."
        )

    # If a future QPA ever grows a Hochschild product surface, FAIL loudly so this
    # honest-scope skip is revisited and a real crosscheck is wired in.
    pytest.fail(
        "QPA now exposes a Hochschild product surface -- wire a real crosscheck "
        f"against the Plan-35 cup/cap/bracket tables. Found: present={present}, "
        f"hochschild_like_names={hochschild_like}."
    )
