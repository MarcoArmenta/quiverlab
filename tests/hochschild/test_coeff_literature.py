"""Gentle HH^1-with-coefficients (arXiv:1811.02211, Chaparro-Schroll-Solotar) +
relative HH^1 (arXiv:2411.03080, Lindell-Rubio y Degrassi).

BLOCKED-until-transcribed: both papers state their explicit HH^1 values through
ribbon-graph / general radical-square-zero combinatorics, NOT a single clean
transcribable integer in the abstract. A direct PDF fetch returns compressed
bytes and the ar5iv HTML conversion of 2411.03080 errors out ("Conversion to HTML
had a Fatal error"), so a verbatim quiver+relations+coefficient+integer+equation
could not be transcribed at implementation time (attempted 2026-08-07 via WebFetch
on ar5iv). Per the plan these pins are xfail(strict=False) until a concrete example
is transcribed VERBATIM (never a fabricated number). The internal identity oracles
(tests/hochschild/test_coeff_identities.py) carry certification meanwhile; the
provenance (both citation keys) is real and surfaced regardless.

To flip a pin live: read the paper's PDF/TeX, fill the exact quiver/relations,
the coefficient bimodule, the integer value, and the equation number into the
docstring, then replace the pytest.fail with the real assertion and drop the xfail.
"""
import pytest
import quiverlab as ql  # noqa: F401  (used once a pin is transcribed)
from quiverlab.hochschild.coefficients import Bimodule  # noqa: F401

pytestmark = pytest.mark.oracle_literature


@pytest.mark.xfail(reason="BLOCKED: transcribe a 1811.02211 gentle HH^1-with-coefficients "
                          "worked example (quiver+relations+M+integer+eq#)", strict=False)
def test_gentle_hh1_with_coefficients():
    # A = ...  # the paper's gentle quiver/relations
    # M = ...  # the paper's coefficient bimodule
    # assert A.hochschild_cohomology(1, coefficients=M).dims[1] == <paper value>  # eq. #
    pytest.fail("gentle HH^1-with-coefficients value not yet transcribed from arXiv:1811.02211")


@pytest.mark.xfail(reason="BLOCKED: transcribe the 2411.03080 radical-square-zero "
                          "relative-HH^1 value (this is the IN-SCOPE vertex-relative "
                          "E=kQ_0 setting P52 ships)", strict=False)
def test_radical_square_zero_relative_hh1():
    # A = ...  # the paper's radical-square-zero quiver/relations
    # assert A.hochschild_cohomology(1, relative_to="vertices").dims[1] == <paper value>  # eq. #
    pytest.fail("radical-square-zero relative-HH^1 value not yet transcribed from arXiv:2411.03080")
