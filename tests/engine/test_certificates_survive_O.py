"""Load-bearing self-certifications must survive ``python -O``.

``python -O`` strips every ``assert`` statement.  The project's design documents
describe a family of internal certificates as "failing loudly" -- d.d = 0,
exactness, corner-typing, cup associativity, vertex-grading -- but when they are
spelled as a bare ``assert`` they VANISH under ``-O`` and a wrong mathematical
answer is returned silently (v1.0.1 wave 2, item 2b).  Every such certificate is
now an explicit ``raise`` (the style of ``resolutions_cs.resolution.assert_dd_zero``
and ``hochschild.koszul_ghms.assert_dd_zero``), which ``-O`` cannot strip.

Two gates: a live ``-O`` subprocess that trips one representative certificate and
requires it to still raise (carrying its message), and a source gate that no
load-bearing math file re-introduces a bare ``assert`` for these certificates.
"""
import ast
import os
import pathlib
import subprocess
import sys
import textwrap

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_SRC = _ROOT / "src"

# The load-bearing math files whose self-certifications were converted from bare
# ``assert`` to an explicit ``raise`` (v1.0.1 wave 2 item 2b).  None of these files
# may carry a bare ``assert``: a certificate spelled as ``assert`` dies under -O.
_NO_BARE_ASSERT = (
    "quiverlab/engine/hh_engine.py",
    "quiverlab/engine/complete_resolution.py",
    "quiverlab/engine/resolutions_bardzell.py",
    "quiverlab/engine/resolutions_minimal.py",
    "quiverlab/hochschild/bv/bracket.py",
    "quiverlab/modules/decompose.py",
    "quiverlab/modules/radtopsoc.py",
    "quiverlab/modules/ar.py",
    "quiverlab/modules/homdims.py",
    "quiverlab/modules/qpa_module.py",
    "quiverlab/derived/_corner.py",
    "quiverlab/resolutions_cs/comparison.py",
    "quiverlab/resolutions_cs/diagonal.py",
)


def test_no_bare_assert_in_load_bearing_math_files():
    offenders = []
    for rel in _NO_BARE_ASSERT:
        src = (_SRC / rel).read_text(encoding="utf-8")
        for node in ast.walk(ast.parse(src)):
            if isinstance(node, ast.Assert):
                offenders.append("%s:%d" % (rel, node.lineno))
    assert not offenders, (
        "bare `assert` in a load-bearing math file -- `python -O` strips it and the "
        "certificate silently vanishes; use an explicit `raise` instead: "
        + "; ".join(offenders))


def test_a_certificate_still_fires_under_dash_O():
    """Trip the AeEngine unit certificate (``hh_engine.Algebra`` needs the unit to
    carry a 1-coordinate for its f-basis change) inside a child running ``python -O``.

    Pre-fix (bare ``assert``): under -O the assert is stripped, so the constructor
    falls through to ``t = int(ts[0])`` and raises IndexError -- the certificate and
    its message are LOST.  Post-fix (explicit ``raise``): the certificate fires as an
    AssertionError carrying its message, even under -O."""
    script = textwrap.dedent(
        """
        import sys
        sys.path.insert(0, %r)
        if __debug__:                       # -O turns __debug__ off; this proves it
            print("NOT-UNDER-O")
            sys.exit(3)
        import numpy as np
        from quiverlab.engine.hh_engine import Algebra
        try:
            Algebra(2, np.zeros((2, 2, 2), dtype=np.int64),
                    np.zeros(2, dtype=np.int64))
        except AssertionError as exc:
            assert "unit must have a coordinate equal to 1" in str(exc), str(exc)
            print("CERT-FIRED")
            sys.exit(0)
        except Exception as exc:            # any OTHER exception = certificate lost
            print("WRONG-EXC:", type(exc).__name__, exc)
            sys.exit(1)
        print("NO-EXC")                     # -O silently accepted a bad unit
        sys.exit(2)
        """
    ) % str(_SRC)
    proc = subprocess.run([sys.executable, "-O", "-c", script],
                          capture_output=True, text=True, encoding="utf-8",
                          env=dict(os.environ, NUMBA_NUM_THREADS="2",
                                   OMP_NUM_THREADS="2"),
                          timeout=180)
    assert proc.returncode == 0, (proc.returncode, proc.stdout, proc.stderr)
    assert "CERT-FIRED" in proc.stdout, (proc.stdout, proc.stderr)
