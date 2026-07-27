"""Structural enforcement of the loud-exactness contract (spec D3, §4.1):
no float/complex literals and no float() calls anywhere under src/.

Beyond bare literals/`float(` calls, the gate also flags the sneakier ways an
IEEE-754 value can enter exact code: `.astype(float)` / `.astype(np.float64)`
dtype casts (positional), `dtype=float` / `dtype='float64'` / `dtype=np.float64`
KEYWORD casts on any call, dtype-LESS numpy array constructors
(`np.zeros/ones/empty/full/eye(...)` with no dtype -- these default to float64),
numpy float/complex linear algebra reached either as `np.linalg.*` OR via a
`from numpy.linalg import inv` / `from numpy import linalg` /
`import numpy.linalg as la` binding, numpy float/complex dtype references
(`np.float64`, `np.floating`, ...), `Decimal` (rounds by design), and
`sympify('0.5')` (a float smuggled in as a string).
Legitimate integer numpy ops (`np.zeros(n, dtype=np.int64)`, `np.eye(m, np.int64)`,
`.astype(np.int64)`, `np.outer`) and integer-returning `math.floor`/`math.ceil`
in exempt glue are deliberately NOT flagged."""
import ast
import pathlib
import re

SRC = pathlib.Path(__file__).resolve().parent.parent / "src" / "quiverlab"

# attr/name prefixes denoting an inexact (float/complex) dtype spelling
_INEXACT_DTYPE = ("float", "complex")
# a decimal-point float literal inside a string (0.5, 1., .5, 3.14e-2) -- but not
# '1/2', 'E(3)', 'sqrt(2)', or a symbol.method access with no digit by the dot
_DECIMAL_RE = re.compile(r"\d+\.\d*|\.\d+")


def _is_inexact_dtype(node) -> bool:
    """True for the `float`/`complex` builtins, an `np.float64`-style attribute,
    or a 'float64'/'complex128' dtype string -- any inexact dtype spelling."""
    if isinstance(node, ast.Name):
        return node.id in ("float", "complex")
    if isinstance(node, ast.Attribute):
        return node.attr.startswith(_INEXACT_DTYPE)
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value.startswith(_INEXACT_DTYPE)
    return False


# numpy array constructors that silently default to float64 when dtype is omitted
_ARRAY_CTORS = ("zeros", "ones", "empty", "full", "eye")


def _names_a_dtype(node) -> bool:
    """True if `node` spells a dtype at all -- exact OR inexact: a builtin
    (`int`/`float`/`complex`/`bool`/`object`), an `np.<dtype>` attribute
    (`np.int64`, `np.float64`, ...), or a dtype string (`'int64'`, `'float64'`).
    Used to decide whether an array constructor got an explicit dtype positionally
    (so `np.eye(m, np.int64)` is NOT mistaken for a dtype-less constructor)."""
    if isinstance(node, ast.Name):
        return node.id in ("int", "uint", "float", "complex", "bool", "object")
    if isinstance(node, ast.Attribute):
        return node.attr.startswith(("int", "uint", "float", "complex", "bool"))
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value.startswith(("int", "uint", "float", "complex", "bool"))
    return False


def _linalg_bound_names(tree) -> set:
    """Names bound to numpy.linalg (or a function pulled out of it) in this module:
    `from numpy.linalg import inv` -> {'inv'}, `from numpy import linalg` ->
    {'linalg'}, `import numpy.linalg as la` -> {'la'}. Each reopens the exact door
    the `np.linalg.*` attribute rule closes; their USE is flagged below. Restricted
    to the `numpy`/`numpy.linalg` modules, so `from quiverlab.fields import linalg`
    (our EXACT linalg) is never captured."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.module == "numpy.linalg":
                for a in node.names:
                    names.add(a.asname or a.name)
            elif node.module == "numpy":
                for a in node.names:
                    if a.name == "linalg":
                        names.add(a.asname or a.name)
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name == "numpy.linalg" and a.asname:
                    names.add(a.asname)
    return names


def _is_sympify(func) -> bool:
    return (
        (isinstance(func, ast.Name) and func.id == "sympify")
        or (isinstance(func, ast.Attribute) and func.attr == "sympify")
    )


def _violations(path: pathlib.Path) -> list[str]:
    # explicit encoding: src/ is UTF-8; the locale default (cp1252 on Windows CI)
    # chokes on the non-ASCII bytes in ported-source comments
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    linalg_names = _linalg_bound_names(tree)
    out = []
    for node in ast.walk(tree):
        # (1) float / complex literals
        if isinstance(node, ast.Constant) and isinstance(node.value, (float, complex)):
            out.append(f"{path.name}:{node.lineno}: literal {node.value!r}")
        # (2) bare float(...) call
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "float"
        ):
            out.append(f"{path.name}:{node.lineno}: float() call")
        # (3) .astype(float) / .astype(np.float64) / .astype('float64') (+complex)
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "astype"
            and node.args
            and _is_inexact_dtype(node.args[0])
        ):
            out.append(f"{path.name}:{node.lineno}: .astype(<inexact dtype>)")
        # (4) np.linalg.* -- float linear algebra (use exact fields.linalg instead)
        if (
            isinstance(node, ast.Attribute)
            and node.attr == "linalg"
            and isinstance(node.value, ast.Name)
            and node.value.id in ("np", "numpy")
        ):
            out.append(f"{path.name}:{node.lineno}: np.linalg.* (use exact fields.linalg)")
        # (5) np.float* / numpy.float* / np.complex* dtype reference
        if (
            isinstance(node, ast.Attribute)
            and node.attr.startswith(_INEXACT_DTYPE)
            and isinstance(node.value, ast.Name)
            and node.value.id in ("np", "numpy")
        ):
            out.append(f"{path.name}:{node.lineno}: {node.value.id}.{node.attr} (inexact dtype)")
        # (6) Decimal (rounds by design -- quiverlab is sympy/int exact)
        if (isinstance(node, ast.Name) and node.id == "Decimal") or (
            isinstance(node, ast.Attribute) and node.attr == "Decimal"
        ):
            out.append(f"{path.name}:{node.lineno}: Decimal")
        # (7) sympify('<decimal>') -- a float smuggled in as a string literal
        if (
            isinstance(node, ast.Call)
            and _is_sympify(node.func)
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
            and _DECIMAL_RE.search(node.args[0].value)
        ):
            out.append(
                f"{path.name}:{node.lineno}: sympify() of a decimal string "
                f"{node.args[0].value!r}"
            )
        # (8) dtype=<inexact> KEYWORD arg on any call -- e.g. np.zeros(3, dtype=float),
        #     np.array([1], dtype='float64'), np.full(s, 0, dtype=np.float64). The
        #     positional .astype(...) rule (3) never sees these.
        if isinstance(node, ast.Call):
            for kw in node.keywords:
                if kw.arg == "dtype" and _is_inexact_dtype(kw.value):
                    out.append(
                        f"{path.name}:{node.lineno}: dtype=<inexact dtype> keyword"
                    )
        # (9) dtype-LESS numpy array constructor (defaults to float64) -- or one
        #     given an inexact positional dtype. np.zeros/ones/empty/full/eye with an
        #     explicit dtype (kwarg OR an int-typed positional like np.eye(m, np.int64))
        #     pass; the codebase writes dtype=np.int64 everywhere, so this forces it.
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in _ARRAY_CTORS
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id in ("np", "numpy")
        ):
            ctor = node.func.attr
            if not any(kw.arg == "dtype" for kw in node.keywords):
                pos_dtypes = [a for a in node.args[1:] if _names_a_dtype(a)]
                if not pos_dtypes:
                    out.append(
                        f"{path.name}:{node.lineno}: np.{ctor}(...) with no explicit "
                        f"dtype (defaults to float64; write dtype=np.int64)"
                    )
                elif any(_is_inexact_dtype(a) for a in pos_dtypes):
                    out.append(
                        f"{path.name}:{node.lineno}: np.{ctor}(<inexact positional dtype>)"
                    )
        # (10) use of a name bound from numpy.linalg (from numpy.linalg import inv;
        #      from numpy import linalg; import numpy.linalg as la)
        if (
            isinstance(node, ast.Name)
            and isinstance(node.ctx, ast.Load)
            and node.id in linalg_names
        ):
            out.append(
                f"{path.name}:{node.lineno}: numpy.linalg via imported name "
                f"{node.id!r} (use exact fields.linalg)"
            )
    return out


def test_package_importable():
    import quiverlab
    assert isinstance(quiverlab.__version__, str)


def test_no_float_literals_or_calls_in_src():
    assert SRC.is_dir(), "src/quiverlab missing"
    bad = [v for f in SRC.rglob("*.py") for v in _violations(f)]
    assert bad == [], "floats are banned in quiverlab core:\n" + "\n".join(bad)


def test_gate_detects_planted_violations(tmp_path):
    bad = tmp_path / "planted.py"
    bad.write_text("x = 0.5\ny = float('2')\nz = 1j\n")
    found = _violations(bad)
    assert len(found) == 3


def test_gate_detects_extended_violations(tmp_path):
    """Each sneaky float-entry pattern the extended gate now covers."""
    bad = tmp_path / "planted_ext.py"
    bad.write_text(
        "import numpy as np\n"
        "import sympy\n"
        "from decimal import Decimal\n"
        "a = M.astype(float)\n"                 # .astype(float)
        "b = M.astype(np.float64)\n"            # .astype(np.float64)
        "c = M.astype('float64')\n"            # .astype('float64')
        "d = np.linalg.inv(M)\n"               # np.linalg.* (+ np.float? no)
        "e = np.float64(1)\n"                  # np.float* dtype ref
        "f = Decimal('1')\n"                   # Decimal
        "g = sympy.sympify('0.5')\n"           # sympify decimal string
        "h = sympify('1.5*x')\n"               # sympify decimal inside expr
    )
    found = _violations(bad)
    joined = "\n".join(found)
    assert any(".astype(<inexact dtype>)" in v for v in found), joined
    assert sum(".astype(<inexact dtype>)" in v for v in found) == 3, joined
    assert any("np.linalg" in v for v in found), joined
    assert any("np.float64 (inexact dtype)" in v for v in found), joined
    assert any(v.endswith("Decimal") for v in found), joined
    assert sum("sympify() of a decimal string" in v for v in found) == 2, joined


def test_gate_detects_dtype_and_linalg_leaks(tmp_path):
    """The critic's three recurrences of the leak the extended gate was built for:
    dtype=<inexact> KEYWORD casts, dtype-less float-defaulting constructors, and
    the `from numpy.linalg import inv` / `from numpy import linalg` re-imports."""
    bad = tmp_path / "planted_kw.py"
    bad.write_text(
        "import numpy as np\n"
        "from numpy.linalg import inv, solve\n"       # binds inv, solve
        "from numpy import linalg\n"                    # binds linalg
        "import numpy.linalg as la\n"                   # binds la
        "a = np.zeros(3, dtype=float)\n"               # dtype=float kwarg
        "b = np.array([1], dtype='float64')\n"        # dtype='float64' kwarg
        "c = np.full((2, 2), 0, dtype=np.float64)\n"  # dtype=np.float64 kwarg
        "d = np.zeros(m)\n"                             # bare -> float64
        "e = np.ones((3, 3))\n"                        # bare -> float64
        "f = np.eye(4)\n"                              # bare -> float64
        "g = np.empty(n)\n"                            # bare -> float64
        "h = np.zeros(3, float)\n"                     # inexact POSITIONAL dtype
        "p = inv(M)\n"                                  # numpy.linalg.inv via import
        "q = solve(M, b)\n"                            # numpy.linalg.solve via import
        "r = linalg.inv(M)\n"                          # from numpy import linalg
        "s = la.inv(M)\n"                              # import numpy.linalg as la
    )
    found = _violations(bad)
    joined = "\n".join(found)
    # dtype=<inexact> keyword: the three kwarg casts (a, b, c)
    assert sum("dtype=<inexact dtype> keyword" in v for v in found) == 3, joined
    # bare constructors d, e, f, g -> four "no explicit dtype" flags
    assert sum("with no explicit dtype (defaults to float64" in v for v in found) == 4, joined
    # the inexact positional dtype (h)
    assert any("inexact positional dtype" in v for v in found), joined
    # every numpy.linalg re-import path (p, q, r, s) is caught at the use site
    assert any("numpy.linalg via imported name 'inv'" in v for v in found), joined
    assert any("numpy.linalg via imported name 'solve'" in v for v in found), joined
    assert any("numpy.linalg via imported name 'linalg'" in v for v in found), joined
    assert any("numpy.linalg via imported name 'la'" in v for v in found), joined


def test_gate_allows_legitimate_integer_ops(tmp_path):
    """Exact/integer numpy + math patterns that must NOT be flagged."""
    ok = tmp_path / "clean.py"
    ok.write_text(
        "import numpy as np\n"
        "import math\n"
        "import sympy\n"
        "from sympy import sympify\n"
        "from quiverlab.fields import linalg\n"        # our EXACT linalg, not numpy's
        "a = np.zeros(3, dtype=np.int64)\n"
        "b = M.astype(np.int64)\n"
        "c = M.astype(int)\n"
        "d = np.outer(u, e_t)\n"
        "e = np.eye(3, dtype=np.int64)\n"
        "e2 = np.eye(m, np.int64)\n"                    # explicit int dtype, positional
        "e3 = np.zeros(3, np.int64)\n"                  # explicit int dtype, positional
        "e4 = np.full((2, 2), 0, dtype=np.int64)\n"    # explicit int dtype, kwarg
        "f = math.floor(x) - math.ceil(y)\n"          # integer-returning glue
        "g = np.flatnonzero(v)\n"                       # not 'float*'
        "h = sympify('E(3)')\n"                         # roots of unity, no decimal
        "i = sympify('1/2 + sqrt(2)')\n"               # rationals/radicals, no decimal
        "j = sympy.floor(z)\n"                          # attr 'floor', not 'float*'
        "k = linalg.solve(A, b)\n"                      # quiverlab.fields.linalg use
        "l = np.zeros_like(w)\n"                        # infers dtype, not a bare ctor
    )
    assert _violations(ok) == []
