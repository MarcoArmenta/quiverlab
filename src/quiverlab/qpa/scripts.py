"""GAP/QPA script builders. Translate a quiverlab Algebra's presentation (its
quiver + relations, over QQ or GF(p)) into QPA constructor calls, and assemble
the enveloping-algebra Hochschild route (HH^n = Ext^n_{A^e}(A,A)) since QPA ships
no HH function (spec §2). Scripts are GAP source strings; no floats.

Scope: cross-check runs on algebras presented over QQ or a prime field GF(p) --
the fields QPA supports exactly. Number-field CC entries and GF(p^n) are out of
the cross-check scope (raise QpaUnavailableError with that reason)."""
from __future__ import annotations


def _gap_field(domain) -> str:
    """QPA base field literal for a quiverlab Domain (QQ or GF(p))."""
    name = getattr(domain, "name", "")
    char = domain.characteristic       # frozen Domain: characteristic is an int attribute
    if char == 0 and name in ("QQ", "Rationals"):
        return "Rationals"
    if char > 0 and getattr(domain, "degree", 1) == 1:      # prime field GF(p)
        return f"GF({char})"
    raise ValueError(
        f"QPA cross-check supports QQ or prime GF(p) only; got domain {name!r} "
        f"(characteristic {char}). Number-field CC and GF(p^n) are out of scope."
    )


def quiver_and_algebra_script(algebra) -> str:
    """Emit GAP source binding `A := kQ/rels` (or `A := kQ` when no relations),
    reconstructing the quiver from algebra.quiver and its relations. Vertices are
    numbered 1..n in quiver order; arrows carry their quiverlab names. VERIFY the
    relation-string translation against QPA's element grammar at execution."""
    Q = algebra.quiver
    verts = list(Q.vertices)
    idx = {v: i + 1 for i, v in enumerate(verts)}               # QPA is 1-based
    arrows = [[idx[Q.source(a)], idx[Q.target(a)], a] for a in Q.arrows]
    arrow_gap = ", ".join(f'[{s}, {t}, "{name}"]' for s, t, name in arrows)
    field = _gap_field(algebra.domain)
    lines = [
        f"Q := Quiver({len(verts)}, [{arrow_gap}]);;",
        f"kQ := PathAlgebra({field}, Q);;",
    ]
    rels = algebra.relations
    if rels:
        # Each relation is a linear combo of parallel paths; render as a QPA element
        # over the generators kQ.<arrow>. VERIFY the exact element syntax at run.
        terms = _relations_to_gap(rels, "kQ")
        lines.append(f"rels := [{terms}];;")
        lines.append("A := kQ/rels;;")
    else:
        lines.append("A := kQ;;")
    return "\n".join(lines)


def _relations_to_gap(relations, kq: str) -> str:
    """Render quiverlab relations (tuples of (coeff, word)) as QPA algebra
    elements `sum coeff * kQ.a1*kQ.a2*...`. Coefficients are exact integers/
    fractions -> GAP integer/rational literals (no floats). VERIFY grammar."""
    out = []
    for rel in relations:
        parts = []
        for coeff, word in rel.terms:              # Relation.terms: ((Fraction, (arrow,...)),...)
            num, den = coeff.numerator, coeff.denominator
            path = "*".join(f"{kq}.{a}" for a in word)
            scal = f"{num}" if den == 1 else f"({num}/{den})"
            parts.append(f"{scal}*{path}")
        out.append(" + ".join(parts))
    return ", ".join(out)


def hochschild_dims_script(algebra, top: int) -> str:
    """Append the enveloping-algebra HH route to the algebra script, binding a GAP
    list `hh := [dim HH^0, ..., dim HH^top]`. HH^n(A) = Ext^n_{A^e}(A,A).

    Dim read: QPA's `ExtAlgebraGenerators(M, n)` returns a list whose FIRST component
    is the list of `dim Ext^i(M, M)` for `i = 0..n` (the standard QPA idiom for reading
    an Ext/HH dimension series). Here `M = AA` is `A` as an `A^e`-module, so
    `ExtAlgebraGenerators(AA, top)[1]` is exactly `[dim HH^0, ..., dim HH^top]`.

    VERIFY AT EXECUTION: (a) `AlgebraAsModuleOverEnvelopingAlgebra` -- the QPA op that
    presents A as a right A^e-module (QPA manual ch.6/8); if the name differs, build the
    bimodule explicitly from `EnvelopingAlgebra` + the regular representation. (b) that
    `ExtAlgebraGenerators(M, n)[1]` is the degreewise-dimension component (GAP is
    1-indexed: `[1]` is the first return value); the `[1,0,0]` fixture is the oracle."""
    base = quiver_and_algebra_script(algebra)
    return base + "\n" + "\n".join([
        "Ae := EnvelopingAlgebra(A);;",
        "AA := AlgebraAsModuleOverEnvelopingAlgebra(A);;    # VERIFY constructor name",
        f"info := ExtAlgebraGenerators(AA, {top});;",
        f"hh := info[1];;                                    # dims Ext^0..Ext^{top}",
    ])


def _gap_scalar(x, field: str) -> str:
    """A single exact matrix entry as a GAP field element (QQ rational or GF(p))."""
    from fractions import Fraction
    if field == "Rationals":
        fr = Fraction(x)
        return f"{fr.numerator}" if fr.denominator == 1 else f"({fr.numerator}/{fr.denominator})"
    # GF(p): x is an int mod p; k*One(GF(p)) is the field element k.1 (no Z(p) logs)
    return f"({int(x)}*One({field}))"


def _gap_matrix(mat, field: str) -> str:
    rows = ["[" + ", ".join(_gap_scalar(x, field) for x in row) + "]" for row in mat]
    return "[" + ", ".join(rows) + "]"


def module_decl(algebra, dimvec_list, arrow_matrices, var: str) -> str:
    """A single ``var := RightModuleOverPathAlgebra(A, dimvec, [[arrow, mat], ...]);;``
    line (assumes the quiver+algebra script `A := ...` is already in the session).
    `arrow_matrices` carries only nonzero arrows (QPA defaults omitted arrows to 0);
    matrices are in QPA's row convention (see modules/qpa_module.graded_form)."""
    field = _gap_field(algebra.domain)
    parts = [f'["{a}", {_gap_matrix(mat, field)}]' for a, mat in arrow_matrices.items()]
    arrows_gap = "[" + ", ".join(parts) + "]"
    return f"{var} := RightModuleOverPathAlgebra(A, {list(dimvec_list)}, {arrows_gap});;"


def ext_algebra_generators_script(algebra, top: int) -> str:
    """Bind ``info := ExtAlgebraGenerators(M, top)`` for ``M = (+) SimpleModules(A)``
    (Plan 27 Yoneda-algebra crosscheck). QPA's ``ExtAlgebraGenerators`` returns a
    three-entry list; the crosscheck reads back two of them, each on its OWN trailing
    statement (``session.run`` returns only the last statement's value):

    * ``info[1]`` = ``[dim Ext^n(M, M)]`` for n = 0..top. Since
      ``Ext^n((+)S_i, (+)S_j) = (+)_{i,j} Ext^n(S_i, S_j)``, this is exactly the TOTAL
      graded dimension ``dim E^n = sum_{i,j} dim Ext^n(S_i, S_j)`` of the Yoneda
      algebra E(A). (Degree 0: ``dim Ext^0(M, M) = dim Hom(M, M) = |Q_0|``.)
    * ``info[2]`` = ``[# minimal E-algebra generators in degree n]`` for n = 0..top,
      where degree 0 counts the ``|Q_0|`` vertex idempotents QPA treats as the
      semisimple base R = k^{Q_0}.

    ``M`` is built in quiver-vertex order (QPA's ``SimpleModules`` ordering)."""
    base = quiver_and_algebra_script(algebra)
    return base + "\n" + "\n".join([
        "M := DirectSumOfQPAModules(SimpleModules(A));;",
        f"info := ExtAlgebraGenerators(M, {top});;",
    ])


def ext_quiver_script(algebra) -> str:
    """Bind ``SS := SimpleModules(A)`` (QPA lists the simples in quiver-vertex order)
    so the caller can read the Ext^1 corner matrix pairwise. Pair each read with
    :func:`ext_quiver_entry` -- one ``ExtOverAlgebra`` eval per corner, honouring the
    one-statement-per-line rule of ``session.run``."""
    return quiver_and_algebra_script(algebra) + "\nSS := SimpleModules(A);;"


def ext_quiver_entry(i: int, j: int) -> str:
    """The single read-back statement giving ``dim Ext^1(S_i, S_j)`` for the corner
    (i, j) (1-based QPA indices, matching :func:`ext_quiver_script`'s ``SS``):
    ``dim Ext^1(S_i, S_j) = Length(ExtOverAlgebra(S_i, S_j)[2])`` (the second entry of
    ``ExtOverAlgebra`` is the basis of Ext^1)."""
    return f"Length(ExtOverAlgebra(SS[{i}], SS[{j}])[2]);"


def quadratic_ideal_script(algebra) -> str:
    """Bind ``isquad := IsQuadraticIdeal(rels)`` -- QPA's quadraticity test on the
    defining ideal I of A = kQ/I (I generated in degree 2). ``quiver_and_algebra_script``
    emits the ``rels`` list only when A HAS relations; a hereditary/semisimple A (I = 0,
    the ``A := kQ`` branch) gets an explicit empty ``rels := []`` -- ``IsQuadraticIdeal([])``
    is ``true`` (vacuously quadratic). Read back ``isquad;``."""
    base = quiver_and_algebra_script(algebra)
    lines = [base]
    if not algebra.relations:
        lines.append("rels := [];;")
    lines.append("isquad := IsQuadraticIdeal(rels);;")
    return "\n".join(lines)


def symmetric_predicates_script(algebra) -> str:
    """Bind ``A := kQ/rels`` (via :func:`quiver_and_algebra_script`) so the caller can
    read QPA's ``IsSymmetricAlgebra(A)`` and ``IsWeaklySymmetricAlgebra(A)`` each on its
    own trailing statement (``session.run`` returns only the last statement's value).
    Plan 29: crosschecks quiverlab's ``is_symmetric`` / ``is_weakly_symmetric`` against
    QPA for algebras that carry a quiver presentation (QQ or prime GF(p))."""
    return quiver_and_algebra_script(algebra)


def module_self_ext_dims_script(algebra, dimvec_M, top: int) -> str:
    """Bind `ext := [dim Ext^0(M,M), ..., dim Ext^top(M,M)]` (self-Ext of one module
    given by its dimension vector) via the SAME idiom `ExtAlgebraGenerators(M, top)[1]`.

    `dimvec_M` is a LIST of ints in `algebra.quiver.vertices` order (the same order
    the quiver builder numbers QPA's 1-based vertices) -- the caller flattens the
    Module.dimension_vector() dict.

    Self-Ext keeps the cross-check on the one confirmed QPA idiom. Distinct-module
    Ext(M,N) (M != N) needs `ExtOverAlgebra` + iterated `NthSyzygy` instead and is a
    flagged post-v1 extension."""
    base = quiver_and_algebra_script(algebra)
    return base + "\n" + "\n".join([
        f"M := RightModuleOverPathAlgebra(A, {list(dimvec_M)}, []);;",
        f"info := ExtAlgebraGenerators(M, {top});;",
        f"ext := info[1];;                                   # dims Ext^0..Ext^{top}(M,M)",
    ])
