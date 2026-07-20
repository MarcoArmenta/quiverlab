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
        f"hh := info[1];;                                    # dims Ext^0..Ext^{top}  [VERIFY component]",
        "Print(hh);",
    ])


def module_self_ext_dims_script(algebra, dimvec_M, top: int) -> str:
    """Bind `ext := [dim Ext^0(M,M), ..., dim Ext^top(M,M)]` (self-Ext of one module
    given by its dimension vector) via the SAME idiom `ExtAlgebraGenerators(M, top)[1]`.

    Self-Ext keeps the cross-check on the one confirmed QPA idiom. Distinct-module
    Ext(M,N) (M != N) needs `ExtOverAlgebra` + iterated `NthSyzygy` instead and is a
    flagged post-v1 extension. VERIFY the `RightModuleOverPathAlgebra` args + the
    `[1]` component read at execution."""
    base = quiver_and_algebra_script(algebra)
    return base + "\n" + "\n".join([
        f"M := RightModuleOverPathAlgebra(A, {list(dimvec_M)}, []);;   # VERIFY args",
        f"info := ExtAlgebraGenerators(M, {top});;",
        f"ext := info[1];;                                   # dims Ext^0..Ext^{top}(M,M)  [VERIFY]",
        "Print(ext);",
    ])
