"""Pure M2 script builders. String generation only -- no subprocess here.

Translation notes (single-vertex kQ/I over GF(p) -> AssociativeAlgebras):
loops become free-algebra generators. A quiverlab presentation is stored as
parsed ``Relation`` objects (``combinat.relations.Relation``: tuples of
``(coeff, word)`` where ``word`` is a tuple of arrow names) -- NOT the raw
strings the user typed -- so relations are RECONSTRUCTED here into M2 element
syntax (which shares ``*`` and ``^``): a maximal run of the same generator of
length >= 3 collapses to ``a^k`` (short repeats read clearer un-collapsed;
both forms are byte-equivalent to M2's free-algebra parser). Multi-vertex
input is refused loudly: M2 has no vertex-idempotent notion (Plan-36 honest
scope). Reads the SAME accessors the QPA builders read
(``algebra.quiver``, ``algebra.relations``, ``algebra.domain.characteristic``).
"""
from __future__ import annotations

from fractions import Fraction

from quiverlab.errors import QuiverlabError

SENTINEL = "<<QL>>"


def _word_to_m2(word) -> str:
    """A path word (tuple of arrow names) as an M2 monomial, collapsing a
    maximal run of one generator of length >= 3 into ``a^k``."""
    parts = []
    i, n = 0, len(word)
    while i < n:
        j = i
        while j < n and word[j] == word[i]:
            j += 1
        run = j - i
        if run >= 3:
            parts.append(f"{word[i]}^{run}")
        else:
            parts.extend([word[i]] * run)
        i = j
    return "*".join(parts)


def _coeff_magnitude_str(coeff) -> str:
    """The (positive) magnitude of an exact rational coefficient as M2 source:
    an integer, or a parenthesised ``num/den`` (M2 reads it in the base field)."""
    fr = Fraction(coeff)
    fr = -fr if fr < 0 else fr
    return f"{fr.numerator}" if fr.denominator == 1 else f"({fr.numerator}/{fr.denominator})"


def _relation_to_m2(rel) -> str:
    """A parsed ``Relation`` as an M2 free-algebra element (the ideal generator).
    Mirrors ``Relation.__repr__`` term ordering, but renders words with
    ``_word_to_m2`` (caret powers) and emits no surrounding whitespace."""
    out = ""
    for coeff, word in rel.terms:
        fr = Fraction(coeff)
        path = _word_to_m2(word)
        term = path if abs(fr) == 1 else f"{_coeff_magnitude_str(fr)}*{path}"
        if out == "":
            out = ("-" + term) if fr < 0 else term
        else:
            out += ("-" if fr < 0 else "+") + term
    return out


def _single_vertex_data(A):
    quiver = A.quiver
    if quiver is None or len(quiver.vertices) != 1:
        n = 0 if quiver is None else len(quiver.vertices)
        raise QuiverlabError(
            "M2 oracle scope is single-vertex algebras only: "
            "AssociativeAlgebras has no quiver/idempotent type "
            f"(got {n} vertices).")
    p = A.domain.characteristic
    if not p:
        raise QuiverlabError("M2 oracle batteries run over GF(p); "
                             "characteristic 0 input not supported here.")
    gens = sorted(quiver.arrows)          # loop names, deterministic order
    rels = [_relation_to_m2(r) for r in A.relations]
    return p, gens, rels


def graded_dims_script(A, top: int) -> str:
    """M2 script printing ``<<QL>> n dim_k(B_n)`` for n = 0..top, where B is
    the degree-``top``-truncated nc Groebner quotient of the free algebra."""
    p, gens, rels = _single_vertex_data(A)
    gen_list = ",".join(gens)
    ideal = ",".join(rels) if rels else ""
    lines = [
        f"kk = ZZ/{p}",
        f"F = kk<|{gen_list}|>",
    ]
    if ideal:
        lines += [
            f"I = ideal {{{ideal}}}",
            f"Igb = NCGB(I, {2 * top + 2})",
            "B = F/I",
        ]
    else:
        lines += ["B = F"]
    lines += [
        f'for n from 0 to {top} do print("{SENTINEL} " | toString n | " " '
        "| toString numgens source ncBasis(n, B))",
    ]
    return "\n".join(lines) + "\n"


def commutative_ext_script(p: int, variables, relations, top: int) -> str:
    """M2 script printing ``<<QL>> n rank(P_n)`` for the minimal graded free
    resolution of the residue field of ZZ/p[variables]/(relations)."""
    vs = ",".join(variables)
    rels = ",".join(r.replace(" ", "") for r in relations)
    return "\n".join([
        "needsPackage \"Complexes\"",
        f"R = ZZ/{p}[{vs}]/({rels})",
        f"C = freeResolution(coker vars R, LengthLimit => {top})",
        f'for n from 0 to {top} do print("{SENTINEL} " | toString n | " " '
        "| toString rank C_n)",
    ]) + "\n"


def parse_sentinels(stdout: str) -> list:
    """Extract the ``<<QL>> n v`` values in degree order; exact ints only."""
    got = {}
    for line in stdout.splitlines():
        line = line.strip()
        if not line.startswith(SENTINEL):
            continue
        _, n_s, v_s = line.split(None, 2)
        got[int(n_s)] = int(v_s)          # int() raises on "1.5" -- wanted
    if sorted(got) != list(range(len(got))):
        raise ValueError(f"non-contiguous degrees in M2 output: {sorted(got)}")
    return [got[n] for n in range(len(got))]
