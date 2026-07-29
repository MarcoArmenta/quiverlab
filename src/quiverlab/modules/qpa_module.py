"""Serialize a quiverlab right-module into QPA's representation form (Plan 23,
item 4 oracle). QPA's ``RightModuleOverPathAlgebra(A, dimvec, [[arrow, mat], ...])``
wants, per arrow ``a: s -> t``, a ``dim_s x dim_t`` matrix acting on ROW vectors
(``m -> m * mat``). Our modules use the column/anti-homomorphism convention
(``m*b = action[b] @ m``, ``m`` a column), so the QPA matrix is the transpose of
the target-in-source-basis block.

The same graded form is the natural serialization for the no-code module schema
(Tier 1b item 3: ``{dims: {v: n_v}, maps: {arrow: [[...]]}}``); :func:`module_blocks`
below emits exactly that shape (Plan 34)."""
from fractions import Fraction

from quiverlab.modules import linalg_mod as lm


def graded_form(M):
    """(dimvec_list, arrow_matrices) in vertex order / QPA's row convention.

    dimvec_list[k] = dim (M e_{v_k}); arrow_matrices[a] = the dim_s x dim_t matrix
    (row convention) for each arrow a with a NONZERO action between positive-dim
    vertices (arrows omitted default to the zero map in QPA)."""
    A = M.algebra
    dom = A.domain
    verts = list(A.quiver.vertices)
    Ev, dv = {}, {}
    for v in verts:
        Pv = M.vertex_projection(v)
        piv = lm.column_space_pivots(Pv, dom)
        Ev[v] = [lm.col(Pv, j) for j in piv]        # a k-basis of M e_v
        dv[v] = len(piv)
    dimvec = [dv[v] for v in verts]
    arrows = {}
    for a, (s, t) in A.quiver.arrows.items():
        if dv[s] == 0 or dv[t] == 0:
            continue
        Es = lm.cols_to_matrix(Ev[s])               # M.dim x dim_s
        Et = lm.cols_to_matrix(Ev[t])               # M.dim x dim_t
        AEs = lm.matmul(M.action[a], Es, dom)       # a applied to each source basis vector
        Ca = lm.solve_columns(Et, AEs, dom)         # coords in the target basis (dim_t x dim_s)
        assert Ca is not None, (
            f"graded_form: arrow {a!r} image escapes the target vertex space "
            "(module is not vertex-graded)")
        Ca_mat = lm.cols_to_matrix(Ca)              # dim_t x dim_s
        if all(dom.is_zero(x) for row in Ca_mat for x in row):
            continue
        arrows[a] = lm.transpose(Ca_mat)            # dim_s x dim_t (QPA row convention)
    return dimvec, arrows


def _json_entry(x, dom=None):
    """A single action-matrix entry as an EXACT JSON scalar: an int when integer,
    else a fraction string like ``"1/2"`` -- exactly the entry grammar the no-code
    module INPUT accepts (schema.py::_valid_entry), so a serialized module round-trips
    straight back into the panel. Never a float.

    Exotic exact elements (e.g. GF(p^n) extension elements, which are internally
    little-endian coefficient TUPLES like ``(0, 1)``) render via the domain's OWN
    human-readable ``to_str`` -- e.g. ``"x^1"`` / ``"1 + x^1"`` -- a display-only
    STRING that is clearly not a plain int/fraction (Plan 34 MAJOR-4: the old
    ``str(x)`` fallback leaked the internal tuple ``"(0, 1)"`` into matrices). Such an
    entry is NOT re-enterable through the input grammar; :func:`module_blocks` flags
    the block ``display_only`` when one appears."""
    if isinstance(x, bool):
        return int(x)
    if isinstance(x, int):
        return int(x)
    if isinstance(x, Fraction):
        return int(x) if x.denominator == 1 else str(x)
    s = str(x)
    try:
        return int(s)
    except ValueError:
        try:
            f = Fraction(s)
        except (ValueError, ZeroDivisionError):
            # Not int/fraction-parseable: let the domain render itself readably
            # (GF(p^n) polynomial-in-generator notation) instead of leaking a tuple.
            if dom is not None:
                to_str = getattr(dom, "to_str", None)
                if callable(to_str):
                    return to_str(x)
            return s
        return int(f) if f.denominator == 1 else str(f)


def _reenterable(val) -> bool:
    """True iff a :func:`_json_entry` result (an int or a string) is re-enterable
    through the no-code module INPUT grammar -- an integer or an ``"a/b"`` fraction.
    The domain's human rendering of an extension-field element (``"x^1"``) is NOT."""
    if isinstance(val, int):
        return True
    try:
        Fraction(val)                            # "3", "-2", "1/2" all pass
        return True
    except (ValueError, ZeroDivisionError, TypeError):
        return False


def module_blocks(M):
    """Serialize module ``M`` as ``{"dims": {str(v): n}, "maps": {arrow: [[..]]}}`` in
    the no-code module INPUT schema (Plan 26/34).

    ``dims`` is the dimension vector (keyed by ``str(vertex)``, sorted). ``maps`` gives,
    per representation-quiver arrow ``a: s -> t`` between POSITIVE-dimensional vertices,
    the ``dim_t x dim_s`` action BLOCK in the COLUMN convention ``Algebra.module``
    consumes (``full[t_rows][s_cols] = block``) -- i.e. the target-in-source-basis
    block WITHOUT the QPA transpose. Zero blocks are kept (an explicit zero matrix at
    that arrow); arrows touching a 0-dimensional vertex are omitted (their block is the
    unique degenerate zero map and carries no information). Entries are exact ints /
    fraction strings, so a computed rad/top/soc feeds straight back into the panel.

    Deterministic (canonical pivots + a unique solve against an independent target
    basis), hence byte-reproducible -- the webapp and the Pyodide GUI both route their
    ``rad_top_soc`` payload through this one function and cannot drift.

    When any entry falls outside the int/fraction input grammar (e.g. a GF(p^n)
    extension-field element, rendered readably but not re-parseable), the block gains
    ``"display_only": true`` (Plan 34 MAJOR-4) so the GUI/webapp can say "display only
    -- not re-enterable" instead of failing when the user feeds it back. The key is
    ADDED ONLY when needed, so int/fraction modules (GF(p), CC-rational) stay
    byte-identical."""
    A = M.algebra
    dom = A.domain
    verts = list(A.quiver.vertices)
    Ev, dv = {}, {}
    for v in verts:
        Pv = M.vertex_projection(v)
        piv = lm.column_space_pivots(Pv, dom)
        Ev[v] = [lm.col(Pv, j) for j in piv]        # a k-basis of M e_v
        dv[v] = len(piv)
    dims = {str(v): dv[v] for v in sorted(verts, key=lambda w: str(w))}
    maps = {}
    display_only = False
    for a, (s, t) in A.quiver.arrows.items():
        if dv[s] == 0 or dv[t] == 0:
            continue
        Es = lm.cols_to_matrix(Ev[s])               # M.dim x dim_s
        Et = lm.cols_to_matrix(Ev[t])               # M.dim x dim_t
        AEs = lm.matmul(M.action[a], Es, dom)       # a applied to each source basis vector
        Ca = lm.solve_columns(Et, AEs, dom)         # coords in the target basis
        assert Ca is not None, (
            f"module_blocks: arrow {a!r} image escapes the target vertex space "
            "(module is not vertex-graded)")
        Ca_mat = lm.cols_to_matrix(Ca)              # dim_t x dim_s (column convention)
        block = []
        for row in Ca_mat:
            out_row = []
            for x in row:
                val = _json_entry(x, dom)
                if not _reenterable(val):
                    display_only = True
                out_row.append(val)
            block.append(out_row)
        maps[a] = block
    result = {"dims": dims, "maps": maps}
    if display_only:
        result["display_only"] = True
    return result


def summand_blocks(M, multiplicity=1):
    """One indecomposable summand of a Krull-Schmidt decomposition, as the report
    and the GUI show it (Marco 2026-07-29).

    A summand isomorphic to a STANDARD indecomposable is NAMED (``standard``:
    ``{"kind": "simple"|"projective"|"injective", "vertex": v}``) and its matrices
    are omitted -- ``S_2`` says everything the reader needs, and printing the action
    of a simple is noise. Every other summand carries its full ``maps`` (the exact
    per-arrow action, in the no-code INPUT schema), because there its dimension
    vector alone does not determine the module.

    Naming is best-effort and never guesses: when ``identify_standard`` cannot
    decide, the summand is simply shown in full."""
    from quiverlab.modules.hom import identify_standard
    out = {"dim_vector": {str(v): int(n) for v, n
                          in sorted(M.dimension_vector().items(), key=lambda kv: str(kv[0]))},
           "multiplicity": int(multiplicity), "indecomposable": True}
    std = identify_standard(M)
    if std is not None:
        out["standard"] = {"kind": std[0], "vertex": str(std[1])}
        return out
    blocks = module_blocks(M)
    out["maps"] = blocks.get("maps", {})
    if blocks.get("display_only"):
        out["display_only"] = True
    return out
