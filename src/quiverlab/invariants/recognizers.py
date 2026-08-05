"""Algebra-class recognizers: is_semisimple .. is_gentle (Plan 38 / C2).

Each takes an Algebra A and returns a bool. Recognizers that need the quiver
presentation (everything except the two Loewy-length ones) refuse loudly on a
structure-constant algebra, mirroring modules.builders._require_provenance.

The special-biserial / string / gentle family turns on ideal membership of
length-2 paths, decided EXACTLY through the Groebner reduction system
(resolutions_cs.build.reduction_system_of): a length-2 path a*b lies in I iff
its normal form is zero. References: Assem-Simson-Skowronski (ASS 2006),
Butler-Ringel (string modules), Assem-Skowronski (gentle algebras)."""
from quiverlab.errors import QuiverlabError


def _require_quiver(A, what):
    if A.quiver is None or A.relations is None:
        raise QuiverlabError(
            f"{what} needs the quiver presentation",
            hint="build the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no quiver")


# -- Loewy-length recognizers (any presentation) ------------------------------
def is_semisimple(A) -> bool:
    """True iff A is semisimple: Loewy length 1 (rad A = 0)."""
    return A.loewy_length() == 1


def is_radical_square_zero(A) -> bool:
    """True iff rad^2 A = 0: Loewy length <= 2."""
    return A.loewy_length() <= 2


# -- structural recognizers (need the quiver) ---------------------------------
def is_hereditary(A) -> bool:
    """True iff A is hereditary (gl.dim <= 1): a path algebra kQ with Q acyclic
    and no relations. An admissible ideal I != 0 forces gl.dim >= 2, so any
    relation makes A non-hereditary."""
    _require_quiver(A, "is_hereditary")
    return (not A.relations) and A.quiver.is_acyclic()


def is_basic(A) -> bool:
    """True for every kQ/I presentation: the vertex idempotents e_v are a
    complete set of primitive orthogonal idempotents whose indecomposable
    projectives e_v A are pairwise non-isomorphic, so A is basic."""
    _require_quiver(A, "is_basic")
    return True


def _in_degree(Q):
    d = {v: 0 for v in Q.vertices}
    for (_s, t) in Q.arrows.values():
        d[t] += 1
    return d


def _out_degree(Q):
    d = {v: 0 for v in Q.vertices}
    for (s, _t) in Q.arrows.values():
        d[s] += 1
    return d


def is_nakayama(A) -> bool:
    """True iff the underlying quiver is a disjoint union of linear A_n paths and
    single oriented cycles: every vertex has in-degree <= 1 and out-degree <= 1
    (each connected component is then uniserial)."""
    _require_quiver(A, "is_nakayama")
    Q = A.quiver
    ind, outd = _in_degree(Q), _out_degree(Q)
    return all(ind[v] <= 1 and outd[v] <= 1 for v in Q.vertices)


def _len2_in_ideal(rs, a, b) -> bool:
    """True iff the length-2 path a*b (first a then b) lies in I: its normal form
    under the reduction system is zero. Assumes a*b is composable
    (target(a) = source(b)) -- the callers only pass composable pairs."""
    nf = rs.normal_form((a, b))
    return all(rs.domain.is_zero(c) for c in nf.values())


def _is_monomial_ideal(A) -> bool:
    """Ideal-invariant monomial test (devil's-advocate fix 2026-08-05): read
    the Groebner reduction system, not the raw generator list -- a rule with
    an empty tail is a monomial rewrite, and the reduced system is canonical
    for the ideal (a non-minimal user presentation must not change verdicts)."""
    from quiverlab.resolutions_cs.build import reduction_system_of
    return all(not r.tail for r in reduction_system_of(A).rules)


def is_special_biserial(A) -> bool:
    """Special biserial (ASS): (1) every vertex has <= 2 arrows in and <= 2 out;
    (2) for every arrow b, at most one arrow c with b*c not in I and at most one
    arrow a with a*b not in I."""
    _require_quiver(A, "is_special_biserial")
    Q = A.quiver
    ind, outd = _in_degree(Q), _out_degree(Q)
    if any(ind[v] > 2 or outd[v] > 2 for v in Q.vertices):
        return False
    from quiverlab.resolutions_cs.build import reduction_system_of
    rs = reduction_system_of(A)
    arrows = list(Q.arrows.items())                 # (name, (source, target))
    for b, (bs, bt) in arrows:
        after = [c for c, (cs, _ct) in arrows
                 if cs == bt and not _len2_in_ideal(rs, b, c)]
        if len(after) > 1:
            return False
        before = [a for a, (_asrc, at) in arrows
                  if at == bs and not _len2_in_ideal(rs, a, b)]
        if len(before) > 1:
            return False
    return True


def is_string(A) -> bool:
    """String algebra: special biserial with a MONOMIAL ideal I."""
    _require_quiver(A, "is_string")
    return is_special_biserial(A) and _is_monomial_ideal(A)


_RECOGNIZER_FLAGS = [
    "is_semisimple", "is_radical_square_zero", "is_hereditary", "is_basic",
    "is_nakayama", "is_special_biserial", "is_string", "is_gentle",
    "is_selfinjective", "is_symmetric",
]


def recognizers_block(A):
    """The no-code ``recognizers`` compute block (Plan 38): the eight structural
    recognizers plus is_selfinjective / is_symmetric, the orientation-blind
    Dynkin/Euclidean type (as a compact string, e.g. "A5" / "~D6"), and the
    finite/tame/wild form type. A flag that REFUSES on this input (e.g. a
    presentation-less algebra) is reported per-flag as ``{"error": msg}``, never
    a silent False and never a crash (the Plan-30 tau-block precedent). SHARED by
    both runners so their blocks are byte-identical; each runner adds
    ``citations`` from ``references``."""
    fns = {
        "is_semisimple": is_semisimple,
        "is_radical_square_zero": is_radical_square_zero,
        "is_hereditary": is_hereditary,
        "is_basic": is_basic,
        "is_nakayama": is_nakayama,
        "is_special_biserial": is_special_biserial,
        "is_string": is_string,
        "is_gentle": is_gentle,
        "is_selfinjective": lambda X: X.is_selfinjective(),
        "is_symmetric": lambda X: X.is_symmetric(),
    }
    flags = {}
    for name in _RECOGNIZER_FLAGS:
        try:
            flags[name] = bool(fns[name](A))
        except Exception as exc:            # honest per-flag error, never a 500
            flags[name] = {"error": str(exc)}
    dynkin = None
    try:
        from quiverlab.invariants.dynkin_type import dynkin_type
        if A.quiver is not None:
            dt = dynkin_type(A.quiver)
            dynkin = None if dt is None else "%s%s" % (dt[0], dt[1])
    except Exception:
        dynkin = None
    try:
        from quiverlab.invariants.forms import form_type
        ftype = form_type(A)
    except Exception:
        ftype = None
    return {"flags": flags, "dynkin_type": dynkin, "form_type": ftype,
            "references": ["assem_book"]}


def is_gentle(A) -> bool:
    """Gentle algebra: a string algebra whose ideal is generated by length-2
    paths, plus the dual biserial condition -- for every arrow b, at most one
    arrow c with b*c IN I and at most one arrow a with a*b IN I."""
    _require_quiver(A, "is_gentle")
    if not is_string(A):
        return False
    Q = A.quiver
    from quiverlab.resolutions_cs.build import reduction_system_of
    rs = reduction_system_of(A)
    # ideal-invariant length-2 test (devil's-advocate fix 2026-08-05): for a
    # monomial ideal the reduced system's leads are its unique minimal
    # monomial generators, so "generated by length-2 paths" == every lead
    # has length 2 -- regardless of how the user wrote the relations.
    if not all(len(r.lead) == 2 for r in rs.rules):
        return False
    arrows = list(Q.arrows.items())
    for b, (bs, bt) in arrows:
        after_in = [c for c, (cs, _ct) in arrows
                    if cs == bt and _len2_in_ideal(rs, b, c)]
        if len(after_in) > 1:
            return False
        before_in = [a for a, (_asrc, at) in arrows
                     if at == bs and _len2_in_ideal(rs, a, b)]
        if len(before_in) > 1:
            return False
    return True
