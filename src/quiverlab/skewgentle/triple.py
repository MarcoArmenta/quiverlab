"""The skew-gentle triple ``(Q, I, Sp)``, the recognizer, and the associated gentle
algebra (Plan 68 / He-Zhou-Zhu 2004.11136 Def 1.1/1.3).

The triple is the user input.  Skew-gentle-ness is decidable: it reduces to gentleness
of the **associated gentle pair** ``(Q^sp, I^sp)`` -- ``Q`` plus one loop ``eps_i`` at
each special ``i``, and ``I^sp = I ∪ {eps_i^2}`` -- plus the structural loop/Sp
conditions.  ``associated_gentle`` materialises the genuinely-gentle algebra
``A^g = kQ^sp/<I ∪ {eps_i^2}>`` (NILPOTENT loops, admissible): the P46 string-
classification substrate and, by HZZ Lemma 1.5, ``dim A^g == dim(skew-gentle A)`` (the
dimension certificate the split constructor checks against).

Convention (HZZ Lemma 1.6): at most one loop per vertex, so a special vertex must not
already carry a loop in ``Q``.  Loop names are identifier-safe (``eps_<v>``), collision-
guarded, and frozen into the dataclass so ``split.py`` / ``modules.py`` read the same
names.  Float-free / field-agnostic (the recognizer is pure combinatorics)."""
from __future__ import annotations

import re
from dataclasses import dataclass

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.recognizers import is_gentle


def _loop_name(v) -> str:
    """An identifier-safe fresh loop name for the special vertex ``v``."""
    return "eps_" + re.sub(r"\W", "_", str(v))


def _relation_tokens(rel: str):
    """The arrow tokens of a relation string (``"a*b"`` -> ``["a", "b"]``).  A skew-
    gentle ``I`` is monomial length-2, so a valid relation yields exactly two tokens;
    a binomial ``"p - q"`` or a length-3 path yields more and is refused upstream."""
    return [t for t in re.split(r"[*+\-\s]+", str(rel)) if t]


@dataclass(frozen=True)
class SkewGentleTriple:
    """A validated skew-gentle triple.

    ``quiver`` is the ORIGINAL ``Q`` (no special loops); ``relations`` are the length-2
    monomial generators of ``I`` on ``Q``; ``special`` is ``Sp ⊆ Q_0``; ``loop_names``
    maps each special vertex to its fresh loop name ``eps_<v>``.
    """
    quiver: Quiver
    relations: tuple
    special: frozenset
    loop_names: dict

    # -- construction --------------------------------------------------------
    @classmethod
    def make(cls, quiver, relations=(), special=()) -> "SkewGentleTriple":
        """Build the triple: freeze the pieces and pick fresh, identifier-safe,
        collision-guarded loop names for each special vertex.  Does NOT validate --
        call ``.validate()`` (or use ``is_skew_gentle_triple``) for the loud checks."""
        special = frozenset(special)
        loop_names: dict = {}
        existing = set(quiver.arrows)
        for i in sorted(special, key=repr):
            base = _loop_name(i)
            nm, k = base, 0
            while nm in existing or nm in loop_names.values():
                k += 1
                nm = f"{base}_{k}"
            loop_names[i] = nm
        return cls(quiver=quiver, relations=tuple(relations), special=special,
                   loop_names=loop_names)

    # -- derived combinatorics ----------------------------------------------
    def q_sp(self) -> Quiver:
        """``Q^sp`` = ``Q`` plus one loop ``eps_i`` at each special vertex."""
        arrows = dict(self.quiver.arrows)
        for i in sorted(self.special, key=repr):
            arrows[self.loop_names[i]] = (i, i)
        return Quiver(list(self.quiver.vertices), arrows)

    def i_sp(self) -> tuple:
        """``I^sp`` = ``I`` plus the NILPOTENT squares ``eps_i*eps_i`` (the admissible
        associated-gentle form; NOT ``eps_i*eps_i - eps_i``, which is non-admissible)."""
        rels = list(self.relations)
        for i in sorted(self.special, key=repr):
            nm = self.loop_names[i]
            rels.append(f"{nm}*{nm}")
        return tuple(rels)

    # -- validation ----------------------------------------------------------
    def validate(self) -> None:
        """Loud on any structural or gentleness failure (``QuiverlabError`` with a
        ``hint``).  Checks, in order: ``Sp ⊆ Q_0``; no special vertex already carries a
        loop in ``Q`` (HZZ Lemma 1.6); every relation is a length-2 monomial path; the
        associated gentle pair ``(Q^sp, I^sp)`` is gentle."""
        Q = self.quiver
        vset = set(Q.vertices)
        for s in sorted(self.special, key=repr):
            if s not in vset:
                raise QuiverlabError(
                    f"SkewGentleTriple: special vertex {s!r} is not a vertex of Q "
                    f"(Q_0 = {sorted(vset, key=repr)})",
                    hint="Sp must be a subset of the quiver's vertices")
        for name, (src, tgt) in Q.arrows.items():
            if src == tgt and src in self.special:
                raise QuiverlabError(
                    f"SkewGentleTriple: special vertex {src!r} already carries a loop "
                    f"{name!r} in Q",
                    hint="a special vertex gets a FRESH loop eps; it may not already "
                         "have one (HZZ Lemma 1.6: at most one loop per vertex)")
        for rel in self.relations:
            toks = _relation_tokens(rel)
            if len(toks) != 2:
                raise QuiverlabError(
                    f"SkewGentleTriple: relation {rel!r} is not a length-2 monomial "
                    f"path (parsed {len(toks)} arrow token(s): {toks})",
                    hint="the ideal I of a skew-gentle triple is generated by length-2 "
                         "monomial paths a*b")
        # the decidable core: the associated gentle pair must be gentle.
        Ag = associated_gentle(self)
        if not is_gentle(Ag):
            raise QuiverlabError(
                "SkewGentleTriple: the associated gentle pair (Q^sp, I^sp) is NOT "
                "gentle, so (Q, I, Sp) is not a skew-gentle triple",
                hint="check the gentle conditions on Q^sp = Q + special loops with "
                     "I^sp = I ∪ {eps_i^2}: <=2 arrows in/out per vertex and the "
                     "branch conditions (HZZ Def 1.1/1.3)")


def associated_gentle(triple: SkewGentleTriple, field=None):
    """The associated gentle algebra ``A^g = kQ^sp/<I ∪ {eps_i^2}>`` (nilpotent loops).

    A genuine gentle algebra (admissible), the P46 string-classification substrate.
    ``dim A^g == dim(skew-gentle A)`` (HZZ Lemma 1.5) -- the dimension certificate the
    split constructor checks against.  Field-agnostic (gentleness is combinatorial); a
    ``field`` is only needed when the caller wants the dim over a specific domain."""
    Q_sp = triple.q_sp()
    I_sp = list(triple.i_sp())
    n_gen = len(list(Q_sp.vertices)) + len(Q_sp.arrows) + len(triple.special)
    bound = 2 * n_gen + 2
    return Q_sp.algebra(relations=I_sp, field=field, degree_bound=bound)


def is_skew_gentle_triple(quiver, relations=(), special=()) -> bool:
    """``True`` iff ``(quiver, relations, special)`` is a valid skew-gentle triple.

    The PRIMARY decidable recognizer: builds the triple, runs ``validate()`` inside a
    ``try``, and returns ``False`` on any ``QuiverlabError`` -- never raises (the
    boolean-recognizer contract shared with ``is_string`` / ``is_gentle``)."""
    try:
        SkewGentleTriple.make(quiver, relations, special).validate()
        return True
    except QuiverlabError:
        return False
