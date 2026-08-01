"""Plan 35 wave 3c -- the classical DICTIONARY: what each (co)homology space MEANS.

Marco: the report must be "complete with all interpretations of the spaces". Every
homological invariant quiverlab computes is a vector space whose ELEMENTS carry standard
representation-theoretic meaning -- Ext classes are exact sequences, HH^1 classes are
derivations, HH^2 classes are infinitesimal deformations, HH_0 is the commutator
quotient, and so on. This module is the SINGLE shared source of that interpretation
prose so the report HTML and both GUI renderers state it identically, never drifting.

Two honesty rules, enforced by construction:

  * CONSTRUCTED vs FRAMING. Where the capture layer records the actual class data (Ext
    Yoneda sequences in ``modules.complex_reps``; the per-degree cochain / cycle reps in
    the HH product + cyclic blocks) the interpretation is backed by explicit objects.
    Where only the dimension is known (the plain HH / cyclic dims tables carry no
    per-class reps yet), the sentence states the STANDARD meaning of the space and is
    plainly a framing sentence -- it never invents element-wise data.
  * Degrees beyond the elementary dictionary (HH^{>=3}, higher Ext / Tor) get the honest
    homological framing (Yoneda / obstruction / derived-functor), explicitly labelled as
    framing, not a fabricated concrete meaning.

Sentences are PLAIN unicode text (light math as ``HH^0``, ``M (x)_A N``, ``0 -> N -> E
-> M -> 0``), so the report's MathML surface and the GUIs' MathJax surface render the
SAME string identically -- the heavy math (the exact sequences, the matrices) goes
through each surface's proper math path, not through this prose. Float-free."""


# --------------------------------------------------------------------------- #
# Ext / Tor
# --------------------------------------------------------------------------- #
def ext_degree(n):
    """The meaning of a class of Ext^n_A(M, N)."""
    if n == 0:
        return ("Ext⁰(M, N) = Homₐ(M, N): its basis classes are the "
                "A-module homomorphisms M → N.")
    if n == 1:
        return ("Each basis class of Ext¹(M, N) is a short exact sequence "
                "0 → N → E → M → 0 (a Baer extension of M by N), up "
                "to equivalence. Below, each extension module E is constructed explicitly "
                "as a pushout and its exactness is verified.")
    return ("Each basis class of Extⁿ(M, N) is an n-fold exact sequence "
            "0 → N → Q → P_{n-2} → ⋯ → P_0 → M → 0 "
            "(Yoneda), up to equivalence. Below, each is spliced explicitly from the "
            "pushout module Q and the minimal resolution of M, and its exactness is "
            "verified at every joint.")


def tor_degree(n):
    """The meaning of a class of Tor_n^A(M, N)."""
    if n == 0:
        return ("Tor₀(M, N) = M ⊗ₐ N, the tensor product itself: the "
                "coequalizer of the two actions M ⊗ A ⊗ N ⇉ M ⊗ N. "
                "Its classes are the cosets m ⊗ n.")
    if n == 1:
        return ("Tor₁(M, N) measures the failure of M (equivalently N) to be flat: "
                "a nonzero class is a syzygy relation among the generators that "
                "− ⊗ N does not see, i.e. an obstruction to flatness.")
    return ("Torₙ(M, N) is the n-th derived functor of − ⊗ₐ N at M: "
            "a higher syzygy / flatness obstruction (homological framing).")


# --------------------------------------------------------------------------- #
# Hochschild cohomology HH^n(A) -- the Gerstenhaber dictionary.
# --------------------------------------------------------------------------- #
def hh_cohomology_degree(n):
    if n == 0:
        return ("HH⁰(A) = Z(A), the CENTRE of A: the classes are the central "
                "elements z (those with za = az for all a in A).")
    if n == 1:
        return ("HH¹(A) = Der(A) / Inn(A), the OUTER DERIVATIONS: each class is a "
                "derivation D : A → A (a k-linear map with the Leibniz rule "
                "D(ab) = D(a) b + a D(b)), determined by its values on the arrow "
                "generators, taken modulo the inner derivations a ↦ ax − xa.")
    if n == 2:
        return ("HH²(A) classifies the INFINITESIMAL DEFORMATIONS of A: each class "
                "is a 2-cocycle μ(a, b) giving a first-order (square-zero) "
                "deformation of the multiplication a * b = ab + t·μ(a, b); the "
                "coboundaries are the trivial (gauge) deformations.")
    return ("HHⁿ(A) controls the higher obstructions to deforming A (its Yoneda "
            "product with HH² carries the obstruction cocycles); homological framing.")


def hh_homology_degree(n):
    if n == 0:
        return ("HH₀(A) = A / [A, A], the COMMUTATOR QUOTIENT: the classes are the "
                "residues of A modulo the subspace spanned by the commutators ab − ba.")
    return ("HHₙ(A) is the n-th Hochschild homology, Tor^{A^e}_n(A, A) -- a "
            "derived-functor / cyclic-theory invariant (homological framing).")


def cyclic_degree(n):
    if n == 0:
        return ("HC₀(A) = A / [A, A]: the same space as HH₀, read as the "
                "TRACE FUNCTIONALS on A (a trace τ with τ(ab) = τ(ba) is "
                "exactly a linear form on A / [A, A]).")
    return ("HCₙ(A) is cyclic homology in degree n, the homology of Connes' (b, B) "
            "total complex -- it packages the S, B, I periodicity of the Hochschild "
            "theory (homological framing).")


# --------------------------------------------------------------------------- #
# Dispatch: (theory, degree) -> sentence. Used by the renderers.
# --------------------------------------------------------------------------- #
_THEORY = {
    "ext": ext_degree,
    "tor": tor_degree,
    "hh_cohomology": hh_cohomology_degree,
    "HH^": hh_cohomology_degree,
    "hh_homology": hh_homology_degree,
    "HH_": hh_homology_degree,
    "cyclic_homology": cyclic_degree,
    "HC_": cyclic_degree,
}


def sentence(theory, n):
    """The interpretation sentence for ``(theory, degree n)``, or ``None`` if the
    theory has no dictionary entry."""
    fn = _THEORY.get(theory)
    return fn(int(n)) if fn is not None else None


# --------------------------------------------------------------------------- #
# HH^1 read-off: a degree-1 cochain rep -> the derivation's value list D(arrow)=value.
# The captured term-sum of a degree-1 HH cocycle is exactly the list of pairs
# ``(word, value)`` with word a single arrow, so the derivation is READ OFF, not
# recomputed.
# --------------------------------------------------------------------------- #
def derivation_values(terms):
    """From a degree-1 HH^1 cochain's labelled term-sum ``[[coeff, word, value], ...]``
    (each ``word`` a single arrow generator), the derivation's action list
    ``D(word) = coeff * value``, grouped by arrow. Presentation only: it relabels the
    already-shipped terms as D(.) values, inventing nothing."""
    by_arrow = {}
    for coeff, word, value in terms:
        piece = value if str(coeff) == "1" else "%s %s" % (coeff, value)
        by_arrow.setdefault(str(word), []).append(piece)
    out = []
    for arrow in sorted(by_arrow):
        out.append("D(%s) = %s" % (arrow, " + ".join(by_arrow[arrow])))
    return out
