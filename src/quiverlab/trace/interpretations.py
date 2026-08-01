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
# Typing statements (Marco 2026-07-31): at the top of every (co)homology section,
# state EXACTLY what object the engine computes and what the tensor/bar-bracket
# notation means -- true to the code (engine.tt_calculus / hochschild.bar cochains
# C^n = Hom_k(Ā^⊗n, A), chains C_n = A ⊗_k Ā^⊗n; resolutions_cs P_n = ⊕ A e_o ⊗ e_t A).
# ``route`` is "bar" (bar / fast / minimal cochain bases) or "cs" (Chouhy-Solotar).
# Plain unicode so the report (MathML surface) and the GUI (MathJax) render identically.
# --------------------------------------------------------------------------- #
_NOTATION_TAIL = (" Here | and ⊗ are tensor products over k, and · inside a word is "
                  "composition of arrows along a path (left to right) -- not a scalar.")


def hh_space_typing(theory, route):
    """The one-paragraph typing statement for an HH cohomology / homology section."""
    coh = theory in ("hh_cohomology", "HH^")
    if route == "cs":
        if coh:
            return ("What the engine computes: Hochschild cohomology as the cohomology "
                    "of Hom_{A^e}(P_•, A), where P_• → A is the Chouhy–Solotar projective "
                    "bimodule resolution with P_n = ⊕_σ A e_{o(σ)} ⊗ e_{t(σ)} A (σ over "
                    "the degree-n ambiguity chains). Degree n collapses to the corner "
                    "space C^n = ⊕_σ e_{o(σ)} A e_{t(σ)}. A basis element v ⊗ p pairs a "
                    "path v ∈ A with a chain word p = a_1·a_2·…." + _NOTATION_TAIL)
        return ("What the engine computes: Hochschild homology as the homology of "
                "A ⊗_{A^e} P_•, where P_• → A is the Chouhy–Solotar projective bimodule "
                "resolution with P_n = ⊕_σ A e_{o(σ)} ⊗ e_{t(σ)} A. Degree n collapses "
                "to the corner space C_n = ⊕_σ e_{t(σ)} A e_{o(σ)}. A basis element "
                "v ⊗ p pairs a path v ∈ A with a chain word p = a_1·a_2·…."
                + _NOTATION_TAIL)
    if coh:
        return ("What the engine computes: Hochschild cohomology as the cohomology of "
                "the normalized bar cochain complex. A degree-n cochain is a k-linear "
                "map C^n = Hom_k(Ā^⊗n, A), where Ā = A/(k·1) is the algebra modulo its "
                "unit (spanned by the arrows and longer paths); by the tensor–hom "
                "adjunction this is Hom_{A^e}(A ⊗ Ā^⊗n ⊗ A, A), a bimodule map out of "
                "the n-th term of the bar resolution of A. A basis functional written "
                "[w_1|…|w_n ↦ v] sends the single basis tensor w_1 ⊗ … ⊗ w_n (each "
                "w_i ∈ Ā) to v ∈ A and every other basis tensor to 0." + _NOTATION_TAIL)
    return ("What the engine computes: Hochschild homology as the homology of the "
            "normalized bar chain complex C_n = A ⊗_k Ā^⊗n (equivalently "
            "A ⊗_{A^e} (A ⊗ Ā^⊗n ⊗ A)). A basis chain written v ⊗ w_1 ⊗ … ⊗ w_n has "
            "v ∈ A and each w_i ∈ Ā." + _NOTATION_TAIL)


def module_reps_label_note(kind):
    """The Ext / Tor class-label explanation (Marco 2026-07-31): what P_v, P_v#k and
    n_{v,j} denote in a term-sum. ``kind`` is "ext" / "tor"."""
    common = ("Notation. In a class term-sum, P_v is the generator of the projective "
              "summand A e_v of the resolution term P_n, and P_v#k the k-th copy of "
              "P_v when the vertex v repeats among the summands of P_n; ")
    if kind == "ext":
        return (common + "n_{v,j} is the j-th basis vector of the vertex-v part N e_v "
                "of N. A basis functional [P_v#k → n_{w,j}] is the A-linear map sending "
                "that generator to n_{w,j} and every other generator to 0.")
    return (common + "n_{v,j} is the j-th basis vector of the vertex-v part e_v N of N. "
            "A basis element P_v#k ⊗ n_{w,j} is the tensor of that generator with "
            "n_{w,j} in P_n ⊗_A N.")


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
# Element-wise read-offs (Plan 35 wave 3d): a captured class's labelled term-sum
# ``[[coeff, word, value], ...]`` (``word`` the ordered list of arrow labels of the
# (co)chain slot, ``value`` a basis label of A) is RELABELLED into the classical
# dictionary object -- the central element (HH^0), the derivation (HH^1), the
# deformation 2-cocycle (HH^2), the commutator residue (HH_0). Presentation only: these
# invent nothing; they read the shipped reps. The words are the SAME labels the
# ``chain_basis`` enumeration lists, so a reader can trace every read-off back to a
# coordinate vector.
# --------------------------------------------------------------------------- #
def _magnitude(coeff, value):
    """One term ``coeff * value`` with a unit coefficient suppressed."""
    return value if str(coeff) == "1" else "%s %s" % (coeff, value)


def _arrow_word(word):
    """A (co)chain slot's word (a list/tuple of arrow labels, or a bare string) rendered
    as one composition label ``a·b·c`` (empty word -> '' at degree 0)."""
    if isinstance(word, (list, tuple)):
        return "·".join(str(w) for w in word)
    return str(word)


def _element_from_terms(terms):
    """A degree-0 (co)chain's element of A as a linear combination of its basis values
    (the words are empty at degree 0)."""
    parts = [_magnitude(coeff, value) for coeff, _word, value in terms]
    return " + ".join(parts) if parts else "0"


def derivation_values(terms):
    """From a degree-1 HH^1 cochain's labelled term-sum ``[[coeff, word, value], ...]``
    (each ``word`` a single arrow generator, as a bare label or a 1-element list), the
    derivation's action list ``D(arrow) = coeff * value``, grouped by arrow. Presentation
    only: it relabels the already-shipped terms as D(.) values, inventing nothing."""
    by_arrow = {}
    for coeff, word, value in terms:
        by_arrow.setdefault(_arrow_word(word), []).append(_magnitude(coeff, value))
    return ["D(%s) = %s" % (a, " + ".join(by_arrow[a])) for a in sorted(by_arrow)]


def central_element(terms):
    """From a degree-0 HH^0 cochain's term-sum, the central element ``z`` of ``Z(A)``."""
    return _element_from_terms(terms)


def commutator_residue(terms):
    """From a degree-0 HH_0 chain's term-sum, the residue in ``A / [A, A]``."""
    return _element_from_terms(terms)


def deformation_cochain(terms):
    """From a degree-2 HH^2 cochain's term-sum ``[[coeff, [a, b], value], ...]``, the
    infinitesimal-deformation 2-cocycle values ``μ(a, b) = coeff * value``, grouped by
    the argument pair (in first appearance order)."""
    by_pair, order = {}, []
    for coeff, word, value in terms:
        pair = _arrow_word(word)
        if pair not in by_pair:
            by_pair[pair] = []
            order.append(pair)
        by_pair[pair].append(_magnitude(coeff, value))
    return ["μ(%s) = %s" % (pair, " + ".join(by_pair[pair])) for pair in order]


# per (theory, degree): the interpretation heading + the read-off builder for one class.
_ELEMENT_READOFF = {
    ("hh_cohomology", 0): (
        "HH⁰ = Z(A): each class is a central element z (za = az for all a in A).",
        lambda t: [central_element(t)]),
    ("hh_cohomology", 1): (
        "HH¹ = Der(A) / Inn(A): each class is an outer derivation, read off as "
        "D(arrow) = value.", derivation_values),
    ("hh_cohomology", 2): (
        "HH² = infinitesimal deformations: each class is the 2-cocycle μ(a, b) of a "
        "first-order deformation a * b = ab + t·μ(a, b).", deformation_cochain),
    ("hh_homology", 0): (
        "HH₀ = A / [A, A]: each class is the residue of an element modulo the "
        "commutators ab − ba.", lambda t: [commutator_residue(t)]),
}
_ELEMENT_READOFF[("HH^", 0)] = _ELEMENT_READOFF[("hh_cohomology", 0)]
_ELEMENT_READOFF[("HH^", 1)] = _ELEMENT_READOFF[("hh_cohomology", 1)]
_ELEMENT_READOFF[("HH^", 2)] = _ELEMENT_READOFF[("hh_cohomology", 2)]
_ELEMENT_READOFF[("HH_", 0)] = _ELEMENT_READOFF[("hh_homology", 0)]


def element_heading(theory, n):
    """The interpretation heading for ``(theory, degree n)``'s element-wise read-off, or
    ``None`` when that space has no element-wise dictionary entry."""
    entry = _ELEMENT_READOFF.get((theory, int(n)))
    return entry[0] if entry is not None else None


def element_readoff(theory, n, terms):
    """The element-wise dictionary read-off of ONE class's labelled term-sum -- a list of
    display strings -- or ``None`` when ``(theory, degree n)`` has no element-wise entry
    (the framing sentence applies instead)."""
    entry = _ELEMENT_READOFF.get((theory, int(n)))
    return entry[1](terms) if entry is not None else None
