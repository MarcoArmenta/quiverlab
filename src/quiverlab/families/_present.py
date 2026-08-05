"""Shared length-lex kernel-enumeration presentation extractor (Plan 44).

Given a quiver ``Q``, a map ``img`` from each arrow of ``Q`` to a coordinate vector in a
structure-constant target algebra ``T`` (words multiply via ``T.multiply``), extract the
ideal ``I = ker(pi: kQ -> T)`` by a length-lex mini-Groebner (the ``TrivialExtension``
idiom, ``families/trivial_extension.py``), build ``kQ/I`` over ``field``, and CERTIFY per
instance by ``dim(kQ/I) == dim_expected`` (widen the enumeration window once). Reused by
the one-point extension, repetitive slices, and the Gabriel-quiver recovery. The
``_relation_string`` / ``_solve_combo`` emitters are the same ones ``trivial_extension``
uses -- presentation-agnostic pure functions."""
from quiverlab.errors import NotFiniteDimensionalError, QuiverlabError
from quiverlab.families.trivial_extension import _relation_string, _solve_combo


def extract_relations(Q, img, T, dom, max_len):
    """Relation strings of ``ker(pi)``: process words by increasing length; keep a
    per-corner reduced set of images of normal words; extend only normal words (on the
    right); a word whose image reduces against prior normal forms of its corner emits the
    parallel relation (a zero image emits a monomial). Deterministic (length, then
    arrow-discovery order)."""
    def pi_word(word):
        acc = img[word[0]]
        for a in word[1:]:
            acc = T.multiply(acc, img[a])
        return acc

    corner_basis = {}                                    # (s_vertex, t_vertex) -> [(word, image)]
    rels = []
    level = [(a,) for a in Q.arrows]                     # length-1 words, arrow order
    length = 1
    while level and length <= max_len:
        nxt = []
        for word in level:
            image = pi_word(word)
            corner = (Q.word_source(word), Q.word_target(word))
            basis = corner_basis.setdefault(corner, [])
            c = _solve_combo([bimg for _bw, bimg in basis], image, dom)
            if c is None:                                # image independent -> new normal form
                basis.append((word, image))
                tv = Q.word_target(word)
                for a in Q.arrows:
                    if Q.source(a) == tv:
                        nxt.append(word + (a,))
            else:                                        # reducible -> emit the relation
                terms = [(dom.one(), word)]
                for kk, ck in enumerate(c):
                    if not dom.is_zero(ck):
                        terms.append((dom.neg(ck), basis[kk][0]))
                rels.append(_relation_string(terms))
        level = nxt
        length += 1
    return rels


def present_from_pi(Q, img, T, dom, dim_expected, base_bound, citations=()):
    """Build ``kQ/ker(pi)`` and certify ``dim == dim_expected`` (widen once). Loud
    ``QuiverlabError`` on a failed certificate; ``NotFiniteDimensionalError`` propagates
    if the presented backbone cannot certify finiteness."""
    B = None
    got = None
    for bound in (base_bound, base_bound + 1):
        rels = extract_relations(Q, img, T, dom, bound)
        try:
            cand = Q.algebra(relations=rels, field=dom)
        except NotFiniteDimensionalError:
            got = "infinite-dimensional"
            continue
        if cand.dim == dim_expected:
            B = cand
            break
        got = cand.dim
    if B is None:
        raise QuiverlabError(
            "presentation certificate failed: dim(kQ/I) must be %d, got %s even after "
            "widening the kernel-enumeration window to length %d"
            % (dim_expected, got, base_bound + 1),
            hint="the length-lex kernel enumeration did not capture every relation of "
                 "this presentation; please report it")
    B._family_citations = tuple(citations)
    return B
