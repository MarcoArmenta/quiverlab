"""Length-lexicographic admissible path order (spec §5, component 3).

Fix the total order on arrows = their INSERTION ORDER in Quiver.arrows: the
arrow at position 0 is smallest. For arrow-words (paths read LEFT TO RIGHT,
Assem-Simson-Skowronski) u, v:

    u < v  iff  len(u) < len(v),
             or len(u) == len(v) and, at the first index i where they differ,
                arrow_index[u[i]] < arrow_index[v[i]].

Empty words (vertices) have length 0 and are the minimum.

Admissibility -- this is a WELL-FOUNDED, TWO-SIDED MULTIPLICATIVE order:

 * Total: any two words compare by (length, then finite lex on ranks).
 * Well-founded: length is a nonnegative integer bounded below by 0; along any
   strictly descending chain the length is non-increasing, so it stabilizes,
   after which the words live in the finite set of words of that fixed length,
   on which lex is a well-order. Hence no infinite descending chain -- so
   reduction by any rule (which replaces a word by strictly smaller words)
   terminates.
 * Two-sided multiplicative: if u < v and w*u*x, w*v*x are both composable
   paths, then w*u*x < w*v*x. Indeed len(w*u*x) - len(w*v*x) = len(u) - len(v)
   <= 0; if the lengths differ the shorter wins directly; if len(u) == len(v)
   the common prefix w matches, the first differing index of u,v reappears
   (shifted by len(w)) as the first differing index of w*u*x, w*v*x with the
   SAME rank comparison, and the common suffix x is irrelevant.

This two-sided compatibility is exactly what makes leading words multiplicative,
so Bergman's Diamond Lemma applies: a reduction system all of whose ambiguities
resolve is confluent, and the irreducible words form a k-basis of kQ/I.

CONVENTION: because paths compose LEFT TO RIGHT, an overlap of leading words
w1, w2 is a word b that is a SUFFIX of w1 and a PREFIX of w2 (w1 = a*b,
w2 = b*c, ambiguity word a*b*c = w1*c = a*w2). Right-to-left treatments (most of
the Groebner/diamond-lemma literature) swap the roles of prefix and suffix; we
write every overlap, inclusion, S-polynomial and reduction in the left-to-right
convention.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PathOrder:
    arrow_index: dict

    def key(self, word):
        """Total-order sort key: (length, ranks). Python tuple order == the order above."""
        return (len(word), tuple(self.arrow_index[a] for a in word))

    def compare(self, u, v):
        ku, kv = self.key(u), self.key(v)
        if ku < kv:
            return -1
        if ku > kv:
            return 1
        return 0

    def leading(self, comb):
        """The largest word appearing in the linear combination comb (word -> coeff),
        or None if comb is empty (the zero element). Callers keep comb free of
        zero coefficients, so every key is a genuine support word."""
        if not comb:
            return None
        return max(comb, key=self.key)

    def sort_words(self, words):
        """Ascending sort by the order (used for deterministic iteration)."""
        return sorted(words, key=self.key)


def path_order(quiver):
    """The admissible order whose arrow ranks are the insertion order of quiver.arrows."""
    return PathOrder(arrow_index={name: i for i, name in enumerate(quiver.arrows)})
