"""Finiteness certificate for a completed reduction system (spec §3.3, §5 c.3).

After completion under a degree bound D, let L = max leading-word length. Every
overlap ambiguity word has length <= 2L-1 and every inclusion word length <= L,
so if 2L-1 <= D EVERY ambiguity was formed and reduced to 0 -- the system is
confluent and the irreducible words are a genuine k-basis of kQ/I. If 2L-1 > D
we refuse to certify (AdmissibilityError: raise the bound). Given a confluent
system we enumerate irreducible words with the Plan-01 automaton
(core.monomial.irreducible_paths): finitely many -> certified finite-dimensional;
a cycle -> NotFiniteDimensionalError naming the offending arrow cycle."""
from quiverlab.core.monomial import irreducible_paths
from quiverlab.errors import AdmissibilityError


def default_degree_bound(relations):
    """A CONSERVATIVE heuristic bound. If completion grows leading words past it,
    the certificate refuses to certify with a loud AdmissibilityError (never a wrong
    answer) telling the user to raise degree_bound; it never silently certifies an
    unconfluent system. Generous enough for every admissible v1 fixture."""
    m = max((r.max_length for r in relations), default=1)
    return max(8, 2 * m + 4)


def check_degree_bound(rules, degree_bound):
    L = max((len(r.lead) for r in rules), default=0)
    need = 2 * L - 1
    if need > degree_bound:
        raise AdmissibilityError(
            f"Groebner completion cannot be certified under degree_bound={degree_bound}: "
            f"leading words reach length {L}, so overlap ambiguities can have length up to "
            f"{need}, which exceeds the bound",
            hint=f"raise degree_bound to at least {need}",
        )


def certified_irreducibles(quiver, rules):
    """Nonempty irreducible words (avoiding every rule's leading word as a factor),
    sorted (len, word) exactly as core.monomial.irreducible_paths -- so the Groebner
    basis indexing agrees elementwise with the Plan-01 monomial route. Raises
    NotFiniteDimensionalError (naming an arrow cycle) when infinitely many exist."""
    forbidden = [r.lead for r in rules]
    return tuple(irreducible_paths(quiver, forbidden))
