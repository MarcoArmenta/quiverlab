"""ConvergenceReport + the standing self-certificate.

The load-bearing arbiter of the whole engine: on EVERY :class:`SpectralSequence`
construction, ``certify_convergence`` asserts

    sum_{p+q=n} dim E_inf^{p,q}  ==  dim H_n(Tot)   for every total degree n

-- a rank identity that must hold by construction (a bounded filtration of a
bounded complex converges strongly, Weibel 5.5.1 / metaplan section-2 SS brief). A
mismatch can only be a page / filtration bookkeeping bug, and it raises
:class:`~quiverlab.errors.QuiverlabError` loudly -- never a silent wrong page. No
floats (``src/`` AST gate)."""
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError


@dataclass
class ConvergenceReport:
    """The convergence data attached to every :class:`SpectralSequence`."""

    e_infinity_page: int          # E_r stabilizes by here (bounded filtration bound)
    degenerates_at: "int | None"  # least r with E_r == E_inf (all higher d_r zero)
    abutment: dict                # {n: dim H_n(Tot)} -- the target of convergence

    def collapse(self):
        """The sequence collapses early (degenerates at ``E_1`` or ``E_2``)."""
        return self.degenerates_at in (1, 2)

    def prose(self):
        """A human sentence for the worked-steps report."""
        degs = sorted(self.abutment)
        abut = ", ".join(f"H_{n}={self.abutment[n]}" for n in degs)
        if self.degenerates_at is None:
            deg = (f"it stabilizes only at the E_{self.e_infinity_page} page "
                   "(the bounded-filtration bound)")
        elif self.degenerates_at == 1:
            deg = "it degenerates at E_1 (E_1 = E_inf; the filtration is trivial)"
        else:
            deg = (f"it degenerates at E_{self.degenerates_at} "
                   "(all higher differentials vanish)")
        return (f"The spectral sequence converges to the homology of the total "
                f"complex ({abut}); {deg}.")


def _page_dims(ss, page):
    return {pq: page.dim(*pq) for pq in ss._candidates}


def certify_convergence(ss):
    """Run the standing self-certificate and return the :class:`ConvergenceReport`.

    ``e_infinity_page = max(width, height) + 1`` (E_r stabilizes by the p-extent, so
    this bound is generous). Asserts the ``E_inf`` totals equal the total-complex
    homology per degree; raises loudly on any mismatch. ``degenerates_at`` is the
    least ``r >= 1`` whose per-cell page dims already equal ``E_inf``'s (page dims
    are monotone non-increasing per cell, so this is exactly the degeneration
    page)."""
    e_inf = max(ss.width, ss.height) + 1
    abutment = dict(ss.total_homology_dims)
    einf_page = ss.page(e_inf)
    totals = {}
    for (p, q) in einf_page.spots:
        totals[p + q] = totals.get(p + q, 0) + einf_page.dim(p, q)
    ns = set(abutment) | set(totals)
    for n in ns:
        if totals.get(n, 0) != abutment.get(n, 0):
            raise QuiverlabError(
                "spectral sequence does not converge to its abutment: at total "
                f"degree {n} the E_inf totals sum to {totals.get(n, 0)} but "
                f"H_{n}(Tot) = {abutment.get(n, 0)} -- a page/filtration "
                "bookkeeping bug",
                hint="the standing self-certificate (E_inf totals == total "
                     "homology) is a rank identity that must hold by construction")
    einf_dims = _page_dims(ss, einf_page)
    degenerates_at = None
    for r in range(1, e_inf + 1):
        if _page_dims(ss, ss.page(r)) == einf_dims:
            degenerates_at = r
            break
    return ConvergenceReport(e_inf, degenerates_at, abutment)
