"""The CS ambiguity S-sequence (arXiv:1406.2300 §3). S_n = 𝒜_{n-1}: S_0=vertices,
S_1=arrows, S_2=tips, S_n(n≥3)=(n-1)-ambiguities. Since ambiguities depend ONLY on the
tips, S_n for A equals Bardzell's AP^n for the tip monomial algebra A_S = kQ/⟨S⟩. We
reuse the ported, bar-oracle-validated
engine.resolutions_bardzell.MonomialPresentation.{associated_paths, left_decomposition}
(Bardzell AP^n = CS 𝒜_{n-1} = quiverlab S_n)."""
from quiverlab.errors import DepthLimitError
from quiverlab.engine.resolutions_bardzell import MonomialPresentation
from quiverlab.resolutions_cs.terms import Chain

_HINT = ("CS ambiguity growth guard; raise max_cells only if the certified range is "
         "genuinely insufficient (SSequence(..., max_degree=D) reuses cached lower degrees)")


class SSequence:
    def __init__(self, rs, max_degree, max_cells=4_000_000):
        self.rs = rs
        self.max_degree = max_degree
        self.max_cells = max_cells
        self._cache = {}
        Q = rs.quiver
        self._vertices = list(Q.vertices)
        arrows = [(name, Q.source(name), Q.target(name)) for name in Q.arrows]  # arrow NAME as id
        self._pres = MonomialPresentation(self._vertices, arrows, list(rs.leading_words()))
        self._maxlen = max(2, max_degree) * max(self._pres.maxrel, 2) + 2       # n-ambiguity length cap

    def tip_presentation(self):
        return self._pres

    def _guard(self, n, count):
        if count > self.max_cells:
            last = max((k for k in self._cache if self._cache.get(k)), default=n - 1)
            raise DepthLimitError(
                f"CS S-sequence at degree {n}: {count} chains exceed max_cells="
                f"{self.max_cells}; certified through degree {last}", hint=_HINT)

    def S(self, n):
        if n in self._cache:
            return self._cache[n]
        if n > self.max_degree:                                   # GO-LOUD (open-item #1): an
            raise DepthLimitError(                                # out-of-range page is a contract
                f"S-sequence degree {n} exceeds the certified range "  # error, not a silent [].
                f"0..{self.max_degree}",
                hint="raise max_degree / top to reach this degree")
        if n < 0:                                                 # negative degrees are empty by
            return []                                             # convention (out of scope here)
        if n == 0:
            out = [Chain((), (), v, v, 0) for v in self._vertices]
        else:
            words = self._pres.associated_paths(n, self._maxlen)   # Bardzell AP^n = S_n
            out = []
            for w in words:
                blocks = tuple(tuple(b) for b in self._pres.left_decomposition(w, n))
                out.append(Chain(tuple(w), blocks,
                                 self._pres.path_src(w), self._pres.path_tgt(w), n))
        self._guard(n, len(out))
        out.sort(key=lambda ch: (self.rs.order.key(ch.word), ch.word))
        self._cache[n] = out
        return out
