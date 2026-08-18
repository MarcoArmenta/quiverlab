"""The classified critical / hypercritical unit-form lists (Plan 62 / R19) -- the
DECISION data for weak nonnegativity, plus the weak-positivity cross-oracle.

* CRITICAL_FORMS: the Euclidean (extended-Dynkin) unit forms ~A_n / ~D_n /
  ~E_{6,7,8}. Each is minimal-not-weakly-positive with an isotropic positive
  radical (q(radical) = 0). Weak positivity itself is decided by Ovsienko's box-6
  (tits.py); this list is only a CROSS-ORACLE there.

* HYPERCRITICAL_FORMS: minimal not-weakly-nonnegative unit forms (von Hohne 1996 /
  de la Pena Banach-26 (1990) / BJP 2019). Each carries a POSITIVE defect vector
  with q(defect) < 0. This is the PRIMARY decision data for weak nonnegativity: a
  form with a connected restriction matching a listed hypercritical form is not
  weakly nonnegative, witnessed by the tabulated defect. The minimal wild TREES
  T_{2,3,7} / T_{2,4,5} / T_{3,3,4} and the 3-Kronecker are encoded here; every
  gram + defect was verified in exact integer arithmetic (2026-08-07/08).

  THE T_{2,3,7} REFUTATION (encoded so no one reintroduces a <= 9-vertex cap): the
  star with arms 1,2,6 is a 10-VARIABLE hypercritical form whose defect
  (12,6,8,4,10,9,7,6,4,2) is SINCERE on all 10 vertices with entries up to 12 --
  no support/box cap reaches it. Hence there is no safe universal support/entry
  bound for weak nonnegativity (the design consequence: list-primary, honest None
  outside coverage).

* HYPERCRITICAL_COVERAGE: the vertex counts for which HYPERCRITICAL_FORMS is
  PROVABLY complete (so a certified weakly-nonnegative True may rest on
  list-completeness). The full von Hohne / de la Pena list (all non-tree
  hypercritical forms) is NOT transcription-verifiable here, so this is left EMPTY
  by design: certified weak-nonnegativity True is emitted through the exact
  positive-semidefinite certificate (PSD => no hypercritical restriction can
  exist, since hypercriticals are indefinite and restrictions of a PSD form are
  PSD -- a genuine, universal completeness certificate), never a guessed box, and
  never a fabricated tame verdict. Outside PSD the verdict is honest None.

Refs: von Hohne Proc. LMS (3) 73 (1996) 47-67; de la Pena Banach Center Publ. 26
(1990) 353-369; Barot-Jimenez-Gonzalez-de la Pena, Springer AA 25 (2019);
Bongartz Math. Ann. 269 (1984).
"""
from __future__ import annotations


# ---------------------------------------------------------------------------
# gram builders (integer symmetric, diagonal 2) from an undirected edge list
# ---------------------------------------------------------------------------
def _gram_from_edges(n, edges):
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        G[i][i] = 2
    for (i, j) in edges:
        G[i][j] -= 1
        G[j][i] -= 1
    return tuple(tuple(row) for row in G)


def _star(arms):
    """Star tree: center 0, arms of the given edge-lengths. Returns (n, gram)."""
    edges = []
    nxt = 1
    for length in arms:
        prev = 0
        for _ in range(length):
            edges.append((prev, nxt))
            prev = nxt
            nxt += 1
    n = nxt
    return n, _gram_from_edges(n, edges)


def _cycle(n):
    return n, _gram_from_edges(n, [(i, (i + 1) % n) for i in range(n)])


def _kron(m):
    # m parallel edges between two vertices -> gram [[2,-m],[-m,2]]
    return 2, ((2, -m), (-m, 2))


def _q(gram, d):
    n = len(gram)
    return sum(gram[i][j] * d[i] * d[j] for i in range(n) for j in range(n)) // 2


# ---------------------------------------------------------------------------
# CRITICAL (Euclidean) forms -- weak-positivity cross-oracle. q(radical) = 0.
# ---------------------------------------------------------------------------
def _crit(name, n, gram, radical):
    assert _q(gram, radical) == 0, (name, "radical not isotropic")
    return {"name": name, "gram": gram, "radical": tuple(radical)}


_A1_n, _A1_G = _kron(2)
_A2_n, _A2_G = _cycle(3)
_D4_n, _D4_G = _star([1, 1, 1, 1])
_E6_n, _E6_G = _star([2, 2, 2])
_E7_n, _E7_G = _star([1, 3, 3])
_E8_n, _E8_G = _star([1, 2, 5])

CRITICAL_FORMS = (
    _crit("~A1", _A1_n, _A1_G, (1, 1)),
    _crit("~A2", _A2_n, _A2_G, (1, 1, 1)),
    _crit("~D4", _D4_n, _D4_G, (2, 1, 1, 1, 1)),
    _crit("~E6", _E6_n, _E6_G, (3, 2, 1, 2, 1, 2, 1)),
    _crit("~E7", _E7_n, _E7_G, (4, 2, 3, 2, 1, 3, 2, 1)),
    _crit("~E8", _E8_n, _E8_G, (6, 3, 4, 2, 5, 4, 3, 2, 1)),
)


# ---------------------------------------------------------------------------
# HYPERCRITICAL forms -- primary weak-nonnegativity decision. q(defect) < 0.
# ---------------------------------------------------------------------------
def _hyp(name, n, gram, defect):
    assert _q(gram, defect) < 0, (name, "defect not negative")
    assert all(x > 0 for x in defect), (name, "defect not sincere/positive")
    return {"name": name, "gram": gram, "defect": tuple(defect)}


_K3_n, _K3_G = _kron(3)
_T334_n, _T334_G = _star([2, 2, 3])
_T245_n, _T245_G = _star([1, 3, 4])
_T237_n, _T237_G = _star([1, 2, 6])

HYPERCRITICAL_FORMS = (
    _hyp("K3", _K3_n, _K3_G, (1, 1)),                     # 3-Kronecker, q = -1
    _hyp("T334", _T334_n, _T334_G, (9, 6, 3, 6, 3, 6, 4, 2)),      # q = -3
    _hyp("T245", _T245_n, _T245_G, (11, 5, 8, 5, 3, 8, 6, 4, 2)),  # q = -2
    _hyp("T237", _T237_n, _T237_G,
         (12, 6, 8, 4, 10, 9, 7, 6, 4, 2)),               # 10-var refutation, q = -1
)

# Provably-complete vertex counts for the ENCODED indefinite hypercritical list.
# Empty by design: the full non-tree hypercritical classification is not
# transcription-verified here, so a certified weakly-nonnegative True rests on the
# exact positive-semidefinite certificate (see the module docstring), never on an
# assumed list completeness. Documented on the verification page (Plan 62 scope).
HYPERCRITICAL_COVERAGE = frozenset()


# ---------------------------------------------------------------------------
# permutation (relabelling) matching of a unit form against the lists
# ---------------------------------------------------------------------------
def _profile(gram, i):
    """The relabelling-invariant signature of vertex i: its diagonal plus the
    sorted multiset of its off-diagonal gram entries."""
    n = len(gram)
    return (gram[i][i], tuple(sorted(gram[i][j] for j in range(n) if j != i)))


def _iso(cand_gram, ref_gram):
    """Return a bijection perm (list: candidate index -> reference index) with
    cand_gram[i][j] == ref_gram[perm[i]][perm[j]] for all i, j, or None. Small n;
    profile-pruned backtracking (fast on the trees/Kroneckers we encode)."""
    n = len(cand_gram)
    if n != len(ref_gram):
        return None
    cand_prof = [_profile(cand_gram, i) for i in range(n)]
    ref_prof = [_profile(ref_gram, j) for j in range(n)]
    if sorted(cand_prof) != sorted(ref_prof):
        return None
    cands_for = {i: [j for j in range(n) if ref_prof[j] == cand_prof[i]]
                 for i in range(n)}
    order = sorted(range(n), key=lambda i: len(cands_for[i]))
    perm = [None] * n
    used = [False] * n

    def bt(k):
        if k == n:
            return True
        i = order[k]
        for j in cands_for[i]:
            if used[j]:
                continue
            ok = True
            for kk in range(k):
                ii = order[kk]
                if cand_gram[i][ii] != ref_gram[j][perm[ii]]:
                    ok = False
                    break
            if not ok:
                continue
            perm[i] = j
            used[j] = True
            if bt(k + 1):
                return True
            used[j] = False
            perm[i] = None
        return False

    return list(perm) if bt(0) else None


def matches_hypercritical(unit_form):
    """Does `unit_form` equal, up to a relabelling of variables, a listed
    hypercritical form? Returns {"name", "defect"} with the defect expressed in the
    caller's index order, or None. This is step (2) of is_weakly_nonnegative."""
    cand = unit_form.gram
    for entry in HYPERCRITICAL_FORMS:
        if len(entry["gram"]) != len(cand):
            continue
        perm = _iso(cand, entry["gram"])          # candidate idx -> reference idx
        if perm is not None:
            defect = tuple(entry["defect"][perm[i]] for i in range(len(cand)))
            return {"name": entry["name"], "defect": defect}
    return None
