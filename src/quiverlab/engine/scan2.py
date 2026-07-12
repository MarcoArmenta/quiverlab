# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""
Search II for a counterexample to Han's conjecture -- a MECHANISM LEDGER.

New relative to Search I (which swept local algebras):
  * new structural classes:  tensor products A (x) B,  trivial extensions T(B)=B|<DB
    (symmetric),  triangular / one-point extensions A[M] (NOT self-injective).
  * a COMPLEXITY diagnostic: the polynomial growth rate of dim HH_n.

Reframing.  Restate Han's conjecture as
      gl.dim A = infinity   ==>   cx(A) >= 1,
where cx(A) := 1 + (degree of polynomial growth of dim_k HH_n(A)) is the complexity of the
Hochschild homology.  A COUNTEREXAMPLE is an algebra of infinite global dimension with cx = 0,
i.e. dim_k HH_n(A) = 0 for all n >> 0 (eventual vanishing).  So measuring cx per family is exactly
measuring the quantity whose positivity is the conjecture, and recording WHICH mechanism forces
cx >= 1 is the ledger we want for a proof.

Builders extend hh_engine.py.  Output: scan2_results.json.
"""

import json
import sys
import numpy as np
from quiverlab.engine.hh_engine import (Algebra, hochschild_homology_dims, check_associative,
                       truncated_polynomial, two_gen_local)


# ----------------------------------------------------------------------
# raw (m, T, unit) extraction and structural combinators
# ----------------------------------------------------------------------
def raw(alg):
    """extract (m, T, unit) from an Algebra (already associative, unit = e_t)."""
    unit = np.zeros(alg.m, dtype=np.int64)
    unit[alg.t] = 1
    return alg.m, alg.T.copy(), unit


def tensor_product(algA, algB, name=None):
    """A (x)_k B.  basis = product basis; (e_i (x) f_j)(e_i' (x) f_j') = (e_i e_i') (x) (f_j f_j')."""
    mA, TA, uA = raw(algA)
    mB, TB, uB = raw(algB)
    m = mA * mB
    T = np.zeros((m, m, m), dtype=np.int64)

    def idx(i, j):
        return i * mB + j
    for i in range(mA):
        for ip in range(mA):
            va = TA[i, ip, :]
            nza = np.nonzero(va)[0]
            if nza.size == 0:
                continue
            for j in range(mB):
                for jp in range(mB):
                    vb = TB[j, jp, :]
                    nzb = np.nonzero(vb)[0]
                    if nzb.size == 0:
                        continue
                    for k in nza:
                        for l in nzb:
                            T[idx(i, j), idx(ip, jp), idx(k, l)] += va[k] * vb[l]
    unit = np.zeros(m, dtype=np.int64)
    for i in np.nonzero(uA)[0]:
        for j in np.nonzero(uB)[0]:
            unit[idx(i, j)] = uA[i] * uB[j]
    return Algebra(m, T, unit, name=name or f"({algA.name}) (x) ({algB.name})")


def trivial_extension(algB, name=None):
    """
    T(B) = B |< DB,  the trivial (square-zero) extension by the dual bimodule DB=Hom_k(B,k).
    Symmetric algebra; always infinite global dimension (when B != 0 is not semisimple-trivial).
    basis: first mB are B; next mB are the dual basis beta_i.
      (b,f)(b',f') = (bb', b.f' + f.b'),   (b.f)(x)=f(xb),  (f.b)(x)=f(bx),  DB.DB=0.
    """
    mB, TB, uB = raw(algB)
    m = 2 * mB
    T = np.zeros((m, m, m), dtype=np.int64)
    # B * B block
    for i in range(mB):
        for j in range(mB):
            T[i, j, :mB] = TB[i, j, :]
    # left action  e_j . beta_i = sum_l TB[l,j,i] beta_l        (index e_j first, beta_i second)
    for j in range(mB):          # b = e_j
        for i in range(mB):      # f = beta_i
            for l in range(mB):
                c = TB[l, j, i]
                if c:
                    T[j, mB + i, mB + l] += c
    # right action  beta_i . e_j = sum_l TB[j,l,i] beta_l
    for i in range(mB):          # f = beta_i
        for j in range(mB):      # b = e_j
            for l in range(mB):
                c = TB[j, l, i]
                if c:
                    T[mB + i, j, mB + l] += c
    # DB * DB = 0  (already)
    unit = np.zeros(m, dtype=np.int64)
    unit[:mB] = uB
    return Algebra(m, T, unit, name=name or f"T({algB.name})")


def triangular_extension(algA, mod_actions, mod_dim, name):
    """
    One-point extension  A[M] = [[A, M],[0, k]]  (upper triangular).  NOT self-injective.
    M is a left A-module of dimension mod_dim given by action matrices mod_actions[i]
    (mod_dim x mod_dim), the left action of basis element e_i of A on M (columns = e_i . m_j).
    The action of the unit of A must be the identity.
    basis: A (mA), then M (mod_dim), then the extra idempotent eps (1).
    products:  e_i.e_k (A);  e_i.m_j (A acts on M);  m_j.eps = m_j;  eps.eps=eps; else 0.
    unit = 1_A + eps.
    """
    mA, TA, uA = raw(algA)
    mM = mod_dim
    m = mA + mM + 1
    eps = mA + mM
    T = np.zeros((m, m, m), dtype=np.int64)
    # A*A
    for i in range(mA):
        for k in range(mA):
            T[i, k, :mA] = TA[i, k, :]
    # e_i . m_j  -> sum_l mod_actions[i][l,j] * m_l   (M block)
    for i in range(mA):
        Ai = mod_actions[i]
        for j in range(mM):
            for l in range(mM):
                c = int(Ai[l, j])
                if c:
                    T[i, mA + j, mA + l] += c
    # m_j . eps = m_j
    for j in range(mM):
        T[mA + j, eps, mA + j] = 1
    # eps.eps = eps
    T[eps, eps, eps] = 1
    unit = np.zeros(m, dtype=np.int64)
    unit[:mA] = uA
    unit[eps] = 1
    return Algebra(m, T, unit, name=name)


def local_3gen_radsq(name="k<x,y,z>/(rad^2)"):
    """3-generated radical-square-zero local algebra: basis {1,x,y,z}, all products of x,y,z = 0.
       dim 4, local, infinite gl.dim, MONOMIAL (a control in a new arity)."""
    m = 4
    T = np.zeros((m, m, m), dtype=np.int64)
    for j in range(m):
        T[0, j, j] = 1
        T[j, 0, j] = 1
    # x,y,z (indices 1,2,3) multiply to 0
    unit = np.zeros(m, dtype=np.int64)
    unit[0] = 1
    return Algebra(m, T, unit, name=name)


# ----------------------------------------------------------------------
# COMPLEXITY diagnostic: growth rate of dim HH_n
# ----------------------------------------------------------------------
def complexity_diagnostic(dims):
    """
    Estimate cx = 1 + (poly growth degree of dim HH_n) over the computed range.
    Returns dict with 'type','complexity', and notes.  cx=0 <=> eventual vanishing.
    Honest: this reflects only the COMPUTED degrees, not a proof of eventual behaviour.
    """
    seq = dims[1:]  # homological degrees n >= 1 (HH_0 is always nonzero / special)
    if len(seq) == 0:
        return {"type": "no_data", "complexity": None}
    if all(d == 0 for d in seq):
        return {"type": "eventually_zero", "complexity": 0}
    if len(seq) < 3:
        return {"type": "insufficient_depth", "complexity": None}
    # eventual vanishing: trailing zeros (the Han counterexample signature). A
    # sequence that decays to zero through a non-constant tail, e.g. [2,1,0,0],
    # is cx 0 and must be flagged -- not mislabelled "growing". Mirrors
    # scan3.complexity_of so the two estimators agree.
    if seq[-1] == 0 and seq[-2] == 0:
        return {"type": "eventually_zero", "complexity": 0}
    last = seq[-4:] if len(seq) >= 4 else seq[-3:]
    if len(set(last)) == 1:
        return {"type": "eventually_constant", "value": int(last[0]), "complexity": 1}
    s = np.array(seq, dtype=np.int64)
    cur = s.copy()
    for k in range(1, 5):
        cur = np.diff(cur)
        if cur.size >= 2 and np.all(cur == 0):
            return {"type": f"polynomial_deg_{k-1}", "complexity": k}
        if cur.size >= 2 and np.all(cur == cur[0]) and cur[0] != 0:
            return {"type": f"polynomial_deg_{k}", "complexity": k + 1}
    # eventual periodicity
    for p in range(2, 5):
        if len(seq) >= 2 * p and all(seq[i] == seq[i - p] for i in range(len(seq) - p, len(seq))):
            if len(set(seq[-p:])) > 1:
                return {"type": f"eventually_periodic_p{p}", "complexity": 1}
    return {"type": "growing_irregular", "complexity": ">=2 (apparent)"}


# ----------------------------------------------------------------------
# analysis with ledger
# ----------------------------------------------------------------------
def analyze2(alg, N, primes, gldim_infinite, self_injective, monomial,
             witness, family, results):
    ok, bad = check_associative(alg)
    if not ok:
        results.append({"name": alg.name, "family": family, "dim": alg.m,
                        "associative": False, "bad_triple": list(bad)})
        print(f"[{family}] {alg.name}: NON-ASSOCIATIVE (bad {bad})")
        return None
    dims = hochschild_homology_dims(alg, N, primes=primes)
    cx = complexity_diagnostic(dims[32003])
    flags = []
    if gldim_infinite and cx.get("complexity") == 0:
        flags = [p for p in primes if complexity_diagnostic(dims[p]).get("complexity") == 0]
    rec = {
        "name": alg.name, "family": family, "dim": alg.m, "associative": True,
        "gldim_infinite": gldim_infinite, "self_injective": self_injective,
        "monomial": monomial, "witness": witness, "max_degree": N,
        "HH": {str(p): dims[p] for p in primes},
        "complexity": cx,
        "counterexample_candidate_primes": flags,
    }
    results.append(rec)
    tag = "   <-- CX=0 CANDIDATE!" if flags else ""
    print(f"[{family}] {alg.name}: HH={dims[32003]}  cx={cx.get('complexity')} ({cx.get('type')}){tag}")
    return rec


# ----------------------------------------------------------------------
# VALIDATION of the new machinery
# ----------------------------------------------------------------------
def kx2():
    return truncated_polynomial(2)


def kA(n):
    """path algebra of linear A_n: 1->2->...->n, dimension n(n+1)/2."""
    # build via paths
    verts = list(range(1, n + 1))
    # basis: e_i (vertices) then paths p_{i->j} for i<j
    paths = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]
    basis = [("v", i) for i in verts] + [("p", i, j) for (i, j) in paths]
    idx = {b: k for k, b in enumerate(basis)}
    m = len(basis)
    T = np.zeros((m, m, m), dtype=np.int64)
    # idempotents
    for i in verts:
        T[idx[("v", i)], idx[("v", i)], idx[("v", i)]] = 1
    # vertex * path (convention: e_source . p = p,  p . e_target = p)
    for (i, j) in paths:
        pij = idx[("p", i, j)]
        T[idx[("v", i)], pij, pij] = 1   # e_source . p = p
        T[pij, idx[("v", j)], pij] = 1   # p . e_target = p
    # path * path: p_{i->j} . p_{j->l} = p_{i->l}   (compose when middle matches)
    for (i, j) in paths:
        for (j2, l) in paths:
            if j == j2:
                T[idx[("p", i, j)], idx[("p", j2, l)], idx[("p", i, l)]] += 1
    unit = np.zeros(m, dtype=np.int64)
    for i in verts:
        unit[idx[("v", i)]] = 1
    return Algebra(m, T, unit, name=f"kA_{n}")


def run_validation(results=None):
    print("=== VALIDATION OF NEW BUILDERS ===\n")
    A = tensor_product(kx2(), kx2(), "kx2 (x) kx2")
    ok, _ = check_associative(A)
    d_tensor = hochschild_homology_dims(A, 6)[32003]
    B = two_gen_local([0,0,0,0],[0,0,0,0],[0,0,0,1],"k[x,y]/(x^2,y^2)")
    d_local = hochschild_homology_dims(B, 6)[32003]
    print(f"tensor kx2(x)kx2 assoc={ok}  HH={d_tensor}")
    print(f"local  k[x,y]/(x^2,y^2)   HH={d_local}")
    print(f"  KUNNETH/tensor cross-check MATCH: {d_tensor == d_local}")
    A2 = kx2(); B2 = kA(2)
    hhA = hochschild_homology_dims(A2, 4)[32003]
    hhB = hochschild_homology_dims(B2, 4)[32003]
    AB = tensor_product(A2, B2, "kx2 (x) kA_2")
    hhAB = hochschild_homology_dims(AB, 4)[32003]
    conv = [sum(hhA[i]*hhB[n-i] for i in range(n+1) if n-i < len(hhB) and i < len(hhA))
            for n in range(5)]
    print(f"  KUNNETH MATCH HH(kx2(x)kA_2)={hhAB} vs conv={conv}: {hhAB[:5]==conv[:5]}")
    TB = trivial_extension(kA(2), "T(kA_2)")
    okT, _ = check_associative(TB)
    print(f"  T(kA_2) dim={TB.m} assoc={okT}")
    act_id = np.array([[1]]); act_x = np.array([[0]])
    Tri = triangular_extension(kx2(), [act_id, act_x], 1, "k[x]/(x^2)[k]")
    okTri, _ = check_associative(Tri)
    print(f"  k[x]/(x^2)[k] dim={Tri.m} assoc={okTri}\n")


def module_simple(a):
    """simple module k for k[x]/(x^a): x acts 0. actions indexed by basis {1,x,...}."""
    acts = [np.array([[0]]) for _ in range(a)]
    acts[0] = np.array([[1]])
    return acts, 1


def module_regular(a):
    """regular module A=k[x]/(x^a) over itself, dim a, basis {1,x,...,x^{a-1}}."""
    acts = []
    for p in range(a):                # action of x^p
        M = np.zeros((a, a), dtype=np.int64)
        for j in range(a):            # x^p . x^j = x^{p+j}
            if p + j < a:
                M[p + j, j] = 1
        acts.append(M)
    return acts, a


def module_semisimple(a, copies):
    """k^copies, x acts 0."""
    acts = [np.zeros((copies, copies), dtype=np.int64) for _ in range(a)]
    acts[0] = np.eye(copies, dtype=np.int64)
    return acts, copies


def run_tensor(results):
    print("\n===== FAMILY: TENSOR PRODUCTS (Kunneth mechanism) =====")
    cases = [
        (tensor_product(kx2(), kA(2), "kx2 (x) kA_2"), 3, "factor kx2 infinite-gldim"),
        (tensor_product(kx2(), kx2(), "kx2 (x) kx2"), 6, "two infinite-gldim factors"),
        (tensor_product(kx2(), truncated_polynomial(3), "kx2 (x) k[x]/(x^3)"), 3, "two infinite factors"),
        (tensor_product(truncated_polynomial(3), kA(2), "k[x]/(x^3) (x) kA_2"), 2, "factor k[x]/x^3 infinite"),
    ]
    for alg, N, note in cases:
        analyze2(alg, N, (32003, 2, 3), True, False, False,
                 "tensor: " + note, "tensor_product", results)


def run_trivext(results):
    print("\n===== FAMILY: TRIVIAL EXTENSIONS T(B) (symmetric; self-injective) =====")
    cases = [
        (trivial_extension(kA(2), "T(kA_2)"), 3),
        (trivial_extension(kx2(), "T(k[x]/(x^2))"), 6),
        (trivial_extension(truncated_polynomial(3), "T(k[x]/(x^3))"), 3),
        (trivial_extension(kA(3), "T(kA_3)"), 2),
    ]
    for alg, N in cases:
        analyze2(alg, N, (32003, 2, 3), True, True, False,
                 "self-injective (symmetric); support varieties apply", "trivial_extension", results)


def run_triang(results):
    print("\n===== FAMILY: TRIANGULAR / ONE-POINT EXTENSIONS A[M] (NOT self-injective) =====")
    a2 = kx2(); a3 = truncated_polynomial(3)
    acts_k2, d = module_simple(2)
    T1 = triangular_extension(a2, acts_k2, d, "k[x]/(x^2)[k]")
    acts_reg2, d = module_regular(2)
    T2 = triangular_extension(a2, acts_reg2, d, "k[x]/(x^2)[A_reg]")
    acts_ss, d = module_semisimple(2, 2)
    T3 = triangular_extension(a2, acts_ss, d, "k[x]/(x^2)[k+k]")
    acts_k3, d = module_simple(3)
    T4 = triangular_extension(a3, acts_k3, d, "k[x]/(x^3)[k]")   # FRONTIER: no 2-truncated cycle
    acts_reg3, d = module_regular(3)
    T5 = triangular_extension(a3, acts_reg3, d, "k[x]/(x^3)[A_reg]")  # FRONTIER, shallow
    spec = [
        (T1, 6, "loop x^2=0 IS a 2-truncated cycle (BHM applies)"),
        (T2, 5, "loop x^2=0 truncated cycle; non-self-injective"),
        (T3, 5, "loop x^2=0 truncated cycle; non-self-injective"),
        (T4, 5, "FRONTIER: loop x^2!=0, NO 2-truncated cycle, non-self-injective"),
        (T5, 2, "FRONTIER: no 2-truncated cycle, non-self-injective (shallow)"),
    ]
    for alg, N, note in spec:
        analyze2(alg, N, (32003, 2, 3), True, False, False,
                 "triangular: " + note, "triangular_extension", results)


def run_local2(results):
    print("\n===== FAMILY: HIGHER LOCAL (new arity / nilpotency) =====")
    analyze2(local_3gen_radsq(), 6, (32003, 2, 3), True, False, True,
             "rad^2=0, 3 generators (monomial)", "local_3gen", results)
    for a, N in [(6, 4), (7, 3)]:
        analyze2(truncated_polynomial(a), N, (32003, 2, 3), True, True, True,
                 "loop, higher nilpotency (self-injective)", "truncated_poly_high", results)


def summarize_and_save(results, outfile):
    assoc = [r for r in results if r.get("associative")]
    inf = [r for r in assoc if r.get("gldim_infinite")]
    cands = [r for r in results if r.get("counterexample_candidate_primes")]
    # complexity tally among infinite-gldim
    from collections import Counter
    cxcount = Counter(str(r["complexity"].get("complexity")) for r in inf)
    print("\n==================== LEDGER SUMMARY ====================")
    print(f"  algebras tested (associative):      {len(assoc)}")
    print(f"  infinite global dimension:          {len(inf)}")
    print(f"  complexity distribution (inf-gldim): {dict(cxcount)}")
    print(f"  cx=0 counterexample candidates:     {len(cands)}")
    if cands:
        for c in cands:
            print(f"     !! {c['name']}")
    else:
        print("  -> every infinite-gldim algebra has complexity >= 1 (HH nonzero in all top degrees)")
    summary = {
        "n_assoc": len(assoc), "n_infinite": len(inf),
        "complexity_distribution": dict(cxcount),
        "n_candidates": len(cands),
        "candidates": [c["name"] for c in cands],
    }
    with open(outfile, "w") as f:
        json.dump({"summary": summary, "results": results}, f, indent=2)
    print(f"  results written to {outfile}")


def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    if phase == "valid":
        run_validation()
        return
    results = []
    if phase in ("tensor", "all"):
        run_tensor(results)
    if phase in ("trivext", "all"):
        run_trivext(results)
    if phase in ("triang", "all"):
        run_triang(results)
    if phase in ("local2", "all"):
        run_local2(results)
    tag = {"tensor": "_T", "trivext": "_X", "triang": "_R", "local2": "_L", "all": ""}.get(phase, "_"+phase)
    if phase == "merge":
        allr = []
        for t in ["_T", "_X", "_R", "_L"]:
            try:
                allr += json.load(open(f"scan2_results{t}.json"))["results"]
            except FileNotFoundError:
                pass
        summarize_and_save(allr, "scan2_results.json")
    else:
        summarize_and_save(results, f"scan2_results{tag}.json")
