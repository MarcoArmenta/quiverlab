# Literature-oracle deep research (2026-07-25)

**Provenance:** Marco's request, 2026-07-25 — "more tests against examples done on
books and the literature … special attention to Cibils, Solotar, Marcos,
Lanzilotta, Suárez-Álvarez, Keller, Rickard, Happel, de la Peña, Lenzing, and
their coauthors." Four parallel deep-research agents (one per author cluster),
hard evidence rule: computed values recorded ONLY from fetched sources
(arXiv PDF/text), with theorem/example numbers; anything weaker is marked
abstract-only / secondary / UNVERIFIED. **These findings become new oracle
batteries via the backlog item added with this document.** Every battery
implemented from here must cite its source in the tests and on
`docs/verification.md` (standing Plan-22 rule).

**Status of the four cluster reports:**
- Solotar cluster — LANDED (below).
- Cibils cluster — LANDED (below).
- Happel / Keller / Rickard — LANDED (below).
- de la Peña / Lenzing / Marcos — LANDED (below).

---

## Cluster report: Solotar school (Solotar, Chouhy, Suárez-Álvarez, Lanzilotta, Redondo, Artenstein, …)

Scope note: quiverlab computes on **finite-dimensional** kQ/I only. Every source
below was fetched and text-extracted (arXiv PDF → pdftotext); theorem/example
numbers quoted. Two flagship Solotar-school topics are **excluded for
infinite-dimensionality** (documented at the end). The already-pinned QCI
homology at a=b=2 is the existing `test_bgms_quantum_ci_homology`; several
findings **extend** it.

### [bergh_erdmann_qci] Bergh & Erdmann (2008) — (Co)homology of quantum complete intersections
- **Citation:** P. A. Bergh, K. Erdmann, "Homology and cohomology of quantum
  complete intersections," Algebra & Number Theory 2 (2008), no. 5, 501–522.
  DOI 10.2140/ant.2008.2.501. arXiv:0709.3029. *(Already `BerghErdmann2008` in
  `references.bib`, but only its a=b=2 homology value is currently pinned.)*
- **Where:** Theorem 3.1 (Hochschild **homology**), Theorem 3.2 (**cohomology**).
- **Algebra:** A = k⟨x,y⟩/(xᵃ, yᵇ, xy − q·yx), a,b > 1 (codim-2 QCI), dim A =
  ab. Stated for **q not a root of unity**.
- **Computed data (verbatim):**
  - **Thm 3.2 (cohomology):** dim HH⁰=2, HH¹=2, HH²=1, HHⁿ=0 for n≥3 — i.e.
    **HH^• = [2,2,1,0,0,…] for every a,b>1**, independent of a,b and (per the
    paper's remark) of characteristic. Total HH^* dimension = 5.
  - **Thm 3.1 (homology):** dim HH₀ = a+b−1; for n≥1, dim HHₙ = a+b (char k
    divides both a,b), = a+b−1 (divides exactly one), = **a+b−2 (divides
    neither)**.
  - Internal check: a=b=2, char 0 → [3,2,2,…] = exactly quiverlab's current pin.
- **Quiverlab test sketch (char 0 / CC only — see caveat):**
  ```python
  A = Quiver([1], {"x":(1,1), "y":(1,1)}).algebra(
        relations=[f"x^{a}", f"y^{b}", "y*x - 2*x*y"], field=CC)   # q=2, not a root of unity
  assert cs_cohomology_dims(A, 6).dims == [2,2,1,0,0,0,0]          # NEW: any a,b>1
  assert cs_homology_dims(A, 6).dims   == [a+b-1] + [a+b-2]*6      # extends a=b=2 pin
  ```
  Points: (a,b) ∈ {(2,2),(2,3),(3,3),(2,4)} → homology [3,2,2,…],[4,3,3,…],
  [5,4,4,…],[5,4,4,…]; cohomology always [2,2,1,0,…]. Bound: degree 6–8.
  `engine="cs"` (non-monomial admissible).
  - **CRITICAL char caveat:** "q not a root of unity" needs an **infinite**
    field → testable **only over CC/char-0**. Over any GF(p), q is always a
    root of unity (different regime, infinite-dim HH^*) — **do not** run these
    vectors over GF(p). The char-p homology branches need an infinite field of
    char p, which quiverlab lacks — not testable, noted honestly.
- **Verification status:** fetched-source (arXiv:0709.3029, Thms 3.1/3.2 verbatim).

### [redondo_roman_triangular] Redondo & Román (2014) — HH of triangular string algebras
- **Citation:** M. J. Redondo, L. Román, "Hochschild cohomology of triangular
  string algebras and its ring structure," J. Pure Appl. Algebra 218 (2014),
  no. 5, 925–936. arXiv:1301.0516.
- **Where:** Theorem 4.3 (general dim HHⁿ), **Example 3** (explicit Aₙ family),
  Theorem 5.2 (ring structure).
- **Algebra (Example 3):** Aₙ = kQ/I, vertices 0..n, **two parallel arrows**
  αᵢ, βᵢ : (i−1)→i; I = ⟨αᵢαᵢ₊₁, βᵢβᵢ₊₁⟩ (length-2 monomial; mixed products
  survive). Monomial, gentle, triangular. Thm 4.3 is a purely combinatorial
  count ⇒ **char-independent**.
- **Computed data (verbatim, Example 3):**
  - **A₁** (2-Kronecker, no relations): **[1,3,0,0,…]**.
  - **A_{2m}**: **[1, 2m, 0, 0, …]** (A₂→[1,2,0,…], A₄→[1,4,0,…]).
  - **A_{2m+1}**: HH⁰=1, HH¹=2m+1, HH^{2m+1}=**2**, else 0 → A₃ →
    **[1,3,0,2,0,…]**, A₅ → **[1,5,0,0,0,2,0,…]** (vanishes at degree 2 then
    **revives** at 2m+1 — a non-monotone discriminator).
  - A₁'s HH¹=3 independently verified via the hereditary Euler characteristic
    (χ = 2 − 4 = −2 ⇒ HH¹ = 3).
  - **Theorem 5.2 (cup):** for any triangular string algebra,
    HHⁿ ∪ HHᵐ = 0 for all n,m>0 — **trivial ring in positive degrees**.
- **Quiverlab test sketch:**
  ```python
  def tri_string(n, field):
      arrows = {}
      for i in range(1, n+1): arrows[f"a{i}"]=(i-1,i); arrows[f"b{i}"]=(i-1,i)
      rels = [f"a{i}*a{i+1}" for i in range(1,n)] + [f"b{i}*b{i+1}" for i in range(1,n)]
      return Quiver(list(range(n+1)), arrows).algebra(relations=rels, field=field)
  # Bardzell/CS cohomology; test over GF(2), GF(3), CC:
  # A1:[1,3,0,0]  A2:[1,2,0,0]  A3:[1,3,0,2,0]  A4:[1,4,0,0,0]  A5:[1,5,0,0,0,2,0]
  # cup: any two HH^{>0} classes cup to 0 (native-cup oracle)
  ```
  Bound: degree ≤ 2m+2 (to catch the revival). **Recommend bar cross-check at
  pin time** for the A₃/A₅ revival values; A₁ already double-verified.
- **Verification status:** fetched-source (arXiv:1301.0516; Thm 4.3, Ex. 3,
  Thm 5.2 verbatim; A₁ independently Euler-verified).

### [redondo_roman_quadratic] Redondo & Román (2018) — Gerstenhaber structure, quadratic string algebras
- **Citation:** M. J. Redondo, L. Román, "Gerstenhaber algebra structure on the
  Hochschild cohomology of quadratic string algebras," Algebr. Represent.
  Theory 21 (2018). DOI 10.1007/s10468-017-9704-1. arXiv:1504.02495.
- **Where:** Thms 3.13, 3.14 (dim HHⁿ, **char-2-dependent**); Cor 3.16
  (recovers Ladkani's gentle formula); Example 3.1 (explicit quiver);
  Thms 4.8, 4.9 (cup/bracket NON-vanishing).
- **Algebra (Example 3.1):** vertices {1,2}, arrows α₁,β₁ : 1→2, α₂ : 2→1;
  I = ⟨α₁α₂, α₂α₁, β₁α₂⟩. **dim A = 6** (basis e₁,e₂,α₁,β₁,α₂,α₂β₁), infinite
  gl.dim.
- **Computed data:** explicit char-2 vs char≠2 branches for dim HHⁿ (char-2
  adds a G-term); Thm 4.8: certain HH∪HH **nonzero** when Gₙ≠∅ (contrast the
  2014 triviality); Thm 4.9: nonzero brackets in char 0. **Honest limitation:**
  the examples are presented via combinatorial sets, not assembled integer
  vectors — reconstructing exact HHⁿ by hand is convention-risky, so use as
  (a) a char-2-vs-char-3 differential test and (b) a cup ≠ 0 predicate for the
  Plan-20 native cup, with a bar cross-check to fix the integers at pin time.
- **Quiverlab test sketch:**
  ```python
  A2 = Quiver([1,2], {"a1":(1,2),"b1":(1,2),"a2":(2,1)}).algebra(
         relations=["a1*a2","a2*a1","b1*a2"], field=GF(2))   # vs field=GF(3)
  ```
- **Verification status:** fetched-source (arXiv:1504.02495; theorem statements
  + Example 3.1 quiver verbatim). Exact integer vectors: needs bar cross-check.

### [alsolotar_toupie] Artenstein, Lanzilotta & Solotar (2020) — HH of toupie algebras
- **Citation:** D. Artenstein, M. Lanzilotta, A. Solotar, "Gerstenhaber
  structure on Hochschild cohomology of toupie algebras," Algebr. Represent.
  Theory 23 (2020), 421–456. DOI 10.1007/s10468-019-09854-y. arXiv:1803.10310.
- **Where:** Def 1 (toupie quiver); Thm 3.3 + Prop 6.3 (HH¹, a-Kronecker
  block); Example 7.4.1 (full worked example). **Field ℂ throughout.**
- **Computed data:**
  - **a-Kronecker Q_a** (0⇒ω, a parallel arrows, no relations): dim HH¹ =
    (a−1)+a(a−1) = **a²−1**; hereditary ⇒ HHⁱ=0 (i≥2). **HH^• = [1, a²−1, 0, …]**
    (a=2→[1,3,0,…] matches Redondo–Román A₁; a=3→[1,8,0,…]; a=4→[1,15,0,…]).
    Char-independent (Euler-verified: χ = 2 − a² ⇒ HH¹ = a²−1).
  - **Example 7.4.1** (13 vertices / 15 arrows; one non-monomial relation +
    monomial branch): **HH^• = [1, 10, 3, 0, 4, 0, …]** + full bracket tables
    and HH¹ ≅ ⟨y₄⟩ ⊕ (sl₂(ℂ) ⋉ (⟨t₁,t₂⟩ ⋉ ⟨z₁₃,z₂₃,z₁₆,z₂₆⟩)).
    **CAVEAT: the quiver is given only as a FIGURE**; the agent's reconstruction
    did not reconcile with the stated vertex/arrow counts — **do not pin
    Example 7.4.1 without consulting the PDF figure**. The a-Kronecker
    sub-result is the safe oracle.
- **Quiverlab test sketch:**
  ```python
  A = Quiver([0,1], {f"a{i}":(0,1) for i in range(a)}).algebra(relations=[], field=CC)
  assert cs_cohomology_dims(A, 4).dims == [1, a*a-1, 0, 0, 0]     # a in {2,3,4}
  ```
- **Verification status:** fetched-source; a-Kronecker Euler-verified;
  Example 7.4.1 builder **unverified (figure missing)**.

### [chouhy_solotar_examples] Chouhy & Solotar (2015) — §7 worked examples (the CS-engine source)
- **Citation:** J. Algebra 432 (2015), 22–61. DOI 10.1016/j.jalgebra.2015.02.019.
  arXiv:1406.2300. *(Already `ChouhySolotar2015`.)*
- **Findings:** §7.1 gives the **explicit bimodule resolution differentials**
  for the dim-4 QCI (a byte-level resolution-term oracle — lower priority, the
  CS engine already structurally validated); §7.2 (Lemmas 7.2/7.3) the general
  xⁿ,yᵐ QCI resolution (supports the general-a,b tests above); the HH values
  themselves defer to Bergh–Erdmann. **EXCLUDED as infinite-dimensional:**
  §7.3 down-up algebras A(α,β,γ); Example 7.0.1 cubic k⟨x,y,z⟩/(xyz−x³−y³−z³).
- **Verification status:** fetched-source.

### [ladkani_gentle] Ladkani (2012) — HH of gentle algebras (ADJACENT, not Solotar cluster)
- **Citation:** S. Ladkani, "Hochschild cohomology of gentle algebras,"
  arXiv:1208.2230 (preprint). Included because Redondo–Román (2018) Cor 3.16
  independently recovers it; the cleanest closed form for the gentle class.
- **Computed data (verbatim):** with the Avella-Alaminos–Geiss invariant φ_Λ:
  dim HH⁰ = 1+φ_Λ(1,0); dim HH¹ = 1−χ(Q)+φ_Λ(1,1) (+φ_Λ(0,1) if char k = 2);
  n≥2: dim HHⁿ = φ_Λ(1,n) + aₙψ_Λ(n) + bₙψ_Λ(n−1), (aₙ,bₙ) =
  (1,0)/(0,1)/(1,1) for [char≠2, n even]/[char≠2, n odd]/[char 2];
  ψ_Λ(n)=Σ_{d|n} φ_Λ(0,d). Strong char-2 discriminator. Caveat: computing φ_Λ
  by hand is error-prone — pin the Redondo–Román integer vectors first.
- **Verification status:** fetched-source. Cluster note: adjacent, not Solotar.

### Exclusions (the infinite-dimensional trap, documented)
- **Reca & Solotar (2018), super Jordan plane** (J. Algebra 507, 120–185,
  arXiv:1707.05345): Nichols algebra, **infinite-dimensional**; excluded.
- **Down-up / quantum generalized Weyl** (CS §7.3; Suárez-Álvarez–Solotar–Vivas
  line): infinite-dimensional; excluded.
- **Han's-conjecture series** (Cibils–Lanzilotta–Marcos–Solotar,
  arXiv:1908.11130, 2101.02597, 2303.17369): structural closure theorems, no
  new small-algebra integer tables beyond the existing zoo scan framing.
- **Igusa–Todorov φ/ψ-dimension** (Lanzilotta et al.): φ-dim is not a quiverlab
  invariant; only bounds gl.dim. No clean mapping.
- **Herscovich–Solotar Yoneda/A∞**: mostly infinite-dim graded; a lead for the
  Plan-27+ Ext-algebra surface, no pinnable f.d. value found.

### BibTeX (fetched-source; first three NEW to references.bib)

```bibtex
@article{RedondoRoman2014,
  author  = {Redondo, Mar\'ia Julia and Rom\'an, Lucrecia},
  title   = {Hochschild cohomology of triangular string algebras and its ring structure},
  journal = {Journal of Pure and Applied Algebra},
  volume  = {218}, number = {5}, year = {2014}, pages = {925--936},
  doi     = {10.1016/j.jpaa.2013.10.011}, note = {arXiv:1301.0516},
}
@article{RedondoRoman2018,
  author  = {Redondo, Mar\'ia Julia and Rom\'an, Lucrecia},
  title   = {Gerstenhaber algebra structure on the {H}ochschild cohomology of quadratic string algebras},
  journal = {Algebras and Representation Theory},
  volume  = {21}, year = {2018}, pages = {61--86},
  doi     = {10.1007/s10468-017-9704-1}, note = {arXiv:1504.02495},
}
@article{ArtensteinLanzilottaSolotar2020,
  author  = {Artenstein, Dalia and Lanzilotta, Marcelo and Solotar, Andrea},
  title   = {Gerstenhaber structure on {H}ochschild cohomology of toupie algebras},
  journal = {Algebras and Representation Theory},
  volume  = {23}, year = {2020}, pages = {421--456},
  doi     = {10.1007/s10468-019-09854-y}, note = {arXiv:1803.10310},
}
@misc{Ladkani2012gentle,
  author = {Ladkani, Sefi},
  title  = {Hochschild cohomology of gentle algebras},
  year   = {2012}, note = {arXiv:1208.2230},
}
```
(`BerghErdmann2008` exists; add its arXiv id `0709.3029` to the note field.)

### Best 5 batteries (agent ranking, highest oracle value per effort)

1. **QCI cohomology [2,2,1,0,…] for all a,b** (BE Thm 3.2, CC only) — new,
   cheap, striking; pins `cs_cohomology_dims` on QCI which the suite never
   checks. Points (a,b) ∈ {(2,2),(2,3),(3,3)}.
2. **QCI homology general a,b [a+b−1, a+b−2, …]** (BE Thm 3.1, CC) — extends
   the single a=b=2 pin to a family.
3. **Triangular-string Aₙ family (all char)** — monomial ⇒ Bardzell + CS;
   A₃ = [1,3,0,2,0] revival discriminator; plus cup-triviality (Thm 5.2).
4. **a-Kronecker [1, a²−1, 0, …]** — first nonzero-HH¹ hereditary oracle
   (complements the existing tree-quiver [1,0,0,…]).
5. **Char-2 discriminator + cup ≠ 0** (RR 2018 + Ladkani) — GF(2)-vs-GF(3)
   branch behavior; native cup pinned against a NON-zero prediction (mirror of
   #3's triviality). Bar cross-check to fix integers.

---

## Cluster report: Cibils school (Cibils, Redondo, Saorín, Marcos; method-descendants Taillefer, Xu–Han–Jiang, Sánchez-Flores, Han)

**Methodological warning from the agent:** automated PDF/HTML summarizers repeatedly
garbled per-degree formulas during this research; the load-bearing
radical-square-zero numbers below were therefore reconstructed from the fetched
cochain complex and **hand-verified** against independently known results
(sl2 for Kronecker; the char-2 anomaly of k[x]/(x²)). Publisher PDFs
(ScienceDirect/Springer/HAL) are bot-walled; ar5iv HTML mirrors + secondary
restatements were the fetched sources. Verification status per entry.

### [cibils1998radsq] Cibils (1998) — HH algebra of radical square zero algebras
- **Citation:** C. Cibils, in *Algebras and Modules II* (Geiranger 1996), CMS
  Conf. Proc. 24, AMS, 1998, 93–101. Framework fetched via Wang arXiv:1511.08348
  (§2–3, "Cibils' Lemma") and Sánchez-Flores arXiv:0711.2810.
- **Framework (fetched verbatim):** for A = kQ/J², the reduced cochain complex is
  `Cⁿ = k(Qₙ//Q₀) ⊕ k(Qₙ//Q₁)` (parallel paths: pairs with same source+target);
  the differential collapses to the two outer terms, splitting as
  `Dₙ: k(Qₙ//Q₀) → k(Q_{n+1}//Q₁)`, so
  `HHⁿ = [k(Qₙ//Q₁)/im D_{n−1}] ⊕ ker Dₙ`.
- **Worked oracles (derived + hand-verified):**
  - `k[x]/(x²)` (one loop): char ≠ 2 → HH⁰=2, HHⁿ=1 (n≥1); **char 2 → HH⁰=2,
    HHⁿ=2 (n≥1)** (differential dies). First-class GF(2)-vs-CC regression pin.
  - Kronecker kK₂ (2 vertices, 2 parallel arrows): HH⁰=1, HH¹=3 (≅ sl₂),
    HHⁿ=0 (n≥2), all char. Euler check 2−4=−2 ✓.
  - Two loops k⟨x,y⟩/(x,y)² over CC (fetched from arXiv:0711.2810):
    HH¹≅gl₂ (dim 4); HH² dim 4; HH³ dim 8; HH⁷ dim 121 (sl₂-module
    decompositions given in the source).
- **Test sketch:** RadicalSquareZero over CC + GF(2) + GF(3), top ≈ 7.
- **Status:** framework fetched-source; loop/Kronecker **derived+hand-verified**;
  Cibils 1998 primary abstract-only (walled).

### [taillefer2001taft] Taillefer (2001) — Cyclic homology of Taft algebras + their Auslander algebras
- **Citation:** R. Taillefer, arXiv:math/0009214 (cf. K-Theory 24 (2001));
  computations via **Cibils' mixed complex**. Field ⊇ ℚ (**char 0 only**).
- **Where:** Thm 2.2 + example; Cor 2.8; Prop 3.6.
- **Algebra:** Taft Λₙ = kΔₙ/mⁿ — the self-injective cyclic Nakayama algebra
  (n vertices, Loewy length n), dim n².
- **Computed data (fetched, ar5iv):** HH₀=kⁿ; HH_p=k^{n−1} for all p≥1.
  HC_{2c}=kⁿ, HC_{2c+1}=k^{n−1} for all c≥0. Auslander algebra ΓΛₙ:
  HH₀=k^{n²}, HH_{≥1}=0; HC_{2p}=k^{n²}, HC_{2p+1}=0.
- **Test sketch:** cyclic self-injective NakayamaAlgebra / cyclic_nakayama(n,n)
  over CC, hochschild_homology + cyclic_homology to depth ≈ 6. Char-0 only.
- **Status:** fetched-source.

### [cibilsmarcosredondosolotar2003] CMRS (2003) — Cohomology of split algebras and trivial extensions
- **Citation:** Glasgow Math. J. 45 (2003) 21–40; arXiv:math/0102194;
  DOI:10.1017/S0017089502008948.
- **Computed data (fetched, ar5iv):** the LES `→ Hⁿ(Λ,M) → HHⁿ(Λ) → Hⁿ(Λ,A) →`;
  for M one-sided projective with M²=0: `Hⁿ(Λ,X) = ⊕_{p+q=n} Ext^q_{A-A}(M^{⊗_A p}, X)`;
  connecting map = cup with 1_M. **HH¹(T(A)) ≠ 0 for EVERY f.d. A** (Z(A) is
  always a summand).
- **Test sketch:** TrivialExtension across the zoo: assert HH¹(T(A)) ≠ 0 (any
  char). Structural, no per-quiver tables — exact numbers from the next entry.
- **Status:** fetched-source.

### [cibilsredondosaorin2004] Cibils–Redondo–Saorín (2004) — HH¹ of the trivial extension of a monomial algebra
- **Citation:** J. Algebra Appl. 3 (2004) 143–159; arXiv:math/0210284.
- **Computed data (fetched, ar5iv):**
  `HH¹(TA) = Z(A) ⊕ HH¹(A) ⊕ HH₁(A)* ⊕ Alt(DA)`;
  char 0: `dim HH¹(A) = s + Σ_{C∈E} w_C − e` (s = nontrivial strong circuits,
  E = efficient circuits, e = |E|, w_C as defined there); char p: e replaced by
  e_{p'} (p'-circuits) — explicit char dependence.
  **Example 2.20:** Q₀ = Z₅ cycle, arrows aᵢ: i→i+1, relations a₄a₃a₂ = 0 =
  a₃a₂a₁ (note the paper's composition order — translate to our left-to-right)
  ⇒ HH¹(A) = 0 in char 0.
- **Test sketch:** build Example 2.20 raw; pin HH¹(A)=0; pin the HH¹(TA)
  decomposition; char-0 vs GF(2)/GF(5) swap on a cyclic quiver.
- **Status:** fetched-source.

### [redondo2008incidence] Redondo (2008) + Cibils (1989) — HH of incidence algebras ≅ simplicial cohomology
- **Citation:** J. London Math. Soc. 77 (2008) 465–480, arXiv:math/0611542;
  original: C. Cibils, J. Pure Appl. Algebra 56 (1989) 221–232.
- **Computed data (fetched):** `HHⁿ(I(Σ)) ≅ SHⁿ(order complex of Σ; k)`;
  vanishing when the order complex is contractible ("no crowns").
- **Test sketch:** IncidenceAlgebra(covers): crown/diamond poset with nerve ≃ S¹
  ⇒ HH⁰=1, HH¹=1, HH^{≥2}=0; contractible nerve ⇒ HH^{≥1}=0. Degree ≤ 3.
  (S¹ pin is DERIVED from the theorem — label as such in the test.)
- **Status:** iso fetched-source; numeric pins derived-from-theorem.

### [xuhanjiang2007truncated] Xu–Han–Jiang (2007) — HH of truncated quiver algebras
- **Citation:** Sci. China Ser. A 50 (2007) 727–736; arXiv:math/0509202.
- **Computed data:** **Thm 3 (clean, fetched): dim_k HH•(A) < ∞ ⇔ gl.dim A < ∞
  ⇔ Q has no oriented cycle** — a boolean oracle across the zoo for kQ/Rᴺ.
  Basic-cycle formula (Thm 1) with a +1 bump exactly when char k | N and
  e | (Ni−N+1): **fetched but transcription-lossy — confirm against the primary
  or our own bar oracle before freezing exact counts.**
- **Status:** Thm 3 fetched-source; exact counts lossy (needs confirmation).

### [cibils1990rigidity] Cibils (1990) — Rigidity of truncated quiver algebras
- Adv. Math. 79 (1990) 18–42. HH²=0 ⇔ rigidity via parallel-path combinatorics.
  **Abstract/citation-only (walled)** — verdict-level use only until read.

### [cibils1990twonilpotent] Cibils (1990) — Cyclic and Hochschild homology of 2-nilpotent algebras
- K-Theory 4 (1990) 131–141. The (b,B) mixed complex specialized to r²=0 —
  the method under Taillefer's computations. **Citation-only**; exact formulas
  not recovered verbatim.

### [cibils1986nocycles] Cibils (1986) — HH_* vanishing for acyclic quivers
- LNM 1177 (1986) 55–59: Q acyclic ⇒ HH_n(A)=0 for n≥1 (degree-0 concentration).
  **Citation-verified** (consistent with fetched Han/XHJ statements).
- Test sketch: acyclic PathAlgebra/IncidenceAlgebra: HH_{≥1}=0, dim HH₀ = |Q₀|
  for A = kQ (paths mod commutators — for kQ/I check the statement scope before
  pinning dim HH₀ beyond hereditary).

### [han2010truncatedcycles] Han (2010) — HH homology dim + truncated oriented cycles
- arXiv:1004.0748 (fetched, ar5iv): a 2-truncated oriented cycle ⇒ hh.dim = ∞ =
  gl.dim; nonzero HH_m at degrees lm−1 (l = cycle length, odd m). Remark 6: a
  2-vertex example with infinite gl.dim and no truncated cycles. Feeds the Han
  bank / batch scans.
- **Status:** fetched-source.

### BibTeX (Cibils cluster)

```bibtex
@incollection{cibils1998radsq,
  author = {Cibils, Claude},
  title = {Hochschild cohomology algebra of radical square zero algebras},
  booktitle = {Algebras and Modules II (Geiranger, 1996)},
  series = {CMS Conf. Proc.}, volume = {24}, pages = {93--101},
  publisher = {Amer. Math. Soc., Providence, RI}, year = {1998}
}
@article{cibils1990rigidity,
  author = {Cibils, Claude}, title = {Rigidity of truncated quiver algebras},
  journal = {Adv. Math.}, volume = {79}, number = {1}, pages = {18--42},
  year = {1990}, doi = {10.1016/0001-8708(90)90057-T}
}
@article{cibils1989incidence,
  author = {Cibils, Claude},
  title = {Cohomology of incidence algebras and simplicial complexes},
  journal = {J. Pure Appl. Algebra}, volume = {56}, number = {3},
  pages = {221--232}, year = {1989}
}
@article{cibils1990twonilpotent,
  author = {Cibils, Claude},
  title = {Cyclic and Hochschild homology of 2-nilpotent algebras},
  journal = {K-Theory}, volume = {4}, number = {2}, pages = {131--141}, year = {1990}
}
@incollection{cibils1986nocycles,
  author = {Cibils, Claude},
  title = {Hochschild homology of an algebra whose quiver has no oriented cycles},
  booktitle = {Representation Theory I (Ottawa, 1984)},
  series = {Lecture Notes in Math.}, volume = {1177}, pages = {55--59},
  publisher = {Springer, Berlin}, year = {1986}
}
@article{cibilsmarcosredondosolotar2003,
  author = {Cibils, Claude and Marcos, Eduardo and Redondo, Mar\'{\i}a Julia and Solotar, Andrea},
  title = {Cohomology of split algebras and of trivial extensions},
  journal = {Glasgow Math. J.}, volume = {45}, number = {1}, pages = {21--40},
  year = {2003}, note = {arXiv:math/0102194}, doi = {10.1017/S0017089502008948}
}
@article{cibilsredondosaorin2004,
  author = {Cibils, Claude and Redondo, Mar\'{\i}a Julia and Saor\'{\i}n, Manuel},
  title = {The first cohomology group of the trivial extension of a monomial algebra},
  journal = {J. Algebra Appl.}, volume = {3}, number = {2}, pages = {143--159},
  year = {2004}, note = {arXiv:math/0210284}
}
@article{redondo2008incidence,
  author = {Redondo, Mar\'{\i}a Julia},
  title = {Hochschild cohomology via incidence algebras},
  journal = {J. London Math. Soc. (2)}, volume = {77}, number = {2},
  pages = {465--480}, year = {2008}, note = {arXiv:math/0611542}
}
@article{taillefer2001taft,
  author = {Taillefer, Rachel},
  title = {Cyclic homology of the Taft algebras and of their Auslander algebras},
  year = {2001}, note = {arXiv:math/0009214; cf. K-Theory 24 (2001)}
}
@article{xuhanjiang2007truncated,
  author = {Xu, Yunge and Han, Yang and Jiang, Wenfeng},
  title = {Hochschild cohomology of truncated quiver algebras},
  journal = {Sci. China Ser. A}, volume = {50}, number = {5}, pages = {727--736},
  year = {2007}, note = {arXiv:math/0509202}
}
@article{sanchezflores2008lie,
  author = {S\'anchez-Flores, Selene},
  title = {The Lie module structure on the Hochschild cohomology groups of monomial algebras with radical square zero},
  journal = {J. Algebra}, year = {2008}, note = {arXiv:0711.2810}
}
@article{han2010truncatedcycles,
  author = {Han, Yang},
  title = {Hochschild homology, global dimension, and truncated oriented cycles},
  year = {2010}, note = {arXiv:1004.0748}
}
```

### Best 5 batteries (agent ranking)

1. **Taft/Nakayama cyclic homology (char 0)** — HH₀=n, HH_{≥1}=n−1; HC even/odd
   = n / n−1; Auslander concentrated in degree 0. Family exists; tests
   cyclic_homology + HH_n at once. Highest value, lowest effort.
2. **Rad²=0 parallel-path HH** — hand-verified char-2 doubling on k[x]/(x²) and
   Kronecker HH¹=3; GF(2)-vs-CC parity regression.
3. **Incidence HHⁿ ≅ SHⁿ(nerve)** — poset topology as a graded oracle.
4. **Trivial-extension HH¹** — Example 2.20 (HH¹=0), the universal
   decomposition, never-vanishing HH¹(TA); char-0/char-p swap.
5. **Truncated finiteness boolean** (dim HH• < ∞ ⇔ acyclic) zoo-wide + the
   basic-cycle char|N bump (pending primary confirmation).

---

## Cluster report: de la Peña / Lenzing / Marcos (+ Green, Martínez-Villa, Barot; method-descendants Cassidy, Leader–Snashall, Etingof–Eu, Hubery, Schremmer)

### Convention mapping (load-bearing, verified)

quiverlab (`invariants/cartan.py`): Cartan `C[i][j] = dim e_i A e_j`; Coxeter
`Φ = −C^{−T} C`; `coxeter_polynomial` = charpoly(Φ) in t. Papers use
`Φ_A = −C^{−1}C^t`. **Verified numerically (3 random triangular Cartans):
`charpoly(−C^{−T}C) ≡ charpoly(−C^{−1}C^t)`** — quiverlab's Coxeter polynomial
equals the papers' exactly (variable rename only), for all triangular algebras
(covers every example below). The agent additionally re-derived, in quiverlab's
convention, the full Dynkin table, the canonical formula, Lehmer's [2,3,7], and
all six Nakayama polynomials — all match.

**Methodology warning:** WebFetch auto-summaries of these PDFs were unreliable
(one stated the OPPOSITE of a theorem); every verbatim value came from
pdftotext on a downloaded PDF. Labels: fetched-source / secondary /
recomputed / UNVERIFIED.

### [LdlP2008-spectral] Lenzing–de la Peña (2008) — Spectral analysis of f.d. algebras and singularities
- **Citation:** EMS Congr. Rep. (ICRA XII) 2008, 541–588. arXiv:0805.1018.
- **Dynkin path algebras (verbatim table; orientation-independent):**
  Aₙ: χ = v_{n+1} = ∏_{d|(n+1), d>1} Φ_d (Coxeter number n+1);
  Dₙ: χ = v₂·v_{2(n−1)}/v_{n−1} = Φ₂·(x^{n−1}+1) (Coxeter number 2(n−1));
  E₆: Φ₃Φ₁₂ (12); E₇: Φ₂Φ₁₈ (18); E₈: Φ₃₀ (30). [v_n := (x^n−1)/(x−1).]
  **Recomputed expansions:** E₆ = t⁶+t⁵−t³+t+1; E₇ = t⁷+t⁶−t⁴−t³+t+1;
  E₈ = t⁸+t⁷−t⁵−t⁴−t³+t+1; A₄ = Φ₅; D₄ = Φ₂²Φ₆; D₅ = Φ₂Φ₈; D₆ = Φ₂²Φ₁₀.
  **⚠ Dₙ caveat:** the pdftotext of the cyclotomic-condition column is
  mis-transcribed ("d|2(n−1), d≠1, d≠n−1" — wrong Φ₂ multiplicity); SHIP the
  v-factorization + small cases, not the condition (correct form should be
  d|2(n−1), d∤(n−1) — check the published table before citing it).
- **Star formula (eq. 1):** χ_{[p₁..p_t]} = (∏v_{p_i})·[(x+1) − x·Σ v_{p_i−1}/v_{p_i}];
  sign of χ(1) ⇔ Dynkin / extended-Dynkin / wild.
- **Cross-invariant identities (both sides quiverlab-computable):**
  - **Prop 1.2 (Happel):** degree-1 coefficient of χ_A = −trace Φ_A =
    Σ_{i≥0} (−1)^i dim HH^i(A). (Verified on A₂.)
  - **Prop 1.5:** χ_A(−1) is a perfect square (triangular A, k alg. closed).
  - **Prop 1.6 (A'Campo):** bipartite hereditary: χ_A(x²) = xⁿ κ_Δ(x+x^{−1}).
  - **B₃ parametric example:** Cartan [[1,a,b],[0,1,c],[0,0,1]] ⇒
    χ = x³+αx²+αx+1, α = abc−a²−b²−c²+3 (recomputed a=b=c=1 → (t+1)(t²+1)).
- **Status:** fetched-source + recomputed (Dₙ column transcription-flagged).

### [LdlP2008-affine] Extended-Dynkin + canonical formula (same paper; corroborated by arXiv:1310.1557 §2.3 and 1310.1910 §2.3–2.6)
- Ã_{p,q}: χ = (x−1)² v_p v_q (the ONLY orientation-dependent case);
  D̃ₙ: (x−1)² v₂² v_{n−2}; Ẽ₆: (x−1)² v₂ v₃²; Ẽ₇: (x−1)² v₂ v₃ v₄;
  Ẽ₈: (x−1)² v₂ v₃ v₅.
- **Canonical algebras (Prop 2.4):** χ_{Λ(p,λ)} = (x−1)² ∏ v_{p_i} — weights
  only, determines them up to order. Recomputed: C(2,2,2)=(x−1)²(x+1)³;
  C(3,3,3)=(x−1)²(x²+x+1)³; C(2,4,4)=(x−1)²(x+1)³(x²+1)²;
  C(2,3,6)=(x−1)²(x+1)²(x²−x+1)(x²+x+1)²; C(2,3,5)=(x−1)²(x+1)(x²+x+1)Φ₅.
- **Status:** fetched-source (3 independent PDFs) + recomputed. Rock-solid.

### [LdlP-cheby] Lenzing–de la Peña (2009) — Chebysheff recursion for Coxeter polynomials
- **Citation:** Linear Algebra Appl. 430 (2009) 947–956. arXiv:math/0611535.
- **Cor 3.2:** f̂_{(p₁..p_t)}(T) = (T+1)(T−1)² ∏v_{p_i}(T) − T·f_{[p₁..p_t]}(T)
  (extended canonical = one-point extension of canonical by an indec.
  projective). **Prop 3.4(d):** f̂_{(..,m+1)} = T·f̂_{(..,m)} − f̂_{(..,m−1)}.
  Prop 3.4(b): f_{(p₁..p_t)}(T²) = χ*_{K₂}(T)·∏ χ*_{[p_i−1]}(T).
- **Status:** fetched-source.

### [dlP-mahler] de la Peña (2013/14) — Mahler measure of Coxeter polynomials ★ spectral.py flagship
- **Citation:** arXiv:1310.1910 (Adv. Math. 2014).
- **§2.5:** the wild hereditary star [2,3,7] (10 vertices) has
  χ = T¹⁰+T⁹−T⁷−T⁶−T⁵−T⁴−T³+T+1 = **Lehmer's polynomial**; spectral radius =
  Mahler measure = μ₀ = 1.176280… (smallest known Salem number). Ordering
  c > ρ_{[2,4,5]} > ρ_{[2,3,m]} > ρ_{[2,3,7]} = μ₀ (m≥8), c = real root of
  T³−T−1. **Recomputed in quiverlab convention: [2,3,7] == Lehmer exactly;
  ρ ordering confirmed** (ρ_{[2,4,5]} ≈ 1.280638, ρ_{[2,3,8]} ≈ 1.230391).
- **§2.4:** A accessible ⇒ HH^i = 0 (i>0), HH⁰ = k. Interlaced towers:
  χ_{A_{s+1}} = (T+1)χ_{A_s} − Tχ_{A_{s−1}}.
- **Status:** fetched-source + recomputed.

### [dlP-cyclo] de la Peña (2013) — Cyclotomic Coxeter polynomials
- **Citation:** arXiv:1310.1557. Corroborates the tables; plus **§5.2 ladder
  R_{2n}** (2×n commutative grid): χ_{2n} = v_{n+1} ⊗ v₃ (root-product tensor);
  R₄ = v₃⊗v₃. "Cyclotomic type" ⇔ Mahler measure 1 — clean `mahler_measure`
  predicate. **Status:** fetched-source.

### [LMR-nakayama] Lenzing–Meltzer–Ruan (2022) — Nakayama algebras and Fuchsian singularities ★ strongest direct battery
- **Citation:** arXiv:2112.15587v2. §2.7 + Prop 6.1.
- **Algebra:** N_n(r) = linear equioriented Aₙ with rad^r = 0 (uniserial
  Nakayama).
- **Verbatim polynomials (Coxeter number in parens) — ALL六 RECOMPUTED AND
  MATCHED under quiverlab's convention:**
  χ(17,8) = (λ+1)(λ¹⁶+λ⁸+1) (24); χ(16,3) = (λ+1)(λ⁶−λ³+1)(λ⁹+1) (18);
  χ(15,6) = χ(15,4) = (λ+1)(λ⁸+λ⁴+1)(λ⁶+1) (12); χ(15,5) = (λ+1)(λ⁴+1)(λ⁵+1)² (40);
  χ(14,7) = (λ+1)(λ⁶−λ³+1)(λ⁷+1) (126). **Prop 6.1:** r≥9 ⇒ χ(r+7,r) =
  (λ+1)(λ⁶−λ³+1)(λ^r+1), Coxeter number lcm(2r,9).
- **Test:** NakayamaAlgebra/truncated linear Aₙ; exact polynomial equality +
  Coxeter number via Φ^m = I. **Status:** fetched-source + recomputed.

### [BdlP-dynkintype] Barot–de la Peña (1999) — Dynkin type of a non-negative unit form
- Expo. Math. 17 (1999) 339–348. Example 3.2: corank-5 D₁₁ unit form with
  Coxeter char. polynomial (T+1)⁶(T−1)⁶(T⁴+T³+T²+T+1) (degree 16). **Quiver
  UNVERIFIED** (bigraph not reconstructable from ASCII art) — polynomial firm,
  hold until the figure is read. **Status:** polynomial fetched; quiver UNVERIFIED.

### [GM-dkoszul] Green–Marcos (2008) + GMMZ (2004) — d-Koszul ★ feeds Plan 27
- **Citations:** arXiv:0812.3408v2; E. Green–E. Marcos–R. Martínez-Villa–P.
  Zhang, "D-Koszul algebras," J. Pure Appl. Algebra 193 (2004) 141–162 Thm 4.1
  (paywalled; verbatim-quoted via two fetched papers — secondary).
- **δ-pattern:** δ(n) = (n/2)d (even), ((n−1)/2)d+1 (odd); d-Koszul ⇔ P^n
  generated in degree δ(n); then Ext_Λ(Λ₀,Λ₀) generated in degrees **0,1,2**.
- **Ext of k[x]/(x^ℓ), ℓ≥3 (Keller/Madsen restatement, fetched):**
  Ext ≅ k[u,v]/(u²), |u|=1, |v|=2, each Ext^i 1-dim, independent of ℓ;
  internal bidegrees u=(1,1), v=(2,d) (the (2,d) half is synthesis from the
  δ-pattern — flagged). Matches Plan 27's k[y,z]/(y²) oracle: adversarially
  good (ranks identical across d, internal grading differs).
- **Status:** fetched-source + secondary (GMMZ Thm 4.1).

### [Cassidy-twodeg] Cassidy (2009) — Ext generated in two degrees (quadratic NON-Koszul witnesses!)
- **Citation:** arXiv:0903.0344v1 (k = CC, connected graded quadratic).
- **Thm 2.7:** for every m≥3 a QUADRATIC algebra C, gl.dim m, Yoneda algebra
  generated in bidegrees (1,1) and (m, m+1); Ext^{ij}=0 for i<j≤m,
  Ext^{m,m+1}≠0 — **m-Koszul in Backelin's sense but NOT Koszul**. Explicit
  presentation: 3m generators, 4+3m relations (m=3 → 10 gens/8 rels).
  **This supplies the quadratic-non-Koszul witness Plan 27's battery lacked**
  (if the presentation is graded-admissible for our engines — check
  finite-dimensionality when building; the agent recorded projective ranks).
- **Status:** fetched-source.

### [LS-stacked] Leader–Snashall (2015) — (D,A)-stacked algebras
- **Citation:** arXiv:1506.01854v2. Def 1.1, Thm 1.4, **Example 1.2**: explicit
  non-monomial (D=4,A=2)-stacked algebra — 7 vertices, 8 arrows (hexagon
  1→2→3→4→5→6→1 plus shortcut 1→7→3), relations
  ⟨(α₁α₂−α₇α₈)α₃α₄, α₃α₄α₅α₆, α₅α₆(α₁α₂−α₇α₈)⟩; minimal resolution
  P^n = e₁Λ⊕e₃Λ⊕e₅Λ (n≥2), generator degrees 0,1,then 2n; Ext algebra
  generated in degrees 0,1,2,3 (Thm 1.4). **Status:** fetched-source.

### [MV-preproj / EtingofEu / Hubery] Preprojective algebras
- **Citations:** Martínez-Villa CMS Conf. Proc. 18 (1996) 487–504 (primary
  UNVERIFIED — not fetchable); Etingof–Eu arXiv:math/0512287 (MRL 14, 2007)
  Thm 3.4.2/Prop 3.2.1; Hubery arXiv:2509.21448 (2025) Thm 7.1/7.5 — both fetched.
- **Non-Dynkin Q:** Π_Q Koszul, matrix Hilbert series (I − Ct + t²)^{−1}
  (C = adjacency of the double). **Dynkin Q:** Π_Q self-injective, Loewy length
  h−1, (h−2,2)-Koszul. **⚠ a WebFetch summary claimed "Koszul iff Dynkin" —
  FALSE**; the verified statements are as here.
- **Test:** PreprojectiveAlgebra: is_selfinjective for Dynkin + Loewy length
  h−1 (A₂→2, D₄→5, E₈→29 — derived, flag); Hilbert-series match to depth for
  non-Dynkin presentations if buildable. **Status:** fetched (EE, Hubery);
  MV primary unverified.

### [Schremmer-Happel] HH of canonical algebras ★
- **Citation:** F. Schremmer, arXiv:2512.08414 (2025), Prop 4.2.8/Cor 4.2.9,
  attributing Happel LNM 1404 / [Hap98 Thm 2.4].
- **Verbatim:** canonical C(p₁..p_t;λ), t≥3: dim HH^i = 1 (i=0), 0 (i=1),
  **t−3** (i=2), 0 else. HH_0 = #vertices = 2+Σ(p_i−1), HH_{≥1} = 0 (acyclic).
  Two-arm hereditary (a₁,a₂): HH¹ ∈ {3,2,1} by weight-1 count.
- **Test:** canonical (2,2,2,2) → HH²=1; (3,3,3),(2,4,4),(2,3,6) → HH²=0;
  (2,2,2,2,2) → HH²=2; cross-check Happel's Prop 1.2 trace identity against
  the degree-1 Coxeter coefficient — a TWO-INVARIANT consistency oracle.
- **Status:** fetched-source (Schremmer; primary Happel LNM 1404).

### [Happel-Cibils-Redondo] Hereditary + rad²=0 HH (via Redondo survey, Resenhas IME-USP 5 (2001), fetched)
- **Hereditary kQ (acyclic, any field):** HH⁰=1; HH¹ = 1 − n + Σ_{α∈Q₁}
  dim e_{t(α)}Ae_{s(α)}; HH^{>1}=0; HH¹=0 ⇔ tree; m-Kronecker HH¹ = m²−1.
  HH_0 = |Q₀|, HH_{≥1} = 0.
- **rad²=0 (Cibils Thm 3.1):** dim HH^n = |Q_n//Q₁| − |Q_{n−1}//Q₀| (n>0),
  HH⁰ = |Q₁//Q₀|+1. **⚠ n=1 discrepancy between Redondo Thm 4.11 (+1) and
  Cibils Thm 3.1 — treat Cibils as authoritative for n≥2 and verify
  HH⁰/HH¹ numerically before pinning.** k[x]/(x²): char≠2 → HH⁰=2, HH^{n>0}=1;
  char 2 → all 2 (agrees with the Cibils-cluster hand-verification). c-crown
  kQ/J² char-placement data (Prop 3.3).
- **Status:** fetched-source (survey); primaries Happel LNM 1404 / Cibils 1998.

### BibTeX (de la Peña / Lenzing / Marcos cluster)

```bibtex
@incollection{LenzingdlPena2008spectral,
  author = {Lenzing, Helmut and de la Pe\~na, Jos\'e Antonio},
  title = {Spectral analysis of finite dimensional algebras and singularities},
  booktitle = {Trends in Representation Theory of Algebras and Related Topics (ICRA XII)},
  series = {EMS Congr. Rep.}, publisher = {Eur. Math. Soc.}, pages = {541--588},
  year = {2008}, note = {arXiv:0805.1018} }
@article{LenzingdlPena2009chebysheff,
  author = {Lenzing, Helmut and de la Pe\~na, Jos\'e Antonio},
  title = {A {C}hebysheff recursion formula for {C}oxeter polynomials},
  journal = {Linear Algebra Appl.}, volume = {430}, number = {4},
  pages = {947--956}, year = {2009}, doi = {10.1016/j.laa.2008.10.003},
  note = {arXiv:math/0611535} }
@article{dlPena2014mahler,
  author = {de la Pe\~na, Jos\'e Antonio},
  title = {On the {M}ahler measure of the {C}oxeter polynomial of an algebra},
  journal = {Adv. Math.}, year = {2014}, note = {arXiv:1310.1910} }
@article{dlPena2013cyclotomic,
  author = {de la Pe\~na, Jos\'e Antonio},
  title = {Algebras whose {C}oxeter polynomials are products of cyclotomic polynomials},
  year = {2013}, note = {arXiv:1310.1557} }
@article{LenzingMeltzerRuan2022nakayama,
  author = {Lenzing, Helmut and Meltzer, Hagen and Ruan, Shiquan},
  title = {Nakayama algebras and {F}uchsian singularities},
  year = {2022}, note = {arXiv:2112.15587} }
@article{BarotdlPena1999dynkintype,
  author = {Barot, Michael and de la Pe\~na, Jos\'e Antonio},
  title = {The {D}ynkin type of a non-negative unit form},
  journal = {Expo. Math.}, volume = {17}, pages = {339--348}, year = {1999} }
@article{GreenMarcos2008dkoszul,
  author = {Green, Edward L. and Marcos, Eduardo N.},
  title = {$d$-{K}oszul algebras, 2-$d$-determined algebras and 2-$d$-{K}oszul algebras},
  year = {2008}, note = {arXiv:0812.3408} }
@article{GreenMarcosMartinezVillaZhang2004Dkoszul,
  author = {Green, Edward L. and Marcos, Eduardo N. and Mart\'inez-Villa, Roberto and Zhang, Pu},
  title = {{D}-{K}oszul algebras}, journal = {J. Pure Appl. Algebra},
  volume = {193}, pages = {141--162}, year = {2004} }
@article{Cassidy2009twodegrees,
  author = {Cassidy, Thomas},
  title = {Quadratic algebras with {E}xt algebras generated in two degrees},
  year = {2009}, note = {arXiv:0903.0344} }
@article{LeaderSnashall2015stacked,
  author = {Leader, Joanna and Snashall, Nicole},
  title = {The {E}xt algebra and a new generalisation of {D}-{K}oszul algebras},
  year = {2015}, note = {arXiv:1506.01854} }
@article{EtingofEu2007preprojective,
  author = {Etingof, Pavel and Eu, Ching-Hwa},
  title = {Koszulity and the {H}ilbert series of preprojective algebras},
  journal = {Math. Res. Lett.}, volume = {14}, year = {2007},
  note = {arXiv:math/0512287} }
@incollection{MartinezVilla1996preprojective,
  author = {Mart\'inez-Villa, Roberto},
  title = {Applications of {K}oszul algebras: the preprojective algebra},
  booktitle = {Representation Theory of Algebras}, series = {CMS Conf. Proc.},
  volume = {18}, pages = {487--504}, year = {1996} }
@book{Happel1989hochschild,
  author = {Happel, Dieter},
  title = {Hochschild cohomology of finite-dimensional algebras},
  series = {Lecture Notes in Math.}, volume = {1404}, publisher = {Springer},
  year = {1989} }
@article{Redondo2001survey,
  author = {Redondo, Mar\'ia Julia},
  title = {Hochschild cohomology: some methods for computations},
  journal = {Resenhas IME-USP}, volume = {5}, year = {2001} }
@article{Schremmer2025wpl,
  author = {Schremmer, Felix},
  title = {Weighted projective lines and {H}ochschild cohomology},
  year = {2025}, note = {arXiv:2512.08414} }
@article{Hubery2025preprojective,
  author = {Hubery, Andrew},
  title = {On the global dimension and {K}oszul property of preprojective algebras},
  year = {2025}, note = {arXiv:2509.21448} }
```

### Best 5 batteries (agent ranking)

1. **Nakayama N_n(r) Coxeter polynomials [LMR-nakayama]** — six exact
   polynomials + the Prop-6.1 family, already reproduced under quiverlab's
   convention; Coxeter number via Φ^m = I. Zero ambiguity, immediate.
2. **Dynkin + canonical Coxeter tables** — triple-corroborated + recomputed;
   the broadest coxeter_polynomial coverage (ship Dₙ v-factorization).
3. **Lehmer [2,3,7] spectral oracle** — flagship spectral_radius /
   mahler_measure pin + strict ρ ordering.
4. **Canonical-algebra HH² = t−3 + Happel trace identity** — two-invariant
   consistency (HH alternating sum = degree-1 Coxeter coefficient).
5. **Ext of k[x]/(x^d) = k[u,v]/(u²) + δ-pattern [GM-dkoszul]** — the Plan-27
   seed; pair with the char-sensitive k[x]/(x²) split for a GF(p)/char-0 probe.
   (Cassidy's quadratic-non-Koszul witnesses close Plan 27's missing-witness gap.)

**Honest-scope notes:** Dₙ cyclotomic condition mis-transcribed (use
v-factorization); Barot–de la Peña D₁₁ quiver unverified (figure); Martínez-
Villa 1996 verified only via reproving papers; k[x]/(x^n) per-degree HH beyond
x² should come from quiverlab's own closed form, not the graded-BV literature.

---

## Cluster report: Happel / Keller / Rickard

**Method:** every entry checked against a fetched source AND, where feasible,
**independently verified by running quiverlab itself** (pure kernel and
GF(32003)) plus **QPA crosschecks** for symmetry claims. Labels:
[ran quiverlab] / [ran QPA]. Springer LNM 1404 and the Happel-1997
ScienceDirect PDF are paywalled → those formulas are secondary-source (open
restatements, arXiv:2509.05135 Cor 4.10) and were re-derived + numerically
confirmed.

**⚠ HEADLINE (actionable, QPA-verified bug):** `is_symmetric` returns
**False** on provably symmetric algebras: (a) multi-vertex symmetric Nakayama
kZ_n/J^L with n | (L−1), n ≥ 2 (Brauer stars — kZ₂/J³, kZ₃/J⁴, kZ₄/J⁵ all QPA
`IsSymmetricAlgebra = true`, quiverlab False; weakly-symmetric detection via
the Nakayama permutation is CORRECT, the downstream "ν inner" search in
`invariants/frobenius.py::is_symmetric_generic` fails only for multi-vertex);
(b) `TrivialExtension(A)` for every A (always symmetric by the classical
theorem; root cause differs — TrivialExtension carries no quiver presentation,
so the path-basis certifier cannot run). Single-vertex k[x]/xᵃ and genuinely
non-symmetric kZ₃/J³ are both correct. Every passing `is_symmetric is True`
assertion in the current suite is single-vertex — an untested gap.

### [happel_hh_hereditary] Happel (1989, LNM 1404) — HH of hereditary kQ
- dim HH⁰ = 1; HH^{i≥2} = 0; **dim HH¹ = 1 − |Q₀| + Σ_{a∈Q₁}
  dim e_{s(a)}Ae_{t(a)}** (general form with relations: arXiv:2509.05135
  Cor 4.10). [ran quiverlab]: m-Kronecker → [1, m²−1, 0, 0] for m = 1,2,3;
  trees → [1,0,0]; the acyclic triangle (1→2, 2→3, 1→3) → [1,2,0].
  The `happel_question` → Happel1989 registry mapping CONFIRMED correct.
- **Status:** secondary-source + numerically confirmed.

### [happel_coxeter_trace] Happel (1997, LAA 258, 169–177) — Coxeter trace = HH Euler characteristic ★
- **tr Φ_A = Σ (−1)^i dim HH^i(A)** for finite gl.dim (Happel's convention).
  **Sign-matched to quiverlab [ran quiverlab]: with our Φ = −C^{−T}C the
  verified identity is `tr(coxeter_matrix) == −Σ(−1)^i dim HH^i`** —
  confirmed on 2/3-Kronecker, A₃, A₄, D₄, the acyclic triangle. PIN THE SIGN
  ON A₃ FIRST (quiverlab tr = −1).
- **Status:** secondary-source + sign-exact numerically confirmed.

### [rickard_brauer_derived] Rickard (1989/1991) — Brauer tree ≅_der Brauer star ★
- Every Brauer tree algebra (e edges, multiplicity m) is derived equivalent to
  the Brauer star = symmetric Nakayama kZ_e/J^{em+1} ⇒ equal HH^n / HH_n /
  HC_n. Star side [ran quiverlab + QPA]: kZ₂/J³, kZ₃/J⁴, kZ₄/J⁵, kZ₅/J⁶ all
  is_selfinjective ✓, QPA symmetric ✓ (quiverlab is_symmetric ✗ — the bug).
  Tree side: build-and-compare recommended (special-biserial relations,
  CS engine), NOT hard-coded relations.
- **Status:** fetched-source + star side verified.

### [keller_derived_invariance] Keller (JPAA 123, 1998) / Rickard (1991) — HH*/HH_*/HC_* derived invariants ★ oracle scheme
- Derived equivalence ⇒ equal HH^n, HH_n, HC_n. **Free self-checking oracle:**
  reflection-equivalent orientations of one graph. [ran quiverlab]: both
  acyclic orientations of the à ₂ triangle give HH* = [1,2,0], HH_* = [3,0,0],
  HC_* = [3,0,3] — matched exactly.
- **Status:** fetched-source + numerically confirmed.

### [armenta_keller_cap] Armenta & Keller (2019, C. R. Acad. Sci.; arXiv:1711.02947) — derived invariance of the CAP product
- The HH*-module structure on HH_* (cap) is a derived invariant — directly
  relevant to quiverlab's native cap (Plans 20/21): compare the induced action
  across a derived-equivalent pair. **Status:** abstract/metadata (statement
  corroborated twice; theorem numbers unconfirmed).

### [happel_trivial_extension] T(A) is symmetric (Happel LMS 119, 1988; Yamagata; ASS2006)
- T(A) = A ⋉ D(A) symmetric for EVERY f.d. A. [ran quiverlab]:
  is_frobenius ✓ but is_symmetric ✗ (no quiver presentation — see headline).
  Suggested fix: give TrivialExtension a double-quiver presentation.
- **Status:** fetched-source theorem + discrepancy reproduced.

### [happel_question_counterexample] BGMS (MRL 12, 2005) / Bergh–Erdmann (2008)
- The QCI A_q (q not a root of unity) has infinite gl.dim with bounded HH^• —
  the negative answer to Happel's question (restated arXiv:2509.05135 Ex 6.2).
  Already wired (`BGMS2005`, `qci_hh_oracle`); ties to the Solotar-cluster
  general-(a,b) extension. **Status:** secondary + already-pinned.

### [euler_char_preprojective] Graded-Euler-characteristic route (arXiv:2606.26255, 2606.15595, 2607.10913)
- Graded refinement of the Happel trace identity + candidate HH/HC tables for
  higher preprojective algebras. **Status: abstract-only — follow-up fetch
  needed before any pin.**

### [happel_ringel_tilted] Happel–Ringel (1982) / Happel–Vossieck (1983) — UNVERIFIED lead
- Canonical sources for concrete tilted-algebra derived pairs (nonzero-HH
  instances of the derived-invariance scheme) and τ/AR data. No machine-usable
  worked example fetched — next-fetch recommendation recorded.

### BibTeX (Happel/Keller/Rickard cluster)

```bibtex
@article{Happel1997,
  author = {Happel, Dieter},
  title = {The trace of the {C}oxeter matrix and {H}ochschild cohomology},
  journal = {Linear Algebra and its Applications}, volume = {258},
  year = {1997}, pages = {169--177}, doi = {10.1016/S0024-3795(96)00195-4} }
@book{Happel1988,
  author = {Happel, Dieter},
  title = {Triangulated Categories in the Representation Theory of Finite Dimensional Algebras},
  series = {London Mathematical Society Lecture Note Series}, volume = {119},
  publisher = {Cambridge University Press}, year = {1988},
  doi = {10.1017/CBO9780511629228} }
@article{Rickard1989stable,
  author = {Rickard, Jeremy}, title = {Derived categories and stable equivalence},
  journal = {Journal of Pure and Applied Algebra}, volume = {61}, number = {3},
  year = {1989}, pages = {303--317}, doi = {10.1016/0022-4049(89)90081-9} }
@article{Rickard1991,
  author = {Rickard, Jeremy}, title = {Derived equivalences as derived functors},
  journal = {Journal of the London Mathematical Society}, volume = {43},
  number = {1}, year = {1991}, pages = {37--48}, doi = {10.1112/jlms/s2-43.1.37} }
@article{Keller1998cyclic,
  author = {Keller, Bernhard},
  title = {Invariance and localization for cyclic homology of {DG} algebras},
  journal = {Journal of Pure and Applied Algebra}, volume = {123},
  number = {1-3}, year = {1998}, pages = {223--273},
  doi = {10.1016/S0022-4049(96)00085-0} }
@article{ArmentaKeller2019,
  author = {Armenta, Marco and Keller, Bernhard},
  title = {Derived invariance of the cap product in {H}ochschild theory},
  journal = {Comptes Rendus Math\'ematique}, year = {2019},
  note = {arXiv:1711.02947} }
@misc{HappelQuestionTau2025,
  title = {Happel's question, {H}an's conjecture and $\tau$-{H}ochschild (co)homology},
  year = {2025}, note = {arXiv:2509.05135} }
@misc{GradedEulerPreproj2026,
  title = {Hochschild (co)homology and cyclic homology via a graded Euler characteristic with applications to higher preprojective algebras},
  year = {2026}, note = {arXiv:2606.26255} }
```

### Best 5 batteries (agent ranking)

1. **Happel-1997 trace identity** — tr(coxeter_matrix) == −Σ(−1)^i dim HH^i:
   cheap, zoo-wide, both sides independent, catches Cartan/Coxeter AND HH bugs.
2. **Derived-invariance across orientations** — self-checking equal
   HH*/HH_*/HC_* (verified on the à ₂ pair).
3. **Hereditary formula** — HH¹ = m²−1 on m-Kronecker + trees [1,0,0].
4. **Brauer-star symmetry oracle** — reproduces the LIVE is_symmetric bug;
   QPA-crosschecked regression battery.
5. **BGMS Happel-question falsifier** — bounded HH^• + infinite gl.dim on
   QuantumCI, tied to the existing qci_hh_oracle.

---

# Consolidated: implementation priorities across all four clusters

Cross-cluster headliners (each battery must cite its bib entry; honest-scope
flags live in the cluster sections):

1. **The is_symmetric bug fix + regression battery** (Happel/Rickard cluster —
   QPA-verified live bug; multi-vertex ν-inner sweep + TrivialExtension
   presentation).
2. **Nakayama Coxeter polynomials** (LMR 2022 — six exact polynomials already
   recomputed under quiverlab's convention + Prop-6.1 family + Coxeter numbers).
3. **Happel-1997 trace identity, zoo-wide** (sign pinned on A₃).
4. **Dynkin/affine/canonical Coxeter tables + Lehmer [2,3,7] spectral pin**
   (ship Dₙ v-factorization; mahler_measure == 1 cyclotomic-type predicate).
5. **Derived-invariance orientation pairs** (HH*/HH_*/HC_*; later the cap
   module structure per Armenta–Keller).
6. **BE QCI cohomology [2,2,1,0,…] all a,b + general-(a,b) homology** (char-0
   only — the GF(p) root-of-unity trap is documented).
7. **Rad²=0 parallel-path HH with the char-2 doubling; triangular-string
   family with the degree-(2m+1) revival; a-Kronecker [1, a²−1, 0…];
   cup-triviality (RR2014) vs cup-nonvanishing (RR2018).**
8. **Taft/Nakayama cyclic homology (char 0)** — first strong HC family oracle.
9. **Canonical-algebra HH² = t−3 cross-linked to the trace identity**
   (two-invariant consistency).
10. **Incidence ≅ nerve cohomology; trivial-extension HH¹ decomposition +
    Example 2.20; truncated finiteness boolean (dim HH• < ∞ ⇔ acyclic).**
11. **Plan-27 feeders:** Ext of k[x]/(x^d) = k[u,v]/(u²) with internal
    bidegrees (1,1),(2,d); Cassidy quadratic-non-Koszul witnesses (gl.dim m,
    generators (1,1)+(m,m+1)) — closes the missing-witness gap; (D,A)-stacked
    Example 1.2; preprojective non-Dynkin Koszul + Dynkin self-injective/LL.
