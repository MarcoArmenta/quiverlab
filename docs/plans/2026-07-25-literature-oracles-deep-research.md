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
- Happel / Keller / Rickard — pending.
- de la Peña / Lenzing / Marcos — pending.

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
