# Computability-expansion deep research — the road to 1.0.0

**Date:** 2026-08-06. **Method:** four author-cluster research scouts (web-verified,
Fable) each shadowed by an adversarial verification critic; all findings below are
**post-adjudication** — every correction the critics established is already applied
to the records. Clusters: (A) Solotar–Suárez-Álvarez (Buenos Aires),
(B) Cibils–Lanzilotta–Marcos–Volkov, (C) Assem–Skowroński–de la Peña,
(D) Keller–Witherspoon–Brüstle–Paquette — plus their coauthor networks
(Redondo, Chouhy, Herscovich, Chaparro, Rossi Bertone, Artenstein, Mata, Barrios,
Le Meur, Saorín, Trepode, Coelho, Liu, Bustamante, Bobiński, Negron, Nguyen,
Shepler, Ocal, Hanson, Tattar, Schroll, Buan, Marsh, Yıldırım, Rock, Mousavand,
Igusa, Todorov, …).

**Verification outcome:** ~90 citations spot-checked by the critics; **zero
fabricated papers**. Corrections applied below: one inverted π₁ oracle (C-scout),
one mis-stated string-algebra criterion (A-scout), a Volkov↔Lambre–Zhou–Zimmermann
reference swap (B-scout), brick-finiteness gates on three τ-tilting-infinite
searches (D-scout), a Plan-40 redundancy re-scope (B-scout), and char-0 /
algebraically-closed-field scope boundaries throughout.

**Release gate (Marco, 2026-08-06):** version **1.0.0 ships when this program is
incorporated** through the full spec → plan → implement → verify cycle. Plans
P51+ are minted from the records below.

**House rules binding on every item:** exact only (floats fail loudly); every
honest-scope boundary below becomes a loud refusal, never a silent gap; every
oracle lands on the verification page; oracle-class markers per Plan 32.

---

## 0. Cross-cluster headlines (double-corroborated)

**H1 — The Gerstenhaber bracket past the bar window** is the single most
strategically valuable item. Two clusters independently converged on it, and the
A-critic verified in code that quiverlab's bracket is bar-window-transported (the
"already have it" dismissal was wrong). Method: Volkov's arbitrary-resolution
formulas (arXiv:1610.05741) / Negron–Witherspoon homotopy liftings (Oke
arXiv:2103.12331; Witherspoon GSM 204). **Prerequisite the critics pinned:** the
resolution must carry a diagonal or a bar-comparison — quiverlab ships BOTH
(Plan-14 comparison, Plan-20 diagonal Δ on CS). Scope: the *new* content is
past-window and off-GF(p); in-window transported bracket is the cross-engine
anchor. → Record R1.

**H2 — Fractional Calabi–Yau dimension of self-injective algebras** — found by
the D-critic as the scout's strongest omission: computable from syzygy/Nakayama
periods, every prerequisite shipped (Π(A_n) known fractionally CY, Frobenius
certifiers, Nakayama automorphism). Low-hanging flagship. → Record R24.

**Deduped:** π₁ appeared in clusters B and C (merged, R11, with the corrected
square oracle); skew-group algebras appeared in B and D (merged, R8).

---

## 1. Hochschild / Tamarkin–Tsygan deep layer

**R1 — Bracket beyond the bar window (homotopy liftings).** [D-scout P1 + A-scout
exclusion reversed; verdict keep-with-corrections]
Object: `[-,-]: HH^p×HH^q → HH^{p+q-1}` natively on the CS/minimal/Bardzell
resolution via homotopy-lifting maps ψ solving a finite linear system (the
`_d_general`-style solve + `reduce_mod_nullspace` canonicalization). Refs:
Volkov arXiv:1610.05741 (Proc. Edinb. Math. Soc.); Oke arXiv:2103.12331;
Witherspoon GSM 204; Negron–Witherspoon A∞-coderivations. Scope: needs the
diagonal/comparison (shipped); new = past-window + any exact field. Oracles:
in-window ≡ transported bracket (cross-engine); antisymmetry+Jacobi (self-cert);
the 2103.12331 Koszul-quiver tables; k[x]/(x^n) zero-entries. Size L. Deps:
resolutions_cs/diagonal.py, hochschild/products.py.

**R2 — BV operator Δ on HH of Frobenius/self-injective algebras.** [B-scout P4,
reference-rebuilt per critic; deepens the backlogged BV item]
Object: Δ: HH^n → HH^{n-1} by transporting Connes B through the σ-twisted
Frobenius duality (design sketch Δ = ∂∘B_σ∘∂⁻¹ — a paraphrase, not a verbatim
formula); bracket recovered from Δ (7-term BV axioms self-certified against the
independent bracket). **Corrected refs:** Lambre–Zhou–Zimmermann,
arXiv:1405.5325, *J. Algebra* 2016 (hypothesis: SEMISIMPLE Nakayama automorphism);
Volkov, arXiv:**1405.5155** (hypothesis: ord(ν) coprime to char k — a DIFFERENT
computable condition; implementation must pick deliberately);
Bian–Itagaki–Kou–Lyu–Zhou arXiv:2603.04834 (2026; self-injective Nakayama,
semisimplicity hypothesis removed — the explicit kZ₁/J^N Δ/bracket pins MUST be
transcribed verbatim from the paper body with equation numbers before any test
uses them; the scout's quoted formulas are unverified). Oracles: Δ²=0,
bracket-from-Δ ≡ independent bracket, symmetric case = Tradler on
k[x,y]/(x²,y²) (HH_• = [4,4,5,6] pin exists), QuantumCI ties. Size M. Deps:
connes_differentials, nakayama_automorphism, is_frobenius, cup/bracket.

**R3 — Tate–Hochschild (singular HH) for self-injective / eventually periodic
Gorenstein.** [D-scout P3; keep]
Object: ĤH^*(A) in negative and positive degrees with cup product, via complete
resolutions (splice minimal A^e with its dual); eventual-periodicity certificate
(invertible homogeneous element). Refs: Keller arXiv:1809.05121; Usui
arXiv:2107.03326; Bergh–Jorgensen arXiv:1109.4019. Scope: eventually-periodic
Gorenstein only (the class where one period makes it finite) — loud refusal
otherwise. Oracles: positive degrees ≡ ordinary HH; k[x]/(x^n) periodic Tate
ring; symmetric Nakayama zoo periods; Tate ≡ ordinary above findim. Size M–L.

**R4 — HH with arbitrary bimodule coefficients + relative HH•(A|B).** [A-scout
P8; keep — verified genuinely absent from the public surface]
Object: HH•(A,M), HH_•(A,M) for any f.d. bimodule M (twisted _1A_σ, D(A),
A/soc…); relative HH over a subalgebra B via the relative bar resolution. Refs:
Chaparro–Schroll–Solotar arXiv:1811.02211 (gentle HH¹ with coefficients);
Lindell–Rubio y Degrassi arXiv:2411.03080. Shallow engine extension (swap the
coefficient in Hom_{A^e}(P•, −)); prerequisite for R2/R8/R9. Size M.

**R5 — Split-algebra / trivial-extension HH long exact sequence.** [B-scout P7;
keep]
Object: HH^*(B⋉M) assembled from HH^*(B) + H^*(B,M) via the LES with
cup-product connecting map; theorem HH¹(B⋉M) ≠ 0 for M ≠ 0 (grading-derivation
witness; state M ≠ 0). Ref: Cibils–Marcos–Redondo–Solotar arXiv:math/0102194.
Turns the existing brute T(B) computation into a structured, self-verifying
decomposition. Oracles: LES-vs-direct degreewise on T(kA_n)/T(kD₄) (in suite);
HH¹≠0 across the zoo. Size M. Deps: R4.

**R6 — Arrow removal/addition HH reductions (Han accelerator).** [A-critic
promotion of a demoted item]
Object: certified finite HH reductions along deleting/adding arrows of kQ/I.
Ref: Cibils–Lanzilotta–Marcos–Solotar arXiv:1812.07655. Pairs with R7's
Jacobi–Zariski machinery; feeds the backlogged Han's campaigns. Size M.

**R7 — Han's-conjecture bounded-extension recognizer + Jacobi–Zariski
sequence.** [B-scout P3; keep-with-corrections]
Object: decide bounded extension B ⊆ A (tensor-nilpotency of A/B over B — an
HONEST CAPPED SEMI-DECISION, powers can grow before dying, no a-priori index
bound; finite pd of A/B over B^e; one-sided projectivity), then transport Han
B ⊨ Han ⇔ A ⊨ Han and expose the Jacobi–Zariski nearly-exact sequence ("exact
twice in three", per arXiv:1908.11130). Refs: CLMS arXiv:2101.02597 (J. Algebra
2022; J-interrupter Ex. 5.3/5.5 as oracles incl. the negative case);
arXiv:1908.11130; arXiv:2409.00945 (recollement approach — AUTHOR ATTRIBUTION
UNVERIFIED, re-check at plan time). Size M–L.

**R8 — Skew group / smash algebras A⋊G as inputs, with HH transfer.** [merged
B-P8 + D-P2; keep]
Object: constructor (A, G-action) → A⋊G (dim |G|·dim A) as a first-class
Algebra; HH via the Ştefan conjugacy-class decomposition
(⊕_{[g]} HH^n(A, {}_gA))^{Z(g)} — needs R4. HARD SCOPE: char k ∤ |G| (Reynolds
averaging); the action is explicit input (never inferred); free-action detection
on the quiver is a finite check enabling the Galois-covering reduction. Refs:
Shepler–Witherspoon (Adv. Math., people.tamu.edu/~sjw/pub/ring.pdf);
arXiv:0911.0938; Cibils–Marcos arXiv:math/0312214 (Proc. AMS 134 (2006) 39–50);
CMRS arXiv:1804.02223. Oracles: trivial G byte-identity; dim law; Z/2 on dual
numbers direct-vs-decomposition; Nakayama as orbit algebras. Size L.

**R9 — Incidence-algebra HH = simplicial cohomology (fast path + oracle
family).** [B-scout P6; keep-with-corrections]
Object: for incidence algebras kP, HH^* ≅ simplicial H^* of the order complex —
CUP-PRODUCT/ring isomorphism (Cibils JPAA 56 (1989) 221–232;
Gerstenhaber–Schack JPAA 30 (1983) 143–156). The "bracket = 0" claim is NOT in
G–S '83 — dropped unless re-sourced (math/0611542 / arXiv:2411.07910 line).
Poset input mode + SNF-fast HH + a second independent oracle for the general
engines (extends the existing B₃ pin to the theorem). Size S–M.

**R10 — GHMS Koszul minimal bimodule resolution (fast HH for Koszul
algebras).** [B-scout P9; keep]
Object: the comultiplicative minimal A^e-resolution of a Koszul algebra built
from the Koszul dual data; fast minimal HH + a third HH oracle class. Ref:
Green–Hartman–Marcos–Solberg, Arch. Math. 85 (2005) 118–127; arXiv:math/0508177.
Oracles: exterior/preprojective HH degreewise vs existing engines; Betti =
Koszul-dual Hilbert coefficients. Size M. Deps: Plan-27 Koszulity + Ext-algebra.

## 2. The Lie / deformation layer on HH

**R11 — Lie-algebra classification of HH¹.** [A-scout P1; keep]
Object: HH¹ = Der/Inn with commutator — computable over ANY exact field
independent of the window-bounded bracket engine (critic-verified). Classify:
derived/lower-central series (any field); radical + Levi + sl₂-count m + toral
rank (HARD char-0 boundary — Levi/Killing fail in char p, and that boundary is
itself mathematically interesting: k[x]/(x^n) gives Jacobson–Witt W₁ when
char|n). Refs: Rubio y Degrassi–Schroll–Solotar arXiv:1903.12145 (Quaest. Math.
2023); Chaparro–Schroll–Solotar arXiv:1811.02211 (J. Algebra 558 (2020));
Eisele–Raedschelders arXiv:1903.07380 (hypothesis: NON-WILD, algebraically
closed, sl₂ needs char ≠ 2); Strametz JAA 5(3) (2006) 245–270; Liu–Xing
arXiv:2306.14372. Oracles: k[x]/(x^n) solvable (char∤n) vs W₁ (char|n);
T(Kronecker) HH¹ ≅ sl₂ (char≠2); the RSS Ext-quiver solvability criterion
(Ext¹(S,S)=0 ∀S and dim Ext¹(S,T)≤1 ⇒ solvable) as a pure quiver test. Size M.

**R12 — HH• as a graded Lie module over HH¹.** [A-scout P2;
keep-with-corrections]
Object: the bracket action [HH¹, HH^n] (degree-1 case = Lie derivative of a
derivation — implementable field-generally WITHOUT the full bracket engine;
state which method), weight/torus decomposition, indecomposable Lie-module
summands (reuse modules/decompose.py; char-0 cleanest). Refs:
Artenstein–Lanzilotta–Solotar arXiv:1803.10310;
Meinel–Nguyen–Pauwels–Redondo–Solotar arXiv:1803.10909 (J. Algebra 580 (2021);
Virasoro-subquotient statement verified verbatim);
Chaparro–Schroll–Solotar–Suárez-Álvarez arXiv:2311.08003 (J. Algebra 708 (2026)
138–231 — journal ref independently confirmed). Size M–L. Deps: R11.

**R13 — L∞ / Maurer–Cartan deformations on Bardzell's complex.** [A-scout P3;
keep-with-corrections]
Object: the ℓₙ brackets on B(A)[1] for monomial/gentle A; MC set = formal
deformations; obstruction [α,α] ∈ HH³; the PRESENTED deformed algebra A_α + its
Ext-algebra (feed back into the engine). Honest scope: char 0; "ℓₙ=0 for n≥5"
is a SUFFICIENT condition and ℓ₄ can be nonzero (honest L∞, not DGLA); MC =
2-cocycles holds in the NILPOTENT regime (gentle, no parallel arrows, no
oriented cycles — Thm 5.4). COST FLAG (critic): quiverlab's Bardzell engine is
GF(p) int64 — this needs a char-0 Bardzell path (CS runs over any Domain and is
the fallback). Refs: Redondo–Rossi Bertone arXiv:2008.08122 (JPAA 226(5) 2022);
Müller–Redondo–Rossi Bertone–Suarez arXiv:2309.02582 (Comm. Alg. 2025);
Redondo–Román–Rossi Bertone–Verdecchia arXiv:2003.10366 (Morita invariance);
Redondo–Román–Rossi Bertone arXiv:2202.01199 (Ext-algebra of deformations —
THREE authors, no Verdecchia); Chouhy arXiv:1708.02933 (degeneration link).
Size L.

## 3. Coverings, fundamental groups, structural recognizers

**R14 — Fundamental group π₁(Q,I) of a presentation + Hurewicz to HH¹.**
[merged B-P5 + C-P1; adjudicated corrections applied]
Object: the COMBINATORIAL presentation group π₁(Q,I) (walks mod
relation-homotopy along a spanning tree) — finite presentation always
computable; abelianization = universal abelian grading group via exact Smith
normal form; Hom(π₁, k⁺) ↪ HH¹ cross-check. HONEST DECIDABILITY (adjudicated):
triviality of an f.p. group is undecidable (Adian–Rabin) — emit π₁^{ab} always,
"simply connected" only via decidable sufficient criteria (separation condition
R16; Le Meur's privileged presentation: char 0 + no double bypasses). The
INTRINSIC π₁ (inverse limit over connected gradings) is NOT bounded-computable
in general — refuse loudly; the char-p oracle π₁(k[x]/(x^p)) = ℤ × C_p
(verified verbatim in arXiv:0906.3069) is the INTRINSIC group, not the
presentation group (a monomial loop relation gives presentation π₁ = ℤ) — use
it only if the grading route is implemented. **CORRECTED ORACLE:** commutative
square WITH the commutativity relation ⇒ π₁ = 1 (the relation identifies the
parallel paths); ℤ only WITHOUT it. Refs: Assem–de la Peña Comm. Alg. 24 (1996)
187–208; Cibils–Redondo–Solotar arXiv:0706.2491, arXiv:1010.6296,
arXiv:0906.3069; Le Meur arXiv:math/0503302; Briggs–Rubio y Degrassi
arXiv:2109.03704 (IMRN 2023: every maximal torus of HH¹ is dual to some π₁ —
verified). Size M.

**R15 — Left/right parts L_A, R_A + support algebras (the substrate).**
[C-scout P2; keep — prerequisite for R17/R18/R21]
Object: L_A/R_A by the closed-under-predecessors pd/id sweep on the knitted AR
quiver (rep-finite scope, loud otherwise); Ext-injectives of add L_A; left/right
support algebras A_λ, A_ρ (theorem: products of tilted algebras — verified by
R17). Refs: Assem–Coelho–Trepode J. Algebra 281 (2004) 518–534;
Assem–Castonguay–Lanzilotta–Vargas arXiv:1102.1188; the "Organising the module
category" survey (São Paulo J. Math. Sci. 2021). Size M.

**R16 — Strongly-simply-connected recognizer (separation condition).** [C-scout
P9; keep]
Object: for triangular A: every convex subcategory satisfies separation (rad P_x
decomposes with supports in distinct components of the non-predecessor
subquiver) — pure finite combinatorics + Krull–Schmidt; witness on failure.
Gate for R19. Refs: Skowroński CMS Conf. Proc. 14 (1993); ASS book;
arXiv:1905.06028. Size M.

**R17 — Tilted-algebra recognizer (Liu–Skowroński faithful section).** [C-scout
P5; keep-with-corrections]
Object: find/refute a faithful section Σ with Hom(X, τY) = 0 (|Σ₀| = n is
IMPLIED, state as consequence); returns slice, tilting module, hereditary type.
Rep-finite via exhaustive section enumeration; the local criterion
(arXiv:1409.2054) is the rep-infinite extension path. Refs: Liu, Arch. Math. 61
(1993) 12–19 (venue verified); Skowroński; arXiv:1409.2054. Oracles: hereditary
⇒ tilted; kZ₃/J² stable tube ⇒ not; cluster-tilted A₃ ⇒ not. Size M.

**R18 — Quasi-tilted / shod / weakly-shod / laura / ada recognizer ladder.**
[C-scout P3+P4+P11; keep, sequenced after R15]
Object: the nested per-instance certificates — quasi-tilted (gl.dim ≤ 2 and
every indec pd≤1 or id≤1), shod, weakly shod (bounded I⇝P paths), laura
(L_A ∪ R_A cofinite; report the finite complement), ada (all projectives+
injectives in L_A ∪ R_A) — and for ada the THEOREM "simply connected ⇔ HH¹ = 0"
(verified in arXiv:1102.1188, algebraically closed k) turning quiverlab's HH¹
into a complete simple-connectedness oracle on that class. Rep-finite scope.
Refs: Happel–Reiten–Smalø Mem. AMS 575; Coelho–Lanzilotta J. Algebra 265
(2003); Assem–Coelho J. Algebra 269 (2003); Smith arXiv:math/0702562;
Bordino–Fernández–Trepode arXiv:1404.5294 (Comm. Alg. 2017). Size M total.

**R19 — Tame/wild certificate via the Tits form (strongly simply connected).**
[C-scout P7; keep-with-corrections]
Object: weak positivity / weak nonnegativity of q_A decided by EXACT bounded
integer search (primary route; the printed hypercritical/critical lists —
de la Peña Banach Center 26 (1990), von Höhne, Unger, and the
Barot–Jiménez-González–de la Peña 2019 book — as transcription-checked
cross-oracles); verdicts: rep-finite ⇔ weakly positive, tame ⇔ weakly
nonnegative. SCOPE: strongly simply connected (gate = R16) over ALGEBRAICALLY
CLOSED fields — the verdict layer is CC-only, the form computations are
field-free. Companion (critic's top find): **Bongartz's criterion** — Math. Ann.
269 (1984) 1–12 — rep-finite ⇔ weakly positive for simply connected, same
machinery. Refs: Brüstle–de la Peña–Skowroński Adv. Math. 226 (2011) 887–951;
arXiv:1905.06028. Oracles: m-Kronecker ladder (1/2/≥3), Dynkin/Euclidean trees,
T_{2,3,7}. Size L.

**R20 — Coxeter spectral analysis.** [C-scout P8; keep-with-corrections]
Object: cyclotomic test (exact ℤ[x] factorization + Φ_n recognition); Mahler
measure and spectral radius as CERTIFIED ALGEBRAIC NUMBERS (minimal polynomial +
rational isolating interval — never a float; sympy CRootOf/Sturm). The Lehmer
dichotomy (M = 1 or ≥ μ₀ ≈ 1.17628) is proven ONLY for restricted classes
(accessible algebras) — an open problem in general; report per-class only.
Refs: de la Peña arXiv:1310.1910, arXiv:1310.1557 (both sole-author, verified);
de la Peña–Takane Arch. Math. 55 (1990) 120–134. Oracles: kA₂ χ = Φ₃;
3-Kronecker χ = x²−7x+1, ρ = (7+3√5)/2; T_{2,3,7} = E₁₀ realizes Lehmer's
degree-10 polynomial (verified). Size S–M.

**R21 — Liu degree theory + AR-component invariants; representation-directed
recognizer.** [C-scout P10+P12; keep]
Object: left/right degrees of irreducible maps, sectional paths,
postprojective/preinjective/regular partition, generalized-standard flags;
directing modules and rep-directed (Γ_A acyclic) — all finite sweeps on the
knitted AR quiver (rep-finite). Refs: Liu JLMS 45 (1992) 32–54, JLMS 47 (1993)
405–416; Ringel LNM 1099; Bongartz CMH 57 (1982). Size M.

**R22 — Cluster-tilted algebras: relation-extension constructor + local-slice
recognizer.** [C-scout P6; keep]
Object: Ĉ = C ⋉ Ext²_C(DC, C) (reuse the Plan-31 ⋉ machinery + Ext engine —
distinct from T(A) = A ⋉ DA); recognizer by local-slice search in Γ_B with
End-reconstruction and Ĉ ≅ B check. Refs: Assem–Brüstle–Schiffler
arXiv:math/0601537 (Bull. LMS 40 (2008) 151–162), arXiv:0707.0038 (J. Algebra
319 (2008) 3464–3479). Oracles: A₃-with-zero-relation → the 3-cycle
cluster-tilted algebra; hereditary C ⇒ Ĉ = C; constructor∘recognizer
round-trip. Size M.

## 4. Homological dimensions (extends Plan 40 — re-scoped per critic)

**R23 — φdim/ψdim as ALGEBRA invariants + the φ-spectrum + LIT certificates.**
[B-scout P1+P2, MAJOR re-scope: quiverlab ALREADY SHIPS φ(M)/ψ(M), the
truncated closed forms, φ=pd, and the self-injective ⇔ φ≡0 oracles (Plan 40)]
Surviving new objects: (a) φdim(A)/ψdim(A) — rep-finite exhaustive sup over the
knitted indecomposables (monotonicity on add makes ⊕-of-all sufficient);
certified lower bounds otherwise; the chain findim ≤ φdim ≤ ψdim ≤ gldim as a
standing self-cert; (b) the φ-spectrum/gaps (Barrios–Mata–Rama
arXiv:1810.12112); (c) LIT-algebra certificates: no known decision procedure in
general (NOT "proven undecidable" — reworded per critic); the decidable
certificate families (φdim ≤ 1, Gorenstein with 𝒟 = Gproj, self-injective,
id(A_A) < ∞) emit a proof-carrying findim bound ψ_𝒟(V) + n + 1. Corrected refs:
Fernandes–Lanzilotta–Mendoza arXiv:1304.0754 (φdim paper — NOT Huard);
Barrios–Lanzilotta–Mata survey arXiv:2310.09283; LIT definitions in
Bravo–Lanzilotta–Mendoza–Vivero arXiv:2002.07866 (JPAA) — arXiv:2105.06273 is
Barrios–Mata "On Lat-Igusa-Todorov algebras"; arXiv:2103.12120, 2311.06148.
Size S–M on top of Plan 40.

**R24 — Fractional Calabi–Yau dimension of self-injective algebras.** [D-critic
omission-find; new record]
Object: the (m, n) with Ω^{m}-shifted Nakayama-twist periodicity certifying A
fractionally CY of dimension m/n in the stable category; every prerequisite
shipped (syzygy periods, ν, Frobenius certifiers, Π(A_n) known fractional CY
values as oracles). Refs to be pinned at plan time from the Keller-school
literature (e.g. the standard fractionally-CY computations for preprojective
and Nakayama algebras). Size S–M.

## 5. τ-tilting, lattices, stability

**R25 — Wall-and-chamber for all ranks via bricks.** [D-scout P7;
keep-with-corrections]
Object: walls D(B) from brick enumeration (exact linear inequalities over
submodule dim-vectors), chambers = g-vector cones. GATE (critic): certified
complete iff τ-tilting-finite/brick-finite — decidable; otherwise a bounded
region with honest truncation. Refs: Brüstle–Smith–Treffinger arXiv:1805.01880
(Adv. Math.); Kaipel–Treffinger arXiv:2302.12699; Asai arXiv:1610.05860.
Oracles: chamber count = # support τ-tilting (cross-engine); kA₂ = 5 chambers /
3 walls (verified — the P₁ wall is a ray); MGS = green paths. Size M.

**R26 — Torsion-lattice congruences, canonical joins, core label order = wide
subcategory poset.** [D-scout P8, reframed: brick LABELS already ship (P45)]
New content only: Con(tors A), the forcing order on bricks, canonical join
representations, core label order ≅ wide subcategories. Refs:
Demonet–Iyama–Reading–Reiten–Thomas arXiv:1711.01785 (Trans. AMS B);
Barnard–Carroll–Zhu. Oracles: join-irreducibles ↔ bricks (count cross-engine);
kA₂ pentagon; kA₃ congruence lattice. Size S–M.

**R27 — τ-exceptional sequences + mutation.** [D-scout P5; keep-with-note]
Object: (signed) τ-exceptional sequences via τ-perpendicular recursion (Jasso
reduction on the shipped τ-tilting engine); bijection with ordered support
τ-tilting. NOTE (critic): mutation transitivity for enumeration proven only in
rank 2 (arXiv:2402.10301) — enumeration for rank ≥ 3 via the bijection, not via
mutation-BFS. Refs: Buan–Marsh J. Algebra 585 (2021) 36–68; arXiv:2211.10428;
Buan–Hanson–Marsh arXiv:2402.10301; Nakayama counting paper
(s10468-021-10060-y) as oracle. Size M.

**R28 — Classical exceptional sequences (hereditary): braid action,
enumeration, c-matrices.** [D-scout P6; keep-with-corrections]
Object: recognizer (Hom/Ext orthogonality), braid mutation σ_i (universal
extension/kernel constructions), enumeration by braid-orbit BFS (transitive for
hereditary — Crawley-Boevey; Ringel), Dynkin closed-form counts as oracles;
c-matrix reading (overlaps existing τ-tilting c-vectors — hereditary-only
scope). CORRECTED refs: cite Crawley-Boevey "Exceptional sequences of
representations of quivers" (Ottawa 1992 proceedings) and Ringel (CMS Conf.
Proc. 14) directly — arXiv:2102.04584 is Alvares–Marcos–Meltzer on weighted
projective lines (mislabeled by the scout, dropped as anchor);
Garver–Igusa–Matherne–Ostroff arXiv:1506.08927; the Dynkin count paper (eudml
284342). Size S–M.

**R29 — τ-cluster morphism category + picture group.** [D-scout P11;
keep-with-corrections]
Object: W(A) from τ-perpendicular subcategories, cube-complex classifying
space, picture-group presentation. GATE: τ-tilting-finite only (infinite
category otherwise — loud refusal). Refs: Hanson–Igusa arXiv:1809.08989 (Comm.
Alg. 2021); Buan–Marsh; Igusa–Todorov signed-exceptional-sequences history
arXiv:2509.10910 (re-slotted as context for R28 too). Oracles: kA₂/kA₃ picture
groups; object count = # wide subcategories (ties R26); Nakayama K(π,1).
Size M.

## 6. Derived / silting

**R30 — Silting: verifier + single mutation + bounded exploration.** [D-scout
P4, MAJOR re-scope per critic]
Object: silting-object VERIFIER in K^b(proj) (Hom(T,T[n≠0]) = 0 + generation —
on the shipped derived stack), silting MUTATION (one approximation triangle),
and the co-t-structure dictionary. HONEST SCOPE: the silting quiver can be
infinite and mutation-transitivity is known only for special classes
(Aihara–Iyama's own statement) — no general BFS enumeration; 2-term slice =
existing τ-tilting (cross-check only, not a deliverable); bounded-radius
exploration with loud truncation; certified-complete only where finiteness is
proven (derived-discrete, local, τ-tilting-finite 2-term). CORRECTED refs:
Aihara–Iyama arXiv:1009.3370 (J. LMS 85 (2012)); **Oppermann**
arXiv:1504.02617 "Quivers for silting mutation"; **Jørgensen** arXiv:1603.09379
"Co-t-structures: the first decade". Oracles: k[x]/(x²) silting = shifts; kA₂
silting quiver (Aihara–Iyama); End(μT) quiver vs Oppermann's rule. Size L.

**R31 — Amiot–Keller cluster categories from (Q,W) (stretch).** [D-critic
omission; exploratory]
Object: the generalized cluster category C_{(Q,W)} for the Jacobi-finite
quivers-with-potential quiverlab already builds (surfaces → gentle Jacobians);
cluster-tilting objects, 2-CY verification. Heavier dg machinery; certified
scope to be established at spec time — candidate for a research-grade plan
after R30. Size XL.

## 7. Gentle / skew-gentle / persistence

**R32 — Skew-gentle algebras: recognizer, classification, τ-tilting.** [D-scout
P10; keep]
Object: skew-gentle triples (Q, I, Sp) recognizer; idempotent-split reduction to
the associated gentle algebra, re-gluing special strings; support τ-tilting via
the orbifold model; brick-finite ⇔ rep-finite certificate. Refs: He–Zhou–Zhu
arXiv:2004.11136; Chen arXiv:2212.06467; Amiot arXiv:2107.02646; Garcia–Lavoué
arXiv:2601.01744 (verified real, Jan 2026). Oracles: Sp = ∅ byte-reduces to the
gentle engine; the 2212.06467 example counts; geometric-vs-engine τ-tilting
cross-check. Size M.

**R33 — Persistence modules over A_n and commutative ladders (TDA bridge).**
[D-scout P9; keep]
Object: barcodes = interval decompositions of A_n/zigzag modules (existing
string/AR machinery); CL(n) = A_n □ A_2 for n ≤ 4 (rep-finite — "length < 5"
verified exact): AR-quiver-indexed generalized persistence diagrams; CL(≥5)
loud rep-infinite refusal. Refs: Escolar–Hiraoka arXiv:1404.7588 (DCG;
Hiraoka-group provenance honestly flagged); Igusa–Rock–Todorov
arXiv:1909.10499; Botnan–Crawley-Boevey. Oracles: A₅ filtration barcode;
Escolar–Hiraoka's CL(3) worked diagram. Size S–M.

**R34 — Homological rep-finite string-algebra test.** [A-scout P6, CORRECTED
statement]
Object: among rep-finite algebras, A is string ⇔ the middle term of every
**extension of indecomposables** (arbitrary Ext¹ classes — NOT just AR
sequences) has ≤ 2 indecomposable summands; finite check over all
Ext¹(X,Y)-middle-terms on the knitted category. Ref: Suárez-Álvarez
arXiv:2105.02948 (Alg. Rep. Theory 26 (2023) 1759–1772); the
Huisgen-Zimmermann–Smalø citation dropped (wrong content — critic). Oracle:
discriminating battery vs is_string across the zoo. Size M.

**R35 — Toupie algebras: recognizer + closed-form HH + sl_a structure.**
[A-scout P5; keep-with-corrections]
Object: trivial graph recognizer (unique source/sink, parallel paths); closed
HH forms as an oracle family; HH¹ ⊇ sl_a identified over ℂ/char 0 ONLY (scope
per critic). Refs: Artenstein–Lanzilotta–Solotar arXiv:1803.10310 (Alg. Rep.
Theory 2020 — cite the journal version); Artenstein thesis (Colibri UdelaR).
Size S–M.

## 8. Koszul beyond quadratic

**R36 — N-Koszul / K₂ / multi-Koszul certifiers.** [A-scout P7;
keep-with-corrections]
Object: generation degrees of Ext•(k,k) off the shipped minimal resolutions;
N-Koszul (single-degree relations + the 2-N alternation pattern — checkable),
K₂ (Cassidy–Shelton, generated in degrees 1,2 — needs an explicit certified
WINDOW, honest inconclusive beyond), multi-Koszul (Herscovich — stated for
connected graded; the f.d. kQ/I transfer needs care at spec time). Quadratic
case = existing Plan-27 (overlap named). Refs: Herscovich arXiv:1305.1678;
Herscovich JPAA 223 (2019) 1054–1072 (year corrected); Chouhy arXiv:1708.02933
(N-Koszul degeneration stability). Size M.

---

## 9. Priority tiers (proposal for the P51+ program)

**Tier α — flagship engine upgrades (double-corroborated / low-hanging):**
R1 (bracket beyond the window — the calculus completion), R4 (HH coefficients —
prerequisite of four other records), R2 (BV — backlog item now executable),
R24 (fractional CY — cheap flagship), R23 (φdim/ψdim on Plan 40).

**Tier β — the recognizer & certificate ladder (C-cluster spine):**
R15 → R16 → {R17, R18, R21} → R19 (+Bongartz) → R20; R14 (π₁) beside them;
R34, R35 as batteries.

**Tier γ — τ-tilting/lattice/derived expansion:** R25, R26, R27, R28, R29;
R30 (silting, re-scoped); R32 (skew-gentle); R33 (persistence/TDA — outreach
value).

**Tier δ — deformation & structure (research-grade):** R11 → R12 → R13;
R3 (Tate), R8 (skew group), R5+R6+R7 (extension/reduction machinery),
R9, R10, R36. R31 last (XL, exploratory).

Each plan minted from a record MUST: re-verify the record's references at spec
time (especially the two flagged attributions R7/2409.00945 and R2's verbatim
formulas), keep the honest-scope boundaries as loud refusals, and land its
oracles on the verification page per the Plan-22/32 standing rules.

---

*Provenance: 4 scout reports + 4 adversarial verification reports, 2026-08-06.
Scout/critic transcripts are session artifacts; this document is the durable
adjudicated record. Version 1.0.0 releases when this program completes
spec → plan → implement → verify (Marco, 2026-08-06).*
