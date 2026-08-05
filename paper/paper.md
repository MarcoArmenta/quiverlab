---
title: "QuiverLab: exact representation theory of quivers with relations in Python"
tags:
  - Python
  - representation theory
  - quivers with relations
  - finite-dimensional algebras
  - Hochschild (co)homology
  - Cup product
  - Gerstenhaber bracket
  - Cap product
  - Auslander–Reiten theory
authors:
  - name: Marco Armenta
    orcid: 0000-0003-3023-7634
    corresponding: true
    email: marco.armenta@usherbrooke.ca
    affiliation: 1
affiliations:
  - index: 1
    name: "Institut quantique, Université de Sherbrooke, Sherbrooke, Québec, Canada"
    ror: 00kybxq39
date: 4 August 2026
bibliography: paper.bib
---

# Summary

`quiverlab` is a pure-Python library for exact computation in the representation
theory of finite-dimensional associative algebras, presented as quivers with
relations, `A = kQ/I`. Given a quiver and a list of relation strings it certifies
that the algebra is finite-dimensional, builds an exact multiplication table, and
computes the working invariants of the field: modules with their Ext and Tor,
minimal projective and injective resolutions, the Auslander–Reiten translates,
Krull–Schmidt decomposition into indecomposables, Yoneda Ext-algebras with a
Koszulity verdict, Hochschild cohomology and homology with their Tamarkin-Tsygan
calculus (cup product, cap action, and bracket), cyclic homology, and Cartan and Coxeter data. Every number is exact. The library works over the rationals, over exact
subfields of the complex numbers (algebraic number fields $\mathbb{Q}(\alpha)$), and over every
finite field `GF(p^n)`, and it fails loudly on any floating-point input rather than
returning an approximation. It is built for algebraists who need not program: a
presentation reaches a certified table in three lines, the same computations run
in a browser with no code, and any computation can emit a human-readable
worked-steps document. Version 0.1.0 ships on PyPI, as a multi-arch container image, and as one-file desktop applications.

# Statement of need

Hochschild cohomology $HH^{\bullet}(A)$, with its Gerstenhaber algebra structure, is a
central invariant in representation theory and deformation theory, and it is hard to
compute by hand past the smallest examples. The same holds for the module-theoretic
and homological data that surround it. Researchers who study quivers with relations,
the standard presentation of finite-dimensional algebras, have had no installable
tool that computes these invariants exactly and without programming overhead.
`quiverlab` is built to be that tool, around two commitments.

The first is exactness. All arithmetic runs over an exact coefficient domain,
rationals, number fields, or finite fields, so there is no floating-point code path
in the algebra, and a static-analysis gate rejects float literals in the source. A
characteristic-`p` phenomenon, such as the collapse of the norm map on `k[x]/(x^p)`,
appears by changing a single argument rather than as rounding noise. The second is
certification per instance. Rather than trust that an algorithm is correct in
general, `quiverlab` checks each answer against explicit conditions and against
independent engines, and refuses loudly when it cannot certify a result.

Two long-term goals organize the project.

- **No code required.** Every computation is reachable without writing code: draw the
  quiver and relations on a browser canvas, specify modules entry by entry in a
  no-code panel, or export a configuration file for a cluster, and read the results as
  rendered mathematics or a printable worked-steps report. Python is a power-user option, never a
  prerequisite.
- **Any computation done in representation theory.** The aim is that whatever a
  representation theorist of finite-dimensional algebras computes in a paper is
  computable here, exactly and with oracle-tested results, or the documentation says
  honestly why not yet; the distance is tracked openly as a coverage program.

# State of the field

The closest existing system is **QPA** [@qpa], a mature package for the computer
algebra system GAP [@gap4]. QPA constructs `kQ/I` by admissible ideals and computes
minimal projective resolutions, module Ext, and much of the Auslander–Reiten
apparatus of the subject; `quiverlab` treats it as its primary external oracle and
reproduces its results wherever the two overlap. QPA also marks the boundary. It
ships no Hochschild cohomology, which must be assembled by hand through the
enveloping algebra, and no cup product or Gerstenhaber bracket; it has no native Tor,
no Koszulity verdict, and no facility for the deep-degree cup and cap products that
`quiverlab` computes past the bar-resolution window, and installing it requires GAP. **SageMath** [@sagemath] provides the free path
algebra but no quotient-by-relations object, and an unreduced bar complex usable only
at toy sizes. **Magma**, **Macaulay2/Singular**, and **QuiverTools** address adjacent
problems (Ext algebras, noncommutative Gröbner bases, moduli of representations) but
none computes finite-dimensional Hochschild theory with its operations. On PyPI there
is nothing of this kind, and we believe `quiverlab` to be the first system to
implement the Chouhy–Solotar resolution [@ChouhySolotar2015] in full; its exact-only kernel, four resolution engines, and no-code surfaces are designed to occupy that space rather than extend any of these systems.

# Software design

`quiverlab` is layered, and the layering is what keeps it exact. A coefficient-domain
protocol carries all arithmetic exactly: rationals by fraction-free elimination,
number fields through `sympy`, and finite fields through a fast integer kernel with a
pure-Python fallback. No engine has a floating-point branch, and the accelerated and
pure kernels are required to agree bit for bit. On top of the domain, a
quiver-with-relations front end runs a noncommutative Gröbner (Buchberger–Mora
overlap) completion with a degree bound and an admissibility certificate: it either
returns a certified finite-dimensional algebra with an irreducible-path basis, or
fails loudly with the offending cycle, never hanging and never guessing.

The homological core offers four interchangeable resolutions of the regular
bimodule: the normalized bar complex [@Hochschild1945], the minimal corner-typed
`A^e` resolution, Bardzell's resolution for monomial algebras [@Bardzell1997], and
the general Chouhy–Solotar resolution [@ChouhySolotar2015], which specializes exactly
to Bardzell in the monomial case. Hochschild dimensions are computed per instance and
gated by three independent checks: $d \circ d = 0$, an order condition, and degreewise
agreement with the bar and minimal engines inside the window where those are
buildable. The Gerstenhaber operations [@Gerstenhaber1963; @NegronWitherspoon2016;
@Volkov2019] are transported to bar cochains inside that window, and, because the cup
and cap products are also computed natively on the small Chouhy–Solotar model through
a comparison-lifted diagonal, past it, where the sign convention is arbitrated by the
Leibniz rule and by agreement with the transported operation rather than assumed.
Cyclic homology rests on Connes' mixed complex [@Connes1985]. Around this core sits a
module-theoretic surface: left and right modules on the same engine, minimal
projective and injective resolutions with Ext and Tor [@GSZ2001], the
Auslander–Reiten translates $\tau$ and $\tau^{-}$ built through the opposite algebra and the
`k`-duality [@ASS2006], Krull–Schmidt decomposition into indecomposables, Yoneda
Ext-algebras with a three-valued Koszulity verdict [@Priddy1970; @Froberg1999;
@PolishchukPositselski2005], Frobenius and symmetric recognition by a
nondegenerate-trace-form certifier [@SkowronskiYamagata2011], and the trivial
extension `T(A)` returned as a genuine quiver presentation [@Happel1988].

The differentiator is not the list of features but how they are trusted. `quiverlab`
is verified against two classes of oracle, documented on a living verification page. The first is theory and literature:
algebras whose invariants the literature or a theorem has already resolved are
constructed, and `quiverlab` must reproduce the published value exactly, over several
characteristics, each source cited at the precision the repository can verify. The
second is cross-engine and external agreement: the independent resolutions must agree
degreewise where they overlap; the accelerated and pure kernels must agree exactly;
and wherever QPA implements a feature, an optional `pip install quiverlab[qpa]`
backend drives GAP to recompute it and refuses to disagree silently. The two classes are
complementary: identical corruption of two internal engines is caught by a
literature pin or by QPA, a mistranscribed literature value by the live
cross-engine agreement. Of the library's 2,772 automated
tests, 1,581 are pinned by at least one such oracle, and the full suite runs on
every change. Where QPA cannot reach, over `CC` and
`GF(p^n)`, for cyclic homology, for the deep cup and cap, and for the Koszul verdict,
the verification page names the theory oracle that stands in, and where nothing yet
certifies a result the library refuses rather than guess.

# Research impact statement

`quiverlab` makes exact Hochschild and module computations routine where they were
manual, so that Gerstenhaber structure, homological finiteness, and Auslander–Reiten
data can be studied across families of algebras and across characteristics rather
than one example at a time. It reproduces published invariants — the Hochschild
dimensions of quantum complete intersections [@BGMS2005; @BerghErdmann2008] and the
vanishing behavior of hereditary algebras [@Happel1989] — and extends them along
parameter sweeps impractical by hand: the quantum-complete-intersection cohomology
$[2,2,1,0,\dots]$ is confirmed independent of the exponents `(a,b)` through `(5,5)` while
the homology $[a{+}b{-}1,\, a{+}b{-}2,\, \dots]$ grows with them, the `m`-Kronecker algebra is
verified to have $HH^1 = m^2{-}1$ [@Happel1989], and Bardzell's minimal resolution
[@Bardzell1997] reaches Hochschild degree 300 on the dimension-220 self-injective
Nakayama algebra $kZ_{20}/J^{11}$, far past the degree at which the bar complex is
buildable. The representative ongoing project is the author's *hanlab* program on Han's
conjecture [@Han2006], that vanishing of high-degree Hochschild homology forces
finite global dimension: `quiverlab` drives its search for counterexamples,
sweeping families of admissible presentations in batch and pushing each algebra's
Hochschild homology to deep degrees on HPC clusters through the library's
checkpoint-and-resume engines and its container/Slurm tier, flagging any candidate
whose homology vanishes while no finite resolution appears. The library is
nonetheless independent of the application that produced it.

Exactness and per-instance certification also change what a computation can be used
for. Because every reported number is checked, a disagreement with the literature is
informative rather than fatal. On Example 2.20 of Cibils–Redondo–Saorín
[@cibilsredondosaorin2004] the validated bar oracle robustly returns `dim HH^1 = 1`,
an explicit surviving oriented cycle where the published statement records `0`;
`quiverlab` pins the verified value over both orientations and several
characteristics and documents the discrepancy. Being exact, pip-installable, and
reproducible, the library lowers the barrier to computational experiments in the
representation theory of finite-dimensional algebras and provides a citable reference
implementation, of the Chouhy–Solotar resolution in particular, against which future
work can be checked.

# AI usage disclosure

`quiverlab` is developed with heavy use of AI, and its verification doctrine is a
direct consequence of that fact. The design and the mathematics are the author's; the
implementation is carried out with Anthropic's Claude Code, whose Claude models
(Opus- and Fable-class, 2026) work as multiple coordinated agents under the
author's direction, with research agents surveying the
literature and paired adversarial critics reviewing each unit of work. This is exactly why every reported invariant is pinned by a literature or theory oracle and, wherever QPA implements it, cross-checked independently: no value `quiverlab` ships rests on unverified generated code. The disagreement recorded above is an
instance of the same discipline: the oracles, not the expectation of the author or of
an assistant, decide what the library reports. The author has reviewed, edited, and
validated all AI-assisted output, and takes full responsibility for the content of
the library and of this paper.

# Acknowledgements

We thank the developers of QPA [@qpa] and GAP [@gap4], and the SymPy and NumPy
projects, on whose tools and ideas this work builds. No external funding supported
this work.

# References
