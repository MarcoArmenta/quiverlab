---
title: "quiverlab: exact Hochschild and representation theory for quivers with relations in Python"
tags:
  - Python
  - representation theory
  - quivers with relations
  - finite-dimensional algebras
  - Hochschild cohomology
  - Gerstenhaber bracket
  - Auslander–Reiten theory
  - exact computation
authors:
  - name: Marco Armenta
    orcid: 0000-0000-0000-0000        # TODO(Marco): replace with the real ORCID before submission
    corresponding: true
    email: drmarcoarmenta@gmail.com
    affiliation: 1
affiliations:
  - index: 1
    name: "Affiliation TODO(Marco)"   # TODO(Marco): set the real affiliation before submission
date: 26 July 2026
bibliography: paper.bib
---

# Summary

`quiverlab` is a pure-Python library for exact computation with finite-dimensional
associative algebras presented as quivers with relations, `A = kQ/I`. Given a quiver
and a list of relation strings it certifies that the algebra is finite-dimensional,
builds an exact multiplication table, and computes a wide range of homological and
representation-theoretic invariants: Hochschild cohomology and homology with their
Gerstenhaber operations (cup product, cap action, and bracket), cyclic homology,
module Ext and Tor, minimal projective and injective resolutions, the
Auslander–Reiten translates, Krull–Schmidt decomposition into indecomposables,
Yoneda Ext-algebras with a Koszulity verdict, and Cartan, Coxeter, and spectral
data. Every number is exact. The library works over the rationals, over exact
subfields of the complex numbers (algebraic number fields `Q(α)`), and over every
finite field `GF(p^n)`, and it fails loudly on any floating-point input rather than
returning an approximation. It is built for research algebraists and does not require
them to program: a presentation reaches a certified table in three lines, the same
computations run in a browser with no code at all, and any computation can emit a
human-readable worked-steps document.

# Statement of need

Hochschild cohomology `HH^\bullet(A)`, with its Gerstenhaber algebra structure, is a
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

Two long-term goals organize the project, and every release is measured against them.

- **No code required.** Every computation is reachable without writing code: draw the
  quiver and relations on a browser canvas, specify modules entry by entry in a
  no-code panel, or export a configuration file for a cluster, and read the results as
  rendered mathematics or a PDF report. Python is a power-user option, never a
  prerequisite.
- **Any computation done in representation theory.** The aim is that whatever a
  representation theorist of finite-dimensional algebras computes in a paper is
  computable here, exactly and with oracle-tested results, or the documentation says
  honestly why not yet. The distance between this goal and the current surface is
  tracked openly as a coverage program.

# State of the field

The closest existing system is **QPA** [@qpa], a mature package for the computer
algebra system GAP [@gap4]. QPA constructs `kQ/I` by admissible ideals and computes
minimal projective resolutions, module Ext, and much of the Auslander–Reiten
apparatus of the subject; `quiverlab` treats it as its primary external oracle and
reproduces its results wherever the two overlap. QPA also marks the boundary. It
ships no Hochschild cohomology, which must be assembled by hand through the
enveloping algebra, and no cup product or Gerstenhaber bracket; it has no native Tor,
no Koszulity verdict, and no facility for the deep-degree cup and cap products that
`quiverlab` computes past the bar-resolution window. Installation requires GAP, a
barrier for users who do not program. **SageMath** [@sagemath] provides the free path
algebra but no quotient-by-relations object, and an unreduced bar complex usable only
at toy sizes. **Magma**, **Macaulay2/Singular**, and **QuiverTools** address adjacent
problems (Ext algebras, noncommutative Gröbner bases, moduli of representations) but
none computes finite-dimensional Hochschild theory with its operations. On PyPI there
is nothing of this kind, and we believe `quiverlab` to be the first system to
implement the Chouhy–Solotar resolution [@ChouhySolotar2015] in full. `quiverlab` is
designed to occupy that space rather than to extend any of these systems: its
exact-only kernel, its four resolution engines, and its no-code surfaces differ in
architecture from each of them.

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
gated by three independent checks: `d∘d = 0`, an order condition, and degreewise
agreement with the bar and minimal engines inside the window where those are
buildable. The Gerstenhaber operations [@Gerstenhaber1963; @NegronWitherspoon2016;
@Volkov2019] are transported to bar cochains inside that window, and, because the cup
and cap products are also computed natively on the small Chouhy–Solotar model through
a comparison-lifted diagonal, past it, where the sign convention is arbitrated by the
Leibniz rule and by agreement with the transported operation rather than assumed.
Cyclic homology rests on Connes' mixed complex [@Connes1985]. Around this core sits a
module-theoretic surface: left and right modules on the same engine, minimal
projective and injective resolutions with Ext and Tor [@GSZ2001], the
Auslander–Reiten translates `τ` and `τ⁻` built through the opposite algebra and the
`k`-duality [@ASS2006], Krull–Schmidt decomposition into indecomposables, Yoneda
Ext-algebras with a three-valued Koszulity verdict [@Priddy1970; @Froberg1999;
@PolishchukPositselski2005], Frobenius and symmetric recognition by a
nondegenerate-trace-form certifier [@SkowronskiYamagata2011], and the trivial
extension `T(A)` returned as a genuine quiver presentation [@Happel1988].

The differentiator is not the list of features but how they are trusted. `quiverlab`
is verified against two classes of oracle, documented on a verification page that
every development plan extends as it ships. The first is theory and literature:
algebras whose invariants the literature or a theorem has already resolved are
constructed, and `quiverlab` must reproduce the published value exactly, over several
characteristics, each source cited at the precision the repository can verify. The
second is cross-engine and external agreement: the independent resolutions must agree
degreewise where they overlap; the accelerated and pure kernels must agree exactly;
and wherever QPA implements a feature, an optional `pip install quiverlab[qpa]`
backend drives GAP to recompute it and refuses to disagree silently. The two classes
are complementary. A regression that corrupted two of the library's own engines
identically would still be caught by a literature pin or by QPA, and a mistranscribed
literature value would be caught by the live cross-engine agreement. Over two thousand
such oracle-pinned tests run on every change. Where QPA cannot reach, over `CC` and
`GF(p^n)`, for cyclic homology, for the deep cup and cap, and for the Koszul verdict,
the verification page names the theory oracle that stands in, and where nothing yet
certifies a result the library refuses rather than guess.

# Research impact statement

`quiverlab` makes exact Hochschild and module computations routine where they were
manual, so that Gerstenhaber structure, homological finiteness, and Auslander–Reiten
data can be studied across families of algebras and across characteristics rather
than one example at a time. It reproduces published invariants, among them the
Hochschild dimensions of quantum complete intersections [@BGMS2005; @BerghErdmann2008]
and the vanishing behavior of hereditary algebras [@Happel1989], and then extends them
along parameter sweeps that would be impractical by hand. The tool grew out of, and
now serves, work on Han's conjecture [@Han2006] and the finiteness questions around
it, but it is independent of the application that produced it.

Exactness and per-instance certification also change what a computation can be used
for. Because every reported number is checked, a disagreement with the literature is
informative rather than fatal. On Example 2.20 of Cibils–Redondo–Saorín
[@cibilsredondosaorin2004] the validated bar oracle robustly returns `dim HH^1 = 1`,
an explicit surviving oriented cycle, where the published statement records `0`;
`quiverlab` pins the verified value, over both orientations and several
characteristics, and documents the discrepancy rather than freezing the printed
number against a live engine disagreement. Being exact, pip-installable, and
reproducible, the library lowers the barrier to computational experiments in the
representation theory of finite-dimensional algebras and provides a citable reference
implementation, of the Chouhy–Solotar resolution in particular, against which future
work can be checked.

# AI usage disclosure

`quiverlab` is developed with heavy use of AI, and its verification doctrine is a
direct consequence of that fact. The design and the mathematics are the author's; the
implementation is carried out by Claude (Anthropic) working as multiple coordinated
agents under the author's direction, with research agents surveying the literature,
paired adversarial critics reviewing each unit of work, and the author verifying the
mathematics line by line. This is exactly why every invariant the library reports is
pinned by a literature or theory oracle and, wherever the GAP package QPA implements
it, cross-checked against an independent computation, so that no value shipped by
`quiverlab` rests on unverified generated code. The disagreement recorded above is an
instance of the same discipline: the oracles, not the expectation of the author or of
an assistant, decide what the library reports. The author takes full responsibility
for the content of the library and of this paper.

# Acknowledgements

We thank the developers of QPA [@qpa] and GAP [@gap4], and the SymPy and NumPy
projects, on whose tools and ideas this work builds.
<!-- TODO(Marco): add financial support / funding acknowledgement before submission -->

# References
