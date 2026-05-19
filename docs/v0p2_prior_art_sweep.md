# v0.2 prior-art second-wave sweep

> Second-wave prior-art search after the first-wave scoping +
> acquisition. Targets the five gaps identified after the first wave:
> Cuyt school, van Iseghem, Aptekarev school, Klein/Bothner/Lisovyy/
> Joshi 2024–2026, and a GitHub code search.
>
> Conducted 2026-05-19, ~30 web fetches, verified citations only.

## Verdict (TL;DR)

**The "first explicit vector / multi-component Padé–Taylor for ODEs
with pole-bridging applied to Painlevé-type systems" claim survives the
second-wave sweep, but with two important qualifications added to the
v0.1 / v0.2 framing.** No paper, code repository, or PhD thesis in any
of the five searched gaps performs the composition of (i) Hermite–Padé
with shared denominator, (ii) Jorba–Zou-style adaptive Taylor stepping,
and (iii) Painlevé / Garnier / Noumi–Yamada targeting. Each gap supplies
algebraic, theoretical, or algorithmic *ingredients* that v0.2 must
cite — none supplies the composition.

The two qualifications:

1. **Name collision and parallel scalar prior art:** Amore 2021
   (arXiv:2111.12198, Physica D 2022) calls its scalar van-der-Pol
   integrator the "Padé–Taylor method (PTM)". The construction (Taylor
   jet + diagonal Padé + residual-based step control) is algorithmically
   adjacent to v0.1's pipeline though scalar, per-component, and without
   FW / GGT robustification or Painlevé pole-field framing. **v0.1's
   docs should acknowledge Amore 2021 as parallel prior art for the
   scalar Padé–Taylor method**; v0.2 should do the same.

2. **Calgo 766 is the missing FORTRAN analogue of GGT 2013 for the
   shared-Q case:** Cabay–Jones–Labahn 1997 (ACM TOMS Calgo 766) is a
   weakly stable FORTRAN implementation of Padé–Hermite and simultaneous
   Padé approximant computation along a Padé table diagonal. It is not
   an ODE solver — but it *is* an oracle the v0.2 shared-Q port can
   mutation-prove against, exactly the way `padeapprox.m` was used in
   v0.1 for the scalar GGT 2013 port. v0.2's ADRs should cite Calgo 766
   alongside Beckermann–Labahn 1992/1994 and van Iseghem 1987.

## Gap 1 — Cuyt school (Antwerp; multivariate Padé)

**Searched:** Cuyt's Wikipedia page, the *Annie@60: A Life in
Approximation* Festschrift, the QD-algorithm-and-multivariate-Padé
Numer. Math. paper (Cuyt 1983), the *Padé Approximants for Operators*
1984 monograph references, the *Multivariate Padé Approximants Revisited*
BIT 1990 paper, the 2008 *Handbook of Continued Fractions for Special
Functions* (Cuyt–Petersen–Verdonk–Waadeland–Jones), and Cuyt's "Padé
approximants in operator theory for the solution of nonlinear
differential and integral equations" (Comput. Math. Appl. 1982) plus
its Cuyt–Van der Cruyssen 1983 system-of-equations companion.

**Found:** Cuyt 1982 is the closest match by *paper title*, but on
inspection (search-result abstract + bibliographic record + Wikipedia
+ Scholarpedia summaries) the *Padé in operator theory* programme is
abstract operator iteration (Halley-method generalisation of Newton,
applied to nonlinear operator equations including ODE IVPs / BVPs
abstractly). It does **not** carry out step-by-step Padé–Taylor
numerical ODE integration with movable-pole step control.

Cuyt's multivariate Padé programme is fully orthogonal to v0.2's
single-z / vector-y setting: multivariate Padé approximates a function
of *several* variables, while v0.2 needs vector Padé approximation in
*one* variable (the analytic-continuation parameter `z`) with multiple
output components.

No PhD student of Cuyt's appears in the searched material to have
written the composition either. The MPF (Multivariate Padé Fraction)
software is unrelated to ODE solving.

**Verdict:** The gap survives — Cuyt's body of work does not contain
"vector Padé–Taylor for ODEs with pole-bridging applied to Painlevé
hierarchies." The metadata stub for Cuyt 1982 is now in
`references/hermite_pade/Cuyt1982_pade_operators_NDE_ComputMathAppl8_metadata.md`.

## Gap 2 — van Iseghem (and the Brezinski / Lille line)

**Searched:** Van Iseghem's 1987 J. Comput. Appl. Math. paper "Vector
orthogonal relations. Vector QD-algorithm", her 2003 *Vector Stieltjes
Continued Fraction and Vector QD Algorithm* paper, her unpublished
1987 Lille PhD thesis (*Approximants de Padé Vectoriels*), the
Brezinski–Van Iseghem 1994 *Handbook of Numerical Analysis* chapter,
and three further Lille-school papers on vector Padé approximants.

**Found:** Van Iseghem 1987 is the *direct* canonical 1987 source for
the algebraic object v0.2 needs — vector Padé approximation in one
variable with one shared denominator across all `d` components. Her
denominator polynomials `P_{s,r}` satisfy vector-orthogonality
recurrences of order `d + 1`, and the vector QD-algorithm is an
explicit alternative backend to the GGT2013-style SVD route. The
Brezinski–Van Iseghem 1994 chapter is the survey-level account
covering scalar, vector, and matrix Padé approximants in the same
algebraic style.

But neither paper applies the construction to numerical ODE solving:
both stay on the pure approximation-theory side (given a power series,
produce the rational approximant). There is no step control, no pole
bridging, no Painlevé / Garnier framing.

**Verdict:** The gap survives. Van Iseghem 1987 supplies the algebraic
foundation; v0.2 must cite it as such, and the design ADR for shared-Q
should follow van Iseghem's recurrence-order-`d+1` recipe (alternative
to or composable with the block-SVD route). Metadata stub:
`references/hermite_pade/VanIseghem1987_vector_QD_JCAM19_metadata.md`.

## Gap 3 — Aptekarev / Steklov school (most load-bearing)

**Searched:** Aptekarev–Kuijlaars 2011 *Russian Math. Surveys* survey
"Hermite–Padé approximations and multiple orthogonal polynomial
ensembles", Aptekarev–Stahl 1992 *Progress in Approximation Theory*
chapter "Asymptotics of Hermite–Padé Polynomials", Aptekarev–Lysov
2010 *Sb. Math.* "Systems of Markov functions generated by graphs",
Stahl's 1986 *Constr. Approx.* paper "Orthogonal polynomials with
complex-valued weight function", Nuttall's conjecture-resolution
literature, Fidalgo–López Lagomasino–Medina 2013 (arXiv:1310.7010),
López Lagomasino–Medina 2014 (arXiv:1406.3737), López Lagomasino
2019 (arXiv:1910.08548), Ikonomov–Kovacheva–Suetin 2016
(arXiv:1603.03314), Ikonomov–Suetin 2026 (arXiv:2605.14760).

**Found:** The convergence-theorem chain v0.2 needs is well-established
in the Aptekarev / Stahl line:

- **Aptekarev–Stahl 1992** (book chapter, paywalled): vector
  potential-theoretic equilibrium-problem framework for Hermite–Padé
  asymptotics. For Markov / Nikishin / Angelesco systems, the
  appropriately normalised Hermite–Padé polynomial zeros distribute
  according to the equilibrium measure, and the type II approximants
  converge uniformly on compact subsets of the complement of the
  contour system. Their zeros identify the poles of the meromorphic
  functions being approximated. This is the convergence theorem the
  v0.2 shared-Q ADR needs to cite.

- **Aptekarev–Kuijlaars 2011** (Russian Math. Surveys, paywalled):
  the canonical modern survey, places Hermite–Padé in the random
  matrix / multiple orthogonal polynomial ensemble framework. Useful
  for ADR context; not the load-bearing theorem source.

- **López Lagomasino 2019 (arXiv:1910.08548)** — fetched. Pedagogical
  survey that re-proves the type II uniform-convergence statement for
  Nikishin systems and cites both Stahl (refs [16, 17, 44, 45] in
  the survey, including the *General Orthogonal Polynomials* CUP 1992
  book) and the Aptekarev–Sorokin equilibrium-problem line. This is the
  citable open-access entry point.

- **Fidalgo–López Lagomasino–Medina 2013 (arXiv:1310.7010)** and
  **López Lagomasino–Medina 2014 (arXiv:1406.3737)** — both fetched.
  These treat the *rational-modification* case (Nikishin system
  perturbed by rational functions, i.e. finitely many added poles)
  which is closer to v0.2's actual setting: a meromorphic ODE solution
  approximated by Hermite–Padé.

- **Ikonomov–Suetin 2026 (arXiv:2605.14760)** — fetched. Most recent
  paper in the line; compares convergence rates of Padé vs rational
  Hermite–Padé for multivalued analytic functions using scalar Green-
  logarithmic potential problems.

None of these papers applies Hermite–Padé to numerical ODE solving
either; the application side is uniformly approximation theory.

**Verdict:** The gap survives. The Aptekarev / Stahl line gives the
v0.2 ADR the convergence-theorem citations it needs (`AptekarevStahl1992`
for the historical formulation, `LopezLagomasino2019` as the
open-access pedagogical entry, `FidalgoLagomasinoMedina2013` and
`LagomasinoMedina2014` for the rational-perturbation case). The
composition with Painlevé / Garnier / Noumi–Yamada ODEs is not in any of
these papers. Metadata stub for Aptekarev–Stahl 1992 at
`references/hermite_pade/AptekarevStahl1992_asymptotics_HP_polynomials_metadata.md`.

## Gap 4 — Klein / Bothner / Lisovyy / Joshi 2024–2026

**Searched:** arXiv author searches for Christian Klein (Dijon),
Thomas Bothner (Bristol), Oleg Lisovyy (Tours), and Nalini Joshi
(Sydney), 2024-01 through 2026-05.

**Found:**

- **Klein 2024–2026:** Klein–Stoilov 2024 (arXiv:2401.04461) on
  *Multi-domain spectral approach to rational-order fractional
  derivatives*; Klein–Stoilov 2024 (arXiv:2401.10141) on *Optimally
  truncated WKB approximation for the 1D stationary Schrödinger
  equation*; Klein–Sparber 2025 (*Ann. Henri Poincaré*) on defocusing
  fractional NLS. **None** on Painlevé equations or multi-component
  Painlevé systems. Klein's Painlevé-numerics line (Grava–Klein 2011,
  Kapaev–Klein–Grava 2015, Klein–Stoilov 2018) appears dormant.

- **Bothner 2024–2026:** Bothner–Wei 2025 (arXiv:2511.18118) and
  Bothner–Jaconelli 2025 (arXiv:2511.23362) on Painlevé asymptotics
  via Riemann–Hilbert; Bothner–Shepherd 2024 (arXiv:2411.18550) on
  random matrix universality. **None** release code; all use the
  established RH-asymptotics methodology.

- **Lisovyy 2024–2026:** Barhoumi–Lisovyy–Miller–Prokhorov 2024
  (SIGMA 20, 019) on Painlevé-III monodromy maps; Iorgov–Iwaki–Lisovyy
  –Zhuravlov 2025 (arXiv:2505.16803) "Many-faced Painlevé I" on
  irregular conformal blocks. **No code releases.** Both papers are
  theoretical (tau function, conformal-block side), not numerical.

- **Joshi 2024–2026:** Joshi 2024 (arXiv:2405.10541) on Segre
  surfaces and Painlevé geometry; Joshi 2024 (arXiv:2410.14327) on
  QRT maps; Joshi 2024 (arXiv:2408.07963) on crystal limit of
  q-difference PVI; Joshi–Roffelsen 2025 (arXiv:2508.18578) and
  Joshi 2026 (arXiv:2602.09298) on arithmetic / finite-characteristic
  dynamics. **None** are vector / multi-component Painlevé numerical
  methods; none release code.

- **Surprise — Adler–Sokolov 2025 (arXiv:2512.18828):** "Vector
  systems of Painlevé type" (December 2025; published 2026, *J. Geom.
  Phys.*) — V.E. Adler & V.V. Sokolov, **not** from any of the four
  groups asked about. The paper applies group reduction to vector
  generalisations of NLS, mKdV, and KdV and obtains ODE systems with
  isomonodromic Lax representations generalising P_I, P_II, P_34, and
  P_IV. It is **pure classification / Lax-pair construction, with no
  numerical method and no code release.** This is exactly the kind of
  recent paper v0.2 docs should cite as the algebraic-motivation
  upstream of multi-component Painlevé: their derivation produces
  systems v0.2's vector solver could integrate. Fetched into
  `references/rh_numerical/AdlerSokolov2025_vector_painleve_type_2512.18828.pdf`.

**Verdict:** The gap survives. No 2024–2026 paper from these groups
releases code or publishes a vector Painlevé numerical method that
scoops v0.2. Adler–Sokolov 2025 is a notable new *target* equation
source for v0.2 — exactly the kind of vector-Painlevé classification
that motivates having a vector solver — and **needs to be added to
PRD §v0.2 north star as a 2025 reinforcement of the "vector Painlevé
hierarchies have classification but no numerical software" thesis.**

## Gap 5 — GitHub code search

**Queries run (via `gh search code` and `gh search repos`):**

```
gh search code "noumi yamada"                  → 0 numerical hits
gh search code "painleve hierarchy"            → 0 production-quality hits;
                                                 a couple of arxiv-metadata mirrors
                                                 and one Mazzocco-Mo PII hierarchy
                                                 paper-supplement reference
gh search code "garnier system" numerical      → 1 todo-list mention in
                                                 openxm-org/OpenXM (computer algebra
                                                 wishlist, no implementation), plus
                                                 arxiv metadata mirrors
gh search code "hermite pade" "differential
   equation"                                   → 1 hit: ore_algebra/guessing.py
                                                 (Mark Kauers; D-finite series
                                                 guessing, scalar Hermite–Padé
                                                 for guessing ODEs from terms —
                                                 NOT a Painlevé / vector ODE solver)
gh search code "tritronquee"                   → multiple hits, all single-component:
                                                 cool-japan/scirs (RK45 PI/PII only),
                                                 bilman/AsymptoticsOfRWIO-Paper-Code
                                                 (P2 paper-supplement Jupyter
                                                 notebook, scalar RH method),
                                                 raeez/chiral-bar-cobar (kappa-
                                                 Painlevé toy code), grantbaker/
                                                 painleve-analysis (Mathematica
                                                 PI–PVI utilities, 2017, no
                                                 hierarchy), dmivilensky/Canonical-
                                                 barriers (PIII D7 only)
gh search code "simultaneous pade" "ode"       → ACM Calgo 766 (FORTRAN,
                                                 vector_pade.f), and
                                                 Beliavsky/Alan-Miller-Fortran
                                                 VEC_PADE.F90 mirror —
                                                 approximation theory, not ODE
gh search code "shared denominator" "pade"     → 0 relevant hits
gh search repos painleve                       → 30 repos; relevant ones:
                                                 JuliaHolomorphic/Painleve.jl
                                                 (Wolfram, 66 KB, last touched
                                                 2021-01-04, 0 stars, 0 forks,
                                                 1 watcher, no releases —
                                                 stub-only, never finished);
                                                 grantbaker/painleve-analysis
                                                 (Mathematica, 2017, PI–PVI
                                                 scalar utilities, never updated);
                                                 rpm4real/painleve-integration
                                                 (paper code, PI only, 2016);
                                                 willyhereman/PainleveTest-
                                                 Mathematica (Painlevé property
                                                 *test*, not solver);
                                                 ericrenone/IMFL-Isomonodromic-
                                                 Frobenius-Learning (2026, ML
                                                 toy demo, not a numerical
                                                 method);
                                                 cool-japan/scirs/scirs2-special
                                                 /src/painleve (RK45 PI/PII only,
                                                 confirmed scalar via raw source)
```

**Verdict — `painleve` topic + code search confirms PRD's claim:**
*GitHub has no production-quality, registered, vector- or
multi-component-Painlevé numerical package.* The most directly
named candidate, `JuliaHolomorphic/Painleve.jl`, is a 2019-created
stub with no releases, no commits since 2021-01-04, 0 stars, mostly
Wolfram code; the most plausible competing scalar numerical Painlevé
code is `cool-japan/scirs`'s `scirs2-special/src/painleve` module
(RK45 PI/PII only, no Padé, no pole-field).

The strongest *adjacent* prior art is `ore_algebra` (Mark Kauers,
SageMath ecosystem): uses scalar Hermite–Padé approximation to *guess*
D-finite differential equations from initial term sequences — a
different problem (recovering the ODE from data) than solving the ODE
forward through movable singularities. v0.2 should cite it in the
"adjacent uses of Hermite–Padé in symbolic / numeric computation" ADR
prose.

## New papers acquired

| Path | Cluster | Gap served | Acquisition |
|---|---|---|---|
| `references/hermite_pade/LopezLagomasino2019_intro_multiple_orthogonal_hermite_pade_1910.08548.pdf` | hermite_pade | Gap 3 | arXiv open PDF, fetched |
| `references/hermite_pade/FidalgoLagomasino_Medina_2013_HP_meromorphic_1310.7010.pdf` | hermite_pade | Gap 3 | arXiv open PDF, fetched |
| `references/hermite_pade/LagomasinoMedina2014_HP_typeI_meromorphic_1406.3737.pdf` | hermite_pade | Gap 3 | arXiv open PDF, fetched |
| `references/hermite_pade/IkonomovSuetin2026_convergence_rational_HP_2605.14760.pdf` | hermite_pade | Gap 3 | arXiv open PDF, fetched |
| `references/rh_numerical/AdlerSokolov2025_vector_painleve_type_2512.18828.pdf` | rh_numerical | Gap 4 (surprise: not from the asked groups but a 2025 vector Painlevé classification paper) | arXiv open PDF, fetched |
| `references/Amore2021_pade_taylor_van_der_pol_PhysicaD_arXiv2111.12198.pdf` | top-level (methodology-adjacent to v0.1) | Cross-cutting: name-collision parallel scalar Padé–Taylor work | arXiv open PDF, fetched |

## New metadata stubs

| Path | Cluster | Gap served | Reason for stub |
|---|---|---|---|
| `references/hermite_pade/CabayJonesLabahn1997_algorithm_766_VECTOR_PADE_ACMTOMS23_metadata.md` | hermite_pade | Gap 2 / Pillar A | ACM TOMS paywalled; FORTRAN source freely available via Calgo archive |
| `references/hermite_pade/Cuyt1982_pade_operators_NDE_ComputMathAppl8_metadata.md` | hermite_pade | Gap 1 | Elsevier paywalled; Wikipedia / Scholarpedia confirm scope |
| `references/hermite_pade/VanIseghem1987_vector_QD_JCAM19_metadata.md` | hermite_pade | Gap 2 | Elsevier paywalled; abstract verified via search |
| `references/hermite_pade/AptekarevStahl1992_asymptotics_HP_polynomials_metadata.md` | hermite_pade | Gap 3 (load-bearing) | Springer book chapter paywalled; convergence-theorem source |
| `references/Amore2021_pade_taylor_van_der_pol_PhysicaD_arXiv2111.12198_metadata.md` | top-level | Cross-cutting | Companion to the PDF — name-collision context |

## Open follow-up questions

1. **Aptekarev–Stahl 1992 chapter acquisition.** The chapter is paywalled
   and not on Stahl's archived personal page (Stahl passed away in
   2013) nor immediately visible on Aptekarev's Keldysh / Steklov
   institutional page. If institutional library access is available,
   fetch the PDF as a v0.2 ADR oracle. Until then, the López Lagomasino
   2019 survey covers the citable convergence statements.

2. **Calgo 766 PDF + source.** The journal paper is paywalled; the
   FORTRAN source is openly redistributable via the ACM Calgo archive
   (https://calgo.acm.org/). If the v0.2 shared-Q port wants a direct
   FORTRAN-vs-Julia mutation-proof oracle (analogous to v0.1's
   `padeapprox.m` round-trip), the source can be wrapped via Julia's
   ccall to LAPACK-style FORTRAN. This is non-blocking but worth a
   beads issue for a v0.2-stage-0 numerical-oracle audit.

3. **`ore_algebra` Hermite–Padé use.** `ore_algebra`'s
   `guessing.py:hermite_pade` is used to recover D-finite differential
   equations from term sequences. v0.2 prose should mention this as
   the *symbolic-computation* counterpart to v0.2's *numerical* use
   (solving forward through movable poles rather than guessing the
   ODE), and acknowledge it as adjacent prior art on the
   "Hermite–Padé applied to ODEs" axis.

4. **PRD §v0.2 north star, "Survey sources" list** (PRD lines 394–414)
   should be updated to add:
   - Amore 2021 (arXiv:2111.12198) — name-collision scalar parallel
     work, no vector / Painlevé scope.
   - Adler–Sokolov 2025 (arXiv:2512.18828) — *Vector systems of
     Painlevé type*, 2025 classification, reinforces the "no
     numerical software exists yet" thesis with new target equations.
   - Van Iseghem 1987 — algebraic foundation for shared-Q vector Padé.
   - Cabay–Jones–Labahn 1997 (Calgo 766) — FORTRAN oracle for
     shared-Q computation.
   - López Lagomasino 2019 — open-access survey citation for the
     Hermite–Padé convergence theorem v0.2 depends on.

5. **Calligraphic v0.1 retroactive update.** If we accept Amore 2021
   as parallel scalar prior art, v0.1's RESEARCH.md and README should
   acknowledge it in a "parallel work" paragraph. This is a v0.1 doc
   maintenance task, not v0.2-blocking.
