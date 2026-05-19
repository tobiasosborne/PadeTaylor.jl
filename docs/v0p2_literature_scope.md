# v0.2 vector Painlevé — literature scope

> Output of the scoping pass that precedes the Stage 0+ deep dive. Each
> row below is a recommendation for the fan-out wave; the fan-out's job
> is to fetch the PDF (or HTML where PDF is unavailable) into
> `references/` following the existing repo convention. Every citation
> has been verified via arXiv abstract pages, Springer/Cambridge DOI
> pages, or author preprint pages between 2026-05-19 02:00 and 04:00.
> Gaps are flagged `UNVERIFIED`.
>
> Two **PRD corrections** surfaced during scoping (PRD lines 402, 403):
>
> 1. `arXiv:0704.2869` (PRD line 402) is by **Yusuke Sasano**, *Studies
>    on the Garnier system in two variables*, not Mazzocco. The
>    Mazzocco Garnier paper actually exists and is `arXiv:math/0106208`
>    (*The geometry of the classical solutions of the Garnier systems*,
>    2001, published IMRN 2002). Both should be acquired; PRD's
>    intended source is Mazzocco's.
> 2. `arXiv:2107.11680` (PRD line 403) is by **Bobrova–Sokolov**, not
>    Adler–Sokolov. The Adler–Sokolov paper at `arXiv:2012.05639`
>    *is* correctly attributed.
>
> The Stage 0+ writer should fix both attributions in PRD as part of
> the deep-dive commit.

## Naming / placement convention (proposed)

The Heun acquisition (`references/heun/` + `references/markdown/heun/`)
sets the cluster precedent. v0.2 has four mathematically distinct
threads, so propose four cluster directories plus the same top-level
placement for the two or three most load-bearing papers in the v0.2
analogue of the FW2011/GGT2013 spot:

- `references/hermite_pade/` — Pillar A. Beckermann–Labahn, Mano–Tsuda,
  Gonnet–Pachón–Trefethen, Nakatsukasa–Sète–Trefethen AAA, multivariate
  Cuyt entries. Top-level pin: **Mano–Tsuda 2017** (the only paper
  that connects Hermite–Padé to Painlevé / Garnier directly).
- `references/noumi_yamada/` — Pillar B. Noumi–Yamada 1998 (the
  founding paper), Matsuda rational-solutions papers (A_4^{(1)},
  A_5^{(1)}), Clarkson–Gómez-Ullate–Grandati–Milson cyclic Maya
  diagrams, Sakai 2001, Aratyn–Gomes–Zimerman, Gómez-Ullate–Grandati–Milson
  2020 survey. Top-level pin: **Noumi–Yamada 1998**.
- `references/painleve_hierarchy/` — Pillar C. Cosgrove 2000-6 (P_I^2
  classification), Kapaev–Klein–Grava 2015, Grava–Klein 2011, Joshi–Mazzocco
  2002, Claeys 2010/2011, Kudryashov, Clarkson–Mansfield Yablonski–
  Vorobiev. Top-level pin: **Kapaev–Klein–Grava 2015** (only published
  pole-field for any hierarchy member).
- `references/garnier/` — Pillar D. Mazzocco 2001, Sasano 2007, Mano–Tsuda
  hypergeometric-Garnier overlap, Calligaris–Mazzocco 2018, Iwasaki–
  Kimura–Shimomura–Yoshida (book; cluster contains a `book_metadata.md`
  noting it is paywalled).
- `references/calogero_moser/` — Pillar F. Airault–McKean–Moser 1977,
  Krichever 1980, Wilson 1993 bispectral.
- `references/rh_numerical/` — Pillar E (adjacent / competitive
  methodology). Olver 2011 Numerische Mathematik framework, Olver 2012
  FoCM Painlevé II, Trogdon–Olver SIAM OT146 (book; metadata only),
  Wechslberger–Bornemann 2012 automatic deformation, Klein–Stoilov
  2018 SIGMA, Claeys–Olver 2011 higher Tracy–Widom, Bornemann 2010,
  Gavrylenko–Lisovyy 2018, Korotkin/Kitaev–Korotkin 1998,
  Bertola–Cafasso 2011, Adler–Sokolov 2020/Bobrova–Sokolov 2021
  (context only).
- **Top-level (alongside FW2011 etc.):** Mano–Tsuda 2017,
  Novokshenov 2014 *Constr. Approx.* (a directly methodology-adjacent
  Padé-for-Painlevé paper that the v0.1 RESEARCH.md missed),
  Kapaev–Klein–Grava 2015. These three are the "v0.2 ground-truth
  pin" analogues of v0.1's FW2011 + GGT2013.

Each cluster directory mirrors `references/heun/`: PDFs flat in the
top of the cluster, with a per-cluster `README.md` describing what's
inside and the per-paper relevance one-liner. Marker conversion outputs
go to `references/markdown/<cluster>/<paper-slug>/`.

## Already in the repo (do not re-fetch)

From `/home/tobias/Projects/PadeTaylor.jl/references/` (top-level v0.1):

- `BerrutTrefethen2004_barycentric_SIAMReview.pdf`
- `FFW2017_painleve_riemann_surfaces_preprint.pdf`
- `FW2011_painleve_methodology_JCP230.pdf`
- `FW2014_second_PII_exploration_FoCM14.pdf`
- `FW2015_imaginary_PII_PhysicaD309.pdf`
- `GGT2013_robust_pade_via_SVD_SIREV55.pdf`
- `JorbaZou2005_taylor_IVP_package_ExpMath14.pdf`
- `Mezzarobba2019_truncation_bounds_dfinite_arXiv1904.pdf`
- `ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.pdf`
- `TrefethenSMIM_2000_book.pdf`
- `WeidemanReddy2000_DMSUITE_ACMTOMS26.pdf`

From `references/heun/`:

- `Fiziev2009_Heun_BH_review_0902.1277.pdf`
- `Hortacsu2018_Heun_applications_physics_1101.0471.pdf`
- `Motygin2015_Heun_evaluation_1506.03848.pdf`
- `Motygin2018_confluent_Heun_numerical_1804.01007.pdf`
- `Motygin2020_regularization_Heun_2010.09053.pdf`
- `Suzuki1998_Teukolsky_Heun_grqc9805064.pdf`
- `references/heun/teukolsky/BertiCardosoStarinets2009_QNM_review_0905.2975.pdf`
- `references/heun/teukolsky/Fiziev2009_Schwarzschild_Heun_0906.5108.pdf`
- `references/heun/dlmf-31/31.{2,3,4,12}.html`

Total existing: 19 PDFs + 4 DLMF HTMLs. Do not refetch any of these.

## Tier 1 — must-have ground truth (12 papers)

| # | Canonical ID | Authors | Year | Title (short) | Venue | Pillar | Relevance (one sentence) | Proposed path | Retrieval URL | Acquisition notes |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `arXiv:1502.06695` / Math. Z. 285 | Mano, Tsuda | 2017 | Hermite–Padé approximation, isomonodromic deformation and hypergeometric integral | Mathematische Zeitschrift 285 (2017) 397–431 | A,D | The single most load-bearing v0.2 paper — directly connects Hermite–Padé approximation to Schlesinger transformations, the sixth Painlevé equation, and the Garnier system, and uses block-Toeplitz determinants throughout, which is exactly the linear-algebra substrate the shared-denominator port needs. | `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf` | https://arxiv.org/abs/1502.06695 | Open arXiv PDF — fetch latest v3 (2016-04-30). |
| 2 | USyd Research Report 2000-6 | Cosgrove | 1999/2000 | Higher-order Painlevé equations in the polynomial class II. Bureau symbol P1 | USyd preprint; later Stud. Appl. Math. 116 (2006) 321–413 | C | Canonical classification of polynomial-class higher-order Painlevé equations (18 equations F-I through F-XVIII), supplying the algebraic forms PadeTaylor needs to build hierarchy-member constructors. | `references/painleve_hierarchy/Cosgrove2000_USyd_2000-6_higher_painleve_polynomial_classII.pdf` | https://ssh.maths.usyd.edu.au/u/ResearchReports/Nonlinear/Cos/2000-6.html | **Page provides DVI/PostScript only, no PDF**; fan-out must `ps2pdf` after fetch. Also cross-check the 2006 *Studies in Applied Mathematics* journal version (DOI 10.1111/j.1467-9590.2006.00346.x) via Wiley. |
| 3 | `arXiv:1306.6161` / Constr. Approx. 41 | Kapaev, Klein, Grava | 2013/2015 | On the tritronquée solutions of P_I^2 | Constr. Approx. 41 (2015) 425–466 | C,E | The only existing published pole-field figure for any Painlevé-hierarchy member; PadeTaylor.jl's PI^(2) acceptance test in PRD line 374–375 cross-validates against this paper specifically. | `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf` | https://arxiv.org/abs/1306.6161 | Open arXiv PDF. |
| 4 | `arXiv:math/9808003` / Funkcial. Ekvac. 41 | Noumi, Yamada | 1998 | Higher order Painlevé equations of type A_l^{(1)} | Funkcial. Ekvac. 41 (1998) 483–503 | B | The founding paper for the A_n^{(1)} Noumi–Yamada hierarchy — defines the systems v0.2 takes as primary target, gives the affine-Weyl symmetry that powers the W(A_n^{(1)}) self-test in PRD line 318. | `references/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41.pdf` | https://arxiv.org/abs/math/9808003 | Open arXiv PDF (math.QA). |
| 5 | Noumi 2004 / AMS Translations 223 | Noumi (M.) | 2004 | Painlevé Equations Through Symmetry | AMS, Translations of Mathematical Monographs vol. 223 | B | The canonical book PRD line 315 names; chapter 4 covers A_n^{(1)} Hamiltonian normal forms, Bäcklund generators, and rational solutions, all needed for the closed-form-family verifier per PRD line 320. | `references/noumi_yamada/Noumi2004_AMS_MMONO223_metadata.md` (book; PDF not free) | https://bookstore.ams.org/mmono-223 | **Book — AMS Translations vol. 223 (paperback ISBN 978-0-8218-3221-9).** No open author copy known; the Japanese 2000 original is also closed. Acquisition route: institutional library scan, AMS member purchase, or used-book buy. Mark v0.2 as *partially* blocked without it (the 1998 paper at #4 covers the affine-Weyl generators; the book adds the rational-solution closed forms). |
| 6 | `arXiv:1611.00958` (note typo correction) / Comm. Math. Phys. 363 | Gavrylenko, Lisovyy | 2016/2018 | Fredholm determinant and Nekrasov sum representations of isomonodromic tau functions | Comm. Math. Phys. 363 (2018) 1–58 | E | The canonical tau-function CFT route — out of scope as a v0.2 method (per PRD line 274–276 and 422), but required reading for the v0.2 documentation to honestly explain why it is *not* the route PadeTaylor takes. Includes the Painlevé VI + Garnier specialisation that scopes the comparison cleanly. | `references/rh_numerical/GavrylenkoLisovyy2018_isomonodromic_tau_CMP363.pdf` | https://arxiv.org/abs/1608.00958 | Open arXiv PDF. Correct ID is `1608.00958`, not `1712.08546` (which is the related Cafasso–Gavrylenko–Lisovyy *Tau functions as Widom constants* — also worth fetching; see Tier 2). |
| 7 | `arXiv:1612.00337` / SISC 40 | Nakatsukasa, Sète, Trefethen | 2018 | The AAA algorithm for rational approximation | SIAM J. Sci. Comput. 40 (2018) A1494–A1522 | A | The Trefethen-school successor to GGT 2013 for rational approximation; AAA generalises to matrix-valued and to "set-AAA" — provides a credible benchmark numerical-Hermite–Padé route to compare the v0.2 shared-Q implementation against. | `references/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf` | https://arxiv.org/abs/1612.00337 | Open arXiv PDF. |
| 8 | `arXiv:1101.2602` / Comm. Math. Phys. 313 | Claeys, Grava | 2011 | The KdV hierarchy: universality and a Painlevé transcendent | Comm. Math. Phys. 313 (2012) 557–595 | C,E | One of the two PRD-named pole-field papers for the PI hierarchy; supplies the connection to the small-dispersion KdV limit that PadeTaylor.jl's Calogero–Moser smoke test (PRD line 361) cross-validates against. | `references/painleve_hierarchy/ClaeysGrava2011_KdV_painleve_CMP313.pdf` | https://arxiv.org/abs/1101.2602 | Open arXiv PDF. |
| 9 | `arXiv:0704.2869` | Sasano | 2007 | Studies on the Garnier system in two variables | arXiv preprint (later refs in J. Math. Soc. Japan) | D | **Correct identification of PRD line 402's arXiv ID** — Sasano's Hamiltonian/symmetry/holomorphy study of the two-variable Garnier system; useful for the v0.2 Garnier-tertiary target's explicit Hamiltonian forms. | `references/garnier/Sasano2007_garnier_two_variables_0704.2869.pdf` | https://arxiv.org/abs/0704.2869 | Open arXiv PDF. |
| 10 | `arXiv:math/0106208` / IMRN 2002 | Mazzocco | 2001/2002 | The geometry of the classical solutions of the Garnier systems | International Mathematics Research Notices 2002 (2002) | D | **This is what PRD line 402 *meant* to cite as the Mazzocco Garnier paper.** Provides the Riemann–Hilbert geometric framework for Garnier classical solutions — the conceptual anchor for the Garnier-tertiary acceptance criteria. | `references/garnier/Mazzocco2002_geometry_garnier_classical_IMRN_math0106208.pdf` | https://arxiv.org/abs/math/0106208 | Open arXiv PDF. |
| 11 | `arXiv:0708.0074` / J. Math. Phys. 53 | Matsuda | 2007/2012 | Rational solutions of the Noumi and Yamada system of type A_4^{(1)} | J. Math. Phys. 53 (2012) 023504 | B | The closed-form rational solutions for A_4^{(1)} — the exact oracle PadeTaylor.jl needs to build a `noumi_yamada_rational` analogue of the v0.1 `pii_rational` verifier (PRD line 352–353). | `references/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53.pdf` | https://arxiv.org/abs/0708.0074 | Open arXiv PDF. |
| 12 | Constr. Approx. 39 | Novokshenov | 2014 | Distributions of Poles to Painlevé Transcendents via Padé Approximations | Constr. Approx. 39 (2014) 85–99 | E,A | **New find missed by v0.1.** Uses a Fair–Luke Padé algorithm to compute pole distributions of PI/PII/PIV transcendents in the complex plane — a *direct methodology-adjacent* paper to PadeTaylor's own algorithm, including for PIV which v0.1 already implements. Mutation-proof opportunity: reproduce Novokshenov's pole pictures with PadeTaylor and confirm agreement, then lift the same method to PJ hierarchies. | `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39.pdf` (top-level — methodology-adjacent) | https://link.springer.com/article/10.1007/s00365-013-9190-6 | Springer paywalled; check arXiv mirror or ResearchGate fallback. If no open copy, mark as institutional-fetch. |

## Tier 2 — strong supporting context (20 papers)

| # | Canonical ID | Authors | Year | Title (short) | Venue | Pillar | Relevance (one sentence) | Proposed path | Retrieval URL | Acquisition notes |
|---|---|---|---|---|---|---|---|---|---|---|
| 13 | Numer. Algorithms 3 | Beckermann, Labahn | 1992 | A uniform approach for Hermite Padé and simultaneous Padé approximants and their matrix-type generalisations | Numerical Algorithms 3 (1992) 45–54 | A | The textbook reduction of Hermite–Padé and simultaneous Padé to a single matrix-Padé framework — provides the theoretical basis for the shared-denominator algorithm v0.2 needs to port from GGT2013. | `references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf` | https://link.springer.com/article/10.1007/BF02141914 | Springer paywalled; **author preprint at https://cs.uwaterloo.ca/~glabahn/Papers/tenerife.pdf** — try this first. |
| 14 | SIMAX 15 | Beckermann, Labahn | 1994 | A uniform approach for the fast computation of matrix-type Padé approximants | SIAM J. Matrix Anal. Appl. 15 (1994) 804–823 | A | Companion to #13 — gives the fast (super-linear) algorithm; relevant once the v0.2 shared-Q port is GREEN and we want to outperform a naive block-SVD. | `references/hermite_pade/BeckermannLabahn1994_fast_matrix_pade_SIMAX15.pdf` | https://epubs.siam.org/doi/abs/10.1137/S0895479892230031 | SIAM paywalled; check Labahn UW page first. |
| 15 | Electron. Trans. Numer. Anal. 38 | Gonnet, Pachón, Trefethen | 2011 | Robust rational interpolation and least-squares | ETNA 38 (2011) 146–167 | A | The Trefethen-school precursor to GGT 2013's SVD-Padé robustification — covers SVD-based Froissart-doublet detection in the rational-interpolation setting, which we'll need to adapt component-by-component. | `references/hermite_pade/GonnetPachonTrefethen2011_robust_rational_interp_ETNA38.pdf` | https://people.maths.ox.ac.uk/trefethen/publication/PDF/2011_141.pdf | Open author preprint at Oxford. |
| 16 | Cambridge UP / Encyc. Math. Appl. 59 | Baker, Graves-Morris | 1996 (2nd ed.) | Padé Approximants | Cambridge University Press, Encyclopedia of Mathematics and its Applications vol. 59 | A | The canonical Padé textbook; PRD line 308 names "chapter 5" — actually chapter on simultaneous Padé tables; cite-by-reference rather than fetch. | `references/hermite_pade/BakerGravesMorris1996_pade_approximants_metadata.md` | https://www.cambridge.org/9780521450072 | **Book — Cambridge UP.** Not openly available; mark as institutional library or used-book. The chapter on simultaneous Padé tables is canonical but not v0.2-blocking (Mano–Tsuda 2017 + Beckermann–Labahn cover the same ground). |
| 17 | `arXiv:0708.2960` / arXiv preprint | Matsuda | 2007 | Rational Solutions of the Noumi and Yamada System of type A_5^{(1)} | arXiv preprint (75 pp.) | B | Closed-form rational solutions for A_5^{(1)} — companion to #11; A_3^{(1)} reduces to PV, A_5^{(1)} is the next nontrivial step beyond the A_4^{(1)} oracle. | `references/noumi_yamada/Matsuda2007_rational_A5_NoumiYamada_0708.2960.pdf` | https://arxiv.org/abs/0708.2960 | Open arXiv PDF. |
| 18 | `arXiv:1811.09274` / Stud. Appl. Math. 144 | Clarkson, Gómez-Ullate, Grandati, Milson | 2018/2020 | Cyclic Maya diagrams and rational solutions of higher order Painlevé systems | Stud. Appl. Math. 144 (2020) 357–397 | B | Unified Maya-diagram construction of rational solutions for the A_{2n} Painlevé / Noumi–Yamada systems; generalised Hermite/Okamoto/Umemura polynomial closed forms — the principled extension of the A_4^{(1)} oracle from #11 to A_{2n}^{(1)} for arbitrary n. | `references/noumi_yamada/Clarkson_etal_2020_cyclic_maya_higher_painleve_StudApplMath144.pdf` | https://arxiv.org/abs/1811.09274 | Open arXiv PDF. |
| 19 | Comm. Math. Phys. 220 | Sakai | 2001 | Rational Surfaces Associated with Affine Root Systems and Geometry of the Painlevé Equations | Comm. Math. Phys. 220 (2001) 165–229 | B,C | The Sakai classification anchor for all (discrete and continuous) Painlevé families via affine-Weyl symmetry; the conceptual foundation under Noumi–Yamada plus the Garnier extension. | `references/noumi_yamada/Sakai2001_rational_surfaces_painleve_CMP220.pdf` | https://link.springer.com/article/10.1007/s002200100446 | Springer paywalled. Author preprint may be on Sakai's page at Kobe Univ; fan-out should check. |
| 20 | `arXiv:math/0212117` / Nonlinearity 16 | Joshi, Mazzocco | 2002/2003 | Existence and uniqueness of tri-tronquée solutions of the second Painlevé hierarchy | Nonlinearity 16 (2003) 427–439 | C | Existence/uniqueness foundation for tritronquée solutions in the PII hierarchy — provides the analogue of the PI tritronquée acceptance test for any future PII^{(n)} extension of v0.2. | `references/painleve_hierarchy/JoshiMazzocco2003_tritronquee_PII_hierarchy_math0212117.pdf` | https://arxiv.org/abs/math/0212117 | Open arXiv PDF. |
| 21 | `arXiv:1107.0214` / Physica D | Claeys | 2011 | Pole-free solutions of the first Painlevé hierarchy and non-generic critical behavior for the KdV equation | arXiv preprint; Physica D | C | Establishes existence of real pole-free solutions across the PI hierarchy — the analogue of the Hastings–McLeod theorem for arbitrary hierarchy members, which gates "is there a smooth answer here?" sanity checks. | `references/painleve_hierarchy/Claeys2011_pole_free_PI_hierarchy_1107.0214.pdf` | https://arxiv.org/abs/1107.0214 | Open arXiv PDF. |
| 22 | `arXiv:1001.2213` / J. Phys. A 43 | Claeys | 2010 | Asymptotics for a special solution to the second member of the Painlevé I hierarchy | J. Phys. A 43 (2010) 434012 | C | Asymptotic structure of the canonical PI^{(2)} smooth solution; provides closed-form behaviour at infinity that PadeTaylor.jl can use to bootstrap initial data for the PI^{(2)} acceptance figure. | `references/painleve_hierarchy/Claeys2010_PI2_asymptotics_JPhysA43.pdf` | https://arxiv.org/abs/1001.2213 | Open arXiv PDF. |
| 23 | `arXiv:1807.04442` / SIGMA 14 | Klein, Stoilov | 2018 | Numerical Approach to Painlevé Transcendents on Unbounded Domains | SIGMA 14 (2018) 068 | E | The Klein-school multidomain spectral method for Painlevé equations on unbounded domains — the most direct competing numerical method we should benchmark v0.2 figures against on PI tritronquée. | `references/rh_numerical/KleinStoilov2018_multidomain_spectral_painleve_SIGMA14.pdf` | https://arxiv.org/abs/1807.04442 | Open arXiv PDF (SIGMA is also open-access). |
| 24 | `arXiv:1111.3527` | Claeys, Olver | 2011 | Numerical study of higher order analogues of the Tracy–Widom distribution | arXiv preprint | E | Sheehan Olver's only published higher-Tracy–Widom paper using RH numerics — context for "no released higher-hierarchy software" and a competing pole-field-adjacent calculation. | `references/rh_numerical/ClaeysOlver2011_higher_tracy_widom_1111.3527.pdf` | https://arxiv.org/abs/1111.3527 | Open arXiv PDF. |
| 25 | `arXiv:0804.2543` / Math. Comp. 79 | Bornemann | 2008/2010 | On the numerical evaluation of Fredholm determinants | Math. Comp. 79 (2010) 871–915 | E | The canonical Nyström-method Fredholm determinant paper — context: the *other* methodology by which Tracy–Widom-type quantities are computed without ever drawing a pole field. | `references/rh_numerical/Bornemann2010_fredholm_determinants_MathComp79.pdf` | https://arxiv.org/abs/0804.2543 | Open arXiv PDF. |
| 26 | `arXiv:math-ph/9810007` / IMRN 1998 | Kitaev, Korotkin | 1998 | On solutions of the Schlesinger equations in terms of Θ-functions | IMRN 1998 no. 17, 877–905 | E | The original theta-function representation of Schlesinger / PVI solutions; context for "what closed-form representations exist that PadeTaylor.jl is not trying to compete with." | `references/rh_numerical/KitaevKorotkin1998_schlesinger_theta_IMRN_mathph9810007.pdf` | https://arxiv.org/abs/math-ph/9810007 | Open arXiv PDF. |
| 27 | `arXiv:1101.3997` / Comm. Math. Phys. 309 | Bertola, Cafasso | 2011/2012 | Fredholm determinants and pole-free solutions to the noncommutative Painlevé II equation | Comm. Math. Phys. 309 (2012) 793–833 | E | Matrix/noncommutative Painlevé II — context for v0.3 (PRD line 271–272 lists matrix Painlevé as v0.3 plausible), not for v0.2; include as a single reference paragraph in Pillar E. | `references/rh_numerical/BertolaCafasso2012_noncommutative_PII_CMP309.pdf` | https://arxiv.org/abs/1101.3997 | Open arXiv PDF. |
| 28 | `arXiv:2012.05639` / TMP 207 | Adler, Sokolov | 2020/2021 | On matrix Painlevé II equations | Theor. Math. Phys. 207 (2021) 560–571 | E | Matrix PII via Painlevé–Kovalevskaya test — context for v0.3, also referenced by PRD line 403. | `references/rh_numerical/AdlerSokolov2021_matrix_PII_TMP207.pdf` | https://arxiv.org/abs/2012.05639 | Open arXiv PDF. |
| 29 | `arXiv:2107.11680` | Bobrova, Sokolov | 2021 | On matrix Painlevé-4 equations. Part 1: Painlevé–Kovalevskaya test | arXiv preprint (math.CA) | E | Matrix PIV — companion of #28; PRD line 403 misattributes to Adler–Sokolov. Same v0.3 disposition. | `references/rh_numerical/BobrovaSokolov2021_matrix_PIV_test_2107.11680.pdf` | https://arxiv.org/abs/2107.11680 | Open arXiv PDF. |
| 30 | `arXiv:1206.2446` / Constr. Approx. 39 | Wechslberger, Bornemann | 2012/2014 | Automatic Deformation of Riemann–Hilbert Problems with Applications to the Painlevé II Transcendents | Constr. Approx. 39 (2014) 137–172 | E | The RH-numerics counterpart to GGT2013 — automatic contour deformation as a graph-shortest-path. Survey-supporting evidence that the RH route, while general, ships no general-purpose code. | `references/rh_numerical/WechslbergerBornemann2014_automatic_RH_painleve_ConstrApprox39.pdf` | https://arxiv.org/abs/1206.2446 | Open arXiv PDF. |
| 31 | CPAM 30 / Wiley | Airault, McKean, Moser | 1977 | Rational and elliptic solutions of the Korteweg–de Vries equation and a related many-body problem | Comm. Pure Appl. Math. 30 (1977) 95–148 | F | The founding Calogero–Moser ↔ KdV pole-soliton paper — gives the closed-form solutions that PadeTaylor's CM smoke test in PRD line 361 cross-validates against. | `references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30.pdf` | https://onlinelibrary.wiley.com/doi/10.1002/cpa.3160300106 | Wiley paywalled; institutional access required. There may be an NYU/Courant preprint mirror — fan-out should check Moser's archive. |
| 32 | Funct. Anal. Appl. 14 / Springer | Krichever | 1980 | Elliptic solutions of the Kadomtsev–Petviashvili equation and integrable systems of particles | Functional Analysis and Its Applications 14 (1980) 282–290 | F | The Krichever pole-soliton correspondence in the elliptic case — required for the KdV-vs-CM cross-validation reference values PRD line 361 cites. | `references/calogero_moser/Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf` | https://link.springer.com/article/10.1007/BF01078304 | Springer paywalled. |

## Tier 3 — nice-to-have / surveys (15 entries)

| # | Canonical ID | Authors | Year | Title (short) | Venue | Pillar | Relevance (one sentence) | Proposed path | Retrieval URL | Acquisition notes |
|---|---|---|---|---|---|---|---|---|---|---|
| 33 | `arXiv:1608.00958` / CMP 365 | Cafasso, Gavrylenko, Lisovyy | 2017/2019 | Tau functions as Widom constants | Comm. Math. Phys. 365 (2019) 741 | E | Companion to #6 — Widom-determinant tau functions; useful background for the "out-of-scope" prose. | `references/rh_numerical/CafassoGavrylenkoLisovyy2019_widom_tau_CMP365.pdf` | https://arxiv.org/abs/1712.08546 | Open arXiv PDF. (PRD line 405 cited this ID as Gavrylenko–Lisovyy; correct citation has Cafasso as first author.) |
| 34 | `arXiv:1010.5725` | Aratyn, Gomes, Zimerman | 2010 | Integrable Origins of Higher Order Painlevé Equations | arXiv preprint (nlin.SI) | B,C | A self-similarity-limit construction of A_n^{(1)}-invariant higher-order Painlevé equations — alternative derivation of the same systems, useful for cross-checking the Hamiltonian forms in #5. | `references/noumi_yamada/AratynGomesZimerman2010_integrable_origins_higher_painleve_1010.5725.pdf` | https://arxiv.org/abs/1010.5725 | Open arXiv PDF. |
| 35 | `arXiv:2009.11668` | Gómez-Ullate, Grandati, Milson | 2020 | Rational solutions of Painlevé systems | Chapter in Nonlinear Systems and Their Remarkable Mathematical Structures vol. 2 | B | Modern survey on rational solutions of Painlevé / Noumi–Yamada systems via dressing chains and Maya diagrams; entry point to the A_n^{(1)} rational-solutions literature. | `references/noumi_yamada/GomezUllateGrandatiMilson2020_rational_painleve_survey_2009.11668.pdf` | https://arxiv.org/abs/2009.11668 | Open arXiv PDF. |
| 36 | Phys. Lett. A | Kudryashov | 1997 | The first and second Painlevé equations of higher order and some relations between them | Phys. Lett. A 224 (1997) 353–360 | C | The Kudryashov hierarchy classification — supplementary to Cosgrove 2000-6, with explicit Bäcklund relations between hierarchy members. | `references/painleve_hierarchy/Kudryashov1997_first_second_painleve_higher_order_PLA224.metadata.md` | https://www.sciencedirect.com/science/article/pii/S0375960197008505 | Elsevier paywalled; institutional access. Probably not v0.2-critical given #2. |
| 37 | `arXiv:1401.1408` / IMRN 2015 | Bertola, Bothner | 2014/2015 | Zeros of large degree Vorob'ev–Yablonski polynomials via a Hankel determinant identity | IMRN 2015 no. 19, 9330 | C | Asymptotic locations of PII rational-solution poles via Hankel determinants — independent oracle for v0.2 PII-rational verification (orthogonal-polynomial path). | `references/painleve_hierarchy/BertolaBothner2015_vorobiev_yablonski_zeros_IMRN.pdf` | https://arxiv.org/abs/1401.1408 | Open arXiv PDF. |
| 38 | `arXiv:1504.00440` / Constr. Approx. 45 | Balogh, Bertola, Bothner | 2015/2017 | Hankel Determinant Approach to Generalized Vorob'ev–Yablonski Polynomials and their Roots | Constr. Approx. 45 (2017) 1 | C | Generalised-Vorob'ev–Yablonski case (PII hierarchy) — same independent-oracle role as #37 but for the hierarchy. | `references/painleve_hierarchy/BaloghBertolaBothner2017_generalized_vorobiev_yablonski_ConstrApprox45.pdf` | https://arxiv.org/abs/1504.00440 | Open arXiv PDF. |
| 39 | `arXiv:1707.05222` / SIGMA 14 | Masoero, Roffelsen | 2017/2018 | Poles of Painlevé IV Rationals and their Distribution | SIGMA 14 (2018) 002 | B | Pole-distribution of generalised Hermite / Okamoto polynomials — direct reference figure for v0.2 to reproduce on PIV-from-A_2^{(1)} as part of the A_2^{(1)} self-validation. | `references/noumi_yamada/MasoeroRoffelsen2018_PIV_rationals_poles_SIGMA14.pdf` | https://arxiv.org/abs/1707.05222 | Open arXiv PDF (SIGMA open-access). |
| 40 | `arXiv:q-alg/9708018` / Nagoya Math. J. 153 | Noumi, Yamada | 1997/1999 | Symmetries in the fourth Painlevé equation and Okamoto polynomials | Nagoya Math. J. 153 (1999) 53–86 | B | The original Schur-function / 3-reduced KP representation of PIV rational solutions — the A_2^{(1)} prototype the v0.2 verifier should reproduce. | `references/noumi_yamada/NoumiYamada1999_PIV_symmetries_okamoto_NagoyaMathJ153.pdf` | https://arxiv.org/abs/q-alg/9708018 | Open arXiv PDF. |
| 41 | `arXiv:1705.03295` / J. Integrable Syst. 3 | Calligaris, Mazzocco | 2017/2018 | Finite orbits of the pure braid group on the monodromy of the 2-variable Garnier system | J. Integrable Syst. 3 (2018) | D | Algebraic-solution classification for 2-variable Garnier — supplies the "do algebraic-solution oracles exist?" answer for the Garnier-tertiary acceptance criteria. | `references/garnier/CalligarisMazzocco2018_finite_orbits_2variable_garnier_JIntegrSyst3.pdf` | https://arxiv.org/abs/1705.03295 | Open arXiv PDF. |
| 42 | `arXiv:1808.09190` / Compositio 156 | Diarra, Loray | 2018/2020 | Classification of algebraic solutions of irregular Garnier systems | Compositio Mathematica 156 (2020) 881–907 | D | Companion to #41 for the irregular Garnier case — useful for the v0.2 documentation to honestly bound which algebraic-solution oracles are available. | `references/garnier/DiarraLoray2020_irregular_garnier_algebraic_Compositio156.pdf` | https://arxiv.org/abs/1808.09190 | Open arXiv PDF. |
| 43 | Vieweg 1991 / Springer | Iwasaki, Kimura, Shimomura, Yoshida | 1991 | From Gauss to Painlevé: A Modern Theory of Special Functions | Vieweg, Aspects of Mathematics vol. 16 | D | The canonical Garnier-systems textbook (Chapters 4–5); foundational for the Garnier-tertiary work but not v0.2-blocking given the Mazzocco / Sasano coverage. | `references/garnier/IwasakiKimuraShimomuraYoshida1991_from_gauss_to_painleve_metadata.md` | https://link.springer.com/book/10.1007/978-3-322-90163-7 | **Book — Vieweg/Springer.** Paywalled; institutional library route. PRD line scoping flag: include metadata-only stub. |
| 44 | `nlin/0306020` | Kimura, Mano | 2003 | Irregular isomonodromic deformations for Garnier systems and Okamoto's canonical transformations | arXiv preprint | D | Irregular-singularity extension of Garnier — bridges into the Painlevé-hierarchy story; supplementary. | `references/garnier/KimuraMano2003_irregular_garnier_nlin0306020.pdf` | https://arxiv.org/abs/nlin/0306020 | Open arXiv PDF. |
| 45 | TaylorIntegration.jl docs + Zenodo records | Pérez-Hernández, Benet | 2019– | TaylorIntegration.jl Julia package | Zenodo DOI 10.5281/zenodo.2562353 (v0.4.1) and successors | G | The Julia-side reference for `jetcoeffs!` and `Vector{Taylor1{T}}` ergonomics — PRD line 311 explicitly names this as the substrate. **No standalone arXiv paper exists.** Defer to code-side reading; archive the Zenodo metadata. | `references/taylor_integration/PerezHernandezBenet_TaylorIntegration_metadata.md` | https://zenodo.org/records/2562353 | No PDF. Cite via Zenodo DOI in v0.2 docs. |
| 46 | `arXiv:1806.08650` | Gavrylenko, Iorgov, Lisovyy | 2018 | On solutions of the Fuji–Suzuki–Tsuda system | arXiv preprint | E | Higher-rank generalisation of the Painlevé VI tau function approach — context for v0.3, not v0.2. | `references/rh_numerical/GavrylenkoIorgovLisovyy2018_FST_system_1806.08650.pdf` | https://arxiv.org/abs/1806.08650 | Open arXiv PDF. |
| 47 | SIAM OT146 (book) | Trogdon, Olver | 2016 | Riemann–Hilbert Problems, Their Numerical Solution, and the Computation of Nonlinear Special Functions | SIAM Other Titles in Applied Mathematics 146 | E | The canonical RH-numerical book; covers RiemannHilbert.jl substrate. **Already named in PRD line 411–412 as a survey source.** | `references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md` | https://epubs.siam.org/doi/book/10.1137/1.9781611974201 | **Book — SIAM.** Paywalled; institutional/SIAM-member access. A sample chapter is openly at https://www.math.uci.edu/~ttrogdon/publications/TrogdonOlver-Sample.pdf — fan-out should fetch the sample plus the front matter for context. |

## Books and other non-PDF resources

| Item | Authors | Acquisition route | v0.2-blocking? |
|---|---|---|---|
| Noumi 2004, *Painlevé Equations Through Symmetry*, AMS MMONO 223 | Noumi (M.) | Library scan, AMS member purchase, used-book | **Partially blocking** — covers rational solutions and Bäcklund generators in book-only form. Mitigation: Noumi–Yamada 1998 (#4) + Matsuda 2007/2012 (#11, #17) + Clarkson et al. 2018 (#18) cover the algorithmically needed content. |
| Baker–Graves-Morris 1996, *Padé Approximants*, Cambridge UP | Baker, Graves-Morris | Library scan, used-book | Non-blocking — Mano–Tsuda 2017 + Beckermann–Labahn cover the simultaneous-Padé table material. |
| Iwasaki–Kimura–Shimomura–Yoshida 1991, *From Gauss to Painlevé*, Vieweg | Iwasaki, Kimura, Shimomura, Yoshida | Library scan, used-book | Non-blocking for v0.2 (Garnier is tertiary). Useful for v0.3. |
| Trogdon–Olver 2016, *Riemann–Hilbert Problems...*, SIAM OT146 | Trogdon, Olver | SIAM member access; sample chapter open | Non-blocking — sample chapter + arXiv:1210.2199 cover what we need for the "competing methodology" narrative in v0.2 ADRs. |
| TaylorIntegration.jl Julia package | Pérez-Hernández, Benet | Code + Zenodo DOI; no paper | Non-blocking — cite as software. |
| Cosgrove USyd 2000-6 preprint | Cosgrove | DVI / PostScript on USyd page → `ps2pdf` | Non-blocking once `ps2pdf` step is done. Cluster placement: `references/painleve_hierarchy/`. |

## Per-pillar narrative

### Pillar A — Hermite–Padé / simultaneous Padé via SVD

The biggest discovery in this scoping pass is **Mano–Tsuda 2017**
(`arXiv:1502.06695`, Math. Z. 285): the paper is a near-perfect match
to PadeTaylor v0.2's needs — it does Hermite–Padé approximation via
block-Toeplitz determinants, connects directly to Schlesinger
transformations and to the Painlevé VI / Garnier systems, and stays
on the *numerical-algebraic* side rather than the deep RH side. It is
the right intellectual heir of GGT 2013 for the multi-numerator case.

Below that, the Beckermann–Labahn pair (#13, #14) gives the abstract
uniform-Padé framework that subsumes simultaneous and Hermite–Padé as
special cases of a single matrix-Padé problem. Their 1992 paper
already covers the conceptual reduction we need; their 1994 SIMAX
paper gives the fast (Newton-style) algorithm we can adapt once the
naive block-SVD port is GREEN. Their preprints are openly available
at the Waterloo author page, so acquisition is easy.

The Trefethen-school numerical-rational-approximation line is
Gonnet–Pachón–Trefethen 2011 (#15, the SVD-Padé robustification
groundwork that GGT2013 built on) and Nakatsukasa–Sète–Trefethen 2018
(#7, AAA — a different, point-wise barycentric algorithm). AAA is not
a direct Hermite–Padé method, but its matrix-valued and set-AAA
extensions are the only Trefethen-school work on multi-output
rational approximation; we should benchmark v0.2 against AAA on a
shared Noumi–Yamada problem to confirm we're at least as accurate.

Baker & Graves-Morris (#16) is named in PRD line 308 but is a
textbook; everything algorithmic that we need is in the Mano–Tsuda
and Beckermann–Labahn papers, so the book itself is non-blocking. If
I had to pick one paper as Pillar A's load-bearing ground truth, it
would be **Mano–Tsuda 2017** without hesitation.

### Pillar B — Noumi–Yamada A_n^{(1)}

The 1998 Noumi–Yamada paper (#4, `arXiv:math/9808003`) is openly
available on arXiv and contains the full A_n^{(1)} system definition
and the affine-Weyl symmetry. That alone unblocks the v0.2 algorithm.
What it does *not* contain is the rational-solution closed forms
needed for the closed-form-family oracle; those live in three places:
- Noumi 2004 book §4 (#5, paywalled),
- Matsuda 2007/2012 papers on A_4^{(1)} (#11) and A_5^{(1)} (#17, openly on arXiv),
- Clarkson–Gómez-Ullate–Grandati–Milson 2018 (#18, openly on arXiv),
  which unifies the construction for arbitrary A_{2n}^{(1)} via cyclic
  Maya diagrams.

If only one of the rational-solutions papers is fetched, fetch **#18
(Clarkson et al. 2018)** — it covers the most ground and is openly
available, so the Noumi 2004 book is reduced to a nice-to-have. The
Masoero–Roffelsen 2018 paper (#39) supplies independent pole-distribution
oracles for A_2^{(1)} (= PIV), which provides the v0.1 ↔ v0.2 cross-check
PRD line 348 requires.

Sakai 2001 (#19) is the conceptual umbrella under which Noumi–Yamada
sits, useful for the v0.2 ADR explaining *why* A_n^{(1)} is the natural
parameterisation. Aratyn–Gomes–Zimerman 2010 (#34) is an independent
construction useful as an oracle for Hamiltonian forms.

### Pillar C — Painlevé hierarchies

Cosgrove 2000-6 (#2) is the canonical classification of fourth-order
polynomial-class Painlevé equations (18 equations) and supplies the
algebraic shape v0.2 needs to construct hierarchy-member problems.
**Acquisition wart:** the USyd page hosts DVI/PostScript only, no PDF
— the fan-out batch dedicated to author-preprint pages must include a
`ps2pdf` step. The 2006 *Studies in Applied Mathematics* version
(DOI 10.1111/j.1467-9590.2006.00346.x) is a fallback through Wiley.

Kapaev–Klein–Grava 2015 (#3, `arXiv:1306.6161`) is the load-bearing
reference for the v0.2 PI^{(2)} acceptance figure — the only published
pole-field figure for any hierarchy member, against which v0.2 must
cross-check. Claeys–Grava 2011 (#8) and Claeys 2010 (#22), 2011 (#21)
fill in the analytic existence/uniqueness story for hierarchy
solutions, particularly the KdV-small-dispersion connection.

Joshi–Mazzocco 2003 (#20) is the PII-hierarchy tritronquée existence
result — useful for a v0.2 extension to PII^{(n)} if we choose to go
there. Kudryashov 1997 (#36) and the Yablonski-Vorobiev hierarchy
papers (#37, #38) are nice-to-have orthogonal oracles.

### Pillar D — Garnier G_n

The PRD-named source `arXiv:0704.2869` is **not by Mazzocco**; it is
Sasano's paper (#9), which is itself useful. The Mazzocco paper PRD
meant to cite is `arXiv:math/0106208` (#10), which gives the
Riemann–Hilbert framework for Garnier classical solutions and is
arguably the load-bearing paper for the tertiary target.

Mano–Tsuda 2017 (#1) also covers Garnier explicitly — another
argument that Pillar A and Pillar D really share one keystone paper.
Iwasaki–Kimura–Shimomura–Yoshida 1991 (#43) is the textbook reference,
paywalled and non-blocking. The Calligaris–Mazzocco 2018 (#41) and
Diarra–Loray 2020 (#42) papers cover the algebraic-solutions side and
are honest oracles for v0.2 testing.

### Pillar E — adjacent numerical methodology

**PRD's "no released software" claim survives the survey.** The
Olver/Trogdon line ships RiemannHilbert.jl and RHPackage, but those
treat PII only (per PRD line 408–410). The Klein-school multidomain
spectral approach (Klein–Stoilov 2018, #23) is the closest competing
numerical method but ships no general code. The Bornemann Fredholm
determinant method (#25), Wechslberger–Bornemann automatic deformation
(#30), and the Gavrylenko–Lisovyy tau-function CFT route (#6, #33)
each address parts of the problem but none provide pole-field maps
for higher-hierarchy or multi-component systems. Bertola–Cafasso 2012
(#27) and Adler–Sokolov 2021 (#28) / Bobrova–Sokolov 2021 (#29) are
matrix-Painlevé context that PRD lists as v0.3.

**Novokshenov 2014 (#12) is the most important new find** —
*Distributions of Poles to Painlevé Transcendents via Padé
Approximations* (Constr. Approx. 39, 85–99) is a directly
methodology-adjacent paper (Fair–Luke Padé applied to PI/PII/PIV pole
distributions) that v0.1 RESEARCH.md missed. PadeTaylor.jl should
explicitly cite Novokshenov as the closest prior art in v0.2 docs, and
reproduce his pictures as a mutation-proof.

### Pillar F — Calogero–Moser pole dynamics

Three papers are sufficient: Airault–McKean–Moser 1977 (#31) for the
KdV-rational/CM correspondence, Krichever 1980 (#32) for the elliptic
case, and Wilson's bispectral work (referenced but not separately
listed — `arXiv:hep-th/9412124` for the bispectral KP / linearisation
connection). All three give closed-form KdV solutions that the v0.2
CM N-particle smoke test (PRD line 361) cross-validates against. AMM
and Krichever are paywalled; institutional fetch.

### Pillar G — vector Taylor jets infrastructure

No PDFs are required. The substrate is the TaylorIntegration.jl and
TaylorSeries.jl Julia packages (#45 metadata). PRD line 311 already
names `jetcoeffs!`. Reading the source is the right level here, not
fetching papers.

## Open questions surfaced during scoping

- **`arXiv:0704.2869` is misattributed in PRD line 402** (Sasano, not
  Mazzocco). The Stage 0+ writer should fix the citation when the v0.2
  PRD or RESEARCH.md is updated → raise as Stage 0+ open question 1.
- **`arXiv:2107.11680` is misattributed in PRD line 403** (Bobrova–Sokolov,
  not Adler–Sokolov). Same fix → raise as Stage 0+ open question 2.
- **Gavrylenko–Lisovyy 2018 ID in PRD line 405 (`arXiv:1712.08546`)
  is wrong:** that ID is the *Cafasso*–Gavrylenko–Lisovyy *Tau functions
  as Widom constants*. The correct ID for the
  Gavrylenko–Lisovyy CMP 363 paper is `arXiv:1608.00958`. Fetch both
  (#6 and #33) → raise as Stage 0+ open question 3.
- **Novokshenov 2014 is missing from PRD's survey** and is the most
  methodology-adjacent paper to PadeTaylor's own algorithm. v0.1
  RESEARCH.md should be retroactively patched to acknowledge this, and
  v0.2 ADRs should explicitly cite Novokshenov as prior art → raise as
  Stage 0+ open question 4.
- **Mano–Tsuda 2017 is missing from PRD's survey** and is the single
  most-load-bearing v0.2 paper. PRD line 308–309 names Baker–Graves-Morris
  and Beckermann–Labahn but misses the paper that *actually* does
  Hermite–Padé for Painlevé/Garnier. → raise as Stage 0+ open question 5.
- **The Cosgrove preprint has DVI/PostScript only.** Adopting it
  requires `ps2pdf` plumbing in the references workflow. Decide whether
  to ship the `.ps`, the `.pdf`, or both → raise as Stage 0+ open
  question 6.
- **Sakai 2001 connection.** Sakai's Painlevé-classification programme
  is the umbrella under which both Noumi–Yamada and Garnier sit. v0.2
  ADRs may want to use Sakai's classification scheme as the type tag
  for `NoumiYamadaProblem` / `GarnierProblem`. → raise as Stage 0+ open
  question 7.

## Fan-out instructions (for the next wave of subagents)

Five parallel batches, each self-contained. Batch sizes are
intentionally uneven: arXiv fetches are cheap and parallelisable, the
DOI/paywall and book batches require human intervention.

### Batch 1 — arXiv PDFs (open, parallelisable, ~22 papers)

A single Sonnet subagent can fetch all of these via
`curl https://arxiv.org/pdf/<ID>.pdf -o <dest>`:

```
arXiv:1502.06695  -> references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf
arXiv:1306.6161   -> references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf
arXiv:1101.2602   -> references/painleve_hierarchy/ClaeysGrava2011_KdV_painleve_CMP313.pdf
arXiv:math/9808003 -> references/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41.pdf
arXiv:1612.00337  -> references/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf
arXiv:0704.2869   -> references/garnier/Sasano2007_garnier_two_variables_0704.2869.pdf
arXiv:math/0106208 -> references/garnier/Mazzocco2002_geometry_garnier_classical_IMRN_math0106208.pdf
arXiv:0708.0074   -> references/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53.pdf
arXiv:0708.2960   -> references/noumi_yamada/Matsuda2007_rational_A5_NoumiYamada_0708.2960.pdf
arXiv:1811.09274  -> references/noumi_yamada/Clarkson_etal_2020_cyclic_maya_higher_painleve_StudApplMath144.pdf
arXiv:math/0212117 -> references/painleve_hierarchy/JoshiMazzocco2003_tritronquee_PII_hierarchy_math0212117.pdf
arXiv:1107.0214   -> references/painleve_hierarchy/Claeys2011_pole_free_PI_hierarchy_1107.0214.pdf
arXiv:1001.2213   -> references/painleve_hierarchy/Claeys2010_PI2_asymptotics_JPhysA43.pdf
arXiv:1807.04442  -> references/rh_numerical/KleinStoilov2018_multidomain_spectral_painleve_SIGMA14.pdf
arXiv:1111.3527   -> references/rh_numerical/ClaeysOlver2011_higher_tracy_widom_1111.3527.pdf
arXiv:0804.2543   -> references/rh_numerical/Bornemann2010_fredholm_determinants_MathComp79.pdf
arXiv:math-ph/9810007 -> references/rh_numerical/KitaevKorotkin1998_schlesinger_theta_IMRN_mathph9810007.pdf
arXiv:1101.3997   -> references/rh_numerical/BertolaCafasso2012_noncommutative_PII_CMP309.pdf
arXiv:2012.05639  -> references/rh_numerical/AdlerSokolov2021_matrix_PII_TMP207.pdf
arXiv:2107.11680  -> references/rh_numerical/BobrovaSokolov2021_matrix_PIV_test_2107.11680.pdf
arXiv:1206.2446   -> references/rh_numerical/WechslbergerBornemann2014_automatic_RH_painleve_ConstrApprox39.pdf
arXiv:1608.00958  -> references/rh_numerical/GavrylenkoLisovyy2018_isomonodromic_tau_CMP363.pdf
arXiv:1712.08546  -> references/rh_numerical/CafassoGavrylenkoLisovyy2019_widom_tau_CMP365.pdf
arXiv:1010.5725   -> references/noumi_yamada/AratynGomesZimerman2010_integrable_origins_higher_painleve_1010.5725.pdf
arXiv:2009.11668  -> references/noumi_yamada/GomezUllateGrandatiMilson2020_rational_painleve_survey_2009.11668.pdf
arXiv:1401.1408   -> references/painleve_hierarchy/BertolaBothner2015_vorobiev_yablonski_zeros_IMRN.pdf
arXiv:1504.00440  -> references/painleve_hierarchy/BaloghBertolaBothner2017_generalized_vorobiev_yablonski_ConstrApprox45.pdf
arXiv:1707.05222  -> references/noumi_yamada/MasoeroRoffelsen2018_PIV_rationals_poles_SIGMA14.pdf
arXiv:q-alg/9708018 -> references/noumi_yamada/NoumiYamada1999_PIV_symmetries_okamoto_NagoyaMathJ153.pdf
arXiv:1705.03295  -> references/garnier/CalligarisMazzocco2018_finite_orbits_2variable_garnier_JIntegrSyst3.pdf
arXiv:1808.09190  -> references/garnier/DiarraLoray2020_irregular_garnier_algebraic_Compositio156.pdf
arXiv:nlin/0306020 -> references/garnier/KimuraMano2003_irregular_garnier_nlin0306020.pdf
arXiv:1806.08650  -> references/rh_numerical/GavrylenkoIorgovLisovyy2018_FST_system_1806.08650.pdf
```

(31 arXiv URLs total; one fan-out agent.) After download, marker
conversion `marker_single <pdf>` → `references/markdown/<cluster>/<slug>/`.

### Batch 2 — Author preprint pages + format conversion (~3 papers)

A single subagent handles the non-PDF cases:

- Cosgrove 2000-6 — `https://ssh.maths.usyd.edu.au/u/ResearchReports/Nonlinear/Cos/2000-6.html` → fetch the gzipped PostScript, gunzip, `ps2pdf`, place at `references/painleve_hierarchy/Cosgrove2000_USyd_2000-6_higher_painleve_polynomial_classII.pdf`. Sanity-check page count > 80.
- Beckermann–Labahn 1992 preprint — `https://cs.uwaterloo.ca/~glabahn/Papers/tenerife.pdf` → `references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf`.
- Gonnet–Pachón–Trefethen 2011 — `https://people.maths.ox.ac.uk/trefethen/publication/PDF/2011_141.pdf` → `references/hermite_pade/GonnetPachonTrefethen2011_robust_rational_interp_ETNA38.pdf`.

### Batch 3 — DOI / publisher paywalled (~5 papers; institutional access needed)

For these, the subagent should first try ResearchGate / Semantic
Scholar mirrors, then arXiv mirrors, and as a last resort flag for
human institutional-library fetch. Each one is verifiably non-arXiv:

- Novokshenov 2014, Constr. Approx. 39 — DOI 10.1007/s00365-013-9190-6.
- Beckermann–Labahn 1994, SIMAX 15 — DOI 10.1137/S0895479892230031 (Labahn UW page first).
- Airault–McKean–Moser 1977, CPAM 30 — DOI 10.1002/cpa.3160300106.
- Krichever 1980, Funct. Anal. Appl. 14 — DOI 10.1007/BF01078304.
- Kudryashov 1997, Phys. Lett. A 224 — DOI 10.1016/S0375-9601(97)00850-5.
- Sakai 2001, CMP 220 — DOI 10.1007/s002200100446 (check author page on Kobe Univ first).

Output: metadata `.md` stub per paper if PDF unobtainable, with the
full bibliographic record and an `acquisition_status: blocked_by_paywall` line.

### Batch 4 — Books and metadata-only entries (~4 books)

Create `<cluster>/<book-slug>_metadata.md` stubs for each book listed
in the "Books and other non-PDF resources" table above. Each stub
contains: full bibliographic record, ISBN, publisher URL, table of
contents (fetch from publisher page), the chapter(s) most relevant to
v0.2, and an `acquisition_status: institutional_library | siam_member
| ams_member` line. No PDFs.

### Batch 5 — Software metadata

Single short batch:

- `references/taylor_integration/PerezHernandezBenet_TaylorIntegration_metadata.md` —
  Zenodo DOI for v0.4.1 (2562353) and current (7953772+), GitHub URL, the
  version range PadeTaylor v0.2 will target.
- (Optional) `references/rh_numerical/RiemannHilbert_jl_metadata.md` —
  github.com/JuliaApproximation/RiemannHilbert.jl, last-commit date,
  PII-only-coverage note from PRD line 408–410.
