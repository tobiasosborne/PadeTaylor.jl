# Heun Function Implementations — Ecosystem Inventory

Compiled 2026-05-18 for PadeTaylor.jl cross-validation oracle selection.  
All cloned repos live under `external/heun-refs/` and are gitignored (reproducible
via shallow clone from the URLs below).

---

## Summary Table

| Short name | Language | Scope | License | Accessibility | Last updated | Cross-val grade |
|---|---|---|---|---|---|---|
| `motygin-Heun` | MATLAB/Octave | General HeunG, both local solutions, analytic continuation | GPL-3.0 | Public, GitHub | 2023-09-23 | **A — PRIMARY** |
| `motygin-confluent-Heun` | MATLAB/Octave | Confluent HeunC/HeunCs, near z=0, z=1, z=∞ | (unlicensed in repo) | Public, GitHub | 2025-10-27 | **A — PRIMARY** |
| `birkandan-integralseries` | Python 3 | General HeunG, integral-series method, real + complex | (unlicensed in repo) | Public, GitHub (alpha) | 2021-11-13 | **B — SECONDARY** |
| Maple built-in | Maple | HeunG, HeunC, HeunD, HeunB, HeunT — all confluences | Proprietary | Paywalled | c. 2005+ (ongoing) | C — reference only |
| Mathematica built-in | Wolfram Language | HeunG, HeunC, HeunD, HeunB, HeunT + derivatives | Proprietary | Paywalled | Introduced v12.1 (2020) | C — reference only |
| `bamber-HeunC-cpp` | C++ / Python ctypes | Confluent HeunC only, real z, double precision | (unlicensed) | Public, GitHub (prototype) | 2020-11-19 | D — not recommended |
| `birkandan-symODE2` | SageMath / Jupyter | Symbolic ODE classification, Heun recognition — no numerics | (unlicensed) | Public, GitHub | 2022-12-25 | D — symbolic only |
| Fiziev Maple worksheets | Maple | Application worksheets (QNMs, astrophysics) | Unknown | **DEAD** — tcpa.uni-sofia.bg/heun/home.html returns HTTP 404 | Last active ~2015 | — |
| mpmath | Python | **No Heun support** | BSD | Public | ongoing | — |
| SymPy | Python | **No numerical Heun** | BSD | Public | ongoing | — |
| SageMath built-in | Python/Cython | **No numerical Heun** (symODE2 adds symbolic recognition only) | GPL | Public | ongoing | — |
| scipy | Python | **No Heun support** (not on any open wishlist) | BSD | Public | ongoing | — |
| Julia ecosystem | Julia | **No registered Heun package** found on JuliaHub or juliapackages.com | — | — | — | — |

---

## Detailed Entries

### 1. motygin-Heun — MATLAB/Octave, General Heun G

**What it does.**  Provides numerical evaluation of both local solutions of the general Heun
equation (the Fuchsian ODE with four regular singular points at 0, 1, a, ∞).  The primary
interfaces are `HeunL0` (first local solution, normalised to 1 at z=0), `HeunS0` (second
local solution), `HeunLmv`/`HeunSmv` (multi-valued versions along arbitrary paths), and
`HeunL`/`HeunS` (accuracy-improved versions near the singularities z=1, z=a, z=∞).  The
algorithm is power-series analytic continuation (Motygin 2015, arXiv:1506.03848): it avoids
direct ODE integration and achieves ~14 digits near well-conditioned points.

**Language / license.**  MATLAB (100%), tested in Octave.  License: GPL-3.0 (in `LICENSE`
file).

**Accessibility.**  Public, GitHub: https://github.com/motygin/Heun_functions  
Cloned to: `external/heun-refs/motygin-Heun/`  
Last commit: 2023-09-23.

**Scope limitations.**  (a) Not defined for γ a non-positive integer (logarithmic cases —
the code detects and adjusts γ by rounding, see `HeunL0.m` lines 35–39, but no second
solution is then returned).  (b) Returns `NaN` with a diagnostic warning string when
accuracy cannot be guaranteed.  (c) The `HeunLmv` multi-valued path integrator has extra
cost.  (d) Limited to double-precision float arithmetic (Octave/MATLAB inherent).

**Cross-validation suitability.**  Primary oracle for the general Heun G function.
Transparent algorithm, returnable error estimate, GPL, thoroughly described in the companion
paper.  Easily invoked from Julia via `run(`octave ...`)` or by porting the algorithm
directly (algorithm is short and clearly commented).

---

### 2. motygin-confluent-Heun — MATLAB/Octave, Confluent Heun C

**What it does.**  Numerical evaluation of the confluent Heun equation
`HeunC''(z) + (γ/z + δ/(z-1) + ε)·HeunC'(z) + (αz - q)/(z(z-1))·HeunC(z) = 0`.
Provides `HeunC0` (first local solution at z=0, value 1), `HeunCs0` (second local
solution), `HeunC`/`HeunCs` with improved accuracy near z=1, and a full complement of
auxiliary routines covering far-field asymptotics (`HeunCinfA`, `HeunCinfB`) and
connection across the branch cut (`HeunCconnect`, `HeunCjoin0infA`, `HeunCjoin0infB`).
Algorithm: power series + asymptotic expansions + analytic continuation (Motygin 2018,
arXiv:1804.01007).

**Language / license.**  MATLAB (100%), Octave tested.  No LICENSE file in the repo; the
paper cites "copyright 2018, license: GNU GPL v3" (visible in source file headers).

**Accessibility.**  Public, GitHub: https://github.com/motygin/confluent_Heun_functions  
Cloned to: `external/heun-refs/motygin-confluent-Heun/`  
Last commit: 2025-10-27 (actively maintained as of 2025).

**Scope limitations.**  (a) Single-confluent Heun only (not double-, bi-, or
tri-confluent).  (b) Branch cut assumed along [1, ∞); evaluation at points on the cut
requires `HeunCconnect`.  (c) Double precision only.  (d) Logarithmic cases (γ = 0, -1,
-2, ...) require separate treatment via `HeunC00`/`HeunCs00`.

**Cross-validation suitability.**  Primary oracle for the confluent Heun C function.  This
is the most actively maintained open-source Heun implementation found in the survey (last
touched October 2025).  The C++ port by Bamber (entry 6) confirms the algorithm is
portable; a Julia port is straightforward.  The companion paper (arXiv:1804.01007) gives
the full algorithmic description.

---

### 3. birkandan-integralseries — Python 3, General Heun G (integral series)

**What it does.**  Implements the unconditionally convergent integral-series representation
of the general Heun function developed by Giscard & Tamar (arXiv:2010.03919), using the
trapezoidal rule for integration.  Two codes: `HeunG.py` (real argument) and
`HeunG_COMPLEX.py` (complex argument).  A third file `num_heun.py` contains shared
subroutines.  The method is especially efficient when evaluating HeunG at many points with
modest accuracy requirements; it is *not* a power-series method and avoids the convergence-
radius restriction of Motygin's approach.

**Language / license.**  Python 3 (NumPy, SciPy, Matplotlib).  No LICENSE file; README
cites CC BY-NC-ND 4.0 on the arXiv preprint — treat code as unspecified / contact authors.

**Accessibility.**  Public, GitHub: https://github.com/tbirkandan/integralseries_heun  
Cloned to: `external/heun-refs/birkandan-integralseries/`  
Last commit: 2021-11-13.  README explicitly warns "alpha version — please visit regularly."

**Scope limitations.**  (a) General HeunG only; no confluent forms.  (b) Alpha-quality
code: no unit tests, no error estimator, no package structure.  (c) Accuracy trades off
against number of quadrature points (configurable but not automatically chosen).
(d) Matches Mathematica's HeunG convention (normalisations may differ from Motygin's).
(e) Not maintained since November 2021.

**Cross-validation suitability.**  Secondary oracle.  Useful as an independent algorithm
cross-check against Motygin (different mathematical basis) for the general Heun G.  Alpha
quality means results should be verified at known special cases before trusting blindly.

---

### 4. Maple built-in (HeunG, HeunC, HeunD, HeunB, HeunT)

**What it does.**  Maple has offered Heun function support since approximately 2005 (Maple
implementation cited in Motygin 2015 as "the only known" package at the time).  Full family
is provided: HeunG (general), HeunC (singly confluent), HeunD (doubly confluent), HeunB
(bi-confluent), HeunT (tri-confluent), plus the `convert/Heun` facility for recognising
Heun solutions of user-supplied ODEs.

**Language / license.**  Maple (proprietary CAS).  Commercial licence required; no public
download.

**Accessibility.**  Paywalled.  Help pages publicly readable at
https://www.maplesoft.com/support/help/maple/view.aspx?path=HeunG  
No code download possible.

**Known limitations / fragility.**  (a) HeunG not defined when γ is a non-positive integer
(logarithmic cases).  (b) HeunC is reported to use numerical integration for |z| > 1
(outside the series-expansion convergence region), making it very slow in that region
(user report, ResearchGate thread, 2018).  (c) FunctionAdvisor metadata is sparse:
branch-cut locations listed as "unknown" as of Maple 2015 (MaplePrimes thread 203384).
(d) Can return infinity or fail to produce a value at regular points (Motygin 2015
arXiv:1506.03848 reports this explicitly as motivation for his alternative algorithm).

**Cross-validation suitability.**  Reference documentation only.  Not accessible for
programmatic use without a Maple licence.  Known fragility in the |z| > 1 regime means
results should not be accepted uncritically even when Maple is available.

---

### 5. Mathematica built-in (HeunG, HeunC, HeunD, HeunB, HeunT + Prime variants)

**What it does.**  Wolfram Language 12.1 (released May 2020) introduced the complete Heun
family: `HeunG`, `HeunC`, `HeunD`, `HeunB`, `HeunT`, plus corresponding `*Prime` variants
for the derivative.  Supports symbolic manipulation (192 equivalent forms via Möbius
transformations, reduction to hypergeometric 2F1, polynomial cases), and arbitrary-
precision numerical evaluation for complex parameters.

**Language / license.**  Wolfram Language (proprietary CAS).  Commercial licence or
Wolfram Cloud access required.

**Accessibility.**  Paywalled.  Documentation at reference.wolfram.com (currently returning
HTTP 403 for unauthenticated access).  Wolfram blog post: https://blog.wolfram.com/2020/05/06/from-sine-to-heun-5-new-functions-for-mathematics-and-physics-in-the-wolfram-language/

**Known limitations / fragility.**  (a) HeunC: "values of HeunC in logarithmic cases
(nonpositive integer γ) are not determined" (official documentation note).  (b) The
Birkandan et al. 2021 paper (arXiv:2106.13729) explicitly benchmarks against Mathematica
HeunG and finds discrepancies in some parameter regimes, suggesting numerical fragility.
(c) Arbitrary-precision mode is slow; double-precision mode has undocumented accuracy
limits.  (d) Version 12.1 introduced the functions; earlier versions (≤12.0) have no Heun
support at all.

**Cross-validation suitability.**  High-prestige reference, but paywalled and has known
fragility.  The Birkandan et al. code was explicitly designed to reproduce Mathematica's
HeunG normalization — so `birkandan-integralseries` serves as an open proxy for Mathematica
HeunG behavior.

---

### 6. bamber-HeunC-cpp — C++ / Python ctypes, Confluent Heun C

**What it does.**  A C++ port of Motygin's confluent-Heun MATLAB code, with a thin Python
`ctypes` wrapper (`HeunC.py`).  Implements `HeunC` (first local solution) only; no second
solution, no connection routines, no asymptotic expansions.  Real argument z (not
complex), double precision.

**Language / license.**  C++ (96%), Python ctypes wrapper.  Jamie Bamber, 2020.
No license file; no paper reference.

**Accessibility.**  Public, GitHub: https://github.com/JamieBamber/Heun_functions_for_Python  
Cloned to: `external/heun-refs/bamber-HeunC-cpp/`  
Last commit: 2020-11-19.  Zero stars/forks; prototype-quality.

**Scope limitations.**  (a) Confluent Heun only, first solution only.  (b) Real argument
only (despite C++ `std::complex` headers being included, the Python interface passes
real-valued components separately).  (c) Ships a precompiled `.so` binary linked against
the author's system libraries — unlikely to load without recompilation.  (d) No tests,
no documentation beyond the README stub.

**Cross-validation suitability.**  Not recommended.  Motygin's original MATLAB code
(entry 2) is superior in every dimension.  The only value of this repo is as confirmation
that Motygin's algorithm is portable to C++.

---

### 7. birkandan-symODE2 — SageMath / Jupyter, Symbolic only

**What it does.**  A SageMath 9.x package (`ode2analyzer.sage`, `hypergeometric_heun.sage`)
that performs symbolic classification of second-order ODEs with polynomial coefficients:
determines singularity structure, characteristic exponents, and recognises when an ODE
reduces to hypergeometric or Heun form.  Outputs closed-form symbolic solutions for the
hypergeometric case; for Heun equations, outputs the equation parameters only — no
numerical evaluation.

**Language / license.**  SageMath / Jupyter Notebook.  Birkandan (Istanbul Technical
University), 2020.  No license file.

**Accessibility.**  Public, GitHub: https://github.com/tbirkandan/symODE2  
Cloned to: `external/heun-refs/birkandan-symODE2/`  
Last commit: 2022-12-25.  Paper: arXiv:2010.01563.

**Scope limitations.**  Purely symbolic; SageMath itself has no numerical Heun evaluation.
Useful for confirming that a given ODE *is* a Heun equation, not for evaluating it.

**Cross-validation suitability.**  Not a numerical oracle.  Could be useful for automated
ODE-to-Heun-parameter conversion in a preprocessing step.

---

### 8. Fiziev Maple worksheets (DEAD)

**What it does (historically).**  Plamen Fiziev (Dept. of Theoretical Physics, Sofia
University; deceased 2017) and Denitsa Staicova maintained a Heun project bibliography
and Maple worksheet archive at http://tcpa.uni-sofia.bg/heun/home.html.  The worksheets
covered confluent Heun applications to black hole quasinormal modes and relativistic wave
equations, using Maple 17.  Fiziev published extensively on Heun functions from 2009–2016
(arXiv: 0904.0245, 1405.6837, 1512.04025, 1606.08539).

**Language / license.**  Maple worksheets (.mw).  License unknown.

**Accessibility.**  **DEAD.**  URL attempted: http://tcpa.uni-sofia.bg/heun/home.html  
Result: HTTP 404 Not Found (verified 2026-05-18).  No archive.org snapshot accessible
from this environment.  arXiv preprints are accessible; they cite the tcpa.uni-sofia.bg
worksheets but do not embed the code.

**Cross-validation suitability.**  Not accessible.  Documented here for completeness and
to prevent future wasted search time.

---

### 9. mpmath — Python arbitrary precision

**What it does.**  mpmath provides arbitrary-precision floating-point arithmetic and a
comprehensive special function library (Bessel, hypergeometric, elliptic, Mathieu, etc.).

**Heun support.**  None.  Confirmed by searching the mpmath GitHub issues and function
registry.  No open feature request for Heun functions found.

**Accessibility.**  Public, GitHub: https://github.com/mpmath/mpmath  BSD licence.

---

### 10. SymPy — Python symbolic CAS

**What it does.**  SymPy provides a Python symbolic mathematics library.

**Heun support.**  No numerical or symbolic Heun function implemented in the main SymPy
codebase (confirmed by GitHub search).  The Birkandan symODE2 package (entry 7) provides
complementary ODE classification in SageMath but not via SymPy.

**Accessibility.**  Public, GitHub: https://github.com/sympy/sympy  BSD licence.

---

### 11. SageMath built-in

**What it does.**  SageMath wraps multiple CASes (Maxima, Pari/GP, etc.) and provides
its own ring/field machinery.

**Heun support.**  No native numerical Heun evaluation.  The symODE2 package (entry 7)
adds Heun equation *recognition* (symbolic only) to SageMath 9.x.

---

### 12. scipy.special — Python scientific computing

**Heun support.**  None.  No open issue or PR in the scipy tracker requesting Heun
functions was found (GitHub search 2026-05-18).  The ENH: "Generalised ufuncs in
special" PR #20320 does not include Heun.

---

### 13. Julia ecosystem

**Heun support.**  No registered package providing Heun function evaluation found on
JuliaHub or juliapackages.com.  `SpecialFunctions.jl` (the standard registry package)
covers Bessel, Airy, gamma, error functions — not Heun.  `FewSpecialFunctions.jl`
(MartinMikkelsen) is a small collection but does not include Heun.  No discourse.julialang.org
thread requesting Heun was found.

**Implication.**  A Julia Heun implementation in PadeTaylor.jl would be the first in the
Julia ecosystem.

---

## Cross-Validation Recommendations

Ranked by practical usefulness as numerical oracles for PadeTaylor.jl work:

**1. motygin-confluent-Heun (PRIMARY — confluent Heun C)**  
GPL-3.0, publicly available, actively maintained (last commit October 2025), full
algorithm documented in arXiv:1804.01007.  Returns value, derivative, error estimate, and
a diagnostic warning string when accuracy fails — exactly the interface needed for
automated cross-validation.  Call from Julia via `run(`octave --no-gui --eval "..."`)` or
port the power-series continuation algorithm directly to Julia (the code is ~1800 lines
of commented MATLAB, clearly structured).

**2. motygin-Heun (PRIMARY — general Heun G)**  
GPL-3.0, publicly available, last commit September 2023, algorithm documented in
arXiv:1506.03848.  Same quality/interface as the confluent package.  Use this for
cross-validating any general-Heun evaluation in PadeTaylor.jl work.

**3. birkandan-integralseries (SECONDARY — general Heun G, independent algorithm)**  
Python 3, public GitHub, alpha quality but based on a mathematically independent method
(integral series vs. power-series continuation).  Useful as a second independent check
when Motygin and Birkandan agree but Mathematica gives a different answer.  Requires
manual installation of numpy/scipy.  Do not use as a sole oracle; verify at known
special cases (confluent Heun → Mathieu, HeunG → 2F1 limits) before trusting.

**Implementations ruled out for oracle use:**  
- Maple/Mathematica: paywalled, known fragility, not automatable.  
- bamber-HeunC-cpp: subset of Motygin with no added value and likely broken .so.  
- birkandan-symODE2: symbolic only.  
- Fiziev worksheets: dead site, Maple required anyway.  
- mpmath/SymPy/scipy/Julia: no Heun support exists.
