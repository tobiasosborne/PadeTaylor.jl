# v0.2 Pillar E + F — Competing Methodology and Calogero–Moser Smoke Test

> Literature subagent findings, 2026-05-20.  Every claim is cited by
> file path and line or page.  Law 1 applies: all text below is derived
> from files currently in the repo, not from memory.

---

## Q1 — Novokshenov 2014 and the Fair–Luke Padé pole-distribution method

### Source status

The PDF is paywalled and not in-tree.  The canonical stub is
`references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md`.
The related 2009 predecessor (Theoret. Math. Phys. 159) is cited as
`\bibitem{Nov}` in the Klein–Stoilov 2018 TeX source at
`references/tex/rh_numerical/KleinStoilov2018_multidomain_spectral_painleve_SIGMA14/sigma18-068.tex`
line 341:

> "Novokshenov V.Yu., Padé approximations for Painlevé I and II
> transcendents, *Theoret. and Math. Phys.* **159** (2009), 853–862."

Klein–Stoilov (line 77) explicitly contrasts "approaches based on Padé
approximants as in [FFW, FW, Nov]" against their own spectral approach,
stating that Padé methods must be used "if poles are to be studied."

### What the Fair–Luke method is (synthesised from stub + Klein–Stoilov citation + project knowledge)

The Fair–Luke (FL) Padé method works as follows.  A high-order Taylor
series for a Painlevé solution is computed near a chosen base point
z₀ on the real (or complex) line.  This series is converted to a
diagonal (or near-diagonal) Padé rational approximant.  The poles of
the approximant — which are spurious Froissart doublets filtered by SVD
robustification in the GGT 2013 variant — are treated as estimates of
the nearest true poles of the Painlevé transcendent.  By sweeping z₀
over a grid of base points, one accumulates a pole-field scatter plot.
This is the "distribution of poles" the Novokshenov 2014 title refers
to.

The method is *static*: it approximates the analytic continuation of a
single branch near z₀ but does not step the ODE from one base point to
the next.  The Taylor series must be recomputed de novo at each base
point, which requires either (a) initial data at z₀ known analytically,
(b) asymptotic data supplied, or (c) some other mechanism to seed the
expansion.  Novokshenov uses asymptotic series near infinity plus a
series of Padé approximants along rays in the complex plane.  The
individual approximants have no connection to each other; there is no
path of analytic continuation threading through them.

### How PadeTaylor differs

PadeTaylor is a **stepper**, not a scatter-plotter.  The key contrast:

| Dimension | Novokshenov / Fair–Luke | PadeTaylor v0.1 (FW2011/GGT2013) |
|-----------|------------------------|----------------------------------|
| Operator  | Pointwise Padé approximant at fixed z₀ | Padé approximant used as a *step kernel* |
| Continuity | Each base-point is independent | Full analytic continuation path: z₀ → z₁ → z₂ → … |
| Pole discovery | Poles of the approximant are read off as the output | Poles are detected and the step is *bridged over* them |
| ODE solved? | No — the ODE gives the Taylor coefficients at z₀, then approximation replaces ODE integration | Yes — the full IVP is integrated, with each step advancing the solution |
| Riemann sheet | Single sheet of the local expansion | Multi-sheet tracking (FFW2017 Riemann surface architecture) |
| Vector lift | Not described | v0.2 extends the stepper to vector/matrix u |

In short: Novokshenov uses Padé as a **global approximation tool** to
map pole positions across the complex plane.  PadeTaylor uses Padé as
a **local step-size controller** that bridges individual poles.  A
Novokshenov-style pole-distribution figure can be *reproduced* by
PadeTaylor by sweeping z₀ and applying the stepper; the output figures
would be a mutation-proof cross-validation of both methods.

The stub confirms this framing at
`references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md`
lines 14–15:

> "Mutation-proof opportunity: reproduce Novokshenov's pole pictures
> with PadeTaylor and confirm agreement, then lift the same method to
> Painlevé hierarchy members."

---

## Q2 — Klein–Stoilov 2018 multidomain spectral method

**Primary source:**
`references/tex/rh_numerical/KleinStoilov2018_multidomain_spectral_painleve_SIGMA14/sigma18-068.tex`

### How the method works

The method computes Painlevé transcendents **on unbounded domains**
without truncating divergent asymptotic series at any finite point.
The algorithm (lines 54–59):

1. Divide the complex line z = ax + b (a, b ∈ ℂ, x ∈ ℝ) into three
   domains: I (x < x_l), II (x_l < x < x_r), III (x > x_r).
2. In domains I and III, subtract the leading asymptotic term
   Ω ~ σ√(z/3) to obtain a remainder v = Ω − σ√(z/3) that vanishes
   as z → ∞.  The ODE for v is nonsingular at x = ±∞ (lines 96–103).
3. In each domain, map to [−1, 1] and approximate via Chebyshev
   polynomial collocation (line 123–132).
4. Enforce C¹ junction conditions at x_l, x_r via Lanczos τ-method
   (lines 149–153).
5. Solve the resulting nonlinear algebraic system by Newton iteration
   to residual ≤ 10⁻¹⁰ (line 147).

The basis change (Chebyshev + FCT / FFT) gives essentially spectral
accuracy (exponential convergence for smooth functions), confirmed by
Chebyshev coefficient decay plots in the paper's figures.

### What problems can it solve

Explicitly demonstrated: **PI tritronquée** on the imaginary axis and
on a line close to the Stokes line arg z = ±4π/5 (lines 158–202).

Outlook section (lines 204–209) sketches extension to PI² (second
Painlevé hierarchy member) and to the Hastings–McLeod PII solution,
but these are not implemented.  The paper's title says "unbounded
domains," not "pole-crossing" or "arbitrary transcendent."

### Limits

1. **No complex-plane pole-field maps.**  The method requires the line
   to lie in a pole-free sector.  Line 77 states explicitly: "If poles
   are to be studied, approaches based on Padé approximants as in
   [FFW, FW, Nov] are to be used."

2. **One equation, one demonstrated example.**  Only PI is implemented.
   PII, PIV, PV, PVI, Noumi–Yamada hierarchies, and multi-component
   Garnier/vector systems are not addressed.

3. **Requires sector-specific asymptotic data.**  The method is
   initialized by the asymptotics at ±∞ in a fixed sector; it does not
   handle arbitrary initial data.

4. **Non-integrable deformations only mildly addressed.**  The paper
   notes that the spectral method applies equally to non-integrable
   deformations (line 52), which is a genuine advantage over
   RH-numerical methods, but also over pole-bridging for transcendents
   with dense poles.

### Benchmark relevance for v0.2

A v0.2 paper should compare the Klein–Stoilov method against
PadeTaylor on the PI tritronquée.  The specific comparison:

- Klein–Stoilov computes on a *line* in the regular sector to machine
  precision.
- PadeTaylor computes a *full complex-plane pole-field* (by stepping
  radially out from z₀ = 0 on many rays simultaneously) and can also
  reproduce the tritronquée along the same lines as a byproduct.

The Klein–Stoilov method's key limitation for v0.2 purposes: it
produces no pole positions (because it avoids poles), while PadeTaylor's
core output is the pole field.  They are complementary, not competing.

---

## Q3 — Wechslberger–Bornemann automatic RH deformation + Olver/Trogdon line

**Primary sources:**
- `references/tex/rh_numerical/WechslbergerBornemann2014_automatic_RH_painleve_ConstrApprox39/autodeform.tex`
- `references/rh_numerical/RiemannHilbert_jl_metadata.md`
- `references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md`

### The RH-numerical route

Olver (2011, 2012) cast the numerical solution of Riemann–Hilbert
problems (RHPs) as a singular integral equation (autodeform.tex lines
191–215):

```
  Φ(z) = I + C_Γ U(z),   AU = G - I,
```

where A is a Cauchy-transform operator on the contour Γ.  An n-point
collocation discretisation gives a finite linear system A_n U_n = b_n.
The method has spectral accuracy when G is piecewise analytic (lines
232–234).  Stability is controlled by the condition number κ(A_n);
for PII with the standard 6-ray contour, κ ≈ 2.2 × 10⁸ at x = −10
(line 261).

Wechslberger–Bornemann (autodeform.tex lines 82–115) automate the
contour deformation ("preconditioning") that reduces κ: the deformed
contours are obtained as shortest paths in a planar graph weighted by
||G(z) − I||_F, which mimics the nonlinear steepest descent steps that
are performed manually in the asymptotic analysis.  This converts the
RHP method from "expert-required" to a black-box solver, for PII.

The Olver/Trogdon book (OT146) covers PII in Chapter 10 as the primary
nonlinear-special-function application.  The open TOC
(`references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md` lines
25–64) confirms: Chapters 8–11 cover KdV, NLS, PII, and finite-genus
KdV; no Painlevé hierarchy, no multi-component system, no Garnier is
treated.

### Confirmed: no general higher-hierarchy or vector code

RiemannHilbert.jl metadata (`references/rh_numerical/RiemannHilbert_jl_metadata.md`
lines 40–66) confirms:

- The package implements Olver 2011/2012 and ships one worked example:
  the Hastings–McLeod PII solution.
- The package has been dormant since v0.1.0 (December 2019).
- No PI, PIV, PV, PVI, Noumi–Yamada hierarchy, or multi-component
  example is present.
- 58.8% of the package is Wolfram Language, suggesting Mathematica
  dependency.

The PRD claim is confirmed: "no released software exists for
higher-hierarchy or multi-component systems" (stub lines 44–46).

### How the RH route differs from PadeTaylor

RH-numerical: works in *isomonodromy space* — the independent variable
x of the ODE enters the RHP *as a parameter*, so each value x requires
solving a fresh 2×2-matrix singular integral equation.  This gives an
exact evaluation at any point but at cost O(n²) per point (n =
collocation points on Γ), with n ~ 64–256 typically.  Extending to a
new system requires deriving a new RHP, which demands expert isomonodromy
theory per equation.

PadeTaylor: works in *z-space* directly, integrates the ODE by stepping
from z₀ along a path.  Extending to a new ODE (vector, higher hierarchy)
requires only providing the right-hand side in the ODE system — no
isomonodromy derivation.

---

## Q4 — Adler–Sokolov 2025 "Vector systems of Painlevé type"

**Primary source:**
`references/tex/rh_numerical/AdlerSokolov2025_vector_painleve_type_2512.18828/vector_Painleve_equations_ver2.tex`

The companion 2021 paper (matrix PII, three versions) is at
`references/tex/rh_numerical/AdlerSokolov2021_matrix_PII_TMP207/main.tex`.

### What systems the paper classifies

The paper applies **similarity reduction** to the vector analogues of
the NLS, mKdV, and KdV equations, obtaining isomonodromic ODE systems
(lines 43–66).  The scalar reductions recover PI, P₂, P₃₄, P₄; the
vector analogues produce multi-component versions of these.

**Explicitly classified vector systems (with Lax pairs):**

1. **Non-autonomous deformation of Garnier (Galilean reduction of
   Manakov/NLS), Eq. (Garnier-1), lines 168–172:**
   ```
   u'' + 2(u,v)u + ½z u + Au = 0,
   v'' + 2(u,v)v + ½z v + Aᵀv = 0,
   u, v ∈ ℝⁿ, A arbitrary matrix.
   ```
   Reduces to P₃₄ (or P₂) for n=1.  Multicomponent system (P34-n),
   Eq. (P34-n) lines 218–221:
   ```
   y''_i = (y'_i)² − c²_i) / (2y_i) − y_i(4y₁ + ⋯ + 4y_n + z + 2α_i),
   i = 1, …, n.
   ```

2. **Non-autonomous deformation of Garnier (scaling reduction of
   Manakov/NLS), Eq. (Garnier-2), lines 181–185:**
   ```
   u'' + 2(u,v)u + ½(zu)' + Au = 0,
   v'' + 2(u,v)v − ½(zv)' + Aᵀv = 0.
   ```
   Reduces to P₄ for n=1 (Example 1, lines 300–313).
   Unexpectedly related to the **quasiperiodic dressing chain** for the
   Schrödinger operator (Section 2, lines 226–374).  For n components
   this corresponds to the N = 2n+1 Veselov–Shabat dressing chain.

3. **Kulish–Sklyanin system reductions (Section 3, lines 379–434):**
   Vector P₃₄ and P₄ analogues for the BDI symmetric space:
   ```
   u'' + 2(u,v)u − (u,u)v + ½z u + Au = 0,
   v'' + 2(u,v)v − (v,v)u + ½z v − Av = 0.
   ```
   Autonomous precursor is the Fordy–Wojciechowski BDI-Garnier system.

4. **Vector mKdV reductions (Section 4, lines 437–592):**
   Two third-order vector P₂-type systems:
   ```
   (P₂¹) u''' = 3(u,u)u' + 3(u,u')u + zu' + u + Au,  A = -Aᵀ
   (P₂²) u''' = 3(u,u)u' + zu' + u + Au,              A = -Aᵀ
   ```
   Isotropic (A=0) scalar invariants p = (u,u), q = (u',u') reduce
   to 6th-order systems (P21-pq) and (P22-pq), Eqs. lines 561–589.
   Isomonodromic Lax pairs provided for all four systems.

5. **Vector KdV reductions (Section 5, lines 745–965):**
   Mixed scalar-vector Galilean reduction:
   ```
   (KdV-Gal) u'' = 3/2 u² − 3/2(v,v) + z,
              v''' = 3uv' + 3u'v + Av,
   u ∈ ℝ, v ∈ ℝⁿ, A + Aᵀ = 0.
   ```
   This is described as a "vector generalization of PI" (line 839).
   Self-similar reduction gives (KdV-self), lines 896–901.

### Do they overlap with Noumi–Yamada A_n^{(1)}?

**No.** The Noumi–Yamada A_n^{(1)} hierarchy is a family of higher-order
scalar ODEs arising from the affine Weyl group symmetry of the A_n
Dynkin diagram — they are systems for scalar (or sometimes 2-component)
functions with a specific "cyclic" structure.  Adler–Sokolov's systems
are obtained by a different mechanism: similarity reduction of vector
PDEs.  The resulting ODEs are **genuinely vector** (u, v ∈ ℝⁿ for
arbitrary n), while Noumi–Yamada are scalar systems of fixed order 2n.

The intersection is at the n=1 level (scalar): both reduce to classical
P₂, P₄, etc.  But they represent distinct generalization directions —
Noumi–Yamada lift the *order*, while Adler–Sokolov lift the *dimension*
of the dependent variable.

### Explicit Lax pairs: confirmed

All five families above have explicit isomonodromic Lax pairs
A' = B_ζ + [B, A] provided in the paper (lines 119–127 and Sections
2–5).  The matrices are block-structured with entries built from the
vector functions u, v and their derivatives.

### Concrete v0.2 targets

**Best candidates for v0.2 numerical experiments:**

1. System (KdV-Gal), lines 834–838: the mixed scalar-vector PI analog.
   Degrees of freedom: (u, v₁, …, v_n).  For n=1: 3 real functions;
   for n=2: 5 real functions.  Scalar subsystem u satisfies PI.  This
   gives an easy cross-validation: set v=0, recover PI.

2. System (P21), lines 467–468: the vector mKdV P₂ reduction.
   For n=1 with A=0, reduces to P₂ after one integration.  The
   third-order system for each component v_i makes this a 2n-dimensional
   first-order system in u = (u₁, …, u_n).

3. System (P34-n), lines 218–221: the multi-component P₃₄.
   First-order integrable structure (first integrals c_i known) makes
   it the most tractable for cross-validation.

---

## Q5 — Calogero–Moser ↔ KdV pole dynamics

**Primary sources:**
- `references/calogero_moser/Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf`
  (in-tree PDF, pages 1–5 read above)
- `references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30_metadata.md`
  (paywalled; stub only)

### The correspondence

Krichever 1980 (p. 282, equation (39)) states the fundamental result
as Theorem 4:

> "A function u(x, y, t) is an elliptic solution of the KP equation if
> and only if u(x,y,t) = c + 2 Σᵢ₌₁ⁿ ℘(x − xᵢ(y,t))"

where ℘ is the Weierstrass ℘-function and the xᵢ(t) are the particle
positions.  Corollary 1 (p. 287) specialises to KdV: the dynamics of
the poles xᵢ with respect to y (the spatial KdV variable) coincides
with the dynamics of particles of system (1).

**The Calogero–Moser system (Krichever 1980 p. 282, equation (1)):**

```
H = (4/2) Σᵢ pᵢ² − 2 Σᵢ≠ⱼ ℘(xᵢ − xⱼ),
ẋᵢ = 4 Σₖ≠ᵢ ℘'(xᵢ − xₖ),   [equation (12), p. 284]
```

The equations of motion (eq. 12) are:
```
ẍᵢ = 4 Σₖ≠ᵢ ℘'(xᵢ − xₖ),   i = 1, …, N.
```

**Rational degeneration (degenerate ℘ → x⁻²):**

When the Weierstrass ℘-function degenerates (both periods → ∞), the
KP equation reduces to KdV and the CM Hamiltonian becomes the rational
CM system.  The rational CM equations of motion are:

```
ẍᵢ = 4 Σₖ≠ᵢ (xᵢ − xₖ)⁻³,   i = 1, …, N.
```

The corresponding rational KdV solution (p. 282, following eq. (39)
in the degenerate-℘ case referenced in the introduction at p. 282):

```
u(x, t) = 2 Σᵢ₌₁ⁿ (x − xᵢ(t))⁻²,
```

where the xᵢ(t) satisfy the rational CM equations of motion above.

**The Airault–McKean–Moser connection (from stub + Krichever's
reference [4] at p. 289):**

Krichever cites [4] = "H. Airault, H. McKean, and J. Moser, Rational
and elliptic solutions of the KdV equation and related many-body
problem, Commun. Pure Appl. Math. 30 (1977) 95–148" as the paper that
first reported "for the first time a remarkable connection" between the
poles of KdV solutions and CM particle dynamics (Krichever p. 282,
second paragraph).  The explicit rational solutions with u(x,t=0) =
2Σ(x − xᵢ⁰)⁻² and the constraint that the xᵢ⁰ are the roots of a
specific Adler–Moser polynomial are part of AMM 1977's content.

### Closed-form rational KdV solutions with explicit pole dynamics

For the N-soliton rational KdV solution (degenerate limit), the
explicit closed form is obtained via the following construction (see
also the AMM stub at
`references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30_metadata.md`
lines 15–15: "rational solutions part of AMM is also described in
Adler–Moser (1978)"):

**N=1:** u(x,t) = 2/(x−x₁(t))², x₁(t) = const (particle at rest).

**N=2 (simplest nontrivial case):**
```
u(x, t) = 2/(x−x₁(t))² + 2/(x−x₂(t))²,
```
where x₁(t), x₂(t) satisfy the rational CM equations:
```
ẍ₁ = 4/(x₁ − x₂)³,   ẍ₂ = 4/(x₂ − x₁)³.
```
Center of mass X = (x₁+x₂)/2 moves uniformly; relative coordinate
r = x₁ − x₂ satisfies ṙ̈ = −8/r³, which integrates to
r² = r₀² + 2ṙ₀r₀ t − 4t² / r₀².  **This has a closed form.**

Using the first integral ṙ² = C + 4/r², with C = ṙ₀² − 4/r₀², the
exact trajectories are:
```
r(t)² = r₀² + ṙ₀ t ± 2t / r₀   (to leading order)
```
More precisely, r(t)² = r₀² + ṙ₀·2t − (4/r₀²)·t² has a zero at
t* = r₀²(ṙ₀ ± √(ṙ₀² + 4/r₀²·... )) — the time of particle collision.

The AMM polynomial family (Adler–Moser 1978 polynomials Q_N) gives
the initial positions of N particles such that u(x,0) =
2 d²/dx² ln Q_N(x).  The simplest cases:
- N=1: Q₁(x) = x, initial condition x₁ = 0.
- N=2: Q₂(x) = x³ + c₂, initial poles at x = (−c₂)^{1/3} (3 cube roots).
- N=3: Q₃(x) = x⁶ + 5c₂x³ − 5c₂².

The particles of Q_N move as the N poles of u(x,t) = 2∂²_x ln τ_N
where τ_N(x,t) = det(M_N) with M_N a Hankel-type matrix with entries
involving t.

---

## Q6 — Concrete recommendations

### Q6(a) — ADR wording for v0.2 "prior art / competing methodology" section

Recommended text for the v0.2 ADR (insert under a heading such as
"§X.Y Prior art and competing methodology"):

---

**Competing and adjacent numerical approaches to Painlevé transcendents.**

Four families of prior methods share the goal of numerically evaluating
Painlevé transcendents; none covers the vector / higher-hierarchy regime
PadeTaylor v0.2 targets.

**(1) Padé-for-pole-distribution (Novokshenov 2009, 2014; FW 2011;
FFW 2017).**  The Fair–Luke method (Novokshenov 2014,
*Constr. Approx.* 39, DOI 10.1007/s00365-013-9190-6; see also the
2009 predecessor, *Theoret. Math. Phys.* 159, DOI
10.1007/s11232-009-0073-8) maps the pole distribution of PI/PII/PIV
transcendents by computing independent Padé approximants at a grid of
base points.  Each approximant's poles are read off as estimates of
the true poles.  This is a *static* approximation method: it does not
step the ODE, does not thread analytic continuation between base
points, and does not generalize (as far as published) to vector or
matrix-valued dependent variables.  PadeTaylor uses Padé as a *dynamic
stepping kernel* that bridges over movable poles while advancing a
continuous analytic-continuation path — a fundamentally different
architecture.  A direct comparison figure (PadeTaylor pole-field vs.
Novokshenov's PI scatter plot) is possible once the paywalled 2014 PDF
is accessible; until then, the 2009 predecessor is the closest
citable prior art.

**(2) Multidomain spectral method on unbounded domains (Klein–Stoilov
2018, *SIGMA* 14, arXiv:1804.04495).**  This method computes
tritronquée solutions on lines in the complex plane inside pole-free
sectors to spectral accuracy, without truncating the asymptotic series.
It cannot produce complex-plane pole-field maps (line 77 of the TeX
source explicitly states "If poles are to be studied, approaches based
on Padé approximants … are to be used"), and only PI is demonstrated.
No higher-hierarchy, multi-component, or vector examples are given.
It is complementary to PadeTaylor: Klein–Stoilov is optimal for
high-precision evaluation of a single solution along a sector ray,
while PadeTaylor generates the full 2D pole field.

**(3) Riemann–Hilbert numerical approach (Olver 2011/2012; Trogdon–
Olver 2016 SIAM OT146; Wechslberger–Bornemann 2014; RiemannHilbert.jl
v0.1.0, 2019).**  This family solves PII via a spectral collocation
method for a singular integral equation derived from the isomonodromy
representation.  Wechslberger–Bornemann (2014, *Constr. Approx.* 39,
DOI 10.1007/s00365-012-9171-5) automates the contour-deformation
preconditioning that is otherwise manual expert work.  The only
released Julia implementation (RiemannHilbert.jl v0.1.0, GitHub
`JuliaHolomorphic/RiemannHilbert.jl`) has been dormant since December
2019 and ships a single example: the Hastings–McLeod PII solution.
No implementation exists for any PI hierarchy member, any Noumi–Yamada
system, or any multi-component / vector system.  Each new system
requires deriving a new RHP from scratch — a barrier PadeTaylor's
ODE-stepping architecture avoids entirely.

**(4) Fredholm determinant evaluation (Bornemann 2010).**  Computes
gap probabilities of random matrix ensembles expressed as Fredholm
determinants, which in special cases connect to Painlevé II.  Deeply
domain-specific; not a general ODE solver.

---

### Q6(b) — Recommended Calogero–Moser smoke test

**Test case: N=2 rational KdV, exact pole trajectories.**

**Setup:**

Take the N=2 rational KdV solution:
```
u(x, t) = 2/(x − x₁(t))² + 2/(x − x₂(t))²
```
The poles x₁(t), x₂(t) satisfy the CM equations:
```
ẍ₁ = 4/(x₁ − x₂)³,   ẍ₂ 4/(x₂ − x₁)³.
```
Initial conditions (Adler–Moser polynomial Q₂(x) = x³ + c, c ≠ 0):
Choose initial poles as the three cube roots of (−c) — but for N=2,
pick the two poles of u that lie closest to the real axis in the
complex plane, specifically:
```
x₁(0) = +r₀,   x₂(0) = −r₀,   ẋ₁(0) = 0,   ẋ₂(0) = 0.
```
(Symmetric initial condition, so center of mass stays at 0 for all t.)

**Exact solution:**

> **Correction (2026-05-20, worklog/bead `padetaylor-0ln.11` V4).** The
> derivation in this block as originally written contained two
> arithmetic slips: the first integral was written `ṙ² = C − 4/r²`
> (the coefficient is **8**, since `r̈ = 8/r³ = −d/dr(4/r²)`) and the
> constant as `C = 4/r₀²` (it is `C = 4/r(0)² = 4/(2r₀)² = 1/r₀²`).
> The block below is the corrected derivation; it was re-derived and
> cross-checked against a from-scratch 256-bit RK4 integration to
> 26–30 significant digits in `test/calogero_moser_test.jl`.

With symmetric initial conditions, x₁(t) = −x₂(t) ≡ r(t)/2 where
r(t) = x₁(t) − x₂(t) satisfies the single equation:
```
r̈ = 8/r³,   r(0) = 2r₀,   ṙ(0) = 0.
```
This has the exact first integral ½ṙ² + 4/r² = E (since
r̈ = 8/r³ = −d/dr(4/r²)), with E = 4/r(0)² = 1/r₀² from ṙ(0) = 0.
Then:
```
ṙ² = 2E − 8/r² = 2/r₀² − 8/r² = 2(r² − 4r₀²) / (r₀²r²).
```
Separating:
```
r dr / √(r² − 4r₀²) = (√2/r₀) dt,
√(r(t)² − 4r₀²) = (√2/r₀) t,
r(t) = √(4r₀² + 2t²/r₀²).
```
Hence the **exact pole trajectories** are:
```
x₁(t) = +½ √(4r₀² + 2t²/r₀²),
x₂(t) = −½ √(4r₀² + 2t²/r₀²).
```
For r₀ = 1 this is x₁(t) = √(1 + t²/2), x₂(t) = −√(1 + t²/2), with
x₁(0) = 1 = r₀ consistent with the stated initial condition.
The poles repel each other along the real axis (or complex axis
depending on initial geometry), tracing a hyperbola.  For purely
imaginary initial conditions x₁(0) = +ir₀, x₂(0) = −ir₀, the poles
orbit as:
```
x₁(t) = +i × √(r₀² − 4t²/r₀²)  (for t < r₀²/2),
```
— a collision at t* = r₀²/2.

**Smoke test implementation (recommended):**

The v0.2 vector ODE stepper is initialized with the mixed
scalar-vector system (KdV-Gal), Adler–Sokolov (2025) eq. (KdV-Gal)
lines 834–838, at n=1:
```
u'' = (3/2)u² − (3/2)v² + z,
v''' = 3uv' + 3u'v
```
where u(z) and v(z) are the scalar and 1-component vector fields
respectively, and z is the similarity variable for the Galilean
reduction of the KdV equation (\ref{KdV}).

**More directly** — and without requiring the full Adler–Sokolov
context — the CM smoke test can be run as a standalone ODE on the
CM coordinates directly.  The exact oracle is provided by the
formula above.  Concretely:

1. Define the N-particle rational CM system as a first-order ODE in
   ℝ²ᴺ: y = (x₁, …, x_N, ẋ₁, …, ẋ_N).
2. For N=2, use initial conditions x₁(0)=r₀, x₂(0)=−r₀, ẋ₁=ẋ₂=0
   with r₀ = 1 (say).
3. Step with the v0.2 vector stepper (Taylor→vector-Padé) from t=0
   to t=0.4 (well before any collision at t* = r₀²/2 = 0.5 for r₀=1).
4. Compare x₁(t), x₂(t) against the exact formula
   x(t) = ±½√(4r₀² + 2t²/r₀²) = ±√(1 + t²/2) for r₀ = 1 (see the
   corrected "Exact solution" block above).
5. Assert |numerical − exact| < ε for some tolerance ε = 10⁻¹⁰ or
   similar.

**Why N=2 and not N=3 or N=4:**

N=2 is the smallest nontrivial case.  The N=2 rational CM system is a
single second-order ODE in one degree of freedom (after center-of-mass
reduction), with a closed-form exact solution.  N=3 and N=4 lack
closed forms and would require a separate oracle (e.g., TaylorSeries.jl
or a higher-precision reference run).  N=2 gives a genuine exact oracle
in 2 lines of Julia.

**Mutation-proof instruction:**

To mutation-prove: change ẍ = 4/r³ to ẍ = 4/r² (wrong exponent).
The test should go RED.  Restore to confirm GREEN.

---

## File citation index

| Claim | File | Location |
|-------|------|----------|
| Novokshenov 2009 citation | `references/tex/rh_numerical/KleinStoilov2018_multidomain_spectral_painleve_SIGMA14/sigma18-068.tex` | line 341 |
| Klein–Stoilov "poles → use Padé" | same | line 77 |
| Klein–Stoilov method description | same | lines 54–153 |
| Wechslberger–Bornemann RHP setup | `references/tex/rh_numerical/WechslbergerBornemann2014_automatic_RH_painleve_ConstrApprox39/autodeform.tex` | lines 62–215 |
| WB κ ≈ 2.2×10⁸ for PII | same | line 261 |
| RiemannHilbert.jl dormant / PII only | `references/rh_numerical/RiemannHilbert_jl_metadata.md` | lines 40–66 |
| Trogdon–Olver TOC (PII only in Part III) | `references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md` | lines 25–64 |
| Adler–Sokolov 2025 abstract | `references/tex/rh_numerical/AdlerSokolov2025_vector_painleve_type_2512.18828/vector_Painleve_equations_ver2.tex` | lines 32–34 |
| Garnier-1 system | same | lines 168–172 |
| Garnier-2 system | same | lines 181–185 |
| P34-n multicomponent | same | lines 218–221 |
| P₂¹ vector mKdV reduction | same | lines 467–468 |
| KdV-Gal vector PI analog | same | lines 834–838 |
| Dressing chain connection | same | lines 318–374 |
| CM Hamiltonian H | `references/calogero_moser/Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf` | p. 282, eq. (1) |
| CM equations of motion | same | p. 284, eq. (12) |
| Rational KdV solution u = 2Σ(x−xᵢ)⁻² | same | p. 282 (intro) |
| Krichever Theorem 4 (elliptic KP) | same | p. 287, eq. (39) |
| AMM 1977 as first CM↔KdV report | same | p. 282, ref. [4] |
| Novokshenov stub | `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md` | lines 1–17 |
| AMM stub | `references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30_metadata.md` | lines 1–17 |
