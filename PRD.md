# PRD: PadeTaylor — General-Purpose Taylor–Padé IVP Solver

> **Status: SKETCH (v0.1).** First draft, written from memory and a single web pass.
> Do **not** start coding from this. **Stage 0 (deep research) is the first step**;
> everything below is provisional and will be replaced by `RESEARCH.md` + `DESIGN.md`.

> **Note (2026-05-19):** v0.1 has shipped — `RESEARCH.md` + `DESIGN.md`
> + 18 ADRs + 51 worklog shards exist, the package handles all six
> scalar Painlevé equations with pole-field maps and Riemann-sheet
> tracking, and the test suite is GREEN at 2447+ assertions. The
> **v0.2 north star — vector / multi-component Painlevé-type systems**
> is appended at the bottom of this file. The v0.1 sketch below
> remains as historical context.

---

## Goal

A general-purpose Taylor–Padé IVP solver for analytic ODEs `y' = f(z, y)`,
`z ∈ ℂ`, supporting:

- Movable poles — pole-immune stepping in the Fornberg–Weideman tradition.
- Arbitrary precision and Float64.
- Complex-analytic problems beyond Painlevé (general `f` built from elementary functions).

Two reference implementations, mutually validating:

1. **Julia**: `PadeTaylor.jl` — arb-prec via `Arblib`, Float64 generic.
2. **Scientist Workbench**: TS-native package under `@cas/`, arb-prec via `@cas/cas-core` primitives.

Goal is for these to be **the reference** — well-specified, well-tested,
the thing others port from. Without 1:1 cross-validation between the two,
neither qualifies.

---

## Stage 0 — Deep research (do this BEFORE planning)

Output: `RESEARCH.md` summarising findings before any design or code.

### Primary literature (read carefully, store PDFs in `references/`)

- **Fornberg & Weideman**, *A numerical methodology for the Painlevé equations*,
  J. Comput. Phys. 230 (2011) 5957–5973. The foundational paper. Algorithmic detail in §3–4.
- **Fornberg & Weideman**, *A computational exploration of the second Painlevé equation*,
  FoCM 14 (2014) 985–1016. Refinements, PII survey.
- **Fasondini, Fornberg, Weideman**, J. Comput. Phys. 344 (2017) 36–50.
  Extension to multi-valued transcendents on Riemann surfaces.
- **Reeger & Fornberg** — PIV survey paper(s).
- **Willers (1974)** — original IVP-with-poles algorithm; FW generalises this.

### Padé robustness (essential)

- **Gonnet, Güttel, Trefethen**, *Robust Padé Approximation via SVD*, SIAM Review 55 (2013).
  This is **the** Padé routine. Read fully before implementing anything.
- Chebfun's `padeapprox` source — reference impl of GGT 2013.
- Hermite–Padé literature (Baker & Graves-Morris) — for v2 / branch-point detection.

### Adjacent Taylor-IVP packages (study designs, do not reinvent)

- **TIDES** (Barrio et al.) — mature high-order Taylor.
- **TaylorIntegration.jl** (Pérez-Hernández, Benet) — Julia; real-axis only, no Padé.
- **TaylorSeries.jl** (Benet, Sanders) — coefficient layer; assess generic-type support
  for `Arb` element type, especially at order 60+.
- **ATOMFT** (Chang & Corliss) — historical reference.
- **Jorba & Zou (2005)** — Taylor IVP package paper, design lessons.

### Validated / arb-prec ODE (adjacent, not duplicates)

- **Mezzarobba** — D-finite / linear ODEs in Arb (`ore_algebra`). Different scope, useful for verified-arithmetic patterns.
- **COSY Infinity** / Taylor models (Makino & Berz) — for verified-bound thinking in v2.

### Source code hunt

- Paper appendices — FW 2011 may include code excerpts.
- Weideman's Stellenbosch page (`appliedmaths.sun.ac.za/~weideman/`) — PFS code does not appear to be posted; confirm.
- arXiv ancillary files for FW papers.
- **Email** Weideman, Fornberg, Fasondini, Reeger. They have working MATLAB and have shared with collaborators historically.

### Ecosystem inventory

- `@cas/cas-core`: does a Taylor / formal power series class exist? Does a bigfloat SVD exist? What's missing?
- Julia: confirm `Arblib.jl` exposes SVD on `ArbMatrix`; confirm `TaylorSeries.jl` is generic enough to use `Arb` element type at order 60+ without breakage.

### Stage 0 deliverable

`RESEARCH.md` containing:

- Confirmed algorithmic spec (Taylor → Padé → step control → optional path network).
- All design choices found in literature with explicit tradeoffs.
- What pieces exist in the ecosystem; what we must build.
- Updated open questions (some below will resolve; new ones will appear).
- Recommendation: proceed to Stage 1, or further research needed.

---

## Stage 1 — Design decisions to lock (after research)

Output: `DESIGN.md` resolving at least:

1. **Padé form**: scalar per component vs Hermite–Padé with shared denominator. (v1 likely scalar; document the v2 path.)
2. **Step control**: coefficient-decay heuristic, Padé-denominator-root distance, or both switchable.
3. **Order strategy**: fixed (probably 30–60) vs adaptive.
4. **Path strategy in v1**: real-axis only with primitive `step_complex` API, vs 2D network from day one. (v1 likely real-axis only.)
5. **API surface**: standalone clean core + thin SciML `CommonSolve` wrapper for Julia.
6. **Verified arithmetic** via Arb balls — v1, v2, or never.
7. **Module boundaries** aligned with the ~200 LOC/file budget.

---

## Algorithmic core (sketch — fill in after research)

```
1. TaylorStep(f, z_n, y_n, N)   → coefficients a_0, ..., a_N
                                   (Taylor-mode AD on operator-overloaded polynomial type)
2. RobustPade(coeffs, [L/M])    → (P, Q)  via Gonnet–Güttel–Trefethen SVD
3. StepControl(coeffs, P, Q, τ) → dz_next
4. Step:  y_{n+1} = (P/Q)(z_n + dz_next)
5. (v2) PathNetwork: graph of steps with consistency checks at meeting points.
```

---

## Implementations

### Julia (`PadeTaylor.jl`)

- Deps: `TaylorSeries.jl`, `Arblib.jl`, `LinearAlgebra` (or generic SVD for `Arb`), `Polynomials.jl` (roots of `Q`).
- Modules (each ≤ ~200 LOC):
  - `RobustPade` — port of GGT 2013.
  - `PadeStepper` — coefficient generation + Padé + step control.
  - `Problems` — IVP problem and solution types.
  - `PathNetwork` (v2).
- Generic over `T <: Number`; `Float64` and `Arb` both first-class.

### Scientist Workbench (TS, monorepo)

- Package: `@cas/pade-taylor` (name TBD).
- **Mirror Julia module structure 1:1.** This is non-negotiable; mirroring is what makes them mutually a reference impl.
- Depends on `@cas/cas-core` for arb-prec primitives, polynomial type, bigfloat SVD.
- **Bigfloat SVD is the heaviest unknown.** If `@cas/cas-core` lacks it: scope as a sub-task and port a known-good algorithm (Demmel–Kahan, or one-sided Jacobi for arb-prec). **Do not hand-roll a fragile version.**

### Cross-validation

Shared test corpus; both impls run it; diff to N digits. Without 1:1 mirroring and shared corpus, neither is a reference impl.

---

## Test corpus

- PI tritronquée (Boutroux IC) — pole locations vs FW 2014 figures.
- PII Hastings–McLeod and Ablowitz–Segur — vs FW 2014 / Reeger–Fornberg.
- Riccati with closed-form meromorphic solution — analytic pole locations as ground truth.
- A Schwarzian / non-Painlevé analytic ODE — proves "general purpose."
- A bog-standard non-stiff problem on the real axis (Lorenz segment, Van der Pol away from stiff regime) — sanity check vs `DifferentialEquations.jl`.
- DLMF / Mathematica `PainleveT` tabulated values where available.

---

## Non-goals (be explicit in docs)

- Stiff systems. PFS is explicit; A-stability is not the concern here. Direct users to Rosenbrock / IMEX from `DifferentialEquations.jl`.
- DAEs.
- Essential singularities (genuinely defeats the method).
- Real-time / embedded.
- Drop-in replacement for SciML's general ODE infrastructure. **USP is correctness through pole fields and arb-prec, not breadth.**

---

## Process / philosophy constraints (per project CLAUDE.md)

- **Beads** for issue tracking. No GitHub Issues.
- **Fail-fast / fail-loud.** No silent fallbacks; numerical breakdown errors explicitly.
- **Ground truth via local PDFs.** FW papers, GGT 2013, etc. in `references/`; string-match against them rather than paraphrasing from memory.
- **No parallel agents** (Julia precompilation cache conflicts).
- **Rigorous reviewer agent** after core changes (Padé routine, stepper, AD layer).
- **~200 LOC per file.** If a module exceeds this, split.

---

## Open questions (provisional)

1. Hermite–Padé in v1 or v2?
2. Step control default: coefficient-decay or Padé-root distance?
3. Verified enclosures via Arb balls — v1, v2, or never?
4. Path-network consistency metric — value, residue, or both, to what tolerance?
5. Naming — `PadeTaylor.jl` sensible, or `MeromorphicIVP.jl` / `Pole field solver.jl` / something else?
6. Licence (MIT/BSD for ecosystem fit, GPL never).

---

## Definition of done (v1)

- [ ] `RESEARCH.md` written and reviewed.
- [ ] `DESIGN.md` written and reviewed.
- [ ] Both impls pass the shared test corpus to declared precision.
- [ ] Both impls reproduce FW 2014 PII pole-field figures qualitatively.
- [ ] Cross-impl diff ≤ 1 ULP at Float64, ≤ N digits at arb-prec for chosen N.
- [ ] Per-file LOC budget respected.
- [ ] Documentation: paper-style spec doc, API reference, minimal tutorial.

---

## Immediate next action for the follow-up agent

**Do not skip Stage 0.** Begin by producing `RESEARCH.md` per the checklist above.
Confirm or refute every provisional claim in this sketch against the literature
and the ecosystem. Return with findings and a recommendation before any code
or detailed design.

---

## v0.2 north star — vector / multi-component Painlevé-type systems

> Added 2026-05-19, post-v0.1.0 ship. Resolves v0.1 Open Question #1
> ("Hermite–Padé in v1 or v2?") with **v2**. Grounded in a survey of
> the released-software landscape conducted 2026-05-19 (four parallel
> research agents; sources listed at the end of this section).

### The opportunity, in one paragraph

A May 2026 survey of the numerical-software landscape for multi-component
Painlevé-type systems — Noumi–Yamada A_n⁽¹⁾ hierarchies, Garnier systems
G_n, Schlesinger isomonodromic deformations, matrix Painlevé equations,
Calogero–Moser pole dynamics, and the Painlevé hierarchies PJ⁽ⁿ⁾ —
returned essentially no released code. GitHub's `painleve` topic page is
empty. DLMF §32.17 cites only paper-bound code for the classical scalar
six. Published per-paper scripts (Grava–Klein 2011, Kapaev–Klein–Grava
2015, Claeys–Olver 2011) were never packaged. The tau-function CFT route
(Gavrylenko–Lisovyy) is the most powerful theoretical alternative for
PVI and Garnier but ships no general code. The Riemann–Hilbert school
(Olver's `RiemannHilbert.jl`) is general-purpose in principle but only
ever shipped a PII example, and is dormant since 2019.

PadeTaylor.jl's core methodology — Taylor → robust Padé → step control,
plus the path-network and edge-gated solve infrastructure — extends to
vector y with one nontrivial design call. Lifting it produces the first
general-purpose released numerical software for any of these systems
and, concretely, **the first complex-plane pole-field maps of any
higher-Painlevé-hierarchy member or any Noumi–Yamada A_n⁽¹⁾ solution
beyond n = 2, 3**.

### Targets, in priority order

**Primary: Noumi–Yamada A_n⁽¹⁾ systems.** Vector `y` of size `n + 1`;
single scalar `z`; Painlevé property by affine-Weyl symmetry. n = 2
reduces to PIV and n = 3 to PV (self-validation against the v0.1 scalar
implementations); n ≥ 4 is novel numerical territory — the pole-field
pictures do not exist in the literature. The path-network, edge-gated
solve, and pole-field extraction transport unchanged because `z` is
still single-complex.

**Secondary: Painlevé hierarchies (PJ⁽²⁾, PJ⁽³⁾).** Scalar high-order
ODEs of order 4, 6, 8 reducing to vector first-order systems. The PI
hierarchy is the most studied; Kapaev–Klein–Grava 2015
(`arXiv:1306.6161`) produced *one* paper of PI⁽²⁾ tritronquée
pole-field plots via ad-hoc methods without releasing code. A reusable
implementation extends every v0.1 capability (named-transcendent
constructors, closed-form-family verifiers, sheet tracking where
relevant) one rung up.

**Tertiary: Garnier G_n with one frozen time.** Vector ODE of dimension
2n in scalar z, monodromy-preserving. Direct heir of v0.1 PVI; the
multivaluedness machinery (`SheetTracker`, `BranchTracker`, sheet-aware
Stage-2) is reusable.

**Out of scope for v0.2:**

- Multi-time evolution (genuine PDE in two or more z variables). Different
  project; the package's whole substrate is single-z.
- Matrix Painlevé in the strict noncommutative-unknown sense. Plausible
  as v0.3 if v0.2 lands cleanly.
- Symbolic Bäcklund-transformation infrastructure. PadeTaylor is a
  numerical solver; the symbolic side belongs in a separate package.
- Tau-function evaluation. Different methodology (CFT series); should
  live in a separate package or an extension if pursued.

### The single load-bearing design call

PRD v0.1 Open Question #1 — *scalar per-component Padé vs Hermite–Padé
with shared denominator* — was deferred to v0.2 and is the v0.2 design
crux. Both options are viable; the trade-off:

- **Per-component Padé.** k separate `(P_i, Q_i)` per step. Trivial port:
  call `robust_pade` once per component over a
  `TaylorSeries.jl::Taylor1{Vector{T}}` jet. Downside: each component
  approximates the same physical pole with its own denominator, so pole
  locations from `extract_poles` fragment k-fold. Step control needs a
  vector norm choice.

- **Simultaneous / Hermite–Padé with shared denominator.** Single `Q`,
  k numerators `P_i`. Right answer for systems whose components share a
  singularity structure — generic for analytic systems, since a
  singularity of the flow appears in every coupled component at once.
  Requires a block-Toeplitz / stacked-SVD rewrite of `RobustPade` — a
  real port, not a wrapper. `PoleField.extract_poles` stays clean:
  roots of the single `Q` are the system's poles.

**Recommendation (provisional, to be confirmed in Stage 0+):** shared
denominator. The whole package premise is that the pole field *is* the
answer; fragmenting it across components defeats the point.

### Stage 0+ — research before v0.2 design

Mirrors v0.1 Stage 0 in shape, narrower in scope. Output: an extension
of `RESEARCH.md` covering:

- **Hermite–Padé via SVD.** Baker & Graves-Morris ch. 5; Beckermann–Labahn
  on simultaneous Padé tables; locate the correct generalisation of
  GGT 2013 Algorithm 2 to the multi-numerator case — do not re-derive.
- **Vector Taylor jets** for ODE systems. `TaylorIntegration.jl::jetcoeffs!`
  is the canonical Julia pattern; benchmark
  `TaylorSeries.jl::Taylor1{Vector{T}}` vs `Vector{Taylor1{T}}` for
  arithmetic ergonomics and type stability.
- **Noumi–Yamada normal forms.** Noumi, *Painlevé Equations Through
  Symmetry*, AMS 2004, §4; Mazzocco surveys. Acquire the canonical
  Hamiltonian forms, the affine-Weyl Bäcklund generators (for self-test:
  W(A_n⁽¹⁾) action permutes computed solutions), and the rational
  closed-form solutions (for a closed-form-family verifier analogous
  to `pii_rational`).
- **PI hierarchy explicit forms.** Cosgrove's USyd classification report
  (2000-6); the Grava–Klein and Kapaev–Klein–Grava papers for the
  pole-field reference picture.

### Stage 1 — v0.2 design decisions to lock

In addition to the shared-vs-per-component Padé crux:

1. **Vector y storage.** `Vector{T}` per node vs `Matrix{T}` of stacked
   jet coefficients vs `Vector{Taylor1{T}}`. Type stability matters.
2. **Step control for vector systems.** Jorba–Zou with `‖c_k‖_∞`,
   `‖c_k‖_2`, or component-wise minimum? Same question for the FFW 2017
   adaptive controller.
3. **`Problems` API.** New `PadeTaylorSystem` type, or generalise
   `PadeTaylorProblem` over scalar vs vector `y`? Backward-compat
   constraint: every v0.1 example must compile and produce
   bit-identical output.
4. **`PathNetwork` Stage-2 interpolation** for vector `y`. Barycentric
   per component vs vector barycentric — same shape, but confirm.
5. **`PoleField` semantics.** Shared denominator → one pole list per
   node; per-component → k lists. The design call above settles this.
6. **`Painleve` layer extensions.** `NoumiYamadaProblem(n; …)`,
   `PainleveHierarchyProblem(family, n; …)`, `GarnierProblem(n; …)` —
   each with named-solution constructors where the literature gives them.

### Test corpus for v0.2

- **Noumi–Yamada n = 2, 3 self-validation.** A_2⁽¹⁾ solutions must
  match PIV from v0.1; A_3⁽¹⁾ must match PV. This is the mutation-proof
  equivalent for the vector lift — disagreement with v0.1 on equations
  both handle is a defect.
- **Noumi–Yamada rational solutions** (Noumi 2004 §4): exact closed-form
  oracle, analogous to `pii_rational`.
- **A_4⁽¹⁾ and A_5⁽¹⁾ pole-field maps:** new mathematics-as-pictures.
  No prior published reference; acceptance is *qualitative* (visual
  coherence + W(A_n⁽¹⁾) symmetry) plus *quantitative* agreement with
  closed-form rational solutions where they exist.
- **PI⁽²⁾ tritronquée pole-field map:** cross-check against
  Kapaev–Klein–Grava 2015 (`arXiv:1306.6161`) pole locations to the
  precision they reported.
- **Calogero–Moser N-particle short-time integration:** vector ODE
  with movable poles in N coordinates; cross-validate against KdV-soliton
  closed forms (Krichever).

### Acceptance criteria for v0.2

- [ ] Shared-denominator `robust_pade` extension passing GGT-style
      mutation-proof tests against an independent Hermite–Padé oracle.
- [ ] Vector-y `Coefficients` + `PadeStepper` + `Problems` cycle GREEN
      on a Noumi–Yamada A_2⁽¹⁾ smoke test, with bit-identical-output
      on the v0.1 PIV self-validation.
- [ ] At least one Noumi–Yamada A_n⁽¹⁾ pole-field figure for n ≥ 4,
      reproducible from the `figures/` project.
- [ ] At least one PI⁽²⁾ tritronquée pole-field figure cross-checked
      against Kapaev–Klein–Grava 2015.
- [ ] All v0.1 tests still GREEN; backward compatibility is non-negotiable.
- [ ] ADRs covering: shared-denominator design, vector-y storage,
      new problem types, vector step-control. The v0.1 ADR ladder
      (0001–0018) sets the bar.

### Why this matters — the survey-derived north star

From the May 2026 survey: **there is no general-purpose public
numerical software for any of these systems.** Not in Julia, Python,
Mathematica, Maple, Sage, or any commercial tool. The unreleased
Fornberg–Weideman MATLAB code is the only prior numerical art for the
scalar Painlevé pole-field methodology, and even it never treated the
higher hierarchies or multi-component systems. Lifting PadeTaylor to
v0.2 is not a port of well-known software — there is no software to
port. It is the construction of the first such tool, and it would
ship the first complex-plane pole-field maps of solutions whose
geometry is currently unknown.

### Survey sources

- GitHub `painleve` topic — empty as of survey date.
- DLMF §32.17 — no higher-hierarchy or multi-component entries.
- Grava & Klein, *KdV hierarchy and a Painlevé transcendent*, `arXiv:1101.2602`.
- Kapaev, Klein & Grava, *On the tritronquée solutions of P_I^2*, `arXiv:1306.6161`.
- Claeys & Olver, *Numerical study of higher Tracy–Widom analogues*, `arXiv:1111.3527`.
- Klein & Stoilov, multidomain spectral for Painlevé on unbounded domains, `arXiv:1807.04442`.
- Mazzocco, *Garnier system in two variables*, `arXiv:0704.2869`.
- Adler & Sokolov, matrix PII / PIV, `arXiv:2012.05639`, `arXiv:2107.11680`.
- Bertola & Cafasso, noncommutative PII, Comm. Math. Phys. 2011.
- Gavrylenko & Lisovyy, Fredholm/Nekrasov tau, Comm. Math. Phys. 2018; `arXiv:1712.08546`.
- Korotkin, Schlesinger via theta functions, `math-ph/9810007`.
- Cosgrove, higher-order Painlevé classification, USyd report 2000-6.
- Olver / Trogdon, `RiemannHilbert.jl` (`github.com/JuliaApproximation/RiemannHilbert.jl`)
  and `RHPackage` (`github.com/dlfivefifty/RHPackage`) — PII only, no
  pole-field maps; dormant since 2019.
- Trogdon & Olver, *Riemann–Hilbert Problems, Their Numerical Solution
  and the Computation of Nonlinear Special Functions*, SIAM OT146, 2016.
- Bornemann, *On the numerical evaluation of Fredholm determinants*,
  `arXiv:0804.2543` — Tracy–Widom on real line, not a Painlevé pole-field tool.

### PRD corrections (2026-05-19, worklog 052 + 053)

Citation-verification during the Stage 0+ scoping pass (one Opus
scout + acquisition fan-out; see `docs/v0p2_literature_scope.md` and
`docs/v0p2_prior_art_sweep.md`) surfaced three misattributions in
the survey sources above. Originals are preserved verbatim to keep
the survey-provenance record honest; corrections are noted here:

- Survey line "Mazzocco, *Garnier system in two variables*,
  `arXiv:0704.2869`" → that arXiv ID is **Sasano**, *Studies on the
  Garnier system in two variables* (2007). The Mazzocco Garnier
  paper does exist and is `arXiv:math/0106208` (*The geometry of the
  classical solutions of the Garnier systems*, IMRN 2002). Both
  Sasano 2007 and Mazzocco 2002 are now in `references/garnier/`.
- Survey line "Adler & Sokolov, matrix PII / PIV, `arXiv:2012.05639`,
  `arXiv:2107.11680`" → `arXiv:2107.11680` is by **Bobrova &
  Sokolov**, *On matrix Painlevé-4 equations. Part 1: Painlevé–
  Kovalevskaya test* (2021). `arXiv:2012.05639` is correctly
  attributed (Adler & Sokolov, matrix PII).
- Survey line "Gavrylenko & Lisovyy, Fredholm/Nekrasov tau, Comm.
  Math. Phys. 2018; `arXiv:1712.08546`" → that arXiv ID is the
  **Cafasso–Gavrylenko–Lisovyy** *Tau functions as Widom constants*
  paper (CMP 365, 2019). The correct Gavrylenko–Lisovyy CMP 363
  (2018) ID is `arXiv:1608.00958`. Both are in
  `references/rh_numerical/`.

### Post-survey additions (2026-05-19, worklog 052 + 053)

The Stage 0+ scoping pass surfaced two papers the original four-agent
survey missed; appended here for completeness:

- **Mano & Tsuda (2017)**, *Hermite–Padé approximation, isomonodromic
  deformation and hypergeometric integral*, Math. Z. 285, 397–431;
  `arXiv:1502.06695`. The single most load-bearing v0.2 paper for
  Pillar A (Hermite–Padé): connects shared-denominator Hermite–Padé
  to Schlesinger transformations and the P_VI / Garnier systems via
  block-Toeplitz determinants. Now in `references/hermite_pade/`.
- **Novokshenov (2014)**, *Distributions of Poles to Painlevé
  Transcendents via Padé Approximations*, Constr. Approx. 39, 85–99;
  DOI `10.1007/s00365-013-9190-6`. The most directly methodology-
  adjacent paper to PadeTaylor's own scalar algorithm: Fair–Luke
  Padé applied to PI/PII/PIV pole distributions. v0.1 RESEARCH.md
  did not cite it; should be reproduced as a mutation-proof of v0.1's
  PIV implementation before the v0.2 vector lift. Springer paywall;
  metadata stub at `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md`.

The second-wave prior-art sweep (worklog 053; five gaps —
Cuyt / Van Iseghem / Aptekarev–Stahl / Klein-Bothner-Lisovyy-Joshi
2024–26 / GitHub code search) further surfaced:

- **Adler & Sokolov (2025)**, *Vector systems of Painlevé type*,
  J. Geom. Phys.; `arXiv:2512.18828`. 2025 classification of vector
  Painlevé systems via group reduction of vector NLS / mKdV / KdV,
  yielding ODE systems with isomonodromic Lax representations
  generalising P_I, P_II, P_34, P_IV. Pure classification, no
  numerics, no code — **independently reinforces the v0.2
  motivation by confirming the target class is alive and
  unaddressed.** Now in `references/rh_numerical/`.
- **Amore (2021)**, *Padé–Taylor method for the van der Pol
  oscillator*, Physica D; `arXiv:2111.12198`. Parallel scalar use
  of the **exact name "Padé–Taylor method"** (PTM). Construction
  is Taylor jet + diagonal Padé + residual step control on the
  scalar van der Pol equation; no FW / GGT / Painlevé framing, no
  shared-Q. Acknowledged in v0.1 README.md §Provenance as parallel
  scalar prior art for the name. Now at top level
  `references/Amore2021_pade_taylor_van_der_pol_PhysicaD_arXiv2111.12198.pdf`.

Convergence-theorem references for the v0.2 shared-Q ADR (Pillar A
analytic underpinning) are tracked in the prior-art sweep report:
Aptekarev–Stahl 1992 (paywalled chapter, metadata-only),
López Lagomasino 2019 (`arXiv:1910.08548`, open-access entry),
Fidalgo–López Lagomasino–Medina 2013 (`arXiv:1310.7010`),
López Lagomasino–Medina 2014 (`arXiv:1406.3737`),
Ikonomov–Suetin 2026 (`arXiv:2605.14760`). All in
`references/hermite_pade/`. Van Iseghem 1987 (algebraic foundation
for shared-Q with denominator-recurrence of order `d + 1`) and
Cabay–Jones–Labahn 1997 (Calgo 766; the FORTRAN shared-Q oracle
v0.2 can mutation-prove against, analogous to v0.1's use of
`padeapprox.m`) are metadata-only stubs in the same cluster.

### Out-of-scope, for the avoidance of doubt

- v0.2 is **single-z, vector-y**. Multi-z (Garnier as PDE, Schlesinger
  as multi-time evolution) is a different project.
- v0.2 is **methodology lift, not performance work**. Performance and
  registration come in v0.3 once the equations are solved correctly.
- v0.2 does **not** pursue tau-function evaluation, quantum / κ-deformed
  Painlevé, or rigorous Arb-ball enclosures (Arb enclosures remain the
  v0.1 Open Question #3, independent of the vector lift).
