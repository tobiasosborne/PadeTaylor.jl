# Changelog

All notable changes to PadeTaylor.jl are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-05-13

First numbered release.  All four architectural layers and five
composition tiers shipped; 1311 / 1311 tests passing.

### Added

**Tier 1 — the four-layer core (Phases Z, 1–6).**

  - `LinAlg` — SVD dispatcher with relative-accuracy guarantees on
    small singular values (LAPACK for `Float64`,
    `GenericLinearAlgebra` one-sided Jacobi for `BigFloat`,
    `Matrix{Arb} → BigFloat` shim via the `PadeTaylorArblibExt`
    extension).
  - `RobustPade.robust_pade` + `PadeApproximant` — GGT 2013
    Algorithm 2 with Chebfun's QR-reweighting (lines 278–280 of
    `padeapprox.m`).
  - `Coefficients.taylor_coefficients_1st` /
    `taylor_coefficients_2nd` — Taylor jet generation via
    `TaylorSeries.jl::Taylor1{T}` (validated for
    `T ∈ {Float64, BigFloat, Arblib.Arb}`).
  - `StepControl.step_jorba_zou` — Jorba–Zou 2005 §3.3.1 eq. 11
    truncation-error step formula (TaylorIntegration.jl-equivalent).
  - `StepControl.step_pade_root` — FW 2011 §3.1 forward-projection
    Padé-denominator-root distance heuristic.
  - `PadeStepper.pade_step!` and `pade_step_with_pade!` — one-step
    orchestrator from `(z, u, u')` to `(z + h, u(z + h), u'(z + h))`.
  - `Problems.PadeTaylorProblem` / `solve_pade` /
    `PadeTaylorSolution` / `taylor_eval` — public driver layer with
    dense output over multi-segment trajectories.

**Tier 1.5 — package extensions (Phases 7–8).**

  - `PadeTaylorCommonSolveExt` — `CommonSolve.jl` adapter exposing
    `PadeTaylorAlg <: CommonSolve.AbstractAlgorithm` with
    `init` / `step!` / `solve!` methods so PadeTaylor.jl participates
    in the wider SciML solver ecosystem.
  - `PadeTaylorArblibExt` — Arbitrary-precision SVD via
    `Arb → BigFloat → GenericLinearAlgebra.svd` (Arblib.jl has no
    native SVD; verified by source inspection).

**Tier 2 — path-network (Phase 10).**

  - `PathNetwork.path_network_solve` + `PathNetworkSolution` — FW 2011
    §3.1 five-direction wedge path-tree solver with Stage-2
    fine-grid barycentric extrapolation from stored Padé
    approximants.  Generic in `T <: AbstractFloat` for `Float64` and
    `Complex{T}`.  Optional Schwarz-reflection symmetry mode
    (`enforce_real_axis_symmetry::Bool` kwarg) for real-coefficient
    real-IC problems.

**Tier 2.5 — pole-field edge detection (Phase 12.5).**

  - `EdgeDetector.laplacian_residual` + `pole_field_mask` — 5-point
    Laplacian classifier per FW 2011 §3.2.2 (threshold on
    `log₁₀|Δu|`).

**Tier 3 — boundary-value composition (Phases 11, 12).**

  - `BVP.bvp_solve` + `BVPSolution` — Chebyshev spectral-collocation
    Newton solver for second-order analytic BVPs on a complex
    segment with Dirichlet BCs (FW 2011 §3.2; Trefethen *SMIM*
    chapters 6 & 13; Berrut–Trefethen barycentric evaluation).
    Step-norm Newton convergence (`eps(T)^(3/4)` default).
  - `Dispatcher.dispatch_solve` + `IVPSegment` / `BVPSegment` /
    `DispatcherSolution` — 1D IVP↔BVP chain composition layer per
    FW 2011 §4.4 with junction derivative-match diagnostics.
  - `LatticeDispatcher.lattice_dispatch_solve` + `LatticeSolution` —
    2D-grid composition with per-row BVP fill on smooth runs
    flanked by IVP cells (FW 2011 line 190).

**Tier 4 — exponential coordinate transforms (Phase 13).**

  - `CoordTransforms.pIII_transformed_rhs` /
    `pV_transformed_rhs` / `pIII_z_to_ζ` / `pIII_ζ_to_z` /
    `pV_z_to_ζ` / `pV_ζ_to_z` — Exponential coordinate maps
    (`z = exp(ζ/2)` for PIII, `z = exp(ζ)` for PV) that remove the
    fixed branch point at `z = 0`.  Helpers-only; integration via
    composition with `PadeTaylorProblem` + `path_network_solve`.

**Tier 5 — Riemann-sheet tracking (Phase 14).**

  - `SheetTracker.pVI_transformed_rhs` — ζ-plane RHS for the sixth
    Painlevé equation (FFW 2017 eq. 3).
  - `SheetTracker.winding_delta` / `accumulate_winding` /
    `sheet_index` — path-side Riemann-sheet bookkeeping primitives
    (signed angle change normalised to `(-π, π]`, cumulative
    winding, sheet-index assignment via `round(total / 2π)`).

**Documentation.**

  - Documenter.jl-built docs site at `docs/build/` (regenerable via
    `julia --project=docs docs/make.jl`).  Sections: Home,
    Architecture, API, Figures.  Local-only per CLAUDE.md Rule 11
    (no `deploydocs`, no CI).
  - Four Architecture Decision Records under `docs/adr/`:
    four-layer architecture (0001), BigFloat-SVD via
    `GenericLinearAlgebra` (0002), Pkg.jl weak-dep extensions (0003),
    path-network architecture (0004).
  - Figure-acceptance catalogue at `docs/figure_catalogue.md`
    covering 79 figures across FW 2011, FW 2014, FW 2015,
    RF 2014, FFW 2017.
  - 18 worklog shards under `docs/worklog/` recording frictions
    surfaced + mutation-proof procedures + algorithmic findings for
    each shipped phase.

### Headline empirical results

  - **FW 2011 Table 5.1 row 1** (long-range integration of the
    equianharmonic Weierstrass ℘-function to `z = 30`):
    **`2.13·10⁻¹⁴` rel-err** in `BigFloat`-256 — beats FW's
    published `8.34·10⁻¹⁴`.
  - **FW 2011 Fig 4.1 step (i)** (tritronquée BVP on `[-20i, +20i]`):
    `u(0)` pinned to `≤ 3.5·10⁻¹³`, `u'(0)` to `≤ 5.3·10⁻¹¹` vs
    FW eq. 4.1 reference values at `N = 240` Chebyshev nodes.
  - **FW 2011 Fig 3.1** (PI tritronquée pole field): qualitative
    PARTIAL reproduction at 25×25 over `[-4, 4]²` — 4-of-5 pole-free
    sectors recovered, conjugate symmetry verified, leading-pole
    magnitude matches Joshi–Kitaev to `≤ 10⁻³`.
  - **Phase-6 pole-bridge demo**: at `z = 1.05` (just past the
    Weierstrass-℘ pole at `z = 1`), the order-30 Padé conversion
    matches the closed-form solution to `3.45·10⁻¹⁰` while plain
    Taylor truncation of the same coefficients diverges to `2.5`
    (relative error) — a `9.86`-orders-of-magnitude gap on
    identical input.

### Cross-validation oracles

The test suite cross-validates against:

  - Mathematica's closed-form `WeierstrassP[z + c₁, {0, c₂}]`;
  - Mathematica's `NDSolve` at `WorkingPrecision = 50`;
  - `mpmath.odefun` at 40 decimal digits (Python);
  - Chebfun's `padeapprox.m` under Octave (Phase 2 oracle);
  - DMSUITE `chebdif` / `chebint` under Octave (Phase 11 BVP oracle);
  - `TaylorIntegration.jl::stepsize` + mpmath + wolframscript
    (Phase 4 three-source pin).

### Project discipline

  - 14 source modules, each under the 200-LOC ceiling (CLAUDE.md
    Rule 6).
  - All load-bearing tests mutation-proven (CLAUDE.md Rule 4):
    perturb impl, confirm RED, restore.
  - Ground-truth-before-code (CLAUDE.md Law 1) enforced via
    line-cited references in commit messages and ADRs to
    `references/markdown/<paper>/<file>.md`.
  - No GitHub CI (CLAUDE.md Rule 11) — quality gates run locally
    via `julia --project=. -e 'using Pkg; Pkg.test()'`.

### Known limitations (v1)

  - **Coordinate-transform tier (Phase 13)** ships RHS factories
    and IC round-trips for PIII / PV.  Non-uniform Stage-1 node
    placement and adaptive Padé `h` are independent follow-ons,
    not bundled.  Suitable for FFW 2017 Fig 1, 4, 5, 6 at one-step
    accuracy.
  - **Sheet-tracking tier (Phase 14)** ships the ζ-plane PVI RHS
    and post-walk winding-number primitives.  Constrained-wedge
    `PathNetwork` routing (that refuses to overstep branch cuts
    during the walk) and sheet-aware Stage-2 evaluation are
    deferred.  Suitable for FFW 2017 Fig 2, 3, 7 at one-step
    accuracy.
  - **`Polynomials.roots` for `Arb` element type** — validated only
    for `Float64` / `Complex{Float64}`.  Arb-precision Padé-root
    step deferred (friction bead `padetaylor-8pi`).
  - **GitHub-pages deployment** — explicitly out of scope per
    CLAUDE.md Rule 11.

### References

  - B. Fornberg & J. A. C. Weideman, *A numerical methodology for
    the Painlevé equations*, J. Comput. Phys. 230 (2011), 5957–5973.
  - P. Gonnet, S. Güttel & L. N. Trefethen, *Robust Padé
    Approximation via SVD*, SIAM Review 55 (2013), 101–117.
  - À. Jorba & M. Zou, *A software package for the numerical
    integration of ODEs by means of high-order Taylor methods*,
    Experimental Mathematics 14 (2005), 99–117.
  - M. Fasondini, B. Fornberg & J. A. C. Weideman, *Methods for
    the computation of the multivalued Painlevé transcendents on
    their Riemann surfaces*, J. Comput. Phys. 344 (2017), 36–50.

[Unreleased]: https://github.com/tobiasosborne/PadeTaylor.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tobiasosborne/PadeTaylor.jl/releases/tag/v0.1.0
