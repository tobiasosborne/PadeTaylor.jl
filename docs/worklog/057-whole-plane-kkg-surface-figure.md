# Worklog 057 — whole-plane KKG 7.4/7.5 P_I^(2) tritronquée surface figure

**Date**: 2026-05-21
**Author**: Claude Opus 4.7 (orchestrator) + serial Opus coding subagents
**Epic**: `padetaylor-0ln` (v0.2 — vector / multi-component Painlevé-type systems)
**Scope**: Orchestrated, serial build of the **headline figure for the v0.2
vector solver** — the P_I^(2) tritronquée V_0(x,0) rendered over the complex
x-plane `[-20,20]²` as Re/Im heatmaps + 3D surfaces, reproducing
Kapaev–Klein–Grava (KKG) 2015 Figs 7.4/7.5. Six phases: F1, F2, the two
Laplace solvers, F3, F4. Continues worklog 056 (which delivered the vector
BVP solver + the V8b pole-field figure).

> **Take-home**: the pole-free sector of the tritronquée is now computed
> THREE independent ways — ray-fan `vector_bvp_solve`, an in-house
> 2D-Chebyshev spectral Laplace solve, and a Gridap FEM Laplace solve — and
> **majority-voted per grid point** (user-directed verification). The three
> agree to ~5·10⁻⁴ across the bulk and the vote *exposes* the inner-arc
> disagreement (where the asymptotic boundary data degrades) rather than
> hiding it. The headline figure ships; its pole-wedge render is a
> documented v1 corner.

## Orchestration method

Parent Opus orchestrator; each phase delegated to one Opus coding subagent
(serial — one Julia process, Rule 7); read-only scoping to a Sonnet recon.
After each phase the orchestrator re-verified the new tests standalone
(Rule 3), confirmed the commit, raised friction beads, closed the bead, then
dispatched the next. Beads created before code (Rule 8).

## What shipped

| Phase | Bead | Deliverable | Commit |
|---|---|---|---|
| F1 | 0ln.32 | `pI2_tritronquee_ic` generalised to any complex point on the V_0 sheet (the asymptotic seed on rays `x = ξe^{iφ}`) | `76e1dad` |
| F2 | 0ln.20 | `VectorPathNetwork` Stage-2 fine-grid fill (`fine_grid` kwarg; split → `VectorPathNetworkStage2.jl`) | `1d3c32b` |
| Laplace (2) | 0ln.35 | `src/Laplace2D.jl` — in-house 2D-Chebyshev tensor-product spectral Laplace solver + **ADR-0024** | `2c184d4` |
| Laplace (3) | 0ln.36 | Gridap FEM Laplace voter — re-homed to the `figures/` project after a package-weak-dep mis-step | `cfe8b86` → `a15a21d` |
| F3 | 0ln.33 | `figures/_kkg_pi2_surface_helpers.jl` — whole-plane compute kernel: triple-method majority vote (sector) + path-network Stage-2 (wedge) + stitch | `7fcbe4b` |
| F4 | 0ln.28 | `figures/kkg_pi2_tritronquee_surface.jl` — the headline figure (Re/Im heatmaps + 3D surfaces) | `71d2dbc` |

Full suite re-verified GREEN at the milestone gate: **4712 / 4712** (14m34s)
after F1+F2+`Laplace2D` landed. F3/F4 are `figures/`-env-only (additive — no
`src/`, no main `test/`), validated by `figures/test_kkg_pi2_surface.jl`
(17888 assertions GREEN) + the rendered PNG.

## The triple-method majority vote (user-directed verification)

The pole-free sector (`arg x ≈ 36°…324°`) solution is computed three ways:
1. **Ray-fan `vector_bvp_solve`** — the P_I^(2) BVP on a fan of radial rays;
   the genuinely independent method (a direct ODE solve).
2. **In-house 2D-Chebyshev Laplace** — `Re u`/`Im u` are harmonic; the
   conformal map `w = log x` sends the annular sector to a rectangle (the
   Laplacian is conformally invariant), solved by a tensor-product Chebyshev
   spectral Laplacian (Trefethen SMIM — KKG's own cited method).
3. **Gridap FEM Laplace** — the same problem, an independent FEM
   discretisation.
Per Cartesian grid point: the **median** of the three (per component). The
agreement/spread map is retained as a diagnostic. Median spread ≈ 5.5·10⁻⁴
in the bulk; ≈ 2.9·10⁻³ at the inner arc — see lesson 57.

## The Gridap dependency mis-step (a recorded lesson)

The Gridap voter was first built as a package weak-dependency extension
`PadeTaylorGridapExt` (commit `cfe8b86`). This **broke the build**:
`Pkg.test()` could no longer resolve its environment — Gridap 0.18 pins
`ForwardDiff 0.10.x`, unsatisfiable against the package's `SpecialFunctions
2.7.2` stack. Corrective commit `a15a21d` backed the weak-dep out entirely
(restoring `Pkg.test()` resolvability) and re-homed Gridap into the separate
`figures/` project (Gridap 0.20.7 — which dropped the ForwardDiff pin —
resolves cleanly there alongside CairoMakie). ADR-0024 carries a Correction
section. The triple-method vote is unchanged; it executes in F3 within the
`figures/` env.

## Hard-won lessons

55. **A heavy package weak-dep can break the core test environment.** Gridap
    (~160 transitive deps) as a `[weakdeps]` entry made `Pkg.test()`
    unresolvable. Figure-only infrastructure belongs in the separate
    `figures/` project — which already absorbs heavy deps (CairoMakie)
    without touching package testability. Justify (Law 2) AND *locate* every
    dependency.
56. **`exp(iπ/3) ≠ −1`.** The orchestrator's F1 brief mis-stated the V_0
    cube-root branch (`θ ∈ [0,2π)`, `θ=π` on the negative axis). The negative
    real axis needs `θ = 3π` so `exp(iθ/3) = exp(iπ) = −1`; the V_0 sheet is
    `θ ∈ (2π,4π]`. The F1 agent caught and corrected the brief (Rule 3 —
    verify even the instructions). Derivatives are then branch-free:
    `Y′ = Y/(3x)`, so the seed components are rational in `Y` and `x`.
57. **A majority vote earns its keep by exposing disagreement.** The three
    sector solvers agree to ~5·10⁻⁴ in the bulk but diverge to ~3·10⁻³ at the
    inner arc — exactly where the `n_terms=2` asymptotic boundary datum is
    only `O(10⁻²)` accurate. The median anchors on the exact ray-fan BVP
    there. A single method would have silently carried the error; the vote
    surfaced it as a spread-map hotspot.
58. **The conformal map collapses the geometry.** `w = log x` sends the 270°
    annular sector to a plain rectangle, and the Laplacian is conformally
    invariant in 2D — so a ~150-LOC in-house rectangle Laplace solver
    (reusing the existing Chebyshev `D₁`/`D₂`) suffices; no FEM package and
    no curved-domain meshing are needed. A heavy general solver was the
    wrong tool for a rectangle-mappable problem.

## Beads — this session

**Created**: `0ln.32` (F1), `0ln.34` (Laplace scoping), `0ln.35` (2D-Cheb),
`0ln.36` (Gridap), `0ln.33` (F3). **Closed**: `0ln.32`, `0ln.20` (F2),
`0ln.35`, `0ln.34` (scoping — superseded by the user's triple-method
decision), `0ln.36`, `0ln.33`, `0ln.28` (F4).
**Open / deferred**: `0ln.19` (V9 — v0.2 docs/release prep), `0ln.21`/`0ln.22`
(V7 deferred sophistications), `0ln.23` (shared-Q wedge selector — would
sharpen the figure's pole wedge), `0ln.27` (ChebUtil de-duplication).

## v1 corners (Rule 9 — documented)

- **`pI2_tritronquee_ic`**: `n_terms=2` is `t=0` only; the `x>0` wedge branch
  is unimplemented (fail-fast).
- **Surface figure**: the pole **wedge render is coarse** — the path-network
  visited hull only spans `|x| ≲ 8` and the grid is 121²; the Stage-2 fill is
  called with `extrapolate=true` beyond the hull (flagged per cell). Forcing
  condition for a crisp wedge: bead `0ln.23` (the shared-Q wedge selector) +
  a wider path-network target fan.
- **Inner-arc boundary data** from the asymptotic series is `O(10⁻²)` at
  `|x|≈2`; the harmonic solve damps it and the ray-fan voter anchors the
  median. Forcing: a pin inside `|x|≲3` needing better than `10⁻²`.
- **Stokes-strip masking**: `±3°` of arc NaN-masked on each `±36°` sector
  edge (the asymptotic seed degrades toward the Stokes lines).
- The ray-fan→Cartesian step is bilinear (Methods 2/3 are the principled
  harmonic fills; the vote reconciles all three).

## Pickup point for the next agent

The vector BVP solver, V8b, and the whole-plane headline figure are all
done. **V9 (`padetaylor-0ln.19`) — v0.2 docs + release prep — is the last
planned epic phase**: README/CHANGELOG/RESEARCH v0.2 sections, ADR review
(ADR-0023 + ADR-0024 are new), HANDOFF refresh, the `padetaylor-0ln` epic
close-out. If a sharper headline figure is wanted before release, do
`0ln.23` (the shared-Q wedge selector) first and re-run F3/F4 with a wider
wedge target fan. The `figures/` project now carries Gridap 0.20.7 — the
surface figure + its test run under `--project=figures`.
