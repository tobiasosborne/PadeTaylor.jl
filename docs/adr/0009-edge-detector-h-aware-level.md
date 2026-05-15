# ADR-0009 — h-aware default level for `EdgeDetector.pole_field_mask`

**Status**: Accepted (2026-05-15) | **Bead**: `padetaylor-f8l`
**Context**: ADR-0004 (path-network architecture) introduced the
EdgeDetector primitive as the gate for `EdgeGatedSolve.edge_gated_pole_field_solve`
and `LatticeDispatcher.lattice_dispatch_solve`.  Both callers expose an
`edge_level` kwarg with a fixed `0.001` default copied from FW 2011
§3.2.2 Fig. 3.3.  The fix proven in worklog 032 is recorded here.

## Context

FW 2011 §3.2.2's discrete-Laplacian edge detector
(`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:202-208`)
thresholds `log₁₀|Δu|` against a "contour level" labelled `0.001` in
Fig. 3.3.  Two observations make a fixed default brittle for a
*general-purpose* solver:

  1. **FW used a single lattice spacing.**  Fig. 3.3 shows a 161×161
     grid over `[-10, 10]²`, implicit `h_grid = 0.125`.  The paper
     does not state how the level should scale with `h`; it just says
     "it is easy to select a contour level (here `0.001`) that gives a
     suitable pole field edge description" (md:208).

  2. **The residual scales as `h²`.**  For a complex-analytic `u`, the
     5-point Laplacian residual is purely the truncation-error term
     `(h²/12)(∂⁴u/∂x⁴ + ∂⁴u/∂y⁴) = h²·u⁗(z)/6` (the analytic
     Laplacian vanishes by Cauchy–Riemann).  Both the smooth-region
     floor AND the pole-annulus residual that the
     [`EdgeGatedSolve`](@ref) flood-fill relies on for cell-to-cell
     connectivity therefore shift by `2·log₁₀(h)` in log space as the
     grid refines.

The combination bites at fine grids: when `h_grid` drops below ≈ 0.25,
the connecting annulus residual falls below `10^0.001 ≈ 1`, the flood-
fill loses connectivity, and `edge_gated_pole_field_solve` stalls at
the seed.  The PI tritronquée at `[-20, 20]²` is the empirical
testbed: at `h = 0.25` (N = 161) the default `level = 0.001` works
and the field grows for 24 passes covering ≈ 18 % of the grid;
at `h = 0.125` (N = 321) it stalls (2 passes, seed only); at
`h = 0.0625` (N = 641) it stalls.  The working levels at the three
spacings are `0.001`, `≈ -0.6`, `≈ -1.0` respectively — matching
`LEVEL0 + 2·log₁₀(h/H0)` with `(H0, LEVEL0) = (0.25, 0.001)` to two
decimals.

The bug is in our *generalisation* of FW's recipe to arbitrary `h`,
not in FW's recipe itself.  We exposed `h` as a free kwarg but kept
the level constant — CLAUDE.md Rule 2 "all bugs are deep", not "pass
edge_level manually."

## Decision

`pole_field_mask`'s default level becomes a sentinel `:auto` that
resolves to

```
level(h) = LEVEL0 + 2·log₁₀(min(h, H0) / H0)
```

with the empirically-calibrated anchor `(H0, LEVEL0) = (0.25, 0.001)`,
clamped at `H0` for `h ≥ H0`.

The clamp is the load-bearing detail.  For `h ≥ H0` the threshold
stays at FW's published `0.001`, which preserves backward compatibility
for the existing FW-figure test grids at `h_grid ∈ {0.333, 0.5, 0.75,
1.0}` — the smooth/annulus gap at those spacings is wide enough that
`0.001` sits comfortably in it, and at the coarsest existing spacing
(`h = 0.75` in `test/edge_gated_solve_test.jl`) the *unclamped* formula
would raise the threshold past the actual pole-annulus residual
(`|Δu| ∝ 1/h³` at near-pole cells, `log₁₀ ≈ 0.38` at `h = 0.75`),
collapsing the mask.  Mutation M3 in worklog 032 demonstrates this
cascade: dropping the `min(h, H0)` clamp makes the unit-test ED.5.1's
clamp assertions fail AND the existing edge-gated-solve tests EG.1.1 +
EG.1.2 fail.

Both wrapper modules thread the new sentinel:

  - `EdgeGatedSolve.edge_gated_pole_field_solve`: `edge_level::Union{Real,Symbol} = :auto`.
  - `LatticeDispatcher.lattice_dispatch_solve`: same; the existing
    `isfinite(edge_level)` fail-fast guard now also accepts `:auto`.

A numeric `level` (or `edge_level`) override still works for callers
who want to reproduce FW's verbatim recipe — `figures/fw2011_fig_3_3.jl`
passes `level = 0.001` explicitly for exactly this reason.

The 1-arg form `pole_field_mask(Δu; level)` has no `h` in scope, so
`:auto` is unresolvable and throws `ArgumentError` per CLAUDE.md
Rule 1.  Pass an explicit numeric `level` or use the 2-arg form.

## Consequences

**Positive**:

  - `edge_gated_pole_field_solve` now grows correctly at arbitrary
    grid spacing.  `examples/tritronquee_3d.jl` ships at N = 641
    (`h_grid = 0.0625`) — 4× the resolution of the v0.1.0 hero — with
    no manual `edge_level` override.
  - Existing FW-figure tests (`fw_fig_41_test.jl`, `phase9_tritronquee_test.jl`,
    `edge_gated_solve_test.jl`, all at `h_grid ≥ H0`) are unaffected
    because the clamp keeps the default at `0.001` in their regime.
  - The new ED.5.1 / ED.5.2 / ED.5.3 testsets (22 assertions) make the
    h-scaling structurally verified; four mutation cycles bite at
    least one of them apiece.

**Negative**:

  - The anchor `(0.25, 0.001)` is empirically pinned from the PI
    tritronquée at `[-20, 20]²`, not derived from FW 2011.  FW's own
    calibration anchor (read off Fig. 3.3) would be `(0.125, 0.001)`
    — at the same `h = 0.125`, our formula gives `level ≈ -0.6`, FW's
    paper says `0.001` works.  The discrepancy is the problem-
    dependent constant in the truncation-error proportionality
    `|Δu| ∝ h²·|u⁗(z)|` — the magnitude of `|u⁗(z)|` in smooth regions
    is solution-dependent, so a single universal calibration anchor
    is not possible.  We pick the anchor that makes the empirical
    PI-tritronquée sweep work; FW's `0.001 at h = 0.125` calibration
    is still recoverable by passing `level = 0.001` explicitly.

  - Calibration is for *PI* (and structurally similar Painlevé
    transcendents).  Problems with very different `|u⁗(z)|`
    distributions might need their own anchor.  Documented in the
    `EdgeDetector` module docstring as the v1 calibration; callers
    can override with an explicit `level`.

## Alternatives considered

**A. Leave the default fixed at 0.001 and document that fine grids
need a manual override.**  Rejected per CLAUDE.md Rule 2 — pushing
the calibration onto every caller is a band-aid, not a fix.  The
existing module docstring already half-acknowledged the issue
(*"calibrated for fine grids"*, *"Grid resolution matters"*); the
honest fix is to make the default track `h`.

**B. Use FW's own calibration `(H0, LEVEL0) = (0.125, 0.001)`.**
Rejected — at the empirical PI testbed, `level = 0.001` at
`h_grid = 0.125` does not work (the bead's three-point sweep
documents this directly).  FW's anchor is for FW's specific
near-tritronquée; ours is for the actual Boutroux tritronquée we
ship as `tritronquee(:I)`.  Both are "correct" in their own setup;
we pick the one that matches our reference solution.

**C. Don't clamp; use the unclamped formula `LEVEL0 + 2·log₁₀(h/H0)`
in both directions.**  Rejected — mutation M3 in worklog 032
demonstrates this cascades into RED on existing tests at `h_grid ∈
{0.5, 0.75}`.  The math supports clamping: in the coarse regime the
smooth/annulus gap is wider than at fine grids, so a lower threshold
suffices; raising the threshold by `2·log₁₀(h/H0)` for `h > H0`
overshoots the actual pole-annulus residual.

**D. Auto-detect the calibration anchor from the input `u_grid`.**
Rejected for v1 — requires histogram analysis of `|Δu|` to find a
bimodal separation, which couples the detector to the solve in ways
the FW reference does not.  Possible future work; would be a
different ADR.
