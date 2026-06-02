# Worklog 070 — steepest-descent wedge ray: circular S¹ distance (C6)

**Date**: 2026-06-02
**Author**: Claude Opus 4.8 (1M) — orchestrated (Opus coding subagent).
**Bead**: `padetaylor-5jvd` (sweep C6, P3 bug, dormant).
**Scope**: fix the `:steepest_descent` wedge-edge selection in
`src/PathNetwork.jl` to use a circular (S¹) angular distance instead of a
linear one.

> **Result**: `_select_candidate`'s `:steepest_descent` branch now selects the
> wedge edge by the true geodesic gap `|rem(θ_sd − off, 2π, RoundNearest)|`.
> New unit test **PN.3.2** pins it against an independent brute-force S¹ oracle
> (two cut-straddling witnesses), **PN.3.3** guards `:min_u` parity.
> `pathnetwork_test.jl` 120/120 GREEN. Mutation-proven.

## Root cause (sweep C6)

`src/PathNetwork.jl:997` selected the wedge edge with `argmin(abs(θ_sd − off))`
— a LINEAR distance on the real line. FW 2011 §5.4.1 (`FW2011_..._JCP230.md:368`)
"the edge of the wedge closest to this steepest-descent direction" is
unambiguously a CIRCULAR S¹ distance. Near `goal_dir ≈ ±π`, `θ_sd` (from
`angle`, range (−π, π]) and an offset `off = goal_dir + θ` straddle the ±π
branch cut: e.g. θ_sd=3, off=−3 are 0.283 apart on S¹ but 6.0 on the real line,
so the linear metric crowns the wrong edge.

**Dormant**: default `step_selection = :min_u` is angle-free; no shipped figure
uses `:steepest_descent`; PN.3.1 only walked right-half-plane targets (goal_dir
near 0), never crossing the cut.

## Fix

`src/PathNetwork.jl` `_select_candidate`:
```julia
twoπ = 2 * T(π)                                   # BigFloat-safe
return argmin(abs(rem(θ_sd - off, twoπ, RoundNearest)) for off in offsets)
```
- `RoundNearest` folds the difference into (−π, π] (the geodesic gap); plain
  `mod2pi` would reintroduce a 0/2π discontinuity.
- The remainder is taken ONLY on the difference `θ_sd − off`. The `offsets`
  generator carries the index→physical-ray correspondence; wrapping each `off`
  individually would break the `argmin` index → `wedge_angles` invariant
  (`:996`). PN.3.2's second witness proves this concretely.

## Test (TDD)

- **PN.3.2** unit-tests `_select_candidate` directly. Witness 1: `goal_dir=3.0`,
  `θ_sd=−3.0`, default wedge `[−π/4,−π/8,0,π/8,π/4]` — linear picks idx 1,
  independent S¹ oracle `min(|d|, 2π−|d|)` picks idx 4. Witness 2:
  `goal_dir=3.13`, `θ_sd=−3.13` — the correct difference-only metric and the
  oracle pick idx 3, but the per-offset-wrap WRONG-fix picks idx 4 (justifying
  the "do not wrap offsets individually" warning). BigFloat repro included.
- **PN.3.3** `:min_u` ⇄ `:steepest_descent` parity on a right-half-plane grid
  (`Δu ≤ 1e-10`) — the fix does not perturb the default path; PN.3.1 unchanged.
- **Mutation-proven**: revert to `abs(θ_sd − off)` ⇒ PN.3.2 RED (idx 1 ≠ 4);
  restored. Wrong-fix (per-offset wrap) caught by witness 2. Procedure recorded
  in the PN.3.2 comment.

## Docs lockstep

`_select_candidate` comment expanded (FW md:368, S¹ rationale); ADR-0004
amendment (2026-06-02) added.

## References

- `references/markdown/FW2011_painleve_methodology_JCP230/...md:368`.
- `src/PathNetwork.jl`; `docs/adr/0004-path-network-architecture.md`;
  worklog 065 (C6).
