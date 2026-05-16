# ADR-0015 — `extrapolate=true` Stage-2 policy (FFW 2017 §2.1.1 alignment)

**Status**: Accepted (2026-05-16) | **Bead**: `padetaylor-afs` | **Worklog**: 045

## Context

Our `PathNetwork.path_network_solve`'s Stage-2 evaluation (worklog
006 / ADR-0004) refuses to evaluate the stored local Padé
approximant at a fine-grid cell when the cell sits more than
`visited_h` from the nearest visited node:

```julia
# src/PathNetwork.jl:596-602
if abs(z_f - z_v) > h_v
    grid_u[i]  = nan_CT       # fail-soft, no silent extrapolation
else
    t = (z_f - z_v) / h_v
    grid_u[i]  = _evaluate_pade(visited_pade[idx_v], t)
end
```

The motivation, per ADR-0004 §"Stage 2" and the PathNetwork.jl
top-of-file docstring (lines 30-32):

> "For each `z_f` in `grid` find the nearest visited `z_v`; if
> `|z_f - z_v| ≤ h`, evaluate `visited_pade[k]` at `t = (z_f - z_v)
> / h`.  Otherwise the slot gets `NaN + NaN·im` — **no silent
> extrapolation**."

This honors CLAUDE.md Rule 1 (fail fast, fail loud): NaN tells the
caller "this cell falls outside the disc of validity of any nearby
local Padé."  Downstream consumers like `PoleField.extract_poles`
and the smoke-test diagnostics rely on the NaN signal to spot
under-resolved regions.

## The divergence from FFW 2017

FFW md:62 specifies Stage 2 differently:

> "Stage 2: The solution is computed at all points on a fine grid.
> This is accomplished as follows. Let ζ_i, 1 ≤ i ≤ N denote the
> points where Padé coefficients are available from Stage 1.  Padé
> steps are taken from each ζ_i to the points on the fine grid to
> which it is the closest point among the ζ_i."

And explicitly md:95 for Figure 1's Stage-2 fill:

> "the Padé coefficients at the 2701 points in the second column
> were used to compute Padé approximations to the solution w at the
> 571×140 = 79940 points of the fine grid."

FFW DO NOT check disc validity — every fine-grid cell gets a Padé
value computed from its nearest visited node, even when the offset
`|z_f - z_i|` exceeds the canonical Padé step magnitude.  The
result: every cell rendered, no white gaps.  Errors past `|t| = 1`
degrade gracefully (the rational function is still well-defined
outside the disc; it just doesn't enjoy the FW 2011 §2.1.4 error
bounds there).

Our every-figure-has-gaps render artefact (~25% of pixels NaN at
the rectangular target grids we ship; 0.9% pathological case in
FFW Fig 2 col 3b alt-sheet panel) is the visible price of choosing
fail-soft over FFW-style extrapolation.

## Decision

Add an opt-in `extrapolate::Bool = false` kwarg to:

  1. `PathNetwork.path_network_solve` — gates the Stage-2 disc-radius
     check.  Default `false` preserves the existing ADR-0004 contract.

  2. `PathNetwork.eval_at_sheet` — the A5 per-point sheet-aware
     accessor; same gate (the disc check is at line 47 of the
     accessor's body).

  3. A new public `PathNetwork.eval_at(sol, z; extrapolate = false)`
     — sheet-blind per-point accessor (the sheet-blind analogue of
     `eval_at_sheet`).  Useful for figure scripts that currently
     roll their own `stage2_eval_blind` helper with the same disc
     check baked in.

When `extrapolate = true`, the disc-radius check is skipped:
every cell / query gets `_evaluate_pade(visited_pade[idx_v], t)`
regardless of how far `t = (z_f - z_v) / h_v` is from 1.

## Why opt-in default-false (vs flipping the default)

  - **2417 existing GREEN assertions** (as of worklog 044) lock in
    the NaN-on-extrapolation behaviour.  Tests like the PathNetwork
    `:fixed` step-size pin assert specific NaN cells; flipping the
    default would either silently regress them or require a
    one-by-one audit.

  - **Downstream consumers** rely on NaN for fail-soft signaling:
    `PoleField.extract_poles` skips NaN nodes implicitly via
    `isfinite(...)` guards; the smoke-test diagnostics (`coverage =
    count(isfinite, ...) / length(...)`) become uninformative when
    every cell is "finite" by extrapolation.

  - **CLAUDE.md Rule 1** is a default-on contract: fail-loud
    behaviour requires an explicit opt-out, not an implicit one.
    Figure scripts choosing extrapolation make the choice visible
    in the call site.

  - **Equivalent expressiveness**: callers wanting FFW-style fills
    pass `extrapolate = true` once at the call site.  No new
    information lost.

## Why a single shared kwarg vs per-method dispatch

The choice between "fail-soft" and "extrapolate" is the SAME
choice across `path_network_solve`, `eval_at_sheet`, and `eval_at`
— it's the Stage-2 disc-radius policy.  A single shared kwarg
name `extrapolate` lets callers think about it as one decision.

(An alternative design would use a `stage2_policy::Symbol` enum
`:fail_soft / :extrapolate` for forward-compatibility with more
policies — e.g. `:linear_extrapolate`, `:bilinear_interpolate`.
Rejected for v1: YAGNI.  Symbols can replace the boolean later
without a breaking change since the bool is opt-in.)

## Consequences

### Test coverage

A new test file `test/extrapolate_test.jl` covers:
  - `extrapolate=false` (default): cell beyond `visited_h` → NaN.
  - `extrapolate=true`: same cell → finite (the Padé extrapolated
    past `|t|=1`).
  - Cell INSIDE `visited_h`: both modes return identical value.
  - `eval_at_sheet` and `eval_at` honour the kwarg.

Target: ~15 GREEN assertions; one mutation (drop the
`extrapolate ||` guard at each of three call sites) bites.

### Documentation

  - This ADR-0015.
  - Worklog 045.
  - PathNetwork.jl Stage-2 docstring updated.
  - eval_at_sheet docstring updated.
  - New eval_at docstring written.
  - Figure-script header docstrings note their use of
    `extrapolate = true`.

### Figure-script changes (the v1 motivation)

Update existing FFW figure scripts to opt into extrapolation:

  - `figures/ffw2017_fig_1.jl` — replace local `stage2_eval` with
    `eval_at(sol, z; extrapolate=true)`, OR pass `extrapolate=true`
    to `path_network_solve` and read `sol.grid_u` directly.
  - `figures/ffw2017_fig_2.jl` — same for the three solves.
  - `figures/ffw2017_fig_4.jl`, `figures/ffw2017_fig_5.jl`,
    `figures/ffw2017_fig_6.jl` — same.

Re-render all 5 to confirm the white gaps close.

### Open follow-ups (v2 work, not in this bead)

  - **FW 2011 figure scripts** (Fig 3.1, 4.1, etc.) — these don't
    discuss extrapolation; FFW md:62's clarification post-dates
    them.  Leave the FW scripts on fail-soft default; they're
    architecturally fine.

  - **`bilinear_interpolate` Stage-2 policy** — Fig 1's existing
    `bilinear_w` helper interpolates between 4 pre-evaluated cells.
    A `:bilinear` `stage2_policy` would do this in the package
    proper, but it requires a regular-grid Stage-2 evaluation
    (which the package's vectorised `path_network_solve(prob, grid)`
    doesn't currently assume).  Separate bead if a caller asks.

  - **Poisson-disk Stage-1 node placement** (`padetaylor-zwh`,
    already filed) — the principled solution to the visual gaps at
    larger Re ζ.  Orthogonal to extrapolate=true.  Both can ship.

## References

  - **FFW 2017 §2.1.1** — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:62` (Stage 2 spec) + md:95 (Figure 1 concrete numerics).
  - **ADR-0004** — `docs/adr/0004-path-network-architecture.md` (the existing fail-soft policy).
  - **ADR-0013** — `docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md` (A4); the visited_sheet field that A5's eval_at_sheet queries.
  - **Worklog 044** — `docs/worklog/044-ffw2017-fig-2.md` (the FFW Fig 2 ship that surfaced the white-gap issue).
  - **`src/PathNetwork.jl`** — the file the new kwargs land in.
