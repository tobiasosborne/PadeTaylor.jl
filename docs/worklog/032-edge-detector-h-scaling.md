# Worklog 032 — `EdgeDetector` h-aware level scaling

**Date**: 2026-05-15
**Author**: Claude Opus
**Bead**: `padetaylor-f8l` — closed.
**Scope**: Fix the `EdgeDetector.pole_field_mask` default-level bug
that blocked any `edge_gated_pole_field_solve` lattice finer than
`h_grid ≈ 0.25` — and ship the 4× resolution `examples/tritronquee_3d.jl`
hero render that motivated it.

> **Take-home**: The 5-point Laplacian residual of an analytic
> function scales as `|Δu| ∝ h²`; in `log₁₀` the smooth-region floor
> and the pole-annulus residual the flood-fill rides on BOTH shift by
> `2·log₁₀(h)` as the grid refines.  A fixed-`0.001` threshold (FW
> 2011's contour at the single resolution they used) over-classifies
> as smooth on fine grids, the flood-fill loses connectivity, and
> region growing stalls at the seed.  Fix: `level = :auto` resolving
> to `LEVEL0 + 2·log₁₀(min(h, H0)/H0)` with anchor `(H0, LEVEL0) =
> (0.25, 0.001)`; clamped at `H0` so coarser grids stay at FW's
> published level and existing tests don't regress.  ADR-0009
> accepted; 1630 → **1652 GREEN** (+ 22 new ED.5.* assertions);
> `examples/tritronquee_3d.png` re-renders at N = 641 / h = 0.0625.

## Ground truth read first — the reconnaissance

Four parallel Sonnet read-only subagents (CLAUDE.md Rule 7 allows
parallel reads) at session open:

  1. **The three affected modules** (`src/EdgeDetector.jl`,
     `src/EdgeGatedSolve.jl`, `src/LatticeDispatcher.jl`) — surface-area
     map: where is the threshold compared, where does `h` flow, where
     does `edge_level` propagate, what fail-fast guards need updating.
     Key finding: `h_grid` is in scope at both wrapper call sites
     (`EdgeGatedSolve:301` and `LatticeDispatcher:195`), so the
     auto-level computation can stay inside `pole_field_mask` rather
     than threaded through each caller.

  2. **FW 2011 §3.2.2 ground truth**.  The agent's verbatim quote
     hunt confirmed:
       - The 5-point Laplacian is the normalised form `(stencil·u)/h²`
         (`FW2011_*.md:204-206` eq. 3.3).
       - The level `0.001` is on `log₁₀|Δu|`, not on bare `|Δu|`
         (md:208).
       - FW's grid is `161×161` over `[-10, 10]²`, implicit
         `h_grid = 0.125` — the paper *never states this number*, must
         be inferred from md:145 + md:153.
       - FW say nothing about h-scaling.  Their `0.001` is a single
         calibration point.  The "proportionality constant in the
         `O(h²)`" phrase (md:206-208) is the truncation prefactor
         (`h²/12)·|u⁗|`, problem-dependent.

  3. **Caller catalogue across the repo**.  Only three call sites
     pass an explicit numeric `level`/`edge_level`:
     `figures/fw2011_fig_3_3.jl:109` (FW reference figure, deliberately
     `0.001`); `test/edge_detector_test.jl:200-201` (API round-trip
     test); `test/edge_detector_test.jl:206-207` (`-3.0`, `3.0`
     deliberate mutation sentinels).  Everything else — every figure
     script, every test exercising `edge_gated_pole_field_solve`,
     `examples/tritronquee_3d.jl` itself — uses the default.  Changing
     the default is therefore a real semantic change at all these call
     sites.

  4. **Existing-test prediction**.  The agent's per-test impact
     analysis surfaced the key risk: tests at `h_grid > H0 = 0.25` —
     PT.1.* (`h_grid ≈ 0.333`), EG.1.* (`h_grid = 0.75`), FF.2.*
     (`h_grid = 0.5`) — would see a *more restrictive* default if the
     formula extrapolated above `H0`.  Verified by reading the test
     files directly afterwards (test agent had conflated the Padé
     step `h = 0.5` with the lattice spacing `h_grid` for some
     entries).

The agent reconnaissance took ~80s wall clock total in parallel; doing
it sequentially or by main-thread `grep` would have been ~10× slower
and surfaced less.

## Design synthesis — the clamp is the load-bearing detail

The bead's empirical sweep (PI tritronquée at `[-20, 20]²`):

  - `h = 0.25` (N = 161): `level = 0.001` works (24 passes, ≈ 18 % field).
  - `h = 0.125` (N = 321): `level = 0.001` STALLS; `level ≈ -0.6` works.
  - `h = 0.0625` (N = 641): `level = 0.001` STALLS; `level ≈ -1.0` works.

These three points match `LEVEL0 + 2·log₁₀(h/H0)` with
`(H0, LEVEL0) = (0.25, 0.001)` to two decimals:
`log₁₀(0.5)·2 = -0.602`; `log₁₀(0.25)·2 = -1.204`.

But the empirical sweep is **silent on `h > H0`**.  Extrapolating the
formula upward gives:

  - `h = 0.333`: `level = +0.250` (vs FW's `0.001`)
  - `h = 0.5`:   `level = +0.603`
  - `h = 0.75`:  `level = +0.955`
  - `h = 1.0`:   `level = +1.204`

A back-of-envelope on the pole-annulus residual at `h = 0.75`:
`|∇²(1/z)| ≈ 2/|z|³`, at a cell `d = h` from a simple pole
`|Δu| ≈ 2/h³ ≈ 4.7`, `log₁₀ ≈ 0.67` — *below* the unclamped threshold
`0.955`.  The unclamped formula would collapse the mask on the existing
edge-gated-solve test grid (`h_grid = 0.75` in `test/edge_gated_solve_test.jl`).

**Decision: clamp the formula at `H0`**:

```julia
level(h) = LEVEL0 + 2·log₁₀(min(h, H0) / H0)
```

For `h ≥ H0` this returns `LEVEL0` (FW's published value).  For `h < H0`
it drops by `2·log₁₀(h/H0)` per the analytical h² scaling.  Backward-
compatible: every existing test at `h_grid ≥ H0` is byte-identical
under the new default.

Mutation M3 (worklog §"Mutation-proof") proves the clamp is load-
bearing: dropping it cascades into RED on EG.1.1 + EG.1.2 at h_grid = 0.75.

## What shipped

### `src/EdgeDetector.jl` — the fix

Module-level calibration constants `_H0 = 0.25`, `_LEVEL0 = 0.001`.
New internal helper `_auto_level(h, ::Type{T}) -> T`.  Default kwarg
on `pole_field_mask` changes:

  - 2-arg form (with `h` in scope): `level::Union{Real,Symbol} = :auto`;
    resolved via `_resolve_level(level, h, T)` which passes numeric
    `Real` through and resolves `:auto` to `_auto_level(h, T)`.
  - 1-arg form (`Δu` only, no `h`): `level::Union{Real,Symbol} = :auto`;
    `:auto` is unresolvable, throws `ArgumentError` per Rule 1.
  - Unknown `Symbol` sentinels throw `ArgumentError` (e.g.
    `pole_field_mask(u, h; level = :bogus)`).

Internal `_mask_from_residual(Δu, level_T)` is the shared resolve →
threshold core; the previous two-method copy-paste collapses to one.

Module docstring expanded with an "Auto-scaled level" section that
explains the `h²` derivation, the FW-vs-PI anchor distinction, and
the clamp.

### `src/EdgeGatedSolve.jl` — propagation

`edge_gated_pole_field_solve(...; edge_level::Union{Real,Symbol} = :auto, ...)`.
Inline docstring updated to point at `EdgeDetector`'s "Auto-scaled
level" section.

### `src/LatticeDispatcher.jl` — propagation + guard widening

Same kwarg type widening.  The existing `isfinite(edge_level)`
fail-fast guard becomes `(edge_level === :auto || (edge_level isa Real
&& isfinite(edge_level)))` — accepts the sentinel, still rejects NaN
/ Inf.

### `docs/adr/0009-edge-detector-h-aware-level.md` (Accepted)

The decision record: context, the `(H0, LEVEL0) = (0.25, 0.001)`
anchor, the clamp, four rejected alternatives.

### `test/edge_detector_test.jl` — three new testsets, 22 assertions

  - **ED.5.1** (10 assertions) — `_auto_level(h, T)` unit test:
    at the anchor returns `LEVEL0`; for `h < H0` drops by
    `2·log₁₀(h/H0)`; for `h ≥ H0` clamps to `LEVEL0`; element-type
    preservation under `BigFloat`-256.

  - **ED.5.2** (9 assertions) — pole-annulus integration test on a
    `u(z) = 1/(z - z₀)` grid at `h ∈ {0.25, 0.125, 0.0625}`: a FIXED
    spatial cell at distance `d = 0.5` from the pole is flagged by
    the new auto default at all three `h`, while the OLD fixed
    `level = 0.001` default misses it at `h = 0.0625` — the
    demonstration of the bug, side by side with the fix.

  - **ED.5.3** (4 assertions) — `:auto` sentinel propagation surface:
    `pole_field_mask(u, h) == pole_field_mask(u, h; level = :auto)`;
    numeric override still works; unknown Symbol fail-loud; 1-arg
    form fail-loud on `:auto`.

  - **ED.2.2 updated** — preserves the original FW Fig 3.3 reproduction
    intent by passing `level = 0.001` *explicitly* (the new auto
    default at `h = 0.05` is `-1.4`, which would flag every interior
    cell of the `1/z` grid and break the test's `count < length/2`
    ceiling).  Comment ties the change to bead `padetaylor-f8l`.

### `examples/tritronquee_3d.jl` — the headline render

`N = 161` → `N = 641` (4× linear resolution, `h_grid = 0.0625`).
Header comment table lists the three useful resolutions; the
`N = 641` line is annotated as "possible only with the auto-scaled
`edge_level`."  No manual `edge_level` override needed in the
script — the new default just works.

## Mutation-proof

Four mutations on `_auto_level` / the default kwarg, verified against
the new ED.5.* tests:

  - **M1** — revert `_auto_level` to `T(_LEVEL0)` (no h-scaling).
    Bites 4 assertions: 3 in ED.5.1 (the h < H0 scaling checks),
    1 in ED.5.2 at `h = 0.0625` (the annulus cell drops below threshold).

  - **M2** — flip sign: `T(_LEVEL0) - 2*log10(min(...)/...)`.  Bites
    5 assertions: 3 in ED.5.1, 2 in ED.5.2 (the sign-flipped
    threshold becomes MORE restrictive at fine grids, opposite of the
    needed correction).

  - **M3** — remove the clamp: `T(_LEVEL0) + 2*log10(T(h)/T(_H0))`
    (no `min`).  Bites 4 ED.5.1 clamp assertions AND cascades to the
    existing `edge_gated_solve_test.jl`: EG.1.1 + EG.1.2 fail at
    `h_grid = 0.75` (the unclamped threshold rises to `0.955` and
    the mask collapses).  The cross-file cascade proves the clamp is
    load-bearing for backward compat — *not* an aesthetic choice.

  - **M4** — revert the 2-arg `pole_field_mask` default from `:auto`
    back to `0.001`.  Bites 2 assertions: ED.5.2 at `h = 0.0625`
    (default no longer h-scales), ED.5.3 (the `pole_field_mask(u, h)
    == pole_field_mask(u, h; level = :auto)` equality fails since the
    default is no longer `:auto`).

Each mutation applied, fast-loop tested via `include("test/edge_detector_test.jl")`
(~1.1 s baseline, ~2 s under mutation), confirmed RED, reverted.
Procedure recorded in `test/edge_detector_test.jl` footer (alongside
the original 2026-05-13 stencil mutations A + B).

## Frictions surfaced

  - **The subagent test-impact prediction conflated Padé `h` with
    lattice `h_grid`.**  Sonnet agent 4's per-test risk analysis used
    "h = 0.5" for EG.1.* (the Padé step) when the actual lattice
    spacing is 0.75 (the `xs = range(-12, 12; length = 33)` step).
    Caught by reading the test files directly afterwards.  Lesson:
    subagent output is a map, not a verdict; the file-level details
    that determine "does this test bite" need first-hand reading
    before designing the fix.  Recorded in CLAUDE.md's broader
    "verify subagent output" Rule 3.

  - **The FW calibration anchor is `(0.125, 0.001)`, not the bead
    author's `(0.25, 0.001)`.**  Sonnet agent 2's read of FW 2011 §3.2.2
    pinned this from md:145 + md:153 (the 161×161 grid over `[-10, 10]²`).
    Two anchors are consistent with the same `h²` scaling law but
    they imply different problem-dependent proportionality constants
    `|u⁗(z)|`.  We pick `(0.25, 0.001)` because *that* is the
    empirically-validated calibration point on the PI tritronquée
    we ship as `tritronquee(:I)`.  FW's own anchor remains recoverable
    by passing `level = 0.001` explicitly (e.g.
    `figures/fw2011_fig_3_3.jl:109`).  Recorded in ADR-0009's
    "Alternatives B" rejection.

  - **The `_auto_level` helper is exported as an internal `_`-prefixed
    name but referenced from the test file via `EdgeDetector._auto_level`.**
    Acceptable for v1 (unit tests for internal helpers are common
    Julia idiom); if a v2 user ever needs to introspect the
    calibration anchor, expose it as a public API instead.

## Hard-won lesson

**A "deep bug" disguised as a kwarg-default cleanup is still a deep
bug.**  The original bead description framed this as "make `level`
h-aware" — sounds like a one-line kwarg change.  But the real work
was:

  1. The h² scaling derivation (FW 2011's eq. 3.3 truncation prefactor
     for an analytic function, in *log* space → additive shift by
     `2·log₁₀(h)`).
  2. The empirical calibration anchor (PI tritronquée at `[-20, 20]²`,
     pinned by the bead's three-point sweep).
  3. The recognition that FW's own anchor is *different* from ours
     (and that this is OK — problem-dependent constant in the
     proportionality, recoverable by explicit override).
  4. The clamp at `H0`, defended by the math (the smooth/annulus gap
     is wide at coarse grids; raising the threshold there overshoots
     the pole-annulus residual at `h = 0.75`).
  5. Mutation M3 cascading into existing tests at `h_grid = 0.75`, as
     a *cross-test* proof that the clamp is load-bearing.

Anything less than (5) would have shipped a fix that broke EG.1.x +
FF.2.x silently.  CLAUDE.md Rule 2 "all bugs are deep" earned its
keep this session.

## What this enables

`examples/tritronquee_3d.jl` at `N = 641` (`h_grid = 0.0625`,
4× the v0.1.0 hero resolution).  The README hero can now be
arbitrarily-resolution as a function of available wall time, not
gated on a calibration-anchor accident.  Future high-resolution
figures (FW Fig 4.3/4.4 quantitative pin on bead `padetaylor-p3l`,
or anything that wants `N > 200` on a `[-10, 10]²`-scale window) are
unblocked.
