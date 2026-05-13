# Worklog 018 — SheetTracker: PVI ζ-plane RHS + winding primitives (bead `padetaylor-grc`)

**Date**: 2026-05-13 (late evening, follow-up to worklog 017)
**Author**: Claude Opus
**Scope**: Phase 14 / Tier-5 v1 — ship the PVI ζ-plane transform per
FFW 2017 eq. (3) (md:144) and winding-number primitives that compute
Riemann-sheet indices from a traversed path.  New module
`src/SheetTracker.jl` (~150 LOC) + 7 testsets / 26 assertions.  Three
mutations bit as predicted.  GREEN at 1285 → 1311.

> Take-home: PVI has TWO fixed branch points (`z=0` AND `z=1`), so no
> single exponential transform removes both (Picard's theorem).  This
> module ships the FFW ζ-plane transform that removes `z = 0` —
> leaving the second branch as a `ζ = 2π·i·k` lattice on the
> imaginary axis — plus path-side winding-number primitives that let
> callers compute sheet indices after a regular PathNetwork walk.
> Constrained-wedge PathNetwork routing (refuses to overstep cuts;
> increments sheet counter on deliberate circumambulation) is the
> production-grade version of FFW 2017 §2.2.2 and remains deferred.

## What changed

  - `src/SheetTracker.jl` (new, ~150 LOC):
    - `pVI_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''`
      Factory returning the ζ-plane P̃_VI RHS closure (FFW eq. 3, md:144).
    - `winding_delta(z_old, z_new, branch) -> Float64`
      Signed angle change normalised to `(-π, π]`.  Caller contract:
      step << distance-to-branch.
    - `accumulate_winding(path, branch) -> Vector{Float64}`
      Cumulative winding along `path` w.r.t. `branch`.  `out[1] == 0`.
    - `sheet_index(total_winding) -> Int`
      `round(Int, total / 2π)`.  FFW2017...md:187-189 sheet convention.

  - `src/PadeTaylor.jl` — `include("SheetTracker.jl")` + re-exports;
    module-header architecture list row 13.

  - **PVI coordinate conversion is mathematically identical to PV's**
    (`u(z) = w(ζ)`, `z = e^ζ`; FFW2017...md:146).  Module docstring
    documents this and instructs callers to reuse `pV_z_to_ζ` /
    `pV_ζ_to_z` from `CoordTransforms`.  No new helper added (would
    be dead-aliased code).

  - `test/sheet_tracker_test.jl` (new, ~200 LOC, 7 testsets / 26
    assertions):
    - **ST.1.1** — PVI ζ-plane RHS at `(ζ=log 2, w=1/2, w'=1)` with
      degenerate parameters `α=β=γ=δ=0`: hand-pinned `w'' = -1`.
    - **ST.1.2** — same point with `α=β=γ=δ=1`: hand-pinned `w'' =
      101/24`.
    - **ST.1.3** — end-to-end PVI direct-vs-transformed agreement at
      `z = 2 → z = 2.05` (away from `z ∈ {0, 1}`): `|Δu| ≤ 1e-10`.
    - **ST.1.4** — winding_delta sign + normalisation tests.
    - **ST.1.5** — accumulate_winding on CCW square (`+2π`), CW
      square (`-2π`), non-enclosing path (`≈ 0`).
    - **ST.1.6** — sheet_index conversion at canonical and edge inputs.
    - **ST.1.7** — branch-point lattice at `ζ = 2π·i·k`: RHS magnitude
      blows up to `> 1e10`; near the branch but not on it, RHS is
      finite + bounded.

  - `docs/figure_catalogue.md §6` row T5 — marked PARTIAL with the
    SheetTracker primitives shipped + deferrals (η-plane PVI eq.,
    constrained-wedge routing).

## Mutation-proof procedure

Three load-bearing mutations applied + restored.

**Mutation O** — in `pVI_transformed_rhs`, swap the `/2` factor on
the `(dw/dζ)²` term to `/3`.  **Verified bite**: ST.1.1 line 35 RED
(hand-pin `-1` fails), ST.1.2 line 51 RED (`101/24` fails), ST.1.3
lines 110-111 RED (end-to-end disagreement).  The RHS hand-pin AND
end-to-end-agreement tests each catch this independently.

**Mutation P** — in `winding_delta`, drop the wrap-to-(-π, π]
normalisation (return raw `Δθ`).  **Verified bite**: ST.1.4 line 136
(near-branch-cut hop's normalised value reads ≈ -2π absolute, fails
`abs(Δθ) < 0.05`); ST.1.5 lines 145 + 150 (closed-loop CCW + CW
accumulations both off by ±2π due to unnormalised intermediate
`Δθ`s); ST.1.6 line 172 (downstream `sheet_index` derives from the
broken `accumulate_winding`).  4 fails total.

**Mutation Q** — in `sheet_index`, swap `round` for `floor`.
**Verified bite**: ST.1.6 line 167 RED — `sheet_index(3π/2)` reads
`0` instead of `+1`.  Other ST.1.6 assertions stay GREEN (their
inputs at integer multiples of `2π` round and floor identically).
One assertion catching the mutation suffices for "load-bearing"; the
test could be strengthened with more sub-2π inputs.

All three mutations restored before commit; full suite GREEN at 1311.

## Frictions surfaced

1. **Initial ST.1.4 test had a wrong direction-assertion.**  My
   first cut asserted `Δθ < 0` for the near-branch-cut step from
   `(-1, +0.01)` to `(-1, -0.01)`, expecting "going clockwise across
   the cut".  But the normalisation reinterprets a `Δθ_raw ≈ -2π`
   step as a SMALL counterclockwise hop (`+0.02`) — that's the
   documented contract: callers must use small steps so single-step
   winding is unambiguous.  My test confused topological direction
   with the algorithm's same-side wrap convention.  Fixed by
   dropping the direction-pin (the magnitude check `|Δθ| < 0.05` is
   the real load-bearing claim).

2. **Initial ST.1.7 expected non-finite RHS at the branch point.**
   At `ζ = 2π·im` exactly in Float64, `exp(ζ) = 1 + O(eps)·im` —
   NOT exactly `1`, because `2π` is a rounded constant in Float64.
   So `(e^ζ - 1)` is `O(eps)`-small but non-zero, and the RHS blows
   up to `~1e30` (huge but finite).  Fixed the test to assert
   magnitude `> 1e10`, not non-finite-ness.  Downstream Padé-Taylor
   will throw on non-finite Taylor coefficients once the magnitude
   exceeds Float64 range — fail-loud still holds at the stepper layer.

## What is NOT shipped (deferral notes)

**η-plane PVI transform** (FFW 2017 eq. 5, md:154).  The second
exponential `ζ = e^η` makes the branch-point-free region `Re η <
log(2π)` more compact at the cost of nested `e^(e^η)` arithmetic.
Useful for FFW Fig 2 first column ONLY; the ζ-plane (eq. 3)
suffices for the other Tier-5 figures.  No bead filed; if a
downstream caller actually needs the η-plane, the transcription
work is well-scoped (~40 LOC + 4 testsets).

**Constrained-wedge PathNetwork routing.**  The production-grade
version of FFW §2.2.2 modifies `path_network_solve`'s wedge selector
to (a) refuse to overstep branch cuts and (b) increment the
sheet-index counter when the walker deliberately circumambulates.
The v1 primitives in this module compute sheet indices AFTER a
walk; they don't constrain wedge selection DURING the walk.  This
is the path-network change called out in the bead description
(`cross_branch=true` kwarg sketch).  Significant change to
PathNetwork.jl; not filed as a bead.

**Sheet-aware barycentric evaluation.**  Once visited-node sheet
indices are computed, Stage-2 evaluation should refuse to evaluate
`grid_z` cells whose sheet differs from the nearest visited node's
sheet (CLAUDE.md Rule 1 fail-fast on cross-sheet queries).  Not
shipped; depends on the routing change above.

## Design decisions

**Module is helpers + ζ-plane only, no driver.**  Same pattern as
CoordTransforms (worklog 017): ship the algebraic primitives,
defer the integration loop.  Callers compose with
`PadeTaylorProblem` + `path_network_solve`.  The η-plane equation
is significantly larger (eq. 5's nested `e^(e^η)` terms) — adding
it doubles the algebra surface area + bug exposure with marginal
v1 value.

**`winding_delta` normalises to `(-π, π]`, not `[-π, π)`.**  The
boundary case `Δθ = +π` (a step that crosses exactly through the
branch cut) stays positive (interpreted as "going counterclockwise
by π"); `Δθ = -π` flips to `+π`.  This is a convention choice;
the test ST.1.4 line 124 verifies it.  Documented in the module
docstring.  Symmetric alternative (`[-π, π)`) is also valid; we
picked the convention closer to `atan2`-style half-open.

**`sheet_index` uses `round`, not `floor` or `ceil`.**  `round` is
symmetric around integer 2π multiples; `floor` would map +π to 0
but -π to -1 (asymmetric), which doesn't match FFW's symmetric
sheet labelling.  Mutation Q confirms `round` is load-bearing.

## Beads

  - `padetaylor-grc` — closed in this session.
  - **Not filed**: η-plane PVI transform, constrained-wedge
    PathNetwork routing, sheet-aware Stage-2 evaluation.  All three
    are clean follow-on work but no in-tree caller is blocked.

## Pointers

  - `src/SheetTracker.jl` — module top docstring documents the
    transforms + primitives + deferrals.
  - `test/sheet_tracker_test.jl` — 7 testsets + mutation-proof
    procedure at end-of-file.
  - `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md`
    :29-37 (PVI equation), :135-160 (ζ + η transforms), :163-189
    (circumambulation + sheet index parametrisation), :191-195
    (Fig 2 caption with parameters + IC).
  - `docs/figure_catalogue.md §5` row T5 — PARTIAL acceptance.
  - `src/CoordTransforms.jl` — PV's `z↔ζ` map reused for PVI.

## Hard-won lessons (for HANDOFF.md §"Hard-won")

29. **"Branch point in Float64" is a magnitude check, not a
    finite-ness check.**  At `ζ = 2π·im` exactly, `exp(ζ)` returns
    `1 + O(eps)·im` (the rounded constant `2π` is non-exact), so
    `(e^ζ - 1)` is tiny-but-non-zero, and the RHS blows up to
    `~1e30` but stays technically finite.  Test the MAGNITUDE
    (`> 1e10`), not `isfinite`.  Downstream Padé-Taylor's
    finite-Taylor-coefficient check is the actual fail-loud point.

30. **`atan2`-style angle differences need (-π, π] normalisation to
    work as winding deltas.**  The raw `angle(z_new) - angle(z_old)`
    has `±2π` discontinuities at the branch cut.  Normalising each
    step to (-π, π] gives the "shortest signed path" interpretation,
    which is correct iff the path step is < π in angular size.
    Callers must respect that contract (e.g., walk with `h <<
    distance-to-branch`).

31. **Topological direction is not the same as the algorithm's
    same-side-wrap convention.**  A step from `(-1, +0.01)` to
    `(-1, -0.01)` is geometrically downward (clockwise around
    origin) but the normalised winding-delta reads `+0.02`
    (counterclockwise small).  The algorithm interprets every
    sub-π step as "going the short way around" the branch — only
    actual circumambulation (loops > π in angular extent built up
    over multiple steps) registers sheet changes.  Don't pin
    direction on individual cross-cut steps in tests.
