# Worklog 030 — `PainleveSolution`: the self-describing solve-output wrapper

**Date**: 2026-05-14 (continues worklog 029)
**Author**: Claude Opus
**Beads**: `padetaylor-p1b`, `padetaylor-26r`, `padetaylor-ylr` — all
closed.  `padetaylor-soi` closed (subsumed by ADR-0007).
**Scope**: Zoom out from FW-figure reproduction and build Painlevé
*infrastructure* — the output-side counterpart to ADR-0006's
`PainleveProblem` builder.

> **Take-home**: `solve_pade` / `path_network_solve` on a
> `PainleveProblem` now return a `PainleveSolution` — a thin wrapper
> carrying the Painlevé identity + coordinate frame around whichever
> raw type the solver produced, with a uniform `z`-frame access
> surface (`sol(z)`, `poles`, `grid_values`, `equation` /
> `parameters` / `solutionname`, `show`).  The `:transformed`
> equations (PIII / PV / PVI) stop being second-class: the wrapper
> carries `from_frame`, so `path_network_solve` now works end-to-end
> for them.  Test suite **1526 → 1590 GREEN**; ADR-0007 accepted.

## Ground truth read first

- `src/Painleve.jl`, `docs/adr/0006-painleve-problem-layer.md` — the
  `PainleveProblem` builder this mirrors.  ADR-0006 refinement #2 (the
  "`:transformed` forwarding throws" decision) is **superseded** by
  ADR-0007.
- `src/Problems.jl:PadeTaylorSolution`, `src/PathNetwork.jl:
  PathNetworkSolution`, `src/PoleField.jl`, plus the `BVPSolution` /
  `DispatcherSolution` / `LatticeSolution` / `EdgeGatedSolution`
  structs — the six unrelated raw solve-output types the wrapper
  unifies.
- `src/CoordTransforms.jl:120-188` — the `pIII_z_to_ζ` / `pV_z_to_ζ`
  transform pairs.  Confirmed at coding time that the coordinate
  component of each map depends *only* on the coordinate argument
  (`ζ = 2 log z` ignores `u`, `u'`), which is what makes the
  `_coord(map, c) = map(c, zero(c), zero(c))[1]` projection idiom
  valid — the same idiom `_build_transformed` already uses for
  `ζ_end`.
- `ext/PadeTaylorArblibExt.jl`, `ext/PadeTaylorCommonSolveExt.jl`,
  `figures/figutil.jl` — the extension pattern and Makie-plotting
  conventions the new `PadeTaylorMakieExt` follows.

## What shipped

### `src/PainleveSolution.jl` (new) — bead `padetaylor-p1b`

`PainleveSolution{S, TF, FF}` — a wrapper struct (`equation`, `params`,
`name`, `frame`, `to_frame`, `from_frame`, `raw`) `include`d into
`module Painleve` so `PainleveProblem` and `PainleveSolution` share one
namespace while each file stays ≤200 LOC (CLAUDE.md Rule 6 is
per-file).  The `z`-frame access surface:

  - **callable `sol(z)`** — `:direct` forwards to `sol.raw`;
    `:transformed` maps `z → ζ`, evaluates, maps `(w,w') → (u,u')`
    back.  Fails loud for a grid-type raw (a `PathNetworkSolution` is
    not a dense interpolant — Rule 1).
  - **`poles(sol)`** — forwards to `PoleField.extract_poles`,
    `from_frame`-mapped for `:transformed`.  Throws (not "returns
    empty") for an unwired raw type.
  - **`grid_values(sol)`** — `(z,u,u')` point values in the `z`-frame;
    for `:transformed` each `(ζ,w,w')` triple is mapped back.  This is
    what makes a `:transformed` path-network result *consumable* (its
    `raw.grid_*` are `ζ`-frame).
  - **`equation` / `parameters` / `solutionname`** accessors;
    **`show`** (multi-line provenance + one-line form).

`name` is `nothing` until the deferred named-transcendent constructors
land — storing the field now keeps that thread purely additive.

### `src/Painleve.jl` — forwarding rewrite

`solve_pade` / `path_network_solve` on a `PainleveProblem` now return a
`PainleveSolution`.  The `:transformed` story (ADR-0007's "partial
flip"):

  - **`path_network_solve` flips fully** — it handles the complex
    `ζ`-domain, so it maps the caller's `z`-frame grid into the
    `ζ`-frame, solves, and the returned wrapper presents `z`-frame
    `poles` / `grid_values`.
  - **`solve_pade` flips partially** — it does fixed-step *real-axis*
    stepping (`state.z < z_end` is undefined for complex `z`), so it
    serves a `:transformed` problem only when the `ζ`-domain is
    real-typed; the common complex-`ζ` case throws, now pointing at
    `path_network_solve` rather than the old "cannot round-trip".

The `_forward_guard` helper (which made *all* `:transformed`
forwarding throw) is deleted.

### `src/PoleField.jl` — bead `padetaylor-26r`

`extract_poles` gained a `PadeTaylorSolution` method.  A single
`solve_pade` trajectory carries the same per-segment Padé store a
path-network does (`sol.pade[k]` centred at `sol.z[k]` with step
`sol.h[k]`), so the extraction is identical — the far-root /
Froissart-residue / greedy-`|t*|`-clustering logic was refactored into
a shared `_extract_poles_core`, and the two public methods are thin
adapters handing it the right three arrays.  The one changed default
is `min_support = 1` for the trajectory method: a single trajectory is
a *chain*, not a *fan*, so a pole is typically bracketed by one or two
consecutive segments — demanding three independent sightings (the
path-network default) would discard every real pole.

### `ext/PadeTaylorMakieExt.jl` (new) — bead `padetaylor-ylr`

A Makie plot recipe, `painleveplot(sol::PainleveSolution)`, loaded when
`Makie` is present.  Presentation-only (ADR-0003): it reads the wrapper
through `grid_values` / `poles` and renders one complex-`z`-plane
`Axis` — a trajectory path for a `PadeTaylorSolution` raw, a coloured
grid scatter for a `PathNetworkSolution` raw, extracted poles overlaid
as red diamonds.  Because it consumes the already-`z`-frame surface, it
is automatically correct for `:transformed` problems.  `Project.toml`
gained `Makie` in `[weakdeps]` / `[extensions]` / `[compat]` (0.24) /
the test target.

### `docs/adr/0007-painleve-solution-wrapper.md` (new, Accepted)

Records the design, the packaging decision (one module / two files),
and the supersession of ADR-0006 refinement #2.  Rejected
alternatives: self-describing the raw types in place, a separate
`PainleveSolutions` module, keeping `:transformed` forwarding throwing,
and a `<: PadeTaylorSolution` subtype.

### Tests — 1526 → 1590 GREEN

  - `test/painleve_solution_test.jl` (new) — `PS.1`–`PS.6`: provenance
    + accessors, the callable (`:direct` transparency, `:transformed`
    IC round-trip, grid-type fail-loud), `poles` forwarding + mapping,
    `show`, `grid_values` (`:direct` passthrough, `:transformed`
    mapping, end-to-end `:transformed` grid round-trip), and the
    unwired-raw fail-loud guards.
  - `test/polefield_test.jl` — `PF.3.1` / `PF.3.2`: extraction from a
    single `solve_pade` trajectory + the `min_support = 1` default
    being load-bearing.
  - `test/ext_makie_test.jl` (new) — `MK.1.*`: `painleveplot` builds a
    titled Figure for both raw types; `show_poles` toggles the overlay.
  - `test/painleve_test.jl` — `PV.3.1` / `PV.4.1` rewritten: forwarding
    now returns a `PainleveSolution` (the old tests read `.grid_u` /
    `.y` straight off the result).

Mutation-proof: ten mutations across the three beads, each confirmed
RED and reverted; procedures in the three test files' footers.  Two
mutation-discovery notes below.

## Frictions

  - **`solve_pade` cannot serve a complex-`ζ` `:transformed` problem.**
    The ADR's first draft said both forwarding methods "no longer
    throw" for `:transformed`.  Caught at design time by reading
    `solve_pade`: it does real-axis stepping (`state.z < z_end`), and a
    `:transformed` problem's `ζ`-domain is in general complex.  The ADR
    was corrected *before* coding to the "partial flip" — `solve_pade`
    still throws for complex `ζ`, now with a `path_network_solve`-
    pointing message.  (Law 1 paid off: the friction surfaced from
    *reading the affected file*, not from a test failure.)

  - **A forwarding-method mutation needs an end-to-end test.**  The
    first draft of `painleve_solution_test.jl` mutation M3 (break the
    `:transformed` grid mapping in `path_network_solve`) claimed it bit
    `PS.5.2` — but `PS.5.2` wraps a *hand-built* raw and so never
    exercises the forwarding method.  Added `PS.5.3`: a real
    `:transformed` PV `path_network_solve` whose invariant is "the
    `z`-frame grid you get back equals the grid you put in".  *Then* M3
    bit.  Lesson: a mutation in a forwarding/glue function must be
    pinned by a test that goes *through* that function, not one that
    constructs the post-glue state directly.

  - **M2 bit two tests, not one.**  Dropping the `_assert_callable`
    guard takes down both `PS.2.2` (a real grid-backed solution) and
    `PS.6.1` (the bogus-raw wrapper) — the footer was corrected to
    record both.  Running the mutation, not predicting it, is the
    discipline (worklog 021 lesson 39).

## What remains genuinely out of scope

  - **Named-transcendent constructors** (`tritronquee(:I)`,
    `hastings_mcleod()`, …) — the `name` field is in place to receive
    them, but each named solution needs its defining ICs/asymptotics
    grounded in the literature.  ADR-0006 "Deferred".
  - **`PainleveProblem` forwarding for the other four solvers**
    (`bvp_solve`, `dispatch_solve`, `lattice_dispatch_solve`,
    `edge_gated_pole_field_solve`) — `PainleveSolution{S}` is already
    generic over their output types; what is missing is the
    `PainleveProblem`-accepting forwarding methods.
  - **A unified `solve(pp; over=…)`** that picks IVP / path-network /
    BVP from the requested domain.

## Hard-won lesson

**The wrapper you build to map frames is also the place the mapping
should happen — so use it.**  ADR-0006 refinement #2 declared
`:transformed` forwarding un-round-trippable because a `ζ`-frame
`PathNetworkSolution`'s *Padé store* has no `z`-frame image.  True —
but it conflated the store with the *point values*, which map back
cleanly.  Once `PainleveSolution` exists and carries `to_frame` /
`from_frame`, the honest move is to do the point-value mapping at the
wrapper boundary (callable, `poles`, `grid_values`) and leave the Padé
store as the documented `ζ`-frame escape hatch in `sol.raw`.  Building
the abstraction and then *not* using it for half the equations would
have shipped the package with PIII / PV / PVI still second-class.
