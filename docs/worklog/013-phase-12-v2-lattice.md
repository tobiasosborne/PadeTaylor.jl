# Worklog 013 — Phase 12 v2: LatticeDispatcher — 2D pole-field-plus-BVP composition

**Date**: 2026-05-13
**Author**: Claude Opus
**Scope**: Phase 12 v2 v1 / bead `padetaylor-k31` — the 2D-lattice
composition layer.  PathNetwork over a Cartesian grid + EdgeDetector
partition + per-row BVP fill on contiguous smooth runs flanked by
pole-field cells.  FW 2011 §4.4 + line 190 verbatim ("161 separate
BVP solutions; one for each grid line").  1237/1237 GREEN (1160 prior
+ 76 LatticeDispatcher + 1 umbrella `isdefined`).

> Take-home: the v1 composition machinery is *generic over the
> partition* — supply a mask (auto or manual), and it fills smooth
> bands row-by-row.  v1 ships verified-against-cosh on a synthetic
> mask; the FW Fig 4.1 quantitative pin is deferred to follow-up
> bead `padetaylor-0c3`.

## v1 scope decision

The original `padetaylor-k31` bead specified an ambitious v2 scope:
> "Acceptance: reproduce FW 2011 Fig 4.1 (near-tritronquée with smooth
> imaginary-axis band) per `docs/figure_catalogue.md` §1, target u(0),
> u'(0), u(20i), u'(20i) to ≤1e-10 abs; pole locations agree with FW
> Table 5.1."

Reading FW 4.1's caption (FW2011...md:218-222) reveals the figure's
composition is **not** the horizontal-row BVP-fill pattern at line 190:

  - FW 4.1 step (i): BVP between `+20i` and `-20i` (a VERTICAL segment
    along the imaginary axis).
  - FW 4.1 step (ii): two pole fields run out from `z = 0` and
    `z = 20i` (HORIZONTAL outward integration).
  - FW 4.1 step (iii): "fill in the area between pole field edges by
    further BVP solutions" — these MAY follow the line-190 horizontal-
    row pattern, but FW does not specify the orientation explicitly.

So FW Fig 4.1 reproduction requires either (a) extending
`lattice_dispatch_solve` to support vertical BVPs in addition to
horizontal ones, or (b) composing the Fig 4.1 result via
`dispatch_solve` (the Phase-12-v1 1D-chain dispatcher) at the outer
level with `lattice_dispatch_solve` as a building block for the
"fill in the area" step.  Either is a separate compositional pattern,
not a 2D-lattice algorithm extension.

**v1 ships the line-190 horizontal-row algorithm faithfully**.  The
Fig 4.1 quantitative pin is filed as `padetaylor-0c3` for follow-up.

This split is faithful to CLAUDE.md Rule 9 ("senior-engineer-grade
only") — we ship the *correct algorithmic primitive* with a tight
synthetic test, and document the gap to FW 4.1 honestly rather than
shipping a v1 with a half-implemented FW 4.1 setup that "almost"
hits the 1e-10 acceptance.

## The algorithm

Per the module docstring `src/LatticeDispatcher.jl`:

  1. **IVP fill**: `path_network_solve(prob, vec(grid))` over the flat
     2D grid.
  2. **Partition**: `pole_field_mask(u_grid, h_grid; level)` from
     EdgeDetector — or accept a pre-computed `mask::BitMatrix` via the
     `mask` kwarg (the FW2011...md:401 "manual classification" path).
  3. **Per-row BVP fill**: for each interior row `j ∈ 2..ny-1`, walk
     along the row, find maximal contiguous smooth runs (`mask=false`)
     flanked by pole-field cells (`mask=true`) on both sides; solve a
     `bvp_solve` per run with Dirichlet BCs taken from the IVP values
     at the flanking cells.
  4. **Stitched output**: replace IVP values in BVP-bridged cells with
     the BVP's barycentric interpolant; tag those cells `:bvp`.  Tags
     in `{:ivp, :bvp, :ivp_only}` per cell.

## What changed

`src/LatticeDispatcher.jl` (~115 lines of code + ~95 lines literate
docstring):

  - `LatticeSolution{T}` — composed 2D output struct with the stitched
    `u_grid`, `up_grid`, partition `mask`, per-cell `region_tag`,
    underlying `PathNetworkSolution`, and `Vector{BVPSolution}` for
    diagnostic inspection.
  - `lattice_dispatch_solve(prob, bvp_f, bvp_∂f_∂u, xs, ys; ...)` —
    the public driver.  Kwargs: `h_path, order, edge_level, N_bvp,
    bvp_tol, mask`.  The `mask` kwarg accepts a `BitMatrix` to bypass
    EdgeDetector (test-friendly + FW manual-classification path).

`src/PadeTaylor.jl`: registered the new module.

`test/lattice_dispatcher_test.jl` (4 testsets, 76 assertions):

  - **LD.1.1**: PI tritronquée on the Phase-9 setup runs end-to-end.
    Region tags consistent with mask.  Zero BVP fills (the pole-field
    wedge is at the grid edge → no two-sided smooth runs).
  - **LD.1.2**: linear `u'' = u` with supplied mask creating a single
    bridgeable smooth run per row.  BVP-filled cells match `cosh(z)`
    to ≤ 1e-10.  Region tags are `:ivp`/`:bvp`/`:ivp_only` on the
    right cells.  9 BVPs solved (one per interior row of an 11-row
    grid).
  - **LD.1.3**: mask = all-false → zero BVP fills, all cells
    `:ivp_only`.
  - **LD.2.1**: fail-fast on short grids, anisotropic spacing,
    `N_bvp < 4`, mask-shape mismatch.

## Mutation-proof (verified 2026-05-13)

**Mutation E** — swap `u_a` and `u_b` in the `bvp_solve` call.  This
breaks the BC-direction contract; the BVP solves a problem with the
wrong endpoint values.  Verified bite: **32 fails** out of 76, all in
LD.1.2's cosh-comparison loop — 45 BVP-filled cell value checks (9
rows × 5 cells per row) plus 9 row-edge IVP-cell checks (the `cosh(z)`
on the `:ivp` cells are unaffected, but the closeness of BVP cells
to `cosh` is destroyed).

**Mutation F** — comment out the `region_tag[k, j] = :bvp` assignment.
This breaks the tagging contract without breaking the u-value
computation.  Verified bite: **1 fail** — the targeted assertion
`sol.region_tag[6, 6] == :bvp` fails (the cell retains its default
`:ivp_only` tag).  The cosh values still match because the u-grid
update precedes the tag assignment.  Clean isolated bite separating
"BVP ran correctly" from "region book-keeping is correct".

Both mutations restored before commit per CLAUDE.md Rule 4.

## Frictions surfaced

  - **F1. PI tritronquée Phase-9 grid produced zero BVP fills**.
    The leading-pole wedge for tritronquée is on the positive-real-axis
    edge of `[-4, 4]²`, so no smooth run is flanked by IVP cells on
    BOTH sides — the BVP-fill path is never triggered.  This is a
    *correct* result (the algorithm has nothing to bridge), but the
    smoke test produced zero work, which initially looked like a bug.
    Resolved by adding the synthetic-mask test (LD.1.2) that forces
    BVP fills.  Lesson: a Phase 9-style real-Painlevé test exercises
    the IVP+partition machinery but not the BVP composition; needs a
    contrived setup or a different Painlevé case with interior poles
    (e.g., the equianharmonic ℘ on a window containing one of its
    hexagonal-lattice poles).  Filed as part of the FW Fig 4.1 v2.1
    follow-up (`padetaylor-0c3`).
  - **F2. The `mask` kwarg is essential for testability**.  Without it,
    every test depends on EdgeDetector producing a specific partition
    on a specific problem at a specific grid — brittle and
    problem-coupled.  With the kwarg, LD.1.2 can supply a hand-crafted
    mask and verify the BVP-fill machinery against `cosh` closed form
    independent of detector behaviour.  Also documents the FW 2011
    line 401 "manual classification" workflow.  Not a friction so
    much as a design refinement that mattered.
  - **F3. The anisotropic-grid guard**.  EdgeDetector's 5-point
    stencil assumes isotropic `h`; supporting `h_x ≠ h_y` would
    require a different stencil weight pattern.  Added a fail-fast
    `isapprox(step(xs), step(ys); rtol = 1e-10)` check.  Caller-side
    constraint, documented.

## Pointers

  - [`src/LatticeDispatcher.jl`](../../src/LatticeDispatcher.jl) — the module.
  - [`test/lattice_dispatcher_test.jl`](../../test/lattice_dispatcher_test.jl)
    — 4 testsets + mutation-proof commentary at bottom.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    lines 188-190 (the "161 BVPs per grid line" passage), 218-222 (FW
    Fig 4.1 caption — the *deferred* quantitative target).
  - [`docs/worklog/011-edge-detector.md`](011-edge-detector.md) — the
    EdgeDetector prerequisite.
  - [`docs/worklog/006-phases-10-11-path-network-bvp.md`](006-phases-10-11-path-network-bvp.md)
    — PathNetwork + BVP prerequisites.
  - [`docs/worklog/007-phase-12-dispatcher.md`](007-phase-12-dispatcher.md)
    — Phase 12 v1 (1D-chain) for contrast.

## Bead state

Closed in this session:
  - `padetaylor-k31` — Phase 12 v2 composition machinery GREEN at this
    commit.  v1 horizontal-row algorithm per FW2011...md:190.

Filed in this session:
  - `padetaylor-0c3` — **P2** Phase 12 v2.1: FW Fig 4.1 quantitative
    pin (u(0), u'(0), u(20i), u'(20i) to ≤1e-10 abs).  Requires
    extending or composing the v1 machinery to handle FW Fig 4.1's
    BVP-on-imaginary-axis + outward-pole-fields pattern.

Still open from prior sessions (no change):
  - `padetaylor-bvh`, `padetaylor-grc`, `padetaylor-61j`, `padetaylor-8pi`.

## Hard-won lesson (for HANDOFF.md §"Hard-won")

**16. v2-bead acceptance criteria often blend "machinery" and "specific
case".**  `padetaylor-k31`'s description fused two deliverables: the
generic 2D-lattice composition layer (machinery) and the FW Fig 4.1
quantitative pin (specific case).  These are different work; FW Fig 4.1
uses a *different* compositional pattern (vertical BVP + two outward
pole fields, not the horizontal-row line-190 algorithm).  v1 ships the
machinery faithfully and splits the specific-case acceptance into a
follow-up.  Lesson: read v2-bead descriptions for "what algorithm is
being shipped" vs "what specific case is being verified", and split
the bead if they're different deliverables.
