# Worklog 014 — PathNetwork shuffle-induced y-asymmetry; tritronquée plot session

**Date**: 2026-05-13 (late evening)
**Author**: Claude Opus
**Scope**: User-driven plotting session that surfaced two structural
bugs in `path_network_solve` when applied to 2D Cartesian grids on a
conjugate-symmetric problem (PI tritronquée).  No source-tree fixes
shipped — both bugs are filed (`padetaylor-dtj`, `padetaylor-0c3`) and
the examples script `examples/tritronquee_3d.jl` carries the
user-space workaround.

> Take-home: `path_network_solve(prob, vec(grid))` over a 2D Cartesian
> lattice for a real-IC + real-coeff ODE does NOT preserve the
> reflection symmetry `u(z̄) = ū(z)` of the analytic solution.  The
> `shuffle(rng, targets)` step (PathNetwork.jl:173) creates an
> asymmetric visited-tree which cascades 4-5 orders of magnitude into
> `|u|` at conjugate-pair cells.  Workaround: walk only `Im(z) > 0`,
> mirror to lower half via `u(z̄) = ū(z)`.  Bit-exact symmetric, ~2×
> faster, applicable only to real-coeff/real-IC problems.

## Session genesis

The user asked for a "3D tritronquée plot" — a routine visualisation
deliverable now that PathNetwork + EdgeDetector + LatticeDispatcher
are all shipped.  Three rounds of scaling (`[-4, 4]² → [-20, 20]² at
41×41 → 121×121 → 501×501`) surfaced increasingly-visible artifacts
that turned out to be structural.

## Bug 1 — shuffle-induced asymmetric visited tree

### Symptoms

At 121×121 over `[-20, 20]²`:
  - Visited-tree node count above real axis (2589) ≠ below (2563) —
    26-node asymmetry from a problem whose analytic solution is
    bit-exact `y → -y` symmetric.
  - Mask under `y → -y` reflection: 2254 mismatches out of 9446
    flagged cells (24%).
  - `max |Δ|u|| over y-reflection`: 14142 at z = 13.33 ± 15.67i; the
    upper cell reads `|u| = 14144`, the lower `|u| = 1.85` — a 4-order
    magnitude divergence at conjugate points.

Verified with four `rng_seed` values (0, 1, 7, 42): all give non-zero
asymmetry, of *different* magnitudes — the asymmetry is shuffle-order-
dependent, not deterministic.

### Mechanism

`path_network_solve` (PathNetwork.jl:173) calls
`targets = shuffle(rng, collect(CT, grid))` before processing.  The
visited-tree is built incrementally — each new target walks from the
already-built tree's nearest visited node.  The wedge walker itself
is symmetry-preserving (the default wedge angles `[-π/4, -π/8, 0,
π/8, π/4]` are symmetric, and complex arithmetic on real-coeff inputs
is conjugation-equivariant), but the *order* in which targets are
processed determines which branches grow first.  A conjugate-pair
target `(z, z̄)` processed at different positions in the shuffled
list sees a different tree state — different "nearest visited" node —
and the walk diverges from there.

Stage 2 then cements the asymmetry: each grid cell evaluates the
canonical Padé at its nearest visited node.  Non-conjugate nearest-
node assignments at `(z, z̄)` propagate the tree asymmetry into per-cell
values.

### Fix (user-space, examples/tritronquee_3d.jl)

Schwarz reflection: for real-coefficient ODEs with real ICs, `u(z̄) =
ū(z)` globally.  Walk only strictly upper-half targets (`Im(z) > 0`),
then populate the lower-half grid via `u_grid[i, N+1-j] =
conj(u_grid[i, j])`.

  - **Bit-exact y-symmetric**: `max |Δ|u||` = 0.0.
  - **~2× faster**: ~17 s vs ~35 s at 501×501 (only half the targets).
  - **Inapplicable for non-real ICs / complex-coeff ODEs**: those need a
    different fix; the mirror trick is exact only when the ODE
    preserves complex conjugation.

### Filed: `padetaylor-dtj` (P2)

Bead description: add `enforce_real_axis_symmetry::Bool = false` kwarg
to `path_network_solve`.  When `true`: filter targets to upper-half +
on-axis, walk those, populate lower-half `grid_u`/`grid_up` via
conjugation.  Opt-in (the kwarg is inapplicable for the general
PathNetwork algorithm; only set by callers who know their ODE
satisfies the symmetry).

## Bug 2 — y = 0 row discontinuity

### Symptoms

After Bug 1's fix, a horizontal stripe of false-positive pole-field
flags appears along `y = 0` from `x ≈ -13` leftward (smooth region,
no actual poles).  Diagnostic at 121×121:

| x | \|u(x, 0)\| | \|u(x, +Δy)\| | jump |
|---|---|---|---|
| -10.0 | 1.99 | 1.11 | 0.88 |
| -8.33 | 2.33 | 1.23 | 1.10 |
| -6.33 | 3.28 | 1.69 | 1.59 |
| -1.0 | 0.42 | 0.43 | 0.008 |
| ... | (smooth near origin) | | |

The y = 0 row's `|u|` values oscillate wildly cell-to-cell while
adjacent rows are smooth.  Far from the origin (in the asymptotic
smooth region), the discrepancy is 0.3–1.6 — comparable to the
function's magnitude.

### Mechanism

With odd `N`, the grid has a cell exactly at `y = 0`.  The wedge
walker walks to that cell along (mostly) the real axis (using wedge
angle `0` since the goal direction is `0`).  The walks to adjacent
cells (`y = ±Δy`) use a slightly off-axis direction and pick wedge
angle `+π/8` or `0` based on a `:min_u` comparison that can flip
cell-to-cell due to floating-point noise.  Different wedge choices →
different visited-node sequences → different accumulated path-
dependent error.

This is FW's "low-level ridges in flat areas" artifact (FW2011...md:208
verbatim: "These irregularities are caused by the fact that different
u entries in (3.3) may have been obtained following entirely different
paths through the complex plane, and then been amplified by the
division by h² in (3.3)").  It is documented behaviour of the
path-network in smooth regions — the IVP solver was designed to be
accurate near pole fields, not in flat regions.

### Fix-attempt 1: reconstruct `u(x, 0)` from `Re(u(x, +Δy))`

By Schwarz reflection, `u(x) ∈ ℝ` for real `x`; Taylor-expand off the
real axis: `u(x + iΔy) = u(x) + iΔy·u'(x) + O(Δy²)`.  So
`Re(u(x, +Δy)) ≈ u(x, 0)` to `O(Δy²) ≈ 0.1` at our grid spacing
(`Δy = 0.333`).  This is smaller than the walker's path-dependent
noise, so it cleans up the smooth-region stripe.

**But broke in the pole-bearing region.**  At positive real `x` where
poles concentrate, `|u|` varies wildly cell-to-cell — the `O(Δy²)`
error of `Re(u(x, +Δy))` reconstruction became visible as a
horizontal line of false-pole flags on the y = 0 row for `x > 2`.

### Fix-attempt 2 (shipped): use even N

With **even N**, no grid cell sits at `y = 0`.  The closest rows are
at `y = ±Δy/2`, both walked by the same upper-half-walk and mirror
pattern, both flanked by similar cells.  No special-case row needed.

The path-dependent ridges in the smooth region don't vanish (they
shift to `y ≈ ±Δy/2` instead of `y = 0`), but they are much smaller
in magnitude (the walker for off-axis cells takes consistent off-axis
paths, not the special-direction real-axis path) and visually almost
invisible at 500×500.

## BVP-cure exploration — does not apply to tritronquée geometry

The user asked whether the residual ridge is curable by the
LatticeDispatcher's per-row BVP fill.  Built `examples/tritronquee_bvp_compose.jl`
to test: 0 BVPs triggered.

**Why 0**: the line-190 horizontal-row algorithm requires smooth runs
*flanked by pole-field cells on BOTH sides*.  Pure tritronquée has a
*single connected pole-bearing region* (3 sectors fused into one
C-shape wrapping the right side of the grid) and a *single
contiguous pole-free wedge* opening to the negative-x grid edge.
Every horizontal row sees pole cells on at most ONE side (right),
with the smooth half extending unboundedly leftward.  Nothing to
bridge.

The FW-prescribed cure for *open* smooth regions like this is the
**FW Fig 4.1 algorithm** (FW2011...md:218-222): vertical BVP between
`+20i` and `-20i` with asymptotic BCs from `u ~ -√(-z/6)` at `±∞i`,
then outward pole-field IVP integrators from the BVP's endpoints
into the pole-bearing sector.  That's a different compositional
pattern than the horizontal-row line-190 algorithm.  Filed as bead
`padetaylor-0c3` for follow-up work.

## What changed in this session

  - `examples/tritronquee_3d.jl` — the public deliverable, now using
    upper-half walk + conjugate mirror + even N.  Reproduces the
    canonical FW Fig 3.1 tritronquée structure at 500×500 over
    `[-20, 20]²` in ~17 s.
  - `examples/tritronquee_bvp_compose.jl` — investigative script
    demonstrating the 0-BVPs-triggered finding for pure tritronquée
    geometry.  Kept in-tree as a reference (it's the "we tried" log
    for the BVP-cure claim).
  - `.gitignore` — `examples/*.png` already ignored from earlier.

## Bead state

Filed in this session:
  - `padetaylor-dtj` **P2** PathNetwork: add `enforce_real_axis_symmetry`
    kwarg.

Still open from prior sessions (no change):
  - `padetaylor-bvh`, `padetaylor-grc`, `padetaylor-61j`, `padetaylor-8pi`,
    `padetaylor-0c3`.

## Hard-won lessons (for HANDOFF.md §"Hard-won")

**17. `shuffle(rng, targets)` in PathNetwork breaks conjugate symmetry
even when the underlying ODE preserves it.**  The shuffle was a FW
2011 line 156 suggestion to avoid pathological orderings, but for
conjugate-symmetric problems the shuffle's random order creates an
asymmetric visited tree whose Stage-2 cascade can blow `|u|` apart by
4-5 orders of magnitude at conjugate-pair cells.  Fix: for real-coeff
real-IC problems, walk only upper-half + mirror.  Bead `padetaylor-dtj`
for the library-level kwarg.

**18. Use even N for 2D Cartesian grids centred on a symmetry axis.**
With odd N, the centre cell of each row (or column) sits exactly on
the axis and gets walked along a special direction (e.g., real-axis
walks for `y = 0` cells use wedge angle exactly `0`).  This produces a
path-dependent discontinuity from off-axis cells.  Even N avoids
the special-case row entirely.

**19. FW Fig 4.1 BVP geometry is fundamentally different from
line-190's horizontal-row algorithm.**  FW Fig 4.1's vertical BVP
between `±20i` with asymptotic BCs serves *open* smooth regions
extending to infinity; line-190's horizontal-row BVP fill serves
*interior* smooth runs bounded by pole fields on both sides.  Our
v1 `lattice_dispatch_solve` implements line-190; FW Fig 4.1 is
`padetaylor-0c3`.  Don't conflate the two.

## Pointers

  - [`examples/tritronquee_3d.jl`](../../examples/tritronquee_3d.jl)
    — the user-space workaround script.
  - [`examples/tritronquee_bvp_compose.jl`](../../examples/tritronquee_bvp_compose.jl)
    — the 0-BVPs investigative script.
  - `src/PathNetwork.jl:173` — the `shuffle(rng, ...)` line.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    lines 152-158 (FW's random-shuffle suggestion), 208 ("low-level
    ridges in flat areas"), 218-222 (FW Fig 4.1 setup).
  - [`docs/worklog/013-phase-12-v2-lattice.md`](013-phase-12-v2-lattice.md)
    — LatticeDispatcher v1 (line-190 algorithm).
